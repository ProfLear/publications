import torch
import numpy as np
import scipy.constants as sc

# Physics Constants as Tensors
mu_0   = sc.physical_constants["vacuum mag. permeability"][0]
beta_e = sc.physical_constants["Bohr magneton"][0]
beta_n = sc.physical_constants["nuclear magneton"][0]
Planck = sc.physical_constants["Planck constant"][0]

def simulate_mims_gpu(
    freq_axis, 
    rs, weights, # Nanoparticle size distribution
    density_obj, rho_bulk, 
    nu_n, g_n, I_scale, ell, g_e, tau, T_m, 
    R_np_md, R_max, target_dr=0.2e-10, 
    dipolar_scale=1.0, contact_scale=1.0, 
    blindspot_offset=0.0, dipolar_np_scale=1.0,
    shape="Lorentzian", eta=0.0, width_scale=1.0,
    device="cuda" if torch.cuda.is_available() else "cpu"
):
    """
    Super-vectorized PyTorch implementation.
    Calculates the ENTIRE size distribution spectrum in one pass.
    """
    # 1. Setup Tensors
    f_axis = torch.as_tensor(freq_axis, device=device, dtype=torch.float32)
    radii = torch.as_tensor(rs, device=device, dtype=torch.float32)      # Shape (S,)
    w_size = torch.as_tensor(weights, device=device, dtype=torch.float32) # Shape (S,)
    
    n_shells = int(np.ceil((R_max - R_np_md) / target_dr))
    dr = (R_max - R_np_md) / n_shells
    
    # Create radial grid: Shape (S, N_shells)
    # Each NP size has its own radial shell midpoints
    shell_idx = torch.arange(n_shells, device=device).float() + 0.5
    r_grid = radii[:, None] + (shell_idx * dr) # Shape (S, N)
    
    # 2. Density Lookup (Vectorized)
    # Assuming density_obj can be represented as an interpolation table on GPU
    # For now, we pre-calculate on CPU or use a simple torch.interp equivalent
    r_md = (r_grid - radii[:, None]) + R_np_md
    # Simple linear interpolation in torch
    xp = torch.tensor(density_obj.bin_centers, device=device, dtype=torch.float32)
    fp = torch.tensor(density_obj.densities, device=device, dtype=torch.float32)
    
    # Custom 1D interpolation for 2D grid
    def torch_interp(x, xp, fp):
        # This is a simplified version; in production use a more robust interp
        idx = torch.bucketize(x, xp) - 1
        idx = torch.clamp(idx, 0, len(xp) - 2)
        x0, x1 = xp[idx], xp[idx+1]
        y0, y1 = fp[idx], fp[idx+1]
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

    rhos = torch_interp(r_md, xp, fp)
    rhos = torch.where(r_md > xp[-1], torch.tensor(rho_bulk, device=device), rhos)
    
    # 3. Physics
    common = mu_0 * g_e * beta_e * g_n * beta_n / Planck
    # I0 depends on the radius of each nanoparticle size: Shape (S, 1)
    I0 = torch.sqrt(1 / (2 * np.pi * ell * (radii**2 + radii * ell + 0.5 * ell**2)))[:, None]
    
    # Orientations: Shape (O,)
    cos_theta = torch.linspace(0, 1, 150, device=device)
    
    # A_contact: Shape (S, N)
    A_contact = contact_scale * (2/3 * common * (I_scale * I0 * torch.exp(-r_grid / ell))**2)
    
    # A_dipolar: Shape (S, N, O)
    r_inv3 = r_grid[..., None]**-3
    angle_term = (3 * cos_theta**2 - 1)
    A_dip_pt = (4*np.pi)**-1 * common * r_inv3 * angle_term
    
    F_dip = (1 - (radii[:, None] / r_grid)**2)**-3
    A_dip = A_dip_pt * ((1 - dipolar_np_scale) + dipolar_np_scale * F_dip[..., None])
    A_dip *= dipolar_scale
    
    # Total A: Shape (S, N, O)
    A = A_contact[..., None] + A_dip
    
    # 4. Mims Filter
    intensity = (1.0 - blindspot_offset) * torch.sin(2 * np.pi * A * tau)**2 + blindspot_offset
    # w: Shape (S, N, O)
    vol_term = (rhos * 4 * np.pi * r_grid**2 * dr)[..., None]
    w = (vol_term / len(cos_theta)) * intensity * w_size[:, None, None]

    # 5. Line Shapes
    mu_plus = nu_n + A / 2
    mu_minus = nu_n - A / 2
    gamma_tm = 1 / (2 * np.pi * T_m) * width_scale
    
    # Shape broadening: (S, N, O)
    shell_gamma = torch.sqrt(gamma_tm**2 + (eta * torch.abs(A)**0.5)**2)
    
    # Final Spectrum Assembly: Shape (Freqs,)
    # We use broadcasting over (S, N, O, Freqs) but carefully to avoid memory blowup
    # Summing iteratively over (S, N, O) is safer
    total_spectrum = torch.zeros_like(f_axis)
    
    # Combine (S, N, O) into a single batch dimension B
    w_flat = w.view(-1, 1)
    mu_p_flat = mu_plus.view(-1, 1)
    mu_m_flat = mu_minus.view(-1, 1)
    gamma_flat = shell_gamma.view(-1, 1)
    
    # Batch processing to fit in GPU memory if needed
    # (S*N*O is roughly 500 * 500 * 150 = 37 million)
    # We can process in chunks
    chunk_size = 100000
    for i in range(0, w_flat.shape[0], chunk_size):
        w_c = w_flat[i:i+chunk_size]
        mu_p_c = mu_p_flat[i:i+chunk_size]
        mu_m_c = mu_m_flat[i:i+chunk_size]
        g_c = gamma_flat[i:i+chunk_size]
        
        if shape[0] == "L":
            g_plus = (w_c / (np.pi * g_c)) * (g_c / ((f_axis - mu_p_c)**2 + g_c**2))
            g_minus = (w_c / (np.pi * g_c)) * (g_c / ((f_axis - mu_m_c)**2 + g_c**2))
        
        total_spectrum += (g_plus + g_minus).sum(dim=0)

    return total_spectrum.cpu().numpy()
