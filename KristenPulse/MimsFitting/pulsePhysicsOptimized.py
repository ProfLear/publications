import scipy.constants as sc
import numpy as np

# constants, all in SI units
mu_0   = sc.physical_constants["vacuum mag. permeability"][0] # N A^-2
beta_e = sc.physical_constants["Bohr magneton"][0]            # J T^-1
beta_n = sc.physical_constants["nuclear magneton"][0]         # J T^-1
Planck = sc.physical_constants["Planck constant"][0]          # J s

def simulate_nanoparticle_mims_spectrum_vectorized(
    freq_axis, 
    R_np_exp, 
    density_obj, 
    rho_bulk, 
    nu_n, 
    g_n, 
    I_scale, 
    ell, 
    g_e, 
    tau, 
    T_m, 
    R_np_md, 
    R_max, 
    target_dr=0.2e-10, 
    dipolar_scale=1, 
    contact_scale=1, 
    blindspot_offset=0, 
    dipolar_np_scale=1, 
    shape="Lorentzian", 
    eta=0, 
    width_scale=1
):
    """
    Fully vectorized version of the Mims ENDOR simulation.
    Eliminates the radial shell loop using NumPy broadcasting.
    """
    # 1. Shell Definition
    n_shells = int(np.ceil((R_max - R_np_md) / target_dr))
    if n_shells <= 0:
        return np.zeros_like(freq_axis)
        
    dr = (R_max - R_np_md) / n_shells
    # Shape (n_shells,)
    r_full = np.linspace(R_np_exp + dr/2, R_np_exp + n_shells*dr - dr/2, n_shells)

    # 2. Setup Grids
    # Shape (n_orientations,)
    cos_theta = np.linspace(0, 1, 150)
    
    # 3. Physics Constants
    common_scalar = mu_0 * g_e * beta_e * g_n * beta_n / Planck
    I0 = np.sqrt(1 / (2 * np.pi * ell * (R_np_exp**2 + R_np_exp * ell + 0.5 * ell**2)))
    
    # 4. Get Density for all shells at once
    # Vectorized get_density (assuming density_obj.get_density can handle arrays or using np.interp)
    # Most MD density lookups use np.interp which is vectorized.
    # We need to ensure mdHelper's get_density is compatible or call it via np.vectorize
    r_md_coords = (r_full - R_np_exp) + R_np_md
    
    # Check if get_density is vectorized. If not, we'll use np.interp directly if possible.
    # For now, let's assume we can get rhos as an array of shape (n_shells,)
    rhos = np.array([density_obj.get_density(r, rho_bulk) for r in r_md_coords])
    
    # Mask for valid density and distance
    mask = (rhos > 0) & (r_full >= (1.78e-9 - R_np_md) + R_np_exp)
    if not np.any(mask):
        return np.zeros_like(freq_axis)
        
    r_v = r_full[mask][:, np.newaxis] # Shape (N, 1)
    rho_v = rhos[mask][:, np.newaxis] # Shape (N, 1)
    ct_v = cos_theta[np.newaxis, :]   # Shape (1, 150)

    # 5. Coupling Calculations (Broadcasting shells and orientations)
    # A_contact: Shape (N, 1)
    A_contact = contact_scale * (2/3 * common_scalar * (I_scale * I0 * np.exp(-r_v / ell))**2)
    
    # A_dipolar: Shape (N, 150)
    A_dipolar_pt = (4*np.pi)**-1 * common_scalar * (r_v**-3) * (3 * ct_v**2 - 1)
    F_dipolar = (1 - (R_np_exp / r_v)**2)**-3
    A_dipolar = A_dipolar_pt * ((1 - dipolar_np_scale) + dipolar_np_scale * F_dipolar)
    A_dipolar *= dipolar_scale
    
    # Total Coupling A: Shape (N, 150)
    A = A_contact + A_dipolar
    
    # 6. Mims Filter & Weighting
    # intensity: Shape (N, 150)
    intensity = (1.0 - blindspot_offset) * np.sin(2 * np.pi * A * tau)**2 + blindspot_offset
    # w: Shape (N, 150)
    w = (rho_v * 4 * np.pi * r_v**2 * dr / len(cos_theta)) * intensity
    
    # 7. Resonance Conditions
    # Shape (N, 150, 1)
    mu_plus = nu_n + A[:, :, np.newaxis] / 2
    mu_minus = nu_n - A[:, :, np.newaxis] / 2
    f_vec = freq_axis[np.newaxis, np.newaxis, :] # Shape (1, 1, n_freqs)
    
    # 8. Line Shapes (Fully Vectorized)
    gamma_tm = 1 / (2 * np.pi * T_m) * width_scale
    
    if shape[0] == "G":
        # shell_sigma: Shape (N, 150, 1)
        shell_sigma = np.sqrt(gamma_tm**2 + (eta * np.abs(A)**0.5)**2)[:, :, np.newaxis]
        norm = shell_sigma * np.sqrt(2 * np.pi)
        g_plus = (w[:, :, np.newaxis] / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2)
        g_minus = (w[:, :, np.newaxis] / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2)
    else: # Lorentzian
        # shell_gamma: Shape (N, 150, 1)
        shell_gamma = np.sqrt(gamma_tm**2 + (eta * np.abs(A)**0.5)**2)[:, :, np.newaxis]
        norm = np.pi * shell_gamma
        g_plus = (w[:, :, np.newaxis] / norm) * (shell_gamma / ((f_vec - mu_plus)**2 + shell_gamma**2))
        g_minus = (w[:, :, np.newaxis] / norm) * (shell_gamma / ((f_vec - mu_minus)**2 + shell_gamma**2))

    # 9. Final Summation
    # Sum over N (shells) and 150 (orientations)
    return np.sum(g_plus + g_minus, axis=(0, 1))
