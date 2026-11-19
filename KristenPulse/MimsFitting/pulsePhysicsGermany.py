# old was 59m 11s
import scipy.constants as sc
import numpy as np
from scipy.signal import fftconvolve

# Constants
mu_0   = sc.physical_constants["vacuum mag. permeability"][0]
beta_e = sc.physical_constants["Bohr magneton"][0]
beta_n = sc.physical_constants["nuclear magneton"][0]
Planck = sc.physical_constants["Planck constant"][0]




def simulate_nanoparticle_mims_spectrum_vectored( # this one is currently working well 2027-07-20
    freq_axis, R_np_exp, density_obj, rho_bulk, nu_n, g_n,
    I_scale, ell, g_e, tau, T_m, R_np_md, R_max, 
    target_dr = 0.2e-10, dipolar_scale = 1, contact_scale = 1, 
    blindspot_offset = 0, dipolar_np_scale = 1, shape = "Lorentzian", 
    eta = 0, width_scale = 1
):
    """
    Accelerated Mims ENDOR simulation preserving original physics equations 
    exactly via 2D vectorization, fractional linear binning, and 1D FFT convolution.
    """
    # --- 1. SETUP 2D COORDINATE GRIDS ---
    n_shells = int(np.ceil((R_max - R_np_md) / target_dr))
    dr = (R_max - R_np_md) / n_shells
    r_axis = np.linspace(R_np_exp + dr/2, R_np_exp + n_shells*dr - dr/2, n_shells)
    cos_theta_axis = np.linspace(0, 1, 150)

    # Pre-calculate density array across shells to drop empty space upfront
    rho_axis = np.array([density_obj.get_density((r - R_np_exp) + R_np_md, rho_bulk) for r in r_axis])
    valid_shells = (rho_axis > 0) & (r_axis >= (1.78e-9 - R_np_md) + R_np_exp)
    
    if not np.any(valid_shells):
        return np.zeros_like(freq_axis)

    r_axis = r_axis[valid_shells]
    rho_filtered = rho_axis[valid_shells, np.newaxis] 

    # Generate the 2D broadcast grids: shape (N_valid_shells, 150)
    r, cos_theta = np.meshgrid(r_axis, cos_theta_axis, indexing='ij')

    df = freq_axis[1] - freq_axis[0]
    f_min = freq_axis[0]
    stick_spectrum = np.zeros_like(freq_axis) 

    # --- 2. ORIGINAL PHYSICS COUPLINGS (Vectorized) ---
    common_scalar = mu_0 * g_e * beta_e * g_n * beta_n / Planck
    f = (1 + 2*R_np_exp/(3*ell))**-1 # for now, we need to just set this. 
    I0 = I_scale * f * np.sqrt(2*f/ell) # reformulation where f indicates fraction of wavefunction outside the particle
    #I0 = I_scale * f * R_np_exp / np.sqrt((R_np_exp**3/3) + (0.5*ell*f**2 * R_np_exp**2)) # uniform density
    #I0 = I_scale * np.sqrt(2/ell) # electron density all at surface.
    #I0 = I_scale * np.sqrt(1/(2*np.pi*ell*((R_np_exp)**2 + (R_np_exp)*ell + 0.5*ell**2))) # erroneous one from before. 

    #print(I0)

    A_contact = (2/3) * common_scalar * (I0 * np.exp(-(r-R_np_exp) / ell)/(r))**2
    A_contact = contact_scale * A_contact

    A_dipolar = (4 * np.pi)**-1 * common_scalar * (r**-3) * (3 * cos_theta**2 - 1)
    F_dipolar = 1#*(1 - (R_np_exp / r)**2)**-3
    A_dipolar = A_dipolar * ((1 - dipolar_np_scale) + dipolar_np_scale * F_dipolar)
    A_dipolar = dipolar_scale * A_dipolar

    A = A_contact + A_dipolar

    # Mims Filter & Weighting
    intensity = (1.0 - blindspot_offset) * np.sin(2 * np.pi * A * tau)**2 + blindspot_offset
    w = (rho_filtered * 4 * np.pi * r**2 * dr / len(cos_theta_axis)) * intensity

    # Resonance split positions
    mu_plus = nu_n + A / 2
    mu_minus = nu_n - A / 2

    # --- 3. FLATTENED FRACTIONAL LINEAR BINNING ---
    idx_plus_exact = ((mu_plus - f_min) / df).ravel()
    idx_minus_exact = ((mu_minus - f_min) / df).ravel()
    w_flat = w.ravel()

    idx_p0 = np.floor(idx_plus_exact).astype(int)
    idx_m0 = np.floor(idx_minus_exact).astype(int)

    frac_p = idx_plus_exact - idx_p0
    frac_m = idx_minus_exact - idx_m0

    valid_p = (idx_p0 >= 0) & (idx_p0 < len(freq_axis) - 1)
    valid_m = (idx_m0 >= 0) & (idx_m0 < len(freq_axis) - 1)

    # Fast bin accumulation
    stick_spectrum += np.bincount(idx_p0[valid_p], weights=w_flat[valid_p] * (1.0 - frac_p[valid_p]), minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_p0[valid_p] + 1, weights=w_flat[valid_p] * frac_p[valid_p], minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_m0[valid_m], weights=w_flat[valid_m] * (1.0 - frac_m[valid_m]), minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_m0[valid_m] + 1, weights=w_flat[valid_m] * frac_m[valid_m], minlength=len(freq_axis))

    # --- 4. CONVOLUTION KERNEL ---
    f_kernel = freq_axis - np.mean(freq_axis)
    gamma_tm = 1 / (2 * np.pi * T_m) * width_scale 

    if shape[0] == "G":
        sigma = np.sqrt((gamma_tm * 1.0)**2) # no longer need the width stretched + (eta * np.mean(np.abs(A))**0.5)**2)
        kernel = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * (f_kernel / sigma)**2)
    elif shape[0] == "L":
        gamma = width_scale * np.sqrt(gamma_tm**2) # no longer need the width stretched+ (eta * np.mean(np.abs(A))**0.5)**2)
        kernel = (1 / (np.pi * gamma)) * (gamma / (f_kernel**2 + gamma**2))
    else:
        return stick_spectrum 

    # Normalize profile and smooth the spectrum
    kernel /= np.sum(kernel)
    return fftconvolve(stick_spectrum, kernel, mode='same')





import numpy as np
from scipy.signal import fftconvolve


def find_R_max(
    R_np_exp, ell, I_scale, T_m, width_scale, contact_scale, dipolar_scale, 
    mu_0, g_e, beta_e, g_n, beta_n, Planck, R_max_cuttoff
):
    # 1. Target threshold based on linewidth
    gamma_hz = (1.0 / (2.0 * np.pi * T_m)) * width_scale  
    threshold = gamma_hz * R_max_cuttoff  # 5% of homogeneous line width

    # 2. Pre-calculate invariant constants OUTSIDE the loop
    common_scalar = mu_0 * g_e * beta_e * g_n * beta_n / Planck
    f = (1.0 + 2.0 * R_np_exp / (3.0 * ell))**-1 
    I0 = I_scale * f * np.sqrt(2.0 * f / ell)

    # 3. Initialize loop variables
    dr = 1e-10  # 1 Angstrom steps
    r = R_np_exp + dr  # Step slightly off surface to avoid F_dipolar singularities
    shift = gamma_hz
    
    # 4. Walk outward until shift < threshold (with a 1000 nm safety valve)
    while shift > threshold and r < 1e-6:
        
        # Fermi contact at distance r
        A_contact_probe = contact_scale * (2.0/3.0) * common_scalar * (I0 * np.exp(-(r - R_np_exp) / ell) / r)**2
        
        # Dipolar at distance r
        T_factor_probe = dipolar_scale * (4.0 * np.pi)**-1 * common_scalar * (r**-3) #*(1 - (R_np_exp / r)**2)**-3 # last term is for spherical shell

        # F_dipolar = (1.0 - (R_np_exp / r)**2)**-3
        # T_factor_probe *= F_dipolar

        # Max coupling magnitude at theta = 0 (dipolar is 2*T)
        shift = A_contact_probe + 2.0 * T_factor_probe

        r += dr

    # 5. Return R_max with a fallback minimum bound
    if r > R_np_exp + 1e-9:
        return r
    else:
        return R_np_exp + 1e-9


def simulate_nanoparticle_mims_spectrum(
    freq_axis, R_np_exp, density_obj, rho_bulk, nu_n, g_n,
    I_scale, ell, g_e, tau, T_m, R_np_md, R_max, #used to be just R_max  
    target_dr = 0.2e-10, dipolar_scale = 1, contact_scale = 1, 
    blindspot_offset = 0, dipolar_np_scale = 1, shape = "Lorentzian", 
    eta = 0, width_scale = 1
):
    """
    Simulates a Mims ENDOR spectrum in the ABSOLUTE laboratory frequency frame.
    freq_axis must span the absolute frequency window (e.g. 45 MHz to 55 MHz).
    This includes the spectra shifting effects of B.
    """
    R_max_cuttoff = R_max # just for now, pass off the value to something else
    # try redefining R_max... 
    R_max = find_R_max(R_np_exp, ell, I_scale, T_m, width_scale, contact_scale, dipolar_scale, 
    mu_0, g_e, beta_e, g_n, beta_n, Planck, R_max_cuttoff)
    #print(R_max)

    # --- 1. SETUP 2D MESHGRIDS ---
    n_shells = int(np.ceil((R_max - R_np_md) / target_dr))
    dr = (R_max - R_np_md) / n_shells
    r_axis = np.linspace(R_np_exp + dr/2, R_np_exp + n_shells*dr - dr/2, n_shells)
    cos_theta_axis = np.linspace(0, 1, 150)

    # Pre-calculate densities
    rho_axis = np.array([density_obj.get_density((r - R_np_exp) + R_np_md, rho_bulk) for r in r_axis])
    valid_shells = (rho_axis > 0) & (r_axis >= (1.78e-9 - R_np_md) + R_np_exp)
    
    if not np.any(valid_shells):
        return np.zeros_like(freq_axis)

    r_axis = r_axis[valid_shells]
    rho_filtered = rho_axis[valid_shells, np.newaxis] 

    r, cos_theta = np.meshgrid(r_axis, cos_theta_axis, indexing='ij')
    sin_theta = np.sqrt(1.0 - cos_theta**2)

    df = freq_axis[1] - freq_axis[0]
    f_min = freq_axis[0]
    stick_spectrum = np.zeros_like(freq_axis) 

    # --- 2. VECTORIZED COUPLINGS WITH PSEUDO-SECULAR B-TERM ---
    common_scalar = mu_0 * g_e * beta_e * g_n * beta_n / Planck
    f = (1 + 2*R_np_exp/(3*ell))**-1 
    I0 = I_scale * f * np.sqrt(2*f/ell)
    
    # Contact and Dipolar
    A_contact = contact_scale * (2/3) * common_scalar * (I0 * np.exp(-(r-R_np_exp) / ell)/(r))**2
    T_factor  = dipolar_scale * (4 * np.pi)**-1 * common_scalar * (r**-3) #*(1 - (R_np_exp / r)**2)**-3 # last term is for spherical shell

    A = A_contact + T_factor * (3 * cos_theta**2 - 1)
    B = 3 * T_factor * sin_theta * cos_theta

    # --- 3. ABSOLUTE SPIN FREQUENCIES (Breit-Rabi) ---
    # Frequencies are kept in the absolute laboratory frame
    nu_alpha_abs = np.sqrt((nu_n + A / 2.0)**2 + (B / 2.0)**2)
    nu_beta_abs  = np.sqrt((nu_n - A / 2.0)**2 + (B / 2.0)**2)

    mu_plus  = nu_alpha_abs
    mu_minus = nu_beta_abs

    #print(f"mu_plus = {mu_plus}")

    # Effective frequency splitting driving the Mims modulation
    delta_nu_eff = np.abs(nu_alpha_abs - nu_beta_abs)

    # Quantization angles rely on absolute Larmor frequency
    eta_alpha = np.arctan2(B, 2.0 * nu_n + A)
    eta_beta  = np.arctan2(B, -2.0 * nu_n + A)
    p_trans   = np.sin((eta_alpha - eta_beta) / 2.0)**2
    #print(f"eta_alpha = {eta_alpha}")

    # Rigorous Mims Filter using effective frequency splitting & transition probability
    mims_filter = (1.0 - blindspot_offset) * np.sin(np.pi * delta_nu_eff * tau)**2 + blindspot_offset
    
    # Base integration weight per orientation/shell
    intensity = mims_filter * p_trans
    #print(f"intensity = {intensity}")
    w = (rho_filtered * 4 * np.pi * r**2 * dr / len(cos_theta_axis)) * intensity

    # --- 4. FLATTENED 2D LINEAR BINNING ---
    # Map absolute frequencies directly onto freq_axis
    idx_plus_exact  = ((mu_plus - f_min) / df).ravel()
    idx_minus_exact = ((mu_minus - f_min) / df).ravel()
    w_flat = w.ravel()
    #print(f"w_flat = {w_flat}")

    # Process mu_plus and mu_minus concurrently
    for idx_exact in (idx_plus_exact, idx_minus_exact):
        idx0 = np.floor(idx_exact).astype(int)
        frac = idx_exact - idx0

        valid = (idx0 >= 0) & (idx0 < len(freq_axis) - 1)
        
        np.add.at(stick_spectrum, idx0[valid], w_flat[valid] * (1.0 - frac[valid]))
        np.add.at(stick_spectrum, idx0[valid] + 1, w_flat[valid] * frac[valid])

    # --- 5. FFT CONVOLUTION KERNEL ---
    # Center kernel at zero offset relative to grid resolution df
    n_pts = len(freq_axis)
    f_kernel = (np.arange(n_pts) - n_pts // 2) * df
    gamma_tm = 1.0 / (2.0 * np.pi * T_m) * width_scale 

    if shape[0].upper() == "G":
        sigma = gamma_tm * width_scale 
        kernel = np.exp(-0.5 * (f_kernel / sigma)**2)
    elif shape[0].upper() == "L":
        gamma = width_scale * gamma_tm
        kernel = gamma / (f_kernel**2 + gamma**2)
    else:
        return stick_spectrum 

    kernel /= np.sum(kernel)

    #print(f"stick spectrum = {stick_spectrum}")
    #print(f"kernel = {kernel}")

    final_spectrum = fftconvolve(stick_spectrum, kernel, mode='same')
    #print(final_spectrum)

    return final_spectrum











def simulate_nanoparticle_mims_spectrum_wB_old(
    freq_axis, R_np_exp, density_obj, rho_bulk, nu_n, g_n,
    I_scale, ell, g_e, tau, T_m, R_np_md, R_max, 
    target_dr = 0.2e-10, dipolar_scale = 1, contact_scale = 1, 
    blindspot_offset = 0, dipolar_np_scale = 1, shape = "Lorentzian", 
    eta = 0, width_scale = 1
):
    """
    Simulates a Mims ENDOR spectrum including the pseudo-secular B-term.
    Fully vectorized over a 2D meshgrid with fractional linear binning.
    """
    # --- 1. SETUP 2D MESHGRIDS ---
    n_shells = int(np.ceil((R_max - R_np_md) / target_dr))
    dr = (R_max - R_np_md) / n_shells
    r_axis = np.linspace(R_np_exp + dr/2, R_np_exp + n_shells*dr - dr/2, n_shells)
    cos_theta_axis = np.linspace(0, 1, 150)

    # Pre-calculate densities to mask out empty physical regions
    rho_axis = np.array([density_obj.get_density((r - R_np_exp) + R_np_md, rho_bulk) for r in r_axis])
    valid_shells = (rho_axis > 0) & (r_axis >= (1.78e-9 - R_np_md) + R_np_exp)
    
    if not np.any(valid_shells):
        return np.zeros_like(freq_axis)

    r_axis = r_axis[valid_shells]
    rho_filtered = rho_axis[valid_shells, np.newaxis] 

    # Create 2D coordinates: shape (N_valid_shells, 150)
    r, cos_theta = np.meshgrid(r_axis, cos_theta_axis, indexing='ij')
    sin_theta = np.sqrt(1.0 - cos_theta**2)

    df = freq_axis[1] - freq_axis[0]
    f_min = freq_axis[0]
    stick_spectrum = np.zeros_like(freq_axis) 

    # --- 2. VECTORIZED COUPLINGS WITH PSEUDO-SECULAR B-TERM ---
    common_scalar = mu_0 * g_e * beta_e * g_n * beta_n / Planck
    f = (1 + 2*R_np_exp/(3*ell))**-1 # for now, we need to just set this. 
    I0 = I_scale * f * np.sqrt(2*f/ell) # reformulation where f indicates fraction of wavefunction outside the particle
    
    # Isotropic contact coupling
    A_contact = (2/3) * common_scalar * (I0 * np.exp(-(r-R_np_exp) / ell)/(r))**2
    A_contact = contact_scale * A_contact

    # Dipolar interaction amplitude (T) including your shell correction factor
    T_factor = (4 * np.pi)**-1 * common_scalar * (r**-3)
    T_factor = dipolar_scale * T_factor

    # Full secular (A) and pseudo-secular (B) expressions
    A = A_contact + T_factor * (3 * cos_theta**2 - 1)
    B = 3 * T_factor * sin_theta * cos_theta

    # Mims Filter & Base Weighting
    # Note: In the presence of strong state mixing (large B), the blindspot modulation 
    # frequency shifts to the true splitting, but sin^2(pi * A * tau) is maintained here
    # as the nominal pulse filter condition.
    intensity = (1.0 - blindspot_offset) * np.sin(np.pi * A * tau)**2 + blindspot_offset
    w = (rho_filtered * 4 * np.pi * r**2 * dr / len(cos_theta_axis)) * intensity

    # --- 3. OFF-LARMOR SPIN FREQUENCIES (Breit-Rabi) ---
    nu_alpha = np.sqrt((nu_n + A / 2)**2 + (B / 2)**2)
    nu_beta  = np.sqrt((nu_n - A / 2)**2 + (B / 2)**2)

    # In ENDOR, the observed nuclear transitions occur at these matched frequencies
    mu_plus  = nu_alpha
    mu_minus = nu_beta

    # --- 4. FLATTENED 2D LINEAR BINNING ---
    idx_plus_exact = ((mu_plus - f_min) / df).ravel()
    idx_minus_exact = ((mu_minus - f_min) / df).ravel()
    w_flat = w.ravel()

    idx_p0 = np.floor(idx_plus_exact).astype(int)
    idx_m0 = np.floor(idx_minus_exact).astype(int)

    frac_p = idx_plus_exact - idx_p0
    frac_m = idx_minus_exact - idx_m0

    valid_p = (idx_p0 >= 0) & (idx_p0 < len(freq_axis) - 1)
    valid_m = (idx_m0 >= 0) & (idx_m0 < len(freq_axis) - 1)

    # Drop frequencies straight into 1D stick spectrum bins
    stick_spectrum += np.bincount(idx_p0[valid_p], weights=w_flat[valid_p] * (1.0 - frac_p[valid_p]), minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_p0[valid_p] + 1, weights=w_flat[valid_p] * frac_p[valid_p], minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_m0[valid_m], weights=w_flat[valid_m] * (1.0 - frac_m[valid_m]), minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_m0[valid_m] + 1, weights=w_flat[valid_m] * frac_m[valid_m], minlength=len(freq_axis))

    # --- 5. FFT CONVOLUTION KERNEL ---
    f_kernel = freq_axis - np.mean(freq_axis)
    gamma_tm = 1 / (2 * np.pi * T_m) * width_scale 

    if shape[0] == "G":
        sigma = gamma_tm * width_scale 
        kernel = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * (f_kernel / sigma)**2)
    elif shape[0] == "L":
        gamma = width_scale * gamma_tm
        kernel = (1 / (np.pi * gamma)) * (gamma / (f_kernel**2 + gamma**2))
    else:
        return stick_spectrum 

    kernel /= np.sum(kernel)
    return fftconvolve(stick_spectrum, kernel, mode='same')






def simulate_nanoparticle_mims_spectrum_full_vector(
    freq_axis, R_np_exp, density_obj, rho_bulk, nu_n, g_n,
    I_scale, ell, g_e, tau, T_m, R_np_md, R_max, 
    target_dr = 0.2e-10, dipolar_scale = 1, contact_scale = 1, 
    blindspot_offset = 0, dipolar_np_scale = 1, shape = "Lorentzian", 
    eta = 0, width_scale = 1
):
    # --- 1. SETUP 2D MESHGRIDS ---
    n_shells = int(np.ceil((R_max - R_np_md) / target_dr))
    dr = (R_max - R_np_md) / n_shells
    r_axis = np.linspace(R_np_exp + dr/2, R_np_exp + n_shells*dr - dr/2, n_shells)
    cos_theta_axis = np.linspace(0, 1, 150)

    # Pre-calculate densities for all shells at once to find valid regions
    rho_axis = np.array([density_obj.get_density((r - R_np_exp) + R_np_md, rho_bulk) for r in r_axis])
    valid_shells = (rho_axis > 0) & (r_axis >= (1.78e-9 - R_np_md) + R_np_exp)
    
    if not np.any(valid_shells):
        return np.zeros_like(freq_axis)

    # Filter our axes upfront to only compute what matters
    r_axis = r_axis[valid_shells]
    rho_filtered = rho_axis[valid_shells, np.newaxis] # Shape (N_valid_shells, 1)

    # Create 2D grids of shape (N_valid_shells, 150)
    r, cos_theta = np.meshgrid(r_axis, cos_theta_axis, indexing='ij')

    df = freq_axis[1] - freq_axis[0]
    f_min = freq_axis[0]
    stick_spectrum = np.zeros_like(freq_axis) 

    # --- 2. FULLY VECTORIZED PHYSICS (No loop!) ---
    common_scalar = mu_0 * g_e * beta_e * g_n * beta_n / Planck
    I0 = I_scale *R_np_exp / np.sqrt((R_np_exp**3/3) + (0.5*ell*I_scale**2 * R_np_exp**2))

    A_contact = (2/3) * common_scalar * (I0 * np.exp(-(r-R_np_exp) / ell)/(r))**2
    A_contact = contact_scale * A_contact

    A_dipolar = (4 * np.pi)**-1 * common_scalar * (r**-3) * (3 * cos_theta**2 - 1)
    A_dipolar = dipolar_scale * A_dipolar

    A = A_contact + A_dipolar

    intensity = (1.0 - blindspot_offset) * np.sin(2 * np.pi * A * tau)**2 + blindspot_offset
    # rho_filtered broadcasts cleanly across the 150 orientation columns
    w = (rho_filtered * 4 * np.pi * r**2 * dr / len(cos_theta_axis)) * intensity

    mu_plus = nu_n + A / 2
    mu_minus = nu_n - A / 2

    # --- 3. FLATTENED 2D LINEAR BINNING ---
    # Flatten the 2D arrays to 1D vectors so np.bincount can process them instantly
    idx_plus_exact = ((mu_plus - f_min) / df).ravel()
    idx_minus_exact = ((mu_minus - f_min) / df).ravel()
    w_flat = w.ravel()

    idx_p0 = np.floor(idx_plus_exact).astype(int)
    idx_m0 = np.floor(idx_minus_exact).astype(int)

    frac_p = idx_plus_exact - idx_p0
    frac_m = idx_minus_exact - idx_m0

    valid_p = (idx_p0 >= 0) & (idx_p0 < len(freq_axis) - 1)
    valid_m = (idx_m0 >= 0) & (idx_m0 < len(freq_axis) - 1)

    # Accumulate all shells and angles in 4 fast array operations
    stick_spectrum += np.bincount(idx_p0[valid_p], weights=w_flat[valid_p] * (1.0 - frac_p[valid_p]), minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_p0[valid_p] + 1, weights=w_flat[valid_p] * frac_p[valid_p], minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_m0[valid_m], weights=w_flat[valid_m] * (1.0 - frac_m[valid_m]), minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_m0[valid_m] + 1, weights=w_flat[valid_m] * frac_m[valid_m], minlength=len(freq_axis))

    # --- 4. CONVOLUTION KERNEL ---
    f_kernel = freq_axis - np.mean(freq_axis)
    gamma_tm = 1 / (2 * np.pi * T_m) * width_scale 

    if shape[0] == "G":
        sigma = gamma_tm * width_scale 
        kernel = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * (f_kernel / sigma)**2)
    elif shape[0] == "L":
        gamma = width_scale * gamma_tm
        kernel = (1 / (np.pi * gamma)) * (gamma / (f_kernel**2 + gamma**2))
    else:
        return stick_spectrum 

    kernel /= np.sum(kernel)
    return fftconvolve(stick_spectrum, kernel, mode='same')












def simulate_nanoparticle_mims_spectrumGermany(
    freq_axis, 
    R_np_exp, 
    density_obj, 
    rho_bulk, 
    nu_n, 
    g_n,
    f, # adjusts the offset of the wavefunction at the surface of the nanoparticle
    ell, 
    g_e, 
    tau, 
    T_m, 
    R_np_md, 
    R_max, 
    target_dr = 0.2e-10, 
    dipolar_scale = 1, 
    contact_scale = 1, 
    blindspot_offset = 0, 
    dipolar_np_scale = 1, 
    shape = "Lorentzian", 
    eta = 0, 
    width_scale = 1
):
    """
    Simulates a Mims ENDOR spectrum using optimized NumPy broadcasting 
    and fast FFT convolution for line broadening.
    """
    # --- 1. SETUP GRIDS ---
    n_shells = int(np.ceil((R_max - R_np_md) / target_dr))
    dr = (R_max - R_np_md) / n_shells
    
    # 1D coordinate arrays
    r_axis = np.linspace(R_np_exp + dr/2, R_np_exp + n_shells*dr - dr/2, n_shells)
    cos_theta_axis = np.linspace(0, 1, 150)
    
    # Create 2D Meshgrids: Shape (n_shells, n_orientations)
    r_grid, cos_theta = np.meshgrid(r_axis, cos_theta_axis, indexing='ij')
    
    # --- 2. VECTORIZED PHYSICS OVER THE 2D MESH ---
    # Filter valid shells upfront based on density and physical limits
    # We evaluate density along the r_axis (1D) first to mask out empty shells
    rho_axis = np.array([density_obj.get_density((r - R_np_exp) + R_np_md, rho_bulk) for r in r_axis])
    valid_mask = (rho_axis > 0) & (r_axis >= (1.78e-9 - R_np_md) + R_np_exp)
    
    if not np.any(valid_mask):
        return np.zeros_like(freq_axis)
        
    # Slice grids to only calculate where physical density exists
    r = r_grid[valid_mask, :]
    cos_theta = cos_theta[valid_mask, :]
    rho = rho_axis[valid_mask, np.newaxis] # broadcast to 2D
    
    common_scalar = mu_0 * g_e * beta_e * g_n * beta_n / Planck
    I0 = f * R_np_exp / np.sqrt((R_np_exp**3/3) + (0.5*ell*f**2 * R_np_exp**2))


    
    # Coupling Calculations
    A_contact = (2/3) * common_scalar * (I0 * np.exp(-r / ell)/(r+R_np_exp))**2
    A_contact = contact_scale * A_contact

    print(np.max(A_contact))

    
    A_dipolar_point = (4 * np.pi)**-1 * common_scalar * (r**-3) * (3 * cos_theta**2 - 1)
    F_dipolar = (1 - (R_np_exp / r)**2)**-3
    A_dipolar = A_dipolar_point * ((1 - dipolar_np_scale) + dipolar_np_scale * F_dipolar)
    A_dipolar = dipolar_scale * A_dipolar
    
    A = A_contact + A_dipolar
    
    # Intensity & Weighting
    intensity = (1.0 - blindspot_offset) * np.sin(2 * np.pi * A * tau)**2 + blindspot_offset
    # Weight per bin: normalized by total orientation count
    w = (rho * 4 * np.pi * r**2 * dr / len(cos_theta_axis)) * intensity
    
    # Resonance Positions
    mu_plus = nu_n + A / 2
    mu_minus = nu_n - A / 2
    
    # --- 3. BUILD THE STICK SPECTRUM (DIGITIZATION) ---
    df = freq_axis[1] - freq_axis[0]
    f_min = freq_axis[0]
    f_max = freq_axis[-1]
    
    stick_spectrum = np.zeros_like(freq_axis)
    
    # Convert physical frequencies to grid indices safely
    idx_plus = np.round((mu_plus - f_min) / df).astype(int)
    idx_minus = np.round((mu_minus - f_min) / df).astype(int)
    
    # Filter indices that fall within the bounds of our freq_axis
    valid_plus = (idx_plus >= 0) & (idx_plus < len(freq_axis))
    valid_minus = (idx_minus >= 0) & (idx_minus < len(freq_axis))
    
    # Populate stick spectrum using np.bincount for rapid accumulation
    stick_spectrum += np.bincount(idx_plus[valid_plus], weights=w[valid_plus], minlength=len(freq_axis))
    stick_spectrum += np.bincount(idx_minus[valid_minus], weights=w[valid_minus], minlength=len(freq_axis))
    
    # --- 4. CONVOLUTION KERNEL ---
    # Center the kernel frequency domain around 0
    f_kernel = freq_axis - np.mean(freq_axis)
    gamma_tm = 1 / (2 * np.pi * T_m) * width_scale
    
    # Note: If strain broadening (eta) is dominant and varies wildly across positions, 
    # an average/representative width is used here.
    mean_A = np.mean(np.abs(A))
    
    if shape[0] == "G":
        sigma = np.sqrt((gamma_tm * 1.0)**2 + (eta * np.mean(np.abs(A))**0.5)**2)
        kernel = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * (f_kernel / sigma)**2)
    elif shape[0] == "L":
        gamma = width_scale * np.sqrt(gamma_tm**2 + (eta * mean_A**0.5)**2)
        kernel = (1 / (np.pi * gamma)) * (gamma / (f_kernel**2 + gamma**2))
    else:
        return stick_spectrum # fallback to sticks if unknown
        
    # Normalize kernel to preserve total spectrum area
    kernel /= np.sum(kernel)
    
    # --- 5. PERFORM FFT CONVOLUTION ---
    total_spectrum = fftconvolve(stick_spectrum, kernel, mode='same')
    
    return total_spectrum