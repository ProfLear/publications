import scipy.constants.physical_constants as pc
import numpy as np

# constants, all in SI units
mu_0   = pc["vacuum mag. permeability"] # J T^-1
beta_e = pc["Bohr magneton"]            # J T^-1
beta_n = pc["nuclear magneton"]         # J T^-1
Planck = pc["Planck constant"]          # J s

def simulate_nanoparticle_mims_spectrum( # this function will return the y-values of a mims spectrum
    freq_axis, # the x-axis in the Mims spectrum

    #properties of the nanoparticle
    R_np, # m; radius of the nanoparticle.

    # properties of the nuclie
    density_obj, # this is an object with attributes related to density of hydrogens over distance
    rho_bulk, # the bulk density for the hydrogens
    nu_n, # Hz,  central frequency of the nuclear frequecy---the frequency about which things split
    g_n, # unitless, nuclear g-factor

    #properties of the electron
    I0, # m^-3/2, Amplitude of the wavefunction at the surface of the nanoparticle. 
    ell, # m, the characteristic spatial distance of the wavefunction
    g_e, # unitless, the isotropic g-value for the metallic electrons

    # properties of the EPR measurement
    tau, # s, the delay time of the Mims pulse sequence <- used to calcluate blind splot
    T_m, # s, the dephasing time for the metallic electron
    sigma, # Hz,  the intrinsic width of the individual features
    eta, # unitless; the amount of broadening strain. Tied to movement of nuclie during experiment

    #properties of the simulation
    R_max, # m; how far out to look for hydrogens when simulating
    target_dr = 0.2 # nm; the target thickness of the various shells
    ):

    """
    Simulates a Mims ENDOR spectrum where the shell thickness (dr) adjusts 
    dynamically to R_max, and density is pulled from a DensityProfile object.

    """

    # --- 1. SHELL DEFINITION ---
    # Target dr defaults to ~0.2 A. We calculate an exact n_shells to keep dr fluid.
    n_shells = int(np.ceil((R_max - R_np) / target_dr)) # scalar unitless.

    # This 'dr' changes smoothly as R_max moves, preventing 'jumpy' fits.
    dr = (R_max - R_np) / n_shells # scalar, nm, 

    # Midpoints of the dynamic shells
    r_full = np.linspace(R_np + dr/2, R_max - dr/2, n_shells) # array, nm

    # --- 2. SETUP GRIDS ---
    df = freq_axis[1] - freq_axis[0] # array, Hz
    cos_theta = np.linspace(0, 1, 150) # array, radians
    total_spectrum = np.zeros_like(freq_axis) # array, units of intensity
    damping = np.exp(-2 * tau / T_m) # scalar, unitless

    # Pre-vectorize the frequency axis for the Gaussian calculation
    f_vec = freq_axis[np.newaxis, :] # convert the freq_axis from an array to a row vector

    # --- 3. INTEGRATION LOOP (Radial) ---
    for r in r_full:

        # Pull density from the object (handles MD, Blending, and Toluene Bulk)
        rho = density_obj.get_density(r, rho_bulk) # scalar, m^-3

        if rho <= 0: continue # if there is no density of nuclei, then we have no need to calculate anything

        # Physics Calculation
        #A_contact = A0 * np.exp(-(max(0, r - R_np)) / ell)
        A_contact = (2*mu_0*g_e*beta_e*g_nucleus*beta_n / (3*Planck))*(I0*exp(-(r - R_np)/ell))**2 # scalar, Hz; coupling constant
        
        #A_dipolar = (79.0 * (g_e / 2.0023)) / (r**3) 
        A_dipolar = (mu_0*g_e*beta_e*g_nucleus*beta_n)/(4*np.pi*Planck*r**3)*(3 * cos_theta**2 - 1) # array, Hz; coupling constant
        
        #A = A_contact + A_dipolar * (3 * cos_theta**2 - 1) # Array coupling constants for 150 angles
        A = A_contact + A_dipolar # array, Hz

        # Mims Filter & Base Weighting
        intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping # scalar, unitless
        w = (rho * 4 * np.pi * r**2 * dr / len(cos_theta)) * intensity # scalar, unitless. accounts for number of protons in the shell

        # --- STRAIN BROADENING (Contact-Only) ---
        # shell_sigma is a single value for this shell since A_contact is isotropic
        shell_sigma = np.sqrt(sigma**2 + (eta * np.abs(A_contact)**0.5)**2)  # scalar, Hz the dependence on A_iso can be changed to reflect different disorders

        # Vectorized Gaussian Calculation for all 150 orientations at once
        A_vec = A[:, np.newaxis] # convert the array to a column vector
        w_vec = w[:, np.newaxis] # convert the array to a column vector
        mu_plus = nu_n + A_vec / 2
        mu_minus = nu_n - A_vec / 2

        # Area normalization: keeps intensity stable as the peak width changes
        norm = shell_sigma * np.sqrt(2 * np.pi)

        g_plus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2) # use column and row vectors to Broadcast
        g_minus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2) # use column and row vectors to Broadcast

        # Sum orientations into the master spectrum
        total_spectrum += np.sum(g_plus + g_minus, axis=0)

    return total_spectrum