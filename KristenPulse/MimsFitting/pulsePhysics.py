import scipy.constants as sc
import numpy as np

# constants, all in SI units
mu_0   = sc.physical_constants["vacuum mag. permeability"][0] # J T^-1
beta_e = sc.physical_constants["Bohr magneton"][0]            # J T^-1
beta_n = sc.physical_constants["nuclear magneton"][0]         # J T^-1
Planck = sc.physical_constants["Planck constant"][0]          # J s

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
    I_scale, # unitless, A scalar to apply to the amplitude of the wavefunction at the surface.  I guess just some wiggle thing here until I understand what is happening
    ell, # m, the characteristic spatial distance of the wavefunction
    g_e, # unitless, the isotropic g-value for the metallic electrons

    # properties of the EPR measurement
    tau, # s, the delay time of the Mims pulse sequence <- used to calcluate blind splot
    T_m, # s, the dephasing time for the metallic electron

    #properties of the simulation
    R_max, # m; how far out to look for hydrogens when simulating
    target_dr = 0.2e-10, # m; the target thickness of the various shells
    dipolar_scale = 1, # float, scale for the relative importance of dipolar coupling.  Set to 0 to turn off coupling mechanism
    contact_scale = 1, # float, scale for the relative importance of contact coupling.  Set to 0 to turn off coupling mechanism

    # properties of the spectrum
    shape = "Lorentzian", # handles what line shape to use, will be either Gaussian, Lorentzian, or Voight
    eta = 0, # unitless; the amount of broadening strain. Tied to movement of nuclie during experiment
    width_scale = 1, # unitless,  scales the intrinsic width of the individual features. Gaussian, we set intrinsic as 1 Hz
    ):

    """
    Simulates a Mims ENDOR spectrum where the shell thickness (dr) adjusts 
    dynamically to R_max, and density is pulled from a DensityProfile object.

    """

    # print(freq_axis)
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

    #start by doing a few calculations of the things that do not depend on r

    common_scalar = mu_0*g_e*beta_e*g_n*beta_n/Planck
    I0 = np.sqrt(1/(2*np.pi*ell*(R_np**2 + R_np*ell + 0.5*ell**2)))
    #print(f"I0 = {I0}")

    #print(common_scalar)

    for r in r_full:  
        #print(density_obj)      
        # Pull density from the object (handles MD, Blending, and Toluene Bulk)
        rho = density_obj.get_density(r, rho_bulk) # scalar, m^-3
        #print(rho)

        if rho <= 0: continue # if there is no density of nuclei, then we have no need to calculate anything
        
        # Physics Calculation
        #A_contact = A0 * np.exp(-(max(0, r - R_np)) / ell)
        A_contact = 2/3*common_scalar*(I_scale*I0*np.exp(-(r - R_np)/ell))**2 *1e2# scalar, Hz; coupling constant
        A_contact = contact_scale * A_contact
        #print(f"contact constant = {A_contact}")
        
        #A_dipolar = (79.0 * (g_e / 2.0023)) / (r**3) 
        A_dipolar = (4*np.pi)**-1*common_scalar*(r**-3)*(3 * cos_theta**2 - 1) # array, Hz; coupling constant
        A_dipolar = dipolar_scale * A_dipolar
        #print(f"dipolar constant = {A_dipolar}")

        #A = A_contact + A_dipolar * (3 * cos_theta**2 - 1) # Array coupling constants for 150 angles
        A = A_contact + A_dipolar # array, Hz
        # print(f"coupling constant = {A}")

        # Mims Filter & Base Weighting
        intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping # scalar, unitless
        w = (rho * 4 * np.pi * r**2 * dr / len(cos_theta)) * intensity # scalar, unitless. accounts for number of protons in the shell

        # calculate the resonance conditions
        # Vectorized Gaussian Calculation for all 150 orientations at once
        A_vec = A[:, np.newaxis] # convert the array to a column vector
        w_vec = w[:, np.newaxis] # convert the array to a column vector
        mu_plus = nu_n + A_vec / 2
        mu_minus = nu_n - A_vec / 2

        # now, create profiles for all the resonance positions calculated above. 
        #print(shape)
        if shape[0] == "G":
            # GAUSSIAN

            # --- STRAIN BROADENING (Contact-Only) ---
            # shell_sigma is a single value for this shell since A_contact is isotropic
            shell_sigma = np.sqrt((width_scale*1)**2 + (eta * np.abs(A_contact)**0.5)**2)  # scalar, Hz the dependence on A_iso can be changed to reflect different disorders

            # Area normalization: keeps intensity stable as the peak width changes
            norm = shell_sigma * np.sqrt(2 * np.pi)

            g_plus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2) # use column and row vectors to Broadcast
            g_minus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2) # use column and row vectors to Broadcast

        elif shape[0] == "L":
            #print("shape is L")
            # LORENTZIAN
            # --- CALCULATE WIDTH FROM DEPHASING TIME ---
            # Conversion from time (s) to frequency width (Hz)
            gamma_tm = 1 / (2 * np.pi * T_m)*width_scale

            #print(gamma_tm)

            # You can still keep eta for strain-based broadening if you suspect 
            # the width increases with the strength of the contact coupling. <-- true??
            shell_gamma = width_scale*np.sqrt(gamma_tm**2 + (eta * np.abs(A_contact)**0.5)**2)

            norm = np.pi * shell_gamma
            g_plus = (w_vec / norm) * (shell_gamma / ((f_vec - mu_plus)**2 + shell_gamma**2))
            g_minus = (w_vec / norm) * (shell_gamma / ((f_vec - mu_minus)**2 + shell_gamma**2))
        elif shape[0] == "V":
            # Voight
            pass # to be written, if needed

        # Sum orientations into the master spectrum
        total_spectrum += np.sum(g_plus + g_minus, axis=0)

        #print(total_spectrum)

    #print(total_spectrum)

    return total_spectrum