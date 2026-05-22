import scipy.constants as sc
import numpy as np

# constants, all in SI units
mu_0   = sc.physical_constants["vacuum mag. permeability"][0] # N A^-2
beta_e = sc.physical_constants["Bohr magneton"][0]            # J T^-1
beta_n = sc.physical_constants["nuclear magneton"][0]         # J T^-1
Planck = sc.physical_constants["Planck constant"][0]          # J s


def simulate_nanoparticle_mims_spectrum( # this function will return the y-values of a mims spectrum
    freq_axis, # the x-axis in the Mims spectrum

    #properties of the experimental nanoparticle
    R_np_exp, # m; radius of the nanoparticle.

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

    #properties of the spectral simulation
    R_np_md, # m, radius of the nanoparticle in the simulation
    R_max, # m; how far out to look for hydrogens when simulating
    target_dr = 0.2e-10, # m; the target thickness of the various shells
    dipolar_scale = 1, # float, scale for the relative importance of dipolar coupling.  Set to 0 to turn off coupling mechanism
    contact_scale = 1, # float, scale for the relative importance of contact coupling.  Set to 0 to turn off coupling mechanism
    blindspot_offset = 0, # float, between 0 and 1, accounts for imperfect blindspot, so ajusts intensity, rasiing the minima from zero
    dipolar_np_scale = 1, # adjusts the distance from the nanoparticle surface that dipolar coupling occurs at

    # properties of the spectrum
    shape = "Lorentzian", # handles what line shape to use, will be either Gaussian, Lorentzian, or Voight
    eta = 0, # unitless; the amount of broadening strain. Tied to movement of nuclie during experiment
    width_scale = 1, # unitless,  scales the intrinsic width of the individual features. Gaussian, we set intrinsic as 1 Hz
    ):

    """
    Simulates a Mims ENDOR spectrum where the shell thickness (dr) adjusts 
    dynamically to R_max, and density is pulled from a DensityProfile object.

    """

    # correct units...
    #ell = ell*1e-9
    #R_max= R_max*1e-9


    # print(freq_axis)
    # --- 1. SHELL DEFINITION ---
    # Target dr defaults to ~0.2 A. We calculate an exact n_shells to keep dr fluid.
    n_shells = int(np.ceil((R_max - R_np_md) / target_dr)) # scalar unitless.
    #print(R_max)
    #print(R_np_md)
    #print(target_dr)
    #print(f"n_shells = {n_shells}")


    # This 'dr' changes smoothly as R_max moves, preventing 'jumpy' fits.
    dr = (R_max - R_np_md) / n_shells # scalar, nm, 

    # Midpoints of the dynamic shells
    r_full = np.linspace(R_np_exp + dr/2, R_np_exp + n_shells*dr - dr/2, n_shells) # array, nm

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
    I0 = np.sqrt(1/(2*np.pi*ell*((R_np_exp)**2 + (R_np_exp)*ell + 0.5*ell**2)))
    #print(f"I0 = {I0}")

    #print(common_scalar)
    p = True # here is just to control when printing happens
    for i, r in enumerate(r_full):  

        #print(density_obj)      
        # Pull density from the object (handles MD, Blending, and Toluene Bulk)
        rho = density_obj.get_density((r - R_np_exp) + R_np_md, rho_bulk) # scalar, m^-3, find the distance from teh NP core, adn then look that far from the simulated one
        #print(rho)

        if rho <= 0 or r < (1.78e-9 - R_np_md) + R_np_exp: continue # if there is no density of nuclei, then we have no need to calculate anything. Also, the second one handles the fact there is a spurious solvent very close to the core in the MD
        
        # Physics Calculation
        #A_contact = A0 * np.exp(-(max(0, r - R_np)) / ell)
        A_contact = 2/3*common_scalar*(I_scale*I0*np.exp(-(r)/ell))**2# scalar, Hz; coupling constant
        A_contact = contact_scale * A_contact

        
        # point diploar with an adjusted shell for the radius
        #A_dipolar = (4*np.pi)**-1*common_scalar*((r-R_np_exp*dipolar_np_scale)**-3)*(3 * cos_theta**2 - 1) # 1.05 array, Hz; coupling constant
        
        # for shell stuff
        A_dipolar = (4*np.pi)**-1*common_scalar*(r**-3)*(3 * cos_theta**2 - 1) # array, Hz; coupling constant. THe point-source equation
        F_dipolar = (1-(R_np_exp/r)**2)**-3
        A_dipolar = A_dipolar*((1-dipolar_np_scale) + dipolar_np_scale*F_dipolar) # multiply the point source by the correction. 

        A_dipolar = dipolar_scale * A_dipolar

        
        #if p and rho > 0 and i < 20:
            #print(f"at distance {r}, rho is {rho}")
            #print(f"contact constant = {A_contact}")
            #print(f"max dipolar constant = {np.max(A_dipolar)}")
            #print("")
            #p = False



        #A = A_contact + A_dipolar * (3 * cos_theta**2 - 1) # Array coupling constants for 150 angles
        A = A_contact + A_dipolar # array, Hz
        # print(f"coupling constant = {A}")

        # Mims Filter & Base Weighting
        intensity = (1.0-blindspot_offset)*np.sin(2*np.pi*A*tau)**2 + blindspot_offset # if using this, do not need the damping calc
        #intensity = ((1 - blindspot_offset) - np.cos(2 * np.pi * A * tau)) * damping + blindspot_offset # array.  Calculate the blindspot, based on distance from resonance (A).
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

            # Conversion from time (s) to frequency intrinsic width (Hz)
            gamma_tm = 1 / (2 * np.pi * T_m) * width_scale # base the width off of the lifetime... just to start

            # --- ANISOTROPIC STRAIN BROADENING (Total Coupling) ---
            # A is a 1D array of length 150. np.abs(A) gives the magnitude at each orientation.
            # We add [:, np.newaxis] to turn it into a (150, 1) column vector.
            shell_sigma = np.sqrt((gamma_tm * 1.0)**2 + (eta * np.abs(A)**0.5)**2)[:, np.newaxis]

            # Area normalization factor is now a (150, 1) column vector
            norm = shell_sigma * np.sqrt(2 * np.pi)

            # Broadcasting math: (150, 1) / (150, 1) * exp( ((1, 500) - (150, 1)) / (150, 1) )
            g_plus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2) 
            g_minus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2) 
        
        elif shape[0] == "L":
            # LORENTZIAN
            # Conversion from time (s) to frequency intrinsic width (Hz)
            gamma_tm = 1 / (2 * np.pi * T_m) * width_scale

            # --- ANISOTROPIC STRAIN BROADENING (Total Coupling) ---
            # We calculate a unique HWHM for each of the 150 orientations and reshape to (150, 1)
            shell_gamma = (width_scale * np.sqrt(gamma_tm**2 + (eta * np.abs(A)**0.5)**2))[:, np.newaxis]

            # Normalization factor is now a (150, 1) column vector
            norm = np.pi * shell_gamma

            # Broadcasting math for Lorentzian profile matrix (150 x 500)
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













#########
# NON SI UNITS!
#########


import scipy.constants as sc
import numpy as np

# constants, all in SI units
mu_0   = sc.physical_constants["vacuum mag. permeability"][0] # N A^-2
beta_e = sc.physical_constants["Bohr magneton"][0]            # J T^-1
beta_n = sc.physical_constants["nuclear magneton"][0]         # J T^-1
Planck = sc.physical_constants["Planck constant"][0]          # J s

def non_SI_simulate_nanoparticle_mims_spectrum( # this function will return the y-values of a mims spectrum
    freq_axis, # Hz --> MHz, the x-axis in the Mims spectrum

    #properties of the nanoparticle
    R_np, # m --> nm; radius of the nanoparticle.

    # properties of the nuclie
    density_obj, # this is an object with attributes related to density of hydrogens over distance
    rho_bulk, # m^-3 --> nm^-3 the bulk density for the hydrogens
    nu_n, # Hz --> MHz,  central frequency of the nuclear frequecy---the frequency about which things split
    g_n, # unitless, nuclear g-factor

    #properties of the electron
    I_scale, # unitless, A scalar to apply to the amplitude of the wavefunction at the surface.  I guess just some wiggle thing here until I understand what is happening
    ell, # m --> nm, the characteristic spatial distance of the wavefunction
    g_e, # unitless, the isotropic g-value for the metallic electrons

    # properties of the EPR measurement
    tau, # s --> micro s, the delay time of the Mims pulse sequence <- used to calcluate blind splot
    T_m, # s --> micro s, the dephasing time for the metallic electron

    #properties of the simulation
    R_max, # m --> nm; how far out to look for hydrogens when simulating
    target_dr = 0.2e-10, # m --> nm; the target thickness of the various shells
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


    # convert incomming units... 
    freq_axis = freq_axis*1e-6 # convert Hz --> MHz
    R_np = R_np*1e9 # convert m --> nm
    #rho_bulk = rho_bulk * 1e-27  # convert m^3 to nm^3 handled in code, when needed
    nu_n = nu_n *1e-6 # convert Hz --> MHz
    ell = ell* 1e9 # convert m to nm
    tau = tau*1e6 # covert s to microseconds
    T_m = T_m*1e6 # convert s --> microseconds
    R_max = R_max *1e9 # convert m to nm
    target_dr = target_dr*1e9 # convert m --> nm





    # print(freq_axis)
    # --- 1. SHELL DEFINITION ---
    # Target dr defaults to ~0.2 A. We calculate an exact n_shells to keep dr fluid.
    n_shells = int(np.ceil((R_max - R_np) / target_dr)) # scalar unitless.

    # This 'dr' changes smoothly as R_max moves, preventing 'jumpy' fits.
    dr = (R_max - R_np) / n_shells # scalar, unit = nm, 

    # Midpoints of the dynamic shells
    r_full = np.linspace(R_np + dr/2, R_max - dr/2, n_shells) # array, nm

    # --- 2. SETUP GRIDS ---
    df = freq_axis[1] - freq_axis[0] # array, MHz
    cos_theta = np.linspace(0, 1, 150) # array, radians
    total_spectrum = np.zeros_like(freq_axis) # array, units of intensity
    damping = np.exp(-2 * tau / T_m) # scalar, unitless

    # Pre-vectorize the frequency axis for the Gaussian calculation
    f_vec = freq_axis[np.newaxis, :] # convert the freq_axis from an array to a row vector



    # --- 3. INTEGRATION LOOP (Radial) ---

    #start by doing a few calculations of the things that do not depend on r

    common_scalar = mu_0*g_e*beta_e*g_n*beta_n/Planck # in units of m^3 Hz... 
    common_scalar = common_scalar*1e27*1e-6# convert to.. nm^3 MHz
    
    I0 = np.sqrt(1/(2*np.pi*ell*(R_np**2 + R_np*ell + 0.5*ell**2))) # in units of... nm^-1.5. This is already correct, since we scalled things already
    #print(f"I0 = {I0}")

    #print(common_scalar)
    p = True # here is just to control when printing happens
    for i, r in enumerate(r_full):  

        #print(density_obj)      
        # Pull density from the object (handles MD, Blending, and Toluene Bulk)
        rho = density_obj.get_density((r-R_np)*1e-9+1.5e-9, rho_bulk)*1e-27 # convert from m^-3 to nm^-3 # scalar, m^-3, find the distance from teh NP core, adn then look that far from the simulated one
        #print(rho)

        if rho <= 0 or r < 1.78: continue # if there is no density of nuclei, then we have no need to calculate anything
        
        # Physics Calculation
        #A_contact = A0 * np.exp(-(max(0, r - R_np)) / ell)
        A_contact = 2/3*common_scalar*(I_scale*I0*np.exp(-(r - R_np)/ell))**2# scalar, Hz; coupling constant
        A_contact = contact_scale * A_contact

        
        #A_dipolar = (79.0 * (g_e / 2.0023)) / (r**3) 
        A_dipolar = (4*np.pi)**-1*common_scalar*(r**-3)*(3 * cos_theta**2 - 1) # array, Hz; coupling constant
        A_dipolar = dipolar_scale * A_dipolar

        
        if p and rho > 0 and i < 20:
            print(f"at distance {r}, rho is {rho}")
            print(f"I0 is {I0}")
            print(f"contact constant = {A_contact}")
            print(f"max dipolar constant = {np.max(A_dipolar)}")
            print("")
            #p = False



        #A = A_contact + A_dipolar * (3 * cos_theta**2 - 1) # Array coupling constants for 150 angles
        A = A_contact + A_dipolar # array, Hz
        # print(f"coupling constant = {A}")

        # Mims Filter & Base Weighting
        intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping # array.  Calculate the blindspot, based on distance from resonance (A).
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