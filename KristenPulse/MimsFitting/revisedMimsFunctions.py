import numpy as np



def simulate_case(md_distances, R_np, params, freq_axis, rho_bulk=0.0453):
    """
    One-stop function to simulate a single ENDOR spectrum.
    """
    # 1. Pre-process MD into Density
    # (Using the logic we developed earlier)
    md_r_centers, md_rho = process_md_to_density(md_distances, R_np)
    
    # 2. Define the Simulation Grid
    # We go from R_np to R_max (calculated from sigma as we discussed)
    R_max = calculate_dynamic_cutoff(params)
    dr = 0.2 # Bin width
    r_grid = np.arange(R_np, R_max + dr, dr)
    
    total_spectrum = np.zeros_like(freq_axis)
    
    # 3. The Integration Loop
    for r in r_grid:
        # Get the correct density (MD or Blended Bulk)
        rho = get_blended_density(r, md_r_centers, md_rho, 
                                  R_stitch=md_r_centers[-1]-2, 
                                  blend_width=5.0, 
                                  rho_bulk=rho_bulk)
        
        n_protons = rho * (4 * np.pi * r**2) * dr
        
        # Calculate the Pake Pattern for this shell
        # This calls our generate_pake_pattern from earlier
        shell_spec = generate_pake_pattern(r, n_protons, params, freq_axis)
        total_spectrum += shell_spec
        
    # 4. Final Broadening
    # Convolve with the Gaussian 'sigma' parameter
    final_spectrum = apply_gaussian_broadening(total_spectrum, params['sigma'])
    
    return final_spectrum




def mims_filter(A, tau):
    """
    Calculates the ideal (undamped) Mims suppression factor.
    
    Parameters:
    A   : ndarray
          Total hyperfine coupling (MHz). This should be A_iso + A_dip.
    tau : float
          The first inter-pulse delay in the Mims sequence (microseconds).
          
    Returns:
    eta : ndarray
          The modulation depth (0 to 2 range).
    """
    # The physical form is 1 - cos(2 * pi * A * tau)
    # A is in MHz, tau is in us, so the product is dimensionless.
    return 1 - np.cos(2 * np.pi * A * tau)

def damp_mims_filter(intensity, tau, T_m):
    """
    Dampens the ENDOR intensity based on the phase memory time.
    
    Parameters:
    intensity : ndarray
                The output from the mims_filter function.
    tau       : float
                The inter-pulse delay (microseconds).
    T_m       : float
                The measured phase memory time (microseconds).
                
    Returns:
    damped_intensity : ndarray
                       The signal scaled by the coherence loss.
    """
    # The signal in a Mims experiment decays as exp(-2 * tau / T_m)
    # This accounts for the loss of electron spin coherence during the tau intervals.
    damping_factor = np.exp(-2 * tau / T_m)
    
    return intensity * damping_factor

def calculate_total_coupling(r, theta, params):
    """
    Calculates the sum of Isotropic (Contact) and Dipolar couplings.
    
    Parameters:
    r      : ndarray
             Distances from the CENTER of the nanoparticle (Angstroms).
    theta  : ndarray or float
             Angle between B0 and the electron-nuclear vector (radians).
    params : dict or LMFit object
             Must contain: 'R_np', 'rho', 'A0', 'g_e', 'g_n'
             
    Returns:
    A_total : ndarray
              Total coupling in MHz.
    """
    # 1. Constants for Dipolar Coupling (T)
    # Prefactor for T in MHz*Angstrom^3 (approx 79.0 for g=2 and 1H)
    # T = (mu0/4pi) * (ge * beta_e * gn * beta_n) / r^3
    # For protons, a common shortcut is: 79.0 * (g_e/2.0023) / r^3
    dipolar_const = 79.0 * (params['g_e'] / 2.0023)
    T = dipolar_const / r**3
    
    # 2. Calculate Dipolar Component at angle theta
    A_dip = T * (3 * np.cos(theta)**2 - 1)
    
    # 3. Calculate Isotropic (Contact) Component
    # We use (r - R_np) to measure distance from the surface
    # np.maximum(0, ...) ensures we don't get negative distances inside the core
    dist_from_surface = np.maximum(0, r - params['R_np'])
    A_iso = params['A0'] * np.exp(-dist_from_surface / params['rho'])
    
    # Total coupling is the sum
    return A_iso + A_dip


def simulate_shell(r, density, dr, params):
    # 1. Convert density to Number of Atoms
    N = density * 4 * np.pi * r**2 * dr
    
    # 2. Get the Pake Pattern (averaging theta from 0 to pi/2)
    # ... (powder average logic here) ...
    
    # 3. Multiply the final Pake pattern by N
    return pake_pattern * N

class NPDensityProfile:
    def __init__(self, surface_distances, R_np, bin_width=0.2, bulk_density=0.11, old_bins = None):
        """
        Processes MD distance lists into a radial density profile.
        
        Parameters:
        surface_distances : ndarray
                            Raw list of proton distances from the NP surface.
        R_np              : float
                            Radius of the metal core (Angstroms).
        bin_width         : float
                            Width of radial bins (Angstroms).
        bulk_density      : float
                            Proton density (H/A^3) beyond the MD box.
        bin_start           : array
                            A pre-existing array of bin minds
        """
        self.R_np = R_np
        self.bulk_density = bulk_density
        
        # 1. Convert to distance from center of the nanoparticle
        r_from_center = surface_distances + R_np
        

        if old_bins is None: 
            # 2. Bin the data
            r_max_md = np.max(r_from_center)
        else:
            r_max_md = np.max(old_bins)
            bin_width = old_bins[1]-old_bins[0]

        bins = np.arange(R_np, r_max_md + bin_width, bin_width)
        counts, edges = np.histogram(r_from_center, bins=bins)
        
        # 3. Calculate shell volumes and density
        r_inner = edges[:-1]
        r_outer = edges[1:]
        volumes = (4/3) * np.pi * (r_outer**3 - r_inner**3)
        self.densities = counts / volumes
        self.bin_centers = (r_inner + r_outer) / 2.0
        
    def get_density(self, r, rho_bulk, r_stitch_center=50, blend_width=3.0):
        """
        Returns a smooth, blended density for any radius r.
        This is the 'Secret Sauce' that makes R_max a smooth fitting parameter.
        """
        if r < self.R_np:
            return 0.0
        
        # MD density at this point
        rho_md = np.interp(r, self.bin_centers, self.densities)
        
        # Blending logic
        r_start = r_stitch_center - blend_width/2
        r_end = r_stitch_center + blend_width/2
        
        if r <= r_start:
            return rho_md
        elif r >= r_end:
            return rho_bulk
        else:
            # Linear blend between MD and Bulk
            w = (r - r_start) / blend_width
            return (1 - w) * rho_md + w * rho_bulk

            

    def get_counts(self, r, dr):
        """Converts density at r back to atom counts for a shell of thickness dr."""
        rho = self.get_density(r)
        return rho * (4 * np.pi * r**2) * dr

    def find_blend_start(self, method = "fit"):
        if method == "fit":
            # fit density to a gassian and use that
            pass
        elif method == "max":
            # just take the maximum value
            pass

        return r_blend_start

    def prepare_blended_density(self, rho_bulk, blend_width=3.0):
        # change this such that you find the bin you need after the blending window, and then end there.  
        # There is no need to go beyond this since we can simply calcluate this on the fly...

        """
        Creates a density profile from R_np to R_max with a tunable blend transition.
        
        Parameters:
        r_centers   : Surface distances from MD histogram.
        densities   : Density values from MD.
        R_np        : Nanoparticle radius.
        R_max       : Integration horizon.
        rho_bulk    : Target bulk density (0.0453 for toluene).
        d_stitch    : The surface distance calculated via the 'Fit-to-Stitch' logic.
        blend_width : Width of the linear transition (Angstroms).
        """
        # 1. Define the grid and the blend boundaries


        r_stitch_center = 50#get_r_start()

        r_stitch_start = r_stitch_center - (blend_width / 2.0)
        r_stitch_end = r_stitch_center + (blend_width / 2.0)
        bin_width = self.bin_centers[1] - self.bin_centers[0]
        bin_start = self.bin_centers[0]

        print(r_stitch_end)

        final_centers = []
        final_rhos = []
        finished = False
        for i, r in enumerate(self.bin_centers):
            if r <= r_stitch_start:
                # Region 1: Pure MD Data
                final_rhos.append(self.densities[i])

            elif r > r_stitch_end:
                # we have reached the end, and we need to just add one more bulk density bin
                final_rhos.append(rho_bulk)
                print("found end!")
                break #exit the control structure
            else:
                # Region 2: Linear Transition
                rho_md = self.densities[i]
                # Calculate weight w (0 at d_start, 1 at d_end)
                w = (r - r_stitch_start) / blend_width
                final_rhos.append((1 - w) * rho_md + w * rho_bulk)

            final_centers.append(bin_start + i*bin_width)
        
        self.bin_centers, self.densities = final_centers, final_rhos
        return

    