import numpy as np

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