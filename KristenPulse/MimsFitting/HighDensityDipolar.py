import marimo

__generated_with = "0.20.4"
app = marimo.App()


@app.cell
def _():
    #import all needed modules
    import marimo as mo
    import numpy as np
    from importlib import reload
    import codechembook.quickPlots as qp
    from scipy.ndimage import gaussian_filter1d
    from plotly.subplots import make_subplots
    from lmfit import Model

    #import modules specifically written for this 
    import pulsePhysics as pp
    import mdHelper as mh

    reload(pp)
    reload(mh)

    return Model, make_subplots, mh, np, pp


@app.cell
def _(np):
    # read spectral data files
    dodecaneH_tolueneH_spectra =  np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)
    dodecaneD_tolueneH_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter="," ,unpack = True, skip_header = 1)
    dodecaneH_tolueneD_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter="," ,unpack = True, skip_header = 1)
    return (
        dodecaneD_tolueneH_spectra,
        dodecaneH_tolueneD_spectra,
        dodecaneH_tolueneH_spectra,
    )


@app.cell
def _(np):
    # read md data files

    # high density files
    hd_toluene_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/toluene_distances_high_density.csv",
                                       delimiter = ",")
    hd_dodecanethiol_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/dodecanol_distances_high_density.csv",
                                       delimiter = ",")
    return hd_dodecanethiol_distances, hd_toluene_distances


@app.cell
def _():
    # read np size data
    return


@app.cell
def _():
    # constants
    R_np_md = 15  #radius in Angstroms for the particle in teh MD simluations
    bulk_toluene_density = 0.0453  # number of hydrogens per cubic angstrom
    bulk_dodecanethiol_density = 0 # the density of hydrogens coming from the ligand at r = infinity
    return R_np_md, bulk_toluene_density


@app.cell
def _(
    R_np_md,
    bulk_toluene_density,
    hd_dodecanethiol_distances,
    hd_toluene_distances,
    make_subplots,
    mh,
    np,
):
    # process md data files

    #toluene
    hd_toluene_profile = mh.NPDensityProfile(hd_toluene_distances, R_np_md, bin_width=0.1)
    hd_toluene_profile.prepare_blended_density(rho_bulk = bulk_toluene_density)

    #ligand
    hd_dodecanethiol_profile = mh.NPDensityProfile(hd_dodecanethiol_distances, R_np_md, old_bins = hd_toluene_profile.bin_centers)


    densityPLot  = make_subplots()

    densityPLot.add_scatter(x = np.array(hd_dodecanethiol_profile.bin_centers)/10, y = hd_dodecanethiol_profile.densities, 
                                  line = dict(color = "#BB6F35"), fill = "tozeroy")
    densityPLot.add_scatter(x = np.array(hd_toluene_profile.bin_centers)/10, y = hd_toluene_profile.densities, 
                           line = dict(color = "#448CC4"), fill = "tozeroy")
    #densityPLot.add_scatter(x = np.array(hd_toluene_profile.bin_centers)/10, y = hd_toluene_profile.densities+hd_dodecanethiol_profile.densities, line = dict(color = "#94A2B1"))

    densityPLot.update_layout(template = "simple_white",showlegend = False, width = 300*2, height = 300*1.4)
    densityPLot.update_xaxes(range = [0, 3.4], title = "radial distance /nm")
    densityPLot.update_yaxes(title = "H density")
    densityPLot.show()
    return hd_dodecanethiol_profile, hd_toluene_profile


@app.cell
def _(
    D_np,
    H_np,
    I1_guess,
    I2_guess,
    I3_guess,
    I_scale_guess,
    Model,
    R_max_guess,
    T_m,
    c1_guess,
    c2_guess,
    c3_guess,
    contact_scale_guess,
    dipolar_scale_guess,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    ell_guess,
    eta_guess,
    g_e,
    g_n,
    len1,
    len2,
    ligand_density_obj,
    np,
    nu_n_guess,
    pp,
    rho_bulk_ligand,
    rho_bulk_solvent,
    shape,
    solvent_density_obj,
    target_dr,
    tau,
    width_scale_guess,
):
    # fit

    # 1 define the function that will be iterated in the fitting
    def triple_mims_simulation(
        xs, # the x-values to calculate intensities for. These are concatenated from the spectra

        # properties of just the final spectra we are fitting to
        len1, len2, # to be used for indices
        I1, I2, I3, # final intensity for each spectra
        c1, c2, c3, # final background for each spectra
    
        # properties of the nanoparticles
        H_np, D_np, # objects containing size info for the protio and duetero ligand protected nanoparticles. This should be paired lists of bin center and frequency.

        # properties of the nuclei
        solvent_density_obj, rho_bulk_solvent,
        ligand_density_obj, rho_bulk_ligand,
        nu_n, # central frequency... if we correct the mims, could be zero
        g_n, # nuclear g-factor, depends on the atom identity

        # properties of the electron
        I_scale, # unitless, a scalar that adjust intenisty, but should default to 1.
        ell, # unit of m, characteristic spatial distance of the decay
        g_e, # electronic g-factor

        #properties of the EPR experiment
        tau, # unit = s, delay time of the pulse sequence
        T_m, # unit = s, coherence lifetime

        # properties of the simulation
        R_max, # unit = m, how far out in the simulation to look 
        target_dr = 0.2E-9, # unit = m, 
        dipolar_scale = 1, 
        contact_scale = 1, 

        # properties of the spectra
        shape = "Lorentzian",
        eta = 0,
        width_scale = 1,
        ):
        # should sigma vary by sample? Or T_M?


        ys = np.array([]) # to hold the ys we calculate, concatenated together, so we can 

        #separate out th xs...
        freq_both =         xs[         :len1     ]
        freq_solvent_only = xs[len1     :len1+len2]
        freq_ligand_only =  xs[len1+len2:         ]

        # for each sample, calculate the contribution from the 

        for freq, I, c, s, l, rs in zip(
            [freq_both, freq_solvent_only, freq_ligand_only],
            [I1, I2, I3],
            [c1, c2, c3],
            [1, 1, 0.01], # solvent weighting, based on deuterium percentage
            [1, 0.02, 1], # ligand weighting, based on dueterium percentage
            [H_np, D_np, H_np],
        ):
            for r, d in zip(rs[0], rs[1]): # go through each particle size 
                y_solvent = pp.simulate_nanoparticle_mims_spectrum(
                    freq,
                
                    r,

                    solvent_density_obj,
                    rho_bulk_solvent,
                    nu_n, 
                    g_n, 

                    I_scale,
                    ell,
                    g_e,

                    tau,
                    T_m,

                    R_max, 
                    target_dr, 
                    dipolar_scale, 
                    contact_scale, 

                    shape,
                    eta, 
                    width_scale
                )

                y_ligand = pp.simulate_nanoparticle_mims_spectrum(
                    freq,
                
                    r,

                    ligand_density_obj,
                    rho_bulk_ligand,
                    nu_n, 
                    g_n, 

                    I_scale,
                    ell,
                    g_e,

                    tau,
                    T_m,

                    R_max, 
                    target_dr, 
                    dipolar_scale, 
                    contact_scale, 

                    shape,
                    eta, 
                    width_scale
                )
    
    
                raw_ys = d*(s*y_solvent + l*y_ligand) # weight contribution of the ligand and solvent by percentage H, and then also weight the total by the density of the nanoparticle size we are looking at
        
            norm_ys = raw_ys/np.max(raw_ys) # normalize all the summed contributions

            scaled_ys = norm_ys*I + c # then adjust by the intensity and background

            ys = np.concatenate([ys, scaled_ys]) # concatenate the final scaled ys for the sample to the total


        return ys

    # second, make a model of this function
    triple_mims_simulation_model = Model(triple_mims_simulation, independent_vars = [
        "xs", "len1", "len2", 
        "H_np", "D_np", 
        "solvent_density_obj", "rho_bulk_solvent", "ligand_density_obj", "rho_bulk_ligand", "g_n", 
        "g_e", 
        "tau", "T_m", 
        "target_dr",
        "shape",])

    # third, make parameters and define them
    triple_mims_simulation_params = triple_mims_simulation_model.make_params()

    triple_mims_simulation_params.add_many(
        ("I1", I1_guess, True, 0, None), 
        ("I2", I2_guess, True, 0, None), 
        ("I3", I3_guess, True, 0, None), # final intensity for each spectra
        ("c1", c1_guess, True, None, None), 
        ("c2", c2_guess, True, None, None),
        ("c3", c3_guess, True, None, None), # final background for each spectra   
    
        ("nu_n", nu_n_guess, True, None, None), # in Hz
    
        ("I_scale", I_scale_guess, False, 0, None), 
        ("ell", ell_guess, True, 0, None), # in meters
    
        ('R_max', R_max_guess, False, 0, None),
    
        ('dipolar_scale', dipolar_scale_guess, False, 0, None),
        ('contact_scale', contact_scale_guess, False, 0, None),    
        ("eta", eta_guess, False, 0, None),
        ("width_scale", width_scale_guess, False, 0, None),
        ) 
    #triple_spectra_model_params['R_max'].diff_step = diff_step=dodecaneH_tolueneH_spectra[0][1] - dodecaneH_tolueneH_spectra[0][0] # Make diff_step for R_max roughly equal to your dr

    # forth, fit the model to the data
    triple_mims_simulation_fit = triple_mims_simulation_model.fit(
        np.concatenate([dodecaneH_tolueneH_spectra[1], dodecaneD_tolueneH_spectra[1], dodecaneH_tolueneD_spectra[1]]), 
        xs = np.concatenate([dodecaneH_tolueneH_spectra[0], dodecaneD_tolueneH_spectra[0], dodecaneH_tolueneD_spectra[0]]),
        params = triple_mims_simulation_params,
    
        len1 = len1,
        len2 = len2,

        H_np = H_np, 
        D_np = D_np,
    
        solvent_density_obj = solvent_density_obj, 
        rho_bulk_solvent = rho_bulk_solvent,                  
        ligand_density_obj = ligand_density_obj, 
        rho_bulk_ligand = rho_bulk_ligand,
        g_n = g_n, 
    
        g_e = g_e,
    
        tau = tau,
        T_m = T_m,
    
        target_dr = target_dr,

        shape = shape, 
    
        iter_cb=my_callback
    )
    return


@app.cell
def _(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    hd_dodecanethiol_profile,
    hd_toluene_profile,
    np,
):
    # VALUES for parameters that are not fit
    len1 = len(dodecaneH_tolueneH_spectra[0])
    len2 = len(dodecaneD_tolueneH_spectra[0])

    H_np = [[1.5e-9],[1.0]] # bin center in m, frequency
    D_np = [[1.5e-9],[1.0]] # bin center in m, frequency

    solvent_density_obj = hd_toluene_profile
    rho_bulk_solvent = 0.0453 # number of hydrogens per cubic angstrom
    ligand_density_obj = hd_dodecanethiol_profile
    rho_bulk_ligand = 0
    g_n = 5.5856946893 # from scipy constants

    g_e = 1.9 # g-value for the nanoparticles

    tau=140e-6 #
    T_m = 140e-6 # in seconds.

    target_dr = 0.1e-9

    shape = "Lorentzian"



    # GUESSES for parameters to be fit
    I1_guess = np.max(dodecaneH_tolueneH_spectra[1]) - np.min(dodecaneH_tolueneH_spectra[1])
    I2_guess = np.max(dodecaneD_tolueneH_spectra[1]) - np.min(dodecaneD_tolueneH_spectra[1])
    I3_guess = np.max(dodecaneH_tolueneD_spectra[1]) - np.min(dodecaneH_tolueneD_spectra[1])

    c1_guess = np.mean(dodecaneH_tolueneH_spectra[1][0:10])
    c2_guess = np.mean(dodecaneD_tolueneH_spectra[1][0:10])
    c3_guess = np.mean(dodecaneH_tolueneD_spectra[1][0:10])

    nu_n_guess = 50e6 # units = s

    I_scale_guess = 1 # unitless

    ell_guess = 1.3e9 # unit = m the decay constant for wavefunction

    R_max_guess = 50e-9 # m, the max radius we look out to

    dipolar_scale_guess = 1 # unitless
    contact_scale_guess = 1 # unitless

    eta_guess = 0 # unitless
    width_scale_guess = 1 # unitless





    return (
        D_np,
        H_np,
        I1_guess,
        I2_guess,
        I3_guess,
        I_scale_guess,
        R_max_guess,
        T_m,
        c1_guess,
        c2_guess,
        c3_guess,
        contact_scale_guess,
        dipolar_scale_guess,
        ell_guess,
        eta_guess,
        g_e,
        g_n,
        len1,
        len2,
        ligand_density_obj,
        nu_n_guess,
        rho_bulk_ligand,
        rho_bulk_solvent,
        shape,
        solvent_density_obj,
        target_dr,
        tau,
        width_scale_guess,
    )


@app.cell
def _():
    # fit and save/display results
    return


@app.function
# random functions

def my_callback(params, iter, resid, *args, **kws):
    # Only execute the print if 'iter' is a multiple of 1000
    if iter % 100 == 0:
        ssr = (resid**2).sum()
        print(f"Iteration: {iter:5d} | Chi-square (SSR): {ssr:.4e}")

    return False


@app.cell
def _():
    # simulation of spectra, just to test...
    return


@app.cell
def _():
    #TESTING

    return


if __name__ == "__main__":
    app.run()
