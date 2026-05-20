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
    import plotly.graph_objects as go

    #import modules specifically written for this 
    import pulsePhysics as pp
    import mdHelper as mh

    reload(pp)
    reload(mh)
    return Model, make_subplots, mh, mo, np, pp, qp, reload


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    To do:
    - [ ] Seems like the errors might be floating point number issues, so what if instead of SI units, we change to... nm and ns? Actually, not totally sure what the best units are
    """)
    return


@app.cell
def _(make_subplots, np):
    # read spectral data files
    dodecaneH_tolueneH_spectra =  np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)
    dodecaneD_tolueneH_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter="," ,unpack = True, skip_header = 1)
    dodecaneH_tolueneD_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter="," ,unpack = True, skip_header = 1)

    # convert from MHz to Hz
    dodecaneH_tolueneH_spectra[0] = dodecaneH_tolueneH_spectra[0]*1e6
    dodecaneD_tolueneH_spectra[0] = dodecaneD_tolueneH_spectra[0]*1e6
    dodecaneH_tolueneD_spectra[0] = dodecaneH_tolueneD_spectra[0]*1e6

    expPlot = make_subplots(rows = 3, cols = 1)
    for _i, _spectrum in enumerate([dodecaneH_tolueneH_spectra, dodecaneD_tolueneH_spectra, dodecaneH_tolueneD_spectra], start = 1):
        expPlot.add_scatter(x = _spectrum[0], y = _spectrum[1], row = _i, col = 1)
    expPlot.update_layout(template = "simple_white")
    expPlot.show()

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
                                       delimiter = ",")*1e-10 # convert to m
    hd_dodecanethiol_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/dodecanol_distances_high_density.csv",
                                       delimiter = ",")*1e-10 # convert to m

    return hd_dodecanethiol_distances, hd_toluene_distances


@app.cell
def _():
    # read np size data
    return


@app.cell
def _():
    # constants
    R_np_md = 15e-10  #radius in m for the particle in teh MD simluations
    bulk_toluene_density = 0.0453e30  # number of hydrogens per cubic m (used to be angtroms)
    bulk_dodecanethiol_density = 0 # the density of hydrogens coming from the ligand at r = infinity
    md_r_step = 0.1e-10 # units = m, 
    return R_np_md, bulk_toluene_density, md_r_step


@app.cell
def _(
    R_np_md,
    bulk_toluene_density,
    hd_dodecanethiol_distances,
    hd_toluene_distances,
    make_subplots,
    md_r_step,
    mh,
    np,
):
    # process md data files

    #toluene
    hd_toluene_profile = mh.NPDensityProfile(hd_toluene_distances, R_np_md, bin_width=md_r_step)
    hd_toluene_profile.prepare_blended_density(rho_bulk = bulk_toluene_density)

    #ligand
    hd_dodecanethiol_profile = mh.NPDensityProfile(hd_dodecanethiol_distances, R_np_md, old_bins = hd_toluene_profile.bin_centers)


    densityPLot  = make_subplots()

    densityPLot.add_scatter(x = np.array(hd_dodecanethiol_profile.bin_centers), y = hd_dodecanethiol_profile.densities, 
                                  line = dict(color = "#BB6F35"), fill = "tozeroy", name = "ligand")
    densityPLot.add_scatter(x = np.array(hd_toluene_profile.bin_centers), y = hd_toluene_profile.densities, 
                           line = dict(color = "#448CC4"), fill = "tozeroy", name = "solvent")
    #densityPLot.add_scatter(x = np.array(hd_toluene_profile.bin_centers)/10, y = hd_toluene_profile.densities+hd_dodecanethiol_profile.densities, line = dict(color = "#94A2B1"))

    densityPLot.update_layout(template = "simple_white",
                              #showlegend = False, 
                              width = 300*2, height = 300*1.4)
    densityPLot.update_xaxes(range = [0, 34e-10], title = "radial distance /nm")
    densityPLot.update_yaxes(title = "H density")
    densityPLot.show()
    return hd_dodecanethiol_profile, hd_toluene_profile


@app.cell
def _(np, pp):
    ### FUNCTION FOR FITTING

    # 1 define the function that will be iterated in the fitting. This is the function that fits all three spectra together. 
    def triple_mims_simulation(
        xs, # the x-values to calculate intensities for. These are concatenated from the spectra

        # properties of just the final spectra we are fitting to
        len1, len2, # to be used for indices
        I1, I2, I3, # final intensity for each spectra
        c1, c2, c3, # final background for each spectra

        # properties of the experimental nanoparticles
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
        R_np_md, # radius of the simulated nanoparticle
        R_max, # unit = m, how far out in the simulation to look 
        target_dr = 0.2E-10, # unit = m, 
        dipolar_scale = 1, 
        contact_scale = 1, 

        # properties of the spectra
        shape = "Lorentzian",
        eta = 0,
        width_scale = 1,
        ):
        # correct units
        I1 = I1*1e6
        I2 = I2*1e6
        I3 = I3*1e6

    
        # should sigma vary by sample? Or T_M?


        ys = np.array([]) 

        # Separate the x-axes
        freq_both =         xs[:len1]
        freq_solvent_only = xs[len1:len1+len2]
        freq_ligand_only =  xs[len1+len2:]

        # Global container to hold unnormalized raw spectra for relative scaling check
        all_raw_samples = []

        # Loop through each of the 3 experimental configurations
        for name, freq, s_weight, l_weight, rs in zip(
            ["both", "solvent only", "ligand only"],
            [freq_both, freq_solvent_only, freq_ligand_only], # frequency info
            [1.0, 1.0, 0.01],  # Solvent H-percentage
            [1.0, 0.02, 1.0],  # Ligand H-percentage
            [H_np, D_np, H_np]  # Nanoparticle size distributions
        ):
            #print(f"sample = {name}")
            # Initialize an empty spectrum for the current sample to sum into
            sample_total_raw = np.zeros_like(freq)

            # Integrate across the entire nanoparticle size distribution
            for r, weight_d in zip(rs[0], rs[1]): # unpack the rs  
            
                # Optimization: Only calculate solvent if it actually contributes
                #print(f"r = {r}")
                if s_weight > 0.0:
                    #print("calculating solvent")
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
                    
                        R_np_md, 
                        R_max, 
                        target_dr, 
                        dipolar_scale, 
                        contact_scale, #contact_scale, # set to 0 to turn off for solvent. 
                    
                        shape, 
                        eta, 
                        width_scale
                        )
                else:
                    y_solvent = np.zeros_like(freq)

                # Optimization: Only calculate ligand if it actually contributes
                if l_weight > 0.00:
                    #print("calculating ligand")
                    y_ligand = pp.simulate_nanoparticle_mims_spectrum(freq, r, ligand_density_obj, rho_bulk_ligand, nu_n, g_n, I_scale, ell, g_e, tau, T_m, R_np_md, R_max, target_dr, dipolar_scale, contact_scale, shape, eta, width_scale)
                else:
                    y_ligand = np.zeros_like(freq)

                # ACCUMULATE: Add this size bin's weighted contribution to the total
                sample_total_raw += weight_d * (s_weight * y_solvent + l_weight * y_ligand)

            # Store the unnormalized integrated spectrum for this sample
            all_raw_samples.append(sample_total_raw)

        # --- SCALING AND CONCATENATION ---
        # Apply intensity scaling and background offsets to each integrated spectrum
        # If your experimental data is already normalized to 1, scaling individually here is fine,
        # but using a global max preserves the true physical intensity ratios.
    
        scaled_both    = (all_raw_samples[0] / np.max(all_raw_samples[0])) * I1 + c1
        scaled_solvent = (all_raw_samples[1] / np.max(all_raw_samples[1])) * I2 + c2
        scaled_ligand  = (all_raw_samples[2] / np.max(all_raw_samples[2])) * I3 + c3

        # Concatenate back into a single vector matching the structure of 'xs'
        ys = np.concatenate([scaled_both, scaled_solvent, scaled_ligand])

        return  ys





    return (triple_mims_simulation,)


@app.cell
def _(
    D_np,
    H_np,
    I1_guess,
    I2_guess,
    I3_guess,
    I_scale_guess,
    R_max_guess,
    R_np_md,
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
    reload,
    rho_bulk_ligand,
    rho_bulk_solvent,
    shape,
    solvent_density_obj,
    target_dr,
    tau,
    triple_mims_simulation,
    width_scale_guess,
):
    reload(pp)
    # simulate, before fit
    sim_spectra = triple_mims_simulation(
        np.concatenate([dodecaneH_tolueneH_spectra[0], dodecaneD_tolueneH_spectra[0], dodecaneH_tolueneD_spectra[0]]), # the x-values to calculate intensities for. These are concatenated from the spectra

        # properties of just the final spectra we are fitting to
        len1, len2, # to be used for indices
        I1_guess, I2_guess, I3_guess, # final intensity for each spectra
        c1_guess, c2_guess, c3_guess, # final background for each spectra

        # properties of the nanoparticles
        H_np, D_np, # objects containing size info for the protio and duetero ligand protected nanoparticles. This should be paired lists of bin center and frequency.

        # properties of the nuclei
        solvent_density_obj, rho_bulk_solvent,
        ligand_density_obj, rho_bulk_ligand,
        nu_n_guess, # central frequency... if we correct the mims, could be zero
        g_n, # nuclear g-factor, depends on the atom identity

        # properties of the electron
        I_scale_guess, # unitless, a scalar that adjust intenisty, but should default to 1.
        ell_guess, # unit of m, characteristic spatial distance of the decay
        g_e, # electronic g-factor

        #properties of the EPR experiment
        tau, # unit = s, delay time of the pulse sequence
        T_m, # unit = s, coherence lifetime

        # properties of the simulation
        R_np_md, 
        R_max_guess, # unit = m, how far out in the simulation to look 
        target_dr, # unit = m, 
        dipolar_scale_guess, 
        contact_scale_guess, 

        # properties of the spectra
        shape,
        eta_guess,
        width_scale_guess,
        )


    return (sim_spectra,)


@app.cell
def _(
    c1_guess,
    c2_guess,
    c3_guess,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    len1,
    len2,
    make_subplots,
    np,
    sim_spectra,
    tau,
):
    sim_compare  = make_subplots(rows = 5, cols = 1)#go.Figure(expPlot)

    # compare raw intensities of the simulated
    sim_compare.add_scatter(x = dodecaneH_tolueneH_spectra[0], y = (sim_spectra[:len1] - c1_guess)/1, row = 1, col = 1, showlegend = False, line = dict(color = "blue"))
    sim_compare.add_scatter(x = dodecaneD_tolueneH_spectra[0], y = (sim_spectra[len1:len1+len2] - c2_guess)/1, row = 1, col = 1, showlegend = False, line = dict(color = "orange"))
    sim_compare.add_scatter(x = dodecaneH_tolueneD_spectra[0], y = (sim_spectra[len1+len2:] - c3_guess)/1, row = 1, col = 1, showlegend = False, line = dict(color = "green"))
    #sim_compare.add_scatter(x = dodecaneH_tolueneH_spectra[0], y = (1 - np.cos(2 * np.pi * dodecaneH_tolueneH_spectra[0] * tau))*np.max(sim_spectra[len1:len1+len2] - c2_guess)/2   , name = "blindspot", line = dict(color = "red"))


    # compare normalized intensities of the simulated
    sim_compare.add_scatter(x = dodecaneH_tolueneH_spectra[0], y = (sim_spectra[:len1] - c1_guess)/np.max(sim_spectra[:len1] - c1_guess), row = 2, col = 1, name = "both sim", line = dict(color = "blue"))
    sim_compare.add_scatter(x = dodecaneD_tolueneH_spectra[0], y = (sim_spectra[len1:len1+len2] - c2_guess)/np.max(sim_spectra[len1:len1+len2] - c2_guess), row = 2, col = 1, name = "solvent only sim", line = dict(color = "orange"))
    sim_compare.add_scatter(x = dodecaneH_tolueneD_spectra[0], y = (sim_spectra[len1+len2:] - c3_guess)/np.max(sim_spectra[len1+len2:] - c3_guess), row = 2, col = 1, name = "ligand only sim", line = dict(color = "green"))
    sim_compare.add_scatter(x = dodecaneH_tolueneH_spectra[0], y = (1 - np.cos(2 * np.pi * dodecaneH_tolueneH_spectra[0] * tau))/2 , row = 2, col = 1, name = "blindspot",line = dict(color = "red"))



    #experimental
    for _i, _spectrum in enumerate([dodecaneH_tolueneH_spectra, dodecaneD_tolueneH_spectra, dodecaneH_tolueneD_spectra], start = 3):
        sim_compare.add_scatter(x = _spectrum[0], y = _spectrum[1], row = _i, col = 1, showlegend = False, line = dict(color = "grey"))

    #simulated
    sim_compare.add_scatter(x = dodecaneH_tolueneH_spectra[0], y = sim_spectra[:len1], row = 3, col = 1, showlegend = False, line = dict(color = "blue"))
    sim_compare.add_scatter(x = dodecaneD_tolueneH_spectra[0], y = sim_spectra[len1:len1+len2], row = 4, col = 1, showlegend = False, line = dict(color = "orange"))
    sim_compare.add_scatter(x = dodecaneH_tolueneD_spectra[0], y = sim_spectra[len1+len2:], row = 5, col = 1, showlegend = False, line = dict(color = "green"))

    sim_compare.update_xaxes(range = [47e6, 53e6])
    sim_compare.update_layout(template = "simple_white", height = 700)
    sim_compare.show()
    return


@app.cell
def _(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    hd_dodecanethiol_profile,
    hd_toluene_profile,
    np,
    pp,
    reload,
):
    reload(pp)

    # VALUES for parameters that are not fit
    len1 = len(dodecaneH_tolueneH_spectra[0])
    len2 = len(dodecaneD_tolueneH_spectra[0])

    H_np = [[1.5e-9],[1.0]] # bin center in m, frequency
    D_np = [[1.5e-9],[1.0]] # bin center in m, frequency

    solvent_density_obj = hd_toluene_profile
    rho_bulk_solvent = 4.528e28 # number of hydrogens per cubic m (used to be cubic angstrom)
    ligand_density_obj = hd_dodecanethiol_profile
    rho_bulk_ligand = 0
    g_n = 5.5856946893 # from scipy constants

    g_e = 2.0 # g-value for the nanoparticles

    tau = 140e-9 # in seconds
    T_m = 674e-9 # in seconds.

    target_dr = 0.1e-9

    shape = "L"

    # GUESSES for parameters to be fit
    I1_guess = (np.max(dodecaneH_tolueneH_spectra[1]) - np.min(dodecaneH_tolueneH_spectra[1]))*1e-6
    I2_guess = (np.max(dodecaneD_tolueneH_spectra[1]) - np.min(dodecaneD_tolueneH_spectra[1]))*1e-6
    I3_guess = (np.max(dodecaneH_tolueneD_spectra[1]) - np.min(dodecaneH_tolueneD_spectra[1]))*1e-6

    c1_guess = np.mean(dodecaneH_tolueneH_spectra[1][0:10])
    c2_guess = np.mean(dodecaneD_tolueneH_spectra[1][0:10])
    c3_guess = np.mean(dodecaneH_tolueneD_spectra[1][0:10])

    nu_n_guess = 50e6 # units = s

    I_scale_guess = 5 # unitless

    ell_guess = 0.7 #e-9 # unit = m the decay constant for wavefunction

    R_max_guess = 20 #e-9 # m, the max radius we look out to

    dipolar_scale_guess = 1 # unitless
    contact_scale_guess = 1 # unitless

    eta_guess = 0 # unitless
    width_scale_guess = 0.8 # unitless


    print(dodecaneH_tolueneH_spectra[0][1] - dodecaneH_tolueneH_spectra[0][0])

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
def fitting(
    D_np,
    H_np,
    I1_guess,
    I2_guess,
    I3_guess,
    I_scale_guess,
    Model,
    R_max_guess,
    R_np_md,
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
    reload,
    rho_bulk_ligand,
    rho_bulk_solvent,
    shape,
    solvent_density_obj,
    target_dr,
    tau,
    triple_mims_simulation,
    width_scale_guess,
):
    reload(pp)
    # fit and save/display results
    def my_callback(params, iter, resid, *args, **kws):
        # Only execute the print if 'iter' is a multiple of 1000
        if iter % 100 == 0:
            ssr = (resid**2).sum()
            print(f"Iteration: {iter:5d} | Chi-square (SSR): {ssr:.4e}")

        return False



    triple_spectra_model = Model(triple_mims_simulation, independent_vars = [
        "xs", 
        "len1", "len2",   
        "H_np", "D_np",
        "solvent_density_obj", "rho_bulk_solvent",
        "ligand_density_obj", "rho_bulk_ligand",
        "g_n",
        "g_e",
        "tau", 
        "T_m", 
        "R_np_md",
        "target_dr", 
        "shape",
        ])


    triple_spectra_model_params = triple_spectra_model.make_params()

    triple_spectra_model_params.add_many(
        ("I1", I1_guess, True, 0, None), 
        ("I2", I2_guess, True, 0, None), 
        ("I3", I3_guess, True, 0, None), # final intensity for each spectra
        ("c1", c1_guess, False, None, None), 
        ("c2", c2_guess, False, None, None),
        ("c3", c3_guess, False, None, None), # final background for each spectra   

        ("nu_n", nu_n_guess, False, None, None), # in Hz

        ("I_scale", I_scale_guess, True, 0, None), 
        ("ell", ell_guess, True, 0, None), # bring in as nm

        ('R_max', R_max_guess, True, 0, None), # bring in as nm

        ('dipolar_scale', dipolar_scale_guess, False, 0, None),
        ('contact_scale', contact_scale_guess, False, 0, None),    
        ("eta", eta_guess, False, 0, None),
        ("width_scale", width_scale_guess, True, 0, None),
        ) 

    # add in parameters for fitting: 
    # R_scale <- parameter to control how far from the nanoparticle surface to adjust when calculating dipolar coupling
    # blind_spot_offset <- how imperfect the blindspot is. 
    # fi <- central resonance for the spectrum.  OR we can correct the spectra to 0? 

    # remove parameters for fitting: 



    triple_spetra_fit = triple_spectra_model.fit(
        np.concatenate([dodecaneH_tolueneH_spectra[1], dodecaneD_tolueneH_spectra[1], dodecaneH_tolueneD_spectra[1]]), 
        xs = np.concatenate([dodecaneH_tolueneH_spectra[0], dodecaneD_tolueneH_spectra[0], dodecaneH_tolueneD_spectra[0]]),
        params = triple_spectra_model_params, 
    
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
        R_np_md = R_np_md,
        target_dr = target_dr, 
        shape = shape,
        iter_cb=my_callback
    )

    return (triple_spetra_fit,)


@app.cell
def _(triple_spetra_fit):
    # fit report

    print(triple_spetra_fit.fit_report())

    return


@app.cell
def _(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    len1,
    len2,
    make_subplots,
    np,
    triple_spetra_fit,
):
    # plot

    fitplot = make_subplots(rows = 4, cols = 1, subplot_titles= ["solvent-H + ligand-H", "solvent-H + ligand-D", "solvent-D + ligand-H"])
    fitplot.update_layout(template = "simple_white", width = 500, height = 500)

    fitplot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneH_spectra[1], line = dict(color = "#ffabff", width = 5), showlegend = False, col = 1, row = 1)
    fitplot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[:len1], line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    fitplot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneD_tolueneH_spectra[1], line = dict(color = "#ffabac", width = 5), showlegend = False, col = 1, row = 2)
    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "black"), showlegend = False, col = 1, row = 2)

    fitplot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneD_spectra[1], line = dict(color = "#ababff", width = 5), showlegend = False, col = 1, row = 3)
    fitplot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= triple_spetra_fit.best_fit[len1+len2:], line = dict(color = "black"), showlegend = False, col = 1, row = 3)

    #add normalized spectra
    fitplot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= (triple_spetra_fit.best_fit[:len1] - triple_spetra_fit.params["c1"].value)/np.max((triple_spetra_fit.best_fit[:len1] - triple_spetra_fit.params["c1"].value)), line = dict(color = "black"), showlegend = False, col = 1, row = 4)
    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= (triple_spetra_fit.best_fit[len1:len1+len2] - triple_spetra_fit.params["c2"].value)/np.max(triple_spetra_fit.best_fit[len1:len1+len2] - triple_spetra_fit.params["c2"].value), line = dict(color = "black"), showlegend = False, col = 1, row = 4)
    fitplot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= (triple_spetra_fit.best_fit[len1+len2:] - triple_spetra_fit.params["c3"].value)/np.max(triple_spetra_fit.best_fit[len1+len2:] - triple_spetra_fit.params["c3"].value), line = dict(color = "black"), showlegend = False, col = 1, row = 4)


    fitplot.update_xaxes(title = "frequency /MHz", row = 3, col = 1)
    fitplot.show()
    return


@app.cell
def _(np, qp):
    #TESTING

    #visualizing the blindplot
    x = np.linspace(-5e6, 5e6, 10000)
    intensity = (1 - np.cos(2 * np.pi * x * 140e-9))

    qp.quickScatter(x = x, y = intensity)
    qp.quickScatter(x = x, y = 0.5* np.sin(2*np.pi*x/2*140e-9)**2)

    return


if __name__ == "__main__":
    app.run()
