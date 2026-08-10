import marimo

__generated_with = "0.20.4"
app = marimo.App(
    width="medium",
    layout_file="layouts/GermanMimsFitting.grid.json",
)


@app.cell
def imports():
    #import all needed modules
    import marimo as mo
    import numpy as np
    from importlib import reload
    import codechembook.quickPlots as qp
    from scipy.ndimage import gaussian_filter1d
    from plotly.subplots import make_subplots
    from lmfit import Model
    import plotly.graph_objects as go
    from datetime import datetime # this is important for naming files

    #import modules specifically written for this 
    import pulsePhysicsGermany as pp
    import mdHelper as mh
    import csv
    import glob


    reload(pp)
    reload(mh)
    return Model, csv, make_subplots, mh, mo, np, pp, reload


@app.cell
def constants():
    # constants <- should I split this up for different functions? Shoudl I make UI elements for each one? That might make it more explicit what is happening. 
    R_np_md = 15e-10  #radius in m for the particle in teh MD simluations
    bulk_toluene_density = 0.0453e30  # number of hydrogens per cubic m (used to be angtroms)
    bulk_dodecanethiol_density = 0 # the density of hydrogens coming from the ligand at r = infinity
    md_r_step = 0.1e-10 # units = m, 
    O_S_difference = 0.08e-9 # units of m
    return O_S_difference, R_np_md, bulk_toluene_density, md_r_step


@app.cell
def process_np_sizes(make_subplots, np):
    # read np size data



    r_Hs = np.genfromtxt("../TEMdata/Feret_All_PdSC12_fit_report.csv", unpack = True, delimiter = ",", usecols = [0], skip_header=1)
    r_Ds = np.genfromtxt("../TEMdata/Feret_All_dPdSC12_lognormal.csv", unpack = True, delimiter = ",", usecols = [0], skip_header=1)

    r_plot = make_subplots(rows = 2, cols = 1)
    r_plot.add_histogram(x = r_Hs, row = 1, col = 1, nbinsx=20)
    r_plot.add_histogram(x = r_Ds, row = 2, col = 1, nbinsx=20)
    r_plot.update_xaxes(range = [0, 10])
    r_plot.show()
    return r_Ds, r_Hs


@app.cell
def process_md_results(
    O_S_difference,
    R_np_md,
    bulk_toluene_density,
    make_subplots,
    md_r_step,
    mh,
    np,
):
    # read md data files

    # high density files
    hd_toluene_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/toluene_distances_high_density.csv",
                                       delimiter = ",")*1e-10 + O_S_difference # convert to m correct for O S distances
    hd_dodecanethiol_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/dodecanol_distances_high_density.csv",
                                       delimiter = ",")*1e-10 + O_S_difference # convert to m correct for O S distances

    # low density files
    ld_toluene_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/toluene_distances_low_density.csv",
                                       delimiter = ",")*1e-10 + O_S_difference # convert to m correct for O S distances
    ld_dodecanethiol_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/dodecanol_distances_low_density.csv",
                                       delimiter = ",")*1e-10 +  O_S_difference# convert to m correct for O S distances

    # process md data files

    # high density
    #toluene
    hd_toluene_profile = mh.NPDensityProfile(hd_toluene_distances, R_np_md, bin_width=md_r_step)
    hd_toluene_profile.prepare_blended_density(rho_bulk = bulk_toluene_density)

    #ligand
    hd_dodecanethiol_profile = mh.NPDensityProfile(hd_dodecanethiol_distances, R_np_md, old_bins = hd_toluene_profile.bin_centers)


    #low density
    #toluene
    ld_toluene_profile = mh.NPDensityProfile(ld_toluene_distances, R_np_md, bin_width=md_r_step)
    ld_toluene_profile.prepare_blended_density(rho_bulk = bulk_toluene_density)

    #ligand
    ld_dodecanethiol_profile = mh.NPDensityProfile(ld_dodecanethiol_distances, R_np_md, old_bins = ld_toluene_profile.bin_centers)



    #
    # Plot
    #
    densityPLot  = make_subplots(rows =2 , cols = 1)

    densityPLot.add_scatter(x = np.array(hd_dodecanethiol_profile.bin_centers)*1e9, y = hd_dodecanethiol_profile.densities/6.022e23*1e-6, 
                                  line = dict(color = "#BB6F35", width = 0), fill = "tozeroy", name = "ligand", fillcolor="#BB6F35")

    #densityPLot.add_scatter(x = np.array(hd_toluene_profile.bin_centers), y = hd_toluene_profile.densities, line = dict(color = "#448CC4"), fill = "tozeroy", name = "solvent")


    densityPLot.add_scatter(x = np.array(hd_toluene_profile.bin_centers)*1e9, y = (hd_toluene_profile.densities+hd_dodecanethiol_profile.densities)/6.022e23*1e-6, line = dict(color = "#448CC4"), fill = "tonexty",)# fillcolor="#448CC4")

    densityPLot.add_scatter(x = np.array(hd_toluene_profile.bin_centers)*1e9, y = (hd_toluene_profile.densities+hd_dodecanethiol_profile.densities)/6.022e23*1e-6, line = dict(color = "black"),)# fillcolor="#448CC4")





    densityPLot.add_scatter(x = np.array(ld_dodecanethiol_profile.bin_centers)*1e9, y = ld_dodecanethiol_profile.densities//6.022e23*1e-6, 
                                  line = dict(color = "#BB6F35", width = 0), fill = "tozeroy", fillcolor="#BB6F35", name = "ligand", row = 2, col = 1)

    #densityPLot.add_scatter(x = np.array(ld_toluene_profile.bin_centers), y = ld_toluene_profile.densities, line = dict(color = "#448CC4"), fill = "tozeroy", name = "solvent", row = 2, col = 1)

    densityPLot.add_scatter(x = np.array(ld_toluene_profile.bin_centers)*1e9, y = (ld_dodecanethiol_profile.densities+ld_toluene_profile.densities)/6.022e23*1e-6, line = dict(color = "#448CC4"), fill = "tonexty", row = 2, col = 1)# fillcolor="#448CC4")

    densityPLot.add_scatter(x = np.array(ld_toluene_profile.bin_centers)*1e9, y = (ld_dodecanethiol_profile.densities+ld_toluene_profile.densities)/6.022e23*1e-6, line = dict(color = "black"),row = 2, col = 1)# fillcolor="#448CC4")

    densityPLot.update_layout(template = "simple_white",
                              #showlegend = False, 
                              width = 300*1.7, height = 300*2.8,
                             showlegend = False)
    densityPLot.update_xaxes(range = [0, 3.4], title = "radial distance /nm")
    densityPLot.update_yaxes(title = "H density / mol⋅mL<sup>-1</sup>", range = [0, 5.5e29/6.022e23*1e-6], showline = False)
    densityPLot.show()
    densityPLot.write_image("densityPlotPub.svg")
    return (
        hd_dodecanethiol_profile,
        hd_toluene_profile,
        ld_dodecanethiol_profile,
        ld_toluene_profile,
    )


@app.cell
def process_exp_spectra(make_subplots, np, nu_n_guess):
    # read spectral data files





    dodecaneH_tolueneH_spectra =  np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)
    dodecaneD_tolueneH_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter="," ,unpack = True, skip_header = 1)
    dodecaneH_tolueneD_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter=",", unpack = True, skip_header = 1)

    # get the number of values in the first two spectra, which is critical for later use in fitting and simulation
    #len1 = len(dodecaneH_tolueneH_spectra[0])
    #len2 = len(dodecaneD_tolueneH_spectra[0])

    # convert from MHz to Hz 
    dodecaneH_tolueneH_spectra[0] = dodecaneH_tolueneH_spectra[0]*1e6 
    dodecaneD_tolueneH_spectra[0] = dodecaneD_tolueneH_spectra[0]*1e6 
    dodecaneH_tolueneD_spectra[0] = dodecaneH_tolueneD_spectra[0]*1e6 


    # baseline correct and normalize...
    normalize = True
    if normalize:
        dodecaneH_tolueneH_spectra[1] = (dodecaneH_tolueneH_spectra[1] - np.mean(dodecaneH_tolueneH_spectra[1][:10]))/np.max((dodecaneH_tolueneH_spectra[1] - np.mean(dodecaneH_tolueneH_spectra[1][:10])))

        dodecaneD_tolueneH_spectra[1] = (dodecaneD_tolueneH_spectra[1] - np.mean(dodecaneD_tolueneH_spectra[1][:10]))/np.max((dodecaneD_tolueneH_spectra[1] - np.mean(dodecaneD_tolueneH_spectra[1][:10])))

        dodecaneH_tolueneD_spectra[1] = (dodecaneH_tolueneD_spectra[1] - np.mean(dodecaneH_tolueneD_spectra[1][:10]))/np.max((dodecaneH_tolueneD_spectra[1] - np.mean(dodecaneH_tolueneD_spectra[1][:10])))

    # can also crop out x-values...
    mask_limit = 5e6
    mask_data = True
    if mask_data:
        mask = (dodecaneH_tolueneH_spectra[0] >= nu_n_guess - mask_limit) & (dodecaneH_tolueneH_spectra[0] <= nu_n_guess + mask_limit)
        trimmed_HH_frequency = dodecaneH_tolueneH_spectra[0][mask]
        trimmed_HH_intensity = dodecaneH_tolueneH_spectra[1][mask]

        mask = (dodecaneD_tolueneH_spectra[0] >= nu_n_guess - mask_limit) & (dodecaneD_tolueneH_spectra[0] <= nu_n_guess + mask_limit)
        trimmed_DH_frequency = dodecaneD_tolueneH_spectra[0][mask]
        trimmed_DH_intensity = dodecaneD_tolueneH_spectra[1][mask]

        mask = (dodecaneH_tolueneD_spectra[0] >= nu_n_guess - mask_limit) & (dodecaneH_tolueneD_spectra[0] <= nu_n_guess + mask_limit)
        trimmed_HD_frequency = dodecaneH_tolueneD_spectra[0][mask]
        trimmed_HD_intensity = dodecaneH_tolueneD_spectra[1][mask]

        dodecaneH_tolueneH_spectra = [trimmed_HH_frequency, trimmed_HH_intensity] 
        dodecaneD_tolueneH_spectra = [trimmed_DH_frequency, trimmed_DH_intensity]
        dodecaneH_tolueneD_spectra = [trimmed_HD_frequency, trimmed_HD_intensity]


    # then when we fit, we can just use c and I of 1

    print(nu_n_guess)
    expPlot = make_subplots(rows = 3, cols = 1)
    for _i, _spectrum in enumerate([dodecaneH_tolueneH_spectra, dodecaneD_tolueneH_spectra, dodecaneH_tolueneD_spectra], start = 1):
        expPlot.add_scatter(x = _spectrum[0]-nu_n_guess, y = _spectrum[1], row = _i, col = 1)
    expPlot.update_layout(template = "simple_white")
    expPlot.show()
    return (
        dodecaneD_tolueneH_spectra,
        dodecaneH_tolueneD_spectra,
        dodecaneH_tolueneH_spectra,
    )


@app.cell
def simulation_conditions():
    '''
    - flag for not having contact
    - flag for not having diplolar
    - flag for how to handle np sizes
    - etc etc.
    '''
    return


@app.cell
def simulation(
    D_np,
    H_np,
    I1_guess,
    I2_guess,
    I3_guess,
    I_scale_guess,
    R_max_guess,
    R_np_md,
    T_m,
    blindspot_offset_guess,
    c1_guess,
    c2_guess,
    c3_guess,
    contact_scale_guess,
    dipolar_np_scale_guess,
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
    '''
    - alter the routine so it uses convolution
    - alter it so it calcluates the contact and dipolar separately, and then adds them together.  Or at least calculates them in different if-then loops
    - 
    '''

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
        1*I_scale_guess, # unitless, a scalar that adjust intenisty, but should default to 1.
        1*ell_guess, # unit of m, characteristic spatial distance of the decay
        g_e, # electronic g-factor

        #properties of the EPR experiment
        tau, # unit = s, delay time of the pulse sequence
        T_m, # unit = s, coherence lifetime

        # properties of the simulation
        R_np_md, 
        1.5*R_max_guess, # unit = m, how far out in the simulation to look 
        target_dr, # unit = m, 
        dipolar_scale_guess, 
        contact_scale_guess, 
        blindspot_offset_guess, 
        dipolar_np_scale_guess, 

        # properties of the spectra
        shape,
        eta_guess,
        width_scale_guess,
        )
    return (sim_spectra,)


@app.cell
def sim_compare(
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

    sim_compare.update_xaxes(
        #range = [-5e6, 5e6],
        )
    sim_compare.update_layout(template = "simple_white", height = 700)
    sim_compare.show()
    return (sim_compare,)


@app.cell
def mims_function(np, pp):
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
        blindspot_offset = 0,
        dipolar_np_scale = 1,

        # properties of the spectra
        shape = "Lorentzian",
        eta = 0,
        width_scale = 1,
        ):


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
            for r, weight_r in zip(rs[0], rs[1]): # unpack the rs  
                if weight_r > 0: # only need to calculate if there are nanoparticles of this size
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
                            blindspot_offset,
                            dipolar_np_scale,

                            shape, 
                            eta, 
                            width_scale
                            )
                    else:
                        y_solvent = np.zeros_like(freq)

                    # Optimization: Only calculate ligand if it actually contributes
                    if l_weight > 0.00:
                        #print("calculating ligand")
                        y_ligand = pp.simulate_nanoparticle_mims_spectrum(freq, r, ligand_density_obj, rho_bulk_ligand, nu_n, g_n, I_scale, ell, g_e, tau, T_m, R_np_md, R_max, target_dr, dipolar_scale, contact_scale, blindspot_offset, dipolar_np_scale, shape, eta, width_scale)
                    else:
                        y_ligand = np.zeros_like(freq)

                    # ACCUMULATE: Add this size bin's weighted contribution to the total
                    sample_total_raw += weight_r * (s_weight * y_solvent + l_weight * y_ligand)

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
def fitting(
    D_np,
    H_np,
    I1_guess,
    I2_guess,
    I3_guess,
    I_scale_guess,
    Model,
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
    eta_guess,
    g_e,
    g_n,
    len1,
    len2,
    ligand_density_obj,
    np,
    pp,
    reload,
    rho_bulk_ligand,
    rho_bulk_solvent,
    shape,
    solvent_density_obj,
    target_dr,
    tau,
    triple_mims_simulation,
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





    for _vary in np.linspace(0, .1, 1):
        print(_vary)

        triple_spectra_model_params = triple_spectra_model.make_params()

        triple_spectra_model_params.add_many(
            ("I1", I1_guess*0.91, True, 0, None), 
            ("I2", I2_guess*0.55, True, 0, None), 
            ("I3", I3_guess*1, True, 0, None), # final intensity for each spectra
            ("c1", c1_guess, False, None, None), 
            ("c2", c2_guess, False, None, None),
            ("c3", c3_guess, False, None, None), # final background for each spectra   

            ("nu_n", 49963664, False, None, None), # in Hz

            ("I_scale", I_scale_guess*9.6, True, 0, None), # for non-B type, try 41
            ("ell", 0.35e-9, True, 0, None), # bring in as nm #1.5*ell_guess

            ('R_max', 0.01, False, 0, None), # bring in as nm

            ('dipolar_scale', dipolar_scale_guess, False, 0, None),
            ('contact_scale', contact_scale_guess, False, 0, None),  
            ('blindspot_offset', 3.98e-4 , False, 0, 1), 
            ('dipolar_np_scale', 0, False, 0, 1), # repurposed for adjusting center frequency of second peak
            ("eta", eta_guess, False, 0, None), # repurposed for adjusting center frequency of third peak
            ("width_scale", 1.13, False, 0, None), # or width_scale_guess  will just scale the width calculated from lifetime
            ) 

        # add in parameters for fitting: 
        # R_scale <- parameter to control how far from the nanoparticle surface to adjust when calculating dipolar coupling
        # blind_spot_offset <- how imperfect the blindspot is. 
        # fi <- central resonance for the spectrum.  OR we can correct the spectra to 0? 

        # remove parameters for fitting: 


        try:
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

            #write_fit_results(f"./playingAround/UniformSphereCorrection.csv", triple_spetra_fit)
            print(triple_spetra_fit.fit_report())
            print(f"{_nm}: {triple_spetra_fit.aic}")
            #write_fit_results(f"R_max_grid_search_{_nm}nm.csv", triple_spetra_fit)
        except: 
            pass
    return (triple_spetra_fit,)


@app.cell
def fit_report(csv):
    def write_fit_results(csv_filename, result):

        with open(csv_filename, 'w', newline='') as f:
            writer = csv.writer(f)

            # 1. Write the Fit Statistics
            writer.writerow(['FIT STATISTICS'])
            writer.writerow(['Data points', getattr(result, 'ndata', 'N/A')])
            writer.writerow(['Variables', getattr(result, 'nvarys', 'N/A')])
            writer.writerow(['Chi-square', getattr(result, 'chisqr', 'N/A')])
            writer.writerow(['Reduced Chi-square', getattr(result, 'redchi', 'N/A')])
            writer.writerow(['Akaike info crit (AIC)', getattr(result, 'aic', 'N/A')])
            writer.writerow(['Bayesian info crit (BIC)', getattr(result, 'bic', 'N/A')])
            if hasattr(result, 'rsquared'):
                writer.writerow(['R-squared', result.rsquared])

            writer.writerow([]) # Empty row for visual separation

            # 2. Write the Fit Parameters Header
            writer.writerow(['FIT PARAMETERS'])
            headers = ['Name', 'Value', 'Standard Error', 'Initial Value', 'Min', 'Max', 'Vary', 'Expression']
            writer.writerow(headers)

            # 3. Extract and Write Parameter Data Rows
            for name, param in result.params.items():
                row = [
                    name,
                    param.value,
                    param.stderr if param.stderr is not None else 'N/A',
                    param.init_value,
                    param.min,
                    param.max,
                    param.vary,
                    param.expr
                ]
                writer.writerow(row)

        return print(f"Fit report saved successfully to {csv_filename}")

    return


@app.cell
def _(
    LDSH_color,
    LHSD_color,
    LHSH_color,
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
    fit_nu = triple_spetra_fit.params["nu_n"].value


    fitplot = make_subplots(rows = 4, cols = 1, subplot_titles= ["solvent-H + ligand-H", "solvent-H + ligand-D", "solvent-D + ligand-H", "comparison of fits"])
    fitplot.update_layout(template = "simple_white", width = 500, height = 800)

    fitplot.add_scatter(x = dodecaneH_tolueneH_spectra[0] - fit_nu, y= dodecaneH_tolueneH_spectra[1], line = dict(color = LHSH_color, width = 5), showlegend = False, col = 1, row = 1)
    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= triple_spetra_fit.best_fit[:len1], line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= dodecaneD_tolueneH_spectra[1], line = dict(color = LDSH_color, width = 5), showlegend = False, col = 1, row = 2)
    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "black"), showlegend = False, col = 1, row = 2)

    fitplot.add_scatter(x = dodecaneH_tolueneD_spectra[0] - fit_nu, y= dodecaneH_tolueneD_spectra[1], line = dict(color = LHSD_color, width = 5), showlegend = False, col = 1, row = 3)
    fitplot.add_scatter(x = dodecaneH_tolueneD_spectra[0] - fit_nu, y= triple_spetra_fit.best_fit[len1+len2:], line = dict(color = "black"), showlegend = False, col = 1, row = 3)

    #add normalized spectra
    fitplot.add_scatter(x = dodecaneH_tolueneH_spectra[0] - fit_nu, y= (triple_spetra_fit.best_fit[:len1] - triple_spetra_fit.params["c1"].value)/np.max((triple_spetra_fit.best_fit[:len1] - triple_spetra_fit.params["c1"].value)), line = dict(color = LHSH_color, width = 3), showlegend = False, col = 1, row = 4)
    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= (triple_spetra_fit.best_fit[len1:len1+len2] - triple_spetra_fit.params["c2"].value)/np.max(triple_spetra_fit.best_fit[len1:len1+len2] - triple_spetra_fit.params["c2"].value), line = dict(color = LDSH_color, width = 3), showlegend = False, col = 1, row = 4)
    fitplot.add_scatter(x = dodecaneH_tolueneD_spectra[0] - fit_nu, y= (triple_spetra_fit.best_fit[len1+len2:] - triple_spetra_fit.params["c3"].value)/np.max(triple_spetra_fit.best_fit[len1+len2:] - triple_spetra_fit.params["c3"].value), line = dict(color = LHSD_color, width = 3), showlegend = False, col = 1, row = 4)


    fitplot.update_xaxes(title = "frequency /MHz", row = 4, col = 1)
    fitplot.show()
    return (fit_nu,)


@app.cell
def _():
    pub_dpi = 300
    main_width = 1.7 *pub_dpi#inch *300 dpi
    main_height = 3.6 *pub_dpi#inch * 300dpi
    main_margin = dict(l = 0, r = 0, t = 0, b = 0)
    LHSH_color = "#715472"
    LHSD_color = "#B66C34"
    LDSH_color = "#448CC4"
    return (
        LDSH_color,
        LHSD_color,
        LHSH_color,
        main_height,
        main_margin,
        main_width,
    )


@app.cell
def fitting_pub_figure(
    LDSH_color,
    LHSD_color,
    LHSH_color,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    fit_nu,
    len1,
    len2,
    main_height,
    main_margin,
    main_width,
    make_subplots,
    triple_spetra_fit,
):



    mainPub = make_subplots(rows = 3, cols = 1, )

    mainPub.update_layout(
        template = "simple_white",
        width = main_width, height = main_height,
        margin = main_margin,
        font = dict(size = 30)
    )

    mainPub.update_xaxes(range = [-4e6, 4e6], linewidth = 4, tickwidth = 4, ticklen = 12)
    mainPub.update_yaxes(range = [-0.2, 1.09], dtick = 1, ticks = "", showticklabels = False, showline = False,  )


    mainPub.add_scatter(x = dodecaneH_tolueneH_spectra[0] - fit_nu, y= dodecaneH_tolueneH_spectra[1] - triple_spetra_fit.best_fit[:len1], line = dict(color = "lightgrey", width = 4), showlegend = False, col = 1, row = 1)
    mainPub.add_scatter(x = dodecaneH_tolueneH_spectra[0] - fit_nu, y= dodecaneH_tolueneH_spectra[1], line = dict(color = LHSH_color, width = 10), showlegend = False, col = 1, row = 1)
    mainPub.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= triple_spetra_fit.best_fit[:len1], line = dict(color = "white", width = 8), showlegend = False, col = 1, row = 1)
    mainPub.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= triple_spetra_fit.best_fit[:len1], line = dict(color = "black", width = 4), showlegend = False, col = 1, row = 1)


    mainPub.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= dodecaneD_tolueneH_spectra[1] - triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "lightgrey", width = 4), showlegend = False, col = 1, row = 2)
    mainPub.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= dodecaneD_tolueneH_spectra[1], line = dict(color = LDSH_color, width = 10), showlegend = False, col = 1, row = 2)
    mainPub.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "white", width = 8), showlegend = False, col = 1, row = 2)
    mainPub.add_scatter(x = dodecaneD_tolueneH_spectra[0] - fit_nu, y= triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "black", width = 4), showlegend = False, col = 1, row = 2)


    mainPub.add_scatter(x = dodecaneH_tolueneH_spectra[0] - fit_nu, y= dodecaneH_tolueneH_spectra[1] - triple_spetra_fit.best_fit[:len1], line = dict(color = "lightgrey", width = 4), showlegend = False, col = 1, row = 3)
    mainPub.add_scatter(x = dodecaneH_tolueneH_spectra[0] - fit_nu, y= dodecaneH_tolueneH_spectra[1] - triple_spetra_fit.best_fit[:len1], line = dict(color = "lightgrey", width = 4), showlegend = False, col = 1, row = 1)
    mainPub.add_scatter(x = dodecaneH_tolueneD_spectra[0] - fit_nu, y= dodecaneH_tolueneD_spectra[1], line = dict(color = LHSD_color, width = 10), showlegend = False, col = 1, row = 3)
    mainPub.add_scatter(x = dodecaneH_tolueneD_spectra[0] - fit_nu, y= triple_spetra_fit.best_fit[len1+len2:], line = dict(color = "white", width = 8), showlegend = False, col = 1, row = 3)
    mainPub.add_scatter(x = dodecaneH_tolueneD_spectra[0] - fit_nu, y= triple_spetra_fit.best_fit[len1+len2:], line = dict(color = "black", width = 4), showlegend = False, col = 1, row = 3)
    mainPub.update_xaxes(row = 3, col = 1, title = "A /Hz")


    mainPub.show()
    return


@app.cell
def ld_binned_contact_dipolar():
    return


@app.cell
def initialize_ui(mo):
    md_selector = mo.ui.radio(label = "MD ligand coverage", options = ["low density", "high density"], value = "high density")
    size_handling_selector = mo.ui.radio(label = "treat size population", options = ["indivdually", "binned", "single mean"], value = "binned")
    coupling_handling_selector = mo.ui.radio(label = "coupling to include", options = ["only dipolar", "only contact", "dipolar + contact"], value = "dipolar + contact")

    normalization_selector = mo.ui.checkbox(label = "normalize spectra?", value = True)

    I1_guess_number = mo.ui.number(label = "I1_guess", value = 1)
    I2_guess_number = mo.ui.number(label = "I2_guess", value = 1)
    I3_guess_number = mo.ui.number(label = "I3_guess", value = 1)

    c1_guess_number = mo.ui.number(label = "c1_guess", value = 0)
    c2_guess_number = mo.ui.number(label = "c2_guess", value = 0)
    c3_guess_number = mo.ui.number(label = "c3_guess", value = 0)
    return (
        coupling_handling_selector,
        md_selector,
        normalization_selector,
        size_handling_selector,
    )


@app.cell
def ui(
    coupling_handling_selector,
    md_selector,
    mo,
    normalization_selector,
    size_handling_selector,
):
    mo.md(f"""
    {md_selector}\n
    {size_handling_selector}\n
    {coupling_handling_selector}\n
    {normalization_selector}\n
    """)
    return


@app.cell
def _(mo, sim_compare):
    mo.ui.plotly(sim_compare)
    return


@app.cell
def _():
    # for correcting resonance positions to 0
    field_1 = 11745
    field_2 = 11729
    nu_G = 0.004262
    nu_n_guess = field_2 * nu_G*1e6  # this is for centered spectra, but if not, use 50e6 # units = s
    return (nu_n_guess,)


@app.cell
def sim_setup(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    hd_dodecanethiol_profile,
    hd_toluene_profile,
    ld_dodecanethiol_profile,
    ld_toluene_profile,
    md_selector,
    normalization_selector,
    np,
    pp,
    r_Ds,
    r_Hs,
    reload,
    size_handling_selector,
):
    reload(pp)

    # VALUES for parameters that are not fit
    len1 = len(dodecaneH_tolueneH_spectra[0])
    len2 = len(dodecaneD_tolueneH_spectra[0])

    # set up the experimental sizes. all_sizes is a flag to use distributions, rather than mean sizes
    if size_handling_selector.value == "binned":
        H_np = [np.array([1.272975688, 1.405767063, 1.538558438, 1.671349814, 1.804141189, 1.936932564, 2.06972394, 2.202515315, 2.33530669, 2.468098066, 2.600889441, 2.733680816, 2.866472192, 2.999263567, 3.132054942, 3.264846318, 3.397637693, 3.530429068, 3.663220444, 3.796011819, 3.928803194, 4.06159457, 4.194385945, 4.32717732, 4.459968696, 4.592760071, 4.725551446, 4.858342822, 4.991134197, 5.123925572, 5.256716947, 5.389508323, 5.522299698, 5.655091073, 5.787882449, 5.920673824, 6.053465199, 6.186256575,])*1e-9/2, [4, 1, 2, 1, 0, 4, 1, 11, 17, 49, 82, 83, 81, 59, 55, 36, 29, 20, 25, 8, 3, 6, 10, 5, 1, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,]]
        D_np = [np.array([1.541959469, 1.728658408, 1.915357347, 2.102056286, 2.288755225, 2.475454163, 2.662153102, 2.848852041, 3.03555098, 3.222249919, 3.408948858, 3.595647796, 3.782346735, 3.969045674, 4.155744613, 4.342443552, 4.52914249, 4.715841429, 4.902540368, 5.089239307,])*1e-9/2, [14, 32, 65, 55, 61, 67, 46, 40, 26, 21, 11, 6, 8, 1, 2, 2, 0, 0, 0, 1,]]#

    elif size_handling_selector.value == "indivdually":
        H_np = [r_Hs*1e-9/2, np.ones_like(r_Hs)]
        D_np = [r_Ds*1e-9/2, np.ones_like(r_Ds)]
    elif size_handling_selector.value == "single mean": 
        H_np = [[2.87e-9/2],[1.0]] # bin center in m, frequency
        D_np = [[2.42e-9/2],[1.0]] # bin center in m, frequency
    else:
        print("make a choice for size handling that exists")


    if md_selector.value == "high density":
        solvent_density_obj = hd_toluene_profile
        ligand_density_obj = hd_dodecanethiol_profile
    elif md_selector.value == "low density": 
        solvent_density_obj = ld_toluene_profile
        ligand_density_obj = ld_dodecanethiol_profile
    else:
        print("make a choice for ligand converage that exists")

    rho_bulk_solvent = 4.528e28 # number of hydrogens per cubic m (used to be cubic angstrom)
    rho_bulk_ligand = 0

    g_n = 5.5856946893 # from scipy constants

    g_e = 2.0 # g-value for the nanoparticles

    tau = 140e-9 # in seconds. From measurements
    T_m = 674e-9 # in seconds. From measurements

    target_dr = 0.1e-9

    shape = "L"

    # GUESSES for parameters to be fit
    norm = True # set equal to true, if the spectra have been purposefully normalized
    if normalization_selector.value:
        I1_guess = 1
        I2_guess = 1
        I3_guess = 1
        c1_guess = 0
        c2_guess = 0
        c3_guess = 0
    else:
        I1_guess = (np.max(dodecaneH_tolueneH_spectra[1]) - np.min(dodecaneH_tolueneH_spectra[1]))
        I2_guess = (np.max(dodecaneD_tolueneH_spectra[1]) - np.min(dodecaneD_tolueneH_spectra[1]))
        I3_guess = (np.max(dodecaneH_tolueneD_spectra[1]) - np.min(dodecaneH_tolueneD_spectra[1]))

        c1_guess = np.mean(dodecaneH_tolueneH_spectra[1][0:10])
        c2_guess = np.mean(dodecaneD_tolueneH_spectra[1][0:10])
        c3_guess = np.mean(dodecaneH_tolueneD_spectra[1][0:10])






    I_scale_guess = 1 # unitless, start with a guess that assumes our simple treatment is correct

    ell_guess = 4.6232e-10 # unit = m the decay constant for wavefunction. Guess comes from some paper Kirsten read.

    R_max_guess = 15e-9 # m, the max radius we look out to.  Seems like we get good fits, no matter what

    dipolar_scale_guess = 1 # unitless
    contact_scale_guess = 1 # unitless

    blindspot_offset_guess = 0.05 # unitless, bounds should be 0 to 1
    dipolar_np_scale_guess = 0.685#/1.04

    eta_guess = 0 # unitless
    width_scale_guess = 0.9 # unitless
    return (
        D_np,
        H_np,
        I1_guess,
        I2_guess,
        I3_guess,
        I_scale_guess,
        R_max_guess,
        T_m,
        blindspot_offset_guess,
        c1_guess,
        c2_guess,
        c3_guess,
        contact_scale_guess,
        dipolar_np_scale_guess,
        dipolar_scale_guess,
        ell_guess,
        eta_guess,
        g_e,
        g_n,
        len1,
        len2,
        ligand_density_obj,
        rho_bulk_ligand,
        rho_bulk_solvent,
        shape,
        solvent_density_obj,
        target_dr,
        tau,
        width_scale_guess,
    )


@app.cell
def specify_files():
    # change the meaning of R_max, so that it is tied to a particular threshold for A.  But what shoudl that be? 
    # then see if we can take reasonable guesses at things...
    # double check how dipolar coupling is calculated.
    return


@app.cell
def _():
    return


@app.cell
def _(R_np_md, T_m, make_subplots, np, r_Hs, triple_spetra_fit):
    import scipy.constants as const


    def find_R_max(
        R_np_exp, ell, I_scale, T_m, width_scale, contact_scale, dipolar_scale, 
        mu_0, g_e, beta_e, g_n, beta_n, Planck, R_max_cuttoff
    ):
        # 1. Target threshold based on linewidth
        gamma_hz = (1.0 / (2.0 * np.pi * T_m)) * width_scale  
        threshold = gamma_hz * R_max_cuttoff  # e.g. 5% of homogeneous line width

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
            T_factor_probe = dipolar_scale * (4.0 * np.pi)**-1 * common_scalar * (r**-3)# *(1 - (R_np_exp / r)**2)**-3 

            # Max coupling magnitude at theta = 0 (dipolar is 2*T)
            shift = A_contact_probe + 2.0 * T_factor_probe

            r += dr

        # 5. Return R_max with a fallback minimum bound
        if r > R_np_exp + 1e-9:
            return r
        else:
            return R_np_exp + 1e-9


    def calculate_hyperfine_vs_r(
        R_np_exp, R_np_md, ell, I_scale, T_m, 
        width_scale=1, contact_scale=1, dipolar_scale=1, 
        target_dr=0.1e-10, R_max_cuttoff=0.05
    ):
        """
        Calculates R_max using the provided logic, then extracts the 
        contact and dipolar hyperfine contributions across the spatial grid.
        """
        # --- Physical Constants ---
        mu_0 = const.mu_0
        g_e = 2.002319 
        beta_e = const.physical_constants['Bohr magneton'][0]
        g_n = 5.585 # Placeholder for 1H
        beta_n = const.physical_constants['nuclear magneton'][0]
        Planck = const.h

        # --- 1. Find R_max ---
        R_max = find_R_max(
            R_np_exp, ell, I_scale, T_m, width_scale, contact_scale, dipolar_scale, 
            mu_0, g_e, beta_e, g_n, beta_n, Planck, R_max_cuttoff
        )
        print(f"Calculated R_max: {R_max * 1e9:.2f} nm")

        # --- 2. Generate the Spatial Grid ---
        n_shells = int(np.ceil((R_max - R_np_md) / target_dr))
        dr = (R_max - R_np_md) / n_shells
        r_axis = np.linspace(R_np_exp + dr/2, R_np_exp + n_shells*dr - dr/2, n_shells)

        # --- 3. Vectorized Calculations ---
        common_scalar = mu_0 * g_e * beta_e * g_n * beta_n / Planck
        f = (1 + 2 * R_np_exp / (3 * ell))**-1 
        I0 = I_scale * f * np.sqrt(2 * f / ell)

        A_contact = contact_scale * (2/3) * common_scalar * (I0 * np.exp(-(r_axis - R_np_exp) / ell) / r_axis)**2

        # Updated T_factor with spherical shell correction
        T_factor = dipolar_scale * (4 * np.pi)**-1 * common_scalar * (r_axis**-3) #* (1 - (R_np_exp / r_axis)**2)**-3

        return r_axis, A_contact, T_factor

    # ==========================================
    # Example Execution and Plotting
    # ==========================================

    params = {
        'R_np_exp': (2.87e-9)/2,  # the mean size for the all H sample
        'R_np_md': R_np_md,
        'ell': triple_spetra_fit.params["ell"].value,       
        'I_scale': triple_spetra_fit.params["I_scale"].value,
        'T_m': T_m,         
        'width_scale': 1,
        'R_max_cuttoff': 0.01  # Represents 5% of linewidth threshold
    }

    # ==========================================
    # Example Execution and Plotting
    # ==========================================


    A_contact_color = "darkcyan"
    T_factor_color = "darkred"

    AcomparePlot = make_subplots(rows = 2, cols = 2,)

    line_syles = ["solid", "dot"]
    A_range = [0.01, 0.01]
    A_total_range = [0.01, 0.01]
    size_list = [np.min(r_Hs*1e-9/2), np.max(r_Hs*1e-9/2)]

    for _i, _r in enumerate(size_list): # this allows cycling through multiple values for r_exp
        params["R_np_exp"] = _r
        r, A_contact, T_factor = calculate_hyperfine_vs_r(**params)

        # Convert units for visualization (meters -> nanometers, Hz -> MHz)
        r_nm = (r-_r) * 1e9
        A_contact_MHz = A_contact / 1e6
        T_factor_MHz = T_factor / 1e6


        if np.log10(np.min([np.min(A_contact_MHz), np.min(T_factor_MHz)])*0.5) < A_range[0]:
            A_range[0] = np.log10(np.min([np.min(A_contact_MHz), np.min(T_factor_MHz)])*0.5)
        if np.log10(np.max([np.max(A_contact_MHz), np.max(T_factor_MHz)])*2) > A_range[1]:
            A_range[1] = np.log10(np.max([np.max(A_contact_MHz), np.max(T_factor_MHz)])*2)
        #A_range = [np.log10(np.min([np.min(A_contact_MHz), np.min(T_factor_MHz)])), np.log10(np.max([np.max(A_contact_MHz), np.max(T_factor_MHz)])*2)]


        if np.log10(np.min(A_contact_MHz + T_factor_MHz)) < A_total_range[0]:
            A_total_range[0] = np.log10(np.min(A_contact_MHz + T_factor_MHz))
        if np.log10(np.max(A_contact_MHz+T_factor_MHz)) > A_total_range[1]:
            A_total_range[1] = np.log10(np.max(A_contact_MHz+T_factor_MHz))
        #A_total_range = [np.log10(np.min(A_contact_MHz + T_factor_MHz)), np.log10(np.max(A_contact_MHz+T_factor_MHz))]

        minorTicks = dict(ticks="outside", nticks = 10, ticklen=2, tickcolor="black", showgrid=False)



        AcomparePlot.update_xaxes(row = 1, col = 1, )#title = "r /nm")
        AcomparePlot.update_yaxes(row = 1, col = 1, type = "log", exponentformat = 'power', title = "A /Mhz", range = A_range, minor = minorTicks, nticks = 10)
        AcomparePlot.add_scatter(row = 1, col = 1, x = r_nm, y = A_contact_MHz, line = dict(color = A_contact_color, dash = line_syles[_i]))
        AcomparePlot.add_scatter(row = 1, col = 1, x = r_nm, y = T_factor_MHz, line = dict(color = T_factor_color, dash = line_syles[_i]))

        AcomparePlot.update_xaxes(row = 2, col = 1, title = "distance from nanoparticle surface /nm")
        AcomparePlot.update_yaxes(row = 2, col = 1, title = "fractional contribution to A<sub>total</sub>", range = [0, 1.01])
        AcomparePlot.add_scatter(row = 2, col = 1, x = r_nm, y = A_contact_MHz/(A_contact_MHz+T_factor_MHz), line = dict(color = A_contact_color, dash = line_syles[_i]))
        AcomparePlot.add_scatter(row = 2, col = 1, x = r_nm, y = T_factor_MHz/(A_contact_MHz+T_factor_MHz), line = dict(color = T_factor_color, dash = line_syles[_i]))



        AcomparePlot.update_xaxes(row = 1, col = 2, type = "log", range = A_total_range, minor = minorTicks)#title = "A_total /Mhz")
        AcomparePlot.update_yaxes(row = 1, col = 2, type = "log", exponentformat = 'power', minor = minorTicks, range = A_range)#title = "A /Mhz")
        AcomparePlot.add_scatter(row = 1, col = 2, x = A_contact_MHz+T_factor_MHz, y = A_contact_MHz, line = dict(color = A_contact_color, dash = line_syles[_i]))
        AcomparePlot.add_scatter(row = 1, col = 2, x = A_contact_MHz+T_factor_MHz, y = T_factor_MHz, line = dict(color = T_factor_color, dash = line_syles[_i]))

        AcomparePlot.update_xaxes(row = 2, col = 2, type = "log", title = "A<sub>total</sub> /Mhz", range = A_total_range, minor = minorTicks)
        AcomparePlot.update_yaxes(row = 2, col = 2, )#title = "fractional contribution to A_total")
        AcomparePlot.add_scatter(row = 2, col = 2, x = A_contact_MHz+T_factor_MHz, y = A_contact_MHz/(A_contact_MHz+T_factor_MHz), line = dict(color = A_contact_color, dash = line_syles[_i]))
        AcomparePlot.add_scatter(row = 2, col = 2, x = A_contact_MHz+T_factor_MHz, y = T_factor_MHz/(A_contact_MHz+T_factor_MHz), line = dict(color = T_factor_color, dash = line_syles[_i]))

        size_labels = ["smallest np", "largest np"]


        AcomparePlot.add_annotation(row = 1, col = 2, x = np.log10(np.max(A_contact_MHz+T_factor_MHz)), y = np.log10(np.max(A_contact_MHz)), text = "A<sub>contact</sub><br>"+f"{size_labels[_i]}", showarrow=False, xanchor = "right", yanchor="bottom", font = dict(color = A_contact_color))
        AcomparePlot.add_annotation(row = 1, col = 2, x = np.log10(np.max(A_contact_MHz+T_factor_MHz)), y = np.log10(np.max(T_factor_MHz)), text = "A<sub>dipolar</sub><br>"+f" {size_labels[_i]}", showarrow=False, xanchor = "center", yanchor="top", font = dict(color = T_factor_color))
        AcomparePlot.update_layout(template = "simple_white", showlegend = False)

    AcomparePlot.show()
    #AcomparePlot.write_image("A_compare_plot.png")
    return


if __name__ == "__main__":
    app.run()
