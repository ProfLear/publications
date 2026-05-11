import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _(
    A0,
    R_max,
    R_np,
    T_m,
    bulk_toluene_density,
    eta,
    g_iso,
    hd_dodecanethiol_profile,
    hd_toluene_profile,
    make_subplots,
    mf,
    np,
    nu_n,
    reload,
    rho_decay,
    sigma,
    simulate_mims_dynamic_contact_strain,
    target_dr,
    tau,
):

    reload(mf)

    print(f"solvent length: {len(hd_toluene_profile.bin_centers)}")
    print(f"ligand length : {len(hd_toluene_profile.densities)}")

    freq_axis = np.linspace(43, 57, 500)
    solventMIMS = simulate_mims_dynamic_contact_strain(
        freq_axis = freq_axis,
        density_obj= hd_toluene_profile,
        R_np = R_np, 
        A0 = A0, 
        rho_decay = rho_decay, 
        g_iso = g_iso, 
        nu_n = nu_n, 
        tau = tau, 
        T_m = T_m, 
        sigma = sigma, 
        rho_bulk = bulk_toluene_density, 
        eta = eta,
        R_max= R_max,
        target_dr = target_dr
        )

    ligandMIMS = simulate_mims_dynamic_contact_strain(
        freq_axis = freq_axis,
        density_obj= hd_dodecanethiol_profile,
        R_np = R_np, 
        A0 = A0, 
        rho_decay = rho_decay, 
        g_iso = g_iso, 
        nu_n = nu_n, 
        tau = tau, 
        T_m = T_m, 
        sigma = sigma,
        rho_bulk = 0, 
        eta = eta,
        R_max = R_max,
        target_dr= target_dr
        )

    ligandAndSolvent = (solventMIMS + ligandMIMS)/np.max(solventMIMS + ligandMIMS)
    solventOnly = (solventMIMS + 0.02*ligandMIMS)/np.max(solventMIMS + 0.02*ligandMIMS)
    ligandOnly = (0.02* solventMIMS + ligandMIMS)/np.max(0.02* solventMIMS + ligandMIMS)

    sim_plot = make_subplots(rows = 3, cols = 1)
    sim_plot.update_layout(template = "simple_white")
    sim_plot.add_scatter(x = freq_axis, 
                         y = ligandAndSolvent, 
                             name = "both", row = 1, col = 1)
    sim_plot.add_scatter(x = freq_axis, 
                         y = solventOnly, 
                             name = "solvent only", row = 2, col = 1)
    sim_plot.add_scatter(x = freq_axis, 
                         y = ligandOnly, 
                             name = "ligand only", row = 3, col = 1)

    sim_plot.show()
    return freq_axis, ligandAndSolvent, ligandOnly, solventOnly


@app.cell
def _(
    freq_axis,
    ligandAndSolvent,
    ligandOnly,
    make_subplots,
    mo,
    np,
    solventOnly,
):
    # import spectral data
    dodecaneH_tolueneH_spectra =  np.genfromtxt("/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/KristenPulse/MIMS/PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)

    dodecaneD_tolueneH_spectra = np.genfromtxt("/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/KristenPulse/MIMS/dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter="," ,unpack = True, skip_header = 1)

    dodecaneH_tolueneD_spectra = np.genfromtxt("/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/KristenPulse/MIMS/Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter="," ,unpack = True, skip_header = 1)

    spectra_Plot = make_subplots(rows = 3, cols = 1)
    spectra_Plot.update_layout(template = "simple_white")
    for _exp, _sim, _row in zip(
        [
            dodecaneH_tolueneH_spectra,
            dodecaneD_tolueneH_spectra,
            dodecaneH_tolueneD_spectra
        ], [
            (ligandAndSolvent)*(np.max(dodecaneH_tolueneH_spectra[1]) - np.min(dodecaneH_tolueneH_spectra[1])) + np.mean(dodecaneH_tolueneH_spectra[1][0:10]) ,
            solventOnly*(np.max(dodecaneD_tolueneH_spectra[1]) - np.min(dodecaneD_tolueneH_spectra[1])) + np.mean(dodecaneD_tolueneH_spectra[1][0:10]),
            ligandOnly*(np.max(dodecaneH_tolueneD_spectra[1]) - np.min(dodecaneH_tolueneD_spectra[1])) + np.mean(dodecaneH_tolueneD_spectra[1][0:10])
        ], 
        range(1,4)):


        spectra_Plot.add_scatter(x = _exp[0], y = _exp[1], 
                             name = "both", row = _row, col = 1)
        spectra_Plot.add_scatter(x = freq_axis, y = _sim, 
                             name = "solvent only", row = _row, col = 1)

    mo.md(f"{mo.ui.plotly(spectra_Plot)}")

    return (
        dodecaneD_tolueneH_spectra,
        dodecaneH_tolueneD_spectra,
        dodecaneH_tolueneH_spectra,
    )


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import mimsFunctions as mf
    from importlib import reload
    import codechembook.quickPlots as qp
    from scipy.ndimage import gaussian_filter1d
    from plotly.subplots import make_subplots
    from lmfit import Model

    return Model, make_subplots, mf, mo, np, qp, reload


@app.cell
def _(R_np, bulk_toluene_density, mf, np, qp, reload):
    # import MD data

    reload(mf)

    # first, let us process the md data. 

    # start with solvent from MD
    hd_toluene_distances = np.genfromtxt("/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/KristenPulse/MDSimulationResults/toluene_distances_high_density.csv",
                                       delimiter = ",")
    hd_toluene_profile = mf.NPDensityProfile(hd_toluene_distances, R_np, bin_width=0.1)

    # then blend to bulk
    hd_toluene_profile.prepare_blended_density(rho_bulk = bulk_toluene_density)

    print(np.max(hd_toluene_profile.bin_centers))

    # then get ligands

    hd_dodecanethiol_distances = np.genfromtxt("/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/KristenPulse/MDSimulationResults/dodecanol_distances_high_density.csv",
                                       delimiter = ",")

    #set up profile to have the same bins as the solvent
    hd_dodecanethiol_profile = mf.NPDensityProfile(hd_dodecanethiol_distances, R_np, old_bins = hd_toluene_profile.bin_centers)



    # then modify the profile for solvent, so that we get to the bulk values. 

    densityPLot  = qp.quickScatter(x = hd_dodecanethiol_profile.bin_centers, y = hd_dodecanethiol_profile.densities, output = None)
    densityPLot.add_scatter(x = hd_toluene_profile.bin_centers, y = hd_toluene_profile.densities)

    densityPLot.show()
    return hd_dodecanethiol_profile, hd_toluene_profile


@app.cell
def _(np):

    def simulate_mims_dynamic_contact_strain(freq_axis, density_obj, R_np, A0, rho_decay, 
                                             g_iso, nu_n, tau, T_m, sigma, rho_bulk, eta, R_max, target_dr = 0.2):
        """
        Simulates a Mims ENDOR spectrum where the shell thickness (dr) adjusts 
        dynamically to R_max, and density is pulled from a DensityProfile object.
        """

        # --- 1. DYNAMIC SHELL DEFINITION ---
        # Target dr is ~0.2 A. We calculate an exact n_shells to keep dr fluid.
        n_shells = int(np.ceil((R_max - R_np) / target_dr))

        # This 'dr' changes smoothly as R_max moves, preventing 'jumpy' fits.
        dr = (R_max - R_np) / n_shells 

        # Midpoints of our dynamic shells
        r_full = np.linspace(R_np + dr/2, R_max - dr/2, n_shells)

        # --- 2. SETUP GRIDS ---
        df = freq_axis[1] - freq_axis[0]
        cos_theta = np.linspace(0, 1, 150) 
        total_spectrum = np.zeros_like(freq_axis)
        damping = np.exp(-2 * tau / T_m)

        # Pre-vectorize the frequency axis for the Gaussian calculation
        f_vec = freq_axis[np.newaxis, :]

        # --- 3. INTEGRATION LOOP (Radial) ---
        for r in r_full:
            # Pull density from the object (handles MD, Blending, and Toluene Bulk)
            rho = density_obj.get_density(r, rho_bulk)

            if rho <= 0: continue

            # Physics Calculation
            T = (79.0 * (g_iso / 2.0023)) / (r**3) 
            A_iso = A0 * np.exp(-(max(0, r - R_np)) / rho_decay)
            A = A_iso + T * (3 * cos_theta**2 - 1) # Array of 150 orientations

            # Mims Filter & Base Weighting
            intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping
            w = (rho * 4 * np.pi * r**2 * dr / len(cos_theta)) * intensity

            # --- HYPERFINE STRAIN BROADENING (Contact-Only) ---
            # shell_sigma is a single value for this shell since A_iso is isotropic
            shell_sigma = np.sqrt(sigma**2 + (eta * np.abs(A_iso)**0.5)**2)  # the dependence on A_iso can be changed to reflect different disorders

            # Vectorized Gaussian Calculation for all 150 orientations at once
            A_vec = A[:, np.newaxis]
            w_vec = w[:, np.newaxis]
            mu_plus = nu_n + A_vec / 2
            mu_minus = nu_n - A_vec / 2

            # Area normalization: keeps intensity stable as the peak width changes
            norm = shell_sigma * np.sqrt(2 * np.pi)

            g_plus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2) # divide by norm to normalize
            g_minus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2) # divide by norm to normalize

            # Sum orientations into the master spectrum
            total_spectrum += np.sum(g_plus + g_minus, axis=0)

        return total_spectrum

    # set up the fit

    def function_for_iso_fit(xs, 
                             len1, len2, # to be used for indices
                             I1, I2, I3, # final intensity for each spectra
                             c1, c2, c3, # final background for each spectra
                             solvent_density_obj, rho_bulk_solvent,
                             ligand_density_obj, rho_bulk_ligand,
                             R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, 
                             eta, R_max, target_dr):
        # should sigma vary by sample? Or T_M?


        ys = np.array([]) # to hold the ys we calculate, concatenated together, so we can 

        #separate out th xs...
        freq_both =         xs[         :len1     ]
        freq_solvent_only = xs[len1     :len1+len2]
        freq_ligand_only =  xs[len1+len2:         ]

        # for each sample, calculate the contribution from the 

        for freq, I, c, s, l in zip(
            [freq_both, freq_solvent_only, freq_ligand_only],
            [I1, I2, I3],
            [c1, c2, c3],
            [1, 1, 0.01], # solvent weighting
            [1, 0.02, 1], # ligand weighting
        ):
            y_solvent = simulate_mims_dynamic_contact_strain(freq, solvent_density_obj, R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk_solvent, eta, R_max, target_dr)
            y_ligand = simulate_mims_dynamic_contact_strain(freq, ligand_density_obj, R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk_ligand, eta, R_max, target_dr)


            raw_ys = s*y_solvent + l*y_ligand

            norm_ys = raw_ys/np.max(raw_ys)

            scaled_ys = norm_ys*I + c

            ys = np.concatenate([ys, scaled_ys])


        return ys

    return function_for_iso_fit, simulate_mims_dynamic_contact_strain


@app.cell
def _(
    A0,
    R_max,
    R_np,
    T_m,
    bulk_dodecanethiol_density,
    bulk_toluene_density,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    eta,
    function_for_iso_fit,
    g_iso,
    hd_dodecanethiol_profile,
    hd_toluene_profile,
    make_subplots,
    np,
    nu_n,
    rho_decay,
    sigma,
    target_dr,
    tau,
):
    #test...

    len1 = len(dodecaneH_tolueneH_spectra[0])
    len2 = len(dodecaneD_tolueneH_spectra[0])

    sim_ys = function_for_iso_fit(xs = np.concatenate([dodecaneH_tolueneH_spectra[0], dodecaneD_tolueneH_spectra[0], dodecaneH_tolueneD_spectra[0]]), 
                             len1 = len1, len2 = len2, # to be used for indices
                             I1 = np.max(dodecaneH_tolueneH_spectra[1]) - np.min(dodecaneH_tolueneH_spectra[1]), I2 = np.max(dodecaneD_tolueneH_spectra[1]) - np.min(dodecaneD_tolueneH_spectra[1]), I3 = np.max(dodecaneH_tolueneD_spectra[1]) - np.min(dodecaneH_tolueneD_spectra[1]), # final intensity for each spectra
                             c1 = np.mean(dodecaneH_tolueneH_spectra[1][0:10]), c2 = np.mean(dodecaneD_tolueneH_spectra[1][0:10]), c3 = np.mean(dodecaneH_tolueneD_spectra[1][0:10]), # final background for each spectra
                             solvent_density_obj= hd_toluene_profile, rho_bulk_solvent = bulk_toluene_density,
                             ligand_density_obj = hd_dodecanethiol_profile, rho_bulk_ligand = bulk_dodecanethiol_density,
                             R_np = R_np, A0 = A0, rho_decay = rho_decay, g_iso = g_iso, nu_n = nu_n, tau = tau, T_m = T_m, sigma = sigma, eta = eta, R_max = R_max, target_dr= target_dr)


    simPlot = make_subplots(rows = 3, cols = 1)
    simPlot.update_layout(template = "simple_white")

    simPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneH_spectra[1], line = dict(color = "lightgrey", width = 10), col = 1, row = 1)
    simPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= sim_ys[:len1], line = dict(color = "red"), col = 1, row = 1)

    simPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneD_tolueneH_spectra[1], line = dict(color = "lightgrey", width = 10), col = 1, row = 2)
    simPlot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= sim_ys[len1:len1+len2], line = dict(color = "red"), col = 1, row = 2)

    simPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneD_spectra[1], line = dict(color = "lightgrey", width = 10), col = 1, row = 3)
    simPlot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= sim_ys[len1+len2:], line = dict(color = "red"), col = 1, row = 3)
    simPlot.show()

    #qp.quickScatter(x = np.concatenate([dodecaneH_tolueneH_spectra[0], dodecaneD_tolueneH_spectra[0], dodecaneH_tolueneD_spectra[0]]), y = sim_ys, mode = "markers")
    return len1, len2


@app.cell
def _():
    # properties of the chemical system
    R_np = 15 # 30 Angstroms
    T_m = 140 # in microseconds. 
    bulk_toluene_density = 0.0453 # number of hydrogens per cubic angstrom
    bulk_dodecanethiol_density = 0 # the density of hydrogens coming from the ligand at r = infinity
    g_iso = 1.9 # g-value for the nanoparticles
    eta = 22 # refects the heterogeneity. value of 0 means no effect

    # to be fit
    nu_n = 50 # the frequency of the hydrogens
    A0 = 6 #
    rho_decay = 1.3 # the decay constant for wavefunction
    sigma = 0.2 # to reflect the dynamics of the particles and resolution of instrument
    R_max = 50 # the max radius we look out to

    # properties of the EPR experiment
    tau=140

    #properties of fit
    target_dr = 0.1

    # to try: 
    '''
    well, I changed the threshold to be 0.1*sigma intead of 0.15* sigma, so may need to go back
    also may need ot think about implimentation of sigma and eta
    may need to try other things.  
        once this works, then I should also try the low denisty stuff
    '''
    return (
        A0,
        R_max,
        R_np,
        T_m,
        bulk_dodecanethiol_density,
        bulk_toluene_density,
        eta,
        g_iso,
        nu_n,
        rho_decay,
        sigma,
        target_dr,
        tau,
    )


@app.cell
def _(
    A0,
    Model,
    R_max,
    R_np,
    T_m,
    bulk_dodecanethiol_density,
    bulk_toluene_density,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    eta,
    function_for_iso_fit,
    g_iso,
    hd_dodecanethiol_profile,
    hd_toluene_profile,
    len1,
    len2,
    np,
    nu_n,
    rho_decay,
    sigma,
    target_dr,
    tau,
):
    def my_callback(params, iter, resid, *args, **kws):
        # Only execute the print if 'iter' is a multiple of 1000
        if iter % 100 == 0:
            ssr = (resid**2).sum()
            print(f"Iteration: {iter:5d} | Chi-square (SSR): {ssr:.4e}")

        return False # Keep fitting




    triple_spectra_model = Model(function_for_iso_fit, independent_vars = ["xs", "len1", "len2", "solvent_density_obj", "rho_bulk_solvent", "ligand_density_obj", "rho_bulk_ligand", "R_np", "g_iso", "tau", "T_m", "target_dr"])
    triple_spectra_model_params = triple_spectra_model.make_params()

    triple_spectra_model_params.add_many(
        #("len1", , False), 
        #("len2", , False), # to be used for indices
        ("I1", np.max(dodecaneH_tolueneH_spectra[1]) - np.min(dodecaneH_tolueneH_spectra[1]), True, 0, None), 
        ("I2", np.max(dodecaneD_tolueneH_spectra[1]) - np.min(dodecaneD_tolueneH_spectra[1]), True, 0, None), 
        ("I3", np.max(dodecaneH_tolueneD_spectra[1]) - np.min(dodecaneH_tolueneD_spectra[1]), True, 0, None), # final intensity for each spectra
        ("c1", np.mean(dodecaneH_tolueneH_spectra[1][0:10]), False, None, None), 
        ("c2", np.mean(dodecaneD_tolueneH_spectra[1][0:10]), False, None, None),
        ("c3", np.mean(dodecaneH_tolueneD_spectra[1][0:10]), False, None, None), # final background for each spectra                  
        #r_grid_solvent, 
        #rho_grid_solvent, 
        #rho_bulk_solvent,                  
        #r_grid_ligand, 
        #rho_grid_ligand, 
        #rho_bulk_ligand,
        #R_np, 
        ("A0", A0, True, 0, None), 
        ("rho_decay", rho_decay, True, 0, None),
        #g_iso, 
        ("nu_n", nu_n, False, None, None),
        #tau, 
        #T_m, 
        ("sigma", sigma, False, 0, None),
        ("eta", eta, True, 0, None),
        ('R_max', R_max, False, R_np, None),
        ) 
    #triple_spectra_model_params['R_max'].diff_step = diff_step=dodecaneH_tolueneH_spectra[0][1] - dodecaneH_tolueneH_spectra[0][0] # Make diff_step for R_max roughly equal to your dr

    triple_spetra_fit = triple_spectra_model.fit(
        np.concatenate([dodecaneH_tolueneH_spectra[1], dodecaneD_tolueneH_spectra[1], dodecaneH_tolueneD_spectra[1]]), 
        xs = np.concatenate([dodecaneH_tolueneH_spectra[0], dodecaneD_tolueneH_spectra[0], dodecaneH_tolueneD_spectra[0]]),
        params = triple_spectra_model_params,
        len1 = len1,
        len2 = len2,
        solvent_density_obj = hd_toluene_profile, 
        rho_bulk_solvent = bulk_toluene_density,                  
        ligand_density_obj = hd_dodecanethiol_profile, 
        rho_bulk_ligand = bulk_dodecanethiol_density,
        R_np = R_np,
        g_iso = g_iso,
        tau = tau,
        T_m = T_m,
        target_dr = target_dr,
        iter_cb=my_callback
    )

    print(triple_spetra_fit.fit_report())
    return my_callback, triple_spectra_model_params, triple_spetra_fit


@app.cell
def _(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    len1,
    len2,
    make_subplots,
    triple_spetra_fit,
):
    fitPlot = make_subplots(rows = 3, cols = 1, subplot_titles= ["solvent and ligand", "solvent only", "ligand only"])
    fitPlot.update_layout(template = "simple_white", width = 500, height = 500)

    fitPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneH_spectra[1], line = dict(color = "#ffabff", width = 10), showlegend = False, col = 1, row = 1)
    fitPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[:len1], line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    fitPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneD_tolueneH_spectra[1], line = dict(color = "#ffabac", width = 10), showlegend = False, col = 1, row = 2)
    fitPlot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "black"), showlegend = False, col = 1, row = 2)

    fitPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneD_spectra[1], line = dict(color = "#ababff", width = 10), showlegend = False, col = 1, row = 3)
    fitPlot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= triple_spetra_fit.best_fit[len1+len2:], line = dict(color = "black"), showlegend = False, col = 1, row = 3)
    fitPlot.show()
    fitPlot.show("fitPlot.png")
    return


@app.cell
def _():
    return


@app.cell
def _(
    Model,
    R_max,
    R_np,
    T_m,
    bulk_dodecanethiol_density,
    bulk_toluene_density,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    eta,
    g_iso,
    hd_dodecanethiol_profile,
    hd_toluene_profile,
    len1,
    len2,
    my_callback,
    np,
    nu_n,
    sigma,
    target_dr,
    tau,
    triple_spectra_model_params,
):

    def simulate_mims_dynamic_dipolar_only(freq_axis, density_obj, R_np, 
                                           g_iso, nu_n, tau, T_m, sigma, rho_bulk, 
                                           eta, R_max, target_dr=0.2):
        # 1. DYNAMIC SHELL DEFINITION
        n_shells = int(np.ceil((R_max - R_np) / target_dr))
        n_shells = max(n_shells, 1) # Guardrail
        dr = (R_max - R_np) / n_shells 
        r_full = np.linspace(R_np + dr/2, R_max - dr/2, n_shells)

        # 2. SETUP
        cos_theta = np.linspace(0, 1, 150) 
        total_spectrum = np.zeros_like(freq_axis)
        damping = np.exp(-2 * tau / T_m)
        f_vec = freq_axis[np.newaxis, :] # Shape (1, 500)

        # 3. INTEGRATION
        for r in r_full:
            rho = density_obj.get_density(r, rho_bulk)
            if rho <= 0: continue

            # Pure Dipolar Physics
            T = (79.0 * (g_iso / 2.0023)) / (r**3) 
            A = T * (3 * cos_theta**2 - 1) # Shape (150,)

            # Mims Filter & Base Weighting
            intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping
            w = (rho * 4 * np.pi * r**2 * dr / len(cos_theta)) * intensity

            # --- ANISOTROPIC DIPOLAR STRAIN ---
            # shell_sigma varies for every orientation in the shell
            shell_sigma = np.sqrt(sigma**2 + (eta * np.abs(A)**0.5)**2)[:, np.newaxis] # Shape (150, 1)
        
            A_vec = A[:, np.newaxis]
            w_vec = w[:, np.newaxis]
            mu_plus = nu_n + A_vec / 2
            mu_minus = nu_n - A_vec / 2
        
            norm = shell_sigma * np.sqrt(2 * np.pi)
        
            # Broadcasting math
            g_plus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2)
            g_minus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2)

            total_spectrum += np.sum(g_plus + g_minus, axis=0)

        return total_spectrum

    # set up the fit

    def function_for_iso_fit_dipolar_only(xs, 
                             len1, len2, # to be used for indices
                             I1, I2, I3, # final intensity for each spectra
                             c1, c2, c3, # final background for each spectra
                             solvent_density_obj, rho_bulk_solvent,
                             ligand_density_obj, rho_bulk_ligand,
                             R_np, g_iso, nu_n, tau, T_m, sigma, 
                             eta, R_max, target_dr):
        # should sigma vary by sample? Or T_M?


        ys = np.array([]) # to hold the ys we calculate, concatenated together, so we can 

        #separate out th xs...
        freq_both =         xs[         :len1     ]
        freq_solvent_only = xs[len1     :len1+len2]
        freq_ligand_only =  xs[len1+len2:         ]

        # for each sample, calculate the contribution from the 

        for freq, I, c, s, l in zip(
            [freq_both, freq_solvent_only, freq_ligand_only],
            [I1, I2, I3],
            [c1, c2, c3],
            [1, 1, 0.01], # solvent weighting
            [1, 0.02, 1], # ligand weighting
        ):
            y_solvent = simulate_mims_dynamic_dipolar_only(freq, solvent_density_obj, R_np, g_iso, nu_n, tau, T_m, sigma, rho_bulk_solvent, eta, R_max, target_dr)
            y_ligand = simulate_mims_dynamic_dipolar_only(freq, ligand_density_obj, R_np, g_iso, nu_n, tau, T_m, sigma, rho_bulk_ligand, eta, R_max, target_dr)


            raw_ys = s*y_solvent + l*y_ligand

            norm_ys = raw_ys/np.max(raw_ys)

            scaled_ys = norm_ys*I + c

            ys = np.concatenate([ys, scaled_ys])


        return ys






    triple_spectra_model_dipolar_only = Model(function_for_iso_fit_dipolar_only, independent_vars = ["xs", "len1", "len2", "solvent_density_obj", "rho_bulk_solvent", "ligand_density_obj", "rho_bulk_ligand", "R_np", "g_iso", "tau", "T_m", "target_dr"])
    triple_spectra_model__diploar_only_params = triple_spectra_model_dipolar_only.make_params()

    triple_spectra_model__diploar_only_params.add_many(
        #("len1", , False), 
        #("len2", , False), # to be used for indices
        ("I1", np.max(dodecaneH_tolueneH_spectra[1]) - np.min(dodecaneH_tolueneH_spectra[1]), True, 0, None), 
        ("I2", np.max(dodecaneD_tolueneH_spectra[1]) - np.min(dodecaneD_tolueneH_spectra[1]), True, 0, None), 
        ("I3", np.max(dodecaneH_tolueneD_spectra[1]) - np.min(dodecaneH_tolueneD_spectra[1]), True, 0, None), # final intensity for each spectra
        ("c1", np.mean(dodecaneH_tolueneH_spectra[1][0:10]), False, None, None), 
        ("c2", np.mean(dodecaneD_tolueneH_spectra[1][0:10]), False, None, None),
        ("c3", np.mean(dodecaneH_tolueneD_spectra[1][0:10]), False, None, None), # final background for each spectra                  
        #r_grid_solvent, 
        #rho_grid_solvent, 
        #rho_bulk_solvent,                  
        #r_grid_ligand, 
        #rho_grid_ligand, 
        #rho_bulk_ligand,
        #R_np, 
        #g_iso, 
        ("nu_n", nu_n, False, None, None),
        #tau, 
        #T_m, 
        ("sigma", sigma, False, 0, None),
        ("eta", eta, True, 0, None),
        ('R_max', R_max, False, R_np, None),
        ) 
    #triple_spectra_model_params['R_max'].diff_step = diff_step=dodecaneH_tolueneH_spectra[0][1] - dodecaneH_tolueneH_spectra[0][0] # Make diff_step for R_max roughly equal to your dr

    triple_spetra_fit_dipolar_only = triple_spectra_model_dipolar_only.fit(
        np.concatenate([dodecaneH_tolueneH_spectra[1], dodecaneD_tolueneH_spectra[1], dodecaneH_tolueneD_spectra[1]]), 
        xs = np.concatenate([dodecaneH_tolueneH_spectra[0], dodecaneD_tolueneH_spectra[0], dodecaneH_tolueneD_spectra[0]]),
        params = triple_spectra_model_params,
        len1 = len1,
        len2 = len2,
        solvent_density_obj = hd_toluene_profile, 
        rho_bulk_solvent = bulk_toluene_density,                  
        ligand_density_obj = hd_dodecanethiol_profile, 
        rho_bulk_ligand = bulk_dodecanethiol_density,
        R_np = R_np,
        g_iso = g_iso,
        tau = tau,
        T_m = T_m,
        target_dr = target_dr,
        iter_cb=my_callback
    )

    print(triple_spetra_fit_dipolar_only.fit_report())
    return (triple_spetra_fit_dipolar_only,)


@app.cell
def _(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    len1,
    len2,
    make_subplots,
    triple_spetra_fit_dipolar_only,
):
    dipolarPlot = make_subplots(rows = 3, cols = 1, subplot_titles= ["solvent and ligand", "solvent only", "ligand only"])
    dipolarPlot.update_layout(template = "simple_white", width = 500, height = 500)

    dipolarPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneH_spectra[1], line = dict(color = "#ffabff", width = 10), showlegend = False, col = 1, row = 1)
    dipolarPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= triple_spetra_fit_dipolar_only.best_fit[:len1], line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    dipolarPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneD_tolueneH_spectra[1], line = dict(color = "#ffabac", width = 10), showlegend = False, col = 1, row = 2)
    dipolarPlot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit_dipolar_only.best_fit[len1:len1+len2], line = dict(color = "black"), showlegend = False, col = 1, row = 2)

    dipolarPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneD_spectra[1], line = dict(color = "#ababff", width = 10), showlegend = False, col = 1, row = 3)
    dipolarPlot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= triple_spetra_fit_dipolar_only.best_fit[len1+len2:], line = dict(color = "black"), showlegend = False, col = 1, row = 3)
    dipolarPlot.show()
    dipolarPlot.write_image("dipolarPlot.png")
    return


@app.cell
def _():
    return


@app.cell
def _(
    A0,
    Model,
    R_max,
    R_np,
    T_m,
    bulk_dodecanethiol_density,
    bulk_toluene_density,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    eta,
    g_iso,
    hd_dodecanethiol_profile,
    hd_toluene_profile,
    len1,
    len2,
    my_callback,
    np,
    nu_n,
    rho_decay,
    sigma,
    target_dr,
    tau,
    triple_spectra_model_params,
):

    def simulate_mims_dynamic_total_strain(freq_axis, density_obj, R_np, A0, rho_decay, 
                                           g_iso, nu_n, tau, T_m, sigma, rho_bulk, eta, R_max, target_dr=0.2):
        """
        Simulates a Mims ENDOR spectrum where:
        1. Shell thickness (dr) adjusts dynamically to R_max for smooth fitting.
        2. Density is pulled from the DensityProfile object (MD + Blending + Bulk).
        3. Broadening (strain) scales with the TOTAL coupling (A), not just A_iso.
        """

        # --- 1. DYNAMIC SHELL DEFINITION ---
        # Ensure R_max is physically outside the NP core to prevent errors
        R_max = max(R_max, R_np + 0.1)
    
        # Calculate exact number of shells to make dr a fluid parameter
        n_shells = int(np.ceil((R_max - R_np) / target_dr))
        n_shells = max(n_shells, 1) # Safety guard
    
        dr = (R_max - R_np) / n_shells 
    
        # Midpoints of our dynamic shells (Shape: n_shells)
        r_full = np.linspace(R_np + dr/2, R_max - dr/2, n_shells)

        # --- 2. SETUP GRIDS & CONSTANTS ---
        df = freq_axis[1] - freq_axis[0]
        cos_theta = np.linspace(0, 1, 150) # 150 orientations for powder average
        total_spectrum = np.zeros_like(freq_axis)
        damping = np.exp(-2 * tau / T_m)
    
        # Pre-vectorize the frequency axis for broadcasting (Shape: 1 x 500)
        f_vec = freq_axis[np.newaxis, :]

        # --- 3. INTEGRATION LOOP (Radial) ---
        for r in r_full:
            # Pull density from the object (handles MD, Blending, and Toluene Bulk)
            rho = density_obj.get_density(r, rho_bulk)
        
            if rho <= 0: continue

            # Physics Calculation: Dipolar (T) and Contact (A_iso)
            T = (79.0 * (g_iso / 2.0023)) / (r**3) 
            A_iso = A0 * np.exp(-(max(0, r - R_np)) / rho_decay)
        
            # A is an array of 150 orientations (Shape: 150,)
            A = A_iso + T * (3 * cos_theta**2 - 1) 

            # Mims Filter & Base Weighting
            # Population weighting: (Shell Volume * Density * Mims Efficiency)
            intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping
            w = (rho * 4 * np.pi * r**2 * dr / len(cos_theta)) * intensity # Shape: (150,)

            # --- ANISOTROPIC HYPERFINE STRAIN ---
            # shell_sigma now varies for EVERY orientation in the shell.
            # We use np.abs(A) so broadening scales with the magnitude of the shift.
            shell_sigma = np.sqrt(sigma**2 + (eta * np.abs(A)**0.5)**2)[:, np.newaxis] # Shape: (150, 1)

            # Prepare vectors for broadcasting
            A_vec = A[:, np.newaxis]   # Shape: (150, 1)
            w_vec = w[:, np.newaxis]   # Shape: (150, 1)
        
            # Calculate branch centers
            mu_plus = nu_n + A_vec / 2
            mu_minus = nu_n - A_vec / 2
        
            # Area normalization: Crucial when sigma varies per orientation
            norm = shell_sigma * np.sqrt(2 * np.pi)

            # Vectorized Gaussian calculation: Result is a (150 x 500) matrix
            g_plus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2)
            g_minus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2)

            # Sum all 150 orientations (rows) into the master spectrum
            total_spectrum += np.sum(g_plus + g_minus, axis=0)

        return total_spectrum

    def function_for_total_strain_fit(xs, 
                             len1, len2, # to be used for indices
                             I1, I2, I3, # final intensity for each spectra
                             c1, c2, c3, # final background for each spectra
                             solvent_density_obj, rho_bulk_solvent,
                             ligand_density_obj, rho_bulk_ligand,
                             R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, 
                             eta, R_max, target_dr):
        # should sigma vary by sample? Or T_M?


        ys = np.array([]) # to hold the ys we calculate, concatenated together, so we can 

        #separate out th xs...
        freq_both =         xs[         :len1     ]
        freq_solvent_only = xs[len1     :len1+len2]
        freq_ligand_only =  xs[len1+len2:         ]

        # for each sample, calculate the contribution from the 

        for freq, I, c, s, l in zip(
            [freq_both, freq_solvent_only, freq_ligand_only],
            [I1, I2, I3],
            [c1, c2, c3],
            [1, 1, 0.01], # solvent weighting
            [1, 0.02, 1], # ligand weighting
        ):
            y_solvent = simulate_mims_dynamic_total_strain(freq, solvent_density_obj, R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk_solvent, eta, R_max, target_dr)
            y_ligand = simulate_mims_dynamic_total_strain(freq, ligand_density_obj, R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk_ligand, eta, R_max, target_dr)


            raw_ys = s*y_solvent + l*y_ligand

            norm_ys = raw_ys/np.max(raw_ys)

            scaled_ys = norm_ys*I + c

            ys = np.concatenate([ys, scaled_ys])


        return ys

    triple_spectra_total_strain_model = Model(function_for_total_strain_fit, independent_vars = ["xs", "len1", "len2", "solvent_density_obj", "rho_bulk_solvent", "ligand_density_obj", "rho_bulk_ligand", "R_np", "g_iso", "tau", "T_m", "target_dr"])
    triple_spectra_model_total_strain_params = triple_spectra_total_strain_model.make_params()

    triple_spectra_model_total_strain_params.add_many(
        #("len1", , False), 
        #("len2", , False), # to be used for indices
        ("I1", np.max(dodecaneH_tolueneH_spectra[1]) - np.min(dodecaneH_tolueneH_spectra[1]), True, 0, None), 
        ("I2", np.max(dodecaneD_tolueneH_spectra[1]) - np.min(dodecaneD_tolueneH_spectra[1]), True, 0, None), 
        ("I3", np.max(dodecaneH_tolueneD_spectra[1]) - np.min(dodecaneH_tolueneD_spectra[1]), True, 0, None), # final intensity for each spectra
        ("c1", np.mean(dodecaneH_tolueneH_spectra[1][0:10]), False, None, None), 
        ("c2", np.mean(dodecaneD_tolueneH_spectra[1][0:10]), False, None, None),
        ("c3", np.mean(dodecaneH_tolueneD_spectra[1][0:10]), False, None, None), # final background for each spectra                  
        #r_grid_solvent, 
        #rho_grid_solvent, 
        #rho_bulk_solvent,                  
        #r_grid_ligand, 
        #rho_grid_ligand, 
        #rho_bulk_ligand,
        #R_np, 
        ("A0", A0, True, 0, None), 
        ("rho_decay", rho_decay, True, 0, None),
        #g_iso, 
        ("nu_n", nu_n, False, None, None),
        #tau, 
        #T_m, 
        ("sigma", sigma, False, 0, None),
        ("eta", eta, True, 0, None),
        ('R_max', R_max, False, R_np, None),
        ) 
    #triple_spectra_model_params['R_max'].diff_step = diff_step=dodecaneH_tolueneH_spectra[0][1] - dodecaneH_tolueneH_spectra[0][0] # Make diff_step for R_max roughly equal to your dr

    triple_spetra_total_strain_fit = triple_spectra_total_strain_model.fit(
        np.concatenate([dodecaneH_tolueneH_spectra[1], dodecaneD_tolueneH_spectra[1], dodecaneH_tolueneD_spectra[1]]), 
        xs = np.concatenate([dodecaneH_tolueneH_spectra[0], dodecaneD_tolueneH_spectra[0], dodecaneH_tolueneD_spectra[0]]),
        params = triple_spectra_model_params,
        len1 = len1,
        len2 = len2,
        solvent_density_obj = hd_toluene_profile, 
        rho_bulk_solvent = bulk_toluene_density,                  
        ligand_density_obj = hd_dodecanethiol_profile, 
        rho_bulk_ligand = bulk_dodecanethiol_density,
        R_np = R_np,
        g_iso = g_iso,
        tau = tau,
        T_m = T_m,
        target_dr = target_dr,
        iter_cb=my_callback
    )

    print(triple_spetra_total_strain_fit.fit_report())
    return (triple_spetra_total_strain_fit,)


@app.cell
def _(triple_spetra_total_strain_fit):
    print(triple_spetra_total_strain_fit.fit_report())
    return


@app.cell
def _(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    len1,
    len2,
    make_subplots,
    triple_spetra_total_strain_fit,
):
    totalPlot = make_subplots(rows = 3, cols = 1, subplot_titles= ["solvent-H + ligand-H", "solvent-H + ligand-D", "solvent-D + ligand-H"])
    totalPlot.update_layout(template = "simple_white", width = 500, height = 500)

    totalPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneH_spectra[1], line = dict(color = "#ffabff", width = 5), showlegend = False, col = 1, row = 1)
    totalPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= triple_spetra_total_strain_fit.best_fit[:len1], line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    totalPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneD_tolueneH_spectra[1], line = dict(color = "#ffabac", width = 5), showlegend = False, col = 1, row = 2)
    totalPlot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_total_strain_fit.best_fit[len1:len1+len2], line = dict(color = "black"), showlegend = False, col = 1, row = 2)

    totalPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneD_spectra[1], line = dict(color = "#ababff", width = 5), showlegend = False, col = 1, row = 3)
    totalPlot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= triple_spetra_total_strain_fit.best_fit[len1+len2:], line = dict(color = "black"), showlegend = False, col = 1, row = 3)
    totalPlot.update_xaxes(title = "frequency /MHz", row = 3, col = 1)
    totalPlot.show()
    totalPlot.write_image("totalFit.png")
    return


@app.cell
def _(triple_spetra_total_strain_fit):
    print(triple_spetra_total_strain_fit.fit_report())
    return


@app.cell
def _(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    len1,
    len2,
    make_subplots,
    triple_spetra_fit,
    triple_spetra_fit_dipolar_only,
    triple_spetra_total_strain_fit,
):
    comparisonPlot = make_subplots(rows = 3, cols = 1, subplot_titles= ["solvent and ligand", "solvent only", "ligand only"])
    comparisonPlot.update_layout(template = "simple_white", width = 900, height = 500)

    comparisonPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[:len1], line = dict(color = "red"), name = "contact broadening + dipolar", showlegend = True, col = 1, row = 1)
    comparisonPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= triple_spetra_fit_dipolar_only.best_fit[:len1], line = dict(color = "blue"), name = "dipolar broadening", showlegend = True, col = 1, row = 1)
    comparisonPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= triple_spetra_total_strain_fit.best_fit[:len1], line = dict(color = "darkcyan"), name = "contact broadening + dipolar broadening", showlegend = False, col = 1, row = 1)

    comparisonPlot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "red"), showlegend = False, col = 1, row = 2)
    comparisonPlot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit_dipolar_only.best_fit[len1:len1+len2], line = dict(color = "blue"), showlegend = False, col = 1, row = 2)
    comparisonPlot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_total_strain_fit.best_fit[len1:len1+len2], line = dict(color = "darkcyan"), showlegend = False, col = 1, row = 2)


    comparisonPlot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= triple_spetra_fit.best_fit[len1+len2:], line = dict(color = "red"), showlegend = False, col = 1, row = 3)
    comparisonPlot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= triple_spetra_fit_dipolar_only.best_fit[len1+len2:], line = dict(color = "blue"), showlegend = False, col = 1, row = 3)
    comparisonPlot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= triple_spetra_total_strain_fit.best_fit[len1+len2:], line = dict(color = "darkcyan"), name = "contact broadening + diploar broadening", showlegend = True, col = 1, row = 3)
    comparisonPlot.update_xaxes(range=[46, 54])
    comparisonPlot.show()
    comparisonPlot.write_image("comparison.png")
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _(np, simulate_mims_component):
    def function_for_iso_fit_for_functions_below(xs, 
                             len1, len2, # to be used for indices
                             I1, I2, I3, # final intensity for each spectra
                             c1, c2, c3, # final background for each spectra
                             r_grid_solvent, rho_grid_solvent, rho_bulk_solvent,
                             r_grid_ligand, rho_grid_ligand, rho_bulk_ligand,
                             R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, 
                             eta, R_max):
        # should sigma vary by sample? Or T_M?


        ys = np.array([]) # to hold the ys we calculate, concatenated together, so we can 

        #separate out th xs...
        freq_both =         xs[         :len1     ]
        freq_solvent_only = xs[len1     :len1+len2]
        freq_ligand_only =  xs[len1+len2:         ]

        # for each sample, calculate the contribution from the 

        for freq, I, c, s, l in zip(
            [freq_both, freq_solvent_only, freq_ligand_only],
            [I1, I2, I3],
            [c1, c2, c3],
            [1, 1, 0.01], # solvent weighting
            [1, 0.02, 1], # ligand weighting
        ):
            y_solvent = simulate_mims_component(freq, r_grid_solvent, s*np.array(rho_grid_solvent), R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk_solvent, eta, R_max)
            y_ligand = simulate_mims_component(freq, r_grid_ligand, l*np.array(rho_grid_ligand), R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk_ligand, eta, R_max)


            raw_ys = y_solvent + y_ligand

            norm_ys = raw_ys/np.max(raw_ys)

            scaled_ys = norm_ys*I + c

            ys = np.concatenate([ys, scaled_ys])


        return ys


    def simulate_mims_component_contact_only_strain(freq_axis, r_grid, rho_grid, R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk, eta, R_max):
        """
        Simulates a Mims ENDOR spectrum with Frequency-Dependent (Hyperfine) Strain.

        Parameters:
        -----------
        ... (all previous parameters) ...
        eta : float
            Hyperfine strain (dimensionless). Scales broadening relative to A_iso.
            Recommended range: 0.0 to 0.2.
        """

        # --- 1. SET THE HORIZON ---
        threshold = 0.15 * sigma
        dr = r_grid[1] - r_grid[0]
        r_dip_limit = (158.0 * (g_iso / 2.0023) / threshold)**(1/3)
        r_con_limit = R_np - rho_decay * np.log(threshold / max(1e-9, A0)) if A0 > threshold else R_np

        # Using the tighter of the physics limits or MD data
        #R_max = max(r_dip_limit, r_con_limit)

        # --- 2. SETUP GRIDS ---
        df = freq_axis[1] - freq_axis[0]
        cos_theta = np.linspace(0, 1, 150) 
        total_spectrum = np.zeros_like(freq_axis)
        damping = np.exp(-2 * tau / T_m)
        r_full = np.arange(R_np, R_max, dr)

        # --- 3. INTEGRATION LOOP ---
        for r in r_full:
            rho = np.interp(r, r_grid, rho_grid) if r <= r_grid[-1] else rho_bulk
            if rho <= 0: continue

            # Physics calculation
            T = (79.0 * (g_iso / 2.0023)) / (r**3) 
            A_iso = A0 * np.exp(-(max(0, r - R_np)) / rho_decay)
            A = A_iso + T * (3 * cos_theta**2 - 1)

            # Mims Filter & Base Weighting
            intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping
            w = (rho * 4 * np.pi * r**2 * dr / len(cos_theta)) * intensity

            # --- NEW: HYPERFINE STRAIN BROADENING ---
            # Calculate a unique sigma for this shell. 
            # As A_iso gets larger (closer to surface), broadening increases.
            shell_sigma = np.sqrt(sigma**2 + (eta * A_iso)**2)

            # Instead of histogramming, we calculate the Gaussian contribution 
            # of each orientation in the shell.
            for i, coupling in enumerate(A):
                if w[i] <= 0: continue

                # Center of the two ENDOR branches
                mu_plus = nu_n + coupling/2
                mu_minus = nu_n - coupling/2

                # Apply Gaussian shape directly
                # (Note: For extreme speed in large fits, you can use a 
                # pre-computed Gaussian look-up table here)
                total_spectrum += w[i] * np.exp(-0.5 * ((freq_axis - mu_plus) / shell_sigma)**2)
                total_spectrum += w[i] * np.exp(-0.5 * ((freq_axis - mu_minus) / shell_sigma)**2)

        # Normalize by shell_sigma to keep area consistent if needed, 
        # but usually total intensity is scaled by a fit parameter anyway.
        return total_spectrum


    def simulate_mims_component_contact_and_dipolar_strain(freq_axis, r_grid, rho_grid, R_np, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk, eta=0.0, R_max=30):
        # --- 1. SETUP ---
        dr = r_grid[1] - r_grid[0]
        df = freq_axis[1] - freq_axis[0]
        cos_theta = np.linspace(0, 1, 150) 
        total_spectrum = np.zeros_like(freq_axis)
        damping = np.exp(-2 * tau / T_m)
        r_full = np.arange(R_np, R_max, dr)

        # Prepare frequency axis for broadcasting (Shape: 1 x 500)
        f_vec = freq_axis[np.newaxis, :]

        # --- 2. INTEGRATION LOOP (Radial) ---
        for r in r_full:
            rho = np.interp(r, r_grid, rho_grid) if r <= r_grid[-1] else rho_bulk
            if rho <= 0: continue

            # Physics
            T = (79.0 * (g_iso / 2.0023)) / (r**3) 
            A_iso = A0 * np.exp(-(max(0, r - R_np)) / rho_decay)
            A = A_iso + T * (3 * cos_theta**2 - 1) # Shape: (150,)

            # Mims Filter & Base Weighting
            intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping
            w = (rho * 4 * np.pi * r**2 * dr / len(cos_theta)) * intensity # Shape: (150,)

            # --- VECTORIZED BROADENING ---
            # 1. Calculate the 150 unique sigmas (Shape: 150 x 1)
            # Use absolute A to ensure strain is positive
            shell_sigma = np.sqrt(sigma**2 + (eta * np.abs(A))**2)[:, np.newaxis]

            # 2. Reshape A and weights for broadcasting (Shape: 150 x 1)
            A_vec = A[:, np.newaxis]
            w_vec = w[:, np.newaxis]

            # 3. Calculate Centers (Shape: 150 x 1)
            mu_plus = nu_n + A_vec / 2
            mu_minus = nu_n - A_vec / 2

            # 4. Matrix Math (Shape: 150 x 500)
            # This calculates all orientations and all frequencies in one step
            # Gaussian = (1 / (s * sqrt(2pi))) * exp(-0.5 * ((f - mu)/s)^2)
            norm = shell_sigma * np.sqrt(2 * np.pi)

            g_plus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2)
            g_minus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2)

            # 5. Sum the 150 orientations (rows) into the total spectrum
            total_spectrum += np.sum(g_plus + g_minus, axis=0)

        return total_spectrum

    def simulate_mims_dynamic_contact_and_dipolar_strain(freq_axis, density_obj, R_np, R_max, A0, rho_decay, 
                              g_iso, nu_n, tau, T_m, sigma, rho_bulk, eta=0.0):
        """
        Simulation that defines bins dynamically based on R_max.
        """
        # 1. Define Dynamic dr
        # We want dr to be ~0.2, but we need an exact integer number of shells
        target_dr = 0.2
        n_shells = int(np.ceil((R_max - R_np) / target_dr))
        # This dr changes slightly as R_max changes, making the function 'smooth'
        dr = (R_max - R_np) / n_shells 

        # 2. Setup Grids
        r_centers = np.linspace(R_np + dr/2, R_max - dr/2, n_shells)
        cos_theta = np.linspace(0, 1, 150)
        total_spectrum = np.zeros_like(freq_axis)
        damping = np.exp(-2 * tau / T_m)
        f_vec = freq_axis[np.newaxis, :]

        # 3. Integration Loop
        for r in r_centers:
            # Use the class method to get density (handles MD, Blending, and Bulk)
            rho = density_obj.get_density(r, rho_bulk)
            if rho <= 0: continue

            # Physics (Contact + Dipolar)
            T = (79.0 * (g_iso / 2.0023)) / (r**3)
            A_iso = A0 * np.exp(-(r - R_np) / rho_decay)
            A = A_iso + T * (3 * cos_theta**2 - 1)

            # Mims Filter & Weight
            intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping
            w = (rho * 4 * np.pi * r**2 * dr / len(cos_theta)) * intensity

            # Vectorized Broadening
            shell_sigma = np.sqrt(sigma**2 + (eta * np.abs(A))**2)[:, np.newaxis]
            mu_plus = nu_n + (A[:, np.newaxis] / 2)
            mu_minus = nu_n - (A[:, np.newaxis] / 2)
            w_vec = w[:, np.newaxis]

            norm = shell_sigma * np.sqrt(2 * np.pi)
            g_plus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2)
            g_minus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2)

            total_spectrum += np.sum(g_plus + g_minus, axis=0)

        return total_spectrum

    return


if __name__ == "__main__":
    app.run()
