import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    the point of this notebook is to work from the Mims with dipolar and see if we can add in the dependnece on particle size...
    """)
    return


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import revisedMimsFunctions as mf
    from importlib import reload
    import codechembook.quickPlots as qp
    from scipy.ndimage import gaussian_filter1d
    from plotly.subplots import make_subplots
    from lmfit import Model

    return Model, make_subplots, mf, mo, np


@app.cell
def _(make_subplots, mo, np):
    # import spectral data
    dodecaneH_tolueneH_spectra =  np.genfromtxt("//Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)

    dodecaneD_tolueneH_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter="," ,unpack = True, skip_header = 1)

    dodecaneH_tolueneD_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter="," ,unpack = True, skip_header = 1)

    names = ["both", "solvent only", "ligand only"]
    exp_plot = make_subplots(rows = 3, cols = 1)
    exp_plot.update_layout(template = "simple_white")
    for _i, _exp, in enumerate(
        [
            dodecaneH_tolueneH_spectra,
            dodecaneD_tolueneH_spectra,
            dodecaneH_tolueneD_spectra
        ], start=0):


        exp_plot.add_scatter(x = _exp[0], y = _exp[1], 
                             name = f"{names[_i]}", row = _i+1, col = 1)

    mo.md(f"{mo.ui.plotly(exp_plot)}")
    return (
        dodecaneD_tolueneH_spectra,
        dodecaneH_tolueneD_spectra,
        dodecaneH_tolueneH_spectra,
    )


@app.cell
def _(R_np_md, bulk_toluene_density, make_subplots, md_bin_size, mf, np):
    # start with solvent from MD
    hd_toluene_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/toluene_distances_high_density.csv",
                                       delimiter = ",")
    hd_toluene_profile = mf.NPDensityProfile(hd_toluene_distances, R_np_md, bin_width=md_bin_size)

    # then blend to bulk
    hd_toluene_profile.prepare_blended_density(rho_bulk = bulk_toluene_density) 

    print(np.max(hd_toluene_profile.bin_centers))

    # then get ligands

    hd_dodecanethiol_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/dodecanol_distances_high_density.csv",
                                       delimiter = ",")

    #set up profile to have the same bins as the solvent
    hd_dodecanethiol_profile = mf.NPDensityProfile(hd_dodecanethiol_distances, R_np_md, old_bins = hd_toluene_profile.bin_centers)



    # then modify the profile for solvent, so that we get to the bulk values. 

    densityPLot  = make_subplots()

    densityPLot.add_scatter(x = np.array(hd_dodecanethiol_profile.bin_centers)/10, y = hd_dodecanethiol_profile.densities, 
                                  line = dict(color = "#BB6F35"), fill = "tozeroy")
    densityPLot.add_scatter(x = np.array(hd_toluene_profile.bin_centers)/10, y = hd_toluene_profile.densities, 
                           line = dict(color = "#448CC4"), fill = "tozeroy")
    #densityPLot.add_scatter(x = np.array(hd_toluene_profile.bin_centers)/10, y = hd_toluene_profile.densities+hd_dodecanethiol_profile.densities, line = dict(color = "#94A2B1"))

    densityPLot.add_scatter(x = [R_np_md/10, R_np_md/10], y = [0, 1], mode = "lines", line = dict(color = "black", dash = "dash"))

    densityPLot.update_layout(template = "simple_white",showlegend = False, width = 300*2, height = 300*1.4)
    densityPLot.update_xaxes(range = [0, 3.4], title = "radial distance from center of particle /nm")
    densityPLot.update_yaxes(title = "H density per A^3")
    densityPLot.show()
    return hd_dodecanethiol_profile, hd_toluene_profile


@app.cell
def single_mims(np):
    # this is the underlying physics for simulating a single Mims spectra

    def simulate_mims_dynamic_contact_strain(freq_axis, density_obj, R_np_md, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk, eta, R_max, target_dr, shape, R_np_exp):
        """
        Simulates a Mims ENDOR spectrum where the shell thickness (dr) adjusts 
        dynamically to R_max, and density is pulled from a DensityProfile object.

        PARAMETERS:
        freq_axis: numpy array of floats, frequency axis, units = MHz
        density_obj: density object, from a special class, contains information about density of hydrogens from MD simulations, units = A^3
        R_np_md: float, the radius of the nanoparticle used in MD simulations, units = A
        A0: float,  "Normalization" factor for the metallic electron wavefunction at a nanoparticle radius of 0A, units = 
        rho_decay: float, the characteristic length scale for decay of the electron wavefunction from the nanoparticle surface, units = A
        g_iso: float, electronic g-factor, unitless
        nu_n: float, the resonance frequency for the hydrogens, units = MHz
        tau: float, delay time for pulse sequence, units = microseconds
        T_m: 
        sigma:
        rho_bulk:
        eta:
        R_max:
        target_dr: 
        shape: 
        R_np_exp: 
        """

        # --- 1. DYNAMIC SHELL DEFINITION ---
        # Target dr is ~0.2 A. We calculate an exact n_shells to keep dr fluid.
        n_shells = int(np.ceil((R_max - R_np_md) / target_dr)) # base number of shells off the md processing

        # This 'dr' changes smoothly as R_max moves, preventing 'jumpy' fits.
        dr = (R_max - R_np_md) / n_shells 

        # Midpoints of our dynamic shells
        # we need both rs for md (to pull density) and experimental (to calcualte coupling etc)
        r_full_md = np.linspace(R_np_md + dr/2, R_np_md + dr*n_shells - dr/2, n_shells)
        r_full_exp = np.linspace(R_np_exp + dr/2, R_np_exp + dr*n_shells - dr/2, n_shells) # but generate the actuall rs, based on experimental values of r. 

        # --- 2. SETUP GRIDS ---
        df = freq_axis[1] - freq_axis[0]
        cos_theta = np.linspace(0, 1, 150) 
        total_spectrum = np.zeros_like(freq_axis)
        damping = np.exp(-2 * tau / T_m)

        # Pre-vectorize the frequency axis for the Gaussian calculation
        f_vec = freq_axis[np.newaxis, :]

        # --- 3. INTEGRATION LOOP (Radial) ---
        for r_md, r_exp in zip(r_full_md, r_full_exp):
            # Pull density from the object (handles MD, Blending, and Toluene Bulk)
        
            rho = density_obj.get_density(r_md, rho_bulk)

            if rho <= 0 or r_md < 17.5: continue

            # Physics Calculation
            T = (79.0 * (g_iso / 2.0023)) / (r_exp**3) 

            scaled_A0 = A0 * np.sqrt(1/(2*np.pi*rho_decay*(r_exp**2 + r_exp*rho_decay + 0.5*rho_decay**2))) # adjust for the fact that different sized nanoparticles will have different intensities of the wavefunction at the surface
        
            A_contact = (scaled_A0 * np.exp(-(max(0, r_exp - R_np_exp)) / rho_decay))**2
            A = A_contact + T * (3 * cos_theta**2 - 1) # Array of 150 orientations

            # Mims Filter & Base Weighting
            intensity = (1 - np.cos(2 * np.pi * A * tau)) * damping
            w = (rho * 4 * np.pi * r_exp**2 * dr / len(cos_theta)) * intensity

            # --- HYPERFINE STRAIN BROADENING (Contact-Only) ---
            # shell_sigma is a single value for this shell since A_iso is isotropic
            shell_sigma = np.sqrt(sigma**2 + (eta * np.abs(A_contact)**0.5)**2)  # the dependence on A_iso can be changed to reflect different disorders

            # Vectorized Gaussian Calculation for all 150 orientations at once
            A_vec = A[:, np.newaxis]
            w_vec = w[:, np.newaxis]
            mu_plus = nu_n + A_vec / 2
            mu_minus = nu_n - A_vec / 2

            # Area normalization: keeps intensity stable as the peak width changes
            if shape[0] == "G": # then this is a gaussian profile
                norm = shell_sigma * np.sqrt(2 * np.pi)
    
                g_plus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2) # divide by norm to normalize
                g_minus = (w_vec / norm) * np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2) # divide by norm to normalize
        
            elif shape[0] == "L": # lorentzian profile
                norm = np.pi * shell_sigma
                g_plus = (w_vec / norm) * (shell_sigma / ((f_vec - mu_plus)**2 + shell_sigma**2))
                g_minus = (w_vec / norm) * (shell_sigma / ((f_vec - mu_minus)**2 + shell_sigma**2))
            
            # Sum orientations into the master spectrum
            total_spectrum += np.sum(g_plus + g_minus, axis=0)

        return total_spectrum

    return (simulate_mims_dynamic_contact_strain,)


@app.cell
def simulate(
    A0,
    R_max,
    R_np_md,
    Rs_np_D,
    Rs_np_H,
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
    make_subplots,
    np,
    nu_n,
    rho_decay,
    sigma,
    simulate_mims_dynamic_contact_strain,
    target_dr,
    tau,
):
    def function_for_iso_fit(xs, 
                             len1, len2, # to be used for indices
                             I1, I2, I3, # final intensity for each spectra
                             c1, c2, c3, # final background for each spectra
                             solvent_density_obj, rho_bulk_solvent,
                             ligand_density_obj, rho_bulk_ligand,
                             R_np_md, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, 
                             eta, R_max, target_dr,
                            shape,
                            Rs_np_H, Rs_np_D, # paired list of nanoparticle radius and density
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
            [1, 1, 0.01], # solvent weighting
            [1, 0.02, 1], # ligand weighting
            [Rs_np_H, Rs_np_D, Rs_np_H],
        ):
            raw_ys = np.zeros_like(freq) # make an empty
            for r, d in zip(rs[0], rs[1]): # for each experimental radius and density pair...
                if d>0:
                    y_solvent = simulate_mims_dynamic_contact_strain(freq, solvent_density_obj, R_np_md, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk_solvent, eta, R_max, target_dr, shape, r)
                    y_ligand = simulate_mims_dynamic_contact_strain(freq, ligand_density_obj, R_np_md, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, rho_bulk_ligand, eta, R_max, target_dr, shape, r)
    
                    raw_ys = raw_ys + d*(s*y_solvent + l*y_ligand)

            norm_ys = raw_ys/np.max(raw_ys)

            scaled_ys = norm_ys*I + c

            ys = np.concatenate([ys, scaled_ys])


        return ys



    # simulate based on this

    sim_ys = function_for_iso_fit(
        xs = np.concatenate([dodecaneH_tolueneH_spectra[0],
                             dodecaneD_tolueneH_spectra[0],
                             dodecaneH_tolueneD_spectra[0]]), 
        len1 = len1, len2 = len2, # to be used for indices
                         
        I1 = np.max(dodecaneH_tolueneH_spectra[1]) - np.min(dodecaneH_tolueneH_spectra[1]), 
        I2 = np.max(dodecaneD_tolueneH_spectra[1]) - np.min(dodecaneD_tolueneH_spectra[1]), 
        I3 = np.max(dodecaneH_tolueneD_spectra[1]) - np.min(dodecaneH_tolueneD_spectra[1]), # final intensity for each spectra
                         
        c1 = np.mean(dodecaneH_tolueneH_spectra[1][0:10]), 
        c2 = np.mean(dodecaneD_tolueneH_spectra[1][0:10]), 
        c3 = np.mean(dodecaneH_tolueneD_spectra[1][0:10]), # final background for each spectra
                         
        solvent_density_obj= hd_toluene_profile, rho_bulk_solvent = bulk_toluene_density,
        ligand_density_obj = hd_dodecanethiol_profile, rho_bulk_ligand = bulk_dodecanethiol_density,
                         
        R_np_md = R_np_md, 
        A0 = A0, 
        rho_decay = rho_decay, 
        g_iso = g_iso, 
        nu_n = nu_n, 
        tau = tau, 
        T_m = T_m, 
        sigma = sigma, 
        eta = eta, 
        R_max = R_max, 
        target_dr= target_dr,
        shape = "G",
        Rs_np_H = Rs_np_H,
        Rs_np_D = Rs_np_D,
    )



    simPlot = make_subplots(rows = 3, cols = 1)
    simPlot.update_layout(template = "simple_white")

    simPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneH_spectra[1], line = dict(color = "lightgrey", width = 10), col = 1, row = 1)
    simPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= sim_ys[:len1], line = dict(color = "red"), col = 1, row = 1)

    simPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneD_tolueneH_spectra[1], line = dict(color = "lightgrey", width = 10), col = 1, row = 2)
    simPlot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= sim_ys[len1:len1+len2], line = dict(color = "red"), col = 1, row = 2)

    simPlot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneD_spectra[1], line = dict(color = "lightgrey", width = 10), col = 1, row = 3)
    simPlot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= sim_ys[len1+len2:], line = dict(color = "red"), col = 1, row = 3)
    simPlot.show()
    return (function_for_iso_fit,)


@app.cell
def _(dodecaneD_tolueneH_spectra, dodecaneH_tolueneH_spectra, np):
    # constants

    md_bin_size = 0.1 # in Angstrom


    len1 = len(dodecaneH_tolueneH_spectra[0])
    len2 = len(dodecaneD_tolueneH_spectra[0])

    R_np_md = 15 # angstrom
    T_m = 640 # in microseconds. 
    bulk_toluene_density = 0.0453 # number of hydrogens per cubic angstrom
    bulk_dodecanethiol_density = 0 # the density of hydrogens coming from the ligand at r = infinity
    g_iso = 1.9 # g-value for the nanoparticles
    eta = 0 # was 22 # refects the heterogeneity. value of 0 means no effect

    # to be fit
    nu_n = 50 # the frequency of the hydrogens
    A0 = 400 #
    rho_decay = 10 # the decay constant for wavefunction, units = A
    sigma = 0.15 # to reflect the dynamics of the particles and resolution of instrument
    R_max = 50 # the max radius we look out to

    # properties of the EPR experiment
    tau = 140

    #properties of fit
    target_dr = 0.1

    Rs_np_H = [
        np.array([1.272975688,
    1.405767063,
    1.538558438,
    1.671349814,
    1.804141189,
    1.936932564,
    2.06972394,
    2.202515315, 
    2.33530669, 
    2.468098066, 
    2.600889441, 
    2.733680816, 
    2.866472192, 
    2.999263567, 
    3.132054942, 
    3.264846318, 
    3.397637693, 
    3.530429068, 
    3.663220444, 
    3.796011819, 
    3.928803194, 
    4.06159457, 
    4.194385945, 
    4.32717732, 
    4.459968696, 
    4.592760071, 
    4.725551446, 
    4.858342822, 
    4.991134197, 
    5.123925572, 
    5.256716947, 
    5.389508323, 
    5.522299698, 
    5.655091073, 
    5.787882449, 
    5.920673824, 
    6.053465199, 
    6.186256575,])*10/2,
        [4, 
    1, 
    2, 
    1, 
    0, 
    4, 
    1, 
    11, 
    17, 
    49, 
    82, 
    83, 
    81, 
    59, 
    55, 
    36, 
    29, 
    20, 
    25, 
    8, 
    3, 
    6, 
    10, 
    5, 
    1, 
    1, 
    2, 
    0, 
    1, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    1,]
    ]


    #[[28.7/2],[1]] # this is the mean size from fitting the TEM
    Rs_np_D = [
        np.array([1.541959469,
    1.728658408,
    1.915357347,
    2.102056286,
    2.288755225,
    2.475454163,
    2.662153102,
    2.848852041,
    3.03555098,
    3.222249919,
    3.408948858,
    3.595647796,
    3.782346735,
    3.969045674,
    4.155744613,
    4.342443552,
    4.52914249,
    4.715841429,
    4.902540368,
    5.089239307,])*10/2,
        [14,
    32,
    65,
    55,
    61,
    67,
    46,
    40,
    26,
    21,
    11,
    6,
    8,
    1,
    2,
    2,
    0,
    0,
    0,
    1,]
    ]#
    [[24.2/2],[1]] # this is the mean size from fitting the TEM

    return (
        A0,
        R_max,
        R_np_md,
        Rs_np_D,
        Rs_np_H,
        T_m,
        bulk_dodecanethiol_density,
        bulk_toluene_density,
        eta,
        g_iso,
        len1,
        len2,
        md_bin_size,
        nu_n,
        rho_decay,
        sigma,
        target_dr,
        tau,
    )


@app.cell
def fit(
    A0,
    Model,
    R_max,
    R_np_md,
    Rs_np_D,
    Rs_np_H,
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

    triple_spectra_model = Model(function_for_iso_fit, independent_vars = ["xs", "len1", "len2", "solvent_density_obj", "rho_bulk_solvent", "ligand_density_obj", "rho_bulk_ligand", "R_np_md", "g_iso", "tau", "T_m", "target_dr", "shape", "Rs_np_H", "Rs_np_D"])
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
        ("sigma", sigma, True, 0, None),
        ("eta", eta, False, 0, None),
        ('R_max', R_max, False, R_np_md, None),
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
        R_np_md = R_np_md,
        g_iso = g_iso,
        tau = tau,
        T_m = T_m,
        target_dr = target_dr,
        shape = "G",
        Rs_np_H = Rs_np_H,
        Rs_np_D = Rs_np_D,
        iter_cb=my_callback
    )

    return (triple_spetra_fit,)


@app.cell
def fit_report(triple_spetra_fit):
    print(triple_spetra_fit.fit_report())
    return


@app.cell
def fit_plot(
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

    return


@app.cell
def saving_fit_info(
    dodecaneH_tolueneH_spectra,
    len1,
    len2,
    make_subplots,
    np,
    triple_spetra_fit,
):
    compare_normalized_fitted_plot = make_subplots()

    compare_normalized_fitted_plot.add_scatter(x = dodecaneH_tolueneH_spectra[0], 
                                               y=(triple_spetra_fit.best_fit[:len1] - triple_spetra_fit.params["c1"].value) /np.max(triple_spetra_fit.best_fit[:len1] - triple_spetra_fit.params["c1"].value), line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    compare_normalized_fitted_plot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y=(triple_spetra_fit.best_fit[len1:len1+len2] - triple_spetra_fit.params["c2"].value) /np.max(triple_spetra_fit.best_fit[len1:len1+len2] - triple_spetra_fit.params["c2"].value), line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    compare_normalized_fitted_plot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y=(triple_spetra_fit.best_fit[len1+len2:] - triple_spetra_fit.params["c3"].value) /np.max(triple_spetra_fit.best_fit[len1+len2:] - triple_spetra_fit.params["c3"].value), line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    compare_normalized_fitted_plot.update_layout(template = "simple_white")
    compare_normalized_fitted_plot.show()
    return


if __name__ == "__main__":
    app.run()
