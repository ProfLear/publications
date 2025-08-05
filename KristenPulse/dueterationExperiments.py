import marimo

__generated_with = "0.13.15"
app = marimo.App(
    width="full",
    layout_file="layouts/dueterationExperiments.grid.json",
)


@app.cell
def _(mo):
    mo.md(r"""# Record of how MIMMS is fit""")
    return


@app.cell
def _():
    # imports and data paths
    import marimo as mo
    import numpy as np
    from pathlib import Path
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import codechembook.quickPlots as qp
    import os
    from lmfit import Model


    try:
        os.chdir(r"C:\Users\benle\Documents\GitHub\publications\KristenPulse")
    except:
        pass


    return Model, Path, go, make_subplots, mo, np


@app.cell
def _(bin_width, np):
    # functions to support the notebook
    def lorentzian(x, A = 1, mu = 1, gamma = 1):
        if isinstance(x, list):
            x = np.array(x)
        return (A / np.pi) * (gamma / ((x - mu)**2 + gamma**2))

    def bin_ligand_and_solvent_together(ligand, solvent, bin_width = 1):
        left_edge = 0
        bin_centers = []
        ligand_counts = []
        solvent_counts = []

        # first, bin the data we have
        while left_edge < max(max(ligand), max(solvent)):
            bin_centers.append(left_edge + bin_width/2)

            ligand_counts.append(np.sum((ligand >= left_edge) & (ligand < left_edge + bin_width)))

            solvent_counts.append(np.sum((solvent >= left_edge) & (solvent < left_edge + bin_width)))

            left_edge = left_edge + bin_width

        return np.array(bin_centers), np.array(ligand_counts), np.array(solvent_counts)

    def extend_bins_based_on_wavefunction(bin_centers, ligand_counts, solvent_counts, alpha = 1, fract_threshold = 0.01):
        extended_bin_centers = []
        extended_sim_solvent = []
        extended_sim_ligand = []
        extended_calc_solvent = []

        latest_fract = 1
        i = 0
        found_ligand = False
        extend = False # flag to keep track of if we have switched to extending the solvent 
        while latest_fract > fract_threshold: # basically go until the contribution from the most recent bin is less than some percent
            # check to see if we should change the extend flag
            if extend:

                new_solvent_value = extended_calc_solvent[-1] * ((extended_bin_centers[-1] + bin_width)**3 - extended_bin_centers[-1]**3)/(extended_bin_centers[-1]**3 - (extended_bin_centers[-1] - bin_width)**3)  

                extended_bin_centers.append(extended_bin_centers[-1] + bin_width)
                extended_calc_solvent.append(new_solvent_value)
                try:
                    extended_sim_ligand.append(ligand_counts[i])
                except:
                    extended_sim_ligand.append(0)
                try:
                    extended_sim_solvent.append(solvent_counts[i])
                except:
                    extended_sim_solvent.append(0)


            else: # we are not yet extending the data
                extended_bin_centers.append(bin_centers[i])
                extended_calc_solvent.append(solvent_counts[i])
                extended_sim_solvent.append(solvent_counts[i])
                extended_sim_ligand.append(ligand_counts[i])

                if ligand_counts[i] > 0:
                    found_ligand = True


                if i == len(bin_centers)-1: # we are are the end of the simulated data, and we MUST switch
                    extend = True
                elif i == len(ligand_counts) and found_ligand: # we have found the end of the ligand, and so we should follow the expected r**2 behavior now
                    extend = True
                elif solvent_counts[i] == max(solvent_counts): # we have found the maximum solvent value, which should then just increase as r**2
                    extend = True

            # now, using the decay function, figure out the fraction that this latest value represents...

            ligand_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_sim_ligand))
            solvent_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_calc_solvent))

            if ligand_contribution + solvent_contribution > 0 and extend: # use of extend ensure we have found the end of the ligand?
                latest_fract = np.exp(-2*alpha * extended_bin_centers[-1]) * extended_calc_solvent[-1] / (ligand_contribution + solvent_contribution)
            else:
                latest_fract = 1
            i = i + 1

        return np.array(extended_bin_centers), np.array(extended_sim_ligand), np.round(np.array(extended_calc_solvent)).astype(int)

    def bin_simulation_distances():
        # bin the distaances into particular widths
        return

    def extend_solvent_bins():
        # work from existing bins and xxx to extend the bins
        return


    return (
        bin_ligand_and_solvent_together,
        extend_bins_based_on_wavefunction,
        lorentzian,
    )


@app.cell
def _(Path, np):
    # import the simulation data
    high_density = {}
    high_density["solvent distances"] = np.genfromtxt(Path(r".\MDSimulationResults\toluene_distances_high_density.csv"), delimiter = ",")
    high_density["ligand distances"] = np.genfromtxt(Path(r".\MDSimulationResults\dodecanol_distances_high_density.csv"), delimiter = ",")

    # then do the low density stuff
    low_density = {}
    low_density["solvent distances"] = np.genfromtxt(Path(r".\MDSimulationResults\toluene_distances_low_density.csv"), delimiter  = ",")
    low_density["ligand distances"] = np.genfromtxt(Path(r".\MDSimulationResults\dodecanol_distances_low_density.csv"), delimiter = ",")


    return high_density, low_density


@app.cell
def _(high_density, low_density, np):
    # import the spectral data

    high_density["experimental spectra"] = {}
    high_density["experimental spectra"]["H25_H8"] = np.genfromtxt(r".\MIMS\PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)
    high_density["experimental spectra"]["D25_H8"] = np.genfromtxt(r".\MIMS\dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter=",", unpack = True, skip_header = 1)
    high_density["experimental spectra"]["H25_D8"] = np.genfromtxt(r".\MIMS\Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter=",", unpack = True, skip_header = 1)

    low_density["experimental spectra"] = {}
    low_density["experimental spectra"]["H25_H8"] = np.genfromtxt(r".\MIMS\PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)
    low_density["experimental spectra"]["D25_H8"] = np.genfromtxt(r".\MIMS\dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter=",", unpack = True, skip_header = 1)
    low_density["experimental spectra"]["H25_D8"] = np.genfromtxt(r".\MIMS\Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter=",", unpack = True, skip_header = 1)

    return


@app.cell
def _(high_density, make_subplots):
    # plot experimental spectrum
    spectral_labels = [
        "ligand hydrogen, solvent hydrogen",
        "ligand dueterium, solvent hydrogen",
        "ligand hydrogen, solvent dueterium"
    ]

    experimental_spectra_plot = make_subplots(rows = 3, cols = 1)
    for _i, _spectrum in enumerate([high_density["experimental spectra"]["H25_H8"], high_density["experimental spectra"]["D25_H8"], high_density["experimental spectra"]["H25_D8"]]):
        experimental_spectra_plot.add_scatter(
            x = _spectrum[0],
            y = _spectrum[1],
            line = dict(color = "lightgrey", width = 6),
            showlegend = False,
        row = _i + 1, col = 1,
        )
        experimental_spectra_plot.add_annotation(
            text = spectral_labels[_i],
            x = _spectrum[0][int(len(_spectrum[0])/2)],
            y = max(_spectrum[1]) + 0.15 * (max(_spectrum[1]) - min(_spectrum[1])), 
            showarrow =False,
            row = _i + 1, col = 1
        )
    experimental_spectra_plot
    return (spectral_labels,)


@app.function
def write_binned_ligand_solvent(centers, ligand, solvent, name = "binned ligand and solvent data"):
    with open(name, "w") as f:
        f.write("bin centers, ligand counts, solvent counts \n")
        for c, l, s in zip(centers, ligand, solvent):
            f.write(f"{c}, {l}, {s}" + "\n")

    print("file written")
    return None


@app.cell
def _(
    bin_ligand_and_solvent_together,
    bin_width,
    high_density,
    low_density,
    make_subplots,
    mo,
):
    # bin the simulation data

    raw_binned_distance_hists = make_subplots(rows = 2, cols = 1, shared_xaxes = True)

    for _i, _dictionary in enumerate([high_density, low_density]):
        _dictionary["distances"] = {}
        _dictionary["distances"]["raw bin centers"], _dictionary["distances"]["raw ligand counts"], _dictionary["distances"]["raw solvent counts"] =  bin_ligand_and_solvent_together(_dictionary["ligand distances"], _dictionary["solvent distances"], bin_width = bin_width)

        raw_binned_distance_hists.add_bar(
            x = _dictionary["distances"]["raw bin centers"], y = _dictionary["distances"]["raw ligand counts"], 
            marker = dict(color = "orange", opacity = 0.33),
            showlegend = False,
            row = _i+1, col = 1,
            )
        raw_binned_distance_hists.add_bar(
            x = _dictionary["distances"]["raw bin centers"], y = _dictionary["distances"]["raw solvent counts"], 
            marker = dict(color = "blue", opacity = 0.33), 
            showlegend = False,
            row = _i+1, col = 1,
            )
    raw_binned_distance_hists.update_xaxes(title = "distance / Angstroms")
    raw_binned_distance_hists.update_layout(bargap = 0, barmode = "overlay")
    #raw_binned_distance_hists.show()

    mo.ui.plotly(raw_binned_distance_hists)

    return (raw_binned_distance_hists,)


@app.cell
def _(
    extend_bins_based_on_wavefunction,
    fract_threshold,
    go,
    high_density,
    high_density_alpha,
    low_density,
    low_density_alpha,
    mo,
    np,
    raw_binned_distance_hists,
):
    # now, extend the data, based on the anticipated contribution, and plot wavefunction and wavefunction squared on histograms

    hist_with_wavefunction_plot = go.Figure(raw_binned_distance_hists)

    for _i, _list in enumerate([[high_density, high_density_alpha], [low_density, low_density_alpha]]):
        #print(len(_list))
        _dictionary, _alpha = _list
        _dictionary["distances"]["extended bin centers"], _dictionary["distances"]["extended ligand counts"], _dictionary["distances"]["extended solvent counts"] =  extend_bins_based_on_wavefunction(_dictionary["distances"]["raw bin centers"], _dictionary["distances"]["raw ligand counts"], _dictionary["distances"]["raw solvent counts"], fract_threshold = fract_threshold, alpha = _alpha)

        hist_with_wavefunction_plot.add_bar(x = _dictionary["distances"]["extended bin centers"], y = _dictionary["distances"]["extended solvent counts"], 
                                            marker = dict(color = "blue", opacity = 0.33),
                                           row = _i + 1, col = 1, 
                                           showlegend = False,
                                           )

        hist_with_wavefunction_plot.add_scatter(x =  _dictionary["distances"]["extended bin centers"], 
                                               y =  np.max(_dictionary["distances"]["extended solvent counts"]) * np.exp(-_alpha*_dictionary["distances"]["extended bin centers"]),
                                                mode = "lines",
                                                line = dict(color = "red"),
                                                showlegend = False,
                                                row = _i +1, col = 1,
                                               )

        hist_with_wavefunction_plot.add_scatter(x =  _dictionary["distances"]["extended bin centers"], 
                                               y =  np.max(_dictionary["distances"]["extended solvent counts"]) * np.exp(-2*_alpha*_dictionary["distances"]["extended bin centers"]),
                                                mode = "lines",
                                                line = dict(color = "red", dash = "dash"),
                                                showlegend = False,
                                                row = _i +1, col = 1,
                                               )

    #hist_with_wavefunction_plot.show()

    mo.ui.plotly(hist_with_wavefunction_plot)


    # write the high density data for kristen's stuff
    write_binned_ligand_solvent(high_density["distances"]["raw bin centers"], high_density["distances"]["extended ligand counts"], high_density["distances"]["extended solvent counts"], name = "binned high density stuff.csv")

    return


@app.cell
def _():
    # settings and parameters
    bin_width = 1
    fract_threshold = 0.01

    high_density_alpha = 0.05
    low_density_alpha = 0.3
    A = 1
    mu = 50
    gamma = 0.1


    return (
        A,
        bin_width,
        fract_threshold,
        gamma,
        high_density_alpha,
        low_density_alpha,
        mu,
    )


@app.cell
def _(Model, model_three_samples):
    # fit data to a model that fits all three spectra at the same time
    # Explicitly list the names of the arguments that are fitting parameters
    triple_sample_param_names = [ # these are the ones that we want to adjust...              
                    "alpha",
                    "A",

                    "gamma", 

                    "I1",
                    "I2",
                    "I3",

                    "mu1",
                    "mu2",
                    "mu3",

                    "a1",
                    "a2",
                    "a3",
                    ]

    triple_sample_model = Model(
        model_three_samples,
        independent_vars=['x'],
        param_names=triple_sample_param_names,
        )



    return (triple_sample_model,)


@app.cell
def _():
    '''
    The fits are pretty sensitive to initial gueeses, and so I should do a grid search on them...
    Turn fitting into a function, fit over a grid of values, and then select the ones with the best AIC and then show them??
    '''
    return


@app.cell
def high_density_fit(
    gamma,
    high_density,
    make_subplots,
    mu,
    np,
    spectral_labels,
    triple_sample_model,
):
    def search_parameter_grid(alpha_endpoints = None, A_endpoints = None, grid_resolution = 10, space = "linear"):
        if alpha_endpoints[0] != alpha_endpoints[-1]:
            #build grid
            if space == "linear":
                alphas = np.linspace(alpha_endpoints[0], alpha_endpoints[1], grid_resolution)
            elif space == "log":
                alphas = np.logspace(alpha_endpoints[0], alpha_endpoints[1], grid_resolution)
            else:
                raise "You need to specify a valid space"
        else:
            alphas = [alpha_endpoints[0]]

        if A_endpoints[0] != A_endpoints[1]:
            if space == "linear":
                As = np.linspace(A_endpoints[0], A_endpoints[1], grid_resolution)
            elif space == "log":
                As = np.logspace(A_endpoints[0], A_endpoints[1], grid_resolution)
            else:
                raise "You need to specify a valid space"
        else:
            As = [A_endpoints[0]]

        collected_fits = {}

        i = 0
        for alpha in alphas:
            for A in As:
                collected_fits[f"alpha = {alpha}, A = {A}"] = {}
                collected_fits[f"alpha = {alpha}, A = {A}"]["fit result"], collected_fits[f"alpha = {alpha}, A = {A}"]["fit plot"] = fit_mimms_spectra(alpha = alpha, A = A)
                i = i+1
                print(f"fit {i} / {len(alphas) * len(As)} completed. alpha = {alpha}, A = {A}")
        return collected_fits







    def fit_mimms_spectra(alpha, A):
        # now, we can fit the high density stuff...
        high_density_triple_fit_params = triple_sample_model.make_params()

        high_density_triple_fit_params.add_many(
            ("alpha",   alpha,  True, 0, None), # decay rate for wavefunction
            ("A",       A,      True, 0, None), # a scalar value that adjust how much the wavefunction overlap leads to coupling

            ("gamma",   gamma,  True, 0, None), # width of lines associated with each contribution

            ("I1",       (max(high_density["experimental spectra"]["H25_H8"][1]) - min(high_density["experimental spectra"]["H25_H8"][1]))/(sum(high_density["distances"]["extended solvent counts"]) + sum(high_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the first spectrum
            ("I2",       (max(high_density["experimental spectra"]["D25_H8"][1]) - min(high_density["experimental spectra"]["D25_H8"][1]))/(sum(high_density["distances"]["extended solvent counts"]))  ,      True, 0, None), # intensity adjuster for the second spectrum
            ("I3",       (max(high_density["experimental spectra"]["H25_D8"][1]) - min(high_density["experimental spectra"]["H25_D8"][1]))/(sum(high_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the third spectrum

            ("mu1",      mu,     True, 48, 52), # center position for the first spectrum
            ("mu2",      mu,     True, 48, 52), # center position for the second spectrum <-- needed?
            ("mu3",      mu,     True, 48, 52), # center position for the second spectrum <--needed?

            ("a1",       min(high_density["experimental spectra"]["H25_H8"][1]),      True, None, None),  # baseline for the 1st spectrum
            ("a2",       min(high_density["experimental spectra"]["D25_H8"][1]),      True, None, None),
            ("a3",       min(high_density["experimental spectra"]["H25_D8"][1]),      True, None, None),

            )


        # now set up to fit them

        high_density_concatenated_xs = np.array([])
        high_density_concatenated_ys = np.array([])
        high_density_spectral_indices = [0]
        for _i, _spectrum in enumerate([high_density["experimental spectra"]["H25_H8"], high_density["experimental spectra"]["D25_H8"], high_density["experimental spectra"]["H25_D8"]]):
            high_density_concatenated_xs = np.concatenate([high_density_concatenated_xs, _spectrum[0]])
            high_density_concatenated_ys = np.concatenate([high_density_concatenated_ys, _spectrum[1]])
            high_density_spectral_indices.append(len(high_density_concatenated_xs)) # append the curent length of the growing array

        high_density_ligand_H_fractions = [1, 0.02, 1]
        high_density_solvent_H_fractions = [1, 1, 0.01]

        high_density_triple_fit_result = triple_sample_model.fit(
            high_density_concatenated_ys,

            x = high_density_concatenated_xs,

            params = high_density_triple_fit_params,

            ligand_binned_centers = high_density["distances"]["extended bin centers"], # to be supplied as not fit
            solvent_binned_centers = high_density["distances"]["extended bin centers"], # to be supplied as not fit

            ligand_binned_counts = high_density["distances"]["extended ligand counts"], # to be supplied as not fit
            solvent_binned_counts = high_density["distances"]["extended solvent counts"], # to be supplied as not fit

            ligand_H_fractions =  high_density_ligand_H_fractions, # to be supplied as not fit
            solvent_H_fractions = high_density_solvent_H_fractions, # to be supplied as not fit

            spectral_indices = high_density_spectral_indices,
            )


        # make a plot
        high_density_triple_fit_plot = make_subplots(rows = 3, cols =1)

        for i in range(3):
            high_density_triple_fit_plot.add_scatter(
                x = high_density_concatenated_xs[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
                y = high_density_concatenated_ys[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
                mode = "lines",
                line = dict(width = 8, color = "grey"),
                row = i + 1, col = 1,
                showlegend = False,
                )
            high_density_triple_fit_plot.add_scatter(
                x = high_density_concatenated_xs[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
                y = high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
                line = dict(width = 4, color = "white"),
                row = i + 1, col = 1,
                showlegend = False,
                )
            high_density_triple_fit_plot.add_scatter(
                x = high_density_concatenated_xs[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
                y = high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
                line = dict(width = 2, color = "red"),
                row = i + 1, col = 1,
                showlegend = False,
                )
            high_density_triple_fit_plot.add_annotation(
                text = spectral_labels[i],
                y = max([
                    max(high_density_concatenated_ys[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]]),
                    max(high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]])
                    ]) + (max([
                        max(high_density_concatenated_ys[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]]),
                        max(high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]])
                        ]) - min([
                            min(high_density_concatenated_ys[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]]),
                            min(high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]])
                            ])) * 0.15,
                showarrow = False,
                row =i + 1, col = 1
                )
        high_density_triple_fit_plot.update_xaxes(    
            title = "magnetic field",
            row = 3, col = 1)

        high_density_triple_fit_plot.update_xaxes(
            range = [46, 54],
            )

        high_density_triple_fit_plot.update_layout(
            template = "simple_white",
            margin = dict(l = 50, r = 10, b = 20, t = 10),
            width = 500, height = 600,
            )

        return high_density_triple_fit_result, high_density_triple_fit_plot


    return (search_parameter_grid,)


@app.cell
def _(search_parameter_grid):
    high_density_results = search_parameter_grid(alpha_endpoints = [0, -5], A_endpoints = [0, -10], grid_resolution = 10, space = "log")


    return (high_density_results,)


@app.cell
def _(go, high_density_results, np):
    AICs = []
    AIC_heats = []
    _i = 0
    temp_row = []
    for _entry in high_density_results:
        if _i%10 ==0 and _i != 0:
            AIC_heats.append(temp_row)
            temp_row = []
        AICs.append(high_density_results[_entry]["fit result"].aic)
        temp_row.append(high_density_results[_entry]["fit result"].aic)
        _i = _i + 1

    print(min(AICs))


    fig = go.Figure(data=go.Heatmap(
                        z = AIC_heats,
                        y = np.linspace(1, 0.01, 10),
                        x = np.linspace(10, 0.1, 10)
                        ))

    fig.update_xaxes(title = "A")
    fig.update_yaxes(title = "alpha")
    fig.show()



    return (AICs,)


@app.cell
def _(AICs, high_density_results):
    for _entry in high_density_results:
        if high_density_results[_entry]["fit result"].aic == min(AICs):
            high_density_results[_entry]["fit plot"].show()
            print(high_density_results[_entry]["fit result"].fit_report())
            high_density_best_fit = high_density_results[_entry]["fit result"]

    return (high_density_best_fit,)


@app.cell
def _(
    A,
    gamma,
    low_density,
    low_density_alpha,
    make_subplots,
    mu,
    np,
    spectral_labels,
    triple_sample_model,
):
    # we can also do a low density thing...

    # now, we can fit the high density stuff...
    low_density_triple_fit_params = triple_sample_model.make_params()

    low_density_triple_fit_params.add_many(
        ("alpha",   low_density_alpha,  True, 0, None), # decay rate for wavefunction
        ("A",       A,      True, 0, None), # a scalar value that adjust how much the wavefunction overlap leads to coupling

        ("gamma",   gamma,  True, 0, None), # width of lines associated with each contribution

        ("I1",       (max(low_density["experimental spectra"]["H25_H8"][1]) - min(low_density["experimental spectra"]["H25_H8"][1]))/(sum(low_density["distances"]["extended solvent counts"]) + sum(low_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the first spectrum
        ("I2",       (max(low_density["experimental spectra"]["D25_H8"][1]) - min(low_density["experimental spectra"]["D25_H8"][1]))/(sum(low_density["distances"]["extended solvent counts"]))  ,      True, 0, None), # intensity adjuster for the second spectrum
        ("I3",       (max(low_density["experimental spectra"]["H25_D8"][1]) - min(low_density["experimental spectra"]["H25_D8"][1]))/(sum(low_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the third spectrum

        ("mu1",      mu,     True, 48, 52), # center position for the first spectrum
        ("mu2",      mu,     True, 48, 52), # center position for the second spectrum <-- needed?
        ("mu3",      mu,     True, 48, 52), # center position for the second spectrum <--needed?

        ("a1",       min(low_density["experimental spectra"]["H25_H8"][1]),      True, None, None),  # baseline for the 1st spectrum
        ("a2",       min(low_density["experimental spectra"]["D25_H8"][1]),      True, None, None),
        ("a3",       min(low_density["experimental spectra"]["H25_D8"][1]),      True, None, None),

        )


    # now set up to fit them

    low_density_concatenated_xs = np.array([])
    low_density_concatenated_ys = np.array([])
    low_density_spectral_indices = [0]
    for _i, _spectrum in enumerate([low_density["experimental spectra"]["H25_H8"], low_density["experimental spectra"]["D25_H8"], low_density["experimental spectra"]["H25_D8"]]):
        low_density_concatenated_xs = np.concatenate([low_density_concatenated_xs, _spectrum[0]])
        low_density_concatenated_ys = np.concatenate([low_density_concatenated_ys, _spectrum[1]])
        low_density_spectral_indices.append(len(low_density_concatenated_xs)) # append the curent length of the growing array

    low_density_ligand_H_fractions = [1, 0.02, 1]
    low_density_solvent_H_fractions = [1, 1, 0.01]

    low_density_triple_fit_result = triple_sample_model.fit(
        low_density_concatenated_ys,

        x = low_density_concatenated_xs,

        params = low_density_triple_fit_params,

        ligand_binned_centers = low_density["distances"]["extended bin centers"], # to be supplied as not fit
        solvent_binned_centers = low_density["distances"]["extended bin centers"], # to be supplied as not fit

        ligand_binned_counts = low_density["distances"]["extended ligand counts"], # to be supplied as not fit
        solvent_binned_counts = low_density["distances"]["extended solvent counts"], # to be supplied as not fit

        ligand_H_fractions =  low_density_ligand_H_fractions, # to be supplied as not fit
        solvent_H_fractions = low_density_solvent_H_fractions, # to be supplied as not fit

        spectral_indices = low_density_spectral_indices,
        )

    print(low_density_triple_fit_result.fit_report())

    #% plot
    triple_fit_plot = make_subplots(rows = 3, cols =1)

    for _i in range(len(low_density_solvent_H_fractions)):
        triple_fit_plot.add_scatter(
            x = low_density_concatenated_xs[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            y = low_density_concatenated_ys[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            mode = "lines",
            line = dict(width = 8, color = "grey"),
            row =_i+ 1, col = 1,
            showlegend = False,
            )
        triple_fit_plot.add_scatter(
            x = low_density_concatenated_xs[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            y = low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            line = dict(width = 4, color = "white"),
            row =_i+ 1, col = 1,
            showlegend = False,
            )
        triple_fit_plot.add_scatter(
            x = low_density_concatenated_xs[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            y = low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            line = dict(width = 2, color = "red"),
            row =_i+ 1, col = 1,
            showlegend = False,
            )
        triple_fit_plot.add_annotation(
            text = spectral_labels[_i],
            y = max([
                max(low_density_concatenated_ys[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]]),
                max(low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]])
                ]) + (max([
                    max(low_density_concatenated_ys[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]]),
                    max(low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]])
                    ]) - min([
                        min(low_density_concatenated_ys[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]]),
                        min(low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]])
                        ])) * 0.15,
            showarrow = False,
            row =_i + 1, col = 1
            )
    triple_fit_plot.update_xaxes(    
        title = "magnetic field",
        row = 3, col = 1)

    triple_fit_plot.update_xaxes(
        range = [46, 54],
        )

    triple_fit_plot.update_layout(
        template = "simple_white",
        margin = dict(l = 50, r = 10, b = 20, t = 10),
        width = 500, height = 600,
        )
    triple_fit_plot
    return (
        low_density_concatenated_xs,
        low_density_concatenated_ys,
        low_density_ligand_H_fractions,
        low_density_solvent_H_fractions,
        low_density_spectral_indices,
        low_density_triple_fit_result,
    )


@app.cell
def _(lorentzian, np):
    # this one is based on a single decay rate (i..e, through space)
    def simulate_pulse_signal(
            x,

            ligand_binned_centers, # to be supplied as not fit
            solvent_binned_centers, # to be supplied as not fit

            ligand_binned_counts, # to be supplied as not fit
            solvent_binned_counts, # to be supplied as not fit

            ligand_H_fraction, # to be supplied as not fit
            solvent_H_fraction, # to be supplied as not fit

            alpha, # decay constant
            A, # scalar value to adjust strength of coupling

            I, # intensity of individual contributions
            mu, # central position, without coupling
            gamma, # width of the lorentzian peak

            a, # baseline correction
            ):

        # for debudgin...
        #print(f"alpha={alpha:.4f}, A={A:.4f}, I={I:.4f}, mu={mu:.4f}, gamma={gamma:.4f}, a={a:.4f}")


        # calculate the psi-squared values at each distance for ligand and solvent
        ligand_psi = np.exp(-alpha * ligand_binned_centers) #included in case we decide we need it
        solvent_psi = np.exp(-alpha * solvent_binned_centers)
        ligand_psi_squareds  = np.exp(-2*alpha * ligand_binned_centers)
        solvent_psi_squareds = np.exp(-2*alpha * solvent_binned_centers)

        # now we are going to sum up contributions
        temp_y_values = np.zeros_like(x)
        for info in [[ligand_psi, ligand_psi_squareds, ligand_binned_counts, ligand_H_fraction], [solvent_psi, solvent_psi_squareds, solvent_binned_counts, solvent_H_fraction]]:
            psis, psi_squareds, counts, H = info # unpack these
            for p, ps, count in zip(psis, psi_squareds, counts):
                temp_y_values = temp_y_values + I*lorentzian(x, mu = mu +   ps*A, gamma = gamma) * H*ps * count
                temp_y_values = temp_y_values + I*lorentzian(x, mu = mu - 1*ps*A, gamma = gamma) * H*ps * count

        return temp_y_values + a



    def model_three_samples(
            x, # x values. These are the concatenated x-values from the experimental spectra we want to fit

            ligand_binned_centers, # to be supplied as not fit
            solvent_binned_centers, # to be supplied as not fit

            ligand_binned_counts, # to be supplied as not fit
            solvent_binned_counts, # to be supplied as not fit

            alpha, # decay constant
            A, # scalar value to adjust strength of coupling

            gamma, # width of the lorentzian peak, shared

            I1, # intensity for 1st spectrum
            I2, # intensity for 2nd spectrum
            I3, # intensity for 3rd spectrum

            mu1, # central position, without coupling, shared
            mu2, # central position, without coupling, shared
            mu3, # central position, without coupling, shared

            a1, # baseline for 1st spectrum
            a2, # baseline for 2nd spectrum
            a3, # baseline for 3rd spectrum

            ligand_H_fractions, # list of ligand H-fractions for each sample
            solvent_H_fractions, # list of solvent H-fractions for each sample
            spectral_indices, # indices associated with the x-values for each spectrum. n+1, where n is the number of samples

            ):

        total_ys = np.array([])

        for I, mu, a, lH, sH, i in zip(
                [I1, I2, I3],
                [mu1, mu2, mu3],
                [a1, a2, a3],
                ligand_H_fractions,
                solvent_H_fractions,
                range(len(ligand_H_fractions)) # get indices for each sample (i.e., for each solvent fraction)
                ):


            part_x = x[spectral_indices[i]: spectral_indices[i+1]] # slice out the x-values associated with this part


            part_y_sim = simulate_pulse_signal(
                    part_x,

                    ligand_binned_centers, # to be supplied as not fit
                    solvent_binned_centers, # to be supplied as not fit

                    ligand_binned_counts, # to be supplied as not fit
                    solvent_binned_counts, # to be supplied as not fit

                    lH, # to be supplied as not fit
                    sH, # to be supplied as not fit

                    alpha, # decay constant
                    A, # scalar value to adjust strength of coupling

                    I, # intensity of individual contributions
                    mu, # central position, without coupling
                    gamma, # width of the lorentzian peak

                    a, # baseline correction
                    )

            total_ys = np.concatenate((total_ys, part_y_sim))
        return total_ys
    return (model_three_samples,)


@app.cell
def _(Model, lorentzian, np):
    # handling through bond and through space

    #need a fucntion that will handle through space and through bond together
    def simulate_pulse_signal_two_exponential(
            x,

            ligand_binned_centers, # to be supplied as not fit
            solvent_binned_centers, # to be supplied as not fit

            ligand_binned_counts, # to be supplied as not fit
            solvent_binned_counts, # to be supplied as not fit

            ligand_H_fraction, # to be supplied as not fit
            solvent_H_fraction, # to be supplied as not fit

            alpha_through_space, # decay constant
            alpha_through_bond, # through bond decay rate
            f_tb, # this is the fraction of the wavefunction that goes through bond.
            A, # scalar value to adjust strength of coupling

            I, # intensity of individual contributions
            mu, # central position, without coupling
            gamma, # width of the lorentzian peak

            a, # baseline correction
            ):


        # for debudgin...
        #print(f"alpha={alpha:.4f}, A={A:.4f}, I={I:.4f}, mu={mu:.4f}, gamma={gamma:.4f}, a={a:.4f}")


        # calculate the psi-squared values at each distance for ligand and solvent
        ligand_psi_through_bond = np.exp(-alpha_through_bond * ligand_binned_centers) 
        ligand_psi_through_space = np.exp(-alpha_through_space * ligand_binned_centers) 
        solvent_psi_through_space = np.exp(-alpha_through_space * solvent_binned_centers)


        liagnd_psi_squareds_through_bond = np.exp(-2*alpha_through_bond * ligand_binned_centers)
        ligand_psi_squareds_through_space  = np.exp(-2*alpha_through_space * ligand_binned_centers)
        solvent_psi_squareds_through_space = np.exp(-2*alpha_through_space * solvent_binned_centers)

        # now we are going to sum up contributions for through space
        temp_y_values_through_space = np.zeros_like(x)
        for info in [
            [ligand_psi_through_space, ligand_psi_squareds_through_space, ligand_binned_counts, ligand_H_fraction], 
            [solvent_psi_through_space, solvent_psi_squareds_through_space, solvent_binned_counts, solvent_H_fraction]
        ]:
            psis, psi_squareds, counts, H = info # unpack these
            for p, ps, count in zip(psis, psi_squareds, counts):
                temp_y_values_through_space = temp_y_values_through_space + I*lorentzian(x, mu = mu +   ps*A, gamma = gamma) * H*ps * count
                temp_y_values_through_space = temp_y_values_through_space + I*lorentzian(x, mu = mu - 1*ps*A, gamma = gamma) * H*ps * count

        # now get the through bond added in
        temp_y_values_through_bond = np.zeros_like(x)
        for info in [
            [ligand_psi_through_bond, liagnd_psi_squareds_through_bond, ligand_binned_counts, ligand_H_fraction], 
        ]:
            psis, psi_squareds, counts, H = info # unpack these
            for p, ps, count in zip(psis, psi_squareds, counts):
                temp_y_values_through_bond = temp_y_values_through_bond + I*lorentzian(x, mu = mu +   ps*A, gamma = gamma) * H*ps * count
                temp_y_values_through_bond = temp_y_values_through_bond + I*lorentzian(x, mu = mu - 1*ps*A, gamma = gamma) * H*ps * count

        return (1-f_tb)*temp_y_values_through_space + f_tb*temp_y_values_through_bond + a

    #function for model
    def model_three_samples_two_exponential(
            x, # x values. These are the concatenated x-values from the experimental spectra we want to fit

            ligand_binned_centers, # to be supplied as not fit
            solvent_binned_centers, # to be supplied as not fit

            ligand_binned_counts, # to be supplied as not fit
            solvent_binned_counts, # to be supplied as not fit

            alpha_through_space, # decay constant
            alpha_through_bond,
            f_tb,
            A, # scalar value to adjust strength of coupling

            gamma, # width of the lorentzian peak, shared

            I1, # intensity for 1st spectrum
            I2, # intensity for 2nd spectrum
            I3, # intensity for 3rd spectrum

            mu1, # central position, without coupling, shared
            mu2, # central position, without coupling, shared
            mu3, # central position, without coupling, shared

            a1, # baseline for 1st spectrum
            a2, # baseline for 2nd spectrum
            a3, # baseline for 3rd spectrum

            ligand_H_fractions, # list of ligand H-fractions for each sample
            solvent_H_fractions, # list of solvent H-fractions for each sample
            spectral_indices, # indices associated with the x-values for each spectrum. n+1, where n is the number of samples

            ):

        total_ys = np.array([])

        for I, mu, a, lH, sH, i in zip(
                [I1, I2, I3],
                [mu1, mu2, mu3],
                [a1, a2, a3],
                ligand_H_fractions,
                solvent_H_fractions,
                range(len(ligand_H_fractions)) # get indices for each sample (i.e., for each solvent fraction)
                ):


            part_x = x[spectral_indices[i]: spectral_indices[i+1]] # slice out the x-values associated with this part


            part_y_sim = simulate_pulse_signal_two_exponential(
                    part_x,

                    ligand_binned_centers, # to be supplied as not fit
                    solvent_binned_centers, # to be supplied as not fit

                    ligand_binned_counts, # to be supplied as not fit
                    solvent_binned_counts, # to be supplied as not fit

                    lH, # to be supplied as not fit
                    sH, # to be supplied as not fit

                    alpha_through_space, # decay constant
                    alpha_through_bond,
                    f_tb,
                    A, # scalar value to adjust strength of coupling

                    I, # intensity of individual contributions
                    mu, # central position, without coupling
                    gamma, # width of the lorentzian peak

                    a, # baseline correction
                    )

            total_ys = np.concatenate((total_ys, part_y_sim))
        return total_ys

    




    # make model

    # Explicitly list the names of the arguments that are fitting parameters
    triple_sample_two_exponential_param_names = [ # these are the ones that we want to adjust...              
                    "alpha_through_space",
                    "alpha_through_bond",
                    "f_tb",
                    "A",

                    "gamma", 

                    "I1",
                    "I2",
                    "I3",

                    "mu1",
                    "mu2",
                    "mu3",

                    "a1",
                    "a2",
                    "a3",
                    ]




    three_samples_two_exponential_model = Model(
        model_three_samples_two_exponential,
        independent_vars=['x'],
        param_names=triple_sample_two_exponential_param_names,
        )


    return (three_samples_two_exponential_model,)


@app.cell
def _(
    A,
    gamma,
    high_density,
    high_density_alpha,
    low_density,
    low_density_concatenated_xs,
    low_density_concatenated_ys,
    low_density_ligand_H_fractions,
    low_density_solvent_H_fractions,
    low_density_spectral_indices,
    mu,
    three_samples_two_exponential_model,
):
    # we can also do a low density thing...

    # now, we can fit the high density stuff...
    high_density_triple_two_exponential_fit_params = three_samples_two_exponential_model.make_params()

    high_density_triple_two_exponential_fit_params.add_many(
        ("alpha_through_space",   high_density_alpha/1,  True, 0, None), # decay rate for wavefunction
        ("alpha_through_bond",    high_density_alpha*1, True, 0, None), # faster through bond decay
        ("f_tb",                  0.5, True, 0, 1),
        ("A",       A,      True, 0, None), # a scalar value that adjust how much the wavefunction overlap leads to coupling

        ("gamma",   gamma/2,  True, 0, None), # width of lines associated with each contribution

        ("I1",       (max(high_density["experimental spectra"]["H25_H8"][1]) - min(high_density["experimental spectra"]["H25_H8"][1]))/(sum(high_density["distances"]["extended solvent counts"]) + sum(high_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the first spectrum
        ("I2",       (max(high_density["experimental spectra"]["D25_H8"][1]) - min(high_density["experimental spectra"]["D25_H8"][1]))/(sum(high_density["distances"]["extended solvent counts"]))  ,      True, 0, None), # intensity adjuster for the second spectrum
        ("I3",       (max(high_density["experimental spectra"]["H25_D8"][1]) - min(high_density["experimental spectra"]["H25_D8"][1]))/(sum(high_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the third spectrum

        ("mu1",      mu,     True, 48, 52), # center position for the first spectrum
        ("mu2",      mu,     True, 48, 52), # center position for the second spectrum <-- needed?
        ("mu3",      mu,     True, 48, 52), # center position for the second spectrum <--needed?

        ("a1",       min(low_density["experimental spectra"]["H25_H8"][1]),      True, None, None),  # baseline for the 1st spectrum
        ("a2",       min(low_density["experimental spectra"]["D25_H8"][1]),      True, None, None),
        ("a3",       min(high_density["experimental spectra"]["H25_D8"][1]),      True, None, None),

        )


    # now set up to fit them
    high_density_triple_two_exponential_fit_result = three_samples_two_exponential_model.fit(
        low_density_concatenated_ys,

        x = low_density_concatenated_xs,

        params = high_density_triple_two_exponential_fit_params,

        ligand_binned_centers = high_density["distances"]["extended bin centers"], # to be supplied as not fit
        solvent_binned_centers = high_density["distances"]["extended bin centers"], # to be supplied as not fit

        ligand_binned_counts = high_density["distances"]["extended ligand counts"], # to be supplied as not fit
        solvent_binned_counts = high_density["distances"]["extended solvent counts"], # to be supplied as not fit

        ligand_H_fractions =  low_density_ligand_H_fractions, # to be supplied as not fit
        solvent_H_fractions = low_density_solvent_H_fractions, # to be supplied as not fit

        spectral_indices = low_density_spectral_indices,
        )

    print(high_density_triple_two_exponential_fit_result.fit_report())
    '''
    #% plot
    triple_fit_plot = make_subplots(rows = 3, cols =1)

    for _i in range(len(low_density_solvent_H_fractions)):
        triple_fit_plot.add_scatter(
            x = low_density_concatenated_xs[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            y = low_density_concatenated_ys[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            mode = "lines",
            line = dict(width = 8, color = "grey"),
            row =_i+ 1, col = 1,
            showlegend = False,
            )
        triple_fit_plot.add_scatter(
            x = low_density_concatenated_xs[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            y = low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            line = dict(width = 4, color = "white"),
            row =_i+ 1, col = 1,
            showlegend = False,
            )
        triple_fit_plot.add_scatter(
            x = low_density_concatenated_xs[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            y = low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]],
            line = dict(width = 2, color = "red"),
            row =_i+ 1, col = 1,
            showlegend = False,
            )
        triple_fit_plot.add_annotation(
            text = spectral_labels[_i],
            y = max([
                max(low_density_concatenated_ys[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]]),
                max(low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]])
                ]) + (max([
                    max(low_density_concatenated_ys[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]]),
                    max(low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]])
                    ]) - min([
                        min(low_density_concatenated_ys[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]]),
                        min(low_density_triple_fit_result.best_fit[low_density_spectral_indices[_i]:low_density_spectral_indices[ _i +1]])
                        ])) * 0.15,
            showarrow = False,
            row =_i + 1, col = 1
            )
    triple_fit_plot.update_xaxes(    
        title = "magnetic field",
        row = 3, col = 1)

    triple_fit_plot.update_xaxes(
        range = [46, 54],
        )

    triple_fit_plot.update_layout(
        template = "simple_white",
        margin = dict(l = 50, r = 10, b = 20, t = 10),
        width = 500, height = 600,
        )
    triple_fit_plot
    '''
    return (high_density_triple_two_exponential_fit_result,)


@app.cell
def _(
    high_density_triple_two_exponential_fit_result,
    low_density_concatenated_xs,
    low_density_concatenated_ys,
    make_subplots,
):
    _t = make_subplots()
    _t.add_scatter(
        x = low_density_concatenated_xs,
        y = low_density_concatenated_ys,
        mode = "markers",
        marker = dict(size = 10),
    )
    _t.add_scatter(
        x = low_density_concatenated_xs,
        y = high_density_triple_two_exponential_fit_result.best_fit,
        mode = "markers",
        marker = dict(size = 5),
    )
    _t.show()
    
    
        #high_density_triple_two_exponential_fit_result.best_fit
    return


@app.cell
def _(
    high_density_triple_two_exponential_fit_result,
    low_density_spectral_indices,
    make_subplots,
):
    def plot_triple_fit(fit, indices):
        # this is just for plotting the fits

        labels = [
            "SC<sub>12</sub>H<sub>25</sub> and C<sub>7</sub>H<sub>8</sub>",
            "SC<sub>12</sub>D<sub>25</sub> and C<sub>7</sub>H<sub>8</sub>",
            "SC<sub>12</sub>H<sub>25</sub> and C<sub>7</sub>D<sub>8</sub>",
            ]

        # Just making some variables for convenience
        # First figure out what the independent variable name(s) is(are)
        independent_vars = fit.model.independent_vars

        # The x data has to be the same for all the independent variables, so
        # so get it from the first one in the list for safety
        xdata = fit.userkws[independent_vars[0]]
        ydata = fit.data


        fig = make_subplots(rows = 3, cols =1)
        for i in range(3):
            fig.add_scatter(
                x = xdata[indices[i]:indices[i+1]],
                y = ydata[indices[i]:indices[i+1]],
                mode = "lines",
                line = dict(color = "grey", width = 10),
                showlegend = False,
                row = i+1, col = 1,
            )
            fig.add_scatter(
                x = xdata[indices[i]:indices[i+1]],
                y = fit.best_fit[indices[i]:indices[i+1]],
                 mode = "lines",
                line = dict(color = "white", width = 4),
                showlegend = False,
                row = i+1, col = 1,
            )
            fig.add_scatter(
                x = xdata[indices[i]:indices[i+1]],
                y = fit.best_fit[indices[i]:indices[i+1]],
                mode = "lines",
                line = dict(color = "red", width = 2),
                showlegend = False,
                row = i+1, col = 1,
            )

            fig.add_annotation(
                text = labels[i],
                x = min(xdata[indices[i]:indices[i+1]]),
                y = max(fit.best_fit[indices[i]:indices[i+1]]),
                xanchor  = "left",
                showarrow = False,
                row = i+1, col = 1,
            )
        

        fig.update_xaxes(title = "magnetic field")
        fig.update_yaxes(title = "intensity")
        fig.update_layout(
            template="simple_white",
            margin = dict(t=5, r = 5)
        )
        
        return fig


    plot_triple_fit(high_density_triple_two_exponential_fit_result, low_density_spectral_indices)
    return (plot_triple_fit,)


@app.cell
def _(
    extend_bins_based_on_wavefunction,
    high_density,
    high_density_best_fit,
    low_density,
    low_density_triple_fit_result,
    np,
):
    # write a function that will accept the fit decays, and then estimate the relative contributions from ligand and solvent

    def calc_ligand_solvent_amounts(alpha, bin_centers, ligand_counts, solvent_counts):

        # get the extended solvent.
        extended_bins, extended_ligand, extended_solvent = extend_bins_based_on_wavefunction(bin_centers, ligand_counts, solvent_counts, alpha = 1, fract_threshold = 0.01)
        # then multiply the psi-squared values times the hydrogen bins and then sum
        ligand = np.sum(np.exp(-2*alpha*extended_bins)*extended_ligand)
        solvent = np.sum(np.exp(-2*alpha*extended_bins)*extended_solvent)

        return ligand, solvent





    high_density_ligand_contribution, high_density_solvent_contribution = calc_ligand_solvent_amounts(
        high_density_best_fit.params["alpha"].value, 
        high_density["distances"]["raw bin centers"], 
        high_density["distances"]["raw ligand counts"], 
        high_density["distances"]["raw solvent counts"],
        )

    low_density_ligand_contribution, low_density_solvent_contribution = calc_ligand_solvent_amounts(
        low_density_triple_fit_result.params["alpha"].value, 
        low_density["distances"]["raw bin centers"], 
        low_density["distances"]["raw ligand counts"], 
        low_density["distances"]["raw solvent counts"],
        )

    print(f"For high density, the ratio of ligand to solvent is {high_density_ligand_contribution/min([high_density_ligand_contribution, high_density_solvent_contribution]):0.2f}:{high_density_solvent_contribution/min([high_density_ligand_contribution, high_density_solvent_contribution]):0.2f}")

    print(f"For low density, the ratio of ligand to solvent is {low_density_ligand_contribution/min([low_density_ligand_contribution, low_density_solvent_contribution]):0.2f}:{low_density_solvent_contribution/min([low_density_ligand_contribution, low_density_solvent_contribution]):0.2f}")
    return (calc_ligand_solvent_amounts,)


@app.cell
def _(Model, calc_ligand_solvent_amounts, high_density):
    # It could be the case that we could FIRST fit data such that we fid alpha that leads to the correct ratio of ligad to solvent values. Once we have that, then we can use this value fixed, and try to fit pulsed signals. 


    def find_ratio_of_ligand_to_solvent(
        alpha, 
        bin_centers, 
        ligand_counts, 
        solvent_counts
        ):

        ligand_sum, solvent_sum = calc_ligand_solvent_amounts(alpha, bin_centers, ligand_counts, solvent_counts)
    
        return [solvent_sum/ligand_sum]

    # Explicitly list the names of the arguments that are fitting parameters
    ligand_to_solvent_param_names = [ # these are the ones that we want to adjust...              
                    "alpha",
                    ]



    ligand_to_solvent_model = Model(
        find_ratio_of_ligand_to_solvent,
        independent_vars = ["bin_centers"],
        param_names = ligand_to_solvent_param_names,
        )

    ligand_to_solvent_params = ligand_to_solvent_model.make_params()
    ligand_to_solvent_params.add_many(
        ("alpha", 0.01, True, 0, None)
    )

    ligand_to_solvent_result = ligand_to_solvent_model.fit(
        [1/3.17], # this is the ratio of solvent/ligand intensity
        params = ligand_to_solvent_params,
        bin_centers = high_density["distances"]["raw bin centers"],
        ligand_counts = high_density["distances"]["raw ligand counts"],
        solvent_counts = high_density["distances"]["raw solvent counts"],
        )

    print(ligand_to_solvent_result.fit_report())
    return (ligand_to_solvent_result,)


@app.cell
def _(
    A,
    gamma,
    ligand_to_solvent_result,
    low_density,
    low_density_concatenated_xs,
    low_density_concatenated_ys,
    low_density_ligand_H_fractions,
    low_density_solvent_H_fractions,
    low_density_spectral_indices,
    mu,
    plot_triple_fit,
    triple_sample_model,
):
    triple_sample_model_defined_alpha_params = triple_sample_model.make_params()
    triple_sample_model_defined_alpha_params.add_many(
        ("alpha",   ligand_to_solvent_result.params["alpha"].value,  False, 0, None), # decay rate for wavefunction
        ("A",       A,      True, 0, None), # a scalar value that adjust how much the wavefunction overlap leads to coupling

        ("gamma",   gamma,  True, 0, None), # width of lines associated with each contribution

        ("I1",       (max(low_density["experimental spectra"]["H25_H8"][1]) - min(low_density["experimental spectra"]["H25_H8"][1]))/(sum(low_density["distances"]["extended solvent counts"]) + sum(low_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the first spectrum
        ("I2",       (max(low_density["experimental spectra"]["D25_H8"][1]) - min(low_density["experimental spectra"]["D25_H8"][1]))/(sum(low_density["distances"]["extended solvent counts"]))  ,      True, 0, None), # intensity adjuster for the second spectrum
        ("I3",       (max(low_density["experimental spectra"]["H25_D8"][1]) - min(low_density["experimental spectra"]["H25_D8"][1]))/(sum(low_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the third spectrum

        ("mu1",      mu,     True, 48, 52), # center position for the first spectrum
        ("mu2",      mu,     True, 48, 52), # center position for the second spectrum <-- needed?
        ("mu3",      mu,     True, 48, 52), # center position for the second spectrum <--needed?

        ("a1",       min(low_density["experimental spectra"]["H25_H8"][1]),      True, None, None),  # baseline for the 1st spectrum
        ("a2",       min(low_density["experimental spectra"]["D25_H8"][1]),      True, None, None),
        ("a3",       min(low_density["experimental spectra"]["H25_D8"][1]),      True, None, None),
    )


    triple_sample_model_defined_alpha_result = triple_sample_model.fit(
        low_density_concatenated_ys,

        x = low_density_concatenated_xs,

        params = triple_sample_model_defined_alpha_params,

        ligand_binned_centers = low_density["distances"]["extended bin centers"], # to be supplied as not fit
        solvent_binned_centers = low_density["distances"]["extended bin centers"], # to be supplied as not fit

        ligand_binned_counts = low_density["distances"]["extended ligand counts"], # to be supplied as not fit
        solvent_binned_counts = low_density["distances"]["extended solvent counts"], # to be supplied as not fit

        ligand_H_fractions =  low_density_ligand_H_fractions, # to be supplied as not fit
        solvent_H_fractions = low_density_solvent_H_fractions, # to be supplied as not fit

        spectral_indices = low_density_spectral_indices,
        )

    print(triple_sample_model_defined_alpha_result.fit_report())

    plot_triple_fit(triple_sample_model_defined_alpha_result, low_density_spectral_indices)
    return


@app.cell
def _():
    # fit a model that fits only the dueterated samples, and then simulates the hydrogen sample from that (allow fitting of baseline and intensity)
    return


if __name__ == "__main__":
    app.run()
