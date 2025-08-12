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
    from LearKristenPulse import extend_bins_based_on_wavefunction, bin_ligand_and_solvent_together


    try:
        os.chdir(r"C:\Users\benle\Documents\GitHub\publications\KristenPulse")
    except:
        pass


    return (
        Model,
        Path,
        bin_ligand_and_solvent_together,
        extend_bins_based_on_wavefunction,
        go,
        make_subplots,
        mo,
        np,
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
    write_binned_ligand_solvent,
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
    bin_width = 0.1
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
def _(Model, model_three_samples_two_exponential):






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
    plot_triple_fit,
):



    plot_triple_fit(high_density_triple_two_exponential_fit_result, low_density_spectral_indices)
    return


@app.cell
def _(
    calc_ligand_solvent_amounts,
    high_density,
    high_density_best_fit,
    low_density,
    low_density_triple_fit_result,
):
    # write a function that will accept the fit decays, and then estimate the relative contributions from ligand and solvent





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
    return


@app.cell
def _(Model, find_ratio_of_ligand_to_solvent, high_density):
    # It could be the case that we could FIRST fit data such that we fid alpha that leads to the correct ratio of ligad to solvent values. Once we have that, then we can use this value fixed, and try to fit pulsed signals. 




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
        [1.5/4.6], # this is the ratio of solvent/ligand intensity
        params = ligand_to_solvent_params,
        bin_centers = high_density["distances"]["raw bin centers"],
        ligand_counts = high_density["distances"]["raw ligand counts"],
        solvent_counts = high_density["distances"]["raw solvent counts"],
        )

    print(ligand_to_solvent_result.fit_report())
    return (ligand_to_solvent_result,)


@app.cell
def _(
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
        ("alpha",   ligand_to_solvent_result.params["alpha"].value/1,  True, 0, None), # decay rate for wavefunction
        ("A",       5,      True, 0.1, None), # a scalar value that adjust how much the wavefunction overlap leads to coupling

        ("gamma",   gamma/5,  True, 0, gamma*5), # width of lines associated with each contribution

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
