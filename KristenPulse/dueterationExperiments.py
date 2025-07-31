import marimo

__generated_with = "0.13.15"
app = marimo.App(
    width="medium",
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
                extended_bin_centers.append(bin_centers[_i])
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

        return np.array(extended_bin_centers), np.array(extended_sim_ligand), np.array(extended_calc_solvent)

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
def _(high_density):
    (max(high_density["experimental spectra"]["H25_H8"][1]) - min(high_density["experimental spectra"]["H25_H8"][1]))/(sum(high_density["distances"]["extended solvent counts"]) + sum(high_density["distances"]["extended ligand counts"]))  

    return


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
def high_density_fit(
    A,
    gamma,
    high_density,
    high_density_alpha,
    make_subplots,
    mu,
    np,
    spectral_labels,
    triple_sample_model,
):
    # now, we can fit the high density stuff...
    high_density_triple_fit_params = triple_sample_model.make_params()

    high_density_triple_fit_params.add_many(
        ("alpha",   high_density_alpha,  True, 0, None), # decay rate for wavefunction
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

    print(high_density_triple_fit_result.fit_report())

    #% plot
    high_density_triple_fit_plot = make_subplots(rows = 3, cols =1)

    for i in range(len(high_density_solvent_H_fractions)):
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
    high_density_triple_fit_plot
    return


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
    return


@app.cell
def _(lorentzian, np):
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
        
        
            part_x = x[spectral_indices[_i]: spectral_indices[i+1]] # slice out the x-values associated with this part
            
        
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
def _():
    # fit a model that fits only the dueterated samples, and then simulates the hydrogen sample from that (allow fitting of baseline and intensity)
    return


if __name__ == "__main__":
    app.run()
