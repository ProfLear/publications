import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    from pathlib import Path
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import codechembook.quickPlots as qp
    from lmfit import Model
    from importlib import reload
    import os

    try:
        os.chdir(r"C:\Users\benle\Documents\GitHub\publications\KristenPulse")
    except:
        pass

    import LearKristenPulse as LKP
    #from LearKristenPulse import extend_bins_based_on_wavefunction, bin_ligand_and_solvent_together

    ligand_H_fractions  = [1, 0.02, 1]
    solvent_H_fractions = [1, 1, 0.01]
    return (
        LKP,
        Path,
        ligand_H_fractions,
        make_subplots,
        mo,
        np,
        reload,
        solvent_H_fractions,
    )


@app.cell
def import_distances(Path, np):
    # import data and bin size data

    high_density = {}
    high_density["solvent distances"] = np.genfromtxt(Path(r".\MDSimulationResults\toluene_distances_high_density.csv"), delimiter = ",")
    high_density["ligand distances"] = np.genfromtxt(Path(r".\MDSimulationResults\dodecanol_distances_high_density.csv"), delimiter = ",")

    # then do the low density stuff
    low_density = {}
    low_density["solvent distances"] = np.genfromtxt(Path(r".\MDSimulationResults\toluene_distances_low_density.csv"), delimiter  = ",")
    low_density["ligand distances"] = np.genfromtxt(Path(r".\MDSimulationResults\dodecanol_distances_low_density.csv"), delimiter = ",")


    return (high_density,)


@app.cell
def bin_distances(
    LKP,
    alpha_input,
    bin_width_input,
    contribution_threshold_input,
    high_density,
    reload,
):
    reload(LKP)

    def bin_and_extend_distances(distances, bin_width, alpha, contribution_threshold=0.01):

        binned_distances = LKP.bin_ligand_and_solvent_together(distances[0], distances[1], bin_width = bin_width)
    
        #print(alpha)
        extended_binned_data = LKP.extend_bins_based_on_wavefunction(binned_distances[0], binned_distances[1], binned_distances[2], alpha = alpha, fract_threshold = contribution_threshold)

    
        return [extended_binned_data[0], extended_binned_data[1], extended_binned_data[2], extended_binned_data[3]]

    extended_binned_data = bin_and_extend_distances([high_density["ligand distances"],high_density["solvent distances"]], bin_width = bin_width_input.value, alpha = alpha_input.value, contribution_threshold = contribution_threshold_input.value)
    return (extended_binned_data,)


@app.cell
def plot_distances_and_wavefunction(make_subplots, np):
    # plot the results of the binned distances and wavefunction
    # show just the wavefunction squared. 

    def generate_binned_distance_plot(binned_distances, alpha):

        colors = ["lightblue", "blue", "darkorange"] # sim solvent, raw solvent, raw ligand
        fig = make_subplots(
            rows = 3, # histogram of distances, heistogram of weighted distances, divided bar chart of ligand/solvent contributions
            cols = 1,
            )
        for i, bd in enumerate(reversed(binned_distances[1:])):
            fig.add_bar( # the binned distances
                x = binned_distances[0],
                y = bd,
                marker = dict(color = colors[i], opacity = 0.66),
                showlegend = False,
                row = 1, col = 1,
                )
    
            fig.add_bar( # binned distances, wieghted by psi-squared
                x = binned_distances[0],
                y = bd*np.exp(-2*alpha*binned_distances[0]), # do the weighting
                marker = dict(color = colors[i], opacity = 0.66),
                showlegend = False,
                row = 2, col = 1,
                )


        fig.add_scatter(# plot psi-squared
            x = binned_distances[0],
            y = np.max(binned_distances[1:])*np.exp(-2*alpha*binned_distances[0]),
            mode="lines",
            line = dict(color = "red",),
            showlegend = False,
            row = 1, col = 1
            )

        summed_ligand = np.sum(binned_distances[1]*np.exp(-2*alpha*binned_distances[0]))
        summed_solvent = np.sum(binned_distances[3]*np.exp(-2*alpha*binned_distances[0]))

        fig.add_bar(# total
            orientation = "h",
            y = [1], x = [1],
            marker = dict(color = colors[0]),
            showlegend = False,
            row = 3, col = 1
            )
        fig.add_bar(# ligand fraction
            orientation = "h",
            y = [1], x = [summed_ligand/(summed_ligand+summed_solvent)],
            marker = dict(color = colors[-1]),
            showlegend = False,
            row = 3, col = 1
            )

        fig.add_annotation(
            text = f"solvent:ligand::{summed_solvent/summed_ligand:.2f}:1",
            x = 0, 
            y = 1,
            xanchor = "left",
            yanchor = "middle",
            showarrow = False,
            row = 3, col = 1,
            )
        fig.add_annotation(
            text = f"ligand:solvent::{summed_ligand/summed_solvent:.2f}:1",
            x = 1, 
            y = 1,
            xanchor = "right",
            yanchor = "middle",
            showarrow = False,
            row = 3, col = 1,
            )
    
        fig.update_xaxes(
            title = "distance from surface / A",
            row = 1, col = 1,
            )
        fig.update_yaxes(
            title = "counts",
            row = 1, col = 1,
            )

        fig.update_xaxes(
            title = "distance from surface / A",
            row = 2, col = 1,
            )
        fig.update_yaxes(
            title = "relative contribution to signal",
            row = 2, col = 1,
            )

        #update axes for the divided bar chart
        fig.update_xaxes(
            title = "fractional contribution to signal",
            range = [-0.01,1.01],
            row = 3, col = 1,
            showgrid = True,
            gridcolor = "black",
            showline = False,
            ticks = "",
            nticks = 11,
            )
        fig.update_yaxes(
            range = [0.25,1.75],
            ticks = "",
            showticklabels = False,
            showline = False,
            row = 3, col = 1,
            )



        fig.update_layout(
            template = "simple_white",
            barmode = "overlay",
            bargap = 0,
            margin = dict(t=5, r = 5),
            width = 600, height=500,
            )

        return fig
    return (generate_binned_distance_plot,)


@app.cell
def distance_pdf(
    alpha_input,
    extended_binned_data,
    generate_binned_distance_plot,
):
    generate_binned_distance_plot(extended_binned_data, alpha = alpha_input.value)
    return


@app.cell
def import_spectra(np):
    # import spectra to fit
    experimental_spectra = {}
    experimental_spectra["H25_H8"] = np.genfromtxt(r".\MIMS\PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)
    experimental_spectra["D25_H8"] = np.genfromtxt(r".\MIMS\dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter=",", unpack = True, skip_header = 1)
    experimental_spectra["H25_D8"] = np.genfromtxt(r".\MIMS\Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter=",", unpack = True, skip_header = 1)
    return (experimental_spectra,)


@app.cell
def _(LKP, experimental_spectra, make_subplots):
    def plot_specra_and_simulation_from_binned_distances(
        exp_spectra, 
        A,
        alpha,
        mu,
        gamma,
    
        I1,
        I2,
        I3, 
    
        a1, 
        a2,
        a3,
    
        binned_distances,
        ligand_H_fractions,
        solvent_H_fractions,
        ):

        labels = [
            "SC<sub>12</sub>H<sub>25</sub> + C<sub>7</sub>H<sub>8</sub>",
            "SC<sub>12</sub><b>D<sub>25</sub></b> + C<sub>7</sub>H<sub>8</sub>",
            "SC<sub>12</sub>H<sub>25</sub> + C<sub>7</sub><b>D<sub>8</sub></b>",
        ]

        experimental_keys = [
            "H25_H8",
            "D25_H8",
            "H25_D8"
        ]

        Is = [I1, I2, I3]
        baselines = [a1, a2, a3]

    
        # get stick spectra...
    
        fig = make_subplots(rows =3, cols = 1)
    
        for i in range(3):
            # plot experimental spectra
            fig.add_scatter(
                x = exp_spectra[experimental_keys[i]][0],
                y = exp_spectra[experimental_keys[i]][1],
                mode = "lines",
                line= dict(color = "grey", width = 10),
                showlegend = False,
                row=i+1, col=1
            )
        
            # add simulated spectra
            simulated_ys = LKP.simulate_pulse_signal(
                experimental_spectra[experimental_keys[i]][0], #experimental xs
        
                binned_distances[0], # bin centers for the ligand
                binned_distances[0], # bin centers for the solvent
        
                binned_distances[1], # binned ligand counts
                binned_distances[3], # binned_solvent counts, 
        
                ligand_H_fractions[i], # to be supplied as not fit
                solvent_H_fractions[i], # to be supplied as not fit
        
                alpha, # decay constant
                A, # scalar value to adjust strength of coupling
        
                Is[i], # intensity of individual contributions
                mu, # central position, without coupling
                gamma, # width of the lorentzian peak
        
                baselines[i], # baseline correction
                )

            fig.add_scatter(
                x = exp_spectra[experimental_keys[i]][0],
                y = simulated_ys,
                mode = "lines",
                line= dict(color = "white", width = 4),
                showlegend = False,
                row=i+1, col=1
            )
            fig.add_scatter(
                x = exp_spectra[experimental_keys[i]][0],
                y = simulated_ys,
                mode = "lines",
                line= dict(color = "red", width = 2),
                showlegend = False,
                row=i+1, col=1
            )

            # add in stick diagram
            '''
            #first generatethe lists for this. 
            stick_xs = []
            stick_ys = [] 
            for x, y in zip(

                ):
                stick_xs.append(x)
                stick_ys.append(0)

                stick_xs.append(x)
                stick_ys.append(y)

                stick_xs.append(None) # create gap in line
                stick_ys.append(None)
            fig.add_scatter(
                x = stick_xs,
                y = stick_ys,
                mode = "lines",
                line = dict(color = "black", width = 2),
                showlegend = False,
                row = i, col = 1,
            )
            '''

        fig.update_layout(
            template="simple_white",
            margin = dict(t = 5, r = 5),
            width = 600, height=500,
            )
    
        return fig

    return (plot_specra_and_simulation_from_binned_distances,)


@app.cell
def _(experimental_spectra, mo, np):
    bin_width_input = mo.ui.number(start = 0.01, step = 0.01, value = 0.1, label = "bin width")

    contribution_threshold_input = mo.ui.number(start = 0.0001, step = 0.0001, value = 0.01, label = "contribution threshold")

    alpha_input = mo.ui.number(start = 1e-12, step = 0.00001, value = 0.1, label="alpha")
    A_input = mo.ui.number(start = 0, step = 0.1, label = "A")

    mu_input = mo.ui.number(value = 50, label = "mu")
    gamma_input = mo.ui.number(start = 0, step =0.1, value = 0.1, label="gamma")



    a_H25_H8_input = mo.ui.number(step = 1, value = np.min(experimental_spectra["H25_H8"][1]), label="a_H25_H8")
    a_D25_H8_input = mo.ui.number(step = 1, value = np.min(experimental_spectra["D25_H8"][1]), label="a_D25_H8")
    a_H25_D8_input = mo.ui.number(step = 1, 
                                  value = np.min(experimental_spectra["H25_D8"][1]), 
                                  label="a_H25_D8")
    return (
        A_input,
        a_D25_H8_input,
        a_H25_D8_input,
        a_H25_H8_input,
        alpha_input,
        bin_width_input,
        contribution_threshold_input,
        gamma_input,
        mu_input,
    )


@app.cell
def _(
    alpha_input,
    experimental_spectra,
    extended_binned_data,
    gamma_input,
    mo,
    np,
):
    # need to set up the I-inputs separately, to avoid a cyclic reference


    I_H25_H8_input = mo.ui.number(start = 0, step = 1, 
                                  value = (np.max(experimental_spectra["H25_H8"][1]) - np.min(experimental_spectra["H25_H8"][1])) / (np.sum(extended_binned_data[3]*np.exp(-2*alpha_input.value * extended_binned_data[0])) + np.sum(extended_binned_data[1]*np.exp(-2*alpha_input.value * extended_binned_data[0]))) *gamma_input.value, 
                                  label="I_H25_H8")

    I_D25_H8_input = mo.ui.number(start = 0, step = 1, value = (np.max(experimental_spectra["D25_H8"][1]) - np.min(experimental_spectra["D25_H8"][1])) / (np.sum(extended_binned_data[3]*np.exp(-2*alpha_input.value * extended_binned_data[0]))) *gamma_input.value, 
                                  label="I_D25_H8")

    I_H25_D8_input = mo.ui.number(start = 0, step = 1, value = (np.max(experimental_spectra["H25_D8"][1]) - np.min(experimental_spectra["H25_D8"][1])) / (np.sum(extended_binned_data[1]*np.exp(-2*alpha_input.value * extended_binned_data[0]))) *gamma_input.value, 
                                  label="I_H25_D8")
    return I_D25_H8_input, I_H25_D8_input, I_H25_H8_input


@app.cell
def _(
    A_input,
    I_D25_H8_input,
    I_H25_D8_input,
    I_H25_H8_input,
    a_D25_H8_input,
    a_H25_D8_input,
    a_H25_H8_input,
    alpha_input,
    bin_width_input,
    contribution_threshold_input,
    gamma_input,
    mo,
    mu_input,
):
    mo.md(
        f"""
    Binning parameters:</br>
    {bin_width_input}</br>
    {contribution_threshold_input}</br>

    </br>
    Spectral parameters:</br>
    {alpha_input}</br>
    {A_input}</br>
    {mu_input}</br>
    {gamma_input}</br>


    {I_H25_H8_input}
    {I_D25_H8_input}
    {I_H25_D8_input}

    {a_H25_H8_input}
    {a_D25_H8_input}
    {a_H25_D8_input}
    """
    )
    return


@app.cell
def _():
    return


@app.cell
def _(
    A_input,
    I_D25_H8_input,
    I_H25_D8_input,
    I_H25_H8_input,
    a_D25_H8_input,
    a_H25_D8_input,
    a_H25_H8_input,
    alpha_input,
    experimental_spectra,
    extended_binned_data,
    gamma_input,
    ligand_H_fractions,
    mu_input,
    plot_specra_and_simulation_from_binned_distances,
    solvent_H_fractions,
):
    plot_specra_and_simulation_from_binned_distances(
        experimental_spectra, 
    
        A_input.value,
        alpha_input.value,
    
        mu_input.value,
        gamma_input.value, 
    
        I_H25_H8_input.value,
        I_D25_H8_input.value,
        I_H25_D8_input.value, 
    
        a_H25_H8_input.value, 
        a_D25_H8_input.value,
        a_H25_D8_input.value,
    
        extended_binned_data,
        ligand_H_fractions,
        solvent_H_fractions,
        )
    return


if __name__ == "__main__":
    app.run()
