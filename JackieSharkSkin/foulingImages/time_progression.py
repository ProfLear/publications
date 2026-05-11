import marimo

__generated_with = "0.21.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    from PIL import Image
    from plotly.subplots import make_subplots 
    import numpy as np
    from pathlib import Path
    from scipy.integrate import odeint
    from lmfit import Parameters, Minimizer, report_fit
    import codechembook.quickPlots as qp
    import plotly.graph_objects as go
    from scipy.stats import norm
    import scipy.stats as stats
    import pandas as pd
    import numpy as np
    import scipy.stats as stats
    import pandas as pd
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    #from statsmodels.stats.multicomp import pairwise_tukeyhsd
    return (
        Image,
        Minimizer,
        Parameters,
        Path,
        go,
        make_subplots,
        norm,
        np,
        odeint,
        pairwise_tukeyhsd,
        pd,
        report_fit,
        stats,
    )


@app.cell
def _(get_class_coverage, np, patterns, samples):
    # washed analysis
    washed_clear_fracs = []
    for _pattern in patterns:
        _temp_frac_clear = []

        for _sample in samples:
            _image_path = f"./retraining_fouling/{_pattern}/washed/{_pattern}_washed{_sample}_segmented.jpg"
            try:
                _temp_frac_clear.append(get_class_coverage(_image_path)[0])
            except:
                print(f"./retraining_fouling/{_pattern}/washed/{_pattern}_washed{_sample}_segmented.jpg")
        _temp_frac_clear = np.array(_temp_frac_clear)
        washed_clear_fracs.append(_temp_frac_clear)
        print(f"{_pattern} washed percent amounts clear are: {_temp_frac_clear*100}. mean = {np.mean(_temp_frac_clear)*100:.1f}%, stderr = {np.std(_temp_frac_clear)/(len(_temp_frac_clear)**0.5)*100:.1f}% clear")
    return (washed_clear_fracs,)


@app.cell
def _(np, pairwise_tukeyhsd, stats, washed_clear_fracs):


    group1, group2, group3 = washed_clear_fracs
    # Perform One-way ANOVA
    washed_f_stat, washed_p_val_anova = stats.f_oneway(group1, group2, group3)

    # Prepare data for Tukey HSD
    data = list(group1) + list(group2) + list(group3)
    labels = (['Group 1'] * len(group1)) + (['Group 2'] * len(group2)) + (['Group 3'] * len(group3))
    washed_tukey = pairwise_tukeyhsd(endog=data, groups=labels, alpha=0.05)

    # Calculate descriptive stats
    means = [np.mean(group1), np.mean(group2), np.mean(group3)]
    stds = [np.std(group1, ddof=1), np.std(group2, ddof=1), np.std(group3, ddof=1)]

    print(f"ANOVA F-statistic: {washed_f_stat}")
    print(f"ANOVA P-value: {washed_p_val_anova}")
    print(f"Means: {means}")
    print(f"Standard Deviations: {stds}")
    print("\nTukey HSD Results:")
    print(washed_tukey)
    return


@app.cell
def _():

    colors = [(7, 115, 177), (231,159, 39), (162, 55, 104)]
    return (colors,)


@app.cell
def _(Path, get_class_coverage):
    patterns = ["flat", "small", "large"]
    samples = ["a", "b", "c", "d"]
    all_coverage = {} # Renamed to avoid collision
    root = "retraining_fouling/" # blank for mine or /jackie_weka

    for _pattern in patterns:
        all_coverage[_pattern] = {}
        for _week in range(1, 13):
            folder = Path(f"./{root}/{_pattern}/week{_week}")
            for _sample in samples:
                image_path = folder / f"{_pattern}_{_week}{_sample}_segmented.jpg"

                try:
                    # Use a unique name for the function output
                    temp_coverage = get_class_coverage(image_path) 
                    if _week not in all_coverage[_pattern]:
                        all_coverage[_pattern][_week] = {}
                    all_coverage[_pattern][_week][_sample] = temp_coverage
                    print(f"{image_path} read!")
                except Exception as e:
                    print(f"{image_path} not found! Error: {e}")

    print(all_coverage)
    return all_coverage, patterns, samples


@app.cell
def _(Image, np):
    def get_class_coverage(image_path):
        img = Image.open(image_path).convert("RGB")
        pixels = img.getdata()

        total_pixels = len(pixels)
        clear_pixels = 0
        thick_pixels = 0
        for pixel in pixels:
            if pixel[0] > 200 and pixel[1] < 10 and pixel[2] < 10:
                clear_pixels +=1
            elif pixel[0] > 100 and pixel[0] < 200 and pixel[1] > 100 and pixel[1] < 200 and pixel[2] > 200: # this is purple
                thick_pixels += 1


        thin_pixels = total_pixels - clear_pixels - thick_pixels

        return np.array([clear_pixels, thin_pixels, thick_pixels])/total_pixels

    return (get_class_coverage,)


@app.cell
def _(all_coverage, colors, make_subplots, patterns):
    timePlot = make_subplots(cols = 1, rows = 3, subplot_titles=["bare", "thin coverage", "thick coverage"])


    for _pattern, _color in zip(patterns, colors):
        xs = [0]
        free_ys_a = [1]
        free_ys_b = [1] 
        free_ys_c = [1]
        free_ys_d = [1]

        thin_ys_a = [0]
        thin_ys_b = [0] 
        thin_ys_c = [0]
        thin_ys_d = [0]

        thick_ys_a = [0]
        thick_ys_b = [0] 
        thick_ys_c = [0]
        thick_ys_d = [0]

        thin_ys = []
        thick_ys = []
        for _week in all_coverage[_pattern]:
            print(_pattern, _week)
            xs.append(int(_week))
            try:
                free_ys_a.append(all_coverage[_pattern][_week]["a"][0])
                thin_ys_a.append(all_coverage[_pattern][_week]["a"][1])
                thick_ys_a.append(all_coverage[_pattern][_week]["a"][2])
            except: 
                pass

            try:
                free_ys_b.append(all_coverage[_pattern][_week]["b"][0])
                thin_ys_b.append(all_coverage[_pattern][_week]["b"][1])
                thick_ys_b.append(all_coverage[_pattern][_week]["b"][2])
            except: 
                pass

            try:
                free_ys_c.append(all_coverage[_pattern][_week]["c"][0])
                thin_ys_c.append(all_coverage[_pattern][_week]["c"][1])
                thick_ys_c.append(all_coverage[_pattern][_week]["c"][2])
            except:
                pass

            try:
                free_ys_d.append(all_coverage[_pattern][_week]["d"][0])
                thin_ys_d.append(all_coverage[_pattern][_week]["d"][1])
                thick_ys_d.append(all_coverage[_pattern][_week]["d"][2])
            except:
                pass

        print(xs)
        timePlot.add_scatter(x = xs, y = free_ys_a, col =1 , row = 1, mode = "lines+markers", marker = dict(color = _color), name = f"{_pattern} a")
        timePlot.add_scatter(x = xs, y = free_ys_b, col =1 , row = 1, mode = "lines+markers", 
                             marker = dict(color = _color), 
                             line = dict(dash = "longdash"),
                             name = f"{_pattern} b")
        timePlot.add_scatter(x = xs, y = free_ys_c, col =1 , row = 1, mode = "lines+markers", marker = dict(color = _color), 
                             line = dict(dash = "dash"),
                             name = f"{_pattern} c")
        timePlot.add_scatter(x = xs, y = free_ys_d, col =1 , row = 1, mode = "lines+markers", marker = dict(color = _color), 
                             line = dict(dash = "dot"),
                             name = f"{_pattern} d")



        timePlot.add_scatter(x = xs, y = thin_ys_a, col =1 , row = 2, mode = "lines+markers", marker = dict(color = _color), name = f"{_pattern} a", showlegend = False)
        timePlot.add_scatter(x = xs, y = thin_ys_b, col =1 , row = 2, mode = "lines+markers", 
                             marker = dict(color = _color), 
                             line = dict(dash = "longdash"),
                             name = f"{_pattern} b", showlegend = False)
        timePlot.add_scatter(x = xs, y = thin_ys_c, col =1 , row = 2, mode = "lines+markers", marker = dict(color = _color), 
                             line = dict(dash = "dash"),
                             name = f"{_pattern} c", showlegend = False)
        timePlot.add_scatter(x = xs, y = thin_ys_d, col =1 , row = 2, mode = "lines+markers", marker = dict(color = _color), 
                             line = dict(dash = "dot"),
                             name = f"{_pattern} d", showlegend = False)



        timePlot.add_scatter(x = xs, y = thick_ys_a, col =1 , row = 3, mode = "lines+markers", marker = dict(color = _color), name = f"{_pattern} a", showlegend = False)
        timePlot.add_scatter(x = xs, y = thick_ys_b, col =1 , row = 3, mode = "lines+markers", 
                             marker = dict(color = _color), 
                             line = dict(dash = "longdash"),
                             name = f"{_pattern} b", showlegend = False)
        timePlot.add_scatter(x = xs, y = thick_ys_c, col =1 , row = 3, mode = "lines+markers", marker = dict(color = _color), 
                             line = dict(dash = "dash"),
                             name = f"{_pattern} c", showlegend = False)
        timePlot.add_scatter(x = xs, y = thick_ys_d, col =1 , row = 3, mode = "lines+markers", marker = dict(color = _color), 
                             line = dict(dash = "dot"),
                             name = f"{_pattern} d", showlegend = False)


    timePlot.update_yaxes(range = [0,1])
    timePlot.update_layout(height = 1800, width = 900)
    timePlot.show()
    return (timePlot,)


@app.cell
def _(colors, make_subplots, np, timePlot):
    ave_timePlot = make_subplots()

    for _i, _sample, _color, _dash in zip(range(3), ["flat", "small", "large"], colors, ["solid", "solid", "solid"]):
        _temp_mean = (np.array(timePlot.data[0+12*_i]['y'])+ 
                np.array(timePlot.data[1+12*_i]['y'])+
                np.array(timePlot.data[2+12*_i]['y'])+
                np.array(timePlot.data[3+12*_i]['y'])
                )/4
        ave_timePlot.add_scatter(x = timePlot.data[0+4*_i]['x'], y = 1 - _temp_mean,
                                error_y = dict(array = (
                                    (np.array(timePlot.data[0+12*_i]['y'])-_temp_mean)**2+ 
                                    (np.array(timePlot.data[1+12*_i]['y'])-_temp_mean)**2+
                                    (np.array(timePlot.data[2+12*_i]['y'])-_temp_mean)**2+
                                    (np.array(timePlot.data[3+12*_i]['y'])-_temp_mean)**2

                                )**0.5
                                              ), 
                                name = _sample, line = dict(color = f"rgb{_color}", dash = _dash))


    ave_timePlot.update_xaxes(title = "weeks", dtick = 4)
    ave_timePlot.update_yaxes(title = "fraction covered")
    ave_timePlot.update_layout(template = "simple_white", paper_bgcolor = "rgba(0,0,0,0)", width = 600, height = 400)
    ave_timePlot.show()
    ave_timePlot.write_image("/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/JackieSharkSkin/fouling/processingReults/timeseries.svg")
    return


@app.cell
def _(
    Minimizer,
    Parameters,
    all_coverage,
    colors,
    objective,
    pairwise_tukeyhsd,
    patterns,
    pd,
    plot_biofilm_plotly,
    report_fit,
    samples,
    stats,
):
    # this will fit each trace individually.abs
    ind_fit_result_dict = {
        'Pattern': ['flat']*4 + ['small']*4 + ['large']*4,
        'k1_rate': [],
        'k2_rate': [] 
    }
    for _pattern, _color in zip(patterns, colors):
        for _sample in samples: # go through a, b, c, d first
            ind_weeks_to_fit = [0]
            ind_total_coverage_to_fit = [0]
            ind_sparse_coverage_to_fit = [0]
            ind_dense_coverage_to_fit = [0]
            for _week in all_coverage[_pattern]:
                try:
                    ind_sparse_coverage_to_fit.append(all_coverage[_pattern][_week][_sample][1])
                    ind_dense_coverage_to_fit.append(all_coverage[_pattern][_week][_sample][2])
                    ind_weeks_to_fit.append(_week)

                except:
                    pass

            # --- Fit ---
            ind_params = Parameters()
            ind_params.add('k1', value=0.5, min=0)  # Growth rate (per week)
            ind_params.add('k2', value=0.1, min=0)  # Transition rate (per week)
            ind_params.add('B0', value=0.01, vary=True) # Initial seeding %
            ind_params.add('C0', value=0.0, vary=False)      # Assume no dense film at start

            ind_minner = Minimizer(objective, ind_params, fcn_args=(ind_weeks_to_fit, ind_sparse_coverage_to_fit, ind_dense_coverage_to_fit))
            ind_result = ind_minner.minimize()
            report_fit(ind_result)

            ind_fit_result_dict["k1_rate"].append(ind_result.params["k1"].value)
            ind_fit_result_dict["k2_rate"].append(ind_result.params["k2"].value)


            ind_fig = plot_biofilm_plotly(ind_result, ind_weeks_to_fit, ind_sparse_coverage_to_fit, ind_dense_coverage_to_fit, _pattern)
            ind_fig.show()

    """
    # 1. Organize your fitted rates into a DataFrame
    # Assume you've collected your k1 values from your 12 fits (3 patterns * 4 replicates)
    data = {
        'Pattern': ['SharkSkin']*4 + ['Smooth']*4 + ['Control']*4,
        'k1_rate': [0.85, 0.82, 0.88, 0.84,  # SharkSkin replicates
                    0.60, 0.62, 0.58, 0.61,  # Smooth replicates
                    0.95, 0.92, 0.98, 0.94]   # Control replicates
    }
    """
    df = pd.DataFrame(ind_fit_result_dict)

    print(df)

    # 2. Perform One-Way ANOVA
    for constant in ["k1_rate", "k2_rate"]:
        f_stat, p_val = stats.f_oneway(
            df[df['Pattern'] == 'flat'][constant],
            df[df['Pattern'] == 'small'][constant],
            df[df['Pattern'] == 'large'][constant]
        )

        print(f"{constant} ANOVA p-value: {p_val:.5f}")

        # 3. Perform Tukey HSD (if ANOVA p-value < 0.05)
        if p_val < 0.05:
            tukey = pairwise_tukeyhsd(endog=df[constant], 
                                      groups=df['Pattern'], 
                                      alpha=0.05)
            print(tukey)
    return


@app.cell
def _(
    Minimizer,
    Parameters,
    all_coverage,
    colors,
    np,
    objective,
    pairwise_z_test,
    patterns,
    plot_biofilm_plotly,
    report_fit,
    samples,
    total_growth_objective,
):
    # this approach fits a line to all the points at once, per pattern. 

    fit_result_dict = {}
    total_fit_result_dict = {}
    for _pattern, _color in zip(patterns, colors):
        weeks_to_fit = [0,0,0,0]
        total_coverage_to_fit = [0, 0, 0, 0]
        sparse_coverage_to_fit = [0, 0, 0, 0]
        dense_coverage_to_fit = [0, 0, 0, 0]

        for _week in all_coverage[_pattern]:
            for _sample in samples: # go through a, b, c, d
                try:
                    sparse_coverage_to_fit.append(all_coverage[_pattern][_week][_sample][1])
                    dense_coverage_to_fit.append(all_coverage[_pattern][_week][_sample][2])
                    weeks_to_fit.append(_week)
                except:
                    pass

        total_coverage_to_fit = np.array(sparse_coverage_to_fit)+np.array(dense_coverage_to_fit) # get the total coverage. 

        # --- Fit ---
        params = Parameters()
        params.add('k1', value=0.5, min=0)  # Growth rate (per week)
        params.add('k2', value=0.1, min=0)  # Transition rate (per week)
        params.add('B0', value=0.01, vary=True) # Initial seeding %
        params.add('C0', value=0.0, vary=False)      # Assume no dense film at start

        minner = Minimizer(objective, params, fcn_args=(weeks_to_fit, sparse_coverage_to_fit, dense_coverage_to_fit))
        result = minner.minimize()
        report_fit(result)

        fit_result_dict[_pattern] = result


        fig = plot_biofilm_plotly(result, weeks_to_fit, sparse_coverage_to_fit, dense_coverage_to_fit, _pattern)
        fig.show()


        total_params = Parameters()
        total_params.add('k_tot', value = 0.5, min = 0)
        total_params.add("T0", value = 0.01, min = 0)

        total_minner = Minimizer(total_growth_objective, total_params, fcn_args=(np.array(weeks_to_fit), total_coverage_to_fit))
        total_result = total_minner.minimize()
        report_fit(total_result)

        total_fit_result_dict[_pattern] = total_result

        total_fig = plot_biofilm_plotly(total_result, weeks_to_fit, total_coverage_to_fit,_pattern)
        total_fig.show()


    print(pairwise_z_test(fit_result_dict, param_name="k1"))
    print(pairwise_z_test(fit_result_dict, param_name="k2"))

    print(pairwise_z_test(total_fit_result_dict, param_name="k_tot"))
    return


@app.cell
def _(go, make_subplots, np, solve_model):

    def plot_biofilm_plotly(result, t_full, data_1, data_2=None, pattern_name=""):
        """
        Handles both Sequential (sparse/dense) and Total (logistic) fits.
        If data_2 is None, it assumes a Total Coverage fit.
        """
        t_full = np.asanyarray(t_full)
        data_1 = np.asanyarray(data_1)

        # Create subplots: Residuals on top, Fit on bottom
        fig = make_subplots(
            rows=2, cols=1, 
            shared_xaxes=True, 
            vertical_spacing=0.07,
            row_heights=[0.3, 0.7],
            subplot_titles=("Residuals", f"Model Fit: {pattern_name}")
        )

        t_plot = np.linspace(t_full.min(), t_full.max(), 150)

        # --- CASE 1: TOTAL COVERAGE MODEL (k_tot) ---
        if 'k_tot' in result.params:
            k_tot = result.params['k_tot'].value
            T0 = result.params['T0'].value
            K = 1.0 # Carrying capacity

            # Analytical Logistic Solution
            fit_curve = K / (1 + ((K - T0) / (T0 + 1e-9)) * np.exp(-k_tot * t_plot))
            fit_at_data = K / (1 + ((K - T0) / (T0 + 1e-9)) * np.exp(-k_tot * t_full))
            residuals = data_1 - fit_at_data

            # Main Plot
            fig.add_trace(go.Scatter(x=t_full, y=data_1, mode='markers', name='Total Exp', 
                                     marker=dict(color='blue', opacity=0.5)), row=2, col=1)
            fig.add_trace(go.Scatter(x=t_plot, y=fit_curve, mode='lines', name='Total Fit', 
                                     line=dict(color='blue', width=3)), row=2, col=1)
            # Residuals
            fig.add_trace(go.Scatter(x=t_full, y=residuals, mode='markers', name='Res.', 
                                     marker=dict(color='blue', symbol='x')), row=1, col=1)

        # --- CASE 2: SEQUENTIAL MODEL (k1, k2) ---
        else:
            k1, k2 = result.params['k1'].value, result.params['k2'].value
            B0, C0 = result.params['B0'].value, result.params['C0'].value
            data_2 = np.asanyarray(data_2)

            # Generate ODE fit lines
            fit_curve = solve_model(t_plot, k1, k2, B0, C0)
            fit_at_data = solve_model(t_full, k1, k2, B0, C0)

            res_sparse = data_1 - fit_at_data[:, 0]
            res_dense = data_2 - fit_at_data[:, 1]

            # Main Plot: Sparse
            fig.add_trace(go.Scatter(x=t_full, y=data_1, mode='markers', name='Sparse Exp', 
                                     marker=dict(color='orange', opacity=0.5)), row=2, col=1)
            fig.add_trace(go.Scatter(x=t_plot, y=fit_curve[:, 0], mode='lines', name='Sparse Fit', 
                                     line=dict(color='orange')), row=2, col=1)
            # Main Plot: Dense
            fig.add_trace(go.Scatter(x=t_full, y=data_2, mode='markers', name='Dense Exp', 
                                     marker=dict(color='green', opacity=0.5)), row=2, col=1)
            fig.add_trace(go.Scatter(x=t_plot, y=fit_curve[:, 1], mode='lines', name='Dense Fit', 
                                     line=dict(color='green')), row=2, col=1)

            # Residuals
            fig.add_trace(go.Scatter(x=t_full, y=res_sparse, mode='markers', name='Sparse Res.', 
                                     marker=dict(color='orange', symbol='circle-open')), row=1, col=1)
            fig.add_trace(go.Scatter(x=t_full, y=res_dense, mode='markers', name='Dense Res.', 
                                     marker=dict(color='green', symbol='circle-open')), row=1, col=1)

        # Common Layout
        fig.add_hline(y=0, line_dash="dash", line_color="black", row=1, col=1)
        fig.update_layout(height=700, template="plotly_white", showlegend=True)
        fig.update_yaxes(title_text="Residuals", row=1, col=1)
        fig.update_yaxes(title_text="% Coverage", row=2, col=1)
        return fig

    return (plot_biofilm_plotly,)


@app.cell
def _(norm, np):
    def pairwise_z_test(results_dict, param_name='k1'):
        """
        results_dict: {'Pattern_A': result_obj, 'Pattern_B': result_obj, ...}
        param_name: The string name of the parameter in lmfit ('k1' or 'k2')
        """
        keys = list(results_dict.keys())
        comparisons = []

        # 3 pairwise comparisons for 3 patterns
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                p1, p2 = keys[i], keys[j]

                # Extract k and SE
                val1 = results_dict[p1].params[param_name].value
                se1 = results_dict[p1].params[param_name].stderr

                val2 = results_dict[p2].params[param_name].value
                se2 = results_dict[p2].params[param_name].stderr

                # Calculate Z and p (two-tailed)
                z_score = np.abs(val1 - val2) / np.sqrt(se1**2 + se2**2)
                p_value = 2 * (1 - norm.cdf(z_score))

                comparisons.append({
                    'Comparison': f"{p1} vs {p2}",
                    'Diff': val1 - val2,
                    'Z-score': z_score,
                    'p-value': p_value,
                    'Significant (adj)': p_value < (0.05 / 3)
                })

        return comparisons

    return (pairwise_z_test,)


@app.cell
def _(np, odeint):
    # Define the model: A (Empty) -> B (Sparse) -> C (Dense)
    def biofilm_ode(y, t, k1, k2):
        # A = 100 - (B + C)
        B, C = y
        A = 1 - (B + C)

        # B grows logistically using available space A
        # B converts to C at rate k2
        dBdt = k1 * B * (A / 1) - k2 * B
        dCdt = k2 * B

        return [dBdt, dCdt]

    def solve_model(t_unique, k1, k2, B0, C0):
        """Solves ODE for unique timepoints only (efficiency)"""
        sol = odeint(biofilm_ode, [B0, C0], t_unique, args=(k1, k2))
        return sol

    def objective(params, t_full, data_B, data_C):
        # Extract the FLOAT values from the Parameter objects
        k1_val = params['k1'].value
        k2_val = params['k2'].value
        B0_val = params['B0'].value
        C0_val = params['C0'].value

        # 1. Get unique timepoints
        t_unique = np.unique(t_full)

        # 2. Pass the FLOAT values to the solver, not the Parameter objects
        sol_unique = solve_model(t_unique, k1_val, k2_val, B0_val, C0_val)

        # 3. Map results back to the full time array
        mapping = {t: i for i, t in enumerate(t_unique)}
        indices = [mapping[t] for t in t_full]

        pred_B = sol_unique[indices, 0]
        pred_C = sol_unique[indices, 1]

        # Calculate residuals
        res_B = pred_B - data_B
        res_C = pred_C - data_C

        return np.concatenate([res_B, res_C])



    def total_growth_objective(params, t_full, total_data):
        '''
        This function models just the total coverage.  So adding together thin (sparse) and thick (dense)
        '''
        k_tot = params['k_tot'].value
        T0 = params['T0'].value

        # Analytical solution to the Logistic Equation:
        # T(t) = K / (1 + ((K - T0)/T0) * exp(-k*t))
        K = 1
        prediction = K / (1 + ((K - T0)/T0) * np.exp(-k_tot * t_full))

        return prediction - total_data

    return objective, solve_model, total_growth_objective


if __name__ == "__main__":
    app.run()
