import marimo

__generated_with = "0.20.4"
app = marimo.App()


@app.cell
def _():
    # import all needed modules
    import marimo as mo
    import numpy as np
    from importlib import reload
    from scipy.ndimage import gaussian_filter1d
    from plotly.subplots import make_subplots
    from lmfit import Model
    import plotly.graph_objects as go
    import torch

    # import optimized modules
    import HighDensityDipolarGPU as hdg
    import mdHelper as mh

    reload(hdg)
    reload(mh)
    return Model, go, hdg, make_subplots, mh, np


@app.cell
def _():
    from datetime import datetime

    # Get today's date
    today = datetime.today()

    # Format it as YYYYMMDD
    formatted_date = today.strftime('%Y%m%d')
    return (formatted_date,)


@app.cell
def _(np):
    # read spectral data files
    field_1 = 11745
    field_2 = 11729
    nu_G = 0.004262

    dodecaneH_tolueneH_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter=",", unpack=True, skip_header=1)
    dodecaneD_tolueneH_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter=",", unpack=True, skip_header=1)
    dodecaneH_tolueneD_spectra = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MimsData/Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter=",", unpack=True, skip_header=1)

    dodecaneH_tolueneH_spectra[0] = (dodecaneH_tolueneH_spectra[0] - field_2 * nu_G) * 1e6
    dodecaneD_tolueneH_spectra[0] = (dodecaneD_tolueneH_spectra[0] - field_2 * nu_G) * 1e6
    dodecaneH_tolueneD_spectra[0] = (dodecaneH_tolueneD_spectra[0] - field_1 * nu_G) * 1e6

    normalize = True
    if normalize:
        for s in [dodecaneH_tolueneH_spectra, dodecaneD_tolueneH_spectra, dodecaneH_tolueneD_spectra]:
            s[1] = (s[1] - np.mean(s[1][:10])) / np.max((s[1] - np.mean(s[1][:10])))

    mask_limit = 5e6
    mask_data = True
    if mask_data:
        res = []
        for s in [dodecaneH_tolueneH_spectra, dodecaneD_tolueneH_spectra, dodecaneH_tolueneD_spectra]:
            mask = (s[0] >= -mask_limit) & (s[0] <= mask_limit)
            res.append([s[0][mask], s[1][mask]])
        dodecaneH_tolueneH_spectra, dodecaneD_tolueneH_spectra, dodecaneH_tolueneD_spectra = res
    return (
        dodecaneD_tolueneH_spectra,
        dodecaneH_tolueneD_spectra,
        dodecaneH_tolueneH_spectra,
    )


@app.cell
def _(np):
    # read md data files
    hd_toluene_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/toluene_distances_high_density.csv", delimiter=",") * 1e-10
    hd_dodecanethiol_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/dodecanol_distances_high_density.csv", delimiter=",") * 1e-10
    ld_toluene_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/toluene_distances_low_density.csv", delimiter=",") * 1e-10
    ld_dodecanethiol_distances = np.genfromtxt("/Users/benjaminlear/GitHub/publications/KristenPulse/MDSimulationResults/dodecanol_distances_low_density.csv", delimiter=",") * 1e-10
    return (
        hd_dodecanethiol_distances,
        hd_toluene_distances,
        ld_dodecanethiol_distances,
        ld_toluene_distances,
    )


@app.cell
def _(np):
    r_Hs = np.genfromtxt("./TEMdata/Feret_All_PdSC12_fit_report.csv", unpack=True, delimiter=",", usecols=[0], skip_header=1)
    r_Ds = np.genfromtxt("./TEMdata/Feret_All_dPdSC12_lognormal.csv", unpack=True, delimiter=",", usecols=[0], skip_header=1)
    return r_Ds, r_Hs


@app.cell
def _():
    R_np_md = 15e-10
    bulk_toluene_density = 0.0453e30
    md_r_step = 0.1e-10
    return R_np_md, bulk_toluene_density, md_r_step


@app.cell
def _(
    R_np_md,
    bulk_toluene_density,
    hd_dodecanethiol_distances,
    hd_toluene_distances,
    ld_dodecanethiol_distances,
    ld_toluene_distances,
    md_r_step,
    mh,
):
    hd_toluene_profile = mh.NPDensityProfile(hd_toluene_distances, R_np_md, bin_width=md_r_step)
    hd_toluene_profile.prepare_blended_density(rho_bulk=bulk_toluene_density)
    hd_dodecanethiol_profile = mh.NPDensityProfile(hd_dodecanethiol_distances, R_np_md, old_bins=hd_toluene_profile.bin_centers)

    ld_toluene_profile = mh.NPDensityProfile(ld_toluene_distances, R_np_md, bin_width=md_r_step)
    ld_toluene_profile.prepare_blended_density(rho_bulk=bulk_toluene_density)
    ld_dodecanethiol_profile = mh.NPDensityProfile(ld_dodecanethiol_distances, R_np_md, old_bins=ld_toluene_profile.bin_centers)
    return ld_dodecanethiol_profile, ld_toluene_profile


@app.cell
def _(hdg, np):
    ### GPU ACCELERATED FUNCTION FOR FITTING

    def triple_mims_simulation_gpu(
        xs, len1, len2, I1, I2, I3, c1, c2, c3,
        H_np, D_np, solvent_density_obj, rho_bulk_solvent,
        ligand_density_obj, rho_bulk_ligand, nu_n, g_n,
        I_scale, ell, g_e, tau, T_m, R_np_md, R_max,
        target_dr=0.2E-10, dipolar_scale=1, contact_scale=1,
        blindspot_offset=0, dipolar_np_scale=1, shape="Lorentzian",
        eta=0, width_scale=1
    ):
        freq_both = xs[:len1]
        freq_solvent_only = xs[len1:len1+len2]
        freq_ligand_only = xs[len1+len2:]

        all_raw_samples = []

        # Configs: (freq, s_weight, l_weight, size_dist)
        configs = [
            (freq_both, 1.0, 1.0, H_np),
            (freq_solvent_only, 1.0, 0.02, D_np),
            (freq_ligand_only, 0.01, 1.0, H_np)
        ]

        for freq, s_w, l_w, rs in configs:
            y_s = np.zeros_like(freq)
            y_l = np.zeros_like(freq)

            if s_w > 0:
                y_s = hdg.simulate_mims_gpu(
                    freq, rs[0], rs[1], solvent_density_obj, rho_bulk_solvent,
                    nu_n, g_n, I_scale, ell, g_e, tau, T_m,
                    R_np_md, R_max, target_dr, dipolar_scale,
                    contact_scale, blindspot_offset, dipolar_np_scale,
                    shape, eta, width_scale,
                    device = "mps" # for apple M chips
                )
            if l_w > 0:
                y_l = hdg.simulate_mims_gpu(
                    freq, rs[0], rs[1], ligand_density_obj, rho_bulk_ligand,
                    nu_n, g_n, I_scale, ell, g_e, tau, T_m,
                    R_np_md, R_max, target_dr, dipolar_scale,
                    contact_scale, blindspot_offset, dipolar_np_scale,
                    shape, eta, width_scale,
                    device = "mps" # for apple M chips
                )
            all_raw_samples.append(s_w * y_s + l_w * y_l)

        scaled_both = (all_raw_samples[0] / np.max(all_raw_samples[0])) * I1 + c1
        scaled_solvent = (all_raw_samples[1] / np.max(all_raw_samples[1])) * I2 + c2
        scaled_ligand = (all_raw_samples[2] / np.max(all_raw_samples[2])) * I3 + c3

        return np.concatenate([scaled_both, scaled_solvent, scaled_ligand])

    return (triple_mims_simulation_gpu,)


@app.cell
def _(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneH_spectra,
    ld_dodecanethiol_profile,
    ld_toluene_profile,
    np,
    r_Ds,
    r_Hs,
):
    density_choice = "LD"

    # VALUES for parameters that are not fit
    len1 = len(dodecaneH_tolueneH_spectra[0])
    len2 = len(dodecaneD_tolueneH_spectra[0])

    sizes = "hist"
    if sizes == "all":
        H_np = [r_Hs*1e-9/2, np.ones_like(r_Hs)]
        D_np = [r_Ds*1e-9/2, np.ones_like(r_Ds)]
    elif sizes == "hist":
        H_np = [np.array([1.272975688, 1.405767063, 1.538558438, 1.671349814, 1.804141189, 1.936932564, 2.06972394, 2.202515315, 2.33530669, 2.468098066, 2.600889441, 2.733680816, 2.866472192, 2.999263567, 3.132054942, 3.264846318, 3.397637693, 3.530429068, 3.663220444, 3.796011819, 3.928803194, 4.06159457, 4.194385945, 4.32717732, 4.459968696, 4.592760071, 4.725551446, 4.858342822, 4.991134197, 5.123925572, 5.256716947, 5.389508323, 5.522299698, 5.655091073, 5.787882449, 5.920673824, 6.053465199, 6.186256575,])*1e-9/2, [4, 1, 2, 1, 0, 4, 1, 11, 17, 49, 82, 83, 81, 59, 55, 36, 29, 20, 25, 8, 3, 6, 10, 5, 1, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,]]
        D_np = [np.array([1.541959469, 1.728658408, 1.915357347, 2.102056286, 2.288755225, 2.475454163, 2.662153102, 2.848852041, 3.03555098, 3.222249919, 3.408948858, 3.595647796, 3.782346735, 3.969045674, 4.155744613, 4.342443552, 4.52914249, 4.715841429, 4.902540368, 5.089239307,])*1e-9/2, [14, 32, 65, 55, 61, 67, 46, 40, 26, 21, 11, 6, 8, 1, 2, 2, 0, 0, 0, 1,]]#
    else:
        H_np = [[2.87e-9/2],[1.0]]
        D_np = [[2.42e-9/2],[1.0]]

    solvent_density_obj = ld_toluene_profile
    ligand_density_obj = ld_dodecanethiol_profile
    rho_bulk_solvent = 4.528e28
    rho_bulk_ligand = 0
    g_n = 5.5856946893
    g_e = 2.0
    tau = 140e-9
    T_m = 674e-9
    target_dr = 0.1e-9
    shape = "L"

    I1_guess, I2_guess, I3_guess = 1, 1, 1
    c1_guess, c2_guess, c3_guess = 0, 0, 0
    nu_n_guess = 0
    I_scale_guess = 13
    ell_guess = 0.5e-9
    R_max_guess = 15e-9
    dipolar_scale_guess = 1
    contact_scale_guess = 1
    blindspot_offset_guess = 0.05
    dipolar_np_scale_guess = 0.685
    eta_guess = 0
    width_scale_guess = 1
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
        density_choice,
        dipolar_np_scale_guess,
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
        sizes,
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
    ligand_density_obj,
    np,
    nu_n_guess,
    rho_bulk_ligand,
    rho_bulk_solvent,
    shape,
    solvent_density_obj,
    target_dr,
    tau,
    triple_mims_simulation_gpu,
    width_scale_guess,
):
    def my_callback(params, iter, resid, *args, **kws):
        if iter % 10 == 0:
            ssr = (resid**2).sum()
            print(f"Iteration: {iter:5d} | Chi-square (SSR): {ssr:.4e}")
        return False

    triple_spectra_model = Model(triple_mims_simulation_gpu, independent_vars=[
        "xs", "len1", "len2", "H_np", "D_np", "solvent_density_obj",
        "rho_bulk_solvent", "ligand_density_obj", "rho_bulk_ligand",
        "g_n", "g_e", "tau", "T_m", "R_np_md", "target_dr", "shape",
    ])

    params = triple_spectra_model.make_params()
    params.add_many(
        ("I1", I1_guess, True, 0, None), ("I2", I2_guess, True, 0, None), ("I3", I3_guess, True, 0, None),
        ("c1", c1_guess, False), ("c2", c2_guess, False), ("c3", c3_guess, False),
        ("nu_n", nu_n_guess, False), ("I_scale", I_scale_guess, True, 0),
        ("ell", ell_guess, True, 0), ("R_max", R_max_guess, False, R_np_md),
        ('dipolar_scale', dipolar_scale_guess, False, 0),
        ('contact_scale', contact_scale_guess, False, 0),
        ('blindspot_offset', blindspot_offset_guess, True, 0, 1),
        ('dipolar_np_scale', dipolar_np_scale_guess, False, 0, 1),
        ("eta", eta_guess, False, 0), ("width_scale", width_scale_guess, True, 0),
    )

    triple_spetra_fit = triple_spectra_model.fit(
        np.concatenate([dodecaneH_tolueneH_spectra[1], dodecaneD_tolueneH_spectra[1], dodecaneH_tolueneD_spectra[1]]),
        xs=np.concatenate([dodecaneH_tolueneH_spectra[0], dodecaneD_tolueneH_spectra[0], dodecaneH_tolueneD_spectra[0]]),
        params=params,
        len1=len(dodecaneH_tolueneH_spectra[0]),
        len2=len(dodecaneD_tolueneH_spectra[0]),
        H_np=H_np, D_np=D_np,
        solvent_density_obj=solvent_density_obj, rho_bulk_solvent=rho_bulk_solvent,
        ligand_density_obj=ligand_density_obj, rho_bulk_ligand=rho_bulk_ligand,
        g_n=g_n, g_e=g_e, tau=tau, T_m=T_m, R_np_md=R_np_md,
        target_dr=target_dr, shape=shape, iter_cb=my_callback
    )
    return (triple_spetra_fit,)


@app.cell
def _(triple_spetra_fit):
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
    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[:len1], line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= dodecaneD_tolueneH_spectra[1], line = dict(color = "#ffabac", width = 5), showlegend = False, col = 1, row = 2)
    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "black"), showlegend = False, col = 1, row = 2)

    fitplot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= dodecaneH_tolueneD_spectra[1], line = dict(color = "#ababff", width = 5), showlegend = False, col = 1, row = 3)
    fitplot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= triple_spetra_fit.best_fit[len1+len2:], line = dict(color = "black"), showlegend = False, col = 1, row = 3)

    #add normalized spectra
    fitplot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= (triple_spetra_fit.best_fit[:len1] - triple_spetra_fit.params["c1"].value)/np.max((triple_spetra_fit.best_fit[:len1] - triple_spetra_fit.params["c1"].value)), line = dict(color = "black"), showlegend = False, col = 1, row = 4)
    fitplot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= (triple_spetra_fit.best_fit[len1:len1+len2] - triple_spetra_fit.params["c2"].value)/np.max(triple_spetra_fit.best_fit[len1:len1+len2] - triple_spetra_fit.params["c2"].value), line = dict(color = "black"), showlegend = False, col = 1, row = 4)
    fitplot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= (triple_spetra_fit.best_fit[len1+len2:] - triple_spetra_fit.params["c3"].value)/np.max(triple_spetra_fit.best_fit[len1+len2:] - triple_spetra_fit.params["c3"].value), line = dict(color = "black"), showlegend = False, col = 1, row = 4)


    fitplot.update_xaxes(title = "frequency /MHz", row = 4, col = 1)
    fitplot.show()
    return


@app.cell
def _(
    density_choice,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    len1,
    len2,
    make_subplots,
    np,
    sizes,
    triple_spetra_fit,
):
    # writing a fit report as a record. 
    report_plot = make_subplots(rows = 3, cols = 1, subplot_titles= ["solvent-H + ligand-H", "solvent-H + ligand-D", "solvent-D + ligand-H"])
    report_plot.update_layout(template = "simple_white", width = 500, height = 500)

    report_plot.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= dodecaneH_tolueneH_spectra[1], line = dict(color = "#ffabff", width = 5), showlegend = False, col = 1, row = 1)
    report_plot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[:len1], line = dict(color = "black"), showlegend = False, col = 1, row = 1)

    report_plot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= dodecaneD_tolueneH_spectra[1], line = dict(color = "#ffabac", width = 5), showlegend = False, col = 1, row = 2)
    report_plot.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "black"), showlegend = False, col = 1, row = 2)

    report_plot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= dodecaneH_tolueneD_spectra[1], line = dict(color = "#ababff", width = 5), showlegend = False, col = 1, row = 3)
    report_plot.add_scatter(x = dodecaneH_tolueneD_spectra[0], y= triple_spetra_fit.best_fit[len1+len2:], line = dict(color = "black"), showlegend = False, col = 1, row = 3)

    report_plot.add_annotation(text = f"{triple_spetra_fit.fit_report()}", showarrow = False, font = dict( color = "rgba(255, 255, 255,0)"), x = 0, y = np.mean(dodecaneH_tolueneH_spectra[1]))

    report_plot.update_xaxes(title = "frequency /Hz", row = 3, col = 1)
    report_plot.show()
    report_plot.write_html(f"{density_choice} {sizes} fit.html")
    return (report_plot,)


@app.cell
def _(
    density_choice,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneH_spectra,
    formatted_date,
    go,
    len1,
    len2,
    report_plot,
    sizes,
    triple_spetra_fit,
):
    pubfit = go.Figure(report_plot)
    pubcolors = ["#715472", "#94B7D4", "#BB6F35"]

    pubfit.update_traces(selector = 0, line = dict(color = "#715472", width = 5))
    pubfit.update_traces(selector = 1, line = dict(color = "white", width = 3))
    pubfit.add_scatter(x = dodecaneH_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[:len1], line = dict(color = "black", width = 1), showlegend = False, col = 1, row = 1)

    pubfit.update_traces(selector = 2, line = dict(color = "#94B7D4", width = 5))
    pubfit.update_traces(selector = 3, line = dict(color = "white", width = 3))
    pubfit.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[len1:len1+len2], line = dict(color = "black", width = 1), showlegend = False, col = 1, row = 2)

    pubfit.update_traces(selector = 4, line = dict(color = "#BB6F35", width = 5))
    pubfit.update_traces(selector = 5, line = dict(color = "white", width = 3))
    pubfit.add_scatter(x = dodecaneD_tolueneH_spectra[0], y= triple_spetra_fit.best_fit[len1+len2:], line = dict(color = "black", width = 1), showlegend = False, col = 1, row = 3)

    pubfit.layout.annotations = [] # clear annotations

    pubfit.update_yaxes(range = [-0.05,1.05], tickvals = [0, 1])
    pubfit.update_yaxes(title = "noramlized intensity", row = 2, col = 1)

    pubfit.update_layout(width = 300*2, height = 300*3)
    pubfit.show()
    pubfit.write_image(f"{formatted_date} {density_choice} {sizes} fit.svg")
    return (pubfit,)


@app.cell
def _(density_choice, formatted_date, make_subplots, pubfit, sizes):
    import copy
    pubfit2 = make_subplots(rows = 3, cols = 1, )

    pubfit2.add_trace(copy.deepcopy(pubfit.data[0]), row = 1, col = 1)
    pubfit2.add_trace(copy.deepcopy(pubfit.data[1]), row = 1, col = 1)
    pubfit2.update_traces(line = dict(width = 12), selector = 0)
    pubfit2.update_traces(line = dict(width = 4), selector = 1)
    #pubfit2.add_trace(copy.deepcopy(pubfit.data[6]), row = 1, col = 1)

    pubfit2.add_scatter(x = pubfit.data[0].x, y = pubfit.data[0].y - pubfit.data[1].y, line = dict(color = "rgba(0, 0, 0, 0.33)"), row = 1, col = 1)



    pubfit2.add_trace(copy.deepcopy(pubfit.data[2]), row = 2, col = 1)
    pubfit2.add_trace(copy.deepcopy(pubfit.data[3]), row = 2, col = 1)
    pubfit2.update_traces(line = dict(width = 12), selector = 3)
    pubfit2.update_traces(line = dict(width = 4), selector = 4)
    #pubfit2.add_trace(copy.deepcopy(pubfit.data[7]), row = 3, col = 1)

    pubfit2.add_scatter(x = pubfit.data[2].x, y = pubfit.data[2].y - pubfit.data[3].y,  line = dict(color = "rgba(0, 0, 0, 0.33)"), row = 2, col = 1)


    pubfit2.add_trace(copy.deepcopy(pubfit.data[4]), row = 3, col = 1)
    pubfit2.add_trace(copy.deepcopy(pubfit.data[5]), row = 3, col = 1)
    pubfit2.update_traces(line = dict(width = 12), selector = 6)
    pubfit2.update_traces(line = dict(width = 4), selector = 7)
    #pubfit2.add_trace(copy.deepcopy(pubfit.data[8]), row = 5, col = 1)

    pubfit2.add_scatter(x = pubfit.data[4].x, y = pubfit.data[4].y - pubfit.data[5].y,  line = dict(color = "rgba(0, 0, 0, 0.33)"), row = 3, col = 1)

    for _r in range(1,4,1):
        pubfit2.update_xaxes(showline = False, row = _r, col = 1)
        pubfit2.update_yaxes(zeroline = True)




    pubfit2.update_xaxes(range = [-4e6, 4e6])
    pubfit2.update_xaxes(title = "A /Hz", row = 3, col = 1)

    pubfit2.update_yaxes(nticks = 2, showline = False)
    pubfit2.update_yaxes(title = "normalized intensity", row = 2, col = 1)

    pubfit2.update_layout(template = "simple_white", showlegend = False, margin = dict(t = 10, b = 50, l = 50, r = 10),plot_bgcolor="#e6e6e6", width = 1.81*300, height = 1.81/1.6*3*300)
    pubfit2.show()

    pubfit2.write_html(f"{formatted_date} {density_choice} {sizes} fit.html")
    pubfit2.write_image(f"{formatted_date} {density_choice} {sizes} fit.svg")
    pubfit2.write_image(f"{formatted_date} {density_choice} {sizes} fit.png")
    return


if __name__ == "__main__":
    app.run()
