"""Comparison of EPR Pulse Inversion Recovery Lifetime Fitting:
Standard Levenberg-Marquardt (lmfit) vs. Variable Projection (VarProPy).

Authors: Benjamin Lear / Antigravity
"""

from __future__ import annotations

import os
import sys
import time
from typing import Any, Dict, List, Tuple
import numpy as np
import scipy.optimize
import lmfit
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from varpropy import ConstantModel, ExponentialModel


def load_bruker_data(dsc_path: str, dta_path: str) -> Dict[str, Any]:
    """Parse Bruker DSC descriptor and binary DTA data."""
    meta: Dict[str, Any] = {}
    with open(dsc_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                meta[parts[0]] = " ".join(parts[1:])

    xpts = int(meta.get("XPTS", 512))
    xmin = float(meta.get("XMIN", 0.0))
    xwid = float(meta.get("XWID", 2555000.0))
    xuni = meta.get("XUNI", "ns").strip("'")

    # Time axis in microseconds
    t_raw = np.linspace(xmin, xmin + xwid, xpts)
    t_us = t_raw / 1000.0 if xuni == "ns" else t_raw

    # Read binary DTA
    raw = np.fromfile(dta_path, dtype=">f8")
    re = raw[0::2]
    im = raw[1::2]
    z = re + 1j * im

    # Automated phase correction (minimizing imaginary variance)
    def phase_loss(phi: float) -> float:
        return float(np.var((z * np.exp(-1j * phi)).imag))

    opt = scipy.optimize.minimize_scalar(phase_loss, bounds=(-np.pi, np.pi), method="bounded")
    phi_opt = float(opt.x)
    z_phased = z * np.exp(-1j * phi_opt)
    # Ensure inversion recovery rises to plateau
    if z_phased.real[-1] < z_phased.real[0]:
        phi_opt += np.pi
        z_phased = z * np.exp(-1j * phi_opt)

    return {
        "t_us": t_us,
        "raw_re": re,
        "raw_im": im,
        "phased_re": z_phased.real,
        "phased_im": z_phased.imag,
        "phi_deg": np.degrees(phi_opt),
        "meta": meta,
    }


def calc_aic(chisqr: float, n: int, k: int) -> float:
    """Akaike Information Criterion (unweighted residual sum of squares)."""
    return float(n * np.log(max(chisqr, 1e-250) / n) + 2 * k)


def calc_bic(chisqr: float, n: int, k: int) -> float:
    """Bayesian Information Criterion."""
    return float(n * np.log(max(chisqr, 1e-250) / n) + k * np.log(n))


def fit_varpro(y: np.ndarray, t: np.ndarray, n_comp: int) -> Any:
    """Fit n-component exponential + baseline using VarProPy."""
    t_span = float(t[-1] - t[0])
    model = ConstantModel()
    for i in range(1, n_comp + 1):
        model = model + ExponentialModel(prefix=f"e{i}_")

    p = model.make_params()
    if n_comp == 1:
        tau_init = [t_span * 0.1]
    elif n_comp == 2:
        tau_init = [t_span * 0.03, t_span * 0.2]
    else:
        tau_init = [t_span * 0.01, t_span * 0.06, t_span * 0.25]

    for i, tau in enumerate(tau_init, 1):
        p[f"e{i}_decay"].value = tau

    t0 = time.perf_counter()
    res = model.fit(y, t=t, params=p)
    dt = time.perf_counter() - t0

    aic = calc_aic(res.chisqr, len(y), res.nvarys)
    bic = calc_bic(res.chisqr, len(y), res.nvarys)
    return res, aic, bic, dt


def fit_lm(y: np.ndarray, t: np.ndarray, n_comp: int) -> Any:
    """Fit n-component exponential + baseline using standard lmfit (Levenberg-Marquardt)."""
    t_span = float(t[-1] - t[0])
    c_init = float(y[-1])
    total_amp = float(y[0] - y[-1])

    if n_comp == 1:
        def ir1(t: np.ndarray, c: float, e1_amplitude: float, e1_decay: float) -> np.ndarray:
            return c + e1_amplitude * np.exp(-t / e1_decay)

        mod = lmfit.Model(ir1, independent_vars=["t"])
        p = mod.make_params(c=c_init, e1_amplitude=total_amp, e1_decay=t_span * 0.1)
        p["e1_decay"].min = 1e-4

    elif n_comp == 2:
        def ir2(
            t: np.ndarray,
            c: float,
            e1_amplitude: float,
            e1_decay: float,
            e2_amplitude: float,
            e2_decay: float,
        ) -> np.ndarray:
            return c + e1_amplitude * np.exp(-t / e1_decay) + e2_amplitude * np.exp(-t / e2_decay)

        mod = lmfit.Model(ir2, independent_vars=["t"])
        p = mod.make_params(
            c=c_init,
            e1_amplitude=total_amp * 0.5,
            e1_decay=t_span * 0.03,
            e2_amplitude=total_amp * 0.5,
            e2_decay=t_span * 0.2,
        )
        p["e1_decay"].min = 1e-4
        p["e2_decay"].min = 1e-4

    elif n_comp == 3:
        def ir3(
            t: np.ndarray,
            c: float,
            e1_amplitude: float,
            e1_decay: float,
            e2_amplitude: float,
            e2_decay: float,
            e3_amplitude: float,
            e3_decay: float,
        ) -> np.ndarray:
            return (
                c
                + e1_amplitude * np.exp(-t / e1_decay)
                + e2_amplitude * np.exp(-t / e2_decay)
                + e3_amplitude * np.exp(-t / e3_decay)
            )

        mod = lmfit.Model(ir3, independent_vars=["t"])
        p = mod.make_params(
            c=c_init,
            e1_amplitude=total_amp * 0.33,
            e1_decay=t_span * 0.01,
            e2_amplitude=total_amp * 0.33,
            e2_decay=t_span * 0.06,
            e3_amplitude=total_amp * 0.33,
            e3_decay=t_span * 0.25,
        )
        p["e1_decay"].min = 1e-4
        p["e2_decay"].min = 1e-4
        p["e3_decay"].min = 1e-4

    t0 = time.perf_counter()
    res = mod.fit(y, params=p, t=t, method="leastsq")
    dt = time.perf_counter() - t0
    return res, res.aic, res.bic, dt


def run_monte_carlo_test(
    y: np.ndarray, t: np.ndarray, n_comp: int = 3, n_trials: int = 50
) -> Dict[str, Any]:
    """Demonstrate basin of attraction and initial guess robustness."""
    rng = np.random.default_rng(42)
    t_span = float(t[-1] - t[0])
    total_amp = float(y[0] - y[-1])
    c_base = float(y[-1])

    vp_model = ConstantModel()
    for i in range(1, n_comp + 1):
        vp_model = vp_model + ExponentialModel(prefix=f"e{i}_")

    # First find target ground truth chisqr via VarPro with default guess
    res_base, _, _, _ = fit_varpro(y, t, n_comp)
    target_chi2 = float(res_base.chisqr)

    # Build LM model
    if n_comp == 2:
        def ir_mc(t: np.ndarray, c: float, e1_amplitude: float, e1_decay: float, e2_amplitude: float, e2_decay: float) -> np.ndarray:
            return c + e1_amplitude * np.exp(-t / e1_decay) + e2_amplitude * np.exp(-t / e2_decay)
        m_lm = lmfit.Model(ir_mc, independent_vars=["t"])
    else:
        def ir_mc(t: np.ndarray, c: float, e1_amplitude: float, e1_decay: float, e2_amplitude: float, e2_decay: float, e3_amplitude: float, e3_decay: float) -> np.ndarray:
            return c + e1_amplitude * np.exp(-t / e1_decay) + e2_amplitude * np.exp(-t / e2_decay) + e3_amplitude * np.exp(-t / e3_decay)
        m_lm = lmfit.Model(ir_mc, independent_vars=["t"])

    vp_evals, vp_success = [], 0
    lm_evals, lm_success = [], 0

    for _ in range(n_trials):
        if n_comp == 2:
            tau_guesses = [
                float(rng.uniform(t_span * 0.005, t_span * 0.1)),
                float(rng.uniform(t_span * 0.08, t_span * 0.5)),
            ]
            amp_guesses = [
                float(rng.uniform(total_amp * 1.5, total_amp * 0.2)),
                float(rng.uniform(total_amp * 1.5, total_amp * 0.2)),
            ]
        else:
            tau_guesses = [
                float(rng.uniform(t_span * 0.002, t_span * 0.04)),
                float(rng.uniform(t_span * 0.03, t_span * 0.15)),
                float(rng.uniform(t_span * 0.12, t_span * 0.6)),
            ]
            amp_guesses = [
                float(rng.uniform(total_amp * 1.2, total_amp * 0.1)),
                float(rng.uniform(total_amp * 1.2, total_amp * 0.1)),
                float(rng.uniform(total_amp * 1.2, total_amp * 0.1)),
            ]
        c_guess = float(rng.uniform(c_base * 0.85, c_base * 1.15))

        # VarPro
        p_vp = vp_model.make_params()
        for i, tau in enumerate(tau_guesses, 1):
            p_vp[f"e{i}_decay"].value = tau
        try:
            r_vp = vp_model.fit(y, t=t, params=p_vp)
            vp_evals.append(r_vp.nfev)
            if r_vp.chisqr <= target_chi2 * 1.05:
                vp_success += 1
        except Exception:
            pass

        # LM
        if n_comp == 2:
            p_lm = m_lm.make_params(
                c=c_guess,
                e1_amplitude=amp_guesses[0],
                e1_decay=tau_guesses[0],
                e2_amplitude=amp_guesses[1],
                e2_decay=tau_guesses[1],
            )
            p_lm["e1_decay"].min = 1e-4
            p_lm["e2_decay"].min = 1e-4
        else:
            p_lm = m_lm.make_params(
                c=c_guess,
                e1_amplitude=amp_guesses[0],
                e1_decay=tau_guesses[0],
                e2_amplitude=amp_guesses[1],
                e2_decay=tau_guesses[1],
                e3_amplitude=amp_guesses[2],
                e3_decay=tau_guesses[2],
            )
            p_lm["e1_decay"].min = 1e-4
            p_lm["e2_decay"].min = 1e-4
            p_lm["e3_decay"].min = 1e-4
        try:
            r_lm = m_lm.fit(y, params=p_lm, t=t, method="leastsq")
            lm_evals.append(r_lm.nfev)
            if r_lm.chisqr <= target_chi2 * 1.05:
                lm_success += 1
        except Exception:
            pass

    return {
        "n_trials": n_trials,
        "n_comp": n_comp,
        "vp_success_rate": vp_success / n_trials * 100.0,
        "vp_avg_evals": float(np.mean(vp_evals)),
        "lm_success_rate": lm_success / n_trials * 100.0,
        "lm_avg_evals": float(np.mean(lm_evals)),
    }


def plot_comparison_html(
    t_us: np.ndarray,
    data_dict: Dict[str, np.ndarray],
    fit_results: Dict[str, Dict[int, Any]],
    output_path: str = "epr_varpro_comparison.html",
) -> None:
    """Generate an interactive Plotly dashboard comparing all fits and residuals."""
    fig = make_subplots(
        rows=3,
        cols=2,
        subplot_titles=(
            "Raw Real Channel & Multi-Exp Fits",
            "Residuals (Raw Real)",
            "Phase-Corrected Real & Fits",
            "Residuals (Phase-Corrected Real)",
            "Raw Imaginary Channel & Fits",
            "Residuals (Raw Imaginary)",
        ),
        vertical_spacing=0.08,
        horizontal_spacing=0.08,
    )

    colors = {1: "#EF553B", 2: "#00CC96", 3: "#AB63FA"}
    curve_mapping = [
        ("Raw Real", "raw_re", 1),
        ("Phase-Corrected Real", "phased_re", 2),
        ("Raw Imaginary", "raw_im", 3),
    ]

    for cname, ckey, row_idx in curve_mapping:
        y = data_dict[ckey]
        # Data trace
        fig.add_trace(
            go.Scatter(
                x=t_us,
                y=y,
                mode="markers",
                name=f"{cname} Data",
                marker=dict(size=3, color="#636EFA", opacity=0.6),
                showlegend=(row_idx == 1),
            ),
            row=row_idx,
            col=1,
        )

        for n_comp in [1, 2, 3]:
            if n_comp in fit_results[cname]:
                r_vp = fit_results[cname][n_comp]["vp"]
                fig.add_trace(
                    go.Scatter(
                        x=t_us,
                        y=r_vp.best_fit,
                        mode="lines",
                        name=f"{n_comp}-Exp (AIC={r_vp.chisqr:.1e})",
                        line=dict(color=colors[n_comp], width=2),
                        showlegend=(row_idx == 1),
                    ),
                    row=row_idx,
                    col=1,
                )
                fig.add_trace(
                    go.Scatter(
                        x=t_us,
                        y=r_vp.residual,
                        mode="lines",
                        name=f"{n_comp}-Exp Residual",
                        line=dict(color=colors[n_comp], width=1.5),
                        showlegend=False,
                    ),
                    row=row_idx,
                    col=2,
                )

    fig.update_layout(
        height=1000,
        title_text="<b>EPR Inversion Recovery Lifetimes: VarPro vs. LM Model Comparison</b>",
        template="plotly_white",
    )
    fig.write_html(output_path)
    print(f"Saved interactive Plotly dashboard to: {output_path}")


def run_analysis(dsc_file: str, dta_file: str, dataset_label: str) -> None:
    print("\n" + "=" * 88)
    print(f"      EPR LIFETIME FITTING: {dataset_label}")
    print("      STANDARD LM (lmfit) vs. VARIABLE PROJECTION (VarProPy)")
    print("=" * 88)

    data = load_bruker_data(dsc_file, dta_file)
    t = data["t_us"]
    print(f"Data Loaded: {len(t)} points, time delay from {t[0]:.2f} to {t[-1]:.2f} \u03bcs")
    print(f"Optimal detector phase correction angle: {data['phi_deg']:.3f}\u00b0")

    curves = {
        "Raw Real": data["raw_re"],
        "Phase-Corrected Real": data["phased_re"],
        "Raw Imaginary": data["raw_im"],
    }

    all_fits: Dict[str, Dict[int, Any]] = {c: {} for c in curves}

    for cname, y in curves.items():
        print("\n" + "#" * 88)
        print(f"  DATASET: {dataset_label} | CURVE: {cname}")
        print("#" * 88)

        for n_comp in [1, 2, 3]:
            # Fit VarPro
            res_vp, aic_vp, bic_vp, dt_vp = fit_varpro(y, t, n_comp)
            # Fit LM
            res_lm, aic_lm, bic_lm, dt_lm = fit_lm(y, t, n_comp)

            all_fits[cname][n_comp] = {
                "vp": res_vp,
                "lm": res_lm,
                "aic_vp": aic_vp,
                "aic_lm": aic_lm,
                "bic_vp": bic_vp,
                "bic_lm": bic_lm,
            }

            print(f"\n--- Model: {n_comp}-Exponential Decay + Baseline ---")
            print(f"{'Statistic':<25} | {'VarProPy (VarPro)':<26} | {'lmfit (Standard LM)':<26}")
            print("-" * 85)
            print(f"{'Chi-square (RSS)':<25} | {res_vp.chisqr:15.5e}            | {res_lm.chisqr:15.5e}")
            print(f"{'Reduced Chi-square':<25} | {res_vp.redchi:15.5e}            | {res_lm.redchi:15.5e}")
            print(f"{'Akaike Info Crit (AIC)':<25} | {aic_vp:15.2f}            | {aic_lm:15.2f}")
            print(f"{'Bayesian Info Crit (BIC)':<25} | {bic_vp:15.2f}            | {bic_lm:15.2f}")
            print(f"{'Function Evaluations':<25} | {res_vp.nfev:15d}            | {res_lm.nfev:15d}")
            print(f"{'Execution Time':<25} | {dt_vp*1000:13.2f} ms          | {dt_lm*1000:13.2f} ms")

            print("\n  Fitted Parameters & Standard Errors:")
            print(f"  {'Parameter':<16} | {'VarPro Value +/- Stderr':<30} | {'LM Value +/- Stderr':<30}")
            print("  " + "-" * 80)
            for pname in res_vp.params.keys():
                vp_p = res_vp.params[pname]
                lm_p = res_lm.params[pname]
                vp_str = f"{vp_p.value:11.4f} +/- {vp_p.stderr:9.4f}" if vp_p.stderr else f"{vp_p.value:11.4f}"
                lm_str = f"{lm_p.value:11.4f} +/- {lm_p.stderr:9.4f}" if lm_p.stderr else f"{lm_p.value:11.4f}"
                pct = f"({vp_p.stderr / abs(vp_p.value) * 100:.1f}%)" if vp_p.stderr and vp_p.value else ""
                print(f"  {pname:<16} | {vp_str:<30} | {lm_str:<30} {pct}")

            print("\n  Parameter Cross-Correlation Matrix (VarPro):")
            vnames = res_vp.var_names
            header = "      " + "".join([f"{v:>13}" for v in vnames])
            print(header)
            for v1 in vnames:
                row = f"  {v1:<8}"
                for v2 in vnames:
                    if v1 == v2:
                        row += f"{'1.000':>13}"
                    elif v2 in res_vp.params[v1].correl:
                        row += f"{res_vp.params[v1].correl[v2]:>13.3f}"
                    else:
                        row += f"{'--':>13}"
                print(row)

    print("\n" + "=" * 88)
    print(f"      ROBUSTNESS & BASIN OF ATTRACTION ({dataset_label}, 50 RANDOMIZED TRIALS)")
    print("=" * 88)
    for n_c in [2, 3]:
        mc = run_monte_carlo_test(data["raw_re"], t, n_comp=n_c, n_trials=50)
        print(f"\n{n_c}-Exponential Model (Randomized Initial Guesses, N=50):")
        print(f"  VarPro: Global Minimum Convergence Rate = {mc['vp_success_rate']:.1f}%  |  Avg Evals = {mc['vp_avg_evals']:.1f}")
        print(f"  LM:     Global Minimum Convergence Rate = {mc['lm_success_rate']:.1f}%  |  Avg Evals = {mc['lm_avg_evals']:.1f}")

    safe_name = dataset_label.replace(" ", "_").lower()
    html_file = f"epr_varpro_{safe_name}.html"
    plot_comparison_html(t, data, all_fits, html_file)
    print(f"\nSaved dashboard to: {html_file}")


def main() -> None:
    datasets = [
        ("PdSC12_Tol_InvRec_tau200_24pi2_10K_9.473067.DSC", "PdSC12_Tol_InvRec_tau200_24pi2_10K_9.473067.DTA", "10K Dataset"),
        ("PdSC12_Tol InvRec tau200 24pi2 40K 9.473067.DSC", "PdSC12_Tol InvRec tau200 24pi2 40K 9.473067.DTA", "40K Dataset (Noisy)"),
    ]

    target = sys.argv[1] if len(sys.argv) > 1 else None
    for dsc, dta, label in datasets:
        if target is None or target.lower() in label.lower() or target.lower() in dsc.lower():
            if os.path.exists(dsc) and os.path.exists(dta):
                run_analysis(dsc, dta, label)


if __name__ == "__main__":
    main()
