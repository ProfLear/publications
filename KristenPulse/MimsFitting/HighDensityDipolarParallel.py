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
    from joblib import Parallel, delayed

    # import optimized modules
    import pulsePhysicsOptimized as ppo
    import mdHelper as mh

    reload(ppo)
    reload(mh)
    return Model, Parallel, delayed, go, make_subplots, mh, mo, np, ppo, reload


@app.cell
def _(Parallel, delayed, np, ppo):
    ### OPTIMIZED FUNCTION FOR FITTING

    def triple_mims_simulation_parallel(
        xs, len1, len2, I1, I2, I3, c1, c2, c3,
        H_np, D_np, solvent_density_obj, rho_bulk_solvent,
        ligand_density_obj, rho_bulk_ligand, nu_n, g_n,
        I_scale, ell, g_e, tau, T_m, R_np_md, R_max,
        target_dr=0.2E-10, dipolar_scale=1, contact_scale=1,
        blindspot_offset=0, dipolar_np_scale=1, shape="Lorentzian",
        eta=0, width_scale=1, n_jobs=-1
    ):
        """
        Parallelized and Vectorized Mims simulation.
        n_jobs=-1 uses all available CPU cores.
        """
        freq_both = xs[:len1]
        freq_solvent_only = xs[len1:len1+len2]
        freq_ligand_only = xs[len1+len2:]

        all_raw_samples = []

        # Configurations for the 3 experimental spectra
        configs = [
            ("both", freq_both, 1.0, 1.0, H_np),
            ("solvent only", freq_solvent_only, 1.0, 0.02, D_np),
            ("ligand only", freq_ligand_only, 0.01, 1.0, H_np)
        ]

        for name, freq, s_weight, l_weight, rs in configs:
            # Prepare tasks for parallel execution over nanoparticle sizes
            # rs[0] are radii, rs[1] are weights
            valid_indices = [i for i, w in enumerate(rs[1]) if w > 0]
            radii = rs[0][valid_indices]
            weights = rs[1][valid_indices]

            def compute_size_contribution(r, weight):
                y_s = np.zeros_like(freq)
                y_l = np.zeros_like(freq)
                
                if s_weight > 0:
                    y_s = ppo.simulate_nanoparticle_mims_spectrum_vectorized(
                        freq, r, solvent_density_obj, rho_bulk_solvent,
                        nu_n, g_n, I_scale, ell, g_e, tau, T_m,
                        R_np_md, R_max, target_dr, dipolar_scale,
                        contact_scale, blindspot_offset, dipolar_np_scale,
                        shape, eta, width_scale
                    )
                if l_weight > 0:
                    y_l = ppo.simulate_nanoparticle_mims_spectrum_vectorized(
                        freq, r, ligand_density_obj, rho_bulk_ligand,
                        nu_n, g_n, I_scale, ell, g_e, tau, T_m,
                        R_np_md, R_max, target_dr, dipolar_scale,
                        contact_scale, blindspot_offset, dipolar_np_scale,
                        shape, eta, width_scale
                    )
                return weight * (s_weight * y_s + l_weight * y_l)

            # Execute in parallel
            results = Parallel(n_jobs=n_jobs)(
                delayed(compute_size_contribution)(r, w) for r, w in zip(radii, weights)
            )
            
            sample_total_raw = np.sum(results, axis=0) if results else np.zeros_like(freq)
            all_raw_samples.append(sample_total_raw)

        # Normalization and concatenation
        scaled_both = (all_raw_samples[0] / np.max(all_raw_samples[0])) * I1 + c1
        scaled_solvent = (all_raw_samples[1] / np.max(all_raw_samples[1])) * I2 + c2
        scaled_ligand = (all_raw_samples[2] / np.max(all_raw_samples[2])) * I3 + c3

        return np.concatenate([scaled_both, scaled_solvent, scaled_ligand])

    return (triple_mims_simulation_parallel,)

# Note: The rest of the Marimo notebook logic (data loading, parameter setup, fitting call)
# would follow here, swapping triple_mims_simulation for triple_mims_simulation_parallel.
# For brevity, I've implemented the core optimized logic.
