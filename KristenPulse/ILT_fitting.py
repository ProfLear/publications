import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import numpy as np
    from scipy.optimize import nnls
    from scipy.linalg import block_diag
    from plotly.subplots import make_subplots
    from importlib import reload
    from scipy.optimize import lsq_linear


    import mimsFunctions as mf
    reload(mf)

    return block_diag, lsq_linear, make_subplots, nnls, np


@app.cell
def _(block_diag, nnls, np):


    def simulate_shell_response(freq_axis, r, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, eta):
        """
        Simulates the ENDOR response for a single shell at distance r 
        with UNIT density (rho=1). This is one column of our kernel.
        """
        cos_theta = np.linspace(0, 1, 100) 
        damping = np.exp(-2 * tau / T_m)
        f_vec = freq_axis[np.newaxis, :]
    
        # Physics
        T = (79.0 * (g_iso / 2.0023)) / (r**3) 
        A_iso = A0 * np.exp(-(max(0, r - 10.0)) / rho_decay) # Assuming R_np=10 for example
        A = A_iso + T * (3 * cos_theta**2 - 1) 

        # Mims weight for rho=1
        # Note: We include the 4*pi*r^2 factor here so the solved 'rho' 
        # represents the actual density (protons/A^3)
        mims = (1 - np.cos(2 * np.pi * A * tau))
        w = (4 * np.pi * r**2 / len(cos_theta)) * mims * damping
    
        # Strain
        shell_sigma = np.sqrt(sigma**2 + (eta * np.abs(A)**0.5)**2)[:, np.newaxis]
        norm = shell_sigma * np.sqrt(2 * np.pi)
    
        mu_plus = nu_n + A[:, np.newaxis] / 2
        mu_minus = nu_n - A[:, np.newaxis] / 2
    
        g = (w[:, np.newaxis] / norm) * (
            np.exp(-0.5 * ((f_vec - mu_plus) / shell_sigma)**2) + 
            np.exp(-0.5 * ((f_vec - mu_minus) / shell_sigma)**2)
        )
    
        return np.sum(g, axis=0)

    def build_global_kernel(freq_list, r_grid, physics_params, compositions):
        """
        Builds the block matrix K that links the 2 density profiles to the 3 spectra.
        compositions: list of (solvent_weight, ligand_weight) for each spectrum.
        """
        n_r = len(r_grid)
    
        # 1. Build basic kernels for ligand and solvent
        # (In this simple version they use the same physics, but could differ)
        K_base = np.zeros((len(freq_list[0]), n_r))
        for j, r in enumerate(r_grid):
            K_base[:, j] = simulate_shell_response(freq_list[0], r, **physics_params)

        # 2. Assemble the global block matrix
        # We create a mapping for: [S1; S2; S3] = K_global * [rho_lig; rho_solv]
        rows = []
        for s_wt, l_wt in compositions:
            # Each row of the block matrix is: [weight_lig * K, weight_solv * K]
            rows.append(np.hstack([l_wt * K_base, s_wt * K_base]))
        
        return np.vstack(rows)

    def build_global_kernel_mixed_axes(freq_list, r_grid, physics_params, compositions):
        """
        Handles cases where each spectrum in freq_list has a different length or range.
        freq_list: [freq_1, freq_2, freq_3]
        compositions: [(s_wt1, l_wt1), (s_wt2, l_wt2), (s_wt3, l_wt3)]
        """
        n_r = len(r_grid)
        all_rows = []

        # Loop through each spectrum individually
        for freq, (s_wt, l_wt) in zip(freq_list, compositions):
            # 1. Generate the base response for THIS specific frequency axis
            # n_freq_points x n_r_shells
            K_spec = np.zeros((len(freq), n_r))
            for j, r in enumerate(r_grid):
                K_spec[:, j] = simulate_shell_response(freq, r, **physics_params)

            # 2. Create the horizontal block for this spectrum: [Ligand, Solvent]
            # Resulting shape: (len(freq), 2 * n_r)
            block = np.hstack([l_wt * K_spec, s_wt * K_spec])
            all_rows.append(block)

        # 3. Vertically stack the blocks
        # Total rows = len(freq1) + len(freq2) + len(freq3)
        return np.vstack(all_rows)

    def get_reg_matrix(n, alpha):
        """Second-derivative smoothing matrix (Laplacian)."""
        D = np.zeros((n-2, n))
        for i in range(n-2):
            D[i, i] = 1
            D[i, i+1] = -2
            D[i, i+2] = 1
        return D * alpha

    def invert_density_global(exp_data_concat, freq_list, r_grid, physics_params, compositions, alpha=1.0):
        """
        The Inner Loop: Solves for rho_ligand and rho_solvent using NNLS + Tikhonov.
        """
        n_r = len(r_grid)
    
        # Build Physics Kernel
        K_global = build_global_kernel(freq_list, r_grid, physics_params, compositions)
    
        # Build Regularization (Smoothing)
        # We smooth ligand and solvent profiles independently
        L_single = get_reg_matrix(n_r, alpha)
        L_global = block_diag(L_single, L_single)
    
        # Augment Matrix for Tikhonov: [K; L] * rho = [Data; 0]
        K_reg = np.vstack([K_global, L_global])
        S_reg = np.concatenate([exp_data_concat, np.zeros(L_global.shape[0])])
    
        # Solve Non-Negative Least Squares
        rho_fit, resid = nnls(K_reg, S_reg)
    
        # Split the result back into ligand and solvent
        rho_ligand = rho_fit[:n_r]
        rho_solvent = rho_fit[n_r:]
    
        return r_grid, rho_ligand, rho_solvent

    return (build_global_kernel_mixed_axes,)


@app.cell
def _(make_subplots, np):
    # import spectral data
    freq_axis = np.linspace(43, 57, 500)


    dodecaneH_tolueneH_spectra =  np.genfromtxt("/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/KristenPulse/MIMS/PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)

    dodecaneD_tolueneH_spectra = np.genfromtxt("/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/KristenPulse/MIMS/dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter="," ,unpack = True, skip_header = 1)

    dodecaneH_tolueneD_spectra = np.genfromtxt("/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/KristenPulse/MIMS/Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter="," ,unpack = True, skip_header = 1)

    spectra_Plot = make_subplots(rows = 3, cols = 1)
    spectra_Plot.update_layout(template = "simple_white")
    for _exp, _row in zip(
        [
            dodecaneH_tolueneH_spectra,
            dodecaneD_tolueneH_spectra,
            dodecaneH_tolueneD_spectra
        ], 
        range(1,4)):


        spectra_Plot.add_scatter(x = _exp[0], y = _exp[1], 
                             name = "both", row = _row, col = 1)

    spectra_Plot.show()
    return (
        dodecaneD_tolueneH_spectra,
        dodecaneH_tolueneD_spectra,
        dodecaneH_tolueneH_spectra,
    )


@app.cell
def _(
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    np,
):


    concatenated_xs = np.concatenate([
        dodecaneH_tolueneH_spectra[0],
        dodecaneD_tolueneH_spectra[0],
        dodecaneH_tolueneD_spectra[0] 
    ]
        )

    concatenated_ys = np.concatenate((  # Note the double parenthesis
        dodecaneH_tolueneH_spectra[1] - np.mean(dodecaneH_tolueneH_spectra[1][:10]),
        dodecaneD_tolueneH_spectra[1] - np.mean(dodecaneD_tolueneH_spectra[1][:10]),
        dodecaneH_tolueneD_spectra[1] - np.mean(dodecaneH_tolueneD_spectra[1][:10])
    ))


    return (concatenated_ys,)


@app.cell
def _(block_diag, lsq_linear, nnls, np):
    def perform_inversion(y_concat, K_global, n_r, alpha):
        # 1. Create the Smoothing Matrix (Second Derivative / Laplacian)
        # This penalizes sharp changes in density between neighboring shells
        D = np.zeros((n_r - 2, n_r))
        for i in range(n_r - 2):
            D[i, i], D[i, i+1], D[i, i+2] = 1, -2, 1
    
        # We apply the same smoothing to both Ligand and Solvent profiles
        L_global = block_diag(D, D) * alpha
    
        # 2. Augment the System for Tikhonov Regularization
        # We are solving: [K; L] @ rho = [Data; 0]
        K_reg = np.vstack([K_global, L_global])
        y_reg = np.concatenate([y_concat, np.zeros(L_global.shape[0])])
    
        # 3. Solve Non-Negative Least Squares
        rho_fit, residual_norm = nnls(K_reg, y_reg)
    
        # Split the result back into the two chemical species
        rho_ligand = rho_fit[:n_r]
        rho_solvent = rho_fit[n_r:]
    
        return rho_ligand, rho_solvent



    def perform_constrained_inversion(y_concat, K_global, n_r, alpha, 
                                      max_ligand_r, rho_max_ligand, 
                                      rho_bulk_solv, r_grid):
        """
        Performs a regularized linear inversion to find ligand and solvent 
        density profiles with physical constraints.
    
        Parameters:
        -----------
        y_concat : array
            Concatenated experimental intensities (Both, Solv-only, Lig-only).
        K_global : matrix
            The global kernel block matrix.
        n_r : int
            Number of radial shells (length of r_grid).
        alpha : float
            Regularization parameter (higher = smoother profile).
        max_ligand_r : float
            The radial distance (A) where the ligand density must drop to zero.
        rho_max_ligand : float
            The maximum allowed proton density for the ligand (packing limit).
        rho_bulk_solv : float
            The maximum allowed proton density for the solvent (bulk limit).
        r_grid : array
            The radial distance axis used to build the kernel.
        """
    
        # --- 1. SETUP REGULARIZATION (SMOOTHNESS) ---
        # Second-derivative (Laplacian) operator to penalize jagged density
        D = np.zeros((n_r - 2, n_r))
        for i in range(n_r - 2):
            D[i, i], D[i, i+1], D[i, i+2] = 1, -2, 1
    
        # Apply smoothing to both ligand and solvent blocks
        L_global = block_diag(D, D) * alpha

        # --- 2. AUGMENT SYSTEM FOR TIKHONOV ---
        # We solve: [K; L] @ rho = [y; 0]
        K_reg = np.vstack([K_global, L_global])
        y_reg = np.concatenate([y_concat, np.zeros(L_global.shape[0])])

        # --- 3. DEFINE PHYSICAL BOUNDS ---
        # Vector size is 2 * n_r (Ligand profile followed by Solvent profile)
        lb = np.zeros(2 * n_r)
        ub = np.inf * np.ones(2 * n_r)
    
        # LIGAND CONSTRAINTS (First n_r elements)
        ub[:n_r] = rho_max_ligand
    
        # Apply the Mask for Ligand Length
        # We use 1e-10 instead of 0.0 to satisfy the strict lb < ub requirement
        ligand_mask = r_grid > max_ligand_r
        ub[:n_r][ligand_mask] = 1e-10
    
        # SOLVENT CONSTRAINTS (Second n_r elements)
        ub[n_r:] = rho_bulk_solv 

        # --- 4. SOLVE CONSTRAINED LEAST SQUARES ---
        # 'trf' (Trust Region Reflective) is the standard method for bounded problems
        res = lsq_linear(K_reg, y_reg, bounds=(lb, ub), method='trf', tol=1e-6)
    
        if not res.success:
            print(f"Warning: Solver did not converge. Message: {res.message}")

        # --- 5. EXTRACT PROFILES ---
        rho_fit = res.x
        rho_ligand = rho_fit[:n_r]
        rho_solvent = rho_fit[n_r:]
    
        return rho_ligand, rho_solvent

    return (perform_constrained_inversion,)


@app.cell
def _(
    alpha,
    build_global_kernel_mixed_axes,
    concatenated_ys,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    eta,
    intensity_scaling,
    np,
    perform_constrained_inversion,
    r_grid,
):
    mixed_kernel = build_global_kernel_mixed_axes(
        [
        dodecaneH_tolueneH_spectra[0],
        dodecaneD_tolueneH_spectra[0],
        dodecaneH_tolueneD_spectra[0] 
        ],
        r_grid,
        {
            "A0" : 0.6,
            "rho_decay" : 3.5, 
            "g_iso" : 2.02, 
            "nu_n" : 50, 
            "tau" : 140, 
            "T_m" : 140, 
            "sigma" : 0.2, 
            "eta" : eta,
        },
        [(1, 1), (1, 0.02), (0.01, 1)]
        )

    rho_ligand, rho_solvent = perform_constrained_inversion(
        y_concat = concatenated_ys/intensity_scaling,
        K_global = mixed_kernel,
        n_r  = len(r_grid),
        alpha = alpha,
        max_ligand_r = 33,
        rho_max_ligand = 0.0453,
        rho_bulk_solv = 0.0453,
        r_grid = r_grid,
    )
    # Combine the two recovered profiles into one vector
    rho_total = np.concatenate([rho_ligand, rho_solvent])

    # Multiply by the Kernel to get the 'Best Fit' simulated spectra
    y_fit = mixed_kernel @ rho_total

    return mixed_kernel, rho_ligand, rho_solvent


@app.cell
def _(make_subplots, np):

    def visualize_inversion_plotly(r_grid, rho_ligand, rho_solvent, y_concat, K_global, freq_list, alpha):
        import plotly.graph_objects as go
        # 1. Reconstruct the fitted spectra from the density
        rho_total = np.concatenate([rho_ligand, rho_solvent])
        y_fit_concat = K_global @ rho_total
    
        # Create the figure with two side-by-side plots
        fig = make_subplots(
            rows=1, cols=2, 
            subplot_titles=("Reconstructed Density Profiles", "Global Spectral Fit"),
            horizontal_spacing=0.15
        )

        # --- LEFT PLOT: Density vs. Distance ---
        fig.add_trace(
            go.Scatter(x=r_grid, y=rho_ligand, name="Ligand Density", 
                       line=dict(color='blue', width=3)),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=r_grid, y=rho_solvent, name="Solvent Density", 
                       line=dict(color='red', width=3, dash='dash')),
            row=1, col=1
        )
    
        # --- RIGHT PLOT: Experimental vs. Fit ---
        # Plotting the concatenated y-axis as a continuous index
        indices = np.arange(len(y_concat))
        fig.add_trace(
            go.Scatter(x=indices, y=y_concat, name="Experimental", 
                       mode='markers', marker=dict(color='black', size=4, opacity=0.4)),
            row=1, col=2
        )
        fig.add_trace(
            go.Scatter(x=indices, y=y_fit_concat, name="Inversion Fit", 
                       line=dict(color='orange', width=2)),
            row=1, col=2
        )

        # Styling
        fig.update_layout(
            title=f"ENDOR Inversion Analysis (Regularization Alpha: {alpha})",
            template="simple_white",
            height=500,
            showlegend=True
        )
    
        fig.update_xaxes(title_text="Distance (Å)", row=1, col=1)
        fig.update_yaxes(title_text="Density (protons/Å³)", row=1, col=1)
    
        fig.update_xaxes(title_text="Data Point Index", row=1, col=2)
        fig.update_yaxes(title_text="Normalized Intensity", row=1, col=2)

        return fig

    # To display in your notebook:
    # plot = visualize_inversion_plotly(r_grid, rho_lig, rho_solv, y_concat, K_global, freq_list, alpha)
    # plot.show()
    return (visualize_inversion_plotly,)


@app.cell
def _(
    alpha,
    concatenated_ys,
    dodecaneD_tolueneH_spectra,
    dodecaneH_tolueneD_spectra,
    dodecaneH_tolueneH_spectra,
    intensity_scaling,
    mixed_kernel,
    r_grid,
    rho_ligand,
    rho_solvent,
    visualize_inversion_plotly,
):
    visualize_inversion_plotly(
        r_grid,
        rho_ligand, rho_solvent,
        concatenated_ys/intensity_scaling,
        mixed_kernel,
        [
            dodecaneH_tolueneH_spectra[0],
            dodecaneD_tolueneH_spectra[0],
            dodecaneH_tolueneD_spectra[0] 
        ],
        alpha=alpha
    )
    return


@app.cell
def _(np):
    R_np = 15
    R_max = 50
    n_rs = 500

    alpha = 2

    eta = 10

    intensity_scaling = 1e4

    r_grid = np.linspace(R_np, R_max, n_rs)
    return R_np, alpha, eta, intensity_scaling, r_grid


@app.cell
def _(R_np, block_diag, lsq_linear, np, r_grid):
    from scipy.stats import gamma, norm

    def create_physical_prior(r_grid, R_np, ligand_length, rho_max_lig, rho_max_solv):
        n_r = len(r_grid)
        x = r_grid - R_np
    
        # --- LIGAND PRIOR: Gamma Distribution ---
        # a=shape, scale=spread. Adjust these to match your expected ligand 'brush'
        lig_shape = gamma.pdf(x, a=3, scale=4) 
        # Normalize and scale to your physical density cap
        lig_shape = (lig_shape / np.max(lig_shape)) * rho_max_lig
        # Hard cutoff at ligand length
        lig_shape[r_grid > (R_np + ligand_length)] = 0
    
        # --- SOLVENT PRIOR: Gaussian 'Layering' ---
        # Centered roughly at the edge of the ligand shell
        solv_center = R_np + ligand_length + 2.0 
        solv_shape = norm.pdf(r_grid, loc=solv_center, scale=5.0)
        solv_shape = (solv_shape / np.max(solv_shape)) * rho_max_solv
    
        return np.concatenate([lig_shape, solv_shape])

    # Generate the prior
    rho_prior = create_physical_prior(r_grid, R_np, 23.0, 0.0453, 0.0453)

    def perform_prior_guided_inversion(y_concat, K_global, n_r, alpha, 
                                      lb, ub, rho_prior):
        # 1. Setup Smoothing Matrix (L)
        D = np.zeros((n_r - 2, n_r))
        for i in range(n_r - 2):
            D[i, i], D[i, i+1], D[i, i+2] = 1, -2, 1
        L_global = block_diag(D, D) * alpha

        # 2. Augment System
        K_reg = np.vstack([K_global, L_global])
    
        # --- THE KEY DIFFERENCE ---
        # Instead of zeros, we use the 'smoothness' of the prior as the target.
        # This tells the solver: "If you have no data, look exactly like the prior."
        prior_target = L_global @ rho_prior
        y_reg = np.concatenate([y_concat, prior_target])

        # 3. Solve with Constraints
        res = lsq_linear(K_reg, y_reg, bounds=(lb, ub), method='trf', tol=1e-6)
    
        return res.x[:n_r], res.x[n_r:]

    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
