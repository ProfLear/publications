import numpy as np
from scipy.optimize import minimize
from scipy.linalg import block_diag
from scipy.stats import lognorm
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os

# --- Configuration & Parameters ---
R_NP = 15.0       # Nanoparticle radius (A)
R_MAX = 70.0      # Maximum radial distance to consider (A)
N_R = 300         # Number of radial shells in the grid

# --- Intensity Scaling Factor ---
# Divide experimental data by this to match the scale of the physical model (protons/A^3)
INTENSITY_SCALING = 2e6

# Physical Parameters for Kernel
PHYSICS_PARAMS = {
    "A0": 8.0,         # Contact hyperfine at surface (MHz)
    "rho_decay": 1.7,   # Decay length for contact hyperfine (A)
    "g_iso": 2.02,      # g-factor
    "nu_n": 50.0,       # Nuclear Larmor frequency (MHz)
    "tau": 0.140,       # Mims tau delay (us)
    "T_m": 1.4,         # Phase memory time (us)
    "sigma": 0.2,       # Inhomogeneous broadening (MHz)
    "eta": 1.3     # Strain parameter
}

# Parameters for Parametric Search
RHO_BULK = 0.0453       # Bulk proton density (protons/A^3)
LIGAND_LENGTH = 20.0    # Expected ligand length for truncation (A)

# Data Files
BASE_PATH = "/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/KristenPulse"
DATA_PATHS = {
    "both": os.path.join(BASE_PATH, "MIMS/PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv"),
    "solvent_h": os.path.join(BASE_PATH, "MIMS/dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv"),
    "ligand_h": os.path.join(BASE_PATH, "MIMS/Mims PdSC12 dTol tau 140 10K 43K scans.csv")
}

# Compositions (Ligand H fraction, Solvent H fraction)
COMPOSITIONS = [
    (1.0, 1.0),    # Both
    (0.02, 1.0),   # Solvent-H (D-ligand)
    (1.0, 0.01)    # Ligand-H (D-solvent)
]

# --- Core Functions ---

def simulate_shell_response(freq_axis, r, A0, rho_decay, g_iso, nu_n, tau, T_m, sigma, eta):
    """
    Simulates the ENDOR response for a single shell at distance r with UNIT density.
    """
    cos_theta = np.linspace(0, 1, 100) 
    damping = np.exp(-2 * tau / T_m)
    f_vec = freq_axis[np.newaxis, :]

    # Physics
    T = (79.0 * (g_iso / 2.0023)) / (r**3) 
    dist_from_surface = max(0, r - R_NP)
    A_iso = A0 * np.exp(-dist_from_surface / rho_decay)
    A = A_iso + T * (3 * cos_theta**2 - 1) 

    # Mims weight (including 4*pi*r^2 so rho is density)
    mims = (1 - np.cos(2 * np.pi * A * tau))
    w = (4 * np.pi * r**2 / len(cos_theta)) * mims * damping

    # Strain-dependent broadening
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
    Builds the block matrix K linking density profiles to multiple spectra.
    """
    n_r = len(r_grid)
    all_rows = []

    for freq, (l_wt, s_wt) in zip(freq_list, compositions):
        K_spec = np.zeros((len(freq), n_r))
        for j, r in enumerate(r_grid):
            K_spec[:, j] = simulate_shell_response(freq, r, **physics_params)

        # Block for this spectrum: [Ligand_contribution, Solvent_contribution]
        block = np.hstack([l_wt * K_spec, s_wt * K_spec])
        all_rows.append(block)

    return np.vstack(all_rows)

def perform_parametric_fit(y_concat, K_global, r_grid, r_np):
    """
    Fits the data using exactly one Log-Normal (Ligand) and one Logistic (Solvent).
    This guarantees a single, smooth distribution with no artifacts.
    """
    n_r = len(r_grid)
    x = r_grid - r_np
    
    def get_rho(params):
        amp_l, s, loc, sc_l, amp_s, k, r0 = params
        
        # SINGLE Log-Normal Distribution
        # loc < 0 ensures non-zero density at the surface (x=0)
        rho_l = amp_l * lognorm.pdf(x, s=s, loc=loc, scale=sc_l)
        
        # Hard physical truncation beyond ligand length (solver only)
        rho_l[x > LIGAND_LENGTH] = 0
        
        # Logistic Solvent Rise
        rho_s = amp_s / (1 + np.exp(-k * (r_grid - r0)))
        
        return np.concatenate([rho_l, rho_s])

    def objective(params):
        rho = get_rho(params)
        y_fit = K_global @ rho
        return np.sum((y_fit - y_concat)**2)

    # Initial Guesses:
    # amp_l, s, loc, sc_l, amp_s, k, r0
    init = [0.1, 0.8, -1.0, 6.0, RHO_BULK, 0.5, r_np + LIGAND_LENGTH]
    
    # Bounding to ensure physical results and surface density
    bnds = [
        (0, None),       # amp_l
        (0.1, 2.0),      # s
        (-5.0, -0.1),    # loc (MUST be negative to ensure intensity at surface)
        (1.0, 20.0),     # sc_l
        (0, RHO_BULK*2), # amp_s
        (0.1, 2.0),      # k
        (r_np, R_MAX)    # r0
    ]
    
    print("Optimizing parametric fit (Single Log-Normal)...")
    res = minimize(objective, init, bounds=bnds, method='L-BFGS-B')
    
    if not res.success:
        print(f"Optimization warning: {res.message}")
        
    best_params = res.x
    rho_fit = get_rho(best_params)
    return rho_fit[:n_r], rho_fit[n_r:], best_params

# --- Main Script ---

if __name__ == "__main__":
    print("Loading data...")
    data = {}
    valid_keys = ["both", "solvent_h", "ligand_h"]
    for key in valid_keys:
        path = DATA_PATHS[key]
        if not os.path.exists(path):
            path = os.path.basename(path)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Could not find data for {key} at {path}")
        data[key] = np.genfromtxt(path, delimiter=",", unpack=True, skip_header=1)

    # Baseline correction and Scaling
    freq_list = [data["both"][0], data["solvent_h"][0], data["ligand_h"][0]]
    y_list = []
    for key in valid_keys:
        y = data[key][1]
        y_corr = (y - np.mean(y[:10])) / INTENSITY_SCALING
        y_list.append(y_corr)
    
    y_concat = np.concatenate(y_list)
    
    # Define grid
    r_grid = np.linspace(R_NP, R_MAX, N_R)
    
    print("Building kernels...")
    K_global = build_global_kernel(freq_list, r_grid, PHYSICS_PARAMS, COMPOSITIONS)
    
    print("Performing parametric inversion...")
    rho_lig_fit, rho_solv_fit, best_p = perform_parametric_fit(y_concat, K_global, r_grid, R_NP)
    
    # --- Visualization ---
    print("Plotting results...")
    rho_total_fit = np.concatenate([rho_lig_fit, rho_solv_fit])
    y_fit = K_global @ rho_total_fit
    
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Single Log-Normal & Logistic Profiles", "Global Spectral Fit"),
        vertical_spacing=0.15
    )
    
    # Density plot
    fig.add_trace(go.Scatter(x=r_grid, y=rho_lig_fit, name="Ligand (Single Log-Normal)", line=dict(color='blue', width=3)), row=1, col=1)
    fig.add_trace(go.Scatter(x=r_grid, y=rho_solv_fit, name="Solvent (Logistic)", line=dict(color='red', width=3)), row=1, col=1)
    
    # Marker for LIGAND_LENGTH
    fig.add_vline(x=R_NP + LIGAND_LENGTH, line_dash="dash", line_color="blue", annotation_text="Max Ligand Length", row=1, col=1)
    
    # Spectral fit plot
    indices = np.arange(len(y_concat))
    fig.add_trace(go.Scatter(x=indices, y=y_concat, name="Experimental (Scaled)", mode='markers', marker=dict(size=2, color='black', opacity=0.5)), row=2, col=1)
    fig.add_trace(go.Scatter(x=indices, y=y_fit, name="Fit", line=dict(color='orange', width=2)), row=2, col=1)
    
    # Dividers for concatenated spectra
    sep1 = len(data["both"][0])
    sep2 = sep1 + len(data["solvent_h"][0])
    for sep in [sep1, sep2]:
        fig.add_vline(x=sep, line_dash="dash", line_color="gray", row=2, col=1)
    
    fig.update_layout(
        height=900, 
        title_text=f"Parametric Fit: s={best_p[1]:.2f}, loc={best_p[2]:.2f}, scale={best_p[3]:.2f}",
        template="simple_white"
    )
    fig.update_xaxes(title_text="Distance (A)", row=1, col=1)
    fig.update_yaxes(title_text="Density (H/A^3)", row=1, col=1)
    fig.update_xaxes(title_text="Data Point Index", row=2, col=1)
    fig.update_yaxes(title_text="Intensity", row=2, col=1)
    
    fig.write_html("regularized_ILT_results.html")
    fig.show()
    print("Clean parametric results saved to regularized_ILT_results.html")
