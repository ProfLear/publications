import numpy as np
import plotly.graph_objects as go

# System constants
R_np = 1.2       # Nanoparticle radius (nm)
k = 2.45         # Core wavevector (nm^-1)
alpha = 5.1      # Decay constant (nm^-1)
A = 1.30         # Inside normalization constant (nm^-1/2)
C = 125       # Outside normalization constant (nm^-1/2)

# Generate radial grid points from 0 to 10 nm
# Avoid absolute 0 to circumvent divide-by-zero, start very close
r_inside = np.linspace(1e-5, R_np, 300)
r_outside = np.linspace(R_np, 10.0, 1000)

# Compute radial wavefunctions R(r)
R_in = A * np.sin(k * r_inside) / r_inside
R_out = C * np.exp(-alpha * r_outside) / r_outside

# Combine grids for a clean continuous visualization
r_total = np.concatenate([r_inside, r_outside])
R_total = np.concatenate([R_in, R_out])

# Build the Plotly figure
fig = go.Figure()

# Plot the continuous radial wavefunction
fig.add_trace(go.Scatter(
    x=r_total, y=R_total,
    mode='lines',
    line=dict(color='#1f77b4', width=3),
    name='R(r) Wavefunction'
))

# Highlight the boundary interface (Pd surface at 1.2 nm)


# Figure Styling
fig.update_layout(
    title=dict(
        text="<b>Radial Wavefunction R(r) of an Electron on a Pd Nanoparticle</b>",
        x=0.5, font=dict(size=16)
    ),
    xaxis=dict(
        title="Radial Distance, r (nm)", 
        range=[0, 4.0],  # Truncated to 4nm to show detail since it's 0 beyond that
        gridcolor='rgba(200, 200, 200, 0.2)'
    ),
    yaxis=dict(
        title="Radial Wavefunction Amplitude, R(r) (nm<sup>-1/2</sup>)",
        gridcolor='rgba(200, 200, 200, 0.2)'
    ),
    template="plotly_white",
    width=800,
    height=450,
    showlegend=False
)

# Display plot or use in a notebook environment
fig.show()
