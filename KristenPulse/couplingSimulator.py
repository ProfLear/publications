import marimo

__generated_with = "0.11.31"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    from plotly.subplots import make_subplots

    return make_subplots, mo, np


@app.cell
def _(make_subplots, np):
    radius = 3
    k = 0.5
    A = 1
    r = np.linspace(radius, 10, 1000)
    psi = A*np.exp(-k*(r-radius))

    fig = make_subplots(rows = 3, cols = 1)
    fig.add_scatter(x = r, y = psi**2, line = dict(color = "magenta"), showlegend = False, row =1, col =1)

    # just constant density of H atoms
    fig.add_scatter(x = r, y = r**2, line = dict(color = "darkcyan"), showlegend = False, row =2, col = 1)




    # if this is just the ligands, standing straight out, then the number of H-atoms is constant with r, and ends at the end of the length...
    n_H_ligand = [10]*500 + [0]*500
    fig.add_scatter(x = r, y = n_H_ligand, line = dict(color = "red"), showlegend = False, row =2, col = 1)

    from codechembook.quickPlots import quickHist
    coupling_bins = quickHist(psi**2*r**2, nbins = 100, output = None)
    coupling_bins2 =  quickHist(psi**2*n_H_ligand, nbins = 100, output = None)


    fig.add_scatter(x = coupling_bins.data[0].x, y = coupling_bins.data[0].y, line = dict(color = "darkcyan"), showlegend = False, row =3, col =1)
    fig.add_scatter(x = -coupling_bins.data[0].x, y = coupling_bins.data[0].y, line = dict(color = "darkcyan"), showlegend = False, row =3, col =1)

    fig.add_scatter(x = coupling_bins2.data[0].x, y = coupling_bins2.data[0].y, line = dict(color = "red"), showlegend = False, row =3, col =1)
    fig.add_scatter(x = -coupling_bins2.data[0].x, y = coupling_bins2.data[0].y, line = dict(color = "red"), showlegend = False, row =3, col =1)

    fig.update_xaxes(range = [-2, 2], row =3, col = 1)

    fig.show()
    return (
        A,
        coupling_bins,
        coupling_bins2,
        fig,
        k,
        n_H_ligand,
        psi,
        quickHist,
        r,
        radius,
    )


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
