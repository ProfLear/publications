import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import os
    from pathlib import Path

    PC_directory = Path(r"C:\Users\benle\Documents\GitHub\Photothermal\zygoprocessing\dataClassApproach")
    os.chdir(PC_directory)
    import profilometryClasses as pc

    data_directory = Path(r"C:\Users\benle\Documents\GitHub\publications\Sandra Textiles\profilometry")
    os.chdir(data_directory)

    return


if __name__ == "__main__":
    app.run()
