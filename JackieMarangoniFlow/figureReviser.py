import marimo

__generated_with = "0.23.6"
app = marimo.App(width="medium")


@app.cell
def _():
    from plotlyRehydrator import load_plotly_html
    import plotly.graph_objects as go

    return (load_plotly_html,)


@app.cell
def _(load_plotly_html):
    oldFig = load_plotly_html("/Users/benjaminlear/GitHub/publications/JackieMarangoniFlow/9401M.html")
    oldFig.show()
    oldFig.layout
    return (oldFig,)


@app.cell
def _(oldFig):
    oldFig.layout
    return


@app.cell
def _(oldFig, revise_figure):
    revise_figure(oldFig,rotate=0, scale  = 3)
    return


@app.cell
def _():
    def revise_figure(oldfig,
                     rotate = 0, scale = 1,
                      colorway = ["orange", "lightgreen"],
                     ):
        import plotly.graph_objects as go
        import numpy as np
        newfig = go.Figure(oldfig)

        newfig.data[0].z = np.rot90(np.array(oldfig.data[0].z), k=rotate/90) # rotate the data
    
        n_xpixels = len(newfig.data[0].z[0])
        n_ypixels = len(newfig.data[0].z)
    
        newfig.layout.xaxis1.range = [0, n_xpixels]
        newfig.layout.xaxis1.showline=False
        newfig.layout.xaxis1.showticklabels=False
        newfig.layout.xaxis1.ticks="outside"
        newfig.layout.xaxis1.tickwidth=4
        newfig.layout.xaxis1.ticklen=8

        newfig.layout.xaxis1.tickvals = [n_xpixels/2 - int(50/scale), n_xpixels/2 + int(50/scale)]
        newfig.layout.xaxis1.title.text = "100 microns"
        newfig.layout.xaxis1.title.standoff = 10


    
        newfig.layout.yaxis1.range = [0, n_ypixels]
        newfig.layout.yaxis1.showline=False
        newfig.layout.yaxis1.showticklabels=False
        newfig.layout.yaxis1.ticks=""
        newfig.layout.yaxis1.title = ""

        if n_xpixels < n_ypixels: # the image is taller than it is wide
            newfig.layout.yaxis1.scaleanchor="x"
            newfig.layout.yaxis1.scaleratio=1
        else:
            newfig.layout.xaxis1.scaleanchor="y"
            newfig.layout.xaxis1.scaleratio=1



        newpaths = []
        oldpaths = []
        for i, s in enumerate(newfig.layout.shapes):
            oldpaths.append(s.path)
            newpaths.append(rotate_svg_path(s.path, angle_degrees = rotate))

        for i, p in enumerate(newpaths):
            newfig.layout.shapes[i].path = p

        print(oldpaths)  
        print(newpaths)
        # annotations
        newfig.layout.annotations = []

        for i, c in enumerate(colorway):
            newfig.layout.shapes[i].fillcolor = c
            newfig.layout.shapes[i].line.color = c

            newfig.data[i+1].line.color = c

        #print(newfig.layout)
        newfig.show()


    def rotate_svg_path(path_str, angle_degrees, center=(0, 0)):
        """
        Rotates an SVG path string by 90, 180, or 270 degrees around a center point.
        """
        import re
        import numpy as np
        # Convert angle to radians
        theta = np.radians(-1*angle_degrees) # the numpy rotation is clockwise for some reason, so need to rotate this the oposite direction than we think
        c, s = np.cos(theta), np.sin(theta)
    
        # 2D Rotation matrix
        R = np.array([[c, -s], 
                      [s,  c]])
    

        # Helper function to rotate a single coordinate match
        def rotate_match(match):
            print(match)
            # Extract the command (M, L, etc.) and the coordinates
            cmd = match.group(1)
            x = float(match.group(2))
            y = float(match.group(3))
        
            # Shift to origin, rotate, and shift back
            xy = np.array([x - cx, y - cy])
            xy_rotated = R @ xy
            new_x = xy_rotated[0] + cx
            new_y = xy_rotated[1] + cy
        
            # Format back to 4 decimal places to keep it clean
            return f"{cmd} {new_x:.4f} {new_y:.4f}"

        # Regex to find SVG commands followed by coordinate pairs (e.g., 'M 1.5 2.3' or 'L -1 4')
        # This handles standard M (move to), L (line to), and Z (close path)

        pattern = r'(\w*)\s*([-\d.]+),([-\d.]+)'
        c_pattern = re.compile(pattern)
        center_points = c_pattern.match(path_str)
        center = [float(center_points.group(2)), float(center_points.group(3))]
        cx, cy = center
        rotated_path = re.sub(pattern, rotate_match, path_str)
        return rotated_path

    return (revise_figure,)


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
