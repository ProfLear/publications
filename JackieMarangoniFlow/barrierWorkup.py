import marimo

__generated_with = "0.23.6"
app = marimo.App(width="medium", layout_file="layouts/barrierWorkup.grid.json")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import importlib
    import supportingFunctions as sf
    from scipy.ndimage import gaussian_filter
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    importlib.reload(sf)
    from scipy.stats import binned_statistic


    return binned_statistic, go, importlib, make_subplots, mo, np, sf


@app.cell
def _(mo):
    # define the file that has the profile data that we wish to analyze
    file_input = mo.ui.text("BarriersBLBottom", label = "filepath")

    #profilometry_file  = "/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/JackieMarangoniFlow/10101TL.txt"
    return (file_input,)


@app.cell
def _(file_input):
    file_input
    return


@app.cell
def filestring(file_input):
    profilometry_file = "/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/JackieMarangoniFlow/" + file_input.value + ".txt"
    print(profilometry_file)
    return (profilometry_file,)


@app.cell
def _(profilometry_file, sf):
    # read in this profile data, being sure to downsample, so that the data is not too large to display.  The downsampling is required in order to display things reasonably. But we can always use the raw for analysis, if we want to preserve data...
    # raw_zs, raw_scales = sf.import_and_downsample(profilometry_file, sampling = 1)
    sampling = 20
    downsampled_zs, downsampled_scales = sf.import_and_downsample(profilometry_file, sampling = sampling)
    return downsampled_scales, downsampled_zs, sampling


@app.cell
def _(downsampled_scales, downsampled_zs, sampling, sf):
    # let us just see what the raw data looks like (downsampled) so that we can get a sense of the form and wave and roughness present. 
    raw_figure = sf.plot_single_areal(downsampled_zs, xy_scales = downsampled_scales*sampling, z_mag = 10)
    raw_figure.show()
    return


@app.cell
def _(
    Path,
    folder,
    np,
    process_zygos_xyz,
    raw_form_wave_rough,
    xyz_average_out_nanvalues,
):
    # export for blender
    def xyz_to_verts_faces (xyz, x_scale = 1, y_scale = 1):
        '''


        Parameters
        ----------
        xyz : 2D array
            the z-values of the surface.
        x_scale : float, optional
            Used to convert from array indices into real space coordinates/distances. The default is 1.
        y_scale : float, optional
            Used to convert from array indices into real space coordinates/distances. The default is 1.

        Returns
        -------
        verts : list of lists
            each sublist is the x, y, and z coordinates.
        faces : list of list
            each sublist is a list of 4 values, which are the indices of the corners of a face. .

        '''
        verts = []
        for i, x in enumerate(xyz): #the row is the x-index
            for j, y in enumerate(x): # walk through the y-indices. 
                verts.append(
                    (
                        (i + 0) * x_scale, 
                        (j + 0) * y_scale, 
                        xyz[i][j])
                    )
        n_x = len(xyz)
        n_y = len(xyz[0])
        faces = []
        for x in range(1, n_x): # this is because we start counting vertices at 1
            for y in range(1, n_y-1): # start at 1 and go second to the end
                faces.append(
                    ((x-1)*n_y + y, # Top left
                     (x-1)*n_y + y + 1, # top right. 
                     x*n_y + y + 1, # bottom right
                     x*n_y + y)   # bottom left
                    )
        #     
        return verts, faces



    def prep_for_blender (file = None, smoothing = None):
        '''
        Writes files that can be read by another script for importing surfaces into blender. 
        This file takes raw data, and then gets form, wave, and rough splines.  Then write 6 files from them.

        Parameters
        ----------
        file : path object, optional
            the file that we want to extract form, wave, and rough and then prepare files for blender. The default is None.
        smoothing : dict, optional
            contains values for form, wave, and rough splines. The default is None.

        Returns
        -------
        None.

        '''
        verts, faces, xyz_array, xy_scale = process_zygos_xyz(folder/file)

        # use a loop to remove NAN values, while the exist
        if np.any(np.isnan(xyz_array)) == False: #there are no nan values to hand
            print("no nan values to handle!")
        while np.any(np.isnan(xyz_array)) == True:
            xyz_array = xyz_average_out_nanvalues(xyz_array)

        #% Set up a plot, so we can determine the correct values of s to use....
        splines, data = raw_form_wave_rough(xyz_array, s_large = smoothing["form"], s_wave = smoothing["wave"], s_rough = smoothing["rough"])

        for name, surface in zip(['1form', '2wave', '3rough', '1raw', '2middle', '3final' ], [splines[0], splines[1], splines[2], data[0], data[1], data[2]]):
            verts, faces = xyz_to_verts_faces(surface)
            write_verts_faces_to_file(file.parent / Path(f"{file.stem}_{name}.csv"), verts = verts, faces = faces)
        print("files written!")


    def write_verts_faces_to_file (filepath, verts = [], faces = []):
        print(f"writing to {filepath}")
        with open(filepath, 'w') as f:
            for coords in [verts, faces]: # go through verts first, then faces
                for row in coords:
                    for element in row:
                        f.write(f'{element}, ')
                    f.write("\n")

    def prep_xyz_for_blender (xyz, filepath):
        verts, faces = xyz_to_verts_faces(xyz)
        write_verts_faces_to_file(filepath, verts = verts, faces = faces)
        print("files written!")

    #prep_xyz_for_blender(downsampled_raw_zs_minus_form, r"")
    return


@app.cell
def _(downsampled_zs, sf):
    custom_scale = [
            [0, 'rgb(3, 52, 78)'],     # Start (Min)
            [0.5, 'rgb(7, 115, 177)'], # Middle
            [1, 'rgb(39, 174, 246)']      # End (Max)
        ]

    downsampled_wave_map = sf.make_profile_image(downsampled_zs, 90, colorscale = custom_scale)
    downsampled_wave_map.update_layout(template = "simple_white", paper_bgcolor = 'rgba(0,0,0,0)', )
    downsampled_wave_map.show()
    #downsampled_wave_map.write_image(r"./test.png")
    #print(f"rms roughness is: {np.sqrt(np.mean((downsampled_raw_zs_minus_form - np.mean(downsampled_raw_zs_minus_form))**2))}")

    len(downsampled_zs[0])
    return (custom_scale,)


@app.cell
def _(importlib, sf):
    importlib.reload(sf)
    return


@app.cell
def _(binned_statistic, np):
    def extract_rect_profile(data, start_pt, end_pt, width):
        """
        Extracts z-heights within a rectangle defined by a line and a width.

        Parameters:
        - data: 2D numpy array (the profilometry map)
        - start_pt: (y1, x1) tuple
        - end_pt: (y2, x2) tuple
        - width: total width of the rectangle (distance perpendicular to line) in units of pixels.
        """
        x1, y1 = start_pt
        x2, y2 = end_pt

        # 1. Create a grid of coordinates
        rows, cols = data.shape
        y, x = np.indices((rows, cols))

        # 2. Vector math
        # Vector from A to B (the line)
        line_vec = np.array([y2 - y1, x2 - x1])
        line_len_sq = np.sum(line_vec**2)
        line_len = np.sqrt(line_len_sq)

        # Vector from A to all points P
        ap_vec_y = y - y1
        ap_vec_x = x - x1

        # 3. Scalar Projection: How far along the line is P?
        # t is the normalized distance (0 to 1) along the segment
        t = (ap_vec_y * line_vec[0] + ap_vec_x * line_vec[1]) / line_len_sq

        # 4. Find distance from the line (Perpendicular)
        # Closest point on the infinite line to P
        closest_y = y1 + t * line_vec[0]
        closest_x = x1 + t * line_vec[1]

        dist_from_line = np.sqrt((y - closest_y)**2 + (x - closest_x)**2)

        # 5. Filter points within the rectangle
        # Must be within [0, 1] of the segment and within width/2 of the line
        mask = (t >= 0) & (t <= 1) & (dist_from_line <= width / 2)

        # 6. Extract values
        z_values = data[mask]
        distances_along_line = t[mask] * line_len

        return distances_along_line, z_values

    def bin_profile_data(distances, z_values, bin_size):
        """
        Bins the point cloud data into a single profile line.

        Parameters:
        - distances: 1D array of distances along the line (from extract_rect_profile)
        - z_values: 1D array of z-heights
        - bin_size: The width of each bin (in the same units as distances)

        Returns:
        - bin_centers: The x-coordinates for the final profile
        - binned_z: The average z-height for each bin
        """
        # Define the bin edges from 0 to the maximum distance
        max_dist = np.max(distances)
        bin_edges = np.arange(0, max_dist + bin_size, bin_size)

        # Calculate the mean z-value for each bin
        # 'statistic' can be 'mean', 'median', 'std', etc.
        binned_z, edges, _ = binned_statistic(distances, z_values, 
                                              statistic='mean', 
                                              bins=bin_edges)

        # Calculate centers of bins for plotting
        bin_centers = (edges[:-1] + edges[1:]) / 2

        return bin_centers, binned_z


    def bin_profile_with_se(distances, z_values, bin_size):
        """
        Bins the data and calculates mean and standard error.
        """
        max_dist = np.max(distances)
        bin_edges = np.arange(0, max_dist + bin_size, bin_size)

        # 1. Calculate Mean
        mu, edges, _ = binned_statistic(distances, z_values, 'mean', bins=bin_edges)

        # 2. Calculate Standard Deviation (sigma)
        sigma, _, _ = binned_statistic(distances, z_values, 'std', bins=bin_edges)

        # 3. Calculate Count (n)
        n, _, _ = binned_statistic(distances, z_values, 'count', bins=bin_edges)

        # 4. Calculate Standard Error
        # Use np.sqrt(n) and handle potential division by zero for empty bins
        se = sigma / np.sqrt(n)

        bin_centers = (edges[:-1] + edges[1:]) / 2

        return bin_centers, mu, se



    return bin_profile_with_se, extract_rect_profile


@app.cell
def _(binned_heights, mo, np):
    mo.md(f"""
    Largest height difference = {(np.max(binned_heights) - np.min(binned_heights))*1E6:.1f}microns
    """)
    return


@app.cell
def _(mo):

    xs_input = mo.ui.text("", label = "x positions")
    ys_input = mo.ui.text("", label = "y positions")
    cutoff_input = mo.ui.number(90, label = "cutoff value")
    return cutoff_input, xs_input, ys_input


@app.cell
def _(xs_input, ys_input):
    try: 
        x_list = list(map(int, xs_input.value.split(",")))
    except:
        x_list = []

    try:
        y_list = list(map(int, ys_input.value.split(",")))
    except:
        y_list = []
    return x_list, y_list


@app.cell
def _(x_list, y_list):
    print(len(x_list))
    print(len(y_list))

    if len(x_list) == len(y_list):
        print("ok")
        for _i in range(0, len(x_list), 2):
            print(_i)
    return


@app.cell
def _(
    bin_profile_with_se,
    custom_scale,
    cutoff_input,
    downsampled_scales,
    downsampled_zs,
    extract_rect_profile,
    file_input,
    go,
    make_subplots,
    mo,
    np,
    x_list,
    y_list,
):
    def get_rectangle_path(start_pt, end_pt, width):
        x1, y1 = start_pt
        x2, y2 = end_pt

        # 1. Vector math for corners
        dx = x2 - x1
        dy = y2 - y1
        length = np.sqrt(dx**2 + dy**2)

        # Unit perpendicular vector
        ux_perp = -dy / length
        uy_perp = dx / length

        # Calculate 4 corners
        c1 = (x1 + (width/2)*ux_perp, y1 + (width/2)*uy_perp)
        c2 = (x1 - (width/2)*ux_perp, y1 - (width/2)*uy_perp)
        c3 = (x2 - (width/2)*ux_perp, y2 - (width/2)*uy_perp)
        c4 = (x2 + (width/2)*ux_perp, y2 + (width/2)*uy_perp)

        # make path
        path = f"M {c1[0]},{c1[1]} L {c2[0]},{c2[1]} L {c3[0]},{c3[1]} L {c4[0]},{c4[1]} Z"

        #print(path)
        return path


    cutoff_input.value

    map_with_profiles = make_subplots(
        rows = int(len(x_list)/2 + 1), cols = 1,
        row_heights= [600] + [100]*int(len(x_list)/2)
                                    )
    map_with_profiles.add_trace(
        go.Heatmap(
            z=downsampled_zs, 
            colorscale=custom_scale,
            showscale = False,
            zmin=np.percentile(downsampled_zs, 100-cutoff_input.value), # this helps with coloration, setting a floor value
            zmax=np.percentile(downsampled_zs, cutoff_input.value), # this helps wiht coloration, setting a ceiling value
            # If you have x and y coordinates, add them here:
            # x=your_x_data,
            # y=your_y_data
        ),
        row=1, col=1
    )
    map_with_profiles.update_xaxes(range = [0, len(downsampled_zs[0])], title = "pixels", row = 1, col = 1)
    map_with_profiles.update_yaxes(range = [0, len(downsampled_zs)], title = "pixels", row = 1, col = 1)

    lineplot = make_subplots(rows = int(len(x_list)/2 + 1), cols = 1)


    #map_with_profile = go.Figure(downsampled_wave_map)

    shapes_list = []
    maxlist = []
    if len(x_list) == len(y_list) and len(x_list) > 0:
        for _i in range(0,len(x_list), 2):
            #print(_i)
            start_x = x_list[_i]
            end_x = x_list[_i+1]
            x_range = end_x - start_x

            start_y = y_list[_i]
            end_y = y_list[_i+1]
            y_range = end_y - start_y

            map_with_profiles.add_annotation(x = start_x + 0.1*x_range, y = start_y + 0.1*y_range, text = f"{int(_i/2 + 1)}i", showarrow=False, font = dict(color = "white"))
            map_with_profiles.add_annotation(x = end_x - 0.1*x_range, y = end_y - 0.1*y_range, text = f"{int(_i/2 + 1)}f", showarrow=False, font = dict(color = "white"))

            shapes_list.append(
            dict(
                type="path",
                path=get_rectangle_path([start_x, start_y], [end_x, end_y], width = 3),
                line=dict(color="white", width=2, dash="solid"),
                fillcolor="rgba(255, 255, 255, 0.66)"
            ),
            # Keep the aspect ratio square so the rectangle doesn't look distorted
        )

            distances, heights = extract_rect_profile(downsampled_zs, [start_x, start_y], [end_x, end_y], width = 20)
            binned_distances, binned_heights, binned_se = bin_profile_with_se(distances, heights,  3)

            # convert units from m to microns
            binned_distances = binned_distances*1E6
            binned_heights = binned_heights * 1e6


            y_offset = np.min(binned_heights-binned_se)
            maxlist.append(np.max(binned_heights-y_offset))

            #lineplot.add_scatter(x = distances, y = heights, mode = "markers")
            #lineplot.add_scatter(x = binned_distances, y = binned_heights+binned_se - y_offset, line = dict(color = 'rgb(39, 174, 246)', width = 0), fill = "tozeroy", fillcolor = 'rgb(39, 174, 246)', showlegend = True, name = f"segment {int(_i/2+1)}")
            #lineplot.add_scatter(x = binned_distances, y = binned_heights-binned_se - y_offset, line = dict(color = 'rgb(39, 174, 246)', width = 0), fill = "tozeroy", fillcolor="white", showlegend = False)
            map_with_profiles.add_scatter(x = binned_distances*downsampled_scales[0], y = binned_heights-y_offset, line = dict(color = 'rgb(7, 115, 177)'), showlegend = False, row = int(_i/2)+2, col = 1)
            map_with_profiles.add_annotation(x = 0, y = 0, text = f"largest height difference = {(np.max(binned_heights) - np.min(binned_heights)):0.1f} microns; stdev = {np.std(binned_heights):.1f} microns", xanchor="left", yanchor = "bottom", showarrow = False, row = int(_i/2)+2, col = 1)
            map_with_profiles.update_yaxes(title = "microns", range = [0, max(maxlist)], row = int(_i/2)+2, col = 1)


    else:
        print("error, x and y lists are not the same length")


    map_with_profiles.update_layout(shapes=shapes_list) # add empty shapes list that we will append

    if len(x_list) > 0:
        map_with_profiles.update_xaxes(title = "microns", row = int(_i/2 + 2), col = 1)


    map_with_profiles.update_layout(
        template = "simple_white",
        width = 600, 
        height = 600 + 100* int(len(x_list)/2),
        yaxis = dict(scaleanchor="x",scaleratio=1),
        title = f"</br>{file_input.value}, </br>x-scale = {downsampled_scales[0]*1e6} microns, </br>y-scale = {downsampled_scales[1]*1e6} microns"
    )



    lineplot.update_layout(template = "simple_white")



    mo.ui.plotly(map_with_profiles)
    return binned_heights, map_with_profiles


@app.cell
def _(cutoff_input, xs_input, ys_input):
    xs_input, ys_input, cutoff_input
    return


@app.cell
def _(mo):
    write_button = mo.ui.run_button(label="Write html file of plots.")
    write_button
    return (write_button,)


@app.cell
def _(file_input, map_with_profiles, mo, write_button):
    mo.stop(not write_button.value)
    map_with_profiles.write_html(f"{file_input.value}.html")
    return


@app.cell
def _():
    #print(map_with_profile.layout.shapes)
    return


@app.cell
def _(downsampled_scales):
    downsampled_scales
    return


@app.cell
def _(write_button):
    write_button.value
    return


if __name__ == "__main__":
    app.run()
