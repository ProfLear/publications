import marimo

__generated_with = "0.19.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import importlib
    import supportingFunctions as sf
    from scipy.ndimage import gaussian_filter
    import plotly.graph_objects as go
    return gaussian_filter, go, importlib, np, sf


@app.cell
def _():
    # define the file that has the profile data that we wish to analyze
    profilometry_file  = "/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/JackieSharkSkin/profilometry/BlankPostBarnacle.txt"
    return (profilometry_file,)


@app.cell
def _(profilometry_file, sf):
    # read in this profile data, being sure to downsample, so that the data is not too large to display.  The downsampling is required in order to display things reasonably. But we can always use the raw for analysis, if we want to preserve data...
    # raw_zs, raw_scales = sf.import_and_downsample(profilometry_file, sampling = 1)
    downsampled_zs, downsampled_scales = sf.import_and_downsample(profilometry_file, sampling = 15)
    return downsampled_scales, downsampled_zs


@app.cell
def _(downsampled_scales, downsampled_zs, gaussian_filter, sf):
    # handle the form correction...
    downsampled_form = gaussian_filter(downsampled_zs, sigma=30)  #50 the form is on a size of 50 pixels in the downsampled data.
    form_figure = sf.plot_single_areal(downsampled_form, xy_scales = downsampled_scales, z_mag = "auto")
    form_figure.show()

    downsampled_raw_zs_minus_form = downsampled_zs - downsampled_form
    downsampled_raw_zs_minus_form_figure = sf.plot_single_areal(downsampled_raw_zs_minus_form, xy_scales = downsampled_scales, z_mag = "auto")

    downsampled_raw_zs_minus_form_figure.show()
    return (downsampled_raw_zs_minus_form,)


@app.cell
def _(downsampled_scales, downsampled_zs, sf):
    # let us just see what the raw data looks like (downsampled) so that we can get a sense of the form and wave and roughness present. 
    raw_figure = sf.plot_single_areal(downsampled_zs, xy_scales = downsampled_scales, z_mag = "auto")
    raw_figure.show()
    return


@app.cell
def _(
    Path,
    downsampled_raw_zs_minus_form,
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

    prep_xyz_for_blender(downsampled_raw_zs_minus_form, r"/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/JackieSharkSkin/profilometry/flatBarnacleBlender.csv")
    return


@app.cell
def _(downsampled_raw_zs_minus_form, np, sf):
    custom_scale = [
            [0, 'rgb(3, 52, 78)'],     # Start (Min)
            [0.5, 'rgb(7, 115, 177)'], # Middle
            [1, 'rgb(39, 174, 246)']      # End (Max)
        ]

    downsampled_wave_map = sf.make_profile_image(downsampled_raw_zs_minus_form, 98, colorscale = custom_scale)
    downsampled_wave_map.update_layout(template = "simple_white", paper_bgcolor = 'rgba(0,0,0,0)')
    downsampled_wave_map.show()
    downsampled_wave_map.write_image(r"/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publicationsData/JackieSharkSkin/profilometry/flatBarnacle.png")
    print(f"rms roughness is: {np.sqrt(np.mean((downsampled_raw_zs_minus_form - np.mean(downsampled_raw_zs_minus_form))**2))}")
    return (downsampled_wave_map,)


@app.cell
def _(downsampled_wave_map, go, sf):
    # from the above, we can also pick points to compare different in heights
    # these are done with pairs of points, x, y, z

    paired_heights = [
        [[237, 559, -90.15247e-6], [162, 589, 49.7457e-6]],
        [[442, 555, -78.74439e-6], [375, 588, 56.19653e-6]],
        [[437, 210, -70.44566e-6], [372, 232, 54.87099e-6]],
        [[542, 35, -92.7806e-6], [463, 56, 43.78794e-6]],
        [[135, 28, -69.0605e-6], [62, 51, 43.50909e-6]],
    ]

    mapped_heights = sf.map_paired_heights(go.Figure(downsampled_wave_map), paired_heights, low_color = "rgb(56, 18, 35)", high_color = "rgb(255, 128, 184)", line_color = "rgb(122, 122, 122)")
    mapped_heights.show()
    mapped_heights.write_image(r"C:\Users\benle\OneDrive - The Pennsylvania State University\publications\JackieSharkSkin\forFigures\largePreMapped.png")
    return


@app.cell
def _(importlib, sf):
    importlib.reload(sf)
    return


if __name__ == "__main__":
    app.run()
