import marimo

__generated_with = "0.17.0"
app = marimo.App(width="medium", layout_file="layouts/exampleWorkup.grid.json")


@app.cell
def _(
    form_checkbox,
    form_figure,
    form_pixel_size_number,
    mo,
    sampling_frequency_number,
    wave_checkbox,
    wave_pixel_size_number,
    z_mag_number,
):
    mo.md(
        f"""
    {sampling_frequency_number} {z_mag_number} {form_pixel_size_number} {form_checkbox} {wave_pixel_size_number} {wave_checkbox}
    {mo.as_html(form_figure)}
    """
    )
    return


@app.cell
def _(import_and_downsample, sampling_frequency_number):
    # import data
    xyz_file = "./MiniPlateSmallPre.txt"

    raw_zs, xy_scales = import_and_downsample(xyz_file, sampling = sampling_frequency_number.value)
    return raw_zs, xy_scales


@app.cell
def _(
    form_pixel_size_number,
    gaussian_filter,
    plot_single_areal,
    raw_zs,
    xy_scales,
    z_mag_number,
):
    # now process, first checking the form
    raw_figure = plot_single_areal(raw_zs, xy_scales = xy_scales, z_mag = z_mag_number.value)

    form = gaussian_filter(raw_zs, sigma=form_pixel_size_number.value)  #50 for small, 90 for large
    form_figure = plot_single_areal(form, xy_scales = xy_scales, z_mag = z_mag_number.value)

    raw_zs_minus_form = raw_zs - form
    raw_zs_minus_form_figure = plot_single_areal(raw_zs_minus_form, xy_scales = xy_scales, z_mag = z_mag_number.value)
    return form_figure, raw_zs_minus_form


@app.cell
def _(histogram_of_heights, raw_zs_minus_form):
    histogram_of_heights(raw_zs_minus_form).show()
    return


@app.cell
def _(make_profile_image, raw_zs_minus_form):
    #show the image that has form yet...
    areal_map = make_profile_image(raw_zs_minus_form, 99.9)
    areal_map.show()
    return


@app.cell
def _(
    form_checkbox,
    gaussian_filter,
    plot_single_areal,
    raw_zs_minus_form,
    wave_pixel_size_number,
    xy_scales,
    z_mag_number,
):
    if form_checkbox.value:
        wave = gaussian_filter(raw_zs_minus_form, sigma=wave_pixel_size_number.value)
        wave_figure = plot_single_areal(wave, xy_scales = xy_scales, z_mag = z_mag_number.value)

        raw_zs_minus_form_minus_wave = raw_zs_minus_form - wave
        raw_zs_minus_form_minus_wave_figure = plot_single_areal(raw_zs_minus_form_minus_wave, xy_scales = xy_scales, z_mag = z_mag_number.value)
    return (raw_zs_minus_form_minus_wave,)


@app.cell
def _(GaussianModel, histogram_of_heights, np, sigma_guess):
    def fit_Gaussians(data, cutoff = 99.9, n_Gaussians = 3, center_guesses = None, reverse_amplitude = False, width_scale = 1):
        hist = histogram_of_heights(data, cutoff = cutoff)
        x_data = hist.data[0].x
        y_data = hist.data[0].y
        model = GaussianModel(prefix = "g1_")
        if n_Gaussians > 1:
            for n in range(2, n_Gaussians+1):
                model = model + GaussianModel(prefix = f"g{n}_")

        #now we have the model built out...
        # so add parameters...
        params = model.make_params()

        if center_guesses == None: # then make guesses
            center_guesses = []
            for n in range(1, n_Gaussians+1):
                center_guesses.append(np.mean(data) + n/n_Gaussians*3*sigma_guess)
        else:
            if len(center_guesses) != n_Gaussians:
                print("you sepcified the wrong number of center guesses")


        amplitude_guesses = [] 
        sigma_guesses = [] 
        for n in range(1, n_Gaussians+1):
            sigma_guess = np.std(data)/2
            amplitude_guesses.append(np.sum(y_data)*sigma_guess/n_Gaussians*n)
            sigma_guesses.append(sigma_guess)

        if reverse_amplitude:
            amplitude_guesses.reverse()


        for n, ag, cg, sg in zip(range(1, n_Gaussians+1), amplitude_guesses, center_guesses, sigma_guesses):
            print(ag, cg, sg)
            params.add_many(
                (f"g{n}_amplitude", ag*width_scale, True, ag/100, None),
                (f"g{n}_center", cg, True, cg - 0.1*(max(center_guesses) - min(center_guesses)), cg + 0.1*(max(center_guesses) - min(center_guesses))),
                (f"g{n}_sigma", sg*width_scale, True, sg/100, 3*np.std(data)),
                )

        result = model.fit(y_data, x = x_data, params = params)

        return result
    return (fit_Gaussians,)


@app.cell
def _(fit_Gaussians, qp, raw_zs_minus_form):
    qp.plotFit(fit_Gaussians(raw_zs_minus_form, cutoff= 99.9, n_Gaussians=4, center_guesses=[-20e-6, -10e-6, 0.3e-6, 7e-6], reverse_amplitude=False, width_scale=0.1), components = True, residual = True, output = None).show()
    return


@app.cell
def _(histogram_of_heights, make_profile_image, raw_zs_minus_form_minus_wave):
    # show the image that is only roughness...

    make_profile_image(raw_zs_minus_form_minus_wave, 98).show()
    histogram_of_heights(raw_zs_minus_form_minus_wave, 98).show()
    return


@app.cell
def _(
    fftconvolve,
    find_peaks_simple,
    form_checkbox,
    go,
    np,
    raw_zs_minus_form,
    wave_checkbox,
    xy_scales,
):
    def analyze_pitch_from_acf(z_data, pixel_size_um=1.0):
        # 1. Calculate Normalized ACF
        data_centered = z_data - np.nanmean(z_data) # make sure the mean is taken out
        data_centered = np.nan_to_num(data_centered) # make sure all nan values gone

        # let's generate a core array, which is the inner 25% of values
        center_data = []
        for i, row in enumerate(data_centered):
            if i > 0.25*len(data_centered) and i < 0.75*len(data_centered):
                temp_row = []
                for j, z in enumerate(row):
                    if j > 0.25*len(row) and j < 0.75*len(row):
                        temp_row.append(z)
                center_data.append(np.array(temp_row))


        #acf_2d = fftconvolve(data_centered, data_centered[::-1, ::-1], mode='same')
        acf_2d = fftconvolve(np.array(center_data), data_centered, mode='same')

        acf_2d /= np.max(abs(acf_2d)) # Normalize 0 to 1

        # 2. Find Peaks
        # min_distance: approximate radius in pixels to ignore noise. 
        # Adjust this if your pitch is very small (< 20 pixels).
        peaks = find_peaks_simple(acf_2d, min_distance=10, threshold_rel=0.3)

        # 3. Analyze Peaks to find Pitch
        center = np.array(acf_2d.shape) // 2

        # Calculate distances from center to all peaks
        distances = []
        valid_peaks = []

        for p in peaks:
            dist = np.linalg.norm(p - center)
            if dist > 0: # Ignore the center peak itself (dist=0)
                distances.append(dist)
                valid_peaks.append(p)

        distances = np.array(distances)
        valid_peaks = np.array(valid_peaks)

        # The Pitch is the minimal non-zero distance (Nearest Neighbor)
        if len(distances) > 0:
            min_idx = np.argmin(distances)
            pitch_pixels = distances[min_idx]
            pitch_um = pitch_pixels * pixel_size_um
            nearest_peak = valid_peaks[min_idx]
        else:
            pitch_pixels = 0
            pitch_um = 0
            nearest_peak = center

        return acf_2d, peaks, pitch_pixels, pitch_um, nearest_peak

    # now for some processing of autocorrelation
    if form_checkbox.value and wave_checkbox.value:
        # --- Usage with Plotly ---

        # Assume 'corrected_data' is your profilometry array
        pixel_size = xy_scales[0] # microns
        acf, peaks, pitch_px, pitch_um, p_near = analyze_pitch_from_acf(raw_zs_minus_form, pixel_size_um=pixel_size)

        print(f"Calculated Pitch: {pitch_px:.2f} pixels ({pitch_um:.2f} um)")

        # Create Plotly Figure
        fig = go.Figure()

        # 1. The Heatmap
        fig.add_trace(go.Heatmap(
            z=acf,
            colorscale='RdBu',
            zmid=0,
            name='ACF'
        ))

        # 2. The Detected Peaks (Red Dots)
        # Note: scatter x is column index, y is row index
        '''
        fig.add_trace(go.Scatter(
            x=peaks[:, 1], 
            y=peaks[:, 0],
            mode='markers',
            marker=dict(color='yellow', size=8, symbol='circle-open', line=dict(width=2)),
            name='Local Maxima'
        ))
        '''

        # 3. The Pitch Vector (Line from Center to Nearest Neighbor)
        '''
        center_y, center_x = acf.shape[0]//2, acf.shape[1]//2
        fig.add_trace(go.Scatter(
            x=[center_x, p_near[1]],
            y=[center_y, p_near[0]],
            mode='lines+markers',
            line=dict(color='green', width=3),
            name='Pitch Vector'
        ))
        '''

        fig.update_layout(
            title=f"Pitch Measurement: {pitch_um*1e6:.2f} microns",
            height=600, width=600,
            yaxis=dict(scaleanchor="x", scaleratio=1) # Keep pixels square
        )
    return (fig,)


@app.cell
def _(fig, mo):
    mo.md(f"""{mo.as_html(fig)}""")
    return


@app.cell
def _(go, make_subplots, maximum_filter, np, qp):
    def histogram_of_heights(data, cutoff=99.9):

        # Apply mask
        valid_values = remove_outliers(data, cutoff = cutoff)

        return qp.quickHist(valid_values.flatten(), output = None)

    def remove_outliers(data, cutoff=100): # default is nothing removed...
        # Create a mask where values are within the range

        p_lower = np.percentile(data, 100-cutoff)
        p_upper = np.percentile(data, cutoff)

        print(f"Lower limit: {p_lower}")
        print(f"Upper limit: {p_upper}")

        mask = (data >= p_lower) & (data <= p_upper)

        # Apply mask
        return data[mask]


    def make_profile_image(zdata, cutoff = 100):
        profile_image = make_subplots()
        profile_image.add_heatmap(z = zdata, colorscale="greys_r", zmin = np.percentile(zdata, 100-cutoff), zmax= np.percentile(zdata, cutoff))
        profile_image.update_xaxes(ticks = "", showticklabels = False)
        profile_image.update_yaxes(ticks = "", showticklabels = False)
        profile_image.update_layout(width = 500, height = 500)
        return profile_image

    def find_peaks_simple(image, min_distance=10, threshold_rel=0.5):
        """
        Finds local maxima without skimage.
        1. Dilates the image (finds max in neighborhood).
        2. Checks where original image equals dilated image.
        """
        # 1. Define neighborhood window size
        size = 2 * min_distance + 1

        # 2. Find local max in that window
        image_max = maximum_filter(image, size=size, mode='constant')

        # 3. Boolean mask where Image == Local Max
        mask = (image == image_max)

        # 4. Filter by threshold (ignore noise peaks)
        threshold = np.max(image) * threshold_rel
        mask = mask & (image > threshold)

        # 5. Get coordinates
        coords = np.argwhere(mask)
        return coords





    def plot_single_areal(areal_zs,  xy_scales = [1,1], z_mag = 1):
        fig = make_subplots(specs=[[{'type': 'surface',}]])
        fig.add_trace(go.Surface(
            z=areal_zs - np.mean(areal_zs),
            x=np.linspace(0, len(areal_zs), len(areal_zs))*xy_scales[0],
            y=np.linspace(0, len(areal_zs[0]), len(areal_zs[0]))*xy_scales[1],
            colorscale = "greys_r",
        ),
            )
        fig.update_layout(autosize=True, 
                          margin=dict(l=65, r=50, b=65, t=90),
                          scene=dict(zaxis=dict(range = [-0.5*len(areal_zs)*xy_scales[1]/z_mag, 0.5*len(areal_zs)*xy_scales[1]/z_mag])),
                          )
        return fig


    def import_and_downsample(xyz_file, sampling = 1):
        with open (xyz_file, "r") as xyz:
            # this block gets the information on the image width and height
            n_headers = 0
            n_data_rows = 0
            n_data_cols = 0
            zs = []
            for i, line in enumerate(xyz, start = 1):
                line = line.strip("\n")  # let us remove newline characters
                if line[0] == "#": # we are in a header
                    n_headers = n_headers + 1
                    if "Width" in line:
                        total_width = float(line.split(" ")[-2])
                        #print(line.split(" "))
                        if line.split(" ")[-1] == "mm":
                            total_width = total_width * 1E-3 # convert from mm to meters
                        else:
                            print("width unit is not recognized")
                    if "Height" in line:
                        total_height = float(line.split(" ")[-2])
                        if line.split(" ")[-1] == "mm":
                            total_height = total_height * 1E-3 # convert from mm to meters
                        else:
                            print("height unit is not recognized")
                    if "Value units" in line:
                        if line.split(" ")[-1] == "m":
                            z_scale = 1 # what we should multiply by to get to meters

                else: # then we are not in a header. This is data. 
                    n_data_rows = n_data_rows + 1
                    if (i - n_headers)%sampling == 0:# this is a row we will sample
                        temp_row = []
                        for j, value in enumerate(line.split("\t")):
                            if j%sampling == 0: # this is a point to sample
                                temp_row.append(float(value)*z_scale)
                                final_sampled_col = j
                        if n_data_cols == 0: # we haven't yet recorded the number of data points
                            n_data_cols = j # so recordit
                        zs.append(np.array(temp_row))
                        final_sampled_row = i

        zs = np.array(zs)
        return zs, [(total_width*final_sampled_col/n_data_cols)/len(zs[0]), (total_height*final_sampled_row/n_data_rows)/len(zs)] # also return x_scale and y_scale in xy_scales
    return (
        find_peaks_simple,
        histogram_of_heights,
        import_and_downsample,
        make_profile_image,
        plot_single_areal,
    )


@app.cell
def _():
    import marimo as mo
    import numpy as np
    from plotly.subplots import make_subplots 
    import plotly.graph_objects as go
    from scipy.ndimage import gaussian_filter
    from scipy.signal import fftconvolve
    from scipy.ndimage import maximum_filter
    import codechembook.quickPlots as qp
    from lmfit import Model
    from lmfit.models import GaussianModel
    import altair as alt
    import pandas as pd


    sampling_frequency_number = mo.ui.number(start = 0, value = 10, label="sampling frequency") # how often to sample a pixel 
    z_mag_number = mo.ui.number(start = 1, value = 1, label = "z-magnification")
    form_pixel_size_number = mo.ui.number(start = 1, value = 50, label = "form size (in pixels)")
    wave_pixel_size_number = mo.ui.number(start = 1, value = 3, label = "wave size (in pixels)")
    form_checkbox = mo.ui.checkbox(value = False, label = "Form accetable?")
    wave_checkbox = mo.ui.checkbox(value = True, label = "Wave acceptable?")
    return (
        GaussianModel,
        fftconvolve,
        form_checkbox,
        form_pixel_size_number,
        gaussian_filter,
        go,
        make_subplots,
        maximum_filter,
        mo,
        np,
        qp,
        sampling_frequency_number,
        wave_checkbox,
        wave_pixel_size_number,
        z_mag_number,
    )


if __name__ == "__main__":
    app.run()
