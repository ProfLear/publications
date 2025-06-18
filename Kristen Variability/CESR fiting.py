import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
    # CESR fitting

    The purpose of this notebook is to explore some of the fitting of the Pd-S-alkane nanoparticles from Kristen
    """
    )
    return


@app.cell
def _():
    # handle all imports here
    import marimo as mo
    import numpy as np
    from plotly.subplots import make_subplots
    from codechembook import quickPlots as qp
    from lmfit import Model
    from pathlib import Path
    from scipy.integrate import cumulative_trapezoid
    from scipy.signal import hilbert
    return (
        Model,
        Path,
        cumulative_trapezoid,
        hilbert,
        make_subplots,
        mo,
        np,
        qp,
    )


@app.cell
def _(make_subplots):
    # for plotting the fit results...

    def plot_fit(fit_result):
        independent_vars = fit_result.model.independent_vars
        plot = make_subplots(rows = 2, cols = 1)
        plot.add_scatter(x = fit_result.userkws[independent_vars[0]], y = fit_result.data, 
                         mode = "lines", 
                         line = dict(color = "darkcyan", width = 5), 
                         showlegend = False
                        )
        plot.add_scatter(x = fit_result.userkws[independent_vars[0]], y = fit_result.best_fit,
                        mode = "lines", 
                         line = dict(color = "black", width = 2), 
                         showlegend = False
                        )

        plot.add_scatter(x = fit_result.userkws[independent_vars[0]], y = fit_result.residual,
                        mode = "lines", 
                         line = dict(color = "black", width = 2), 
                         showlegend = False,
                         row = 2, col = 1,
                        )

        ymax = max([
            abs(max(fit_result.data)), 
            abs(min(fit_result.data)),
            abs(max(fit_result.best_fit)), 
            abs(min(fit_result.best_fit)),
            ])

        plot.update_xaxes(showline = False)
        plot.update_yaxes(range = [-1.05*ymax, 1.05*ymax], zeroline = True)
        plot.update_layout(template = "simple_white", title = f"{fit_result.model.name}: AIC:{int(fit_result.aic)}")
        plot.show()


    return (plot_fit,)


@app.cell
def _(cumulative_trapezoid, np):
    # data handlers

    def double_integral(x, y):
        """
        Computes the double integral of a given spectrum (y vs. x).

        This is done by performing a cumulative trapezoidal integration twice.

        Args:
            x (np.ndarray): The array of x-values (e.g., wavenumber, frequency).
                            Must be the same size as y.
            y (np.ndarray): The array of y-values (e.g., absorbance, intensity).
                            Must be the same size as x.

        Returns:
            tuple[np.ndarray, np.ndarray]: A tuple containing:
                - first_integral (np.ndarray): The result of the first cumulative integration.
                - second_integral (np.ndarray): The result of the second integration.
                Both arrays have the same size as the input x.
        """
        # Ensure inputs are numpy arrays
        x = np.asarray(x)
        y = np.asarray(y)

        if x.shape != y.shape:
            raise ValueError("Input arrays x and y must have the same shape.")

        # 1. Compute the first cumulative integral.
        # We use initial=0 to make the output array the same size as the input.
        # This means the integral starts at 0.
        first_integral = cumulative_trapezoid(y, x, initial=0)

        # 2. Compute the second integral by integrating the first result.
        # Again, use initial=0 to maintain the array size.
        second_integral = cumulative_trapezoid(first_integral, x, initial=0)

        return first_integral, second_integral



    def get_instrument_ID(file):
        with open(file, "r") as f:
            for line in f:
                if line[:2] == "##": # this determines the instrument being used...
                    instrument_ID = "##" # will use this to help parse things

                elif ";" in line:
                    instrument_ID = ";"

                break

        return instrument_ID


    # import EPR
    def import_epr_data(file):
        B_fields, intensities = [], []

        instrument_ID = get_instrument_ID(file)
        print(instrument_ID)
        if instrument_ID ==";":
            with open(file, "r") as f:
                for line in f:
                    if "Frequency;" in line:
                        instrument_frequency = float(line.split(";")[1])*10**9
                    if line[0].isdigit():
                        B, I = line.split(";")
                        B_fields.append(float(B))
                        intensities.append(float(I))

        if instrument_ID == "##":

            frequency_string = str(file.stem).split("_")[-1]
            if frequency_string[-3]=="k":
                instrument_frequency = int(frequency_string[:-3])*1000 # convert to Hz
            with open(file, "r") as f:
                for line in f:
                    #print(line)
                    if line[0].isdigit():
                        B, I = line.split("   ") #separator is 3 spaces...
                        B_fields.append(float(B))
                        intensities.append(float(I))


        return [np.array(B_fields), np.array(intensities)], instrument_frequency


    def get_lineshape_points(xydata):
        # we will get the positive peak, negative peak, and zero crossing...

        index_max = np.argmax(xydata[1])
        index_min = np.argmin(xydata[1])

        # a better solution might be to calculate the derivative for the y-data, and then find the max that is between the min and max arguments
        index_zero = np.argmin(abs(xydata[1]-np.mean(xydata[1])))

        return [index_max, index_zero, index_min]

    #trim EPR data
    def trim_data(x, y, limits):

        index1 = np.argmin(abs(x-limits[0]))
        index2 = np.argmin(abs(x-limits[1]))

        return [x[min([index1, index2]):max([index1, index2])], y[min([index1, index2]):max([index1, index2])]]


    # import sizes

    def import_sizes(file, col = 0, delimiter = None):
        sizes = []
        with open(file, "r") as f:
            for row in f:
                if delimiter:
                    size_string = row.split(delimiter)[col]
                else: # means there is only a single column
                    size_string = row


                if size_string[0].isdigit():
                    sizes.append(float(size_string))

        return sizes

    return (
        double_integral,
        get_lineshape_points,
        import_epr_data,
        import_sizes,
        trim_data,
    )


@app.cell
def _(np):
    # LINESHAPES

    # lorentzian

    def lorentzian(x, I, mu, gamma):
        '''
        sigma = 1/2 FWHM
        '''
        return (I/np.pi)*gamma*((x-mu)**2 + gamma**2)**-1

    # voight
    from scipy.special import voigt_profile

    def voigt(x, I, mu, sigma, gamma):
        """
        Calculates a Voigt profile.

        This function serves as a user-friendly wrapper for scipy.special.voigt_profile,
        mapping common spectroscopic parameters to the function's arguments.

        The Voigt profile is a convolution of a Gaussian and a Lorentzian profile.

        Args:
            x (np.ndarray): Array of independent variable values (e.g., frequency,
                          wavenumber, or wavelength).
            amplitude (float): The integrated area of the profile. A value of 1.0
                             normalizes the profile to unit area.
            center (float): The center position (mean) of the profile (x₀).
            sigma (float): The standard deviation of the Gaussian component.
            gamma (float): The half-width at half-maximum (HWHM) of the
                         Lorentzian component.

        Returns:
            np.ndarray: The calculated Voigt profile values corresponding to each x.
        """
        # voigt_profile requires sigma and gamma, and returns a profile with unit area.
        # We can scale it by the desired amplitude.
        # Note: The amplitude here corresponds to the *area* of the profile.
        # If you want to control the peak height instead, the scaling is more complex.
        profile = voigt_profile(x - mu, sigma, gamma)

        # The amplitude scales the area of the profile.
        # To make the amplitude control the peak height instead of area, you would
        # need to divide by the peak height of the unit-area profile.
        # For now, this function treats amplitude as the integrated area.

        return I * profile


    return lorentzian, voigt


@app.cell
def _(hilbert, np):
    # modifiers of lineshapes

    def is_uniform(x, atol=1e-8):
        """
        Checks if a 1D array is uniformly spaced within a given tolerance.

        Args:
            x (np.ndarray): The array to check.
            atol (float): The absolute tolerance for checking floating point equality.

        Returns:
            bool: True if the array is uniformly spaced, False otherwise.
        """
        if len(x) < 2:
            return True
        # Calculate the differences between consecutive elements
        diffs = np.diff(x)
        # Check if all differences are close to the first difference
        return np.allclose(diffs, diffs[0], atol=atol)

    def hilbert_transform(x, y, tolerance=1e-8):
        """
        Computes the Hilbert transform of a signal (y vs. x), automatically
        handling both uniform and non-uniform spacing of the x-axis.

        - If x is uniform, it computes the transform directly.
        - If x is non-uniform, it resamples to a uniform grid, computes the
          transform, and then resamples the result back to the original x-grid.

        Args:
            x (np.ndarray): The x-axis data.
            y (np.ndarray): The y-axis data (signal).
            tolerance (float): The tolerance used to check for x-axis uniformity.

        Returns:
            tuple[np.ndarray, np.ndarray]: A tuple containing:
                - hilbert_transform (np.ndarray): The Hilbert transform corresponding
                                                  to the original x-axis.
                - amplitude_envelope (np.ndarray): The amplitude envelope corresponding
                                                   to the original x-axis.
        """
        x = np.asarray(x)
        y = np.asarray(y)

        if x.shape != y.shape:
            raise ValueError("Input arrays x and y must have the same shape.")

        # 1. Check if the x-axis is uniform
        if is_uniform(x, atol=tolerance):
            #print("-> Data is uniform. Applying Hilbert transform directly.")
            analytic_signal = hilbert(y)

            hilbert_transform = np.imag(analytic_signal)
            amplitude_envelope = np.abs(analytic_signal)

        else:
            #print("-> Data is non-uniform. Resampling before and after transform.")

            # --- Step A: Resample to a uniform grid ---
            num_points = len(x)
            x_uniform = np.linspace(x.min(), x.max(), num_points)
            y_uniform = np.interp(x_uniform, x, y) # Note the argument order!

            # --- Step B: Compute Hilbert transform on the uniform data ---
            analytic_uniform = hilbert(y_uniform)
            hilbert_uniform = np.imag(analytic_uniform)
            envelope_uniform = np.abs(analytic_uniform)

            # --- Step C: Resample the results BACK to the original x-grid ---
            hilbert_transform = np.interp(x, x_uniform, hilbert_uniform)
            amplitude_envelope = np.interp(x, x_uniform, envelope_uniform)

        return hilbert_transform, amplitude_envelope
    return (hilbert_transform,)


@app.cell
def _():
    # powder pattern generator


    return


@app.cell
def _():
    # accounting for population

    return


@app.cell
def _(np):
    # size dependencies

    # g-factor
    def get_resonance_position_from_size(particle_radius, mu_0, mu_volume, mu_surface, atom_radius= 0, fractional_coverage = 1):

        particle_volume = 4/3*np.pi*particle_radius**3
        shell_volume    = particle_volume-(particle_radius - 2*atom_radius)**3
        atom_volume     = 4/3*np.pi*atom_radius**3
    
    
        #return mu_0 - mu_volume/(particle_volume/atom_volume)**3 + mu_surface*fractional_coverage * (shell_volume/atom_volume)/(particle_volume/atom_volume) 
        return mu_0 + mu_volume*(atom_radius/particle_radius)**3 + mu_surface*fractional_coverage * (particle_radius**3 - (particle_radius-2*atom_radius)**3)/(particle_radius**3) 

    # relaxation
    def get_linewidth_from_size(size, w_0, w_matrix, offset):
    
        return w_0 + w_matrix/(size-offset)**2
    return get_linewidth_from_size, get_resonance_position_from_size


@app.cell
def _(mo):
    mo.md(r"""# Fitting models""")
    return


@app.cell
def _(
    Path,
    get_lineshape_points,
    import_epr_data,
    import_sizes,
    qp,
    trim_data,
):
    # import the data we will be using
    epr_file = Path(r"C:\Users\benle\Downloads\20240325_125158009_240325_PdSC12tol_6K_0dB_wide.csv")

    size_file = Path(r"C:\Users\benle\Downloads\Feret.csv")


    exp_data, instrument_frequency = import_epr_data(epr_file)

    sizes = import_sizes(size_file)

    exp_plot = qp.quickScatter(x = exp_data[0], y = exp_data[1])


    trimmed_exp_data = trim_data(exp_data[0], exp_data[1], [285, 407])
    trimmed_exp_plot = qp.quickScatter(x = trimmed_exp_data[0], y = trimmed_exp_data[1], output = None)
    exp_lineshape_indicies = get_lineshape_points(trimmed_exp_data)

    for i in exp_lineshape_indicies:
        trimmed_exp_plot.add_scatter(x = [trimmed_exp_data[0][i], trimmed_exp_data[0][i]], y = [trimmed_exp_data[1][exp_lineshape_indicies[0]], trimmed_exp_data[1][exp_lineshape_indicies[-1]]], mode = "lines")

    trimmed_exp_plot.update_layout(title = f"instrument frequency = {instrument_frequency} Hz")
    trimmed_exp_plot.update_traces(showlegend = False)
    trimmed_exp_plot.show()
    return exp_lineshape_indicies, sizes, trimmed_exp_data


@app.function
def linear_background(x, a, b):
    return a + b*x


@app.cell
def _(
    Model,
    double_integral,
    exp_lineshape_indicies,
    lorentzian,
    np,
    plot_fit,
    trimmed_exp_data,
):
    # Model1: Lorentzian 
    def lorentzian_epr(x, I, mu, gamma, a, b):
        absorbance = lorentzian(x, I, mu, gamma) 
        epr = np.gradient(absorbance, x)

        return epr + linear_background(x, a, b)

    Model1 = Model(lorentzian_epr)
    Model1_params = Model1.make_params()
    Model1_params.add_many(
        ("I", double_integral(trimmed_exp_data[0], trimmed_exp_data[1]-np.mean(trimmed_exp_data[1]))[1][-1], True, 0, None),
        ("mu", trimmed_exp_data[0][exp_lineshape_indicies[1]], True, 0, None),
        ("gamma", abs(trimmed_exp_data[0][exp_lineshape_indicies[0]] - trimmed_exp_data[0][exp_lineshape_indicies[-1]]), True, 0, None),
        ("a", np.mean(trimmed_exp_data[1]) , True, None, None),
        ("b", 0, True, None, None),
    )

    Model1_result = Model1.fit(trimmed_exp_data[1], params = Model1_params, x = trimmed_exp_data[0])

    print(Model1_result.fit_report())

    plot_fit(Model1_result)

    return (Model1_result,)


@app.cell
def _(Model1_result):
    Model1_result.params["I"].value
    return


@app.cell
def _(Model, Model1_result, np, plot_fit, trimmed_exp_data, voigt):
    # Model 2: Voight
    def voigt_epr(x, I, mu, sigma, gamma, a, b):
        absorbance = voigt(x, I, mu, sigma, gamma) 
        epr = np.gradient(absorbance, x)

        return epr + linear_background(x, a, b)

    Model2 = Model(voigt_epr)
    Model2_params = Model2.make_params()
    Model2_params.add_many(
        ("I",     Model1_result.params["I"].value, True, 0, None),
        ("mu",    Model1_result.params["mu"].value, True, 0, None),
        ("sigma", Model1_result.params["gamma"].value, True, 0, None),
        ("gamma", Model1_result.params["gamma"].value, True, 0, None),
        ("a",     Model1_result.params["a"].value , True, None, None),
        ("b",     Model1_result.params["b"].value, True, None, None),
    )

    Model2_result = Model2.fit(trimmed_exp_data[1], params = Model2_params, x = trimmed_exp_data[0])
    print(Model2_result.aic)
    print(Model2_result.fit_report())
    plot_fit(Model2_result)
    return Model2_params, Model2_result


@app.cell
def _(
    Model,
    Model1_result,
    hilbert_transform,
    lorentzian,
    np,
    plot_fit,
    trimmed_exp_data,
):
    # Model 3: Lorentzian + dispersion
    def lorentzian_dispersion_epr(x, I, D, mu, gamma, a, b):
        # 1. Calculate the dynamic absorption component based on current parameters
        absorption = lorentzian(x, I, mu, gamma) 

        # 2. Calculate the corresponding dynamic dispersion from THAT absorption
        #    We use the 'smart' function that handles uniform/non-uniform x
        #    The [0] selects the Hilbert transform (dispersive part).
        dispersion, _ = hilbert_transform(x, absorption) 

        # 3. Now mix the two components, which are a physically linked pair
        sum_signal = (1 - D) * absorption + D * dispersion

        # 4. Take the derivative for the EPR signal
        epr = np.gradient(sum_signal, x)

        # 5. Add the linear background
        return epr + linear_background(x, a, b)

    Model3 = Model(lorentzian_dispersion_epr)
    Model3_params = Model3.make_params()
    Model3_params.add_many(
        ("I",     Model1_result.params["I"].value, True, 0, None),
        ("D",     0.5, True, 0, 1),
        ("mu",    Model1_result.params["mu"].value, True, 0, None),
        ("gamma", Model1_result.params["gamma"].value, True, 0, None),
        ("a",     Model1_result.params["a"].value, True, None, None),
        ("b",     Model1_result.params["b"].value, True, None, None),
    )

    Model3_result = Model3.fit(trimmed_exp_data[1], params = Model3_params, x = trimmed_exp_data[0])

    print(Model3_result.aic)
    print(Model3_result.fit_report())
    plot_fit(Model3_result)
    return (Model3_result,)


@app.cell
def _(
    Model,
    Model2_result,
    Model3_result,
    hilbert_transform,
    np,
    plot_fit,
    trimmed_exp_data,
    voigt,
):
    # Model 4: Voight + dispersion
    def voigt_dispersion_epr(x, I, D, mu, sigma, gamma, a, b):
        # 1. Calculate the dynamic absorption component based on current parameters
        absorption = voigt(x, I, mu, sigma, gamma) 

        # 2. Calculate the corresponding dynamic dispersion from THAT absorption
        #    We use the 'smart' function that handles uniform/non-uniform x
        #    The [0] selects the Hilbert transform (dispersive part).
        dispersion, _ = hilbert_transform(x, absorption) 

        # 3. Now mix the two components, which are a physically linked pair
        sum_signal = (1 - D) * absorption + D * dispersion

        # 4. Take the derivative for the EPR signal
        epr = np.gradient(sum_signal, x)

        # 5. Add the linear background
        return epr + linear_background(x, a, b)

    Model4 = Model(voigt_dispersion_epr)
    Model4_params = Model4.make_params()
    Model4_params.add_many(
        ("I",     Model3_result.params["I"].value, True, 0, None),
        ("D",     Model3_result.params["D"].value, True, 0, 1),
        ("mu",    Model3_result.params["mu"].value, True, 0, None),
        ("sigma", Model2_result.params["sigma"].value, True, 0, None),
        ("gamma", Model2_result.params["gamma"].value, True, 0, None),
        ("a",     Model3_result.params["a"].value , True, None, None),
        ("b",     Model3_result.params["b"].value, True, None, None),
    )

    Model4_result = Model4.fit(trimmed_exp_data[1], params = Model4_params, x = trimmed_exp_data[0])

    print(Model4_result.aic)
    print(Model4_result.fit_report())
    plot_fit(Model4_result)

    return


@app.cell
def _():
    # Model 5: Powder Lorentzian
    return


@app.cell
def _():
    # Model 6: Powder Voight
    return


@app.cell
def _(
    Model,
    Model1_result,
    exp_lineshape_indicies,
    get_linewidth_from_size,
    get_resonance_position_from_size,
    lorentzian,
    np,
    plot_fit,
    sizes,
    trimmed_exp_data,
):
    # Model 7: Lorentzian size effects
    def lorentzian_sizes_epr(x, sizes, I, mu_0, mu_volume, mu_surface, gamma_0, gamma_volume, gamma_offset, a, b):
        summed_epr = np.zeros_like(x)
        for s in sizes:
            mu = get_resonance_position_from_size(s/2, mu_0, mu_volume, mu_surface, atom_radius = 0.14) # divide by 2 to get the radius
            gamma = get_linewidth_from_size(s/2, gamma_0, gamma_volume, gamma_offset)
            absorbance = lorentzian(x, 1, mu, gamma) 
            summed_epr = summed_epr + absorbance

        epr = np.gradient(I*summed_epr, x)

        return epr + linear_background(x, a, b)

    Model7 = Model(lorentzian_sizes_epr, independent_vars=['x', 'sizes'])
    Model7_params = Model7.make_params()
    Model7_params.add_many(
        ("I",            Model1_result.params["I"].value/len(sizes), True, 0, None),
        ("mu_0",         trimmed_exp_data[0][exp_lineshape_indicies[2]], True, 0, None), # start at one end of the feature
        ("mu_volume",    0, True, None, None),
        ("mu_surface",   0, True, None, None),
        ("gamma_0",      Model1_result.params["gamma"].value/len(sizes)**0.5, True, 0, None),
        ("gamma_volume", 0, False, None, None),
        ("gamma_offset", 0, False, 0, None),
        ("a",            Model1_result.params["a"].value, True, None, None),
        ("b",            Model1_result.params["b"].value, True, None, None),
    )

    Model7_result = Model7.fit(trimmed_exp_data[1], params = Model7_params, x = trimmed_exp_data[0], sizes = sizes)
    print(Model7_result.aic)
    print(Model7_result.fit_report())

    plot_fit(Model7_result)
    return (Model7_result,)


@app.cell
def _(Model7_result):
    Model7_result.params["mu_volume"].value
    return


@app.cell
def _(
    Model,
    Model2_params,
    Model7_result,
    get_linewidth_from_size,
    get_resonance_position_from_size,
    np,
    plot_fit,
    sizes,
    trimmed_exp_data,
    voigt,
):
    # Model 8: Voight size effects
    # ? Maybe only a single value for sigma? Like the degree of heterogenous broadening will be the same for everything?
    def voight_sizes_epr(x, sizes, I, mu_0, mu_volume, mu_surface, gamma_0, gamma_volume, gamma_offset, sigma_0, sigma_volume, sigma_surface, a, b):
        summed_epr = np.zeros_like(x)
        for s in sizes:
            mu = get_resonance_position_from_size(s/2, mu_0, mu_volume, mu_surface, atom_radius = 0.14)
            gamma = get_linewidth_from_size(s/2, gamma_0, gamma_volume, gamma_offset)
            sigma = get_resonance_position_from_size(s/2, sigma_0, sigma_volume, sigma_surface, atom_radius = 0.14)
            absorbance = voigt(x, 1, mu, sigma, gamma) 
            summed_epr = summed_epr + absorbance

        epr = np.gradient(I*summed_epr, x)

        return epr + linear_background(x, a, b)

    Model8 = Model(voight_sizes_epr, independent_vars=['x', 'sizes'])
    Model8_params = Model8.make_params()
    Model8_params.add_many(
        ("I",            Model7_result.params["I"].value, True, 0, None),
        ("mu_0",         Model7_result.params["mu_0"].value, True, 0, None), # start at one end of the feature
        ("mu_volume",    Model7_result.params["mu_volume"].value, True, None, None),
        ("mu_surface",   Model7_result.params["mu_surface"].value, True, None, None),
        ("sigma_0",      Model7_result.params["gamma_0"].value*Model2_params["sigma"].value/Model2_params["gamma"].value, True, 0, None),
        ("sigma_volume",  0 , False, None, None),
        ("sigma_surface", 0 , False, None, None),
        ("gamma_0",       Model7_result.params["gamma_0"].value, True, 0, None),
        ("gamma_volume",  0, False, None, None),
        ("gamma_offset",  0, False, 0, None),
        ("a",             Model7_result.params["a"].value, True, None, None),
        ("b",             Model7_result.params["b"].value, True, None, None),
    )

    Model8_result = Model8.fit(trimmed_exp_data[1], params = Model8_params, x = trimmed_exp_data[0], sizes = sizes)
    print(Model8_result.aic)
    print(Model8_result.fit_report())
    plot_fit(Model8_result)
    return (Model8_result,)


@app.cell
def _(
    Model,
    Model7_result,
    get_linewidth_from_size,
    get_resonance_position_from_size,
    hilbert_transform,
    lorentzian,
    np,
    plot_fit,
    sizes,
    trimmed_exp_data,
):
    # Model 9: Size dependent dysonian
    # ? Maybe only a single value for sigma? Like the degree of heterogenous broadening will be the same for everything?
    def lorentzian_dispersion_sizes_epr(x, sizes, I, mu_0, mu_volume, mu_surface, gamma_0, gamma_volume, gamma_offset, D_0, D_volume, D_surface, a, b):
        summed_epr = np.zeros_like(x)
        for s in sizes:
            mu = get_resonance_position_from_size(s, mu_0, mu_volume, mu_surface, offset = 0.28)
            gamma = get_linewidth_from_size(s, gamma_0, gamma_volume, gamma_offset)
            D = get_resonance_position_from_size(s, D_0, D_volume, D_surface)
            absorption = lorentzian(x, 1, mu, gamma) 
            dispersion, _ = hilbert_transform(x, absorption) 

            # 3. Now mix the two components, which are a physically linked pair
            sum_signal = (1 - D) * absorption + D * dispersion
            summed_epr = summed_epr + sum_signal

        epr = np.gradient(I*summed_epr, x)

        return epr + linear_background(x, a, b)

    Model9 = Model(lorentzian_dispersion_sizes_epr, independent_vars=['x', 'sizes'])
    Model9_params = Model9.make_params()
    Model9_params.add_many(
        ("I", Model7_result.params["I"].value, True, 0, None),
        ("mu_0", Model7_result.params["mu_0"].value, True, 0, None), # start at one end of the feature
        ("mu_volume", Model7_result.params["mu_volume"].value, True, None, None),
        ("mu_surface", Model7_result.params["mu_surface"].value, True, None, None),
        ("D_0",        0.5, True, 0, None),
        ("D_volume",   0, False, None, None),
        ("D_surface",  0, False, None, None),
        ("gamma_0",    Model7_result.params["gamma_0"].value, True, 0, None),
        ("gamma_volume", Model7_result.params["gamma_volume"].value, False, None, None),
        ("gamma_offset", Model7_result.params["gamma_offset"].value, True, None, None),
        ("a", Model7_result.params["a"].value, True, None, None),
        ("b", Model7_result.params["b"].value, True, None, None),
    )

    Model9_result = Model9.fit(trimmed_exp_data[1], params = Model9_params, x = trimmed_exp_data[0], sizes = sizes)
    print(Model9_result.aic)
    print(Model9_result.fit_report())

    plot_fit(Model9_result)
    return (Model9_result,)


@app.cell
def _(
    Model,
    Model8_result,
    Model9_result,
    get_linewidth_from_size,
    get_resonance_position_from_size,
    hilbert_transform,
    np,
    plot_fit,
    sizes,
    trimmed_exp_data,
    voigt,
):
    # Model 9: Size dependent dysonian
    # ? Maybe only a single value for sigma? Like the degree of heterogenous broadening will be the same for everything?
    def voigt_dispersion_sizes_epr(x, sizes, I, mu_0, mu_volume, mu_surface, sigma, gamma_0, gamma_volume, gamma_offset, D_0, D_volume, D_surface, a, b):
        summed_epr = np.zeros_like(x)
        for s in sizes:
            mu = get_resonance_position_from_size(s, mu_0, mu_volume, mu_surface)
            gamma = get_linewidth_from_size(s, gamma_0, gamma_volume, gamma_offset)
            D = get_resonance_position_from_size(s, D_0, D_volume, D_surface)
            absorption = voigt(x, 1, mu, sigma, gamma) 
            dispersion, _ = hilbert_transform(x, absorption) 

            # 3. Now mix the two components, which are a physically linked pair
            sum_signal = (1 - D) * absorption + D * dispersion
            summed_epr = summed_epr + sum_signal

        epr = np.gradient(I*summed_epr, x)

        return epr + linear_background(x, a, b)

    Model10 = Model(voigt_dispersion_sizes_epr, independent_vars=['x', 'sizes'])
    Model10_params = Model10.make_params()
    Model10_params.add_many(
        ("I", Model8_result.params["I"].value, True, 0, None),
        ("mu_0", Model8_result.params["mu_0"].value, True, 0, None), # start at one end of the feature
        ("mu_volume", Model8_result.params["mu_volume"].value, True, None, None),
        ("mu_surface", Model8_result.params["mu_surface"].value, True, None, None),
        ("D_0",        Model9_result.params["D_0"].value, True, 0, None),
        ("D_volume",   Model9_result.params["D_volume"], False, None, None),
        ("D_surface",  Model9_result.params["D_surface"], False, None, None),
        ("sigma", Model8_result.params["sigma_0"], False, None, None),
        ("gamma_0",    Model8_result.params["gamma_0"].value, True, 0, None),
        ("gamma_volume", Model8_result.params["gamma_volume"].value, False, None, None),
        ("gamma_offset", Model8_result.params["gamma_offset"].value, True, None, None),
        ("a", Model8_result.params["a"].value, True, None, None),
        ("b", Model8_result.params["b"].value, True, None, None),
    )

    Model10_result = Model10.fit(trimmed_exp_data[1], params = Model10_params, x = trimmed_exp_data[0], sizes = sizes)
    print(Model10_result.aic)
    print(Model10_result.fit_report())

    plot_fit(Model10_result)
    return


@app.function
def EPR_fitting_with_sizes(exp_x, exp_y, frequency, sizes):
    # first, fit the spectrum to a single Lorentzian

    # then fit a single Voigt 

    # then fit a voigt with size effects but no gamma or sigma size dependence, dividing the gamma and sigma by the length of the sizes...
    
    return


if __name__ == "__main__":
    app.run()
