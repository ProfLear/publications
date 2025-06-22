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
    from scipy.constants import physical_constants

    # constants we will need, assumes the magnetic field is in T
    h = physical_constants["Planck constant"][0] # J s
    u_B = physical_constants["Bohr magneton"][0] # J T-1
    g_e = physical_constants["electron g factor"] # no units

    n_mesh_points = 500 # number of points in the sphere we want to sample

    # Model 7: Lorentzian size effects
    g_metal = 2.15
    return (
        Model,
        Path,
        cumulative_trapezoid,
        g_metal,
        h,
        hilbert,
        make_subplots,
        mo,
        n_mesh_points,
        np,
        qp,
        u_B,
    )


@app.function
def convert_B_to_g():
    return


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

        plot.update_xaxes(title = "magnetic field /T", showline = False)
        plot.update_yaxes(title = "intensity", range = [-1.05*ymax, 1.05*ymax], zeroline = True)
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
        # need to import magnetic field in units of T
        B_fields, intensities = [], []

        instrument_ID = get_instrument_ID(file)
        print(instrument_ID)
        if instrument_ID ==";":
            with open(file, "r") as f:
                for line in f:
                    if "Frequency;" in line:
                        instrument_frequency = float(line.split(";")[1])*10**9 # convert to Hz
                    if line[0].isdigit():
                        B, I = line.split(";")
                        B_fields.append(float(B)/1) # convert to mT
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
                        B_fields.append(float(B)/10) # convert to mT
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
def _(n_mesh_points, np):
    # powder pattern generator stuff

    def sphere_sampling(mesh_sphere) :
        import numpy as np
        gr = (1 + 5**0.5)/2  #golden ratio
    
        #this is one approach, that uses the appraoch from LATTICE 3 in Extreme learning: 
        #http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/
    
        #we are going to generate a paired list of points. That is, two lists, where common idexes represent a pair of values. 

        #start with the i = 0 points.  This will be theta = 0, phi = 0...
        # we are starting with x and y, because the appraoch first generates the uniform sampling of a rectangle, and then converts this to a sphere. 
        # so we start in cartesian space...
    
        #for speed, we are going to work with a list first, and then convert this to an array once we have the full length...
        x = [0]
        y = [0]

        #then move on to i = 1 through i = mesh_sphere-1
        #i is both an index and the way we calculate the values
        i = 1

        while i < mesh_sphere :
            x.append((i + 6) / (mesh_sphere + 11))
            y.append(i / gr)
            i += 1

        #finally, the value for i = mesh_sphere

        x.append(1)
        y.append(0)
    
        #then we convert to an array to make math easy... 
        x = np.array(x)
        y = np.array(y)
    
        y = y % 1  #we want to convert to the unit square, and we do this by using the modulo operator
    
        # then we want to convert the unit square to equally spaced polar coordinates. 
        thetas = np.arccos(2*x-1)
        phis = 2*np.pi*y
    
        #convert the array back to list, so that sorting through it is easier...
        thetas = np.array(thetas).tolist()
        phis = np.array(phis).tolist()

        #print(thetas)

        #########
        # REDUCE THE SAMPLING TO A QUADRANT
        ########
        # the above gives a uniform sampling of a sphere, but we really only need 1/8 of a sphere. We need theta 0 -> pi/2 and phi 0 -> pi/2  Thus, the simulation will be 8x faster!  So, let us trim this...

        #below is one way to do this.  There is CERTAINLY a more elegant way to do this...
    
        #first, collect the thetas that are correct...
        # we are going to step through the theta vector, and pull out all thetas with the right value, as well as their paired value in phi..
        tp = 0
        temp_thetas1 = []
        temp_phis1   = []
        while tp < len(thetas) :
            if thetas[tp] <= np.pi/2 :
                temp_thetas1.append(thetas[tp])
                temp_phis1.append(phis[tp])
            tp += 1  #advance the counter
      
    
        #then, collect the phis that are correct...
        tp = 0
        temp_thetas2 = []
        temp_phis2 = []
        while tp < len(temp_phis1) :
            if temp_phis1[tp] <= np.pi/2 :
                temp_thetas2.append(temp_thetas1[tp])
                temp_phis2.append(temp_phis1[tp])
            tp += 1  #advance the counter
    
    

        thetas = temp_thetas2
        phis = temp_phis2
    
        #lets keep the output as lists for right now.  I think that is good for speeding up the latter portions...
        # we can always convert them to arrays later, if we want. 
        return [thetas, phis]



    thetas_phis = sphere_sampling(n_mesh_points*8)

    def get_values_from_orientation(x_value, y_value, z_value, orientations):
        thetas, phis = orientations # unpack this
        return  (x_value**2 * (np.sin(thetas) * np.cos(phis))**2 + y_value**2 * (np.sin(thetas) * np.sin(phis))**2 + z_value**2 * np.cos(thetas)**2 )**0.5


    return get_values_from_orientation, thetas_phis


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
    epr_file = Path(r"C:\Users\benle\Documents\GitHub\publications\Kristen Variability\20240325_125158009_240325_PdSC12tol_6K_0dB_wide.csv")

    size_file = Path(r"C:\Users\benle\Documents\GitHub\publications\Kristen Variability\20240325_125158009_240325_PdSC12tol_6K_0dB_wide SIZES.csv")


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
    return (
        exp_lineshape_indicies,
        instrument_frequency,
        sizes,
        trimmed_exp_data,
    )


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
        ("gamma", Model2_result.params["gamma"].value/10, True, 0, None),# start with this as a minor component
        ("a",     Model3_result.params["a"].value , True, None, None),
        ("b",     Model3_result.params["b"].value, True, None, None),
    )

    Model4_result = Model4.fit(trimmed_exp_data[1], params = Model4_params, x = trimmed_exp_data[0])

    print(Model4_result.aic)
    print(Model4_result.fit_report())
    plot_fit(Model4_result)

    return


@app.cell
def _(
    Model,
    Model1_result,
    exp_lineshape_indicies,
    get_values_from_orientation,
    lorentzian,
    np,
    plot_fit,
    thetas_phis,
    trimmed_exp_data,
):
    # Model 5: Powder Lorentzian
    def lorentzian_powder_epr(x, I, Bx, By, Bz, gamma_x, gamma_y, gamma_z, a, b, thetas_phis):

        summed_absorptions = np.zeros_like(x)
    
        # first, get the vectors of the g-value, and gamma
        B_resonances = get_values_from_orientation(Bx, By, Bz, thetas_phis)
        gammas = get_values_from_orientation(gamma_x, gamma_y, gamma_z, thetas_phis)

        # then get the mean position in B-field
        #B_resonances = 1000* h*instrument_frequency / (gs * u_B) # get resonance position in terms of mT

        for B, g in zip(B_resonances, gammas):
            summed_absorptions = summed_absorptions+lorentzian(x, 1, B, g)

        epr = np.gradient(I*summed_absorptions, x)

        return epr+linear_background(x, a, b)

    Model5 = Model(lorentzian_powder_epr, independent_vars = ["x", "thetas_phis"])
    Model5_params = Model5.make_params()
    Model5_params.add_many(
        ("I",       Model1_result.params["I"].value/len(thetas_phis), True, 0, None),
        ("Bx",      trimmed_exp_data[0][exp_lineshape_indicies[0]], True, trimmed_exp_data[0][0], trimmed_exp_data[0][-1]),
        ("By",      trimmed_exp_data[0][exp_lineshape_indicies[1]], True, trimmed_exp_data[0][0], trimmed_exp_data[0][-1]),
        ("Bz",      trimmed_exp_data[0][exp_lineshape_indicies[2]], True, trimmed_exp_data[0][0], trimmed_exp_data[0][-1]),
        ("gamma_x", Model1_result.params["gamma"].value/3      , True, 0, None),
        ("gamma_y", Model1_result.params["gamma"].value/3      , True, 0, None),
        ("gamma_z", Model1_result.params["gamma"].value/3      , True, 0, None),
        ("a",       Model1_result.params["a"].value, True, None, None),
        ("b",       Model1_result.params["b"].value, True, None, None),
        )

    Model5_result = Model5.fit(trimmed_exp_data[1], params = Model5_params, x = trimmed_exp_data[0], thetas_phis = thetas_phis)
    print(Model5_result.aic)
    print(Model5_result.fit_report())

    plot_fit(Model5_result)
    return


@app.cell
def _(
    Model,
    Model1_result,
    Model2_result,
    exp_lineshape_indicies,
    get_values_from_orientation,
    np,
    plot_fit,
    thetas_phis,
    trimmed_exp_data,
    voigt,
):
    # Model 6: Powder Voight

    def voigt_powder_epr(x, I, Bx, By, Bz, sigma_x, sigma_y, sigma_z, gamma_x, gamma_y, gamma_z, a, b, thetas_phis):

        summed_absorptions = np.zeros_like(x)
    
        # first, get the vectors of the g-value, and gamma
        B_resonances = get_values_from_orientation(Bx, By, Bz, thetas_phis)
        sigmas = get_values_from_orientation(sigma_x, sigma_y, sigma_z, thetas_phis)
        gammas = get_values_from_orientation(gamma_x, gamma_y, gamma_z, thetas_phis)

        # then get the mean position in B-field
        #B_resonances = 1000* h*instrument_frequency / (gs * u_B) # get resonance position in terms of mT

        for B, s, g in zip(B_resonances, sigmas, gammas):
            summed_absorptions = summed_absorptions+voigt(x, 1, B, s, g)

        epr = np.gradient(I*summed_absorptions, x)

        return epr+linear_background(x, a, b)

    Model6 = Model(voigt_powder_epr, independent_vars = ["x", "thetas_phis"])
    Model6_params = Model6.make_params()
    Model6_params.add_many(
        ("I",       Model1_result.params["I"].value/len(thetas_phis), True, 0, None),
        ("Bx",      trimmed_exp_data[0][exp_lineshape_indicies[0]],   True, trimmed_exp_data[0][0], trimmed_exp_data[0][-1]),
        ("By",      trimmed_exp_data[0][exp_lineshape_indicies[1]],   True, trimmed_exp_data[0][0], trimmed_exp_data[0][-1]),
        ("Bz",      trimmed_exp_data[0][exp_lineshape_indicies[2]],   True, trimmed_exp_data[0][0], trimmed_exp_data[0][-1]),
        ("sigma_x", Model2_result.params["sigma"].value/3,            True, 0, None),
        ("sigma_y", Model2_result.params["sigma"].value/3,            True, 0, None),
        ("sigma_z", Model2_result.params["sigma"].value/3,            True, 0, None),
        ("gamma_x", Model2_result.params["gamma"].value/3,            True, 0, None),
        ("gamma_y", Model2_result.params["gamma"].value/3,            True, 0, None),
        ("gamma_z", Model2_result.params["gamma"].value/3,            True, 0, None),
        ("a",       Model1_result.params["a"].value,                  True, None, None),
        ("b",       Model1_result.params["b"].value,                  True, None, None),
        )

    Model6_result = Model6.fit(trimmed_exp_data[1], params = Model6_params, x = trimmed_exp_data[0], thetas_phis = thetas_phis)
    print(Model6_result.aic)
    print(Model6_result.fit_report())

    plot_fit(Model6_result)
    return (Model6_result,)


@app.cell
def _(Model6_result, h, instrument_frequency, u_B):
    gx = h*instrument_frequency/(u_B * Model6_result.params["Bx"]/1000)
    gy = h*instrument_frequency/(u_B * Model6_result.params["By"]/1000)
    gz = h*instrument_frequency/(u_B * Model6_result.params["Bz"]/1000)

    print(f"gx = {gx:.4f}")
    print(f"gy = {gy:.4f}")
    print(f"gz = {gz:.4f}")
    print(f"<g> = {(gx+gy+gz)/3:.4f}")
    return


@app.cell
def _(
    Model,
    Model1_result,
    g_metal,
    get_linewidth_from_size,
    get_resonance_position_from_size,
    h,
    instrument_frequency,
    lorentzian,
    np,
    plot_fit,
    sizes,
    trimmed_exp_data,
    u_B,
):

    def lorentzian_sizes_epr(x, sizes, I, mu_0, mu_volume, mu_surface, gamma_0, gamma_volume, gamma_offset, a, b):
        summed_epr = np.zeros_like(x)
        for s in sizes:
            mu = get_resonance_position_from_size(s/2, mu_0, mu_volume, mu_surface, atom_radius = 0.14) # divide size by 2 to get the radius
            gamma = get_linewidth_from_size(s/2, gamma_0, gamma_volume, gamma_offset)
            absorbance = lorentzian(x, 1, mu, gamma) 
            summed_epr = summed_epr + absorbance

        epr = np.gradient(I*summed_epr, x)

        return epr + linear_background(x, a, b)

    Model7 = Model(lorentzian_sizes_epr, independent_vars=['x', 'sizes'])
    Model7_params = Model7.make_params()
    Model7_params.add_many(
        ("I",            Model1_result.params["I"].value/len(sizes), True, 0, None),
        ("mu_0",         1000* h*instrument_frequency / (u_B*g_metal), True, 0, None), #calculate the resonance position from g_metal, and don't change it 1000* h*instrument_frequency / (u_B*g_metal) # start at one end of the feature trimmed_exp_data[0][exp_lineshape_indicies[0]]
        ("mu_volume",    12737, True, None, None),
        ("mu_surface",   -16, True, None, None),
        ("gamma_0",      9.6, True, 0, None),
        ("gamma_volume", 0, False, None, None),
        ("gamma_offset", 0, False, 0, None),
        ("a",            -11031, True, None, None),
        ("b",            32, True, None, None),
    )

    Model7_result = Model7.fit(trimmed_exp_data[1], params = Model7_params, x = trimmed_exp_data[0], sizes = sizes)
    print(Model7_result.aic)
    print(Model7_result.fit_report())

    plot_fit(Model7_result)
    return Model7_params, Model7_result


@app.cell(disabled=True)
def _(
    Model1_result,
    Model7_params,
    exp_lineshape_indicies,
    sizes,
    trimmed_exp_data,
):
    Model7_params.add_many(
        ("I",            Model1_result.params["I"].value/len(sizes), True, 0, None),
        ("mu_0",         trimmed_exp_data[0][exp_lineshape_indicies[0]], True, 0, None), #calculate the resonance position from g_metal, and don't change it 1000* h*instrument_frequency / (u_B*g_metal) # start at one end of the feature trimmed_exp_data[0][exp_lineshape_indicies[0]]
        ("mu_volume",    0, True, None, None),
        ("mu_surface",   0, True, None, None),
        ("gamma_0",      Model1_result.params["gamma"].value/len(sizes)**0.5, True, 0, None),
        ("gamma_volume", 0, False, None, None),
        ("gamma_offset", 0, False, 0, None),
        ("a",            Model1_result.params["a"].value, True, None, None),
        ("b",            Model1_result.params["b"].value, True, None, None),
    )
    return


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
    Model8_result,
    get_resonance_position_from_size,
    h,
    instrument_frequency,
    np,
    qp,
    sizes,
    u_B,
):
    resonance_positions = []
    g_values = []
    for s in sizes:
        resonance_position = get_resonance_position_from_size(s/2, Model8_result.params["mu_0"], Model8_result.params["mu_volume"], Model8_result.params["mu_surface"], atom_radius = 0.14)
        g_values.append(h*instrument_frequency / (u_B*resonance_position/1000)) # factor of 1000 to convert to T
        resonance_positions.append(resonance_position)

    resonance_positions = np.array(resonance_positions)

    g_values = np.array(g_values)
    log_g_mean = np.mean(np.log(g_values))
    log_g_variance = np.std(np.log(g_values))**2

    log_size_mean = np.mean(np.log(sizes))
    log_size_variance = np.std(np.log(sizes))**2

    lognormal_size_mean = np.exp(log_size_mean + log_size_variance/2)
    lognormal_size_mean_resonance = get_resonance_position_from_size(lognormal_size_mean, Model8_result.params["mu_0"], Model8_result.params["mu_volume"], Model8_result.params["mu_surface"], atom_radius = 0.14)

    print(instrument_frequency)
    print(f"g_0 = {h*instrument_frequency/(u_B * Model8_result.params["mu_0"]/1000)}")
    print(f"g_volume = {h*instrument_frequency/(u_B * Model8_result.params["mu_volume"]/1000)}")
    print(f"g_surface = {h*instrument_frequency/(u_B * Model8_result.params["mu_surface"]/1000)}")
    print(f"mean g-value = {np.mean(g_values):.4f}")
    print(f"g-value for mean size = {h*instrument_frequency/(u_B * lognormal_size_mean_resonance/1000):.4f}")

    g_hist = qp.quickHist(x = g_values)
    return


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
