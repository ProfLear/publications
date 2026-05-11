import numpy as np


# functions to support the notebook
def lorentzian(x, A = 1, mu = 1, gamma = 1):
    if isinstance(x, list):
        x = np.array(x)
    return (A / np.pi) * (gamma / ((x - mu)**2 + gamma**2))

def bin_ligand_and_solvent_together(ligand, solvent, bin_width = 1):
    import numpy as np

    left_edge = 0
    bin_centers = []
    ligand_counts = []
    solvent_counts = []

    # first, bin the data we have
    while left_edge < max(max(ligand), max(solvent)):
        bin_centers.append(left_edge + bin_width/2)

        ligand_counts.append(np.sum((ligand >= left_edge) & (ligand < left_edge + bin_width)))

        solvent_counts.append(np.sum((solvent >= left_edge) & (solvent < left_edge + bin_width)))

        left_edge = left_edge + bin_width

    return np.array(bin_centers), np.array(ligand_counts), np.array(solvent_counts)

def extend_bins_based_on_wavefunction(bin_centers, ligand_counts, solvent_counts, alpha = 1, fract_threshold = 0.01):
    extended_bin_centers = []
    extended_sim_solvent = []
    extended_sim_ligand = []
    extended_calc_solvent = []

    #print(f"alpha = {alpha}")

    bin_width = bin_centers[1] - bin_centers[0] # determine the bin width being used 

    latest_fract = 1
    i = 0
    found_ligand = False
    extend = False # flag to keep track of if we have switched to extending the solvent 
    while latest_fract > fract_threshold: # basically go until the contribution from the most recent bin is less than some percent
        # check to see if we should change the extend flag
        if extend:

            new_solvent_value = extended_calc_solvent[-1] * ((extended_bin_centers[-1] + bin_width)**3 - extended_bin_centers[-1]**3)/(extended_bin_centers[-1]**3 - (extended_bin_centers[-1] - bin_width)**3)  

            extended_bin_centers.append(extended_bin_centers[-1] + bin_width)
            extended_calc_solvent.append(new_solvent_value)
            try:
                extended_sim_ligand.append(ligand_counts[i])
            except:
                extended_sim_ligand.append(0)
            try:
                extended_sim_solvent.append(solvent_counts[i])
            except:
                extended_sim_solvent.append(0)


        else: # we are not yet extending the data
            extended_bin_centers.append(bin_centers[i])
            extended_calc_solvent.append(solvent_counts[i])
            extended_sim_solvent.append(solvent_counts[i])
            extended_sim_ligand.append(ligand_counts[i])

            if ligand_counts[i] > 0:
                found_ligand = True


            if i == len(bin_centers)-1: # we are are the end of the simulated data, and we MUST switch
                extend = True
            elif i == len(ligand_counts) and found_ligand: # we have found the end of the ligand, and so we should follow the expected r**2 behavior now
                extend = True
            elif solvent_counts[i] == max(solvent_counts): # we have found the maximum solvent value, which should then just increase as r**2
                extend = True


        # now, using the decay function, figure out the fraction that this latest value represents...

        ligand_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_sim_ligand))
        solvent_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_calc_solvent))

        if ligand_contribution + solvent_contribution > 0 and extend: # use of extend ensure we have found the end of the ligand?
            latest_fract = np.exp(-2*alpha * extended_bin_centers[-1]) * extended_calc_solvent[-1] / (ligand_contribution + solvent_contribution)
        else:
            latest_fract = 1
        #print(extend)
        i = i + 1

    return np.array(extended_bin_centers), np.array(extended_sim_ligand), np.array(extended_sim_solvent), np.round(np.array(extended_calc_solvent)).astype(int)


def search_parameter_grid(alpha_endpoints = None, A_endpoints = None, grid_resolution = 10, space = "linear"):
    if alpha_endpoints[0] != alpha_endpoints[-1]:
        #build grid
        if space == "linear":
            alphas = np.linspace(alpha_endpoints[0], alpha_endpoints[1], grid_resolution)
        elif space == "log":
            alphas = np.logspace(alpha_endpoints[0], alpha_endpoints[1], grid_resolution)
        else:
            raise "You need to specify a valid space"
    else:
        alphas = [alpha_endpoints[0]]

    if A_endpoints[0] != A_endpoints[1]:
        if space == "linear":
            As = np.linspace(A_endpoints[0], A_endpoints[1], grid_resolution)
        elif space == "log":
            As = np.logspace(A_endpoints[0], A_endpoints[1], grid_resolution)
        else:
            raise "You need to specify a valid space"
    else:
        As = [A_endpoints[0]]

    collected_fits = {}

    i = 0
    for alpha in alphas:
        for A in As:
            collected_fits[f"alpha = {alpha}, A = {A}"] = {}
            collected_fits[f"alpha = {alpha}, A = {A}"]["fit result"], collected_fits[f"alpha = {alpha}, A = {A}"]["fit plot"] = fit_mimms_spectra(alpha = alpha, A = A)
            i = i+1
            print(f"fit {i} / {len(alphas) * len(As)} completed. alpha = {alpha}, A = {A}")
    return collected_fits







def fit_mimms_spectra(alpha, A):
    # now, we can fit the high density stuff...
    high_density_triple_fit_params = triple_sample_model.make_params()

    high_density_triple_fit_params.add_many(
        ("alpha",   alpha,  True, 0, None), # decay rate for wavefunction
        ("A",       A,      True, 0, None), # a scalar value that adjust how much the wavefunction overlap leads to coupling

        ("gamma",   gamma,  True, 0, None), # width of lines associated with each contribution

        ("I1",       (max(high_density["experimental spectra"]["H25_H8"][1]) - min(high_density["experimental spectra"]["H25_H8"][1]))/(sum(high_density["distances"]["extended solvent counts"]) + sum(high_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the first spectrum
        ("I2",       (max(high_density["experimental spectra"]["D25_H8"][1]) - min(high_density["experimental spectra"]["D25_H8"][1]))/(sum(high_density["distances"]["extended solvent counts"]))  ,      True, 0, None), # intensity adjuster for the second spectrum
        ("I3",       (max(high_density["experimental spectra"]["H25_D8"][1]) - min(high_density["experimental spectra"]["H25_D8"][1]))/(sum(high_density["distances"]["extended ligand counts"]))  ,      True, 0, None), # intensity adjuster for the third spectrum

        ("mu1",      mu,     True, 48, 52), # center position for the first spectrum
        ("mu2",      mu,     True, 48, 52), # center position for the second spectrum <-- needed?
        ("mu3",      mu,     True, 48, 52), # center position for the second spectrum <--needed?

        ("a1",       min(high_density["experimental spectra"]["H25_H8"][1]),      True, None, None),  # baseline for the 1st spectrum
        ("a2",       min(high_density["experimental spectra"]["D25_H8"][1]),      True, None, None),
        ("a3",       min(high_density["experimental spectra"]["H25_D8"][1]),      True, None, None),

        )


    # now set up to fit them

    high_density_concatenated_xs = np.array([])
    high_density_concatenated_ys = np.array([])
    high_density_spectral_indices = [0]
    for _i, _spectrum in enumerate([high_density["experimental spectra"]["H25_H8"], high_density["experimental spectra"]["D25_H8"], high_density["experimental spectra"]["H25_D8"]]):
        high_density_concatenated_xs = np.concatenate([high_density_concatenated_xs, _spectrum[0]])
        high_density_concatenated_ys = np.concatenate([high_density_concatenated_ys, _spectrum[1]])
        high_density_spectral_indices.append(len(high_density_concatenated_xs)) # append the curent length of the growing array

    high_density_ligand_H_fractions = [1, 0.02, 1]
    high_density_solvent_H_fractions = [1, 1, 0.01]

    high_density_triple_fit_result = triple_sample_model.fit(
        high_density_concatenated_ys,

        x = high_density_concatenated_xs,

        params = high_density_triple_fit_params,

        ligand_binned_centers = high_density["distances"]["extended bin centers"], # to be supplied as not fit
        solvent_binned_centers = high_density["distances"]["extended bin centers"], # to be supplied as not fit

        ligand_binned_counts = high_density["distances"]["extended ligand counts"], # to be supplied as not fit
        solvent_binned_counts = high_density["distances"]["extended solvent counts"], # to be supplied as not fit

        ligand_H_fractions =  high_density_ligand_H_fractions, # to be supplied as not fit
        solvent_H_fractions = high_density_solvent_H_fractions, # to be supplied as not fit

        spectral_indices = high_density_spectral_indices,
        )


    # make a plot
    high_density_triple_fit_plot = make_subplots(rows = 3, cols =1)

    for i in range(3):
        high_density_triple_fit_plot.add_scatter(
            x = high_density_concatenated_xs[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
            y = high_density_concatenated_ys[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
            mode = "lines",
            line = dict(width = 8, color = "grey"),
            row = i + 1, col = 1,
            showlegend = False,
            )
        high_density_triple_fit_plot.add_scatter(
            x = high_density_concatenated_xs[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
            y = high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
            line = dict(width = 4, color = "white"),
            row = i + 1, col = 1,
            showlegend = False,
            )
        high_density_triple_fit_plot.add_scatter(
            x = high_density_concatenated_xs[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
            y = high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]],
            line = dict(width = 2, color = "red"),
            row = i + 1, col = 1,
            showlegend = False,
            )
        high_density_triple_fit_plot.add_annotation(
            text = spectral_labels[i],
            y = max([
                max(high_density_concatenated_ys[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]]),
                max(high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]])
                ]) + (max([
                    max(high_density_concatenated_ys[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]]),
                    max(high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]])
                    ]) - min([
                        min(high_density_concatenated_ys[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]]),
                        min(high_density_triple_fit_result.best_fit[high_density_spectral_indices[i]:high_density_spectral_indices[i+1]])
                        ])) * 0.15,
            showarrow = False,
            row =i + 1, col = 1
            )
    high_density_triple_fit_plot.update_xaxes(    
        title = "magnetic field",
        row = 3, col = 1)

    high_density_triple_fit_plot.update_xaxes(
        range = [46, 54],
        )

    high_density_triple_fit_plot.update_layout(
        template = "simple_white",
        margin = dict(l = 50, r = 10, b = 20, t = 10),
        width = 500, height = 600,
        )

    return high_density_triple_fit_result, high_density_triple_fit_plot


def simulate_pulse_signal(
    x,

    ligand_binned_centers, # to be supplied as not fit
    solvent_binned_centers, # to be supplied as not fit

    ligand_binned_counts, # to be supplied as not fit
    solvent_binned_counts, # to be supplied as not fit

    ligand_H_fraction, # to be supplied as not fit
    solvent_H_fraction, # to be supplied as not fit

    alpha, # decay constant
    A, # scalar value to adjust strength of coupling

    I, # intensity of individual contributions
    mu, # central position, without coupling
    gamma, # width of the lorentzian peak

    a, # baseline correction
    ):

    # for debudgin...
    #print(f"alpha={alpha:.4f}, A={A:.4f}, I={I:.4f}, mu={mu:.4f}, gamma={gamma:.4f}, a={a:.4f}")


    # calculate the psi-squared values at each distance for ligand and solvent
    ligand_psi = np.exp(-alpha * ligand_binned_centers) #included in case we decide we need it
    solvent_psi = np.exp(-alpha * solvent_binned_centers)
    ligand_psi_squareds  = np.exp(-2*alpha * ligand_binned_centers)
    solvent_psi_squareds = np.exp(-2*alpha * solvent_binned_centers)

    # now we are going to sum up contributions
    temp_y_values = np.zeros_like(x)
    for info in [[ligand_psi, ligand_psi_squareds, ligand_binned_counts, ligand_H_fraction], [solvent_psi, solvent_psi_squareds, solvent_binned_counts, solvent_H_fraction]]:
        psis, psi_squareds, counts, H = info # unpack these
        for p, ps, count in zip(psis, psi_squareds, counts):
            temp_y_values = temp_y_values + I*lorentzian(x, mu = mu +   ps*A, gamma = gamma) * H*ps * count
            temp_y_values = temp_y_values + I*lorentzian(x, mu = mu - 1*ps*A, gamma = gamma) * H*ps * count

    return temp_y_values + a



def model_three_samples(
        x, # x values. These are the concatenated x-values from the experimental spectra we want to fit

        ligand_binned_centers, # to be supplied as not fit
        solvent_binned_centers, # to be supplied as not fit

        ligand_binned_counts, # to be supplied as not fit
        solvent_binned_counts, # to be supplied as not fit

        alpha, # decay constant
        A, # scalar value to adjust strength of coupling

        gamma, # width of the lorentzian peak, shared

        I1, # intensity for 1st spectrum
        I2, # intensity for 2nd spectrum
        I3, # intensity for 3rd spectrum

        mu1, # central position, without coupling, shared
        mu2, # central position, without coupling, shared
        mu3, # central position, without coupling, shared

        a1, # baseline for 1st spectrum
        a2, # baseline for 2nd spectrum
        a3, # baseline for 3rd spectrum

        ligand_H_fractions, # list of ligand H-fractions for each sample
        solvent_H_fractions, # list of solvent H-fractions for each sample
        spectral_indices, # indices associated with the x-values for each spectrum. n+1, where n is the number of samples

        ):

    total_ys = np.array([])

    for I, mu, a, lH, sH, i in zip(
            [I1, I2, I3],
            [mu1, mu2, mu3],
            [a1, a2, a3],
            ligand_H_fractions,
            solvent_H_fractions,
            range(len(ligand_H_fractions)) # get indices for each sample (i.e., for each solvent fraction)
            ):


        part_x = x[spectral_indices[i]: spectral_indices[i+1]] # slice out the x-values associated with this part


        part_y_sim = simulate_pulse_signal(
                part_x,

                ligand_binned_centers, # to be supplied as not fit
                solvent_binned_centers, # to be supplied as not fit

                ligand_binned_counts, # to be supplied as not fit
                solvent_binned_counts, # to be supplied as not fit

                lH, # to be supplied as not fit
                sH, # to be supplied as not fit

                alpha, # decay constant
                A, # scalar value to adjust strength of coupling

                I, # intensity of individual contributions
                mu, # central position, without coupling
                gamma, # width of the lorentzian peak

                a, # baseline correction
                )

        total_ys = np.concatenate((total_ys, part_y_sim))
    return total_ys


    # handling through bond and through space

#need a fucntion that will handle through space and through bond together
def simulate_pulse_signal_two_exponential(
        x,

        ligand_binned_centers, # to be supplied as not fit
        solvent_binned_centers, # to be supplied as not fit

        ligand_binned_counts, # to be supplied as not fit
        solvent_binned_counts, # to be supplied as not fit

        ligand_H_fraction, # to be supplied as not fit
        solvent_H_fraction, # to be supplied as not fit

        alpha_through_space, # decay constant
        alpha_through_bond, # through bond decay rate
        f_tb, # this is the fraction of the wavefunction that goes through bond.
        A, # scalar value to adjust strength of coupling

        I, # intensity of individual contributions
        mu, # central position, without coupling
        gamma, # width of the lorentzian peak

        a, # baseline correction
        ):


    # for debudgin...
    #print(f"alpha={alpha:.4f}, A={A:.4f}, I={I:.4f}, mu={mu:.4f}, gamma={gamma:.4f}, a={a:.4f}")


    # calculate the psi-squared values at each distance for ligand and solvent
    ligand_psi_through_bond = np.exp(-alpha_through_bond * ligand_binned_centers) 
    ligand_psi_through_space = np.exp(-alpha_through_space * ligand_binned_centers) 
    solvent_psi_through_space = np.exp(-alpha_through_space * solvent_binned_centers)


    liagnd_psi_squareds_through_bond = np.exp(-2*alpha_through_bond * ligand_binned_centers)
    ligand_psi_squareds_through_space  = np.exp(-2*alpha_through_space * ligand_binned_centers)
    solvent_psi_squareds_through_space = np.exp(-2*alpha_through_space * solvent_binned_centers)

    # now we are going to sum up contributions for through space
    temp_y_values_through_space = np.zeros_like(x)
    for info in [
        [ligand_psi_through_space, ligand_psi_squareds_through_space, ligand_binned_counts, ligand_H_fraction], 
        [solvent_psi_through_space, solvent_psi_squareds_through_space, solvent_binned_counts, solvent_H_fraction]
    ]:
        psis, psi_squareds, counts, H = info # unpack these
        for p, ps, count in zip(psis, psi_squareds, counts):
            temp_y_values_through_space = temp_y_values_through_space + I*lorentzian(x, mu = mu +   ps*A, gamma = gamma) * H*ps * count
            temp_y_values_through_space = temp_y_values_through_space + I*lorentzian(x, mu = mu - 1*ps*A, gamma = gamma) * H*ps * count

    # now get the through bond added in
    temp_y_values_through_bond = np.zeros_like(x)
    for info in [
        [ligand_psi_through_bond, liagnd_psi_squareds_through_bond, ligand_binned_counts, ligand_H_fraction], 
    ]:
        psis, psi_squareds, counts, H = info # unpack these
        for p, ps, count in zip(psis, psi_squareds, counts):
            temp_y_values_through_bond = temp_y_values_through_bond + I*lorentzian(x, mu = mu +   ps*A, gamma = gamma) * H*ps * count
            temp_y_values_through_bond = temp_y_values_through_bond + I*lorentzian(x, mu = mu - 1*ps*A, gamma = gamma) * H*ps * count

    return (1-f_tb)*temp_y_values_through_space + f_tb*temp_y_values_through_bond + a

#function for model
def model_three_samples_two_exponential(
        x, # x values. These are the concatenated x-values from the experimental spectra we want to fit

        ligand_binned_centers, # to be supplied as not fit
        solvent_binned_centers, # to be supplied as not fit

        ligand_binned_counts, # to be supplied as not fit
        solvent_binned_counts, # to be supplied as not fit

        alpha_through_space, # decay constant
        alpha_through_bond,
        f_tb,
        A, # scalar value to adjust strength of coupling

        gamma, # width of the lorentzian peak, shared

        I1, # intensity for 1st spectrum
        I2, # intensity for 2nd spectrum
        I3, # intensity for 3rd spectrum

        mu1, # central position, without coupling, shared
        mu2, # central position, without coupling, shared
        mu3, # central position, without coupling, shared

        a1, # baseline for 1st spectrum
        a2, # baseline for 2nd spectrum
        a3, # baseline for 3rd spectrum

        ligand_H_fractions, # list of ligand H-fractions for each sample
        solvent_H_fractions, # list of solvent H-fractions for each sample
        spectral_indices, # indices associated with the x-values for each spectrum. n+1, where n is the number of samples

        ):

    total_ys = np.array([])

    for I, mu, a, lH, sH, i in zip(
            [I1, I2, I3],
            [mu1, mu2, mu3],
            [a1, a2, a3],
            ligand_H_fractions,
            solvent_H_fractions,
            range(len(ligand_H_fractions)) # get indices for each sample (i.e., for each solvent fraction)
            ):


        part_x = x[spectral_indices[i]: spectral_indices[i+1]] # slice out the x-values associated with this part


        part_y_sim = simulate_pulse_signal_two_exponential(
                part_x,

                ligand_binned_centers, # to be supplied as not fit
                solvent_binned_centers, # to be supplied as not fit

                ligand_binned_counts, # to be supplied as not fit
                solvent_binned_counts, # to be supplied as not fit

                lH, # to be supplied as not fit
                sH, # to be supplied as not fit

                alpha_through_space, # decay constant
                alpha_through_bond,
                f_tb,
                A, # scalar value to adjust strength of coupling

                I, # intensity of individual contributions
                mu, # central position, without coupling
                gamma, # width of the lorentzian peak

                a, # baseline correction
                )

        total_ys = np.concatenate((total_ys, part_y_sim))
    return total_ys


def plot_triple_fit(fit, indices):
    # this is just for plotting the fits

    labels = [
        "SC<sub>12</sub>H<sub>25</sub> and C<sub>7</sub>H<sub>8</sub>",
        "SC<sub>12</sub>D<sub>25</sub> and C<sub>7</sub>H<sub>8</sub>",
        "SC<sub>12</sub>H<sub>25</sub> and C<sub>7</sub>D<sub>8</sub>",
        ]

    # Just making some variables for convenience
    # First figure out what the independent variable name(s) is(are)
    independent_vars = fit.model.independent_vars

    # The x data has to be the same for all the independent variables, so
    # so get it from the first one in the list for safety
    xdata = fit.userkws[independent_vars[0]]
    ydata = fit.data


    fig = make_subplots(rows = 3, cols =1)
    for i in range(3):
        fig.add_scatter(
            x = xdata[indices[i]:indices[i+1]],
            y = ydata[indices[i]:indices[i+1]],
            mode = "lines",
            line = dict(color = "grey", width = 10),
            showlegend = False,
            row = i+1, col = 1,
        )
        fig.add_scatter(
            x = xdata[indices[i]:indices[i+1]],
            y = fit.best_fit[indices[i]:indices[i+1]],
             mode = "lines",
            line = dict(color = "white", width = 4),
            showlegend = False,
            row = i+1, col = 1,
        )
        fig.add_scatter(
            x = xdata[indices[i]:indices[i+1]],
            y = fit.best_fit[indices[i]:indices[i+1]],
            mode = "lines",
            line = dict(color = "red", width = 2),
            showlegend = False,
            row = i+1, col = 1,
        )

        fig.add_annotation(
            text = labels[i],
            x = min(xdata[indices[i]:indices[i+1]]),
            y = max(fit.best_fit[indices[i]:indices[i+1]]),
            xanchor  = "left",
            showarrow = False,
            row = i+1, col = 1,
        )


    fig.update_xaxes(title = "magnetic field")
    fig.update_yaxes(title = "intensity")
    fig.update_layout(
        template="simple_white",
        margin = dict(t=5, r = 5)
    )

    return fig



def calc_ligand_solvent_amounts(alpha, bin_centers, ligand_counts, solvent_counts):

    # get the extended solvent.
    extended_bins, extended_ligand, extended_solvent = extend_bins_based_on_wavefunction(bin_centers, ligand_counts, solvent_counts, alpha = 1, fract_threshold = 0.01)
    # then multiply the psi-squared values times the hydrogen bins and then sum
    ligand = np.sum(np.exp(-2*alpha*extended_bins)*extended_ligand)
    solvent = np.sum(np.exp(-2*alpha*extended_bins)*extended_solvent)

    return ligand, solvent

def find_ratio_of_ligand_to_solvent(
    alpha, 
    bin_centers, 
    ligand_counts, 
    solvent_counts
    ):

    ligand_sum, solvent_sum = calc_ligand_solvent_amounts(alpha, bin_centers, ligand_counts, solvent_counts)

    return [solvent_sum/ligand_sum]



def write_binned_ligand_solvent(centers, ligand, solvent, name = "binned ligand and solvent data"):
    with open(name, "w") as f:
        f.write("bin centers, ligand counts, solvent counts \n")
        for c, l, s in zip(centers, ligand, solvent):
            f.write(f"{c}, {l}, {s}" + "\n")

    print("file written")
    return None