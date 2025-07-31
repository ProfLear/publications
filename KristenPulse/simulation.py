# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 20:30:01 2025

@author: benle
"""

import numpy as np
from pathlib import Path
from plotly.subplots import make_subplots
import codechembook.quickPlots as qp
import os
from lmfit import Model


try:
    os.chdir(r"C:\Users\benle\Documents\GitHub\publications\KristenPulse")
except:
    pass


#%% Define functions for simulating and fitting spectra

def lorentzian(x, A = 1, mu = 1, gamma = 1):
    return (A / np.pi) * (gamma / ((x - mu)**2 + gamma**2))

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
    ligand_psi_squareds  = np.exp(-2*alpha * ligand_binned_centers)
    solvent_psi_squareds = np.exp(-2*alpha * solvent_binned_centers)
    
    # now we are going to sum up contributions
    temp_y_values = np.zeros_like(x)
    for info in [[ligand_psi_squareds, ligand_binned_counts, ligand_H_fraction], [solvent_psi_squareds, solvent_binned_counts, solvent_H_fraction]]:
        psi_squareds, counts, H = info # unpack these
        for ps, count in zip(psi_squareds, counts):
            temp_y_values = temp_y_values + I*lorentzian(x, mu = mu +   ps*A, gamma = gamma) * H*ps * count
            temp_y_values = temp_y_values + I*lorentzian(x, mu = mu - 1*ps*A, gamma = gamma) * H*ps * count
             
    return temp_y_values + a



#%% import both distance and spectral data

H25_H8_HD = {}
D25_H8_HD = {}
H25_D8_HD = {}

# import distance data
high_density_toluene_file = Path(r".\MDSimulationResults\toluene_distances_high_density.csv")
high_density_ligand_file = Path(r".\MDSimulationResults\dodecanol_distances_high_density.csv")
low_density_toluene_file = Path(r".\MDSimulationResults\toluene_distances_low_density.csv")
low_density_ligand_file = Path(r".\MDSimulationResults\dodecanol_distances_low_density.csv")

high_density_toluene_distanaces = np.genfromtxt(high_density_toluene_file, delimiter = ",")
high_density_ligand_distanaces = np.genfromtxt(high_density_ligand_file, delimiter = ",")
low_density_toluene_distanaces = np.genfromtxt(low_density_toluene_file, delimiter = ",")
low_density_ligand_distanaces = np.genfromtxt(low_density_ligand_file, delimiter = ",")

#% load the experimental data

H25_H8_HD["experimental spectrum"] = np.genfromtxt(r".\MIMS\PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)
D25_H8_HD["experimental spectrum"] = np.genfromtxt(r".\MIMS\dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter=",", unpack = True, skip_header = 1)
H25_D8_HD["experimental spectrum"] = np.genfromtxt(r".\MIMS\Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter=",", unpack = True, skip_header = 1)



#%%

high_density_ligand_distance_hist = qp.quickHist(x = high_density_ligand_distanaces, width = 1)
high_density_toluene_distance_hist= qp.quickHist(x = high_density_toluene_distanaces, width = 1)
low_density_ligand_distance_hist= qp.quickHist(x = low_density_ligand_distanaces, width = 1)
low_density_toluene_distance_hist= qp.quickHist(x = low_density_toluene_distanaces, width = 1)


#%% plot together
bothhist = make_subplots(rows =2, cols = 1)
bothhist.add_bar(
    x = high_density_ligand_distance_hist.data[0].x, 
    y = high_density_ligand_distance_hist.data[0].y,
    row = 1, col =1)
bothhist.add_bar(
    x = high_density_toluene_distance_hist.data[0].x, 
    y = high_density_toluene_distance_hist.data[0].y,
    row = 2, col =1)
for i in [1,2]:
    bothhist.add_scatter(x =[max(high_density_ligand_distanaces)]*2,
                         y = [0, 500], mode = "lines",
                         row = i, col = 1)
bothhist.update_xaxes(range = [0, 40])
bothhist.show("png")

#%% now simulate






H25_H8_HD = {
    "ligand"  : {"distances": high_density_ligand_distanaces,  "fraction hydrogen" : 1},
    "solvent" : {"distances": high_density_toluene_distanaces, "fraction hydrogen" : 1},
    "experimental spectrum" : None,
    }
D25_H8_HD = {
    "ligand"  : {"distances": high_density_ligand_distanaces,  "fraction hydrogen" : 0.012},
    "solvent" : {"distances": high_density_toluene_distanaces, "fraction hydrogen" : 1},
    "experimental spectrum" : None,
    }
H25_D8_HD = {
    "ligand"  : {"distances": high_density_ligand_distanaces,  "fraction hydrogen" : 1},
    "solvent" : {"distances": high_density_toluene_distanaces, "fraction hydrogen" : 0.01},
    "experimental spectrum" : None,
    }

#print(f"Type of H25_H8_HD['ligand']['fraction hydrogen']: {type(H25_H8_HD['ligand']['fraction hydrogen'])}, Value: {H25_H8_HD['ligand']['fraction hydrogen']}")
#print(f"Type of H25_H8_HD['solvent']['fraction hydrogen']: {type(H25_H8_HD['solvent']['fraction hydrogen'])}, Value: {H25_H8_HD['solvent']['fraction hydrogen']}")




qp.quickScatter(x = H25_H8_HD["experimental spectrum"][0], y = H25_H8_HD["experimental spectrum"][1])






#%% now fit individually



#%
# numbers after the "#" were for the non-psi^2 weigted intensities
dc = 0.0017#0.15
I =  5e2#4.e-5 * 0.75e7
mu = 49.95
cs = 2e-2#5
g = .15#0.2
a = -35.3e6

H25_D8_HD_simulated = simulate_pulse_signal(
    H25_D8_HD["ligand"]["distances"],
    H25_D8_HD["solvent"]["distances"],
    H25_D8_HD["ligand"]["fraction hydrogen"],
    H25_D8_HD["solvent"]["fraction hydrogen"],
    decay_constant = dc,
    I = I, 
    mu = mu, 
    coupling_scaler = cs,
    gamma = g,
    a = a,
    x = H25_D8_HD["experimental spectrum"][0]
    )

qp.quickScatter(
    x = [H25_D8_HD_simulated[0], H25_D8_HD["experimental spectrum"][0]], 
    y = [H25_D8_HD_simulated[1], H25_D8_HD["experimental spectrum"][1]],
    )

#%%




# Explicitly list the names of the arguments that are fitting parameters
fit_parameter_names = ['decay_constant', 'I', 'mu', 'coupling_scaler', 'gamma', 'a']

# Create the Model, telling it which arguments are parameters and which are independent
singleModel = Model(model_single_sample,
                    independent_vars=['x'],
                    param_names=fit_parameter_names)



singleModelParams = singleModel.make_params()
singleModelParams.add_many(
    ("decay_constant",  dc, True, 0, None),
    ("I",               I,  True, 0, None),
    ("mu",              mu, True, None, None), 
    ("coupling_scaler", cs, True, 0, None),
    ("gamma",           g,  False, 1e-9, None),
    ("a",               a,  True, None, None),
    )

#norm_y = (H25_H8_HD["experimental spectrum"][1]-a)/max(H25_H8_HD["experimental spectrum"][1]-a)
H25_D8_HD_result = singleModel.fit(
    H25_D8_HD["experimental spectrum"][1], 
    params = singleModelParams,
    x = H25_D8_HD["experimental spectrum"][0],
    ligand_distances   = H25_D8_HD["ligand"]["distances"], 
    solvent_distances  = H25_D8_HD["solvent"]["distances"],
    ligand_H_fraction  = H25_D8_HD["ligand"]["fraction hydrogen"],
    solvent_H_fraction = H25_D8_HD["solvent"]["fraction hydrogen"],
    )

print("for H25_D8_HD_result:")
print(H25_D8_HD_result.fit_report())
print("\n")
qp.plotFit(H25_D8_HD_result, residual = True)
#%%

# numbers after the "#" were for the non-psi^2 weigted intensities
dc = 0.0017#0.15
I =  2e2#4.e-5 * 0.75e7
mu = 49.95
cs = 2e-2#5
g = .15#0.2
a = -32.3e6

D25_H8_HD_simulated = simulate_pulse_signal(
    D25_H8_HD["ligand"]["distances"],
    D25_H8_HD["solvent"]["distances"],
    D25_H8_HD["ligand"]["fraction hydrogen"],
    D25_H8_HD["solvent"]["fraction hydrogen"],
    decay_constant = dc,
    I = I, 
    mu = mu, 
    coupling_scaler = cs,
    gamma = g,
    a = a,
    x = D25_H8_HD["experimental spectrum"][0]
    )

qp.quickScatter(
    x = [D25_H8_HD_simulated[0], D25_H8_HD["experimental spectrum"][0]], 
    y = [D25_H8_HD_simulated[1], D25_H8_HD["experimental spectrum"][1]],
    )
#%%

singleModelParams.add_many(
    ("decay_constant",  dc, True, 0, None),
    ("I",               I,  True, 0, None),
    ("mu",              mu, True, None, None), 
    ("coupling_scaler", cs, True, 0, None),
    ("gamma",           g,  True, 1e-9, None),
    ("a",               a,  True, None, None),
    )

D25_H8_HD_result = singleModel.fit(
    D25_H8_HD["experimental spectrum"][1], 
    params = singleModelParams,
    x = D25_H8_HD["experimental spectrum"][0],
    ligand_distances   = D25_H8_HD["ligand"]["distances"], 
    solvent_distances  = D25_H8_HD["solvent"]["distances"],
    ligand_H_fraction  = D25_H8_HD["ligand"]["fraction hydrogen"],
    solvent_H_fraction = D25_H8_HD["solvent"]["fraction hydrogen"],
    )

print("for D25_H8_HD_result:")
print(D25_H8_HD_result.fit_report())
print("\n")
qp.plotFit(D25_H8_HD_result, residual = True)

#%%

H25_H8_HD_result = singleModel.fit(
    H25_H8_HD["experimental spectrum"][1], 
    params = singleModelParams,
    x = H25_H8_HD["experimental spectrum"][0],
    ligand_distances   = H25_H8_HD["ligand"]["distances"], 
    solvent_distances  = H25_H8_HD["solvent"]["distances"],
    ligand_H_fraction  = H25_H8_HD["ligand"]["fraction hydrogen"],
    solvent_H_fraction = H25_H8_HD["solvent"]["fraction hydrogen"],
    )

print("for H25_D8_HD_result:")
print(H25_D8_HD_result.fit_report())
print("\n")
qp.plotFit(H25_D8_HD_result, residual = True)


#%% now, supply all the y-values in a single place, and fit that...







# Explicitly list the names of the arguments that are fitting parameters
triple_fit_parameter_names = [
                                'decay_constant', 'mu', 'coupling_scaler', 'gamma',
                                'I1', 'I2', 'I3',
                                'a1', 'a2', 'a3',
                              ]

# Create the Model, telling it which arguments are parameters and which are independent
tripleModel = Model(model_three_samples,
                    independent_vars=['x'],
                    param_names=triple_fit_parameter_names)



tripleModelParams = tripleModel.make_params()
tripleModelParams.add_many(
    ("decay_constant",  dc, True, 0, None),
    ("mu",              mu, True, None, None), 
    ("coupling_scaler", cs, True, 0, None),
    ("gamma",           g,  True, 1e-9, None),
    ("I1", H25_H8_HD_result.params["I"].value,  True, 0, None),
    ("I2", D25_H8_HD_result.params["I"].value,  True, 0, None),
    ("I3", H25_D8_HD_result.params["I"].value,  True, 0, None),
    ("a1", H25_H8_HD_result.params["a"].value,  True, None, None),
    ("a2", D25_H8_HD_result.params["a"].value,  True, None, None),
    ("a3", H25_D8_HD_result.params["a"].value,  True, None, None),
    )

x_values_list = [H25_H8_HD["experimental spectrum"][0], D25_H8_HD["experimental spectrum"][0], H25_D8_HD["experimental spectrum"][0]]
y_values_list = [H25_H8_HD["experimental spectrum"][1], D25_H8_HD["experimental spectrum"][1], H25_D8_HD["experimental spectrum"][1]]

trimmed_x_values_list = []
trimmed_y_values_list = []
xlimits = [-np.inf, np.inf]
for x, y in zip(x_values_list, y_values_list):
    mask = (x > xlimits[0]) & (x < xlimits[1])
    trimmed_x_values_list.append(x[mask])
    trimmed_y_values_list.append(y[mask])
x_values = np.concatenate(trimmed_x_values_list)
y_values = np.concatenate(trimmed_y_values_list)

tripleFitResult = tripleModel.fit(
    y_values, 
    params = tripleModelParams,
    x = x_values,
    ligand_distances   = H25_H8_HD["ligand"]["distances"], 
    solvent_distances  = H25_H8_HD["solvent"]["distances"],
    ligand_H_fractions = [1, 0.012, 1],
    solvent_H_fractions = [1, 1, 0.01],
    spectral_indices = [
        0, 
        len(trimmed_x_values_list[0]),
        len(trimmed_x_values_list[0]) + len(trimmed_x_values_list[1]),
        len(trimmed_x_values_list[0]) + len(trimmed_x_values_list[1]) + len(trimmed_x_values_list[2])
                         ],
    )

print(tripleFitResult.fit_report())
#qp.plotFit(tripleFitResult, residual = True)

#%% work up the triple fit...

complete_result = make_subplots(rows = 3, cols = 1,
                                subplot_titles = [
                                    "all hydrogen",
                                    "dueterated ligand",
                                    "dueterated solvent",
                                    ])

experimental_xs = [H25_H8_HD["experimental spectrum"][0], D25_H8_HD["experimental spectrum"][0], H25_D8_HD["experimental spectrum"][0]]
experimental_ys = [H25_H8_HD["experimental spectrum"][1], D25_H8_HD["experimental spectrum"][1], H25_D8_HD["experimental spectrum"][1]]
ligand_H_fractions = [1, 0.012, 1]
solvent_H_fractions = [1, 1, 0.01]
annotations = [
    "all hydrogen",
    "dueterated ligand",
    "dueterated solvent",
    ]

'''
ligand_distances, 
solvent_distances,

ligand_H_fraction = 1,
solvent_H_fraction = 1,

decay_constant = 1,
coupling_scaler = 1,

I = 1,
mu = 0,
gamma = 1, 
a = 0, 
x = None,
'''

for i in range(3):
    sim_x, sim_y = simulate_pulse_signal(
        ligand_distances = H25_H8_HD["ligand"]["distances"], 
        solvent_distances = H25_H8_HD["ligand"]["distances"],
        
        ligand_H_fraction = ligand_H_fractions[i],
        solvent_H_fraction = solvent_H_fractions[i], 

        decay_constant = tripleFitResult.params["decay_constant"].value, 
        coupling_scaler = tripleFitResult.params["coupling_scaler"].value,
        
        I = tripleFitResult.params[f"I{i+1}"].value,
        mu = tripleFitResult.params["mu"].value, 
        gamma = tripleFitResult.params["gamma"].value, 

        a = 0,
        x = trimmed_x_values_list[i],

        )
    
    
    complete_result.add_scatter( # data
                                x = trimmed_x_values_list[i],
                                y = trimmed_y_values_list[i] - tripleFitResult.params[f"a{i+1}"].value,
                                line = dict(width = 8, color = "grey"),
                                showlegend = False,
                                row=i+1, col = 1,
        )
    
    complete_result.add_scatter( # fit
                                x = sim_x,
                                y = sim_y,
                                line = dict(width = 2, color = "white"),
                                showlegend = False,
                                row=i+1, col = 1,
        )
    complete_result.add_scatter( # fit
                                x = sim_x,
                                y = sim_y,
                                line = dict(width = 1, color = "red"),
                                showlegend = False,
                                row=i+1, col = 1,
        )
#%
complete_result.update_layout(template = "simple_white",
                              width = 1.5*300, height = 1.5*300,
                              margin = dict(l = 50, r = 5, t = 20, b = 20))
complete_result.show('png')
    


#%% 

print(f"Type of dc: {type(dc)}, Value: {dc}")
print(f"Type of I: {type(I)}, Value: {I}")
print(f"Type of mu: {type(mu)}, Value: {mu}")
print(f"Type of cs: {type(cs)}, Value: {cs}")
print(f"Type of g: {type(g)}, Value: {g}")
print(f"Type of a (initial value): {type(0)}, Value: {0}") # 'a' is initialized with 0, so it's not 'a'

singleModelParams.add_many(
    ("decay_constant",  dc, True, 0, None),
    ("I",               I,  True, 0, None),
    ("mu",              mu, True, None, None),
    ("coupling_scaler", cs, True, 0, None),
    ("gamma",           g,  True, 0, None),
    ("a",               0,  True, None, None),
)

#%%

simulate_pulse_signal(
    H25_H8_HD,
    x = H25_H8_HD["experimental spectrum"][0],
    decay_constant = dc,
    I =I,
    mu = mu, 
    coupling_scaler = cs,
    gamma = g,
    a = 0
    )



#%% This set of functions is working, but needs to be re-written to be better...

def simulate_pulse_signal(experiment_dictionary, decay_constant = 1, gamma = 1, coupling_scaler = 1, a = 0, mu = 0, I = 1, x = None, ):
     
    experiment_dictionary["ligand"]["couplings"] = calcluate_psi_squareds(experiment_dictionary["ligand"]["distances"], decay_constant = decay_constant)
    experiment_dictionary["solvent"]["couplings"] = calcluate_psi_squareds(experiment_dictionary["solvent"]["distances"], decay_constant = decay_constant)
   
    max_x = max([np.max(experiment_dictionary["ligand"]["couplings"]), np.max(experiment_dictionary["solvent"]["couplings"])]) + 6*gamma
    if x is not None:
        x_values = x.copy()
    else:
        x_values = np.linspace(mu-max_x, mu+max_x, 1000)
    
    for source in ["ligand", "solvent"]:
        temp_y_values = np.zeros_like(x_values)
        #print(source)
        #print(experiment[source]["couplings"])
        for p in experiment_dictionary[source]["couplings"]:
            temp_y_values = temp_y_values + I*lorentzian(x_values, mu = mu + p*coupling_scaler, gamma = gamma)
            temp_y_values = temp_y_values + I*lorentzian(x_values, mu = mu - 1*p*coupling_scaler, gamma = gamma)
        
        experiment_dictionary[source]["spectrum contribution"] = temp_y_values.copy()*experiment_dictionary[source]["fraction hydrogen"]
    experiment_dictionary["x values"] = x_values   
    experiment_dictionary['total spectrum'] = experiment_dictionary["ligand"]["spectrum contribution"] + experiment_dictionary["solvent"]["spectrum contribution"] + a
    
    return experiment_dictionary['total spectrum'] # return just the y-values, for fitting

simulate_pulse_signal(H25_H8_HD, x = H25_H8_HD["experimental spectrum"][0], decay_constant = dc, mu = mu, I = I, gamma = g, coupling_scaler = cs,)
simulate_pulse_signal(D25_H8_HD, x = H25_H8_HD["experimental spectrum"][0], decay_constant = dc, mu = mu, I = I, gamma = g, coupling_scaler = cs,)
simulate_pulse_signal(H25_D8_HD, x = H25_H8_HD["experimental spectrum"][0], decay_constant = dc, mu = mu, I = I, gamma = g, coupling_scaler = cs,)

total_sim = qp.quickScatter(
    x = [H25_H8_HD["x values"],       H25_D8_HD["x values"],       D25_H8_HD["x values"]], 
    y = [H25_H8_HD['total spectrum'], H25_D8_HD['total spectrum'], D25_H8_HD['total spectrum']], 
    name = ["all", "ligand", "solvent"],
    output = None)  

total_sim.update_xaxes(title = "coupling value", )
total_sim.update_yaxes(title = "relative intensity")
total_sim.update_layout(margin = dict(l = 0, r=0, t = 0, b = 0))
total_sim.show('png')

comp_sim = qp.quickScatter(
    x = [H25_H8_HD["x values"],       H25_H8_HD["experimental spectrum"][0]],
    y = [H25_H8_HD['total spectrum'], (H25_H8_HD["experimental spectrum"][1]-a)/max(H25_H8_HD["experimental spectrum"][1]-a)]
    )



#%% going to make an attempt at fitting using histogram weights.
# numbers after the "#" were for the non-psi^2 weigted intensities
alpha = 0.02
I =  1e2#4.e-5 * 0.75e7
mu = 49.95
A = 1 
gamma = .25#0.2
a = -46.8e6

#%

fract_threshold = 0.01
bin_width = 0.9
#need to set up a common axis, and bin data for it. 
left_edge = 0
bin_centers = []
binned_sim_solvent = []
binned_sim_ligand = []

# first, bin the data we have
while left_edge < max(max(high_density_ligand_distanaces), max(high_density_toluene_distanaces)):
    bin_centers.append(left_edge + bin_width/2)
    binned_sim_ligand.append(np.sum((high_density_ligand_distanaces >= left_edge) & (high_density_ligand_distanaces < left_edge + bin_width)))
    binned_sim_solvent.append(np.sum((high_density_toluene_distanaces >= left_edge) & (high_density_toluene_distanaces < left_edge + bin_width)))
    left_edge = left_edge + bin_width
    
# then, we will need to add more bins, in order to account for the missing solvent.  
# but we need to determine how much of the "new" data to have..
# that wil be based on the wavefunction decay
extended_sim_solvent = []
extended_sim_ligand = []
extended_calc_solvent = []
extended_bin_centers = []
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
            extended_sim_ligand.append(binned_sim_ligand[i])
        except:
            extended_sim_ligand.append(0)
        try:
            extended_sim_solvent.append(binned_sim_solvent[i])
        except:
            extended_sim_solvent.append(0)
        
        
    else: # we are not yet extending the data
        extended_bin_centers.append(bin_centers[i])
        extended_calc_solvent.append(binned_sim_solvent[i])
        extended_sim_solvent.append(binned_sim_solvent[i])
        extended_sim_ligand.append(binned_sim_ligand[i])
        
        if binned_sim_ligand[i] > 0:
            found_ligand = True
        
        
        if i == len(bin_centers)-1: # we are are the end of the simulated data, and we MUST switch
            extend = True
        elif binned_sim_ligand[i] == 0 and found_ligand: # we have found the end of the ligand, and so we should follow the expected r**2 behavior now
            extend = True
        elif binned_sim_solvent[i] == max(binned_sim_solvent): # we have found the maximum solvent value, which should then just increase as r**2
            extend = True
    
    # now, using the decay function, figure out the fraction that this latest value represents...
    
    ligand_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_sim_ligand))
    solvent_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_calc_solvent))
    
    if ligand_contribution + solvent_contribution > 0 and extend: # use of extend ensure we have found the end of the ligand?
        latest_fract = np.exp(-2*alpha * extended_bin_centers[-1]) * extended_calc_solvent[-1] / (ligand_contribution + solvent_contribution)
    else:
        latest_fract = 1
    i = i + 1
#%
extended_bin_centers = np.array(extended_bin_centers)
extended_sim_ligand = np.array(extended_sim_ligand)
extended_sim_solvent = np.array(extended_sim_solvent)
extended_calc_solvent = np.array(extended_calc_solvent)


#% calculate contributions...
ligand_contributions = np.exp(-2*alpha * extended_bin_centers) * extended_sim_ligand
solvent_contributions = np.exp(-2*alpha * extended_bin_centers) * extended_calc_solvent

total_contribution = sum(ligand_contributions) + sum(solvent_contributions)


report_hist = make_subplots(rows = 2, cols = 1)
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_sim_solvent,
    marker = dict(color = "darkgreen", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_calc_solvent,
    marker = dict(color = "green", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_sim_ligand,
    marker = dict(color = "orange", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_scatter(
    x = extended_bin_centers, 
    y = np.exp(-alpha * extended_bin_centers) * max([max(extended_calc_solvent), max(extended_sim_ligand)]),
    mode = "lines",
    line = dict(color = "red",),
    showlegend = False, 
    row = 1, col = 1,
    )
report_hist.add_scatter(
    x = extended_bin_centers, 
    y = np.exp(-2*alpha * extended_bin_centers) * max([max(extended_calc_solvent), max(extended_sim_ligand)]),
    mode = "lines",
    line = dict(color = "red", dash = "dot"),
    showlegend = False, 
    row = 1, col = 1,
    )
report_hist.update_yaxes(
    title = "total hydrogens",
    range = [0, max([max(extended_calc_solvent), max(extended_sim_ligand)])],
    row = 1, col = 1,
    )


report_hist.add_bar(
    x = extended_bin_centers,
    y = solvent_contributions/total_contribution, # + ligand_contributions,
    marker = dict(color = "green", opacity = 0.66),
    showlegend = False,
    row = 2, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = ligand_contributions/total_contribution,
    marker = dict(color = "orange", opacity = 0.66),
    showlegend = False,
    row = 2, col = 1,
    )
report_hist.update_yaxes(
    title = "fraction contribution",
    row = 2, col = 1,
    )

report_hist.update_xaxes(
    title = "distance from closest Pd atom /Å",
    range = [0, extended_bin_centers[-1] + bin_width/2],
    )
report_hist.update_layout(
    template = "simple_white",
    barmode = "overlay",
    bargap = 0,
    width = 500, height = 400,
    margin = dict(l = 50, r = 10, t = 10, b = 25)
    )
report_hist.show("png")


#% Now try to simulate the data
sim_y = simulate_pulse_signal(
        H25_H8_HD["experimental spectrum"][0],
        
        extended_bin_centers, 
        extended_bin_centers,
        
        extended_sim_ligand,
        extended_calc_solvent,
        
        ligand_H_fraction = 1,
        solvent_H_fraction = 1,
        
        alpha = alpha,
        A = A,
        I = I,
        mu = mu,
        
        gamma = gamma, 
        a = a, 
        )


H25_H8_HD["simulated spectrum"] = [H25_H8_HD["experimental spectrum"][0], sim_y]

#%
spectra_sim_plot=make_subplots(rows = 3, cols = 1)

spectra_sim_plot.add_scatter(
    x = H25_H8_HD["experimental spectrum"][0],
    y = H25_H8_HD["experimental spectrum"][1],
    mode = "lines",
    line = dict(width = 6, color = "lightgrey"),
    showlegend = False,
    row = 1, col = 1
    )
spectra_sim_plot.add_scatter(
    x = H25_H8_HD["simulated spectrum"][0],
    y = H25_H8_HD["simulated spectrum"][1],
    mode = "lines",
    line = dict(width = 2, color = "red"),
    showlegend = False,
    row = 1, col = 1
    )


spectra_sim_plot.update_layout(
    template = "simple_white",
    )

spectra_sim_plot.show("png")


#% then try to fit the all hydrogenated sample

# Explicitly list the names of the arguments that are fitting parameters
single_sample_param_names = [ # these are the ones that we want to adjust...              
                "alpha",
                "A",
                
                "I",
                "mu",
                "gamma", 
                
                "a", 
                ]

single_sample_model = Model(
    simulate_pulse_signal,
    independent_vars=['x'],
    param_names=single_sample_param_names,
    )

single_sample_model_params = single_sample_model.make_params()

single_sample_model_params.add_many(
    ("alpha",   alpha,  True, 0, None),
    ("A",       A,      True, 0, None),
    
    ("I",       I,      True, 0, None),
    ("mu",      mu,     True, None, None),
    ("gamma",   gamma,  True, 0, None),
    
    ("a",       0,      True, None, None),
    )

single_sample_model_result = single_sample_model.fit(
    H25_H8_HD["experimental spectrum"][1],
    
    x = H25_H8_HD["experimental spectrum"][0],
    
    params = single_sample_model_params,
    
    ligand_binned_centers = extended_bin_centers, # to be supplied as not fit
    solvent_binned_centers = extended_bin_centers, # to be supplied as not fit
    
    ligand_binned_counts = extended_sim_ligand, # to be supplied as not fit
    solvent_binned_counts = extended_calc_solvent, # to be supplied as not fit
    
    ligand_H_fraction = 1, # to be supplied as not fit
    solvent_H_fraction = 1, # to be supplied as not fit

    )

print(single_sample_model_result.fit_report())


spectra_fit_plot=make_subplots(rows = 1, cols = 1)

spectra_fit_plot.add_scatter(
    x = H25_H8_HD["experimental spectrum"][0],
    y = H25_H8_HD["experimental spectrum"][1],
    mode = "lines",
    line = dict(width = 8, color = "grey"),
    showlegend = False,
    row = 1, col = 1
    )
spectra_fit_plot.add_scatter(
    x = H25_H8_HD["experimental spectrum"][0],
    y = single_sample_model_result.best_fit,
    mode = "lines",
    line = dict(width = 4, color = "white"),
    showlegend = False,
    row = 1, col = 1
    )
spectra_fit_plot.add_scatter(
    x = H25_H8_HD["experimental spectrum"][0],
    y = single_sample_model_result.best_fit,
    mode = "lines",
    line = dict(width = 2, color = "red"),
    showlegend = False,
    row = 1, col = 1
    )


spectra_fit_plot.update_layout(
    template = "simple_white",
    width = 500, height = 300,
    margin = dict(l = 50, b = 20, r = 10, t = 10)
    )

spectra_fit_plot.show("png")


#%% Now let us fit three at once

I1 = I
I2 = I
I3 = I

mu1 = 50
mu2 = 50
mu3 = 50

a1 = H25_H8_HD["experimental spectrum"][1][-1]
a2 = D25_H8_HD["experimental spectrum"][1][-1]
a3 = H25_D8_HD["experimental spectrum"][1][-1]

def model_three_samples(
        x,
        
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
                        
        ligand_H_fractions,
        solvent_H_fractions,
        spectral_indices,
        
        ):
    
    total_ys = np.array([])
    
    for I, mu, a, lH, sH, i in zip(
            [I1, I2, I3],
            [mu1, mu2, mu3],
            [a1, a2, a3],
            ligand_H_fractions,
            solvent_H_fractions,
            range(len(ligand_H_fractions))
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



#% then simulate them first for guesses

# H25_H8, D25_H8, H25_D8
concatenated_xs = list(H25_H8_HD["experimental spectrum"][0]) + list(D25_H8_HD["experimental spectrum"][0]) + list(H25_D8_HD["experimental spectrum"][0])
concatenated_ys = list(H25_H8_HD["experimental spectrum"][1]) + list(D25_H8_HD["experimental spectrum"][1]) + list(H25_D8_HD["experimental spectrum"][1])

ligand_H_fractions = [1, 0, 1]
solvent_H_fractions = [1, 1, 0]
spectral_indices = [
    0,
    len(H25_H8_HD["experimental spectrum"][0]),
    len(H25_H8_HD["experimental spectrum"][0]) + len(D25_H8_HD["experimental spectrum"][0]), 
    len(H25_H8_HD["experimental spectrum"][0]) + len(D25_H8_HD["experimental spectrum"][0]) + len(H25_D8_HD["experimental spectrum"][0]), 
    ]

simulated_three_samples = model_three_samples(
        concatenated_xs,
        
        extended_bin_centers, # to be supplied as not fit
        extended_bin_centers, # to be supplied as not fit
        
        extended_sim_ligand, # to be supplied as not fit
        extended_calc_solvent, # to be supplied as not fit
        
        alpha, # decay constant
        A, # scalar value to adjust strength of coupling
        
        gamma, # width of the lorentzian peak, shared
                
        I1, # intensity for 1st spectrum
        I2, # intensity for 2nd spectrum
        I3, # intensity for 3rd spectrum
        
        mu1, 
        mu2,
        mu3,
        
        a1, # baseline for 1st spectrum
        a2, # baseline for 2nd spectrum
        a3, # baseline for 3rd spectrum
                        
        ligand_H_fractions,
        solvent_H_fractions,
        spectral_indices,      
        )

triple_sim_plot = make_subplots(rows = 3, cols =1)

for i in range(len(ligand_H_fractions)):
    triple_sim_plot.add_scatter(
        x = concatenated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = concatenated_ys[spectral_indices[i]:spectral_indices[i+1]],
        mode = "lines",
        line = dict(width = 8, color = "grey"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    triple_sim_plot.add_scatter(
        x = concatenated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = simulated_three_samples[spectral_indices[i]:spectral_indices[i+1]],
        line = dict(width = 4, color = "white"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    triple_sim_plot.add_scatter(
        x = concatenated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = simulated_three_samples[spectral_indices[i]:spectral_indices[i+1]],
        line = dict(width = 2, color = "red"),
        row = i + 1, col = 1,
        showlegend = False,
        )

triple_sim_plot.update_layout(template = "simple_white")
triple_sim_plot.show("png")
    
    

#%% high density ligand coverage. 

# Explicitly list the names of the arguments that are fitting parameters
triple_sample_param_names = [ # these are the ones that we want to adjust...              
                "alpha",
                "A",
                
                "gamma", 
                
                "I1",
                "I2",
                "I3",
                
                "mu1",
                "mu2",
                "mu3",
                
                "a1",
                "a2",
                "a3",
                ]

triple_sample_model = Model(
    model_three_samples,
    independent_vars=['x'],
    param_names=triple_sample_param_names,
    )






triple_sample_model_params = triple_sample_model.make_params()

triple_sample_model_params.add_many(
    ("alpha",   alpha,  True, 0, None),
    ("A",       A,      True, 0, None),
    
    ("gamma",   gamma,  True, 0, None),
    
    ("I1",       I1,      True, 0, None),
    ("I2",       I2,      True, 0, None),
    ("I3",       I3,      True, 0, None),
    
    ("mu1",      mu,     True, 48, 52),
    ("mu2",      mu,     True, 48, 52),
    ("mu3",      mu,     True, 48, 52),

    ("a1",       a1,      True, None, None),
    ("a2",       a2,      True, None, None),
    ("a3",       a3,      True, None, None),

    )

triple_sample_model_result = triple_sample_model.fit(
    concatenated_ys,
    
    x = concatenated_xs,
    
    params = triple_sample_model_params,
    
    ligand_binned_centers = extended_bin_centers, # to be supplied as not fit
    solvent_binned_centers = extended_bin_centers, # to be supplied as not fit
    
    ligand_binned_counts = extended_sim_ligand, # to be supplied as not fit
    solvent_binned_counts = extended_calc_solvent, # to be supplied as not fit
    
    ligand_H_fractions =  [1, 0, 1], # to be supplied as not fit
    solvent_H_fractions = [1, 1, 0], # to be supplied as not fit
    
    spectral_indices = spectral_indices,
    )

print(triple_sample_model_result.fit_report())

#% plot
names = [
    "all hydrogenated",
    "only hydrogenated solvent",
    "only hydrogenated ligand",
    ]

triple_fit_plot = make_subplots(rows = 3, cols =1)

for i in range(len(ligand_H_fractions)):
    triple_fit_plot.add_scatter(
        x = concatenated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = concatenated_ys[spectral_indices[i]:spectral_indices[i+1]],
        mode = "lines",
        line = dict(width = 8, color = "grey"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    triple_fit_plot.add_scatter(
        x = concatenated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]],
        line = dict(width = 4, color = "white"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    triple_fit_plot.add_scatter(
        x = concatenated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]],
        line = dict(width = 2, color = "red"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    triple_fit_plot.add_annotation(
        text = names[i],
        y = max([
            max(concatenated_ys[spectral_indices[i]:spectral_indices[i+1]]),
            max(triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]])
            ]) + (max([
                max(concatenated_ys[spectral_indices[i]:spectral_indices[i+1]]),
                max(triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]])
                ]) - min([
                    min(concatenated_ys[spectral_indices[i]:spectral_indices[i+1]]),
                    min(triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]])
                    ])) * 0.15,
        showarrow = False,
        row =i + 1, col = 1
        )
triple_fit_plot.update_xaxes(    
    title = "magnetic field",
    row = 3, col = 1)

triple_fit_plot.update_xaxes(
    range = [46, 54],
    )

triple_fit_plot.update_layout(
    template = "simple_white",
    margin = dict(l = 50, r = 10, b = 20, t = 10),
    width = 500, height = 600,
    )
triple_fit_plot.show("png")



#%% Finally, histograms and psi squared based onthe fit...

alpha = triple_sample_model_result.params["alpha"].value
A = triple_sample_model_result.params["A"].value

fract_threshold = 0.01
bin_width = 0.9
#need to set up a common axis, and bin data for it. 
left_edge = 0
bin_centers = []
binned_sim_solvent = []
binned_sim_ligand = []

# first, bin the data we have
while left_edge < max(max(high_density_ligand_distanaces), max(high_density_toluene_distanaces)):
    bin_centers.append(left_edge + bin_width/2)
    binned_sim_ligand.append(np.sum((high_density_ligand_distanaces >= left_edge) & (high_density_ligand_distanaces < left_edge + bin_width)))
    binned_sim_solvent.append(np.sum((high_density_toluene_distanaces >= left_edge) & (high_density_toluene_distanaces < left_edge + bin_width)))
    left_edge = left_edge + bin_width
    
# then, we will need to add more bins, in order to account for the missing solvent.  
# but we need to determine how much of the "new" data to have..
# that wil be based on the wavefunction decay
extended_sim_solvent = []
extended_sim_ligand = []
extended_calc_solvent = []
extended_bin_centers = []
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
            extended_sim_ligand.append(binned_sim_ligand[i])
        except:
            extended_sim_ligand.append(0)
        try:
            extended_sim_solvent.append(binned_sim_solvent[i])
        except:
            extended_sim_solvent.append(0)
        
        
    else: # we are not yet extending the data
        extended_bin_centers.append(bin_centers[i])
        extended_calc_solvent.append(binned_sim_solvent[i])
        extended_sim_solvent.append(binned_sim_solvent[i])
        extended_sim_ligand.append(binned_sim_ligand[i])
        
        if binned_sim_ligand[i] > 0:
            found_ligand = True
        
        
        if i == len(bin_centers)-1: # we are are the end of the simulated data, and we MUST switch
            extend = True
        elif binned_sim_ligand[i] == 0 and found_ligand: # we have found the end of the ligand, and so we should follow the expected r**2 behavior now
            extend = True
        elif binned_sim_solvent[i] == max(binned_sim_solvent): # we have found the maximum solvent value, which should then just increase as r**2
            extend = True
    
    # now, using the decay function, figure out the fraction that this latest value represents...
    
    ligand_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_sim_ligand))
    solvent_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_calc_solvent))
    
    if ligand_contribution + solvent_contribution > 0 and extend: # use of extend ensure we have found the end of the ligand?
        latest_fract = np.exp(-2*alpha * extended_bin_centers[-1]) * extended_calc_solvent[-1] / (ligand_contribution + solvent_contribution)
    else:
        latest_fract = 1
    i = i + 1
#%
extended_bin_centers = np.array(extended_bin_centers)
extended_sim_ligand = np.array(extended_sim_ligand)
extended_sim_solvent = np.array(extended_sim_solvent)
extended_calc_solvent = np.array(extended_calc_solvent)


#% calculate contributions...
ligand_contributions = np.exp(-2*alpha * extended_bin_centers) * extended_sim_ligand
solvent_contributions = np.exp(-2*alpha * extended_bin_centers) * extended_calc_solvent

total_contribution = sum(ligand_contributions) + sum(solvent_contributions)


report_hist = make_subplots(rows = 2, cols = 1)
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_sim_solvent,
    marker = dict(color = "darkgreen", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_calc_solvent,
    marker = dict(color = "green", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_sim_ligand,
    marker = dict(color = "orange", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_scatter(
    x = extended_bin_centers, 
    y = np.exp(-alpha * extended_bin_centers) * max([max(extended_calc_solvent), max(extended_sim_ligand)]),
    mode = "lines",
    line = dict(color = "red",),
    showlegend = False, 
    row = 1, col = 1,
    )
report_hist.add_scatter(
    x = extended_bin_centers, 
    y = np.exp(-2*alpha * extended_bin_centers) * max([max(extended_calc_solvent), max(extended_sim_ligand)]),
    mode = "lines",
    line = dict(color = "red", dash = "dot"),
    showlegend = False, 
    row = 1, col = 1,
    )
report_hist.update_yaxes(
    title = "total hydrogens",
    range = [0, max([max(extended_calc_solvent), max(extended_sim_ligand)])],
    row = 1, col = 1,
    )


report_hist.add_bar(
    x = extended_bin_centers,
    y = solvent_contributions/total_contribution, # + ligand_contributions,
    marker = dict(color = "green", opacity = 0.66),
    showlegend = False,
    row = 2, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = ligand_contributions/total_contribution,
    marker = dict(color = "orange", opacity = 0.66),
    showlegend = False,
    row = 2, col = 1,
    )
report_hist.update_yaxes(
    title = "fraction contribution",
    row = 2, col = 1,
    )

report_hist.update_xaxes(
    title = "distance from closest Pd atom /Å",
    range = [0, extended_bin_centers[-1] + bin_width/2],
    )
report_hist.update_layout(
    template = "simple_white",
    barmode = "overlay",
    bargap = 0,
    width = 500, height = 400,
    margin = dict(l = 50, r = 10, t = 10, b = 25)
    )
report_hist.show("png")





#%%

#%% low density ligand coverage. 
# first, set up the bin centers and the bin counts...

bin_centers =[]
binned_sim_ligand = []
binned_sim_solvent = []
left_edge = 0
# first, bin the data we have
while left_edge < max(max(low_density_ligand_distanaces), max(low_density_toluene_distanaces)):
    bin_centers.append(left_edge + bin_width/2)
    binned_sim_ligand.append(np.sum((low_density_ligand_distanaces >= left_edge) & (low_density_ligand_distanaces < left_edge + bin_width)))
    binned_sim_solvent.append(np.sum((low_density_toluene_distanaces >= left_edge) & (low_density_toluene_distanaces < left_edge + bin_width)))
    left_edge = left_edge + bin_width

test = make_subplots()
test.add_bar(x = bin_centers, y = binned_sim_ligand)
test.show("png")
#%

# then, we will need to add more bins, in order to account for the missing solvent.  
# but we need to determine how much of the "new" data to have..
# that wil be based on the wavefunction decay
extended_sim_solvent = []
extended_sim_ligand = []
extended_calc_solvent = []
extended_bin_centers = []
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
            extended_sim_ligand.append(binned_sim_ligand[i])
        except:
            extended_sim_ligand.append(0)
        try:
            extended_sim_solvent.append(binned_sim_solvent[i])
        except:
            extended_sim_solvent.append(0)
        
        
    else: # we are not yet extending the data
        extended_bin_centers.append(bin_centers[i])
        extended_calc_solvent.append(binned_sim_solvent[i])
        extended_sim_solvent.append(binned_sim_solvent[i])
        extended_sim_ligand.append(binned_sim_ligand[i])
        
        if binned_sim_ligand[i] > 0:
            found_ligand = True
        
        
        if i == len(bin_centers)-1: # we are are the end of the simulated data, and we MUST switch
            extend = True
        elif binned_sim_ligand[i] == 0 and found_ligand: # we have found the end of the ligand, and so we should follow the expected r**2 behavior now
            extend = True
        elif binned_sim_solvent[i] == max(binned_sim_solvent): # we have found the maximum solvent value, which should then just increase as r**2
            extend = True
    
    # now, using the decay function, figure out the fraction that this latest value represents...
    
    ligand_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_sim_ligand))
    solvent_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_calc_solvent))
    
    if ligand_contribution + solvent_contribution > 0 and extend: # use of extend ensure we have found the end of the ligand?
        latest_fract = np.exp(-2*alpha * extended_bin_centers[-1]) * extended_calc_solvent[-1] / (ligand_contribution + solvent_contribution)
    else:
        latest_fract = 1
    i = i + 1
#%
extended_bin_centers = np.array(extended_bin_centers)
extended_sim_ligand = np.array(extended_sim_ligand)
extended_sim_solvent = np.array(extended_sim_solvent)
extended_calc_solvent = np.array(extended_calc_solvent)

#% calculate contributions...
ligand_contributions = np.exp(-2*alpha * extended_bin_centers) * extended_sim_ligand
solvent_contributions = np.exp(-2*alpha * extended_bin_centers) * extended_calc_solvent

total_contribution = sum(ligand_contributions) + sum(solvent_contributions)


report_hist = make_subplots(rows = 2, cols = 1)
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_sim_solvent,
    marker = dict(color = "darkgreen", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_calc_solvent,
    marker = dict(color = "green", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_sim_ligand,
    marker = dict(color = "orange", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_scatter(
    x = extended_bin_centers, 
    y = np.exp(-alpha * extended_bin_centers) * max([max(extended_calc_solvent), max(extended_sim_ligand)]),
    mode = "lines",
    line = dict(color = "red",),
    showlegend = False, 
    row = 1, col = 1,
    )
report_hist.add_scatter(
    x = extended_bin_centers, 
    y = np.exp(-2*alpha * extended_bin_centers) * max([max(extended_calc_solvent), max(extended_sim_ligand)]),
    mode = "lines",
    line = dict(color = "red", dash = "dot"),
    showlegend = False, 
    row = 1, col = 1,
    )
report_hist.update_yaxes(
    title = "total hydrogens",
    range = [0, max([max(extended_calc_solvent), max(extended_sim_ligand)])],
    row = 1, col = 1,
    )


report_hist.add_bar(
    x = extended_bin_centers,
    y = solvent_contributions/total_contribution, # + ligand_contributions,
    marker = dict(color = "green", opacity = 0.66),
    showlegend = False,
    row = 2, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = ligand_contributions/total_contribution,
    marker = dict(color = "orange", opacity = 0.66),
    showlegend = False,
    row = 2, col = 1,
    )
report_hist.update_yaxes(
    title = "fraction contribution",
    row = 2, col = 1,
    )

report_hist.update_xaxes(
    title = "distance from closest Pd atom /Å",
    range = [0, extended_bin_centers[-1] + bin_width/2],
    )
report_hist.update_layout(
    template = "simple_white",
    barmode = "overlay",
    bargap = 0,
    width = 500, height = 400,
    margin = dict(l = 50, r = 10, t = 10, b = 25)
    )
report_hist.show("png")


#%% then fit

triple_sample_model_result = triple_sample_model.fit(
    concatenated_ys,
    
    x = concatenated_xs,
    
    params = triple_sample_model_params,
    
    ligand_binned_centers = extended_bin_centers, # to be supplied as not fit
    solvent_binned_centers = extended_bin_centers, # to be supplied as not fit
    
    ligand_binned_counts = extended_sim_ligand, # to be supplied as not fit
    solvent_binned_counts = extended_calc_solvent, # to be supplied as not fit
    
    ligand_H_fractions =  [1, 0, 1], # to be supplied as not fit
    solvent_H_fractions = [1, 1, 0], # to be supplied as not fit
    
    spectral_indices = spectral_indices,
    )

print(triple_sample_model_result.fit_report())

#% plot
names = [
    "all hydrogenated",
    "only hydrogenated solvent",
    "only hydrogenated ligand",
    ]

triple_fit_plot = make_subplots(rows = 3, cols =1)

for i in range(len(ligand_H_fractions)):
    triple_fit_plot.add_scatter(
        x = concatenated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = concatenated_ys[spectral_indices[i]:spectral_indices[i+1]],
        mode = "lines",
        line = dict(width = 8, color = "grey"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    triple_fit_plot.add_scatter(
        x = concatenated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]],
        line = dict(width = 4, color = "white"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    triple_fit_plot.add_scatter(
        x = concatenated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]],
        line = dict(width = 2, color = "red"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    triple_fit_plot.add_annotation(
        text = names[i],
        y = max([
            max(concatenated_ys[spectral_indices[i]:spectral_indices[i+1]]),
            max(triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]])
            ]) + (max([
                max(concatenated_ys[spectral_indices[i]:spectral_indices[i+1]]),
                max(triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]])
                ]) - min([
                    min(concatenated_ys[spectral_indices[i]:spectral_indices[i+1]]),
                    min(triple_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]])
                    ])) * 0.15,
        showarrow = False,
        row =i + 1, col = 1
        )
triple_fit_plot.update_xaxes(    
    title = "magnetic field",
    row = 3, col = 1)

triple_fit_plot.update_xaxes(
    range = [46, 54],
    )

triple_fit_plot.update_layout(
    template = "simple_white",
    margin = dict(l = 50, r = 10, b = 20, t = 10),
    width = 500, height = 600,
    )
triple_fit_plot.show("png")


#%% Finally, histograms and psi squared based onthe fit...

alpha = triple_sample_model_result.params["alpha"].value
A = triple_sample_model_result.params["A"].value

fract_threshold = 0.01
bin_width = 0.9
#need to set up a common axis, and bin data for it. 
left_edge = 0
bin_centers = []
binned_sim_solvent = []
binned_sim_ligand = []

# first, bin the data we have
while left_edge < max(max(high_density_ligand_distanaces), max(high_density_toluene_distanaces)):
    bin_centers.append(left_edge + bin_width/2)
    binned_sim_ligand.append(np.sum((high_density_ligand_distanaces >= left_edge) & (high_density_ligand_distanaces < left_edge + bin_width)))
    binned_sim_solvent.append(np.sum((high_density_toluene_distanaces >= left_edge) & (high_density_toluene_distanaces < left_edge + bin_width)))
    left_edge = left_edge + bin_width
    
# then, we will need to add more bins, in order to account for the missing solvent.  
# but we need to determine how much of the "new" data to have..
# that wil be based on the wavefunction decay
extended_sim_solvent = []
extended_sim_ligand = []
extended_calc_solvent = []
extended_bin_centers = []
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
            extended_sim_ligand.append(binned_sim_ligand[i])
        except:
            extended_sim_ligand.append(0)
        try:
            extended_sim_solvent.append(binned_sim_solvent[i])
        except:
            extended_sim_solvent.append(0)
        
        
    else: # we are not yet extending the data
        extended_bin_centers.append(bin_centers[i])
        extended_calc_solvent.append(binned_sim_solvent[i])
        extended_sim_solvent.append(binned_sim_solvent[i])
        extended_sim_ligand.append(binned_sim_ligand[i])
        
        if binned_sim_ligand[i] > 0:
            found_ligand = True
        
        
        if i == len(bin_centers)-1: # we are are the end of the simulated data, and we MUST switch
            extend = True
        elif binned_sim_ligand[i] == 0 and found_ligand: # we have found the end of the ligand, and so we should follow the expected r**2 behavior now
            extend = True
        elif binned_sim_solvent[i] == max(binned_sim_solvent): # we have found the maximum solvent value, which should then just increase as r**2
            extend = True
    
    # now, using the decay function, figure out the fraction that this latest value represents...
    
    ligand_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_sim_ligand))
    solvent_contribution = sum(np.exp(-2*alpha * np.array(extended_bin_centers)) * np.array(extended_calc_solvent))
    
    if ligand_contribution + solvent_contribution > 0 and extend: # use of extend ensure we have found the end of the ligand?
        latest_fract = np.exp(-2*alpha * extended_bin_centers[-1]) * extended_calc_solvent[-1] / (ligand_contribution + solvent_contribution)
    else:
        latest_fract = 1
    i = i + 1
#%
extended_bin_centers = np.array(extended_bin_centers)
extended_sim_ligand = np.array(extended_sim_ligand)
extended_sim_solvent = np.array(extended_sim_solvent)
extended_calc_solvent = np.array(extended_calc_solvent)


#% calculate contributions...
ligand_contributions = np.exp(-2*alpha * extended_bin_centers) * extended_sim_ligand
solvent_contributions = np.exp(-2*alpha * extended_bin_centers) * extended_calc_solvent

total_contribution = sum(ligand_contributions) + sum(solvent_contributions)


report_hist = make_subplots(rows = 2, cols = 1)
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_sim_solvent,
    marker = dict(color = "darkgreen", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_calc_solvent,
    marker = dict(color = "green", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = extended_sim_ligand,
    marker = dict(color = "orange", opacity = 0.66),
    showlegend = False,
    row = 1, col = 1,
    )
report_hist.add_scatter(
    x = extended_bin_centers, 
    y = np.exp(-alpha * extended_bin_centers) * max([max(extended_calc_solvent), max(extended_sim_ligand)]),
    mode = "lines",
    line = dict(color = "red",),
    showlegend = False, 
    row = 1, col = 1,
    )
report_hist.add_scatter(
    x = extended_bin_centers, 
    y = np.exp(-2*alpha * extended_bin_centers) * max([max(extended_calc_solvent), max(extended_sim_ligand)]),
    mode = "lines",
    line = dict(color = "red", dash = "dot"),
    showlegend = False, 
    row = 1, col = 1,
    )
report_hist.update_yaxes(
    title = "total hydrogens",
    range = [0, max([max(extended_calc_solvent), max(extended_sim_ligand)])],
    row = 1, col = 1,
    )


report_hist.add_bar(
    x = extended_bin_centers,
    y = solvent_contributions/total_contribution, # + ligand_contributions,
    marker = dict(color = "green", opacity = 0.66),
    showlegend = False,
    row = 2, col = 1,
    )
report_hist.add_bar(
    x = extended_bin_centers,
    y = ligand_contributions/total_contribution,
    marker = dict(color = "orange", opacity = 0.66),
    showlegend = False,
    row = 2, col = 1,
    )
report_hist.update_yaxes(
    title = "fraction contribution",
    row = 2, col = 1,
    )

report_hist.update_xaxes(
    title = "distance from closest Pd atom /Å",
    range = [0, extended_bin_centers[-1] + bin_width/2],
    )
report_hist.update_layout(
    template = "simple_white",
    barmode = "overlay",
    bargap = 0,
    width = 500, height = 400,
    margin = dict(l = 50, r = 10, t = 10, b = 25)
    )
report_hist.show("png")


#%% Dual dueterated fitting
# here, we only fit the two dueterated samples, so that we don't overly bias towards the ligands. 

def model_two_samples(
        x,
        
        ligand_binned_centers, # to be supplied as not fit
        solvent_binned_centers, # to be supplied as not fit
        
        ligand_binned_counts, # to be supplied as not fit
        solvent_binned_counts, # to be supplied as not fit

        alpha, # decay constant
        A, # scalar value to adjust strength of coupling
        
        gamma, # width of the lorentzian peak, shared
                
        I1, # intensity for 1st spectrum
        I2, # intensity for 2nd spectrum
        
        mu1, # central position, without coupling, shared
        mu2, # central position, without coupling, shared
        
        a1, # baseline for 1st spectrum
        a2, # baseline for 2nd spectrum
                        
        ligand_H_fractions,
        solvent_H_fractions,
        spectral_indices,
        
        ):
    
    total_ys = np.array([])
    
    for I, mu, a, lH, sH, i in zip(
            [I1,  I2],
            [mu1, mu2],
            [a1,  a2],
            ligand_H_fractions,
            solvent_H_fractions,
            range(len(ligand_H_fractions))
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





# Explicitly list the names of the arguments that are fitting parameters
double_sample_param_names = [ # these are the ones that we want to adjust...              
                "alpha",
                "A",
                
                "gamma", 
                
                "I1",
                "I2",
                
                "mu1",
                "mu2",
                
                "a1",
                "a2",
                ]

double_sample_model = Model(
    model_two_samples,
    independent_vars=['x'],
    param_names=double_sample_param_names,
    )






double_sample_model_params = double_sample_model.make_params()

double_sample_model_params.add_many(
    ("alpha",   alpha,  True, 0, None),
    ("A",       A,      True, 0, None),
    
    ("gamma",   gamma,  True, 0, None),
    
    ("I1",       I2,      True, 0, None),
    ("I2",       I3,      True, 0, None),
    
    ("mu1",      mu,     True, None, None),
    ("mu2",      mu,     True, None, None),

    ("a1",       a2,      True, None, None),
    ("a2",       a3,      True, None, None),
    )



# build up the xs and ys



dueterated_xs = list(D25_H8_HD["experimental spectrum"][0]) + list(H25_D8_HD["experimental spectrum"][0])
dueterated_ys = list(D25_H8_HD["experimental spectrum"][1]) + list(H25_D8_HD["experimental spectrum"][1])
spectral_indices = [
    0,
    len(D25_H8_HD["experimental spectrum"][0]), 
    len(D25_H8_HD["experimental spectrum"][0]) + len(H25_D8_HD["experimental spectrum"][0]), 
    ]

double_sample_model_result = double_sample_model.fit(
    dueterated_ys,
    
    x = dueterated_xs,
    
    params = double_sample_model_params,
    
    ligand_binned_centers = extended_bin_centers, # to be supplied as not fit
    solvent_binned_centers = extended_bin_centers, # to be supplied as not fit
    
    ligand_binned_counts = extended_sim_ligand, # to be supplied as not fit
    solvent_binned_counts = extended_calc_solvent, # to be supplied as not fit
    
    ligand_H_fractions =  [0, 1], # to be supplied as not fit
    solvent_H_fractions = [1, 0], # to be supplied as not fit
    
    spectral_indices = spectral_indices,
    )

print(double_sample_model_result.fit_report())

#%% plot
names = [
    "only hydrogenated solvent",
    "only hydrogenated ligand",
    ]

double_fit_plot = make_subplots(rows = 2, cols =1)

for i in range(len([0,1])):
    double_fit_plot.add_scatter(
        x = dueterated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = dueterated_ys[spectral_indices[i]:spectral_indices[i+1]],
        mode = "lines",
        line = dict(width = 8, color = "grey"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    double_fit_plot.add_scatter(
        x = dueterated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = double_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]],
        line = dict(width = 4, color = "white"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    double_fit_plot.add_scatter(
        x = dueterated_xs[spectral_indices[i]:spectral_indices[i+1]],
        y = double_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]],
        line = dict(width = 2, color = "red"),
        row = i + 1, col = 1,
        showlegend = False,
        )
    double_fit_plot.add_annotation(
        text = names[i],
        y = max([
            max(concatenated_ys[spectral_indices[i]:spectral_indices[i+1]]),
            max(double_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]])
            ]) + (max([
                max(concatenated_ys[spectral_indices[i]:spectral_indices[i+1]]),
                max(double_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]])
                ]) - min([
                    min(concatenated_ys[spectral_indices[i]:spectral_indices[i+1]]),
                    min(double_sample_model_result.best_fit[spectral_indices[i]:spectral_indices[i+1]])
                    ])) * 0.15,
        showarrow = False,
        row =i + 1, col = 1
        )
double_fit_plot.update_xaxes(    
    title = "magnetic field",
    row = 3, col = 1)

double_fit_plot.update_xaxes(
    range = [46, 54],
    )

double_fit_plot.update_layout(
    template = "simple_white",
    margin = dict(l = 50, r = 10, b = 20, t = 10),
    width = 500, height = 600*2/3,
    )
double_fit_plot.show("png")
