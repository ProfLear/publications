# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 20:30:01 2025

@author: benle
"""

import numpy as np
from pathlib import Path
from plotly.subplots import make_subplots
import codechembook.quickPlots as qp



# import distance data
high_density_toluene_file = Path(r".\MDSimulationResults\toluene_distances_high_density.csv")
high_density_ligand_file = Path(r".\MDSimulationResults\dodecanol_distances_high_density.csv")
low_density_toluene_file = Path(r".\MDSimulationResults\toluene_distances_low_density.csv")
low_density_ligand_file = Path(r".\MDSimulationResults\dodecanol_distances_low_density.csv")

high_density_toluene_distanaces = np.genfromtxt(high_density_toluene_file, delimiter = ",")
high_density_ligand_distanaces = np.genfromtxt(high_density_ligand_file, delimiter = ",")
low_density_toluene_distanaces = np.genfromtxt(low_density_toluene_file, delimiter = ",")
low_density_ligand_distanaces = np.genfromtxt(low_density_ligand_file, delimiter = ",")


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


def lorentzian(x, A = 1, mu = 1, gamma = 1):
    return (A / np.pi) * (gamma / ((x - mu)**2 + gamma**2))

def calcluate_psi_squareds(distances, decay_constant = 1):
    return np.exp(-1*decay_constant * distances)**2



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



#% load the experimental data
H25_H8_HD["experimental spectrum"] = np.genfromtxt(r".\MIMS\PdSC12_Tol MIMS tau140 11730G 10K 9.466690 Q250207_03.csv", delimiter="," ,unpack = True, skip_header = 1)
D25_H8_HD["experimental spectrum"] = np.genfromtxt(r".\MIMS\dPdSC12_Tol MIMS tau140 11729G 10K 9.517308 Q250325_04_export.csv", delimiter=",", unpack = True, skip_header = 1)
H25_D8_HD["experimental spectrum"] = np.genfromtxt(r".\MIMS\Mims PdSC12 dTol tau 140 10K 43K scans.csv", delimiter=",", unpack = True, skip_header = 1)

qp.quickScatter(x = H25_H8_HD["experimental spectrum"][0], y = H25_H8_HD["experimental spectrum"][1])






#%% now fit individually

from lmfit import Model

def simulate_pulse_signal(ligand_distances, 
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
                          ):
    
    ligand_distances = np.atleast_1d(ligand_distances)
    solvent_distances = np.atleast_1d(solvent_distances)
     
    ligand_couplings  = calcluate_psi_squareds(ligand_distances, decay_constant = decay_constant)
    solvent_couplings = calcluate_psi_squareds(solvent_distances, decay_constant = decay_constant)
   
    if x is not None: # so, if we provided x-values, then use those.
        x_values = x.copy()
    else: # ifnot, then use something "reasonable"
        max_x = max([np.max(ligand_couplings), np.max(solvent_couplings)]) + 6*gamma
        x_values = np.linspace(mu-max_x, mu+max_x, 1000)
    
    temp_y_values = np.zeros_like(x_values)
    for coupling_pair in [[ligand_couplings, ligand_H_fraction], [solvent_couplings, solvent_H_fraction]]:
        couplings, H = coupling_pair
        for p in couplings:
            temp_y_values = temp_y_values + I*lorentzian(x_values, mu = mu + p*coupling_scaler, gamma = gamma) * H*p
            temp_y_values = temp_y_values + I*lorentzian(x_values, mu = mu - 1*p*coupling_scaler, gamma = gamma) * H*p
             
    return x_values, temp_y_values + a


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

def model_single_sample(x, 
                        ligand_distances, 
                        solvent_distances,
                        ligand_H_fraction,
                        solvent_H_fraction, 
                        decay_constant, 
                        I, 
                        mu, 
                        coupling_scaler, 
                        gamma, 
                        a,
                        ):
    
    x_sim, y_sim = simulate_pulse_signal(
        ligand_distances, 
        solvent_distances,
        ligand_H_fraction,
        solvent_H_fraction, 
        x = x,
        decay_constant = decay_constant, 
        I = I,
        mu = mu, 
        coupling_scaler = coupling_scaler,
        gamma = gamma, 
        a = a,       
        )
    #print(y_sim)
    return y_sim


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


def model_three_samples(x, 
                        ligand_distances, 
                        solvent_distances,
                        decay_constant, 
                        mu, 
                        coupling_scaler, 
                        gamma, 
                        I1,
                        I2,
                        I3,
                        a1,
                        a2,
                        a3,
                        ligand_H_fractions,
                        solvent_H_fractions,
                        spectral_indices,
                        ):
    
    total_ys = np.array([])
    
    for I, a, lH, sH, i in zip(
            [I1, I2, I3],
            [a1, a2, a3],
            ligand_H_fractions,
            solvent_H_fractions,
            range(len(ligand_H_fractions))
            ):
        
        
        part_x = x[spectral_indices[i]: spectral_indices[i+1]] # slice out the x-values associated with this part
        
            
        
        _, part_y_sim = simulate_pulse_signal(
            ligand_distances, 
            solvent_distances,
            lH,
            sH, 
            x = part_x,
            decay_constant = decay_constant, 
            I = I,
            mu = mu, 
            coupling_scaler = coupling_scaler,
            gamma = gamma, 
            a = a,       
            )
        
        total_ys = np.concatenate((total_ys, part_y_sim))
    return total_ys




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
