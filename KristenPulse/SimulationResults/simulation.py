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
high_density_toluene_file = Path(r"C:\Users\benle\Documents\GitHub\publications\KristenPulse\SimulationResults\toluene_distances_high_density.csv")
high_density_ligand_file = Path(r"C:\Users\benle\Documents\GitHub\publications\KristenPulse\SimulationResults\dodecanol_distances_high_density.csv")
low_density_toluene_file = Path(r"C:\Users\benle\Documents\GitHub\publications\KristenPulse\SimulationResults\toluene_distances_low_density.csv")
low_density_ligand_file = Path(r"C:\Users\benle\Documents\GitHub\publications\KristenPulse\SimulationResults\dodecanol_distances_low_density.csv")

high_density_toluene_distanaces = np.genfromtxt(high_density_toluene_file, delimiter = ",")
high_density_ligand_distanaces = np.genfromtxt(high_density_ligand_file, delimiter = ",")
low_density_toluene_distanaces = np.genfromtxt(low_density_toluene_file, delimiter = ",")
low_density_ligand_distanaces = np.genfromtxt(low_density_ligand_file, delimiter = ",")


high_density_ligand_distance_hist = qp.quickHist(x = high_density_ligand_distanaces)
high_density_toluene_distance_hist= qp.quickHist(x = high_density_toluene_distanaces)
low_density_ligand_distance_hist= qp.quickHist(x = low_density_ligand_distanaces)
low_density_toluene_distance_hist= qp.quickHist(x = low_density_toluene_distanaces)


#%% import mims data


#%% now simulate


def lorentzian(x, A = 1, mu = 1, gamma = 1):
    return (A / np.pi) * (gamma / ((x - mu)**2 + gamma**2))

def calcluate_psi_squareds(distances, decay_constant = 1):
    return np.exp(-1*decay_constant * distances)**2

def simulate_pulse_signal(experiment, decay_constant = 1, gamma = 1, coupling_scaler = 1):
    experiment["ligand"]["couplings"] = calcluate_psi_squareds(experiment["ligand"]["distances"], decay_constant = decay_constant)
    experiment["solvent"]["couplings"] = calcluate_psi_squareds(experiment["solvent"]["distances"], decay_constant = decay_constant)
   
    max_x = max([np.max(experiment["ligand"]["couplings"]), np.max(experiment["solvent"]["couplings"])]) + 6*gamma
    
    x_values = np.linspace(-max_x, max_x, 1000)
    for source in experiment:
        temp_y_values = np.zeros_like(x_values)
        for p in experiment[source]["couplings"]:
            temp_y_values = temp_y_values + lorentzian(x_values, mu = p*coupling_scaler, gamma = gamma)
            temp_y_values = temp_y_values + lorentzian(x_values, mu = -1*p*coupling_scaler, gamma = gamma)
        
        experiment[source]["spectrum contribution"] = temp_y_values.copy()*experiment[source]["fraction hydrogen"]
    experiment["x values"] = x_values   
    experiment['total spectrum'] = experiment["ligand"]["spectrum contribution"] + experiment["solvent"]["spectrum contribution"]
    
    return x_values, experiment['total spectrum'] # return as tuple of lists

dc = 0.18
g = 0.008
l = 0.2

H25_H8_HD = {
    "ligand"  : {"distances": high_density_ligand_distanaces,  "fraction hydrogen" : 1},
    "solvent" : {"distances": high_density_toluene_distanaces, "fraction hydrogen" : 1},
    }
D25_H8_HD = {
    "ligand"  : {"distances": high_density_ligand_distanaces,  "fraction hydrogen" : 0.012},
    "solvent" : {"distances": high_density_toluene_distanaces, "fraction hydrogen" : 1},
    }
H25_D8_HD = {
    "ligand"  : {"distances": high_density_ligand_distanaces,  "fraction hydrogen" : 1},
    "solvent" : {"distances": high_density_toluene_distanaces, "fraction hydrogen" : 0.01},
    }

simulate_pulse_signal(H25_H8_HD, decay_constant = dc, gamma = g, coupling_scaler = l)
simulate_pulse_signal(D25_H8_HD, decay_constant = dc, gamma = g, coupling_scaler = l)
simulate_pulse_signal(H25_D8_HD, decay_constant = dc, gamma = g, coupling_scaler = l)

total_sim = qp.quickScatter(
    x = [H25_H8_HD["x values"],       H25_D8_HD["x values"],       D25_H8_HD["x values"]], 
    y = [H25_H8_HD['total spectrum'], H25_D8_HD['total spectrum'], D25_H8_HD['total spectrum']], 
    name = ["all", "ligand", "solvent"],
    output = None)  

total_sim.update_xaxes(title = "coupling value", range = [-0.1, 0.1])
total_sim.update_yaxes(title = "relative intensity")
total_sim.update_layout(margin = dict(l = 0, r=0, t = 0, b = 0))
total_sim.show('png')



#%% now fit