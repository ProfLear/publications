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

high_density_toluene_distanaces = np.genfromtxt(high_density_toluene_file, delimiter = ",")
high_density_ligand_distanaces = np.genfromtxt(high_density_ligand_file, delimiter = ",")

high_density_ligand_distance_hist = qp.quickHist(x = high_density_ligand_distanaces)
high_density_toluene_distance_hist= qp.quickHist(x = high_density_toluene_distanaces)


# now simulate

def lorentzian(x, A = 1, mu = 1, gamma = 1):
    return (A / np.pi) * (gamma / ((x - mu)**2 + gamma**2))


def simulate_pulse_signal(distances, decay_constant = 1, gamma = 1, coupling_scaler = 1):
    psi_squareds = []
    for distance in distances:  
        psi_squareds.append((np.exp(-1*decay_constant * distance))**2)
        
    max_x = np.max(psi_squareds) + 6*gamma
    
    x_values = np.linspace(-max_x, max_x, 1000)
    y_values = np.zeros_like(x_values)
    for p in psi_squareds:
        y_values = y_values + lorentzian(x_values, mu = p*coupling_scaler, gamma = gamma)
        y_values = y_values + lorentzian(x_values, mu = -1*p*coupling_scaler, gamma = gamma)
        
    
    return x_values, y_values # return as tuple of lists

dc = 0.18
g = 0.008
l = 0.2
x_solvent, y_solvent = simulate_pulse_signal(
    high_density_toluene_distanaces,
    decay_constant = dc,
    gamma = g,
    )

x_ligand, y_ligand = simulate_pulse_signal(
    high_density_ligand_distanaces,
    decay_constant = dc,
    gamma = g
    )
   
total_sim = qp.quickScatter(x = [x_ligand, x_solvent], y = [y_ligand, y_solvent], output = None)  
solvent_sim = qp.quickScatter(x = x_solvent, y = y_solvent, output = None) 

total_sim.update_xaxes(range = [-l, l])
total_sim.show("png")