# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 21:43:45 2025

@author: benle
"""

import numpy as np
from scipy.optimize import lsq_linear
from plotly.subplots import make_subplots
from pathlib import Path


# Simulate some data: sum of exponentials with noise
def simulate_data(t, taus, amps, noise_level=0.01):
    y = np.zeros_like(t)
    for tau, amp in zip(taus, amps):
        y += amp * np.exp(-t / tau)
    noise = noise_level * np.random.normal(size=len(t))
    return y + noise

# Inverse Laplace via discretized tau space with optional regularization
def inverse_laplace_fit(t, f, tau_grid, regularization=0.0):
    K = np.exp(-np.outer(t, 1 / tau_grid))

    if regularization > 0:
        ATA = K.T @ K + regularization * np.eye(K.shape[1])
        ATy = K.T @ f
        # Use regularized design matrix in lsq_linear as workaround
        K_reg = np.vstack([K, np.sqrt(regularization) * np.eye(K.shape[1])])
        f_reg = np.concatenate([f, np.zeros(K.shape[1])])
        result = lsq_linear(K_reg, f_reg, bounds=(0, np.inf))
    else:
        result = lsq_linear(K, f, bounds=(0, np.inf))
    
    return result.x

# AIC calculation
def compute_aic(f_data, f_fit, weights, threshold=1e-6):
    n = len(f_data)
    rss = np.sum((f_data - f_fit) ** 2)
    k = np.count_nonzero(weights > threshold)
    aic = 2 * k + n * np.log(rss / n)
    return aic

folder = Path("C:/Users/benle/Downloads")
file = Path("PdSC12_Tol 2pESEEM tau180 40K 9.458519 Q250212_20.csv")
# Load real data
t, f = np.genfromtxt(folder/file,
    unpack=True, delimiter=",")

# AIC model selection loop
grid_sizes = range(1, 100, 2)
regularization = 1e-2
aic_values = []
fit_results = []
threshold = 1e-6

for n_grid in grid_sizes:
    tau_grid = np.logspace(0, 6, n_grid)
    weights = inverse_laplace_fit(t, f, tau_grid, regularization=regularization)
    K = np.exp(-np.outer(t, 1 / tau_grid))
    f_fit = K @ weights
    aic = compute_aic(f, f_fit, weights, threshold = threshold)
    aic_values.append(aic)
    fit_results.append((tau_grid, weights, f_fit))

# Find best model
best_idx = np.argmin(aic_values)
best_tau, best_weights, best_f_fit = fit_results[best_idx]

#%%
fig = make_subplots(rows = 3, cols = 1)
fig.add_scatter(x = list(grid_sizes), y = aic_values, 
                mode = "lines+markers", line = dict(color = "grey"), marker = dict(color = "grey"),
                showlegend = False, row = 1, col = 1)
fig.add_scatter(x = [list(grid_sizes)[best_idx]], y = [aic_values[best_idx]], mode = "markers", showlegend = False, marker = dict(color = "red"), row = 1, col = 1)
fig.add_annotation(text = f"n = {list(grid_sizes)[best_idx]}", x = list(grid_sizes)[best_idx], y = aic_values[best_idx], arrowcolor = "red", font = dict(color = "red"), row = 1, col = 1)

fig.add_scatter(x = t, y = f, mode = "lines", showlegend = False, line = dict(width = 3, color = "grey"), row = 2, col = 1)
fig.add_scatter(x = t, y = best_f_fit, mode = "lines", line = dict(color = "red",), showlegend = False, row = 2, col = 1)



fig.add_scatter(y = best_weights, x = best_tau, 
                mode = "lines", line = dict(color = "grey", shape = "hvh"),
                showlegend = False, row =3, col = 1)
fig.add_scatter(y = best_weights[best_weights > threshold], x = best_tau[best_weights > threshold], 
                mode = "markers", marker = dict(color = "red"),
                showlegend = False, row =3, col = 1)
laplace_string = "I(t) ~ "
for weight, tau in zip(best_weights[best_weights > threshold], best_tau[best_weights > threshold]):
    laplace_string = laplace_string + f"{weight:.0f}exp(-{tau:.0f}t) + "
    '''
    fig.add_annotation(text = f"</br>A ={weight:.0f}</br>k = {tau:.0f}",
                    x = np.log10(tau), y = weight, 
                    #xref = "x3", yref = "y3",
                    arrowcolor = "red", font = dict(color = "red"),
                    row = 3, col = 1)
    '''
fig.add_annotation(text = laplace_string[:-2],
                   font = dict(color = "red"),
                   xanchor = "right",
                   x = max(t), y = max(f),
                   showarrow = False,
                   row = 2, col = 1)

fig.update_xaxes(title = "number of grid points", row = 1, col = 1)
fig.update_yaxes(title = "AIC", row = 1, col = 1)

fig.update_xaxes(title = "time /ns", row = 2, col = 1)
fig.update_yaxes(title = "intensity", row = 2, col = 1)

fig.update_xaxes(title = "tau", type = "log", row = 3, col = 1)
fig.update_yaxes(title = "weight", range = [0, max(best_weights[best_weights>threshold])*1.05], row = 3, col = 1)
fig.update_layout(title = file.stem,
    template = "simple_white",
    height = 800, width = 800)
fig.show("png")


#%% might be intereseting to then use this as initial guesses for a more conventional fit?
from lmfit import Model

def zero_background(x):
    return 0 

def exp_decay(x, A, tau):
    return np.array(A*np.exp(-x/tau))


'''
thinking of ways to try to process the results.
1. Can just pass them on
2. Can look for all indices with weights above the threshold, and then average out all contiguous units
3. Can order by weights and systematically drop the lowest weights one after another
4. COmbine 2 and 3... take all contiguous points and calculate their collective weights. Then average out the ones wiht the lowest weights, one after one another. 
'''
#%% 4
def extract_contiguous_regions(mask, weights, taus):
    c_weights = []
    c_taus= []
    start = None

    for i, val in enumerate(mask):
        if val:
            if start is None:
                start = i  # Start of new region
        elif start is not None:
            # End of region
            c_weights.append(list(weights[start:i]))
            c_taus.append(list(taus[start:i]))
            #regions.append((weights[start:i], taus[start:i]))
            start = None

    # Handle case where mask ends in True
    if start is not None:
        c_weights.append(weights[start:])
        c_taus.append(taus[start:])
        #regions.append((values[start:], taus[start:]))

    return [c_weights, c_taus]

contiguous_weights, contiguous_taus = extract_contiguous_regions(best_weights > threshold, best_weights, best_tau)

def sum_inside_contiguous_region(regions):
    sums = []
    for region in regions:
        total = 0
        for value in region:
            total =+ value
        sums.append(total)
    return sums

contiguous_sums = sum_inside_contiguous_region(contiguous_weights)
print(contiguous_sums)

sorted_result = sorted(zip(contiguous_sums, contiguous_weights, contiguous_taus), reverse=True)

#now, go through this in order (smallest to largest) and, each time there are two values (or more) average them.

#start by fitting the un-avarged stuff

# then go thorugh and fit each subsequent modification
for sw, st in zip(sorted_weights, sorted_taus):
    if len(sw) > 1:
        # calculate weighted avarge of tau
        # sum up the weights
        
        # fit with modified parameters
        # record AIC

#%%


def identify_contiguous_regions():
    
    return regions

#3 Order all invididual ones by weight, and then start trimming off the end...
sorted_weights = sorted(best_weights)
sorted_taus = # sort same as done for "best_weights

fits = []
for i in range(len(sorted_weights)):
    fits.append(fit_decay(t, f, sorted_weights[:-i], sorted_taus[:-i]))

current_AIC = fits[0].aic
best_fit = fits[0]
for fit in fits:
    if fit.aic < curent_AIC:
        best_fit = fit
        
      
    


#4
    

decay_model  = Model(zero_background)
decay_model_params = decay_model.make_params()
for i, weight, tau in zip(range(1, len(best_weights[best_weights > threshold])+2)[:-1], best_weights[best_weights > threshold][:-1], best_tau[best_weights > threshold][:-1]):
    exp_component = Model(exp_decay, prefix = f"c{i}_")
    exp_component_params = exp_component.make_params()
    decay_model_params.update(exp_component.make_params(
        A =   dict(value=weight, vary = True, min=0,),
        tau = dict(value=tau,    vary = True, min=0,),
        ),
        )
    decay_model = decay_model + exp_component

init = decay_model.eval(decay_model_params, x=t)
decay_fit = decay_model.fit(f, x = t, params = decay_model_params)
print(decay_fit.fit_report())
from codechembook.quickPlots import plotFit

fig = plotFit(decay_fit)
fig.add_scatter(x=t, y = init)
fig.show("png")
