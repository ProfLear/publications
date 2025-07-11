# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 20:40:20 2025

@author: benle
"""

import numpy as np
import profilometryClasses as pc
from pathlib import Path
import codechembook.quickPlots as qp

#%% Paper stuff
paper_roughnesses = []
paper_all_heights = []
for i in range(1, 6):
    paper_file = Path(f"C:/Users/benle/Documents/GitHub/publications/Sandra Textiles/profilometry/10x_05xzoom_paper_0{i}.xyz")
    
    paperAreal = pc.makeSample(file = paper_file, instrument = "zygos")
    #paperAreal.raw.averaged.result.render()
    
    paperAreal.raw.averaged.result.fitRectbiSpline(name = "form")
    
    #paperAreal.raw.averaged.result.form.residual.render(zscale = 0.01)
    paperAreal.raw.averaged.result.form.residual.getArealRoughness()
    paper_roughnesses.append(paperAreal.raw.averaged.result.form.residual.rms_roughness)
    #% add all zs...
    for z in np.array(paperAreal.raw.averaged.result.form.residual.zs).flatten():
        paper_all_heights.append(z)
#%%
paper_roughnesses = np.array(paper_roughnesses)
print(f"For the uncoated paper, mean roughness is {np.mean(paper_roughnesses):.2f} with a standard error of {np.std(paper_roughnesses)/5:.2f}.")

paper_hist = qp.quickHist(x = paper_all_heights, width = 1, output = None)

#%% photothermal stuff <-- problem file 5 is the best to render
photothermal_roughnesses = []
photothermal_all_heights = []
for i in range(1, 6):
    photothermal_file = Path(f"C:/Users/benle/Documents/GitHub/publications/Sandra Textiles/profilometry/10x_05xzoom_laser_0{i}.xyz")
    
    photothermalAreal = pc.makeSample(file = photothermal_file, instrument = "zygos")
    photothermalAreal.raw.averaged.result.fitRectbiSpline(name = "form")
    
    #photothermalAreal.raw.averaged.result.form.residual.render(zscale = 0.1)
    photothermalAreal.raw.averaged.result.form.residual.getArealRoughness()
    
    photothermal_roughnesses.append(photothermalAreal.raw.averaged.result.form.residual.rms_roughness)
    #% add all zs...
    for z in np.array(photothermalAreal.raw.averaged.result.form.residual.zs).flatten():
        photothermal_all_heights.append(z)
photothermal_roughnesses = np.array(photothermal_roughnesses)
print(f"For the photothermally coated paper, mean roughness is {np.mean(photothermal_roughnesses):.2f} with a standard error of {np.std(photothermal_roughnesses)/5:.2f}.")

photothermal_hist = qp.quickHist(x = photothermal_all_heights, width = 1, output = None)

#%%
from plotly.subplots import make_subplots

all_hist = make_subplots(rows = 2, cols = 1)
all_hist.add_bar(x = paper_hist.data[0].x, y = paper_hist.data[0].y, 
                 marker = dict(color = "#0773b1"),
                 row = 1, col = 1, name = "uncoated paper",
                 showlegend = False,
                 )
all_hist.add_bar(x = photothermal_hist.data[0].x, y = photothermal_hist.data[0].y, 
                 marker = dict(color = "#a33769"),
                 row = 2, col = 1, name = "photothermally cured",
                 showlegend = False,
                 )

all_hist.add_annotation(text = "uncoated paper",
                        font = dict(color = "#0773b1"),
                        row = 1, col = 1,
                        x = 0, y = 600000,
                        xanchor = "left",
                        showarrow = False,
                        )
all_hist.add_annotation(text = "photothermally-coated",
                        font = dict(color = "#a33769"),
                        row = 2, col = 1,
                        x = 0, y = 450000,
                        xanchor = "left",
                        showarrow = False,
                        )

all_hist.update_xaxes(range = [-30, 20], title = "heights /microns", zeroline = True)
all_hist.update_yaxes(range = [0, max([max(photothermal_hist.data[0].y), max(paper_hist.data[0].y)])*1.1], title= "counts", showline=False,
                      nticks = 4)
all_hist.update_layout(template = "simple_white", bargap = 0,
                       margin = dict(l = 0, r = 0, t = 50, b = 0),
                       height = 300*2, width = 3.3*300,
                       font = dict(size = 20),
                       )
all_hist.show("png")

#%%
np.sum((paper_all_heights - np.mean(paper_all_heights))**2/len(paper_all_heights))**0.5
np.sum((photothermal_all_heights - np.mean(photothermal_all_heights))**2/len(photothermal_all_heights))**0.5

#%% render things with a shared scale...
import profilometryClasses as pc

paperAreal.raw.averaged.result.form.residual.render(zscale = 500, zmode = "absolute")
photothermalAreal.raw.averaged.result.form.residual.render(zscale = 500, zmode = "absolute")


#%%

def rms_roughness(array):
    return np.sum(np.mean((array-np.mean(array))**2))**0.5