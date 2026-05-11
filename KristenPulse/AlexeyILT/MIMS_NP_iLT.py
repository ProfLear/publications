# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 14:56:06 2026

@author: Alexey
"""

import brukerread as br
import numpy as np
from scipy.optimize import nnls
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

##############################################################################
# Settings
##############################################################################

fpath = ".//"
fname = []; title = []; tau_ns  = []; field = []; lwcolor = []; scale = [];
# fname.append("Q250207_03.DSC");title.append("NA/NA"  );
# tau_ns.append(140); field.append(11730); lwcolor.append('darkgray')
fname.append("Q241204_07.DSC");title.append("nat.abund"  ); scale.append(1);
tau_ns.append(140); field.append(11745); lwcolor.append('chocolate')
fname.append("Q250207_07.DSC");title.append("d-NP"  ); scale.append(3);
tau_ns.append(140); field.append(11729); lwcolor.append('olivedrab')
fname.append("Q241204_08.DSC");title.append("d-Tol"); scale.append(1);
tau_ns.append(140); field.append(11745); lwcolor.append('slateblue')

r_start = 2
r_end = 20 #was 20 to start


input_dict = {
    'fpath':fpath, 
    'fname':None,
    'field':None,
    'tau_ns':None,
    'r_start':r_start,  # [Angst]
    'r_end': r_end,    # [Angst]
    'nPoints': (r_end - r_start)*3,  # # must be less or equal to experimental number of points, or will be reduced
    'r_scale':'log', #  'log' (log10 base) or 'linear'
    ### things to experimentwith vvvvv
    'BL_damp':0.1, # 0 - full blindspot, 1 - no blindspot
    'Fermi_Contact_a': 1, # [Angst] exponential factor on how Aiso depends on distance, probably <1 if to judge from s-type orbital 
    'Fermi_Contact_A0':1500, # [MHz] preexponential factor for Aiso at r=0, based on Stoll's Hdot study ... probably does not matter as isotropic goes down quick 
    'Dipolar_A0_fudge':1, # Fudge factor for correcting dipolar contributios 1==no fudging .. realistically, it is just a scale for the effective distance
    'unit_lw':0.055, # [MHz] Unitary liinewidth for Pake pattern = sigma in gaussian
    'Regularisation':50,
    }


###############################################################################
# Generate ENDOR spectrum based on 
#    - Pake pattern with added isotropic shifts
#    - Mims ENDOR blindspots with a blinespot dampenning
# nu - RF Frequency - Larmor (MHz)
# T  - Magnitude of dipol-dipol HFC (MHz)
# Aiso - isotropic HFC
# tau_us - tau for calculating blind spot (us) Note that typical measurement uses [ns] !
# BL_damp - How strict the blindspot. 0 - 100% suppression, 1 - no suppression
# lw - unitary linewidth [MHz]
def Pake(nu, T, Aiso, tau_us, BL_damp=0, lw=0):
    # Powder average: sample over theta
    theta = np.linspace(0, np.pi/2, 10000)
    weights = np.sin(theta)  # powder averaging weight
    
    Spec = np.zeros_like(nu)
    
    for ms in [+0.5, -0.5]:
        # A(theta) = T * (1 - 3*cos^2(theta))
        A = T * (1 - 3 * np.cos(theta)**2)+Aiso
        nu_res = ms * A  # resonance frequency for each theta at zero Larmor

        # Bin contributions into spectrum
        # histogram does the binning natively
        I_contrib, _ = np.histogram(nu_res, bins=len(nu),
                                range=(nu[0], nu[-1]), 
                                weights=weights)
        Spec += I_contrib
    Spec*=(1.0-BL_damp)*np.sin(2*np.pi*nu*tau_us)**2 + BL_damp
    # make line broadenning if lw is not 0
    if lw>0:
        nu_kernel = np.linspace(-nu[-1]/2, nu[-1]/2, len(nu))
        kernel = np.exp(-nu_kernel**2 / (2 * lw**2))
        kernel /= kernel.sum()
        I_fft      = np.fft.rfft(Spec)
        kernel_fft = np.fft.rfft(kernel)
        Spec    = np.fft.fftshift(np.fft.irfft(I_fft * kernel_fft, n=len(Spec)))
    return Spec

###############################################################################
# Main Function that does it all. 
def Process_file(input_dict):
    nMean = 40  # number of points at the beginning and the end of the data set to use for baseline
    filename = input_dict['fpath']+input_dict['fname'] #
    Field = input_dict['field'] #11730 # [G]
    tau_ns = input_dict['tau_ns'] #140  # [ns] tau in Mims
    
    
    BL_damp = input_dict['BL_damp'] #0.05
    r_start = input_dict['r_start'] #2.5 # [Angst]
    r_end   = input_dict['r_end'] #15 # [Angst]
    nPoints = input_dict['nPoints'] #500 # must be less or equal to experimental number of points

    r_scale = input_dict['r_scale'] #'log' # 'log' (log10 base) or 'linear'

    Fermi_Contact_a = input_dict['Fermi_Contact_a'] #0.5  # [Angst] exponential factor on how Aiso depends on distance, probably <1 if to judge from s-type orbital 
    Fermi_Contact_A0   = input_dict['Fermi_Contact_A0'] #1500  # [MHz] preexponential factor for Aiso at r=0 ... probably negligiblle based on Stoll's Hdot study
    # A0​=(2μ0/3h)⋅ge⋅μB⋅gn⋅μn​⋅ρ​(r​) 

    Dipolar_A0_fudge         = input_dict['Dipolar_A0_fudge'] # 1 # Fudge factor for correcting dipolar contributios 1==no fudging

    unit_lw = input_dict['unit_lw'] #0.06 # [MHz] Unitary liinewidth for Pake pattern = sigma in gaussian

    Regularisation =input_dict['Regularisation'] #500
    
    
    nu_G = 0.004262 # MHz/G
    T_1A = 78.973 # MHz at 1A for 1H
    
    ##############################################################################
    # Load File
    ##############################################################################
    
    ax, y, dsc = br.brukerread(filename)
    nu_n = nu_G*Field
    
    xFreq = np.array(ax['x']-nu_n) # make a copy, just in case ... python crap
    
    yReal = np.array(y.real, dtype=np.float64) # make a copy, just in case we need original ... python crap

    my = (np.mean(yReal[0:nMean])+np.mean(yReal[-nMean:]))/2.0  # average of first and last 10 points is the BL
    yReal/=-my
    yReal+=1.0
    
    ##############################################################################
    # Biuld kernel
    ##############################################################################
    tau_us = tau_ns*1e-3
    if nPoints>len(xFreq):
        nPoints=len(xFreq)
        print(f'nPoints was adjusted to nPoints={len(xFreq)}')
    if r_scale == 'log':
        xR = np.logspace(np.log10(r_start), np.log10(r_end), nPoints)
    else:
        xR = np.linspace(r_start, r_end, nPoints)
    
    Kernel = np.zeros((len(xFreq), len(xR)))
    
    T = Dipolar_A0_fudge*T_1A/xR**3
    Aiso = Fermi_Contact_A0*np.exp(-2*xR/Fermi_Contact_a)
    
    for ii, rR in enumerate(xR):
        Kernel[:, ii] = Pake(xFreq, T[ii], Aiso[ii], tau_us, BL_damp, unit_lw)
    
    ##############################################################################
    # Simple SVD iLT
    ##############################################################################
    U, S, Vt = np.linalg.svd(Kernel)
    
    k =len(xR)
    
    coeffs1 = (U[:, :k].T @ yReal) / S[:k]
    svd_Y = U[:, :k] @ (S[:k] * coeffs1)
    
    filters2 = S[:k]  / ((S[:k] )**2 + Regularisation**2)
    coeffs2 = filters2 * (U[:, :k].T @ yReal)
    reg_Y = U[:, :k] @ (S[:k] * coeffs2)
    
    dist_svd = Vt[:k, :].T @ coeffs1
    dist_reg = Vt[:k, :].T @ coeffs2
    
    # filters = sigma / (sigma**2 + lam**2)
    #     coeffs = filters * (U.T @ b)
    #     return Vt.T @ coeffs
    
    ##############################################################################
    # Non-negative Tikhonov (NNLS) iLT
    ##############################################################################
    
    # Build augmented system: [K; lambda*I] p = [y; 0]
    # This encodes the Tikhonov penalty while enforcing p >= 0 via nnls
    
    K_aug = np.vstack([Kernel, Regularisation * np.eye(len(xR))])
    y_aug = np.concatenate([yReal, np.zeros(len(xR))])
    
    dist_nn, residual = nnls(K_aug, y_aug)
    
    # Reconstruct the fit in frequency domain
    nn_Y = Kernel @ dist_nn


    return xFreq, yReal, svd_Y, reg_Y, nn_Y, xR, dist_nn, dist_reg, T, Aiso


##############################################################################
# Prepare to plot results
##############################################################################
fig = plt.figure(1, figsize=(12, 4))
fig.clf()
axs = fig.subplot_mosaic([['ENDOR', 'P(r)'], ['ENDOR', 'As']])

shY= 0.2

for ii, _ in enumerate(fname):
    input_dict['fname']  = fname[ii] #
    input_dict['field']  = field[ii]
    input_dict['tau_ns'] = tau_ns[ii]

    svdlab = ''
    tikhlab = ''
    tikhlabdiff = ''
    blspotlab = ''
    if ii==len(fname)-1:
        svdlab = 'unconstr. Tikhonov (SVD)'
        tikhlab = '$Reg_{NNLS}$ (Tikhonov NNLS fit)'
        tikhlabdiff = 'exp-$Reg_{NNLS}$'
        blspotlab = 'Mims Blindspot'
        
    xFreq, yReal, svd_Y, reg_Y, nn_Y, xR, dist_nn, dist_reg, T, Aiso = Process_file(input_dict)
    
    # --- Left panel: frequency domain ---
    if scale[ii] != 1:
        axs['ENDOR'].plot(xFreq, scale[ii]*yReal+shY*ii,  label=f'{fname[ii]} (x{scale[ii]})', color=lwcolor[ii], linewidth=3)
    else:
        axs['ENDOR'].plot(xFreq, scale[ii]*yReal+shY*ii,  label=f'{fname[ii]}', color=lwcolor[ii], linewidth=3)
    
    # axs['ENDOR'].plot(xFreq, scale[ii]*svd_Y+shY*ii,  label=svdlab, color='gray', linewidth=0.5)
    axs['ENDOR'].plot(xFreq, scale[ii]*reg_Y+shY*ii,  label=svdlab, color='black', linewidth=0.5)
    axs['ENDOR'].plot(xFreq, scale[ii]*nn_Y+shY*ii,  label= tikhlab, color='maroon', linewidth=2, linestyle=':')
    axs['ENDOR'].plot(xFreq, scale[ii]*(yReal-nn_Y)+shY*ii,  label=tikhlabdiff, color='turquoise', linewidth=1, linestyle=':')

    # axs['ENDOR'].plot(xFreq, Kernel[:,-1],  label='Tikhonov non_neg fit', color='k', linewidth=3, linestyle='-')
    BL = (1.0-input_dict['BL_damp'])*np.sin(2*np.pi*xFreq*input_dict['tau_ns']*1e-3)**2 + input_dict['BL_damp']
    axs['ENDOR'].plot(xFreq, BL*shY +shY*ii, label=blspotlab, color='gray', linewidth=0.5, linestyle='--')
    
    # --- Right panel: distance domain ---
    # axs['P(r)'].plot(xR, dist_svd, label='SVD distribution')
    axs['P(r)'].plot(xR, scale[ii]*dist_nn, label=f'{title[ii]} NNLS ', color=lwcolor[ii], linewidth=2, linestyle='-', marker='o')

    axs['P(r)'].plot(xR, scale[ii]*dist_reg, label=svdlab, color='gray', linewidth=0.5)
    axs['P(r)'].grid(axis='both', color='0.95', which='both')
    

    
axs['ENDOR'].set_xlabel(r'$\nu_{RF}$-$\nu_L$, MHz')
axs['ENDOR'].set_ylabel('ENDOR intensity, arb.un.')
axs['ENDOR'].set_title('Mims ENDOR')
axs['ENDOR'].legend()


axs['P(r)'].set_xlabel('Effective distance, Å')
axs['P(r)'].set_ylabel('P(r), arb.un.')
axs['P(r)'].set_title(r'Distance Distribution of ${^1}$H')
axs['P(r)'].set_xscale('log')
# axs['P(r)'].xaxis.set_minor_locator(ticker.LogLocator(base=10, numticks=50))
axs['P(r)'].legend()
    

axs['As'].plot(xR, T, label='T', color='black', linewidth=1, marker='.')
axs['As'].plot(xR, Aiso, label='Aiso', color='blue', linewidth=1, marker='.')
axs['As'].set_xscale('log')
axs['As'].set_yscale('log')
axs['As'].set_ylim(1e-4, 30)
axs['As'].grid(axis='both', color='0.95', which='both')
axs['As'].legend()
axs['As'].set_xlabel('Effective distance, Å')
axs['As'].set_ylabel('HF coupling constant, MHz')
axs['As'].set_title(r'Distance Distribution of HF couplings (modelled)')
axs['As'].set_title('')
    
plt.tight_layout()
plt.show()
# ###############################################################################
# # Test Pake (comment when not needed)
# sig1 = Pake(xFreq, 5, 0, 0.0)
# sig2 = Pake(xFreq, 5, 0, 0.1)

# plt.figure(figsize=(8, 4))
# plt.plot(xFreq, yReal)
# plt.plot(xFreq, sig1)
# plt.plot(xFreq, sig2)
# plt.tight_layout()
# plt.show()
# ###############################################################################