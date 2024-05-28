# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 19:13:13 2024

@author: Lukas
"""

from d2_plotter import D2_Plotter
from scan_plotter import Timing_scan
from integrator import Integrator
import pandas as pd
from analysis_pipeline import data_analysis
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 20
plt.rcParams['lines.markersize'] = 6
plt.rcParams['lines.linewidth'] = 2

def baseline_als(y, lam, p, niter=10):
    """
    Function for baseline fitting
    
    Parameters
    ----------
    y: Array
        Intensity data  
    lam: float
        Smoothness of the curve
    p: float
        Parameter for y-position of baseline relative to data;
        0.5 => Line is in the middle                
    """
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

values_scans = []

scan_numbers = [4141, 4142, 4113]
experiments = ["nai3moltiming6", "nai3moltiming6", "nai3moltiming6"]

for i,j in enumerate(list(zip(scan_numbers, experiments))):
    scan = j[0]
    experiment = j[1]
    roi = (1455,126,60,21)
    
    # Script loading
    fast_scan = Timing_scan()
    fast_scan_detector = D2_Plotter()
    fast_scan.scan_number = scan
    fast_scan.experiment = experiment
    fast_scan.roi = roi
    
    fast_scan_integrator = Integrator()
    fast_scan_integrator.image_number = 0
    fast_scan_integrator.experiment = experiment
    fast_scan_integrator.scan_number = scan
    fast_scan_integrator.xlims = [roi[0], roi[0]+roi[2]]
    fast_scan_integrator.ylims = [roi[1], roi[1]+roi[3]]
    
    fast_scan_detector.scan_number = scan
    fast_scan_detector.experiment = experiment
    fast_scan_detector.rois = {"roi": roi}
    
    fast_scan()
    fast_scan.plot_intensity()
    fast_scan.plot_fitted_params_x()
    fast_scan_detector()
    fast_scan_integrator()
    
    values_scans.append(fast_scan.output_intensity())

fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(values_scans[0][0]), np.array(values_scans[0][1]), yerr = np.array(values_scans[0][2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(values_scans[0][0]), baseline_als(np.array(values_scans[0][1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(np.array(values_scans[0][0]), np.array(values_scans[0][1])-np.array(values_scans[0][2]), np.array(values_scans[0][1])+np.array(values_scans[0][2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
# ax1.axvline(x = 0.6, color = "red")
# ax1.axvline(x = 3.5, color = "red")
ax1.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(values_scans[1][0]), np.array(values_scans[1][1]), yerr = np.array(values_scans[1][2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(values_scans[1][0]), baseline_als(np.array(values_scans[1][1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(np.array(values_scans[1][0]), np.array(values_scans[1][1])-np.array(values_scans[1][2]), np.array(values_scans[1][1])+np.array(values_scans[1][2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
# ax1.axvline(x = -1.2, color = "red")
# ax1.axvline(x = 0.9, color = "red", linestyle = "--")
# ax1.axvline(x = 3.2, color = "red")
ax1.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

fig, (ax1, ax2) = plt.subplots(2, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(values_scans[0][0]), np.array(values_scans[0][1]), yerr = np.array(values_scans[0][2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(values_scans[0][0]), baseline_als(np.array(values_scans[0][1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(np.array(values_scans[0][0]), np.array(values_scans[0][1])-np.array(values_scans[0][2]), np.array(values_scans[0][1])+np.array(values_scans[0][2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
# ax1.axvline(x = 0.6, color = "red")
# ax1.axvline(x = 3.5, color = "red")
ax2.errorbar(np.array(values_scans[1][0]), np.array(values_scans[1][1]), yerr = np.array(values_scans[1][2]), color='#1c9e15', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax2.plot(np.array(values_scans[1][0]), baseline_als(np.array(values_scans[1][1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax2.fill_between(np.array(values_scans[1][0]), np.array(values_scans[1][1])-np.array(values_scans[1][2]), np.array(values_scans[1][1])+np.array(values_scans[1][2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
ax2.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

shift = 0
fig, (ax1, ax2) = plt.subplots(2, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(values_scans[0][0]), np.array(values_scans[0][1]), yerr = np.array(values_scans[0][2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(values_scans[0][0]), baseline_als(np.array(values_scans[0][1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(np.array(values_scans[0][0]), np.array(values_scans[0][1])-np.array(values_scans[0][2]), np.array(values_scans[0][1])+np.array(values_scans[0][2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
ax2.errorbar(np.array(values_scans[1][0])+shift, np.array(values_scans[1][1]), yerr = np.array(values_scans[1][2]), color='#1c9e15', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax2.plot(np.array(values_scans[1][0])+shift, baseline_als(np.array(values_scans[1][1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax2.fill_between(np.array(values_scans[1][0])+shift, np.array(values_scans[1][1])-np.array(values_scans[1][2]), np.array(values_scans[1][1])+np.array(values_scans[1][2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
ax2.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

fig =  plt.figure(dpi = 300)
for i in range(0, len(scan_numbers)):
    if i == 0:
        plt.subplot(len(scan_numbers), 1, i+1)
        ax = fig.gca()
    else:
        plt.subplot(len(scan_numbers), 1, i+1, sharex = ax)
        ax = fig.gca()
    ax.errorbar(np.array(values_scans[i][0]), np.array(values_scans[i][1]), yerr = np.array(values_scans[i][2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
    ax.plot(np.array(values_scans[i][0]), baseline_als(np.array(values_scans[i][1]),0.2, 0.5), linestyle = "--", color='#11570d')
    ax.fill_between(np.array(values_scans[i][0]), np.array(values_scans[i][1])-np.array(values_scans[i][2]), np.array(values_scans[i][1])+np.array(values_scans[i][2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
ax.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

def concatenate(params1: list, params2: list) -> list:
    params_list = []
    params1 = list([list(map(lambda x: round(x, 2), params1[0].tolist())), params1[1].tolist(), params1[2].tolist()])
    params2 = list([list(map(lambda x: round(x, 2), params2[0].tolist())), params2[1].tolist(), params2[2].tolist()])
    for x in params1[0]:
        if (x not in params2[0]):
            params_list.append([x, params1[1][params1[0].index(x)], params1[2][params1[0].index(x)]])
        else:
            params_list.append([x, (params1[1][params1[0].index(x)] + params2[1][params2[0].index(x)])/2,  (params1[2][params1[0].index(x)] + params2[2][params2[0].index(x)])/2])
    for x in params2[0]:
        if x not in params1[0]:
            params_list.append([x, params2[1][params2[0].index(x)], params2[2][params2[0].index(x)]])
    param_list_sorted = sorted(params_list, key = lambda x: x[0])
    print(param_list_sorted)
    return [[x[0] for x in param_list_sorted], [x[1] for x in param_list_sorted], [x[2] for x in param_list_sorted]]

sumed_values = concatenate(values_scans[0], [values_scans[1][0]+shift, values_scans[1][1], values_scans[1][2]])
fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(sumed_values[0]), np.array(sumed_values[1]), yerr = np.array(sumed_values[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(sumed_values[0]), baseline_als(np.array(sumed_values[1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(sumed_values[0], np.array(sumed_values[1])-np.array(sumed_values[2]), np.array(sumed_values[1])+np.array(sumed_values[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
# ax1.axvline(x = 0.6, color = "red")
# ax1.axvline(x = 3.5, color = "red")
ax1.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
