# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 11:40:31 2024

@author: Lukas
"""

from d2_plotter import D2_Plotter
from scan_plotter import Timing_scan
from integrator import Integrator
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

scan_numbers = [4190, 4202, 4204]
experiments = ["nai3moltiming7", "nai3moltiming7", "nai3moltiming7"]

params = []

for i,j in enumerate(list(zip(scan_numbers, experiments))):
    scan = j[0]
    experiment = j[1]
    roi = (1455,126,60,21)
    
    print(scan)
    
    # Script loading
    fast_scan = Timing_scan()
    fast_scan.scan_number = scan
    fast_scan.experiment = experiment
    fast_scan.roi = roi
    fast_scan.save_directory = "H:/supex623/supex623/LTP_2023_11/05mol_fast_scans_2/"
    
    fast_scan()
    fast_scan.plot_intensity()
    params.append(fast_scan.output_intensity())
    # fast_scan_detector()
    # fast_scan_integrator()
    
cm = 1/2.54
fig, axs = plt.subplots(1, 3, dpi = 500, sharey= True, figsize = (14.85*cm, 5.65*cm))
axs[1].errorbar(np.array(params[0][0]), np.array(params[0][1]), yerr = np.array(params[0][2]), color='#11570d', ls='--', marker='o', mec='#424242',mfc='white', mew=0.7, markersize = 2)

axs[1].set_xlabel("time (s)")

axs[1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1].yaxis.set_minor_locator(AutoMinorLocator())

axs[1].tick_params(which='both', width=1)
axs[1].tick_params(which='minor', length=2)
axs[1].tick_params(which='both', direction="in")
axs[1].tick_params(which='both', top = True, right = True)
axs[1].axvline(x = 20, color = "green")
axs[1].axvline(x = 60, color = "red")
axs[1].text(110,0.49, "(b)")

axs[0].errorbar(np.array(params[1][0]), np.array(params[1][1]), yerr = np.array(params[1][2]), color='#11570d', ls='--', marker='o', mec='#424242',mfc='white', mew=0.7, markersize = 2)
axs[0].set_xlabel("time (s)")
axs[0].set_ylabel("intensity (a.U.)")
axs[0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0].yaxis.set_minor_locator(AutoMinorLocator())

axs[0].tick_params(which='both', width=1)
axs[0].tick_params(which='minor', length=2)
axs[0].tick_params(which='both', direction="in")
axs[0].tick_params(which='both', top = True, right = True)
axs[0].axvline(x = 20, color = "green")
axs[0].axvline(x = 60, color = "red")
axs[0].text(110,0.49, "(a)")

axs[2].errorbar(np.array(params[2][0]), np.array(params[2][1]), yerr = np.array(params[2][2]), color='#11570d', ls='--', marker='o', mec='#424242',mfc='white', mew=0.7, markersize = 2)
axs[2].set_xlabel("time (s)")

axs[2].xaxis.set_minor_locator(AutoMinorLocator())
axs[2].yaxis.set_minor_locator(AutoMinorLocator())

axs[2].tick_params(which='both', width=1)
axs[2].tick_params(which='minor', length=2)
axs[2].tick_params(which='both', direction="in")
axs[2].tick_params(which='both', top = True, right = True)
axs[2].axvline(x = 20, color = "green")
axs[2].axvline(x = 60, color = "red")
axs[2].text(110,0.49, "(c)")