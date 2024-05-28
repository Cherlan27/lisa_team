# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 11:54:25 2024

@author: Lukas
"""

from d2_plotter import D2_Plotter
from scan_plotter import Timing_scan
from integrator import Integrator
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


dataframe = pd.read_excel("H:/supex623/supex623/LTP_2023_11/3mol_fast_time.xlsx", sheet_name= "05mol")
scan_numbers = [4141]
experiments = ["nai3moltiming6"]

def gaussian_double(data, b0, amp1, pos1, amp2, pos2, sig):
    return b0 + amp1* np.exp(-(data-pos1)**2/(2*sig**2)) + amp2* np.exp(-(data-pos2)**2/(2*sig**2))

values = []

for i,j in enumerate(list(zip(scan_numbers, experiments))):
    scan = j[0]
    experiment = j[1]
    roi = (1455,126,60,21)
    print(scan)
    fast_scan_integrator = Integrator()
    fast_scan_integrator.image_number = 0
    fast_scan_integrator.experiment = experiment
    fast_scan_integrator.scan_number = scan
    fast_scan_integrator.xlims = [roi[0], roi[0]+roi[2]]
    fast_scan_integrator.ylims = [roi[1], roi[1]+roi[3]]
    fast_scan_integrator()
    values.append(fast_scan_integrator.output_integration_values_y())
    
    for i in range(1, int(fast_scan_integrator.scan_cmd.split(" ")[-2])+1):
        print(i)
        fast_scan_integrator.image_number = i
        fast_scan_integrator()
        values.append(fast_scan_integrator.output_integration_values_y())


range_scan = np.linspace(float(fast_scan_integrator.scan_cmd.split(" ")[-4]), float(fast_scan_integrator.scan_cmd.split(" ")[-3]), int(fast_scan_integrator.scan_cmd.split(" ")[-2])+1)

fig = plt.figure(dpi = 300, figsize = (8, 20))
ax = fig.gca()
for i,scan in enumerate(values):
        plt.plot(scan[0], scan[1] + 2000*i, color = "blue")
        time_pos = range_scan[i]*10**(9)
        plt.text(1460, 2000*i + 300, "{0:.2f} ns".format(time_pos))
ax.set_ylabel("Intensity")
ax.set_xlabel("Pixel")
plt.show()

fit_value_list = []
err_value_list = []

fig = plt.figure(dpi = 300, figsize = (8, 20))
ax = fig.gca()
for i,scan in enumerate(values):
        plt.plot(scan[0], scan[1] + 2000*i, color = "blue")
        fit_values, corr_values = curve_fit(gaussian_double, scan[0], scan[1], p0=[0,2000, 1482, 2000, 1491, 4])
        err_values = np.sqrt(np.diag(abs(corr_values)))
        fit_value_list.append(fit_values)
        err_value_list.append(err_values)
        plt.plot(scan[0], gaussian_double(scan[0], *fit_values) + 2000*i, color = "green", linestyle = "--")
        time_pos = range_scan[i]*10**(9)
        plt.text(1460, 2000*i + 300, "{0:.2f} ns".format(time_pos))
plt.axvline(x = fit_value_list[0][2], color = "orange")
plt.axvline(x = fit_value_list[0][4], color = "orange")
ax.set_ylabel("Intensity")
ax.set_xlabel("Pixel")
plt.show()

# Amp and pos
fig, ax = plt.subplots(3, 2, dpi = 300, figsize = (10,10))
ax[0,0].errorbar(range_scan*10**(9), [x[1] for x in fit_value_list], yerr =  [x[1] for x in err_value_list], marker = "o", linestyle = "--", markersize = 3)
ax[0,1].errorbar(range_scan*10**(9), [x[3] for x in fit_value_list], yerr =  [x[3] for x in err_value_list], marker = "o", linestyle = "--", markersize = 3)
ax[1,0].errorbar(range_scan*10**(9), [x[2] for x in fit_value_list], yerr =  [x[2] for x in err_value_list], marker = "o", linestyle = "--", markersize = 3)
ax[1,1].errorbar(range_scan*10**(9), [x[4] for x in fit_value_list], yerr =  [x[4] for x in err_value_list], marker = "o", linestyle = "--", markersize = 3)
ax[2,0].errorbar(range_scan*10**(9), [x[5] for x in fit_value_list], yerr =  [x[5] for x in err_value_list], marker = "o", linestyle = "--", markersize = 3)
ax[2,1].errorbar(range_scan*10**(9), [x[5] for x in fit_value_list], yerr =  [x[5] for x in err_value_list], marker = "o", linestyle = "--", markersize = 3)
ax[0,0].set_title("Peak 1 Amplitude")
ax[0,1].set_title("Peak 2 Amplitude")
ax[1,0].set_title("Peak 1 Position")
ax[1,1].set_title("Peak 2 Position")
ax[2,0].set_title("Common Peak Sigma")
ax[2,1].set_title("Common Peak Sigma")
ax[2,0].set_xlabel("time in ns")
ax[2,1].set_xlabel("time in ns")
ax[0,0].set_ylabel("Amplitude")
ax[1,0].set_ylabel("Position")
ax[2,0].set_ylabel("Sigma")

plt.tight_layout()
plt.show()

