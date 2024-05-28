# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 14:12:59 2024

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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.optimize import curve_fit
from scipy.special import erf

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
    

data = pd.read_excel("H:/supex623/supex623/LTP_2023_11/3mol_fast_time.xlsx", sheet_name= "3mol")
data["Timeranges_no"] = [[x.split(" ")  [0], x.split(" ")[2]] if type(x) is str else np.nan for x in data["Timescale"] ]
data["Timerange_start"] = [list(map(float, (map(lambda x: str.replace(x, ",", "."), x))))[0] if type(x) is list else np.nan for x in data["Timeranges_no"]]
data["Timerange_end"] = [list(map(float, (map(lambda x: str.replace(x, ",", "."), x))))[1] if type(x) is list else np.nan for x in data["Timeranges_no"]]

# All scans in the range of -0.5 to 0.5 ns
scan_numbers = data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142)]["Scan_number"].values.tolist()
experiments =  data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142)]["experiment"].values.tolist()
params = []
first_run = True

for scan_id,experiment_id in list(zip(scan_numbers, experiments)):
    roi = (1455,126,60,21)
    
    # Script loading
    fast_scan = Timing_scan()
    fast_scan.scan_number = scan_id
    fast_scan.experiment = experiment_id
    fast_scan.roi = roi
    fast_scan()
    fast_scan.plot_intensity()
    print(scan_id)
    params.append(fast_scan.output_intensity())
    if first_run == True:
        params_sum = params[0]
        first_run = False
    else:
        concated_params = concatenate(params_sum, params[-1])
        

def exp_invfunc( xData, c, tau, b0):
    """
    Exponential model function (Version 1)
    
    Parameters
    ----------
    xData: Array
        x-axis data 
    c: float
        Factor in front of exp
    tau: float
        Relaxation time factor
    b0: float
        x-axis starting parameter   
        
    Result
    ----------
    Array
        y-values for the exponential function
    """
    return c * (1 - np.exp(-xData / tau)) + b0

def exp_invfunc2(xData, c, tau, b0):
    """
    Exponential model function (Version 2)
    
    Parameters
    ----------
    xData: Array
        x-axis data 
    c: float
        Factor in front of exp
    tau: float
        Relaxation time factor
    b0: float
        x-axis starting parameter  

    Result
    ----------
    Array
        y-values for the exponential function            
    """
    return ((c * (1 - np.exp(-xData / tau))) - b0)*(-1)
        
    
def error(xData, amp, mu, sigma, b0):
    """
    Exponential model function (Version 2)
    
    Parameters
    ----------
    xData: Array
        x-axis data 
    c: float
        Factor in front of exp
    tau: float
        Relaxation time factor
    b0: float
        x-axis starting parameter  

    Result
    ----------
    Array
        y-values for the exponential function            
    """
    return amp*erf((xData - mu)/(sigma*np.sqrt(2))) + b0

def error_2(xData, amp, mu, sigma):
    """
    Exponential model function (Version 2)
    
    Parameters
    ----------
    xData: Array
        x-axis data 
    c: float
        Factor in front of exp
    tau: float
        Relaxation time factor
    b0: float
        x-axis starting parameter  

    Result
    ----------
    Array
        y-values for the exponential function            
    """
    return amp*erf((xData - mu)/(sigma*np.sqrt(2))) + 9.205 # 9.25


plt.rcParams["font.size"] = 8

cm = 1/2.54
fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12.85*cm, 5.65*cm), dpi = 600, sharex=True)
ax1.errorbar(np.array(params[0][0]), np.array(params[0][1])*1000, yerr = np.array(params[0][2])*1000, color='#11570d', ls='--', marker='o', mec='#424242',mfc='white', mew=0.7)
ax1.fill_between(np.array(params[0][0]), np.array(params[0][1])*1000-np.array(params[0][2])*1000, np.array(params[0][1])*1000+np.array(params[0][2])*1000, facecolor='#11570d', alpha=0.3, edgecolor='none')
ax1.set_ylabel("intensity (a.U.)")
ax1.set_xlabel("time (ns)")

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.tick_params(which='both', width=1)
ax1.tick_params(which='minor', length=2)
ax1.tick_params(which='both', direction="in")
ax1.tick_params(which='both', top = True, right = True)
ax1.set_xlim([-4, 10])
ax1.set_ylim([9.1, 9.7])

ax2.errorbar(np.array(params[1][0]), np.array(params[1][1])*1000, yerr = np.array(params[1][2])*1000, color='#11570d', ls='--', marker='o', mec='#424242',mfc='white', mew=0.7)
ax2.fill_between(np.array(params[1][0]), np.array(params[1][1])*1000-np.array(params[1][2])*1000, np.array(params[1][1])*1000+np.array(params[1][2])*1000, facecolor='#11570d', alpha=0.3, edgecolor='none')
ax2.set_xlabel("time (ns)")

ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

ax2.tick_params(which='both', width=1)
ax2.tick_params(which='minor', length=2)
ax2.tick_params(which='both', direction="in")
ax2.tick_params(which='both', top = True, right = True)
ax2.set_yticklabels([])

ax2.set_xlim([-4, 10])
ax2.set_ylim([9.1, 9.7])

plt.subplots_adjust(wspace=0.1, hspace=0)
plt.show()


the_shfit = 0.4
cm = 1/2.54
final_1 = 8
final_2 = 30

time_axis = np.array(concated_params[0]) + the_shfit

params_01, params_corr_01 = curve_fit(error, np.array(time_axis[0:9]), np.array(concated_params[1][0:9])*1000, p0 = [-0.2, 0,0.4,9.2])
err_01 = np.sqrt(np.diag(params_corr_01))
params_02, params_corr_02 = curve_fit(error_2, np.array(time_axis[final_1:final_2]), np.array(concated_params[1][final_1:final_2])*1000, p0 = [0.2, 0.6, 4], bounds = ([-1, 0.5, 0],[1, 5, 6]), maxfev = 50000)
err_02 = np.sqrt(np.diag(params_corr_02))

fig, (ax1) = plt.subplots(1, 1, figsize = (7.25*cm, 4.65*cm), dpi = 600, sharex=True)
ax1.errorbar(np.array(time_axis), np.array(concated_params[1])*1000, yerr = np.array(concated_params[2])*1000, color='#11570d', ls='--', marker='o', mec='#424242',mfc='white', mew=0.7)

ax1.errorbar(np.array(time_axis[0:9]), error(np.array(time_axis[0:9]), *params_01), color = "red")
ax1.errorbar(np.array(time_axis[final_1:final_2]), error_2(np.array(time_axis[final_1:final_2]), *params_02), color = "red")


ax1.fill_between(np.array(time_axis), np.array(concated_params[1])*1000-np.array(concated_params[2])*1000, np.array(concated_params[1])*1000+np.array(concated_params[2])*1000, facecolor='#11570d', alpha=0.3, edgecolor='none')
ax1.set_ylabel("intensity (a.U.)")
ax1.set_xlabel("time (ns)")

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.tick_params(which='both', width=1)
ax1.tick_params(which='minor', length=2)
ax1.tick_params(which='both', direction="in")
ax1.tick_params(which='both', top = True, right = True)
#ax1.set_ylim([9.1, 9.7])

plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

print(f"Zeitfaktor:  {params_01[2]} +/- {err_01[2]} ns")
print(f"Zeitfaktor:  {params_02[2]} +/- {err_02[2]} ns")
print(params_02)


the_shfit = 0.55
cm = 1/2.54
final_1 = 8
final_2 = 20

time_axis = np.array(concated_params[0]) + the_shfit

params_01, params_corr_01 = curve_fit(error, np.array(time_axis[0:9]), np.array(concated_params[1][0:9])*1000, p0 = [-0.2, 0,0.4,9.2])
err_01 = np.sqrt(np.diag(params_corr_01))
params_02, params_corr_02 = curve_fit(error_2, np.array(time_axis[final_1:final_2]), np.array(concated_params[1][final_1:final_2])*1000, p0 = [0.2, 0.6, 4], bounds = ([-1, 0.5, 0],[1, 5, 6]), maxfev = 50000)
err_02 = np.sqrt(np.diag(params_corr_02))

fig, (ax1) = plt.subplots(1, 1, figsize = (7.25*cm, 5.65*cm), dpi = 600, sharex=True)
ax1.errorbar(np.array(time_axis), np.array(concated_params[1])*1000, yerr = np.array(concated_params[2])*1000, color='#11570d', ls='--', marker='o', mec='#424242',mfc='white', mew=0.7)

ax1.errorbar(np.array(time_axis[0:9]), error(np.array(time_axis[0:9]), *params_01), color = "red")
ax1.errorbar(np.array(time_axis[final_1:final_2]), error_2(np.array(time_axis[final_1:final_2]), *params_02), color = "red")


ax1.fill_between(np.array(time_axis), np.array(concated_params[1])*1000-np.array(concated_params[2])*1000, np.array(concated_params[1])*1000+np.array(concated_params[2])*1000, facecolor='#11570d', alpha=0.3, edgecolor='none')
ax1.set_ylabel("intensity (a.U.)")
ax1.set_xlabel("time (ns)")

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.tick_params(which='both', width=1)
ax1.tick_params(which='minor', length=2)
ax1.tick_params(which='both', direction="in")
ax1.tick_params(which='both', top = True, right = True)
#ax1.set_ylim([9.1, 9.7])

ax1.legend(["Fit"])

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

print(f"Zeitfaktor:  {params_01[2]} +/- {err_01[2]} ns")
print(f"Zeitfaktor:  {params_02[2]} +/- {err_02[2]} ns")
print(params_02)