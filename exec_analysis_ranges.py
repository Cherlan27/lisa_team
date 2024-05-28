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
data["Timeranges_no"] = [[x.split(" ")[0], x.split(" ")[2]] if type(x) is str else np.nan for x in data["Timescale"] ]
data["Timerange_start"] = [list(map(float, (map(lambda x: str.replace(x, ",", "."), x))))[0] if type(x) is list else np.nan for x in data["Timeranges_no"]]
data["Timerange_end"] = [list(map(float, (map(lambda x: str.replace(x, ",", "."), x))))[1] if type(x) is list else np.nan for x in data["Timeranges_no"]]

# All scans in the range of -0.5 to 0.5 ns
scan_numbers = data[(data["Timerange_start"] == -0.5) & (data["Timerange_end"] == 0.5) & (data["Usable"] == "y")]["Scan_number"].values.tolist()
experiments =  data[(data["Timerange_start"] == -0.5) & (data["Timerange_end"] == 0.5) & (data["Usable"] == "y")]["experiment"].values.tolist()
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
        
fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(concated_params[0]), np.array(concated_params[1]), yerr = np.array(concated_params[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(concated_params[0]), baseline_als(np.array(concated_params[1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(np.array(concated_params[0]), np.array(concated_params[1])-np.array(concated_params[2]), np.array(concated_params[1])+np.array(concated_params[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
ax1.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()


# All scans in the range of -0.5 to 0.5 ns
data = pd.read_excel("H:/supex623/supex623/LTP_2023_11/3mol_fast_time.xlsx", sheet_name= "3mol")
scan_numbers = data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142) | (data["Scan_number"] == 4149) | (data["Scan_number"] == 4171)| (data["Scan_number"] == 4172) | (data["Scan_number"] == 4178)]["Scan_number"].values.tolist()
experiments =  data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142) | (data["Scan_number"] == 4149) | (data["Scan_number"] == 4171)| (data["Scan_number"] == 4172) | (data["Scan_number"] == 4178)]["experiment"].values.tolist()
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
        
fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(concated_params[0]), np.array(concated_params[1]), yerr = np.array(concated_params[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(concated_params[0]), baseline_als(np.array(concated_params[1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(np.array(concated_params[0]), np.array(concated_params[1])-np.array(concated_params[2]), np.array(concated_params[1])+np.array(concated_params[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
ax1.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()



data = pd.read_excel("H:/supex623/supex623/LTP_2023_11/3mol_fast_time.xlsx", sheet_name= "5mol")
data["Timeranges_no"] = [[x.split(" ")  [0], x.split(" ")[2]] if type(x) is str else np.nan for x in data["Timescale"] ]
data["Timerange_start"] = [list(map(float, (map(lambda x: str.replace(x, ",", "."), x))))[0] if type(x) is list else np.nan for x in data["Timeranges_no"]]
data["Timerange_end"] = [list(map(float, (map(lambda x: str.replace(x, ",", "."), x))))[1] if type(x) is list else np.nan for x in data["Timeranges_no"]]

# All scans in the range of -0.5 to 0.5 ns
scan_numbers = data[(data["Scan_number"] == 2541) | (data["Scan_number"] == 2542) | (data["Scan_number"] == 2549) | (data["Scan_number"] == 2550)| (data["Scan_number"] == 2554) | (data["Scan_number"] == 2555) | (data["Scan_number"] == 2561) | (data["Scan_number"] == 2562)]["Scan_number"].values.tolist()
experiments =  data[(data["Scan_number"] == 2541) | (data["Scan_number"] == 2542) | (data["Scan_number"] == 2549) | (data["Scan_number"] == 2550)| (data["Scan_number"] == 2554) | (data["Scan_number"] == 2555) | (data["Scan_number"] == 2561) | (data["Scan_number"] == 2562)]["experiment"].values.tolist()
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
        
fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(concated_params[0]), np.array(concated_params[1]), yerr = np.array(concated_params[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(concated_params[0]), baseline_als(np.array(concated_params[1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(np.array(concated_params[0]), np.array(concated_params[1])-np.array(concated_params[2]), np.array(concated_params[1])+np.array(concated_params[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
ax1.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

data = pd.read_excel("H:/supex623/supex623/LTP_2023_11/3mol_fast_time.xlsx", sheet_name= "05mol")
data["Timeranges_no"] = [[x.split(" ")  [0], x.split(" ")[2]] if type(x) is str else np.nan for x in data["Timescale"] ]
data["Timerange_start"] = [list(map(float, (map(lambda x: str.replace(x, ",", "."), x))))[0] if type(x) is list else np.nan for x in data["Timeranges_no"]]
data["Timerange_end"] = [list(map(float, (map(lambda x: str.replace(x, ",", "."), x))))[1] if type(x) is list else np.nan for x in data["Timeranges_no"]]

# All scans in the range of -0.5 to 0.5 ns
scan_numbers = data[(data["Scan_number"] == 4241) | (data["Scan_number"] == 4242)]["Scan_number"].values.tolist()
experiments =  data[(data["Scan_number"] == 4241) | (data["Scan_number"] == 4242)]["experiment"].values.tolist()
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
        
fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(concated_params[0]), np.array(concated_params[1]), yerr = np.array(concated_params[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(concated_params[0]), baseline_als(np.array(concated_params[1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(np.array(concated_params[0]), np.array(concated_params[1])-np.array(concated_params[2]), np.array(concated_params[1])+np.array(concated_params[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
ax1.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

# All scans in the range of -0.5 to 0.5 ns
scan_numbers = data[(data["Scan_number"] == 4160) | (data["Scan_number"] == 4161)]["Scan_number"].values.tolist()
experiments =  data[(data["Scan_number"] == 4160) | (data["Scan_number"] == 4161)]["experiment"].values.tolist()
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
        
fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
ax1.errorbar(np.array(concated_params[0]), np.array(concated_params[1]), yerr = np.array(concated_params[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
ax1.plot(np.array(concated_params[0]), baseline_als(np.array(concated_params[1]),0.2, 0.5), linestyle = "--", color='#11570d')
ax1.fill_between(np.array(concated_params[0]), np.array(concated_params[1])-np.array(concated_params[2]), np.array(concated_params[1])+np.array(concated_params[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
ax1.set_xlabel("t in ns")
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()


# All scans in the range of -5 to 10 ns
# scan_numbers = data[(data["Timerange_start"] == -5) & (data["Timerange_end"] == 10) & (data["Usable"] == "y")]["Scan_number"].values.tolist()
# experiments =  data[(data["Timerange_start"] == -5) & (data["Timerange_end"] == 10) & (data["Usable"] == "y")]["experiment"].values.tolist()

# params = []
# first_run = True

# for scan_id,experiment_id in list(zip(scan_numbers, experiments)):
#     roi = (1455,126,60,21)
    
#     # Script loading
#     fast_scan = Timing_scan()
#     fast_scan.scan_number = scan_id
#     fast_scan.experiment = experiment_id
#     fast_scan.roi = roi
#     fast_scan()
#     fast_scan.plot_intensity()
#     print(scan_id)
#     params.append(fast_scan.output_intensity())
#     if first_run == True:
#         params_sum = params[0]
#         first_run = False
#     else:
#         concated_params = concatenate(params_sum, params[-1])
        
# fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
# ax1.errorbar(np.array(concated_params[0]), np.array(concated_params[1]), yerr = np.array(concated_params[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
# ax1.plot(np.array(concated_params[0]), baseline_als(np.array(concated_params[1]),0.2, 0.5), linestyle = "--", color='#11570d')
# ax1.fill_between(np.array(concated_params[0]), np.array(concated_params[1])-np.array(concated_params[2]), np.array(concated_params[1])+np.array(concated_params[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
# ax1.set_xlabel("t in ns")
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.show()


# All scans in the range of -5 to 10 ns
# scan_numbers = data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142)]["Scan_number"].values.tolist()
# experiments =  data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142)]["experiment"].values.tolist()

# params = []
# first_run = True

# for scan_id,experiment_id in list(zip(scan_numbers, experiments)):
#     roi = (1455,126,60,21)
    
#     # Script loading
#     fast_scan = Timing_scan()
#     fast_scan.scan_number = scan_id
#     fast_scan.experiment = experiment_id
#     fast_scan.roi = roi
#     fast_scan()
#     fast_scan.plot_intensity()
#     print(scan_id)
#     params.append(fast_scan.output_intensity())
#     if first_run == True:
#         params_sum = params[0]
#         first_run = False
#     else:
#         concated_params = concatenate(params_sum, params[-1])
        
# fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
# ax1.errorbar(np.array(concated_params[0]), np.array(concated_params[1]), yerr = np.array(concated_params[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
# ax1.plot(np.array(concated_params[0]), baseline_als(np.array(concated_params[1]),0.2, 0.5), linestyle = "--", color='#11570d')
# ax1.fill_between(np.array(concated_params[0]), np.array(concated_params[1])-np.array(concated_params[2]), np.array(concated_params[1])+np.array(concated_params[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
# ax1.set_xlabel("t in ns")
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.show()


# All scans in the range of -5 to 10 ns
# scan_numbers = data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142) | (data["Scan_number"] == 4173) | (data["Scan_number"] == 4171)]["Scan_number"].values.tolist()
# experiments =  data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142) | (data["Scan_number"] == 4173) | (data["Scan_number"] == 4171)]["experiment"].values.tolist()

# params = []
# first_run = True

# for scan_id,experiment_id in list(zip(scan_numbers, experiments)):
#     roi = (1455,126,60,21)
    
#     # Script loading
#     fast_scan = Timing_scan()
#     fast_scan.scan_number = scan_id
#     fast_scan.experiment = experiment_id
#     fast_scan.roi = roi
#     fast_scan()
#     fast_scan.plot_intensity()
#     print(scan_id)
#     params.append(fast_scan.output_intensity())
#     if first_run == True:
#         params_sum = params[0]
#         first_run = False
#     else:
#         concated_params = concatenate(params_sum, params[-1])
        
# fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
# ax1.errorbar(np.array(concated_params[0]), np.array(concated_params[1]), yerr = np.array(concated_params[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
# ax1.plot(np.array(concated_params[0]), baseline_als(np.array(concated_params[1]),0.2, 0.5), linestyle = "--", color='#11570d')
# ax1.fill_between(np.array(concated_params[0]), np.array(concated_params[1])-np.array(concated_params[2]), np.array(concated_params[1])+np.array(concated_params[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
# ax1.set_xlabel("t in ns")
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.show()



# All scans in the range of -5 to 10 ns
# scan_numbers = data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142) | (data["Scan_number"] == 4113)]["Scan_number"].values.tolist()
# experiments =  data[(data["Scan_number"] == 4141) | (data["Scan_number"] == 4142) | (data["Scan_number"] == 4113)]["experiment"].values.tolist()

# params = []
# first_run = True

# for scan_id,experiment_id in list(zip(scan_numbers, experiments)):
#     roi = (1455,126,60,21)
    
#     # Script loading
#     fast_scan = Timing_scan()
#     fast_scan.scan_number = scan_id
#     fast_scan.experiment = experiment_id
#     fast_scan.roi = roi
#     fast_scan()
#     fast_scan.plot_intensity()
#     print(scan_id)
#     params.append(fast_scan.output_intensity())
#     if first_run == True:
#         params_sum = params[0]
#         first_run = False
#     else:
#         concated_params = concatenate(params_sum, params[-1])
        
# fig, (ax1) = plt.subplots(1, 1, dpi = 300, sharex=True)
# ax1.errorbar(np.array(concated_params[0]), np.array(concated_params[1]), yerr = np.array(concated_params[2]), color='#11570d', ls='None', marker='o', mec='#424242',mfc='white', mew=1.2)
# ax1.plot(np.array(concated_params[0]), baseline_als(np.array(concated_params[1]),0.2, 0.5), linestyle = "--", color='#11570d')
# ax1.fill_between(np.array(concated_params[0]), np.array(concated_params[1])-np.array(concated_params[2]), np.array(concated_params[1])+np.array(concated_params[2]), facecolor='#11570d', alpha=0.3, edgecolor='none')
# ax1.set_xlabel("t in ns")
# plt.subplots_adjust(wspace=0, hspace=0)
plt.show()