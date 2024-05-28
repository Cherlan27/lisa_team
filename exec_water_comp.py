# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 16:01:43 2024

@author: Lukas
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 15:01:16 2024

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

xrr_analysis = data_analysis("h2o", range(1320, 1328))
refnx_xrr1, ml_xrr1 = xrr_analysis()
xrr_analysis = data_analysis("nai1mol", range(3069, 3077))
refnx_xrr2, ml_xrr2 = xrr_analysis()
# xrr_analysis = data_analysis("nai3molxrr", range(3324, 3333))
# refnx_xrr3, ml_xrr3 = xrr_analysis()
# xrr_analysis = data_analysis("nai5molxrr", range(3227, 3236))
# refnx_xrr4, ml_xrr4 = xrr_analysis()

cm = 1/2.54
plt.rcParams["font.size"] = 8

fig = plt.figure(figsize = (8.85*cm, 4.65*cm), dpi = 600)
ax = fig.gca()
plt.tight_layout()
ax.errorbar(refnx_xrr2[0], refnx_xrr2[1], color = "#228b22", marker = "o",mfc= "none", linestyle = "", markersize = 2)
ax.plot(refnx_xrr2[0], refnx_xrr2[3], color = "#228b22", marker = "", linestyle = "-")
ax.errorbar(refnx_xrr1[0], refnx_xrr1[1], color = "#ffa54f", marker = "o",mfc= "none", linestyle = "", markersize = 2)
ax.plot(refnx_xrr1[0], refnx_xrr1[3], color = "#ffa54f", marker = "", linestyle = "-")
# ax.errorbar(refnx_xrr3[0], refnx_xrr3[1], refnx_xrr3[2], color = "#FA5858", marker = "o",mfc='none', linestyle = "")
# ax.plot(refnx_xrr3[0], refnx_xrr3[3], color = "#FA5858", marker = "", linestyle = "-")
# ax.errorbar(refnx_xrr4[0], refnx_xrr4[1], refnx_xrr4[2], color = "#DF0101", marker = "o",mfc='none', linestyle = "")
# ax.plot(refnx_xrr4[0], refnx_xrr4[3], color = "#DF0101", marker = "", linestyle = "-")
ax.legend(["Good alignment","Misalignment"])
plt.xlabel(r'$q_z$ ($\mathring{A}^{-1}$)')
plt.ylabel(r'R/R$_F$')
plt.yscale("log")
plt.show()
