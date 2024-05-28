# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 15:00:28 2024

@author: Lukas
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['text.usetex'] = True


file = "H:/supex623/supex623/LTP_2023_11"

file_5mol = "{0}/5mol_reflectivity_fit_params.txt".format(file)
file_3mol = "{0}/3mol_reflectivity_fit_params.txt".format(file)
file_1mol = "{0}/1mol_reflectivity_fit_params.txt".format(file)
file_05mol = "{0}/0p5mol_reflectivity_fit_params.txt".format(file)
file_01mol = "{0}/0p1mol_reflectivity_fit_params.txt".format(file)

pd_5mol = pd.read_csv(file_5mol, sep = "\t")
pd_3mol = pd.read_csv(file_3mol, sep = "\t")
pd_1mol = pd.read_csv(file_1mol, sep = "\t")
pd_05mol = pd.read_csv(file_05mol, sep = "\t")
pd_01mol = pd.read_csv(file_01mol, sep = "\t")

print(pd_5mol)

pp_off_err = [pd_01mol["sld_error"][0], pd_05mol["sld_error"][0], pd_1mol["sld_error"][0], pd_3mol["sld_error"][0], pd_5mol["sld_error"][0]]
pp_900_err = [pd_01mol["sld_error"][1], pd_05mol["sld_error"][1], pd_1mol["sld_error"][1], pd_3mol["sld_error"][1], pd_5mol["sld_error"][1]]
pp_1200_err = [pd_01mol["sld_error"][3], pd_05mol["sld_error"][3], pd_1mol["sld_error"][2], pd_3mol["sld_error"][3], pd_5mol["sld_error"][2]]
pp_1500_err = [pd_1mol["sld_error"][3], pd_5mol["sld_error"][3]]

pp_off = [pd_01mol["sld"][0], pd_05mol["sld"][0], pd_1mol["sld"][0], pd_3mol["sld"][0], pd_5mol["sld"][0]]
pp_900 = [pd_01mol["sld"][1], pd_05mol["sld"][1], pd_1mol["sld"][1], pd_3mol["sld"][1], pd_5mol["sld"][1]]
pp_1200 = [pd_01mol["sld"][3], pd_05mol["sld"][3], pd_1mol["sld"][2], pd_3mol["sld"][3], pd_5mol["sld"][2]]
pp_1500 = [pd_1mol["sld"][3], pd_5mol["sld"][3]]

conc = [0.1, 0.5, 1, 3, 5]
conc_1500 = [1, 5]

plt.rcParams['font.size'] = 16.5
fig = plt.figure(dpi = 300)
ax = fig.gca()
ax.errorbar(conc, pp_off, yerr = pp_off_err, color = "blue", marker = "o", linestyle = "")
ax.errorbar(conc, pp_1200, yerr = pp_1200_err, color = "red", marker = "o", linestyle = "")
ax.errorbar(conc_1500, pp_1500, yerr = pp_1500_err, color = "#AF2103", marker = "o", linestyle = "")
ax.set_xlabel("NaI concentration in Mol/l")
ax.set_ylabel(r"SLD in $10^{-6}$ $\mathring{A}^{-2}$")
ax.legend(["Laser off", r"$E_{puls}$: 1.7 $\mu$J", r"$E_{puls}$: 3.55 $\mu$J"])
plt.show()