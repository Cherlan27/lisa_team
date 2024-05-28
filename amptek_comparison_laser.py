# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 16:53:54 2024

@author: Lukas
"""

from amptek_plotter import Amptek_plotter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import preprocessing as pre

amptek = Amptek_plotter()
amptek.experiment = "nai3mol"
amptek.scan_number_off = [2846]
amptek.scan_number_on = [2852]
amptek()
amptek.plot_signal()
amptek.plot_fitting()

signal_3mol_900 = amptek.get_signal()
signal_fit_3mol_900 = amptek.get_fitting()

amptek = Amptek_plotter()
amptek.experiment = "nai3mol"
amptek.scan_number_off = [2846]
amptek.scan_number_on = [2861]
amptek()
amptek.plot_signal()
amptek.plot_fitting()

signal_3mol_1200 = amptek.get_signal()
signal_fit_3mol_1200 = amptek.get_fitting()

amptek = Amptek_plotter()
amptek.experiment = "nai3mol"
amptek.scan_number_off = [2846]
amptek.scan_number_on = [2867]
amptek()
amptek.plot_signal()
amptek.plot_fitting()

signal_3mol_1500 = amptek.get_signal()
signal_fit_3mol_1500 = amptek.get_fitting()

fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize = (6,12), sharex = True)
fig.tight_layout()
plt.subplots_adjust(hspace=0)
ax1.plot(signal_fit_3mol_900[0], signal_fit_3mol_900[2], marker = "o", color = "blue", linestyle = "--", markerfacecolor = "white")
ax1.plot(signal_fit_3mol_900[1], signal_fit_3mol_900[3], marker = "o", color = "red", linestyle = "--", markerfacecolor = "white")
ax2.plot(signal_fit_3mol_1200[0], signal_fit_3mol_1200[2], marker = "o", color = "blue", linestyle = "--", markerfacecolor = "white")
ax2.plot(signal_fit_3mol_1200[1], signal_fit_3mol_1200[3], marker = "o", color = "red", linestyle = "--", markerfacecolor = "white")
ax3.plot(signal_fit_3mol_1500[0], signal_fit_3mol_1500[2], marker = "o", color = "blue", linestyle = "--", markerfacecolor = "white")
ax3.plot(signal_fit_3mol_1500[1], signal_fit_3mol_1500[3], marker = "o", color = "red", linestyle = "--", markerfacecolor = "white")
ax1.set_ylabel("Intensity 900")
ax3.set_xlabel("Energy (in keV)")

plt.plot()