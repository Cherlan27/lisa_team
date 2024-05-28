# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 15:48:51 2024

@author: Lukas
"""

from amptek_plotter import Amptek_plotter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import preprocessing as pre

plt.rcParams['text.usetex'] = True


amptek = Amptek_plotter()
amptek.experiment = "nai0p1mol"
amptek.scan_number_off = [3896]
amptek.scan_number_on = [3917]
amptek()
amptek.plot_signal()
amptek.plot_fitting()

signal_01mol = amptek.get_signal()
signal_fit_01mol = amptek.get_fitting()

amptek = Amptek_plotter()
amptek.experiment = "nai0p5mol"
amptek.scan_number_off = [3844]
amptek.scan_number_on = [3857]
amptek()
amptek.plot_signal()
amptek.plot_fitting()

signal_05mol = amptek.get_signal()
signal_fit_05mol = amptek.get_fitting()

amptek = Amptek_plotter()
amptek.experiment = "nai1mol"
amptek.scan_number_off = [3013]
amptek.scan_number_on = [3026]
amptek()
amptek.plot_signal()
amptek.plot_fitting()

signal_1mol = amptek.get_signal()
signal_fit_1mol = amptek.get_fitting()

amptek = Amptek_plotter()
amptek.experiment = "nai3mol"
amptek.scan_number_off = [2846]
amptek.scan_number_on = [2862]
amptek()
amptek.plot_signal()
amptek.plot_fitting()

signal_3mol = amptek.get_signal()
signal_fit_3mol = amptek.get_fitting()

amptek = Amptek_plotter()
amptek.experiment = "nai5molfluorescence"
amptek.scan_number_off = [3991]
amptek.scan_number_on = [4011]
amptek()
amptek.plot_signal()
amptek.plot_fitting()

signal_5mol = amptek.get_signal()
signal_fit_5mol = amptek.get_fitting()

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5,1, figsize = (4,12), sharex = True)
fig.tight_layout()
plt.subplots_adjust(hspace=0)
ax1.plot(signal_fit_01mol[0], signal_fit_01mol[2], marker = "o", color = "blue", linestyle = "--", markerfacecolor = "white")
ax1.plot(signal_fit_01mol[1], signal_fit_01mol[3], marker = "o", color = "red", linestyle = "--", markerfacecolor = "white")
ax1.set_ylim([0.1,3.5])
ax2.plot(signal_fit_05mol[0], signal_fit_05mol[2], marker = "o", color = "blue", linestyle = "--", markerfacecolor = "white")
ax2.plot(signal_fit_05mol[1], signal_fit_05mol[3], marker = "o", color = "red", linestyle = "--", markerfacecolor = "white")
ax2.set_ylim([0.1,3.5])
ax3.plot(signal_fit_1mol[0], signal_fit_1mol[2], marker = "o", color = "blue", linestyle = "--", markerfacecolor = "white")
ax3.plot(signal_fit_1mol[1], signal_fit_1mol[3], marker = "o", color = "red", linestyle = "--", markerfacecolor = "white")
ax3.set_ylim([0.1,3.5])
ax4.plot(signal_fit_3mol[0], signal_fit_3mol[2], marker = "o", color = "blue", linestyle = "--", markerfacecolor = "white")
ax4.plot(signal_fit_3mol[1], signal_fit_3mol[3], marker = "o", color = "red", linestyle = "--", markerfacecolor = "white")
ax4.set_ylim([0.1,3.5])
ax5.plot(signal_fit_5mol[0], signal_fit_5mol[2], marker = "o", color = "blue", linestyle = "--", markerfacecolor = "white")
ax5.plot(signal_fit_5mol[1], signal_fit_5mol[3], marker = "o", color = "red", linestyle = "--", markerfacecolor = "white")
ax5.set_ylim([0.1,3.5])
ax3.set_ylabel("Intensity in arb. u.")
ax5.set_xlabel(r"$q_z$ in $\mathring{A}^{-1}$")

plt.plot()

plt.rcParams['font.size'] = 24
height_plot = 0.75
fig = plt.figure(dpi = 300, figsize = (8, 20))
ax = fig.gca()
for i in range(len(signal_3mol[1])):
    plt.plot(signal_3mol[0], signal_3mol[1][i] + height_plot*i, color = "blue")
    plt.plot(signal_3mol[0], signal_3mol[2][i] + height_plot*i, color = "red")
ax.set_ylabel("Intensity in arb. u.")
ax.set_xlabel("Energy in keV")
ax.set_xlim([3,4.5])
plt.show()


plt.rcParams["font.size"] = 8

cm = 1/2.54
fig = plt.figure(figsize = (7.25*cm, 4.65*cm), dpi = 600)

ax1 = fig.gca()
ax1.plot(signal_fit_3mol[0], signal_fit_3mol[2], marker = "o", color = "blue", linestyle = "--", markerfacecolor = "white", markersize = 4.5)
ax1.plot(signal_fit_3mol[1], signal_fit_3mol[3], marker = "o", color = "red", linestyle = "--", markerfacecolor = "white", markersize = 4.5)
ax1.set_ylabel("intensity (a. U.)")
ax1.set_xlabel(r"$q_z$ ($\mathring{A}^{-1}$)")
plt.tight_layout()
plt.show()