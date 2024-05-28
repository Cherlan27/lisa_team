# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:20:28 2024

@author: Lukas
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 14:22:53 2024

@author: Lukas
"""

from d2_plotter import D2_Plotter
from scan_plotter import Timing_scan
from integrator import Integrator
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

x_axis_sums = []
y_axis_sums = []
imgs = []

scan_numbers = [1477,1477,1555]
image_number = [0,15,15]
experiment = ["align_18keV_2", "align_18keV_2", "align_18keV_2"]
data_directory = "K:/SYNCHROTRON/Murphy/2019-10_P08_11006759_giri_warias/raw_data/raw"

for i, (j, k, l) in enumerate(zip(scan_numbers, image_number, experiment)):
    roi = (1434,345,50,151)
    
    print(i, k)
    fast_scan_integrator = Integrator()
    fast_scan_integrator.image_number = k
    fast_scan_integrator.experiment = l
    fast_scan_integrator.data_directory = data_directory
    fast_scan_integrator.intensity_scale = (1, 60)
    fast_scan_integrator.scan_number = j
    fast_scan_integrator.xlims = [roi[0], roi[0]+roi[2]]
    fast_scan_integrator.ylims = [roi[1], roi[1]+roi[3]]

    fast_scan_integrator()
    imgs.append(fast_scan_integrator.output_img())
    x_axis_sums.append(fast_scan_integrator.output_integration_values_x())
    y_axis_sums.append(fast_scan_integrator.output_integration_values_y())

for i,j in enumerate([1]):
    scan = 3976
    experiment = "nai3moltiming5"
    image = j
    roi = (1455,126,60,21)
    
    print(scan)
    
    fast_scan_integrator = Integrator()
    fast_scan_integrator.image_number = image
    fast_scan_integrator.experiment = experiment
    fast_scan_integrator.intensity_scale = (1, 60)
    fast_scan_integrator.scan_number = scan
    fast_scan_integrator.xlims = [roi[0], roi[0]+roi[2]]
    fast_scan_integrator.ylims = [roi[1], roi[1]+roi[3]]

    fast_scan_integrator()
    imgs.append(fast_scan_integrator.output_img())
    x_axis_sums.append(fast_scan_integrator.output_integration_values_x())
    y_axis_sums.append(fast_scan_integrator.output_integration_values_y())
    
    
rotate_detector = True
intensity_scale = (1, 1000)
log_scale = True
xdim = "pixels"
ydim = "pixels"
Q = None,
colormap = "Blues"

y_axis_diff = []
y_axis_diff.append(np.array(list(y_axis_sums[0][0])) - y_axis_sums[0][0][25])
y_axis_diff.append(np.array(list(y_axis_sums[1][0])) - y_axis_sums[0][0][25])
y_axis_diff.append(np.array(list(y_axis_sums[2][0])) - y_axis_sums[0][0][25])
y_axis_diff.append(np.array(list(y_axis_sums[3][0])) - y_axis_sums[3][0][34])

x_axis_diff = []
x_axis_diff.append(np.array(list(x_axis_sums[0][0])) - x_axis_sums[0][0][8])
x_axis_diff.append(np.array(list(x_axis_sums[1][0])) - x_axis_sums[0][0][8])
x_axis_diff.append(np.array(list(x_axis_sums[2][0])) - x_axis_sums[0][0][8])
x_axis_diff.append(np.array(list(x_axis_sums[3][0])) - x_axis_sums[3][0][9])

roi = (1434,345,50,12)
xlims = [roi[0], roi[0]+roi[2]]
ylims = [roi[1], roi[1]+roi[3]]

    
mosaic = [['.', 'A panel',],['B panel','C panel']]
kw = dict(layout='constrained')

plt.rcParams["font.size"] = 10
plt.rcParams['lines.markersize'] = 1

cm = 1/2.54
fig = plt.figure(layout="constrained", dpi = 500, figsize=(8.85*cm,8.5*cm))
first_pic, second_pic, third_pic, fourth_pic = fig.subfigures(nrows=1, ncols=4)
fig.tight_layout()

axd_01 = first_pic.subplot_mosaic(mosaic, width_ratios = [0.25,1], height_ratios = [0.1,1])

ax_x_integrate = axd_01['A panel']
ax_x_integrate.set_xticklabels("")
ax_x_integrate.set_yticklabels("")
ax_x_integrate.tick_params(which='both', top = True, left = False, right = False, bottom = True, direction = "in", width=1)

ax_y_integrate = axd_01['B panel']
ax_y_integrate.set_xticklabels("")
ax_y_integrate.set_ylabel(r"$\Delta$z$_{det}$ (pixel)")
ax_y_integrate.tick_params(which='both', bottom = False, left = True, right = True, direction = "in", width=1)

ax_detector = axd_01['C panel']
ax_detector.set_yticklabels("")
#ax_detector.set_xlabel(r"$\Delta$y$_{det}$ (pixel)")
ax_detector.tick_params(which='both', top = True, left = True, right = True, bottom = True, direction = "in", width=1)

smin, smax = intensity_scale
if log_scale:
    data = np.log10(imgs[0])
    vmin = max(np.log10(smin), 0)
    vmax = np.log10(smax)
else:
    data = imgs[0]
    vmin, vmax = smin, smax

if rotate_detector == False:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[0].shape[1]),
                                  np.arange(imgs[0].shape[0]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[0].shape[1]),
                                  np.arange(imgs[0].shape[0]))
    # plot data
    ax_detector.pcolormesh(xdata - x_axis_sums[0][0][8],  ydata - y_axis_sums[0][0][25] , data,  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim([-8, 8])
        ax_detector.set_ylim([-25, 25])
    
else:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[0].shape[0]),
                                  np.arange(imgs[0].shape[1]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[0].shape[0]),
                                  np.arange(imgs[0].shape[1]))
    # plot data
    ax_detector.pcolormesh(xdata - x_axis_sums[0][0][8], ydata - y_axis_sums[0][0][25], np.transpose(data),  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim([-8, 8])
        ax_detector.set_ylim([-25, 25])

if rotate_detector == True:
    ax_x_integrate.plot(x_axis_diff[0], x_axis_sums[0][1], marker = 'o', linestyle = "None")
    ax_x_integrate.set_xlim([-8, 8])
    ax_x_integrate.text(-13, 2000, "a)")
    ax_x_integrate.text(13, 2000, "b)")
    ax_x_integrate.text(39, 2000, "c)")
    ax_x_integrate.text(64.5, 2000, "d)")

    ax_y_integrate.plot(y_axis_sums[0][1], y_axis_diff[0], marker = 'o', linestyle = "None")
    ax_y_integrate.invert_xaxis()
    ax_y_integrate.set_ylim([-25, 25])
# else:
#     ax_x_integrate.plot(y_axis_sums[1][0], y_axis_sums[1][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_x_integrate.set_xlim(fast_scan_integrator.xlims)

#     ax_y_integrate.invert_yaxis()
#     ax_y_integrate.plot(x_axis_sums[0][0], x_axis_sums[0][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_y_integrate.set_ylim(fast_scan_integrator.ylims)
plt.subplots_adjust(wspace=0, hspace=0)


axd_02 = second_pic.subplot_mosaic(mosaic, width_ratios = [0.25,1], height_ratios = [0.1,1])

ax_x_integrate = axd_02['A panel']
ax_x_integrate.set_xticklabels("")
ax_x_integrate.set_yticklabels("")
ax_x_integrate.tick_params(which='both', top = True, left = False, right = False, bottom = True, direction = "in", width=1)


ax_y_integrate = axd_02['B panel']
ax_y_integrate.set_xticklabels("")
ax_y_integrate.set_yticklabels("")
ax_y_integrate.tick_params(which='both', bottom = False, left = True, right = True, direction = "in", width=1)

ax_detector = axd_02['C panel']
ax_detector.set_yticklabels("")
#ax_detector.set_xlabel(r"$\Delta$y$_{det}$ (pixel)")
ax_detector.tick_params(which='both', top = True, left = True, right = True, bottom = True, direction = "in", width=1)

smin, smax = intensity_scale
if log_scale:
    data = np.log10(imgs[1])
    vmin = max(np.log10(smin), 0)
    vmax = np.log10(smax)
else:
    data = imgs[1]
    vmin, vmax = smin, smax

if rotate_detector == False:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[1].shape[1]),
                                  np.arange(imgs[1].shape[0]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[1].shape[1]),
                                  np.arange(imgs[1].shape[0]))
    # plot data
    ax_detector.pcolormesh(xdata - x_axis_sums[0][0][8], ydata - y_axis_sums[0][0][25], data,  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim([-8, 8])
        ax_detector.set_ylim([-25, 25])
    
else:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[1].shape[0]),
                                  np.arange(imgs[1].shape[1]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[1].shape[0]),
                                  np.arange(imgs[1].shape[1]))
    # plot data
    ax_detector.pcolormesh(xdata - x_axis_sums[0][0][8], ydata - y_axis_sums[0][0][25], np.transpose(data),  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim([-8, 8])
        ax_detector.set_ylim([-25, 25])

if rotate_detector == True:
    ax_x_integrate.plot(x_axis_diff[1], x_axis_sums[1][1], marker = 'o', linestyle = "None")
    ax_x_integrate.set_xlim([-8, 8])

    ax_y_integrate.plot(y_axis_sums[1][1], y_axis_diff[1], marker = 'o', linestyle = "None")
    ax_y_integrate.invert_xaxis()
    ax_y_integrate.set_ylim([-25, 25])
# else:
#     ax_x_integrate.plot(y_axis_sums[1][0], y_axis_sums[1][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_x_integrate.set_xlim(fast_scan_integrator.xlims)

#     ax_y_integrate.invert_yaxis()
#     ax_y_integrate.plot(x_axis_sums[0][0], x_axis_sums[0][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_y_integrate.set_ylim(fast_scan_integrator.ylims)
plt.subplots_adjust(wspace=0, hspace=0)

axd_03 = third_pic.subplot_mosaic(mosaic, width_ratios = [0.25,1], height_ratios = [0.1,1])

ax_x_integrate = axd_03['A panel']
ax_x_integrate.set_xticklabels("")
ax_x_integrate.set_yticklabels("")
ax_x_integrate.tick_params(which='both', top = True, left = False, right = False, bottom = True, direction = "in", width=1)

ax_y_integrate = axd_03['B panel']
ax_y_integrate.set_xticklabels("")
ax_y_integrate.set_yticklabels("")
ax_y_integrate.tick_params(which='both', bottom = False, left = True, right = True, direction = "in", width=1)

ax_detector = axd_03['C panel']
ax_detector.set_yticklabels("")
#ax_detector.set_xlabel(r"$\Delta$y$_{det}$ (pixel)")
ax_detector.tick_params(which='both', top = True, left = True, right = True, bottom = True, direction = "in", width=1)

smin, smax = intensity_scale
if log_scale:
    data = np.log10(imgs[2])
    vmin = max(np.log10(smin), 0)
    vmax = np.log10(smax)
else:
    data = imgs[2]
    vmin, vmax = smin, smax

if rotate_detector == False:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[2].shape[1]),
                                  np.arange(imgs[2].shape[0]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[2].shape[1]),
                                  np.arange(imgs[2].shape[0]))
    # plot data
    ax_detector.pcolormesh(xdata - x_axis_sums[0][0][8], ydata - y_axis_sums[0][0][25], data,  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim([-8, 8])
        ax_detector.set_ylim([-25, 25])
    
else:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[2].shape[0]),
                                  np.arange(imgs[2].shape[1]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[2].shape[0]),
                                  np.arange(imgs[2].shape[1]))
    # plot data
    ax_detector.pcolormesh(xdata - x_axis_sums[0][0][8], ydata - y_axis_sums[0][0][25], np.transpose(data),  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim([-8, 8])
        ax_detector.set_ylim([-25, 25])

if rotate_detector == True:
    ax_x_integrate.plot(x_axis_diff[2], x_axis_sums[2][1], marker = 'o', linestyle = "None")
    ax_x_integrate.set_xlim([-8, 8])

    ax_y_integrate.plot(y_axis_sums[2][1], y_axis_diff[2], marker = 'o', linestyle = "None")
    ax_y_integrate.invert_xaxis()
    ax_y_integrate.set_ylim([-25, 25])
# else:
#     ax_x_integrate.plot(y_axis_sums[1][0], y_axis_sums[1][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_x_integrate.set_xlim(fast_scan_integrator.xlims)

#     ax_y_integrate.invert_yaxis()
#     ax_y_integrate.plot(x_axis_sums[0][0], x_axis_sums[0][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_y_integrate.set_ylim(fast_scan_integrator.ylims)
plt.subplots_adjust(wspace=0, hspace=0)


roi = (1455,126,60,21)
xlims = [roi[0], roi[0]+roi[2]]
ylims = [roi[1], roi[1]+roi[3]]

axd_04 = fourth_pic.subplot_mosaic(mosaic, width_ratios = [0.25,1], height_ratios = [0.1,1])

ax_x_integrate = axd_04['A panel']
ax_x_integrate.set_xticklabels("")
ax_x_integrate.set_yticklabels("")
ax_x_integrate.tick_params(which='both', top = True, left = False, right = False, bottom = True, direction = "in", width=1)

ax_y_integrate = axd_04['B panel']
ax_y_integrate.set_xticklabels("")
ax_y_integrate.set_yticklabels("")
ax_y_integrate.tick_params(which='both', bottom = False, left = True, right = True, direction = "in", width=1)

ax_detector = axd_04['C panel']
ax_detector.set_yticklabels("")
#.set_xlabel(r"$\Delta$y$_{det}$ (pixel)")
ax_detector.tick_params(which='both', top = True, left = True, right = True, bottom = True, direction = "in", width=1)

smin, smax = intensity_scale
if log_scale:
    data = np.log10(imgs[3])
    vmin = max(np.log10(smin), 0)
    vmax = np.log10(smax)
else:
    data = imgs[3]
    vmin, vmax = smin, smax

if rotate_detector == False:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[3].shape[1]),
                                  np.arange(imgs[3].shape[0]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[3].shape[1]),
                                  np.arange(imgs[3].shape[0]))
    # plot data
    ax_detector.pcolormesh(xdata, ydata, data,  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim([-8, 8])
        ax_detector.set_ylim([-25, 25])
    
else:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[3].shape[0]),
                                  np.arange(imgs[3].shape[1]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[3].shape[0]),
                                  np.arange(imgs[3].shape[1]))
    # plot data
    ax_detector.pcolormesh(xdata - x_axis_sums[3][0][9], ydata - y_axis_sums[3][0][34], np.transpose(data),  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim([-8, 8])
        ax_detector.set_ylim([-25, 25])

if rotate_detector == True:
    ax_x_integrate.plot(x_axis_diff[3], x_axis_sums[3][1], marker = 'o', linestyle = "None")
    ax_x_integrate.set_xlim([-8, 8])

    ax_y_integrate.plot(y_axis_sums[3][1], y_axis_diff[3], marker = 'o', linestyle = "None")
    ax_y_integrate.invert_xaxis()
    ax_y_integrate.set_ylim([-25, 25])
# else:
#     ax_x_integrate.plot(y_axis_sums[1][0], y_axis_sums[1][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_x_integrate.set_xlim(fast_scan_integrator.xlims)

#     ax_y_integrate.invert_yaxis()
#     ax_y_integrate.plot(x_axis_sums[0][0], x_axis_sums[0][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_y_integrate.set_ylim(fast_scan_integrator.ylims)
fig.supxlabel(r"$\Delta$y$_{det}$ (pixel)", fontsize = 10)
plt.subplots_adjust(wspace=0, hspace=0)

plt.show()
