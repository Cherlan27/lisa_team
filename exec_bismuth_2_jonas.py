# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 09:20:50 2024

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

scan_numbers = [2025]
experiments = ["bismuth"]

x_axis_sums = []
y_axis_sums = []
imgs = []

for i,j in enumerate([1,50]):
    scan = 2025
    experiment = "bismuth"
    image = j
    roi = (1050, 0, 200, 400)
    
    print(scan)
    
    fast_scan_integrator = Integrator()
    fast_scan_integrator.image_number = image
    fast_scan_integrator.experiment = experiment
    fast_scan_integrator.intensity_scale = (1, 100)
    fast_scan_integrator.scan_number = scan
    fast_scan_integrator.xlims = [roi[0], roi[0]+roi[2]]
    fast_scan_integrator.ylims = [roi[1], roi[1]+roi[3]]

    fast_scan_integrator()
    imgs.append(fast_scan_integrator.output_img())
    x_axis_sums.append(fast_scan_integrator.output_integration_values_x())
    y_axis_sums.append(fast_scan_integrator.output_integration_values_y())
    


rotate_detector = True
intensity_scale = (1, 40)
log_scale = True
xdim = "pixels"
ydim = "pixels"
Q = None,
colormap = "Blues"
xlims = fast_scan_integrator.xlims
ylims = fast_scan_integrator.ylims

    
mosaic = [['B panel','C panel']]
mosaic2 = [['A panel',"."]]
kw = dict(layout='constrained')

plt.rcParams["font.size"] = 8
plt.rcParams['lines.markersize'] = 1

zoll=0.393791
fig = plt.figure(layout="constrained", dpi = 500, figsize=(8*zoll,2))
left, right, sum_up = fig.subfigures(nrows=1, ncols=3)
fig.tight_layout()

axd_01 = left.subplot_mosaic(mosaic, width_ratios = [0.25,1])

# ax_x_integrate = axd_01['A panel']
# ax_x_integrate.set_xticklabels("")
# ax_x_integrate.set_yticklabels("")
# ax_x_integrate.tick_params(which='both', top = True, left = False, right = False, bottom = True, direction = "in", width=1)
# ax_x_integrate.set_xticks([0,100,200])


ax_y_integrate = axd_01['B panel']
ax_y_integrate.set_xticklabels("")
ax_y_integrate.tick_params(direction = "in")
ax_y_integrate.set_ylabel(r"z$_{det}$ (pixel)")
ax_y_integrate.tick_params(which='both', bottom = False, left = True, right = True, direction = "in", width=1)

ax_detector = axd_01['C panel']
ax_detector.set_yticklabels("")
ax_detector.set_xlabel(r"y$_{det}$ (pixel)")
ax_detector.tick_params(which='both', top = True, left = True, right = True, bottom = True, direction = "in", width=1)
ax_detector.set_xticks([0,200,400])

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
    ax_detector.pcolormesh(xdata, ydata, data,  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim(xlims)
        ax_detector.set_ylim(ylims)
    
else:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[0].shape[0]),
                                  np.arange(imgs[0].shape[1]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[0].shape[0]),
                                  np.arange(imgs[0].shape[1]))
    # plot data
    ax_detector.pcolormesh(xdata, ydata, np.transpose(data),  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim(ylims)
        ax_detector.set_ylim(xlims)

if rotate_detector == True:
    # ax_x_integrate.plot(x_axis_sums[0][0], x_axis_sums[0][1], marker = 'o', linestyle = "None")
    # ax_x_integrate.set_xlim(fast_scan_integrator.ylims)

    ax_y_integrate.plot(y_axis_sums[0][1], y_axis_sums[0][0], marker = 'o', linestyle = "None", markersize = 0.5)
    ax_y_integrate.invert_xaxis()
    ax_y_integrate.set_ylim(fast_scan_integrator.xlims)
# else:
#     ax_x_integrate.plot(y_axis_sums[1][0], y_axis_sums[1][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_x_integrate.set_xlim(fast_scan_integrator.xlims)

#     ax_y_integrate.invert_yaxis()
#     ax_y_integrate.plot(x_axis_sums[0][0], x_axis_sums[0][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_y_integrate.set_ylim(fast_scan_integrator.ylims)
plt.subplots_adjust(wspace=0, hspace=0)


axd_02 = right.subplot_mosaic(mosaic, width_ratios = [0.25,1])

# ax_x_integrate = axd_02['A panel']
# ax_x_integrate.set_xticklabels("")
# ax_x_integrate.set_yticklabels("")
# ax_x_integrate.tick_params(which='both', top = True, left = False, right = False, bottom = True, direction = "in", width=1)
# ax_x_integrate.set_xticks([0,100,200])

ax_y_integrate = axd_02['B panel']
ax_y_integrate.set_xticklabels("")
ax_y_integrate.set_yticklabels("")
ax_y_integrate.tick_params(which='both', bottom = False, left = True, right = True, direction = "in", width=1)

ax_detector = axd_02['C panel']
ax_detector.set_yticklabels("")
ax_detector.set_xlabel(r"y$_{det}$ (pixel)")
ax_detector.tick_params(which='both', top = True, left = True, right = True, bottom = True, direction = "in", width=1)
ax_detector.set_xticks([0,200,400])

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
    ax_detector.pcolormesh(xdata, ydata, data,  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim(xlims)
        ax_detector.set_ylim(ylims)
    
else:
    # make xdata
    if xdim == "pixels" or Q is None:
        xdata, _ = np.meshgrid(np.arange(imgs[1].shape[0]),
                                  np.arange(imgs[1].shape[1]))
    if ydim == "pixels" or Q is None:
        _, ydata = np.meshgrid(np.arange(imgs[1].shape[0]),
                                  np.arange(imgs[1].shape[1]))
    # plot data
    ax_detector.pcolormesh(xdata, ydata, np.transpose(data),  cmap=colormap, 
              vmin=vmin, vmax=vmax)
    
    if xlims != None and ylims != None:
        ax_detector.set_xlim(ylims)
        ax_detector.set_ylim(xlims)

if rotate_detector == True:
    # ax_x_integrate.plot(x_axis_sums[1][0], x_axis_sums[1][1], marker = 'o', linestyle = "None", color = "green")
    # ax_x_integrate.set_xlim(fast_scan_integrator.ylims)

    ax_y_integrate.plot(y_axis_sums[1][1], y_axis_sums[1][0], marker = 'o', linestyle = "None", color = "orange", markersize = 0.5)
    ax_y_integrate.invert_xaxis()
    ax_y_integrate.set_ylim(fast_scan_integrator.xlims)
# else:
#     ax_x_integrate.plot(y_axis_sums[1][0], y_axis_sums[1][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_x_integrate.set_xlim(fast_scan_integrator.xlims)

#     ax_y_integrate.invert_yaxis()
#     ax_y_integrate.plot(x_axis_sums[0][0], x_axis_sums[0][1], marker = 'o', markersize = 4, linestyle = "None")
#     ax_y_integrate.set_ylim(fast_scan_integrator.ylims)

axd_02 = sum_up.subplot_mosaic(mosaic2, width_ratios = [0.5,0.5])
ax_y_integrate = axd_02['A panel']
ax_y_integrate.plot(y_axis_sums[0][1], y_axis_sums[0][0], marker = 'o', linestyle = "None", markersize = 0.5)
ax_y_integrate.plot(y_axis_sums[1][1], y_axis_sums[1][0], marker = 'o', linestyle = "None", color = "orange", markersize = 0.5)
ax_y_integrate.tick_params(which='both', top = False, left = True, right = True, bottom = False, direction = "in", width=1)
ax_y_integrate.invert_xaxis()
ax_y_integrate.set_yticklabels("")
ax_y_integrate.set_xticklabels("")
ax_y_integrate.tick_params(direction = "in")
ax_y_integrate.set_ylim(fast_scan_integrator.xlims)

plt.savefig("C:/Users/Lukas/Desktop/jonas_paper/test.png",bbox_inches='tight')
plt.show()


    
# for i,j in enumerate([1,50]):
#     scan = 2025
#     experiment = "bismuth"
#     image = j
#     roi = (1100, 0, 100, 250)
    
#     print(scan)
    
#     fast_scan_detector = D2_Plotter()
#     fast_scan_detector.scan_number = scan
#     fast_scan_detector.image_number = image
#     fast_scan_detector.experiment = experiment
#     fast_scan_detector.rois = {"roi": roi}
#     fast_scan_detector.xlims = [0,1500]
#     fast_scan_detector.ylims = [0,500]
    
#     fast_scan_detector()