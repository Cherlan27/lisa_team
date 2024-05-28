# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 18:29:05 2024

@author: Lukas
"""

import h5py
import numpy
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Rectangle
import hashlib
import os
from PIL import Image
from eigene.p08_detector_read import p08_detector_read
from eigene.fio_reader import read

class D2_Plotter:
    def __init__(self):
        self.flatfield_directory = "./eigene/Module_2017-004_GaAs_MoFluor_Flatfielddata.tif"
        self.pixel_mask_directory = "./eigene/Module_2017-004_GaAs_mask.tif"
        self.use_flatfield = True
        self.use_maskfield = True
        
        self.data_directory = "K:/SYNCHROTRON/Murphy/2023-11_P08_11016542_hoevelmann_petersdorf/raw/"
        self.save_directory = "H:/supex623/supex623/LTP_2023_11/0_5mol_fast_scans/"
        
        self.saving_plot = True
        
        self.scan_number = 2034
        self.experiment = "water"
        self.detector = "lambda"

        self.image_number =  0
        self.log_scale = True
        self.intensity_scale = (1, 1000)
        self.colormap = "viridis"
        self.xdim = "pixels"
        self.ydim = "pixels"
        self.Q = None,
        self.rois = {#"roi": (1470,360,60,20),
                #"roi": (14, 85, 60, 22),
                #"roi": (35, 85, 22, 15),
                "roi": (1450,125,75,20),
                #"roi2_old": (1506, 185, 18, 15),
                #"roi1": (1485,180,60,22),
                #"roi_1": (140,85,60,22),
                #"roi": (1450,354,90,20),
                }
        self.xlims = [1400,1550]
        self.ylims = [100,210]
    
    def __call__(self):
        self.reading_file()
        self.apply_flatfield()
        self.apply_pixel_mask()
        self.plot_detector()
    
    def reading_file(self):
        self.img = p08_detector_read(self.data_directory, self.experiment, self.scan_number, self.detector)()[self.image_number]
        self.fio_filename = "{0}/{1}_{2:05}.fio".format(self.data_directory, self.experiment, self.scan_number)
        self.header, self.column_names, self.data2, self.scan_cmd = read(self.fio_filename)
        
        if self.data2 == {}:
            self.x_axis = "burst"
            self.x_axis_data = numpy.linspace(1, int(self.scan_cmd.split()[1]), int(float(self.scan_cmd.split()[1])/float(self.scan_cmd.split()[2])))
        else:
            self.x_axis = list(self.data2.keys())[0]
            self.x_axis_data = (self.data2[self.x_axis])
        
    def apply_flatfield(self):
        if self.use_flatfield == True:
            flatfield = numpy.array(Image.open(self.flatfield_directory))
            self.img = self.img / flatfield
            
            
    def apply_pixel_mask(self):
        if self.use_maskfield == True:
            self.img_mask = numpy.array(Image.open(self.pixel_mask_directory))  # h5py.File(self.pixel_mask_directory, "r")
            self.mask = numpy.zeros_like(self.img_mask)
            self.mask[(self.img_mask > 1)] = 1
            self.mask = (self.mask == 1)
            self.mask_value = 0
            
            for k,j in enumerate(zip(numpy.where(self.mask != False)[0], numpy.where(self.mask != False)[1])):
                self.mask_intens_weighted = numpy.zeros((2,3,3))
                
                if j[0]-1 > 0 and j[1]-1 > 0:
                    self.mask_intens_weighted[0,0,0] = self.img[j[0]-1, j[1]-1]
                    self.mask_intens_weighted[1,0,0] = numpy.sqrt(2)/2
                if j[0]-1 > 0:
                    self.mask_intens_weighted[0,1,0] = self.img[j[0]-1, j[1]]
                    self.mask_intens_weighted[1,1,0] = 1
                if j[0]-1 > 0 and j[1]+1 < 1554:
                    self.mask_intens_weighted[0,2,0] = self.img[j[0]-1, j[1]+1]
                    self.mask_intens_weighted[1,2,0] = numpy.sqrt(2)/2
                
                if j[1]-1 > 0:
                    self.mask_intens_weighted[0,0,0] = self.img[j[0], j[1]-1]
                    self.mask_intens_weighted[1,0,0] = 1
                if j[1]+1 > 1554:
                    self.mask_intens_weighted[0,2,0] = self.img[j[0], j[1]+1]
                    self.mask_intens_weighted[1,2,0] = 1
                    
                if j[0]+1 > 516 and j[1]-1 > 0:
                    self.mask_intens_weighted[0,0,0] = self.img[j[0]+1, j[1]-1]
                    self.mask_intens_weighted[1,0,0] = numpy.sqrt(2)/2
                if j[0]+1 > 516:
                    self.mask_intens_weighted[0,1,0] = self.img[j[0]+1, j[1]]
                    self.mask_intens_weighted[1,1,0] = 1
                if j[0]+1 > 516 and j[1]-1 < 1554:
                    self.mask_intens_weighted[0,2,0] = self.img[j[0]+1, j[1]+1]
                    self.mask_intens_weighted[1,2,0] = numpy.sqrt(2)/2
                
                temp_masked_intensity = numpy.multiply(self.mask_intens_weighted[0], self.mask_intens_weighted[1]).sum()/self.mask_intens_weighted[1].sum()
                self.img[j[0],j[1]] = temp_masked_intensity
    
    def plot_detector(self):
        rcParams["figure.figsize"] = 8, 5
        rcParams["font.size"] = 20
        rcParams["text.usetex"] = False
        rcParams["text.latex.preamble"] = r"\usepackage{sfmath}"
        fig = plt.figure(dpi = 300)
        fig.patch.set_color("white")
        ax = fig.gca()
        
        AXIS_LABELS = dict(pixels="pixels", qx=r"q$_x$ [\AA$^{-1}$]",
                           qy=r"q$_y$ [\AA$^{-1}$]", qz=r"q$_z$ [\AA$^{-1}$]")

        smin, smax = self.intensity_scale

        if self.log_scale:
            self.data = numpy.log10(self.img)
            vmin = max(numpy.log10(smin), 0)
            vmax = numpy.log10(smax)
        else:
            self.data = self.img
            vmin, vmax = smin, smax

        # make xdata
        if self.xdim == "pixels" or self.Q is None:
            self.xdata, _ = numpy.meshgrid(numpy.arange(self.img.shape[1]),
                                      numpy.arange(self.img.shape[0]))
        elif self.xdim == "qx":
            self.xdata = self.Q[:,:,0] * 1e-10
        elif self.xdim == "qy":
            self.xdata = self.Q[:,:,1] * 1e-10
        elif self.xdim == "qz":
            self.xdata = self.Q[:,:,2] * 1e-10
        # make ydata                            
        if self.ydim == "pixels" or self.Q is None:
            _, self.ydata = numpy.meshgrid(numpy.arange(self.img.shape[1]),
                                      numpy.arange(self.img.shape[0]))
        elif self.ydim == "qx":
            self.ydata = self.Q[:,:,0] * 1e-10
        elif self.ydim == "qy":
            self.ydata = self.Q[:,:,1] * 1e-10
        elif self.ydim == "qz":
            self.ydata = self.Q[:,:,2] * 1e-10
        # plot data
        ax.pcolormesh(self.xdata, self.ydata, self.data,  cmap=self.colormap, 
                  vmin=vmin, vmax=vmax)
        # plot rois
        for self.roi_name in self.rois:
            roi = self.rois[self.roi_name]
            
            color = "#" + hashlib.md5(self.roi_name.encode('utf-8')).hexdigest()[:6]
            # top line
            x = self.xdata[roi[1], roi[0]:(roi[0]+roi[2]+1)]
            y = self.ydata[roi[1], roi[0]:(roi[0]+roi[2]+1)]
            ax.plot(x, y, color=color)
            # bottom line
            x = self.xdata[roi[1]+roi[3], roi[0]:(roi[0]+roi[2]+1)]
            y = self.ydata[roi[1]+roi[3], roi[0]:(roi[0]+roi[2]+1)]
            ax.plot(x, y, color=color)
            # left line
            x = self.xdata[roi[1]:(roi[1]+roi[3]+1), roi[0]]
            y = self.ydata[roi[1]:(roi[1]+roi[3]+1), roi[0]]
            ax.plot(x, y, color=color)
            # right line
            x = self.xdata[roi[1]:(roi[1]+roi[3]+1), roi[0]+roi[2]]
            y = self.ydata[roi[1]:(roi[1]+roi[3]+1), roi[0]+roi[2]]
            ax.plot(x, y, color=color)
            # roi name
            ax.annotate(self.roi_name, (self.xdata[roi[1], roi[0]], self.ydata[roi[1], roi[0]]),
                        color=color, ha="left", va="bottom", size=12)

        # configure plot
        if self.xdim == "pixels" or self.Q is None:
            ax.xaxis.set_minor_locator(MultipleLocator(20))
        if self.ydim == "pixels" or self.Q is None:

            ax.yaxis.set_minor_locator(MultipleLocator(20))
        ax.set_aspect("auto")
        ax.set_xlabel(AXIS_LABELS[self.xdim] if self.Q is not None else "x-pixel", fontsize = 20)
        ax.set_ylabel(AXIS_LABELS[self.ydim] if self.Q is not None else "y-pixel", fontsize = 20)
        
        if self.xlims != None and self.ylims != None:
            ax.set_xlim(self.xlims)
            ax.set_ylim(self.ylims)
        plt.subplots_adjust(bottom=0.12)
        
        plt.show()
        self.save_plot(fig, "detector_image")
    
    def save_plot(self, fig, plot_name):
        if self.saving_plot == True:
            if not os.path.exists(self.save_directory):
                os.mkdir(self.save_directory)
            out_filename_fits = self.save_directory + "scan_" + str(self.experiment) + "_" + str(self.x_axis) + "_" + str(self.scan_number) + "_" + str(plot_name) + ".png"
            fig.savefig(out_filename_fits)