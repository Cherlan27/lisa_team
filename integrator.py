# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:22:40 2024

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

class Integrator:
    def __init__(self):
        self.flatfield_directory = "./eigene/Module_2017-004_GaAs_MoFluor_Flatfielddata.tif"
        self.pixel_mask_directory = "./eigene/Module_2017-004_GaAs_mask.tif"
        self.use_flatfield = True
        self.use_maskfield = True
        
        self.data_directory = "K:/SYNCHROTRON/Murphy/2023-11_P08_11016542_hoevelmann_petersdorf/raw/"
        self.save_directory = "H:/supex623/supex623/LTP_2023_11/0_5mol_fast_scans/"
        
        self.saving_plot = True
        
        self.scan_number = 3227
        self.experiment = "nai5molxrr"
        self.detector = "lambda"

        self.image_number =  0
        self.log_scale = True
        self.intensity_scale = (1, 20000)
        self.colormap = "Blues" # "viridis"
        self.xdim = "pixels"
        self.ydim = "pixels"
        self.Q = None,
        self.xlims = [1450,1530]
        self.ylims =  [125,145]
        
        self.rotate_detector = True
        
    def __call__(self):
        self.reading_file()
        self.apply_flatfield()
        self.apply_pixel_mask()
        self.plot_detector()
        self.plot_detector_data()
        self.integrating_signal()
        self.plot_integration()

    def reading_file(self):
        self.img = p08_detector_read(self.data_directory, self.experiment, self.scan_number, self.detector)()[self.image_number]
        try: 
            self.fio_filename = "{0}/{1}_{2:05}.fio".format(self.data_directory, self.experiment, self.scan_number)
            self.header, self.column_names, self.data2, self.scan_cmd = read(self.fio_filename)
        except FileNotFoundError:
            self.fio_filename = "{0}/{1}_{2:05}.log".format(self.data_directory, self.experiment, self.scan_number)
            self.header, self.column_names, self.data2, self.scan_cmd = read(self.fio_filename)
        finally:
            print("File is not existing!")
        
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

    def integrating_signal(self):
            self.x_axes_sum = numpy.ones((self.ylims[1]-self.ylims[0]+1))
            self.x_axes_pixels = range(self.ylims[0], self.ylims[1]+1)
            
            for i,x in  enumerate(range(self.ylims[0], self.ylims[1]+1)) :
                self.x_axes_sum[i] = self.img[x,self.xlims[0]:self.xlims[1]].sum()
                
            self.y_axes_sum = numpy.ones((self.xlims[1]-self.xlims[0]+1))
            self.y_axes_pixels = range(self.xlims[0], self.xlims[1]+1)
            
            for i,y in  enumerate(range(self.xlims[0], self.xlims[1]+1)) :
                self.y_axes_sum[i] = self.img[self.ylims[0]:self.ylims[1],y].sum()

    def plot_detector(self):
        mosaic = [['.', 'A panel',],['B panel','C panel']]
        
        kw = dict(layout='constrained')
        self.fig, self.axs = plt.subplot_mosaic(mosaic, width_ratios = [0.25,1], height_ratios = [0.1,1], figsize = (3,10), **kw)
        self.fig.tight_layout()
        
        self.ax_x_integrate = self.axs['A panel']
        self.ax_x_integrate.set_xticklabels("")
        self.ax_x_integrate.set_yticklabels("")
        self.ax_x_integrate.tick_params(direction = "in")
        
        self.ax_y_integrate = self.axs['B panel']
        self.ax_y_integrate.set_xticklabels("")
        self.ax_y_integrate.tick_params(direction = "in")
        
        self.ax_detector = self.axs['C panel']
        self.ax_detector.set_yticklabels("")
        self.ax_detector.tick_params(direction = "in")
        
        #self.fig.text(0.5, -0.05, 'y', ha='center')
        #self.fig.text(-0.05, 0.5, 'x', va='center', rotation='vertical')
        plt.subplots_adjust(wspace=0, hspace=0)
        
    def plot_detector_data(self):
        smin, smax = self.intensity_scale

        if self.log_scale:
            self.data = numpy.log10(self.img)
            vmin = max(numpy.log10(smin), 0)
            vmax = numpy.log10(smax)
        else:
            self.data = self.img
            vmin, vmax = smin, smax

        if self.rotate_detector == False:
            # make xdata
            if self.xdim == "pixels" or self.Q is None:
                self.xdata, _ = numpy.meshgrid(numpy.arange(self.img.shape[1]),
                                          numpy.arange(self.img.shape[0]))
            if self.ydim == "pixels" or self.Q is None:
                _, self.ydata = numpy.meshgrid(numpy.arange(self.img.shape[1]),
                                          numpy.arange(self.img.shape[0]))
            # plot data
            self.ax_detector.pcolormesh(self.xdata, self.ydata, self.data,  cmap=self.colormap, 
                      vmin=vmin, vmax=vmax)
            
            if self.xlims != None and self.ylims != None:
                self.ax_detector.set_xlim(self.xlims)
                self.ax_detector.set_ylim(self.ylims)
        else:
            # make xdata
            if self.xdim == "pixels" or self.Q is None:
                self.xdata, _ = numpy.meshgrid(numpy.arange(self.img.shape[0]),
                                          numpy.arange(self.img.shape[1]))
            if self.ydim == "pixels" or self.Q is None:
                _, self.ydata = numpy.meshgrid(numpy.arange(self.img.shape[0]),
                                          numpy.arange(self.img.shape[1]))
            # plot data
            self.ax_detector.pcolormesh(self.xdata, self.ydata, numpy.transpose(self.data),  cmap=self.colormap, 
                      vmin=vmin, vmax=vmax)
            
            if self.xlims != None and self.ylims != None:
                self.ax_detector.set_xlim(self.ylims)
                self.ax_detector.set_ylim(self.xlims)
    
    def plot_integration(self): 
        if self.rotate_detector == True:
            self.ax_x_integrate.plot(self.x_axes_pixels, self.x_axes_sum, marker = 'o', markersize = 4, linestyle = "None")
            self.ax_x_integrate.set_xlim(self.ylims)
        
            self.ax_y_integrate.plot(self.y_axes_sum, self.y_axes_pixels, marker = 'o', markersize = 4, linestyle = "None")
            self.ax_y_integrate.invert_xaxis()
            self.ax_y_integrate.set_ylim(self.xlims)
        else:
            self.ax_x_integrate.plot(self.y_axes_pixels, self.y_axes_sum, marker = 'o', markersize = 4, linestyle = "None")
            self.ax_y_integrate.invert_yaxis()
            self.ax_x_integrate.set_xlim(self.xlims)
        
            self.ax_y_integrate.plot(self.x_axes_sum, self.x_axes_pixels, marker = 'o', markersize = 4, linestyle = "None")
            self.ax_y_integrate.set_ylim(self.ylims)
        
        plt.show()
        self.save_plot(self.fig, "detector_integration")
        
    def output_integration_values_y(self):
        return self.y_axes_pixels, self.y_axes_sum
    
    def output_integration_values_x(self):
        return self.x_axes_pixels, self.x_axes_sum
    
    def output_img(self):
        return self.img
        
    def save_plot(self, fig, plot_name):
        if self.saving_plot == True:
            if not os.path.exists(self.save_directory):
                os.mkdir(self.save_directory)
            out_filename_fits = self.save_directory + "scan_" + str(self.experiment) + "_" + str(self.x_axis) + "_" + str(self.scan_number) + "_" + str(plot_name) + ".png"
            fig.savefig(out_filename_fits)
    