# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 14:03:07 2024

@author: Lukas
"""

from eigene.fio_reader import read
import numpy
import h5py
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator
import matplotlib.lines as mlines
import pandas
import os
from eigene.abs_overlap_fit_poly import Absorber
import sys
from PIL import Image
from eigene.p08_detector_read import p08_detector_read
import matplotlib.gridspec as gridspec
from scipy import optimize
from typing import List

class Amptek_plotter():
    def __init__(self):
        self.data_directory = "K:/SYNCHROTRON/Murphy/2023-11_P08_11016542_hoevelmann_petersdorf/raw"
        self.save_directory = "H:/supex623/supex623/LTP_2023_11/Amptek/"
        self.scan_number_off = [2845]
        self.scan_number_on = [2866]
        self.experiment = "nai3mol"
        self.first_image_number = 1
        self.last_image_number = "all" # "all"
        self.peak_number = 2
        self.height_plot = 5000
        self.normalization = "gaussian" #heigth, gaussian, auto
        self.checking_plot = True
        self.save_plot = True
        self.xlim = [5,11]
        
    def __call__(self):
        self.calibrate_energy()
        self.get_flourescence_off()
        self.get_flourescence_on()
        self.normalisation()
        self.fitting_off()
        self.fitting_on()

    
    def calibrate_energy(self):
        self.energy_d = (14.96-11.92)/(1683-1343)
        self.energy_spectrum = numpy.linspace(14.96+self.energy_d*(0-1683), 14.96+self.energy_d*(2047-1683), 2048)
        
    def get_flourescence_off(self):
        self.spectrum_off = numpy.ones(shape = (len(self.scan_number_off), len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_off[0]))), 2048))
        self.q_off = numpy.ones(shape = (len(self.scan_number_off), len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_off[0])))))
        
        self.params_max_off = numpy.ones(shape = (len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_on[0]))), 3))
        
        for i,scan in enumerate(self.scan_number_off):
            
            if self.last_image_number == "all":
                self.last_image_number = len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, scan)))
            
            self.fio_filename = "{0}/{1}_{2:05}.fio".format(self.data_directory, self.experiment, scan)
            self.header_gen, self.column_names_gen, self.data_gen, self.scan_cmd_gen = read(self.fio_filename)
            
            self.q_off[i] = self.data_gen["qz"]
            
            for j,scan_point in enumerate(range(self.first_image_number, self.last_image_number + 1)):
                self.amptek_filename_off = "{0}/{1}_{2:05}/{1}_{2:05}_mca_s{3:01}.fio".format(self.data_directory, self.experiment, scan, scan_point)
                self.header_off, self.column_names_off, self.data_off, self.scan_cmd_off = read(self.amptek_filename_off)
                self.sample_time_off = self.header_off["Sample_time"]
                self.spectrum_off[i,j] = self.data_off["amptek_spectrum"]/self.sample_time_off/self.data_gen["ion2"][j]
                
                if self.normalisation == "gaussian":
                    initial = [500, self.energy_spectrum[1910],0.05]
                    self.params_max_off[j], _ = optimize.curve_fit(self.gaussian, self.energy_spectrum, self.spectrum_off[i,j], p0=initial, maxfev = 200000)
                else:
                    self.params_max_off[j,0] = self.spectrum_off[i,j,1890:2047].max()
                
    def get_flourescence_on(self):
        self.spectrum_on = numpy.ones(shape = (len(self.scan_number_on), len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_on[0]))), 2048))
        self.q_on = numpy.ones(shape = (len(self.scan_number_on), len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_on[0])))))
        
        self.params_max_on = numpy.ones(shape = (len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_on[0]))), 3))
        
        for i,scan in enumerate(self.scan_number_on):
            
            if self.last_image_number == "all":
                self.last_image_number = len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, scan)))
            
            self.fio_filename = "{0}/{1}_{2:05}.fio".format(self.data_directory, self.experiment, scan)
            self.header_gen, self.column_names_gen, self.data_gen, self.scan_cmd_gen = read(self.fio_filename)
            
            self.q_on[i] = self.data_gen["qz"]
            
            for j,scan_point in enumerate(range(self.first_image_number, self.last_image_number + 1)):
                self.amptek_filename_on = "{0}/{1}_{2:05}/{1}_{2:05}_mca_s{3:01}.fio".format(self.data_directory, self.experiment, scan, scan_point)
                self.header_on, self.column_names_on, self.data_on, self.scan_cmd_on = read(self.amptek_filename_on)
                self.sample_time_on = self.header_on["Sample_time"]
                self.spectrum_on[i,j] = self.data_on["amptek_spectrum"]/self.sample_time_on/self.data_gen["ion2"][j]
                
                if self.normalisation == "gaussian":
                    initial = [500, self.energy_spectrum[1910],0.05]
                    self.params_max_on[j], _ = optimize.curve_fit(self.gaussian, self.energy_spectrum, self.spectrum_on[i,j], p0=initial, maxfev = 200000)
                else:
                    self.params_max_on[j,0] = self.spectrum_on[i,j,1890:2047].max()
            
    def gaussian(self, xdata: List[float], amp: float, pos: float, sigma: float, *args, **kwargs) -> List[float]:
        return amp * numpy.exp(-(xdata-pos)**2/(2*sigma**2))
    
    def gaussian_2peaks(self, xdata, bg0, bg1, amp1, pos1, sigma1, amp2, pos2, sigma2):
        return bg0 + bg1 * xdata + amp1 * numpy.exp(-(xdata-pos1)**2/(2*sigma1**2)) + amp2 * numpy.exp(-(xdata-pos2)**2/(2*sigma2**2))
    
    def fitting_off(self):
        self.total_spectrum_off = self.spectrum_off.sum(axis=0)
        self.params_off = numpy.ones(shape = (len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_off[0]))), 8))
        self.params_off_err = numpy.ones(shape = (len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_off[0]))), 8))
        
        
        for j,scan_point in enumerate(range(self.first_image_number, self.last_image_number + 1)):
            if self.peak_number == 2:
                initial = [0, 0,2, self.energy_spectrum[421],0.05, 2, self.energy_spectrum[450], 0.05]
                self.params_off[j], params_off_cov = optimize.curve_fit(self.gaussian_2peaks, self.energy_spectrum, self.total_spectrum_off[j], p0=initial, maxfev=200000)
                self.params_off_err[j] = numpy.sqrt(numpy.diag(abs(params_off_cov)))
                
            if self.checking_plot == True:
                plt.plot(self.energy_spectrum,self.total_spectrum_off[j])
                plt.plot(self.energy_spectrum,self.gaussian_2peaks(self.energy_spectrum, *self.params_off[j]))
                #plt.plot(self.energy_spectrum,self.gaussian(self.energy_spectrum, *self.params_max_off[j]))
                plt.show()
                plt.close()
        
    def fitting_on(self):
        self.total_spectrum_on = self.spectrum_on.sum(axis=0)
        self.params_on = numpy.ones(shape = (len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_on[0]))), 8))
        self.params_on_err = numpy.ones(shape = (len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, self.scan_number_on[0]))), 8))
        
        
        for j,scan_point in enumerate(range(self.first_image_number, self.last_image_number + 1)):
            if self.peak_number == 2:
                initial = [0, 0,2, self.energy_spectrum[421],0.05, 2000, self.energy_spectrum[450], 0.05]
                self.params_on[j], params_on_cov = optimize.curve_fit(self.gaussian_2peaks, self.energy_spectrum, self.total_spectrum_on[j], p0=initial, maxfev=200000)
                self.params_on_err[j] = numpy.sqrt(numpy.diag(abs(params_on_cov)))
                
                
        
            if self.checking_plot == True:
                plt.plot(self.energy_spectrum,self.total_spectrum_on[j])
                plt.plot(self.energy_spectrum,self.gaussian_2peaks(self.energy_spectrum, *self.params_on[j]))
                plt.show()
                plt.close()
    
    def normalisation(self):
        self.spectrum_on = self.spectrum_on/numpy.expand_dims(self.params_max_on[:,0], axis=(1))
        self.spectrum_off = self.spectrum_off/numpy.expand_dims(self.params_max_off[:,0], axis=(1))
    
    def plot_signal(self):
        fig = plt.figure(dpi = 300, figsize = (8, 20))
        ax = fig.gca()
        for i,scan in enumerate(self.scan_number_off):
            for i in range(len(os.listdir("{0}/{1}_{2:05}/".format(self.data_directory, self.experiment, scan)))):
                plt.plot(self.energy_spectrum, self.spectrum_off[0][i] + self.height_plot*i, color = "blue")
                plt.plot(self.energy_spectrum, self.spectrum_on[0][i] + self.height_plot*i, color = "red")
                text_q = str(self.q_off[0,i])[0:6] + "(1/A)"
                plt.text(2.8, self.height_plot*i + 500, text_q)
        ax.set_ylabel("Intensity (arb. u.)")
        ax.set_xlabel("Energy (keV)")
        ax.set_xlim([2.5,5])
        plt.show()
        
    def plot_fitting(self):
        fig = plt.figure(dpi = 300)
        ax = fig.gca()
        plt.plot(self.q_off[0], self.params_off[:,2], color = "blue", linestyle = "--", marker = "o")
        plt.plot(self.q_on[0], self.params_on[:,2], color = "red", linestyle = "--", marker = "o")
        ax.set_ylabel("Intensity (arb. u.)")
        ax.set_xlabel("qz (1/A)")
        plt.show()
    
    def get_signal(self):
        return self.energy_spectrum, self.spectrum_off[0], self.spectrum_on[0]
    
    def get_fitting(self):
        return self.q_off[0], self.q_on[0], self.params_off[:,2], self.params_on[:,2]
        
        
