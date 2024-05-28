# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 10:23:14 2024

@author: Lukas
"""

import numpy
from scipy.special import wofz
from eigene.fio_reader import read
import h5py
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from PIL import Image
from eigene.p08_detector_read import p08_detector_read
import os

# Scan fitter

class Timing_scan:
    def __init__(self):
        self.qz_range_peaksplit = 0.4
        self.peaksplit_detection = False
        self.peak_splitting = False
        self.peak_split_detected = False
        
        self._rel_peak_height = 0.25
        self._peak_px_distance = 2
        
        self.fitfunction = "gaussian"
        self.normalization = "apd"
        self.laser_on = None
        self.laser_off = None
        self.peak_detector = "norm"
        
        self.use_flatfield = True
        self.use_maskfield = True
        self.flatfield_directory = "./eigene/Module_2017-004_GaAs_MoFluor_Flatfielddata.tif"
        self.pixel_mask_directory = "./eigene/Module_2017-004_GaAs_mask.tif"
        
        self.data_directory = "K:/SYNCHROTRON/Murphy/2023-11_P08_11016542_hoevelmann_petersdorf/raw"
        self.save_directory = "H:/supex623/supex623/LTP_2023_11/3mol_fast_scans_2/"
        
        self.saving_plot = True
        
        self.experiment = "nai3moltiming5"
        self.detector = "lambda"
        self.roi = (1455,126,60,21)
        self.roi_offset = 20
        self.scan_number = 3976
        
        self.ration = 0.45226130653266333
        
    def __call__(self):
        self.reading_file()
        self.get_x_axis()
        self.define_fit_variables()
        self.get_intensity()
        self.fit_axis_x()
        self.prepare_axis()


    def gaussian(self, xdata, bg0, bg1, amp, pos, sigma):
        return bg0 + bg1 * xdata + amp * numpy.exp(-(xdata-pos)**2/(2*sigma**2))
    
    def gaussian_double(self, xdata, bg0, bg1, amp1, pos1, sigma1, amp2, pos2, sigma2):
                    return bg0 + bg1 * xdata + (amp1 * numpy.exp(-(xdata-pos1)**2/(2*sigma1**2)) + amp2 * numpy.exp(-(xdata-pos2)**2/(2*sigma2**2)))
        
    def lorentzian(self, xdata, bg0, bg1, amp, pos, sigma):
        return bg0 + bg1 * xdata +  amp * 1/(1+((xdata-pos)/sigma)**2)
    
    def lorentzian_double(self, xdata, bg0, bg1, amp1, pos1, gamma1, amp2, pos2, gamma2):
                    return bg0 + bg1 * xdata + (amp1 * gamma1/(numpy.pi*((xdata-pos1)**2+gamma1**2)) + amp2 * gamma2/(numpy.pi*((xdata-pos2)**2+gamma2**2)))
    
    def voigt(self, xdata, bg0, bg1, amp, pos, sigma, gamma):
        w = wofz((xdata - pos +1j*gamma)/(sigma*numpy.sqrt(2)))
        return bg0 + bg1 * xdata + amp * w.real/(sigma*numpy.sqrt(2*numpy.pi))*1/max(self.voigt2(xdata, bg0, bg1,pos,sigma,gamma))
        
    def voigt2(self, xdata, bg0, bg1, pos, sigma, gamma):
        w = wofz((xdata - pos +1j*gamma)/(sigma*numpy.sqrt(2)))
        return  w.real/(sigma*numpy.sqrt(2*numpy.pi))    
    
    def voigt_ratio(self, xdata, bg0, bg1, amp, pos, sigma):
        w = wofz((xdata - pos +1j*sigma*self.ratio)/(sigma*numpy.sqrt(2)))
        return bg0 + bg1 * xdata + amp * w.real/(sigma*numpy.sqrt(2*numpy.pi))*1/max(self.voigt2(xdata, bg0, bg1,pos,sigma,sigma*self.ratio))
    
    def voigt_ratio_double(self, xdata, bg0, bg1, amp1, pos1, sigma1, amp2, pos2, sigma2):
        w1 = wofz((xdata - pos1 +1j*sigma1*self.ratio)/(sigma1*numpy.sqrt(2)))
        voi1 = bg0 + bg1 * xdata + amp1 * w1.real/(sigma1*numpy.sqrt(2*numpy.pi))*1/max(self.voigt2(xdata, bg0, bg1,pos1,sigma1,sigma1*self.ratio))  
        w2 = wofz((xdata - pos2 +1j*sigma2*self.ratio)/(sigma2*numpy.sqrt(2)))
        voi2 = bg0 + bg1 * xdata + amp2 * w2.real/(sigma1*numpy.sqrt(2*numpy.pi))*1/max(self.voigt2(xdata, bg0, bg1,pos2,sigma2,sigma2*self.ratio))  
        return voi1 + voi2
    
    def reading_file(self):
        fio_filename = "{0}/{1}_{2:05}.fio".format(self.data_directory, self.experiment, self.scan_number)
        self.header, self.column_names, self.data, self.scan_cmd = read(fio_filename)
        filename = "{0}/{1}_{2:05}/{3}/{4}_{5:05}_{6:05}.nxs".format(self.data_directory, self.experiment, self.scan_number, self.detector, self.experiment, self.scan_number, 0)
        self.det_file = h5py.File(filename)
        
        if self.data == {}:
            self.monitor = p08_detector_read(self.data_directory, self.experiment, self.scan_number, "pilc")()[self.normalization]
        else:
            self.monitor = self.data[self.normalization]
        
    def get_x_axis(self):
        if self.data == {}:
            self.x_axis = "burst"
            self.x_axis_data = numpy.linspace(float(self.scan_cmd.split()[2]), int(self.scan_cmd.split()[1])*float(self.scan_cmd.split()[2]), int(float(self.scan_cmd.split()[1])))
        else:
            self.x_axis = list(self.data.keys())[0]
            self.x_axis_data = (self.data[self.x_axis])
            
    def apply_flatfield(self):
        if self.use_flatfield == True:
            flatfield = numpy.array(Image.open(self.flatfield_directory))
            self.data = self.data / flatfield
            
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
                    self.mask_intens_weighted[0,0,0] = self.data[j[0]-1, j[1]-1]
                    self.mask_intens_weighted[1,0,0] = numpy.sqrt(2)/2
                if j[0]-1 > 0:
                    self.mask_intens_weighted[0,1,0] = self.data[j[0]-1, j[1]]
                    self.mask_intens_weighted[1,1,0] = 1
                if j[0]-1 > 0 and j[1]+1 < 1554:
                    self.mask_intens_weighted[0,2,0] = self.data[j[0]-1, j[1]+1]
                    self.mask_intens_weighted[1,2,0] = numpy.sqrt(2)/2
                
                if j[1]-1 > 0:
                    self.mask_intens_weighted[0,0,0] = self.data[j[0], j[1]-1]
                    self.mask_intens_weighted[1,0,0] = 1
                if j[1]+1 > 1554:
                    self.mask_intens_weighted[0,2,0] = self.data[j[0], j[1]+1]
                    self.mask_intens_weighted[1,2,0] = 1
                    
                if j[0]+1 > 516 and j[1]-1 > 0:
                    self.mask_intens_weighted[0,0,0] = self.data[j[0]+1, j[1]-1]
                    self.mask_intens_weighted[1,0,0] = numpy.sqrt(2)/2
                if j[0]+1 > 516:
                    self.mask_intens_weighted[0,1,0] = self.data[j[0]+1, j[1]]
                    self.mask_intens_weighted[1,1,0] = 1
                if j[0]+1 > 516 and j[1]-1 < 1554:
                    self.mask_intens_weighted[0,2,0] = self.data[j[0]+1, j[1]+1]
                    self.mask_intens_weighted[1,2,0] = numpy.sqrt(2)/2
                
                temp_masked_intensity = numpy.multiply(self.mask_intens_weighted[0], self.mask_intens_weighted[1]).sum()/self.mask_intens_weighted[1].sum()
                self.data[j[0],j[1]] = temp_masked_intensity
   
    def define_fit_variables(self):
        if self.fitfunction == "gaussian" or self.fitfunction == "lorentzian":
            num_fit_params = 5
        elif self.fitfunction == "voigt_ratio" or self.fitfunction == "voigt":
            num_fit_params = 6
            
        self.fit_x_1 = numpy.ones((num_fit_params, len(self.det_file["/entry/instrument/detector/data"])))
        self.fit_x_2 = numpy.ones((num_fit_params, len(self.det_file["/entry/instrument/detector/data"])))
        self.fit_err_x_1 = numpy.ones((num_fit_params, len(self.det_file["/entry/instrument/detector/data"])))
        self.fit_err_x_2 = numpy.ones((num_fit_params, len(self.det_file["/entry/instrument/detector/data"])))
        
        self.fit_y_1 = numpy.ones((num_fit_params, len(self.det_file["/entry/instrument/detector/data"])))
        self.fit_y_2 = numpy.ones((num_fit_params, len(self.det_file["/entry/instrument/detector/data"])))
        self.fit_err_y_1 = numpy.ones((num_fit_params, len(self.det_file["/entry/instrument/detector/data"])))
        self.fit_err_y_2 = numpy.ones((num_fit_params, len(self.det_file["/entry/instrument/detector/data"])))
        
        self.intensity = numpy.ones((2, len(self.det_file["/entry/instrument/detector/data"])))

    def get_intensity(self):
        for n in range(len(self.det_file["/entry/instrument/detector/data"])):
            
            self.data = numpy.array(self.det_file["/entry/instrument/detector/data"][n]) 
            self.apply_flatfield()
            self.apply_pixel_mask()
    
            self.roi_data = self.data[self.roi[1]:(self.roi[1]+self.roi[3]), self.roi[0]:(self.roi[0]+self.roi[2])]
            self.x_pixels = numpy.arange(self.roi[2]) + self.roi[0]
            self.y_pixels = numpy.arange(self.roi[3]) + self.roi[1]
            self.x_data = self.roi_data.sum(axis=0)
            self.y_data = self.roi_data.sum(axis=1)
            
            self.int_bg0 = self.data[(self.roi[1]+self.roi[3]+self.roi_offset):(self.roi[1]+2*self.roi[3]+self.roi_offset),
                                              (self.roi[0]):(self.roi[0]+self.roi[2])].sum()
            self.int_bg1 = self.data[(self.roi[1]-self.roi[3]-self.roi_offset):(self.roi[1]-self.roi_offset),
                                              (self.roi[0]):(self.roi[0]+self.roi[2])].sum()
            
            self.int_specular = self.roi_data.sum()
            self.intensity[0, n] = (self.int_specular - (self.int_bg0 + self.int_bg1) / 2.0) / self.monitor[n]
            if self.normalization == "apd" or self.normalization == "apd2":
                self.intensity[1, n] = ((numpy.sqrt(self.int_specular) + (numpy.sqrt(self.int_bg0) + numpy.sqrt(self.int_bg1)) / 2.0) /self.monitor[n]
                                             + abs (numpy.sqrt(self.int_specular - (self.int_bg0 + self.int_bg1) / 2.0) / self.monitor[n]))
            elif self.normalization == "ion" or self.normalization == "ion2":
                self.intensity[1, n] = ((numpy.sqrt(self.int_specular) + (numpy.sqrt(self.int_bg0) + numpy.sqrt(self.int_bg1)) / 2.0) /self.monitor[n]
                                             + abs (0.1 * (self.int_specular - (self.int_bg0 + self.int_bg1) / 2.0) / self.monitor[n]))
            else:
                self.intensity[1, n] = (numpy.sqrt(self.int_specular) + (numpy.sqrt(self.int_bg0) + numpy.sqrt(self.int_bg1)) / 2.0) /self.monitor[n]
     
    def detect_peak_split(self):
        pass
    
    def fit_axis_x(self):
        for n in range(len(self.det_file["/entry/instrument/detector/data"])):
            self.data = numpy.array(self.det_file["/entry/instrument/detector/data"][n]) 
            self.apply_flatfield()
            self.roi_data = self.data[self.roi[1]:(self.roi[1]+self.roi[3]), self.roi[0]:(self.roi[0]+self.roi[2])]
            self.x_pixels = numpy.arange(self.roi[2]) + self.roi[0]
            self.x_data = self.roi_data.sum(axis=0)
            
            underline = 0
            slope = 0
            height = abs(self.x_data.max())
            pos = numpy.where(self.x_data == self.x_data.max())
            pos = pos[0][0] + self.x_pixels[0]
            sigma = 2
            if self.fitfunction == "gaussian":
                x_pinit = [underline, slope, height, pos, sigma]
                x_popt_l, x_pcov_l = curve_fit(self.gaussian, self.x_pixels, self.x_data, x_pinit, maxfev=200000)
            elif self.fitfunction == "lorentzian":
                x_pinit = [underline, slope, height, pos, sigma]
                x_popt_l, x_pcov_l = curve_fit(self.lorentzian, self.x_pixels, self.x_data, x_pinit, maxfev=200000)
            elif self.fitfunction == "voigt":
                x_pinit = [underline, slope, height, pos, sigma, sigma]
                x_popt_l, x_pcov_l = curve_fit(self.voigt, self.x_pixels, self.x_data, x_pinit, maxfev=500000)
                self.fit_x_1[5, n] = x_popt_l[5]
                self.fit_err_x_1[5, n] = numpy.sqrt(numpy.diag(abs(x_pcov_l)))[5]
            elif self.fitfunction == "voigt_ratio":
                boundaries = ([min(self.x_data)*0.2,-0.000001,height*0.7, pos-2, 0],[min(self.x_data)*1.2,0.0000005,height*1.2, pos+3, 20])
                x_pinit = [underline, slope, height, pos, sigma]
                x_popt_l, x_pcov_l = curve_fit(self.voigt_ratio, self.x_pixels, self.x_data, p0=x_pinit,bounds=boundaries, maxfev=500000)
                self.fit_x_1[5, n] = x_popt_l[4]
                self.fit_err_x_1[5, n] = numpy.sqrt(numpy.diag(abs(x_pcov_l)))[4]
            for k in range(5):
                self.fit_x_1[k, n] = x_popt_l[k]
                self.fit_err_x_1[k, n] = numpy.sqrt(numpy.diag(abs(x_pcov_l)))[k]
            
    def prepare_axis(self):
        if self.x_axis == "pp_delay":
            self.x_axis_data = self.x_axis_data * 1e09
            self.x_axis_label = "t in ns"
        elif self.x_axis == "delayg":
            self.x_axis_data = self.x_axis_data * 1e09
            self.x_axis_label = "t in ns"
        elif self.x_axis == "lmy":
            self.x_axis_label = "lmy (arb. u.)"
        elif self.x_axis == "dumm":
            self.x_axis_label = "t in s"
        elif self.x_axis == "mdds":
           self.x_axis_data = (self.x_axis_data - self.x_axis_data[int(len(self.x_axis_data)/2)]) / 360
           self.x_axis_label = "t in ns"
        elif self.x_axis == "bcl":
           self.x_axis_data = (self.x_axis_data - self.x_axis_data[int(len(self.x_axis_data)/2)])
           self.x_axis_label = "t in ns"
        elif self.x_axis == "burst":
            self.x_axis_label = "t in s"
        else :
            self.x_axis_label = "(arb. u.)"
    
    def output_intensity(self):
        return self.x_axis_data, self.intensity[0, :], self.intensity[1, :]
    
    def plot_intensity(self):
        fig = plt.figure(figsize=(18,11), dpi = 300)
        fig.patch.set_color('white')
        plt.title(self.x_axis + " scan")
        ax = fig.gca()
        ax.errorbar(self.x_axis_data, self.intensity[0, :], yerr = self.intensity[1, :],
                      color='#ABABAB', ls='--', marker='o', mec='#424242',
                      mfc='white', mew=1.2)
        ax.set_xlabel(self.x_axis_label)
        ax.set_ylabel('Intensity in arb. u.')
        if self.laser_on != None and self.laser_off != None:
            ax.axvline(x=self.laser_on, color='r')
            if len(self.scan_cmd) > 4:
                ax.axvline(x=self.laser_off, color='b')
        self.save_plot(fig, "intensity")
    
    def plot_fitted_params_x(self):
        fig = plt.figure(figsize=(18,11), dpi = 300)
        ay_pos = fig.add_subplot(221)
        plt.subplots_adjust(wspace = .15)
        plt.subplots_adjust(hspace = .12)
        ay_pos.errorbar(self.x_axis_data, self.fit_x_1[3], yerr=self.fit_err_x_1[3],
                      color='#ABABAB', ls='none', marker='o', mec='#424242',
                      mfc='white', mew=1.2)
        ay_pos.set_ylabel('x position')
        ay_pos.set_title("Peak Position")
        if self.laser_on != None and self.laser_off != None:
            ay_pos.axvline(x=self.laser_on, color='r')
            if len(self.scan_cmd) > 4:
                ay_pos.axvline(x=self.laser_off, color='b')
        
        
        ay_amp = fig.add_subplot(222)
        ay_amp.errorbar(self.x_axis_data, self.fit_x_1[2], yerr=self.fit_err_x_1[2],
                      color='#ABABAB', ls='none', marker='o', mec='#424242',
                      mfc='white', mew=1.2)
        ay_amp.set_ylabel('x amp')
        ay_amp.set_title("Peak Amplitude")
        if self.laser_on != None and self.laser_off != None:
            ay_amp.axvline(x=self.laser_on, color='r')
            if len(self.scan_cmd) > 4:
                ay_amp.axvline(x=self.laser_off, color='b')
        
        
        ay_fwhm = fig.add_subplot(223)
        ay_fwhm.errorbar(self.x_axis_data, self.fit_x_1[4], yerr=self.fit_err_x_1[4],
                      color='#ABABAB', ls='none', marker='o', mec='#424242',
                      mfc='white', mew=1.2)
        ay_fwhm.set_xlabel(self.x_axis_label)
        ay_fwhm.set_ylabel('x fwhm')
        ay_fwhm.set_title("Peak FWHM")
        if self.laser_on != None and self.laser_off != None:
            ay_fwhm.axvline(x=self.laser_on, color='r')
            if len(self.scan_cmd) > 4:
                ay_fwhm.axvline(x=self.laser_off, color='b')
        
        ay_int = fig.add_subplot(224)
        ay_int.errorbar(self.x_axis_data, self.intensity[0, :], yerr = self.intensity[1, :],
                      color='#ABABAB', ls='--', marker='o', mec='#424242',
                      mfc='white', mew=1.2)
        ay_int.set_ylabel("Int (a.u.)")
        ay_int.set_xlabel(self.x_axis_label)
        ay_int.set_title("Intensity Sum")
        if self.laser_on != None and self.laser_off != None:
            ay_int.axvline(x=self.laser_on, color='r')
            if len(self.scan_cmd) > 4:
                ay_int.axvline(x=self.laser_off, color='b')
        plt.tight_layout()
        plt.suptitle(self.x_axis + " scan")
        plt.show()
        self.save_plot(fig, "fit_params")
        
    def save_plot(self, fig, plot_name):
        if self.saving_plot == True:
            if not os.path.exists(self.save_directory):
                os.mkdir(self.save_directory)
            out_filename_fits = self.save_directory + "scan_" + str(self.experiment) + "_" + str(self.x_axis) + "_" + str(self.scan_number) + "_" + str(plot_name) + ".png"
            fig.savefig(out_filename_fits)