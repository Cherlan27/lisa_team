# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 13:38:08 2022

@author: Petersdorf
"""

# -*- coding: utf-8 -*-
"""
@author: Petersdorf
"""

import numpy
import matplotlib.pyplot as plt
from lmfit import Model
from scipy import sparse
from scipy.sparse.linalg import spsolve

class fitting_bursts():
    """
    Fitting class for burst scan data
    
    Parameters
    ----------
    x_data: Array
        x-axis data  
    intensity: Array
        y-axis intensity data
    frequency: int
        Frequency of the burst scan data
    lim_1: int
        Time value - point of time when laser on
    lim_2: int
        Time value - point of time when first exponential function ends and baseline starts
    lim_3: int
        Time value - point of time when laser off
    lim_4: int
        Time value - point of time when second exponential function ends and baseline starts
    """
        
    def __init__(self, x_data, intensity, frequency, lims, smoothness, baseline_plot):

        self.x_data = x_data
        self.intensity = intensity
        self.frequency = frequency
        
        
        if len(lims) == 4:
            self.exp_num = 2
            self.lim_1 = int(lims[0]*self.frequency)
            self.lim_2 = int(lims[1]*self.frequency)
            self.lim_3 = int(lims[2]*self.frequency)
            self.lim_4 = int(lims[3]*self.frequency)
            self.lim_5 = len(intensity)
        elif len(lims) == 5:
            self.exp_num = 2.5
            self.lim_1 = int(lims[0]*self.frequency)
            self.lim_2 = int(lims[1]*self.frequency)
            self.lim_3 = int(lims[2]*self.frequency)
            self.lim_4 = int(lims[3]*self.frequency)
            self.lim_5 = int(lims[4]*self.frequency)
            self.lim_6 = len(intensity)
        elif len(lims) == 6:
            self.exp_num = 3
            self.lim_1 = int(lims[0]*self.frequency)
            self.lim_2 = int(lims[1]*self.frequency)
            self.lim_3 = int(lims[2]*self.frequency)
            self.lim_4 = int(lims[3]*self.frequency)
            self.lim_5 = int(lims[4]*self.frequency)
            self.lim_6 = int(lims[5]*self.frequency)
            self.lim_7 = len(intensity)
            
        self.baseline_values = []
        self.result_exp = []
        self.baseline_valeus2 = []
        self.result_exp2 = []
        self.baseline_values3 = []
        
        self.fit_report = ""
        self.fit_report2 = ""
        
        self.smoothness = smoothness
        self.baseline_plot = baseline_plot
    
    def baseline_als(self, y, lam, p, niter=10):
        """
        Function for baseline fitting
        
        Parameters
        ----------
        y: Array
            Intensity data  
        lam: float
            Smoothness of the curve
        p: float
            Parameter for y-position of baseline relative to data;
            0.5 => Line is in the middle                
        """
        L = len(y)
        D = sparse.csc_matrix(numpy.diff(numpy.eye(L), 2))
        w = numpy.ones(L)
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = spsolve(Z, w*y)
            w = p * (y > z) + (1-p) * (y < z)
        return z
    
    def exp_invfunc(self, xData, c, tau, b0):
        """
        Exponential model function (Version 1)
        
        Parameters
        ----------
        xData: Array
            x-axis data 
        c: float
            Factor in front of exp
        tau: float
            Relaxation time factor
        b0: float
            x-axis starting parameter   
            
        Result
        ----------
        Array
            y-values for the exponential function
        """
        return c * (1 - numpy.exp(-xData / tau)) + b0

    def exp_invfunc2(self, xData, c, tau, b0):
        """
        Exponential model function (Version 2)
        
        Parameters
        ----------
        xData: Array
            x-axis data 
        c: float
            Factor in front of exp
        tau: float
            Relaxation time factor
        b0: float
            x-axis starting parameter  

        Result
        ----------
        Array
            y-values for the exponential function            
        """
        return ((c * (1 - numpy.exp(-xData / tau))) - b0)*(-1)
        

    def first_exp(self, x_data, intensity, y_start, y_limit, c_value, tau_value, baseline_values):
        """
        First expontential fit
            - Fit function is based on lmfit (lmfit can have fixed parameter for fitting)
            - Tries fit with a fixed start parameter
            - If fit can't be finished because fixation, the fit interval of starting parameter is increased
        
        Parameters
        ----------
        intensity: Array
            intensity data
        y_start: int
            Number in array - Start of the exponential function
        y_limit: int
            Number in array - End of the exponential function
        c: float
            Factor in front of exp
        tau: float
            Relaxation time factor
        baseline_values: Array
            Values of the previous baseline
        
        Result
        ----------
        Array
            y-values for the exponential function    
        """
        exp_func_model2 = Model(self.exp_invfunc2)
        fitted = False
        n = 1
        while fitted == False:
            try:
                params = exp_func_model2.make_params(c=c_value, tau=tau_value, b0 = baseline_values[-1])
                params['b0'].vary = False
                
                params['tau'].min = 0.02
                params['tau'].max = 10950
                
                result = exp_func_model2.fit(intensity[y_start:y_limit], params, xData=numpy.linspace(0,(y_limit-y_start),(y_limit-y_start))/self.frequency)
                fitted = True
            except:
                try:
                    params = exp_func_model2.make_params(c=c_value, tau=tau_value, b0 = baseline_values[-1])
                    params['b0'].vary = True
                    print(baseline_values[-1]*(1-0.01*n))
                    params['b0'].min = baseline_values[-1]*(1-0.01*n)
                    params['b0'].max = baseline_values[-1]*(1+0.01*n)
                    
                    params['tau'].min = 0.02
                    params['tau'].max = 10950
                    
                    result = exp_func_model2.fit(intensity[y_start:y_limit], params, xData=numpy.linspace(0,(y_limit-y_start),(y_limit-y_start))/self.frequency)
                    fitted = True
                except:
                    print("Change of starting parameter, b0.min:" + str(1-0.01*n))
            n = n + 1        
        return result
    
    
    def fit(self):
        """
        Actual fit function - Uses the baseline and exponential models with class parameters
        """
        #baseline_values = self.first_baseline(self.intensity, self.lim_1)
        baseline_values = self.baseline_als(self.intensity[0:self.lim_1], self.smoothness, 0.5)
        self.baseline_values = baseline_values
        
        baseline_values2 = self.baseline_als(self.intensity[self.lim_2:self.lim_3], self.smoothness, 0.5)
        self.baseline_values2 = baseline_values2
        
        baseline_values3 = self.baseline_als(self.intensity[self.lim_4:self.lim_5], self.smoothness, 0.5)
        self.baseline_values3 = baseline_values3
        
        result_exp = self.first_exp(self.x_data, self.intensity, self.lim_1, self.lim_2, 0.005, 5, baseline_values)
        self.result_exp = result_exp.best_fit
        self.fit_report = result_exp.fit_report()
    
        result_exp2 = self.first_exp(self.x_data, self.intensity, self.lim_3, self.lim_4, 0.005, 5, baseline_values2)
        self.result_exp2 = result_exp2.best_fit
        self.fit_report2 = result_exp2.fit_report()
        
        if self.exp_num > 2:
            if self.exp_num == 3:
                baseline_values4 = self.baseline_als(self.intensity[self.lim_6:self.lim_7], self.smoothness, 0.5)
                self.baseline_values4 = baseline_values4
            
                result_exp3 = self.first_exp(self.x_data, self.intensity, self.lim_5, self.lim_6, 0.005, 5, baseline_values3)
                
            if self.exp_num == 2.5:
                print(self.result_exp2[-1])
                result_exp3 = self.first_exp(self.x_data, self.intensity, self.lim_5, self.lim_6, 0.005, 5, self.result_exp2)
                
            self.result_exp3 = result_exp3.best_fit
            self.fit_report3 = result_exp3.fit_report()
    
    def fit_results(self):
        """
        Returns
        -------
        Array with two elements
            It contains the fit reports with the fitted parameters
        """
        if self.exp_num == 2:
            return self.fit_report, self.fit_report2
        elif self.exp_num == 2.5 or self.exp_num  == 3:
            return self.fit_report, self.fit_report2, self.fit_report3
        
    def plot(self):
        """
        Plot function of the fitted function
        """
        plt.plot(self.x_data[self.lim_1:self.lim_2], self.result_exp, color = "orange", linewidth = "4.5")
        plt.plot(self.x_data[self.lim_3:self.lim_4], self.result_exp2, color = "orange", linewidth = "4.5")
        if self.exp_num == 3 or self.exp_num == 2.5:
            plt.plot(self.x_data[self.lim_5:self.lim_6], self.result_exp3, color = "orange", linewidth = "4.5")
        if self.baseline_plot == True:
            plt.plot(self.x_data[0:self.lim_1], self.baseline_values, color = "orange", linewidth = "4.5")
            plt.plot(self.x_data[self.lim_2:self.lim_3], self.baseline_values2, color = "orange", linewidth = "4.5")
            plt.plot(self.x_data[self.lim_4:self.lim_5], self.baseline_values3, color = "orange", linewidth = "4.5")
            if self.exp_num == 3:
                plt.plot(self.x_data[self.lim_6:self.lim_7], self.baseline_values4, color = "orange", linewidth = "4.5")
