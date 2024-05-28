# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 03:24:26 2023

@author: petersdorf
"""

7#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy
from scipy.optimize import curve_fit
            
            
class Absorber(object):
    
    def __init__(self, values={x: 10.0 for x in range(1, 10)}):
        super(Absorber, self).__init__()
        self._values = values
        self._data = []
        
    def __call__(self, n):
        if n == 0:
            return 1.0
        else:
            return self._values[n] * self.__call__(n-1)
            
    def add_dataset(self, abs_value, qz, intensity):
        """
        Add dataset for absorber factor determination from overlaps.
        """
        self._data.append((int(abs_value+0.05), qz, intensity))
    
    def function1(self, x_1, a_1, b_1, c_1): # not all parameters are used here, c is shared
            return a_1+b_1*x_1+c_1*x_1**2

    def function2(self, x_2, a_1, b_1, c_1, u): # not all parameters are used here, c is shared
            return (a_1+b_1*x_2+c_1*x_2**2)*u
    
    def combinedFunction(self, comboData, a, b, c, u):
        # single data reference passed in, extract separate data
        extract1 = comboData[:len(self.x1)] # first data
        extract2 = comboData[len(self.x1):] # second data

        result1 = self.function1(extract1, a, b, c)
        result2 = self.function2(extract2, a, b, c, u)

        return numpy.append(result1, result2)
    
    def calculate_from_overlaps(self):
        temp = {}
        zero_test = 0
        # compare each dataset with others with abs_value-1
        for n in range(len(self._data)):
            abs_value, qz, intensity = self._data[n]
            for m in range(len(self._data)):
                abs_value_a, qz_a, intensity_a = self._data[m]
                if abs_value_a == abs_value - 1 and zero_test == 0:
                    
                    fit_mask = numpy.where(qz >= qz_a[0])[0]
                    if fit_mask.size == 1:
                        fit_mask = numpy.append(numpy.array([fit_mask[0]-1]), fit_mask)
                    fit_mask_a = numpy.where(qz_a <= qz[len(qz)-1])[0]
                    if fit_mask_a.size == 1:
                        fit_mask_a = numpy.append(fit_mask_a, numpy.array([1]))
                    
                    self.x1 = numpy.array(qz)[fit_mask]
                    self.x2 = numpy.array(qz_a)[fit_mask_a]
                    self.y1 = numpy.array(intensity)[fit_mask]
                    self.y2 = numpy.array(intensity_a)[fit_mask_a]
                    
                    comboY = numpy.append(self.y1, self.y2)
                    comboX = numpy.append(self.x1, self.x2)
                    
                    initialParameters = numpy.array([1.0, 1.0, 1.0, 1.0])
                    fittedParameters, pcov = curve_fit(self.combinedFunction, comboX, comboY, initialParameters)
                    
                    if abs_value not in temp:
                        temp[abs_value] = []
                    temp[abs_value].append(fittedParameters[-1])
                    
                    if abs_value_a == 0:
                        zero_test =+ 1
        
        result = {n: temp[n][0] for n in temp.keys()}
        print(result)
        self._values.update(result)