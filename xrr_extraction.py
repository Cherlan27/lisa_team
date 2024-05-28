# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 11:53:48 2023

@author: petersdorf
"""

from make_xrr import make_xrr

raw_file = "./raw"
processed_file = "./processed/reflectivity"
scan_number = range(10, 18+1)
experiment_sample = "test"
save_name_xrr = "{0}_{1:05}".format(experiment_sample, scan_number[0])


xrr_scan = make_xrr()
xrr_scan.scan_numbers = scan_number
xrr_scan.experiment = experiment_sample
xrr_scan.roi = (14, 70, 40, 22)
xrr_scan.data_directory = raw_file
xrr_scan.out_data_directory = processed_file
xrr_scan()
xrr_scan.saving_xrr_data(save_name_xrr)