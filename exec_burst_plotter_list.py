# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 15:19:25 2024

@author: Lukas
"""

from d2_plotter import D2_Plotter
from scan_plotter import Timing_scan
from integrator import Integrator
import pandas as pd
import numpy as np

scan_numbers = range(4208, 4227+1)
experiments = ["nai3moltiming7" for i in scan_numbers]

for i,j in enumerate(list(zip(scan_numbers, experiments))):
    scan = j[0]
    experiment = j[1]
    roi = (1455,126,60,21)
    
    print(scan)
    
    # Script loading
    fast_scan = Timing_scan()
    fast_scan.scan_number = scan
    fast_scan.experiment = experiment
    fast_scan.roi = roi
    fast_scan.save_directory = "H:/supex623/supex623/LTP_2023_11/05mol_fast_scans_2/"
    
    # fast_scan_integrator = Integrator()
    # fast_scan_integrator.image_number = 0
    # fast_scan_integrator.experiment = experiment
    # fast_scan_integrator.scan_number = scan
    # fast_scan_integrator.xlims = [roi[0], roi[0]+roi[2]]
    # fast_scan_integrator.ylims = [roi[1], roi[1]+roi[3]]
    
    # fast_scan_detector = D2_Plotter()
    # fast_scan_detector.scan_number = scan
    # fast_scan_detector.experiment = experiment
    # fast_scan_detector.rois = {"roi": roi}
    
    fast_scan()
    fast_scan.plot_intensity()
    fast_scan.plot_fitted_params_x()
    # fast_scan_detector()
    # fast_scan_integrator()