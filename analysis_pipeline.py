# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 14:18:26 2023

@author: petersdorf
"""

from refnx_models import zero_layer_refnx, lipid_bilayer_refnx, salt_onelayer_refnx
from refnx_analyser import refnx_analyser
from make_xrr import make_xrr
from predicter_xrr_ml import prediction_sample
import pandas
import time 
import re
import os

class data_analysis():
    def __init__(self, scan_name, scan_number):
        self.scan_name = scan_name
        self.scan_number = scan_number
    # Defintion of sample name
    
    def __call__(self):
        scan_no = self.scan_number[0]
        scan_no = "{:05d}".format(scan_no)
        file_name = "xrr_" + self.scan_name + "_" + str(scan_no)
        
        # Preprocessing of the experimental scan data
        
        xrr_scan = make_xrr()
        
        xrr_scan.data_directory = "K:/SYNCHROTRON/Murphy/2023-11_P08_11016542_hoevelmann_petersdorf/raw/"
        xrr_scan.out_data_directory = "H:/supex623/supex623/LTP_2023_11/reflectivity/"
        
        xrr_scan.experiment = self.scan_name
        xrr_scan.out_experiment = self.scan_name
        xrr_scan.scan_numbers = self.scan_number
        xrr_scan.nom_scan_numbers = self.scan_number
        xrr_scan()
        xrr_scan.saving_xrr_data(file_name)
        xrr_scan.plot_xrr()
        xrr_scan.plot_xrr_crit()
        
        ### -------------------------------------------------------------------------------------------------------- ###
        # mlreflect
        if True:
            mlreflect_pred_layer00 = prediction_sample(file_name, self.scan_number, output_file_end = "layer0" , pathfile = "trained_0_layer_v2.h5")
            mlreflect_pred_layer00()
    
            mlrefelct_file_0layer_name = "H:/supex623/supex623/LTP_2023_11/reflectivity/" + file_name + "_layer0_xrr_fitparams.dat"
            mlreflect_file_0layer_fitparams = pandas.read_csv(mlrefelct_file_0layer_name, sep="\t")
            
            sample0_h2o_roughness = float(mlreflect_file_0layer_fitparams["H2O_roughness"])
            sample0_h2o_sld = float(mlreflect_file_0layer_fitparams["H2O_sld"])
                
            # mlreflect_pred_layer01 = prediction_sample(file_name, self.scan_number, output_file_end = "layer1" , pathfile = "trained_1_layer_v2.h5")
            # mlreflect_pred_layer01()
    
            # mlrefelct_file_1layer_name = "H:/supex623/supex623/LTP_2023_11/reflectivity/" + file_name + "_layer1_xrr_fitparams.dat"
            # mlreflect_file_1layer_fitparams = pandas.read_csv(mlrefelct_file_1layer_name, sep="\t")
    
            
            # sample1_layer_thickness = float(mlreflect_file_1layer_fitparams["Layer1_thickness"])
            # sample1_h2o_roughness = float(mlreflect_file_1layer_fitparams["H2O_roughness"])
            # sample1_layer_roughness = float(mlreflect_file_1layer_fitparams["Layer1_roughness"])
            # sample1_h2o_sld = float(mlreflect_file_1layer_fitparams["H2O_sld"])
            # sample1_layer1_sld = float(mlreflect_file_1layer_fitparams["Layer1_sld"])
            
            ### -------------------------------------------------------------------------------------------------------- ###
            # Sample paraemters defintion - Layer 0
            
            # Define the sample with the initinal values
            layer0_sample = zero_layer_refnx(roughness = sample0_h2o_roughness, sld = sample0_h2o_sld)
            
            # Define the boundaries of the fit values
            # If boundaries are set, the values are set for the fit to vary
            layer0_sample.set_roughness_bounds_h2o(bounds = (sample0_h2o_roughness*(1-0.1),sample0_h2o_roughness*(1+0.1)))
            layer0_sample.set_sld_bounds_h2o(bounds = (sample0_h2o_sld*(1-0.1),sample0_h2o_sld*(1+0.1)))
            
            # Define the model (background and scaling) and set their fit boundaries
            # If boundaries are set, the values are set for the fit to vary
            layer0_sample.set_model(bkg= 10e-10)
            layer0_sample.set_bkg_bounds_model(bounds = (1e-12, 5e-10))
            #h2o_sample.set_scale_bounds_model(bounds = (1,1.1))
            
            # Define the density slicing and the interface types
            layer0_sample.set_slicing(0.1)
            layer0_sample.set_interface(1, "Tanh")
            
            ### -------------------------------------------------------------------------------------------------------- ###
            # # Sample paraemters defintion - Layer 1
            
            # # Define the sample with the initinal values
            # layer1_sample = salt_onelayer_refnx(roughness = sample1_h2o_roughness, sld = sample1_h2o_sld)
            
            # # Define the boundaries of the fit values
            # # If boundaries are set, the values are set for the fit to vary
            # layer1_sample.set_roughness_bounds_h2o(bounds = (sample1_h2o_roughness*(1-0.1),sample1_h2o_roughness*(1+0.1)))
            # layer1_sample.set_sld_bounds_h2o(bounds = (sample1_h2o_sld*(1-0.1),sample1_h2o_sld*(1+0.1)))
    
            # layer1_sample.set_thickness_bounds_layer(bounds = (sample1_layer_thickness*(1-0.1),sample1_layer_thickness*(1+0.1)))
            # layer1_sample.set_roughness_bounds_layer(bounds = (sample1_layer_roughness*(1-0.1),sample1_layer_roughness*(1+0.1)))
            # layer1_sample.set_sld_bounds_layer(bounds = (sample1_layer1_sld*(1-0.1),sample1_layer1_sld*(1+0.1)))
            
            # # Define the model (background and scaling) and set their fit boundaries
            # # If boundaries are set, the values are set for the fit to vary
            # layer1_sample.set_model(bkg= 10e-10)
            # layer1_sample.set_bkg_bounds_model(bounds = (1e-12, 5e-10))
            # #h2o_sample.set_scale_bounds_model(bounds = (1,1.1))
            
            # # Define the density slicing and the interface types
            # layer1_sample.set_slicing(0.1)
            # layer1_sample.set_interface(1, "Tanh")
            # layer1_sample.set_interface(2, "Tanh")
            
            ### -------------------------------------------------------------------------------------------------------- ###
            # Refnx fit procedure, load file and model - Layer 0
            
            analysis = refnx_analyser(file_name, file_ending = "layer0")
            analysis.load_file("H:/supex623/supex623/LTP_2023_11/reflectivity/" + file_name + ".dat")
            
            # Load the model and define the fitting range (mask) of the experimental data
            analysis.load_model(layer0_sample)
            analysis.prepare_data(mask = [0.05, 0.85])
            analysis.fit()
            
            # Plot the reflectivity, the reflectivity devied by Fresnel reflectivity and the SLD profile
            analysis.save_plots("H:/supex623/supex623/LTP_2023_11/reflectivity")
            analysis.plot_refl(fresnel_roughness = 0)
            analysis.plot_refl_fr(fresnel_roughness = 0)
            analysis.plot_sld()
            
            # Print the fitting results and save the data
            analysis.get_results()
            analysis.save_data("H:/supex623/supex623/LTP_2023_11/reflectivity")
    
            ### -------------------------------------------------------------------------------------------------------- ###
            # Refnx fit procedure, load file and model - Layer 1
            
            # analysis = refnx_analyser(file_name, file_ending = "layer1")
            # analysis.load_file("H:/supex623/supex623/LTP_2023_11/reflectivity/" + file_name + ".dat")
            
            # # Load the model and define the fitting range (mask) of the experimental data
            # analysis.load_model(layer1_sample)
            # analysis.prepare_data(mask = [0.05, 0.85])
            # analysis.fit()
            
            # # Plot the reflectivity, the reflectivity devied by Fresnel reflectivity and the SLD profile
            # analysis.save_plots("H:/supex623/supex623/LTP_2023_11/reflectivity")
            # analysis.plot_refl(fresnel_roughness = 2.4)
            # analysis.plot_refl_fr(fresnel_roughness = 2.4)
            # analysis.plot_sld()
            
            # # Print the fitting results and save the data
            # analysis.get_results()
            # analysis.save_data("H:/supex623/supex623/LTP_2023_11/reflectivity")
        else:
            print("Could not be fitted")
        
        
        ### -------------------------------------------------------------------------------------------------------- ###
        # Load analysed data and model
        #analysis = refnx_analyser()
        #analysis.load_data("/mnt/Analyse", "h2o_sample")
        #analysis.plot_refl(fresnel_roughness = 2.4)
        #analysis.plot_refl_fr(fresnel_roughness = 2.4)
        #analysis.plot_sld()
        #analysis.get_results()
        return analysis.get_reflectiviy_data(), mlreflect_pred_layer00.get_reflectiviy_data_ml()

# def analysis_trigger(scan_information):
#     print(scan_information)
#     scan_directory = scan_information[0].replace(scan_information[0].split("/")[-1], "")
#     scan_name = scan_information[0].split("/")[-1].split("_")[0]
#     scan_number_start = scan_information[0].split("/")[-1].split("_")[-1]
#     scan_number_end = scan_information[-1].split("/")[-1].split("_")[-1]
#     scan_number = range(int(scan_number_start), int(scan_number_end))
#     scan_save_name = scan_information[0].split("/")[-1]
    
#     print(scan_directory)
#     print(scan_name)
#     print(scan_number)
    
#     data_analysis(scan_directory, scan_name, scan_number, scan_save_name)
    
#     time.sleep(2)

# def searcher(filename, processed_file, raw_file):
#     while True:
#         xrr = []
#         f = open(filename, "r")
#         scans = f.read().split("\n")
#         for k in range(len(scans)):
#             #xrr.append(scans[k].split(",")[0].split("/")[-1][0:-1])
#             xrr.append(scans[k].replace("[", "").replace("]", "").replace("'", "").split(","))
#         processed_file_list = [f for f in os.listdir(processed_file)]
#         raw_file_list = [f for f in os.listdir(raw_file)]
#         for j in range(len(scans)):
#             data_saved = "xrr_" + xrr[j][0].split("/")[-1] + ".dat"
#             data_raw = xrr[j][0].split("/")[-1]
#             if data_saved not in processed_file_list and data_raw in raw_file_list:
#                 analysis_trigger(xrr[j])
#         time.sleep(5)
