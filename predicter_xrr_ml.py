"""
Created on Mon Feb 13 15:56:03 2023

@author: petersdorf
"""

import mlreflect
from mlreflect.utils import check_gpu
import matplotlib.pyplot as plt
import numpy as np
from mlreflect.data_generation import Layer, Substrate, AmbientLayer, MultilayerStructure
from mlreflect.training import Trainer
from mlreflect.data_generation import ReflectivityGenerator
from mlreflect.curve_fitter import CurveFitter
import pandas as pd
from mlreflect.models import TrainedModel
import pandas
import os
from scipy import sparse
from scipy.sparse.linalg import spsolve

class prediction_sample():
    def __init__(self, exp_file, scan_number, output_file_end = "layer0" , pathfile = "trained_0_layer.h5", save_directory = "H:/supex623/supex623/LTP_2023_11/reflectivity/"):
        df = pandas.read_csv("H:/supex623/supex623/LTP_2023_11/reflectivity/" + exp_file + ".dat", sep = "\t")
        self.qz = df["//qz"]
        self.file_name = exp_file
        self.inties = abs(df["intensity_normalized"])
        self.inties_e = df["e_intensity_normalized"]
        self.path_file = pathfile
        self.scan_number = scan_number
        self.save_directory = save_directory
        self.output_file_end = output_file_end
        print(self.qz)

    def fresnel(self, qz, qc = 0.0216, roughness=0):
        """
        Fresnel helper function for plotting
        Input:
            qz: wave vector values (list)
            qc: critical angle (float)
            roughness: surface roughness (float)
        """
        return (np.exp(-qz**2 * roughness**2) *
                abs((qz - np.sqrt((qz**2 - qc**2) + 0j)) /
                (qz + np.sqrt((qz**2 - qc**2)+0j)))**2)
    
    def __call__(self):
        model = TrainedModel()
        model.from_file(self.path_file)
        curve_fitter = CurveFitter(model)
        
        num_fraction_bound = len(model.sample.layers)*3 + 2
        fraction_bounds = tuple(list(np.ones(num_fraction_bound) * 0.5))
        try:
            experimental_fit_output = curve_fitter.fit_curve(self.inties, self.qz, polish = True, optimize_q = True, fraction_bounds = fraction_bounds)
        except RuntimeError:
            experimental_fit_output = curve_fitter.fit_curve(self.inties[0:-5], self.qz[0:-5], polish = True, optimize_q = True, fraction_bounds = fraction_bounds)
        pred_experimental_reflectivity = experimental_fit_output["predicted_reflectivity"]
        self.ml_reflectivity_result = pred_experimental_reflectivity[0]
        pred_experimental_test_labels = experimental_fit_output["predicted_parameters"]
        
        fig = plt.figure(dpi = 300)
        #plt.title("XRR - scan no. " + str(self.scan_number))
        plt.semilogy(self.qz, self.inties, 'o', markerfacecolor = "white", markeredgecolor = "blue", label = "Experiment")
        plt.semilogy(self.qz, pred_experimental_reflectivity[0], label = "Prediction", color="red")
        
        plt.legend()
        plt.xlabel(r"q$_z$ in $\mathring{A}^{-1}$")
        plt.ylabel("Reflectivity (norm)")
        plt.show()
        if os.path.exists(self.save_directory):
            save_directory_comp = self.save_directory + self.file_name + "_" + self.output_file_end + ".png"
        else:
            os.makedirs(self.save_directory)
            save_directory_comp = self.save_directory + self.file_name + "_" + self.output_file_end + ".png"
        plt.savefig(save_directory_comp)
        
        fig2 = plt.figure(dpi = 300)
        #plt.title("XRR - scan no. " + str(self.scan_number))
        plt.semilogy(self.qz, self.inties/self.fresnel(self.qz), 'o', markerfacecolor = "white", markeredgecolor = "blue", label = "Experiment")
        plt.semilogy(self.qz, pred_experimental_reflectivity[0]/self.fresnel(self.qz), label = "Prediction", color="red")
        plt.ylim([0.05, 2])
        plt.legend()
        plt.xlabel(r"$q_z$ in $\mathring{A}^{-1}$")
        plt.ylabel(r"R/R$_F$")
        plt.show()
        
        if os.path.exists(self.save_directory):
            save_directory_comp = self.save_directory + self.file_name + "_" + self.output_file_end + ".png"
        else:
            os.makedirs(self.save_directory)
            save_directory_comp = self.save_directory + self.file_name + "_" + self.output_file_end + ".png"
        plt.savefig(save_directory_comp)
        
        out_filename = self.save_directory + self.file_name + "_" + self.output_file_end + "_xrr_data.dat"
        df = pandas.DataFrame()
        df["qz"] = self.qz
        df["intensity_normalized"] = self.inties
        df["e_intensity_normalized"] = self.inties_e
        df["intensity_fit"] = pred_experimental_reflectivity[0]
        df.to_csv(out_filename, sep="\t", index=False)
        
        out_filename = self.save_directory + self.file_name + "_" + self.output_file_end + "_xrr_fitparams.dat"
        df = pandas.DataFrame()
        print(pred_experimental_test_labels)
        for keys, values in pred_experimental_test_labels.items():
            df[keys] = [float(str(values).split("\t")[0].split("\n")[0].split(" ")[4])]
        df.to_csv(out_filename, sep="\t", index=False)
        
        print("Reflectivity scan: " + str(self.scan_number))
        for keys, values in pred_experimental_test_labels.items():
            if str(keys).find("Name") == -1:
                print(str(keys) + ": " + str(values).split("\t")[0].split("\n")[0].split(" ")[4])

    def get_reflectiviy_data_ml(self):
        return [self.qz, self.inties, self.ml_reflectivity_result[0]]