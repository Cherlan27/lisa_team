# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:00:55 2023

@author: petersdorf
"""

import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pickle
from matplotlib import rcParams

import refnx
from refnx.dataset import ReflectDataset, Data1D
from refnx.analysis import Transform, CurveFitter, Objective, Model, Parameter
from refnx.reflect import SLD, Slab, ReflectModel
from refnx.reflect import Linear, Tanh, Interface, Erf, Sinusoidal, Exponential
from eigene.refnx_transform_fresnel import Transform_Fresnel

class refnx_analyser():
    def __init__(self, analyse_name = "Dummy", file_ending = "layer0"):
        self.save_plot = False
        self.analyse_name = analyse_name
        self.file_ending = file_ending
    
    def load_data(self, file, name):
        """
        Load the analysed data (objective, structure and data)
        Input:
            file: File directory where the objective, strcture and data file is saved
            name: Name of the data files
        """
        objective_file = file + "/" + name + "_object.data"
        structure_file = file + "/" + name + "_structure.data"
        data_file = file + "/" + name + "_data.data"
        
        with open(objective_file, 'rb') as file:
            self.objective = pickle.load(file)
        
        with open(structure_file, 'rb') as file:
            self.structure = pickle.load(file)
        
        with open(data_file, 'rb') as file:
            self.data = pickle.load(file)

        self.data_unpolished = self.data
        
        self.x_unpolished = self.data_unpolished.x
        self.y_err_unpolished = self.data_unpolished.y_err
        self.y_unpolished = self.data_unpolished.y
        self.qc_fit = np.sqrt(16*np.pi*self.objective.parameters[1][-1][1][0].value*10**(-6))
        self.load_type = "Data"

    def load_file(self, file):
        """ 
        Load the preprocessed reflectivity dat file into the class
        Input:
            file: dat file as saved in the make_xrr file with columns for qz, intensity and intensity error
        """
        self.file = file
        self.data_unpolished = ReflectDataset(self.file)
        self.x_unpolished = self.data_unpolished.x
        self.y_err_unpolished = self.data_unpolished.y_err
        self.y_unpolished = self.data_unpolished.y
        self.load_type = "File"
        
    def prepare_data(self, mask = [0.06, 0.8]):
        """ 
        Prepare data with specific mask - fitting range of the experimental data
        Input:
            mask: List of two values - range of qz values which should be considered for fitting
        """
        self.masky = np.logical_and(self.data_unpolished.x < mask[1], self.data_unpolished.x > mask[0])
        
        if self.load_type == "File":
            self.data = ReflectDataset(self.file, mask = self.masky)
        elif self.load_type == "Data":
            self.data = ReflectDataset(self.data , mask = self.masky)
        
    def load_model(self, model):
        """ 
        Load the model into the class
        Input:
            model: refnx reflectivity model 
        """
        self.model = model.get_model()
        self.structure = model.get_structure()
    
    def fit(self, fit_function = 'differential_evolution', transform = "fresnel", use_weights = True, qc = 0.0216, roughness = 0):
        """ 
        Fit of the sample system
        Input:
            fit_function: optimization function (string)
            transform: transformation function for transforming the y-values for fitting (str)
            use_weigths: Use weigths for the fit in dependance of the y error (boolean)
            qc: critical angle (float)
            roughness: roughness value of the surface (float)
        """
        self.objective = Objective(self.model, self.data, transform=Transform_Fresnel(transform, use_weights = use_weights, qc = qc, roughness = roughness))
        self.fitter = CurveFitter(self.objective)
        self.fitter.fit(fit_function)
        self.qc_fit = np.sqrt(16*np.pi*self.objective.parameters[1][-1][1][0].value*10**(-6))
        
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
        
    def plot_sld(self):
        """
        Plot the electron density of the fitted sample system
        """
        self.fig3 = plt.figure(dpi = 300)
        plt.plot(*self.structure.sld_profile())
        plt.ylabel('SLD in $10^{-6} \AA^{-2}$')
        plt.xlabel(r'z in $\AA$')
        if self.save_plot == True:
            plt.savefig(self.save_file + "/" + self.analyse_name + "_" + self.file_ending + "_sld.png")
        plt.show()
        plt.close()
        
    def plot_refl(self, fresnel_roughness = 0):
        """
        Plot the reflectivity with the fresnel reflectivity and the model fit
        """
        self.fig = plt.figure(dpi = 300)
        ax = self.fig.gca()
        ax.set_yscale("log")
        self.x,self.y,self.y_err,self.model_fit = (self.data.x, self.data.y, self.data.y_err, self.objective.model(self.data_unpolished.x))
        
        ax.errorbar(self.x_unpolished, self.fresnel(self.x_unpolished, qc = self.qc_fit, roughness = fresnel_roughness), c="#424242", ls="--")
        ax.errorbar(self.x_unpolished,self.y_unpolished,self.y_err_unpolished,color="blue",marker="o", markerfacecolor = "None")
        ax.errorbar(self.x_unpolished,self.model_fit, color="red")
        
        ax.legend(["Fresnel", "Mes. points", "Fit"])
        plt.xlabel(r'q$_z$ in $\mathring{A}^{-1}$')
        plt.ylabel(r'Reflectivity (norm)')
        
        if self.save_plot == True:
            plt.savefig(self.save_file + "/" + self.analyse_name + "_" + self.file_ending + "_refl.png")
        plt.show()
        plt.close()
        
    def plot_refl_fr(self, fresnel_roughness = 0):
        """
        Plot the R/R_F reflectivity with the the model fit
        """
        self.fig2 = plt.figure(dpi = 300)
        ax = self.fig2.gca()
        self.x,self.y,self.y_err,self.model_fit = (self.data.x, self.data.y, self.data.y_err, self.objective.model(self.data_unpolished.x))

        ax.errorbar(self.x_unpolished,self.y_unpolished/self.fresnel(self.x_unpolished, qc = self.qc_fit, roughness = fresnel_roughness),self.y_err_unpolished/self.fresnel(self.x_unpolished, qc = self.qc_fit, roughness = fresnel_roughness),color="blue",marker="o", markerfacecolor = "None", linestyle = "None")
        ax.errorbar(self.x_unpolished,self.model_fit/self.fresnel(self.x_unpolished, qc = self.qc_fit, roughness = fresnel_roughness), color="red")
        
        ax.legend(["Mes. points", "Fit"])
        plt.xlabel(r'$q_z$ in $\mathring{A}^{-1}$')
        plt.ylabel(r'R/R$_F$')
        plt.ylim([0.05, 2])
        if self.save_plot == True:
            plt.savefig(self.save_file + "/" + self.analyse_name + "_" + self.file_ending + "_refl_fr.png")
        plt.yscale("log")
        plt.show()
        plt.close()
        
    def get_reflectiviy_data(self, fresnel_roughness = 0):
        return [self.x_unpolished, self.y_unpolished/self.fresnel(self.x_unpolished, qc = self.qc_fit, roughness = fresnel_roughness), self.y_err_unpolished/self.fresnel(self.x_unpolished, qc = self.qc_fit, roughness = fresnel_roughness), self.model_fit/self.fresnel(self.x_unpolished, qc = self.qc_fit, roughness = 0)]
        
    def save_plots(self, file):
        """
        Function to activate saving plots
        Input:
            file: file directory for saving plots (str)
        """
        self.save_plot = True
        self.save_file = file
    
    def get_fitted_data(self):
        """
        Get the fitted data
        Output:
            self.x_unpolished: x-values - wave vector qz
            self.model_fit: y-values - intensity
        """
        return self.x_unpolished, self.model_fit
    
    def get_results(self):
        """
        Print the structure with fitted and non-fitted parameters
        """
        for i in range(len(self.objective.parameters[0])):
            print(str(self.objective.parameters[0][i]) + "\n")
        for i in range(len(self.objective.parameters[1])):
            for j in range(len(self.objective.parameters[1][i])):
                print(str(self.objective.parameters[1][i][j]) + "\n")
        
    def save_data(self, file):
        """
        Save all data (objective, strcutre, data) in a specific file
        input:
            file: file directory
            name: name for saving the file
        """
        objective_file = file + "/" + self.analyse_name + "_" + self.file_ending +  "_object.data"
        structure_file = file + "/" + self.analyse_name + "_" + self.file_ending + "_structure.data"
        data_file = file + "/" + self.analyse_name +"_" + self.file_ending +  "_data.data"
        
        filehandler = open(objective_file, 'wb')
        pickle.dump(self.objective, filehandler)
        filehandler.close()
        
        filehandler = open(structure_file, 'wb')
        pickle.dump(self.structure, filehandler)
        filehandler.close()
        
        filehandler = open(data_file, 'wb')
        pickle.dump(self.data_unpolished, filehandler)
        filehandler.close()
        
        parameter_file = file + "/" + self.analyse_name + "_" + self.file_ending +  "_parameter.txt"
        f = open(parameter_file,'w')
        for i in range(len(self.objective.parameters[0])):
            f.write(str(self.objective.parameters[0][i]) + "\n")
        for i in range(len(self.objective.parameters[1])):
            for j in range(len(self.objective.parameters[1][i])):
                f.write(str(self.objective.parameters[1][i][j]) + "\n")
            
        f.write("Calculation of the critical angle based on SLD: " + str(self.qc_fit) + "\n")
        f.close()
    
