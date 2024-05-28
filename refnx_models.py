# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 17:03:36 2023

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

class zero_layer_refnx():
    def __init__(self, roughness = 2.4, sld = 9.8):
        self.air = SLD(0, name='air')
        self.roughness_h2o = roughness
        self.sld_h2o = sld
        
        self.h2o_bulk_sld = SLD(self.sld_h2o, name='h2o_bulk')
        self.h2o_bulk = self.h2o_bulk_sld(0, self.roughness_h2o)
        
    def set_structure(self):
        """ 
        Define the structure of the sample based on slabs devided by interfaces |
        """
        self.structure = self.air | self.h2o_bulk
    
    def set_model(self, bkg = 2e-11, dq = 5.0, scale = 1):
        """ 
        Define the model based on the structure, the noise background bkg, dq and
        the scaling factor of the reflectivity curve
        """
        self.set_structure()
        self.model = ReflectModel(self.structure, bkg=bkg, dq=dq, scale=scale)
        self.model.scale.setp(bounds=(0.7, 1.3), vary=False)
        self.model.bkg.setp(bounds=(1e-17, 9e-8), vary=False)
        
    def set_scale_bounds_model(self, bounds = (0.7,1.3), vary = True):
        """ 
        Change the bounds for the scaling reflectivity factor and if the value should be fitted
        """
        self.model.scale.setp(bounds= bounds, vary= vary)
        
    def set_bkg_bounds_model(self, bounds = (1e-17, 9e-8), vary = True):
        """ 
        Change the bounds for the backgorund of the reflectivity and if the value should be fitted
        """
        self.model.bkg.setp(bounds= bounds, vary= vary)
        
    def set_sld_h2o(self, sld):
        """ 
        Change the initial value for the scattering length density of the water bulk
        """
        self.sld_h2o = sld
        self.h2o_bulk_sld = SLD(sld, name='h2o_bulk')
        self.h2o_bulk = self.h2o_bulk_sld(0, self.roughness_h2o)
    
    def set_roughness_h2o(self, roughness):
        """ 
        Change the initial value for the roughness interface of the water bulk
        """
        self.roughness_h2o = roughness
        self.h2o_bulk = self.h2o_bulk_sld(0, roughness)
    
    def set_roughness_bounds_h2o(self, bounds = (2.2, 3.8), vary = True):
        """ 
        Change the bounds for the roughness interface of the water bulk and if the value should be fitted
        """
        self.h2o_bulk.rough.setp(bounds = bounds, vary = vary)
        
    def set_sld_bounds_h2o(self, bounds = (8.0, 15.5), vary = True):
        """ 
        Change the bounds for the scattering length density of the water bulk and if the value should be fitted
        """
        self.h2o_bulk.sld.real.setp(bounds = bounds, vary = vary)
        
    def set_interface(self, number, interface_type):
        """ 
        Change the interface type
        Input:
            number: Position of the interface
            interface_type: Type of the interface
        """
        if interface_type == "Tanh":
            interface_type_func = Tanh()
        elif interface_type == "Linear":
            interface_type_func = Linear()
        elif interface_type == "Erf":
            interface_type_func = Erf()
        elif interface_type == "Sinusoidal":
            interface_type_func = Sinusoidal()
        elif interface_type == "Exponential":
            interface_type_func = Exponential()
        self.structure[number].interfaces= interface_type_func
    
    def set_slicing(self, parameter):
        """ 
        Change the slicing of the slabs - The smaller the number the higher the number of the density slicing
        """
        self.structure.contract = parameter
    
    def get_model(self):
        """ 
        Get the model of the sample system
        """
        return self.model
    
    # def save_data(self, file, name):
    #     """
    #     Save all data (model and structure) in a specific file
    #     input:
    #         file: file directory
    #         name: name for saving the file
    #     """
    #     structure_file = file + "/" + name + "_structure_sample.data"
    #     model_file = file + "/" + name + "_model_sample.data"
        
    #     filehandler = open(structure_file, 'wb')
    #     pickle.dump(self.structure, filehandler)
    #     filehandler.close()
        
    #     filehandler = open(model_file, 'wb')
    #     pickle.dump(self.model, filehandler)
    #     filehandler.close()
        
    # def load_data(self, file, name):
    #     """
    #     Load the analysed data (objective, structure and data)
    #     Input:
    #         file: File directory where the structure and model file is saved
    #         name: Name of the data files
    #     """
    #     structure_file = file + "/" + name + "_structure_sample.data"
    #     model_file = file + "/" + name + "_model_sample.data"
        
    #     with open(structure_file, 'rb') as file:
    #         self.structure = pickle.load(file)
        
    #     with open(model_file, 'rb') as file:
    #         self.model = pickle.load(file)
            
    def get_structure(self):
        """ 
        Get the structure of the sample system
        """
        return self.structure
    
    
class lipid_bilayer_refnx(zero_layer_refnx):
    def __init__(self,
                 roughness = 2.5,
                 sld = 9.8,
                 roughness_head_value = 2.4,
                 thickness_head_value = 5,
                 sld_head_value = 10 ,
                 roughness_tail_value = 2.4,
                 thickness_tail_value = 10,
                 sld_tail_value = 10,
                 ):
        super().__init__(roughness = roughness, sld = sld)
        
        self.roughness_head_value = roughness_head_value
        self.thickness_head_value = thickness_head_value
        self.sld_head_value = sld_head_value
        
        self.roughness_tail_value = roughness_tail_value
        self.thickness_tail_value = thickness_tail_value
        self.sld_tail_value = sld_tail_value
        
        self.head_sld = SLD(self.sld_head_value, name='head')
        self.head_group = self.head_sld(self.thickness_head_value, self.roughness_head_value)
        
        self.tail_sld = SLD(self.sld_tail_value, name='tail')
        self.tail_group = self.tail_sld(self.thickness_tail_value, self.roughness_tail_value)
        
    
    def set_sld_head(self, sld):
        """ 
        Change the initial value for the scattering length density of the head group
        """
        self.sld_head_value = sld
        self.head_sld = SLD(sld, name='h2o_bulk')
        self.head_group = self.head_sld(self.thickness_head_value, self.roughness_head_value)
    
    def set_roughness_head(self, roughness):
        """ 
        Change the initial value for the roughness of the head group
        """
        self.roughness_h2o = roughness
        self.head_group = self.head_sld(self.thickness_head_value, roughness)
        
    def set_thickness_head(self, thickness):
        """ 
        Change the initial value for the thickness of the head group
        """
        self.thickness_head_value = thickness
        self.h2o_bulk = self.h2o_bulk_sld(thickness, self.roughness_head_value)
        
    def set_sld_bounds_head(self, bounds = (8.0, 15.5), vary = True):
        """ 
        Change the bounds for the scattering length density of the head group and if the value should be fitted
        """
        self.head_group.sld.real.setp(bounds = bounds, vary = vary)
        
    def set_roughness_bounds_head(self, bounds = (2.2, 3.8), vary = True):
        """ 
        Change the bounds for the roughness of the head group and if the value should be fitted
        """
        self.head_group.rough.setp(bounds = bounds, vary = vary)
        
        
    def set_thickness_bounds_head(self, bounds = (2.2, 3.8), vary = True):
        """ 
        Change the bounds for the thickness of the head group and if the value should be fitted
        """
        self.head_group.rough.setp(bounds = bounds, vary = vary)
        
    def set_sld_tail(self, sld):
        """ 
        Change the initial value for the scattering length density of the tail group
        """
        self.sld_tail_value = sld
        self.tail_sld = SLD(sld, name='h2o_bulk')
        self.tail_group = self.tail.sld(self.thickness_tail_value, self.roughness_tail_value)
    
    def set_roughness_tail(self, roughness):
        """ 
        Change the initial value for the roughness of the tail group
        """
        self.roughness_h2o = roughness
        self.tail_group = self.tail_sld(self.thickness_tail_value, roughness)
        
    def set_thickness_tail(self, thickness):
        """ 
        Change the initial value for the thickness of the tail group
        """
        self.thickness_tail_value = thickness
        self.tail_group = self.tail_sld(thickness, self.roughness_tail_value)
        
    def set_sld_bounds_tail(self, bounds = (8.0, 15.5), vary = True):
        """ 
        Change the bounds for the scattering length density of the tail group and if the value should be fitted
        """
        self.tail_group.sld.real.setp(bounds = bounds, vary = vary)
        
    def set_roughness_bounds_tail(self, bounds = (2.2, 3.8), vary = True):
        """ 
        Change the bounds for the roughness of the tail group and if the value should be fitted
        """
        self.tail_group.rough.setp(bounds = bounds, vary = vary)
        
        
    def set_thickness_bounds_tail(self, bounds = (2.2, 3.8), vary = True):
        """ 
        Change the bounds for the thickness of the tail group and if the value should be fitted
        """
        self.tail_group.rough.setp(bounds = bounds, vary = vary)

    def set_structure(self):
        self.structure = self.air | self.tail_group | self.head_group | self.h2o_bulk

        
class salt_onelayer_refnx(zero_layer_refnx):
    def __init__(self,
                 roughness = 2.5,
                 sld = 9.8,
                 roughness_layer_value = 2.4,
                 thickness_layer_value = 5,
                 sld_layer_value = 10 ,
                 ):
        super().__init__(roughness = roughness, sld = sld)
        
        self.roughness_layer_value = roughness_layer_value
        self.thickness_layer_value = thickness_layer_value
        self.sld_layer_value = sld_layer_value
        
        self.layer_sld = SLD(self.sld_layer_value, name='layer')
        self.layer = self.layer_sld(self.thickness_layer_value, self.roughness_layer_value)
        
    
    def set_sld_layer(self, sld):
        """ 
        Change the initial value for the scattering length density of the salt layer
        """
        self.sld_layer_value = sld
        self.layer_sld = SLD(sld, name='h2o_bulk')
        self.layer = self.layer_sld(self.thickness_layer_value, self.roughness_layer_value)
    
    def set_roughness_layer(self, roughness):
        """ 
        Change the initial value for the roughness of the salt layer
        """
        self.roughness_layer_value = roughness
        self.layer = self.layer_sld(self.thickness_layer_value, roughness)
        
    def set_thickness_layer(self, thickness):
        """ 
        Change the initial value for the thickness of the salt layer
        """
        self.thickness_layer_value = thickness
        self.layer = self.layer_sld(thickness, self.roughness_layer_value)
        
    def set_sld_bounds_layer(self, bounds = (8.0, 15.5), vary = True):
        """ 
        Change the bounds for the scattering length density of the salt layer and if the value should be fitted
        """
        self.layer.sld.real.setp(bounds = bounds, vary = vary)
        
    def set_roughness_bounds_layer(self, bounds = (2.2, 3.8), vary = True):
        """ 
        Change the bounds for the roughness interface of the salt layer  and if the value should be fitted
        """
        self.layer.rough.setp(bounds = bounds, vary = vary)
        
        
    def set_thickness_bounds_layer(self, bounds = (2.2, 3.8), vary = True):
        """ 
        Change the bounds for the thickness of the interface of the salt layer  and if the value should be fitted
        """
        self.layer.rough.setp(bounds = bounds, vary = vary)
        

    def set_structure(self):
        self.structure = self.air |  self.layer | self.h2o_bulk

        
