# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 14:34:50 2023

@author: petersdorf
"""

from mlreflect.data_generation import Layer, Substrate, AmbientLayer, MultilayerStructure

class structure_creater():
    def __init__(self, number_layer, thickness, roughness, sld):
        self.thickness = thickness
        self.roughness = roughness
        self.sld = sld        
        self.substrate = Substrate('H2O', self.roughness[0], self.sld[0]) 
        self.ambient = AmbientLayer('ambient', 0)
        self.sample = MultilayerStructure()
        self.sample.set_substrate(self.substrate)
        self.sample.set_ambient_layer(self.ambient)
        self.add_layering(number_layer)
    
    def add_layering(self, number_layer):
        for i in range(number_layer):
            j = i + 1
            str_name = "Layer" + str(j)
            layer_gen = Layer(str_name, self.thickness[j], self.roughness[j], self.sld[j])
            self.sample.add_layer(layer_gen)
        
    def get_sample(self):
        return self.sample
        
        
        
        
