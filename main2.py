# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:36:58 2023

@author: petersdorf
"""

from training_water import water_training

train_h2o = water_training(1, [(0, 1), (0, 15)], [(1,8.2), (1,8.2)], [(7, 18), (7, 18)])
train_h2o.print_structure()
train_h2o.training_generator(2**19)
train_h2o.get_loss(100)
train_h2o.save_loss(0)
train_h2o.fit_model()
train_h2o.save_fit(0)
#train_h2o.get_prediction_tests()

train_h2o.saving_model("trained_1_layer_v2.h5")

#train_h2o = water_training(1, (0, 20), (1.2,4.5), (7, 15))
#train_h2o.print_structure()#
#train_h2o.training_generator(2**18)
#train_h2o.get_loss(200)
#train_h2o.save_loss(1)
#train_h2o.fit_model()
#train_h2o.save_fit(1)
#train_h2o.get_prediction_tests()

#train_h2o.saving_model("trained_1_layer.h5")

#train_h2o = water_training(2, (0, 10), (1.6,3.5), (7, 15))
#train_h2o.print_structure()#
#train_h2o.training_generator(2**18)
#train_h2o.get_loss(100)
#train_h2o.save_loss(2)
#train_h2o.fit_model()
#train_h2o.save_fit(2)
#train_h2o.get_prediction_tests()

#train_h2o.saving_model("trained_2_layer.h5")

