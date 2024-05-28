# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 17:28:02 2023

@author: petersdorf

Test generator for document writing Scan as csv

"""

import time
import pandas as pd
import os.path



def scanning(filename, scan_number):
    df = pd.DataFrame()
    df["Type"] = ["Scan"]
    df["Status"] = [str(scan_number)]
    file_exists = os.path.isfile(filename)
    if file_exists:
        df.to_csv(filename, index = False, header = None, sep = "\t", mode = "a")
    else:
        df.to_csv(filename, index = False, header = True, sep = "\t", mode = "a")
    
def start_scan(filename):
    df = pd.DataFrame()
    df["Type"] = ["Start"]
    df["Status"] = ["Not analysed"]
    file_exists = os.path.isfile(filename)
    if file_exists:
        df.to_csv(filename, index = False, header = None, sep = "\t", mode = "a")
    else:
        df.to_csv(filename, index = False, header = True, sep = "\t", mode = "a")
    
def end_scan(filename):
    df = pd.DataFrame()
    df["Type"] = ["End"]
    df["Status"] = ["Not analysed"]
    file_exists = os.path.isfile(filename)
    if file_exists:
        df.to_csv(filename, index = False, header = None, sep = "\t", mode = "a")
    else:
        df.to_csv(filename, index = False, header = True, sep = "\t", mode = "a")


scan_number = 14

while True:
    time.sleep(2)
    start_scan("xrr_list.dat")
    for i in range(5):
        time.sleep(3)
        scanning("xrr_list.dat",scan_number)
        scan_number = scan_number + 1
    time.sleep(2)
    end_scan("xrr_list.dat")
    
