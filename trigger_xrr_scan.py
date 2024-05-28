# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 17:42:43 2023

@author: petersdorf

Check Scans und trigger ML

"""

import numpy as np
import pandas as pd
import time

def searcher(filename):
    scan_numbers = []
    scan_start_analysing = False
    scan_end_analysing = False
    while True:
        df = pd.read_csv(filename, sep="\t")
        for i in range(len(df["Type"])):
            if df["Type"][i] == "Start" and df["Status"][i] == "Not analysed":
                pos_start_value = i
                scan_start_analysing = True
            if df["Type"][i] == "End" and df["Status"][i] == "Not analysed":
                 pos_end_value = i
                 scan_end_analysing = True
            if scan_start_analysing == True and scan_end_analysing == True:
                for j in range(pos_start_value+1, pos_end_value):
                    scan_numbers.append(int(df["Status"][j]))
                scan_status_pos = [pos_start_value, pos_end_value]
                
                mlreflect_trigger(scan_numbers, scan_status_pos)
                
                df = pd.read_csv(filename, sep="\t")
                df.loc[pos_start_value, "Status"] = "Analysed"
                df.loc[pos_end_value, "Status"] = "Analysed"
                df.to_csv(filename, header = True, index = False, sep="\t") 
                
                scan_numbers = []
                scan_start_analysing = False
                scan_end_analysing = False
                
def mlreflect_trigger(scan_numbers, scan_status_pos):
    print(scan_numbers)
    print(scan_status_pos)
    
    time.sleep(5)
    # Extraction of data
    # mlreflect prediction
    
    

searcher("xrr_list.dat")

