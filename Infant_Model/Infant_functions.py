"""
Created on Tue April 20 2020

@author: magdalenarussell
"""
import numpy as np
import sys
import pandas as pd

# WHAT ARE THESE??
# constants
h = 10**(-8)
err = .002

# Load patient data (Assuming APDs are already calculated 
# (but I could work that calculation into this script eventually!))

def load_p_data(patient_names = 'all', filepath = '../_ignore/AllRunsAvg.csv', timescale = 'years'):
    """Load patient data from .csv file.
    Keyword arguments:
    patient_names -- the patient numbers for analysis (either 'all' or a specific patient (i.e. 'pt93'))
    filepath -- the path to the .csv file 
    timescale -- the scale for time (either 'year' or 'monthes')
    """
    data = pd.read_csv(filepath)
    data = data[data.AvgAPD1.notnull()] # Select only rows with average values for APD...
    #if patient_names != 'all':
    #    data = data[data["Sample"] == patient_names]
    data_all = {}
    samples = list(data.Sample.unique())
    for samp in samples:
        v = data[data["Sample"] == samp]
        if timescale == 'years':
            time = list(v['ActualTOI (year)'])
            actualTI = time - 0.125
        elif timescale == 'monthes':
            time = list(v['ActualTOI (Month)'])
            actualTI = time - 1.5
        vload = list(v['VL'])
        apd1 = list(v['AvgAPD1'])
        apd5 = list(v['AvgAPD5'])
        apd10 = list(v['AvgAPD10'])
        pos_start = list(v['HXB2nt_start']) 
        pos_end = list(v['HXB2nt_end'])  
        fragment = list(v['Fragment'])
        actualTI = 
        data_all[samp] = (time, actualTI, vload, apd1, apd5, apd10, pos_start,pos_end, fragment)
    data_all['pat_names'] = samples
    return data_all


