#________________________________________________________________________________________________
#   PROGRAM NAME:  Control.py
#   PURPOSE:       Read in Config.txt File and Prepare Data for Composition Analysis
#   CREATED:       Jackson Glass 2024
#
#   LABORATORY FOR ASTROPHYSICS AND SURFACE PHYSICS
#   UNIVERSITY OF VIRGINIA
#________________________________________________________________________________________________

#------------------------------------------------------------------------------------------------
#               IMPORT STATEMENTS
#------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy.signal import argrelmin, argrelmax

#------------------------------------------------------------------------------------------------
#               READ IN CONFIG TEXT FILE
#------------------------------------------------------------------------------------------------

f = open("xpsConfig.txt", "r")
lines = f.readlines()
data_lines = []

for line in lines:
    if ">>>" in line:
        data_lines.append(line[3:])

#------------------------------------------------------------------------------------------------
#               READ AND RETURN CONTROL DATA
#------------------------------------------------------------------------------------------------

# Determine filename
def Return_Filename():
    return data_lines[0].strip()

# Create elements list
def Return_Elements():
    return (data_lines[1].strip().split(','))[::2]

# Create peak centers list
def Return_Peak_Centers():
    return np.array((data_lines[1].strip().split(','))[1::2]).astype('float')

# Create CRSF values list
def Return_CRSF():
    return data_lines[4].strip().split(',')

#------------------------------------------------------------------------------------------------
#               UPDATE MANUALLY SET REGION CSV's
#------------------------------------------------------------------------------------------------

# Update Edited Regions from xpsConfig.txt
def Write_Edited_Regions():
    return 0