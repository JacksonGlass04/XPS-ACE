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
import os 
import sys
import argparse

#------------------------------------------------------------------------------------------------
#               READ IN CONFIG TEXT FILE
#------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()


parser.add_argument("-custom_region", "--CustomRegion",
                    default=False,
                    help='This flag allows the user to set custom XPS regions using xpsConfig.txt, use \'-custom_region True\' to save regions')

args = parser.parse_args()

# f = open("./xpsConfig.txt", "r")
f = open(os.path.join(sys.path[0], "xpsConfig.txt"), "r")
lines = f.readlines()
data_lines = []

for line in lines:
    if ">>>" in line:
        data_lines.append(line[3:])

Input_df = pd.read_csv("Data_Table.csv")
Input_df = Input_df.to_numpy()

#------------------------------------------------------------------------------------------------
#               READ AND RETURN CONTROL DATA
#------------------------------------------------------------------------------------------------

##### DEPRECATED METHODS FOR DATA READ IN #####

# # Determine filename
# def Return_Filename():
#     return data_lines[0].strip()

# # Create elements list
# def Return_Elements():
#     return (data_lines[1].strip().split(','))[::2]

# # Create peak centers list
# def Return_Peak_Centers():
#     return np.array((data_lines[1].strip().split(','))[1::2]).astype('float')

##### CURRENT METHODS FOR DATA READ IN #####

# Determine filename
def Return_Filename():
    return data_lines[0].strip()


# # Create elements list
# def Return_Elements():
#     return (data_lines[1].strip().split(','))[::2]

# Create elements list
def Return_Elements():
    return Input_df[:,0]

# Create peak centers list
def Return_Peak_Centers():
    return np.array(Input_df[:,1])

# Obtain CRSF values from csv file
def Return_CRSF():
    return Input_df[:,2]

#------------------------------------------------------------------------------------------------
#               UPDATE MANUALLY SET REGION CSV's
#------------------------------------------------------------------------------------------------

# Update Edited Regions from xpsConfig.txt
def Write_Edited_Regions():

    # If we wish the change the region, we overwrite the file 'Filename_Transition_Low/HighBE_Region'
    Changes_list = data_lines[3].replace(" ", "").replace("\n", "").split(',')
    Transitions = Changes_list[::4]
    LowHigh = Changes_list[1::4]
    StartingBE = np.array(Changes_list[2::4],dtype=float)
    EndingBE = np.array(Changes_list[3::4],dtype=float)

    df = pd.read_csv(Return_Filename(), skiprows=4, header=None, error_bad_lines = False)
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.to_numpy()
    df = df.T
    df[0] = np.flip(df[0])
    df[1] = np.flip(df[1])
    BEstep = np.abs(np.round(df[0,0] - df[0,1],3))

    # Read in changes from xpsConfig.txt, separate in to groups of 4
    for i in range(len(Transitions)):

        # For each group of 4, generate the file that needs replacing
        fname = f'{Return_Filename()[:-4]}_{Transitions[i]}_{LowHigh[i]}BE_Region'

        # Generate the correct list of binding energies
        BE_list = np.arange(StartingBE[i],EndingBE[i]+BEstep,BEstep)

        # Write the list to the correct file
        np.savetxt(f'./{Return_Filename()[:-4]}_Regions/{fname}',BE_list)

if(args.CustomRegion):

    Write_Edited_Regions()
