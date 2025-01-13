#________________________________________________________________________________________________
#   PROGRAM NAME:  Composition.py
#   PURPOSE:       Calculate Atomic % and Area Under Curve of XPS Data
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
import Control
import Background
import os
import sys

#------------------------------------------------------------------------------------------------
#               GLOBAL CONSTANTS TAKEN FROM Control.py
#------------------------------------------------------------------------------------------------

filename = Control.Return_Filename()

# df = pd.read_csv(filename, skiprows=4, header=None)
df = pd.read_csv(os.path.join(sys.path[0], filename), skiprows=4, header=None)
df = df.apply(pd.to_numeric, errors='coerce')
df = df.to_numpy()
df = df.T
df[0] = np.flip(df[0])
df[1] = np.flip(df[1])

# Max value
maxBE = max(df[0])
# Min value
minBE = min(df[0])
# Length of list
lenBE = len(df[0])
# Step size
BEstep = np.abs(np.round(df[0,0] - df[0,1],3))

# Indexing function
def Index(val):
    ind = np.where(df[0] == val)[0][0]
    return int(ind)

peak_centers = Control.Return_Peak_Centers()
elements = Control.Return_Elements()

#------------------------------------------------------------------------------------------------
#               IMPORT LEFT AND RIGHT REGION
#------------------------------------------------------------------------------------------------

def Import_Region_Data(pk):

    current_directory = os.getcwd()
    directory_reg = os.path.join(current_directory, f'{filename[:-4]}_Regions')

    L_Reg = np.loadtxt(f'{directory_reg}/{filename[:-4]}_{elements[pk].strip()}_LowBE_Region')
    R_Reg = np.loadtxt(f'{directory_reg}/{filename[:-4]}_{elements[pk].strip()}_HighBE_Region')

    return L_Reg, R_Reg

L_Reg, R_Reg = Import_Region_Data(0)

#------------------------------------------------------------------------------------------------
#               LOOK OVER ALL POSSIBLE ENDPOINTS
#------------------------------------------------------------------------------------------------

def Calc_Area_Array(pk,low_reg,high_reg):
    # Loop over low_reg and high_reg
    m = 0
    AUC_array = np.zeros(len(high_reg)*len(low_reg)) #### THESE MAY NEED TO BE FLIPPED, CHECK LATER
    for Emin in low_reg:
        for Emax in high_reg:
            # j = 0
            # BG = 0

            # Curve = j - BG

            # AUC = np.trapz(Curve,np.flip(df[Index(Emax):Index(Emin)+1,0]),dx=BEstep)
            # print(Emin,Emax)

            upperBEdist = Emax - peak_centers[pk]
            lowerBEdist = peak_centers[pk] - Emin
            AUC = Background.area_btwn_spect(pk, upperBEdist, lowerBEdist)

            AUC_array[m] = AUC
            m += 1

    # REMOVE invalid values from array (negative)
    AUC_array = AUC_array[AUC_array > 0]

    return AUC_array

Calc_Area_Array(0,L_Reg,R_Reg)

#------------------------------------------------------------------------------------------------
#               SUBSET 1 SIGMA REGION FROM AREA DISTRIBUTION
#------------------------------------------------------------------------------------------------

def Subset_Area_Histogram(AUC_array):

    mu = np.mean(AUC_array)
    sigma = np.std(AUC_array)

    # Lower bound area
    if (mu - sigma) < np.min(AUC_array):
        Lbound = np.min(AUC_array)
    
    else:
        Lbound = mu - sigma

    # Upper bound area
    if (mu + sigma) > np.max(AUC_array):
        Ubound = np.max(AUC_array)
    
    else:
        Ubound = mu + sigma

    return Lbound, Ubound

#------------------------------------------------------------------------------------------------
#               CALC MIN/MAX FOR ALL ELEMENTS
#------------------------------------------------------------------------------------------------

def Calc_Area_Extrema_Array():

    # Create empty array for values, of the form [Lower_Bound_List,Upper_Bound_List]
    Extrema_Array = np.zeros([len(elements),len(elements)])

    # Loop over all peaks
    for pk in range(len(peak_centers)):

        # Get Area Array
        L_Reg, R_Reg = Import_Region_Data(pk)
        Area_Array = Calc_Area_Array(pk, L_Reg, R_Reg)

        # Get upper and lower bound
        Lbound, Ubound = Subset_Area_Histogram(Area_Array)

        # Add bounds to array
        Extrema_Array[pk,0] = Lbound
        Extrema_Array[pk,1] = Ubound


    Extrema_df = pd.DataFrame({"Transition":elements,"MinAUC":Extrema_Array[:,0],"MaxAUC":Extrema_Array[:,1]})

    return Extrema_df

#------------------------------------------------------------------------------------------------
#               At%, BASED ON EXTREMA
#------------------------------------------------------------------------------------------------

#### Using df from AUC_Extrema, calculate the min and max composition for each element
    ## The max composition will draw from MaxAUC, and force all others to draw from MinAUC
    ## The min composition does the opposite of this

def CalculateErrorBars():

    Extrema_df = Calc_Area_Extrema_Array()
    # print(AUC_df)

    ## Import CRSF
    CRSF = Control.Return_CRSF()
    CRSF = np.array(CRSF,dtype=float)

    MaxComposition = np.zeros(len(elements))
    MinComposition = np.zeros(len(elements))

    ## Loop over all transitions
    for i in range(len(elements)):
        #### Compute maximum %
        ## Find max AUC for element i
        Element_Max_AUC = Extrema_df.iloc[i,2] * (1/CRSF[i])
        # print(Element_Max_AUC)

        ## Find min AUC for all other elements
        ## Turn min AUCs in to np array and drop index i
        ## Divide min AUCs by their CRSF_values
        Other_Elements_Min = (np.delete(np.array(Extrema_df["MinAUC"]),i) / np.delete(np.copy(CRSF),i))

        ## Find %s from 
        Normalization = Element_Max_AUC + np.sum(Other_Elements_Min)
        MaxComposition[i] = Element_Max_AUC / Normalization

        #### Compute minumum %
        ## Find min AUC for element i
        Element_Min_AUC = Extrema_df.iloc[i,1] / CRSF[i]

        ## Find max AUC for all other elements
        ## Turn max AUCs in to np array and drop index i
        ## Divide max AUCs by their CRSF_values
        Other_Elements_Max = (np.delete(np.array(Extrema_df["MaxAUC"]),i) / np.delete(np.copy(CRSF),i))

        ## Find %s from 
        Normalization = Element_Min_AUC + np.sum(Other_Elements_Max)
        MinComposition[i] = Element_Min_AUC / Normalization

    # Save results to a csv file
    ACE_df = pd.DataFrame({'Element':elements,'Min':MinComposition, 'Max':MaxComposition})
    ACE_df.to_csv('ACE.csv',index=False)

    return MaxComposition, MinComposition

#------------------------------------------------------------------------------------------------
#               MAIN FUNCTION
#------------------------------------------------------------------------------------------------

def Main():

    print(CalculateErrorBars())

Main()

