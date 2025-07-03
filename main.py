#________________________________________________________________________________________________
#   PROGRAM NAME:  Region.py
#   PURPOSE:       Calculate Region to be used for background subtraction in XPS
#   CREATED:       Jackson Glass 2024
#
#   LABORATORY FOR ASTROPHYSICS AND SURFACE PHYSICS
#   UNIVERSITY OF VIRGINIA
#________________________________________________________________________________________________

#------------------------------------------------------------------------------------------------
#               IMPORT STATEMENTS
#------------------------------------------------------------------------------------------------

import os
import control
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy.signal import argrelmin, argrelmax


#------------------------------------------------------------------------------------------------
#               GLOBAL CONSTANTS TAKEN FROM control.py
#------------------------------------------------------------------------------------------------

filename = control.Return_Filename()
outdir = 'regions'

df = pd.read_csv(os.path.join("input", filename),skiprows=4,header=None,on_bad_lines='skip')
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

peak_centers = control.Return_Peak_Centers()
elements = control.Return_Elements()

#------------------------------------------------------------------------------------------------
#               SMOOTHING/MINIMA
#------------------------------------------------------------------------------------------------

def smooth_signal(data):
    wl = 7
    po = 2
    return savgol_filter(data, window_length=wl, polyorder=po)

def Minima(pk):
    # Distance from the peak, in terms of Binding Energy
    deV = 12

    # Smooth the region of interest
    id_min = control.NearestIdx(df[0], peak_centers[pk]-deV)
    id_max= control.NearestIdx(df[0], peak_centers[pk]+deV)

    x = df[0,id_min:id_max]
    j = df[1,id_min:id_max]
    
    j_smooth = smooth_signal(j)

    # Split in to right data and exclude peak +-2eV
    half_len = int(len(j_smooth) / 2)
    n_steps = int(2 / BEstep)

    x_low = x[0:half_len-n_steps+1]
    x_high = x[half_len+n_steps:]

    lowSpec = j_smooth[0:half_len-n_steps+1]
    highSpec = j_smooth[half_len+n_steps:]

    # Find minima -- special case for no minima
    try:
        # lowMin = np.min(lowSpec)
        lowMinInd = argrelmin(lowSpec)[0][0]
        lowMin = lowSpec[lowMinInd]
        lowBE = x_low[np.where(lowSpec == lowMin)[0][0]]

    except:
        lowMin = 0
        lowBE = 0

    try:
        highMinInd = argrelmin(highSpec)[0][0]
        highMin = highSpec[highMinInd]
        highBE = x_high[np.where(highSpec == highMin)[0][0]]

    except:
        highMin = 0
        highBE = 0

    return lowMin, highMin, lowBE, highBE

#------------------------------------------------------------------------------------------------
#               HALF WIDTH HALF MAX
#------------------------------------------------------------------------------------------------

def HWHM(pk):
    # We only find the HWHM towards the low BE side, and use this for both sides
    deV_low = 8

    id_min = control.NearestIdx(df[0], peak_centers[pk]-deV_low)
    id_max= control.NearestIdx(df[0], peak_centers[pk])

    x = df[0,id_min:control.NearestIdx(df[0], id_max)]
    j = df[1,id_min:control.NearestIdx(df[1], id_max)]

    #### First find the peak maximum, subtract off the minima
    #       Get the minima from the above function
    #       If the minima does not exist, use np.min on the spectrum
    if Minima(pk)[0] != 0:
        j_subtracted = j - Minima(pk)[0]
    elif Minima(pk)[0] == 0:
        j_subtracted = j - np.min(j)
    
    # Calculate the half max by dividing the max by 2
    Half_Max = np.max(j_subtracted) / 2

    # Find the value closest to the half max
    L_index = np.argmin(np.abs(j - Half_Max))

    plt.plot(j_subtracted)

    
    return x[L_index]

#------------------------------------------------------------------------------------------------
#               INFLECTION POINTS FROM DERIVATIVE
#------------------------------------------------------------------------------------------------

def DerivativePoints(pk):
    # Distance from the peak, in terms of Binding Energy
    up = 10

    # Smooth the region of interest
    id_min = control.NearestIdx(df[0], peak_centers[pk]-up)
    id_max= control.NearestIdx(df[0], peak_centers[pk]+up)

    x = df[0,id_min:id_max]
    j = df[1,id_min:id_max]

    
    j_smooth = smooth_signal(j)

    # Calculate the derivative of the smoothed function
    # djdE = np.gradient(j_smooth)

    # TRYING SOMETHING HERE -- SMOOTH THE DERIVATIVE AFTER np.gradient
    djdE = smooth_signal(np.gradient(j))
    # djdE = smooth_signal(djdE)

    wl = 9
    po = 2
    djdE = savgol_filter(np.gradient(j), window_length=wl, polyorder=po)

    # Split in to right data and exclude peak +-2eV
    half_len = int(len(j_smooth) / 2)
    n_steps = int(1)

    x_low = x[0:half_len-n_steps+1]
    x_high = x[half_len:]

    lowSpec = np.abs(djdE[0:half_len-n_steps+1])
    highSpec = np.abs(djdE[half_len:])

    # Find the inflection points where |djdE| == 0
    try:
        lowMinInd = argrelmax(lowSpec)[0][-1]
        lowMin = lowSpec[lowMinInd]
        lowBE = x_low[np.where(lowSpec == lowMin)[0][0]]
    except:
        lowMin = 0
        lowBE = 0

    try:
        highMinInd = argrelmax(highSpec)[0][0]
        highMin = highSpec[highMinInd]
        highBE = x_high[np.where(highSpec == highMin)[0][0]]

    except:
        highMin = 0
        highBE = 0


    return lowBE, highBE

# Function which finds the nth derivative inflection
def DerivativeInflections(n,pk):

    # Peak center
    peak_center = peak_centers[pk]

    # Difference between center and 1st inflection
    lowBE, highBE = DerivativePoints(pk)

    lowDiff = peak_center - lowBE
    highDiff = peak_center - highBE

    lowInflect = peak_center - n*lowDiff
    highInflect = peak_center - n*highDiff
    
    # Returns left and right locations
    return lowInflect, highInflect


#------------------------------------------------------------------------------------------------
#               SUBSETTING BACKGROUND REGION
#------------------------------------------------------------------------------------------------

# First, collect minima and inflection points
def SubsetLeftRegion(pk):
    minimaBE = Minima(pk)[2]
    Inflec3 = DerivativeInflections(3,pk)[0]
    # Check if minima is further than 3 nInf
    # if ...
    if minimaBE < Inflec3:
        # Calculate 5 nInf
        Inflec5 = DerivativeInflections(5,pk)[0]
        # Return region between 3 nInf, 5nInf
        return np.arange(Inflec5,Inflec3+BEstep,BEstep)

    # If min closer than 3 nInf
    Inflec1 = DerivativeInflections(1,pk)[0]
    Inflec2 = DerivativeInflections(2,pk)[0]

    # Create full list of (2,3,1)nInf --> Ordering to preserve 2nInf preference
    Inflections = np.array([Inflec2,Inflec3,Inflec1])

    # Find closest point to min, favoring 2nInf
    diffArr = np.abs(minimaBE - Inflections)
    Inflec_Index = np.argmin(diffArr)

    ChosenInflec = Inflections[Inflec_Index]
    Separation = np.abs(ChosenInflec-minimaBE)
    halfSep = Separation / 2
    halfSepMod = halfSep % BEstep

    # If Inflec closer to pk than Min
    if (ChosenInflec > minimaBE):
        
        # Return the region from M-Inflec_Dist to I+Inflec+Dist
        return np.arange(minimaBE-Separation-halfSepMod,ChosenInflec+Separation+BEstep,BEstep) - Separation/2

    # If Min closer to pk than Inflec
    if (ChosenInflec < minimaBE):

        # Return the region from I-Inflec_Dist to M+Inflec+Dist
        return np.arange(ChosenInflec-Separation-halfSepMod,minimaBE+Separation+BEstep,BEstep) + Separation/2
    
    if (ChosenInflec == minimaBE):

        return np.arange(ChosenInflec-BEstep,ChosenInflec+2*BEstep,BEstep)

    return 0

# First, collect minima and inflection points
def SubsetRightRegion(pk):
    minimaBE = Minima(pk)[3]
    Inflec3 = DerivativeInflections(3,pk)[1]
    # Check if minima is further than 3 nInf
    if minimaBE > Inflec3:
        # Calculate 5 nInf
        Inflec5 = DerivativeInflections(5,pk)[1]
        # Return region between 3 nInf, 5nInf
        return np.arange(Inflec3,Inflec5+BEstep,BEstep)

    # If min closer than 3 nInf
    Inflec1 = DerivativeInflections(1,pk)[1]
    Inflec2 = DerivativeInflections(2,pk)[1]

    # Create full list of (2,3,1)nInf --> Ordering to preserve 2nInf preference
    Inflections = np.array([Inflec2,Inflec3,Inflec1])

    # Find closest point to min, favoring 2nInf
    diffArr = np.abs(minimaBE - Inflections)
    Inflec_Index = np.argmin(diffArr)

    ChosenInflec = Inflections[Inflec_Index]
    Separation = np.abs(ChosenInflec-minimaBE)
    halfSep = Separation / 2
    halfSepMod = halfSep % BEstep
    
    if (ChosenInflec > minimaBE):
        
        # Return the region from M-Inflec_Dist to I+Inflec+Dist
        return np.arange(minimaBE-Separation-halfSepMod,ChosenInflec+Separation+BEstep,BEstep) - halfSep

    # If Min closer to pk than Inflec
    if (ChosenInflec < minimaBE):

        # Return the region from I-Inflec_Dist to M+Inflec+Dist 
        return np.arange(ChosenInflec-Separation-halfSepMod,minimaBE+Separation+BEstep,BEstep) + halfSep
    
    if (ChosenInflec == minimaBE):

        return np.arange(ChosenInflec-BEstep,ChosenInflec+2*BEstep,BEstep)

    return 0

#------------------------------------------------------------------------------------------------
#               PLOTTING REGION
#------------------------------------------------------------------------------------------------

def Plot_Region(pk):

    _, _, lowBE, highBE = Minima(pk)

    lowBEder, highBEder = DerivativePoints(pk)

    low_der2, high_der2 = DerivativeInflections(2,pk)

    low_der3, high_der3 = DerivativeInflections(3,pk)

    plt.axvline(lowBEder,color='navy',linestyle='dashed',label='Peak Inflections')
    plt.axvline(highBEder,color='navy',linestyle='dashed')

    plt.axvline(low_der2,color='navy',linestyle='dashed')
    plt.axvline(high_der2,color='navy',linestyle='dashed')

    plt.axvline(low_der3,color='navy',linestyle='dashed')
    plt.axvline(high_der3,color='navy',linestyle='dashed')

    plt.axvline(peak_centers[pk],linestyle='dotted',color='black')

    if (lowBE != 0):
        plt.axvline(lowBE,color='red',linestyle='dashed')
    if (highBE != 0):
        plt.axvline(highBE,color='red',linestyle='dashed',label='Minima')

    L_Reg = SubsetLeftRegion(pk)
    plt.axvspan(L_Reg[0],L_Reg[-1],alpha=0.5)

    R_Reg = SubsetRightRegion(pk)
    plt.axvspan(R_Reg[0],R_Reg[-1],alpha=0.5,label='Endpoint Region')

    up = 10

    id_min = control.NearestIdx(df[0], peak_centers[pk]-up)
    id_max= control.NearestIdx(df[0], peak_centers[pk]+up)

    x = df[0,id_min:id_max]
    j = df[1,id_min:id_max]

    plt.title(f'{elements[pk]} - {filename[:-4]} - Step = {BEstep}')
    plt.xlabel('Binding Energy')
    plt.ylabel('Counts')

    plt.plot(x,j,color='black')
    # plt.plot(x,np.gradient(j)+np.min(j),color='black')

    plt.legend(loc = 'lower right')
    plt.grid()

    bottom, top = plt.ylim()
    y_max = 0.8*(top-bottom)+bottom

    plt.text(L_Reg[0], y_max, f'{L_Reg[0]:.2f}',fontsize=14,ha='right',rotation='vertical')
    plt.text(L_Reg[-1], y_max, f'{L_Reg[-1]:.2f}',fontsize=14,rotation='vertical')

    plt.text(R_Reg[0], y_max, f'{R_Reg[0]:.2f}',fontsize=14,ha='right',rotation='vertical')
    plt.text(R_Reg[-1], y_max, f'{R_Reg[-1]:.2f}',fontsize=14,rotation='vertical')

    plt.savefig(f'{outdir}/{elements[pk].strip()}_{filename[:-4]}_Region.png')

    np.savetxt(f'{outdir}/{filename[:-4]}_{elements[pk].strip()}_LowBE_Region',L_Reg, fmt='%.1f')
    np.savetxt(f'{outdir}/{filename[:-4]}_{elements[pk].strip()}_HighBE_Region',R_Reg, fmt='%.1f')
    
    plt.tight_layout()

    plt.show()

current_directory = os.getcwd()
for folder in ['regions', 'output']:
    directory = os.path.join(current_directory, folder)
    if not os.path.exists(directory):
        os.makedirs(directory)

def Main():

    for i in range(len(peak_centers)):
        Plot_Region(i)

Main()

