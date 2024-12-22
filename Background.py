#________________________________________________________________________________________________
#   PROGRAM NAME:  Background.py
#   PURPOSE:       Calculate and Subtract off Shirely Background in XPS Spectra Given Endpoints
#   CREATED:       Jackson Glass 2024
#
#   LABORATORY FOR ASTROPHYSICS AND SURFACE PHYSICS
#   UNIVERSITY OF VIRGINIA
#________________________________________________________________________________________________

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy.signal import argrelmin, argrelmax
import Control
import os

filename = Control.Return_Filename()

df = pd.read_csv(filename,skiprows=4,header=None)
df = df.apply(pd.to_numeric, errors='coerce')
df = df.to_numpy()
# df = df.T
# df[0] = np.flip(df[0])
# df[1] = np.flip(df[1])

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
    ind = np.where(df[:,0] == val)[0][0]
    return int(ind)

peak_centers = Control.Return_Peak_Centers()
elements = Control.Return_Elements()

def ShirleyBG(df, Emin,Emax):

    ### Subset the spectrum between Emin and Emax
        ### With how we are reading in the data, j[0] occurs at Emax, j[-1] occurs at Emin

    j = np.copy(df[Index(Emax):Index(Emin)+1,1])

    ### Point Averaging, take j[0] and j[-1] and 
    num_avg_adj = 2
    j[0] = np.average(df[Index(Emax)-num_avg_adj:Index(Emax)+num_avg_adj+1,1])
    j[-1] = np.average(df[Index(Emin)-num_avg_adj:Index(Emin)+num_avg_adj+1,1])

    ### If the intensity at the minimum binding energy is greater than the intensity at the maximum, a straight line is returned
    if df[Index(Emax),1] < df[Index(Emin),1]:
        ## Calculate slope
        # slope = np.abs(df[Index(Emin),1] - df[Index(Emax),1]) / np.abs(df[Index(Emin),0] - df[Index(Emax),0])
        slope = np.abs(j[-1] - j[0]) / np.abs(df[Index(Emin),0] - df[Index(Emax),0])
        # line = slope * array + j(Emax)
        BE_array = np.arange(0,np.abs(df[Index(Emin),0] - df[Index(Emax),0])+BEstep,step=BEstep)
        return (slope * BE_array + (df[Index(Emax),1]))

    ### Create unchanging j0 array with 0 subtraction
    j0 = np.copy(j)

    Bg = np.zeros(len(j))

    k = (j0[0] - j0[-1]) / np.trapz(j0 - j0[-1])

    Eflip = (df[Index(Emax):Index(Emin)+1,0])
    i = 0

    for E in Eflip:
            
        h = np.arange(E,Emax+BEstep,step=BEstep)
        print(h)
        integrand = j[0:len(h)] - j0[-1]
        print(integrand)

        Bg[-len(h)] = np.trapz(integrand)
        # print(np.trapz(integrand))

        i += 1

    Bg = np.flip(np.copy(Bg))
        
    k = (j0[0] - j0[-1]) / np.trapz(j - j[-1])
    j = (j0 - k*Bg + j0[-1])

    BgTemp = (k*Bg + j0[-1])
    diff = -(BgTemp[0] + np.flip(BgTemp)[0])
    # return np.flip(-BgTemp - diff)
    return -BgTemp - diff

print(ShirleyBG(df, peak_centers[1]-8,peak_centers[1]+6))
