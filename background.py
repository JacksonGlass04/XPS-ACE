#________________________________________________________________________________________________
#   PROGRAM NAME:  Background.py
#   PURPOSE:       Calculate and Subtract off Shirely Background in XPS Spectra Given Endpoints
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
from scipy.integrate import simpson
import control
import os

#------------------------------------------------------------------------------------------------
#               GLOBAL CONSTANTS TAKEN FROM control.py
#------------------------------------------------------------------------------------------------

#filename = control.Return_Filename()
#
#df = pd.read_csv(filename,skiprows=4,header=None, on_bad_lines='skip')
#df = df.apply(pd.to_numeric, errors='coerce')
#df = df.to_numpy()
#
#df = df.T
#
#df[0] = np.flip(df[0])
#df[1] = np.flip(df[1])
#
## Max value
#maxBE = max(df[0])
## Min value
#minBE = min(df[0])
## Length of list
#lenBE = len(df[0])
## Step size
#BEstep = np.abs(np.round(df[0,0] - df[0,1],3))

peak_centers = control.Return_Peak_Centers()
elements = control.Return_Elements()

#------------------------------------------------------------------------------------------------
#               FUNCTIONS TO CALC k_n and B_n(E)
#------------------------------------------------------------------------------------------------

# Func to calculate new k value
def calc_kn(j_min, j_max, j, B):
    diff = (j_min - j_max)
    area = simpson((j - B),dx=0.5)
    return (diff / area)

# Func to get new background
def calc_bg(kn, j, E0, bg_in):
    # Filling Bn(E) array via loop
    bg_out = np.zeros(len(bg_in))
    E_max = E0[0]

    j_max = j[0]
    j_min = j[-1]

    i = 0
    for E in E0:
        # Integrate from E to E_max
        # Subset j and bg_in to go from (E, E_max)
        j_int = j[:i]
        bg_in_int = bg_in[:i]

        if i == 0:
            bg_val = 0
        else:
            bg_val = kn * simpson((j_int - bg_in_int),dx=0.5)

        bg_out[i] = bg_val
        i += 1

    bg_out += j_max
    bg_out[0] = j_max

    diff = j_min - bg_out[0]

    return bg_out


#------------------------------------------------------------------------------------------------
#               LOOP FOR ITERATIVE BACKGROUND
#------------------------------------------------------------------------------------------------

# Given a peak location and 2 endpoints, return the background
def iterative_shirley(df, pk, upperBEdist, lowerBEdist):

    Emax = peak_centers[pk] + upperBEdist
    Emin = peak_centers[pk] - lowerBEdist

    id_max = control.NearestIdx(df[0], Emax)
    id_min = control.NearestIdx(df[0], Emin)

    E0 = np.copy(df[0, id_min:id_max+1])
    j0 = np.copy(df[1, id_min:id_max+1])


    j_max = j0[-1]
    j_min = j0[0]

    B0 = np.array([j_min]*len(j0))

    n_iterations = 7

    B = B0

    for i in range(n_iterations):

        # Get new value for k
        k = calc_kn(j_min, j_max, j0, B)

        # Get new background
        # plt.plot(E0,B)
        # plt.plot(E0,j0)
        B = calc_bg(k, j0, E0, B)

    return B

#------------------------------------------------------------------------------------------------
#               FUNCS. FOR AREA UNDER CURVE
#------------------------------------------------------------------------------------------------

# Given a peak location and 2 endpoints, calculate the area between the spectrum and background
def area_btwn_spect(df, pk, upperBEdist, lowerBEdist):

    # Get spectrum values

    Emax = peak_centers[pk] + upperBEdist
    Emin = peak_centers[pk] - lowerBEdist


    id_max = control.NearestIdx(df[0], Emax)
    id_min = control.NearestIdx(df[0], Emin)

    j0 = np.copy(df[1, id_min:id_max+1])

    # Calculate background
    B = iterative_shirley(df, pk, upperBEdist, lowerBEdist)

    Area = simpson((j0 - B),dx=0.5)

    return Area

