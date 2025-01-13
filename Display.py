#________________________________________________________________________________________________
#   PROGRAM NAME:  Display.py
#   PURPOSE:       Display Error Bars of XPS Atomic Composition
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
import os
import sys

#------------------------------------------------------------------------------------------------
#               READING IN DATA
#------------------------------------------------------------------------------------------------

data = pd.read_csv("ACE.csv")
data['Min'] = 100*data['Min']
data['Max'] = 100*data['Max']

filename = Control.Return_Filename()

#------------------------------------------------------------------------------------------------
#               PLOTTING DATA
#------------------------------------------------------------------------------------------------

x = np.arange(0,data.shape[0])
plt.scatter(x,data['Min'],marker='_',color='black')
plt.scatter(x,data['Max'],marker='_',color='black')
plt.vlines(x,data['Min'],data['Max'],color='black')
plt.xticks(np.arange(0,data.shape[0]),data['Element'])
plt.title(filename[:-4])
plt.ylabel('At %')
plt.xlabel('Element')
plt.savefig('ErrobarPlot.png',dpi=300)
plt.show()