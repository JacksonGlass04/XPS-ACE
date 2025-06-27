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

import os
import argparse
import control
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
 
parser = argparse.ArgumentParser()


parser.add_argument("-fform", "--fileformat",
                    default="png",
                    help='Exports plot in other than default png format, such as svg, pdf, jpeg')

parser.add_argument("-scale", "--plotscale",
                    default="linear",
                    help='Determines the scale of the plot, takes the arguments \'log\' or \'linear\'')

args = parser.parse_args()

casename = control.Return_Casename()
outdir = 'output'

#------------------------------------------------------------------------------------------------
#               READING IN DATA
#------------------------------------------------------------------------------------------------

data = pd.read_csv(f"{casename}_comp.csv")
data['Min'] = data['Min']
data['Max'] = data['Max']

#------------------------------------------------------------------------------------------------
#               PLOTTING DATA
#------------------------------------------------------------------------------------------------

x = np.arange(0,data.shape[0])
plt.scatter(x,data['Min'],marker='_',color='black')
plt.scatter(x,data['Max'],marker='_',color='black')
plt.vlines(x,data['Min'],data['Max'],color='black')
plt.xticks(np.arange(0,data.shape[0]),data['Element'])
plt.yscale(args.plotscale)
plt.title(f'Atomic% {args.plotscale}')
plt.ylabel('At %')
plt.xlabel('Element')
plt.savefig(os.path.join(outdir,f'{casename}_error_{args.plotscale}.{args.fileformat}'),dpi=300)
plt.show()