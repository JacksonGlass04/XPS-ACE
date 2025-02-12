## BACKGROUND INFORMATION ON THIS PROGRAM

This program was written as part of the Virginia Initiative on Cosmic Origins (VICO) Summer 2024 program, by Jackson Glass (contact at jacksoncglass@gmail.com), for the Laboratory for Astrophysics and Surface Physics at UVA as part of the OSIRIS-REx analysis. 

This program aims to calculate error bars on the atomic composition of an xps sample, specifically the error which is introduced by changing the variable endpoints of the iterative Shirley background, which is commonly used in XPS analysis.


## EXPLANATION OF PROGRAMS AND FILES

- xpsConfig.txt -- This is the config file from which the user sets up their analysis preferences and provides data

- Background.py -- A program for calculating the iterative shirley background of a sample

- Control.py -- A program for reading in data and distributing it among other programs

- Region.py -- A program for finding appropriate endpoint regions for each element being analyzed

- Composition.py -- A program for calculating the composition error bars based on all sets of endpoints

- Display.py -- From the calculated error bars, this program creates a plot of each element's error bar



## HOW TO RUN THE ERROR ANALYSIS PROGRAM
STEP 1: In Multipak, load your data, and select elements you wish to analyze.

#### STEP 2: From multipak, export the spectrum counts, as well as corrected RSF values

    - To export spectrum counts: Multipak has an option to export the spectrum counts as ASCII, this will be saved by default as .csv
    - To export CRSF values: Generate the atomic composition summary table (saved as summary.txt), and reference the table CRSF values

#### STEP 3: Open xpsConfig.txt, we will now set the data to be analyzed, as well as elements and their peak locations

    - Follow instructions in xpsConfig.txt to set both the data to be analyzed, as well as the transitions and their peak centers

#### STEP 4: Run Region.py from the console, i.e. >>> python Region.py

    - This will create and display chosen regions of 'acceptable' endpoints. 
    - If you wish to edit these regions, consult the xpsConfig.txt section on custom regions

#### STEP 5: Run Composition.py from the console, i.e. >>> python Composition.py

    - This will use the chosen regions to calculate atomic percent error bars
    - The results will be printed to the terminal, and saved as a csv in ACE.csv

#### STEP 6: (Optional) To create a plot of the errorbars, run >>> python Display.py

    - This plot can be log or linear, with user-chosen image extension
    - Run 'python Display.py -h' to see all options



## OVERVIEW OF ERROR CALCULATION METHODS

This program calculates the error in atomic composition as a function of the choice of background endpoints across the sample. To calculate the error, we first must find a method for isolating possible locations where background endpoints could be placed.

-Where could they be placed -- Derivatives

A 'good' choice of background endpoint is subjective, but has some common characteristics. A good endpoint occurs near a flatter region, not on the slope of the peak, and generally this region should be a ways away from the inflection of the peak. To model this, we calculate the minima of a smoothed version of the peak, as well as the distance from the peak center to it's inflection (call it say, deltaE). We then generate a region which is symmetric about the distance between the smoothed minima and the nearest multiple of the distance deltaE from the peak. If there is no minima close to the peak, then the region is chosen between 3-5 deltaE.

In the event that this region selection is not optimal, we allow the user to manually override with their own choice of background region.

To turn these sets of endpoints in to error bars, we must calculate the area between the spectrum and background for each choice of background, saving the area to an array. We take the 1-sigma region of this area distribution for each element, and compute the extremal atomic percent with the remaining values. The minimum-% error bar is found by minimizing the area of element X, while maximizing all other element areas. A similar process is done for the maximum-% error bar.


## TROUBLE SHOOTING GUIDE / POTENTIAL ERRORS

- If you run in to an issue you cannot resolve - contact me at jacksoncglass@gmail.com

- Potential common errors include:

    1) Incorrectly specifying parameters or values in xpsConfig.txt
    2) Accidentally accessing data that does not exist (this is more of a problem for high-res spectra)
    3) Forgetting to update CRSF values for your specific XPS setup



## LIST OF DEPENDENCIES

- numpy
- matplotlib.pylot
- scipy.signal
- scipy.integrate
- pandas
- argparse
- sys
- os

