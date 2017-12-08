#! 
# Description:
#************* 
# Python script containing various functions for plotting, integrating etc.
#
# Created by:  Aylwin Iantchenko (07-12-2017)
##########################################################################################################################################
# Import packages
##########################################################################################################################################
import h5py
from scipy.integrate import simps
import numpy as np
from scipy import integrate
from math import sin,cos,radians
import matplotlib.pyplot as plt
##########################################################################################################################################
# Initialize constants and arrays
##########################################################################################################################################
# For plotting
defCol = ['blue','green','red','cyan','purple','yellow','black'\
,'grey','blue','green','red','cyan','purple','yellow','black','grey'] # color
LineType = ['-','--','-o','--p','-^','-s','-d','--o','--p','--^','--s',\
'--d','-o','--p','-^','-s','-d','--o','--p','--^','--s','--d']# Line type


##########################################################################################################################################
# Other useful functions
##########################################################################################################################################  
def readH5(path,varNames,blnScalar):
# function to read *.h5 files
 f = h5py.File(path, 'r')

 data = [] # initiate variable
 for var in varNames:
  if blnScalar:
   dataset = f[var]
   data.append(dataset[()])

  else:
  # Get the data
   data.append(list(f[var]))
 
 return tuple(data)  

def plotter(x,y,Set,colSet):
# Function to simplify plotting
 plt.figure(Set[0])
 ax=plt.gca()

 plt.plot(x,y,Set[1],color=colSet,lw=2,ms = 4,label=Set[2]) # plot things

 ax.set_xlabel(Set[3])
 ax.set_ylabel(Set[4])
 plt.legend()

def getTrigFuns(theta):
# Function to generate sine and cosine
# theta should be in radians
 sineVal,cosineVal = [],[] # initiate list
 
 for val in theta: sineVal.append(sin(val)) # get sine
 for val in theta: cosineVal.append(cos(val)) # get cosine
 return sineVal,cosineVal

def OnedInt(y,x):
 # performs a one dimensional integral over y(x)
 result = np.trapz(y,x = x)
 return result

def SortValsY(x,y):
# Function to sort y and x after inceasing x
 ToSort = zip(x,y)
 ToSort = sorted(ToSort)
 x = [point[0] for point in ToSort]
 y = [point[1] for point in ToSort]
 return y

def SortVals(x,y):
# Function to sort y and x after inceasing x
 ToSort = zip(x,y)
 ToSort = sorted(ToSort)
 x = [point[0] for point in ToSort]
 y = [point[1] for point in ToSort]
 return y,x
