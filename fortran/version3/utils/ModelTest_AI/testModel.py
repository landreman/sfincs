#! 
# Description:
#************* 
# Python script to compare output from SFINCS with theoretical model.
# The script calculates the sine and cosine parts of the perturbed impurity density
# and compares the output with the formulas for n_sine,n_cosine in Section IV.A of the following paper:
# http://aip.scitation.org/doi/abs/10.1063/1.873593

# Input:
#************* 
# 1. (Required) paths to folder containing SFINCS output files (ending on .h5). Multiple entries should be separated by comma. 
# 2. (Required) x-variable to use in plotting. You can use any expression that can be executed in the funtion "getPert" inside modelFunctions.py. Examples: epsilont
# mHats[0], Lperp, Lperp*epsilont/2.0 ...
# 3. (Optional) labels to use in plot. If nothing is given the output will be labelled as data0,data1, etc.  Multiple entries should be separated by comma. 
#
# example:
# >>> python testModel.py /sfincsOut_1,/sfincsOut_2 Lperp label_1,label_2
# This will use all the *.h5 output located in folders  /sfincsOut_1 and /sfincsOut_2, and label them as label1 and label2 respectively. The
# result is plotted vs Lperp.
# 
# Created by:  Aylwin Iantchenko (07-12-2017)
##########################################################################################################################################
# Import packages
##########################################################################################################################################
from sys import argv,exit
from modelFunctions import *
import h5py
import glob
##########################################################################################################################################
# Read input
##########################################################################################################################################

# Read  path
try:
 path = argv[1].split(',')
except:
 exit('%%%%%%%%%%% testModel.py:: Please provide path to folder. Aborting ... %%%%%%%%%%%')

# Read  x-variable name
try:
 varPlot = argv[2]
except:
 exit('%%%%%%%%%%% testModel.py:: Please provide x-variable expression to use in plot. Examples are: Lperp,1.0/Lperp,epsilont,nHats[1] %%%%%%%%%%%')

# Read labels
try:
 labels = argv[3].split(',')
except:
 labels = []
 [labels.append('data' + str(idx)) for idx in range(len(path))]
 print "%%%%%%%%%%% testModel.py:: No labels given. Data will be labelled as data0,data1 etc %%%%%%%%%%%"

##########################################################################################################################################
# Main section: Get all output
##########################################################################################################################################

for idxPath in range(len(path)):


 filenames = glob.glob(path[idxPath]+'*.h5')
  
 Theory_ns,Theory_nc = [],[]
 SFINCS_ns,SFINCS_nc = [],[]
 xAll = []
 for IDXname in range(len(filenames)):

# Get theoretical prediction
  nsTheory, ncTheory, x = getPert(filenames[IDXname],False,[], [], [],[], [], [],[], [], [], [], [],[],[],[],[],[],[],['ns','nc',varPlot],[])

# Get SFINCS output
  nsSFINCS, ncSFINCS = getDensityParts(filenames[IDXname],[],[])

# Save values
  SFINCS_ns.append(nsSFINCS)
  SFINCS_nc.append(ncSFINCS)
  Theory_ns.append(nsTheory)
  Theory_nc.append(ncTheory)
  xAll.append(x)


##########################################################################################################################################
# Plot data
##########################################################################################################################################
# Sort all values
 SFINCS_ns = SortValsY(xAll,SFINCS_ns)
 SFINCS_nc = SortValsY(xAll,SFINCS_nc)
 Theory_nc = SortValsY(xAll,Theory_nc)
 Theory_ns,xAll = SortVals(xAll,Theory_ns)

# plot SFINCS
 plotter(xAll,SFINCS_ns,[1,LineType[idxPath],'nsSFINCS_' + labels[idxPath],varPlot,'Val'],defCol[0])
 plotter(xAll,SFINCS_nc,[2,LineType[idxPath],'ncSFINCS_' + labels[idxPath],varPlot,'Val'],defCol[1])

# plot model
 plotter(xAll,Theory_ns,[1,LineType[idxPath],'nsTheory_' + labels[idxPath],varPlot,'Val'],defCol[2])
 plotter(xAll,Theory_nc,[2,LineType[idxPath],'ncTheory_' + labels[idxPath],varPlot,'Val'],defCol[3])

# Show plot
plt.show()

