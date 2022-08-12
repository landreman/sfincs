#!/usr/bin/env python

# This python script plots the output of a SFINCS scan for scanType = 3
# This is a scan in which any 1 numeric input parameter is varied.

outputFilename = "sfincsOutput.h5"

#import matplotlib.pyplot as plt
import matplotlib
import h5py
import numpy
import inspect, math, os
import pickle
import sys

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsScanPlot_3"
print ("This is "+ this_filename)

makePDF = False
for arg in sys.argv:
   if arg.lower()=='pdf':
      makePDF = True

matplotlib.rcParams.update({'font.size': 8})

if makePDF:
   matplotlib.use('PDF')
   # Change default font size
   font = {'size':6}
   matplotlib.rc('font', **font)
   matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=0.5)
   matplotlib.rc('axes',linewidth=0.7)

import matplotlib.pyplot as plt
#######################

numRuns = 0
scanVariableValues = []
outputs = []

## Load some other required subroutines:
#execfile(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/sfincsScan_common")
exec(open(os.path.dirname(os.path.abspath(__file__))+"/sfincsScan_common").read())

# Determine which variable was scanned by looking at input.namelist in the main directory:
inputFilename = "input.namelist"
with open(inputFilename, 'r') as f:
    inputFile = f.readlines()
scanVariable = readScanVariable("scanVariable","string")

# Was the scan linear or logarithmic?
scanVariableScale = readScanVariable("scanVariableScale","string")
if scanVariableScale == "log" or scanVariableScale == "logarithmic":
   scanScale = "log"
elif scanVariableScale == "lin" or scanVariableScale == "linear":
   scanScale = "linear"
else:
   print ("ERROR: You must set scanVariableScale to either linear or log.")
   exit(1)

def uniq(seq):
   checked = []
   for e in seq:
       if e not in checked:
           checked.append(e)
   return checked

# Get a list of the subdirectories:                                                                        
directories = filter(os.path.isdir, os.listdir("."))

atLeastOneDirectorySucceeded = False

for directory in directories:
    filename = directory+"/"+outputFilename
    if not os.path.isfile(filename):
        print ("Directory "+directory+" does not have a "+outputFilename+" file (yet).")
        continue

    try:
        f = h5py.File(filename,'r')
    except IOError:
        print ("Unable to open "+filename+" even though this file exists.")
        continue

    try:
        # Try reading a field that should definitely be present in the output file for any run that completed.
        dummy = f["FSABFlow"][()]
    except:
        print ("Unable to read "+filename+" even though this file exists.")
        continue

    try:
        finished = f["finished"]
    except KeyError:
        print ("Run in directory "+directory+" does not appear to have finished.")
        continue

    print ("Processing directory "+directory)

    # The expression [()] converts from an h5py dataset to a numpy ndarray:
    integerToRepresentTrue = (f["integerToRepresentTrue"][()])
    inputRadialCoordinateForGradients_new = f["inputRadialCoordinateForGradients"][()]
    RHSMode_new = f["RHSMode"][()]
    Nspecies_new = f["Nspecies"][()]
    Zs = f["Zs"][()]
    includePhi1_new = (f["includePhi1"][()] == integerToRepresentTrue)
    if numRuns == 0:
       inputRadialCoordinateForGradients = inputRadialCoordinateForGradients_new
       RHSMode = RHSMode_new
       Nspecies = Nspecies_new
       includePhi1 = includePhi1_new
    else:
       if inputRadialCoordinateForGradients != inputRadialCoordinateForGradients_new:
          print ("Error! inputRadialCoordinateForGradients is not consistent among runs.")
          exit(1)
       if RHSMode != RHSMode_new:
          print ("Error! RHSMode is not consistent among runs.")
          exit(1)
       if Nspecies != Nspecies_new:
          print ("Error! Nspecies is not consistent among runs.")
          exit(1)
       if includePhi1 != includePhi1_new:
          print ("Error! includePhi1 is not consistent among runs.")
          exit(1)

    if RHSMode != 1 and RHSMode != 2 and RHSMode != 3:
        print ("Error! sfincsScanPlot is not yet set up for RHSMode = "+str(RHSMode))
        exit(1)

    # Read the variable that was scanned from input.namelist rather than from the hdf5 file.
    # This is done in case the variable is not saved in the hdf5 file, and also so we do not have
    # to worry about lower/upper case.
    with open(directory+"/"+inputFilename, 'r') as localInputFile:
       inputFile = localInputFile.readlines()
    scanVariableValues.append(readVariable(scanVariable,"float"))

    FSABFlow = f["FSABFlow"][()]
    FSABjHat = f["FSABjHat"][()]
    particleFlux_vm = f["particleFlux_vm_psiHat"][()]
    heatFlux_vm = f["heatFlux_vm_psiHat"][()]
    elapsedTime = f["elapsed time (s)"][()]
    try:
       sources = f["sources"][()]
    except:
       sources = numpy.array([0,0])
    atLeastOneDirectorySucceeded = True
    if RHSMode>1:
       transportMatrix = f["transportMatrix"][()]
    if includePhi1:
        particleFlux_vd = f["particleFlux_vd_psiHat"][()]
        heatFlux_vd = f["heatFlux_vd_psiHat"][()]

    if RHSMode > 1:
       results = []
       for icol in range(transportMatrix.shape[1]):
          for irow in range(transportMatrix.shape[0]):
             results.append(transportMatrix[irow,icol])
       outputs.append(results)
    else:
       # RHSMode = 1
       if Nspecies==1:
          if includePhi1:
             outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],particleFlux_vd[0,0],heatFlux_vm[0,0],heatFlux_vd[0,0]])
          else:
             outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],heatFlux_vm[0,0]])
             #outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],heatFlux_vm[0,0],elapsedTime])
       else:
          results = []
          iteration = -1
          for ispecies in range(Nspecies):
             if includePhi1:
                results.append(FSABFlow[ispecies,iteration])
                results.append(particleFlux_vm[ispecies,iteration])
                results.append(particleFlux_vd[ispecies,iteration])
                results.append(heatFlux_vm[ispecies,iteration])
                results.append(heatFlux_vd[ispecies,iteration])
             else:
                results.append(FSABFlow[ispecies,iteration])
                results.append(particleFlux_vm[ispecies,iteration])
                results.append(heatFlux_vm[ispecies,iteration])
                results.append(sources[0,ispecies,iteration])
                results.append(sources[1,ispecies,iteration])
          results.append(FSABjHat[iteration])
          outputs.append(results)

    numRuns += 1

    print ("Successfully read run in directory "+directory)

if not atLeastOneDirectorySucceeded:
   print ("Error! There do not seem to be any completed sfincs jobs in subdirectories of this directory.")
   exit(1)


# Sort:
scanVariableValues_sorted = sorted(scanVariableValues)
outputs_sorted = []
for scanVariableValue in scanVariableValues_sorted:
   outputs_sorted.append(outputs[scanVariableValues.index(scanVariableValue)])
 
outputs_array = numpy.array(outputs_sorted)

if RHSMode > 1:
   yAxisLabels=[]
   numQuantities = transportMatrix.shape[0]*transportMatrix.shape[1]
   for irow in range(transportMatrix.shape[0]):
      for icol in range(transportMatrix.shape[1]):
         yAxisLabels.append("L"+str(irow+1)+str(icol+1))
else:
   if Nspecies==1:
      if includePhi1:
         numQuantities = 5
         yAxisLabels=["FSABFlow","particleFlux\nvm_psiHat","particleFlux\nvd_psiHat","heatFlux\nvm_psiHat","heatFlux\nvd_psiHat"]
      else:
         numQuantities = 3
         yAxisLabels=["FSABFlow","particleFlux\nvm_psiHat","heatFlux\nvm_psiHat","elapsed time (s)"]
   else:
      if includePhi1:
         numQuantities = Nspecies*5+1
         yAxisLabels=[]
         for i in range(1,Nspecies+1):
            yAxisLabels.append("FSABFlow (species "+str(i)+")")
            yAxisLabels.append("particleFlux vm_psiHat (species "+str(i)+")")
            yAxisLabels.append("particleFlux vd_psiHat (species "+str(i)+")")
            yAxisLabels.append("heatFlux vm_psiHat (species "+str(i)+")")
            yAxisLabels.append("heatFlux vd_psiHat (species "+str(i)+")")
         yAxisLabels.append("FSABjHat")
        
      else:
         numQuantities = Nspecies*5+1
         yAxisLabels=[]
         for i in range(1,Nspecies+1):
            yAxisLabels.append("FSABFlow (species "+str(i)+")")
            yAxisLabels.append("particleFlux_vm_psiHat (species "+str(i)+")")
            yAxisLabels.append("heatFlux_vm_psiHat (species "+str(i)+")")
            yAxisLabels.append("sources(1) (species "+str(i)+")")
            yAxisLabels.append("sources(2) (species "+str(i)+")")
         yAxisLabels.append("FSABjHat")
            


# ***************************************************
# Now make the plot
# ***************************************************

fig = plt.figure()
fig.patch.set_facecolor('white')

numCols = math.ceil(math.sqrt(numQuantities*1.0))
numRows = math.ceil(numQuantities*1.0/numCols)

xdata = []
ydata = []
xlabels = []
ylabels = []
xscales = []
yscales = []
ymins = []
ymaxs = []
linespec = '.-'

for iQuantity in range(numQuantities):
   plt.subplot(numRows,numCols,iQuantity+1)
   xdata.append(scanVariableValues_sorted)
   if "FSA" in yAxisLabels[iQuantity]:
      # Use linear scale for flows
      ydata.append(outputs_array[:,iQuantity])
      yscales.append('linear')
   else:
      # Use log scale for fluxes
      ydata.append(abs(outputs_array[:,iQuantity]))
      yscales.append('log')
   
   xlabels.append(scanVariable)
   ylabels.append(yAxisLabels[iQuantity])
   if scanScale=="log":
      xscales.append('log')
   else:
      xscales.append('linear') 

   plt.plot(xdata[-1],ydata[-1],linespec)
   plt.xscale(xscales[-1])
   plt.yscale(yscales[-1])
   plt.xlabel(xlabels[-1])
   plt.ylabel(ylabels[-1])
   ymin,ymax = plt.ylim()
   ymins.append(ymin)
   ymaxs.append(ymax)

outputFile = open('sfincsScan.dat','wb')
scanType=3
data = {'scanType':scanType, 'numQuantities':numQuantities, 'numRows':numRows,'numCols':numCols,
        'xdata':xdata, 'ydata':ydata, 'xlabels':xlabels, 'ylabels':ylabels,
        'xscales':xscales, 'yscales':yscales, 'ymins':ymins, 'ymaxs':ymaxs,
        'linespec':linespec}
# pickle.dump(scanType,outputFile)
# pickle.dump(numQuantities,outputFile)
# pickle.dump(numRows,outputFile)
# pickle.dump(numCols,outputFile)
# pickle.dump(xdata,outputFile)
# pickle.dump(ydata,outputFile)
# pickle.dump(xlabels,outputFile)
# pickle.dump(ylabels,outputFile)
# pickle.dump(xscales,outputFile)
# pickle.dump(yscales,outputFile)
# pickle.dump(ymins,outputFile)
# pickle.dump(ymaxs,outputFile)
pickle.dump(data,outputFile)
outputFile.close()

#titleString = "Plot generated by "+ inspect.getfile(inspect.currentframe()) + "\nRun in "+os.getcwd()
titleString = "Plot generated by "+ this_filename + "\nRun in "+os.getcwd()
ax = fig.add_axes([0,0,1,1], frameon=False)
ax.text(0.5,0.99,titleString,horizontalalignment='center',verticalalignment='top')

# If this script was called with any command-line arguments, then do not actually display the plot:
if len(sys.argv) == 1:
    plt.show()


if makePDF:
   if len(sys.argv)>2 : #Use the substituted name as file name
      print ("Writing plot to " + os.getcwd() + "/" + sys.argv[2] + ".pdf.")
      plt.savefig(sys.argv[2] + ".pdf")
   else: #Use script name as file name
      print ("Writing plot to " + os.getcwd() + "/" + os.path.basename(__file__) + ".pdf.")
      plt.savefig(os.path.basename(__file__) + ".pdf")
#######################

