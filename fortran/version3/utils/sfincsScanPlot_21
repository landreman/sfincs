#!/usr/bin/env python

# This python script plots the output of a SFINCS scan for scanType = 3
# This is a scan in which any 1 numeric input parameter is varied.

outputFilename = "sfincsOutput.h5"

import matplotlib
#import matplotlib.pyplot as plt
import h5py
import numpy
import inspect, math, os
import pickle
import sys
from subprocess import call
from builtins import input

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsScanPlot_21"
print ("This is "+ this_filename)

makePDF = False
for arg in sys.argv:
   if arg.lower()=='pdf':
      makePDF = True

if makePDF:
   matplotlib.use('PDF')
   # Change default font size                                                                                                                               
   font = {'size':6}
   matplotlib.rc('font', **font)
   matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=0.5)
   matplotlib.rc('axes',linewidth=0.7)

import matplotlib.pyplot as plt

numRuns = 0
scanVariableValues = []
outputs = []

xArrayPosition = -1
yArrayPosition = -1

def print_options() :
   # Get a list of the subdirectories:                                                                                                                                        
   directories = filter(os.path.isdir, os.listdir("."))
   atLeastOneOutputExists = False
   
   for directory in directories:
      filename = directory+"/"+outputFilename
      if not os.path.isfile(filename):
         continue
      
      try:
         f = h5py.File(filename,'r')
      except IOError:
         continue
      
      try:
         # Try reading a field that should definitely be present in the output file for any run that completed.                                                              \
         dummy = f["FSABFlow"][()]
      except:
         continue
      
      try:
         finished = f["finished"]
      except KeyError:
         continue
      #if reaching this far the run is ok                                                                                                                                     
      call(["h5dump", "-n", filename])
      atLeastOneOutputExists = True
      break
   
   
   if not atLeastOneOutputExists:
      print ("Error! Cannot read any sfincs output file in subdirectories of this directory.")
      exit(1)




while True:
    #xVariableName = raw_input("Which variable on x-axis? [to quit write quit; to get a list of datasets (quantities to plot) in the output .h5 files write options] ")
    xVariableName = input("Which variable on x-axis? [to quit write quit; to get a list of datasets (quantities to plot) in the output .h5 files write options] ")
    if xVariableName.lower() == "quit":
        exit(0)
    elif xVariableName.lower() == "options":
       print_options()
    else :
        break


if len(str(xVariableName).split("_")) != 1 :
   try:
      xArrayPosition = int(str(xVariableName).split("_")[-1])
      xVariableName = '_'.join(str(xVariableName).split("_")[:-1])
      if xArrayPosition < 1 :
         print ("WARNING: Index has to be larger than 0, using 1 instead.")
         xArrayPosition = 1
   except ValueError:
      xArrayPosition = -1

    
while True:
    #yVariableName = raw_input("Which variable on y-axis? [to quit write quit; to get a list of datasets (quantities to plot) in the output .h5 files write options] ")
    yVariableName = input("Which variable on y-axis? [to quit write quit; to get a list of datasets (quantities to plot) in the output .h5 files write options] ")
    if yVariableName.lower() == "quit":
        exit(0)
    elif yVariableName.lower() == "options":
       print_options()
    else :
        break

if len(str(yVariableName).split("_")) != 1 :
   try:
      yArrayPosition = int(str(yVariableName).split("_")[-1])
      yVariableName = '_'.join(str(yVariableName).split("_")[:-1])
      if yArrayPosition < 1:
         print ("WARNING: Index has to be larger than 0, using 1 instead.")
         yArrayPosition = 1
   except ValueError:
      yArrayPosition = -1

xArrayPositionString = ""
yArrayPositionString = ""


if xArrayPosition != -1 :
   xArrayPositionString = "[" + str(xArrayPosition) + "]"
if yArrayPosition != -1 :
   yArrayPositionString = "[" + str(yArrayPosition) + "]"

print ("Plot " + yVariableName + yArrayPositionString + " as a function of " + xVariableName + xArrayPositionString + ".")


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
       # Try reading the x variable
       if xArrayPosition != -1:
          xVariable = (f[xVariableName][()])[xArrayPosition - 1]
       else :
          xVariable = f[xVariableName][()]

       if f["includePhi1"][()] == f["integerToRepresentTrue"][()] and xVariable.size > 1:
          xVariable = xVariable[-1]

    except:
       print ("Unable to read value of " + xVariableName + xArrayPositionString + " in " + filename + ".")
       continue
    print ("Read " + xVariableName + xArrayPositionString + " = " + str(xVariable) + " from " + filename + ".") 
    try:
       # Try reading the y variable
       if yArrayPosition != -1:
          yVariable = (f[yVariableName][()])[yArrayPosition - 1]
       else :
          yVariable = f[yVariableName][()]

       if f["includePhi1"][()] == f["integerToRepresentTrue"][()] and yVariable.size > 1:
          yVariable = yVariable[-1]

    except:
       print ("Unable to read value of " + yVariableName + yArrayPositionString + " in " + filename + ".")
       continue
    print ("Read " + yVariableName + yArrayPositionString + " = " + str(yVariable) + " from " + filename + ".")
    
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

    scanVariableValues.append(xVariable)
    outputs.append(yVariable)
    
    atLeastOneDirectorySucceeded = True

    numRuns += 1

    print ("Successfully read run in directory "+directory)

if not atLeastOneDirectorySucceeded:
   print ("Error! There do not seem to be any completed sfincs jobs in subdirectories of this directory or none of the output files contain both " + xVariableName + xArrayPositionString +  " and " + yVariableName + yArrayPositionString + ".")
   exit(1)


# Sort:
scanVariableValues_sorted = sorted(scanVariableValues)
outputs_sorted = []
for scanVariableValue in scanVariableValues_sorted:
   outputs_sorted.append(outputs[scanVariableValues.index(scanVariableValue)])
 
outputs_array = numpy.array(outputs_sorted)

xAxisLabels=[]
yAxisLabels=[]

numQuantities = 1
xAxisLabels.append(xVariableName + xArrayPositionString) 
yAxisLabels.append(yVariableName + yArrayPositionString)


logXaxis = False
logYaxis = False

symlogXaxis = False
symlogYaxis = False

#inputXscale = raw_input("Scale on x-axis? [linear (default) / log / symlog] ")
inputXscale = input("Scale on x-axis? [linear (default) / log / symlog] ")
if inputXscale.lower() == "log":
    logXaxis = True
elif inputXscale.lower() == "symlog":
   symlogXaxis = True

#inputYscale = raw_input("Scale on y-axis? [linear (default) / log / symlog] ")
inputYscale = input("Scale on y-axis? [linear (default) / log / symlog] ")
if inputYscale.lower() == "log":
    logYaxis = True
elif inputYscale.lower() == "symlog":
   symlogYaxis = True


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

plt.subplot(numRows,numCols,1)
xdata = scanVariableValues_sorted

#if not logYaxis :
#   ydata = outputs_array
#   yscales.append('linear')
#else:
#   ydata = outputs_array
#   yscales.append('log')

ydata = outputs_array
if logYaxis :
   yscales.append('log')
elif symlogYaxis :
   yscales.append('symlog')
else :
   yscales.append('linear') 
   

xlabels.append(xAxisLabels[0])
ylabels.append(yAxisLabels[0]) 

if logXaxis : 
   xscales.append('log') 
elif symlogXaxis :
   xscales.append('symlog')
else:
   xscales.append('linear') 
   
plt.plot(xdata,ydata,linespec)
try: 
   if symlogXaxis :
      plt.xscale(xscales[-1], linthreshx=numpy.amin(numpy.absolute(xdata)))
   else :
      plt.xscale(xscales[-1])
except: 
   plt.xscale('symlog')  

try: 
   if symlogYaxis :
      plt.yscale(yscales[-1], linthreshy=numpy.amin(numpy.absolute(ydata)))
   else :
      plt.yscale(yscales[-1])
except: 
   plt.yscale('symlog') 


plt.xlabel(xlabels[-1])
plt.ylabel(ylabels[-1])
ymin,ymax = plt.ylim()
ymins.append(ymin)
ymaxs.append(ymax)

outputFile = open('sfincsScan.dat','wb')
scanType=21
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

