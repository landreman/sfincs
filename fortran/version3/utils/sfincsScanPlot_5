#!/usr/bin/env python

# This python script plots the output of a SFINCS scan for scanType = 5.
# This is a scan of radius, with a scan of Er at each radius.

#import matplotlib.pyplot as plt
import matplotlib
import h5py
import numpy
import inspect, math, os
import subprocess
import pickle
import sys

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsScanPlot_5"
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

Nradii = 0
radii_wish = []
radii_actual = []
radiiWithDuplicates = []
ydata = []
root_types = []
FSABHat2s = []

# Load some other required subroutines:
#execfile(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/sfincsScan_common")
exec(open(os.path.dirname(os.path.abspath(__file__))+"/sfincsScan_common").read())

# Determine which variable was scanned by looking at input.namelist in the main directory:
with open(inputFilename, 'r') as f:
    inputFile = f.readlines()

## inputRadialCoordinate = readVariable("inputRadialCoordinate","int") 
inputRadialCoordinate = readVariable("inputRadialCoordinate","int", False) 
if inputRadialCoordinate == None :
   inputRadialCoordinate = readDefault("inputRadialCoordinate","int")

if inputRadialCoordinate==0:
    radiusName = "psiHat"
elif inputRadialCoordinate==1:
    radiusName = "psiN"
elif inputRadialCoordinate==2:
    radiusName = "rHat"
elif inputRadialCoordinate==3:
    radiusName = "rN"
else:
    print ("Error! Invalid inputRadialCoordinate.")
    exit(1)

#######################
inputRadialCoordinateForGradients = readVariable("inputRadialCoordinateForGradients","int", False) 
if inputRadialCoordinateForGradients == None :
   inputRadialCoordinateForGradients = readDefault("inputRadialCoordinateForGradients","int")

if inputRadialCoordinateForGradients==0:
    generalErName = "dPhiHatdpsiHat"
elif inputRadialCoordinateForGradients==1:
    generalErName = "dPhiHatdpsiN"
elif inputRadialCoordinateForGradients==2:
    generalErName = "dPhiHatdrHat"
elif inputRadialCoordinateForGradients==3:
    generalErName = "dPhiHatdrN"
elif inputRadialCoordinateForGradients==4:
    generalErName = "Er"
else:
    print ("Error! Invalid inputRadialCoordinateForGradients.")
    exit(1)

originalDirectory = os.getcwd()

# Get a list of the subdirectories:                                                                        
directories = sorted(filter(os.path.isdir, os.listdir(".")))

isThisTheFirstRadius = True
for directory in directories:
    try:
        print ("*************************************************")
        print ("Processing directory "+directory)
        print ("*************************************************")

        fullDirectory = originalDirectory + "/" + directory

        os.chdir(fullDirectory)

        # Call sfincsScan in the subdirectory.
        # The "noplot" below is just an arbitrary command-line argument, which suppresses plotting:
        #submitCommand = "sfincsScanPlot noplot"
        #submitCommand = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/sfincsScanPlot_2 noplot"
        submitCommand = os.path.dirname(os.path.abspath(__file__))+"/sfincsScanPlot_2 noplot"
        subprocess.call(submitCommand.split(" "))
        os.chdir(originalDirectory)

        # If the following file does not exist, it may mean there was no solution for ambipolarity.        
        filename = "ambipolarSolutions.dat"
        fullFilename = fullDirectory+"/"+filename
        if not os.path.isfile(fullFilename):
            print ("Directory "+directory+" does not have a "+filename+" file. This may mean a solution has not (yet) been found for ambipolarity.")
            continue

        inputFile = open(fullFilename,'rb')
        Nradii += 1
        data = pickle.load(inputFile)
        inputFile.close()

    except:
        os.chdir(originalDirectory)
        raise

    radii_wish.append(data["radius_wish"])
    radii_actual.append(data["radius_actual"])
    # Add a copy of radius for each root
    radiiWithDuplicates += ([radii_actual[-1]] * len(data["roots"]))

    if isThisTheFirstRadius:
        isThisTheFirstRadius = False
        # We will plot all the quantities plotted by sfincsScanPlot_2, plus extras: the density and temperature of each species, and the ambipolar Er
        oldNumQuantities = data["numQuantities"]
        nHats = data["nHats"]
        Nspecies = len(nHats)
        numQuantities = oldNumQuantities + 1 + 2*Nspecies
        ydata = [[]]*numQuantities

    nHats = data["nHats"]
    THats = data["THats"]
    FSABHat2s.append(["FSABHat2"])
    for i in range(Nspecies):
       numRoots = len(data["roots"])
       ydata[i] = ydata[i] + ([nHats[i]]*numRoots)
       ydata[i+Nspecies] = ydata[i+Nspecies] + ([THats[i]]*numRoots)
    ydata[2*Nspecies] = ydata[2*Nspecies] + list(data["roots"])
    root_types += data["root_types"]
    latestOutputs = data["outputs_ambipolar"]
    for i in range(oldNumQuantities):
        ydata[i+1+2*Nspecies] = ydata[i+1+2*Nspecies] + list(latestOutputs[i])

    print ("*************************************************")
    print ("Finished with directory "+directory)
    print ("*************************************************")

if Nradii<1:
    print ("Unable to read data for an ambipolar solution at any radius. Stopping.")
    exit(0)

print ("Successfully read data for the ambipolar solutions at "+str(Nradii)+" radii.")
#if Nradii<2:
#    print "Since fewer than 2 radii are available, a plot to summarize results vs radius will not be generated."
#    exit(0)

ylabels = []
for i in range(Nspecies):
   ylabels.append('nHat, species '+str(i+1))
for i in range(Nspecies):
   ylabels.append('THat, species '+str(i+1))

ylabels += ["Ambipolar "+generalErName] + data["ylabels"]

print ("root_types:",root_types)

# ***************************************************
# Plot actual vs "wish" radii
# ***************************************************

if False:
   fig = plt.figure(2)
   fig.patch.set_facecolor('white')

   maxRadius = max([max(radii_wish), max(radii_actual), 0])
   minRadius = min([min(radii_wish), min(radii_actual), 0])

   plt.plot([minRadius, maxRadius], [minRadius, maxRadius],':k')
   plt.plot(radii_wish, radii_actual,'.')
   plt.xlabel(radiusName+"_wish")
   plt.ylabel("Actual "+radiusName+" used")
   plt.title("All points should lie exactly on the dashed line!")
   # In matplotlib plots, it is sometimes hard to see the first or last points, so add some margin:
   difference = maxRadius - minRadius
   maxRadius += 0.05*difference
   minRadius -= 0.05*difference
   plt.xlim(minRadius,maxRadius)
   plt.ylim(minRadius,maxRadius)

# ***************************************************
# Sort data by increasing radius
# ***************************************************

perm = numpy.argsort(radii_actual)
radii = list(numpy.array(radii_actual)[perm])
FSABHat2s = list(numpy.array(FSABHat2s)[perm])

perm = numpy.argsort(radiiWithDuplicates)
radiiWithDuplicates = list(numpy.array(radiiWithDuplicates)[perm])
root_types = numpy.array(root_types)[perm]
for iQuantity in range(numQuantities):
   ydata[iQuantity] = list(numpy.array(ydata[iQuantity])[perm])

# ***************************************************
# Now make the plot
# ***************************************************

fig = plt.figure(figsize=(16,8.5))
fig.patch.set_facecolor('white')

numCols = math.ceil(math.sqrt(numQuantities*1.0))
numRows = math.ceil(numQuantities*1.0/numCols)

# We have already populated ydata and ylabels. The rest of the following arrays still need to be populated.
xdata = []
xlabels = []
xscales = []
yscales = []
ymins = []
ymaxs = []
linespec = '.-'

for iQuantity in range(numQuantities):
   plt.subplot(numRows,numCols,iQuantity+1)
   xdata.append(radiiWithDuplicates)
   xlabels.append(radiusName)
   xscales.append('linear') 
   yscales.append('linear') 

   if iQuantity < 2*Nspecies:
      # Density and temperature
      linespec='.-'
      plt.plot(xdata[iQuantity],ydata[iQuantity],linespec)
   else:
      linespec='.'
      x = numpy.array(xdata[iQuantity])
      y = numpy.array(ydata[iQuantity])
      mask = numpy.array([t=='ion' for t in root_types])
      plt.plot(x[mask],y[mask],linespec,color='r',label='ion root')
      mask = numpy.array([t=='unstable' for t in root_types])
      plt.plot(x[mask],y[mask],linespec,color='g',label='unstable root')
      mask = numpy.array([t=='electron' for t in root_types])
      plt.plot(x[mask],y[mask],linespec,color='b',label='electron root')
      mask = numpy.array([t=='unknown' for t in root_types])
      plt.plot(x[mask],y[mask],linespec,color='orange',label='unknown root')
   if iQuantity==2*Nspecies:
      # Plot of Er
      plt.legend(loc=0,fontsize=7)

   plt.xscale(xscales[iQuantity])
   plt.yscale(yscales[iQuantity])
   plt.xlabel(xlabels[iQuantity])
   plt.ylabel(ylabels[iQuantity])
   ymin,ymax = plt.ylim()
   ymins.append(ymin)
   ymaxs.append(ymax)

outputFile = open('sfincsScan.dat','wb')
scanType=5
data = {'scanType':scanType, 'numQuantities':numQuantities, 'numRows':numRows,'numCols':numCols,
        'xdata':xdata, 'ydata':ydata, 'xlabels':xlabels, 'ylabels':ylabels,
        'xscales':xscales, 'yscales':yscales, 'ymins':ymins, 'ymaxs':ymaxs,
        'linespec':linespec, 'radii':radii, 'FSABHat2':FSABHat2s, 'Nspecies':Nspecies, 'root_types':root_types}
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
plt.tight_layout()

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
