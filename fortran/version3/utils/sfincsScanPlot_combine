#!/usr/bin/env python

# This python script combines several the results of multiple sfincsScans,
# plotting the results on the same axes.
# You can either call this script directly, or else call sfincsScanPlot with arguments.
# In the latter case, sfincsScanPlot will just call this script with the same arguments.

import matplotlib
import matplotlib.pyplot as plt
import h5py
import numpy
import inspect, os
import pickle
import sys, subprocess

print ("This is "+ inspect.getfile(inspect.currentframe()))

if len(sys.argv)<2:
    print ("Error! You must call this script with command-line arguments to indicate the directories to plot.")
    exit(1)

matplotlib.rcParams.update({'font.size': 8})

scanTypes = []
numQuantities = []
numRows = []
numCols = []
xdata = []
ydata = []
xlabels = []
ylabels = []
xscales = []
yscales = []
ymins = []
ymaxs = []
linespecs = []
numScans = 0

directories = []

originalDirectory = os.getcwd()

for whichFile in range(1,len(sys.argv)):
    try:
        filename = sys.argv[whichFile]+'/sfincsScan.dat'
        print ("*************************************************")
        print ("Processing directory "+sys.argv[whichFile])
        print ("*************************************************")

        os.chdir(sys.argv[whichFile])
        directories.append(os.getcwd())
        this_script = os.path.abspath(__file__)
        # The "1" below is just an arbitrary command-line argument, which suppresses plotting:
        submitCommand = this_script[:-8] + " noplot"
        #subprocess.call(submitCommand.split(" "))
        os.chdir(originalDirectory)

        inputFile = open(filename,'rb')

        numScans += 1

        # scanTypes.append(pickle.load(inputFile))
        # numQuantities.append(pickle.load(inputFile))
        # numRows.append(pickle.load(inputFile))
        # numCols.append(pickle.load(inputFile))
        # xdata.append(pickle.load(inputFile))
        # ydata.append(pickle.load(inputFile))
        # xlabels.append(pickle.load(inputFile))
        # ylabels.append(pickle.load(inputFile))
        # xscales.append(pickle.load(inputFile))
        # yscales.append(pickle.load(inputFile))
        # ymins.append(pickle.load(inputFile))
        # ymaxs.append(pickle.load(inputFile))

        data = pickle.load(inputFile)
        inputFile.close()

        scanTypes.append(data['scanType'])
        numQuantities.append(data['numQuantities'])
        numRows.append(data['numRows'])
        numCols.append(data['numCols'])
        xdata.append(data['xdata'])
        ydata.append(data['ydata'])
        xlabels.append(data['xlabels'])
        ylabels.append(data['ylabels'])
        xscales.append(data['xscales'])
        yscales.append(data['yscales'])
        ymins.append(data['ymins'])
        ymaxs.append(data['ymaxs'])
        linespecs.append(data['linespec'])

    except:
        os.chdir(originalDirectory)
        raise

    if whichFile>1:
        if scanTypes[-1] != scanTypes[0]:
            print ("ERROR: scanType is not consistent among the directories you requested.")
            exit(1)
        if numQuantities[-1] != numQuantities[0]:
            print ("ERROR: numQuantities is not consistent among the directories you requested.")
            exit(1)
        if numRows[-1] != numRows[0]:
            print ("ERROR: numRows is not consistent among the directories you requested.")
            exit(1)
        if numCols[-1] != numCols[0]:
            print ("ERROR: numCols is not consistent among the directories you requested.")
            exit(1)
        for i in range(numQuantities[0]):
            if xlabels[-1][i] != xlabels[0][i]:
                print ("ERROR: xlabels are not consistent among the directories you requested.")
                exit(1)

# ***************************************************
# Now make the plot
# ***************************************************

fig = plt.figure(figsize=(16,8.5))
fig.patch.set_facecolor('white')

# Set the default color cycle:
#matplotlib.rcParams['axes.color_cycle'] = ['r', 'b', 'c','y','g','k']

for iQuantity in range(numQuantities[0]):
   plt.subplot(numRows[0],numCols[0],iQuantity+1)
   for whichScan in range(numScans):
       plt.plot(xdata[whichScan][iQuantity],ydata[whichScan][iQuantity],linespecs[0],label=directories[whichScan])
   plt.xscale(xscales[0][iQuantity])
   plt.yscale(yscales[0][iQuantity])
   plt.xlabel(xlabels[0][iQuantity])
   plt.ylabel(ylabels[0][iQuantity])
   if len(xlabels[0][iQuantity])==0:
       plt.gca().axes.xaxis.set_ticklabels([])
   if len(ylabels[0][iQuantity])==0:
       plt.gca().axes.yaxis.set_ticklabels([])

   ymin_combined = ymins[0][iQuantity]
   ymax_combined = ymaxs[0][iQuantity]
   for whichScan in range(1,numScans):
       ymin_combined = min(ymin_combined, ymins[whichScan][iQuantity])
       ymax_combined = max(ymax_combined, ymaxs[whichScan][iQuantity])
   plt.ylim(ymin_combined, ymax_combined)
   if iQuantity==0:
       plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, borderaxespad=0., prop={'size':9})

titleString = "Plot generated by "+ inspect.getfile(inspect.currentframe())
plt.figtext(0.5,0.01,titleString,horizontalalignment='center',verticalalignment='bottom')

plt.show()

