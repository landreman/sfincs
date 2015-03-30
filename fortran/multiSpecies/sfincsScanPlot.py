#!/usr/bin/env python

# This python script plots the output of a SFINCS convergence scan.

outputFilename = "sfincsOutput.h5"

hideRepeatedScales = True
#hideRepeatedScales = False

import matplotlib.pyplot as plt
import h5py
import numpy
import inspect, math, os

print "This is "+ inspect.getfile(inspect.currentframe())

numRuns = 0
Nthetas = []
Nzetas = []
Nxis = []
NLs = []
Nxs = []
NxPotentialsPerVths = []
xMaxs = []
solverTolerances = []

outputs = []

baseCaseIndex = -1

def uniq(seq):
   checked = []
   for e in seq:
       if e not in checked:
           checked.append(e)
   return checked

# Get a list of the subdirectories:                                                                        
directories = filter(os.path.isdir, os.listdir("."))

for directory in directories:
    filename = directory+"/"+outputFilename
    if not os.path.isfile(filename):
        print "Directory "+directory+" does not have a "+outputFilename+" file (yet)."
        continue

    try:
        f = h5py.File(filename,'r')
    except IOError:
        print "Unable to open "+filename+" even though this file exists."
        continue

    try:
        # Try reading a field that should definitely be present in the output file for any run that completed.
        dummy = f["/run  1/FSABFlow"][()]
    except:
        print "Unable to read "+filename+" even though this file exists."
        continue

    print "Processing directory "+directory

    # The expression [()] converts from an h5py dataset to a numpy ndarray:
    integerToRepresentTrue = (f["/run  1/integerToRepresentTrue"][()])
    Nthetas.append(f["/run  1/Ntheta"][()])
    Nzetas.append(f["/run  1/Nzeta"][()])
    Nxis.append(f["/run  1/Nxi"][()])
    NLs.append(f["/run  1/NL"][()])
    Nxs.append(f["/run  1/Nx"][()])
    NxPotentialsPerVths.append(f["/run  1/NxPotentialsPerVth"][()])
    xMaxs.append(f["/run  1/xMax"][()])
    solverTolerances.append(f["/run  1/solverTolerance"][()])

    results = f["/run  1/ArrayFirstSpeciesParticleFluxCoefficients"][()]
    outputs.append(results)
    if directory=="00":
        baseCaseIndex=numRuns
    numRuns += 1

    print "Successfully read run in directory "+directory
 

if baseCaseIndex<0:
    print "Error! No baseCase run available."
    exit(1)

numRows = 3
yAxisLabels = ["Lzz11","Lzi11","Lzz12+Lzi12"]
            
parametersToVary = []
abscissae = []
convergeds = []
quantities = []

# Check whether Ntheta was scanned:
data = Nthetas
values = uniq(data)
values.sort()
if len(values)>1:
    # Yes, a scan was performed in this variable.
    parametersToVary.append("Ntheta")
    abscissae.append(values)
    # Ensure only 1 value is repeated:
    numRepeatedValues = 0
    for value in values:
        if data.count(value)>1:
            numRepeatedValues += 1
            convergedValue = value
    if numRepeatedValues>1:
        print "Error! Multiple values of "+parametersToVary[-1]+" were repeated."
        exit(1)
    if numRepeatedValues==0:
       convergedValue = data[baseCaseIndex]
    thisScanResults = []
    for value in values:
        if value==convergedValue:
            index = baseCaseIndex
        else:
            index = data.index(value)
        thisScanResults.append(outputs[index])
    print "A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue)
    print "Values of the scanned parameter: ",abscissae[-1]
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

# Check whether Nzeta was scanned:
data = Nzetas
values = uniq(data)
values.sort()
if len(values)>1:
    # Yes, a scan was performed in this variable.
    parametersToVary.append("Nzeta")
    abscissae.append(values)
    # Ensure only 1 value is repeated:
    numRepeatedValues = 0
    for value in values:
        if data.count(value)>1:
            numRepeatedValues += 1
            convergedValue = value
    if numRepeatedValues>1:
        print "Error! Multiple values of "+parametersToVary[-1]+" were repeated."
        exit(1)
    if numRepeatedValues==0:
       convergedValue = data[baseCaseIndex]
    thisScanResults = []
    for value in values:
        if value==convergedValue:
            index = baseCaseIndex
        else:
            index = data.index(value)
        thisScanResults.append(outputs[index])
    print "A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue)
    print "Values of the scanned parameter: ",abscissae[-1]
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

# Check whether Nxi was scanned:
data = Nxis
values = uniq(data)
values.sort()
if len(values)>1:
    # Yes, a scan was performed in this variable.
    parametersToVary.append("Nxi")
    abscissae.append(values)
    # Ensure only 1 value is repeated:
    numRepeatedValues = 0
    for value in values:
        if data.count(value)>1:
            numRepeatedValues += 1
            convergedValue = value
    if numRepeatedValues>1:
        print "Error! Multiple values of "+parametersToVary[-1]+" were repeated."
        exit(1)
    if numRepeatedValues==0:
       convergedValue = data[baseCaseIndex]
    thisScanResults = []
    for value in values:
        if value==convergedValue:
            index = baseCaseIndex
        else:
            index = data.index(value)
        thisScanResults.append(outputs[index])
    print "A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue)
    print "Values of the scanned parameter: ",abscissae[-1]
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

# Check whether NL was scanned:
data = NLs
values = uniq(data)
values.sort()
if len(values)>1:
    # Yes, a scan was performed in this variable.
    parametersToVary.append("NL")
    abscissae.append(values)
    # Ensure only 1 value is repeated:
    numRepeatedValues = 0
    for value in values:
        if data.count(value)>1:
            numRepeatedValues += 1
            convergedValue = value
    if numRepeatedValues>1:
        print "Error! Multiple values of "+parametersToVary[-1]+" were repeated."
        exit(1)
    if numRepeatedValues==0:
       convergedValue = data[baseCaseIndex]
    thisScanResults = []
    for value in values:
        if value==convergedValue:
            index = baseCaseIndex
        else:
            index = data.index(value)
        thisScanResults.append(outputs[index])
    print "A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue)
    print "Values of the scanned parameter: ",abscissae[-1]
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

# Check whether Nx was scanned:
data = Nxs
values = uniq(data)
values.sort()
if len(values)>1:
    # Yes, a scan was performed in this variable.
    parametersToVary.append("Nx")
    abscissae.append(values)
    # Ensure only 1 value is repeated:
    numRepeatedValues = 0
    for value in values:
        if data.count(value)>1:
            numRepeatedValues += 1
            convergedValue = value
    if numRepeatedValues>1:
        print "Error! Multiple values of "+parametersToVary[-1]+" were repeated."
        exit(1)
    if numRepeatedValues==0:
       convergedValue = data[baseCaseIndex]
    thisScanResults = [] 
    for value in values:
       if value==convergedValue:
          index = baseCaseIndex
       else:
          index = data.index(value)
       thisScanResults.append(outputs[index])
    print "A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue)
    print "Values of the scanned parameter: ",abscissae[-1]
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

# Check whether NxPotentialsPerVth was scanned:
data = NxPotentialsPerVths
values = uniq(data)
values.sort()
if len(values)>1:
    # Yes, a scan was performed in this variable.
    parametersToVary.append("NxPotentialsPerVth")
    abscissae.append(values)
    # Ensure only 1 value is repeated:
    numRepeatedValues = 0
    for value in values:
        if data.count(value)>1:
            numRepeatedValues += 1
            convergedValue = value
    if numRepeatedValues>1:
        print "Error! Multiple values of "+parametersToVary[-1]+" were repeated."
        exit(1)
    if numRepeatedValues==0:
       convergedValue = data[baseCaseIndex]
    thisScanResults = []
    for value in values:
        if value==convergedValue:
            index = baseCaseIndex
        else:
            index = data.index(value)
        thisScanResults.append(outputs[index])
    print "A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue)
    print "Values of the scanned parameter: ",abscissae[-1]
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

# Check whether xMax was scanned:
data = xMaxs
values = uniq(data)
values.sort()
if len(values)>1:
    # Yes, a scan was performed in this variable.
    parametersToVary.append("xMax")
    abscissae.append(values)
    # Ensure only 1 value is repeated:
    numRepeatedValues = 0
    for value in values:
        if data.count(value)>1:
            numRepeatedValues += 1
            convergedValue = value
    if numRepeatedValues>1:
        print "Error! Multiple values of "+parametersToVary[-1]+" were repeated."
        exit(1)
    if numRepeatedValues==0:
       convergedValue = data[baseCaseIndex]
    thisScanResults = []
    for value in values:
        if value==convergedValue:
            index = baseCaseIndex
        else:
            index = data.index(value)
        thisScanResults.append(outputs[index])
    print "A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue)
    print "Values of the scanned parameter: ",abscissae[-1]
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

# Check whether solverTolerances was scanned:
data = [-math.log10(x) for x in solverTolerances]
values = uniq(data)
values.sort()
if len(values)>1:
    # Yes, a scan was performed in this variable.
    parametersToVary.append("-log10(solverTolerance)")
    abscissae.append(values)
    # Ensure only 1 value is repeated:
    numRepeatedValues = 0
    for value in values:
        if data.count(value)>1:
            numRepeatedValues += 1
            convergedValue = value
    if numRepeatedValues>1:
        print "Error! Multiple values of "+parametersToVary[-1]+" were repeated."
        exit(1)
    if numRepeatedValues==0:
       convergedValue = data[baseCaseIndex]
    thisScanResults = []
    for value in values:
        if value==convergedValue:
            index = baseCaseIndex
        else:
            index = data.index(value)
        thisScanResults.append(outputs[index])
    print "A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue)
    print "Values of the scanned parameter: ",abscissae[-1]
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))


# ***************************************************
# Now make the plot
# ***************************************************

fig = plt.figure()
fig.patch.set_facecolor('white')
#ax = plt.axes([0,0,1,1],axisbg='w')

numParameters = len(parametersToVary)
numQuantities = numRows
numCols = int(numParameters)
numRows = int(numRows)

maxs = []
mins = []
for iQuantity in range(numQuantities):
    thisMax = quantities[0][:,iQuantity].max()
    thisMin = quantities[0][:,iQuantity].min()
    if numRuns>1:
        for iParameter in range(numParameters-1):
            thisMax = max(thisMax,quantities[iParameter+1][:,iQuantity].max())
            thisMin = min(thisMin,quantities[iParameter+1][:,iQuantity].min())
    if thisMin >= thisMax:
        thisMin -= 0.5
        thisMax += 0.5
    maxs.append(thisMax)
    mins.append(thisMin)

for iQuantity in range(numQuantities):
    for iParameter in range(numParameters):
        plt.subplot(numRows,numCols,iParameter+1+(iQuantity)*numParameters)
        plt.plot(abscissae[iParameter],quantities[iParameter][:,iQuantity],'.-')
        if iQuantity==numQuantities-1:
            plt.xlabel(parametersToVary[iParameter])
        else:
            if hideRepeatedScales:
                plt.gca().axes.xaxis.set_ticklabels([])
        if iParameter==0:
            plt.ylabel(yAxisLabels[iQuantity])
        else:
            if hideRepeatedScales:
                plt.gca().axes.yaxis.set_ticklabels([])
        plt.ylim(mins[iQuantity],maxs[iQuantity])
        plt.plot([convergeds[iParameter],convergeds[iParameter]], [mins[iQuantity],maxs[iQuantity]], 'r')


titleString = "Plot generated by "+ inspect.getfile(inspect.currentframe()) + "\nRun in "+os.getcwd()
ax = fig.add_axes([0,0,1,1], frameon=False)
ax.text(0.5,0.99,titleString,horizontalalignment='center',verticalalignment='top')

plt.show()

