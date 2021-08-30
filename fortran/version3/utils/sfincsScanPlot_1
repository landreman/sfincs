#!/usr/bin/env python

# This python script plots the output of a SFINCS convergence scan.

from builtins import input

outputFilename = "sfincsOutput.h5"

hideRepeatedScales = True
#hideRepeatedScales = False

import matplotlib
#import matplotlib.pyplot as plt
import h5py
import numpy
import inspect, math, os
import pickle
import sys
#from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib.offsetbox import AnchoredText

import matplotlib.ticker as ticker

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsScanPlot_1"
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

numRuns = 0
Nthetas = []
Nzetas = []
Nxis = []
NLs = []
Nxs = []
NxPotentialsPerVths = []
xMaxs = []
solverTolerances = []
Zs = []

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

    if math.isnan(dummy[0,0]):
       print ("Run in directory "+directory+" has NaNs, so skipping it.")
       continue

    integerToRepresentTrue = (f["integerToRepresentTrue"][()])
    includePhi1_new = (f["includePhi1"][()] == integerToRepresentTrue)

    if includePhi1_new:
       try:
          if f["didNonlinearCalculationConverge"][()] != integerToRepresentTrue:
             print ("The nonlinear iteration in directory "+directory+" does not appear to have converged.")
             continue
       except:
          print ("WARNING! Cannot read didNonlinearCalculationConverge in directory "+directory+" even though this is a nonlinear calculation.")
          print ("You should manually check if the nonlinear iteration has converged!")

    print ("Processing directory "+directory)

    # The expression [()] converts from an h5py dataset to a numpy ndarray:
    #integerToRepresentTrue = (f["integerToRepresentTrue"][()])
    Nthetas.append(f["Ntheta"][()])
    Nzetas.append(f["Nzeta"][()])
    Nxis.append(f["Nxi"][()])
    NLs.append(f["NL"][()])
    Nxs.append(f["Nx"][()])
    NxPotentialsPerVths.append(f["NxPotentialsPerVth"][()])
    xMaxs.append(f["xMax"][()])
    solverTolerances.append(f["solverTolerance"][()])
    Zs.append(f["Zs"][()])
    
    RHSMode_new = f["RHSMode"][()]
    Nspecies_new = f["Nspecies"][()]
    #includePhi1_new = (f["includePhi1"][()] == integerToRepresentTrue)
    if numRuns == 0:
        RHSMode = RHSMode_new
        Nspecies = Nspecies_new
        includePhi1 = includePhi1_new
    else:
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

    FSABFlow = f["FSABFlow"][()]
    FSABjHat = f["FSABjHat"][()]
    particleFlux_vm = f["particleFlux_vm_psiHat"][()]
    heatFlux_vm = f["heatFlux_vm_psiHat"][()]
    elapsedTime = f["elapsed time (s)"][()]
    try:
       sources = f["sources"][()]
    except:
       sources = numpy.array([0,0])
    radCurr_vm = 0
    radCurr_vd = 0

    if RHSMode>1:
       transportMatrix = f["transportMatrix"][()]
    if includePhi1:
        particleFlux_vd = f["particleFlux_vd_psiHat"][()]
        heatFlux_vd = f["heatFlux_vd_psiHat"][()]
        Phi1atOrigin = (f["Phi1Hat"][()])[0,0,:]

    if RHSMode > 1:
       results = []
       for icol in range(transportMatrix.shape[1]):
          for irow in range(transportMatrix.shape[0]):
             results.append(transportMatrix[irow,icol])
       outputs.append(results)
    else:
       if Nspecies==1:
          if includePhi1:
             outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],particleFlux_vd[0,0],heatFlux_vm[0,0],heatFlux_vd[0,0], Phi1atOrigin[0]])
          else:
             if sources.shape[0]>1:
                try:
                   outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],heatFlux_vm[0,0],sources[0,0,0],sources[1,0,0]])
                except:
                   outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],heatFlux_vm[0,0],sources[0],sources[1]])
             else:
                outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],heatFlux_vm[0,0],sources[0,0,0]])
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
                radCurr_vm += particleFlux_vm[ispecies,iteration]*Zs[-1][ispecies]
                radCurr_vd += particleFlux_vd[ispecies,iteration]*Zs[-1][ispecies]
             else:
                results.append(FSABFlow[ispecies,iteration])
                results.append(particleFlux_vm[ispecies,iteration])
                results.append(heatFlux_vm[ispecies,iteration])
                radCurr_vm += particleFlux_vm[ispecies,iteration]*Zs[-1][ispecies]
          results.append(FSABjHat[iteration])
          results.append(radCurr_vm)

          if includePhi1:
             results.append(radCurr_vd) 
             results.append(Phi1atOrigin[iteration])

          outputs.append(results)

    if directory=="baseCase":
        baseCaseIndex=numRuns
    numRuns += 1

    print ("Successfully read run in directory "+directory)
 

if baseCaseIndex<0:
    print ("Error! No baseCase run available.")
    exit(1)

if RHSMode>1:
   yAxisLabels=[]
   for icol in range(transportMatrix.shape[1]):
      for irow in range(transportMatrix.shape[0]):
         yAxisLabels.append("L"+str(irow+1)+str(icol+1))
else:
   if Nspecies==1:
      if includePhi1:
         yAxisLabels=["FSABFlow","particleFlux\nvm_psiHat","particleFlux\nvd_psiHat","heatFlux\nvm_psiHat","heatFlux\nvd_psiHat", "Phi1Hat\n[0,0]"]
      else:
         if sources.shape[0]>1:
            yAxisLabels=["FSABFlow","particleFlux\nvm_psiHat","heatFlux\nvm_psiHat","particle source",'heat source']
         else:
            yAxisLabels=["FSABFlow","particleFlux\nvm_psiHat","heatFlux\nvm_psiHat","source"]
   else:
      if includePhi1:
         yAxisLabels=[]
         for i in range(1,Nspecies+1):
            yAxisLabels.append("FSABFlow\n(species "+str(i)+")")
            yAxisLabels.append("particleFlux\nvm_psiHat\n(species "+str(i)+")")
            yAxisLabels.append("particleFlux\nvd_psiHat\n(species "+str(i)+")")
            yAxisLabels.append("heatFlux\nvm_psiHat\n(species "+str(i)+")")
            yAxisLabels.append("heatFlux\nvd_psiHat\n(species "+str(i)+")")
         yAxisLabels.append("FSABjHat")
         yAxisLabels.append("J_net\nvm_psiHat")
         yAxisLabels.append("J_net\nvd_psiHat")
         yAxisLabels.append("Phi1Hat\n[0,0]")
        
      else:
         yAxisLabels=[]
         for i in range(1,Nspecies+1):
            yAxisLabels.append("FSABFlow\n(species "+str(i)+")")
            yAxisLabels.append("particleFlux\nvm_psiHat\n(species "+str(i)+")")
            yAxisLabels.append("heatFlux\nvm_psiHat\n(species "+str(i)+")")
         yAxisLabels.append("FSABjHat")
         yAxisLabels.append("Radial Current")



##SINCE Phi1Hat CAN BE LARGE WE ONLY PLOT IT FOR THE MIN, MAX AND BASE CASE OF EACH SCAN
PlotPhi1Bool = False
if includePhi1: 
   #YesOrNo = raw_input("Do you want to make a plot of Phi1? [n (default) / y]: ")
   YesOrNo = input("Do you want to make a plot of Phi1? [n (default) / y]: ")
   if YesOrNo == "y" or YesOrNo == "Y" or YesOrNo.lower() == "yes":
      PlotPhi1Bool = True

if PlotPhi1Bool and includePhi1:
    if not os.path.isfile("baseCase/" + outputFilename):
        print ("Directory " + "baseCase" + " does not have a " + outputFilename + " file (yet). Will not plot Phi1!")
        PlotPhi1Bool = False
    try:
        f = h5py.File("baseCase/" + outputFilename, 'r')
    except IOError:
        print ("Unable to open " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
        PlotPhi1Bool = False

    try:
        Phi1BaseCase = (f["Phi1Hat"][()])[:,:,-1]
        thetaBaseCase = f["theta"][()]
        zetaBaseCase = f["zeta"][()]
    except:
        print ("Unable to read Phi1Hat from " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
        PlotPhi1Bool = False
    Phi1_Array = []
    Phi1_Theta_Array = []
    Phi1_Zeta_Array = []
    Phi1_Labels_Array = []
#########################
            
numRows = len(yAxisLabels)
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
        print ("Error! Multiple values of "+parametersToVary[-1]+" were repeated.")
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
    print ("A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue))
    print ("Values of the scanned parameter: ",abscissae[-1])
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

    if PlotPhi1Bool and includePhi1:
       if (abscissae[-1])[0] == convergedValue:
          minScanValue = "baseCase"
       else :
          minScanValue = "Ntheta" + str((abscissae[-1])[0])

       if (abscissae[-1])[-1] == convergedValue:
          maxScanValue = "baseCase"
       else :
          maxScanValue = "Ntheta" + str((abscissae[-1])[-1])

       Phi1_Labels_Array.append(minScanValue)
       Phi1_Labels_Array.append(maxScanValue)
       try:
          fmin = h5py.File(minScanValue + "/" + outputFilename, 'r')
          fmax = h5py.File(maxScanValue + "/" + outputFilename, 'r')
       except IOError:
          print ("Unable to open " + minScanValue + "/" + outputFilename + " or " + maxScanValue + "/" + outputFilename + ". Will not plot Phi1!")
          PlotPhi1Bool = False
          
       try:
          Phi1_Array.append( (fmin["Phi1Hat"][()])[:,:,-1] )
          Phi1_Array.append( (fmax["Phi1Hat"][()])[:,:,-1] )
          Phi1_Theta_Array.append( fmin["theta"][()] )
          Phi1_Theta_Array.append( fmax["theta"][()] )
          Phi1_Zeta_Array.append( fmin["zeta"][()] )
          Phi1_Zeta_Array.append( fmax["zeta"][()] )
       except:
          print ("Unable to read Phi1Hat from " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
          PlotPhi1Bool = False
    ########################

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
        print ("Error! Multiple values of "+parametersToVary[-1]+" were repeated.")
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
    print ("A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue))
    print ("Values of the scanned parameter: ",abscissae[-1])
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

    if PlotPhi1Bool and includePhi1:
       if (abscissae[-1])[0] == convergedValue:
          minScanValue = "baseCase"
       else :
          minScanValue = "Nzeta" + str((abscissae[-1])[0])

       if (abscissae[-1])[-1] == convergedValue:
          maxScanValue = "baseCase"
       else :
          maxScanValue = "Nzeta" + str((abscissae[-1])[-1])

       Phi1_Labels_Array.append(minScanValue)
       Phi1_Labels_Array.append(maxScanValue)
       try:
          fmin = h5py.File(minScanValue + "/" + outputFilename, 'r')
          fmax = h5py.File(maxScanValue + "/" + outputFilename, 'r')
       except IOError:
          print ("Unable to open " + minScanValue + "/" + outputFilename + " or " + maxScanValue + "/" + outputFilename + ". Will not plot Phi1!")
          PlotPhi1Bool = False
          
       try:
          Phi1_Array.append( (fmin["Phi1Hat"][()])[:,:,-1] )
          Phi1_Array.append( (fmax["Phi1Hat"][()])[:,:,-1] )
          Phi1_Theta_Array.append( fmin["theta"][()] )
          Phi1_Theta_Array.append( fmax["theta"][()] )
          Phi1_Zeta_Array.append( fmin["zeta"][()] )
          Phi1_Zeta_Array.append( fmax["zeta"][()] )
       except:
          print ("Unable to read Phi1Hat from " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
          PlotPhi1Bool = False
    ########################

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
        print ("Error! Multiple values of "+parametersToVary[-1]+" were repeated.")
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
    print ("A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue))
    print ("Values of the scanned parameter: ",abscissae[-1])
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

    if PlotPhi1Bool and includePhi1:
       if (abscissae[-1])[0] == convergedValue:
          minScanValue = "baseCase"
       else :
          minScanValue = "Nxi" + str((abscissae[-1])[0])

       if (abscissae[-1])[-1] == convergedValue:
          maxScanValue = "baseCase"
       else :
          maxScanValue = "Nxi" + str((abscissae[-1])[-1])

       Phi1_Labels_Array.append(minScanValue)
       Phi1_Labels_Array.append(maxScanValue)
       try:
          fmin = h5py.File(minScanValue + "/" + outputFilename, 'r')
          fmax = h5py.File(maxScanValue + "/" + outputFilename, 'r')
       except IOError:
          print ("Unable to open " + minScanValue + "/" + outputFilename + " or " + maxScanValue + "/" + outputFilename + ". Will not plot Phi1!")
          PlotPhi1Bool = False
          
       try:
          Phi1_Array.append( (fmin["Phi1Hat"][()])[:,:,-1] )
          Phi1_Array.append( (fmax["Phi1Hat"][()])[:,:,-1] )
          Phi1_Theta_Array.append( fmin["theta"][()] )
          Phi1_Theta_Array.append( fmax["theta"][()] )
          Phi1_Zeta_Array.append( fmin["zeta"][()] )
          Phi1_Zeta_Array.append( fmax["zeta"][()] )
       except:
          print ("Unable to read Phi1Hat from " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
          PlotPhi1Bool = False
    ########################

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
        print ("Error! Multiple values of "+parametersToVary[-1]+" were repeated.")
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
    print ("A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue))
    print ("Values of the scanned parameter: ",abscissae[-1])
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

    if PlotPhi1Bool and includePhi1:
       if (abscissae[-1])[0] == convergedValue:
          minScanValue = "baseCase"
       else :
          minScanValue = "NL" + str((abscissae[-1])[0])

       if (abscissae[-1])[-1] == convergedValue:
          maxScanValue = "baseCase"
       else :
          maxScanValue = "NL" + str((abscissae[-1])[-1])

       Phi1_Labels_Array.append(minScanValue)
       Phi1_Labels_Array.append(maxScanValue)
       try:
          fmin = h5py.File(minScanValue + "/" + outputFilename, 'r')
          fmax = h5py.File(maxScanValue + "/" + outputFilename, 'r')
       except IOError:
          print ("Unable to open " + minScanValue + "/" + outputFilename + " or " + maxScanValue + "/" + outputFilename + ". Will not plot Phi1!")
          PlotPhi1Bool = False
          
       try:
          Phi1_Array.append( (fmin["Phi1Hat"][()])[:,:,-1] )
          Phi1_Array.append( (fmax["Phi1Hat"][()])[:,:,-1] )
          Phi1_Theta_Array.append( fmin["theta"][()] )
          Phi1_Theta_Array.append( fmax["theta"][()] )
          Phi1_Zeta_Array.append( fmin["zeta"][()] )
          Phi1_Zeta_Array.append( fmax["zeta"][()] )
       except:
          print ("Unable to read Phi1Hat from " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
          PlotPhi1Bool = False
    ########################

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
        print ("Error! Multiple values of "+parametersToVary[-1]+" were repeated.")
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
    print ("A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue))
    print ("Values of the scanned parameter: ",abscissae[-1])
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))


    if PlotPhi1Bool and includePhi1:
       if (abscissae[-1])[0] == convergedValue:
          minScanValue = "baseCase"
       else :
          minScanValue = "Nx" + str((abscissae[-1])[0])

       if (abscissae[-1])[-1] == convergedValue:
          maxScanValue = "baseCase"
       else :
          maxScanValue = "Nx" + str((abscissae[-1])[-1])

       Phi1_Labels_Array.append(minScanValue)
       Phi1_Labels_Array.append(maxScanValue)
       try:
          fmin = h5py.File(minScanValue + "/" + outputFilename, 'r')
          fmax = h5py.File(maxScanValue + "/" + outputFilename, 'r')
       except IOError:
          print ("Unable to open " + minScanValue + "/" + outputFilename + " or " + maxScanValue + "/" + outputFilename + ". Will not plot Phi1!")
          PlotPhi1Bool = False
          
       try:
          Phi1_Array.append( (fmin["Phi1Hat"][()])[:,:,-1] )
          Phi1_Array.append( (fmax["Phi1Hat"][()])[:,:,-1] )
          Phi1_Theta_Array.append( fmin["theta"][()] )
          Phi1_Theta_Array.append( fmax["theta"][()] )
          Phi1_Zeta_Array.append( fmin["zeta"][()] )
          Phi1_Zeta_Array.append( fmax["zeta"][()] )
       except:
          print ("Unable to read Phi1Hat from " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
          PlotPhi1Bool = False
    ########################

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
        print ("Error! Multiple values of "+parametersToVary[-1]+" were repeated.")
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
    print ("A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue))
    print ("Values of the scanned parameter: ",abscissae[-1])
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

    if PlotPhi1Bool and includePhi1:
       if (abscissae[-1])[0] == convergedValue:
          minScanValue = "baseCase"
       else :
          minScanValue = "NxPotentialsPerVth" + "{:.4g}".format((abscissae[-1])[0])

       if (abscissae[-1])[-1] == convergedValue:
          maxScanValue = "baseCase"
       else :
          maxScanValue = "NxPotentialsPerVth" + "{:.4g}".format((abscissae[-1])[-1])

       Phi1_Labels_Array.append(minScanValue)
       Phi1_Labels_Array.append(maxScanValue)
       try:
          fmin = h5py.File(minScanValue + "/" + outputFilename, 'r')
          fmax = h5py.File(maxScanValue + "/" + outputFilename, 'r')
       except IOError:
          print ("Unable to open " + minScanValue + "/" + outputFilename + " or " + maxScanValue + "/" + outputFilename + ". Will not plot Phi1!")
          PlotPhi1Bool = False
          
       try:
          Phi1_Array.append( (fmin["Phi1Hat"][()])[:,:,-1] )
          Phi1_Array.append( (fmax["Phi1Hat"][()])[:,:,-1] )
          Phi1_Theta_Array.append( fmin["theta"][()] )
          Phi1_Theta_Array.append( fmax["theta"][()] )
          Phi1_Zeta_Array.append( fmin["zeta"][()] )
          Phi1_Zeta_Array.append( fmax["zeta"][()] )
       except:
          print ("Unable to read Phi1Hat from " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
          PlotPhi1Bool = False
    ########################

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
        print ("Error! Multiple values of "+parametersToVary[-1]+" were repeated.")
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
    print ("A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue))
    print ("Values of the scanned parameter: ",abscissae[-1])
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

    if PlotPhi1Bool and includePhi1:
       if (abscissae[-1])[0] == convergedValue:
          minScanValue = "baseCase"
       else :
          minScanValue = "xMax" + "{:.4g}".format((abscissae[-1])[0])

       if (abscissae[-1])[-1] == convergedValue:
          maxScanValue = "baseCase"
       else :
          maxScanValue = "xMax" + "{:.4g}".format((abscissae[-1])[-1])

       Phi1_Labels_Array.append(minScanValue)
       Phi1_Labels_Array.append(maxScanValue)
       try:
          fmin = h5py.File(minScanValue + "/" + outputFilename, 'r')
          fmax = h5py.File(maxScanValue + "/" + outputFilename, 'r')
       except IOError:
          print ("Unable to open " + minScanValue + "/" + outputFilename + " or " + maxScanValue + "/" + outputFilename + ". Will not plot Phi1!")
          PlotPhi1Bool = False
          
       try:
          Phi1_Array.append( (fmin["Phi1Hat"][()])[:,:,-1] )
          Phi1_Array.append( (fmax["Phi1Hat"][()])[:,:,-1] )
          Phi1_Theta_Array.append( fmin["theta"][()] )
          Phi1_Theta_Array.append( fmax["theta"][()] )
          Phi1_Zeta_Array.append( fmin["zeta"][()] )
          Phi1_Zeta_Array.append( fmax["zeta"][()] )
       except:
          print ("Unable to read Phi1Hat from " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
          PlotPhi1Bool = False
    ########################

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
        print ("Error! Multiple values of "+parametersToVary[-1]+" were repeated.")
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
    print ("A scan was performed in "+parametersToVary[-1]+", and the repeated value was "+str(convergedValue))
    print ("Values of the scanned parameter: ",abscissae[-1])
    convergeds.append(convergedValue)
    quantities.append(numpy.array(thisScanResults))

    if PlotPhi1Bool and includePhi1:
       if (abscissae[-1])[0] == convergedValue:
          minScanValue = "baseCase"
       else :
          minScanValue = "solverTolerance" + "{:.4g}".format(math.pow(10, -(abscissae[-1])[0]))

       if (abscissae[-1])[-1] == convergedValue:
          maxScanValue = "baseCase"
       else :
          maxScanValue = "solverTolerance" + "{:.4g}".format(math.pow(10, -(abscissae[-1])[-1]))

       Phi1_Labels_Array.append(minScanValue)
       Phi1_Labels_Array.append(maxScanValue)
       try:
          fmin = h5py.File(minScanValue + "/" + outputFilename, 'r')
          fmax = h5py.File(maxScanValue + "/" + outputFilename, 'r')
       except IOError:
          print ("Unable to open " + minScanValue + "/" + outputFilename + " or " + maxScanValue + "/" + outputFilename + ". Will not plot Phi1!")
          PlotPhi1Bool = False
          
       try:
          Phi1_Array.append( (fmin["Phi1Hat"][()])[:,:,-1] )
          Phi1_Array.append( (fmax["Phi1Hat"][()])[:,:,-1] )
          Phi1_Theta_Array.append( fmin["theta"][()] )
          Phi1_Theta_Array.append( fmax["theta"][()] )
          Phi1_Zeta_Array.append( fmin["zeta"][()] )
          Phi1_Zeta_Array.append( fmax["zeta"][()] )
       except:
          print ("Unable to read Phi1Hat from " + "baseCase/" + outputFilename + " even though this file exists. Will not plot Phi1!")
          PlotPhi1Bool = False
    ########################


# ***************************************************
# Now make the plot
# ***************************************************

fig = plt.figure(figsize=(16,8.5))
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

    expand = (thisMax-thisMin)*0.1
    thisMax = thisMax + expand
    thisMin = thisMin - expand
    maxs.append(thisMax)
    mins.append(thisMin)

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
    for iParameter in range(numParameters):
        ax = plt.subplot(numRows,numCols,iParameter+1+(iQuantity)*numParameters)
        xdata.append(abscissae[iParameter])
        ydata.append(quantities[iParameter][:,iQuantity])
        plt.plot(xdata[-1],ydata[-1],linespec)

        # Add annotation showing the relative spread in the result:
        avg = abs(max(ydata[-1])+min(ydata[-1]))/2
        if avg>0:
           relVariation = (max(ydata[-1])-min(ydata[-1]))/avg
        else:
           relVariation = 0
        label = 'spread='+("%.1f" % (relVariation*100))+'%'
        if min(xdata[-1])<convergeds[iParameter] and max(xdata[-1]) > convergeds[iParameter]:
           ydata2 = []
           for i in range(len(xdata[-1])):
              if xdata[-1][i] >= convergeds[iParameter]:
                 ydata2.append(ydata[-1][i])
           avg = abs(max(ydata2)+min(ydata2))/2
           if avg>0:
              relVariation = (max(ydata2)-min(ydata2))/avg
           else:
              relVariation = 0
           label = label+', '+("%.1f" % (relVariation*100))+'% above base'
        at = AnchoredText(label,loc=1,frameon=False)
        #at = AnchoredText(label,loc=1,frameon=False,prop=dict(size='x-small'))
        ax.add_artist(at)

        if iQuantity==numQuantities-1:
            xlabels.append(parametersToVary[iParameter])
            plt.xlabel(xlabels[-1])
        else:
            xlabels.append('')
            if hideRepeatedScales:
                plt.gca().axes.xaxis.set_ticklabels([])
        if iParameter==0:
            ylabels.append(yAxisLabels[iQuantity])
            plt.ylabel(yAxisLabels[iQuantity])
        else:
            ylabels.append('')
            if hideRepeatedScales:
                plt.gca().axes.yaxis.set_ticklabels([])
        plt.ylim(mins[iQuantity],maxs[iQuantity])
        ymins.append(mins[iQuantity])
        ymaxs.append(maxs[iQuantity])
        xscales.append('linear')
        yscales.append('linear')
        plt.plot([convergeds[iParameter],convergeds[iParameter]], [mins[iQuantity],maxs[iQuantity]], 'r')


#titleString = "Plot generated by "+ inspect.getfile(inspect.currentframe()) + "\nRun in "+os.getcwd()
titleString = "Plot generated by "+ this_filename + "\nRun in "+os.getcwd()
ax = fig.add_axes([0,0,1,1], frameon=False)
ax.text(0.5,0.99,titleString,horizontalalignment='center',verticalalignment='top')

outputFile = open('sfincsScan.dat','wb')
scanType=1
numQuantitiesForPickle = numRows*numCols
data = {'scanType':scanType, 'numQuantities':numQuantitiesForPickle, 'numRows':numRows,'numCols':numCols,
        'xdata':xdata, 'ydata':ydata, 'xlabels':xlabels, 'ylabels':ylabels,
        'xscales':xscales, 'yscales':yscales, 'ymins':ymins, 'ymaxs':ymaxs,
        'linespec':linespec}
# pickle.dump(scanType,outputFile)
# pickle.dump(numParameters*numQuantities,outputFile)
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


# *******************************
# MAKE PLOT OF Phi1 IF APPLICABLE
# *******************************

if PlotPhi1Bool and includePhi1:

   def fmt_cbar(x, pos):
      if x == 0.0:
         return r'${}$'.format(x)
      a, b = '{:.1e}'.format(x).split('e')
      b = int(b)
      return r'${} \cdot 10^{{{}}}$'.format(a, b)

   def fmt_xy_axis(x, pos):
      return r'${}$'.format(x)

   numRowsPlot2 = 3
   numColsPlot2 = numParameters
   numContours = 100
   numLevels = 5

   fig2 = plt.figure(2)
   fig2.patch.set_facecolor('white')


#   print numpy.amax(Phi1BaseCase)
#   print numpy.amin(Phi1BaseCase)
#   print numParameters
#   print len(Phi1_Array)
#   print len(Phi1_Theta_Array)
#   print len(Phi1_Zeta_Array)
#   print len(Phi1_Labels_Array)

   for i in range(0, len(Phi1_Array), 2):

      ##Add min Plot##
      ################

      delta = (numpy.amax(Phi1_Array[i]) - numpy.amin(Phi1_Array[i])) / numLevels
      ContourLevels = numpy.arange(numpy.amin(Phi1_Array[i]), numpy.amax(Phi1_Array[i]) + delta/2.0, delta)
      
      ax = plt.subplot(numRowsPlot2, numColsPlot2, i/2 + 1)
      Phi1Plot = plt.contourf(Phi1_Zeta_Array[i], Phi1_Theta_Array[i], Phi1_Array[i].transpose(), numContours)
      Phi1Plot2 = plt.contour(Phi1Plot,levels=ContourLevels, colors='k', hold='on')
      plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
      plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')
      plt.gca().axes.xaxis.set_label_coords(0.5,-0.08)
      plt.gca().axes.yaxis.set_label_coords(-0.12,0.5)
      ax.xaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis))
      ax.xaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))
      ax.yaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis))
      ax.yaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))
      
      cbar = plt.colorbar(Phi1Plot, format=ticker.FuncFormatter(fmt_cbar))#, ticks=ContourLevels)
      cbar.ax.set_ylabel(r'$\hat{\Phi}_1$', rotation=0, labelpad=10)
      
      plt.clabel(Phi1Plot2, fmt=ticker.FuncFormatter(fmt_cbar), colors='k', fontsize=10, inline=False)
      ax.text(.5,.9,Phi1_Labels_Array[i],
               horizontalalignment='center',
               transform=ax.transAxes, fontsize=12)

      ##Add base case plot##
      ######################
      if (i  == 2*int(math.floor(numColsPlot2/2))):
         delta = (numpy.amax(Phi1BaseCase) - numpy.amin(Phi1BaseCase)) / numLevels
         ContourLevels = numpy.arange(numpy.amin(Phi1BaseCase), numpy.amax(Phi1BaseCase) + delta/2.0, delta)
         
         ax = plt.subplot(numRowsPlot2,numColsPlot2, i/2 + 1 + numColsPlot2)
         Phi1Plot = plt.contourf(zetaBaseCase, thetaBaseCase, Phi1BaseCase.transpose(), numContours)
         Phi1Plot2 = plt.contour(Phi1Plot,levels=ContourLevels, colors='k', hold='on')
         plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
         plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')
         plt.gca().axes.xaxis.set_label_coords(0.5,-0.08)
         plt.gca().axes.yaxis.set_label_coords(-0.12,0.5)
         #ax.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
         #ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))         
         #ax.ticklabel_format(fmt = ticker.FuncFormatter(fmt_xy_axis))
         ax.xaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis))
         ax.xaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))
         ax.yaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis))
         ax.yaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))

#         plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
         
         cbar = plt.colorbar(Phi1Plot, format=ticker.FuncFormatter(fmt_cbar))#, ticks=ContourLevels)
         cbar.ax.set_ylabel(r'$\hat{\Phi}_1$', rotation=0, labelpad=10)
         
#         fmtPlot = ticker.LogFormatterMathtext()
#         fmtPlot.create_dummy_axis()
      #      fmtPlot = "%.4g"
         
         plt.clabel(Phi1Plot2, fmt=ticker.FuncFormatter(fmt_cbar), colors='k', fontsize=10, inline=False)
         ax.text(.5,.9,'baseCase',
                  horizontalalignment='center',
                  transform=ax.transAxes, fontsize=12)
         
      ##Add max Plot##
      ################

      delta = (numpy.amax(Phi1_Array[i+1]) - numpy.amin(Phi1_Array[i+1])) / numLevels
      ContourLevels = numpy.arange(numpy.amin(Phi1_Array[i+1]), numpy.amax(Phi1_Array[i+1]) + delta/2.0, delta)
      
      ax = plt.subplot(numRowsPlot2, numColsPlot2, i/2 + 1 + 2*numColsPlot2)
      Phi1Plot = plt.contourf(Phi1_Zeta_Array[i+1], Phi1_Theta_Array[i+1], Phi1_Array[i+1].transpose(), numContours)
      Phi1Plot2 = plt.contour(Phi1Plot,levels=ContourLevels, colors='k', hold='on')
      plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
      plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')
      plt.gca().axes.xaxis.set_label_coords(0.5,-0.08)
      plt.gca().axes.yaxis.set_label_coords(-0.12,0.5)
      ax.xaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis))
      ax.xaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))
      ax.yaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis))
      ax.yaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))
      
      cbar = plt.colorbar(Phi1Plot, format=ticker.FuncFormatter(fmt_cbar))#, ticks=ContourLevels)
      cbar.ax.set_ylabel(r'$\hat{\Phi}_1$', rotation=0, labelpad=10)
      
      plt.clabel(Phi1Plot2, fmt=ticker.FuncFormatter(fmt_cbar), colors='k', fontsize=10, inline=False)
      ax.text(.5,.9,Phi1_Labels_Array[i+1],
               horizontalalignment='center',
               transform=ax.transAxes, fontsize=12)
      
########################   


# If this script was called with any command-line arguments, then do not actually display the plot:
if len(sys.argv) == 1:
   plt.show()

if makePDF:
   fig = plt.figure(1)
   ##plt.savefig('sfincsScanPlot_1.pdf', orientation = 'landscape', papertype='letter')
   if len(sys.argv)>2 : #Use the substituted name as file name                                                                                      
      print ("Writing plot to " + os.getcwd() + "/" + sys.argv[2] + ".pdf.")
      plt.savefig(sys.argv[2] + ".pdf", orientation = 'landscape', papertype='letter')
   else: #Use script name as file name                                                                                                              
      print ("Writing plot to " + os.getcwd() + "/" + os.path.basename(__file__) + ".pdf.")
      plt.savefig(os.path.basename(__file__) + ".pdf", orientation = 'landscape', papertype='letter')
   #######################
 
   #plt.savefig('sfincsScanPlot_1.pdf', orientation = 'landscape', papertype='letter',bbox_inches='tight',pad_inches=0)

   if PlotPhi1Bool and includePhi1:
      fig2 = plt.figure(2)
      fig2.set_size_inches(numColsPlot2*3, numRowsPlot2*2.5)
      if len(sys.argv)>2 : #Use the substituted name as file name                                                                                      
         print ("Writing plot to " + os.getcwd() + "/" + sys.argv[2] + "_Phi1" + ".pdf.")
         plt.savefig(sys.argv[2] + "_Phi1" + ".pdf", orientation = 'landscape', papertype='letter')
      else: #Use script name as file name                                                                                                              
         print ("Writing plot to " + os.getcwd() + "/" + os.path.basename(__file__) + "_Phi1" + ".pdf.")
         plt.savefig(os.path.basename(__file__) + "_Phi1" + ".pdf", orientation = 'landscape', papertype='letter')
