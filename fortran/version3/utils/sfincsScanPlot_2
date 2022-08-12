#!/usr/bin/env python

# This python script plots the output of a SFINCS E_r scan.

import time
start_time = time.time()

inputFilename = "input.namelist"
outputFilename = "sfincsOutput.h5"

import matplotlib
#import matplotlib.pyplot as plt
import h5py
import numpy
import inspect, math, os
import sys
from timeit import default_timer as timer

from scipy.interpolate import interp1d
try:
    # Older versions of scipy do not have this interpolation routine
    from scipy.interpolate import PchipInterpolator
    whichInterpolator = 0
except:
    whichInterpolator = 1

from scipy.optimize import brentq
import pickle
from sys import argv

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsScanPlot_2"
print ("This is "+ this_filename)
print ("Time for imports: ", time.time() - start_time)

# If this script was called with any command-line arguments, then do not actually display the plot:
display_plot = (len(argv) == 1)

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

# Load some other required subroutines:
startTime = timer()
#execfile(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/sfincsScan_common")
exec(open(os.path.dirname(os.path.abspath(__file__))+"/sfincsScan_common").read())
endTime = timer()
print ("Time for exec in sfincsScanPlot_2:",endTime-startTime)

numRuns = 0
Ers = []
outputs = []
radialCurrents = []

# In case this script is being called within a radial scan, determine the radius. To do so,
# first determine which radial variable was used by looking at input.namelist in the main directory:
with open(inputFilename, 'r') as f:
    inputFile = f.readlines()

print ("Time after input file is opened: ", time.time() - start_time)

#inputRadialCoordinate = readVariable("inputRadialCoordinate","int")
inputRadialCoordinate = readVariable("inputRadialCoordinate","int", False)

if inputRadialCoordinate == None:
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
# Finally get the radius requested:
##radius_wish = readVariable(radiusName+"_wish","float")
radius_wish = readVariable(radiusName+"_wish","float", False)                                                                                                                           
if radius_wish == None:
    radius_wish = readDefault(radiusName+"_wish","float")

# Get a list of the subdirectories:                                                                        
print ("Time before listdir: ", time.time() - start_time)
directories = filter(os.path.isdir, os.listdir("."))
print ("Time after listdir: ", time.time() - start_time)

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

    ##check for nonlinear convergence##
    Check_integerToRepresentTrue = (f["integerToRepresentTrue"][()])
    Check_includePhi1 = (f["includePhi1"][()] == Check_integerToRepresentTrue)
    if Check_includePhi1 :
        try:
            if (f["didNonlinearCalculationConverge"][()] != Check_integerToRepresentTrue):
                print ("The nonlinear solver in directory "+directory+" does not appear to have converged.")
                continue
        except KeyError:
            print ("The nonlinear solver in directory "+directory+" does not appear to have finished.")
            continue
    ##################################################################

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

    if inputRadialCoordinateForGradients==0:
       varName = "dPhiHatdpsiHat"
    elif inputRadialCoordinateForGradients==1:
       varName = "dPhiHatdpsiN"
    elif inputRadialCoordinateForGradients==2:
       varName = "dPhiHatdrHat"
    elif inputRadialCoordinateForGradients==3:
       varName = "dPhiHatdrN"
    elif inputRadialCoordinateForGradients==4:
       varName = "Er"
    else:
       print ("Error! Invalid inputRadialCoordinateForGradients.")
       exit(1)

    if RHSMode != 1 and RHSMode != 2 and RHSMode != 3:
        print ("Error! sfincsScanPlot is not yet set up for RHSMode = "+str(RHSMode))
        exit(1)

    doAmbipolaritySolve = (RHSMode==1)

    radius_actual = f[radiusName][()] # Thus, radius_actual corresponds to the last file read.
    Ers.append(f[varName][()])
    FSABFlow = f["FSABFlow"][()]
    FSABjHatOverRootFSAB2 = f["FSABjHatOverRootFSAB2"][()]
    FSABjHat = f["FSABjHat"][()]
    particleFlux_vm = f["particleFlux_vm_rHat"][()]
    heatFlux_vm = f["heatFlux_vm_rHat"][()]
    elapsedTime = f["elapsed time (s)"][()]
    try:
        sources = f["sources"][()]
    except:
        sources = numpy.array([0,0])
    nHats = f["nHats"][()]
    THats = f["THats"][()]
    atLeastOneDirectorySucceeded = True
    if RHSMode>1:
       transportMatrix = f["transportMatrix"][()]
    if includePhi1:
        particleFlux_vd = f["particleFlux_vd_rHat"][()]
        #heatFlux_vd = f["heatFlux_vd_rHat"][()]
        heatFlux_withoutPhi1 = f["heatFlux_withoutPhi1_rHat"][()]

    if RHSMode > 1:
       results = []
       for icol in range(transportMatrix.shape[1]):
          for irow in range(transportMatrix.shape[0]):
             results.append(transportMatrix[irow,icol])
       outputs.append(results)
       radialCurrents.append(0)
    else:
       # RHSMode = 1
       if Nspecies==1:
          if includePhi1:
             radialCurrents.append(Zs[0]*particleFlux_vd[0,-1])
             #outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],particleFlux_vd[0,0],heatFlux_vm[0,0],heatFlux_vd[0,0]])
             outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],particleFlux_vd[0,0],heatFlux_vm[0,0],heatFlux_withoutPhi1[0,0]])
          else:
             radialCurrents.append(Zs[0]*particleFlux_vm[0,-1])
             try:
                 outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],heatFlux_vm[0,0],sources[0,0,0],sources[1,0,0]])
             except:
                 outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],heatFlux_vm[0,0],sources[0],sources[1]])
             #outputs.append([FSABFlow[0,0],particleFlux_vm[0,0],heatFlux_vm[0,0],elapsedTime])
       else:
          results = []
          iteration = -1
          for ispecies in range(Nspecies):
             if includePhi1:
                results.append(FSABFlow[ispecies,iteration])
                #results.append(particleFlux_vm[ispecies,iteration])
                results.append(particleFlux_vd[ispecies,iteration])
                #results.append(heatFlux_vm[ispecies,iteration])
                #results.append(heatFlux_vd[ispecies,iteration])
                results.append(heatFlux_withoutPhi1[ispecies,iteration])
             else:
                results.append(FSABFlow[ispecies,iteration])
                results.append(particleFlux_vm[ispecies,iteration])
                results.append(heatFlux_vm[ispecies,iteration])
          #results.append(FSABjHatOverRootFSAB2[iteration])
          results.append(FSABjHat[iteration])
          if includePhi1:
             radialCurrents.append(sum(Zs*particleFlux_vd[:,-1]))
          else:
             radialCurrents.append(sum(Zs*particleFlux_vm[:,-1]))
          results.append(radialCurrents[-1])
          outputs.append(results)

    numRuns += 1

    print ("Successfully read run in directory "+directory)

if not atLeastOneDirectorySucceeded:
   print ("Error! There do not seem to be any completed sfincs jobs in subdirectories of this directory.")
   exit(1)

conversionFactorToddpsiHat = f['dPhiHatdpsiHat'][()] / f[varName][()]
conversionFactorToddpsiN   = f['dPhiHatdpsiN'][()]   / f[varName][()]
conversionFactorToddrHat   = f['dPhiHatdrHat'][()]   / f[varName][()]
conversionFactorToddrN     = f['dPhiHatdrN'][()]     / f[varName][()]
conversionFactorToEr       = f['Er'][()]             / f[varName][()]


# Sort by Er:
Ers_sorted = sorted(Ers)
outputs_sorted = []
radialCurrents_sorted = []
print ("Here comes radialCurrents")
print (radialCurrents)
for Er in Ers_sorted:
   outputs_sorted.append(outputs[Ers.index(Er)])
   radialCurrents_sorted.append(radialCurrents[Ers.index(Er)])
 
outputs_array = numpy.array(outputs_sorted)

radialCurrentLabel = "radial current"

if RHSMode > 1:
   yAxisLabels=[]
   for irow in range(transportMatrix.shape[0]):
      for icol in range(transportMatrix.shape[1]):
         yAxisLabels.append("L"+str(irow+1)+str(icol+1))
else:
   if Nspecies==1:
      if includePhi1:
         #yAxisLabels=["FSABFlow","particleFlux\nvm_rHat","particleFlux\nvd_rHat","heatFlux\nvm_rHat","heatFlux\nvd_rHat"]
         yAxisLabels=["FSABFlow","particleFlux\nvm_rHat","particleFlux\nvd_rHat","heatFlux\nvm_rHat","heatFlux\nwithoutPhi1_rHat"]
      else:
         #yAxisLabels=["FSABFlow","particleFlux\nvm_rHat","heatFlux\nvm_rHat","elapsed time (s)"]
         yAxisLabels=["FSABFlow","particleFlux\nvm_rHat","heatFlux\nvm_rHat","source 1","source 2"]
   else:
      if includePhi1:
         yAxisLabels=[]
         for i in range(1,Nspecies+1):
            yAxisLabels.append("FSABFlow (species "+str(i)+")")
            yAxisLabels.append("particleFlux rHat (species "+str(i)+")")
            #yAxisLabels.append("particleFlux vm_rHat (species "+str(i)+")")
            #yAxisLabels.append("particleFlux vd_rHat (species "+str(i)+")")
            #yAxisLabels.append("heatFlux vm_rHat (species "+str(i)+")")
            #yAxisLabels.append("heatFlux vd_rHat (species "+str(i)+")")
            #yAxisLabels.append("heatFlux withoutPhi1_rHat (species "+str(i)+")")
            yAxisLabels.append("heatFlux rHat (species "+str(i)+")")
##         yAxisLabels.append("FSABjHatOverRootFSAB2")
         yAxisLabels.append("FSABjHat")
         yAxisLabels.append(radialCurrentLabel)
        
      else:
         yAxisLabels=[]
         for i in range(1,Nspecies+1):
            yAxisLabels.append("FSABFlow (species "+str(i)+")")
            #yAxisLabels.append("particleFlux_vm_rHat (species "+str(i)+")")
            #yAxisLabels.append("heatFlux_vm_rHat (species "+str(i)+")")
            yAxisLabels.append("particleFlux rHat (species "+str(i)+")")
            yAxisLabels.append("heatFlux rHat (species "+str(i)+")")
         #yAxisLabels.append("FSABjHatOverRootFSAB2")
         yAxisLabels.append("FSABjHat")
         yAxisLabels.append(radialCurrentLabel)
            

numQuantities = len(yAxisLabels)

# ***************************************************
# Solve for the ambipolar E_r
# ***************************************************

# The PchipInterpolator routine seems to return an error when there are only 2 elements to the input array.
if len(Ers_sorted)<3:
    whichInterpolator = 1

ambipolaritySolveSucceeded = False
if doAmbipolaritySolve:
   if max(radialCurrents)>0 and min(radialCurrents)<0:
      try:
         # We have runs that straddle j_psi = 0, so we can find at least 1 ambipolar solution.
         print ("Solving for E_r.")
         # Note: For interpolation it does not work well to use scipy's "interp1d" function, since this tends to give
         # large overshoots for the typical data we get in an Er scan (with sharp structure near Er=0).
         if whichInterpolator==0:
             ##Added a try/except because PchipInterpolator seem to raise an error if there are only two points in the interpolation data
             try:
                 interpolator = PchipInterpolator(Ers_sorted, radialCurrents_sorted)
             except:
                 from scipy.interpolate import interp1d
                 whichInterpolator = 1
                 interpolator = interp1d(Ers_sorted, radialCurrents_sorted)
         else:
             interpolator = interp1d(Ers_sorted, radialCurrents_sorted)
         # In the next line, 500 is arbitrary, and could be replaced by another large integer.
         Er_fine = numpy.linspace(min(Ers), max(Ers), num=500)
         radialCurrent_fine = interpolator(Er_fine)
         # Approximately find points where the sign of radial current flips:
         positiveCurrent = (radialCurrent_fine>0)
         signFlips = (positiveCurrent[:-1] != positiveCurrent[1:])
         numRoots = sum(1 for x in signFlips if x)
         print ("Number of roots found for E_r: ",numRoots)
         roots = []
         #for i in range(numRoots):
         for index,value in enumerate(signFlips):
            if value:
               roots.append(brentq(interpolator,Er_fine[index],Er_fine[index+1]))
         # Convert standard array to numpy array:
         roots = numpy.sort(numpy.array(roots))

         # Classify roots as ion, unstable, or electron root:
         if len(roots)==1:
             if roots[0]*conversionFactorToEr>0:
                 root_types=['electron']
             else:
                 root_types=['ion']
         elif len(roots)==3:
             root_types=['ion','unstable','electron']
         else:
             # If the number of roots is not 1 or 3, then flag the case as complicated:
             root_types=['unknown']*len(roots)

         print ("Roots found for the radial derivative of the electrostatic potential, for various radial coordinates:")
         print ("d PhiHat / d psiHat = ",roots * conversionFactorToddpsiHat)
         print ("d PhiHat / d psiN   = ",roots * conversionFactorToddpsiN)
         print ("d PhiHat / d rHat   = ",roots * conversionFactorToddrHat)
         print ("d PhiHat / d rN     = ",roots * conversionFactorToddrN)
         print ("Er                  = ",roots * conversionFactorToEr)

         print ("Root types:",root_types)
         # Interpolate the other results onto the ambipolar Er(s):
         print ("**** Output quantities at these ambipolar E_r value(s) ****")
         outputs_ambipolar = []
         for iQuantity in range(numQuantities):
            if whichInterpolator==0:
                quantityInterpolator = PchipInterpolator(Ers_sorted, outputs_array[:,iQuantity])
            else:
                quantityInterpolator = interp1d(Ers_sorted, outputs_array[:,iQuantity])
            thisOutput_ambipolar = quantityInterpolator(roots)
            outputs_ambipolar.append(thisOutput_ambipolar)
            print (yAxisLabels[iQuantity],": ",thisOutput_ambipolar)

         ambipolaritySolveSucceeded = True
      except Exception as e:
         print ("For some reason, the ambipolarity solve failed. Here is the error message:")
         print (e)
         # Eventually remove this next line:
         raise
   else:
      print ("All runs have the same sign of radial current, so we cannot solve for the ambipolar E_r.")
else:
   print ("Not doing an ambipolarity solve since RHSMode > 1.")

print ("Time at end of ambipolarity calculations: ", time.time() - start_time)

# ***************************************************
# Organize all the data that will be plotted, or saved in the pickle file
# ***************************************************

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
   xdata.append(Ers_sorted)
   ydata.append(outputs_array[:,iQuantity])
   xlabels.append(varName)
   ylabels.append(yAxisLabels[iQuantity])
   ymins.append(-numpy.inf)
   ymaxs.append(numpy.inf)
   xscales.append('linear')
   yscales.append('linear')

# ***************************************************
# Now make the plot
# ***************************************************

if display_plot:
    fig = plt.figure(figsize=(16,8.5))
    fig.patch.set_facecolor('white')

    for iQuantity in range(numQuantities):
       plt.subplot(numRows,numCols,iQuantity+1)
       plt.plot(xdata[iQuantity],ydata[iQuantity],linespec)
       if ambipolaritySolveSucceeded:
          for i in range(len(roots)):
             plt.plot(roots[i],outputs_ambipolar[iQuantity][i],'xr')
       plt.xlabel(xlabels[iQuantity])
       plt.ylabel(ylabels[iQuantity])
       ymin,ymax = plt.ylim()
       if yAxisLabels[iQuantity]==radialCurrentLabel:
          plt.plot([min(Ers),max(Ers)],[0,0],':k')
          if ambipolaritySolveSucceeded:
             plt.plot(Er_fine,radialCurrent_fine,'-g')   

    plt.tight_layout()
    #titleString = "Plot generated by "+ inspect.getfile(inspect.currentframe()) + "\nRun in "+os.getcwd()
    titleString = "Plot generated by "+ this_filename + "\nRun in "+os.getcwd()
    ax = fig.add_axes([0,0,1,1], frameon=False)
    ax.text(0.5,0.99,titleString,horizontalalignment='center',verticalalignment='top')
print ("Time at end of plotting: ", time.time() - start_time)

# Save info which may be used by sfincsScan_combine:
scanType=2
data = {'scanType':scanType, 'numQuantities':numQuantities, 'numRows':numRows,'numCols':numCols,
        'xdata':xdata, 'ydata':ydata, 'xlabels':xlabels, 'ylabels':ylabels,
        'xscales':xscales, 'yscales':yscales, 'ymins':ymins, 'ymaxs':ymaxs,
        'linespec':linespec, 'nHats':nHats, 'THats':THats}
outputFile = open('sfincsScan.dat','wb')
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
print ("Time at end of first dat file close: ", time.time() - start_time)

filename = 'ambipolarSolutions.dat'
try:
   # Remove any file which might already exist.
   os.remove(filename)
except:
   # An exception will occur if the file does not exist, but we do not care.
   pass

if ambipolaritySolveSucceeded:
   # Save info which may be used by sfincsScanPlot_5 (radial scan):
   data.update({'roots':roots, 'outputs_ambipolar':outputs_ambipolar, 'radius_wish':radius_wish, 'radius_actual':radius_actual, 'root_types':root_types, 'nHats':nHats, 'THats':THats})
   outputFile = open(filename,'wb')
   pickle.dump(data,outputFile)
   outputFile.close()
print ("Time at end of second dat file close: ", time.time() - start_time)

if display_plot:
   plt.show()

if makePDF:
   if len(sys.argv)>2 : #Use the substituted name as file name                                                                                      
      print ("Writing plot to " + os.getcwd() + "/" + sys.argv[2] + ".pdf.")
      plt.savefig(sys.argv[2] + ".pdf")
   else: #Use script name as file name                                                                                                              
      print ("Writing plot to " + os.getcwd() + "/" + os.path.basename(__file__) + ".pdf.")
      plt.savefig(os.path.basename(__file__) + ".pdf")



print ("Time to run sfincsScanPlot_2: ", time.time() - start_time)
