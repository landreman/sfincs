#!/usr/bin/env python

# This script will not work if called directly.
# From the command line, you should call sfincsScan instead.

import os, inspect
from builtins import input

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsScan_3"
print ("This is "+ this_filename)
print ("Beginning a scan of one specified parameter.")

scanVariable = readScanVariable("scanVariable","string").lower()

# Validate:
warningText = "WARNING: While "+scanVariable+" is an input namelist parameter with a numerical value, it is probably a mistake to scan it."
ErWarning = "WARNING: For scans of the radial electric field, it is preferable to use scanType=2, in which case a solution for ambipolarity will be attempted."
resolutionWarning = "WARNING: For scans of resolution parameters, it is preferable to use scanType=1, so you can scan multiple resolution parameters concurrently."
if scanVariable == "rhsmode":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "inputradialcoordinate":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "inputradialcoordinateforgradients":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "psihat_wish":
    mustBeInteger = False
elif scanVariable == "psin_wish":
    mustBeInteger = False
elif scanVariable == "rhat_wish":
    mustBeInteger = False
elif scanVariable == "rn_wish":
    mustBeInteger = False
elif scanVariable == "geometryscheme":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "b0overbbar":
    mustBeInteger = False
elif scanVariable == "ghat":
    mustBeInteger = False
elif scanVariable == "ihat":
    mustBeInteger = False
elif scanVariable == "iota":
    mustBeInteger = False
elif scanVariable == "epsilon_t":
    mustBeInteger = False
elif scanVariable == "epsilon_h":
    mustBeInteger = False
elif scanVariable == "helicity_l":
    mustBeInteger = True
elif scanVariable == "helicity_n":
    mustBeInteger = True
elif scanVariable == "epsilon_antisymm":
    mustBeInteger = False
elif scanVariable == "helicity_antisymm_l":
    mustBeInteger = True
elif scanVariable == "helicity_antisymm_n":
    mustBeInteger = True
elif scanVariable == "psiahat":
    mustBeInteger = False
elif scanVariable == "ahat":
    mustBeInteger = False
elif scanVariable == "vmecradialoption":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "min_bmn_to_load":
    mustBeInteger = False
elif scanVariable == "zs":
    mustBeInteger = False
elif scanVariable == "mhats":
    mustBeInteger = False
elif scanVariable == "nhats":
    mustBeInteger = False
elif scanVariable == "thats":
    mustBeInteger = False
elif scanVariable == "dnhatdpsihats":
    mustBeInteger = False
elif scanVariable == "dthatdpsihats":
    mustBeInteger = False
elif scanVariable == "dnhatdpsins":
    mustBeInteger = False
elif scanVariable == "dthatdpsins":
    mustBeInteger = False
elif scanVariable == "dnhatdrhats":
    mustBeInteger = False
elif scanVariable == "dthatdrhats":
    mustBeInteger = False
elif scanVariable == "dnhatdrns":
    mustBeInteger = False
elif scanVariable == "dthatdrns":
    mustBeInteger = False
elif scanVariable == "delta":
    mustBeInteger = False
elif scanVariable == "alpha":
    mustBeInteger = False
elif scanVariable == "eparallelhat":
    mustBeInteger = False
elif scanVariable == "dphihatdpsihat":
    print (ErWarning)
    mustBeInteger = False
elif scanVariable == "dphihatdpsin":
    print (ErWarning)
    mustBeInteger = False
elif scanVariable == "dphihatdrhat":
    print (ErWarning)
    mustBeInteger = False
elif scanVariable == "dphihatdrn":
    print (ErWarning)
    mustBeInteger = False
elif scanVariable == "nu_n":
    mustBeInteger = False
elif scanVariable == "krook":
    mustBeInteger = False
elif scanVariable == "nuprime":
    mustBeInteger = False
elif scanVariable == "estar":
    mustBeInteger = False
elif scanVariable == "collisionoperator":
    mustBeInteger = True
elif scanVariable == "constraintscheme":
    mustBeInteger = True
elif scanVariable == "magneticdriftscheme":
    mustBeInteger = True
elif scanVariable == "ntheta":
    print (resolutionWarning)
    mustBeInteger = True
elif scanVariable == "nzeta":
    print (resolutionWarning)
    mustBeInteger = True
elif scanVariable == "nxi":
    print (resolutionWarning)
    mustBeInteger = True
elif scanVariable == "nl":
    print (resolutionWarning)
    mustBeInteger = True
elif scanVariable == "nx":
    print (resolutionWarning)
    mustBeInteger = True
elif scanVariable == "nxpotentialspervth":
    print (resolutionWarning)
    mustBeInteger = False
elif scanVariable == "xmax":
    print (resolutionWarning)
    mustBeInteger = False
elif scanVariable == "solvertolerance":
    print (resolutionWarning)
    mustBeInteger = False
elif scanVariable == "thetaderivativescheme":
    mustBeInteger = True
elif scanVariable == "zetaderivativescheme":
    mustBeInteger = True
elif scanVariable == "xgridscheme":
    mustBeInteger = True
elif scanVariable == "xgrid_k":
    mustBeInteger = True
elif scanVariable == "xpotentialsgridscheme":
    mustBeInteger = True
elif scanVariable == "exbderivativescheme":
    mustBeInteger = True
elif scanVariable == "magneticdriftderivativescheme":
    mustBeInteger = True
elif scanVariable == "whichparallelsolvertofactorpreconditioner":
    mustBeInteger = True
elif scanVariable == "petscpreallocationstrategy":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "preconditioner_species":
    mustBeInteger = True
elif scanVariable == "preconditioner_x":
    mustBeInteger = True
elif scanVariable == "preconditioner_x_min_l":
    mustBeInteger = True
elif scanVariable == "preconditioner_theta":
    mustBeInteger = True
elif scanVariable == "preconditioner_theta_min_l":
    mustBeInteger = True
elif scanVariable == "preconditioner_zeta":
    mustBeInteger = True
elif scanVariable == "preconditioner_zeta_min_l":
    mustBeInteger = True
elif scanVariable == "preconditioner_xi":
    mustBeInteger = True
elif scanVariable == "export_f_theta_option":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "export_f_zeta_option":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "export_f_theta":
    print (warningText)
    mustBeInteger = False
elif scanVariable == "export_f_zeta":
    print (warningText)
    mustBeInteger = False
elif scanVariable == "export_f_xi_option":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "export_f_xi":
    print (warningText)
    mustBeInteger = False
elif scanVariable == "export_f_x_option":
    print (warningText)
    mustBeInteger = True
elif scanVariable == "export_f_x":
    print (warningText)
    mustBeInteger = False
elif scanVariable == "ripplescale":
    mustBeInteger = False
else:
    print ("Error! The selected scanVariable is either not a valid input namelist parameter or it is one that does not take a numerical value.")
    exit(1)

print ("The variable that will be scanned is "+scanVariable)

if mustBeInteger:
    print ("Only integer values of this variable are permitted.")
    scanVarTypeStr = "int"
else:
    print ("Non-integer values of this variable are permitted.")
    scanVarTypeStr = "float"

scanVariableMin = readScanVariable("scanVariableMin",scanVarTypeStr)
scanVariableMax = readScanVariable("scanVariableMax",scanVarTypeStr)

try:
    dummy = readVariable(scanVariable,scanVarTypeStr)
except:
    print ("The variable you wish to scan must be explicitly assigned some value in the appropriate namelist,")
    print ("even though the value will be ignored in the scan.")
    exit(1)

scanVariableN = readScanVariable("scanVariableN","int")
scanVariableScale = readScanVariable("scanVariableScale","string")
if scanVariableScale == "log" or scanVariableScale == "logarithmic":
    if scanVariableMin * scanVariableMax < 0:
        print ("ERROR: scanVariableMin and scanVariableMax must have the same sign for a logarithmic scan.")
        exit(1)

    scanVariableValues = logspace(abs(scanVariableMin), abs(scanVariableMax), scanVariableN)
    if scanVariableMin<0:
        scanVariableValues[:] = [-x for x in scanVariableValues]

elif scanVariableScale == "lin" or scanVariableScale == "linear":
    scanVariableValues = linspace(scanVariableMin, scanVariableMax, scanVariableN)
else:
    print ("ERROR: You must set scanVariableScale to either linear or log.")
    exit(1)


if mustBeInteger:
    scanVariableValues[:] = [int(x) for x in scanVariableValues]
    # Remove duplicates:
    scanVariableValues = list(set(scanVariableValues))

print ("Here are the values of "+scanVariable+" we will use for this scan:")
print (scanVariableValues)

directories = [scanVariable+"_"+"{:.4g}".format(scanVariableValue) for scanVariableValue in scanVariableValues]

scanVariableValues_copy = list(scanVariableValues)
directories_copy = list(directories)

# See if any runs with the same description already exist.
# This happens if you re-run sfincsScan more than once in the same directory.
for i in range(len(scanVariableValues_copy)):
    directory = directories_copy[i]
    if os.path.exists(directory):
        print ("Warning: directory "+directory+" already exists, so skipping this run.")
        scanVariableValues.remove(scanVariableValues_copy[i])
        directories.remove(directory)

print ()
print ("Here are the directories that will be created:")
print (directories)

while True:
    #proceed=raw_input("Should I go ahead and launch these "+str(len(scanVariableValues))+" jobs? [y/n] ")
    proceed=input("Should I go ahead and launch these "+str(len(scanVariableValues))+" jobs? [y/n] ")
    if proceed=="y" or proceed=="n":
        break
    print ("You must enter either y or n.")

if proceed=="n":
    exit(0)
print ("launching jobs...")

# Read in the job.sfincsScan file:
with open(jobFilename, 'r') as f:
    jobFile = f.readlines()

for runNum in range(len(scanVariableValues)):
    directory = directories[runNum]
    print ("Beginning to handle job "+str(runNum+1)+" of "+str(len(scanVariableValues))+": "+directory)

    # To be extra safe, check again to see if the directory exists.
    if os.path.exists(directory):
        print ("Warning: directory "+directory+" already exists.")
        i = -1
        while True:
            i += 1
            directory2 = directory+"_"+str(i)
            if not os.path.exists(directory2):
                break
        directory = directory2
    os.makedirs(directory)
    os.chdir(directory)

    # Copy the job.sfincsScan file:
    thisJobFile = list(jobFile)
    # This next function is defined separately for each system in sfincsScan
    nameJobFile(thisJobFile,directory)
    f = open(jobFilename,"w")
    f.writelines(thisJobFile)
    f.close()

    # Now copy the input.namelist file:
    f = open(filename,"w")
    for line in inputFile:
        # This next line works because previously we have ensured the scan variable explicitly appears
        # in a namelist.  Otherwise we would have had to possibly add the variable to the appropriate namelist.
        if namelistLineContains(line,scanVariable):
            line = "  "+scanVariable+" = "+str(scanVariableValues[runNum])+" ! Set by sfincsScan_3.\n"
        f.write(line)
    f.close()

    # Submit the sfincs job:
    try:
        # We need to include .split(" ") to separate the command-line arguments into an array of strings.   
        # I'm not sure why python requires this. 
        submissionResult = subprocess.call(submitCommand.split(" "))
        #submissionResult=0
    except:
        print ("ERROR! Unable to submit run "+directory+" for some reason.")
        raise
    else:
        if submissionResult==0:
            print ("No errors submitting job "+directory)
        else:
            print ("Nonzero exit code returned when trying to submit job "+directory)

    os.chdir("..")
