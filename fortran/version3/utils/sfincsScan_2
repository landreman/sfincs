#!/usr/bin/env python

# This script will not work if called directly.
# From the command line, you should call sfincsScan instead.

import os, inspect
from builtins import input

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsScan_2"
print ("This is "+ this_filename)
print ("Beginning a scan of E_r.")

## inputRadialCoordinateForGradients = readVariable("inputRadialCoordinateForGradients","int")

inputRadialCoordinateForGradients = readVariable("inputRadialCoordinateForGradients","int", False)
if inputRadialCoordinateForGradients == None :
    inputRadialCoordinateForGradients = readDefault("inputRadialCoordinateForGradients","int")

NErs = readScanVariable("NErs","int")
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

generalEr_min = readScanVariable(varName+"Min","float")
generalEr_max = readScanVariable(varName+"Max","float")

generalErs = linspace(generalEr_max, generalEr_min, NErs)

print ("Radial electric field is specified using "+varName)
print ("Here are the "+str(len(generalErs))+" values we will use for this quantity:")
print (generalErs)

directories = [varName+"{:.4g}".format(generalEr) for generalEr in generalErs]
generalErs_copy = list(generalErs)
directories_copy = list(directories)

# See if any runs with the same description already exist.
# This happens if you re-run sfincsScan more than once in the same directory.
for i in range(len(generalErs_copy)):
    directory = directories_copy[i]
    if os.path.exists(directory):
        print ("Warning: directory "+directory+" already exists, so skipping this run.")
        generalErs.remove(generalErs_copy[i])
        directories.remove(directory)

print ()
print ("Here are the directories that will be created:")
print (directories)

if waitBeforeSubmitting:
    while True:
        #proceed=raw_input("Should I go ahead and launch these "+str(len(generalErs))+" jobs? [y/n] ")
        proceed=input("Should I go ahead and launch these "+str(len(generalErs))+" jobs? [y/n] ")
        if proceed=="y" or proceed=="n":
            break
        print ("You must enter either y or n.")

    if proceed=="n":
        exit(0)

print ("launching jobs...")

# Read in the job.sfincsScan file:
with open(jobFilename, 'r') as f:
    jobFile = f.readlines()

for runNum in range(len(generalErs)):
    directory = directories[runNum]
    print ("Beginning to handle job "+str(runNum+1)+" of "+str(len(generalErs))+": "+directory)

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

    #################################################

    f = open(filename,"w")
    for line in inputFile:
        if namelistLineContains(line,varName):
##            line = "  "+varName+" = "+str(Ers[runNum])+" ! Set by sfincsScan_2.\n" 
            continue

        if line.strip().find("&physicsParameters") == 0 :
            line += "  "+varName+" = "+str(generalErs[runNum])+" ! Set by sfincsScan_2.\n"

        f.write(line)
    f.close()

    #################################################

    # Submit the sfincs job:
    try:
        # We need to include .split(" ") to separate the command-line arguments into an array of strings.   
        # I'm not sure why python requires this. 
        submissionResult = subprocess.call(submitCommand.split(" "))
    except:
        print ("ERROR! Unable to submit run "+directory+" for some reason.")
        raise
    else:
        if submissionResult==0:
            print ("No errors submitting job "+directory)
        else:
            print ("Nonzero exit code returned when trying to submit job "+directory)

    os.chdir("..")
