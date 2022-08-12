#!/usr/bin/env python

# This script will not work if called directly.
# From the command line, you should call sfincsScan instead.

import os, inspect
from builtins import input

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsScan_1"
print ("This is "+ this_filename)
print ("Beginning a convergence scan.")

#####################################################################################################################3

#try:
Ntheta = readVariable("Ntheta","int",required=False)

if Ntheta != None:
    baseValuePresent = True
else:
    Ntheta = 15
    baseValuePresent = False
##except:
##    # Default
##  Ntheta = 15
##  baseValuePresent = False

#try:
NthetaMaxFactor = readScanVariable("NthetaMaxFactor","float",required=False)
NthetaMinFactor = readScanVariable("NthetaMinFactor","float",required=False)
NthetaNumRuns = readScanVariable("NthetaNumRuns","int",required=False)

if NthetaMaxFactor != None and NthetaMinFactor != None and NthetaNumRuns != None:
    scanPresent = True
#except:
else:
    NthetaMaxFactor = 1
    NthetaMinFactor = 1
    NthetaNumRuns = 0
    scanPresent = False

if scanPresent and (not baseValuePresent):
    print ("Error! If scanning Ntheta, you must set the base value for Ntheta.")
    exit(1)

if Ntheta<0:
    print ("Error! Ntheta must be positive")
    exit(1)
if NthetaMaxFactor<0:
    print ("Error! NthetaMaxFactor must be positive")
    exit(1)
if NthetaMinFactor<0:
    print ("Error! NthetaMinFactor must be positive")
    exit(1)
Nthetas = logspace_odd(NthetaMinFactor*Ntheta,NthetaMaxFactor*Ntheta,NthetaNumRuns)
try:
    Nthetas.remove(Ntheta)
except ValueError:
    pass
print ("Nthetas:",Nthetas)

#####################################################################################################################3

#try:
Nzeta = readVariable("Nzeta","int",required=False)
if Nzeta != None:
    baseValuePresent = True
#except:
else:
    # Default
    Nzeta = 15
    baseValuePresent = False

#try:
NzetaMaxFactor = readScanVariable("NzetaMaxFactor","float",required=False)
NzetaMinFactor = readScanVariable("NzetaMinFactor","float",required=False)
NzetaNumRuns = readScanVariable("NzetaNumRuns","int",required=False)
if NzetaMaxFactor != None and NzetaMinFactor != None and NzetaNumRuns != None:
    scanPresent = True
#except:
else:
    NzetaMaxFactor = 1
    NzetaMinFactor = 1
    NzetaNumRuns = 0
    scanPresent = False

if scanPresent and (not baseValuePresent):
    print ("Error! If scanning Nzeta, you must set the base value for Nzeta.")
    exit(1)

if Nzeta<0:
    print ("Error! Nzeta must be positive")
    exit(1)
if NzetaMaxFactor<0:
    print ("Error! NzetaMaxFactor must be positive")
    exit(1)
if NzetaMinFactor<0:
    print ("Error! NzetaMinFactor must be positive")
    exit(1)
Nzetas = logspace_odd(NzetaMinFactor*Nzeta,NzetaMaxFactor*Nzeta,NzetaNumRuns)
try:
    Nzetas.remove(Nzeta)
except ValueError:
    pass
print ("Nzetas:",Nzetas)

#####################################################################################################################3

#try:
Nxi = readVariable("Nxi","int",required=False)
if Nxi != None:
    baseValuePresent = True
#except:
else:
    # Default
    Nxi = 16
    baseValuePresent = False

#try:
NxiMaxFactor = readScanVariable("NxiMaxFactor","float",required=False)
NxiMinFactor = readScanVariable("NxiMinFactor","float",required=False)
NxiNumRuns = readScanVariable("NxiNumRuns","int",required=False)
if NxiMaxFactor != None and NxiMinFactor != None and NxiNumRuns != None:
    scanPresent = True
#except:
else:
    NxiMaxFactor = 1
    NxiMinFactor = 1
    NxiNumRuns = 0
    scanPresent = False

if scanPresent and (not baseValuePresent):
    print ("Error! If scanning Nxi you must set the base value for Nxi.")
    exit(1)

if Nxi<0:
    print ("Error! Nxi must be positive")
    exit(1)
if NxiMaxFactor<0:
    print ("Error! NxiMaxFactor must be positive")
    exit(1)
if NxiMinFactor<0:
    print ("Error! NxiMinFactor must be positive")
    exit(1)
Nxis = logspace_int(NxiMinFactor*Nxi,NxiMaxFactor*Nxi,NxiNumRuns)
try:
    Nxis.remove(Nxi)
except ValueError:
    pass
print ("Nxis:",Nxis)

#####################################################################################################################3

#try:
NL = readVariable("NL","int",required=False)
if NL != None:
    baseValuePresent = True
#except:
else:
    # Default
    NL = 4
    baseValuePresent = False

#try:
NLMaxFactor = readScanVariable("NLMaxFactor","float",required=False)
NLMinFactor = readScanVariable("NLMinFactor","float",required=False)
NLNumRuns = readScanVariable("NLNumRuns","int",required=False)
if NLMaxFactor != None and NLMinFactor != None and NLNumRuns != None:
    scanPresent = True
#except:
else:
    NLMaxFactor = 1
    NLMinFactor = 1
    NLNumRuns = 0
    scanPresent = False

if scanPresent and (not baseValuePresent):
    print ("Error! If scanning NL, you must set the base value for NL.")
    exit(1)

if NL<0:
    print ("Error! NL must be positive")
    exit(1)
if NLMaxFactor<0:
    print ("Error! NLMaxFactor must be positive")
    exit(1)
if NLMinFactor<0:
    print ("Error! NLMinFactor must be positive")
    exit(1)
NLs = logspace_int(NLMinFactor*NL,NLMaxFactor*NL,NLNumRuns)
try:
    NLs.remove(NL)
except ValueError:
    pass
print ("NLs:",NLs)

#####################################################################################################################3

#try:
Nx = readVariable("Nx","int",required=False)
if Nx != None:
    baseValuePresent = True
#except:
else:
    # Default
    Nx = 5
    baseValuePresent = False

#try:
NxMaxFactor = readScanVariable("NxMaxFactor","float",required=False)
NxMinFactor = readScanVariable("NxMinFactor","float",required=False)
NxNumRuns = readScanVariable("NxNumRuns","int",required=False)
if NxMaxFactor != None and NxMinFactor != None and NxNumRuns != None:
    scanPresent = True
#except:
else:
    NxMaxFactor = 1
    NxMinFactor = 1
    NxNumRuns = 0
    scanPresent = False

if scanPresent and (not baseValuePresent):
    print ("Error! If scanning Nx, you must set the base value for Nx.")
    exit(1)

if Nx<0:
    print ("Error! Nx must be positive")
    exit(1)
if NxMaxFactor<0:
    print ("Error! NxMaxFactor must be positive")
    exit(1)
if NxMinFactor<0:
    print ("Error! NxMinFactor must be positive")
    exit(1)
Nxs = logspace_int(NxMinFactor*Nx,NxMaxFactor*Nx,NxNumRuns)
try:
    Nxs.remove(Nx)
except ValueError:
    pass
print ("Nxs:",Nxs)

#####################################################################################################################3

#try:
NxPotentialsPerVth = readVariable("NxPotentialsPerVth","float",required=False)
if NxPotentialsPerVth != None:
    baseValuePresent = True
#except:
else:
    # Default
    NxPotentialsPerVth = 40.0
    baseValuePresent = False

#try:
NxPotentialsPerVthMaxFactor = readScanVariable("NxPotentialsPerVthMaxFactor","float",required=False)
NxPotentialsPerVthMinFactor = readScanVariable("NxPotentialsPerVthMinFactor","float",required=False)
NxPotentialsPerVthNumRuns = readScanVariable("NxPotentialsPerVthNumRuns","int",required=False)
if NxPotentialsPerVthMaxFactor != None and NxPotentialsPerVthMinFactor != None and NxPotentialsPerVthNumRuns != None:
    scanPresent = True
#except:
else:
    NxPotentialsPerVthMaxFactor = 1
    NxPotentialsPerVthMinFactor = 1
    NxPotentialsPerVthNumRuns = 0
    scanPresent = False

if scanPresent and (not baseValuePresent):
    print ("Error! If scanning NxPotentialsPerVth, you must set the base value for NxPotentialsPerVth.")
    exit(1)

if NxPotentialsPerVth<0:
    print ("Error! NxPotentialsPerVth must be positive")
    exit(1)
if NxPotentialsPerVthMaxFactor<0:
    print ("Error! NxPotentialsPerVthMaxFactor must be positive")
    exit(1)
if NxPotentialsPerVthMinFactor<0:
    print ("Error! NxPotentialsPerVthMinFactor must be positive")
    exit(1)
NxPotentialsPerVths = logspace(NxPotentialsPerVthMinFactor*NxPotentialsPerVth,NxPotentialsPerVthMaxFactor*NxPotentialsPerVth,NxPotentialsPerVthNumRuns)
for value in NxPotentialsPerVths:
    if abs(value-NxPotentialsPerVth)/value < 1e-6:
        NxPotentialsPerVths.remove(value)
print ("NxPotentialsPerVths:",NxPotentialsPerVths)

#####################################################################################################################3

#try:
xMax = readVariable("xMax","float",required=False)
if xMax != None:
    baseValuePresent = True
#except:
else:
    # Default
    xMax = 5.0
    baseValuePresent = False

#try:
xMaxMaxFactor = readScanVariable("xMaxMaxFactor","float",required=False)
xMaxMinFactor = readScanVariable("xMaxMinFactor","float",required=False)
xMaxNumRuns = readScanVariable("xMaxNumRuns","int",required=False)
if xMaxMaxFactor != None and xMaxMinFactor != None and xMaxNumRuns != None:
    scanPresent = True
#except:
else:
    xMaxMaxFactor = 1
    xMaxMinFactor = 1
    xMaxNumRuns = 0
    scanPresent = False

if scanPresent and (not baseValuePresent):
    print ("Error! If scanning xMax, you must set the base value for xMax.")
    exit(1)

if xMax<0:
    print ("Error! xMax must be positive")
    exit(1)
if xMaxMaxFactor<0:
    print ("Error! xMaxMaxFactor must be positive")
    exit(1)
if xMaxMinFactor<0:
    print ("Error! xMaxMinFactor must be positive")
    exit(1)
xMaxs = logspace(xMaxMinFactor*xMax,xMaxMaxFactor*xMax,xMaxNumRuns)
for value in xMaxs:
    if abs(value-xMax)/value < 1e-6:
        xMaxs.remove(value)
print ("xMaxs:",xMaxs)

#####################################################################################################################3

#try:
solverTolerance = readVariable("solverTolerance","float",required=False)
if solverTolerance != None:
    baseValuePresent = True
#except:
else:
    # Default
    solverTolerance = 1e-6
    baseValuePresent = False

#try:
solverToleranceMaxFactor = readScanVariable("solverToleranceMaxFactor","float",required=False)
solverToleranceMinFactor = readScanVariable("solverToleranceMinFactor","float",required=False)
solverToleranceNumRuns = readScanVariable("solverToleranceNumRuns","int",required=False)
if solverToleranceMaxFactor != None and solverToleranceMinFactor != None and solverToleranceNumRuns != None:
    scanPresent = True
#except:
else:
    solverToleranceMaxFactor = 1
    solverToleranceMinFactor = 1
    solverToleranceNumRuns = 0
    scanPresent = False

if scanPresent and (not baseValuePresent):
    print ("Error! If scanning solverTolerance, you must set the base value for solverTolerance.")
    exit(1)

if solverTolerance<0:
    print ("Error! solverTolerance must be positive")
    exit(1)
if solverToleranceMaxFactor<0:
    print ("Error! solverToleranceMaxFactor must be positive")
    exit(1)
if solverToleranceMinFactor<0:
    print ("Error! solverToleranceMinFactor must be positive")
    exit(1)
solverTolerances = logspace(solverToleranceMinFactor*solverTolerance,solverToleranceMaxFactor*solverTolerance,solverToleranceNumRuns)
for value in solverTolerances:
    if abs((math.log(value)-math.log(solverTolerance))/math.log(value)) < 1e-6:
        solverTolerances.remove(value)
print ("solverTolerances:",solverTolerances)

#####################################################################################################################3

numRunsInScan = 1+len(Nthetas)+len(Nzetas)+len(Nxis)+len(NLs)+len(Nxs)+len(NxPotentialsPerVths)+len(xMaxs)+len(solverTolerances)

baseCase = [Ntheta,Nzeta,Nxi,NL,Nx,NxPotentialsPerVth,xMax,solverTolerance]

parametersForScan = []
for i in range(numRunsInScan):
    parametersForScan.append(list(baseCase))

currentIndex = 1
descriptions = ["baseCase"]

for i in range(len(Nthetas)):
    parametersForScan[currentIndex][0] = Nthetas[i]
    descriptions.append("Ntheta" + str(Nthetas[i]))
    currentIndex += 1

for i in range(len(Nzetas)):
    parametersForScan[currentIndex][1] = Nzetas[i]
    descriptions.append("Nzeta" + str(Nzetas[i]))
    currentIndex += 1

for i in range(len(Nxis)):
    parametersForScan[currentIndex][2] = Nxis[i]
    descriptions.append("Nxi" + str(Nxis[i]))
    currentIndex += 1

for i in range(len(NLs)):
    parametersForScan[currentIndex][3] = NLs[i]
    descriptions.append("NL" + str(NLs[i]))
    currentIndex += 1

for i in range(len(Nxs)):
    parametersForScan[currentIndex][4] = Nxs[i]
    descriptions.append("Nx" + str(Nxs[i]))
    currentIndex += 1

for i in range(len(NxPotentialsPerVths)):
    parametersForScan[currentIndex][5] = NxPotentialsPerVths[i]
    descriptions.append("NxPotentialsPerVth" + "{0:.4g}".format(NxPotentialsPerVths[i]))
    currentIndex += 1

for i in range(len(xMaxs)):
    parametersForScan[currentIndex][6] = xMaxs[i]
    descriptions.append("xMax" + "{0:.4g}".format(xMaxs[i]))
    currentIndex += 1

for i in range(len(solverTolerances)):
    parametersForScan[currentIndex][7] = solverTolerances[i]
    descriptions.append("solverTolerance" + "{0:.4g}".format(solverTolerances[i]))
    currentIndex += 1

if currentIndex != numRunsInScan:
    print ("Error! Something went wrong.")
    exit(1)

if len(parametersForScan) != len(descriptions):
    print ("Error! Something went wrong.")
    exit(1)

# See if any runs with the same description already exist.
# This happens if you re-run sfincsScan more than once in the same directory.
runNum = 0
while runNum < numRunsInScan:
    directory = descriptions[runNum]
    if os.path.exists(directory):
        print ("Warning: directory "+directory+" already exists, so skipping this run.")
        numRunsInScan -= 1
        descriptions.pop(runNum)
        parametersForScan.pop(runNum)
        runNum -= 1
    runNum += 1
        

print ()
print ("Here are the parameters for the "+str(numRunsInScan)+" runs we will launch:")
print ("[Ntheta, Nzeta, Nxi, NL, Nx, NxPotentialsPerVth, xMax, solverTolerance]")
print ("-----------------------------------------------------------------------")
for line in parametersForScan:
    print (line)

print ()
print ("Here are the directories that will be created:")
print (descriptions)

if waitBeforeSubmitting:
    while True:
        #proceed=raw_input("Should I go ahead and launch these "+str(numRunsInScan)+" jobs? [y/n] ")
        proceed=input("Should I go ahead and launch these "+str(numRunsInScan)+" jobs? [y/n] ")
        if proceed=="y" or proceed=="n":
            break
        print ("You must enter either y or n.")

    if proceed=="n":
        exit(0)

print ("launching jobs...")

# Read in the job.sfincsScan file:
with open(jobFilename, 'r') as f:
    jobFile = f.readlines()

for runNum in range(numRunsInScan):
    directory = descriptions[runNum]
    print ("Beginning to handle job "+str(runNum+1)+" of "+str(numRunsInScan)+": "+directory)

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
        if namelistLineContains(line,"Ntheta"):
            line = "  Ntheta = "+str(parametersForScan[runNum][0])+" ! Set by sfincsScan_1.\n"
        if namelistLineContains(line,"Nzeta"):
            line = "  Nzeta = "+str(parametersForScan[runNum][1])+" ! Set by sfincsScan_1.\n"
        if namelistLineContains(line,"Nxi"):
            line = "  Nxi = "+str(parametersForScan[runNum][2])+" ! Set by sfincsScan_1.\n"
        if namelistLineContains(line,"NL"):
            line = "  NL = "+str(parametersForScan[runNum][3])+" ! Set by sfincsScan_1.\n"
        if namelistLineContains(line,"Nx"):
            line = "  Nx = "+str(parametersForScan[runNum][4])+" ! Set by sfincsScan_1.\n"
        if namelistLineContains(line,"NxPotentialsPerVth"):
            line = "  NxPotentialsPerVth = "+str(parametersForScan[runNum][5])+" ! Set by sfincsScan_1.\n"
        if namelistLineContains(line,"xMax"):
            line = "  xMax = "+str(parametersForScan[runNum][6])+" ! Set by sfincsScan_1.\n"
        if namelistLineContains(line,"solverTolerance"):
            line = "  solverTolerance = "+str(parametersForScan[runNum][7])+" ! Set by sfincsScan_1.\n"
        f.write(line)
    f.close()

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
