import os, inspect

print "This is "+ inspect.getfile(inspect.currentframe())
print "Beginning a convergence scan."

Ntheta = readVariable("Ntheta","int")
NthetaMaxFactor = readScanVariable("NthetaMaxFactor","float")
NthetaMinFactor = readScanVariable("NthetaMinFactor","float")
NthetaNumRuns = readScanVariable("NthetaNumRuns","int")
if Ntheta<0:
    print "Error! Ntheta must be positive"
    exit(1)
if NthetaMaxFactor<0:
    print "Error! NthetaMaxFactor must be positive"
    exit(1)
if NthetaMinFactor<0:
    print "Error! NthetaMinFactor must be positive"
    exit(1)
Nthetas = logspace_odd(NthetaMinFactor*Ntheta,NthetaMaxFactor*Ntheta,NthetaNumRuns)
try:
    Nthetas.remove(Ntheta)
except ValueError:
    pass
print "Nthetas:",Nthetas

Nzeta = readVariable("Nzeta","int")
NzetaMaxFactor = readScanVariable("NzetaMaxFactor","float")
NzetaMinFactor = readScanVariable("NzetaMinFactor","float")
NzetaNumRuns = readScanVariable("NzetaNumRuns","int")
if Nzeta<0:
    print "Error! Nzeta must be positive"
    exit(1)
if NzetaMaxFactor<0:
    print "Error! NzetaMaxFactor must be positive"
    exit(1)
if NzetaMinFactor<0:
    print "Error! NzetaMinFactor must be positive"
    exit(1)
Nzetas = logspace_odd(NzetaMinFactor*Nzeta,NzetaMaxFactor*Nzeta,NzetaNumRuns)
try:
    Nzetas.remove(Nzeta)
except ValueError:
    pass
print "Nzetas:",Nzetas

Nxi = readVariable("Nxi","int")
NxiMaxFactor = readScanVariable("NxiMaxFactor","float")
NxiMinFactor = readScanVariable("NxiMinFactor","float")
NxiNumRuns = readScanVariable("NxiNumRuns","int")
if Nxi<0:
    print "Error! Nxi must be positive"
    exit(1)
if NxiMaxFactor<0:
    print "Error! NxiMaxFactor must be positive"
    exit(1)
if NxiMinFactor<0:
    print "Error! NxiMinFactor must be positive"
    exit(1)
Nxis = logspace_int(NxiMinFactor*Nxi,NxiMaxFactor*Nxi,NxiNumRuns)
try:
    Nxis.remove(Nxi)
except ValueError:
    pass
print "Nxis:",Nxis

NL = readVariable("NL","int")
NLMaxFactor = readScanVariable("NLMaxFactor","float")
NLMinFactor = readScanVariable("NLMinFactor","float")
NLNumRuns = readScanVariable("NLNumRuns","int")
if NL<0:
    print "Error! NL must be positive"
    exit(1)
if NLMaxFactor<0:
    print "Error! NLMaxFactor must be positive"
    exit(1)
if NLMinFactor<0:
    print "Error! NLMinFactor must be positive"
    exit(1)
NLs = logspace_int(NLMinFactor*NL,NLMaxFactor*NL,NLNumRuns)
try:
    NLs.remove(NL)
except ValueError:
    pass
print "NLs:",NLs

Nx = readVariable("Nx","int")
NxMaxFactor = readScanVariable("NxMaxFactor","float")
NxMinFactor = readScanVariable("NxMinFactor","float")
NxNumRuns = readScanVariable("NxNumRuns","int")
if Nx<0:
    print "Error! Nx must be positive"
    exit(1)
if NxMaxFactor<0:
    print "Error! NxMaxFactor must be positive"
    exit(1)
if NxMinFactor<0:
    print "Error! NxMinFactor must be positive"
    exit(1)
Nxs = logspace_int(NxMinFactor*Nx,NxMaxFactor*Nx,NxNumRuns)
try:
    Nxs.remove(Nx)
except ValueError:
    pass
print "Nxs:",Nxs

numRunsInScan = 1+len(Nthetas)+len(Nzetas)+len(Nxis)+len(NLs)+len(Nxs)

baseCase = [Ntheta,Nzeta,Nxi,NL,Nx]

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

if currentIndex != numRunsInScan:
    print "Error! Something went wrong."
    exit(1)

if len(parametersForScan) != len(descriptions):
    print "Error! Something went wrong."
    exit(1)

# See if any runs with the same description already exist.
# This happens if you re-run sfincsScan more than once in the same directory.
runNum = 0
while runNum < numRunsInScan:
    directory = descriptions[runNum]
    if os.path.exists(directory):
        print "Warning: directory "+directory+" already exists, so skipping this run."
        numRunsInScan -= 1
        descriptions.pop(runNum)
        parametersForScan.pop(runNum)
        runNum -= 1
    runNum += 1
        

print
print "Here are the parameters for the "+str(numRunsInScan)+" runs we will launch:"
print "[Ntheta, Nzeta, Nxi, NL, Nx]"
print "----------------------------"
for line in parametersForScan:
    print line

print
print "Here are the directories that will be created:"
print descriptions

while True:
    proceed=raw_input("Should I go ahead and launch these jobs? [y/n] ")
    if proceed=="y" or proceed=="n":
        break
    print "You must enter either y or n."

if proceed=="n":
    exit(0)
print "launching jobs..."

# Read in the job.sfincsScan file:
with open(jobFilename, 'r') as f:
    jobFile = f.readlines()

for runNum in range(numRunsInScan):
    directory = descriptions[runNum]
    print "Beginning to handle job "+directory

    # To be extra safe, check again to see if the directory exists.
    if os.path.exists(directory):
        print "Warning: directory "+directory+" already exists."
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
    # This next function is defined separately for each system in sfincsScan.py
    nameJobFile(thisJobFile,directory)
    f = open(jobFilename,"w")
    f.writelines(thisJobFile)
    f.close()

    # Now copy the input.namelist file:
    f = open(filename,"w")
    for line in inputFile:
        if namelistLineContains(line,"Ntheta"):
            line = "  Ntheta = "+str(parametersForScan[runNum][0])+" ! Set by sfincsScan.\n"
        if namelistLineContains(line,"Nzeta"):
            line = "  Nzeta = "+str(parametersForScan[runNum][1])+" ! Set by sfincsScan.\n"
        if namelistLineContains(line,"Nxi"):
            line = "  Nxi = "+str(parametersForScan[runNum][2])+" ! Set by sfincsScan.\n"
        if namelistLineContains(line,"NL"):
            line = "  NL = "+str(parametersForScan[runNum][3])+" ! Set by sfincsScan.\n"
        if namelistLineContains(line,"Nx"):
            line = "  Nx = "+str(parametersForScan[runNum][4])+" ! Set by sfincsScan.\n"
        f.write(line)
    f.close()

    # Submit the sfincs job:
    try:
        # We need to include .split(" ") to separate the command-line arguments into an array of strings.   
        # I'm not sure why python requires this. 
        submissionResult = subprocess.call(submitCommand.split(" "))
    except:
        print "ERROR! Unable to submit run "+directory+" for some reason."
        raise
    else:
        if submissionResult==0:
            print "No errors submitting job "+directory
        else:
            print "Nonzero exit code returned when trying to submit job "+directory

    os.chdir("..")
