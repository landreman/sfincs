#!/usr/bin/env python

# This script will not work if called directly.
# From the command line, you should call sfincsScan instead.
#
#
# This script launches a number of runs with parameters listed in the file runspec.dat.
# Such a file contains first a number of comment lines, beginning with ! (or % or #). 
# The last line beginning with ! contains the names of the variables listed on the lines
# below it. Each data line corresponds to one sfincs run. The variables not listed
# are left unchanged in the input.namelist. 
#
# On the hydra system the job file variables node, tasks_per_node, 
# ConsumableCpus and wall_clock_limit can also set in runspec.dat.
#
# Here is an example of a runspec.dat file:

#! This is a very short and ridiculous runspec.dat example file
#! The parameter 'wall_clock_limit' is specific to the LoadLeveler
#! used on the hydra cluster
#!
#! rN_wish THats_1  THats_2  dPhiHatdpsiN wall_clock_limit
#      0.5     1.0       1.0          0.0               30
#      0.2     1.5       1.3          0.7               60
#      0.2     2.0       1.7         2.04               90
 

import os, inspect
from builtins import input

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsScan_22"
print ("This is "+ this_filename)
print ("Beginning a runspec.dat scan.")

################## Subroutine to read the runspec file #################################
def readRunspec(therunspecFilename):
    def numberatend(strng):
        thenumb=0
        if strng[-2]=="_":
            try:
                thenumb=int(strng[-1])
            except:
                thenumb=0
        return thenumb

    if not os.path.isfile(therunspecFilename):
        print ("Error! The file "+therunspecFilename+" must be present in the directory from which you call sfincsScan, or the full path must be given.")
        exit(1)
    # Load the runspec file:
    with open(therunspecFilename, 'r') as rsf:
        runspecLines = rsf.readlines()
    
    lind=0
    while runspecLines[lind][0]=="!" or runspecLines[lind][0]=="%" or runspecLines[lind][0]=="#":
        lind+=1
    runparamnames=runspecLines[lind-1][2:].split()
    noElems=len(runparamnames)
    datalines=runspecLines[lind:]
    while len(datalines[-1]) <= 1 :
        datalines=datalines[0:-1]     #remove trailing empty lines
    runparams=[[0 for j in range(noElems)] for i in range(len(datalines))]
    for datalineind in range(len(datalines)):
        thesplittedstring=datalines[datalineind].split()
        if not len(thesplittedstring)==noElems:
            print ("Error! In the runspec file, the number of variable names does not match the number of variables!")
            exit(1)
        for parind in range(noElems):
           thestring=thesplittedstring[parind]
           thestring=thestring.replace('d','e').replace('D','e')
           if ("." in thestring) or ("e" in thestring):
               runparams[datalineind][parind]=float(thestring)
           else:
               runparams[datalineind][parind]=int(thestring)
    
    rsf.close()

    runnumbersatend=[0 for j in range(noElems)]
    runvectorlength=[0 for j in range(noElems)]
    for parind in range(noElems):
        runnumbersatend[parind]=numberatend(runparamnames[parind])
        if runnumbersatend[parind]>0:
            for subind in range(runnumbersatend[parind]):
                runvectorlength[parind+1-runnumbersatend[parind]+subind]=runnumbersatend[parind]

    return runparamnames, runparams, runnumbersatend, runvectorlength

##############################################################################

##try:
##    runspecFilename = readScanVariable("runSpecFile","string",required=False,stringValueCaseSensitive=True)
runspecFilename = readScanVariable("runSpecFile","string",required=False,stringValueCaseSensitive=True)
#except:
if runspecFilename == None :
    # Default
    runspecFilename = "runspec.dat"


runparamnames, runparams, runnumbersatend, runvectorlength=readRunspec(runspecFilename)

numErscans=len(runparams)

##Check that variables are defined in SFINCS input file
for paramname in runparamnames:
    if len(str(paramname).split("_")) != 1 :
        try: ##Check if the last part of the variable name can be transformed into an integer
            dummy = int(str(paramname).split("_")[-1]) #If this fails an exception is raised, then do nothing to paramname 
            paramname = '_'.join(str(paramname).split("_")[:-1])
        except :
            pass
    testReadName = readVariable(paramname, "string", False)
    testReadScanName = readScanVariable(paramname, "string", False)
    if testReadName == None and testReadScanName == None:
        print ("The variable " + paramname + " you wish to scan must be explicitly assigned some value in the appropriate namelist,")
        print ("even though the value will be ignored in the scan.")
        exit(1)
#######################

# First, print the data from runspec.dat to the screen
print (' '+' '.join('%9s' % v for v in runparamnames))
for lind in range(numErscans):
    oneline=""
    for parind in range(len(runparamnames)):
        theformat= "%"+str(max(9,len(runparamnames[parind])))+".3e"
        oneline=oneline+' '+theformat % runparams[lind][parind]
    print (oneline)


while True:
    #proceed=raw_input("Should I go ahead and launch these "+str(numErscans)+" scans? [y/n] ")
    proceed=input("Should I go ahead and launch these "+str(numErscans)+" scans? [y/n] ")
    if proceed=="y" or proceed=="n":
        break
    print ("You must enter either y or n.")

if proceed=="n":
    exit(0)
print ("launching scans...")


# Read in the job.sfincsScan file:
with open(jobFilename, 'r') as f:
    jobFile = f.readlines()

for runNum in range(numErscans):
   dirNum = runNum-1
   while True:
      dirNum+=1
      dirName=str(dirNum)
      if dirNum<10:
         dirName = "0" + dirName
      if not os.path.exists(dirName):
          break
   os.mkdir(dirName)
   os.chdir(dirName)

   jobName="Sfnx."+dirName
   print ("Beginning to handle Er scan "+str(runNum+1)+" of "+str(numErscans)+": "+dirName)

   # Copy the job.sfincsScan file:
   thisJobFile = list(jobFile)
   # This next function is defined separately for each system in sfincsScan
   nameJobFile(thisJobFile,jobName)
   if sfincsSystem=="hydra":
       for parind in range(len(runparamnames)):
           if runparamnames[parind]=="node":
               for line in thisJobFile:
                   if "# @ node =" in line:
                       thisJobFile[thisJobFile.index(line)]="# @ node = "+str(runparams[runNum][parind])+"\n"
           if runparamnames[parind]=="tasks_per_node":
               for line in thisJobFile:
                   if "# @ tasks_per_node =" in line:
                       thisJobFile[thisJobFile.index(line)]="# @ tasks_per_node = "+str(runparams[runNum][parind])+"\n"
           if runparamnames[parind]=="ConsumableCpus":
               for line in thisJobFile:
                   if "# @ resources = ConsumableCpus" in line:
                       thisJobFile[thisJobFile.index(line)]="# @ resources = ConsumableCpus("+str(runparams[runNum][parind])+")\n"
           if runparamnames[parind]=="wall_clock_limit":
               hours=runparams[runNum][parind]//60
               minutes=runparams[runNum][parind]%60
               if hours<10:
                  hour_str="0"+str(hours)
               if minutes<10:
                  minute_str="0"+str(minutes)
               for line in thisJobFile:
                   if "# @ wall_clock_limit =" in line:
                       thisJobFile[thisJobFile.index(line)]="# @ wall_clock_limit = "+hour_str+":"+minute_str+":00\n"

   f = open(jobFilename,"w")
   f.writelines(thisJobFile)
   f.close()

   # Now copy the input.namelist file:
   f = open(filename,"w")
   for line in inputFile:
       
       if line[-22:]==" ! Set by sfincsScan.\n":
           line=line[:-22]+"\n" #Remove this old scan note from the file. Only indicate the present scan variables.
       if namelistLineContainsSS(line,'scanType'):
           line="!ss scanType = 2  ! set by sfincsScan_22.\n"
       for parind in range(len(runparamnames)):
           if runnumbersatend[parind]==0:               
               if namelistLineContains(line,runparamnames[parind]):
                   line = "  "+runparamnames[parind]+" = "+str(runparams[runNum][parind])+" ! Set by sfincsScan_22.\n"
               if namelistLineContainsSS(line,runparamnames[parind]):
                   line = "!ss "+runparamnames[parind]+" = "+str(runparams[runNum][parind])+" ! Set by sfincsScan_22.\n"
           else:
               if namelistLineContains(line,runparamnames[parind][:-2]):
                   if runnumbersatend[parind]==1:
                       str_stack="  "+runparamnames[parind][:-2]+" = "+str(runparams[runNum][parind])
                   else:
                       str_stack=str_stack+"  "+str(runparams[runNum][parind])
                   if runnumbersatend[parind]==runvectorlength[parind]:
                       line = str_stack+" ! Set by sfincsScan_22.\n"
       f.write(line)
   f.close()

   #################################################
   # Submit the Er scan:
   #submitCommand = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/sfincsScan -f"
   submitCommand = os.path.dirname(os.path.abspath(__file__))+"/sfincsScan -f" 
   #############################
   print ('submitCommand = '+submitCommand)
   try:
       # We need to include .split(" ") to separate the command-line arguments into an array of strings.   
       # I'm not sure why python requires this. 
       submissionResult = subprocess.call(submitCommand.split(" "))
       #submissionResult=0
   except:
       print ("ERROR! Unable to submit run "+dirName+" for some reason.")
       raise
   else:
       if submissionResult==0:
           print ("No errors submitting job "+dirName)
       else:
           print ("Nonzero exit code returned when trying to submit job "+dirName)

   os.chdir("..")

