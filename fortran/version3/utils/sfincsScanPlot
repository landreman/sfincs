#!/usr/bin/env python

# This python script plots various types of scan in which the sfincs executable is run
# multiple times.  All this script does is check input.namelist to read scanType, and then
# call the appropriate script sfincsScanPlot_X.  You are also welcome to run
# sfincsScanPlot_X directly.

import os, inspect, sys, subprocess, time

print ("This is "+ inspect.getfile(inspect.currentframe()))

# If called with any command-line arguments other than 'noplot' or 'pdf', then redirect to sfincsScanPlot_combine.
#if len(sys.argv)>1 and (sys.argv[1] != "noplot"):
if len(sys.argv)>1 and (sys.argv[1] != "noplot") and (sys.argv[1].lower() != "pdf"):
    scriptName =  os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/sfincsScanPlot_combine"
    command = sys.argv
    command[0] = scriptName
    start_time = time.time()
    subprocess.call(command)
    print ("Time for subprocess call to sfincsScanPlot_combine: ", time.time() - start_time)
    exit(0)

inputFilename = "input.namelist"
#commentCode = "!ss"

if not os.path.isfile(inputFilename):
    print ("Error! The file "+inputFilename+" must be present in the directory from which you call sfincsScanPlot.")
    exit(1)

# Load some other required subroutines:
start_time = time.time()
#execfile(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/sfincsScan_common")
exec(open(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/sfincsScan_common").read())
print ("Time for exec call to sfincsScanPlot_common: ", time.time() - start_time)

# Load the input file:
with open(inputFilename, 'r') as f:
    inputFile = f.readlines()

# def readScanVariable(varName, intOrFloat):
#     # This subroutine reads the special scan commands in the input.namelist that are hidden from fortran:
#     varName = varName.lower()
#     returnValue = None
    # numValidLines = 0
    # for line in inputFile:
    #     line2 = line.strip().lower()
    #     # We need enough characters for the comment code, varName, =, and value: 
    #     if len(line2)<len(commentCode)+3:
    #         continue

    #     if not line2[:len(commentCode)]==commentCode:
    #         continue

    #     line3 = line2[len(commentCode):].strip()

    #     if len(line3) < len(varName)+2:
    #         continue

    #     if not line3[:len(varName)]==varName:
    #         continue

    #     line4 = line3[len(varName):].strip()

    #     if not line4[0] =="=":
    #         continue

    #     line5 = line4[1:].strip();
        
    #     if intOrFloat=="int":
    #         try:
    #             returnValue = int(line5)
    #             numValidLines += 1
    #         except:
    #             print "Warning! I found a definition for the variable "+varName+" in "+inputFilename+" but I was unable to parse the line to get an integer."
    #             print "Here is the line in question:"
    #             print line
    #     elif intOrFloat=="float":
    #         try:
    #             returnValue = float(line5)
    #             numValidLines += 1
    #         except:
    #             print "Warning! I found a definition for the variable "+varName+" in "+inputFilename+" but I was unable to parse the line to get a float."
    #             print "Here is the line in question:"
    #             print line

    # if returnValue==None:
    #     print "Error! Unable to find a valid setting for the scan variable "+varName+" in "+inputFilename+"."
    #     print "A definition should have the following form:"
    #     if intOrFloat=="int":
    #         print commentCode+" "+varName+" = 1"
    #     elif intOrFloat=="float":
    #         print commentCode+" "+varName+" = 1.5"
    #     exit(1)

    # if numValidLines > 1:
    #     print "Warning! More than 1 valid definition was found for the variable "+varName+". The last one will be used."

    # print "Read "+varName+" = "+str(returnValue)
    # return returnValue

scanType = readScanVariable("scanType","int")

scriptName =  os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/sfincsScanPlot_"+str(scanType) 

if not os.path.isfile(scriptName):
    print ("Error! The file "+scriptName+" does not exist, meaning that plots for scanType = "+str(scanType)+ " are not yet supported.")
    exit(1)

if len(sys.argv)>1 and (sys.argv[1].lower() == "pdf"):
    try:
        command = sys.argv
        command[0] = scriptName
        start_time = time.time()
        subprocess.call(command)
        print ("Time for subprocess call to sfincsScanPlot_x: ", time.time() - start_time)
        exit(0)
    except IOError:
        print ("Unable to run "+scriptName+" even though the file exists.")
        raise
##

try:
    #execfile(scriptName)
    exec(open(scriptName).read())
except IOError:
    print ("Unable to run "+scriptName+" even though the file exists.")
    raise


print ("Good bye!")
