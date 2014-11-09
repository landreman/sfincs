#!/usr/bin/env python

# This script is designed to be called by "make test" in the parent directory,
# not to be called directly.  The reason is that there are several system-dependent
# variables used which are set in the makefiles.

import os
import subprocess
import time

def verifyVariableExists(str):
    try:
        temp = os.environ[str]
    except:
        print "Error!  Variable "+str+" is not set.  This error may be caused by calling runExamples.py directly rather than by calling 'make test'."
        raise

    return temp

isABatchSystemUsed = verifyVariableExists("SFINCS_IS_A_BATCH_SYSTEM_USED")
runLargeExamplesStr = verifyVariableExists("SFINCS_RUN_LARGE_EXAMPLES")
commandToSubmitJob = verifyVariableExists("SFINCS_COMMAND_TO_SUBMIT_JOB")
retestStr = verifyVariableExists("SFINCS_RETEST")

if retestStr=="yes":
    retest=True
elif retestStr=="no":
    retest=False
else:
    print "Error! SFINCS_RETEST must be either 'yes' or 'no'. There is likely an error in the main makefile."
    exit(1)

if runLargeExamplesStr=="yes":
    runLargeExamples = True
elif runLargeExamplesStr=="no":
    runLargeExamples = False
else:
    print "Error! SFINCS_RUN_LARGE_EXAMPLES must be either 'yes' or 'no'.  There is likely an error in the makefile for your system."
    exit(1)

# Try executing h5dump with no arguments, just to see if h5dump is available:
        
try:
    p = subprocess.Popen("h5dump", stdout=subprocess.PIPE)
except:
    print "Error! Use of 'make test' for SFINCS requires the utility h5dump, which is usually available whenever the HDF5 library is available. However, an exception was raised when trying to execute h5dump just now."
    raise

wereThereAnyErrors = False
examplesWithErrors = []

# Get a list of the subdirectories:
subdirectories = filter(os.path.isdir, os.listdir("."))

# From this list, keep only subdirectories that have a "tests.py" file:
# For convenience, run small examples before large ones. Hence, we loop twice.
examplesToRun = []
for subdirectory in subdirectories:
    hasSmall = os.path.isfile(subdirectory+"/tests_small.py")
    hasLarge = os.path.isfile(subdirectory+"/tests_large.py")
    if hasSmall and hasLarge:
        print "Error! The subdirectory examples/"+subdirectory+" has both a tests_small.py and tests_large.py file. It must have one or the other or neither, but not both."
        exit(1)
    if hasSmall:
        examplesToRun.append(subdirectory)
if runLargeExamples:
    for subdirectory in subdirectories:
        hasLarge = os.path.isfile(subdirectory+"/tests_large.py")
        if hasLarge:
            examplesToRun.append(subdirectory)

if len(examplesToRun)==0:
    if runLargeExamples:
        print "There are no examples in the examples/ directory with a tests_small.py file or a tests_large.py file.  Therefore it is not possible to run any tests."
    else:
        print "There are no examples in the examples/ directory with a tests_small.py file.  Therefore it is not possible to run any tests."
    exit(1)

if isABatchSystemUsed == "yes":
    print
    print "Note: in order to make this testing program platform-independent, we test for whether jobs have completed"
    print "by trying to read sfincsOutput.h5 using h5dump rather than by directly checking the batch queue.  Therefore"
    print "if any of the example jobs crash before opening sfincsOutput.h5, or without cleanly writing sfincsOutput.h5,"
    print "this testing program will not realize the jobs have crashed."
    print 

if runLargeExamples:
    print "Based on the subdirectories of examples/ that contain a tests_small.py or tests_large.py file, the following examples will be used as tests:"
else:
    print "Based on the subdirectories of examples/ that contain a tests_small.py file, the following examples will be used as tests:"

for example in examplesToRun:
    print "   " + example

if isABatchSystemUsed == "no":
    for subdirectory in examplesToRun:
        print "Preparing to check example: "+subdirectory
        try:
            os.chdir(subdirectory)
        except:
            print "Error occurred when trying to change directory to "+subdirectory
            raise

        print "Moved to working directory "+os.getcwd()

        if not retest:
            try:
                os.remove("sfincsOutput.h5")
            except:
                pass
            # If sfincsOutput.h5 does not exist, there will be an exception, but we can safely ignore it.

            print "Lanching SFINCS..."
            try:
                # Next we launch SFINCS.
                # We need to include .split(" ") to separate the command-line arguments into an array of strings. 
                # I'm not sure why python requires this.
                subprocess.call(commandToSubmitJob.split(" "))
            except:
                print "An error occurred when attempting to launch sfincs."
                print "This is likely due to an error in the SFINCS_COMMAND_TO_SUBMIT_JOB parameter in the makefile for your system."
                print "The present value of SFINCS_COMMAND_TO_SUBMIT_JOB is: " + os.environ["SFINCS_COMMAND_TO_SUBMIT_JOB"]
                raise

        print "About to run tests on output."

        hasSmall = os.path.isfile("tests_small.py")
        if hasSmall:
            filename = "tests_small.py"
        else:
            filename = "tests_large.py"

        try:
            testResults = subprocess.call("./"+filename)
        except:
            print "An error occurred when attempting to run "+filename+" in the following directory:"
            print(os.getcwd)
            raise

        if testResults > 0:
            wereThereAnyErrors = True
            examplesWithErrors.append(subdirectory)

        # Step back one directory
        os.chdir("..")

    print "-----------------------------------------------"
    print "Done with tests."
    print "Examples attempted:"
    for subdirectory in examplesToRun:
        print "  " + subdirectory

elif isABatchSystemUsed == "yes":
    examplesSubmitted = []
    examplesNotSubmitted = []
    devnull = open(os.devnull,'w')

    if retest:
        examplesSubmitted = examplesToRun
    else:
        for subdirectory in examplesToRun:
            print "Preparing to submit example: "+subdirectory
            try:
                os.chdir(subdirectory)
            except:
                print "Error occurred when trying to change directory to "+subdirectory
                raise
            
            print "Moved to working directory "+os.getcwd()
            
            try:
                os.remove("sfincsOutput.h5")
            except:
                pass
            # If sfincsOutput.h5 does not exist, there will be an exception, but we can safely ignore it.

            print "Submitting job for example "+subdirectory+"..."
            try:
                # Next we launch SFINCS.
                # We need to include .split(" ") to separate the command-line arguments into an array of strings. 
                # I'm not sure why python requires this.
                submissionResult = subprocess.call(commandToSubmitJob.split(" "))
            except:
                examplesNotSubmitted.append(subdirectory)
                print "Unable to submit example "+subdirectory+" for some reason. The problem may be SFINCS_COMMAND_TO_SUBMIT_JOB for your system. Skipping this example."
            else:
                if submissionResult==0:
                    examplesSubmitted.append(subdirectory)
                    print "No errors submitting example "+subdirectory+"."
                else:
                    examplesNotSubmitted.append(subdirectory)
                    print "Nonzero exit code returned when trying to submit example "+subdirectory+". Skipping this example."

                # Step back one directory
                os.chdir("..")

    if len(examplesSubmitted)==0:
        print "Unable to submit any of the examples."
        exit(1)

    #status = {subdirectory:"unprocessed" for subdirectory in examplesSubmitted}
    status = {}
    for subdirectory in examplesSubmitted:
        status[subdirectory] = "unprocessed"

    keepGoing = True
    # At this point we have submitted all the examples we were going to try submitting.
    # Periodically ping each example to see if its sfincsOutput.h5 output is readable.
    while keepGoing:
        time.sleep(4)
        print "Checking whether jobs have finished (by trying to read sfincsOutput.h5 files.) Press Ctrl-C to quit."
        keepGoing = False
        for subdirectory in examplesSubmitted:
            filename = subdirectory+"/sfincsOutput.h5"
            if status[subdirectory]=="completedWithErrors":
                print " - Example "+subdirectory+" completed with at least one test failed."
            elif status[subdirectory]=="completedWithoutErrors":
                print " + Example "+subdirectory+" completed and passed all tests."
            elif not os.path.isfile(filename):
                print "   Example "+subdirectory+" has not yet started."
                keepGoing = True
            else:
                try:
                    result=subprocess.call(["h5dump",filename],stdout=devnull,stderr=devnull)
                except:
                    print "Error occurred attempting to run h5dump"
                    raise

                if result != 0:
                    print "   Example "+subdirectory+" has started but I cannot yet read its sfincsOutput.h5."
                    keepGoing = True
                else:
                    print "   Example "+subdirectory+" has completed. Running tests..."
                    os.chdir(subdirectory)

                    hasSmall = os.path.isfile("tests_small.py")
                    if hasSmall:
                        testFilename = "tests_small.py"
                    else:
                        testFilename = "tests_large.py"
                        
                    try:
                        testResults = subprocess.call("./"+testFilename)
                    except:
                        print "An error occurred when attempting to run "+testFilename+" in the following directory:"
                        print(os.getcwd)
                        raise

                    if testResults > 0:
                        wereThereAnyErrors = True
                        examplesWithErrors.append(subdirectory)
                        status[subdirectory] = "completedWithErrors"
                    else:
                        status[subdirectory] = "completedWithoutErrors"

                    # Step back one directory
                    os.chdir("..")
        print

    print "-----------------------------------------------"
    print "Done with tests."
    print "Examples attempted:"
    for subdirectory in examplesSubmitted:
        print "  " + subdirectory

else:
    print "Error! Environment variable SFINCS_IS_A_BATCH_SYSTEM_USED must be either 'yes' or 'no'.  There is probably a mistake in your system's makefile"
    exit(1)

print
# Report whether any tests failed.
if wereThereAnyErrors:
    print "AT LEAST ONE TEST WAS FAILED."
    print "Examples which failed:"
    for x in examplesWithErrors:
        print "   "+x
else:
    print "ALL TESTS THAT WERE RUN WERE PASSED SUCCESSFULLY."

print
