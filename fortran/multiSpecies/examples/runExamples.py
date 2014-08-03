#!/usr/bin/env python

# This script is designed to be called by "make test" in the parent directory,
# not to be called directly.  The reason is that there are several system-dependent
# environment variables used which are set in the makefiles.

import os
import subprocess

try:
    isABatchSystemUsed = os.environ["SFINCS_IS_A_BATCH_SYSTEM_USED"]
except:
    print "Error!  Environment variable SFINCS_IS_A_BATCH_SYSTEM_USED is not set.  This error may be caused by calling runExamples.py directly rather than by calling 'make test'."
    raise

try:
    temp = os.environ["SFINCS_COMMAND_TO_SUBMIT_JOB"]
except:
    print "Error!  Environment variable SFINCS_COMMAND_TO_SUBMIT_JOB is not set.  This error may be caused by calling runExamples.py directly rather than by calling 'make test'."
    raise

wereThereAnyErrors = False
examplesWithErrors = []

# Get a list of the subdirectories:
subdirectories = filter(os.path.isdir, os.listdir("."))

# From this list, keep only subdirectories that have a "tests.py" file:
examplesToRun = []
for subdirectory in subdirectories:
    if os.path.isfile(subdirectory+"/tests.py"):
        examplesToRun.append(subdirectory)

if len(examplesToRun)==0:
    print "There are no examples in the examples/ directory with a tests.py file.  Therefore it is not possible to run any tests."
    exit(1)

print "Based on the subdirectories of examples/ that contain a tests.py file, the following examples may be used as tests:"
for example in examplesToRun:
    print example

if isABatchSystemUsed == "no":
    for subdirectory in examplesToRun:
        print "Preparing to run example: "+subdirectory
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
        # If sfincsOutput.h5 does not exist, there would be an exception, but we can safely ignore it.

        print "Lanching SFINCS..."
        try:
            # Next we launch SFINCS.
            # We need to include .split(" ") to separate the command-line arguments into an array of strings. 
            # I'm not sure why python requires this.
            subprocess.call(os.environ["SFINCS_COMMAND_TO_SUBMIT_JOB"].split(" "))
        except:
            print "An error occurred when attempting to launch sfincs."
            print "This is likely due to an error in the SFINCS_COMMAND_TO_SUBMIT_JOB parameter in the makefile for your system."
            print "The present value of SFINCS_COMMAND_TO_SUBMIT_JOB is: " + os.environ["SFINCS_COMMAND_TO_SUBMIT_JOB"]
            raise

        print "About to run tests on output."

        try:
            testResults = subprocess.call("./tests.py")
        except:
            print "An error occurred when attempting to run tests.py in the following directory:"
            print(os.getcwd)
            raise

        if testResults > 0:
            wereThereAnyErrors = True
            examplesWithErrors.append(subdirectory)

elif isABatchSystemUsed == "yes":
    pass
else:
    print "Error! Environment variable SFINCS_IS_A_BATCH_SYSTEM_USED must be either 'yes' or 'no'.  There is probably a mistake in your system's makefile"
    exit(1)

# If this system does not use a queue manager
#   For each subdirectory
#     If not runLargeExamples,
#       Compute matrix size for this example
#       If matrix is large, continue
#     delete a previous sfincsOutput.h5 file if one exists.
#     Launch job, waiting for it to finish.
#     testResults = subprocess.call(subdirectory+"/tests.py")

# else (i.e. system does use a queue manager)
#   For each subdirectory
#     If not runLargeExamples,
#       Compute matrix size for this example
#       If matrix is large, continue
#     delete a previous sfincsOutput.h5 file if one exists.
#     Submit job but do not wait for it to finish.
#   Periodically ping runs to see if they are done
#     If a run finishes, call tests.py

# Report whether any tests failed.
if wereThereAnyErrors:
    print "-----------------------------------------------"
    print "AT LEAST ONE TEST WAS FAILED."
    print "Examples which failed:"
    for x in examplesWithErrors:
        print "   "+x
else:
    print "-----------------------------------------------"
    print "ALL TESTS COMPLETED SUCCESSFULLY."
