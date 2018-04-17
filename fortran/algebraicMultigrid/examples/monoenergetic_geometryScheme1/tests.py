#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("transportMatrix[0,0;;;]", -0.0929182, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,1;;;]", -1.08278, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,0;;;]", -1.08408, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,1;;;]", 438.826, desiredTolerance)

exit(numFailures > 0)
