#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("transportMatrix[0,0;;;]", -9.43187887182127862E-008, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,1;;;]", 3.87601E-005, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,0;;;]", 4.9935E-006, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,1;;;]", 5.1597834378562917, desiredTolerance)

exit(numFailures > 0)
