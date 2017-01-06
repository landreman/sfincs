#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("transportMatrix[0,0;;;]", -0.000629935, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,1;;;]", -0.00199572, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,2;;;]", -0.06232, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,0;;;]", -0.00199797, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,1;;;]", -0.010742, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,2;;;]", -0.172505, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,0;;;]", -0.0623828, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,1;;;]", -0.17279, desiredTolerance)

numFailures += shouldBe("transportMatrix[2,2;;;]", 35.1156, desiredTolerance)

exit(numFailures > 0)
