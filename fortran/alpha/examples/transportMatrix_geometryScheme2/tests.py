#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("transportMatrix[0,0;;;]", -0.0114476, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,1;;;]", -0.04203, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,2;;;]",  0.027763, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,0;;;]", -0.0420358, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,1;;;]", -0.31792, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,2;;;]", -0.0226085, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,0;;;]", 0.0277744, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,1;;;]", -0.0225649, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,2;;;]", 25.7517, desiredTolerance)

exit(numFailures > 0)
