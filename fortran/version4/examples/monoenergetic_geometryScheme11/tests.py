#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("transportMatrix[0,0;;;]", -2.8156407763945913E-003, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,1;;;]", -0.00107726, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,0;;;]", -0.00106362, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,1;;;]", 0.885378, desiredTolerance)

exit(numFailures > 0)
