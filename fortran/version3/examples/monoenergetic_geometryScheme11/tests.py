#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("transportMatrix[0,0;;;]", -2.84704283288131658E-003, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,1;;;]", -1.07314943693241670E-003, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,0;;;]", -1.07288314523113848E-003, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,1;;;]", 0.88506147786933709, desiredTolerance)

exit(numFailures > 0)
