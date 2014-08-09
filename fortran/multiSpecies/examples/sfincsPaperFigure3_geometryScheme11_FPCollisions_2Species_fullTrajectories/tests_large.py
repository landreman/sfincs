#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

species = 0
numFailures += shouldBe("FSABFlow", species, -9.193502564773663E-003, desiredTolerance)
numFailures += shouldBe("particleFlux", species, 1.078924587797651E-006, desiredTolerance)
numFailures += shouldBe("heatFlux", species, 2.337320329279748E-006, desiredTolerance)

species = 1
numFailures += shouldBe("FSABFlow", species, -1.441006015653415E-003, desiredTolerance)
numFailures += shouldBe("particleFlux", species, 5.355677905882999E-008, desiredTolerance)
numFailures += shouldBe("heatFlux", species, 1.075782052474404E-007, desiredTolerance)

numFailures += shouldBe("FSABjHat", 0, -7.752496549120249E-003, desiredTolerance)

exit(numFailures > 0)
