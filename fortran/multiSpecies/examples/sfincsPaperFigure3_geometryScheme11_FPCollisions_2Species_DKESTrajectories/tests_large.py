#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

species = 0
numFailures += shouldBe("FSABFlow", species, -8.923796180722230E-003, desiredTolerance)
numFailures += shouldBe("particleFlux", species, 1.110019018426244E-006, desiredTolerance)
numFailures += shouldBe("heatFlux", species, 2.419159180052273E-006, desiredTolerance)

species = 1
numFailures += shouldBe("FSABFlow", species, -1.270506048510871E-003, desiredTolerance)
numFailures += shouldBe("particleFlux", species, 5.357759543503480E-008, desiredTolerance)
numFailures += shouldBe("heatFlux", species, 1.076434969068456E-007, desiredTolerance)

numFailures += shouldBe("FSABjHat", 0, -7.653290132211359E-003, desiredTolerance)

exit(numFailures > 0)
