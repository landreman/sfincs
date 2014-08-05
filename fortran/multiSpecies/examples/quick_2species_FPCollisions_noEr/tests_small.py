#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

species = 0
numFailures += shouldBe("FSABFlow", species, -9.459196929759253E-004, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -8.827673309080052E-009, desiredTolerance)
numFailures += shouldBe("heatFlux", species, -8.418539619831971E-009, desiredTolerance)

species = 1
numFailures += shouldBe("FSABFlow", species, -1.202828424761789E-005, desiredTolerance)
numFailures += shouldBe("particleFlux", species, 9.012860128087754E-011, desiredTolerance)
numFailures += shouldBe("heatFlux", species, 2.007511670631882E-011, desiredTolerance)

numFailures += shouldBe("FSABjHat", 0, -1.018089398461633E-003, desiredTolerance)

exit(numFailures > 0)
