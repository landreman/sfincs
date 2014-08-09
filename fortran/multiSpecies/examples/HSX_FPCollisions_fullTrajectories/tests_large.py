#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

species = 0
numFailures += shouldBe("FSABFlow", species, -15.4685095524092, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -2.145587846683186E-006, desiredTolerance)
numFailures += shouldBe("heatFlux", species, -1.801092374929527E-004, desiredTolerance)

species = 1
numFailures += shouldBe("FSABFlow", species, 36.5147587636708, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -2.081637115299071E-005, desiredTolerance)
numFailures += shouldBe("heatFlux", species, -4.445459659182958E-002, desiredTolerance)

numFailures += shouldBe("FSABjHat", 0, -51.9832683160801, desiredTolerance)

exit(numFailures > 0)
