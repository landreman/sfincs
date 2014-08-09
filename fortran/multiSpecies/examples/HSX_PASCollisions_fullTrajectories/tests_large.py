#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

species = 0
numFailures += shouldBe("FSABFlow", species, -0.944579828254162, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -1.129709459372926E-005, desiredTolerance)
numFailures += shouldBe("heatFlux", species, -8.596590745223173E-004, desiredTolerance)

species = 1
numFailures += shouldBe("FSABFlow", species, 58.0722332960741, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -2.074186516830158E-005, desiredTolerance)
numFailures += shouldBe("heatFlux", species, -4.428462393819304E-002, desiredTolerance)

numFailures += shouldBe("FSABjHat", 0, -59.0168131243283, desiredTolerance)

exit(numFailures > 0)
