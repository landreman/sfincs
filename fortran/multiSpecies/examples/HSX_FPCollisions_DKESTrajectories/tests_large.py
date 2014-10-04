#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.006

numFailures = 0

species = 0
numFailures += shouldBe("FSABFlow", species, -13.0866765009603, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -2.205311016736238E-006, desiredTolerance)
numFailures += shouldBe("heatFlux", species, -1.492152938311221E-004, desiredTolerance)

species = 1
numFailures += shouldBe("FSABFlow", species, 38.2501896319117, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -2.104401838496566E-005, desiredTolerance)
numFailures += shouldBe("heatFlux", species, -4.498651646222836E-002, desiredTolerance)

numFailures += shouldBe("FSABjHat", 0, -51.3368661328720, desiredTolerance)

exit(numFailures > 0)
