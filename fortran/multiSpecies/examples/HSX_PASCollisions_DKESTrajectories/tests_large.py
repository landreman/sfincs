#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

species = 0
numFailures += shouldBe("FSABFlow", species, -0.770757651696965, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -1.161190184424680E-005, desiredTolerance)
numFailures += shouldBe("heatFlux", species, -8.819637499826962E-004, desiredTolerance)

species = 1
numFailures += shouldBe("FSABFlow", species, 58.1306611699325, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -2.098719213759062E-005, desiredTolerance)
numFailures += shouldBe("heatFlux", species, -4.482791404112483E-002, desiredTolerance)

numFailures += shouldBe("FSABjHat", 0, -58.9014188216295, desiredTolerance)

exit(numFailures > 0)
