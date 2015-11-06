#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

species = 0
numFailures += shouldBe("FSABFlow", species, -0.17791531603134328, desiredTolerance)
numFailures += shouldBe("particleFlux", species, -4.02276909126271142E-008, desiredTolerance)
numFailures += shouldBe("heatFlux", species, 1.30783676291250351E-007, desiredTolerance)


exit(numFailures > 0)
