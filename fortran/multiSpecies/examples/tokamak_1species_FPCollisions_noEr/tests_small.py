#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

species = 0
numFailures += shouldBe("FSABFlow", species, 3.31053540633345494E-002, desiredTolerance)
numFailures += shouldBe("heatFlux", species, 9.80304762108666200E-008, desiredTolerance)


exit(numFailures > 0)
