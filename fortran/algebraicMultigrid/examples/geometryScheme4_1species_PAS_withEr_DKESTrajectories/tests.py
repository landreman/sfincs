#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.00157974, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -7.53922e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -2.54483e-07, desiredTolerance)

exit(numFailures > 0)
