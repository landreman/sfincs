#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("FSABFlow[0,0;;;]", -0.178341, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -3.74943e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", 1.31912e-07, desiredTolerance)

exit(numFailures > 0)
