#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.00402022, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -5.99922e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -2.14282e-07, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", 0.00379729, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  -4.23325e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  -9.20425e-09, desiredTolerance)

exit(numFailures > 0)
