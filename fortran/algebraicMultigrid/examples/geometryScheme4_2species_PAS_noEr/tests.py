#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.00163061, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -1.45143e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -3.85562e-07, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", -0.000208742, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  -2.68622e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  -7.01907e-09, desiredTolerance)

exit(numFailures > 0)
