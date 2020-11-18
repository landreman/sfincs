#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.0180645, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -6.66677e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -2.75351e-05, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", -0.00520509, desiredTolerance) 
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", -3.12404e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", -1.87371e-05, desiredTolerance)


# Species 3
numFailures += shouldBe("FSABFlow[2,0;;;]", -7.86796e-06, desiredTolerance) 
numFailures += shouldBe("particleFlux_vm_psiHat[2,0;;;]", 2.06407e-10, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[2,0;;;]", 1.29136e-08, desiredTolerance)


exit(numFailures > 0)
