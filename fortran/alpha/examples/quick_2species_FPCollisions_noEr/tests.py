#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Check the first species:
numFailures += shouldBe("FSABFlow[0,0;;;]", -0.0134817, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -1.36021e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -3.38547e-08, desiredTolerance)

# Check the second species:
numFailures += shouldBe("FSABFlow[1,0;;;]", -0.000166248, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", 1.03276e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", 2.9169e-10, desiredTolerance)

numFailures += shouldBe("FSABjHat[0;;;]", -0.0144792, desiredTolerance)

exit(numFailures > 0)
