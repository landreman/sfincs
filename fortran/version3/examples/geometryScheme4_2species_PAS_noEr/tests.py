#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.00166075, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -1.45834e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -3.87758e-07, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", -0.000214135, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  -2.6945e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  -7.05096e-09, desiredTolerance)

exit(numFailures > 0)
