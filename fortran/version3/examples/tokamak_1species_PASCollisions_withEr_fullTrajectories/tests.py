#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("FSABFlow[0,0;;;]", -0.0466595, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -1.5659e-06, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -6.64492e-07, desiredTolerance)

exit(numFailures > 0)
