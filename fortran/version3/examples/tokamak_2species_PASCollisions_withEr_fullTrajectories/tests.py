#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("FSABFlow[0,0;;;]", -0.0654328, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -9.21578e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -7.12684e-07, desiredTolerance)

numFailures += shouldBe("FSABFlow[1,0;;;]", -0.0665217, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", 4.25463e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", 3.82137e-08, desiredTolerance)

exit(numFailures > 0)
