#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("FSABFlow[0,0;;;]", 1.31689894494766199E-002, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", 2.43400486962836789E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", 1.78992270500802452E-007, desiredTolerance)

exit(numFailures > 0)
