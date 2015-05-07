#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", -7.742907217793762E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]",  -8.063408399439195E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -4.337518167984696E-006, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", -1.034603693826060E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  -1.171846913091040E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  -6.254944896891315E-008, desiredTolerance)

exit(numFailures > 0)
