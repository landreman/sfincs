#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Check the first species:
numFailures += shouldBe("FSABFlow[0,0;;;]", 15.4685095524092, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -2.145587846683186E-006, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -1.801092374929527E-004, desiredTolerance)

# Check the second species
numFailures += shouldBe("FSABFlow[1,0;;;]", -36.5147587636708, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", -2.081637115299071E-005, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", -4.445459659182958E-002, desiredTolerance)

numFailures += shouldBe("FSABjHat[0;;;]", 51.9832683160801, desiredTolerance)

exit(numFailures > 0)
