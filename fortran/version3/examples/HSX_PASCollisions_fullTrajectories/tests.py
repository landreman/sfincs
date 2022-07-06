#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.01

numFailures = 0

# Check the first species:
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.944579828254162, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", 1.129709459372926E-005, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", 8.596590745223173E-004, desiredTolerance)

# Check the second species:
numFailures += shouldBe("FSABFlow[1,0;;;]", -58.0722332960741, desiredTolerance)
#numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", 2.074186516830158E-005, desiredTolerance) # Value for Nxi_for_x_option=0
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", 2.05506e-05, desiredTolerance)
#numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", 4.428462393819304E-002, desiredTolerance) # Value for Nxi_for_x_option=0 
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", 0.0440612, desiredTolerance)

numFailures += shouldBe("FSABjHat[0;;;]", 59.0168131243283, desiredTolerance)

exit(numFailures > 0)
