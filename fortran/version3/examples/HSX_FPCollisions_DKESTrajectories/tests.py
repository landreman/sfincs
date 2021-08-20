#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.006

numFailures = 0

# Check the first species:
numFailures += shouldBe("FSABFlow[0,0;;;]", 13.0866765009603, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", 2.205311016736238E-006, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", 1.492152938311221E-004, desiredTolerance)

# Check the second species:
numFailures += shouldBe("FSABFlow[1,0;;;]", -38.2501896319117, desiredTolerance)
#numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  2.104401838496566E-005, desiredTolerance) # Value for Nxi_for_x_option=0
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  2.08554e-05, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", 4.498651646222836E-002, desiredTolerance)

numFailures += shouldBe("FSABjHat[0;;;]", 51.3368661328720, desiredTolerance)

exit(numFailures > 0)
