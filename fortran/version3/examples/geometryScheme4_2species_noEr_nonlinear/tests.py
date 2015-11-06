#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# First, perform tests on the first iteration, i.e. the equivalent linear calculation:

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.00221551, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[0,0;;;]", -2.44012e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,0;;;]",  -2.44012e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]",  -1.13488e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,0;;;]", -1.3789e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,0;;;]",  -1.3789e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[0,0;;;]", -1.83009e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,0;;;]",  -1.83009e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -3.31401e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,0;;;]",  -3.49702e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,0;;;]", -3.49702e-07, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", 0.00176907, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[1,0;;;]", -2.44012e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[1,0;;;]",  -2.44012e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  2.16605e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[1,0;;;]", -2.74069e-09, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[1,0;;;]",  -2.74069e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[1,0;;;]", -1.83009e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[1,0;;;]",  -1.83009e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  2.35903e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[1,0;;;]",  5.28943e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[1,0;;;]", 5.28943e-09, desiredTolerance)

# Next, perform tests on the nonlinear solution:

# Species 1
numFailures += shouldBe("FSABFlow[0,2;;;]", 0.00221542, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[0,2;;;]", -2.44018e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,2;;;]",  -2.44017e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,2;;;]",  -1.13491e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,2;;;]", -1.37893e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,2;;;]",  -1.37893e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[0,2;;;]", -1.83013e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,2;;;]",  -1.8302e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,2;;;]",  -3.31408e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,2;;;]",  -3.4971e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,2;;;]", -3.49709e-07, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,2;;;]",  0.00176902, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[1,2;;;]", -2.44018e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[1,2;;;]",  -2.44017e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,2;;;]",  2.1661e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[1,2;;;]", -2.74074e-09, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[1,2;;;]",  -2.74066e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[1,2;;;]", -1.83013e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[1,2;;;]",  -1.83012e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,2;;;]",  2.35911e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[1,2;;;]",  5.28982e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[1,2;;;]", 5.28971e-09, desiredTolerance)

exit(numFailures > 0)
