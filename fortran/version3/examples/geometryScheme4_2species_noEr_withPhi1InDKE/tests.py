#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

# First, perform tests on the first iteration, i.e. the equivalent linear calculation:

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 2.21104e-03, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[0,0;;;]", -2.43920e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,0;;;]", -2.43936e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]",  -1.13436e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,0;;;]", -1.37828e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,0;;;]", -1.37829e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[0,0;;;]", -1.82940e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,0;;;]",  -1.83328e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -3.31268e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,0;;;]",  -3.49601e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,0;;;]", -3.49562e-07, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", 1.76459e-03, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[1,0;;;]", -2.43922e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[1,0;;;]",  -2.43936e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  2.16504e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[1,0;;;]", -2.74175e-09, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[1,0;;;]",  -2.74316e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[1,0;;;]", -1.82941e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[1,0;;;]",  -1.82992e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  2.35765e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[1,0;;;]",  5.27729e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[1,0;;;]", 5.28241e-09, desiredTolerance)

# Next, perform tests on the nonlinear solution:

# Species 1
numFailures += shouldBe("FSABFlow[0,2;;;]", 2.21155e-03, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[0,2;;;]", -2.43923e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,2;;;]",  -2.43937e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,2;;;]",  -1.13439e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,2;;;]", -1.37832e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,2;;;]", -1.37833e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[0,2;;;]", -1.82942e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,2;;;]",  -1.83337e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,2;;;]",  -3.31276e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,2;;;]",  -3.49610e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,2;;;]", -3.49571e-07, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,2;;;]",  1.76525e-03, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[1,2;;;]", -2.43924e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[1,2;;;]",  -2.43937e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,2;;;]",  2.16507e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[1,2;;;]", -2.74171e-09, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[1,2;;;]",  -2.74308e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[1,2;;;]", -1.82943e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[1,2;;;]",  -1.82993e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,2;;;]",  2.35769e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[1,2;;;]",  5.27760e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[1,2;;;]", 5.28265e-09, desiredTolerance)

exit(numFailures > 0)
