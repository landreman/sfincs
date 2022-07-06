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
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.00221551, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[0,0;;;]", -2.44012e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,0;;;]", -2.44012e-08 , desiredTolerance)
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

exit(numFailures > 0)
