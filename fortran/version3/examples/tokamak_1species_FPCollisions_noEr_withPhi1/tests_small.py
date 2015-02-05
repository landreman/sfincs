#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("FSABFlow[0,0;;;]", 3.31053540633345494E-002, desiredTolerance)

numFailures += shouldBe("particleFlux_vE0_psiHat[0,0;;;]", -2.35529e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,0;;;]",  -2.35529e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]",   2.35599e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,0;;;]",  6.93197e-11, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,0;;;]",   6.93197e-11, desiredTolerance)

numFailures += shouldBe("heatFlux_vE0_psiHat[0,0;;;]", -8.83235e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,0;;;]",  -8.83235e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",   2.45283e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,0;;;]",  1.56959e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,0;;;]",   1.56959e-07, desiredTolerance)

exit(numFailures > 0)
