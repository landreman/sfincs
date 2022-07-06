#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("FSABFlow[0,1;;;]", 1.316898944950224e-02, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[0,1;;;]", 9.53017e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,1;;;]",  9.53017e-08, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,1;;;]",  1.48084e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,1;;;]", 2.43386e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,1;;;]",  2.43386e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[0,1;;;]", 3.57381e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,1;;;]",  3.57381e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,1;;;]", 1.19420e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,1;;;]",  1.55158e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,1;;;]", 1.55158e-07, desiredTolerance)

exit(numFailures > 0)
