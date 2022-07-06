#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

# Since the particle flux is ~ 0, use a wider relative tolerance for it:
easyTolerance = 0.5

numFailures = 0

# Tests for the first iteration, i.e. the equivalent linear run:

numFailures += shouldBe("FSABFlow[0,0;;;]", 3.311251583332576E-002, desiredTolerance)

numFailures += shouldBe("particleFlux_vE0_psiHat[0,0;;;]", -1.17795e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,0;;;]",  -1.17795e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]",   1.17827e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,0;;;]", 3.23113e-11, easyTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,0;;;]",  3.13815e-11, easyTolerance)

numFailures += shouldBe("heatFlux_vE0_psiHat[0,0;;;]", -4.41729e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,0;;;]",  -4.41060e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",   1.71681e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,0;;;]",  1.27508e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,0;;;]",  1.27575e-07, desiredTolerance)

# Tests after solving the nonlinear equation:

numFailures += shouldBe("FSABFlow[0,2;;;]", 3.310715818526402e-02, desiredTolerance)

numFailures += shouldBe("particleFlux_vE0_psiHat[0,2;;;]",  -1.17720e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,2;;;]",  -1.17720e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,2;;;]",   1.17752e-07, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,2;;;]",  3.20756e-11, easyTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,2;;;]",   3.23295e-11, easyTolerance)

numFailures += shouldBe("heatFlux_vE0_psiHat[0,2;;;]", -4.41451e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,2;;;]",  -4.41470e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,2;;;]",   1.71620e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,2;;;]",  1.27475e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,2;;;]",   1.27473e-07, desiredTolerance)


exit(numFailures > 0)
