#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Tests for the first iteration, i.e. the equivalent linear run:

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

# Tests after solving the nonlinear equation:

numFailures += shouldBe("FSABFlow[0,2;;;]", 3.31922449972840428E-002, desiredTolerance)

numFailures += shouldBe("particleFlux_vE0_psiHat[0,2;;;]", -2.34987895320895338E-007, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,2;;;]",  -2.34987895320895391E-007, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,2;;;]",   2.35056764754528482E-007 , desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,2;;;]",  6.88694336331440671E-011, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,2;;;]",   6.88694336330911275E-011, desiredTolerance)

numFailures += shouldBe("heatFlux_vE0_psiHat[0,2;;;]", -8.81204607453358409E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,2;;;]",  -8.83959592137427189E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,2;;;]",   2.44972680382320321E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,2;;;]",  1.56852219636984480E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,2;;;]",   1.56576721168577616E-007, desiredTolerance)


exit(numFailures > 0)
