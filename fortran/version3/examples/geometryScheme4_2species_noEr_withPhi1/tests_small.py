#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 4.336504239369484E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[0,0;;;]", -4.775967182538037E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,0;;;]",  -4.775967182678429E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]",  -2.221356667912423E-007, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,0;;;]", -2.698953386166227E-007, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,0;;;]",  -2.698953386180266E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[0,0;;;]", -3.581975386903532E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,0;;;]",  -3.581975386103843E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -6.486635327060895E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,0;;;]",  -6.844832865671280E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,0;;;]", -6.844832865751248E-007, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", 3.462923598912980E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[1,0;;;]", -4.775967182538040E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[1,0;;;]",  -4.775967182534334E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  4.239701830565301E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[1,0;;;]", -5.362653519727391E-009, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[1,0;;;]",  -5.362653519690333E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[1,0;;;]", -3.581975386903532E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[1,0;;;]",  -3.581975386103843E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  4.617413382709388E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[1,0;;;]",  1.035437995817237E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[1,0;;;]", 1.035437995805857E-008, desiredTolerance)

exit(numFailures > 0)
