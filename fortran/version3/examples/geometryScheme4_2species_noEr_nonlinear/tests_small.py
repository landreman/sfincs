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

# Next, perform tests on the nonlinear solution:

# Species 1
numFailures += shouldBe("FSABFlow[0,2;;;]", 4.33578528498198851E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[0,2;;;]", -4.77644868055974515E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,2;;;]",  -4.77638647580263808E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,2;;;]",  -2.22155107125627676E-007, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[0,2;;;]", -2.69919593931225141E-007, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,2;;;]",  -2.69918971883654043E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[0,2;;;]", -3.58233651041980920E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,2;;;]",  -3.58285952216835333E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,2;;;]",  -6.48714892897319018E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,2;;;]",  -6.84543488119002584E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[0,2;;;]", -6.84538258001517084E-007, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,2;;;]", 3.46244718852103776E-003 , desiredTolerance)
numFailures += shouldBe("particleFlux_vE0_psiHat[1,2;;;]", -4.77644868055974714E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[1,2;;;]",  -4.77638647580264602E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,2;;;]",  4.24014357371686910E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vd1_psiHat[1,2;;;]", -5.36305106842878037E-009, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[1,2;;;]",  -5.36242902085776917E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vE0_psiHat[1,2;;;]", -3.58233651041981118E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[1,2;;;]",  -3.58224814125042240E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,2;;;]",  4.61798847231940205E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[1,2;;;]",  1.03574033106897966E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vd1_psiHat[1,2;;;]", 1.03565196189959087E-008, desiredTolerance)

exit(numFailures > 0)
