#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 3.94745358179843242E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -6.91502213934273899E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,0;;;]", 7.93408957611179042E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -2.28094633391901390E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,0;;;]",  5.93648871465064875E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_withoutPhi1_psiHat[0,0;;;]",  -2.18166464527565219E-007, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", 3.75563367693844971E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  -1.21915380196313510E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[1,0;;;]", 7.93408957611178877E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  -1.91753662031250112E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[1,0;;;]",  5.94010581069900424E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_withoutPhi1_psiHat[1,0;;;]",  -9.24719733878882335E-009, desiredTolerance)

exit(numFailures > 0)
