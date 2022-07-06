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
numFailures += shouldBe("FSABFlow[0,1;;;]", 3.975949427891570E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,1;;;]", -6.903048168949459E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[0,1;;;]", 8.029232903602056E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,1;;;]", -2.278881358190897E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[0,1;;;]",  6.026691785508448E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_withoutPhi1_psiHat[0,1;;;]", -2.178417007800692E-007, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,1;;;]", 3.775807768573963E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,1;;;]", -1.228316424896260E-008, desiredTolerance)
numFailures += shouldBe("particleFlux_vE_psiHat[1,1;;;]", 8.029232903622479E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,1;;;]",  -1.928912254970882E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vE_psiHat[1,1;;;]",  6.011305636773584E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_withoutPhi1_psiHat[1,1;;;]", -9.243435350916246E-009, desiredTolerance)

exit(numFailures > 0)
