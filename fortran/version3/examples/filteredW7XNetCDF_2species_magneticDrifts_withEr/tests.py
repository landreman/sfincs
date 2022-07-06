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
numFailures += shouldBe("FSABFlow[0,0;;;]", -5.099594120038621E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -2.481551056456413E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -1.089413032476085E-006, desiredTolerance)

# Species 2
##numFailures += shouldBe("FSABFlow[1,0;;;]", -1.294425955700637E-003, desiredTolerance) ##Commented by AM 2018-03
numFailures += shouldBe("FSABFlow[1,0;;;]", -1.30125E-003, 10.0*desiredTolerance) ##Modified by AM 2018-08
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", -5.390144269114580E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", -3.060810653487707E-008, desiredTolerance)

exit(numFailures > 0)
