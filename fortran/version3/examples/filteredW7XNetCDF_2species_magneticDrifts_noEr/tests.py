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
numFailures += shouldBe("FSABFlow[0,0;;;]", -3.829380699567663E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -3.046431720535173E-007, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -1.462349694439989E-006, desiredTolerance)

# Species 2
##numFailures += shouldBe("FSABFlow[1,0;;;]", -3.594147540383717E-004, desiredTolerance) ##Commented by AM 2018-03
numFailures += shouldBe("FSABFlow[1,0;;;]", -3.64566E-004, 20.0*desiredTolerance) ##Modified by AM 2018-08
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", -6.257597194755813E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", -3.411992535031735E-008, desiredTolerance)

exit(numFailures > 0)
