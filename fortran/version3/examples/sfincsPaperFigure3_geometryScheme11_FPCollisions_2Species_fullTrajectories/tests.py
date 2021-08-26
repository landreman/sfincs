#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.03

numFailures = 0

# Check the first species:
numFailures += shouldBe("FSABFlow[0,0;;;]", -9.193502564773663E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", 1.078924587797651E-006, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", 2.337320329279748E-006, desiredTolerance)

# Check the second species:
numFailures += shouldBe("FSABFlow[1,0;;;]", -1.441006015653415E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", 5.355677905882999E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", 1.075782052474404E-007, desiredTolerance)

numFailures += shouldBe("FSABjHat[0;;;]", -7.752496549120249E-003, desiredTolerance)

exit(numFailures > 0)
