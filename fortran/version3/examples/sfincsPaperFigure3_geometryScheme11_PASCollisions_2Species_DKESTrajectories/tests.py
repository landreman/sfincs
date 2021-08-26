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
numFailures += shouldBe("FSABFlow[0,0;;;]", -6.982140236880133E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", 1.116826148932465E-006, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", 2.382995682114280E-006, desiredTolerance)

# Check the second species:
numFailures += shouldBe("FSABFlow[1,0;;;]", 7.812106245409801E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", 5.480946026266353E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", 1.133473454775712E-007, desiredTolerance)

numFailures += shouldBe("FSABjHat[0;;;]", -1.479424648228993E-002, desiredTolerance)

exit(numFailures > 0)
