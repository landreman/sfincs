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
numFailures += shouldBe("FSABFlow[0,0;;;]", 3.97594692502329251E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -6.09853558592433280E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -2.17831731531552077E-007, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", 3.77579981496331714E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  -4.23844614210280866E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  -9.23322435605986916E-009, desiredTolerance)

exit(numFailures > 0)
