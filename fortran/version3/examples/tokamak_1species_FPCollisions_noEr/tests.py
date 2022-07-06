#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.003

numFailures = 0

numFailures += shouldBe("FSABFlow[0,0;;;]", 3.31053540633345494E-002, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", 9.80304762108666200E-008, desiredTolerance)

exit(numFailures > 0)
