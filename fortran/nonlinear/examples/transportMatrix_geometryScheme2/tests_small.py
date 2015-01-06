#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

# Check the first species:
#numFailures += shouldBe("FSABFlow[0,0;;;]", -9.459196929759253E-004, desiredTolerance)
#numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -8.827673309080052E-009, desiredTolerance)
#numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -8.418539619831971E-009, desiredTolerance)
#numFailures += shouldBe("particleFlux_vd_psiHat[0,0;;;]", -8.827673309080052E-009, desiredTolerance)
#numFailures += shouldBe("heatFlux_vd_psiHat[0,0;;;]", -8.418539619831971E-009, desiredTolerance)


exit(numFailures > 0)
