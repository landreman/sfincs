#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

## Old values for Nxi_for_x_option=0:
## Species 1
#numFailures += shouldBe("FSABFlow[0,0;;;]", -0.00417179, desiredTolerance)
#numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]",  -4.34448e-07, desiredTolerance)
#numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -2.33701e-06, desiredTolerance)
#
## Species 2
#numFailures += shouldBe("FSABFlow[1,0;;;]", -0.000557432, desiredTolerance)
#numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  -6.31378e-09, desiredTolerance)
#numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  -3.3701e-08, desiredTolerance)

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", -0.00417896, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -4.36428e-07 , desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -2.34336e-06, desiredTolerance)

# Species 2
##numFailures += shouldBe("FSABFlow[1,0;;;]", -0.00056099, desiredTolerance) ##Commented by AM 2018-03
numFailures += shouldBe("FSABFlow[1,0;;;]", -0.000561962, 5.0*desiredTolerance) ##Modified by AM 2018-08
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", -6.31647e-09, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  -3.37095e-08, desiredTolerance)

exit(numFailures > 0)
