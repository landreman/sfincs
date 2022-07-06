#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.002

numFailures = 0

# Check the first species:
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.770757651696965, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", 1.161190184424680E-005, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", 8.819637499826962E-004, desiredTolerance)

# Check the second species:
#numFailures += shouldBe("FSABFlow[1,0;;;]", -58.1306611699325, desiredTolerance) ##Commented by AM 2018-11
numFailures += shouldBe("FSABFlow[1,0;;;]", -58.1306611699325, 10.0*desiredTolerance) ##Added by AM 2018-11
#numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", 2.098719213759062E-005, desiredTolerance) # Value for Nxi_for_x_option=0
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", 2.07956e-05, desiredTolerance)
#numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",4.482791404112483E-002, desiredTolerance) # Value for Nxi_for_x_option=0
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",0.0446038, desiredTolerance)

#numFailures += shouldBe("FSABjHat[0;;;]", 58.9014188216295, desiredTolerance) ##Commented by AM 2018-11
numFailures += shouldBe("FSABjHat[0;;;]", 58.9014188216295, 10.0*desiredTolerance) ##Added by AM 2018-11

exit(numFailures > 0)
