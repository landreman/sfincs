#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

# Check the first species:
numFailures += shouldBe("FSABFlow[0,0;;;]", 308.77570247458709, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -1.19377717309288199E-004, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]", -5.46764373674423623E-005, desiredTolerance)

# Check the second species:
numFailures += shouldBe("FSABFlow[1,0;;;]", 4.0387888836037282, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]", -1.48789578198661934E-006, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]", -1.35574234718297328E-006, desiredTolerance)

numFailures += shouldBe("FSABjHat[0;;;]", 333.00843577620947, desiredTolerance)

exit(numFailures > 0)
