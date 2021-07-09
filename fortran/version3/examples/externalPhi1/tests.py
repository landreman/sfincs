#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

exec(compile(open('../testsCommon.py', "rb").read(), '../testsCommon.py', 'exec'))

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 0.0269953, desiredTolerance)
numFailures += shouldBe("particleFlux_vd_psiHat[0,0;;;]", -6.69433e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[0,0;;;]", -2.76758e-05, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", 0.00826827, desiredTolerance) 
numFailures += shouldBe("particleFlux_vd_psiHat[1,0;;;]", -3.18062e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[1,0;;;]", -1.88472e-05, desiredTolerance)


# Species 3
numFailures += shouldBe("FSABFlow[2,0;;;]", 2.58643e-07, desiredTolerance) 
numFailures += shouldBe("particleFlux_vd_psiHat[2,0;;;]", 3.5988e-12, desiredTolerance)
numFailures += shouldBe("heatFlux_vd_psiHat[2,0;;;]", 6.82697e-11, desiredTolerance)


exit(numFailures > 0)
