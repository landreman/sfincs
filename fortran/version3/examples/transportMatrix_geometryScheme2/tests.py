#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

# The values below are for Nxi=14, which it turns out was a bit too low for convergence to 3%:
#numFailures += shouldBe("transportMatrix[0,0;;;]", -0.0113502, desiredTolerance)
#numFailures += shouldBe("transportMatrix[0,1;;;]",  -0.0414953, desiredTolerance)
#numFailures += shouldBe("transportMatrix[0,2;;;]",  0.0276516, desiredTolerance)
#numFailures += shouldBe("transportMatrix[1,0;;;]",   -0.041501, desiredTolerance)
#numFailures += shouldBe("transportMatrix[1,1;;;]",   -0.314625, desiredTolerance)
##numFailures += shouldBe("transportMatrix[1,2;;;]",   -0.0239723, desiredTolerance) # Value for Nxi_for_x_option=0
#numFailures += shouldBe("transportMatrix[1,2;;;]",   -0.0238869, desiredTolerance)
#numFailures += shouldBe("transportMatrix[2,0;;;]",   0.0276755, desiredTolerance)
##numFailures += shouldBe("transportMatrix[2,1;;;]",   -0.0238248, desiredTolerance) # Value for Nxi_for_x_option=0
#numFailures += shouldBe("transportMatrix[2,1;;;]",   -0.0237431, desiredTolerance)
#numFailures += shouldBe("transportMatrix[2,2;;;]",   25.9315, desiredTolerance)

# The values below are for Nxi=16, in which case the outputs are converged to 3%:
numFailures += shouldBe("transportMatrix[0,0;;;]", -0.0113763, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,1;;;]", -0.0416693, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,2;;;]",  0.0276852, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,0;;;]", -0.0416749, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,1;;;]", -0.315883, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,2;;;]", -0.0231003, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,0;;;]",  0.0277078, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,1;;;]", -0.0229643, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,2;;;]", 25.8266, desiredTolerance)

exit(numFailures > 0)
