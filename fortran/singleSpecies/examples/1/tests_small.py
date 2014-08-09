#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("transportMatrix", 0, -0.0113502, desiredTolerance)
numFailures += shouldBe("transportMatrix", 1, -0.0414953, desiredTolerance)
numFailures += shouldBe("transportMatrix", 2, 0.0276516, desiredTolerance)
numFailures += shouldBe("transportMatrix", 3,  -0.041501, desiredTolerance)
numFailures += shouldBe("transportMatrix", 4,  -0.314625, desiredTolerance)
numFailures += shouldBe("transportMatrix", 5,  -0.0239723, desiredTolerance)
numFailures += shouldBe("transportMatrix", 6,  0.0276755, desiredTolerance)
numFailures += shouldBe("transportMatrix", 7,  -0.0238248, desiredTolerance)
numFailures += shouldBe("transportMatrix", 8,  25.9315, desiredTolerance)

exit(numFailures > 0)
