#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("transportMatrix[0,0;;;]", -0.000630586, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,1;;;]", -0.00199776, desiredTolerance)

# This next line changed sign with Hakan's sign corrections:
numFailures += shouldBe("transportMatrix[0,2;;;]", -0.0620745, desiredTolerance)

numFailures += shouldBe("transportMatrix[1,0;;;]", -0.00199966, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,1;;;]", -0.0107542, desiredTolerance)

# These next 3 lines changed sign with Hakan's sign corrections:
numFailures += shouldBe("transportMatrix[1,2;;;]", -0.171768, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,0;;;]", -0.0621547, desiredTolerance)
numFailures += shouldBe("transportMatrix[2,1;;;]", -0.172102, desiredTolerance)

numFailures += shouldBe("transportMatrix[2,2;;;]", 35.0602, desiredTolerance)

exit(numFailures > 0)
