#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("transportMatrix[0,0;;;]", -2.84548264465028485E-003, desiredTolerance)
numFailures += shouldBe("transportMatrix[0,1;;;]", -1.06835854356196296E-003, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,0;;;]", -1.07164063354908494E-003, desiredTolerance)
numFailures += shouldBe("transportMatrix[1,1;;;]", 0.88506114839264771, desiredTolerance)

exit(numFailures > 0)
