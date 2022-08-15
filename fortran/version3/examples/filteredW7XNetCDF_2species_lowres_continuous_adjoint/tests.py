#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("dPhidPsidLambda[0,0,0;;;]",1.56229, desiredTolerance)
numFailures += shouldBe("dPhidPsidLambda[0,1,0;;;]",-0.0193202, desiredTolerance)
numFailures += shouldBe("dPhidPsidLambda[0,2,0;;;]",-0.0401214, desiredTolerance)
numFailures += shouldBe("dPhidPsidLambda[0,3,0;;;]",0.54171, desiredTolerance)
numFailures += shouldBe("dPhidPsidLambda[1,0,0;;;]",-50.0455, desiredTolerance)
numFailures += shouldBe("dPhidPsidLambda[2,0,0;;;]",-100.937, desiredTolerance)
numFailures += shouldBe("dPhidPsidLambda[3,0,0;;;]",-27.8693, desiredTolerance)
numFailures += shouldBe("dPhidPsidLambda[4,0,0;;;]",30.5579, desiredTolerance)


# Species 1
numFailures += shouldBe("dHeatFluxdLambda[0,1,0,0;;;]", 6.08663e-09, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,2,0,0;;;]", -7.07285e-06, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,3,0,0;;;]", 2.54599e-09, desiredTolerance)

# Species 2
numFailures += shouldBe("dParallelFlowdLambda[0,0,1,0;;;]", 0.000514809, desiredTolerance)

exit(numFailures > 0)
