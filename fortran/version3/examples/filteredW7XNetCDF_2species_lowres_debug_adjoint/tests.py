#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0



# Species 1
numFailures += shouldBe("dHeatFluxdLambda[0,0,0,0;;;]", -2.91416e-07, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,1,0,0;;;]", 1.45988e-08, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,2,0,0;;;]", -1.63331e-08, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,3,0,0;;;]", 4.79669e-09, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[4,0,0,0;;;]", -2.74338e-06, desiredTolerance)

numFailures += shouldBe("dParallelFlowdLambda[0,0,0,0;;;]", -0.000397396, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,1,0,0;;;]", 0.0046219, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,2,0,0;;;]", -8.46991e-05, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,3,0,0;;;]", -0.00975908, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[3,0,0,0;;;]", -0.0520897, desiredTolerance)


numFailures += shouldBe("dParticleFluxdLambda[0,0,0,0;;;]", -7.09016e-08, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,1,0,0;;;]", 3.74166e-09, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,2,0,0;;;]", -4.25987e-09, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,3,0,0;;;]", 3.95029e-08, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[2,0,0,0;;;]", 1.5141e-07, desiredTolerance)

# Species 2
numFailures += shouldBe("dHeatFluxdLambda[0,0,1,0;;;]", -9.26131e-09, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,1,1,0;;;]", 5.45938e-10, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,2,1,0;;;]", -5.11362e-10, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,3,1,0;;;]", 4.79669e-09, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[4,0,1,0;;;]", -9.29347e-08, desiredTolerance)

numFailures += shouldBe("dParallelFlowdLambda[0,0,1,0;;;]", -0.000332624, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,1,1,0;;;]", 0.00170519, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,2,1,0;;;]", -4.18839e-05, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,3,1,0;;;]", -0.00334459, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[3,0,1,0;;;]", -0.0279159, desiredTolerance)

numFailures += shouldBe("dParticleFluxdLambda[0,0,1,0;;;]", -2.18328e-09, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,1,1,0;;;]", 1.66491e-10, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,2,1,0;;;]", -1.27847e-10, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,3,1,0;;;]", 1.12286e-09, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[2,0,1,0;;;]", 1.31782e-08, desiredTolerance)

exit(numFailures > 0)
