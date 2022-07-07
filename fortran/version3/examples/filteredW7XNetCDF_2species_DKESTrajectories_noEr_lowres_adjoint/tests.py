#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

exec(open('../testsCommon.py').read())


desiredTolerance = 0.001

numFailures = 0



# Species 1
numFailures += shouldBe("dHeatFluxdLambda[0,0,0,0;;;]", -2.4742e-07, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,1,0,0;;;]", 1.22259e-08, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,2,0,0;;;]", -1.37189e-08, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,3,0,0;;;]", 1.35085e-07, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[4,0,0,0;;;]", -2.42081e-06, desiredTolerance)

numFailures += shouldBe("dParallelFlowdLambda[0,0,0,0;;;]", -4.32542e-05, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,1,0,0;;;]", 0.00345475, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,2,0,0;;;]", -6.38307e-05, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,3,0,0;;;]", -0.00737457, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[3,0,0,0;;;]", -0.0386586, desiredTolerance)


numFailures += shouldBe("dParticleFluxdLambda[0,0,0,0;;;]", -5.79136e-08, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,1,0,0;;;]", 2.95706e-09, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,2,0,0;;;]", -3.36903e-09, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,3,0,0;;;]", 3.16934e-08, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[2,0,0,0;;;]", 2.87945e-07, desiredTolerance)

# Species 2
numFailures += shouldBe("dHeatFluxdLambda[0,0,1,0;;;]", -1.07825e-08, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,1,1,0;;;]", 6.2772e-10, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,2,1,0;;;]", -6.03516e-10, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[0,3,1,0;;;]", 5.623e-09, desiredTolerance)
numFailures += shouldBe("dHeatFluxdLambda[4,0,1,0;;;]", -1.06544e-07, desiredTolerance)

numFailures += shouldBe("dParallelFlowdLambda[0,0,1,0;;;]", -9.47518e-05, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,1,1,0;;;]", 0.000871583, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,2,1,0;;;]", -2.52776e-05, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[0,3,1,0;;;]", -0.00168835, desiredTolerance)
numFailures += shouldBe("dParallelFlowdLambda[3,0,1,0;;;]", -0.0205893, desiredTolerance)

numFailures += shouldBe("dParticleFluxdLambda[0,0,1,0;;;]", -2.6892e-09, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,1,1,0;;;]", 1.99397e-10, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,2,1,0;;;]", -1.63386e-10, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[0,3,1,0;;;]", 1.42444e-09, desiredTolerance)
numFailures += shouldBe("dParticleFluxdLambda[2,0,1,0;;;]", 9.07324e-09, desiredTolerance)

exit(numFailures > 0)
