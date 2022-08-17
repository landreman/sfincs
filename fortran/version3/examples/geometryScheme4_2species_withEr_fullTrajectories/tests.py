#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

#execfile('../testsCommon.py')
exec(open('../testsCommon.py').read())

desiredTolerance = 0.001

numFailures = 0

# Species 1
numFailures += shouldBe("FSABFlow[0,0;;;]", 3.97594692502329251E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[0,0;;;]", -6.09853558592433280E-008, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[0,0;;;]",  -2.17831731531552077E-007, desiredTolerance)

# Species 2
numFailures += shouldBe("FSABFlow[1,0;;;]", 3.77579981496331714E-003, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat[1,0;;;]",  -4.23844614210280866E-009, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat[1,0;;;]",  -9.23322435605986916E-009, desiredTolerance)

# before integral outputs
numFailures += shouldBe("NTVBeforeSurfaceIntegral[0,1,0,0;;;]", 1.97402e-06, desiredTolerance)
numFailures += shouldBe("NTVBeforeSurfaceIntegral[1,10,1,0;;;]", 4.57305e-07, desiredTolerance)
numFailures += shouldBe("heatFluxBeforeSurfaceIntegral_vE[1,10,1,0;;;]", 0.0, desiredTolerance)
numFailures += shouldBe("heatFluxBeforeSurfaceIntegral_vE[0,6,0,0;;;]", 0.0, desiredTolerance)
numFailures += shouldBe("heatFluxBeforeSurfaceIntegral_vm[0,4,1,0;;;]", 6.88662e-06, desiredTolerance)
numFailures += shouldBe("heatFluxBeforeSurfaceIntegral_vm[1,8,0,0;;;]", 3.17866e-06, desiredTolerance)
numFailures += shouldBe("heatFluxBeforeSurfaceIntegral_vm0[0,9,1,0;;;]", -6.88738e-06, desiredTolerance)
numFailures += shouldBe("heatFluxBeforeSurfaceIntegral_vm0[1,3,0,0;;;]", -7.79995e-06, desiredTolerance)
numFailures += shouldBe("momentumFluxBeforeSurfaceIntegral_vE[1,6,1,0;;;]", 0.0, desiredTolerance)
numFailures += shouldBe("momentumFluxBeforeSurfaceIntegral_vE[0,4,0,0;;;]", 0.0, desiredTolerance)
numFailures += shouldBe("momentumFluxBeforeSurfaceIntegral_vm[0,7,1,0;;;]", -1.45942e-11, desiredTolerance)
numFailures += shouldBe("momentumFluxBeforeSurfaceIntegral_vm[1,9,0,0;;;]", 4.51749e-08, desiredTolerance)
numFailures += shouldBe("momentumFluxBeforeSurfaceIntegral_vm0[1,2,0,0;;;]", 0.0, desiredTolerance)
numFailures += shouldBe("particleFluxBeforeSurfaceIntegral_vE[0,3,1,0;;;]", 0.0, desiredTolerance)
numFailures += shouldBe("particleFluxBeforeSurfaceIntegral_vE0[1,12,0,0;;;]", 0.0, desiredTolerance)
numFailures += shouldBe("particleFluxBeforeSurfaceIntegral_vm[0,12,1,0;;;]", -3.43298e-06, desiredTolerance)
numFailures += shouldBe("particleFluxBeforeSurfaceIntegral_vm[1,3,0,0;;;]", -6.24729e-06, desiredTolerance)
numFailures += shouldBe("particleFluxBeforeSurfaceIntegral_vm0[0,11,1,0;;;]", -5.72834e-06, desiredTolerance)
numFailures += shouldBe("particleFluxBeforeSurfaceIntegral_vm0[0,8,0,0;;;]", 3.64947e-06, desiredTolerance)

# vs x outputs
numFailures += shouldBe("FSABFlow_vs_x[4,0,0;;;]", 0.000461825, desiredTolerance)
numFailures += shouldBe("FSABFlow_vs_x[0,1,0;;;]", 1.49619e-07, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat_vs_x[2,0,0;;;]", 1.19665e-08, desiredTolerance)
numFailures += shouldBe("heatFlux_vm_psiHat_vs_x[3,1,0;;;]", -3.83757e-09, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat_vs_x[0,0,0;;;]", 3.05541e-12, desiredTolerance)
numFailures += shouldBe("particleFlux_vm_psiHat_vs_x[2,1,0;;;]", -3.69496e-10, desiredTolerance)





exit(numFailures > 0)
