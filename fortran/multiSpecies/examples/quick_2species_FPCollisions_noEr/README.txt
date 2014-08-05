This example is a very low resolution case in the model W7-X field (geometryScheme=4) which should be able to fit on a single processor. This may be a good example to try first.


Expected results:
----------------------------------------------------
Species 1:
 FSABFlow:                 -9.46E-004
 particleFlux:             -8.83E-009
 heatFlux:                 -8.42E-009

Species 2:
 FSABFlow:                 -1.20E-005
 particleFlux:              9.01E-011
 heatFlux:                  2.01E-011

FSABjHat (bootstrap current):   -1.02E-003

Don't pay much attention to FSADensityPerturbation, FSAPressurePerturbation, and momentumFlux. These quantities converge to 0, so the precise value you obtain is extremely sensitive to roundoff error etc and may not be the same from run to run.
