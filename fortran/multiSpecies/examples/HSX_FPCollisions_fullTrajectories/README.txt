Example:
HSX, 2 species (hydrogen + electrons), r/a=0.22
E_r included. This run is not at the ambipolar E_r.

Normalizations are:
nBar = 10^18 / m^3
TBar = 1 eV
BBar = 1 T
RBar = 1 m
mBar = proton mass

Expected results:
----------------------------------------------------
Species 1:
 FSABFlow:                  -15.5
 particleFlux:             -2.15E-006
 heatFlux:                 -1.80E-004

Species 2:
 FSABFlow:                   36.5
 particleFlux:             -2.08E-005
 heatFlux:                 -4.45E-002

FSABjHat (bootstrap current):   -51.98

Don't pay much attention to FSADensityPerturbation, FSAPressurePerturbation, and momentumFlux. These quantities converge to 0, so the precise value you obtain is extremely sensitive to roundoff error etc and may not be the same from run to run.
