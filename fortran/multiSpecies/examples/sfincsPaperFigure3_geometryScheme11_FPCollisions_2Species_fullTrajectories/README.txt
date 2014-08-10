The conditions for this example are those of figure 3 in the paper on SFINCS:
Landreman, Smith, Mollen, and Helander, Physics of Plasmas 21, 042503 (2014).

Instead of a full scan of E_r as in the figure, here we consider just
a single value for E_r = -9 kV/m, corresponding to dPhiHatdpsiN = 2.5.

On Apr 18, 2014, this job took ~3 minutes on Edison at NERSC.

Expected results
----------------------------------------------------
Species 1:
 FSADensityPerturbation:    < 1E-010 (Precise value doesn't matter.)
 FSABFlow:                  9.19E-003
 FSAPressurePerturbation:   < 1E-011 (Precise value doesn't matter.)
 particleFlux:              -1.08E-006
 momentumFlux:              < 1E-008 (Precise value doesn't matter.)
 heatFlux:                  -2.34E-006
Species 2:
 FSADensityPerturbation:    < 1E-012 (Precise value doesn't matter.)
 FSABFlow:                  1.44E-003
 FSAPressurePerturbation:   < 1E-012 (Precise value doesn't matter.)
 particleFlux:              -5.36E-008
 momentumFlux:              < 1E-011 (Precise value doesn't matter.)
 heatFlux:                  -1.08E-007
FSABjHat (bootstrap current):  7.75E-003
