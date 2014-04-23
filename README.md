sfincs
======

SFINCS: the Stellarator Fokker-Planck Iterative Neoclassical Conservative Solver.

SFINCS is a code designed to compute neoclassical effects in nonaxisymmetric toroidal plasmas, such as stellarators and perturbed tokamaks.  The code solves a linear drift-kinetic equation for the distribution function.  In addition to the neoclassical fluxes, flows, and bootstrap current, one can also obtain other moments such as the density variation on a flux surface.

Documentation is available in several places. The primary published reference is [Landreman, Smith, Mollen, and Helander, Physics of Plasmas 21, 042503 (2014)](https://github.com/landreman/sfincs/blob/master/doc/LandremanSmithMollenHelander_2014_PoP_v21_p042503_SFINCS.pdf?raw=true). Further details of the velocity-space discretization and Fokker-Planck collision operator may be found in [Landreman and Ernst, Journal of Computational Physics 243, 130 (2013)](). Also see the [technical documentation](https://github.com/landreman/sfincs/blob/master/doc/20131219-01 Technical documentation for SFINCS with multiple species.pdf?raw=true) and other files in the `/doc/` folder of the repository, as well as [the wiki.](https://github.com/landreman/sfincs/wiki) 

SFINCS is a descendant of the axisymmetric codes [here](https://github.com/landreman/tokamakDriftKineticEquationSolver).  SFINCS is also closely related to the radially global tokamak code [PERFECT](https://github.com/landreman/perfect).  All of these codes have two independent velocity-space variables: speed and pitch angle.  SFINCS and PERFECT both have two independent spatial variables: poloidal and toroidal angles in SFINCS, and poloidal and radial coordinates in PERFECT.

Versions of SFINCS in both matlab (serial) and fortran (parallel) are available in this repository. The matlab and fortran versions are entirely independent of each other, and they have nearly identical functionality and variable names.  Finer resolution can be used in the fortan version because it is parallelized and therefore able to access more memory.  To accurately resolve the distribution function for experimentally relevant magnetic fields and collisionality, you will likely need the fortran version to access enough memory.

For historical reasons, both single-species and multi-species versions of SFINCS are available here. The single-species version is useful for computing the ion transport matrix. For other applications, the multi-species version is recommended.

The fortran versions rely on the [PETSc library](http://www.mcs.anl.gov/petsc/) and use the HDF5 library for output. The makefiles for the fortran versions of SFINCS may need to be adapted to link to these libraries on your computing system.

For inquiries, contact Matt Landreman at matt dot landreman at gmail dot com.
