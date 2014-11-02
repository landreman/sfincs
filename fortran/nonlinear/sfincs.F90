! Main program

#include <finclude/petscsysdef.h>
#include "PETScVersions.F90"

program sfincs

  use globalVariables
  use writeHDF5Output
  use readInput
  use petscsysdef
  use solver

  implicit none

  PetscErrorCode ierr
  PetscLogDouble :: startTime, time1

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)
  masterProc = (myRank==0)

  ! In the future, if we want to divide the processors into sub-communicators, this next line would change:
  MPIComm = PETSC_COMM_WORLD

  if (masterProc) then
     print *,"****************************************************************************"
     print *,"SFINCS: Stellarator Fokker-Plank Iterative Neoclassical Conservative Solver"
     print *,"Grid in theta and zeta, Legendre modal in xi, polynomial spectral collocation in x."
     print *,"Nonlinear version."
#if defined(PETSC_USE_REAL_SINGLE)
     print *,"Using single precision."
#else
     print *,"Using double precision."
#endif
     if (numProcs==1) then
        print *,"Serial job (1 process) detected."
     else
        print "(a, i4, a)", " Parallel job (",numProcs," processes) detected."
     end if
  end if

  call PetscTime(time1, ierr)
  startTime = time1

  call readNamelistInput()

  ! If running with >1 proc,
  ! make sure either superlu_dist or mumps is installed, and pick which one
  ! of these packages to use:
  call chooseParallelDirectSolver()

  if (masterProc) then
     print *,"---- Physics parameters: ----"
     print *,"Number of particle species = ", Nspecies
     print *,"Delta (rho* at reference parameters)          = ", Delta
     print *,"alpha (e Phi / T at reference parameters)     = ", alpha
     print *,"nu_n (collisionality at reference parameters) = ", nu_n
     print *,"dPhiHatdpsiN (radial electric field)          = ", dPhiHatdpsiN
     if (nonlinear) then
        print *,"Nonlinear run"
     else
        print *,"Linear run"
     end if
  end if

  ! Before investing the time in solving the system, make sure
  ! it is possible to at least open the output file.
  call openOutputFile()

  ! Change this next subroutine?
  !call allocateArraysForSingleRun()

  ! Do various calculations that will not need to be repeated at each
  ! iteration, such as setting up the coordinate grids and evaluating
  ! the magnetic field and its derivatives on the spatial grid.
  call createGrids()

  ! Solve the main system, either linear or nonlinear.
  ! This step takes more time than everything else combined.
  ! Diagnostics should be computed within the solver, for 2 reasons:
  !   1. There might be >1 RHS
  !   2. If doing a nonlinear run, we should also save linear results, which we get for free.
  call mainSolverLoop()

  ! Build the HDF5 output file:
  call writeOutputFile()

  call PetscFinalize(ierr)

end program sfincs
