! Main program

#include <finclude/petscsysdef.h>
#include <petscversion.h>

! For PETSc versions prior to 3.4, the PetscTime subroutine was called PetscGetTime.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 4))
#define PetscTime PetscGetTime
#endif
!Hereafter in this code, use PetscTime.

program sfincs

  use globalVariables
  use printToStdout
  use writeHDF5Output
  use readInput
  use scan
  use petscsysdef

  implicit none

  PetscErrorCode ierr
  PetscLogDouble :: startTime, time1

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)
  masterProc = (myRank==0)

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
  call allocateArraysForSingleRun()

  ! Do various calculations that will not need to be repeated at each
  ! iteration, such as setting up the coordinate grids and evaluating
  ! the magnetic field and its derivatives on the spatial grid.
  call createGrids()

  ! Solve the main system, either linear or nonlinear.
  ! This step takes more time than everything else combined.
  ! Diagnostics should be computed within the solver, for 2 reasons:
  !   1. There might be >1 RHS
  !   2. If doing a nonlinear run, we should also save linear results, which we get for free.
  call solveSystem()

  ! Build the HDF5 output file:
  call writeOutputFile(1)

  call PetscFinalize(ierr)

end program sfincs
