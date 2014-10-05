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
  integer :: runNum

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)
  ! This program uses 1-based indices to number the MPI processes, consistent with the Fortran
  ! convention for array indices, not the 0-based indexing used by C and MPI.
  myRank = myRank + 1
  masterProc = (myRank==1)

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

  call chooseParallelDirectSolver()

  call openOutputFile()

  call allocateArraysForSingleRun()

  call createGrids()

  call solveDKE()

  call writeRunToOutputFile(1)

  call PetscFinalize(ierr)

end program sfincs
