! Main program

#include <finclude/petscsysdef.h>
#include "PETScVersions.F90"

program sfincs

  use globalVariables
  use writeHDF5Output
  use readInput
  use petscsysdef
  use solver
  use geometry
  use radialCoordinates

  implicit none

  PetscErrorCode ierr
  PetscLogDouble :: startTime, time1

!!$  ! Begin testing Chebyshev stuff
!!$  integer :: i
!!$  PetscScalar :: xMin
!!$  PetscScalar, allocatable :: interp(:,:)
!!$
!!$  Nx=5
!!$  xMin = -5.3d+0
!!$  xMax = -2.1d+0
!!$  NxPotentials = 5
!!$  allocate(ddx(Nx,Nx))
!!$  allocate(x(Nx))
!!$  allocate(xWeights(Nx))
!!$  allocate(xPotentials(NxPotentials))
!!$  allocate(interp(NxPotentials,Nx))
!!$  xPotentials = 0.0d+0
!!$  xPotentials(1) = -2.3d+0
!!$  xPotentials(2) = -3.0d+0
!!$  xPotentials(3) = -4.4d+0
!!$  xPotentials(4) = -2.1d+0
!!$  xPotentials(5) = -3.1d+0
!!$  call ChebyshevGrid(Nx,xMin,xMax,x,xWeights,ddx)
!!$  call ChebyshevInterpolationMatrix(Nx,NxPotentials,x,xPotentials,interp)
!!$  print *,"x:",x
!!$  print *,"xWeights:",xWeights
!!$  print *,"Here comes ddx:"
!!$  do i=1,Nx
!!$     print *,ddx(i,:)
!!$  end do
!!$  print *,"Here comes interp:"
!!$  do i=1,NxPotentials
!!$     print *,interp(i,:)
!!$  end do
!!$  stop
!!$  ! End testing Chebyshev stuff

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)
  masterProc = (myRank==0)

  ! In the future, if we want to divide the processors into sub-communicators, this next line would change:
  MPIComm = PETSC_COMM_WORLD

  if (masterProc) then
     print *,"****************************************************************************"
     print *,"SFINCS: Stellarator Fokker-Plank Iterative Neoclassical Conservative Solver"
     print *,"Version 3"
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
  call validateInput()

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
     if (nonlinear) then
        print *,"Nonlinear run"
     else
        print *,"Linear run"
     end if
  end if

  ! Initialize NPeriods, psiAHat, and aHat.  We need to know NPeriods before
  ! we can initialize the zeta grid.
  call initializeGeometry()

  ! Do various calculations that will not need to be repeated at each
  ! iteration, such as setting up the coordinate grids and evaluating
  ! the magnetic field and its derivatives on the spatial grid.
  call createGrids()

  if (RHSMode==3) then
     ! Monoenergetic coefficient computation.
     ! Overwrite nu_n and dPhiHatd* using nuPrime and EStar.

     nu_n = nuPrime * B0OverBBar / (GHat + iota * IHat)
     dPhiHatdpsiHat = 2 / (alpha * Delta) * EStar * iota * B0OverBBar / GHat
  end if

  ! For input quantities that depend on the radial coordinate, pick out the values for the selected
  ! radial coordinate, and use these values to over-write values for the other radial coordinates.
  call setInputRadialCoordinate()

  ! Create HDF5 data structures, and save the quantities that will not change
  ! at each iteration of the solver (i.e. save all quantities except diagnostics.)
  call initializeOutputFile()

  ! Solve the main system, either linear or nonlinear.
  ! This step takes more time than everything else combined.
  call mainSolverLoop()

  call finalizeHDF5()
  call PetscFinalize(ierr)

  if (masterProc) then
     print *,"Goodbye!"
  end if

end program sfincs
