module scan

  use globalVariables
  use petscsysdef

  implicit none

#include <finclude/petscsysdef.h>

  integer, allocatable :: NthetasForScan(:), NzetasForScan(:), NxisForScan(:), &
       NLsForScan(:), NxsForScan(:), constraintSchemesForScan(:)
  PetscScalar, allocatable :: NxPotentialsPerVthsForScan(:), xMaxsForScan(:)
  PetscScalar, allocatable :: solverTolerancesForScan(:)
  integer :: numRunsInScan, minScanUnit, maxScanUnit

contains

  ! ----------------------------------------------------------

  subroutine allocateArraysForSingleRun()

    implicit none

    call setConstraintScheme()

    allocate(NthetasForScan(1))
    allocate(NzetasForScan(1))
    allocate(NxisForScan(1))
    allocate(NLsForScan(1))
    allocate(NxsForScan(1))
    allocate(NxPotentialsPerVthsForScan(1))
    allocate(xMaxsForScan(1))
    allocate(constraintSchemesForScan(1))

    NthetasForScan(1) = Ntheta
    NzetasForScan(1) = Nzeta
    NxisForScan(1) = Nxi
    NLsForScan(1) = NL
    NxsForScan(1) = Nx
    NxPotentialsPerVthsForScan(1) = NxPotentialsPerVth
    xMaxsForScan(1) = xMax
    constraintSchemesForScan(1) = constraintScheme

    numRunsInScan = 1

  end subroutine allocateArraysForSingleRun

  ! ----------------------------------------------------------

  subroutine processConvergenceScanParameters()

    implicit none

    integer, allocatable :: Nthetas(:), Nzetas(:), Nxis(:), NLs(:), Nxs(:)
    PetscScalar, allocatable :: NxPotentialsPerVths(:), xMaxs(:), tempArray(:)
    PetscScalar, allocatable :: solverTolerances(:)
    integer :: i, j, k, currentIndex
    PetscScalar :: verySmall = 1d-5


    call setConstraintScheme()

    allocate(tempArray(NthetaNumRuns))
    allocate(NThetas(NthetaNumRuns))
    call logspace(Ntheta*NthetaMinFactor, Ntheta*NthetaMaxFactor, NthetaNumRuns, tempArray)
    Nthetas = nint(tempArray)
    deallocate(tempArray)

    allocate(tempArray(NzetaNumRuns))
    allocate(Nzetas(NzetaNumRuns))
    call logspace(Nzeta*NzetaMinFactor, Nzeta*NzetaMaxFactor, NzetaNumRuns, tempArray)
    Nzetas = nint(tempArray)
    deallocate(tempArray)

    allocate(tempArray(NxiNumRuns))
    allocate(Nxis(NxiNumRuns))
    call logspace(Nxi*NxiMinFactor, Nxi*NxiMaxFactor, NxiNumRuns, tempArray)
    Nxis = nint(tempArray)
    deallocate(tempArray)

    allocate(tempArray(NLNumRuns))
    allocate(NLs(NLNumRuns))
    call logspace(NL*NLMinFactor, NL*NLMaxFactor, NLNumRuns, tempArray)
    NLs = nint(tempArray)
    deallocate(tempArray)

    allocate(tempArray(NxNumRuns))
    allocate(Nxs(NxNumRuns))
    call logspace(Nx*NxMinFactor, Nx*NxMaxFactor, NxNumRuns, tempArray)
    Nxs = nint(tempArray)
    deallocate(tempArray)

    allocate(NxPotentialsPerVths(NxPotentialsPerVthNumRuns))
    call logspace(NxPotentialsPerVth*NxMinFactor, NxPotentialsPerVth*NxMaxFactor, &
         NxPotentialsPerVthNumRuns, NxPotentialsPerVths)

    allocate(xMaxs(xMaxNumRuns))
    call logspace(xMax*xMaxMinFactor, xMax*xMaxMaxFactor, xMaxNumRuns, xMaxs)

    allocate(solverTolerances(solverToleranceNumRuns))
    call logspace(solverTolerance*solverToleranceMinFactor, solverTolerance*solverToleranceMaxFactor, &
         solverToleranceNumRuns, solverTolerances)

    numRunsInScan = 1 + NthetaNumRuns + NzetaNumRuns + NxiNumRuns &
         + NLNumRuns + NxNumRuns + NxPotentialsPerVthNumRuns + xMaxNumRuns + solverToleranceNumRuns

    allocate(NthetasForScan(numRunsInScan))
    allocate(NzetasForScan(numRunsInScan))
    allocate(NxisForScan(numRunsInScan))
    allocate(NLsForScan(numRunsInScan))
    allocate(NxsForScan(numRunsInScan))
    allocate(NxPotentialsPerVthsForScan(numRunsInScan))
    allocate(xMaxsForScan(numRunsInScan))
    allocate(solverTolerancesForScan(numRunsInScan))
    
    NthetasForScan = Ntheta
    NzetasForScan = Nzeta
    NxisForScan = Nxi
    NLsForScan = NL
    NxsForScan = Nx
    NxPotentialsPerVthsForScan = NxPotentialsPerVth
    xMaxsForScan = xMax
    solverTolerancesForScan = solverTolerance

    currentIndex=2

    do i=1, NThetaNumRuns
       NthetasForScan(currentIndex) = Nthetas(i)
       currentIndex = currentIndex + 1
    end do

    do i=1, NzetaNumRuns
       NzetasForScan(currentIndex) = Nzetas(i)
       currentIndex = currentIndex + 1
    end do

    do i=1, NxiNumRuns
       NxisForScan(currentIndex) = Nxis(i)
       currentIndex = currentIndex + 1
    end do

    do i=1, NLNumRuns
       NLsForScan(currentIndex) = NLs(i)
       currentIndex = currentIndex + 1
    end do

    do i=1, NxNumRuns
       NxsForScan(currentIndex) = Nxs(i)
       currentIndex = currentIndex + 1
    end do

    do i=1, NxPotentialsPerVthNumRuns
       NxPotentialsPerVthsForScan(currentIndex) = NxPotentialsPerVths(i)
       currentIndex = currentIndex + 1
    end do

    do i=1, xMaxNumRuns
       xMaxsForScan(currentIndex) = xMaxs(i)
       currentIndex = currentIndex + 1
    end do

    do i=1, solverToleranceNumRuns
       solverTolerancesForScan(currentIndex) = solverTolerances(i)
       currentIndex = currentIndex + 1
    end do

    if (currentIndex /= numRunsInScan+1) then
       print *,"Error - something went wrong:"
       print *,"  currentIndex:",currentIndex
       print *,"  numRunsInScan:",numRunsInScan
       stop
    end if

    if (forceOddNthetaAndNzeta) then
       do i=1,numRunsInScan
          if (mod(NthetasForScan(i), 2) == 0) then
             NthetasForScan(i) = NthetasForScan(i) + 1
          end if
          if (mod(NzetasForScan(i), 2) == 0) then
             NzetasForScan(i) = NzetasForScan(i) + 1
          end if
       end do
    end if

    i=1
    do
       if (i .ge. numRunsInScan) then
          exit
       end if
       j = i + 1
       do
          if (j > numRunsInScan) then
             exit
          end if
          if (  &
               (abs(NthetasForScan(i)-NthetasForScan(j)) < verySmall) .and. &
               (abs(NzetasForScan(i)-NzetasForScan(j)) < verySmall) .and. &
               (abs(NxisForScan(i)-NxisForScan(j)) < verySmall) .and. &
               (abs(NLsForScan(i)-NLsForScan(j)) < verySmall) .and. &
               (abs(NxsForScan(i)-NxsForScan(j)) < verySmall) .and. &
               (abs(NxPotentialsPerVthsForScan(i)-NxPotentialsPerVthsForScan(j)) < verySmall) .and. &
               (abs(log(solverTolerancesForScan(i))-log(solverTolerancesForScan(j))) < verySmall) .and. &
               (abs(xMaxsForScan(i)-xMaxsForScan(j)) < verySmall) ) then

             ! Item j is a duplicate, so remove it:
             do k=j+1, numRunsInScan
                NthetasForScan(k-1) = NthetasForScan(k)
                NzetasForScan(k-1) = NzetasForScan(k)
                NxisForScan(k-1) = NxisForScan(k)
                NLsForScan(k-1) = NLsForScan(k)
                NxsForScan(k-1) = NxsForScan(k)
                NxPotentialsPerVthsForScan(k-1) = NxPotentialsPerVthsForScan(k)
                solverTolerancesForScan(k-1) = solverTolerancesForScan(k)
                xMaxsForScan(k-1) = xMaxsForScan(k)
             end do
             numRunsInScan = numRunsInScan - 1
          else
             j = j + 1
          end if
       end do
       i = i + 1
    end do

    allocate(constraintSchemesForScan(numRunsInScan))
    constraintSchemesForScan = constraintScheme

  end subroutine processConvergenceScanParameters

  !---------------------------------------------------------------------------------------------

  subroutine prepareEStarScan()

    implicit none

    numRunsInScan = NEStar

    call setConstraintScheme()

    allocate(NthetasForScan(numRunsInScan))
    allocate(NzetasForScan(numRunsInScan))
    allocate(NxisForScan(numRunsInScan))
    allocate(NLsForScan(numRunsInScan))
    allocate(NxsForScan(numRunsInScan))
    allocate(NxPotentialsPerVthsForScan(numRunsInScan))
    allocate(xMaxsForScan(numRunsInScan))
    allocate(constraintSchemesForScan(numRunsInScan))

    NthetasForScan = Ntheta
    NzetasForScan = Nzeta
    NxisForScan = Nxi
    NLsForScan = NL
    NxsForScan = Nx
    NxPotentialsPerVthsForScan = NxPotentialsPerVth
    xMaxsForScan = xMax
    constraintSchemesForScan = constraintScheme

    allocate(EStars(NEStar))
    call logspace(EStarMin, EStarMax, NEStar, EStars)

  end subroutine prepareEStarScan

  !---------------------------------------------------------------------------------------------

  subroutine linspace(min, max, N, array)
    ! This subroutine acts just like the Matlab function of the same name.

    implicit none

    PetscScalar, intent(in) :: min, max
    integer, intent(in) :: N
    PetscScalar :: array(:)
    integer :: i

    if (N < 0) then
       print *,"Error! 'N' must be at least 0."
       stop
    end if

    do i=1,N
       array(i) = (i-one)*(max-min)/(N-one) + min
    end do

  end subroutine linspace

  !---------------------------------------------------------------------------------------------

  subroutine logspace(min, max, N, array)
    ! NOTE: this function works differently than the MATLAB logspace function
    ! in that it returns values from min to max rather than from 10^min to 10^max.

    implicit none

    PetscScalar, intent(in) :: min, max
    integer, intent(in) :: N
    PetscScalar :: array(:)
    integer :: i

    if (min<0) then
       print *,"Error! 'min' must be >0."
       stop
    end if

    if (max<0) then
       print *,"Error! 'max' must be >0."
       stop
    end if

    call linspace(log(min), log(max), N, array)
    do i=1,N
       array(i) = exp(array(i))
    end do

  end subroutine logspace

  !---------------------------------------------------------------------------------------------

  subroutine setMPICommunicatorsForScan()
    ! This subroutine divides the PETSC_COMM_WORLD MPI communicator into smaller communicators.
    ! Each sub-group of processes can then work independently on different runs which are needed
    ! for a parameter scan.

    implicit none

    integer :: i, j, numProcsForThisComm
    integer :: currentCommunicator, currentMinUnit, firstIndex, lastIndex
    integer, dimension(:), allocatable :: commMaxProcs, procsToInclude
    integer :: convergenceScanUnits
    PetscErrorCode ierr
    MPI_Group :: mpi_group_world
    MPI_Group, dimension(:), allocatable :: mpi_groups_scan
    MPI_Comm, dimension(:), allocatable :: mpi_comms_scan

    ! In this code, the rank of each procedure is considered a 1-based index (consistent
    ! with Fortran convention) rather than a 0-based index (used in MPI calls).

    ! For now, each run required for the parameter scan will be considered a "unit".
    ! In the future, I might want to group several runs into each unit for better efficiency:
    ! some units would consist of 1 or 2 expensive runs, while other units would consist
    ! of many inexpensive runs.
    convergenceScanUnits = numRunsInScan

    if (.not. parallelizeOverScan) then
       numCommunicators=1
       MPIComm = PETSC_COMM_WORLD
       myRankInSubComm = myRank
       numProcsInSubComm = numProcs
       masterProcInSubComm = masterProc
       myCommunicatorIndex = 1
       minScanUnit = 1
       maxScanUnit = numRunsInScan
    else



       numCommunicators = min(numProcs, convergenceScanUnits)

       allocate(minUnits(numProcs))
       allocate(maxUnits(numProcs))
       allocate(procsToInclude(numProcs))
       allocate(mpi_comms_scan(numCommunicators))
       allocate(mpi_groups_scan(numCommunicators))
       allocate(commMinProcs(numCommunicators))
       allocate(commMaxProcs(numCommunicators))

       do i=1,numProcs
          minUnits(i) = floor(convergenceScanUnits * (1.d+0) * (i-1) / numProcs) + 1
          maxUnits(i) = floor(convergenceScanUnits * (1.d+0) * i / numProcs)
          maxUnits(i) = max(maxUnits(i), minUnits(i))
       end do
       minScanUnit = minUnits(myRank)
       maxScanUnit = maxUnits(myRank)

       if (masterProc) then
          do i=1,numProcs
             print "(a, i4, a, i3, a, i3)", "MPI proc ",i," is responsible for scan units ",&
                  minUnits(i)," to ",maxUnits(i)
          end do
       end if

       commMinProcs(1) = 1
       commMaxProcs(numCommunicators) = numProcs
       currentCommunicator = 1
       currentMinUnit = 1
       myCommunicatorIndex = -1
       do i=2,numProcs
          if (minUnits(i) /= currentMinUnit) then
             currentMinUnit = minUnits(i)
             commMaxProcs(currentCommunicator) = i-1
             currentCommunicator = currentCommunicator + 1
             commMinProcs(currentCommunicator) = i
          end if
          if (myRank == i) then
             myCommunicatorIndex = currentCommunicator
          end if
       end do

       if (myRank == 1) then
          myCommunicatorIndex = 1
       end if
       if (myCommunicatorIndex == -1) then
          print "(a, i4, a)","Error! Somehow, myCommunicatorIndex for proc ",myRank," did not get assigned."
          stop
       end if

       if (currentCommunicator /= numCommunicators) then
          if (masterProc) then
             print *,"Error! Something went wrong in assigning processors to communicators."
          end if
          stop
       end if

       if (masterProc) then
          print "(a, i4, a)","Creating ",numCommunicators, &
               " MPI communicators for parallelizing the parameter scan."
          do i=1,numCommunicators
             print "(a, i4, a, i4, a, i4, a, i3, a, i3)", "Communicator ",i," consists of procs ", &
                  commMinProcs(i)," through ",commMaxProcs(i), " and will handle scan units ", &
                  minUnits(commMinProcs(i))," through ",maxUnits(commMinProcs(i))
          end do
       end if

       call MPI_Comm_group(PETSC_COMM_WORLD, mpi_group_world, ierr)

       do i=1,numCommunicators
          numProcsForThisComm = commMaxProcs(i)-commMinProcs(i) + 1
          do j=1,numProcsForThisComm
             procsToInclude(j) = commMinProcs(i) + j - 2
             ! Above, we needed to subtract an additional 1 to convert from Fortran 1-based indices
             ! to MPI 0-based indices.
          end do
          call MPI_Group_incl(mpi_group_world, numProcsForThisComm, procsToInclude(1:numProcsForThisComm), mpi_groups_scan(i), ierr)
          call MPI_Comm_create(PETSC_COMM_WORLD, mpi_groups_scan(i), mpi_comms_scan(i), ierr)
       end do

       ! Next, set the MPI communicator that this process will use for the rest of the program execution:
       MPIComm = mpi_comms_scan(myCommunicatorIndex)
       CHKERRQ(ierr)

       call MPI_COMM_SIZE(MPIComm, numProcsInSubComm, ierr)
       call MPI_COMM_RANK(MPIComm, myRankInSubComm, ierr)
       myRankInSubComm = myRankInSubComm + 1
       masterProcInSubComm = (myRankInSubComm == 1)
    end if
  end subroutine setMPICommunicatorsForScan

  ! ------------------------------------------------------------------------

  subroutine chooseParallelDirectSolver()

    implicit none

    isAParallelDirectSolverInstalled = .false.

    if ((whichParallelSolverToFactorPreconditioner<1) .or. (whichParallelSolverToFactorPreconditioner>2)) then
       print *,"Error! Invalid setting for whichParallelSolverToFactorPreconditioner"
       stop
    end if

#ifdef PETSC_HAVE_MUMPS
    isAParallelDirectSolverInstalled = .true.
    if (masterProc) then
       print *,"mumps detected"
    end if
#else
    whichParallelSolverToFactorPreconditioner = 2
    if (masterProc) then
       print *,"mumps not detected"
    end if
#endif

#ifdef PETSC_HAVE_SUPERLU_DIST
    isAParallelDirectSolverInstalled = .true.
    if (masterProc) then
       print *,"superlu_dist detected"
    end if
#else
    if (masterProc) then
       print *,"superlu_dist not detected"
    end if
    if (whichParallelSolverToFactorPreconditioner==2) then
       whichParallelSolverToFactorPreconditioner = 1
    end if
#endif

  end subroutine chooseParallelDirectSolver



end module scan

