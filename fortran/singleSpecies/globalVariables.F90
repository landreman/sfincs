module globalVariables

  implicit none

#include <finclude/petscsysdef.h>

  integer, parameter :: integerToRepresentTrue  =  1
  integer, parameter :: integerToRepresentFalse = -1

  PetscScalar :: one = 1., oneHalf = 0.5d+0
  PetscScalar :: zero = 0., two = 2., three = 3., four = 4., five = 5.

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0

  ! ********************************************************
  ! ********************************************************
  !
  ! Options for program flow control:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: programMode, RHSMode, numRHSs
  logical :: saveMatlabOutput, saveMatricesAndVectorsInBinary
  integer :: outputScheme
  character(len=200) :: MatlabOutputFilename
  character(len=200) :: binaryOutputFilename
  character(len=200) :: outputFilename
  logical :: solveSystem

  ! ********************************************************
  ! ********************************************************
  !
  ! Geometry input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: geometryScheme
  PetscScalar ::  GHat, IHat, iota, B0OverBBar
  PetscScalar :: epsilon_t, epsilon_h, epsilon_antisymm
  integer :: NPeriods, helicity_l, helicity_n, helicity_antisymm_l, helicity_antisymm_n
  character(len=200) :: fort996boozer_file, JGboozer_file, JGboozer_file_NonStelSym
  PetscScalar :: normradius_wish, min_Bmn_to_load

  ! ********************************************************
  ! ********************************************************
  !
  ! Physics input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  PetscScalar :: Delta
  PetscScalar :: omega
  PetscScalar :: psiAHat

  PetscScalar :: nuN, nuPrime

  integer :: speciesMode
  ! 0 = ions
  ! 1 = electrons

  PetscScalar :: THat
  PetscScalar :: nHat
  PetscScalar :: EHat
  PetscScalar :: dnHatdpsi
  PetscScalar :: dTHatdpsi
  PetscScalar :: dPhiHatdpsi, EStar, EStarMin, EStarMax
  integer :: NEStar
  PetscScalar, allocatable, dimension(:) :: EStars

  integer :: collisionOperator
  ! 0 = Full linearized Fokker-Planck operator
  ! 1 = pitch-angle scattering with no momentum-conserving term
  ! 2 = pitch-angle scattering with model momentum conserving term

  logical :: includeXDotTerm
  logical :: includeElectricFieldTermInXiDot
  logical :: useDKESExBDrift, include_fDivVE_term


  ! ********************************************************
  ! ********************************************************
  !
  ! Numerical input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: thetaDerivativeScheme
  ! 0 = spectral colocation
  ! 1 = 2nd order finite differences
  ! 2 = 4th order dinite differences

  integer :: Ntheta = 10
  PetscScalar :: NthetaMinFactor = 0.8d+0, NthetaMaxFactor=2d+0
  integer :: NthetaNumRuns = 2

  integer :: Nzeta
  PetscScalar :: NzetaMinFactor, NzetaMaxFactor
  integer :: NzetaNumRuns = 2

  integer :: Nxi = 13
  PetscScalar :: NxiMinFactor = 0.8d+0, NxiMaxFactor=2d+0
  integer :: NxiNumRuns = 2

  integer :: NL = 4
  PetscScalar :: NLMinFactor = 0.8d+0, NLMaxFactor=2d+0
  integer :: NLNumRuns = 2

  integer :: Nx 
  PetscScalar :: NxMinFactor, NxMaxFactor
  integer :: NxNumRuns 

  PetscScalar  :: NxPotentialsPerVth = 2d+0
  PetscScalar :: NxPotentialsPerVthMinFactor = 0.8d+0, NxPotentialsPerVthMaxFactor=2d+0
  integer :: NxPotentialsPerVthNumRuns = 2

  PetscScalar :: xMax
  PetscScalar :: xMaxMinFactor, xMaxMaxFactor
  integer :: xMaxNumRuns

  PetscScalar :: solverTolerance = 1d-6
  PetscScalar :: solverToleranceMinFactor = 1d-1
  PetscScalar :: solverToleranceMaxFactor = 1
  integer :: solverToleranceNumRuns = 0

  logical :: forceOddNthetaAndNzeta = .true.
  ! If forceOddNthetaAndNzeta is set to true, 1 is added to Ntheta any time a run is attempted with even Ntheta,
  ! and 1 is added to Nzeta any time a run is attempted with even Nzeta.
  ! This can be useful because the iterative solvers sometimes do not work with even Ntheta.

  logical :: useIterativeSolver = .false.

  integer :: whichParallelSolverToFactorPreconditioner = 1
  ! Options for whichParallelSolverToFactorPreconditioner:
  ! 1 = use mumps if it is detected, otherwise use superlu_dist
  ! 2 = force use of superlu_dist, if it is available

  logical :: isAParallelDirectSolverInstalled = .false.

  PetscScalar :: multipleOfNuPrimeToAddOnDiagonalOfPreconditioner = 1d+0

  integer :: preconditioner_x, preconditioner_x_min_L, preconditioner_zeta
  integer :: preconditioner_theta, preconditioner_xi

  integer :: constraintScheme

  integer :: PETSCPreallocationStrategy

  ! ********************************************************
  !
  !  Outputs and numerical data which will be saved in the output file
  !
  ! ********************************************************

  PetscScalar, dimension(:), allocatable :: theta, zeta, sources
  PetscScalar, dimension(:,:), allocatable :: BHat, dBHatdtheta, dBHatdzeta
  PetscScalar, dimension(:,:), allocatable :: flow, densityPerturbation, pressurePerturbation
  PetscScalar, dimension(:,:), allocatable :: particleFluxBeforeSurfaceIntegral
  PetscScalar, dimension(:,:), allocatable :: momentumFluxBeforeSurfaceIntegral
  PetscScalar, dimension(:,:), allocatable :: heatFluxBeforeSurfaceIntegral
  PetscScalar, dimension(:,:), allocatable :: NTVBeforeSurfaceIntegral
  PetscScalar :: FSADensityPerturbation, FSAFlow, FSAPressurePerturbation
  PetscScalar :: particleFlux, momentumFlux, heatFlux, NTV, VPrimeHat, FSABHat2

  PetscLogDouble :: elapsedTime
  integer :: didItConverge
  PetscScalar :: transportMatrix(3,3)

  ! ********************************************************
  !
  !  Variables related to parallelization:
  !
  ! ********************************************************

  logical :: parallelizeOverScan

  integer :: numProcs, myRank 
  logical :: masterProc
  ! The above quantities refer to the PETSC_COMM_WORLD communicator, not to the smaller communicators
  ! used for parameter scans.

  ! The quantities below refer to the sub-communicator:
  integer :: numCommunicators
  integer, dimension(:), allocatable :: commMinProcs
  integer, dimension(:), allocatable :: minUnits, maxUnits
  MPI_Comm :: MPIComm
  integer :: myRankInSubComm, numProcsInSubComm
  logical :: masterProcInSubComm
  integer :: myCommunicatorIndex

contains

  ! ------------------------------------------------------------------------

  subroutine setConstraintScheme()

    implicit none

    if (constraintScheme < 0) then
       if (collisionOperator == 0) then
          constraintScheme = 1
       else
          constraintScheme = 2
       end if
    end if

  end subroutine setConstraintScheme

  ! ------------------------------------------------------------------------

  subroutine validateInput()

    implicit none

    ! Validate some input quantities:

    if (preconditioner_theta < 0 .or. preconditioner_theta > 1) then
       print *,"Error! preconditioner_theta must be 0 or 1."
       stop
    end if
    if (preconditioner_zeta < 0 .or. preconditioner_zeta > 1) then
       print *,"Error! preconditioner_zeta must be 0 or 1."
       stop
    end if
    if (preconditioner_xi < 0 .or. preconditioner_xi > 1) then
       print *,"Error! preconditioner_xi must be 0 or 1."
       stop
    end if
    if (preconditioner_x < 0 .or. preconditioner_x > 4) then
       print *,"Error! preconditioner_x must be in the range [0, 4]."
       stop
    end if

    if (collisionOperator < 0 .or. collisionOperator > 2) then
       print *,"Error! collisionOperator must be 0, 1, or 2."
       stop
    end if

    if (constraintScheme < 0 .or. constraintScheme > 2) then
       print *,"Error! constraintScheme must be 0, 1, or 2."
       stop
    end if

  end subroutine validateInput

  ! ------------------------------------------------------------------------

  subroutine deallocateArrays()

    implicit none

    deallocate(zeta)
    deallocate(theta)
    deallocate(BHat)
    deallocate(dBHatdtheta)
    deallocate(dBHatdzeta)
   
    if (masterProcInSubComm) then
       deallocate(densityPerturbation)
       deallocate(pressurePerturbation)
       deallocate(flow)
       deallocate(particleFluxBeforeSurfaceIntegral)
       deallocate(momentumFluxBeforeSurfaceIntegral)
       deallocate(heatFluxBeforeSurfaceIntegral)
       deallocate(NTVBeforeSurfaceIntegral)
       if (constraintScheme > 0) then
          deallocate(sources)
       end if
    end if

  end subroutine deallocateArrays

end module globalVariables

