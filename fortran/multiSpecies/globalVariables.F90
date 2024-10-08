module globalVariables

  implicit none

#include <finclude/petscsysdef.h>

  character(len=50), parameter :: inputFilename = "input.namelist"

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
  PetscScalar :: normradius_wish, min_Bmn_to_load, dGdpHat

  ! ********************************************************
  ! ********************************************************
  !
  ! Species quantities
  !
  ! ********************************************************
  ! ********************************************************

  integer, parameter :: maxNumSpecies = 100
  integer, parameter :: speciesNotInitialized = -9999
  integer :: NSpecies
  PetscScalar, dimension(maxNumSpecies) :: Zs, mHats, NHats, THats, dNHatdpsiNs, dTHatdpsiNs

  ! ********************************************************
  ! ********************************************************
  !
  ! Physics input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  PetscScalar :: Delta
  PetscScalar :: alpha
  PetscScalar :: psiAHat
  PetscScalar :: nu_n

  integer :: speciesMode
  ! 0 = ions
  ! 1 = electrons

  PetscScalar :: EParallelHat
  PetscScalar :: dPhiHatdpsiN, dPhiHatdpsiN_min, dPhiHatdpsiN_max
  integer :: NErs
  PetscScalar, allocatable, dimension(:) :: dPhiHatdpsiNs

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
  ! 0 = spectral collocation
  ! 1 = 2nd order finite differences
  ! 2 = 4th order dinite differences

  integer :: xGridScheme = 1, xPotentialsGridScheme = 2
  logical :: pointAtX0

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
  ! This can be useful because the iterative solvers sometimes do not work with even Ntheta or Nzeta.
  ! This parameter should be true unless you know what you are doing.

  logical :: useIterativeSolver = .false.

  integer :: whichParallelSolverToFactorPreconditioner = 1
  ! Options for whichParallelSolverToFactorPreconditioner:
  ! 1 = use mumps if it is detected, otherwise use superlu_dist
  ! 2 = force use of superlu_dist, if it is available

  logical :: isAParallelDirectSolverInstalled = .false.

  PetscScalar :: multipleOfNuPrimeToAddOnDiagonalOfPreconditioner = 1d+0

  integer :: preconditioner_x, preconditioner_x_min_L, preconditioner_zeta
  integer :: preconditioner_theta, preconditioner_xi, preconditioner_species
  integer :: preconditioner_theta_min_L = 2
  integer :: preconditioner_zeta_min_L = 2

  integer :: constraintScheme

  integer :: PETSCPreallocationStrategy

  ! ********************************************************
  !
  !  Outputs and numerical data which will be saved in the output file
  !
  ! ********************************************************

  integer :: matrixSize
  PetscScalar :: normradius
  PetscScalar, dimension(:), allocatable :: theta, zeta
  PetscScalar, dimension(:,:), allocatable :: BHat, dBHatdtheta, dBHatdzeta, sources, jHat, Phi1Hat
  PetscScalar, dimension(:,:,:), allocatable :: flow, densityPerturbation, pressurePerturbation
  PetscScalar, dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral
  PetscScalar, dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral
  PetscScalar, dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral
  PetscScalar, dimension(:,:,:), allocatable :: NTVBeforeSurfaceIntegral
  PetscScalar, dimension(:,:), allocatable :: NTVKernel
  PetscScalar, dimension(:), allocatable :: FSADensityPerturbation, FSABFlow, FSAPressurePerturbation
  PetscScalar, dimension(:), allocatable :: particleFlux, momentumFlux, heatFlux, NTV
  PetscScalar :: VPrimeHat, FSABHat2, FSABjHat

  PetscLogDouble :: elapsedTime
  integer :: didItConverge
  PetscScalar :: transportMatrix(3,3)

  !!Added by AM 2014-09!!
  PetscScalar :: ArrayFirstSpeciesParticleFluxCoefficients(3)
  !!Added by AM 2015-05
  PetscScalar :: ArrayFirstSpeciesHeatFluxCoefficients(3)
  PetscScalar :: ArraySecondSpeciesParticleFluxCoefficients(3)
  PetscScalar :: ArraySecondSpeciesHeatFluxCoefficients(3)
  !!!!!!!!!!!!!!!!!!!!!!!

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

    if (preconditioner_theta < 0 .or. preconditioner_theta > 3) then
       print *,"Error! preconditioner_theta must be 0, 1, 2, or 3."
       stop
    end if
    if (preconditioner_zeta < 0 .or. preconditioner_zeta > 3) then
       print *,"Error! preconditioner_zeta must be 0, 1, 2, or 3."
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
    deallocate(NTVKernel)
   
    if (masterProcInSubComm) then
       deallocate(FSADensityPerturbation)
       deallocate(FSABFlow)
       deallocate(FSAPressurePerturbation)
       deallocate(particleFlux)
       deallocate(momentumFlux)
       deallocate(heatFlux)
       deallocate(NTV)

       deallocate(densityPerturbation)
       deallocate(flow)
       deallocate(pressurePerturbation)
       deallocate(particleFluxBeforeSurfaceIntegral)
       deallocate(momentumFluxBeforeSurfaceIntegral)
       deallocate(heatFluxBeforeSurfaceIntegral)
       deallocate(NTVBeforeSurfaceIntegral)

       deallocate(jHat)
       deallocate(Phi1Hat)

       if (constraintScheme > 0) then
          deallocate(sources)
       end if
    end if

  end subroutine deallocateArrays

end module globalVariables

