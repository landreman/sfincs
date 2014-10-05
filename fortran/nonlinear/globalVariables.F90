module globalVariables

  implicit none

#include <finclude/petscsysdef.h>

  character(len=50), parameter :: inputFilename = "input.namelist"

  integer, parameter :: integerToRepresentTrue  =  1
  integer, parameter :: integerToRepresentFalse = -1

  PetscScalar, parameter :: one = 1., oneHalf = 0.5d+0
  PetscScalar, parameter :: zero = 0., two = 2., three = 3., four = 4., five = 5.

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0

  ! ********************************************************
  ! ********************************************************
  !
  ! Options for program flow control:
  !
  ! ********************************************************
  ! ********************************************************

  logical :: saveMatlabOutput, saveMatricesAndVectorsInBinary
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
  PetscScalar ::  GHat, IHat, iota, B0OverBBar, psiAHat
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

  logical :: nonlinear
  PetscScalar :: Delta
  PetscScalar :: alpha
  PetscScalar :: nu_n

  PetscScalar :: EParallelHat
  PetscScalar :: dPhiHatdpsiN

  integer :: collisionOperator
  ! 0 = Full linearized Fokker-Planck operator
  ! 1 = pitch-angle scattering with no momentum-conserving term

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
  integer :: zetaDerivativeScheme
  ! 0 = spectral collocation
  ! 1 = 2nd order finite differences
  ! 2 = 4th order dinite differences

  integer :: Ntheta
  integer :: Nzeta
  integer :: Nxi
  integer :: NL
  integer :: Nx 
  PetscScalar  :: NxPotentialsPerVth
  PetscScalar :: xMax
  PetscScalar :: solverTolerance

  logical :: forceOddNthetaAndNzeta
  ! If forceOddNthetaAndNzeta is set to true, 1 is added to Ntheta any time a run is attempted with even Ntheta,
  ! and 1 is added to Nzeta any time a run is attempted with even Nzeta.
  ! This can be useful because the iterative solvers sometimes do not work with even Ntheta or Nzeta.
  ! This parameter should be true unless you know what you are doing.

  logical :: useIterativeSolver

  integer :: whichParallelSolverToFactorPreconditioner
  ! Options for whichParallelSolverToFactorPreconditioner:
  ! 1 = use mumps if it is detected, otherwise use superlu_dist
  ! 2 = force use of superlu_dist, if it is available

  integer :: preconditioner_x, preconditioner_x_min_L, preconditioner_zeta
  integer :: preconditioner_theta, preconditioner_xi, preconditioner_species

  integer :: constraintScheme

  integer :: PETSCPreallocationStrategy

  ! ********************************************************
  !
  !  Other variables that are used by multiple subroutines
  !
  ! ********************************************************

  integer :: matrixSize, NxPotentials
  PetscScalar :: normradius
  PetscScalar, dimension(:), allocatable :: theta, zeta, x
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights
  PetscScalar, dimension(:,:), allocatable :: ddtheta, ddzeta
  PetscScalar, dimension(:), allocatable :: xWeights, xPotentials
  PetscScalar, dimension(:), allocatable :: x2, expx2
  PetscScalar, dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
  PetscScalar, dimension(:,:), allocatable :: ddx_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddtheta_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddzeta_preconditioner
  PetscScalar, dimension(:,:), allocatable :: regridPolynomialToUniform

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
  integer :: didLinearCalculationConverge, didNonlinearCalculationConverge

  ! ********************************************************
  !
  !  Variables related to parallelization:
  !
  ! ********************************************************

  integer :: numProcs, myRank 
  logical :: masterProc
  integer :: izetaMin, izetaMax, localNzeta
  logical :: procThatHandlesConstraints

end module globalVariables

