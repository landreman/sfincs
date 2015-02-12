module globalVariables

  implicit none

#include <finclude/petscvecdef.h>

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
  ! General options:
  !
  ! ********************************************************
  ! ********************************************************

  logical :: saveMatlabOutput, saveMatricesAndVectorsInBinary
  character(len=200) :: MatlabOutputFilename
  character(len=200) :: binaryOutputFilename
  character(len=200) :: outputFilename
  logical :: solveSystem
  integer :: RHSMode, whichRHS

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
  character(len=200) :: equilibriumFile
  PetscScalar :: min_Bmn_to_load, dGdpHat, aHat
  PetscScalar :: psiHat_wish, psiN_wish, rHat_wish, rN_wish
  PetscScalar :: psiHat, psiN, rHat, rN
  integer :: inputRadialCoordinate = 2
  integer :: inputRadialCoordinateForGradients = 2
  logical :: force0RadialCurrentInEquilibrium = .true.
  integer :: VMECRadialOption = 0

  ! ********************************************************
  ! ********************************************************
  !
  ! Species quantities
  !
  ! ********************************************************
  ! ********************************************************

  integer, parameter :: maxNumSpecies = 100
  integer, parameter :: speciesNotInitialized = -9999
  integer :: NSpecies = -9999
  PetscScalar, dimension(maxNumSpecies) :: Zs = 1, mHats = 1.0, NHats = 1.0, THats = 1.0
  PetscScalar, dimension(maxNumSpecies) :: dNHatdpsiHats = 0, dTHatdpsiHats = 0
  PetscScalar, dimension(maxNumSpecies) :: dNHatdpsiNs   = 0, dTHatdpsiNs   = 0
  PetscScalar, dimension(maxNumSpecies) :: dNHatdrHats   = 0, dTHatdrHats   = 0
  PetscScalar, dimension(maxNumSpecies) :: dNHatdrNs     = 0, dTHatdrNs     = 0

  ! ********************************************************
  ! ********************************************************
  !
  ! Physics input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  logical :: nonlinear = .false.
  PetscScalar :: Delta
  PetscScalar :: alpha
  PetscScalar :: nu_n

  PetscScalar :: EParallelHat = 0
  PetscScalar :: dPhiHatdpsiHat = 0, dPhiHatdpsiN = 0, dPhiHatdrHat = 0, dPhiHatdrN = 0

  integer :: collisionOperator = 0
  ! 0 = Full linearized Fokker-Planck operator
  ! 1 = pitch-angle scattering with no momentum-conserving term

  logical :: includeXDotTerm = .true.
  logical :: includeElectricFieldTermInXiDot = .true.
  logical :: useDKESExBDrift = .true.
  logical :: include_fDivVE_term = .true.
  logical :: includeTemperatureEquilibrationTerm = .true.
  logical :: includePhi1 = .true.

  PetscScalar :: nuPrime = 0, EStar = 0

  ! ********************************************************
  ! ********************************************************
  !
  ! Numerical input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: thetaDerivativeScheme = 2
  integer :: zetaDerivativeScheme = 2
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

  logical :: useIterativeLinearSolver

  integer :: whichParallelSolverToFactorPreconditioner
  ! Options for whichParallelSolverToFactorPreconditioner:
  ! 1 = use mumps if it is detected, otherwise use superlu_dist
  ! 2 = force use of superlu_dist, if it is available

  integer :: preconditioner_x, preconditioner_x_min_L, preconditioner_zeta
  integer :: preconditioner_theta, preconditioner_xi, preconditioner_species
  logical :: reusePreconditioner

  integer :: constraintScheme

  integer :: PETSCPreallocationStrategy

  ! ********************************************************
  !
  !  Other variables that are used by multiple subroutines
  !
  ! ********************************************************

  integer :: matrixSize, NxPotentials
  PetscScalar, dimension(:), allocatable :: theta, zeta, x
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights
  PetscScalar, dimension(:,:), allocatable :: ddtheta, ddzeta
  PetscScalar, dimension(:), allocatable :: xWeights, xPotentials
  PetscScalar :: maxxPotentials, zetaMax, xMaxNotTooSmall
  PetscScalar, dimension(:), allocatable :: x2, expx2
  PetscScalar, dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
  PetscScalar, dimension(:,:), allocatable :: ddx_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddtheta_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddzeta_preconditioner
  PetscScalar, dimension(:,:), allocatable :: regridPolynomialToUniform

  integer, parameter :: COORDINATE_SYSTEM_UNINITIALIZED = -1
  integer, parameter :: COORDINATE_SYSTEM_BOOZER = 1
  integer, parameter :: COORDINATE_SYSTEM_VMEC = 2
  integer :: coordinateSystem = COORDINATE_SYSTEM_UNINITIALIZED

  PetscScalar, dimension(:,:), allocatable :: BHat, dBHatdtheta, dBHatdzeta, dBHatdpsiHat, DHat
  PetscScalar, dimension(:,:), allocatable :: BHat_sub_psi, dBHat_sub_psi_dtheta, dBHat_sub_psi_dzeta
  PetscScalar, dimension(:,:), allocatable :: BHat_sub_theta, dBHat_sub_theta_dzeta, dBHat_sub_theta_dpsiHat
  PetscScalar, dimension(:,:), allocatable :: BHat_sub_zeta, dBHat_sub_zeta_dtheta, dBHat_sub_zeta_dpsiHat
  PetscScalar, dimension(:,:), allocatable :: BHat_sup_theta, dBHat_sup_theta_dzeta, dBHat_sup_theta_dpsiHat
  PetscScalar, dimension(:,:), allocatable :: BHat_sup_zeta, dBHat_sup_zeta_dtheta, dBHat_sup_zeta_dpsiHat
  PetscScalar, dimension(:,:), allocatable :: sources, jHat, Phi1Hat, dPhi1Hatdtheta, dPhi1Hatdzeta
  PetscScalar, dimension(:,:,:), allocatable :: densityPerturbation, totalDensity
  PetscScalar, dimension(:,:,:), allocatable :: pressurePerturbation, totalPressure
  PetscScalar, dimension(:,:,:), allocatable :: flow, velocityUsingFSADensity, velocityUsingTotalDensity
  PetscScalar, dimension(:,:,:), allocatable :: MachUsingFSAThermalSpeed
  PetscScalar, dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vm0
  PetscScalar, dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vm
  PetscScalar, dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vE0
  PetscScalar, dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vE
  PetscScalar, dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vm0
  PetscScalar, dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vm
  PetscScalar, dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vE0
  PetscScalar, dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vE
  PetscScalar, dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vm0
  PetscScalar, dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vm
  PetscScalar, dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vE0
  PetscScalar, dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vE
  PetscScalar, dimension(:,:,:), allocatable :: NTVBeforeSurfaceIntegral
  PetscScalar, dimension(:,:), allocatable :: NTVKernel
  PetscScalar, dimension(:), allocatable :: FSADensityPerturbation, FSAPressurePerturbation
  PetscScalar, dimension(:), allocatable :: FSABFlow, FSABVelocityUsingFSADensity
  PetscScalar, dimension(:), allocatable :: FSABVelocityUsingFSADensityOverB0
  PetscScalar, dimension(:), allocatable :: FSABVelocityUsingFSADensityOverRootFSAB2

  PetscScalar, dimension(:), allocatable :: particleFlux_vm0_psiHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vm0_psiN
  PetscScalar, dimension(:), allocatable :: particleFlux_vm0_rHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vm0_rN

  PetscScalar, dimension(:), allocatable :: particleFlux_vm_psiHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vm_psiN
  PetscScalar, dimension(:), allocatable :: particleFlux_vm_rHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vm_rN

  PetscScalar, dimension(:), allocatable :: particleFlux_vE0_psiHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vE0_psiN
  PetscScalar, dimension(:), allocatable :: particleFlux_vE0_rHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vE0_rN

  PetscScalar, dimension(:), allocatable :: particleFlux_vE_psiHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vE_psiN
  PetscScalar, dimension(:), allocatable :: particleFlux_vE_rHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vE_rN

  PetscScalar, dimension(:), allocatable :: particleFlux_vd1_psiHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vd1_psiN
  PetscScalar, dimension(:), allocatable :: particleFlux_vd1_rHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vd1_rN

  PetscScalar, dimension(:), allocatable :: particleFlux_vd_psiHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vd_psiN
  PetscScalar, dimension(:), allocatable :: particleFlux_vd_rHat
  PetscScalar, dimension(:), allocatable :: particleFlux_vd_rN

  PetscScalar, dimension(:), allocatable :: momentumFlux_vm0_psiHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vm0_psiN
  PetscScalar, dimension(:), allocatable :: momentumFlux_vm0_rHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vm0_rN

  PetscScalar, dimension(:), allocatable :: momentumFlux_vm_psiHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vm_psiN
  PetscScalar, dimension(:), allocatable :: momentumFlux_vm_rHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vm_rN

  PetscScalar, dimension(:), allocatable :: momentumFlux_vE0_psiHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vE0_psiN
  PetscScalar, dimension(:), allocatable :: momentumFlux_vE0_rHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vE0_rN

  PetscScalar, dimension(:), allocatable :: momentumFlux_vE_psiHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vE_psiN
  PetscScalar, dimension(:), allocatable :: momentumFlux_vE_rHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vE_rN

  PetscScalar, dimension(:), allocatable :: momentumFlux_vd1_psiHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vd1_psiN
  PetscScalar, dimension(:), allocatable :: momentumFlux_vd1_rHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vd1_rN

  PetscScalar, dimension(:), allocatable :: momentumFlux_vd_psiHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vd_psiN
  PetscScalar, dimension(:), allocatable :: momentumFlux_vd_rHat
  PetscScalar, dimension(:), allocatable :: momentumFlux_vd_rN

  PetscScalar, dimension(:), allocatable :: heatFlux_vm0_psiHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vm0_psiN
  PetscScalar, dimension(:), allocatable :: heatFlux_vm0_rHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vm0_rN

  PetscScalar, dimension(:), allocatable :: heatFlux_vm_psiHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vm_psiN
  PetscScalar, dimension(:), allocatable :: heatFlux_vm_rHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vm_rN

  PetscScalar, dimension(:), allocatable :: heatFlux_vE0_psiHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vE0_rHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vE0_rN
  PetscScalar, dimension(:), allocatable :: heatFlux_vE0_psiN

  PetscScalar, dimension(:), allocatable :: heatFlux_vE_psiHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vE_rHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vE_rN
  PetscScalar, dimension(:), allocatable :: heatFlux_vE_psiN

  PetscScalar, dimension(:), allocatable :: heatFlux_vd1_psiHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vd1_psiN
  PetscScalar, dimension(:), allocatable :: heatFlux_vd1_rHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vd1_rN

  PetscScalar, dimension(:), allocatable :: heatFlux_vd_psiHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vd_psiN
  PetscScalar, dimension(:), allocatable :: heatFlux_vd_rHat
  PetscScalar, dimension(:), allocatable :: heatFlux_vd_rN

  PetscScalar, dimension(:), allocatable :: heatFlux_withoutPhi1_psiHat
  PetscScalar, dimension(:), allocatable :: heatFlux_withoutPhi1_psiN
  PetscScalar, dimension(:), allocatable :: heatFlux_withoutPhi1_rHat
  PetscScalar, dimension(:), allocatable :: heatFlux_withoutPhi1_rN

  PetscScalar, dimension(:), allocatable :: NTV
  PetscScalar :: VPrimeHat, FSABHat2, FSABjHat, FSABjHatOverB0, FSABjHatOverRootFSAB2
  PetscScalar :: lambda

  PetscScalar :: ddpsiN2ddpsiHat, ddrHat2ddpsiHat, ddrN2ddpsiHat
  PetscScalar :: ddpsiHat2ddpsiN, ddpsiHat2ddrHat, ddpsiHat2ddrN

  PetscLogDouble :: elapsedTime
  integer :: didLinearCalculationConverge, didNonlinearCalculationConverge
  integer :: iterationForMatrixOutput
  logical :: firstMatrixCreation

  integer :: transportMatrixSize = 3
  PetscScalar, dimension(:,:), allocatable :: transportMatrix

  Vec :: f0

  ! ********************************************************
  !
  !  Variables related to parallelization:
  !
  ! ********************************************************

  MPI_Comm :: MPIComm
  integer :: numProcs, myRank 
  logical :: masterProc
  ! Eventually, keep only the zeta variables, and drop the theta variables
  integer :: ithetaMin, ithetaMax, localNtheta
  integer :: izetaMin, izetaMax, localNzeta
  logical :: procThatHandlesConstraints

end module globalVariables

