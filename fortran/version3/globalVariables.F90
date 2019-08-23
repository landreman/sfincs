module globalVariables

#include "PETScVersions.F90"

  implicit none


  character(len=200) :: inputFilename = "input.namelist"

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

  logical :: saveMatlabOutput = .false., saveMatricesAndVectorsInBinary = .false.
  character(len=200) :: MatlabOutputFilename = "sfincsMatrices"
  character(len=200) :: binaryOutputFilename = "sfincsBinary"
  character(len=200) :: outputFilename = "sfincsOutput.h5"
  character(len=200) :: externalPhi1Filename = "externalPhi1.h5" !!Added by AM 2018-12
  logical :: solveSystem = .true.
  integer :: RHSMode = 1, whichRHS
  logical :: isAParallelDirectSolverInstalled

  ! Quantities related to ambipolar solve
  logical :: ambipolarSolve = .false.
  integer :: ambipolarSolveOption = 2
  integer :: NEr_ambipolarSolve = 20
  PetscScalar :: Er_search_tolerance_dx = 1.d-8
  PetscScalar :: Er_search_tolerance_f = 1.d-10
  PetscScalar :: Er_min = -100
  PetscScalar :: Er_max = 100

  ! ********************************************************
  ! ********************************************************
  !
  ! Geometry input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: geometryScheme = 1
  PetscScalar ::  GHat = 3.7481d+0, IHat = 0.0, iota = 0.4542d+0, B0OverBBar = 1.0, psiAHat = 0.15596d+0
  PetscScalar ::  diotadpsiHat = 0.0, pPrimeHat = 0.0 !!Added by HS 2016-09-19
  PetscScalar :: epsilon_t = -0.07053d+0, epsilon_h = 0.05067d+0, epsilon_antisymm = 0.0
  integer :: NPeriods = 0, helicity_l = 2, helicity_n = 10, helicity_antisymm_l = 1, helicity_antisymm_n = 0
  character(len=200) :: equilibriumFile = ""
  character(len=200) :: EParallelHatSpec_bcdatFile = ""
  PetscScalar :: min_Bmn_to_load = 0.0, dGdpHat, aHat = 0.5585d+0
  PetscScalar :: psiHat_wish = -1, psiN_wish = 0.25, rHat_wish = -1, rN_wish = 0.5
  PetscScalar :: psiHat, psiN, rHat, rN
  integer :: inputRadialCoordinate = 3
  integer :: inputRadialCoordinateForGradients = 4
  logical :: force0RadialCurrentInEquilibrium = .true.
  integer :: VMECRadialOption = 1
  integer :: VMEC_Nyquist_option = 1
  PetscScalar :: rippleScale = 1

  ! For use in Boozer bmnc stellopt
  integer, parameter :: mmax_boozer = 40
  integer, parameter :: nmax_boozer = 40
  PetscScalar, dimension(0:mmax_boozer,-nmax_boozer:nmax_boozer) :: boozer_bmnc
  PetscScalar, dimension(0:mmax_boozer,-nmax_boozer:nmax_boozer) :: boozer_bmns

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
  ! For all the array variables declared in the next 5 lines, default values
  ! are set in readInput.F90
  PetscScalar, dimension(maxNumSpecies) :: Zs, mHats, NHats, THats
  PetscScalar, dimension(maxNumSpecies) :: dNHatdpsiHats, dTHatdpsiHats
  PetscScalar, dimension(maxNumSpecies) :: dNHatdpsiNs,   dTHatdpsiNs
  PetscScalar, dimension(maxNumSpecies) :: dNHatdrHats,   dTHatdrHats
  PetscScalar, dimension(maxNumSpecies) :: dNHatdrNs,     dTHatdrNs
  PetscScalar, dimension(maxNumSpecies) :: EParallelHatSpec = 0 !!Added by HS 2017-02

  !!Added by AM 2016-02!!
  PetscScalar :: adiabaticZ = -1.0, adiabaticMHat = 5.446170214d-4, adiabaticNHat = 1.0, adiabaticTHat = 1.0
  logical :: withAdiabatic = .false.

  !!Added by HM 2017-09!!
  PetscScalar :: NBIspecZ = 1.0,  NBIspecNHat = 0.0
  logical :: withNBIspec = .false.
  !!!!!!!!!!!!!!!!!!!!!!!

  ! ********************************************************
  ! ********************************************************
  !
  ! Sensitivity options:
  !
  ! ********************************************************
  ! ********************************************************

  logical :: debugAdjoint = .false.
  PetscScalar :: deltaLambda = 1.d-4 ! for finite diff testing
  logical :: adjointBootstrapOption = .false.
  logical :: adjointRadialCurrentOption = .false.
  logical :: adjointTotalHeatFluxOption = .false.
  ! These are initialized to .false. in readInput
  logical, dimension(:), allocatable :: adjointHeatFluxOption, adjointParticleFluxOption, adjointParallelFlowOption
  logical :: discreteAdjointOption = .true.


  ! ********************************************************
  ! ********************************************************
  !
  ! Physics input parameters:
  !
  ! ********************************************************
  ! ********************************************************

!!  logical :: nonlinear = .false.
 !!Commented by AM 2016-02
  PetscScalar :: Delta = 4.5694d-3
  PetscScalar :: alpha = 1.0
  PetscScalar :: nu_n = 8.330d-3

  PetscScalar :: EParallelHat = 0
  PetscScalar :: dPhiHatdpsiHat = 0, dPhiHatdpsiN = 0, dPhiHatdrHat = 0, dPhiHatdrN = 0, Er = 0

  integer :: collisionOperator = 0
  ! 0 = Full linearized Fokker-Planck operator
  ! 1 = pitch-angle scattering with no momentum-conserving term

  logical :: includeXDotTerm = .true.
  logical :: includeElectricFieldTermInXiDot = .true.
  logical :: useDKESExBDrift = .false.
  logical :: include_fDivVE_term = .false.
  logical :: includeTemperatureEquilibrationTerm = .false.
  logical :: includePhi1 = .false.
  logical :: includePhi1InCollisionOperator = .false.
!!  logical :: includeRadialExBDrive = .false. !!Commented by AM 2016-02
  logical :: includePhi1InKineticEquation = .true. !!Added by AM 2016-03
  logical :: readExternalPhi1 = .false. !!Added by AM 2018-12

  PetscScalar :: nuPrime = 1, EStar = 0

  integer :: magneticDriftScheme = 0

  integer :: quasineutralityOption = 1

  PetscScalar :: Krook = 0

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

  integer :: ExBDerivativeSchemeTheta = 0
  integer :: ExBDerivativeSchemeZeta = 0
  integer :: magneticDriftDerivativeScheme = 3
  integer :: xDotDerivativeScheme = 0

  integer :: xPotentialsGridScheme = 2, xPotentialsInterpolationScheme
  integer :: xGridScheme = 5, xInterpolationScheme
  logical :: pointAtX0

  integer :: Ntheta = 15
  integer :: Nzeta = 15
  integer :: Nxi = 16
  integer :: NL = 4
  integer :: Nx = 5
  PetscScalar  :: NxPotentialsPerVth = 40.0
  PetscScalar :: xMax = 5.0
  PetscScalar :: solverTolerance = 1d-6

  logical :: forceOddNthetaAndNzeta=.true.
  ! If forceOddNthetaAndNzeta is set to true, 1 is added to Ntheta any time a run is attempted with even Ntheta,
  ! and 1 is added to Nzeta any time a run is attempted with even Nzeta.
  ! This can be useful because the iterative solvers sometimes do not work with even Ntheta or Nzeta.
  ! This parameter should be true unless you know what you are doing.

  logical :: useIterativeLinearSolver=.true.

  integer :: whichParallelSolverToFactorPreconditioner = 1
  ! Options for whichParallelSolverToFactorPreconditioner:
  ! 1 = use mumps if it is detected, otherwise use superlu_dist
  ! 2 = force use of superlu_dist, if it is available

  integer :: preconditioner_x=1, preconditioner_x_min_L=0, preconditioner_zeta=0
  integer :: preconditioner_theta=0, preconditioner_xi=1, preconditioner_species=1
  integer :: preconditioner_theta_min_L=0, preconditioner_zeta_min_L=0
  logical :: reusePreconditioner=.true.
  integer :: preconditioner_magnetic_drifts_max_L = 2

  integer :: constraintScheme=-1

  integer :: PETSCPreallocationStrategy=1
  integer :: Nxi_for_x_option = 1

  ! ********************************************************
  !
  !  Other variables that are used by multiple subroutines
  !
  ! ********************************************************

  integer, dimension(:), allocatable :: Nxi_for_x, min_x_for_L
  integer :: matrixSize, NxPotentials
  PetscScalar, dimension(:), allocatable :: theta, zeta, x, x_plus1
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights
  PetscScalar, dimension(:,:), allocatable :: ddtheta, ddzeta
  PetscScalar, dimension(:,:), allocatable :: ddtheta_ExB_plus, ddtheta_ExB_minus
  PetscScalar, dimension(:,:), allocatable :: ddzeta_ExB_plus, ddzeta_ExB_minus
  PetscScalar, dimension(:,:), allocatable :: ddtheta_magneticDrift_plus, ddtheta_magneticDrift_minus
  PetscScalar, dimension(:,:), allocatable :: ddzeta_magneticDrift_plus, ddzeta_magneticDrift_minus
  PetscScalar, dimension(:), allocatable :: xWeights, xPotentials
  PetscScalar :: maxxPotentials, zetaMax, xMaxNotTooSmall
  PetscScalar, dimension(:), allocatable :: x2, expx2
  PetscScalar, dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
  PetscScalar, dimension(:,:), allocatable :: ddx_xDot_plus, ddx_xDot_preconditioner_plus
  PetscScalar, dimension(:,:), allocatable :: ddx_xDot_minus, ddx_xDot_preconditioner_minus
  PetscScalar, dimension(:,:), allocatable :: ddx_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddtheta_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddzeta_preconditioner
  PetscScalar, dimension(:,:), allocatable :: interpolateXToXPotentials

  PetscScalar, dimension(:,:,:,:,:), allocatable :: RosenbluthPotentialTerms

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
  PetscScalar, dimension(:,:), allocatable :: BDotCurlB, uHat, gradpsidotgradB_overgpsipsi, gpsipsi
  PetscScalar, dimension(:,:), allocatable :: sources, jHat, Phi1Hat, dPhi1Hatdtheta, dPhi1Hatdzeta, bcdata
  PetscScalar, dimension(:,:,:), allocatable :: densityPerturbation, totalDensity
  PetscScalar, dimension(:,:,:), allocatable :: pressurePerturbation, totalPressure, pressureAnisotropy
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

  PetscScalar, dimension(:), allocatable :: classicalParticleFlux_psiHat
  PetscScalar, dimension(:), allocatable :: classicalParticleFlux_psiN
  PetscScalar, dimension(:), allocatable :: classicalParticleFlux_rHat
  PetscScalar, dimension(:), allocatable :: classicalParticleFlux_rN
  PetscScalar, dimension(:), allocatable :: classicalParticleFluxNoPhi1_psiHat
  PetscScalar, dimension(:), allocatable :: classicalParticleFluxNoPhi1_psiN
  PetscScalar, dimension(:), allocatable :: classicalParticleFluxNoPhi1_rHat
  PetscScalar, dimension(:), allocatable :: classicalParticleFluxNoPhi1_rN
  PetscScalar, dimension(:), allocatable :: classicalHeatFlux_psiHat
  PetscScalar, dimension(:), allocatable :: classicalHeatFlux_psiN
  PetscScalar, dimension(:), allocatable :: classicalHeatFlux_rHat
  PetscScalar, dimension(:), allocatable :: classicalHeatFlux_rN
  PetscScalar, dimension(:), allocatable :: classicalHeatFluxNoPhi1_psiHat
  PetscScalar, dimension(:), allocatable :: classicalHeatFluxNoPhi1_psiN
  PetscScalar, dimension(:), allocatable :: classicalHeatFluxNoPhi1_rHat
  PetscScalar, dimension(:), allocatable :: classicalHeatFluxNoPhi1_rN

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

  PetscScalar, dimension(:,:), allocatable :: particleFlux_vm_psiHat_vs_x
  PetscScalar, dimension(:,:), allocatable :: heatFlux_vm_psiHat_vs_x
  PetscScalar, dimension(:,:), allocatable :: FSABFlow_vs_x

  PetscScalar, dimension(:), allocatable :: NTV
  PetscScalar :: VPrimeHat, FSABHat2, FSABHat, FSABjHat, FSABjHatOverB0, FSABjHatOverRootFSAB2
  PetscScalar :: lambda

  PetscScalar :: ddpsiN2ddpsiHat, ddrHat2ddpsiHat, ddrN2ddpsiHat
  PetscScalar :: ddpsiHat2ddpsiN, ddpsiHat2ddrHat, ddpsiHat2ddrN

  PetscLogDouble :: elapsedTime
  integer :: didLinearCalculationConverge, didNonlinearCalculationConverge
  integer :: iterationForMatrixOutput, iterationForResidualOutput = 0, iterationForStateVectorOutput = 0
  logical :: firstMatrixCreation

  integer :: transportMatrixSize = 3
  PetscScalar, dimension(:,:), allocatable :: transportMatrix

  Vec :: f0

! ********************************************************
  !
  !  Variables related to sensitivity:
  !
  ! ********************************************************

  !> @var ns Values of n to use for Fourier derivatives.
  !> @var ms Values of m to sue for Fourier derivatives.
  integer, dimension(:), allocatable :: ns_sensitivity, ms_sensitivity
  integer :: nMaxAdjoint = 0, mMaxAdjoint = 0
  integer :: nMinAdjoint = 0, mMinAdjoint = 0
  integer :: NModesAdjoint = 0, NLambdas = 4
  PetscScalar, dimension(:,:), allocatable :: dRadialCurrentdLambda,dTotalHeatFluxdLambda,dBootstrapdLambda
  PetscScalar, dimension(:,:,:), allocatable :: dParticleFluxdLambda, dHeatFluxdLambda, dParallelFlowdLambda
  PetscScalar, dimension(:,:), allocatable :: dPhidPsidLambda
  ! Below are quantities related to finite diff testing
  PetscScalar, dimension(:,:), allocatable :: dRadialCurrentdLambda_finitediff,dTotalHeatFluxdLambda_finitediff,dBootstrapdLambda_finitediff
  PetscScalar, dimension(:,:,:), allocatable :: dParticleFluxdLambda_finitediff, dHeatFluxdLambda_finitediff, dParallelFlowdLambda_finitediff
  PetscScalar, dimension(:,:,:), allocatable :: particleFluxPercentError, heatFluxPercentError, parallelFlowPercentError
  PetscScalar, dimension(:,:), allocatable :: radialCurrentPercentError, totalHeatFluxPercentError, bootstrapPercentError
  PetscScalar, dimension(:), allocatable :: bmnc_init
  PetscScalar, dimension(:,:), allocatable :: DHat_init, BHat_init, dBHatdtheta_init, dBHatdzeta_init, BHat_sup_theta_init, dBHat_sup_theta_dzeta_init, BHat_sup_zeta_init, dBHat_sup_zeta_dtheta_init, BHat_sub_theta_init, dBHat_sub_theta_dzeta_init, BHat_sub_zeta_init, dBHat_sub_zeta_dtheta_init
  PetscScalar :: GHat_init, IHat_init, iota_init
  PetscScalar, dimension(:,:), allocatable :: dPhidPsidLambda_finitediff, dPhidPsiPercentError
  ! Related to ambipolar solve
  PetscScalar :: dRadialCurrentdEr

  ! ********************************************************
  !
  !  Variables related to parallelization:
  !
  ! ********************************************************

  MPI_Comm :: MPIComm
  integer :: numProcs, myRank 
  logical :: masterProc

  integer :: ithetaMin, ithetaMax, localNtheta
  integer :: izetaMin, izetaMax, localNzeta
  logical :: procThatHandlesConstraints

end module globalVariables

