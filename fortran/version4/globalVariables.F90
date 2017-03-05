module globalVariables

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

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

  logical :: saveMatlabOutput = .false., saveMatricesAndVectorsInBinary = .false.
  character(len=200) :: MatlabOutputFilename = "sfincsMatrices"
  character(len=200) :: binaryOutputFilename = "sfincsBinary"
  character(len=200) :: outputFilename = "sfincsOutput.h5"
  logical :: solveSystem = .true.
  integer :: RHSMode = 1, whichRHS
  logical :: isAParallelDirectSolverInstalled

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
  PetscScalar :: min_Bmn_to_load = 0.0, dGdpHat, aHat = 0.5585d+0
  PetscScalar :: psiHat_wish = -1, psiN_wish = 0.25, rHat_wish = -1, rN_wish = 0.5
  PetscScalar :: psiHat, psiN, rHat, rN
  integer :: inputRadialCoordinate = 3
  integer :: inputRadialCoordinateForGradients = 4
  logical :: force0RadialCurrentInEquilibrium = .true.
  integer :: VMECRadialOption = 1
  PetscScalar :: rippleScale = 1

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

  !!Added by AM 2016-02!!
  PetscScalar :: adiabaticZ = -1.0, adiabaticMHat = 5.446170214d-4, adiabaticNHat = 1.0, adiabaticTHat = 1.0
  logical :: withAdiabatic = .false.
  !!!!!!!!!!!!!!!!!!!!!!!


  ! ********************************************************
  ! ********************************************************
  !
  ! Physics input parameters:
  !
  ! ********************************************************
  ! ********************************************************

!!  logical :: nonlinear = .false. !!Commented by AM 2016-02
  PetscScalar :: Delta = 4.5694d-3
  PetscScalar :: gamma = 1.0
  PetscScalar :: nu_n = 8.330d-3

  PetscScalar :: EParallelHat = 0
  PetscScalar :: dPhiHatdpsiHat = 0, dPhiHatdpsiN = 0, dPhiHatdrHat = 0, dPhiHatdrN = 0, Er = 0

  integer :: collisionOperator = 0
  ! 0 = Full linearized Fokker-Planck operator
  ! 1 = pitch-angle scattering with no momentum-conserving term

  logical :: includeXDotTerm = .true.
  logical :: includeElectricFieldTermInXiDot = .true.
  integer :: ExB_option = 1
  logical :: includeTemperatureEquilibrationTerm = .false.
  logical :: includePhi1 = .false.
!!  logical :: includeRadialExBDrive = .false. !!Commented by AM 2016-02
  logical :: includePhi1InKineticEquation = .true. !!Added by AM 2016-03

  PetscScalar :: nuPrime = 1, EStar = 0

  integer :: magneticDriftScheme = 0

  integer :: quasineutralityOption = 1

  PetscScalar :: ln_Lambda = 16
  !!!!!!!!!!!!!!!!!!!!!!!

  ! ********************************************************
  ! ********************************************************
  !
  ! Numerical input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: theta_derivative_option = 8
  integer :: zeta_derivative_option = 8
  integer :: xi_derivative_option = 8
  integer :: pitch_angle_scattering_option = 3
  integer :: xi_quadrature_option = 3

  !integer :: xDotDerivativeScheme = 0

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

  logical :: useIterativeLinearSolver=.true.

  integer :: whichParallelSolverToFactorPreconditioner = 1
  ! Options for whichParallelSolverToFactorPreconditioner:
  ! 1 = use mumps if it is detected, otherwise use superlu_dist
  ! 2 = force use of superlu_dist, if it is available

  integer :: preconditioner_x=1
  integer :: preconditioner_zeta_derivative_option=4
  integer :: preconditioner_theta_derivative_option = 4
  integer :: preconditioner_xi_derivative_option = 4
  integer :: preconditioner_species=1
  integer :: preconditioner_pitch_angle_scattering_option=2
  integer :: preconditioner_field_term_xi_option = 1 ! Should be 1 or 2 eventually.
  logical :: reusePreconditioner=.true.

  integer :: constraintScheme=-1

  integer :: PETSCPreallocationStrategy=1
  integer :: Nxi_for_x_option = 0
  integer :: x_scaling_option = 1
  integer :: spatial_scaling_option = 1
  logical :: fieldsplit

  ! ********************************************************
  !
  !  Other variables that are used by multiple subroutines
  !
  ! ********************************************************

  integer, dimension(:), allocatable :: Nxi_for_x, min_x_for_L
  integer :: matrixSize, NxPotentials
  PetscScalar, dimension(:), allocatable :: theta, zeta, x, x_plus1, xi
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights, xiWeights
  PetscScalar, dimension(:,:), allocatable, target :: ddtheta_plus, ddtheta_minus, ddtheta_plus_preconditioner, ddtheta_minus_preconditioner
  PetscScalar, dimension(:,:), allocatable, target :: ddzeta_plus, ddzeta_minus, ddzeta_plus_preconditioner, ddzeta_minus_preconditioner
  PetscScalar, dimension(:,:), allocatable, target :: ddxi_plus, ddxi_minus, ddxi_plus_preconditioner, ddxi_minus_preconditioner
  PetscScalar, dimension(:,:), allocatable, target :: pitch_angle_scattering_operator, pitch_angle_scattering_operator_preconditioner
  PetscScalar, dimension(:), allocatable :: xWeights, xPotentials
  PetscScalar :: maxxPotentials, zetaMax, xMaxNotTooSmall
  PetscScalar, dimension(:), allocatable :: x2, expx2
  PetscScalar, dimension(:,:), allocatable :: d2dx2, ddxPotentials, d2dx2Potentials
  PetscScalar, dimension(:,:), allocatable, target :: ddx, ddx_preconditioner
  PetscScalar, dimension(:,:), allocatable :: interpolateXToXPotentials

  PetscScalar, dimension(:,:,:,:,:), allocatable :: RosenbluthPotentialTerms
  PetscScalar, dimension(:,:), allocatable :: Legendre_polynomials
  PetscScalar, dimension(:,:,:), allocatable :: Legendre_projection

  integer, parameter :: COORDINATE_SYSTEM_UNINITIALIZED = -1
  integer, parameter :: COORDINATE_SYSTEM_BOOZER = 1
  integer, parameter :: COORDINATE_SYSTEM_VMEC = 2
  integer :: coordinateSystem = COORDINATE_SYSTEM_UNINITIALIZED

  PetscScalar, dimension(:,:), allocatable :: BHat, dBHatdtheta, dBHatdzeta, dBHatdpsiHat, sqrt_g, spatial_scaling, x_scaling
  PetscScalar, dimension(:,:), allocatable :: BHat_sub_psi, dBHat_sub_psi_dtheta, dBHat_sub_psi_dzeta
  PetscScalar, dimension(:,:), allocatable :: BHat_sub_theta, dBHat_sub_theta_dzeta, dBHat_sub_theta_dpsiHat
  PetscScalar, dimension(:,:), allocatable :: BHat_sub_zeta, dBHat_sub_zeta_dtheta, dBHat_sub_zeta_dpsiHat
  PetscScalar, dimension(:,:), allocatable :: BHat_sup_theta, dBHat_sup_theta_dzeta, dBHat_sup_theta_dpsiHat
  PetscScalar, dimension(:,:), allocatable :: BHat_sup_zeta, dBHat_sup_zeta_dtheta, dBHat_sup_zeta_dpsiHat
  PetscScalar, dimension(:,:), allocatable :: BDotCurlB, uHat, gradpsidotgradB_overgpsipsi
  PetscScalar, dimension(:,:), allocatable :: sources, jHat, Phi1Hat, dPhi1Hatdtheta, dPhi1Hatdzeta
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
  PetscScalar :: VPrimeHat, FSABHat2, FSABjHat, FSABjHatOverB0, FSABjHatOverRootFSAB2
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
  !  Variables related to parallelization:
  !
  ! ********************************************************

  MPI_Comm :: MPIComm
  integer :: numProcs, myRank 
  logical :: masterProc

  integer :: ithetaMin, ithetaMax, localNtheta
  integer :: izetaMin, izetaMax, localNzeta
  logical :: procThatHandlesConstraints

  Mat :: Mat_for_Jacobian

end module globalVariables

