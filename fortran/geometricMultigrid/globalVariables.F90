module globalVariables

  use kinds

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

  character(len=50), parameter :: inputFilename = "input.namelist"

  integer, parameter :: integerToRepresentTrue  =  1
  integer, parameter :: integerToRepresentFalse = -1

  real(prec), parameter :: one = 1., oneHalf = 0.5d+0
  real(prec), parameter :: zero = 0., two = 2., three = 3., four = 4., five = 5.

  real(prec), parameter :: pi = 3.14159265358979d+0
  real(prec), parameter :: sqrtpi = 1.77245385090552d+0

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
  real(prec) ::  GHat = 3.7481d+0, IHat = 0.0, iota = 0.4542d+0, B0OverBBar = 1.0, psiAHat = 0.15596d+0
  real(prec) ::  diotadpsiHat = 0.0, pPrimeHat = 0.0 !!Added by HS 2016-09-19
  real(prec) :: epsilon_t = -0.07053d+0, epsilon_h = 0.05067d+0, epsilon_antisymm = 0.0
  integer :: NPeriods = 0, helicity_l = 2, helicity_n = 10, helicity_antisymm_l = 1, helicity_antisymm_n = 0
  character(len=200) :: equilibriumFile = ""
  real(prec) :: min_Bmn_to_load = 0.0, dGdpHat, aHat = 0.5585d+0
  real(prec) :: psiHat_wish = -1, psiN_wish = 0.25, rHat_wish = -1, rN_wish = 0.5
  real(prec) :: psiHat, psiN, rHat, rN
  integer :: inputRadialCoordinate = 3
  integer :: inputRadialCoordinateForGradients = 4
  logical :: force0RadialCurrentInEquilibrium = .true.
  integer :: VMECRadialOption = 1
  real(prec) :: rippleScale = 1

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
  real(prec), dimension(maxNumSpecies) :: Zs, mHats, NHats, THats
  real(prec), dimension(maxNumSpecies) :: dNHatdpsiHats, dTHatdpsiHats
  real(prec), dimension(maxNumSpecies) :: dNHatdpsiNs,   dTHatdpsiNs
  real(prec), dimension(maxNumSpecies) :: dNHatdrHats,   dTHatdrHats
  real(prec), dimension(maxNumSpecies) :: dNHatdrNs,     dTHatdrNs

  !!Added by AM 2016-02!!
  real(prec) :: adiabaticZ = -1.0, adiabaticMHat = 5.446170214d-4, adiabaticNHat = 1.0, adiabaticTHat = 1.0
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
  real(prec) :: Delta = 4.5694d-3
  real(prec) :: gamma = 1.0
  real(prec) :: nu_n = 8.330d-3

  real(prec) :: EParallelHat = 0
  real(prec) :: dPhiHatdpsiHat = 0, dPhiHatdpsiN = 0, dPhiHatdrHat = 0, dPhiHatdrN = 0, Er = 0

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

  real(prec) :: nuPrime = 1, EStar = 0

  integer :: magneticDriftScheme = 0

  integer :: quasineutralityOption = 1

  real(prec) :: ln_Lambda = 16
  !!!!!!!!!!!!!!!!!!!!!!!

  ! ********************************************************
  ! ********************************************************
  !
  ! Numerical input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: theta_derivative_option = 10
  integer :: zeta_derivative_option = 8
  integer :: xi_derivative_option = 8
  integer :: pitch_angle_scattering_option = 3
  integer :: xi_quadrature_option = 3
  real(prec) :: theta_upwinding_factor = 0.2
  real(prec) :: zeta_upwinding_factor = 0.0

  !integer :: xDotDerivativeScheme = 0

  integer :: xPotentialsGridScheme = 2, xPotentialsInterpolationScheme
  integer :: xGridScheme = 5, xInterpolationScheme
  logical :: pointAtX0

  integer :: Ntheta = 15
  integer :: Nzeta = 15
  integer :: Nxi = 16
  integer :: NL = 4
  integer :: Nx = 5
  real(prec)  :: NxPotentialsPerVth = 40.0
  real(prec) :: xMax = 5.0
  real(prec) :: solverTolerance = 1d-6

  logical :: useIterativeLinearSolver=.true.

  integer :: whichParallelSolverToFactorPreconditioner = 1
  ! Options for whichParallelSolverToFactorPreconditioner:
  ! 1 = use mumps if it is detected, otherwise use superlu_dist
  ! 2 = force use of superlu_dist, if it is available

  integer :: preconditioning_option = 1
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
  integer :: spatial_scaling_option = 3
  integer :: constraint_scaling_option = 1
  integer :: null_space_option=0

  ! ********************************************************
  !
  !  Other variables that are used by multiple subroutines
  !
  ! ********************************************************

  integer :: matrixSize, NxPotentials
  real(prec), dimension(:), allocatable :: theta, zeta, x, x_plus1, xi
  real(prec), dimension(:), allocatable :: thetaWeights, zetaWeights, xiWeights
!  real(prec), dimension(:,:), allocatable, target :: ddtheta_plus, ddtheta_minus, ddtheta_plus_preconditioner, ddtheta_minus_preconditioner
!  real(prec), dimension(:,:), allocatable, target :: ddzeta_plus, ddzeta_minus, ddzeta_plus_preconditioner, ddzeta_minus_preconditioner
!  real(prec), dimension(:,:), allocatable, target :: ddxi_plus, ddxi_minus, ddxi_plus_preconditioner, ddxi_minus_preconditioner
!  real(prec), dimension(:,:), allocatable, target :: pitch_angle_scattering_operator, pitch_angle_scattering_operator_preconditioner
  real(prec), dimension(:), allocatable :: xWeights, xPotentials
  real(prec) :: maxxPotentials, zetaMax, xMaxNotTooSmall
  real(prec), dimension(:), allocatable :: x2, expx2
  real(prec), dimension(:,:), allocatable :: d2dx2, ddxPotentials, d2dx2Potentials
  real(prec), dimension(:,:), allocatable, target :: ddx, ddx_preconditioner
  real(prec), dimension(:,:), allocatable :: interpolateXToXPotentials

  ! This next array gets used in apply_dense_terms, so let's use type PetscScalar instead of real(prec).
  PetscScalar, dimension(:,:,:,:,:), allocatable :: RosenbluthPotentialTerms

  integer, parameter :: COORDINATE_SYSTEM_UNINITIALIZED = -1
  integer, parameter :: COORDINATE_SYSTEM_BOOZER = 1
  integer, parameter :: COORDINATE_SYSTEM_VMEC = 2
  integer :: coordinateSystem = COORDINATE_SYSTEM_UNINITIALIZED

  real(prec), dimension(:,:), allocatable :: x_scaling
  real(prec), dimension(:,:), allocatable :: sources, jHat, Phi1Hat, dPhi1Hatdtheta, dPhi1Hatdzeta
  real(prec), dimension(:,:,:), allocatable :: densityPerturbation, totalDensity
  real(prec), dimension(:,:,:), allocatable :: pressurePerturbation, totalPressure, pressureAnisotropy
  real(prec), dimension(:,:,:), allocatable :: flow, velocityUsingFSADensity, velocityUsingTotalDensity
  real(prec), dimension(:,:,:), allocatable :: MachUsingFSAThermalSpeed
  real(prec), dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vm0
  real(prec), dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vm
  real(prec), dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vE0
  real(prec), dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vE
  real(prec), dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vm0
  real(prec), dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vm
  real(prec), dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vE0
  real(prec), dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vE
  real(prec), dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vm0
  real(prec), dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vm
  real(prec), dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vE0
  real(prec), dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vE
  real(prec), dimension(:,:,:), allocatable :: NTVBeforeSurfaceIntegral
  real(prec), dimension(:,:), allocatable :: NTVKernel
  real(prec), dimension(:), allocatable :: FSADensityPerturbation, FSAPressurePerturbation
  real(prec), dimension(:), allocatable :: FSABFlow, FSABVelocityUsingFSADensity
  real(prec), dimension(:), allocatable :: FSABVelocityUsingFSADensityOverB0
  real(prec), dimension(:), allocatable :: FSABVelocityUsingFSADensityOverRootFSAB2

  real(prec), dimension(:), allocatable :: particleFlux_vm0_psiHat
  real(prec), dimension(:), allocatable :: particleFlux_vm0_psiN
  real(prec), dimension(:), allocatable :: particleFlux_vm0_rHat
  real(prec), dimension(:), allocatable :: particleFlux_vm0_rN

  real(prec), dimension(:), allocatable :: particleFlux_vm_psiHat
  real(prec), dimension(:), allocatable :: particleFlux_vm_psiN
  real(prec), dimension(:), allocatable :: particleFlux_vm_rHat
  real(prec), dimension(:), allocatable :: particleFlux_vm_rN

  real(prec), dimension(:), allocatable :: particleFlux_vE0_psiHat
  real(prec), dimension(:), allocatable :: particleFlux_vE0_psiN
  real(prec), dimension(:), allocatable :: particleFlux_vE0_rHat
  real(prec), dimension(:), allocatable :: particleFlux_vE0_rN

  real(prec), dimension(:), allocatable :: particleFlux_vE_psiHat
  real(prec), dimension(:), allocatable :: particleFlux_vE_psiN
  real(prec), dimension(:), allocatable :: particleFlux_vE_rHat
  real(prec), dimension(:), allocatable :: particleFlux_vE_rN

  real(prec), dimension(:), allocatable :: particleFlux_vd1_psiHat
  real(prec), dimension(:), allocatable :: particleFlux_vd1_psiN
  real(prec), dimension(:), allocatable :: particleFlux_vd1_rHat
  real(prec), dimension(:), allocatable :: particleFlux_vd1_rN

  real(prec), dimension(:), allocatable :: particleFlux_vd_psiHat
  real(prec), dimension(:), allocatable :: particleFlux_vd_psiN
  real(prec), dimension(:), allocatable :: particleFlux_vd_rHat
  real(prec), dimension(:), allocatable :: particleFlux_vd_rN

  real(prec), dimension(:), allocatable :: momentumFlux_vm0_psiHat
  real(prec), dimension(:), allocatable :: momentumFlux_vm0_psiN
  real(prec), dimension(:), allocatable :: momentumFlux_vm0_rHat
  real(prec), dimension(:), allocatable :: momentumFlux_vm0_rN

  real(prec), dimension(:), allocatable :: momentumFlux_vm_psiHat
  real(prec), dimension(:), allocatable :: momentumFlux_vm_psiN
  real(prec), dimension(:), allocatable :: momentumFlux_vm_rHat
  real(prec), dimension(:), allocatable :: momentumFlux_vm_rN

  real(prec), dimension(:), allocatable :: momentumFlux_vE0_psiHat
  real(prec), dimension(:), allocatable :: momentumFlux_vE0_psiN
  real(prec), dimension(:), allocatable :: momentumFlux_vE0_rHat
  real(prec), dimension(:), allocatable :: momentumFlux_vE0_rN

  real(prec), dimension(:), allocatable :: momentumFlux_vE_psiHat
  real(prec), dimension(:), allocatable :: momentumFlux_vE_psiN
  real(prec), dimension(:), allocatable :: momentumFlux_vE_rHat
  real(prec), dimension(:), allocatable :: momentumFlux_vE_rN

  real(prec), dimension(:), allocatable :: momentumFlux_vd1_psiHat
  real(prec), dimension(:), allocatable :: momentumFlux_vd1_psiN
  real(prec), dimension(:), allocatable :: momentumFlux_vd1_rHat
  real(prec), dimension(:), allocatable :: momentumFlux_vd1_rN

  real(prec), dimension(:), allocatable :: momentumFlux_vd_psiHat
  real(prec), dimension(:), allocatable :: momentumFlux_vd_psiN
  real(prec), dimension(:), allocatable :: momentumFlux_vd_rHat
  real(prec), dimension(:), allocatable :: momentumFlux_vd_rN

  real(prec), dimension(:), allocatable :: heatFlux_vm0_psiHat
  real(prec), dimension(:), allocatable :: heatFlux_vm0_psiN
  real(prec), dimension(:), allocatable :: heatFlux_vm0_rHat
  real(prec), dimension(:), allocatable :: heatFlux_vm0_rN

  real(prec), dimension(:), allocatable :: heatFlux_vm_psiHat
  real(prec), dimension(:), allocatable :: heatFlux_vm_psiN
  real(prec), dimension(:), allocatable :: heatFlux_vm_rHat
  real(prec), dimension(:), allocatable :: heatFlux_vm_rN

  real(prec), dimension(:), allocatable :: heatFlux_vE0_psiHat
  real(prec), dimension(:), allocatable :: heatFlux_vE0_rHat
  real(prec), dimension(:), allocatable :: heatFlux_vE0_rN
  real(prec), dimension(:), allocatable :: heatFlux_vE0_psiN

  real(prec), dimension(:), allocatable :: heatFlux_vE_psiHat
  real(prec), dimension(:), allocatable :: heatFlux_vE_rHat
  real(prec), dimension(:), allocatable :: heatFlux_vE_rN
  real(prec), dimension(:), allocatable :: heatFlux_vE_psiN

  real(prec), dimension(:), allocatable :: heatFlux_vd1_psiHat
  real(prec), dimension(:), allocatable :: heatFlux_vd1_psiN
  real(prec), dimension(:), allocatable :: heatFlux_vd1_rHat
  real(prec), dimension(:), allocatable :: heatFlux_vd1_rN

  real(prec), dimension(:), allocatable :: heatFlux_vd_psiHat
  real(prec), dimension(:), allocatable :: heatFlux_vd_psiN
  real(prec), dimension(:), allocatable :: heatFlux_vd_rHat
  real(prec), dimension(:), allocatable :: heatFlux_vd_rN

  real(prec), dimension(:), allocatable :: heatFlux_withoutPhi1_psiHat
  real(prec), dimension(:), allocatable :: heatFlux_withoutPhi1_psiN
  real(prec), dimension(:), allocatable :: heatFlux_withoutPhi1_rHat
  real(prec), dimension(:), allocatable :: heatFlux_withoutPhi1_rN

  real(prec), dimension(:,:), allocatable :: particleFlux_vm_psiHat_vs_x
  real(prec), dimension(:,:), allocatable :: heatFlux_vm_psiHat_vs_x
  real(prec), dimension(:,:), allocatable :: FSABFlow_vs_x

  real(prec), dimension(:), allocatable :: NTV
  real(prec) :: VPrimeHat, FSABHat2, FSABjHat, FSABjHatOverB0, FSABjHatOverRootFSAB2
  real(prec) :: lambda

  real(prec) :: ddpsiN2ddpsiHat, ddrHat2ddpsiHat, ddrN2ddpsiHat
  real(prec) :: ddpsiHat2ddpsiN, ddpsiHat2ddrHat, ddpsiHat2ddrN

  PetscLogDouble :: sfincs_start_time, iteration_start_time
  PetscLogDouble :: time_for_iteration
  integer :: number_of_KSP_iterations
  integer :: didLinearCalculationConverge, didNonlinearCalculationConverge
  integer :: iterationForMatrixOutput, iterationForResidualOutput = 0, iterationForStateVectorOutput = 0
  logical :: firstMatrixCreation

  integer :: transportMatrixSize = 3
  real(prec), dimension(:,:), allocatable :: transportMatrix

  Vec :: f0

  KSP :: inner_KSP
  PC :: inner_preconditioner

  ! ********************************************************
  !
  !  Variables related to parallelization:
  !
  ! ********************************************************

  MPI_Comm :: MPIComm
  integer :: numProcs, myRank 
  logical :: masterProc
  logical :: procThatHandlesConstraints

  ! ********************************************************
  !
  ! Variables related to multigrid
  !
  ! ********************************************************

  type :: multigrid_level
     integer :: Ntheta, Nzeta, Nxi, matrixSize

     real(prec), allocatable, dimension(:) :: theta, zeta, xi, y
     real(prec), dimension(:), allocatable :: thetaWeights, zetaWeights, xiWeights

     real(prec), dimension(:,:), allocatable :: ddtheta_plus
     real(prec), dimension(:,:), allocatable :: ddtheta_minus
     real(prec), dimension(:,:), allocatable :: ddtheta_plus_preconditioner
     real(prec), dimension(:,:), allocatable :: ddtheta_minus_preconditioner
     
     real(prec), dimension(:,:), allocatable :: ddzeta_plus
     real(prec), dimension(:,:), allocatable :: ddzeta_minus
     real(prec), dimension(:,:), allocatable :: ddzeta_plus_preconditioner
     real(prec), dimension(:,:), allocatable :: ddzeta_minus_preconditioner

     real(prec), dimension(:,:), allocatable :: ddxi_plus
     real(prec), dimension(:,:), allocatable :: ddxi_minus
     real(prec), dimension(:,:), allocatable :: ddxi_plus_preconditioner
     real(prec), dimension(:,:), allocatable :: ddxi_minus_preconditioner

     real(prec), dimension(:,:), allocatable :: pitch_angle_scattering_operator
     real(prec), dimension(:,:), allocatable :: pitch_angle_scattering_operator_preconditioner

     real(prec), dimension(:,:), allocatable :: BHat, dBHatdtheta, dBHatdzeta, dBHatdpsiHat, sqrt_g
     real(prec), dimension(:,:), allocatable :: BHat_sub_psi, dBHat_sub_psi_dtheta, dBHat_sub_psi_dzeta
     real(prec), dimension(:,:), allocatable :: BHat_sub_theta, dBHat_sub_theta_dzeta, dBHat_sub_theta_dpsiHat
     real(prec), dimension(:,:), allocatable :: BHat_sub_zeta, dBHat_sub_zeta_dtheta, dBHat_sub_zeta_dpsiHat
     real(prec), dimension(:,:), allocatable :: BHat_sup_theta, dBHat_sup_theta_dzeta, dBHat_sup_theta_dpsiHat
     real(prec), dimension(:,:), allocatable :: BHat_sup_zeta, dBHat_sup_zeta_dtheta, dBHat_sup_zeta_dpsiHat
     real(prec), dimension(:,:), allocatable :: BDotCurlB, uHat, gradpsidotgradB_overgpsipsi
     ! This next array gets used in apply_dense_terms, so let's use type PetscScalar instead of real(prec).
     PetscScalar, dimension(:,:,:), allocatable :: Legendre_projection

     integer :: ithetaMin, ithetaMax
     integer :: izetaMin, izetaMax

     integer, dimension(:), allocatable :: Nxi_for_x, min_x_for_L
     integer :: DKE_size
     integer, dimension(:), allocatable :: first_index_for_x
     real(prec), dimension(:,:), allocatable :: spatial_scaling

     Mat :: low_order_matrix, high_order_matrix, smoothing_matrix

  end type multigrid_level

  type (multigrid_level), allocatable, dimension(:), target :: levels

  Mat, allocatable, dimension(:) :: multigrid_prolongation_matrices, multigrid_restriction_matrices
  integer :: smoothing_option = 1
  integer :: coarsen_option = 1
  integer :: defect_option = 2
  logical :: coarsen_theta = .false.
  logical :: coarsen_zeta = .true.
  logical :: coarsen_xi = .false.
  integer :: N_levels
  integer, allocatable, dimension(:) :: Ntheta_levels, Nzeta_levels, Nxi_levels
  integer :: Ntheta_min = 7
  integer :: Nzeta_min = 7
  integer :: Nxi_min = 9

end module globalVariables

