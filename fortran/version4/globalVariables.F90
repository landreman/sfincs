module globalVariables

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

  ! Double precision, for calculations that are always done in double even when PETSc is compiled with single-precision
  integer, parameter :: prec = SELECTED_REAL_KIND(12,100)

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

  real(prec) :: adiabaticZ = -1.0, adiabaticMHat = 5.446170214d-4, adiabaticNHat = 1.0, adiabaticTHat = 1.0
  logical :: withAdiabatic = .false.


  ! ********************************************************
  ! ********************************************************
  !
  ! Physics input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  real(prec) :: Delta = 4.5694d-3
  real(prec) :: alpha = 1.0
  real(prec) :: nu_n = 8.330d-3

  real(prec) :: EParallelHat = 0
  real(prec) :: dPhiHatdpsiHat = 0, dPhiHatdpsiN = 0, dPhiHatdrHat = 0, dPhiHatdrN = 0, Er = 0

  integer :: collisionOperator = 0
  ! 0 = Full linearized Fokker-Planck operator
  ! 1 = pitch-angle scattering with no momentum-conserving term

  logical :: includeXDotTerm = .true.
  logical :: includeElectricFieldTermInXiDot = .true.
  logical :: useDKESExBDrift = .false.
  logical :: includeTemperatureEquilibrationTerm = .false.
  logical :: includePhi1 = .false.
  logical :: includePhi1InKineticEquation = .true. !!Added by AM 2016-03

  real(prec) :: nuPrime = 1, EStar = 0

  integer :: magneticDriftScheme = 0

  !!Added by AM 2016-02!!
  integer :: quasineutralityOption = 1
  !!!!!!!!!!!!!!!!!!!!!!!

  ! ********************************************************
  ! ********************************************************
  !
  ! Numerical input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: xPotentialsGridScheme = 2, xPotentialsInterpolationScheme
  integer :: xGridScheme = 5, xInterpolationScheme
  logical :: pointAtX0

  integer :: NFourier = 10
  integer :: NFourier2=0
  integer :: mmax=32, nmax=32
  integer :: Ntheta, Nzeta
  integer :: Nxi = 16
  integer :: NL = 4
  integer :: Nx = 5
  real(prec)  :: NxPotentialsPerVth = 40.0
  real(prec) :: xMax = 5.0
  real(prec) :: solverTolerance = 1d-6
  real(prec) :: Fourier_threshold = 1d-12


  integer :: whichParallelSolverToFactorPreconditioner = 1
  ! Options for whichParallelSolverToFactorPreconditioner:
  ! 1 = use mumps if it is detected, otherwise use superlu_dist
  ! 2 = force use of superlu_dist, if it is available

  integer :: preconditioner_x=1, preconditioner_x_min_L=0
  integer :: preconditioner_Fourier=0, preconditioner_xi=1, preconditioner_species=1
  integer :: preconditioner_Fourier_min_L=0
  real(prec) :: preconditioner_Fourier_threshold=0.03
  logical :: reusePreconditioner=.true.
  integer :: preconditioner_magnetic_drifts_max_L = 2
  logical :: preconditioner_drop_xDot = .true.
  logical :: preconditioner_drop_xiDot = .true.
  integer :: preconditioner_Fourier_max_nnz_per_row = 6

  integer :: constraintScheme=-1

  integer :: PETSCPreallocationStrategy=1
  integer :: Nxi_vs_x_option = 1

  ! ********************************************************
  !
  !  Other variables that are used by multiple subroutines
  !
  ! ********************************************************

  integer, dimension(:), allocatable :: Nxi_for_x, min_x_for_L
  integer, dimension(:), allocatable :: xm, xn
  real(prec), dimension(:,:), allocatable :: predictedAmplitudes
  integer :: matrixSize, NxPotentials
  real(prec), dimension(:), allocatable :: theta, zeta, x, x_plus1
  real(prec) :: thetaWeight, zetaWeight
  real(prec), dimension(:,:), allocatable :: ddtheta, ddzeta
  real(prec), dimension(:), allocatable :: xWeights, xPotentials
  real(prec) :: maxxPotentials, zetaMax, xMaxNotTooSmall
  real(prec), dimension(:), allocatable :: x2, expx2
  real(prec), dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
  real(prec), dimension(:,:), allocatable :: ddx_preconditioner
  real(prec), dimension(:,:), allocatable :: interpolateXToXPotentials

  real(prec), dimension(:,:,:,:,:), allocatable :: Rosenbluth_H
  real(prec), dimension(:,:,:,:,:), allocatable :: Rosenbluth_dHdxb
  real(prec), dimension(:,:,:,:,:), allocatable :: Rosenbluth_d2Gdxb2

  integer, parameter :: COORDINATE_SYSTEM_UNINITIALIZED = -1
  integer, parameter :: COORDINATE_SYSTEM_BOOZER = 1
  integer, parameter :: COORDINATE_SYSTEM_VMEC = 2
  integer :: coordinateSystem = COORDINATE_SYSTEM_UNINITIALIZED

  real(prec), dimension(:,:), allocatable :: BHat, dBHatdtheta, dBHatdzeta, dBHatdpsiHat, DHat
  real(prec), dimension(:,:), allocatable :: BHat_sub_psi, dBHat_sub_psi_dtheta, dBHat_sub_psi_dzeta
  real(prec), dimension(:,:), allocatable :: BHat_sub_theta, dBHat_sub_theta_dzeta, dBHat_sub_theta_dpsiHat
  real(prec), dimension(:,:), allocatable :: BHat_sub_zeta, dBHat_sub_zeta_dtheta, dBHat_sub_zeta_dpsiHat
  real(prec), dimension(:,:), allocatable :: BHat_sup_theta, dBHat_sup_theta_dzeta, dBHat_sup_theta_dpsiHat
  real(prec), dimension(:,:), allocatable :: BHat_sup_zeta, dBHat_sup_zeta_dtheta, dBHat_sup_zeta_dpsiHat
  real(prec), dimension(:,:), allocatable :: BDotCurlB, uHat

  real(prec), dimension(:,:), allocatable :: sources, jHat_realSpace
  real(prec), dimension(:), allocatable :: jHat_Fourier
  real(prec), dimension(:), allocatable :: Phi1Hat_Fourier, dPhi1Hatdtheta_Fourier, dPhi1Hatdzeta_Fourier
  real(prec), dimension(:,:), allocatable :: Phi1Hat_realSpace, dPhi1Hatdtheta_realSpace, dPhi1Hatdzeta_realSpace

  real(prec), dimension(:,:,:), allocatable :: densityNonadiabaticPerturbation_realSpace, totalDensity_realSpace
  real(prec), dimension(:,:,:), allocatable :: pressureNonadiabaticPerturbation_realSpace, totalPressure_realSpace
  real(prec), dimension(:,:,:), allocatable :: pressureAnisotropy_realSpace
  real(prec), dimension(:,:,:), allocatable :: flow_realSpace, velocityUsingFSADensity_realSpace
  real(prec), dimension(:,:,:), allocatable :: MachUsingFSAThermalSpeed_realSpace
  real(prec), dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vm0_realSpace
  real(prec), dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vm_realSpace
  real(prec), dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vE0_realSpace
  real(prec), dimension(:,:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vE_realSpace
  real(prec), dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vm0_realSpace
  real(prec), dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vm_realSpace
  real(prec), dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vE0_realSpace
  real(prec), dimension(:,:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vE_realSpace
  real(prec), dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vm0_realSpace
  real(prec), dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vm_realSpace
  real(prec), dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vE0_realSpace
  real(prec), dimension(:,:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vE_realSpace

  real(prec), dimension(:,:), allocatable :: densityNonadiabaticPerturbation_Fourier, totalDensity_Fourier
  real(prec), dimension(:,:), allocatable :: pressureNonadiabaticPerturbation_Fourier, totalPressure_Fourier
  real(prec), dimension(:,:), allocatable :: pressureAnisotropy_Fourier
  real(prec), dimension(:,:), allocatable :: flow_Fourier, velocityUsingFSADensity_Fourier
  real(prec), dimension(:,:), allocatable :: MachUsingFSAThermalSpeed_Fourier
  real(prec), dimension(:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vm0_Fourier
  real(prec), dimension(:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vm_Fourier
  real(prec), dimension(:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vE0_Fourier
  real(prec), dimension(:,:), allocatable :: particleFluxBeforeSurfaceIntegral_vE_Fourier
  real(prec), dimension(:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vm0_Fourier
  real(prec), dimension(:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vm_Fourier
  real(prec), dimension(:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vE0_Fourier
  real(prec), dimension(:,:), allocatable :: momentumFluxBeforeSurfaceIntegral_vE_Fourier
  real(prec), dimension(:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vm0_Fourier
  real(prec), dimension(:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vm_Fourier
  real(prec), dimension(:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vE0_Fourier
  real(prec), dimension(:,:), allocatable :: heatFluxBeforeSurfaceIntegral_vE_Fourier

  real(prec), dimension(:,:,:), allocatable :: NTVBeforeSurfaceIntegral
  real(prec), dimension(:,:), allocatable :: NTVKernel

  real(prec), dimension(:), allocatable :: FSADensityNonadiabaticPerturbation, FSAPressureNonadiabaticPerturbation
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

  real(prec) :: elapsedTime
  PetscLogDouble :: startTime
  integer :: didLinearCalculationConverge, didNonlinearCalculationConverge
  integer :: iterationForMatrixOutput, iterationForResidualOutput = 0, iterationForStateVectorOutput = 0
  logical :: firstMatrixCreation

  integer :: transportMatrixSize = 3
  real(prec), dimension(:,:), allocatable :: transportMatrix

  Vec :: f0

  real(prec), dimension(:,:), allocatable :: FourierMatrix_streaming
  real(prec), dimension(:,:), allocatable :: FourierMatrix_ExB
  real(prec), dimension(:,:), allocatable :: FourierMatrix_mirror
  real(prec), dimension(:,:), allocatable :: FourierMatrix_xiDot
  real(prec), dimension(:,:), allocatable :: FourierMatrix_xDot


  ! ********************************************************
  !
  !  Variables related to parallelization:
  !
  ! ********************************************************

  MPI_Comm :: MPIComm
  integer :: numProcs, myRank 
  logical :: masterProc

  !integer :: iFourierMin, iFourierMax, localNFourier
  integer :: LMin, LMax, localNxi
  logical :: procThatHandlesConstraints
  Mat :: MatForJacobian

end module globalVariables

