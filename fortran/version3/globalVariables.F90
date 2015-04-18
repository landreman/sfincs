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
  logical :: isAParallelDirectSolverInstalled

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

  integer :: xPotentialsGridScheme = 2, xPotentialsInterpolationScheme
  integer :: xGridScheme = 1, xInterpolationScheme
  logical :: pointAtX0

  integer :: Ntheta
  integer :: Nzeta
  integer :: Nxi
  integer :: NL
  integer :: Nx 
  PetscScalar  :: NxPotentialsPerVth
  PetscScalar :: xMax
  PetscScalar :: solverTolerance

  logical :: forceOddNthetaAndNzeta=.true.
  ! If forceOddNthetaAndNzeta is set to true, 1 is added to Ntheta any time a run is attempted with even Ntheta,
  ! and 1 is added to Nzeta any time a run is attempted with even Nzeta.
  ! This can be useful because the iterative solvers sometimes do not work with even Ntheta or Nzeta.
  ! This parameter should be true unless you know what you are doing.

  logical :: useIterativeLinearSolver=.true.

  integer :: whichParallelSolverToFactorPreconditioner
  ! Options for whichParallelSolverToFactorPreconditioner:
  ! 1 = use mumps if it is detected, otherwise use superlu_dist
  ! 2 = force use of superlu_dist, if it is available

  integer :: preconditioner_x=1, preconditioner_x_min_L=0, preconditioner_zeta=0
  integer :: preconditioner_theta=0, preconditioner_xi=0, preconditioner_species=1
  integer :: preconditioner_theta_min_L=0, preconditioner_zeta_min_L=0
  logical :: reusePreconditioner=.true.

  integer :: constraintScheme=-1

  integer :: PETSCPreallocationStrategy

  ! ********************************************************
  !
  !  Other variables that are used by multiple subroutines
  !
  ! ********************************************************

  integer :: matrixSize, NxPotentials
  PetscScalar, dimension(:), allocatable :: theta, zeta, x, x_plus1
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights
  PetscScalar, dimension(:,:), allocatable :: ddtheta, ddzeta
  PetscScalar, dimension(:), allocatable :: xWeights, xPotentials
  PetscScalar :: maxxPotentials, zetaMax, xMaxNotTooSmall
  PetscScalar, dimension(:), allocatable :: x2, expx2
  PetscScalar, dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
  PetscScalar, dimension(:,:), allocatable :: ddx_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddtheta_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddzeta_preconditioner
  PetscScalar, dimension(:,:), allocatable :: interpolateXToXPotentials

  PetscScalar, dimension(:,:), allocatable :: ddtheta_upwind_diagonal
  PetscScalar, dimension(:,:), allocatable :: ddtheta_upwind_offDiagonal
  PetscScalar, dimension(:,:), allocatable :: ddzeta_upwind_diagonal
  PetscScalar, dimension(:,:), allocatable :: ddzeta_upwind_offDiagonal
  ! First 100 values computed using m20150417_01_testingLegendreIntegrals.m
  PetscScalar, dimension(100) :: LegendreIntegrals_diagonal = &
     (/                0.5,&
                      0.25,&
                     0.125,&
                   0.09375,&
                 0.0703125,&
                0.05859375,&
               0.048828125,&
            0.042724609375,&
         0.037384033203125,&
        0.0336456298828125,&
        0.0302810668945312,&
        0.0277576446533203,&
        0.0254445075988769,&
        0.0236270427703857,&
        0.0219393968582153,&
        0.0205681845545769,&
        0.0192826730199158,&
        0.0182114134076983,&
        0.0171996682183817,&
        0.0163396848074626,&
        0.0155227005670895,&
        0.0148171232685854,&
        0.0141436176654679,&
        0.0135543002627401,&
        0.0129895377517926,&
        0.0124899401459544,&
        0.0120095578326485,&
         0.011580645052911,&
        0.0111670505867356,&
        0.0107948155671778,&
        0.0104349883816052,&
          0.01010889499468,&
       0.00979299202609626,&
       0.00950496284885815,&
       0.00922540511800936,&
       0.00896914386473133,&
       0.00872000097959992,&
       0.00849052726961045,&
       0.00826709234146281,&
       0.00806041503292624,&
       0.00785890465710308,&
       0.00767178787955299,&
       0.00748912626337315,&
       0.00731891884829652,&
       0.00715257978356249,&
       0.00699708891870245,&
       0.00684497829003499,&
       0.00670237457565927,&
       0.00656274177199971,&
        0.0064314869365597,&
       0.00630285719782852,&
       0.00618164840556258,&
       0.00606277055160946,&
       0.00595049702287595,&
       0.00584030263356343,&
       0.00573601151510694,&
       0.00563358273805146,&
        0.0055364520011885,&
       0.00544099593220249,&
       0.00535031266666578,&
       0.00526114078888802,&
       0.00517628367938983,&
       0.00509279523294806,&
       0.00501322030743325,&
        0.0049348887401296,&
       0.00486011769861249,&
       0.00478647955166382,&
        0.0047160901464923,&
       0.00464673587963213,&
        0.0045803539384945,&
       0.00451492031080174,&
       0.00445221308426282,&
       0.00439037679142584,&
       0.00433104737532549,&
       0.00427251970809136,&
       0.00421630234351121,&
       0.00416082468109659,&
       0.00410748077492868,&
        0.0040548207649937,&
       0.00400413550543128,&
       0.00395408381161339,&
       0.00390586327732542,&
       0.00385823079833364,&
       0.00381229947930587,&
       0.00376691496169507,&
       0.00372311362493116,&
       0.00367982160603661,&
       0.00363800545142257,&
       0.00359666448038368,&
       0.00355670154171275,&
       0.00351718263569372,&
       0.00347895238965357,&
       0.00344113768976604,&
       0.00340452984200257,&
       0.00336831143942807,&
       0.00333322486193403,&
       0.00329850376962222,&
       0.00326484556789138,&
       0.00323153081719861,&
       0.00319921550902662 /)


  PetscScalar, dimension(:,:,:,:,:), allocatable :: Rosenbluth_H
  PetscScalar, dimension(:,:,:,:,:), allocatable :: Rosenbluth_dHdxb
  PetscScalar, dimension(:,:,:,:,:), allocatable :: Rosenbluth_d2Gdxb2

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

