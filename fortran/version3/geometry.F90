! Note: all the "wish" radial coordinates (psiHat_wish, psiN_wish, rHat_wish, and rN_wish)
! are consistently set to an equivalent radius before any subroutines in this module are run.  
! However, the computeBHat subroutine in this module must set the final rN.  The other 3
! final radial coordinates (psiHat, psiN, and rHat) will be set from rN later on.

! Reading in the magnetic geometry takes a negligible amount of time, so the code here
! just has every processor independently open and read the input file.  An alternative
! would be to have one processor read the file and then broadcast the geometry info to the others.
! Since computing time is not an issue for this part of the code, there is little advantage to the
! latter approach.

module geometry

#include "PETScVersions.F90"

  use globalVariables
  use radialCoordinates
  use read_wout_mod, only: read_wout_file, Aminor, phi, nfp, ns, xm, xn, xm_nyq, xn_nyq, mpol, ntor, mnmax, mnmax_nyq, &
       lasym, presf, phip, iotas, bmnc, bmns, gmnc, gmns, &
       bsubumnc, bsubumns, bsubvmnc, bsubvmns, bsubsmnc, bsubsmns, bsupumnc, bsupumns, bsupvmnc, bsupvmns, rmnc, rmns, zmnc, zmns

  implicit none

contains

  ! -----------------------------------------------------------------------------------

  subroutine initializeGeometry()
    ! For each geometryScheme, this subroutine must set the following variables:
    !   NPeriods
    !   psiAHat (if the value in input.namelist is to be over-written.)
    !   aHat (if the value in input.namelist is to be over-written.)
    ! 
    ! This subroutine sets NPeriods, which is the number of identical toroidal segments
    ! in the stellarator (e.g. 5 for W7-X, 10 for LHD, 4 for HSX.)
    ! Also, if psiAHat and/or aHat from input.namelist are going to be over-written,
    ! this is done now.  (We need to set psiAHat and aHat early on so the "wish" radius
    ! can be set from any of the radial coordinates.)

    implicit none

    integer :: fileUnit, didFileAccessWork, Turkin_sign, ierr, iopen
    character(len=200) :: lineOfFile
    integer, dimension(4) :: headerIntegers
    PetscScalar, dimension(3) :: headerReals
    integer :: tag, dummy(1), i
    integer :: status(MPI_STATUS_SIZE)

    select case (geometryScheme)
    case (1)
       NPeriods = max(1, helicity_n)

    case (2)
       NPeriods = 10
       aHat = 0.5585d+0 ! (meters)
       psiAHat = (aHat ** 2) / two
       rN_wish = 0.5d+0
       inputRadialCoordinate = 3

       if (masterProc) then
          print *,"---------------------------------------------------------"
          print *,"Since geometryScheme=2, we will ignore the *_wish parameters and use the flux surface rN = 0.5."
       end if

    case (3)
       NPeriods = 10
       aHat = 0.5400d+0 ! (meters)
       psiAHat = (aHat ** 2) / two
       rN_wish = 0.5d+0
       inputRadialCoordinate = 3

       if (masterProc) then
          print *,"---------------------------------------------------------"
          print *,"Since geometryScheme=3, we will ignore the *_wish parameters and use the flux surface rN = 0.5."
       end if

    case (4)
       NPeriods = 5
       aHat = 0.5109d+0 ! (meters)
       psiAHat = -0.384935d+0 ! Tesla * meters^2 / radian
       rN_wish = 0.5d+0
       inputRadialCoordinate = 3

       if (masterProc) then
          print *,"---------------------------------------------------------"
          print *,"Since geometryScheme=4, we will ignore the *_wish parameters and use the flux surface rN = 0.5."
       end if

    case (5)
       ! Read VMEC file, defining the effective minor radius aHat to be the quantity called Aminor_p in vmec's wout file, which is called just Aminor in read_wout_mod.F.
       ! Libstell does not allow multiple procs to open an ASCII wout file simultaneously, so have each proc read it 1 at a time.
       if (masterProc) then
          call read_wout_file(equilibriumfile, ierr, iopen)
          if (iopen .ne. 0) then
             print *, 'Error opening wout file:',equilibriumFile
             stop
          end if
          if (ierr .ne. 0) then
             print *, 'Error reading wout file:',equilibriumFile
             stop
          end if

          tag=0
          dummy=0
          do i = 1,numProcs-1
             ! Ping each proc 1 at a time by sending a dummy value:
             call MPI_SEND(dummy,1,MPI_INT,i,tag,MPIComm,ierr)
             ! Wait for a value to be returned before continuing.
             call MPI_RECV(dummy,1,MPI_INT,i,MPI_ANY_TAG,MPIComm,status,ierr)
          end do
       else
          ! First, wait for the dummy message from proc 0:
          call MPI_RECV(dummy,1,MPI_INT,0,MPI_ANY_TAG,MPIComm,status,ierr)
          call read_wout_file(equilibriumfile, ierr, iopen)
          if (iopen .ne. 0) then
             print *, 'Error opening wout file:',equilibriumFile
             stop
          end if
          if (ierr .ne. 0) then
             print *, 'Error reading wout file:',equilibriumFile
             stop
          end if
          ! Send a value back to proc 0 to tell it to continue.
          call MPI_SEND(dummy,1,MPI_INT,0,tag,MPIComm,ierr)
       end if

       NPeriods = nfp
       psiAHat = phi(ns)/(2*pi)
       aHat = Aminor

!!$    case (6)
!!$       ! Read Erika Strumberger's NEMEC format used at IPP.
!!$       call read_NEMEC(equilibriumFile)
!!$       NPeriods = vmec%nfp
!!$       psiAHat = vmec%phi(vmec%ns)/(2*pi)
!!$       ! Aminor_p is not set in this format, so we must specify a value in the input namelist.
!!$       if (aHat-0.5585d+0 < 1e-4 .and. masterProc) then
!!$          print *,"###############################################################################################"
!!$          print *,"WARNING: It appears you have not explicitly set aHat to a value other than the default."
!!$          print *,"This is probably wrong. For geometryScheme=6, aHat is not set using the NEMEC equilibrium file."
!!$          print *,"###############################################################################################"
!!$       end if

    case (10)
       print *,"Error! This geometryScheme has not been implemented yet."

    case (11)
       fileUnit = 11
       open(unit=fileUnit, file=equilibriumFile, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file."
          stop
       else
          Turkin_sign = 1
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             ! Skip lines that begin with "CC":
             if (lineOfFile(1:2) /= "CC") then
                exit
             else
                if (index(lineOfFile,"CStconfig")>0) then
                   Turkin_sign = -1 ! This file was saved by Yuriy Turkin.
                                    ! He has then changed a sign, which I here change back to Joachim Geiger's convention
                end if                
             end if
          end do

          ! Read header line:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) headerIntegers, headerReals
          if (didFileAccessWork /= 0) then
             print *,"Unable to read header from the magnetic equilibrium file ",equilibriumFile
             stop
          end if

          NPeriods = headerIntegers(4)
          psiAHat  = headerReals(1)/2/pi !Convert the flux from Tm^2 to Tm^2/rad
          aHat     = headerReals(2)      !minor radius in meters

       end if

       close(unit = fileUnit)
       if (masterProc) then
          print *,"Successfully opened magnetic equilibrium file ",trim(equilibriumFile),".  Nperiods = ",Nperiods
       end if

       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       psiAHat=psiAHat*(-1)*Turkin_sign           !toroidal direction sign switch, and possibly Yuriy Turkin's sign

    case (12)
       fileUnit = 11
       open(unit=fileUnit, file=equilibriumFile, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file."
          stop
       else
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             ! Skip lines that begin with "CC":
             if (lineOfFile(1:2) /= "CC") exit
          end do

          ! Read header line:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) headerIntegers, headerReals
          if (didFileAccessWork /= 0) then
             print *,"Unable to read header from the magnetic equilibrium file ",equilibriumFile
             stop
          end if

          NPeriods = headerIntegers(4)
          psiAHat  = headerReals(1)/2/pi !Convert the flux from Tm^2 to Tm^2/rad
          aHat     = headerReals(2)      !minor radius in meters

       end if

       close(unit = fileUnit)
       if (masterProc) then
          print *,"Successfully opened magnetic equilibrium file ",trim(equilibriumFile),".  Nperiods = ",Nperiods
       end if

       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       psiAHat=psiAHat*(-1)           !toroidal direction sign switch    

    case default
       print *,"Error! Invalid setting for geometryScheme."
       stop
    end select

  end subroutine initializeGeometry

  ! -----------------------------------------------------------------------------------

  subroutine computeBHat()
! This subroutine initializes a bunch of variables, including a bunch of (theta,zeta) arrays
! that store |B| and the components of \vect{B}, as well as derivatives of these quantities.
!
! For Boozer coordinates, you must set the following variables:
!  BHat, dBHatdpsiHat, dBHatdtheta, dBHatdzeta,
!  IHat, GHat, iota,
!  dIHatdpsiHat, dGHatdpsiHat, diotadpsiHat,
!  BHat_sub_psi, dBHat_sub_psi_dtheta, dBHat_sub_psi_dzeta.
! Then call setBoozerCoordinates(), which is a subroutine that sets the other arrays.
!
! For coordinates other than Boozer coordinates, you must set all the following arrays directly:
!  BHat, dBHatdtheta, dBHatdzeta, dBHatdpsiHat, DHat
!  BHat_sub_psi, dBHat_sub_psi_dtheta, dBHat_sub_psi_dzeta
!  BHat_sub_theta, dBHat_sub_theta_dzeta, dBHat_sub_theta_dpsiHat
!  BHat_sub_zeta, dBHat_sub_zeta_dtheta, dBHat_sub_zeta_dpsiHat
!  BHat_sup_theta, dBHat_sup_theta_dzeta, dBHat_sup_theta_dpsiHat
!  BHat_sup_zeta, dBHat_sup_zeta_dtheta, dBHat_sup_zeta_dpsiHat


! Note: all the "wish" radial coordinates (psiHat_wish, psiN_wish, rHat_wish, and rN_wish)
! are consistently set to an equivalent radius before any subroutines in this module are run.  
! However, the computeBHat subroutine in this module must set the final rN.  The other 3
! final radial coordinates (psiHat, psiN, and rHat) will be set from rN later on.

    implicit none

    integer :: itheta,izeta
    PetscScalar :: DHat11, BHat_sub_theta11, BHat_sub_zeta11

    ! Using the selected radial coordinate, set input quantities for the other radial coordinates:
    call setInputRadialCoordinateWish()
    ! Note that this call only sets the "wish" radial coordinates, not the final radial coordinates
    ! or the input gradients. These quantities will be set later as we load the magnetic
    ! geometry, in case the final radial coordinate is different from the "wish" values.

    ! Set the radius to a silly value here to make sure the proper value is set eventually:
    rN = -9999

    select case (geometryScheme)
    case (1,2,3,4,11,12)
       coordinateSystem = COORDINATE_SYSTEM_BOOZER
       call computeBHat_Boozer()
    case (5,6,7)
       coordinateSystem = COORDINATE_SYSTEM_VMEC
       call computeBHat_VMEC()
    case default
       print *,"Error! Invalid setting for geometryScheme."
       stop
    end select

    BDotCurlB = DHat * (  BHat_sub_theta * dBHat_sub_psi_dzeta &
                        - BHat_sub_theta * dBHat_sub_zeta_dpsiHat &
                        + BHat_sub_zeta * dBHat_sub_theta_dpsiHat &
                        - BHat_sub_zeta * dBHat_sub_psi_dtheta)

    if (.not. force0RadialCurrentInEquilibrium) then
       BDotCurlB = BDotCurlB + DHat * BHat_sub_psi * (dBHat_sub_zeta_dtheta - dBHat_sub_theta_dzeta)
    end if

    ! Validate geometry arrays:
    DHat11 = DHat(1,1)
    BHat_sub_theta11 = BHat_sub_theta(1,1)
    BHat_sub_zeta11 = BHat_sub_zeta(1,1)
    do itheta=1,Ntheta
       do izeta=1,Nzeta
          if (BHat(itheta,izeta) <= 0) then
             print *,"ERROR! BHat is not everywhere positive!"
             stop
          end if
          if (DHat(itheta,izeta)*DHat11 <= 0) then
             print *,"ERROR! DHat does not have the same sign everywhere!"
             stop
          end if
          if (ExBDerivativeSchemeZeta > 0 .and. BHat_sub_theta(itheta,izeta)*BHat_sub_theta11 < -1d-15) then
             print *,"ERROR! ExBDerivativeSchemeZeta>0 assumes BHat_sub_theta has the same sign everywhere,"
             print *,"       but the sign of BHat_sub_theta is not the same everywhere."
             stop
          end if
          if (ExBDerivativeSchemeTheta > 0 .and. BHat_sub_zeta(itheta,izeta)*BHat_sub_zeta11 < -1d-15) then
             print *,"ERROR! ExBDerivativeSchemeTheta>0 assumes BHat_sub_zeta has the same sign everywhere,"
             print *,"       but the sign of BHat_sub_zeta is not the same everywhere."
             stop
          end if
       end do
    end do
  end subroutine computeBHat

  ! ---------------------------------------------------------------------------------------
          
  subroutine computeBHat_Boozer

    use transforms
    ! Note that the BHarmonics_amplitudes in this subroutine are normalized by B0, not by BBar!

    implicit none

    integer :: itheta, izeta, NHarmonics, NHarmonicsL, NHarmonicsH, i, m, n
    integer, dimension(:), allocatable :: BHarmonics_l, BHarmonics_n
    integer, dimension(:), allocatable :: BHarmonics_lL, BHarmonics_nL, BHarmonics_lH, BHarmonics_nH
    PetscScalar, dimension(:), allocatable :: BHarmonics_amplitudes
    PetscScalar, dimension(:), allocatable :: BHarmonics_amplitudesL, BHarmonics_amplitudesH
    PetscScalar, dimension(:), allocatable :: RHarmonics_L, RHarmonics_H
    PetscScalar, dimension(:), allocatable :: ZHarmonics_L, ZHarmonics_H
    PetscScalar, dimension(:), allocatable :: dzHarmonics_L, dzHarmonics_H
    logical, dimension(:), allocatable :: BHarmonics_parity
    logical, dimension(:), allocatable :: BHarmonics_parityL, BHarmonics_parityH
    PetscScalar, dimension(:,:), allocatable :: hHat, duHatdtheta, duHatdzeta
    !!!!!!!
    ! MFM
    logical :: include_zero, sine_term
    integer :: nzeta_fft=200, ntheta_fft=200
    PetscScalar :: zeta_fft_max, dzeta_fft, dtheta_fft
    PetscScalar, dimension(200) :: theta_fft, zeta_fft
    PetscScalar, dimension(:,:), allocatable :: hHat_temp, hHat_star, BHat_temp, BHat_temp2
    PetscScalar, dimension(:), allocatable :: hHat_temp_k
    PetscScalar, dimension(:,:), allocatable :: hHat_tempL, hHat_starL, BHat_tempL, BHat_temp2L
    PetscScalar, dimension(:), allocatable :: hHat_temp_kL
    PetscScalar, dimension(:,:), allocatable :: hHat_tempH, hHat_starH, BHat_tempH, BHat_temp2H
    PetscScalar, dimension(:), allocatable :: hHat_temp_kH
    PetscScalar, dimension(:), allocatable :: BHarmonics_amplitudes_temp, BHarmonics_amplitudesL_temp, BHarmonics_amplitudesH_temp
    !!!!!!!
    PetscScalar :: R0
    PetscScalar, dimension(:,:), allocatable :: BHatL, dBHatdthetaL, dBHatdzetaL
    PetscScalar, dimension(:,:), allocatable :: BHatH, dBHatdthetaH, dBHatdzetaH
    PetscScalar, dimension(:,:), allocatable :: RHat,  dRHatdtheta,  dRHatdzeta,  d2RHatdtheta2,  d2RHatdzeta2,  d2RHatdthetadzeta
    PetscScalar, dimension(:,:), allocatable :: RHatL,dRHatdthetaL,dRHatdzetaL,d2RHatdtheta2L, d2RHatdzeta2L,d2RHatdthetadzetaL
    PetscScalar, dimension(:,:), allocatable :: RHatH,dRHatdthetaH,dRHatdzetaH,d2RHatdtheta2H,d2RHatdzeta2H,d2RHatdthetadzetaH
    PetscScalar, dimension(:,:), allocatable :: ZHat,  dZHatdtheta,  dZHatdzeta,  d2ZHatdtheta2,  d2ZHatdzeta2, d2ZHatdthetadzeta
    PetscScalar, dimension(:,:), allocatable :: ZHatL, dZHatdthetaL, dZHatdzetaL, d2ZHatdtheta2L, d2ZHatdzeta2L,d2ZHatdthetadzetaL
    PetscScalar, dimension(:,:), allocatable :: ZHatH, dZHatdthetaH, dZHatdzetaH, d2ZHatdtheta2H, d2ZHatdzeta2H,d2ZHatdthetadzetaH
    PetscScalar, dimension(:,:), allocatable :: Dz,  dDzdtheta,  dDzdzeta,  d2Dzdtheta2,  d2Dzdzeta2,  d2Dzdthetadzeta
    PetscScalar, dimension(:,:), allocatable :: DzL, dDzdthetaL, dDzdzetaL, d2Dzdtheta2L, d2Dzdzeta2L, d2DzdthetadzetaL
    PetscScalar, dimension(:,:), allocatable :: DzH, dDzdthetaH, dDzdzetaH, d2Dzdtheta2H, d2Dzdzeta2H, d2DzdthetadzetaH
    PetscScalar, dimension(:,:), allocatable :: geomang, dgeomangdtheta, dgeomangdzeta, d2geomangdtheta2, d2geomangdzeta2, d2geomangdthetadzeta
    PetscScalar, dimension(:,:), allocatable :: dXdtheta, dXdzeta, dYdtheta, dYdzeta
    PetscScalar, dimension(:,:), allocatable :: d2Xdtheta2, d2Xdthetadzeta, d2Xdzeta2, d2Ydtheta2, d2Ydthetadzeta, d2Ydzeta2
    PetscScalar, dimension(:,:), allocatable :: gradpsiX, gradpsiY, gradpsiZ, CX, CY, CZ
    
    integer :: fileUnit, didFileAccessWork
    character(len=200) :: lineOfFile
    integer, dimension(4) :: headerIntegers
    PetscScalar, dimension(3) :: headerReals
    PetscScalar, dimension(6) :: surfHeader
    PetscScalar, dimension(4) :: dataNumbers
    PetscScalar, dimension(8) :: data8Numbers
    integer, dimension(2) :: dataIntegers
    integer :: no_of_modes_old, no_of_modes_new, modeind, numB0s, startn, stopn
    PetscScalar :: iota_old, iota_new, GHat_old, GHat_new, IHat_old, IHat_new
    PetscScalar :: pPrimeHat_old, pPrimeHat_new, invFSA_BHat2
    logical :: end_of_file, proceed, include_mn, nearbyRadiiGiven, nonStelSym
    integer, parameter :: max_no_of_modes = 10000
    integer, dimension(max_no_of_modes) :: modesm_old, modesm_new, modesn_old, modesn_new
    PetscScalar, dimension(max_no_of_modes) :: modesb_old, modesb_new, modesR_old, modesR_new
    PetscScalar, dimension(max_no_of_modes) :: modesZ_old, modesZ_new, modesDz_old, modesDz_new
    PetscScalar :: rN_old,  rN_new, B0_old, B0_new, B0OverBBarL, B0OverBBarH
    PetscScalar :: R0_old, R0_new, R0L, R0H
    PetscScalar :: hHatHarmonics_amplitude, uHatHarmonics_amplitude
    PetscScalar :: dBHat_sub_psi_dthetaHarmonics_amplitude, dBHat_sub_psi_dzetaHarmonics_amplitude
    PetscScalar :: DeltapsiHat !, diotadpsiHat moved to global variables 2016-09-15 HS
    PetscScalar :: RadialWeight = 1.0 ! weight of closest surface with rN<=rN_wish

    ! For the BHarmonics_parity array, 
    ! true indicates the contribution to B(theta,zeta) has the form
    ! cos(l * theta - n * zeta)
    ! while false indicates the contribution to B(theta,zeta) has the form
    ! sin(l * theta - n * zeta)

    ! Initialise some quantities which will otherwise only be calculated
    ! in geometryScheme 11 and 12 
    dBHatdpsiHat = 0
    BHat_sub_psi = 0
    dBHat_sub_psi_dtheta = 0
    dBHat_sub_psi_dzeta = 0

    dBHat_sub_theta_dpsiHat = 0
    dBHat_sub_theta_dzeta = 0 !Always zero in Boozer coords.

    dBHat_sub_zeta_dpsiHat = 0
    dBHat_sub_zeta_dtheta = 0 !Always zero in Boozer coords.

    dBHat_sup_theta_dpsiHat = 0
    dBHat_sup_theta_dzeta = 0

    dBHat_sup_zeta_dpsiHat = 0
    dBHat_sup_zeta_dtheta = 0
    diotadpsiHat = 0

    nearbyRadiiGiven = .false. !Will be the case for all geometrySchemes except 11 and 12

    select case (geometryScheme)
    case (1)
       ! Three-helicity model:
       ! B = B0 * [1 + epsilon_t * cos(theta) + epsilon_h * cos(helicity_l * theta - helicity_n * zeta) ...
       !             + epsilon_antisymm * sin(helicity_antisymm_l * theta - helicity_antisymm_n * zeta) ...


       NHarmonics = 3
       allocate(BHarmonics_l(NHarmonics))
       allocate(BHarmonics_n(NHarmonics))
       allocate(BHarmonics_amplitudes(NHarmonics))
       allocate(BHarmonics_parity(NHarmonics))
       BHarmonics_parity = .true.

       i = 1
       BHarmonics_l(i) = 1
       BHarmonics_n(i) = 0
       BHarmonics_amplitudes(i) = epsilon_t

       i = 2
       BHarmonics_l(i) = helicity_l
       if (helicity_n == 0) then
          BHarmonics_n(i) = 0
       else
          BHarmonics_n(i) = 1
       end if
       BHarmonics_amplitudes(i) = epsilon_h

       i = 3
       BHarmonics_parity(i) = .false.
       BHarmonics_l(i) = helicity_antisymm_l
       if (helicity_n == 0) then
          BHarmonics_n(i) = helicity_antisymm_n
       else
          BHarmonics_n(i) = helicity_antisymm_n / helicity_n
       end if
       BHarmonics_amplitudes(i) = epsilon_antisymm

       if (helicity_n == 0) then
          if (helicity_antisymm_n .ne. 0) then
             print *,"WARNING: Typically, helicity_antisymm_n should be an integer multiple of helicity_n (possibly zero)."
          end if
       else
          if (mod(helicity_antisymm_n, helicity_n) .ne. 0) then
             print *,"WARNING: Typically, helicity_antisymm_n should be an integer multiple of helicity_n (possibly zero)."
          end if
       end if

       dGdpHat = 0 !Not implemented as an input for this case yet, could be put in namelist input if needed
       !rN = -1 !dummy
       rN = rN_wish
       BHarmonics_amplitudes=BHarmonics_amplitudes*B0OverBBar

    case (2)
       ! A three-harmonic approximation of the LHD standard configuration.
       ! Values taken from Table 1 of
       ! Beidler et al, Nuclear Fusion 51, 076001 (2011).

       iota = 0.4542d+0

       NHarmonics = 3
       allocate(BHarmonics_l(NHarmonics))
       allocate(BHarmonics_n(NHarmonics))
       allocate(BHarmonics_amplitudes(NHarmonics))
       allocate(BHarmonics_parity(NHarmonics))
       BHarmonics_parity = .true.

       i = 1
       BHarmonics_l(i) = 1
       BHarmonics_n(i) = 0
       BHarmonics_amplitudes(i) = -0.07053d+0

       i = 2
       BHarmonics_l(i) = 2
       BHarmonics_n(i) = 1
       BHarmonics_amplitudes(i) = 0.05067d+0

       i = 3
       BHarmonics_l(i) = 1
       BHarmonics_n(i) = 1
       BHarmonics_amplitudes(i) = -0.01476d+0
                    
       B0OverBBar = 1.0d+0  ! (Tesla)
       R0 = 3.7481d+0 ! (meters)
       GHat = B0OverBBar * R0
       IHat = 0
       dGdpHat = 0
       !rN = -1 !dummy
       rN = rN_wish
       BHarmonics_amplitudes=BHarmonics_amplitudes*B0OverBBar
                    
    case (3)
       ! A four-harmonic approximation of the LHD inward-shifted configuration.
       ! Values taken from Table 1 of
       ! Beidler et al, Nuclear Fusion 51, 076001 (2011).

       iota = 0.4692d+0

       NHarmonics = 4
       allocate(BHarmonics_l(NHarmonics))
       allocate(BHarmonics_n(NHarmonics))
       allocate(BHarmonics_amplitudes(NHarmonics))
       allocate(BHarmonics_parity(NHarmonics))
       BHarmonics_parity = .true.

       i = 1
       BHarmonics_l(i) = 1
       BHarmonics_n(i) = 0
       BHarmonics_amplitudes(i) = -0.05927d+0

       i = 2
       BHarmonics_l(i) = 2
       BHarmonics_n(i) = 1
       BHarmonics_amplitudes(i) = 0.05267d+0

       i = 3
       BHarmonics_l(i) = 1
       BHarmonics_n(i) = 1
       BHarmonics_amplitudes(i) = -0.04956d+0
                    
       i = 4
       BHarmonics_l(i) = 0
       BHarmonics_n(i) = 1
       BHarmonics_amplitudes(i) = 0.01045d+0
                    
       B0OverBBar = 1.0d+0  ! (Tesla)
       R0 = 3.6024d+0 ! (meters)
       GHat = B0OverBBar * R0
       IHat = 0
       dGdpHat = 0
       !normradius = -1 !dummy
       rN = rN_wish
       BHarmonics_amplitudes=BHarmonics_amplitudes*B0OverBBar

    case (4)
       ! A three-harmonic approximation of the W7-X standard configuration.
       ! Values taken from Table 1 of
       ! Beidler et al, Nuclear Fusion 51, 076001 (2011).

       iota = 0.8700d+0

       NHarmonics = 3
       allocate(BHarmonics_l(NHarmonics))
       allocate(BHarmonics_n(NHarmonics))
       allocate(BHarmonics_amplitudes(NHarmonics))
       allocate(BHarmonics_parity(NHarmonics))
       BHarmonics_parity = .true.

       i = 1
       BHarmonics_l(i) = 0
       BHarmonics_n(i) = 1
       BHarmonics_amplitudes(i) = 0.04645d+0

       i = 2
       BHarmonics_l(i) = 1
       BHarmonics_n(i) = 1
       BHarmonics_amplitudes(i) = -0.04351d+0

       i = 3
       BHarmonics_l(i) = 1
       BHarmonics_n(i) = 0
       BHarmonics_amplitudes(i) = -0.01902d+0
                    
       B0OverBBar = 3.089d+0  ! (Tesla)
       R0 = 5.5267d+0 ! (meters)
       GHat = -17.885d+0
       IHat = 0
       dGdpHat = 0
       !normradius = -1 !dummy
       rN = rN_wish
       BHarmonics_amplitudes=BHarmonics_amplitudes*B0OverBBar

    case (11)
       ! Read Boozer coordinate file in .bc format used at IPP Greifswald

       fileUnit = 11
       open(unit=fileUnit, file=equilibriumFile, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file ",equilibriumFile
          stop
       else
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             !print *,lineOfFile
             if (lineOfFile(1:2) /= "CC") exit
          end do

          ! Read header line:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) headerIntegers, headerReals
          !print *,'headerIntegers, headerReals='
          !print *,headerIntegers
          !print *,headerReals
          if (didFileAccessWork /= 0) then
             print *,"Unable to read header from the magnetic equilibrium file ",equilibriumFile
             stop
          end if

          !NPeriods = headerIntegers(4)
          !psiAHat  = headerReals(1)/2/pi !Convert the flux from Tm^2 to Tm^2/rad
          !aHat     = headerReals(2)      !minor radius in meters

          end_of_file = .false.

          rN_old = 0
          no_of_modes_old = 0
          modesm_old = 0
          modesn_old = 0
          modesb_old = 0
          modesR_old = 0
          modesZ_old = 0
          modesDz_old = 0
          iota_old = 0
          GHat_old = 0
          IHat_old = 0
          B0_old = 0
          R0_old = 0
          pPrimeHat_old = 0

          rN_new = 0
          no_of_modes_new = 0
          modesm_new = 0
          modesn_new = 0
          modesb_new = 0
          modesR_new = 0
          modesZ_new = 0
          modesDz_new = 0
          iota_new = 0
          GHat_new = 0
          IHat_new = 0
          B0_new = 0
          R0_new = 0
          pPrimeHat_new = 0

          ! Skip a line
          read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
          !print *,'s line='
          !print *,lineOfFile

          do 
             if ((rN_new .ge. rN_wish) .or. end_of_file) exit

             rN_old = rN_new
             no_of_modes_old = no_of_modes_new
             modesm_old = modesm_new
             modesn_old = modesn_new
             modesb_old = modesb_new
             modesR_old = modesR_new
             modesZ_old = modesZ_new
             modesDz_old = modesDz_new
             iota_old = iota_new
             GHat_old = GHat_new
             IHat_old = IHat_new
             B0_old = B0_new
             R0_old = R0_new
             pPrimeHat_old = pPrimeHat_new
             numB0s = 0

             if (.not.(index(lineOfFile,"[A]") > 0)) then !units were not on this line, so skip next
                ! Skip a line:
                read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
                !print *,'skip unit [A] line='
                !print *,lineOfFile
             end if

             ! Read the header for the magnetic surface:
             read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) surfHeader
             !print *,'surfHeader='
             !print *,surfHeader

             rN_new = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
             iota_new = surfHeader(2)
             ! Note that G and I have a minus sign in the following two lines
             ! because Ampere's law comes with a minus sign in the left-handed
             ! (r,pol,tor) system.
             GHat_new = -surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
             IHat_new = -surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
             pPrimeHat_new = surfheader(5)/psiAHat*(4*pi*1e-7)   ! dpdpsi=pPrimeHat/mu_0

             ! Skip units line:
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             proceed = .true.
             modeind = 0
             !print *,'s line='
             !print *,lineOfFile


             do
                if (.not. proceed) exit

                read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
                if (didFileAccessWork /= 0) then
                   proceed = .false.
                   end_of_file = .true.
                else if (index(lineOfFile,"s") > 0) then
                   ! Next flux surface has been reached
                   proceed = .false.
                   !print *,'s line='
                   !print *,lineOfFile
                else
                  read(unit=lineOfFile, fmt=*) dataIntegers, dataNumbers
                  !print *,'dataIntegers='
                  !print *,dataIntegers
                  !print *,'dataNumbers='
                  !print *,dataNumbers
                  if (dataIntegers(1) == 0 .and. dataIntegers(2) == 0) then
                      B0_new = dataNumbers(4)
                      R0_new = dataNumbers(1)
                      numB0s = numB0s + 1
                   else if (abs(dataNumbers(4)) > min_Bmn_to_load) then
                      modeind = modeind + 1
                      if (modeind > max_no_of_modes) then
                         print *,"The value of max_no_of_modes in geometry.F90 was insufficient."
                         print *,"Either increase this value and recompile, or else increase min_Bmn_to_load."
                         stop
                      end if
                      modesm_new(modeind) = dataIntegers(1)
                      modesn_new(modeind) = dataIntegers(2)
                      modesR_new(modeind) = dataNumbers(1)
                      modesZ_new(modeind) = dataNumbers(2)
                      modesDz_new(modeind) = dataNumbers(3)
                      modesb_new(modeind) = dataNumbers(4)
                   end if
                end if
             end do
             if (numB0s == 0) then
                print *,"Error: no (0,0) mode found in magnetic equilibrium file ",equilibriumFile
             else if (numB0s > 1) then
                print *,"Error: more than 1 (0,0) mode found in magnetic equilibrium file ",equilibriumFile
             end if
             no_of_modes_new = modeind
          end do

       end if

       close(unit = fileUnit)
       if (masterProc) then
          print *,"Successfully read magnetic equilibrium from file ",trim(equilibriumFile)
       end if

       DeltapsiHat = psiAHat * (rN_new*rN_new-rN_old*rN_old)
       nearbyRadiiGiven = .true.

       if (VMECRadialOption == 1) then !Choose the nearest flux surface available
          if (abs(rN_old - rN_wish) < abs(rN_new - rN_wish)) then
             RadialWeight = 1.0
             rN = rN_old
          else
             RadialWeight = 0.0
             rN = rN_new
          end if
       else !Linear interpolation in s=rN^2
          RadialWeight = (rN_new*rN_new-rN_wish*rN_wish) / (rN_new*rN_new-rN_old*rN_old)
          rN   = rN_wish
       end if
       iota = iota_old*RadialWeight+iota_new*(1.0-RadialWeight)
       GHat = GHat_old*RadialWeight+GHat_new*(1.0-RadialWeight)
       IHat = IHat_old*RadialWeight+IHat_new*(1.0-RadialWeight)
       B0OverBBar = B0_old*RadialWeight+B0_new*(1.0-RadialWeight)
       R0         = R0_old*RadialWeight+R0_new*(1.0-RadialWeight)
       pPrimeHat = pPrimeHat_old*RadialWeight+pPrimeHat_new*(1.0-RadialWeight)

       B0OverBBarL=B0_old
!       print *,"B0OverBBar = ",B0OverBBar
!       print *,"B0OverBBarL = ",B0OverBBarL 
!       print *,"rN_old = ",rN_old
       B0OverBBarH=B0_new
!       print *,"B0OverBBarH = ",B0OverBBarH 
!       print *,"rN_new = ",rN_new
       R0L=R0_old
       R0H=R0_new
       NHarmonicsL = no_of_modes_old
       NHarmonicsH = no_of_modes_new
       allocate(BHarmonics_lL(NHarmonicsL))
       allocate(BHarmonics_nL(NHarmonicsL))
       allocate(BHarmonics_amplitudesL(NHarmonicsL))
       allocate(RHarmonics_L(NHarmonicsL))
       allocate(ZHarmonics_L(NHarmonicsL))
       allocate(DzHarmonics_L(NHarmonicsL))
       allocate(BHarmonics_parityL(NHarmonicsL))
       allocate(BHarmonics_lH(NHarmonicsH))
       allocate(BHarmonics_nH(NHarmonicsH))
       allocate(BHarmonics_amplitudesH(NHarmonicsH))
       allocate(RHarmonics_H(NHarmonicsH))
       allocate(ZHarmonics_H(NHarmonicsH))
       allocate(DzHarmonics_H(NHarmonicsH))
       allocate(BHarmonics_parityH(NHarmonicsH))
       BHarmonics_lL = modesm_old(1:NHarmonicsL)
       BHarmonics_nL = modesn_old(1:NHarmonicsL)
       BHarmonics_amplitudesL = modesb_old(1:NHarmonicsL)
       RHarmonics_L = modesR_old(1:NHarmonicsL)
       ZHarmonics_L = modesZ_old(1:NHarmonicsL)
       DzHarmonics_L = modesDz_old(1:NHarmonicsL)*2*pi/NPeriods
       BHarmonics_lH = modesm_new(1:NHarmonicsH)
       BHarmonics_nH = modesn_new(1:NHarmonicsH)
       BHarmonics_amplitudesH = modesb_new(1:NHarmonicsH)
       RHarmonics_H  = modesR_new(1:NHarmonicsH)
       ZHarmonics_H  = modesZ_new(1:NHarmonicsH)
       DzHarmonics_H = modesDz_new(1:NHarmonicsH)*2*pi/NPeriods
       BHarmonics_parityL = .true.
       BHarmonics_parityH = .true.

       !dGdpHat=(GHat_new-GHat_old)/(rN_new*rN_new-rN_old*rN_old)/pPrimeHat

       !Unnnecessary step:
       !BHarmonics_amplitudes = BHarmonics_amplitudes / B0OverBBar

       if (GHat*psiAHat>0) then
          !Note that GHat and psiAHat already have the opposite sign to the corresponding quantities in the .bc file
          !Therefore, the flip is performed if they have the same sign here.
          !print *,"This is a stellarator symmetric file from Joachim Geiger. It will now be turned 180 degrees around a horizontal axis <=> flip the sign of G and I, so that it matches the sign of its total toroidal flux."
          GHat    =-GHat
          GHat_new=-GHat_new
          GHat_old=-GHat_old
          IHat    =-IHat
          IHat_new=-IHat_new
          IHat_old=-IHat_old
          !dGdpHat=-dGdpHat
       end if
       
       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       !(The toroidal direction sign switch psiAHat=psiAHat*(-1) was already made in the initializeGeometry routine)
       GHat     = GHat*(-1)                   !toroidal direction sign switch
       GHat_new = GHat_new*(-1)               !toroidal direction sign switch
       GHat_old = GHat_old*(-1)               !toroidal direction sign switch
       iota     = iota*(-1)                   !toroidal direction sign switch
       iota_new = iota_new*(-1)               !toroidal direction sign switch
       iota_old = iota_old*(-1)               !toroidal direction sign switch
       if (.not. nearbyRadiiGiven) then
          BHarmonics_n=BHarmonics_n*(-1) !toroidal direction sign switch
       else
          BHarmonics_nL=BHarmonics_nL*(-1) !toroidal direction sign switch
          BHarmonics_nH=BHarmonics_nH*(-1) !toroidal direction sign switch
       end if
       DzHarmonics_L = DzHarmonics_L*(-1) !toroidal direction sign switch 
       DzHarmonics_H = DzHarmonics_H*(-1) !toroidal direction sign switch 
       
       dBHat_sub_zeta_dpsiHat = (GHat_new-GHat_old)/DeltapsiHat
       dBHat_sub_theta_dpsiHat =(IHat_new-IHat_old)/DeltapsiHat
       diotadpsiHat= (iota_new-iota_old)/DeltapsiHat

    case (12)
       ! Read Boozer coordinate file in a generalisation of the .bc format used at IPP Greifswald for non-stellarator symmetric equilibria 

       fileUnit = 11
       open(unit=fileUnit, file=equilibriumFile, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file ",equilibriumFile
          stop
       else
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             if (lineOfFile(1:2) /= "CC") exit
          end do

          ! Read header line:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) headerIntegers, headerReals
          if (didFileAccessWork /= 0) then
             print *,"Unable to read header from the magnetic equilibrium file ",equilibriumFile
             stop
          end if

          !NPeriods = headerIntegers(4)
          !psiAHat  = headerReals(1)/2/pi !Convert the flux from Tm^2 to Tm^2/rad
          !aHat     = headerReals(2)      !minor radius in meters

          end_of_file = .false.

          rN_old = 0
          no_of_modes_old = 0
          modesm_old = 0
          modesn_old = 0
          modesb_old = 0
          modesR_old = 0
          modesZ_old = 0
          modesDz_old = 0
          iota_old = 0
          GHat_old = 0
          IHat_old = 0
          B0_old = 0
          R0_old = 0
          pPrimeHat_old = 0

          rN_new = 0
          no_of_modes_new = 0
          modesm_new = 0
          modesn_new = 0
          modesb_new = 0
          modesR_new = 0
          modesZ_new = 0
          modesDz_new = 0
          iota_new = 0
          GHat_new = 0
          IHat_new = 0
          B0_new = 0
          R0_new = 0
          pPrimeHat_new = 0

          ! Skip a line
          read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile

          do 
             if ((rN_new .ge. rN_wish) .or. end_of_file) exit

             rN_old = rN_new
             no_of_modes_old = no_of_modes_new
             modesm_old = modesm_new
             modesn_old = modesn_new
             modesb_old = modesb_new
             modesR_old = modesR_new
             modesZ_old = modesZ_new
             modesDz_old = modesDz_new
             iota_old = iota_new
             GHat_old = GHat_new
             IHat_old = IHat_new
             B0_old = B0_new
             R0_old = R0_new
             pPrimeHat_old = pPrimeHat_new
             numB0s = 0

             ! Skip a line:
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             ! Read the header for the magnetic surface:
             read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) surfHeader

             rN_new = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
             iota_new = surfHeader(2)
             ! Note that G and I has a minus sign in the following two lines
             ! because Ampere's law comes with a minus sign in the left-handed
             ! (r,pol,tor) system.
             GHat_new = -surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
             IHat_new = -surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
             pPrimeHat_new = surfheader(5)/psiAHat*(4*pi*1e-7)   ! dpdpsi=pPrimeHat/mu_0

             ! Skip units line:
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             proceed = .true.
             modeind = 0
             do
                if (.not. proceed) exit

                read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
                if (didFileAccessWork /= 0) then
                   proceed = .false.
                   end_of_file = .true.
                else if (index(lineOfFile,"s") > 0) then
                   ! Next flux surface has been reached
                   proceed = .false.
                else
                   read(unit=lineOfFile, fmt=*) dataIntegers, data8Numbers
                   if (dataIntegers(1) == 0 .and. dataIntegers(2) == 0) then
                      B0_new = data8Numbers(7)
                      R0_new = data8Numbers(1)
                      numB0s = numB0s + 1
                   else if (abs(data8Numbers(7)) > min_Bmn_to_load) then
                      if (modeind + 2 > max_no_of_modes) then
                         print *,"The value of max_no_of_modes in geometry.F90 was insufficient."
                         print *,"Either increase this value and recompile, or else increase min_Bmn_to_load."
                         stop
                      end if
                      modeind = modeind + 1
                      modesm_new(modeind) = dataIntegers(1)
                      modesn_new(modeind) = dataIntegers(2)
                      modesR_new(modeind) = data8Numbers(1) !Cosinus component
                      modesZ_new(modeind) = data8Numbers(4) !Sinus component
                      modesDz_new(modeind)= data8Numbers(6) !Sinus component
                      modesb_new(modeind) = data8Numbers(7) !Cosinus component
                      modeind = modeind + 1
                      modesm_new(modeind) = dataIntegers(1)
                      modesn_new(modeind) = dataIntegers(2)
                      modesR_new(modeind) = data8Numbers(2) !Sinus component
                      modesZ_new(modeind) = data8Numbers(3) !Cosinus component
                      modesDz_new(modeind)= data8Numbers(5) !Cosinus component
                      modesb_new(modeind) = data8Numbers(8) !Sinus component
                   end if
                end if
             end do
             if (numB0s == 0) then
                print *,"Error: no (0,0) mode found in magnetic equilibrium file ",equilibriumFile
             else if (numB0s > 1) then
                print *,"Error: more than 1 (0,0) mode found in magnetic equilibrium file ",equilibriumFile
             end if
             no_of_modes_new = modeind
          end do

       end if

       close(unit = fileUnit)
       if (masterProc) then
          print *,"Successfully read magnetic equilibrium from file ",trim(equilibriumFile)
       end if

       DeltapsiHat = psiAHat * (rN_new*rN_new-rN_old*rN_old)
       nearbyRadiiGiven = .true.

       if (VMECRadialOption == 1) then !Choose the nearest flux surface available
          if (abs(rN_old - rN_wish) < abs(rN_new - rN_wish)) then
             RadialWeight = 1.0
             rN = rN_old
          else
             RadialWeight = 0.0
             rN = rN_new
          end if
       else !Linear interpolation in s=rN^2
          RadialWeight = (rN_new*rN_new-rN_wish*rN_wish) / (rN_new*rN_new-rN_old*rN_old)
          rN   = rN_wish
       end if
       iota = iota_old*RadialWeight+iota_new*(1.0-RadialWeight)
       GHat = GHat_old*RadialWeight+GHat_new*(1.0-RadialWeight)
       IHat = IHat_old*RadialWeight+IHat_new*(1.0-RadialWeight)
       B0OverBBar = B0_old*RadialWeight+B0_new*(1.0-RadialWeight)
       R0         = R0_old*RadialWeight+R0_new*(1.0-RadialWeight)
       pPrimeHat = pPrimeHat_old*RadialWeight+pPrimeHat_new*(1.0-RadialWeight)
       
       B0OverBBarL=B0_old
       B0OverBBarH=B0_new
       R0L=R0_old
       R0H=R0_new
       NHarmonicsL = no_of_modes_old
       NHarmonicsH = no_of_modes_new
       allocate(BHarmonics_lL(NHarmonicsL))
       allocate(BHarmonics_nL(NHarmonicsL))
       allocate(BHarmonics_amplitudesL(NHarmonicsL))
       allocate(RHarmonics_L(NHarmonicsL))
       allocate(ZHarmonics_L(NHarmonicsL))
       allocate(DzHarmonics_L(NHarmonicsL))
       allocate(BHarmonics_parityL(NHarmonicsL))
       allocate(BHarmonics_lH(NHarmonicsH))
       allocate(BHarmonics_nH(NHarmonicsH))
       allocate(BHarmonics_amplitudesH(NHarmonicsH))
       allocate(RHarmonics_H(NHarmonicsH))
       allocate(ZHarmonics_H(NHarmonicsH))
       allocate(DzHarmonics_H(NHarmonicsH))
       allocate(BHarmonics_parityH(NHarmonicsH))
       BHarmonics_lL = modesm_old(1:NHarmonicsL)
       BHarmonics_nL = modesn_old(1:NHarmonicsL)
       BHarmonics_amplitudesL = modesb_old(1:NHarmonicsL)
       RHarmonics_L = modesR_old(1:NHarmonicsL)
       ZHarmonics_L = modesZ_old(1:NHarmonicsL)
       DzHarmonics_L = modesDz_old(1:NHarmonicsL)*2*pi/NPeriods
       BHarmonics_lH = modesm_new(1:NHarmonicsH)
       BHarmonics_nH = modesn_new(1:NHarmonicsH)
       BHarmonics_amplitudesH = modesb_new(1:NHarmonicsH)
       RHarmonics_H = modesR_new(1:NHarmonicsH)
       ZHarmonics_H = modesZ_new(1:NHarmonicsH)
       DzHarmonics_H = modesDz_new(1:NHarmonicsH)*2*pi/NPeriods
       do i = 0, NHarmonicsL/2-1
          BHarmonics_parityL(2*i+1)=.true.
          BHarmonics_parityL(2*i+2)=.false.
       end do
       do i = 0, NHarmonicsH/2-1
          BHarmonics_parityH(2*i+1)=.true.
          BHarmonics_parityH(2*i+2)=.false.
       end do

       !This unnecessary scaling has has been removed !!!!!!!
       !BHarmonics_amplitudes = BHarmonics_amplitudes / B0OverBBar

       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       !(The toroidal direction sign switch psiAHat=psiAHat*(-1) was already made in the initializeGeometry routine)
       GHat     = GHat*(-1)                   !toroidal direction sign switch
       GHat_new = GHat_new*(-1)               !toroidal direction sign switch
       GHat_old = GHat_old*(-1)               !toroidal direction sign switch
       iota     = iota*(-1)                   !toroidal direction sign switch
       iota_new = iota_new*(-1)               !toroidal direction sign switch
       iota_old = iota_old*(-1)               !toroidal direction sign switch
       if (.not. nearbyRadiiGiven) then
          BHarmonics_n=BHarmonics_n*(-1) !toroidal direction sign switch
       else
          BHarmonics_nL=BHarmonics_nL*(-1) !toroidal direction sign switch
          BHarmonics_nH=BHarmonics_nH*(-1) !toroidal direction sign switch
       end if
       DzHarmonics_L = DzHarmonics_L*(-1) !toroidal direction sign switch 
       DzHarmonics_H = DzHarmonics_H*(-1) !toroidal direction sign switch 

       dBHat_sub_zeta_dpsiHat = (GHat_new-GHat_old)/DeltapsiHat
       dBHat_sub_theta_dpsiHat =(IHat_new-IHat_old)/DeltapsiHat
       diotadpsiHat= (iota_new-iota_old)/DeltapsiHat

    case default
       print *,"Error! Invalid geometryScheme"
       stop
    end select


    if (.not. nearbyRadiiGiven) then
       ! Initialize arrays:
       BHat = B0OverBBar ! This includes the (0,0) component.
       dBHatdtheta = 0
       dBHatdzeta = 0

       !I do not Bother to calculate Sugama's drift for geometryScheme=1,2,3,4.
!!$       RHat = 0 
!!$       dRHatdtheta      = 0
!!$       dRHatdzeta       = 0
!!$       d2RHatdtheta2    = 0
!!$       d2RHatdzeta2     = 0
!!$       d2RHatdthetadzeta= 0
!!$       ZHatL = 0 
!!$       dZHatdtheta      = 0
!!$       dZHatdzeta       = 0
!!$       d2ZHatdtheta2    = 0
!!$       d2ZHatdzeta2     = 0
!!$       d2ZHatdthetadzeta= 0
!!$       DzL = 0 
!!$       dDzdtheta      = 0
!!$       dDzdzeta       = 0
!!$       d2Dzdtheta2    = 0
!!$       d2Dzdzeta2     = 0
!!$       d2Dzdthetadzeta= 0

       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MFM scaling of non-quasisymmetric terms
       if (symm_breaking == 0) then
          print *,"There will be no scaling of symmetry-breaking terms"

       else if (symm_breaking == 1) then

          if (symm_type == 0) then
             print *,"The non-quasi-AXI-symmetric Boozer harmonics of |B| will be scaled"
             do i = 1, NHarmonics
                if (BHarmonics_n(i) == 0) then
                   BHarmonics_amplitudes(i) = BHarmonics_amplitudes(i)
                else
                   BHarmonics_amplitudes(i) = BHarmonics_amplitudes(i)*epsilon_symmbreak
                end if
             end do

          else if (symm_type == 1) then
             print *,"The non-quasi-HELICALLY-symmetric Boozer harmonics of |B| will be scaled"
             do i = 1, NHarmonics
                if (qhs_poloidal*(-1)*BHarmonics_n(i) == qhs_toroidal*BHarmonics_l(i)) then
                   BHarmonics_amplitudes(i) = BHarmonics_amplitudes(i)
                else
                   BHarmonics_amplitudes(i) = BHarmonics_amplitudes(i)*epsilon_symmbreak
                end if
             end do

          else
             print *,"The type of quasisymmetry-breaking must take integer values of either 0 (QAS) or 1 (QHS) ... Exiting"
             stop
          end if
       else if (symm_breaking == 2) then
          zeta_fft_max = 2*pi/NPeriods
          dzeta_fft = zeta_fft_max/nzeta_fft
          dtheta_fft = 2.*pi/ntheta_fft
          allocate(hHat_temp(ntheta_fft,nzeta_fft))
          allocate(hHat_temp_k(NHarmonics*2+1))
          allocate(hHat_star(ntheta_fft,nzeta_fft))
          allocate(BHat_temp(ntheta_fft,nzeta_fft))
          allocate(BHat_temp2(ntheta_fft,nzeta_fft))
          allocate(BHarmonics_amplitudes_temp(NHarmonics*2+1))
          hHat_temp = 0
          hHat_temp_k = 0
          hHat_star = 0 !1.0/(B0OverBBar*B0OverBBar)
          BHat_temp = 0
          BHat_temp2 = 0
          include_zero = 0
          sine_term = 0
          BHat = 0

          zeta_fft(0) = 0
          theta_fft(0) = 0
          do i = 1, nzeta_fft
             zeta_fft(i) = dzeta_fft + zeta_fft(i-1)
          end do
          do i = 1, ntheta_fft
             theta_fft(i) = dtheta_fft + theta_fft(i-1)
          end do

          call inverse_b(BHarmonics_amplitudes, BHat_temp, BHarmonics_n, BHarmonics_l, NHarmonics, include_zero, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)
          BHat_temp = BHat_temp + B0OverBBar
          include_zero = 1
          sine_term = 1

          do itheta = 1, ntheta_fft
             do izeta = 1, nzeta_fft
                hHat_temp(itheta,izeta) = 1.0/(BHat_temp(itheta,izeta)*BHat_temp(itheta,izeta))
             end do
          end do
          call forward_b(hHat_temp, hHat_temp_k, BHarmonics_n, BHarmonics_l, NHarmonics, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)

          if (symm_type == 1) then ! QHS      
             do i = 2, NHarmonics+1 ! The (0,0) mode is not be scaled
                if (qhs_poloidal*(-1)*BHarmonics_n(i-1) == qhs_toroidal*BHarmonics_l(i-1)) then
                   hHat_temp_k(i) = hHat_temp_k(i)
                else
                   hHat_temp_k(i) = epsilon_symmbreak*hHat_temp_k(i)
                end if
             end do
          else
             do i = 2, NHarmonics+1 ! The (0,0) mode is not be scaled
                if (BHarmonics_n(i-1)==0) then ! if mode is QAS
                   hHat_temp_k(i) = hHat_temp_k(i)
                else
                   hHat_temp_k(i) = epsilon_symmbreak*hHat_temp_k(i)
                end if
             end do
          end if
          
          call inverse_b(hHat_temp_k, hHat_star, BHarmonics_n, BHarmonics_l, NHarmonics, include_zero, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)
          do itheta = 1, ntheta_fft
             do izeta = 1, nzeta_fft
                BHat_temp2(itheta,izeta) = sqrt(1.0/hHat_star(itheta,izeta))
             end do
          end do
          BHarmonics_amplitudes_temp = 0
          BHarmonics_amplitudes = 0
          sine_term = 0
          call forward_b(BHat_temp2, BHarmonics_amplitudes_temp, BHarmonics_n, BHarmonics_l, NHarmonics, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)

          BHat = BHarmonics_amplitudes_temp(1) ! (0,0) mode amplitude
          do i = 2, NHarmonics+1
             BHarmonics_amplitudes(i-1) = BHarmonics_amplitudes_temp(i)
          end do

          deallocate(hHat_temp)
          deallocate(hHat_temp_k)
          deallocate(hHat_star)
          deallocate(BHat_temp)
          deallocate(BHat_temp2)
          deallocate(BHarmonics_amplitudes_temp)
       else
          print *,"You have selected an invalid integer for the method of scaling symmetry-breaking terms. Input must be 0 (no scaling), 1 (scaling through |B|), or 2 (scaling through invB^2) ... Exiting"
          stop
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do i = 1, NHarmonics
          if (BHarmonics_parity(i)) then   ! The cosine components of BHat
             include_mn = .false.
             if ((abs(BHarmonics_n(i))<=int(Nzeta/2.0)).and.(BHarmonics_l(i)<=int(Nzeta/2.0))) then
                include_mn = .true.
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   BHat(itheta,:) = BHat(itheta,:) + BHarmonics_amplitudes(i) * &
                   cos(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)
                   
                   dBHatdtheta(itheta,:) = dBHatdtheta(itheta,:) - BHarmonics_amplitudes(i) * BHarmonics_l(i) * &
                   sin(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)
                      
                   dBHatdzeta(itheta,:) = dBHatdzeta(itheta,:) + BHarmonics_amplitudes(i) * Nperiods * BHarmonics_n(i) * &
                   sin(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)
                end do
             end if
          else  ! The sine components of BHat
             include_mn=.false.
             if ((abs(BHarmonics_n(i))<=int(Nzeta/2.0)).and.(BHarmonics_l(i)<=int(Nzeta/2.0))) then
                include_mn=.true.
             end if
             if (BHarmonics_l(i)==0 .or. real(BHarmonics_l(i))==Ntheta/2.0) then
                if (BHarmonics_n(i)==0 .or. abs(real(BHarmonics_n(i)))==Nzeta/2.0 ) then
                   include_mn=.false.
                end if
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   BHat(itheta,:) = BHat(itheta,:) + BHarmonics_amplitudes(i) * &
                   sin(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)
                      
                   dBHatdtheta(itheta,:) = dBHatdtheta(itheta,:) + BHarmonics_amplitudes(i) * BHarmonics_l(i) * &
                   cos(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)
                      
                   dBHatdzeta(itheta,:) = dBHatdzeta(itheta,:) - BHarmonics_amplitudes(i) * Nperiods * BHarmonics_n(i) * &
                   cos(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)
                   
                end do
             end if
          end if
       end do
    else !Two nearby radii L and H are given 
       allocate(BHatL(Ntheta,Nzeta))
       allocate(dBHatdthetaL(Ntheta,Nzeta))
       allocate(dBHatdzetaL(Ntheta,Nzeta))

       allocate(RHatL(Ntheta,Nzeta))
       allocate(dRHatdthetaL(Ntheta,Nzeta))
       allocate(dRHatdzetaL(Ntheta,Nzeta))
       allocate(d2RHatdtheta2L(Ntheta,Nzeta))
       allocate(d2RHatdzeta2L(Ntheta,Nzeta))
       allocate(d2RHatdthetadzetaL(Ntheta,Nzeta))

       allocate(ZHatL(Ntheta,Nzeta))
       allocate(dZHatdthetaL(Ntheta,Nzeta))
       allocate(dZHatdzetaL(Ntheta,Nzeta))
       allocate(d2ZHatdtheta2L(Ntheta,Nzeta))
       allocate(d2ZHatdzeta2L(Ntheta,Nzeta))
       allocate(d2ZHatdthetadzetaL(Ntheta,Nzeta))

       allocate(DzL(Ntheta,Nzeta))
       allocate(dDzdthetaL(Ntheta,Nzeta))
       allocate(dDzdzetaL(Ntheta,Nzeta))
       allocate(d2Dzdtheta2L(Ntheta,Nzeta))
       allocate(d2Dzdzeta2L(Ntheta,Nzeta))
       allocate(d2DzdthetadzetaL(Ntheta,Nzeta))

       BHatL = B0OverBBarL ! This includes the (0,0) component.
       dBHatdthetaL = 0
       dBHatdzetaL = 0

       RHatL = R0L ! This includes the (0,0) component.
       dRHatdthetaL      = 0
       dRHatdzetaL       = 0
       d2RHatdtheta2L    = 0
       d2RHatdzeta2L     = 0
       d2RHatdthetadzetaL= 0
       ZHatL = 0 ! This includes the (0,0) component.
       dZHatdthetaL      = 0
       dZHatdzetaL       = 0
       d2ZHatdtheta2L    = 0
       d2ZHatdzeta2L     = 0
       d2ZHatdthetadzetaL= 0
       DzL = 0 ! This includes the (0,0) component.
       dDzdthetaL      = 0
       dDzdzetaL       = 0
       d2Dzdtheta2L    = 0
       d2Dzdzeta2L     = 0
       d2DzdthetadzetaL= 0
       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MFM scaling of non-quasisymmetric terms
       if (symm_breaking == 0) then
          print *,"There will be no scaling of symmetry-breaking terms"

       else if (symm_breaking == 1) then

          if (symm_type == 0) then
             print *,"The non-quasi-AXI-symmetric Boozer harmonics of |B| will be scaled"
             do i = 1, NHarmonicsL
                if (BHarmonics_nL(i) == 0) then
                   BHarmonics_amplitudesL(i) = BHarmonics_amplitudesL(i)
                else
                   BHarmonics_amplitudesL(i) = BHarmonics_amplitudesL(i)*epsilon_symmbreak
                end if
             end do

          else if (symm_type == 1) then
             print *,"The non-quasi-HELICALLY-symmetric Boozer harmonics of |B| will be scaled"
             do i = 1, NHarmonicsL
                if (qhs_poloidal*(-1)*BHarmonics_nL(i) == qhs_toroidal*BHarmonics_lL(i)) then
                   BHarmonics_amplitudesL(i) = BHarmonics_amplitudesL(i)
                else
                   BHarmonics_amplitudesL(i) = BHarmonics_amplitudesL(i)*epsilon_symmbreak
                end if
             end do

          else
             print *,"The type of quasisymmetry-breaking must take integer values of either 0 (QAS) or 1 (QHS) ... Exiting"
             stop
          end if
       else if (symm_breaking == 2) then
          zeta_fft_max = 2*pi/NPeriods
          dzeta_fft = zeta_fft_max/nzeta_fft
          dtheta_fft = 2.*pi/ntheta_fft
          allocate(hHat_tempL(ntheta_fft,nzeta_fft))
          allocate(hHat_temp_kL(NHarmonicsL*2+1))
          allocate(hHat_starL(ntheta_fft,nzeta_fft))
          allocate(BHat_tempL(ntheta_fft,nzeta_fft))
          allocate(BHat_temp2L(ntheta_fft,nzeta_fft))
          allocate(BHarmonics_amplitudesL_temp(NHarmonicsL*2+1))
          hHat_tempL = 0
          hHat_temp_kL = 0
          hHat_starL = 0
          BHat_tempL = 0
          BHat_temp2L = 0
          include_zero = 0
          sine_term = 0
          BHatL = 0

          zeta_fft(0) = 0
          theta_fft(0) = 0
          do i = 1, nzeta_fft
             zeta_fft(i) = dzeta_fft + zeta_fft(i-1)
          end do
          do i = 1, ntheta_fft
             theta_fft(i) = dtheta_fft + theta_fft(i-1)
          end do

          call inverse_b(BHarmonics_amplitudesL, BHat_tempL, BHarmonics_nL, BHarmonics_lL, NHarmonicsL, include_zero, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)
          BHat_tempL = BHat_tempL + B0OverBBarL
          include_zero = 1
          sine_term = 1

          do itheta = 1, ntheta_fft
             do izeta = 1, nzeta_fft
                hHat_tempL(itheta,izeta) = 1.0/(BHat_tempL(itheta,izeta)*BHat_tempL(itheta,izeta))
!                if (zeta(izeta) == 0) then
!                   print *,"BHat_tempL = ",BHat_tempL(itheta,izeta)
!                   print *,"hHat_temp = ",hHat_tempL(itheta,izeta)
!                end if

             end do
          end do
          call forward_b(hHat_tempL, hHat_temp_kL, BHarmonics_nL, BHarmonics_lL, NHarmonicsL, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)
!          do i = 2, NHarmonicsL+1
!          end do

          if (symm_type == 1) then ! QHS      
             do i = 2, NHarmonicsL+1 ! The (0,0) mode is not scaled
                if (qhs_poloidal*(-1)*BHarmonics_nL(i-1) == qhs_toroidal*BHarmonics_lL(i-1)) then
!                   print *,"BHarmonics_nL = ",BHarmonics_nL(i-1),"BHarmonics_lL = ",BHarmonics_lL(i-1)
                   hHat_temp_kL(i) = hHat_temp_kL(i)
                else
                   hHat_temp_kL(i) = epsilon_symmbreak*hHat_temp_kL(i)
                end if
             end do
          else
             do i = 2, NHarmonicsL+1 ! The (0,0) mode is not scaled
!                print *,"hHat_temp_kL = ",hHat_temp_kL(i)
                if (BHarmonics_nL(i-1)==0) then ! if mode is QAS
                   hHat_temp_kL(i) = hHat_temp_kL(i)
!                   print *,"n = ",BHarmonics_nL(i-1)
!                   print *,"m = ",BHarmonics_lL(i-1)
!                   print *,"hHat_temp_kL(i) ", hHat_temp_kL(i)
                else
                   hHat_temp_kL(i) = epsilon_symmbreak*hHat_temp_kL(i)
!                   print *,"n = ",BHarmonics_nL(i-1)
!                   print *,"m = ",BHarmonics_lL(i-1)
!                   print *,"hHat_temp_kL(i) ", hHat_temp_kL(i)
                end if
             end do
          end if

          call inverse_b(hHat_temp_kL, hHat_starL, BHarmonics_nL, BHarmonics_lL, NHarmonicsL, include_zero, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)

          do itheta = 1, ntheta_fft
!             print *,"theta = ",theta(itheta)
!             print *,"hHat_star = ",hHat_starL(itheta,1)
             do izeta = 1, nzeta_fft
!                print *,"zeta = ",zeta(izeta)
                BHat_temp2L(itheta,izeta) = sqrt(1.0/hHat_starL(itheta,izeta))
!                if (izeta == 5) then
                !   print *,"hHat_starL(itheta,1) =",hHat_starL(itheta,izeta)
                !   print *,"BHat_temp2L(itheta,1) = ",BHat_temp2L(itheta,izeta)
!                end if
!                print *,"hHat_star = ",hHat_starL(itheta,izeta)
!                print *,"BHat_temp2 = ",BHat_temp2L(itheta,izeta)
             end do
          end do
          BHarmonics_amplitudesL_temp = 0
          BHarmonics_amplitudesL = 0
          sine_term = 0

          call forward_b(BHat_temp2L, BHarmonics_amplitudesL_temp, BHarmonics_nL, BHarmonics_lL, NHarmonicsL, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)

          BHatL = BHarmonics_amplitudesL_temp(1) ! (0,0) mode amplitude                  
          do i = 2, NHarmonicsL+1
             if (abs(BHarmonics_amplitudesL_temp(i)) < 1.e-10) then
                BHarmonics_amplitudesL(i-1) = 0
             else
                BHarmonics_amplitudesL(i-1) = BHarmonics_amplitudesL_temp(i)
             end if
!             print *,"n = ",BHarmonics_nL(i-1)
!             print *,"m = ",BHarmonics_lL(i-1)
!             print *,"BHarmonics_amplitudesL = ",BHarmonics_amplitudesL_temp(i)
         end do

          deallocate(hHat_tempL)
          deallocate(hHat_temp_kL)
          deallocate(hHat_starL)
          deallocate(BHat_tempL)
          deallocate(BHat_temp2L)
          deallocate(BHarmonics_amplitudesL_temp)
       else
          print *,"You have selected an invalid integer for the method of scaling symmetry-breaking terms. Input must be 0 (no scaling), 1 (scaling through |B|), or 2 (scaling through invB^2) ... Exiting"
          stop
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       do i = 1, NHarmonicsL
!          print *,"BHarmonicsL(i) = ",BHarmonics_amplitudesL(i)
          if (BHarmonics_parityL(i)) then   ! The cosine components of BHat
             include_mn = .false.
             if ((abs(BHarmonics_nL(i))<=int(Nzeta/2.0)).and.(BHarmonics_lL(i)<=int(Ntheta/2.0))) then
                include_mn = .true.
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   BHatL(itheta,:) = BHatL(itheta,:) + BHarmonics_amplitudesL(i) * &
                   cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)
                      
                   dBHatdthetaL(itheta,:) = dBHatdthetaL(itheta,:) - BHarmonics_amplitudesL(i) * BHarmonics_lL(i) * &
                   sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)
                      
                   dBHatdzetaL(itheta,:) = dBHatdzetaL(itheta,:) + BHarmonics_amplitudesL(i) * Nperiods * BHarmonics_nL(i) * &
                   sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)
                                       
                   !The following are only there to calculate Sugama's magnetic drift
                   RHatL(itheta,:) = RHatL(itheta,:) + RHarmonics_L(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dRHatdthetaL(itheta,:) = dRHatdthetaL(itheta,:) - RHarmonics_L(i) * BHarmonics_lL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dRHatdzetaL(itheta,:) = dRHatdzetaL(itheta,:) + RHarmonics_L(i) * Nperiods * BHarmonics_nL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2RHatdtheta2L(itheta,:) = d2RHatdtheta2L(itheta,:) - RHarmonics_L(i) * BHarmonics_lL(i)**2 * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2RHatdzeta2L(itheta,:) = d2RHatdzeta2L(itheta,:) - RHarmonics_L(i) * (Nperiods * BHarmonics_nL(i))**2 * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2RHatdthetadzetaL(itheta,:) = d2RHatdthetadzetaL(itheta,:) + RHarmonics_L(i) * BHarmonics_lL(i)*Nperiods*BHarmonics_nL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   ZHatL(itheta,:) = ZHatL(itheta,:) + ZHarmonics_L(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dZHatdthetaL(itheta,:) = dZHatdthetaL(itheta,:) + ZHarmonics_L(i) * BHarmonics_lL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dZHatdzetaL(itheta,:) = dZHatdzetaL(itheta,:) - ZHarmonics_L(i) * Nperiods * BHarmonics_nL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2ZHatdtheta2L(itheta,:) = d2ZHatdtheta2L(itheta,:) - ZHarmonics_L(i) * BHarmonics_lL(i)**2 * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2ZHatdzeta2L(itheta,:) = d2ZHatdzeta2L(itheta,:) - ZHarmonics_L(i) * (Nperiods * BHarmonics_nL(i))**2 * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2ZHatdthetadzetaL(itheta,:) = d2ZHatdthetadzetaL(itheta,:) + ZHarmonics_L(i) * BHarmonics_lL(i)*Nperiods*BHarmonics_nL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   DzL(itheta,:) = DzL(itheta,:) + DzHarmonics_L(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dDzdthetaL(itheta,:) = dDzdthetaL(itheta,:) + DzHarmonics_L(i) * BHarmonics_lL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dDzdzetaL(itheta,:) = dDzdzetaL(itheta,:) - DzHarmonics_L(i) * Nperiods * BHarmonics_nL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2Dzdtheta2L(itheta,:) = d2Dzdtheta2L(itheta,:) - DzHarmonics_L(i) * BHarmonics_lL(i)**2 * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2Dzdzeta2L(itheta,:) = d2Dzdzeta2L(itheta,:) - DzHarmonics_L(i) * (Nperiods * BHarmonics_nL(i))**2 * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2DzdthetadzetaL(itheta,:) = d2DzdthetadzetaL(itheta,:) + DzHarmonics_L(i) * BHarmonics_lL(i)*Nperiods*BHarmonics_nL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                end do
             end if
          else  ! The sine components of BHat
             include_mn=.false.
             if ((abs(BHarmonics_nL(i))<=int(Nzeta/2.0)).and.(BHarmonics_lL(i)<=int(Ntheta/2.0))) then
                include_mn=.true.
             end if
             if (BHarmonics_lL(i)==0 .or. real(BHarmonics_lL(i))==Ntheta/2.0) then
                if (BHarmonics_nL(i)==0 .or. abs(real(BHarmonics_nL(i)))==Nzeta/2.0 ) then
                   include_mn=.false.
                end if
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   BHatL(itheta,:) = BHatL(itheta,:) + BHarmonics_amplitudesL(i) * &
                   sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)
                   
                   dBHatdthetaL(itheta,:) = dBHatdthetaL(itheta,:) + BHarmonics_amplitudesL(i) * BHarmonics_lL(i) * &
                   cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)
                   
                   dBHatdzetaL(itheta,:) = dBHatdzetaL(itheta,:) - BHarmonics_amplitudesL(i) * Nperiods * BHarmonics_nL(i) * &
                   cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)
                      
                   !The following are only there to calculate Sugama's magnetic drift
                   RHatL(itheta,:) = RHatL(itheta,:) + RHarmonics_L(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dRHatdthetaL(itheta,:) = dRHatdthetaL(itheta,:) + RHarmonics_L(i) * BHarmonics_lL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dRHatdzetaL(itheta,:) = dRHatdzetaL(itheta,:) - RHarmonics_L(i) * Nperiods * BHarmonics_nL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2RHatdtheta2L(itheta,:) = d2RHatdtheta2L(itheta,:) - RHarmonics_L(i) * BHarmonics_lL(i)**2 * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2RHatdzeta2L(itheta,:) = d2RHatdzeta2L(itheta,:) - RHarmonics_L(i) * (Nperiods * BHarmonics_nL(i))**2 * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2RHatdthetadzetaL(itheta,:) = d2RHatdthetadzetaL(itheta,:) + RHarmonics_L(i) * BHarmonics_lL(i)*Nperiods*BHarmonics_nL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   ZHatL(itheta,:) = ZHatL(itheta,:) + ZHarmonics_L(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dZHatdthetaL(itheta,:) = dZHatdthetaL(itheta,:) - ZHarmonics_L(i) * BHarmonics_lL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dZHatdzetaL(itheta,:) = dZHatdzetaL(itheta,:) + ZHarmonics_L(i) * Nperiods * BHarmonics_nL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2ZHatdtheta2L(itheta,:) = d2ZHatdtheta2L(itheta,:) - ZHarmonics_L(i) * BHarmonics_lL(i)**2 * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2ZHatdzeta2L(itheta,:) = d2ZHatdzeta2L(itheta,:) - ZHarmonics_L(i) * (Nperiods * BHarmonics_nL(i))**2 * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2ZHatdthetadzetaL(itheta,:) = d2ZHatdthetadzetaL(itheta,:) + ZHarmonics_L(i) * BHarmonics_lL(i)*Nperiods*BHarmonics_nL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   DzL(itheta,:) = DzL(itheta,:) + DzHarmonics_L(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dDzdthetaL(itheta,:) = dDzdthetaL(itheta,:) - DzHarmonics_L(i) * BHarmonics_lL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   dDzdzetaL(itheta,:) = dDzdzetaL(itheta,:) + DzHarmonics_L(i) * Nperiods * BHarmonics_nL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2Dzdtheta2L(itheta,:) = d2Dzdtheta2L(itheta,:) - DzHarmonics_L(i) * BHarmonics_lL(i)**2 * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2Dzdzeta2L(itheta,:) = d2Dzdzeta2L(itheta,:) - DzHarmonics_L(i) * (Nperiods * BHarmonics_nL(i))**2 * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)

                   d2DzdthetadzetaL(itheta,:) = d2DzdthetadzetaL(itheta,:) + DzHarmonics_L(i) * BHarmonics_lL(i)*Nperiods*BHarmonics_nL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)


                end do
             end if
          end if
       end do
       allocate(BHatH(Ntheta,Nzeta))
       allocate(dBHatdthetaH(Ntheta,Nzeta))
       allocate(dBHatdzetaH(Ntheta,Nzeta))

       allocate(RHatH(Ntheta,Nzeta))
       allocate(dRHatdthetaH(Ntheta,Nzeta))
       allocate(dRHatdzetaH(Ntheta,Nzeta))
       allocate(d2RHatdtheta2H(Ntheta,Nzeta))
       allocate(d2RHatdzeta2H(Ntheta,Nzeta))
       allocate(d2RHatdthetadzetaH(Ntheta,Nzeta))

       allocate(ZHatH(Ntheta,Nzeta))
       allocate(dZHatdthetaH(Ntheta,Nzeta))
       allocate(dZHatdzetaH(Ntheta,Nzeta))
       allocate(d2ZHatdtheta2H(Ntheta,Nzeta))
       allocate(d2ZHatdzeta2H(Ntheta,Nzeta))
       allocate(d2ZHatdthetadzetaH(Ntheta,Nzeta))

       allocate(DzH(Ntheta,Nzeta))
       allocate(dDzdthetaH(Ntheta,Nzeta))
       allocate(dDzdzetaH(Ntheta,Nzeta))
       allocate(d2Dzdtheta2H(Ntheta,Nzeta))
       allocate(d2Dzdzeta2H(Ntheta,Nzeta))
       allocate(d2DzdthetadzetaH(Ntheta,Nzeta))

       BHatH = B0OverBBarH ! This includes the (0,0) component.
       dBHatdthetaH = 0
       dBHatdzetaH = 0

       RHatH = R0H ! This includes the (0,0) component.
       dRHatdthetaH      = 0
       dRHatdzetaH       = 0
       d2RHatdtheta2H    = 0
       d2RHatdzeta2H     = 0
       d2RHatdthetadzetaH= 0
       ZHatH = 0 ! This includes the (0,0) component.
       dZHatdthetaH      = 0
       dZHatdzetaH       = 0
       d2ZHatdtheta2H    = 0
       d2ZHatdzeta2H     = 0
       d2ZHatdthetadzetaH= 0
       DzH = 0 ! This includes the (0,0) component.
       dDzdthetaH      = 0
       dDzdzetaH       = 0
       d2Dzdtheta2H    = 0
       d2Dzdzeta2H     = 0
       d2DzdthetadzetaH= 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MFM scaling of non-quasisymmetric terms
       if (symm_breaking == 0) then
          print *,"There will be no scaling of symmetry-breaking terms"

       else if (symm_breaking == 1) then

          if (symm_type == 0) then
             print *,"The non-quasi-AXI-symmetric Boozer harmonics of |B| will be scaled"
             do i = 1, NHarmonicsH
                if (BHarmonics_nH(i) == 0) then
                   BHarmonics_amplitudesH(i) = BHarmonics_amplitudesH(i)
                else
                   BHarmonics_amplitudesH(i) = BHarmonics_amplitudesH(i)*epsilon_symmbreak
                end if
             end do

          else if (symm_type == 1) then
             print *,"The non-quasi-HELICALLY-symmetric Boozer harmonics of |B| will be scaled"
             do i = 1, NHarmonicsH
                if (qhs_poloidal*(-1)*BHarmonics_nH(i) == qhs_toroidal*BHarmonics_lH(i)) then
                   BHarmonics_amplitudesH(i) = BHarmonics_amplitudesH(i)
                else
                   BHarmonics_amplitudesH(i) = BHarmonics_amplitudesH(i)*epsilon_symmbreak
                end if
             end do

          else
             print *,"The type of quasisymmetry-breaking must take integer values of either 0 (QAS) or 1 (QHS) ... Exiting"
             stop
          end if
       else if (symm_breaking == 2) then
          zeta_fft_max = 2*pi/NPeriods
          dzeta_fft = zeta_fft_max/nzeta_fft
          dtheta_fft = 2.*pi/ntheta_fft
          allocate(hHat_tempH(ntheta_fft,nzeta_fft))
          allocate(hHat_temp_kH(NHarmonicsH*2+1))
          allocate(hHat_starH(ntheta_fft,nzeta_fft))
          allocate(BHat_tempH(ntheta_fft,nzeta_fft))
          allocate(BHat_temp2H(ntheta_fft,nzeta_fft))
          allocate(BHarmonics_amplitudesH_temp(NHarmonicsH*2+1))
          hHat_tempH = 0
          hHat_temp_kH = 0
          hHat_starH = 0 !1.0 / (B0OverBBar*B0OverBBar)
          BHat_tempH = 0
          BHat_temp2H = 0
          include_zero = 0
          sine_term = 0
          BHatH = 0

          zeta_fft(0) = 0
          theta_fft(0) = 0
          do i = 1, nzeta_fft
             zeta_fft(i) = dzeta_fft + zeta_fft(i-1)
          end do
          do i = 1, ntheta_fft
             theta_fft(i) = dtheta_fft + theta_fft(i-1)
          end do

          call inverse_b(BHarmonics_amplitudesH, BHat_tempH, BHarmonics_nH, BHarmonics_lH, NHarmonicsH, include_zero, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)
          BHat_tempH = BHat_tempH + B0OverBBarH
          include_zero = 1
          sine_term = 1

          do itheta = 1, ntheta_fft
             do izeta = 1, nzeta_fft
                hHat_tempH(itheta,izeta) = 1.0/(BHat_tempH(itheta,izeta)*BHat_tempH(itheta,izeta))
             end do
          end do
          call forward_b(hHat_tempH, hHat_temp_kH, BHarmonics_nH, BHarmonics_lH, NHarmonicsH, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)

          if (symm_type == 1) then ! QHS      
             do i = 2, NHarmonicsH+1 ! The (0,0) mode is not scaled
                if (qhs_poloidal*(-1)*BHarmonics_nH(i-1) == qhs_toroidal*BHarmonics_lH(i-1)) then
                   hHat_temp_kH(i) = hHat_temp_kH(i)
                else
                   hHat_temp_kH(i) = epsilon_symmbreak*hHat_temp_kH(i)
                end if
             end do
          else
             do i = 2, NHarmonicsH+1 ! The (0,0) mode is not scaled
                if (BHarmonics_nH(i-1)==0) then ! if mode is QAS
                   hHat_temp_kH(i) = hHat_temp_kH(i)
                else
                   hHat_temp_kH(i) = epsilon_symmbreak*hHat_temp_kH(i)
                end if
             end do
          end if
          
          call inverse_b(hHat_temp_kH, hHat_starH, BHarmonics_nH, BHarmonics_lH, NHarmonicsH, include_zero, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)
          do itheta = 1, ntheta_fft
             do izeta = 1, nzeta_fft
                BHat_temp2H(itheta,izeta) = sqrt(1.0/hHat_starH(itheta,izeta))
             end do
          end do

          BHarmonics_amplitudesH_temp = 0
          BHarmonics_amplitudesH = 0
          sine_term = 0
          call forward_b(BHat_temp2H, BHarmonics_amplitudesH_temp, BHarmonics_nH, BHarmonics_lH, NHarmonicsH, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)
          
          BHatH = BHarmonics_amplitudesH_temp(1) ! (0,0) mode amplitude
          do i = 2, NHarmonicsH+1
             if (abs(BHarmonics_amplitudesH_temp(i)) < 1.e-10) then
                BHarmonics_amplitudesH(i-1) = 0
             else
                BHarmonics_amplitudesH(i-1) = BHarmonics_amplitudesH_temp(i)
             end if
          end do

          deallocate(hHat_tempH)
          deallocate(hHat_temp_kH)
          deallocate(hHat_starH)
          deallocate(BHat_tempH)
          deallocate(BHat_temp2H)
          deallocate(BHarmonics_amplitudesH_temp)
       else
          print *,"You have selected an invalid integer for the method of scaling symmetry-breaking terms. Input must be 0 (no scaling), 1 (scaling through |B|), or 2 (scaling through invB^2) ... Exiting"
          stop
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do i = 1, NHarmonicsH
!          print *,"BHarmonicsH(i) = ",BHarmonics_amplitudesH(i)
          if (BHarmonics_parityH(i)) then   ! The cosine components of BHat
             include_mn = .false.
             if ((abs(BHarmonics_nH(i))<=int(Nzeta/2.0)).and.(BHarmonics_lH(i)<=int(Ntheta/2.0))) then
                include_mn = .true.
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   BHatH(itheta,:) = BHatH(itheta,:) + BHarmonics_amplitudesH(i) * &
                   cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)
                      
                   dBHatdthetaH(itheta,:) = dBHatdthetaH(itheta,:) - BHarmonics_amplitudesH(i) * BHarmonics_lH(i) * &
                   sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)
                      
                   dBHatdzetaH(itheta,:) = dBHatdzetaH(itheta,:) + BHarmonics_amplitudesH(i) * Nperiods * BHarmonics_nH(i) * &
                   sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   !The following are only there to calculate Sugama's magnetic drift
                   RHatH(itheta,:) = RHatH(itheta,:) + RHarmonics_H(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dRHatdthetaH(itheta,:) = dRHatdthetaH(itheta,:) - RHarmonics_H(i) * BHarmonics_lH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dRHatdzetaH(itheta,:) = dRHatdzetaH(itheta,:) + RHarmonics_H(i) * Nperiods * BHarmonics_nH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2RHatdtheta2H(itheta,:) = d2RHatdtheta2H(itheta,:) - RHarmonics_H(i) * BHarmonics_lH(i)**2 * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2RHatdzeta2H(itheta,:) = d2RHatdzeta2H(itheta,:) - RHarmonics_H(i) * (Nperiods * BHarmonics_nH(i))**2 * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2RHatdthetadzetaH(itheta,:) = d2RHatdthetadzetaH(itheta,:) + RHarmonics_H(i) * BHarmonics_lH(i)*Nperiods*BHarmonics_nH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   ZHatH(itheta,:) = ZHatH(itheta,:) + ZHarmonics_H(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dZHatdthetaH(itheta,:) = dZHatdthetaH(itheta,:) + ZHarmonics_H(i) * BHarmonics_lH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dZHatdzetaH(itheta,:) = dZHatdzetaH(itheta,:) - ZHarmonics_H(i) * Nperiods * BHarmonics_nH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2ZHatdtheta2H(itheta,:) = d2ZHatdtheta2H(itheta,:) - ZHarmonics_H(i) * BHarmonics_lH(i)**2 * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2ZHatdzeta2H(itheta,:) = d2ZHatdzeta2H(itheta,:) - ZHarmonics_H(i) * (Nperiods * BHarmonics_nH(i))**2 * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2ZHatdthetadzetaH(itheta,:) = d2ZHatdthetadzetaH(itheta,:) + ZHarmonics_H(i) * BHarmonics_lH(i)*Nperiods*BHarmonics_nH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   DzH(itheta,:) = DzH(itheta,:) + DzHarmonics_H(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dDzdthetaH(itheta,:) = dDzdthetaH(itheta,:) + DzHarmonics_H(i) * BHarmonics_lH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dDzdzetaH(itheta,:) = dDzdzetaH(itheta,:) - DzHarmonics_H(i) * Nperiods * BHarmonics_nH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2Dzdtheta2H(itheta,:) = d2Dzdtheta2H(itheta,:) - DzHarmonics_H(i) * BHarmonics_lH(i)**2 * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2Dzdzeta2H(itheta,:) = d2Dzdzeta2H(itheta,:) - DzHarmonics_H(i) * (Nperiods * BHarmonics_nH(i))**2 * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2DzdthetadzetaH(itheta,:) = d2DzdthetadzetaH(itheta,:) + DzHarmonics_H(i) * BHarmonics_lH(i)*Nperiods*BHarmonics_nH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)


                end do
             end if
          else  ! The sine components of BHat
             include_mn=.false.
             if ((abs(BHarmonics_nH(i))<=int(Nzeta/2.0)).and.(BHarmonics_lH(i)<=int(Ntheta/2.0))) then
                include_mn=.true.
             end if
             if (BHarmonics_lH(i)==0 .or. real(BHarmonics_lH(i))==Ntheta/2.0) then
                if (BHarmonics_nH(i)==0 .or. abs(real(BHarmonics_nH(i)))==Nzeta/2.0 ) then
                   include_mn=.false.
                end if
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   BHatH(itheta,:) = BHatH(itheta,:) + BHarmonics_amplitudesH(i) * &
                   sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)
                      
                   dBHatdthetaH(itheta,:) = dBHatdthetaH(itheta,:) + BHarmonics_amplitudesH(i) * BHarmonics_lH(i) * &
                   cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)
                      
                   dBHatdzetaH(itheta,:) = dBHatdzetaH(itheta,:) - BHarmonics_amplitudesH(i) * Nperiods * BHarmonics_nH(i) * &
                   cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   !The following are only there to calculate Sugama's magnetic drift
                   RHatH(itheta,:) = RHatH(itheta,:) + RHarmonics_H(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dRHatdthetaH(itheta,:) = dRHatdthetaH(itheta,:) + RHarmonics_H(i) * BHarmonics_lH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dRHatdzetaH(itheta,:) = dRHatdzetaH(itheta,:) - RHarmonics_H(i) * Nperiods * BHarmonics_nH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2RHatdtheta2H(itheta,:) = d2RHatdtheta2H(itheta,:) - RHarmonics_H(i) * BHarmonics_lH(i)**2 * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2RHatdzeta2H(itheta,:) = d2RHatdzeta2H(itheta,:) - RHarmonics_H(i) * (Nperiods * BHarmonics_nH(i))**2 * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2RHatdthetadzetaH(itheta,:) = d2RHatdthetadzetaH(itheta,:) + RHarmonics_H(i) * BHarmonics_lH(i)*Nperiods*BHarmonics_nH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   ZHatH(itheta,:) = ZHatH(itheta,:) + ZHarmonics_H(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dZHatdthetaH(itheta,:) = dZHatdthetaH(itheta,:) - ZHarmonics_H(i) * BHarmonics_lH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dZHatdzetaH(itheta,:) = dZHatdzetaH(itheta,:) + ZHarmonics_H(i) * Nperiods * BHarmonics_nH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2ZHatdtheta2H(itheta,:) = d2ZHatdtheta2H(itheta,:) - ZHarmonics_H(i) * BHarmonics_lH(i)**2 * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2ZHatdzeta2H(itheta,:) = d2ZHatdzeta2H(itheta,:) - ZHarmonics_H(i) * (Nperiods * BHarmonics_nH(i))**2 * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2ZHatdthetadzetaH(itheta,:) = d2ZHatdthetadzetaH(itheta,:) + ZHarmonics_H(i) * BHarmonics_lH(i)*Nperiods*BHarmonics_nH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   DzH(itheta,:) = DzH(itheta,:) + DzHarmonics_H(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dDzdthetaH(itheta,:) = dDzdthetaH(itheta,:) - DzHarmonics_H(i) * BHarmonics_lH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   dDzdzetaH(itheta,:) = dDzdzetaH(itheta,:) + DzHarmonics_H(i) * Nperiods * BHarmonics_nH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2Dzdtheta2H(itheta,:) = d2Dzdtheta2H(itheta,:) - DzHarmonics_H(i) * BHarmonics_lH(i)**2 * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2Dzdzeta2H(itheta,:) = d2Dzdzeta2H(itheta,:) - DzHarmonics_H(i) * (Nperiods * BHarmonics_nH(i))**2 * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                   d2DzdthetadzetaH(itheta,:) = d2DzdthetadzetaH(itheta,:) + DzHarmonics_H(i) * BHarmonics_lH(i)*Nperiods*BHarmonics_nH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)

                end do
             end if
          end if
       end do
       do itheta = 1,Ntheta
          BHat(itheta,:) = BHatL(itheta,:)*RadialWeight + BHatH(itheta,:)*(1.0-RadialWeight)
          dBHatdtheta(itheta,:) = dBHatdthetaL(itheta,:)*RadialWeight + dBHatdthetaH(itheta,:)*(1.0-RadialWeight)
          dBHatdzeta(itheta,:)  =  dBHatdzetaL(itheta,:)*RadialWeight +  dBHatdzetaH(itheta,:)*(1.0-RadialWeight)
          dBHatdpsiHat(itheta,:)= (BHatH(itheta,:)-BHatL(itheta,:)) / DeltapsiHat

       end do
       
       if (nearbyRadiiGiven) then !implemented only for geometryScheme 11 and 12 so far
          allocate(RHat(Ntheta,Nzeta))
          allocate(dRHatdtheta(Ntheta,Nzeta))
          allocate(dRHatdzeta(Ntheta,Nzeta))
          allocate(d2RHatdtheta2(Ntheta,Nzeta))
          allocate(d2RHatdzeta2(Ntheta,Nzeta))
          allocate(d2RHatdthetadzeta(Ntheta,Nzeta))
          allocate(ZHat(Ntheta,Nzeta))
          allocate(dZHatdtheta(Ntheta,Nzeta))
          allocate(dZHatdzeta(Ntheta,Nzeta))
          allocate(d2ZHatdtheta2(Ntheta,Nzeta))
          allocate(d2ZHatdzeta2(Ntheta,Nzeta))
          allocate(d2ZHatdthetadzeta(Ntheta,Nzeta))
          allocate(Dz(Ntheta,Nzeta))
          allocate(dDzdtheta(Ntheta,Nzeta))
          allocate(dDzdzeta(Ntheta,Nzeta))
          allocate(d2Dzdtheta2(Ntheta,Nzeta))
          allocate(d2Dzdzeta2(Ntheta,Nzeta))
          allocate(d2Dzdthetadzeta(Ntheta,Nzeta))
          do itheta = 1,Ntheta
             RHat(itheta,:) = RHatL(itheta,:)*RadialWeight + RHatH(itheta,:)*(1.0-RadialWeight)
             dRHatdtheta(itheta,:) = dRHatdthetaL(itheta,:)*RadialWeight + dRHatdthetaH(itheta,:)*(1.0-RadialWeight)
             dRHatdzeta(itheta,:)  =  dRHatdzetaL(itheta,:)*RadialWeight +  dRHatdzetaH(itheta,:)*(1.0-RadialWeight)
             d2RHatdtheta2(itheta,:) = d2RHatdtheta2L(itheta,:)*RadialWeight + d2RHatdtheta2H(itheta,:)*(1.0-RadialWeight)
             d2RHatdzeta2(itheta,:)  =  d2RHatdzeta2L(itheta,:)*RadialWeight +  d2RHatdzeta2H(itheta,:)*(1.0-RadialWeight)
             d2RHatdthetadzeta(itheta,:) = d2RHatdthetadzetaL(itheta,:)*RadialWeight + d2RHatdthetadzetaH(itheta,:)*(1.0-RadialWeight)

             ZHat(itheta,:) = ZHatL(itheta,:)*RadialWeight + ZHatH(itheta,:)*(1.0-RadialWeight)
             dZHatdtheta(itheta,:) = dZHatdthetaL(itheta,:)*RadialWeight + dZHatdthetaH(itheta,:)*(1.0-RadialWeight)
             dZHatdzeta(itheta,:)  =  dZHatdzetaL(itheta,:)*RadialWeight +  dZHatdzetaH(itheta,:)*(1.0-RadialWeight)
             d2ZHatdtheta2(itheta,:) = d2ZHatdtheta2L(itheta,:)*RadialWeight + d2ZHatdtheta2H(itheta,:)*(1.0-RadialWeight)
             d2ZHatdzeta2(itheta,:)  =  d2ZHatdzeta2L(itheta,:)*RadialWeight +  d2ZHatdzeta2H(itheta,:)*(1.0-RadialWeight)
             d2ZHatdthetadzeta(itheta,:) = d2ZHatdthetadzetaL(itheta,:)*RadialWeight + d2ZHatdthetadzetaH(itheta,:)*(1.0-RadialWeight)

             Dz(itheta,:) = DzL(itheta,:)*RadialWeight + DzH(itheta,:)*(1.0-RadialWeight)
             dDzdtheta(itheta,:) = dDzdthetaL(itheta,:)*RadialWeight + dDzdthetaH(itheta,:)*(1.0-RadialWeight)
             dDzdzeta(itheta,:)  =  dDzdzetaL(itheta,:)*RadialWeight +  dDzdzetaH(itheta,:)*(1.0-RadialWeight)
             d2Dzdtheta2(itheta,:) = d2Dzdtheta2L(itheta,:)*RadialWeight + d2Dzdtheta2H(itheta,:)*(1.0-RadialWeight)
             d2Dzdzeta2(itheta,:)  =  d2Dzdzeta2L(itheta,:)*RadialWeight +  d2Dzdzeta2H(itheta,:)*(1.0-RadialWeight)
             d2Dzdthetadzeta(itheta,:) = d2DzdthetadzetaL(itheta,:)*RadialWeight + d2DzdthetadzetaH(itheta,:)*(1.0-RadialWeight)
          end do
       end if
    end if
    
    ! ---------------------------------------------------------------------------------------
    ! Calculate parallel current u from harmonics of 1/B^2. Used in NTV calculation.
    ! \nabla_\parallel u = (2/B^4) \nabla B \times \vector{B} \cdot \iota \nabla \psi 
    ! ---------------------------------------------------------------------------------------
    allocate(hHat(Ntheta,Nzeta))
    allocate(duHatdtheta(Ntheta,Nzeta))
    allocate(duHatdzeta(Ntheta,Nzeta))
    
    uHat = 0
    duHatdtheta = 0
    duHatdzeta = 0
    BHat_sub_psi = 0
    dBHat_sub_psi_dtheta = 0
    dBHat_sub_psi_dzeta = 0
    invFSA_BHat2 = 0
    do itheta = 1,Ntheta
       do izeta = 1,Nzeta
          hHat(itheta,izeta) = 1.0 / (BHat(itheta,izeta)*BHat(itheta,izeta))
          invFSA_BHat2 = invFSA_BHat2 + hHat(itheta,izeta)/Ntheta/Nzeta
       end do
    end do
    
    if (.not. nearbyRadiiGiven) then
       nonStelSym = any(.not. BHarmonics_parity)
    else
       nonStelSym = any(.not. BHarmonics_parityL) .or. any(.not. BHarmonics_parityH)
    end if
    if (nonStelSym) then !sine components exist
       do m = 0,int(Ntheta/2.0) !Nyquist max freq.
          if (m == 0) then
             startn=1
          else if (real(m)==Ntheta/2.0) then
             startn=0
          else if (real(int(Nzeta/2.0))==Nzeta/2.0) then
             startn=-int(Nzeta/2.0)+1
          else
             startn=-int(Nzeta/2.0)
          end if
          stopn=int(Nzeta/2.0)
          do n = startn,stopn 
             !cos
             hHatHarmonics_amplitude = 0
             do itheta = 1,Ntheta
                if ((m == 0 .and. real(n)==Nzeta/2.0) .or. (real(m)==Ntheta/2.0 .and. n==0) .or. &
                     (real(m)==Ntheta/2.0 .and. real(n)==Nzeta/2.0)) then
                   hHatHarmonics_amplitude = hHatHarmonics_amplitude + 1.0/(Ntheta*Nzeta) * &
                        dot_product(cos(m * theta(itheta)  - n * NPeriods * zeta), hHat(itheta,:))
                else
                   hHatHarmonics_amplitude = hHatHarmonics_amplitude + 2.0/(Ntheta*Nzeta) * &
                        dot_product(cos(m * theta(itheta)  - n * NPeriods * zeta), hHat(itheta,:))
                end if
             end do
             uHatHarmonics_amplitude = &
                  iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude
             dBHat_sub_psi_dthetaHarmonics_amplitude = &
                  -pPrimeHat/iota*(uHatHarmonics_amplitude-iota*IHat*hHatHarmonics_amplitude)
             dBHat_sub_psi_dzetaHarmonics_amplitude = &
                  pPrimeHat*(uHatHarmonics_amplitude+GHat*hHatHarmonics_amplitude)

             do itheta = 1,Ntheta
                uHat(itheta,:) = uHat(itheta,:) &
                     + uHatHarmonics_amplitude * cos(m * theta(itheta) - n * NPeriods * zeta)
                duHatdtheta(itheta,:) = duHatdtheta(itheta,:) &
                     - uHatHarmonics_amplitude * m * sin(m * theta(itheta) - n * NPeriods * zeta)
                duHatdzeta(itheta,:) = duHatdzeta(itheta,:) &
                     + uHatHarmonics_amplitude * n * NPeriods * sin(m * theta(itheta) - n * NPeriods * zeta)
                dBHat_sub_psi_dtheta(itheta,:) = dBHat_sub_psi_dtheta(itheta,:) &
		     + dBHat_sub_psi_dthetaHarmonics_amplitude * cos(m * theta(itheta) - n * NPeriods * zeta)
                dBHat_sub_psi_dzeta(itheta,:) = dBHat_sub_psi_dzeta(itheta,:) &
		     + dBHat_sub_psi_dzetaHarmonics_amplitude * cos(m * theta(itheta) - n * NPeriods * zeta)
             end do
             if (n==0) then
                do itheta = 1,Ntheta
                   BHat_sub_psi(itheta,:) = BHat_sub_psi(itheta,:) &
                        + dBHat_sub_psi_dthetaHarmonics_amplitude / m &
                        * sin(m * theta(itheta) - n * NPeriods * zeta)
                end do
             else
                do itheta = 1,Ntheta
                   BHat_sub_psi(itheta,:) = BHat_sub_psi(itheta,:) &
                        - dBHat_sub_psi_dzetaHarmonics_amplitude / n / NPeriods  &
                        * sin(m * theta(itheta) - n * NPeriods * zeta)
                end do
             end if

             !sin
             hHatHarmonics_amplitude = 0
             if (.not.((m == 0 .and. real(n)==Nzeta/2.0) .or. (real(m)==Ntheta/2.0 .and. n==0) .or. &
                     (real(m)==Ntheta/2.0 .and. real(n)==Nzeta/2.0))) then
                do itheta = 1,Ntheta
                   hHatHarmonics_amplitude = hHatHarmonics_amplitude + 2.0/(Ntheta*Nzeta) * &
                        dot_product(sin(m * theta(itheta)  - n * NPeriods * zeta), hHat(itheta,:))
                end do
             end if
             uHatHarmonics_amplitude = &
                  iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude
             dBHat_sub_psi_dthetaHarmonics_amplitude = &
                  -pPrimeHat/iota*(uHatHarmonics_amplitude-iota*IHat*hHatHarmonics_amplitude)
             dBHat_sub_psi_dzetaHarmonics_amplitude = &
                  pPrimeHat*(uHatHarmonics_amplitude+GHat*hHatHarmonics_amplitude)
             do itheta = 1,Ntheta
                uHat(itheta,:) = uHat(itheta,:) &
                     + uHatHarmonics_amplitude * sin(m * theta(itheta) - n * NPeriods * zeta)
                duHatdtheta(itheta,:) = duHatdtheta(itheta,:) &
                     + uHatHarmonics_amplitude * m * cos(m * theta(itheta) - n * NPeriods * zeta)
                duHatdzeta(itheta,:) = duHatdzeta(itheta,:) &
                     - uHatHarmonics_amplitude * n * NPeriods * cos(m * theta(itheta) - n * NPeriods * zeta)
                dBHat_sub_psi_dtheta(itheta,:) = dBHat_sub_psi_dtheta(itheta,:) &
		     + dBHat_sub_psi_dthetaHarmonics_amplitude * sin(m * theta(itheta) - n * NPeriods * zeta)
                dBHat_sub_psi_dzeta(itheta,:) = dBHat_sub_psi_dzeta(itheta,:) &
		     + dBHat_sub_psi_dzetaHarmonics_amplitude * sin(m * theta(itheta) - n * NPeriods * zeta)
             end do
             if (n==0) then
                do itheta = 1,Ntheta
                   BHat_sub_psi(itheta,:) = BHat_sub_psi(itheta,:) &
                        - dBHat_sub_psi_dthetaHarmonics_amplitude / m &
                        * cos(m * theta(itheta) - n * NPeriods * zeta)
                end do
             else
                do itheta = 1,Ntheta
                   BHat_sub_psi(itheta,:) = BHat_sub_psi(itheta,:) &
                        + dBHat_sub_psi_dzetaHarmonics_amplitude / n / NPeriods &
                        * cos(m * theta(itheta) - n * NPeriods * zeta)
                end do
             end if
          end do
       end do
    else !only cosinus components
       do m = 0,int(Ntheta/2.0) !Nyquist max freq.
          if (m == 0) then
             startn=1
          else if (real(m)==Ntheta/2.0) then
             startn=0
          else if (real(int(Nzeta/2.0))==Nzeta/2.0) then
             startn=-int(Nzeta/2.0)+1
          else
             startn=-int(Nzeta/2.0)
          end if
          stopn=int(Nzeta/2.0)
          do n = startn,stopn 
             !cos
             hHatHarmonics_amplitude = 0
             do itheta = 1,Ntheta
                if ((m == 0 .and. real(n)==Nzeta/2.0) .or. (real(m)==Ntheta/2.0 .and. n==0) .or. &
                     (real(m)==Ntheta/2.0 .and. real(n)==Nzeta/2.0)) then
                   hHatHarmonics_amplitude = hHatHarmonics_amplitude + 1.0/(Ntheta*Nzeta) * &
                        dot_product(cos(m * theta(itheta)  - n * NPeriods * zeta), hHat(itheta,:))
                else
                   hHatHarmonics_amplitude = hHatHarmonics_amplitude + 2.0/(Ntheta*Nzeta) * &
                        dot_product(cos(m * theta(itheta)  - n * NPeriods * zeta), hHat(itheta,:))
                end if
             end do
             uHatHarmonics_amplitude = &
                  iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude
             dBHat_sub_psi_dthetaHarmonics_amplitude = &
                  -pPrimeHat/iota*(uHatHarmonics_amplitude-iota*IHat*hHatHarmonics_amplitude)
             dBHat_sub_psi_dzetaHarmonics_amplitude = &
                  pPrimeHat*(uHatHarmonics_amplitude+GHat*hHatHarmonics_amplitude)
             do itheta = 1,Ntheta
                uHat(itheta,:) = uHat(itheta,:) &
                     + uHatHarmonics_amplitude * cos(m * theta(itheta) - n * NPeriods * zeta)
                duHatdtheta(itheta,:) = duHatdtheta(itheta,:) &
                     - uHatHarmonics_amplitude * m * sin(m * theta(itheta) - n * NPeriods * zeta)
                duHatdzeta(itheta,:) = duHatdzeta(itheta,:) &
                     + uHatHarmonics_amplitude * n * NPeriods * sin(m * theta(itheta) - n * NPeriods * zeta)
                dBHat_sub_psi_dtheta(itheta,:) = dBHat_sub_psi_dtheta(itheta,:) &
		     + dBHat_sub_psi_dthetaHarmonics_amplitude * cos(m * theta(itheta) - n * NPeriods * zeta)
                dBHat_sub_psi_dzeta(itheta,:) = dBHat_sub_psi_dzeta(itheta,:) &
		     + dBHat_sub_psi_dzetaHarmonics_amplitude * cos(m * theta(itheta) - n * NPeriods * zeta)
             end do
             if (n==0) then
                do itheta = 1,Ntheta
                   BHat_sub_psi(itheta,:) = BHat_sub_psi(itheta,:) &
                        + dBHat_sub_psi_dthetaHarmonics_amplitude / m &
                        * sin(m * theta(itheta) - n * NPeriods * zeta)
                end do
             else
                do itheta = 1,Ntheta
                   BHat_sub_psi(itheta,:) = BHat_sub_psi(itheta,:) &
                        - dBHat_sub_psi_dzetaHarmonics_amplitude / n / NPeriods  &
                        * sin(m * theta(itheta) - n * NPeriods * zeta)
                end do
             end if
          end do
       end do
    end if
   
    ! This method is also right, but I prefer the above one to make sure that the 00 comp. is 0
    !do itheta = 1,Ntheta
    !   dBHat_sub_psi_dtheta(itheta,:) = -pPrimeHat/iota*(uHat(itheta,:)-iota*IHat*(hHat(itheta,:)-1/FSABHat2))
    !   dBHat_sub_psi_dzeta(itheta,:)  =  pPrimeHat*(uHat(itheta,:)+GHat*(hHat(itheta,:)-1/FSABHat2))
    !end do

    do itheta = 1,Ntheta
       do izeta = 1,Nzeta
          NTVKernel(itheta,izeta) = 2.0/5.0 / BHat(itheta,izeta) * ( &
               (uHat(itheta,izeta) - GHat*invFSA_BHat2) * (iota * dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta)) &
               + iota * hHat(itheta,izeta) * (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta)) )
       end do
    end do

    ! ---------------------------------------------------------------------------------------
    ! Calculate the normal curvature, which is needed in Sugama's magnetic drift
    ! ---------------------------------------------------------------------------------------
    if (.not. nearbyRadiiGiven) then
       gradpsidotgradB_overgpsipsi=0 !implemented only for geometryScheme 11 and 12 so far
    else  
       allocate(geomang(Ntheta,Nzeta))
       allocate(dgeomangdtheta(Ntheta,Nzeta))
       allocate(dgeomangdzeta(Ntheta,Nzeta))
       allocate(d2geomangdtheta2(Ntheta,Nzeta))
       allocate(d2geomangdzeta2(Ntheta,Nzeta))
       allocate(d2geomangdthetadzeta(Ntheta,Nzeta))
       allocate(dXdtheta(Ntheta,Nzeta))
       allocate(dXdzeta(Ntheta,Nzeta))
       allocate(dYdtheta(Ntheta,Nzeta))
       allocate(dYdzeta(Ntheta,Nzeta))
       allocate(d2Xdtheta2(Ntheta,Nzeta))
       allocate(d2Xdthetadzeta(Ntheta,Nzeta))
       allocate(d2Xdzeta2(Ntheta,Nzeta))
       allocate(d2Ydtheta2(Ntheta,Nzeta))
       allocate(d2Ydthetadzeta(Ntheta,Nzeta))
       allocate(d2Ydzeta2(Ntheta,Nzeta))
       allocate(gradpsiX(Ntheta,Nzeta))
       allocate(gradpsiY(Ntheta,Nzeta))
       allocate(gradpsiZ(Ntheta,Nzeta))
       !allocate(gpsipsi(Ntheta,Nzeta)) Has been upgraded to a global variable 2018.02.21/HS
       allocate(CX(Ntheta,Nzeta))
       allocate(CY(Ntheta,Nzeta))
       allocate(CZ(Ntheta,Nzeta))

       do itheta = 1,Ntheta
          geomang(itheta,:) = Dz(itheta,:) - zeta !geometric toroidal angle, (R, geomang, Z) is right handed
          dgeomangdtheta(itheta,:)      = dDzdtheta(itheta,:)
          dgeomangdzeta(itheta,:)       = dDzdzeta(itheta,:) - 1
          d2geomangdtheta2(itheta,:)    = d2Dzdtheta2(itheta,:)
          d2geomangdzeta2(itheta,:)     = d2Dzdzeta2(itheta,:)
          d2geomangdthetadzeta(itheta,:)= d2Dzdthetadzeta(itheta,:)

          dXdtheta(itheta,:)=dRHatdtheta(itheta,:)*cos(geomang(itheta,:))-RHat(itheta,:)*dgeomangdtheta(itheta,:)*sin(geomang(itheta,:))
          dXdzeta(itheta,:) =dRHatdzeta(itheta,:) *cos(geomang(itheta,:))-RHat(itheta,:)*dgeomangdzeta(itheta,:) *sin(geomang(itheta,:))
          dYdtheta(itheta,:)=dRHatdtheta(itheta,:)*sin(geomang(itheta,:))+RHat(itheta,:)*dgeomangdtheta(itheta,:)*cos(geomang(itheta,:))
          dYdzeta(itheta,:) =dRHatdzeta(itheta,:) *sin(geomang(itheta,:))+RHat(itheta,:)*dgeomangdzeta(itheta,:) *cos(geomang(itheta,:))

          d2Xdtheta2(itheta,:)=d2RHatdtheta2(itheta,:)*cos(geomang(itheta,:)) &
          -2*dRHatdtheta(itheta,:)*dgeomangdtheta(itheta,:)*sin(geomang(itheta,:))&
          -RHat(itheta,:)*d2geomangdtheta2(itheta,:)*sin(geomang(itheta,:)) &
          -RHat(itheta,:)*dgeomangdtheta(itheta,:)**2*cos(geomang(itheta,:))
          d2Xdthetadzeta(itheta,:)=d2RHatdthetadzeta(itheta,:)*cos(geomang(itheta,:)) &
          -(dRHatdtheta(itheta,:)*dgeomangdzeta(itheta,:)+dRHatdzeta(itheta,:)*dgeomangdtheta(itheta,:))*sin(geomang(itheta,:)) &
          -RHat(itheta,:)*d2geomangdthetadzeta(itheta,:)*sin(geomang(itheta,:)) &
          -RHat(itheta,:)*dgeomangdtheta(itheta,:)*dgeomangdzeta(itheta,:)*cos(geomang(itheta,:))
          d2Xdzeta2(itheta,:)=d2RHatdzeta2(itheta,:)*cos(geomang(itheta,:)) &
          -2*dRHatdzeta(itheta,:)*dgeomangdzeta(itheta,:)*sin(geomang(itheta,:))&
          -RHat(itheta,:)*d2geomangdzeta2(itheta,:)*sin(geomang(itheta,:)) &
          -RHat(itheta,:)*dgeomangdzeta(itheta,:)**2*cos(geomang(itheta,:))

          d2Ydtheta2(itheta,:)=d2RHatdtheta2(itheta,:)*sin(geomang(itheta,:)) &
          +2*dRHatdtheta(itheta,:)*dgeomangdtheta(itheta,:)*cos(geomang(itheta,:))&
          +RHat(itheta,:)*d2geomangdtheta2(itheta,:)*cos(geomang(itheta,:)) &
          -RHat(itheta,:)*dgeomangdtheta(itheta,:)**2*sin(geomang(itheta,:))
          d2Ydthetadzeta(itheta,:)=d2RHatdthetadzeta(itheta,:)*sin(geomang(itheta,:)) &
          +(dRHatdtheta(itheta,:)*dgeomangdzeta(itheta,:)+dRHatdzeta(itheta,:)*dgeomangdtheta(itheta,:))*cos(geomang(itheta,:)) &
          +RHat(itheta,:)*d2geomangdthetadzeta(itheta,:)*cos(geomang(itheta,:)) &
          -RHat(itheta,:)*dgeomangdtheta(itheta,:)*dgeomangdzeta(itheta,:)*sin(geomang(itheta,:))
          d2Ydzeta2(itheta,:)=d2RHatdzeta2(itheta,:)*sin(geomang(itheta,:)) &
          +2*dRHatdzeta(itheta,:)*dgeomangdzeta(itheta,:)*cos(geomang(itheta,:))&
          +RHat(itheta,:)*d2geomangdzeta2(itheta,:)*cos(geomang(itheta,:)) &
          -RHat(itheta,:)*dgeomangdzeta(itheta,:)**2*sin(geomang(itheta,:))

          gradpsiX(itheta,:)=BHat(itheta,:)*BHat(itheta,:)/(GHat+iota*IHat)* &
               (dYdtheta(itheta,:)*dZHatdzeta(itheta,:)-dZHatdtheta(itheta,:)*dYdzeta(itheta,:))
          gradpsiY(itheta,:)=BHat(itheta,:)*BHat(itheta,:)/(GHat+iota*IHat)* &
               (dZHatdtheta(itheta,:)*dXdzeta(itheta,:)-dXdtheta(itheta,:)*dZHatdzeta(itheta,:))
          gradpsiZ(itheta,:)=BHat(itheta,:)*BHat(itheta,:)/(GHat+iota*IHat)* &
               (dXdtheta(itheta,:)*dYdzeta(itheta,:)-dYdtheta(itheta,:)*dXdzeta(itheta,:))
          gpsipsi(itheta,:)=gradpsiX(itheta,:)*gradpsiX(itheta,:)+&
                            gradpsiY(itheta,:)*gradpsiY(itheta,:)+&
                            gradpsiZ(itheta,:)*gradpsiZ(itheta,:)

          CX(itheta,:)=(d2Xdzeta2(itheta,:)+2*iota*d2Xdthetadzeta(itheta,:)+iota**2*d2Xdtheta2(itheta,:))&
               *(BHat(itheta,:)**2/(GHat+iota*IHat))**2
          CY(itheta,:)=(d2Ydzeta2(itheta,:)+2*iota*d2Ydthetadzeta(itheta,:)+iota**2*d2Ydtheta2(itheta,:))&
               *(BHat(itheta,:)**2/(GHat+iota*IHat))**2
          CZ(itheta,:)=(d2ZHatdzeta2(itheta,:)+2*iota*d2ZHatdthetadzeta(itheta,:)+iota**2*d2ZHatdtheta2(itheta,:))&
               *(BHat(itheta,:)**2/(GHat+iota*IHat))**2

          !normal_curvature(itheta,:)=1/(BHat(itheta,:)**2*sqrt(gpsipsi(itheta,:))) * &
          !                           CX(itheta,:)*gradpsi_X(itheta,:)+ &
          !                           CY(itheta,:)*gradpsi_Y(itheta,:)+ &
          !                           CZ(itheta,:)*gradpsi_Z(itheta,:)
          gradpsidotgradB_overgpsipsi(itheta,:) = (CX(itheta,:)*gradpsiX(itheta,:)+ &
                                                   CY(itheta,:)*gradpsiY(itheta,:)+ &
                                                   CZ(itheta,:)*gradpsiZ(itheta,:)) &
                                                   /(BHat(itheta,:)*gpsipsi(itheta,:)) &
                                                  - pPrimeHat/BHat(itheta,:)
       end do
       
       deallocate(geomang)
       deallocate(dgeomangdtheta)
       deallocate(dgeomangdzeta)
       deallocate(d2geomangdtheta2)
       deallocate(d2geomangdzeta2)
       deallocate(d2geomangdthetadzeta)
       deallocate(dXdtheta)
       deallocate(dXdzeta)
       deallocate(dYdtheta)
       deallocate(dYdzeta)
       deallocate(d2Xdtheta2)
       deallocate(d2Xdthetadzeta)
       deallocate(d2Xdzeta2)
       deallocate(d2Ydtheta2)
       deallocate(d2Ydthetadzeta)
       deallocate(d2Ydzeta2)
       deallocate(gradpsiX)
       deallocate(gradpsiY)
       deallocate(gradpsiZ)
       !deallocate(gpsipsi)
       deallocate(CX)
       deallocate(CY)
       deallocate(CZ)
    end if
    
    if (.not. nearbyRadiiGiven) then
       deallocate(BHarmonics_l)
       deallocate(BHarmonics_n)
       deallocate(BHarmonics_amplitudes)
       deallocate(BHarmonics_parity)
    else
       deallocate(BHarmonics_lL)
       deallocate(BHarmonics_nL)
       deallocate(BHarmonics_amplitudesL)
       deallocate(BHarmonics_parityL)
       deallocate(BHarmonics_lH)
       deallocate(BHarmonics_nH)
       deallocate(BHarmonics_amplitudesH)
       deallocate(BHarmonics_parityH)
       deallocate(RHatL)
       deallocate(dRHatdthetaL)
       deallocate(dRHatdzetaL)
       deallocate(d2RHatdtheta2L)
       deallocate(d2RHatdzeta2L)
       deallocate(d2RHatdthetadzetaL)
       deallocate(ZHatL)
       deallocate(dZHatdthetaL)
       deallocate(dZHatdzetaL)
       deallocate(d2ZHatdtheta2L)
       deallocate(d2ZHatdzeta2L)
       deallocate(d2ZHatdthetadzetaL)
       deallocate(DzL)
       deallocate(dDzdthetaL)
       deallocate(dDzdzetaL)
       deallocate(d2Dzdtheta2L)
       deallocate(d2Dzdzeta2L)
       deallocate(d2DzdthetadzetaL)
       deallocate(RHatH)
       deallocate(dRHatdthetaH)
       deallocate(dRHatdzetaH)
       deallocate(d2RHatdtheta2H)
       deallocate(d2RHatdzeta2H)
       deallocate(d2RHatdthetadzetaH)
       deallocate(ZHatH)
       deallocate(dZHatdthetaH)
       deallocate(dZHatdzetaH)
       deallocate(d2ZHatdtheta2H)
       deallocate(d2ZHatdzeta2H)
       deallocate(d2ZHatdthetadzetaH)
       deallocate(DzH)
       deallocate(dDzdthetaH)
       deallocate(dDzdzetaH)
       deallocate(d2Dzdtheta2H)
       deallocate(d2Dzdzeta2H)
       deallocate(d2DzdthetadzetaH)
    end if

    ! Set the Jacobian and various other components of B:

    DHat = BHat * BHat / (GHat + iota * IHat)
    BHat_sup_theta = iota * DHat
    BHat_sup_zeta = DHat
    BHat_sub_theta = IHat
    BHat_sub_zeta = GHat

    if (nearbyRadiiGiven) then
       dBHat_sup_zeta_dpsiHat = 2.0 * BHat *dBHatdpsiHat / (GHat + iota * IHat) &
            -(dBHat_sub_zeta_dpsiHat + iota*dBHat_sub_theta_dpsiHat + diotadpsiHat*IHat) &
            / (GHat + iota * IHat) / (GHat + iota * IHat)
       dBHat_sup_zeta_dtheta = 2.0 * BHat *dBHatdtheta / (GHat + iota * IHat)
       
       dBHat_sup_theta_dpsiHat = iota * dBHat_sup_zeta_dpsiHat + diotadpsiHat* DHat
       dBHat_sup_theta_dzeta = iota * 2.0 * BHat *dBHatdzeta / (GHat + iota * IHat)
    end if

    !possible double-check
    !dBHatdpsiHat= sqrt(hHat)/2.0*( dBHat_sup_theta_dpsiHat*BHat_sub_theta &
    !                              +dBHat_sub_theta_dpsiHat*BHat_sup_theta &
    !                              +dBHat_sup_zeta_dpsiHat *BHat_sub_zeta &
    !                              +dBHat_sub_zeta_dpsiHat *BHat_sup_zeta)
  end subroutine computeBHat_Boozer

  ! -----------------------------------------------------------------------------------------

  subroutine computeBHat_VMEC

    implicit none

    integer :: vmecRadialIndex_full(2)
    integer :: vmecRadialIndex_half(2)
    PetscScalar :: vmecRadialWeight_full(2)
    PetscScalar :: vmecRadialWeight_half(2)
    PetscScalar, dimension(:), allocatable :: dr2, psiN_full, psiN_half
    PetscScalar, dimension(:), allocatable :: vmec_dBHatdpsiHat, vmec_dBHat_sub_theta_dpsiHat, vmec_dBHat_sub_zeta_dpsiHat
    PetscScalar, dimension(:), allocatable :: vmec_dRdpsiHat, vmec_dZdpsiHat, vmec_dpdpsiHat
    integer :: i, j, index, isurf, itheta, izeta, m, n, imn, imn_nyq
    PetscScalar :: min_dr2, angle, sin_angle, cos_angle, b, b00, temp, dphi, dpsi
    integer :: numSymmetricModesIncluded, numAntisymmetricModesIncluded
    PetscScalar :: scaleFactor
    PetscScalar, dimension(:,:), allocatable :: R, dRdtheta, dRdzeta, dRdpsiHat, dZdtheta, dZdzeta, dZdpsiHat
    PetscScalar, dimension(:,:), allocatable :: dXdtheta, dXdzeta, dXdpsiHat, dYdtheta, dYdzeta, dYdpsiHat
    PetscScalar, dimension(:,:), allocatable :: g_sub_theta_theta, g_sub_theta_zeta, g_sub_zeta_zeta, g_sub_psi_theta, g_sub_psi_zeta, g_sub_psi_psi
    logical :: non_Nyquist_mode_available, found_imn


    ! This subroutine is written so that only psiN_wish is used, not the other *_wish quantities.

    if (masterProc) then
       print *,"Reading VMEC geometry from file ",trim(equilibriumFile)
    end if

    allocate(vmec_dBHatdpsiHat(ns))
    allocate(vmec_dBHat_sub_theta_dpsiHat(ns))
    allocate(vmec_dBHat_sub_zeta_dpsiHat(ns))
    allocate(vmec_dRdpsiHat(ns))
    allocate(vmec_dZdpsiHat(ns))
    allocate(vmec_dpdpsiHat(ns))

    allocate(R(Ntheta,Nzeta))
    allocate(dRdtheta(Ntheta,Nzeta))
    allocate(dRdzeta(Ntheta,Nzeta))
    allocate(dRdpsiHat(Ntheta,Nzeta))

    allocate(dXdtheta(Ntheta,Nzeta))
    allocate(dXdzeta(Ntheta,Nzeta))
    allocate(dXdpsiHat(Ntheta,Nzeta))

    allocate(dYdtheta(Ntheta,Nzeta))
    allocate(dYdzeta(Ntheta,Nzeta))
    allocate(dYdpsiHat(Ntheta,Nzeta))

    allocate(dZdtheta(Ntheta,Nzeta))
    allocate(dZdzeta(Ntheta,Nzeta))
    allocate(dZdpsiHat(Ntheta,Nzeta))

    allocate(g_sub_theta_theta(Ntheta,Nzeta))
    allocate(g_sub_theta_zeta(Ntheta,Nzeta))
    allocate(g_sub_zeta_zeta(Ntheta,Nzeta))
    allocate(g_sub_psi_theta(Ntheta,Nzeta))
    allocate(g_sub_psi_zeta(Ntheta,Nzeta))
    allocate(g_sub_psi_psi(Ntheta,Nzeta))

    ! There is a bug in libstell read_wout_file for ASCII-format wout files, in which the xm_nyq and xn_nyq arrays are sometimes
    ! not populated. The next few lines here provide a workaround:
    if (maxval(abs(xm_nyq)) < 1 .and. maxval(abs(xn_nyq)) < 1) then
       if (mnmax_nyq == mnmax) then
          if (masterProc) print *,"xm_nyq and xn_nyq arrays are not populated in the wout file. Using xm and xn instead."
          xm_nyq = xm
          xn_nyq = xn
       else
          if (masterProc) print *,"Error! xm_nyq and xn_nyq arrays are not populated in the wout file, and mnmax_nyq != mnmax."
          stop
       end if
    end if

    ! --------------------------------------------------------------------------------
    ! Do some sanity checking to ensure the VMEC arrays have some expected properties.
    ! --------------------------------------------------------------------------------

    if (abs(phi(1)) > 1d-14) then
       if (masterProc) then
          print *,"Error! VMEC phi array does not begin with 0."
       end if
       stop
    end if

    dphi = phi(2) - phi(1)
    do j=3,ns
       if (abs(phi(j)-phi(j-1)-dphi) > 1d-11) then
          if (masterProc) then
             print *,"Error! VMEC phi array is not uniformly spaced."
          end if
          stop
       end if
    end do

    ! The variable called 'phips' in the wout file is called just 'phip' in read_wout_mod.F.
    ! phips is on the half-mesh, so skip first point.
    do j=2,ns
       if (abs(phip(j)+phi(ns)/(2*pi)) > 1d-11) then
          if (masterProc) then
             print *,"Error! VMEC phips array is not constant and equal to -phi(ns)/(2*pi)."
          end if
          stop
       end if
    end do

    ! The first mode in the m and n arrays should be m=n=0:
    if (xm(1) .ne. 0) stop "First element of xm in the wout file should be 0."
    if (xn(1) .ne. 0) stop "First element of xn in the wout file should be 0."
    if (xm_nyq(1) .ne. 0) stop "First element of xm_nyq in the wout file should be 0."
    if (xn_nyq(1) .ne. 0) stop "First element of xn_nyq in the wout file should be 0."

    ! --------------------------------------------------------------------------------
    ! End of sanity checks.
    ! --------------------------------------------------------------------------------

    allocate(psiN_full(ns))
    psiN_full = phi / (2*pi*psiAHat)

    ! Build an array of the half grid
    allocate(psiN_half(ns-1))
    do i = 1,ns-1
       psiN_half(i) = (psiN_full(i) + psiN_full(i+1))/two
    end do

    ! --------------------------------------------------------------------------------
    ! Now choose the "actual" radius to use, based on psiN_wish:
    ! --------------------------------------------------------------------------------

    ! VMECRadialOption
    ! 0 = use exact radius requested.
    ! 1 = use nearest value of the VMEC half grid.
    ! 2 = use nearest value of the VMEC full grid.  I'm not sure why you would ever want to choose this option,
    !     but I've implemented it for completeness.
    
    select case (VMECRadialOption)
    case (0)
       ! Use exact radius requested.
       psiN = psiN_wish

    case (1)
       ! Use nearest value of the VMEC half grid

       ! Compute differences
       allocate(dr2(ns-1))
       dr2 = (psiN_half - psiN_wish) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,ns-1
          if (dr2(j)<min_dr2) then
             index = j
             min_dr2 = dr2(j)
          end if
       end do

       psiN = psiN_half(index)
       deallocate(dr2)

    case (2)
       ! Use nearest value of the VMEC full grid

       ! Compute differences
       allocate(dr2(ns))
       dr2 = (psiN_full - psiN_wish) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,ns
          if (dr2(j)<min_dr2) then
             index = j
             min_dr2 = dr2(j)
          end if
       end do

       psiN = psiN_full(index)
       deallocate(dr2)

    case default
       if (masterProc) then
          print *,"Error! Invalid VMECRadialOption"
       end if
       stop
    end select

    ! --------------------------------------------------------------------------------
    ! Done choosing the actual radius to use.
    ! --------------------------------------------------------------------------------

    ! In general, we get quantities for SFINCS by linear interpolation, taking a weighted average of the quantity from
    ! 2 surfaces in the VMEC file. Sometimes the weights are 0 and 1, i.e. no interpolation is needed.

    ! For any VMEC quantity Q on the full grid, the value used in SFINCS will be
    !  Q_sfincs = Q(vmecRadialIndex_full(1))*vmecRadialWeight_full(1) + Q(vmecRadialIndex_full(2))*vmecRadialWeight_full(2)

    ! For any VMEC quantity Q on the half grid, the value used in SFINCS will be
    !  Q_sfincs = Q(vmecRadialIndex_half(1))*vmecRadialWeight_half(1) + Q(vmecRadialIndex_half(2))*vmecRadialWeight_half(2)

    ! In VMEC, quantities on the half grid have the same number of array elements (ns) as quantities on the full grid,
    ! but the first array element is 0.



    ! Handle quantities for the full grid
    if (psiN>1) then
       if (masterProc) then
          print *,"Error! psiN cannot be >1"
       end if
       stop
    elseif (psiN<0) then
       if (masterProc) then
          print *,"Error! psiN cannot be <0"
       end if
       stop
    elseif (psiN==1) then
       vmecRadialIndex_full(1) = ns-1
       vmecRadialIndex_full(2) = ns
       vmecRadialWeight_full(1) = zero
    else
       ! psiN is >= 0 and <1
       ! This is the most common case.
       vmecRadialIndex_full(1) = floor(psiN*(ns-1))+1
       vmecRadialIndex_full(2) = vmecRadialIndex_full(1) + 1
       vmecRadialWeight_full(1) = vmecRadialIndex_full(1) - psiN*(ns-one)
    end if
    vmecRadialWeight_full(2) = one - vmecRadialWeight_full(1)
    
    ! Handle quantities for the half grid
    if (psiN < psiN_half(1)) then
       if (masterProc) then
          print *,"Warning: extrapolating beyond the end of VMEC's half grid."
          print *,"(Extrapolating towards the magnetic axis.)"
       end if
       ! We start at element 2 since element 1 is always 0 for quantities on the half grid.
       vmecRadialIndex_half(1) = 2
       vmecRadialIndex_half(2) = 3
       vmecRadialWeight_half(1) = (psiN_half(2) - psiN) / (psiN_half(2) - psiN_half(1))

    elseif (psiN > psiN_half(ns-1)) then
       if (masterProc) then
          print *,"Warning: extrapolating beyond the end of VMEC's half grid."
          print *,"(Extrapolating towards the last closed flux surface.)"
       end if
       vmecRadialIndex_half(1) = ns-1
       vmecRadialIndex_half(2) = ns
       vmecRadialWeight_half(1) = (psiN_half(ns-1) - psiN) &
            / (psiN_half(ns-1) - psiN_half(ns-2))

    elseif (psiN == psiN_half(ns-1)) then
       ! We are exactly at the last point of the half grid
       vmecRadialIndex_half(1) = ns-1
       vmecRadialIndex_half(2) = ns
       vmecRadialWeight_half(1) = zero
    else
       ! psiN is inside the half grid.
       ! This is the most common case.
       vmecRadialIndex_half(1) = floor(psiN*(ns-1) + 0.5d+0)+1
       if (vmecRadialIndex_half(1) < 2) then
          ! This can occur sometimes due to roundoff error.
          vmecRadialIndex_half(1) = 2
       end if
       vmecRadialIndex_half(2) = vmecRadialIndex_half(1) + 1
       vmecRadialWeight_half(1) = vmecRadialIndex_half(1) - psiN*(ns-one) - (0.5d+0)
    end if
    vmecRadialWeight_half(2) = one-vmecRadialWeight_half(1)

    if (masterProc) then
!!$       print *,"vmecRadialIndex_full:",vmecRadialIndex_full
!!$       print *,"vmecRadialWeight_full:",vmecRadialWeight_full
!!$       print *,"vmecRadialIndex_half:",vmecRadialIndex_half
!!$       print *,"vmecRadialWeight_half:",vmecRadialWeight_half
       if (abs(vmecRadialWeight_half(1)) < 1e-14) then
          print "(a,i3,a,i3,a)"," Using radial index ",vmecRadialIndex_half(2)," of ",ns," from vmec's half mesh."
       elseif (abs(vmecRadialWeight_half(2)) < 1e-14) then
          print "(a,i3,a,i3,a)"," Using radial index ",vmecRadialIndex_half(1)," of ",ns," from vmec's half mesh."
       else
          print "(a,i3,a,i3,a,i3,a)", " Interpolating using radial indices ",vmecRadialIndex_half(1)," and ",vmecRadialIndex_half(2),&
               " of ",ns," from vmec's half mesh."
          print "(a,f17.14,a,f17.14)", " Weights for half mesh = ",vmecRadialWeight_half(1)," and ",vmecRadialWeight_half(2)
          print "(a,i3,a,i3,a,i3,a)", " Interpolating using radial indices ",vmecRadialIndex_full(1)," and ",vmecRadialIndex_full(2),&
               " of ",ns," from vmec's full mesh."
          print "(a,f17.14,a,f17.14)", " Weights for full mesh = ",vmecRadialWeight_full(1)," and ",vmecRadialWeight_full(2)
       end if
    end if

    iota = iotas(vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
         + iotas(vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)

    dpsi = phi(2)/(2*pi)
    vmec_dpdpsiHat = 0
    vmec_dpdpsiHat(2:ns) = (presf(2:ns) - presf(1:ns-1)) / dpsi
    pPrimeHat = (4*pi*1.0d-7) * ( &
         vmec_dpdpsiHat(vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
         + vmec_dpdpsiHat(vmecRadialIndex_half(2)) * vmecRadialWeight_half(2))

    BHat = zero
    DHat = zero
    dBHatdtheta = zero
    dBHatdzeta = zero
    dBHatdpsiHat = zero

    BHat_sub_psi = zero
    dBHat_sub_psi_dtheta = zero
    dBHat_sub_psi_dzeta = zero

    BHat_sub_theta = zero
    dBHat_sub_theta_dpsiHat = zero
    dBHat_sub_theta_dzeta = zero

    BHat_sub_zeta = zero
    dBHat_sub_zeta_dpsiHat = zero
    dBHat_sub_zeta_dtheta = zero

    BHat_sup_theta = zero
    BHat_sup_zeta = zero

    ! First, get the (m=0,n=0) component of |B|, which will be used for testing whether
    ! other harmonics of |B| are large enough to include. Note that bmnc(1,:) represents the m=n=0 component.
    b00 = bmnc(1,vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
               + bmnc(1,vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)

    ! --------------------------------------------------------------------------------
    ! At last, we are now ready to
    ! loop over all the VMEC Fourier modes to build the SFINCS arrays.
    ! --------------------------------------------------------------------------------

    numSymmetricModesIncluded = 0
    numAntisymmetricModesIncluded = 0

    ! We take the following approach:
    ! Include a given (m,n) mode of |B|, B sub u, B sup u, etc
    ! if and only if that (m,n) mode of |B| (normalized to B00) is > min_Bmn_to_load.
    ! I.e., the size of a |B| harmonic controls not only whether that Fourier mode of |B|
    ! is included, but also whether that Fourier mode of the other fields like B sub u is included.
    do imn_nyq = 1, mnmax_nyq ! All the quantities we need except R and Z use the _nyq mode numbers.
       m = xm_nyq(imn_nyq)
       n = xn_nyq(imn_nyq)/nfp

       if (abs(m) >= mpol .or. abs(n) > ntor) then
          if (VMEC_Nyquist_option==1) then 
             !if (masterProc) print *,"Skipping m=",m,"n=",n," b/c Nyquist."
             cycle
          end if

          non_Nyquist_mode_available = .false.
       else
          non_Nyquist_mode_available = .true.
          ! Find the imn in the non-Nyquist arrays that corresponds to the same m and n.
          found_imn = .false.
          do imn = 1,mnmax
             if (xm(imn)==m .and. xn(imn)==n*nfp) then
                found_imn = .true.
                exit
             end if
          end do
          if ((xm(imn) .ne. m) .or. (xn(imn) .ne. n*nfp)) stop "Something went wrong!"
          if (.not. found_imn) stop "Error! imn could not be found matching the given imn_nyq."
       end if

       ! -----------------------------------------------------
       ! First, consider just the stellarator-symmetric terms:
       ! -----------------------------------------------------
       
       b = bmnc(imn_nyq,vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
            + bmnc(imn_nyq,vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)
       
       ! Set scaleFactor to rippleScale for non-axisymmetric or non-quasisymmetric modes
       scaleFactor = setScaleFactor(n,m)
       b = b*scaleFactor
	
       if (abs(b/b00) >= min_Bmn_to_load) then
          ! This (m,n) mode is sufficiently large to include.
          !if (masterProc) then
          !   print *,"Including mode with m = ",m,", n = ",n
          !end if
          numSymmetricModesIncluded = numSymmetricModesIncluded + 1

          ! Evaluate the radial derivatives we will need:
          dpsi = phi(2)/(2*pi)  ! Doesn't need to be in the loops, but here for convenience.

          ! B, B_sub_theta, and B_sub_zeta are on the half mesh, so their radial derivatives are on the full mesh.
          ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.

          vmec_dBHatdpsiHat(2:ns-1) = (bmnc(imn_nyq,3:ns) - bmnc(imn_nyq,2:ns-1)) / dpsi
          ! Simplistic "extrapolation" at the endpoints:
          vmec_dBHatdpsiHat(1) = vmec_dBHatdpsiHat(2)
          vmec_dBHatdpsiHat(ns) = vmec_dBHatdpsiHat(ns-1)

          vmec_dBHat_sub_theta_dpsiHat(2:ns-1) = (bsubumnc(imn_nyq,3:ns) - bsubumnc(imn_nyq,2:ns-1)) / dpsi
          vmec_dBHat_sub_theta_dpsiHat(1) = vmec_dBHat_sub_theta_dpsiHat(2)
          vmec_dBHat_sub_theta_dpsiHat(ns) = vmec_dBHat_sub_theta_dpsiHat(ns-1)

          vmec_dBHat_sub_zeta_dpsiHat(2:ns-1) = (bsubvmnc(imn_nyq,3:ns) - bsubvmnc(imn_nyq,2:ns-1)) / dpsi
          vmec_dBHat_sub_zeta_dpsiHat(1) = vmec_dBHat_sub_zeta_dpsiHat(2)
          vmec_dBHat_sub_zeta_dpsiHat(ns) = vmec_dBHat_sub_zeta_dpsiHat(ns-1)

          if (non_Nyquist_mode_available) then
             vmec_dRdpsiHat(2:ns) = (rmnc(imn,2:ns) - rmnc(imn,1:ns-1)) / dpsi
             vmec_dRdpsiHat(1) = 0

             vmec_dZdpsiHat(2:ns) = (zmns(imn,2:ns) - zmns(imn,1:ns-1)) / dpsi
             vmec_dZdpsiHat(1) = 0
          else
             vmec_dRdpsiHat = 0
             vmec_dZdpsiHat = 0
          end if

          ! End of evaluating radial derivatives.
             
          do itheta = 1,Ntheta
             do izeta = 1,Nzeta
                angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                cos_angle = cos(angle)
                sin_angle = sin(angle)
                
                BHat(itheta,izeta) = BHat(itheta,izeta) + b * cos_angle
                
                dbHatdtheta(itheta,izeta) = dBHatdtheta(itheta,izeta) - m * b * sin_angle
                
                dbHatdzeta(itheta,izeta) = dBHatdzeta(itheta,izeta) + n * NPeriods * b * sin_angle
                
                do isurf = 1,2
                   ! Handle Jacobian:
                   ! SFINCS's DHat is the INVERSE Jacobian of the (psiHat, theta, zeta) coordinates.
                   ! VMEC's gmnc and gmns are the Jacobian of the (psiN, theta, zeta) coordinates.
                   ! Because one uses psiHat and the other uses psiN, we need a factor of psiAHat for conversion.
                   ! We will also set DHat = 1 / DHat at the end of this loop.
                   temp = gmnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf) / (psiAHat)
                   temp = temp*scaleFactor
                   DHat(itheta,izeta) = DHat(itheta,izeta) + temp * cos_angle

                   ! Handle B sup theta:
                   ! Note that VMEC's bsupumnc and bsupumns are exactly the same as SFINCS's BHat_sup_theta, with no conversion factors of 2pi needed.
                   temp = bsupumnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                   temp = temp*scaleFactor
                   BHat_sup_theta(itheta,izeta) = BHat_sup_theta(itheta,izeta) + temp * cos_angle
                   dBHat_sup_theta_dzeta(itheta,izeta) = dBHat_sup_theta_dzeta(itheta,izeta) + n * NPeriods * temp * sin_angle

                   ! Handle B sup zeta:
                   ! Note that VMEC's bsupvmnc and bsupvmns are exactly the same as SFINCS's BHat_sup_zeta, with no conversion factors of 2pi or Nperiods needed.
                   temp = bsupvmnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                   temp = temp*scaleFactor
                   BHat_sup_zeta(itheta,izeta) = BHat_sup_zeta(itheta,izeta) + temp * cos_angle
                   dBHat_sup_zeta_dtheta(itheta,izeta) = dBHat_sup_zeta_dtheta(itheta,izeta) - m * temp * sin_angle

                   ! Handle B sub theta:
                   ! Note that VMEC's bsubumnc and bsubumns are exactly the same as SFINCS's BHat_sub_theta, with no conversion factors of 2pi needed.
                   temp = bsubumnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                   temp = temp*scaleFactor
                   BHat_sub_theta(itheta,izeta) = BHat_sub_theta(itheta,izeta) + temp * cos_angle
                   dBHat_sub_theta_dzeta(itheta,izeta) = dBHat_sub_theta_dzeta(itheta,izeta) + n * NPeriods * temp * sin_angle

                   ! Handle B sub zeta:
                   ! Note that VMEC's bsubvmnc and bsubvmns are exactly the same as SFINCS's BHat_sub_zeta, with no conversion factors of 2pi needed.
                   temp = bsubvmnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                   temp = temp*scaleFactor
                   BHat_sub_zeta(itheta,izeta) = BHat_sub_zeta(itheta,izeta) + temp * cos_angle
                   dBHat_sub_zeta_dtheta(itheta,izeta) = dBHat_sub_zeta_dtheta(itheta,izeta) - m * temp * sin_angle

                   ! Handle B sub psi.
                   ! Unlike the other components of B, this one is on the full mesh.
                   ! Notice B_psi = B_s * (d s / d psi), and (d s / d psi) = 1 / psiAHat
                   temp = bsubsmns(imn_nyq,vmecRadialIndex_full(isurf)) / psiAHat * vmecRadialWeight_full(isurf)
                   temp = temp*scaleFactor
                   BHat_sub_psi(itheta,izeta) = BHat_sub_psi(itheta,izeta) + temp * sin_angle
                   dBHat_sub_psi_dtheta(itheta,izeta) = dBHat_sub_psi_dtheta(itheta,izeta) + m * temp * cos_angle
                   dBHat_sub_psi_dzeta(itheta,izeta)  = dBHat_sub_psi_dzeta(itheta,izeta) - n * NPeriods * temp * cos_angle

                   ! Handle dBHatdpsiHat.
                   ! Since bmnc is on the half mesh, its radial derivative is on the full mesh.
                   temp = vmec_dBHatdpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                   temp = temp*scaleFactor
                   dBHatdpsiHat(itheta,izeta) = dBHatdpsiHat(itheta,izeta) + temp * cos_angle

                   ! Handle dBHat_sub_theta_dpsiHat.
                   ! Since bsubumnc is on the half mesh, its radial derivative is on the full mesh.
                   temp = vmec_dBHat_sub_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                   temp = temp*scaleFactor
                   dBHat_sub_theta_dpsiHat(itheta,izeta) = dBHat_sub_theta_dpsiHat(itheta,izeta) + temp * cos_angle

                   ! Handle dBHat_sub_zeta_dpsiHat.
                   ! Since bsubvmnc is on the half mesh, its radial derivative is on the full mesh.
                   temp = vmec_dBHat_sub_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                   temp = temp*scaleFactor
                   dBHat_sub_zeta_dpsiHat(itheta,izeta) = dBHat_sub_zeta_dpsiHat(itheta,izeta) + temp * cos_angle

                   ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                   if (non_Nyquist_mode_available) then

                      ! Handle R, which is on the full mesh
                      temp = rmnc(imn,vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scaleFactor
                      R(itheta,izeta) = R(itheta,izeta) + temp * cos_angle
                      dRdtheta(itheta,izeta) = dRdtheta(itheta,izeta) - temp * m * sin_angle
                      dRdzeta(itheta,izeta)  = dRdzeta(itheta,izeta)  + temp * n * Nperiods * sin_angle

                      ! Handle Z, which is on the full mesh
                      temp = zmns(imn,vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scaleFactor
                      !Z(itheta,izeta) = Z(itheta,izeta) + temp * sin_angle  ! We don't actually need Z itself, only derivatives of Z.
                      dZdtheta(itheta,izeta) = dZdtheta(itheta,izeta) + temp * m * cos_angle
                      dZdzeta(itheta,izeta)  = dZdzeta(itheta,izeta)  - temp * n * Nperiods * cos_angle

                      ! Handle dRdpsiHat.
                      ! Since R is on the full mesh, its radial derivative is on the half mesh.
                      temp = vmec_dRdpsiHat(vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scaleFactor
                      dRdpsiHat(itheta,izeta) = dRdpsiHat(itheta,izeta) + temp * cos_angle

                      ! Handle dZdpsiHat.
                      ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                      temp = vmec_dZdpsiHat(vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scaleFactor
                      dZdpsiHat(itheta,izeta) = dZdpsiHat(itheta,izeta) + temp * sin_angle

                   end if
                end do
             end do
          end do
       else
          !if (masterProc) then
          !   print *,"NOT including mode with m = ",m,", n = ",n
          !end if
       end if

       ! -----------------------------------------------------
       ! Now consider the stellarator-asymmetric terms.
       ! NOTE: This functionality has not been tested as thoroughly !!!
       ! -----------------------------------------------------

       if (lasym) then

          b = bmns(imn_nyq,vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
               + bmns(imn_nyq,vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)

          ! Set scaleFactor to rippleScale for non-axisymmetric or non-quasisymmetric modes
          scaleFactor = setScaleFactor(n,m)
          b = b*scaleFactor
          
          if (abs(b/b00) >= min_Bmn_to_load) then
             ! This (m,n) mode is sufficiently large to include.
             !if (masterProc) then
             !   print *,"Including stellarator-asymmetric mode with m = ",m,", n = ",n
             !end if
             numAntisymmetricModesIncluded = numAntisymmetricModesIncluded + 1
             
             ! Evaluate the radial derivatives we will need:
             dpsi = phi(2)/(2*pi)  ! Doesn't need to be in the loops, but here for convenience.
             
             ! B, B_sub_theta, and B_sub_zeta are on the half mesh, so their radial derivatives are on the full mesh.
             ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.
             
             vmec_dBHatdpsiHat(2:ns-1) = (bmns(imn_nyq,3:ns) - bmns(imn_nyq,2:ns-1)) / dpsi
             ! Simplistic "extrapolation" at the endpoints:
             vmec_dBHatdpsiHat(1) = vmec_dBHatdpsiHat(2)
             vmec_dBHatdpsiHat(ns) = vmec_dBHatdpsiHat(ns-1)
             
             vmec_dBHat_sub_theta_dpsiHat(2:ns-1) = (bsubumns(imn_nyq,3:ns) - bsubumns(imn_nyq,2:ns-1)) / dpsi
             vmec_dBHat_sub_theta_dpsiHat(1) = vmec_dBHat_sub_theta_dpsiHat(2)
             vmec_dBHat_sub_theta_dpsiHat(ns) = vmec_dBHat_sub_theta_dpsiHat(ns-1)
             
             vmec_dBHat_sub_zeta_dpsiHat(2:ns-1) = (bsubvmns(imn_nyq,3:ns) - bsubvmns(imn_nyq,2:ns-1)) / dpsi
             vmec_dBHat_sub_zeta_dpsiHat(1) = vmec_dBHat_sub_zeta_dpsiHat(2)
             vmec_dBHat_sub_zeta_dpsiHat(ns) = vmec_dBHat_sub_zeta_dpsiHat(ns-1)
             
             if (non_Nyquist_mode_available) then
                vmec_dRdpsiHat(2:ns) = (rmns(imn,2:ns) - rmns(imn,1:ns-1)) / dpsi
                vmec_dRdpsiHat(1) = 0
             
                vmec_dZdpsiHat(2:ns) = (zmnc(imn,2:ns) - zmnc(imn,1:ns-1)) / dpsi
                vmec_dZdpsiHat(1) = 0
             else
                vmec_dRdpsiHat = 0
                vmec_dZdpsiHat = 0
             end if
             
             ! End of evaluating radial derivatives.
             
             do itheta = 1,Ntheta
                do izeta = 1,Nzeta
                   angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                   cos_angle = cos(angle)
                   sin_angle = sin(angle)
                      
                   BHat(itheta,izeta) = BHat(itheta,izeta) + b * sin_angle
                   dbHatdtheta(itheta,izeta) = dBHatdtheta(itheta,izeta) + m * b * cos_angle
                   dbHatdzeta(itheta,izeta) = dBHatdzeta(itheta,izeta) - n * NPeriods * b * cos_angle
                   
                   do isurf = 1,2
                      ! Handle Jacobian:
                      ! SFINCS's DHat is the INVERSE Jacobian of the (psiHat, theta, zeta) coordinates.
                      ! VMEC's gmnc and gmns are the Jacobian of the (psiN, theta, zeta) coordinates.
                      ! Because one uses psiHat and the other uses psiN, we need a factor of psiAHat for conversion.
                      ! We will also set DHat = 1 / DHat at the end of this loop.
                      temp = gmns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf) / (psiAHat)
                      temp = temp*scaleFactor
                      DHat(itheta,izeta) = DHat(itheta,izeta) + temp * sin_angle
                         
                      ! Handle B sup theta:
                      ! Note that VMEC's bsupumnc and bsupumns are exactly the same as SFINCS's BHat_sup_theta, with no conversion factors of 2pi needed.
                      temp = bsupumns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scaleFactor	
                      BHat_sup_theta(itheta,izeta) = BHat_sup_theta(itheta,izeta) + temp * sin_angle
                      dBHat_sup_theta_dzeta(itheta,izeta) = dBHat_sup_theta_dzeta(itheta,izeta) - n * NPeriods * temp * cos_angle
                      
                      ! Handle B sup zeta:
                      ! Note that VMEC's bsupvmnc and bsupvmns are exactly the same as SFINCS's BHat_sup_zeta, with no conversion factors of 2pi or Nperiods needed.
                      temp = bsupvmns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scaleFactor
                      BHat_sup_zeta(itheta,izeta) = BHat_sup_zeta(itheta,izeta) + temp * sin_angle
                      dBHat_sup_zeta_dtheta(itheta,izeta) = dBHat_sup_zeta_dtheta(itheta,izeta) + m * temp * cos_angle
                      
                      ! Handle B sub theta:
                      ! Note that VMEC's bsubumnc and bsubumns are exactly the same as SFINCS's BHat_sub_theta, with no conversion factors of 2pi needed.
                      temp = bsubumns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scaleFactor
                      BHat_sub_theta(itheta,izeta) = BHat_sub_theta(itheta,izeta) + temp * sin_angle
                      dBHat_sub_theta_dzeta(itheta,izeta) = dBHat_sub_theta_dzeta(itheta,izeta) - n * NPeriods * temp * cos_angle
                      
                      ! Handle B sub zeta:
                      ! Note that VMEC's bsubvmnc and bsubvmns are exactly the same as SFINCS's BHat_sub_zeta, with no conversion factors of 2pi needed.
                      temp = bsubvmns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scaleFactor
                      BHat_sub_zeta(itheta,izeta) = BHat_sub_zeta(itheta,izeta) + temp * sin_angle
                      dBHat_sub_zeta_dtheta(itheta,izeta) = dBHat_sub_zeta_dtheta(itheta,izeta) + m * temp * cos_angle
                      
                      ! Handle B sub psi.
                      ! Unlike the other components of B, this one is on the full mesh.
                      ! Notice B_psi = B_s * (d s / d psi), and (d s / d psi) = 1 / psiAHat
                      temp = bsubsmnc(imn_nyq,vmecRadialIndex_full(isurf)) / psiAHat * vmecRadialWeight_full(isurf)
                      temp = temp*scaleFactor
                      BHat_sub_psi(itheta,izeta) = BHat_sub_psi(itheta,izeta) + temp * cos_angle
                      dBHat_sub_psi_dtheta(itheta,izeta) = dBHat_sub_psi_dtheta(itheta,izeta) - m * temp * sin_angle
                      dBHat_sub_psi_dzeta(itheta,izeta)  = dBHat_sub_psi_dzeta(itheta,izeta) + n * NPeriods * temp * sin_angle
                      
                      ! Handle dBHatdpsiHat.
                      ! Since bmns is on the half mesh, its radial derivative is on the full mesh.
                      temp = vmec_dBHatdpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scaleFactor
                      dBHatdpsiHat(itheta,izeta) = dBHatdpsiHat(itheta,izeta) + temp * sin_angle
                      
                      ! Handle dBHat_sub_theta_dpsiHat.
                      ! Since bsubumns is on the half mesh, its radial derivative is on the full mesh.
                      temp = vmec_dBHat_sub_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scaleFactor
                      dBHat_sub_theta_dpsiHat(itheta,izeta) = dBHat_sub_theta_dpsiHat(itheta,izeta) + temp * sin_angle
                      
                      ! Handle dBHat_sub_zeta_dpsiHat.
                      ! Since bsubvmns is on the half mesh, its radial derivative is on the full mesh.
                      temp = vmec_dBHat_sub_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scaleFactor
                      dBHat_sub_zeta_dpsiHat(itheta,izeta) = dBHat_sub_zeta_dpsiHat(itheta,izeta) + temp * sin_angle
                      
                      ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                      if (non_Nyquist_mode_available) then

                         ! Handle R, which is on the full mesh
                         temp = rmns(imn,vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         temp = temp*scaleFactor
                         R(itheta,izeta) = R(itheta,izeta) + temp * sin_angle
                         dRdtheta(itheta,izeta) = dRdtheta(itheta,izeta) + temp * m * cos_angle
                         dRdzeta(itheta,izeta)  = dRdzeta(itheta,izeta)  - temp * n * Nperiods * cos_angle

                         ! Handle Z, which is on the full mesh
                         temp = zmnc(imn,vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         temp = temp*scaleFactor
                         ! Z(itheta,izeta) = Z(itheta,izeta) + temp * cos_angle   ! We don't actually need Z itself, only derivatives of Z.
                         dZdtheta(itheta,izeta) = dZdtheta(itheta,izeta) - temp * m * sin_angle
                         dZdzeta(itheta,izeta)  = dZdzeta(itheta,izeta)  + temp * n * Nperiods * sin_angle

                         ! Handle dRdpsiHat.
                         ! Since R is on the full mesh, its radial derivative is on the half mesh.
                         temp = vmec_dRdpsiHat(vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         temp = temp*scaleFactor
                         dRdpsiHat(itheta,izeta) = dRdpsiHat(itheta,izeta) + temp * sin_angle
                         
                         ! Handle dZdpsiHat.
                         ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                         temp = vmec_dZdpsiHat(vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         temp = temp*scaleFactor
                         dZdpsiHat(itheta,izeta) = dZdpsiHat(itheta,izeta) + temp * cos_angle
                      end if
                   end do
                end do
             end do
          else
             !if (masterProc) then
             !   print *,"NOT including mode with m = ",m,", n = ",n
             !end if
          end if
       end if
    end do

    if (masterProc) then
       print "(a,i4,a,i4,a)"," Including ",numSymmetricModesIncluded," of ",(2*ntor+1)*mpol," stellarator-symmetric modes from the VMEC file."
       if (lasym) then
          print "(a,i4,a,i4,a)"," Including ",numAntisymmetricModesIncluded," of ",(2*ntor+1)*mpol," stellarator-antisymmetric modes from the VMEC file."
       else
          print *,"Equilibrium is stellarator-symmetric."
       end if
    end if

    ! Convert Jacobian to inverse Jacobian:
    DHat = one / DHat

    do izeta = 1,Nzeta
       cos_angle = cos(zeta(izeta))
       sin_angle = sin(zeta(izeta))

       ! X = R * cos(zeta)
       dXdtheta(:,izeta) = dRdtheta(:,izeta) * cos_angle
       dXdzeta(:,izeta) = dRdzeta(:,izeta) * cos_angle - R(:,izeta) * sin_angle
       dXdpsiHat(:,izeta) = dRdpsiHat(:,izeta) * cos_angle

       ! Y = R * sin(zeta)
       dYdtheta(:,izeta) = dRdtheta(:,izeta) * sin_angle
       dYdzeta(:,izeta) = dRdzeta(:,izeta) * sin_angle + R(:,izeta) * cos_angle
       dYdpsiHat(:,izeta) = dRdpsiHat(:,izeta) * sin_angle
    end do

    g_sub_theta_theta = dXdtheta*dXdtheta  + dYdtheta*dYdtheta  + dZdtheta*dZdtheta
    g_sub_theta_zeta  = dXdtheta*dXdzeta   + dYdtheta*dYdzeta   + dZdtheta*dZdzeta
    g_sub_zeta_zeta   = dXdzeta *dXdzeta   + dYdzeta *dYdzeta   + dZdzeta *dZdzeta
    g_sub_psi_theta   = dXdpsiHat*dXdtheta + dYdpsiHat*dYdtheta + dZdpsiHat*dZdtheta
    g_sub_psi_zeta    = dXdpsiHat*dXdzeta  + dYdpsiHat*dYdzeta  + dZdpsiHat*dZdzeta
    g_sub_psi_psi     = dXdpsiHat*dXdpsiHat+ dYdpsiHat*dYdpsiHat+ dZdpsiHat*dZdpsiHat

    gradpsidotgradB_overgpsipsi = dBHatdpsiHat &
         + ((g_sub_theta_zeta*g_sub_psi_zeta-g_sub_psi_theta*g_sub_zeta_zeta)*dBHatdtheta &
         + (g_sub_psi_theta*g_sub_theta_zeta-g_sub_theta_theta*g_sub_psi_zeta)*dBHatdzeta) &
         / (g_sub_theta_theta*g_sub_zeta_zeta - g_sub_theta_zeta*g_sub_theta_zeta)

    gpsipsi = 1.0/(g_sub_psi_psi +&
                     (g_sub_psi_theta*(g_sub_theta_zeta*g_sub_psi_zeta-g_sub_psi_theta*g_sub_zeta_zeta) &
                      +g_sub_psi_zeta*(g_sub_psi_theta*g_sub_theta_zeta-g_sub_theta_theta*g_sub_psi_zeta)) &
                     /(g_sub_theta_theta*g_sub_zeta_zeta - g_sub_theta_zeta*g_sub_theta_zeta))
    
    ! These next lines should be replaced eventually with a proper calculation:
    dBHat_sup_theta_dpsiHat = 0
    dBHat_sup_zeta_dpsiHat = 0

    deallocate(psiN_full)
    deallocate(psiN_half)
    deallocate(vmec_dBHatdpsiHat)
    deallocate(vmec_dBHat_sub_theta_dpsiHat)
    deallocate(vmec_dBHat_sub_zeta_dpsiHat)
    deallocate(vmec_dRdpsiHat)
    deallocate(vmec_dZdpsiHat)
    deallocate(vmec_dpdpsiHat)
    deallocate(R,dRdtheta,dRdzeta,dRdpsiHat)
    deallocate(dXdtheta,dXdzeta,dXdpsiHat)
    deallocate(dYdtheta,dYdzeta,dYdpsiHat)
    deallocate(dZdtheta,dZdzeta,dZdpsiHat)
    deallocate(g_sub_theta_theta,g_sub_theta_zeta,g_sub_zeta_zeta,g_sub_psi_theta,g_sub_psi_zeta,g_sub_psi_psi)


    rN = sqrt(psiN)

    if (masterProc) then
       print *,"Successfully read VMEC geometry file ",trim(equilibriumFile)
    end if

  end subroutine computeBHat_VMEC

  ! ------------------------------------------------------------------------------------------

  subroutine computeBIntegrals

    implicit none

    integer :: itheta, izeta

    ! This subroutine computes VPrimeHat, FSABHat2, and (if needed) B0OverBBar, GHat, and IHat.

    VPrimeHat = 0
    FSABHat2 = 0
    FSABHat = 0
    FSAInvB2_minusInvFSAB2 = 0
    do itheta=1,Ntheta
       do izeta=1,Nzeta
          VPrimeHat = VPrimeHat + thetaWeights(itheta) * zetaWeights(izeta) / DHat(itheta,izeta)
          FSABHat2 = FSABHat2 + thetaWeights(itheta) * zetaWeights(izeta) &
               * BHat(itheta,izeta) * BHat(itheta,izeta) / DHat(itheta,izeta)
          FSABHat  = FSABHat + thetaWeights(itheta) * zetaWeights(izeta) &
               * BHat(itheta,izeta) / DHat(itheta,izeta)
       end do
    end do

    FSABHat2 = FSABHat2 / VPrimeHat
    FSABHat  = FSABHat  / VPrimeHat

    ! MFM 03/25/19
    do itheta=1,Ntheta
       do izeta=1,Nzeta
          FSAInvB2_minusInvFSAB2 = FSAInvB2_minusInvFSAB2 + thetaWeights(itheta) * zetaWeights(izeta) &
               * ((1/(BHat(itheta,izeta)*BHat(itheta,izeta))) - (1/FSABHat2)) / DHat(itheta,izeta)
       end do
    end do
          
    FSAInvB2_minusInvFSAB2 = FSAInvB2_minusInvFSAB2 / VPrimeHat

    if (coordinateSystem .ne. COORDINATE_SYSTEM_BOOZER) then
       ! Compute B0, the (m=0,n=0) Boozer harmonic.
       ! We compute it using B0 = <B^3> / <B^2>, where the right hand side
       ! can be computed in any coordinate system.
       B0OverBBar = 0
       do itheta=1,Ntheta
          do izeta=1,Nzeta
             B0OverBBar = B0OverBBar + thetaWeights(itheta) * zetaWeights(izeta) &
                  * (BHat(itheta,izeta) ** 3) / DHat(itheta,izeta)
          end do
       end do

       B0OverBBar = B0OverBBar / VPrimeHat / FSABHat2

       ! Compute G and H for the case of non-Boozer coordinates.
       ! We don't actually need G or H for anything, but they are not much work to compute,
       ! so let's compute them in case it turns out to be convenient.
       !
       ! (We could do this computation on a finer (theta,zeta) grid than the grid
       ! used for the kinetic calculation. But let's do this on the same grid for 2 reasons:
       ! (1) simplicity, and (2) since the trapezoid rule is spectrally
       ! accurate on a uniform periodic grid, so a fine grid is not required.)
       GHat = dot_product(thetaWeights, matmul(BHat_sub_zeta,  zetaWeights)) / (4*pi*pi)
       IHat = dot_product(thetaWeights, matmul(BHat_sub_theta, zetaWeights)) / (4*pi*pi)

    end if

    if (masterProc) then
       print *,"---- Geometry parameters: ----"
       print *,"Geometry scheme = ", geometryScheme
       print *,"psiAHat (Normalized toroidal flux at the last closed flux surface) = ", psiAHat
       print *,"aHat (Radius of the last closed flux surface in units of RHat) = ", aHat
       if (geometryScheme==1) then
          print *,"epsilon_t = ", epsilon_t
          print *,"epsilon_h = ", epsilon_h
          print *,"epsilon_antisymm = ", epsilon_antisymm
       end if
       print *,"GHat (Boozer component multiplying grad zeta) = ", GHat
       print *,"IHat (Boozer component multiplying grad theta) = ", IHat
       print *,"iota (Rotational transform) = ", iota
    end if

  end subroutine computeBIntegrals

! Set scale to value of rippleScale for non-axisymmetric or non-quasisymmetric components
  function setScaleFactor(n,m) result(scale)
    integer :: n
    integer :: m
    PetscScalar :: scale 
    scale = 1
    if (helicity_n == 0) then
      if (n /=0 .or. m/=helicity_l) then
        scale = rippleScale
      end if
    else
      if ((n /= 0) .and. (helicity_l/helicity_n) /= (m/n)) then
        scale = rippleScale
      else if (n == 0) then
        scale = rippleScale
      end if
    end if
  end function setScaleFactor

  ! ---------------------------------------------------------------------------------------
          
  subroutine load_bcdat_file

    implicit none

    integer :: itheta, izeta, NHarmonics, NHarmonicsL, NHarmonicsH, i, m, n
    integer, dimension(:), allocatable :: BHarmonics_l, BHarmonics_n
    integer, dimension(:), allocatable :: BHarmonics_lL, BHarmonics_nL, BHarmonics_lH, BHarmonics_nH
    PetscScalar, dimension(:), allocatable :: BHarmonics_amplitudes
    PetscScalar, dimension(:), allocatable :: BHarmonics_amplitudesL, BHarmonics_amplitudesH
    logical, dimension(:), allocatable :: BHarmonics_parity
    logical, dimension(:), allocatable :: BHarmonics_parityL, BHarmonics_parityH
    PetscScalar, dimension(:,:), allocatable :: bcdataL
    PetscScalar, dimension(:,:), allocatable :: bcdataH

    character(len=200) :: bcdatFile
    integer :: fileUnit, didFileAccessWork
    character(len=200) :: lineOfFile
    integer, dimension(4) :: headerIntegers
    PetscScalar, dimension(3) :: headerReals
    integer :: bcStelSym, dataStelSym
    PetscScalar, dimension(6) :: surfHeader
    PetscScalar :: dataNumbers
    PetscScalar, dimension(2) :: data2Numbers
    integer, dimension(2) :: dataIntegers
    integer :: no_of_modes_old, no_of_modes_new, modeind, numB0s, startn, stopn
    PetscScalar :: iota_old, iota_new, GHat_old, GHat_new, IHat_old, IHat_new
    PetscScalar :: pPrimeHat_old, pPrimeHat_new, invFSA_BHat2
    logical :: end_of_file, proceed, include_mn, nearbyRadiiGiven, nonStelSym
    integer, parameter :: max_no_of_modes = 10000
    integer, dimension(max_no_of_modes) :: modesm_old, modesm_new, modesn_old, modesn_new
    PetscScalar, dimension(max_no_of_modes) :: modesb_old, modesb_new
    PetscScalar :: rN_old,  rN_new, B0_old, B0_new, bcdata00L, bcdata00H
    PetscScalar :: DeltapsiHat !, diotadpsiHat moved to global variables 2016-09-15 HS
    PetscScalar :: RadialWeight = 1.0 ! weight of closest surface with rN<=rN_wish
    integer :: bcdat_NPeriods
    PetscScalar :: bcdat_psiAHat, bcdat_aHat, bcdat_iota 
    PetscScalar :: bcdat_GHat, bcdat_IHat, bcdata00
    ! For the BHarmonics_parity array, 
    ! true indicates the contribution to B(theta,zeta) has the form
    ! cos(l * theta - n * zeta)
    ! while false indicates the contribution to B(theta,zeta) has the form
    ! sin(l * theta - n * zeta)

    !nearbyRadiiGiven = .false. !Will be the case for all geometrySchemes except 11 and 12

    bcdatFile=EParallelHatSpec_bcdatFile !I am only using load_bcdat_file for this purpose, but it could be generalised.
    
    fileUnit = 11
    open(unit=fileUnit, file=bcdatFile, action="read", status="old", iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
       print *,"Unable to open the bcdat file ",bcdatFile
       stop
    end if
    do
       read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
       !print *,lineOfFile
       if (lineOfFile(1:2) /= "CC") exit
    end do

    ! Read the text line " Stellarator symmetry of bc and bcdat ... "
    ! read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
    !print *, lineOfFile
    ! Read stellarator symmetry flags:
    read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) bcStelSym,dataStelSym
    if (didFileAccessWork /= 0) then
       print *,"Unable to read StelSym header from the bcdat file ",bcdatFile
       stop
    end if

    ! Read the text line " m0b  n0b nsurf nper flux/[Tm^2]  ... "
    read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
    ! Read header line:
    read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) headerIntegers, headerReals
    !print *,'headerIntegers, headerReals='
    !print *,headerIntegers
    !print *,headerReals
    if (didFileAccessWork /= 0) then
       print *,"Unable to read header from the bcdat file ",bcdatFile
       stop
    end if
    bcdat_NPeriods = headerIntegers(4)
    bcdat_psiAHat  = headerReals(1)/2/pi !Convert the flux from Tm^2 to Tm^2/rad
    bcdat_aHat     = headerReals(2)      !minor radius in meters

    if ((geometryScheme /= 11) .and. (geometryScheme /= 12)) then
       print *,"Loading bcdat files is only compatible with geometryScheme 11 and 12!"
       stop
    end if
    if (bcdat_NPeriods /= NPeriods) then !.or. (bcdat_psiAHat /= psiAHat) .or. (bcdat_aHat /= aHat)) then
       print *,"The bcdat file header is inconsistent with the bc file header!"
       print *,"bcdat: NPeriods = ", bcdat_NPeriods
       print *,"bc   : NPeriods = ", NPeriods
       print *,"bcdat: psiAHat = ", bcdat_psiAHat
       print *,"bc   : psiAHat = ", psiAHat
       print *,"bcdat: aHat    = ", bcdat_NPeriods
       print *,"bc   : aHat    = ", aHat 
       stop
    end if
    
    end_of_file = .false.

    rN_old = 0
    no_of_modes_old = 0
    modesm_old = 0
    modesn_old = 0
    modesb_old = 0
    iota_old = 0
    GHat_old = 0
    IHat_old = 0
    B0_old = 0
    pPrimeHat_old = 0

    rN_new = 0
    no_of_modes_new = 0
    modesm_new = 0
    modesn_new = 0
    modesb_new = 0
    iota_new = 0
    GHat_new = 0
    IHat_new = 0
    B0_new = 0
    pPrimeHat_new = 0
    
    select case (dataStelSym)
    case (1)! Read a bcdat file with stellarator symmetric data

       ! Skip a line
       read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
       !print *,'s line='
       !print *,lineOfFile

       do 
          if ((rN_new .ge. rN_wish) .or. end_of_file) exit

          rN_old = rN_new
          no_of_modes_old = no_of_modes_new
          modesm_old = modesm_new
          modesn_old = modesn_new
          modesb_old = modesb_new
          iota_old = iota_new
          GHat_old = GHat_new
          IHat_old = IHat_new
          B0_old = B0_new
          pPrimeHat_old = pPrimeHat_new
          numB0s = 0

          if (.not.(index(lineOfFile,"[A]") > 0)) then !units were not on this line, so skip next
             ! Skip a line:
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             !print *,'skip unit [A] line='
             !print *,lineOfFile
          end if

          ! Read the header for the magnetic surface:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) surfHeader
          !print *,'surfHeader='
          !print *,surfHeader

          rN_new = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
          iota_new = surfHeader(2)
          ! Note that G and I have a minus sign in the following two lines
          ! because Ampere's law comes with a minus sign in the left-handed
          ! (r,pol,tor) system.
          GHat_new = -surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
          IHat_new = -surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
          pPrimeHat_new = surfheader(5)/psiAHat*(4*pi*1e-7)   ! dpdpsi=pPrimeHat/mu_0

          ! Skip units line:
          read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
          proceed = .true.
          modeind = 0
          !print *,'s line='
          !print *,lineOfFile


          do
             if (.not. proceed) exit

             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             if (didFileAccessWork /= 0) then
                proceed = .false.
                end_of_file = .true.
             else if (index(lineOfFile,"s") > 0) then
                ! Next flux surface has been reached
                proceed = .false.
                !print *,'s line='
                !print *,lineOfFile
             else
                read(unit=lineOfFile, fmt=*) dataIntegers, dataNumbers
                !print *,'dataIntegers='
                !print *,dataIntegers
                !print *,'dataNumbers='
                !print *,dataNumbers
                if (dataIntegers(1) == 0 .and. dataIntegers(2) == 0) then
                   B0_new = dataNumbers
                   numB0s = numB0s + 1
                else !!! if (abs(dataNumbers) > min_bcdata_to_load) then
                   modeind = modeind + 1
                   if (modeind > max_no_of_modes) then
                      print *,"The value of max_no_of_modes in geometry.F90 was insufficient."
                      print *,"Increase this value and recompile."
                      stop
                   end if
                   modesm_new(modeind) = dataIntegers(1)
                   modesn_new(modeind) = dataIntegers(2)
                   modesb_new(modeind) = dataNumbers
                end if
             end if
          end do
          if (numB0s == 0) then
             print *,"Warning: no (0,0) mode found in bcdat file ",bcdatFile
          else if (numB0s > 1) then
             print *,"Error: more than 1 (0,0) mode found in bcdat file ",bcdatFile
          end if
          no_of_modes_new = modeind
       end do

       
       close(unit = fileUnit)
       if (masterProc) then
          print *,"Successfully read data from bcdat file ",trim(bcdatFile)
       end if

       DeltapsiHat = psiAHat * (rN_new*rN_new-rN_old*rN_old)
       nearbyRadiiGiven = .true.


       if (VMECRadialOption == 1) then !Choose the nearest flux surface available
          if (abs(rN_old - rN_wish) < abs(rN_new - rN_wish)) then
             RadialWeight = 1.0
             rN = rN_old
          else
             RadialWeight = 0.0
             rN = rN_new
          end if
       else !Linear interpolation in s=rN^2
          RadialWeight = (rN_new*rN_new-rN_wish*rN_wish) / (rN_new*rN_new-rN_old*rN_old)
          rN   = rN_wish
       end if
       bcdat_iota = iota_old*RadialWeight+iota_new*(1.0-RadialWeight)
       bcdat_GHat = GHat_old*RadialWeight+GHat_new*(1.0-RadialWeight)
       bcdat_IHat = IHat_old*RadialWeight+IHat_new*(1.0-RadialWeight)
       bcdata00 = B0_old*RadialWeight+B0_new*(1.0-RadialWeight)
       pPrimeHat = pPrimeHat_old*RadialWeight+pPrimeHat_new*(1.0-RadialWeight)

       bcdata00L=B0_old
       bcdata00H=B0_new
       NHarmonicsL = no_of_modes_old
       NHarmonicsH = no_of_modes_new
       allocate(BHarmonics_lL(NHarmonicsL))
       allocate(BHarmonics_nL(NHarmonicsL))
       allocate(BHarmonics_amplitudesL(NHarmonicsL))
       allocate(BHarmonics_parityL(NHarmonicsL))
       allocate(BHarmonics_lH(NHarmonicsH))
       allocate(BHarmonics_nH(NHarmonicsH))
       allocate(BHarmonics_amplitudesH(NHarmonicsH))
       allocate(BHarmonics_parityH(NHarmonicsH))
       BHarmonics_lL = modesm_old(1:NHarmonicsL)
       BHarmonics_nL = modesn_old(1:NHarmonicsL)
       BHarmonics_amplitudesL = modesb_old(1:NHarmonicsL)
       BHarmonics_lH = modesm_new(1:NHarmonicsH)
       BHarmonics_nH = modesn_new(1:NHarmonicsH)
       BHarmonics_amplitudesH = modesb_new(1:NHarmonicsH)
       BHarmonics_parityL = .true.
       BHarmonics_parityH = .true.


       if (GHat*psiAHat>0) then
          !Note that GHat and psiAHat already have the opposite sign to the corresponding quantities in the .bc file
          !Therefore, the flip is performed if they have the same sign here.
          !print *,"This is a stellarator symmetric file from Joachim Geiger. It will now be turned 180 degrees around a horizontal axis <=> flip the sign of G and I, so that it matches the sign of its total toroidal flux."
          bcdat_GHat    =-bcdat_GHat
          GHat_new=-GHat_new
          GHat_old=-GHat_old
          bcdat_IHat    =-bcdat_IHat
          IHat_new=-IHat_new
          IHat_old=-IHat_old
       end if

       if (.not. nearbyRadiiGiven) then
          BHarmonics_n=BHarmonics_n*(-1) !toroidal direction sign switch
       else
          BHarmonics_nL=BHarmonics_nL*(-1) !toroidal direction sign switch
          BHarmonics_nH=BHarmonics_nH*(-1) !toroidal direction sign switch
       end if

       !Here I may implement to double-check whether bcdat_GHat=GHat, bcdat_IHat=IHat, bcdat_iota=iota ...
       !Not implemented yet


    case (0)! Read a bcdat file with NON-stellarator-symmetric data

       ! Skip a line
       read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile

       do 
          if ((rN_new .ge. rN_wish) .or. end_of_file) exit

          rN_old = rN_new
          no_of_modes_old = no_of_modes_new
          modesm_old = modesm_new
          modesn_old = modesn_new
          modesb_old = modesb_new
          iota_old = iota_new
          GHat_old = GHat_new
          IHat_old = IHat_new
          B0_old = B0_new
          pPrimeHat_old = pPrimeHat_new
          numB0s = 0

          ! Skip a line:
          read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
          ! Read the header for the magnetic surface:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) surfHeader

          rN_new = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
          iota_new = surfHeader(2)
          ! Note that G and I has a minus sign in the following two lines
          ! because Ampere's law comes with a minus sign in the left-handed
          ! (r,pol,tor) system.
          GHat_new = -surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
          IHat_new = -surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
          pPrimeHat_new = surfheader(5)/psiAHat*(4*pi*1e-7)   ! dpdpsi=pPrimeHat/mu_0

          ! Skip units line:
          read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
          proceed = .true.
          modeind = 0
          do
             if (.not. proceed) exit

             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             if (didFileAccessWork /= 0) then
                proceed = .false.
                end_of_file = .true.
             else if (index(lineOfFile,"s") > 0) then
                ! Next flux surface has been reached
                proceed = .false.
             else
                read(unit=lineOfFile, fmt=*) dataIntegers, data2Numbers
                if (dataIntegers(1) == 0 .and. dataIntegers(2) == 0) then
                   B0_new = data2Numbers(1)
                   numB0s = numB0s + 1
                else !!!!if (abs(data2Numbers(1)) > min_bcdata_to_load) then
                   if (modeind + 2 > max_no_of_modes) then
                      print *,"The value of max_no_of_modes in geometry.F90 was insufficient."
                      print *,"Increase this value and recompile."
                      stop
                   end if
                   modeind = modeind + 1
                   modesm_new(modeind) = dataIntegers(1)
                   modesn_new(modeind) = dataIntegers(2)
                   modesb_new(modeind) = data2Numbers(1) !Cosinus component
                   modeind = modeind + 1
                   modesm_new(modeind) = dataIntegers(1)
                   modesn_new(modeind) = dataIntegers(2)
                   modesb_new(modeind) = data2Numbers(2) !Sinus component
                end if
             end if
          end do
          if (numB0s == 0) then
             print *,"Warning: no (0,0) mode found in the bcdat file ",bcdatFile
          else if (numB0s > 1) then
             print *,"Error: more than 1 (0,0) mode found in the bcdat file ",bcdatFile
          end if
          no_of_modes_new = modeind
       end do
       
       close(unit = fileUnit)
       if (masterProc) then
          print *,"Successfully read data from bcdat file ",trim(bcdatFile)
       end if

       DeltapsiHat = psiAHat * (rN_new*rN_new-rN_old*rN_old)
       nearbyRadiiGiven = .true.

       if (VMECRadialOption == 1) then !Choose the nearest flux surface available
          if (abs(rN_old - rN_wish) < abs(rN_new - rN_wish)) then
             RadialWeight = 1.0
             rN = rN_old
          else
             RadialWeight = 0.0
             rN = rN_new
          end if
       else !Linear interpolation in s=rN^2
          RadialWeight = (rN_new*rN_new-rN_wish*rN_wish) / (rN_new*rN_new-rN_old*rN_old)
          rN   = rN_wish
       end if
       bcdat_iota = iota_old*RadialWeight+iota_new*(1.0-RadialWeight)
       bcdat_GHat = GHat_old*RadialWeight+GHat_new*(1.0-RadialWeight)
       bcdat_IHat = IHat_old*RadialWeight+IHat_new*(1.0-RadialWeight)
       bcdata00 = B0_old*RadialWeight+B0_new*(1.0-RadialWeight)
       pPrimeHat = pPrimeHat_old*RadialWeight+pPrimeHat_new*(1.0-RadialWeight)
       
       bcdata00L=B0_old
       bcdata00H=B0_new
       NHarmonicsL = no_of_modes_old
       NHarmonicsH = no_of_modes_new
       allocate(BHarmonics_lL(NHarmonicsL))
       allocate(BHarmonics_nL(NHarmonicsL))
       allocate(BHarmonics_amplitudesL(NHarmonicsL))
       allocate(BHarmonics_parityL(NHarmonicsL))
       allocate(BHarmonics_lH(NHarmonicsH))
       allocate(BHarmonics_nH(NHarmonicsH))
       allocate(BHarmonics_amplitudesH(NHarmonicsH))
       allocate(BHarmonics_parityH(NHarmonicsH))
       BHarmonics_lL = modesm_old(1:NHarmonicsL)
       BHarmonics_nL = modesn_old(1:NHarmonicsL)
       BHarmonics_amplitudesL = modesb_old(1:NHarmonicsL)
       BHarmonics_lH = modesm_new(1:NHarmonicsH)
       BHarmonics_nH = modesn_new(1:NHarmonicsH)
       BHarmonics_amplitudesH = modesb_new(1:NHarmonicsH)
       do i = 0, NHarmonicsL/2-1
          BHarmonics_parityL(2*i+1)=.true.
          BHarmonics_parityL(2*i+2)=.false.
       end do
       do i = 0, NHarmonicsH/2-1
          BHarmonics_parityH(2*i+1)=.true.
          BHarmonics_parityH(2*i+2)=.false.
       end do

       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       !(The toroidal direction sign switch psiAHat=psiAHat*(-1) was already made in the initializeGeometry routine)
       bcdat_GHat     = bcdat_GHat*(-1)                   !toroidal direction sign switch
       GHat_new = GHat_new*(-1)               !toroidal direction sign switch
       GHat_old = GHat_old*(-1)               !toroidal direction sign switch
       bcdat_iota     = bcdat_iota*(-1)                   !toroidal direction sign switch
       iota_new = iota_new*(-1)               !toroidal direction sign switch
       iota_old = iota_old*(-1)               !toroidal direction sign switch
       if (.not. nearbyRadiiGiven) then
          BHarmonics_n=BHarmonics_n*(-1) !toroidal direction sign switch
       else
          BHarmonics_nL=BHarmonics_nL*(-1) !toroidal direction sign switch
          BHarmonics_nH=BHarmonics_nH*(-1) !toroidal direction sign switch
       end if

       !Here I may implement to double-check whether bcdat_GHat=GHat, bcdat_IHat=IHat, bcdat_iota=iota ...
       !Not implemented yet

    case default
       print *,"Error! Invalid symmetry flag in bcdat file!"
       stop
    end select

    if (.not. nearbyRadiiGiven) then
       ! Initialize arrays:
       bcdata = bcdata00 ! This includes the (0,0) component.
       
       do i = 1, NHarmonics
          if (BHarmonics_parity(i)) then   ! The cosine components of bcdata
             include_mn = .false.
             if ((abs(BHarmonics_n(i))<=int(Nzeta/2.0)).and.(BHarmonics_l(i)<=int(Nzeta/2.0))) then
                include_mn = .true.
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   bcdata(itheta,:) = bcdata(itheta,:) + BHarmonics_amplitudes(i) * &
                        cos(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)
                end do
             end if
          else  ! The sine components of bcdata
             include_mn=.false.
             if ((abs(BHarmonics_n(i))<=int(Nzeta/2.0)).and.(BHarmonics_l(i)<=int(Nzeta/2.0))) then
                include_mn=.true.
             end if
             if (BHarmonics_l(i)==0 .or. real(BHarmonics_l(i))==Ntheta/2.0) then
                if (BHarmonics_n(i)==0 .or. abs(real(BHarmonics_n(i)))==Nzeta/2.0 ) then
                   include_mn=.false.
                end if
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   bcdata(itheta,:) = bcdata(itheta,:) + BHarmonics_amplitudes(i) * &
                        sin(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)
                end do
             end if
          end if
       end do
    else !Two nearby radii L and H are given 
       allocate(bcdataL(Ntheta,Nzeta))

       bcdataL = bcdata00L ! This includes the (0,0) component.
       
       do i = 1, NHarmonicsL
          if (BHarmonics_parityL(i)) then   ! The cosine components of bcdata
             include_mn = .false.
             if ((abs(BHarmonics_nL(i))<=int(Nzeta/2.0)).and.(BHarmonics_lL(i)<=int(Ntheta/2.0))) then
                include_mn = .true.
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   bcdataL(itheta,:) = bcdataL(itheta,:) + BHarmonics_amplitudesL(i) * &
                        cos(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)
                end do
             end if
          else  ! The sine components of bcdata
             include_mn=.false.
             if ((abs(BHarmonics_nL(i))<=int(Nzeta/2.0)).and.(BHarmonics_lL(i)<=int(Ntheta/2.0))) then
                include_mn=.true.
             end if
             if (BHarmonics_lL(i)==0 .or. real(BHarmonics_lL(i))==Ntheta/2.0) then
                if (BHarmonics_nL(i)==0 .or. abs(real(BHarmonics_nL(i)))==Nzeta/2.0 ) then
                   include_mn=.false.
                end if
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   bcdataL(itheta,:) = bcdataL(itheta,:) + BHarmonics_amplitudesL(i) * &
                        sin(BHarmonics_lL(i) * theta(itheta) - NPeriods * BHarmonics_nL(i) * zeta)
                end do
             end if
          end if
       end do
       allocate(bcdataH(Ntheta,Nzeta))

       bcdataH = bcdata00H ! This includes the (0,0) component.

       do i = 1, NHarmonicsH
          if (BHarmonics_parityH(i)) then   ! The cosine components of bcdata
             include_mn = .false.
             if ((abs(BHarmonics_nH(i))<=int(Nzeta/2.0)).and.(BHarmonics_lH(i)<=int(Ntheta/2.0))) then
                include_mn = .true.
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   bcdataH(itheta,:) = bcdataH(itheta,:) + BHarmonics_amplitudesH(i) * &
                        cos(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)
                end do
             end if
          else  ! The sine components of bcdata
             include_mn=.false.
             if ((abs(BHarmonics_nH(i))<=int(Nzeta/2.0)).and.(BHarmonics_lH(i)<=int(Ntheta/2.0))) then
                include_mn=.true.
             end if
             if (BHarmonics_lH(i)==0 .or. real(BHarmonics_lH(i))==Ntheta/2.0) then
                if (BHarmonics_nH(i)==0 .or. abs(real(BHarmonics_nH(i)))==Nzeta/2.0 ) then
                   include_mn=.false.
                end if
             end if
             if (Nzeta==1) then
                include_mn = .true.
             end if
             if (include_mn) then
                do itheta = 1,Ntheta
                   bcdataH(itheta,:) = bcdataH(itheta,:) + BHarmonics_amplitudesH(i) * &
                        sin(BHarmonics_lH(i) * theta(itheta) - NPeriods * BHarmonics_nH(i) * zeta)
                end do
             end if
          end if
       end do
       do itheta = 1,Ntheta
          bcdata(itheta,:) = bcdataL(itheta,:)*RadialWeight + bcdataH(itheta,:)*(1.0-RadialWeight)
       end do      

    end if

    if (.not. nearbyRadiiGiven) then
       deallocate(BHarmonics_l)
       deallocate(BHarmonics_n)
       deallocate(BHarmonics_amplitudes)
       deallocate(BHarmonics_parity)
    else
       deallocate(BHarmonics_lL)
       deallocate(BHarmonics_nL)
       deallocate(BHarmonics_amplitudesL)
       deallocate(BHarmonics_parityL)
       deallocate(BHarmonics_lH)
       deallocate(BHarmonics_nH)
       deallocate(BHarmonics_amplitudesH)
       deallocate(BHarmonics_parityH)
    end if


  end subroutine load_bcdat_file

  ! -----------------------------------------------------------------------------------------

end module geometry

