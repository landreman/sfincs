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

  use globalVariables
  use radialCoordinates
  use petscsysdef
  use readVMEC

  implicit none

#include <finclude/petscsysdef.h>

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

    integer :: fileUnit, didFileAccessWork
    character(len=200) :: lineOfFile
    integer, dimension(4) :: headerIntegers
    PetscScalar, dimension(3) :: headerReals

    select case (geometryScheme)
    case (1)
       NPeriods = max(1, helicity_n)

    case (2)
       NPeriods = 10
       aHat = 0.5585d+0 ! (meters)
       psiAHat = (aHat ** 2) / two
       rN_wish = 0.5
       inputRadialCoordinate = 3

       if (masterProc) then
          print *,"---------------------------------------------------------"
          print *,"Since geometryScheme=2, we will ignore the *_wish parameters and use the flux surface rN = 0.5."
       end if

    case (3)
       NPeriods = 10
       aHat = 0.5400d+0 ! (meters)
       psiAHat = (aHat ** 2) / two
       rN_wish = 0.5
       inputRadialCoordinate = 3

       if (masterProc) then
          print *,"---------------------------------------------------------"
          print *,"Since geometryScheme=3, we will ignore the *_wish parameters and use the flux surface rN = 0.5."
       end if

    case (4)
       NPeriods = 5
       aHat = 0.5109d+0 ! (meters)
       psiAHat = -0.384935d+0 ! Tesla * meters^2 / radian
       rN_wish = 0.5
       inputRadialCoordinate = 3

       if (masterProc) then
          print *,"---------------------------------------------------------"
          print *,"Since geometryScheme=4, we will ignore the *_wish parameters and use the flux surface rN = 0.5."
       end if

    case (5)
       ! Read VMEC file, defining the effective minor radius aHat to be VMEC's Aminor_p
       call read_VMEC(equilibriumfile)
       NPeriods = vmec%nfp
       psiAHat = vmec%phi(vmec%ns)/(2*pi)
       aHat = vmec%Aminor_p

    case (6)
       ! Read VMEC file, using an effective minor radius aHat which is NOT VMEC's Aminor_p,
       ! but instead a definition used at IPP
       call read_VMEC(equilibriumFile)
       NPeriods = vmec%nfp
       psiAHat = vmec%phi(vmec%ns)/(2*pi)

       ! This next line is not correct and should be fixed!!!
       aHat = vmec%Aminor_p

    case (10)
       print *,"Error! This geometryScheme has not been implemented yet."

    case (11)
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
    case (5,6)
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
    
  end subroutine computeBHat

  ! ---------------------------------------------------------------------------------------
          
  subroutine computeBHat_Boozer

    ! Note that the BHarmonics_amplitudes in this subroutine are normalized by B0, not by BBar!

    implicit none

    integer :: itheta, izeta, NHarmonics, i, m, n
    integer, dimension(:), allocatable :: BHarmonics_l, BHarmonics_n
    PetscScalar, dimension(:), allocatable :: BHarmonics_amplitudes
    logical, dimension(:), allocatable :: BHarmonics_parity
    PetscScalar, dimension(:,:), allocatable :: hHat, uHat, duHatdtheta, duHatdzeta
    PetscScalar :: R0
    
    integer :: fileUnit, didFileAccessWork
    character(len=200) :: lineOfFile
    integer, dimension(4) :: headerIntegers
    PetscScalar, dimension(3) :: headerReals
    PetscScalar, dimension(6) :: surfHeader
    PetscScalar, dimension(4) :: dataNumbers
    PetscScalar, dimension(8) :: data8Numbers
    integer, dimension(2) :: dataIntegers
    integer :: no_of_modes_old, no_of_modes_new, modeind, numB0s, startn, stopn
    PetscScalar :: iota_old, iota_new, G_old, G_new, I_old, I_new
    PetscScalar :: pPrimeHat, pPrimeHat_old, pPrimeHat_new
    logical :: end_of_file, proceed
    integer, parameter :: max_no_of_modes = 10000
    integer, dimension(max_no_of_modes) :: modesm_old, modesm_new, modesn_old, modesn_new
    PetscScalar, dimension(max_no_of_modes) :: modesb_old, modesb_new
    PetscScalar :: rN_old,  rN_new, B0_old, B0_new
    PetscScalar :: hHatHarmonics_amplitude, uHatHarmonics_amplitude


    ! For the BHarmonics_parity array, 
    ! true indicates the contribution to B(theta,zeta) has the form
    ! cos(l * theta - n * zeta)
    ! while false indicates the contribution to B(theta,zeta) has the form
    ! sin(l * theta - n * zeta)

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

       BHat_sub_psi = 0
       dBHat_sub_psi_dtheta = 0
       dBHat_sub_psi_dzeta = 0

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
                    
       BHat_sub_psi = 0
       dBHat_sub_psi_dtheta = 0
       dBHat_sub_psi_dzeta = 0

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

       BHat_sub_psi = 0
       dBHat_sub_psi_dtheta = 0
       dBHat_sub_psi_dzeta = 0

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

       BHat_sub_psi = 0
       dBHat_sub_psi_dtheta = 0
       dBHat_sub_psi_dzeta = 0

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
          iota_old = 0
          G_old = 0
          I_old = 0
          B0_old = 0
          pPrimeHat_old = 0

          rN_new = 0
          no_of_modes_new = 0
          modesm_new = 0
          modesn_new = 0
          modesb_new = 0
          iota_new = 0
          G_new = 0
          I_new = 0
          B0_new = 0
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
             iota_old = iota_new
             G_old = G_new
             I_old = I_new
             B0_old = B0_new
             pPrimeHat_old = pPrimeHat_new
             numB0s = 0

             ! Skip a line:
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             ! Read the header for the magnetic surface:
             read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) surfHeader

             rN_new = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
             iota_new = surfHeader(2)
             ! Note that G and I have a minus sign in the following two lines
             ! because Ampere's law comes with a minus sign in the left-handed
             ! (r,pol,tor) system.
             G_new = -surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
             I_new = -surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
             pPrimeHat_new = surfheader(5)*(4*pi*1e-7)       ! p=pHat \bar{B}^2 / \mu_0

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
                   read(unit=lineOfFile, fmt=*) dataIntegers, dataNumbers
                   if (dataIntegers(1) == 0 .and. dataIntegers(2) == 0) then
                      B0_new = dataNumbers(4)
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

       if (abs(rN_old - rN_wish) < abs(rN_new - rN_wish)) then
          iota = iota_old
          GHat = G_old
          IHat = I_old
          rN = rN_old
          B0OverBBar = B0_old
          NHarmonics = no_of_modes_old
          pPrimeHat=pPrimeHat_old
          allocate(BHarmonics_l(NHarmonics))
          allocate(BHarmonics_n(NHarmonics))
          allocate(BHarmonics_amplitudes(NHarmonics))
          allocate(BHarmonics_parity(NHarmonics))
          BHarmonics_parity = .true.
          BHarmonics_l = modesm_old(1:NHarmonics)
          BHarmonics_n = modesn_old(1:NHarmonics)
          BHarmonics_amplitudes = modesb_old(1:NHarmonics)
       else
          iota = iota_new
          GHat = G_new
          IHat = I_new
          rN = rN_new
          B0OverBBar = B0_new
          NHarmonics = no_of_modes_new
          pPrimeHat=pPrimeHat_new
          allocate(BHarmonics_l(NHarmonics))
          allocate(BHarmonics_n(NHarmonics))
          allocate(BHarmonics_amplitudes(NHarmonics))
          allocate(BHarmonics_parity(NHarmonics))
          BHarmonics_parity = .true.
          BHarmonics_l = modesm_new(1:NHarmonics)
          BHarmonics_n = modesn_new(1:NHarmonics)
          BHarmonics_amplitudes = modesb_new(1:NHarmonics)
       end if
       dGdpHat=(G_new-G_old)/(rN_new*rN_new-rN_old*rN_old)/pPrimeHat

       BHarmonics_amplitudes = BHarmonics_amplitudes / B0OverBBar

       ! These next lines could be replaced with the actual values from the equilibrium:
       BHat_sub_psi = 0
       dBHat_sub_psi_dtheta = 0
       dBHat_sub_psi_dzeta = 0

       if (GHat*psiAHat<0) then
          !print *,"This is a stellarator symmetric file from Joachim Geiger. It will now be turned 180 degrees around a horizontal axis <=> flip the sign of G and I, so that it matches the sign of its total toroidal flux."
          GHat=-GHat
          IHat=-IHat
          dGdpHat=-dGdpHat
       end if
       
       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       psiAHat=psiAHat*(-1)           !toroidal direction sign switch
       GHat = GHat*(-1)               !toroidal direction sign switch
       iota = iota*(-1)               !toroidal direction sign switch
       BHarmonics_n=BHarmonics_n*(-1) !toroidal direction sign switch

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
          iota_old = 0
          G_old = 0
          I_old = 0
          B0_old = 0
          pPrimeHat_old = 0

          rN_new = 0
          no_of_modes_new = 0
          modesm_new = 0
          modesn_new = 0
          modesb_new = 0
          iota_new = 0
          G_new = 0
          I_new = 0
          B0_new = 0
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
             iota_old = iota_new
             G_old = G_new
             I_old = I_new
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
             G_new = -surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
             I_new = -surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
             pPrimeHat_new = surfheader(5)*(4*pi*1e-7)       ! p=pHat \bar{B}^2 / \mu_0

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
                      modesb_new(modeind) = data8Numbers(7) !Cosinus component
                      modeind = modeind + 1
                      modesm_new(modeind) = dataIntegers(1)
                      modesn_new(modeind) = dataIntegers(2)
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

       if (abs(rN_old - rN_wish) < abs(rN_new - rN_wish)) then
          iota = iota_old
          GHat = G_old
          IHat = I_old
          rN = rN_old
          B0OverBBar = B0_old
          NHarmonics = no_of_modes_old
          pPrimeHat=pPrimeHat_old
          allocate(BHarmonics_l(NHarmonics))
          allocate(BHarmonics_n(NHarmonics))
          allocate(BHarmonics_amplitudes(NHarmonics))
          BHarmonics_l = modesm_old(1:NHarmonics)
          BHarmonics_n = modesn_old(1:NHarmonics)
          BHarmonics_amplitudes = modesb_old(1:NHarmonics)
       else
          iota = iota_new
          GHat = G_new
          IHat = I_new
          rN = rN_new
          B0OverBBar = B0_new
          NHarmonics = no_of_modes_new
          pPrimeHat=pPrimeHat_new
          allocate(BHarmonics_l(NHarmonics))
          allocate(BHarmonics_n(NHarmonics))
          allocate(BHarmonics_amplitudes(NHarmonics))
          BHarmonics_l = modesm_new(1:NHarmonics)
          BHarmonics_n = modesn_new(1:NHarmonics)
          BHarmonics_amplitudes = modesb_new(1:NHarmonics)
       end if
       dGdpHat=(G_new-G_old)/(rN_new*rN_new-rN_old*rN_old)/pPrimeHat

       allocate(BHarmonics_parity(NHarmonics))
       do i = 0, NHarmonics/2-1
          BHarmonics_parity(2*i+1)=.true.
          BHarmonics_parity(2*i+2)=.false.
       end do

       BHarmonics_amplitudes = BHarmonics_amplitudes / B0OverBBar

       ! These next 3 lines could be replaced with the actual data from the geometry input file.
       BHat_sub_psi = 0
       dBHat_sub_psi_dtheta = 0
       dBHat_sub_psi_dzeta = 0

       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       psiAHat=psiAHat*(-1)           !toroidal direction sign switch
       GHat = GHat*(-1)               !toroidal direction sign switch
       iota = iota*(-1)               !toroidal direction sign switch
       BHarmonics_n=BHarmonics_n*(-1) !toroidal direction sign switch

    case default
       print *,"Error! Invalid geometryScheme"
       stop
    end select

    ! Initialize arrays:
    BHat = B0OverBBar ! This includes the (0,0) component.
    dBHatdtheta = 0
    dBHatdzeta = 0

    do i = 1, NHarmonics
       if (BHarmonics_parity(i)) then   ! The cosine components of BHat
          do itheta = 1,Ntheta
             BHat(itheta,:) = BHat(itheta,:) + B0OverBBar * BHarmonics_amplitudes(i) * &
                  cos(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)

             dBHatdtheta(itheta,:) = dBHatdtheta(itheta,:) - B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_l(i) * &
                  sin(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)

             dBHatdzeta(itheta,:) = dBHatdzeta(itheta,:) + B0OverBBar * BHarmonics_amplitudes(i) * Nperiods * BHarmonics_n(i) * &
                  sin(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)

          end do
       else  ! The sine components of BHat
          do itheta = 1,Ntheta
             BHat(itheta,:) = BHat(itheta,:) + B0OverBBar * BHarmonics_amplitudes(i) * &
                  sin(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)

             dBHatdtheta(itheta,:) = dBHatdtheta(itheta,:) + B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_l(i) * &
                  cos(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)

             dBHatdzeta(itheta,:) = dBHatdzeta(itheta,:) - B0OverBBar * BHarmonics_amplitudes(i) * Nperiods * BHarmonics_n(i) * &
                  cos(BHarmonics_l(i) * theta(itheta) - NPeriods * BHarmonics_n(i) * zeta)

          end do
       end if
    end do
    ! ---------------------------------------------------------------------------------------
    ! Calculate parallel current u from harmonics of 1/B^2. Used in NTV calculation.
    ! \nabla_\parallel u = (2/B^4) \nabla B \times \vector{B} \cdot \iota \nabla \psi 
    ! ---------------------------------------------------------------------------------------
    allocate(hHat(Ntheta,Nzeta))
    allocate(uHat(Ntheta,Nzeta))
    allocate(duHatdtheta(Ntheta,Nzeta))
    allocate(duHatdzeta(Ntheta,Nzeta))
    
    uHat = 0
    duHatdtheta = 0
    duHatdzeta = 0
    do itheta = 1,Ntheta
       do izeta = 1,Nzeta
          hHat(itheta,izeta) = 1.0 / (BHat(itheta,izeta)*BHat(itheta,izeta))
       end do
    end do
    
    if (any(.not. BHarmonics_parity)) then !sine components exist
       do m = 0,int(Ntheta/2.0-1) !Nyquist max freq.
          if (m == 0) then
             startn=1
          else
             startn=-int(Nzeta/2.0)
          end if
          stopn=int(Nzeta/2.0-1)
          do n = startn,stopn 
             !cos
             hHatHarmonics_amplitude = 0
             do itheta = 1,Ntheta
                hHatHarmonics_amplitude = hHatHarmonics_amplitude + 2.0/(Ntheta*Nzeta) * &
                     dot_product(cos(m * theta(itheta)  - n * NPeriods * zeta), hHat(itheta,:))
             end do
             uHatHarmonics_amplitude = &
                  iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude
             do itheta = 1,Ntheta
                uHat(itheta,:) = uHat(itheta,:) &
                     + uHatHarmonics_amplitude * cos(m * theta(itheta) - n * NPeriods * zeta)
                duHatdtheta(itheta,:) = duHatdtheta(itheta,:) &
                     - uHatHarmonics_amplitude * m * sin(m * theta(itheta) - n * NPeriods * zeta)
                duHatdzeta(itheta,:) = duHatdzeta(itheta,:) &
                     + uHatHarmonics_amplitude * n * NPeriods * sin(m * theta(itheta) - n * NPeriods * zeta)
             end do

             !sin
             hHatHarmonics_amplitude = 0
             do itheta = 1,Ntheta
                hHatHarmonics_amplitude = hHatHarmonics_amplitude + 2.0/(Ntheta*Nzeta) * &
                     dot_product(sin(m * theta(itheta)  - n * NPeriods * zeta), hHat(itheta,:))
             end do
             uHatHarmonics_amplitude = &
                  iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude
             do itheta = 1,Ntheta
                uHat(itheta,:) = uHat(itheta,:) &
                     + uHatHarmonics_amplitude * sin(m * theta(itheta) - n * NPeriods * zeta)
                duHatdtheta(itheta,:) = duHatdtheta(itheta,:) &
                     + uHatHarmonics_amplitude * m * cos(m * theta(itheta) - n * NPeriods * zeta)
                duHatdzeta(itheta,:) = duHatdzeta(itheta,:) &
                     - uHatHarmonics_amplitude * n * NPeriods * cos(m * theta(itheta) - n * NPeriods * zeta)
             end do
          end do
       end do
    else !only cosinus components
       do m = 0,int(Ntheta/2.0-1) !Nyquist max freq.
          if (m == 0) then
             startn=1
          else
             startn=-int(Nzeta/2.0)
          end if
          stopn=int(Nzeta/2.0-1)
          do n = startn,stopn 
             !cos
             hHatHarmonics_amplitude = 0
             do itheta = 1,Ntheta
                hHatHarmonics_amplitude = hHatHarmonics_amplitude + 2.0/(Ntheta*Nzeta) * &
                     dot_product(cos(m * theta(itheta)  - n * NPeriods * zeta), hHat(itheta,:))
             end do
             uHatHarmonics_amplitude = &
                  iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude
             do itheta = 1,Ntheta
                uHat(itheta,:) = uHat(itheta,:) &
                     + uHatHarmonics_amplitude * cos(m * theta(itheta) - n * NPeriods * zeta)
                duHatdtheta(itheta,:) = duHatdtheta(itheta,:) &
                     - uHatHarmonics_amplitude * m * sin(m * theta(itheta) - n * NPeriods * zeta)
                duHatdzeta(itheta,:) = duHatdzeta(itheta,:) &
                     + uHatHarmonics_amplitude * n * NPeriods * sin(m * theta(itheta) - n * NPeriods * zeta)
             end do
          end do
       end do
    end if

    do itheta = 1,Ntheta
       do izeta = 1,Nzeta
          NTVkernel(itheta,izeta) = 2.0/5.0 * ( &
               dGdpHat / BHat(itheta,izeta) * (iota * dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta)) &
               + 1.0/2.0 * (iota * (duHatdtheta(itheta,izeta) &
                                + uHat(itheta,izeta) * 2.0/BHat(itheta,izeta) * dBHatdtheta(itheta,izeta)) &
                        + duHatdzeta(itheta,izeta) & 
                        + uHat(itheta,izeta) * 2.0/BHat(itheta,izeta) * dBHatdzeta(itheta,izeta)) )
       end do
    end do

    deallocate(BHarmonics_l)
    deallocate(BHarmonics_n)
    deallocate(BHarmonics_amplitudes)
    deallocate(BHarmonics_parity)


    ! Set the Jacobian and various other components of B:

    DHat = BHat * BHat / (GHat + iota * IHat)
    BHat_sup_theta = iota * DHat
    BHat_sup_zeta = DHat
    BHat_sub_theta = IHat
    BHat_sub_zeta = GHat

    ! Eventually we could replace the next lines with a proper calculation:

    dBHat_sub_theta_dpsiHat = 0
    dBHat_sub_theta_dzeta = 0

    dBHat_sub_zeta_dpsiHat = 0
    dBHat_sub_zeta_dtheta = 0

    dBHat_sup_theta_dpsiHat = 0
    dBHat_sup_theta_dzeta = 0

    dBHat_sup_zeta_dpsiHat = 0
    dBHat_sup_zeta_dtheta = 0

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
    integer :: i, j, index, isurf, itheta, izeta, m, n
    PetscScalar :: min_dr2, angle, sin_angle, cos_angle, b, b00, temp, dphi, dpsi
    integer :: numSymmetricModesIncluded, numAntisymmetricModesIncluded

    ! This subroutine is written so that only psiN_wish is used, not the other *_wish quantities.

    if (masterProc) then
       print *,"Reading VMEC geometry from file ",trim(equilibriumFile)
    end if

    allocate(vmec_dBHatdpsiHat(vmec%ns))
    allocate(vmec_dBHat_sub_theta_dpsiHat(vmec%ns))
    allocate(vmec_dBHat_sub_zeta_dpsiHat(vmec%ns))

    ! --------------------------------------------------------------------------------
    ! Do some sanity checking to ensure the VMEC arrays have some expected properties.
    ! --------------------------------------------------------------------------------

    if (abs(vmec%phi(1)) > 1d-14) then
       if (masterProc) then
          print *,"Error! VMEC phi array does not begin with 0."
       end if
       stop
    end if

    dphi = vmec%phi(2) - vmec%phi(1)
    do j=3,vmec%ns
       if (abs(vmec%phi(j)-vmec%phi(j-1)-dphi) > 1d-14) then
          if (masterProc) then
             print *,"Error! VMEC phi array is not uniformly spaced."
          end if
          stop
       end if
    end do

    ! phips is on the half-mesh, so skip first point.
    do j=2,vmec%ns
       if (abs(vmec%phips(j)+vmec%phi(vmec%ns)/(2*pi)) > 1d-14) then
          if (masterProc) then
             print *,"Error! VMEC phips array is not constant and equal to -phi(ns)/(2*pi)."
          end if
          stop
       end if
    end do

    ! --------------------------------------------------------------------------------
    ! End of sanity checks.
    ! --------------------------------------------------------------------------------

    allocate(psiN_full(vmec%ns))
    psiN_full = vmec%phi / (2*pi*psiAHat)

    ! Build an array of the half grid
    allocate(psiN_half(vmec%ns-1))
    do i = 1,vmec%ns-1
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
       allocate(dr2(vmec%ns-1))
       dr2 = (psiN_half - psiN_wish) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,vmec%ns-1
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
       allocate(dr2(vmec%ns))
       dr2 = (psiN_full - psiN_wish) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,vmec%ns
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

    ! In VMEC, quantities on the half grid have the same number of array elements (vmec%ns) as quantities on the full grid,
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
       vmecRadialIndex_full(1) = vmec%ns-1
       vmecRadialIndex_full(2) = vmec%ns
       vmecRadialWeight_full(1) = zero
    else
       ! psiN is >= 0 and <1
       ! This is the most common case.
       vmecRadialIndex_full(1) = floor(psiN*(vmec%ns-1))+1
       vmecRadialIndex_full(2) = vmecRadialIndex_full(1) + 1
       vmecRadialWeight_full(1) = vmecRadialIndex_full(1) - psiN*(vmec%ns-one)
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

    elseif (psiN > psiN_half(vmec%ns-1)) then
       if (masterProc) then
          print *,"Warning: extrapolating beyond the end of VMEC's half grid."
          print *,"(Extrapolating towards the last closed flux surface.)"
       end if
       vmecRadialIndex_half(1) = vmec%ns-1
       vmecRadialIndex_half(2) = vmec%ns
       vmecRadialWeight_half(1) = (psiN_half(vmec%ns-1) - psiN) &
            / (psiN_half(vmec%ns-1) - psiN_half(vmec%ns-2))

    elseif (psiN == psiN_half(vmec%ns-1)) then
       ! We are exactly at the last point of the half grid
       vmecRadialIndex_half(1) = vmec%ns-1
       vmecRadialIndex_half(2) = vmec%ns
       vmecRadialWeight_half(1) = zero
    else
       ! psiN is inside the half grid.
       ! This is the most common case.
       vmecRadialIndex_half(1) = floor(psiN*(vmec%ns-1) + 0.5)+1
       if (vmecRadialIndex_half(1) < 2) then
          ! This can occur sometimes due to roundoff error.
          vmecRadialIndex_half(1) = 2
       end if
       vmecRadialIndex_half(2) = vmecRadialIndex_half(1) + 1
       vmecRadialWeight_half(1) = vmecRadialIndex_half(1) - psiN*(vmec%ns-one) - (0.5d+0)
    end if
    vmecRadialWeight_half(2) = one-vmecRadialWeight_half(1)

    if (masterProc) then
!!$       print *,"vmecRadialIndex_full:",vmecRadialIndex_full
!!$       print *,"vmecRadialWeight_full:",vmecRadialWeight_full
!!$       print *,"vmecRadialIndex_half:",vmecRadialIndex_half
!!$       print *,"vmecRadialWeight_half:",vmecRadialWeight_half
       if (abs(vmecRadialWeight_half(1)) < 1e-14) then
          print "(a,i3,a,i3,a)"," Using radial index ",vmecRadialIndex_half(2)," of ",vmec%ns," from vmec's half mesh."
       elseif (abs(vmecRadialWeight_half(2)) < 1e-14) then
          print "(a,i3,a,i3,a)"," Using radial index ",vmecRadialIndex_half(1)," of ",vmec%ns," from vmec's half mesh."
       else
          print "(a,i3,a,i3,a,i3,a)", " Interpolating using radial indices ",vmecRadialIndex_half(1)," and ",vmecRadialIndex_half(2),&
               " of ",vmec%ns," from vmec's half mesh."
       end if
    end if

    iota = vmec%iotas(vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
         + vmec%iotas(vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)

    BHat = zero
    DHat = zero
    dBHatdtheta = zero
    dBHatdzeta = zero

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
    ! other harmonics of |B| are large enough to include.
    m = 0
    n = 0
    b00 = vmec%bmnc(n,m,vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
               + vmec%bmnc(n,m,vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)

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
    do n = -vmec%ntor, vmec%ntor
       do m = 0, vmec%mpol-1
          ! -----------------------------------------------------
          ! First, consider just the stellarator-symmetric terms:
          ! -----------------------------------------------------

          b = vmec%bmnc(n,m,vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
               + vmec%bmnc(n,m,vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)

          if (abs(b/b00) >= min_Bmn_to_load) then
             ! This (m,n) mode is sufficiently large to include.
             !if (masterProc) then
             !   print *,"Including mode with m = ",m,", n = ",n
             !end if
             numSymmetricModesIncluded = numSymmetricModesIncluded + 1

             ! Evaluate the radial derivatives we will need:
             dpsi = vmec%phi(2)/(2*pi)  ! Doesn't need to be in the loops, but here for convenience.

             vmec_dBHatdpsiHat(2:vmec%ns-1) = (vmec%bmnc(n,m,3:vmec%ns) - vmec%bmnc(n,m,2:vmec%ns-1)) / dpsi
             ! Simplistic "extrapolation" at the endpoints:
             vmec_dBHatdpsiHat(1) = vmec_dBHatdpsiHat(2)
             vmec_dBHatdpsiHat(vmec%ns) = vmec_dBHatdpsiHat(vmec%ns-1)

             vmec_dBHat_sub_theta_dpsiHat(2:vmec%ns-1) = (vmec%bsubumnc(n,m,3:vmec%ns) - vmec%bsubumnc(n,m,2:vmec%ns-1)) / dpsi
             vmec_dBHat_sub_theta_dpsiHat(1) = vmec_dBHat_sub_theta_dpsiHat(2)
             vmec_dBHat_sub_theta_dpsiHat(vmec%ns) = vmec_dBHat_sub_theta_dpsiHat(vmec%ns-1)

             vmec_dBHat_sub_zeta_dpsiHat(2:vmec%ns-1) = (vmec%bsubvmnc(n,m,3:vmec%ns) - vmec%bsubvmnc(n,m,2:vmec%ns-1)) / dpsi
             vmec_dBHat_sub_zeta_dpsiHat(1) = vmec_dBHat_sub_zeta_dpsiHat(2)
             vmec_dBHat_sub_zeta_dpsiHat(vmec%ns) = vmec_dBHat_sub_zeta_dpsiHat(vmec%ns-1)

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
                      temp = vmec%gmnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf) / (psiAHat)
                      DHat(itheta,izeta) = DHat(itheta,izeta) + temp * cos_angle

                      ! Handle B sup theta:
                      ! Note that VMEC's bsupumnc and bsupumns are exactly the same as SFINCS's BHat_sup_theta, with no conversion factors of 2pi needed.
                      temp = vmec%bsupumnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      BHat_sup_theta(itheta,izeta) = BHat_sup_theta(itheta,izeta) + temp * cos_angle
                      dBHat_sup_theta_dzeta(itheta,izeta) = dBHat_sup_theta_dzeta(itheta,izeta) + n * NPeriods * temp * sin_angle

                      ! Handle B sup zeta:
                      ! Note that VMEC's bsupvmnc and bsupvmns are exactly the same as SFINCS's BHat_sup_zeta, with no conversion factors of 2pi or Nperiods needed.
                      temp = vmec%bsupvmnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      BHat_sup_zeta(itheta,izeta) = BHat_sup_zeta(itheta,izeta) + temp * cos_angle
                      dBHat_sup_zeta_dtheta(itheta,izeta) = dBHat_sup_zeta_dtheta(itheta,izeta) - m * temp * sin_angle

                      ! Handle B sub theta:
                      ! Note that VMEC's bsubumnc and bsubumns are exactly the same as SFINCS's BHat_sub_theta, with no conversion factors of 2pi needed.
                      temp = vmec%bsubumnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      BHat_sub_theta(itheta,izeta) = BHat_sub_theta(itheta,izeta) + temp * cos_angle
                      dBHat_sub_theta_dzeta(itheta,izeta) = dBHat_sub_theta_dzeta(itheta,izeta) + n * NPeriods * temp * sin_angle

                      ! Handle B sub zeta:
                      ! Note that VMEC's bsubvmnc and bsubvmns are exactly the same as SFINCS's BHat_sub_zeta, with no conversion factors of 2pi needed.
                      temp = vmec%bsubvmnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      BHat_sub_zeta(itheta,izeta) = BHat_sub_zeta(itheta,izeta) + temp * cos_angle
                      dBHat_sub_zeta_dtheta(itheta,izeta) = dBHat_sub_zeta_dtheta(itheta,izeta) - m * temp * sin_angle

                      ! Handle B sub psi.
                      ! Unlike the other components of B, this one is on the full mesh.
                      ! Notice B_psi = B_s * (d s / d psi), and (d s / d psi) = 1 / psiAHat
!!$                      print *,"UUU"
!!$                      print *,"itheta=",itheta
!!$                      print *,"izeta=",izeta
!!$                      print *,"temp=",temp
!!$                      print *,"size(BHat_sub_psi,1)",size(BHat_sub_psi,1)
!!$                      print *,"size(BHat_sub_psi,2)",size(BHat_sub_psi,2)
!!$                      print *,"size(vmec%bsubsmns,1)",size(vmec%bsubsmns,1)
!!$                      print *,"size(vmec%bsubsmns,2)",size(vmec%bsubsmns,2)
!!$                      print *,"size(vmec%bsubsmns,3)",size(vmec%bsubsmns,3)
!!$                      print *,"size(vmec%bsubsmnc,1)",size(vmec%bsubsmnc,1)
!!$                      print *,"size(vmec%bsubsmnc,2)",size(vmec%bsubsmnc,2)
!!$                      print *,"size(vmec%bsubsmnc,3)",size(vmec%bsubsmnc,3)
!!$                      print *,"size(vmec%bsubumnc,1)",size(vmec%bsubumnc,1)
!!$                      print *,"size(vmec%bsubumnc,2)",size(vmec%bsubumnc,2)
!!$                      print *,"size(vmec%bsubumnc,3)",size(vmec%bsubumnc,3)
!!$                      print *,"size(vmec%bmnc,1)",size(vmec%bmnc,1)
!!$                      print *,"size(vmec%bmnc,2)",size(vmec%bmnc,2)
!!$                      print *,"size(vmec%bmnc,3)",size(vmec%bmnc,3)
!!$                      print *,"BHat_sub_psi(itheta,izeta)=",BHat_sub_psi(itheta,izeta)
!!$                      print *,"isurf=",isurf
!!$                      print *,"vmecRadialIndex_full(isurf)",vmecRadialIndex_full(isurf)
                      temp = vmec%bsubsmns(n,m,vmecRadialIndex_full(isurf)) / psiAHat * vmecRadialWeight_full(isurf)
                      BHat_sub_psi(itheta,izeta) = BHat_sub_psi(itheta,izeta) + temp * sin_angle
                      dBHat_sub_psi_dtheta(itheta,izeta) = dBHat_sub_psi_dtheta(itheta,izeta) + m * temp * cos_angle
                      dBHat_sub_psi_dzeta(itheta,izeta)  = dBHat_sub_psi_dzeta(itheta,izeta) - n * NPeriods * temp * cos_angle

                      ! Handle dBHatdpsiHat.
                      ! Since bmnc is on the half mesh, its radial derivative is on the full mesh.
                      temp = vmec_dBHatdpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      dBHatdpsiHat(itheta,izeta) = dBHatdpsiHat(itheta,izeta) + temp * cos_angle

                      ! Handle dBHat_sub_theta_dpsiHat.
                      ! Since bsubumnc is on the half mesh, its radial derivative is on the full mesh.
                      temp = vmec_dBHat_sub_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      dBHat_sub_theta_dpsiHat(itheta,izeta) = dBHat_sub_theta_dpsiHat(itheta,izeta) + temp * cos_angle

                      ! Handle dBHat_sub_zeta_dpsiHat.
                      ! Since bsubvmnc is on the half mesh, its radial derivative is on the full mesh.
                      temp = vmec_dBHat_sub_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      dBHat_sub_zeta_dpsiHat(itheta,izeta) = dBHat_sub_zeta_dpsiHat(itheta,izeta) + temp * cos_angle

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
          ! NOTE: This functionality has not been tested!!!
          ! -----------------------------------------------------

          if (vmec%iasym > 0) then

             b = vmec%bmns(n,m,vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
                  + vmec%bmns(n,m,vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)
             
             if (abs(b/b00) >= min_Bmn_to_load) then
                ! This (m,n) mode is sufficiently large to include.
                !if (masterProc) then
                !   print *,"Including stellarator-asymmetric mode with m = ",m,", n = ",n
                !end if
                numAntisymmetricModesIncluded = numAntisymmetricModesIncluded + 1

                ! Evaluate the radial derivatives we will need:
                dpsi = vmec%phi(2)/(2*pi)  ! Doesn't need to be in the loops, but here for convenience.

                vmec_dBHatdpsiHat(2:vmec%ns-1) = (vmec%bmns(n,m,3:vmec%ns) - vmec%bmns(n,m,2:vmec%ns-1)) / dpsi
                ! Simplistic "extrapolation" at the endpoints:
                vmec_dBHatdpsiHat(1) = vmec_dBHatdpsiHat(2)
                vmec_dBHatdpsiHat(vmec%ns) = vmec_dBHatdpsiHat(vmec%ns-1)

                vmec_dBHat_sub_theta_dpsiHat(2:vmec%ns-1) = (vmec%bsubumns(n,m,3:vmec%ns) - vmec%bsubumns(n,m,2:vmec%ns-1)) / dpsi
                vmec_dBHat_sub_theta_dpsiHat(1) = vmec_dBHat_sub_theta_dpsiHat(2)
                vmec_dBHat_sub_theta_dpsiHat(vmec%ns) = vmec_dBHat_sub_theta_dpsiHat(vmec%ns-1)

                vmec_dBHat_sub_zeta_dpsiHat(2:vmec%ns-1) = (vmec%bsubvmns(n,m,3:vmec%ns) - vmec%bsubvmns(n,m,2:vmec%ns-1)) / dpsi
                vmec_dBHat_sub_zeta_dpsiHat(1) = vmec_dBHat_sub_zeta_dpsiHat(2)
                vmec_dBHat_sub_zeta_dpsiHat(vmec%ns) = vmec_dBHat_sub_zeta_dpsiHat(vmec%ns-1)

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
                         temp = vmec%gmns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf) / (psiAHat)
                         DHat(itheta,izeta) = DHat(itheta,izeta) + temp * sin_angle
                         
                         ! Handle B sup theta:
                         ! Note that VMEC's bsupumnc and bsupumns are exactly the same as SFINCS's BHat_sup_theta, with no conversion factors of 2pi needed.
                         temp = vmec%bsupumns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         BHat_sup_theta(itheta,izeta) = BHat_sup_theta(itheta,izeta) + temp * sin_angle
                         dBHat_sup_theta_dzeta(itheta,izeta) = dBHat_sup_theta_dzeta(itheta,izeta) - n * NPeriods * temp * cos_angle
                         
                         ! Handle B sup zeta:
                         ! Note that VMEC's bsupvmnc and bsupvmns are exactly the same as SFINCS's BHat_sup_zeta, with no conversion factors of 2pi or Nperiods needed.
                         temp = vmec%bsupvmns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         BHat_sup_zeta(itheta,izeta) = BHat_sup_zeta(itheta,izeta) + temp * sin_angle
                         dBHat_sup_zeta_dtheta(itheta,izeta) = dBHat_sup_zeta_dtheta(itheta,izeta) + m * temp * cos_angle
                         
                         ! Handle B sub theta:
                         ! Note that VMEC's bsubumnc and bsubumns are exactly the same as SFINCS's BHat_sub_theta, with no conversion factors of 2pi needed.
                         temp = vmec%bsubumns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         BHat_sub_theta(itheta,izeta) = BHat_sub_theta(itheta,izeta) + temp * sin_angle
                         dBHat_sub_theta_dzeta(itheta,izeta) = dBHat_sub_theta_dzeta(itheta,izeta) - n * NPeriods * temp * cos_angle
                         
                         ! Handle B sub zeta:
                         ! Note that VMEC's bsubvmnc and bsubvmns are exactly the same as SFINCS's BHat_sub_zeta, with no conversion factors of 2pi needed.
                         temp = vmec%bsubvmns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         BHat_sub_zeta(itheta,izeta) = BHat_sub_zeta(itheta,izeta) + temp * sin_angle
                         dBHat_sub_zeta_dtheta(itheta,izeta) = dBHat_sub_zeta_dtheta(itheta,izeta) + m * temp * cos_angle

                         ! Handle B sub psi.
                         ! Unlike the other components of B, this one is on the full mesh.
                         ! Notice B_psi = B_s * (d s / d psi), and (d s / d psi) = 1 / psiAHat
                         temp = vmec%bsubsmnc(n,m,vmecRadialIndex_full(isurf)) / psiAHat * vmecRadialWeight_full(isurf)
                         BHat_sub_psi(itheta,izeta) = BHat_sub_psi(itheta,izeta) + temp * cos_angle
                         dBHat_sub_psi_dtheta(itheta,izeta) = dBHat_sub_psi_dtheta(itheta,izeta) - m * temp * sin_angle
                         dBHat_sub_psi_dzeta(itheta,izeta)  = dBHat_sub_psi_dzeta(itheta,izeta) + n * NPeriods * temp * sin_angle

                         ! Handle dBHatdpsiHat.
                         ! Since bmns is on the half mesh, its radial derivative is on the full mesh.
                         temp = vmec_dBHatdpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         dBHatdpsiHat(itheta,izeta) = dBHatdpsiHat(itheta,izeta) + temp * sin_angle

                         ! Handle dBHat_sub_theta_dpsiHat.
                         ! Since bsubumns is on the half mesh, its radial derivative is on the full mesh.
                         temp = vmec_dBHat_sub_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         dBHat_sub_theta_dpsiHat(itheta,izeta) = dBHat_sub_theta_dpsiHat(itheta,izeta) + temp * sin_angle

                         ! Handle dBHat_sub_zeta_dpsiHat.
                         ! Since bsubvmns is on the half mesh, its radial derivative is on the full mesh.
                         temp = vmec_dBHat_sub_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         dBHat_sub_zeta_dpsiHat(itheta,izeta) = dBHat_sub_zeta_dpsiHat(itheta,izeta) + temp * sin_angle

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
    end do

    if (masterProc) then
       print "(a,i4,a,i4,a)"," Including ",numSymmetricModesIncluded," of ",(2*vmec%ntor+1)*vmec%mpol," stellarator-symmetric modes from the VMEC file."
       if (vmec%iasym > 0) then
          print "(a,i4,a,i4,a)"," Including ",numAntisymmetricModesIncluded," of ",(2*vmec%ntor+1)*vmec%mpol," stellarator-antisymmetric modes from the VMEC file."
       else
          print *,"Equilibrium is stellarator-symmetric."
       end if
    end if

    ! Convert Jacobian to inverse Jacobian:
    DHat = one / DHat

    ! These next lines should be replaced eventually with a proper calculation:
    dBHat_sup_theta_dpsiHat = 0
    dBHat_sup_zeta_dpsiHat = 0

    deallocate(psiN_full)
    deallocate(psiN_half)
    deallocate(vmec_dBHatdpsiHat)
    deallocate(vmec_dBHat_sub_theta_dpsiHat)
    deallocate(vmec_dBHat_sub_zeta_dpsiHat)

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
    do itheta=1,Ntheta
       do izeta=1,Nzeta
          VPrimeHat = VPrimeHat + thetaWeights(itheta) * zetaWeights(izeta) / DHat(itheta,izeta)
          FSABHat2 = FSABHat2 + thetaWeights(itheta) * zetaWeights(izeta) &
               * BHat(itheta,izeta) * BHat(itheta,izeta) / DHat(itheta,izeta)
       end do
    end do

    FSABHat2 = FSABHat2 / VPrimeHat

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

  end subroutine computeBIntegrals


end module geometry

