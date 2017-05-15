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

  use kinds
  use globalVariables
  use radialCoordinates
  use readVMEC

  implicit none

  private

  public :: initialize_geometry, evaluate_B_on_grids

  integer :: vmecRadialIndex_full(2)
  integer :: vmecRadialIndex_half(2)
  real(prec) :: vmecRadialWeight_full(2)
  real(prec) :: vmecRadialWeight_half(2)

  integer, parameter :: max_num_modes = 10000
  integer, dimension(2) :: num_modes
  integer, dimension(max_num_modes,2) :: modes_l, modes_n
  logical, dimension(max_num_modes,2) :: modes_parity
  real(prec), dimension(max_num_modes,2) :: modes_B, modes_R, modes_Z, modes_delta_zeta
  real(prec), dimension(2) :: Boozer_radial_weights

contains

  ! -----------------------------------------------------------------------------------

  subroutine initialize_geometry()
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
    real(prec), dimension(3) :: headerReals

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

    case (5,7)
       ! Read VMEC file, defining the effective minor radius aHat to be VMEC's Aminor_p
       ! geometryScheme=7 is a hack for testing NEMEC_compute_missing_fields.f90
       call read_VMEC(equilibriumfile,geometryScheme==7)
       NPeriods = vmec%nfp
       psiAHat = vmec%phi(vmec%ns)/(2*pi)
       aHat = vmec%Aminor_p

    case (6)
       ! Read Erika Strumberger's NEMEC format used at IPP.
       call read_NEMEC(equilibriumFile)
       NPeriods = vmec%nfp
       psiAHat = vmec%phi(vmec%ns)/(2*pi)
       ! Aminor_p is not set in this format, so we must specify a value in the input namelist.
       if (aHat-0.5585d+0 < 1e-4 .and. masterProc) then
          print *,"###############################################################################################"
          print *,"WARNING: It appears you have not explicitly set aHat to a value other than the default."
          print *,"This is probably wrong. For geometryScheme=6, aHat is not set using the NEMEC equilibrium file."
          print *,"###############################################################################################"
       end if

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Done reading psiAHat, Nperiods, and aHat.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       call load_B_Fourier_modes_Boozer()
    case (5,6,7)
       coordinateSystem = COORDINATE_SYSTEM_VMEC
       call load_B_Fourier_modes_VMEC()
    case default
       print *,"Error! Invalid setting for geometryScheme."
       stop
    end select

  end subroutine initialize_geometry

  ! -----------------------------------------------------------------------------------

  subroutine evaluate_B_on_grids()
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
    !  BHat, dBHatdtheta, dBHatdzeta, dBHatdpsiHat, sqrt_g
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

    integer :: itheta, izeta, level
    real(prec) :: sqrt_g_11

    select case (coordinateSystem)
    case (COORDINATE_SYSTEM_BOOZER)
       call evaluate_B_on_grids_Boozer()
    case (COORDINATE_SYSTEM_VMEC)
       call evaluate_B_on_grids_VMEC()
    case default
       print *,"Error! Invalid setting for coordinateSystem."
       stop
    end select

    sqrt_g_11 = levels(1)%sqrt_g(1,1)

    do level = 1,N_levels
       levels(level)%BDotCurlB = (1/levels(level)%sqrt_g) * (  levels(level)%BHat_sub_theta * levels(level)%dBHat_sub_psi_dzeta &
            - levels(level)%BHat_sub_theta * levels(level)%dBHat_sub_zeta_dpsiHat &
            + levels(level)%BHat_sub_zeta * levels(level)%dBHat_sub_theta_dpsiHat &
            - levels(level)%BHat_sub_zeta * levels(level)%dBHat_sub_psi_dtheta)

       if (.not. force0RadialCurrentInEquilibrium) then
          levels(level)%BDotCurlB = levels(level)%BDotCurlB + levels(level)%BHat_sub_psi / levels(level)%sqrt_g * (levels(level)%dBHat_sub_zeta_dtheta - levels(level)%dBHat_sub_theta_dzeta)
       end if

       ! Validate geometry arrays:
       do itheta=1,levels(level)%Ntheta
          do izeta=1,levels(level)%Nzeta
             if (levels(level)%BHat(itheta,izeta) <= 0) then
                print *,"ERROR! BHat is not everywhere positive!"
                stop
             end if
             if (levels(level)%sqrt_g(itheta,izeta)*sqrt_g_11 <= 0) then
                print *,"ERROR! sqrt_g does not have the same sign everywhere!"
                stop
             end if
          end do
       end do
    end do

  end subroutine evaluate_B_on_grids

  ! ---------------------------------------------------------------------------------------

  subroutine load_B_Fourier_modes_Boozer

    ! Note that the BHarmonics_amplitudes in this subroutine are normalized by B0, not by BBar!

    implicit none

    integer :: j
    integer :: fileUnit, didFileAccessWork
    character(len=200) :: lineOfFile
    integer, dimension(4) :: headerIntegers
    real(prec), dimension(3) :: headerReals
    real(prec), dimension(6) :: surfHeader
    real(prec), dimension(4) :: dataNumbers
    real(prec), dimension(8) :: data8Numbers
    integer, dimension(2) :: dataIntegers
    integer :: modeind, numB0s
    logical :: end_of_file, proceed, include_mn, nearbyRadiiGiven, nonStelSym
    real(prec) :: DeltapsiHat

    ! For the BHarmonics_parity array, 
    ! true indicates the contribution to B(theta,zeta) has the form
    ! cos(l * theta - n * zeta)
    ! while false indicates the contribution to B(theta,zeta) has the form
    ! sin(l * theta - n * zeta)

    diotadpsiHat = 0

    modes_parity = .true.
    num_modes = 0
    Boozer_radial_weights = zero
    Boozer_radial_weights(1) = one

    modes_l = 0
    modes_n = 0
    modes_B = 0
    modes_R = 0
    modes_Z = 0
    modes_delta_zeta = 0

    GHat_surfaces = 0
    IHat_surfaces = 0
    iota_surfaces = 0
    pPrimeHat_surfaces = 0
    rN_surfaces = 0

    select case (geometryScheme)
    case (1)
       ! Three-helicity model:
       ! B = B0 * [1 + epsilon_t * cos(theta) + epsilon_h * cos(helicity_l * theta - helicity_n * zeta) ...
       !             + epsilon_antisymm * sin(helicity_antisymm_l * theta - helicity_antisymm_n * zeta) ...

       num_modes(1) = 3

       i = 1
       modes_l(i,1) = 1
       modes_n(i,1) = 0
       modes_B(i,1) = epsilon_t

       i = 2
       modes_l(i,1) = helicity_l
       if (helicity_n == 0) then
          modes_n(i,1) = 0
       else
          modes_n(i,1) = 1
       end if
       modes_B(i,1) = epsilon_h

       i = 3
       modes_parity(i,1) = .false.
       modes_l(i,1) = helicity_antisymm_l
       if (helicity_n == 0) then
          modes_n(i,1) = helicity_antisymm_n
       else
          modes_n(i,1) = helicity_antisymm_n / helicity_n
       end if
       modes_B(i,1) = epsilon_antisymm

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
       modes_B = modes_B * B0OverBBar

    case (2)
       ! A three-harmonic approximation of the LHD standard configuration.
       ! Values taken from Table 1 of
       ! Beidler et al, Nuclear Fusion 51, 076001 (2011).

       iota = 0.4542d+0

       num_modes(1) = 3

       i = 1
       modes_l(i,1) = 1
       modes_n(i,1) = 0
       modes_B(i,1) = -0.07053d+0

       i = 2
       modes_l(i,1) = 2
       modes_n(i,1) = 1
       modes_B(i,1) = 0.05067d+0

       i = 3
       modes_l(i,1) = 1
       modes_n(i,1) = 1
       modes_B(i,1) = -0.01476d+0

       B0OverBBar = 1.0d+0  ! (Tesla)
       R0 = 3.7481d+0 ! (meters)
       GHat = B0OverBBar * R0
       IHat = 0
       dGdpHat = 0
       !rN = -1 !dummy
       rN = rN_wish
       modes_B = modes_B * B0OverBBar

    case (3)
       ! A four-harmonic approximation of the LHD inward-shifted configuration.
       ! Values taken from Table 1 of
       ! Beidler et al, Nuclear Fusion 51, 076001 (2011).

       iota = 0.4692d+0

       num_modes(1) = 4

       i = 1
       modes_l(i,1) = 1
       modes_n(i,1) = 0
       modes_B(i,1) = -0.05927d+0

       i = 2
       modes_l(i,1) = 2
       modes_n(i,1) = 1
       modes_B(i,1) = 0.05267d+0

       i = 3
       modes_l(i,1) = 1
       modes_n(i,1) = 1
       modes_B(i,1) = -0.04956d+0

       i = 4
       modes_l(i,1) = 0
       modes_n(i,1) = 1
       modes_B(i,1) = 0.01045d+0

       B0OverBBar = 1.0d+0  ! (Tesla)
       R0 = 3.6024d+0 ! (meters)
       GHat = B0OverBBar * R0
       IHat = 0
       dGdpHat = 0
       !normradius = -1 !dummy
       rN = rN_wish
       modes_B = modes_B * B0OverBBar

    case (4)
       ! A three-harmonic approximation of the W7-X standard configuration.
       ! Values taken from Table 1 of
       ! Beidler et al, Nuclear Fusion 51, 076001 (2011).

       iota = 0.8700d+0

       num_modes(1) = 3

       i = 1
       modes_l(i,1) = 0
       modes_n(i,1) = 1
       modes_B(i,1) = 0.04645d+0

       i = 2
       modes_l(i,1) = 1
       modes_n(i,1) = 1
       modes_B(i,1) = -0.04351d+0

       i = 3
       modes_l(i,1) = 1
       modes_n(i,1) = 0
       modes_B(i,1) = -0.01902d+0

       B0OverBBar = 3.089d+0  ! (Tesla)
       R0 = 5.5267d+0 ! (meters)
       GHat = -17.885d+0
       IHat = 0
       dGdpHat = 0
       !normradius = -1 !dummy
       rN = rN_wish
       modes_B = modes_B * B0OverBBar

    case (11, 12)
       ! Read Boozer coordinate file in .bc format used at IPP Greifswald

       fileUnit = 11
       open(unit=fileUnit, file=equilibriumFile, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file ",equilibriumFile
          stop
       end if

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

       rN_surfaces = 0
       num_modes = 0
       iota_surfaces = 0
       GHat_surfaces = 0
       IHat_surfaces = 0
       B0 = 0
       R0 = 0
       pPrimeHat = 0
       num_surfaces_read = 0

       ! Skip a line
       read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile

       do 
          if ((rN_surfaces(1) .ge. rN_wish) .or. end_of_file) exit

          ! Set new surface 3 equal to previous surface 2,
          ! then set new surface 2 equal to previous surface 1:
          do j = 1,1 ! Eventually I may want to change this to include 3 surfaces
             modes_l(1:num_modes(j),j+1) = modes_l(1:num_modes(j),j)
             modes_n(1:num_modes(j),j+1) = modes_n(1:num_modes(j),j)
             modes_B(1:num_modes(j),j+1) = modes_B(1:num_modes(j),j)
             modes_R(1:num_modes(j),j+1) = modes_R(1:num_modes(j),j)
             modes_Z(1:num_modes(j),j+1) = modes_Z(1:num_modes(j),j)
             modes_zeta_shift(1:num_modes(j),j+1) = modes_zeta_shift(1:num_modes(j),j)
             iota_surfaces(j+1) = iota_surfaces(j)
             GHat_surfaces(j+1) = GHat_surfaces(j)
             IHat_surfaces(j+1) = IHat_surfaces(j)
             pPrimeHat_surfaces(j+1) = pPrimeHat_surfaces(j)
             rN_surfaces(j+1) = rN_surfaces(j)
             num_modes_(j+1) = num_modes(j)
          end do
          B0 = 0
          R0 = 0
          numB0s = 0

          ! Skip a line:
          read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
          ! Read the header for the magnetic surface:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) surfHeader

          rN_surfaces(1) = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
          iota_surfaces(1) = surfHeader(2)
          ! Note that G and I have a minus sign in the following two lines
          ! because Ampere's law comes with a minus sign in the left-handed
          ! (r,pol,tor) system.
          GHat_surfaces(1) = -surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
          IHat_surfaces(1) = -surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
          pPrimeHat_surfaces(1) = surfheader(5)/psiAHat*(4*pi*1e-7)   ! dpdpsi=pPrimeHat/mu_0

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
                   if (geometryScheme==11) then
                      B0 = dataNumbers(4)
                   else
                      B0 = dataNumbers(7)
                   end if
                   R0 = dataNumbers(1)
                   numB0s = numB0s + 1
                else if (abs(dataNumbers(4)) > min_Bmn_to_load) then
                   if (modeind + 2 > max_num_modes) then
                      print *,"The value of max_no_of_modes in geometry.F90 was insufficient."
                      print *,"Either increase this value and recompile, or else increase min_Bmn_to_load."
                      stop
                   end if
                   if (geometryScheme==11) then
                      modeind = modeind + 1
                      modes_l(modeind,1) = dataIntegers(1)
                      modes_n(modeind,1) = dataIntegers(2)
                      modes_R(modeind,1) = dataNumbers(1)
                      modes_Z(modeind,1) = dataNumbers(2)
                      modes_delta_zeta(modeind,1) = dataNumbers(3)
                      modes_B(modeind,1) = dataNumbers(4)
                   else
                      ! geometryScheme==12
                      modeind = modeind + 1
                      modes_l(modeind) = dataIntegers(1)
                      modes_n(modeind) = dataIntegers(2)
                      modes_R(modeind) = data8Numbers(1) !Cosinus component
                      modes_Z(modeind) = data8Numbers(4) !Sinus component
                      modes_delta_zeta(modeind)= data8Numbers(6) !Sinus component
                      modes_B(modeind) = data8Numbers(7) !Cosinus component
                      modeind = modeind + 1
                      modes_l(modeind) = dataIntegers(1)
                      modes_n(modeind) = dataIntegers(2)
                      modes_R(modeind) = data8Numbers(2) !Sinus component
                      modes_Z(modeind) = data8Numbers(3) !Cosinus component
                      modes_delta_zeta(modeind)= data8Numbers(5) !Cosinus component
                      modes_B(modeind) = data8Numbers(8) !Sinus component
                      modes_parity(modeind) = .false.
                   end if
                end if
             end if
          end do
          if (numB0s == 0) then
             print *,"Error: no (0,0) mode found in magnetic equilibrium file ",equilibriumFile
          else if (numB0s > 1) then
             print *,"Error: more than 1 (0,0) mode found in magnetic equilibrium file ",equilibriumFile
          end if
          num_modes(1) = modeind
          num_surfaces_read = num_surface_read + 1
       end do


       close(unit = fileUnit)
       if (masterProc) print *,"Successfully read magnetic equilibrium from file ",trim(equilibriumFile)

       DeltapsiHat = psiAHat * (rN_surface(1)*rN_surface(1) - rN_surface(2)*rN_surface(2))

       if (num_surfaces_read==0) then
          if (masterProc) print *,"Error! No surfaces found in equilibrium file ",trim(equilibriumFile)
          stop
       elseif (num_surfaces_read==1) then
          if (masterProc) print *,"Only 1 magnetic surface found in the equilibrium file."
          Boozer_radial_weights(1) = one
          Boozer_radial_weights(2) = zero
       else
          ! At least 2 surfaces found in the equilibrium file
          if (VMECRadialOption == 1) then !Choose the nearest flux surface available
             if (abs(rN_surfaces(1) - rN_wish) < abs(rN_surfaces(2) - rN_wish)) then
                Boozer_radial_weights(1) = one
                Boozer_radial_weights(2) = zero
                rN = rN_surfaces(1)
             else
                Boozer_radial_weights(1) = zero
                Boozer_radial_weights(2) = one
                rN = rN_surfaces(2)
             end if
          else !Linear interpolation in s=rN^2
             Boozer_radial_weights(2) = (rN_surfaces(1)**2 - rN_wish**2) / (rN_surfaces(1)**2 - rN_surfaces(2)**2)
             Boozer_radial_weights(1) = one - Boozer_radial_weights(1)
             rN   = rN_wish
          end if
       end if
       if (masterProc) then
          print *,"Boozer_radial_weights:",Boozer_radial_weights
          print *,"num_modes:",num_modes
       end if

       iota = dot_product(Boozer_radial_weights, iota_surfaces)
       GHat = dot_product(Boozer_radial_weights, GHat_surfaces)
       IHat = dot_product(Boozer_radial_weights, IHat_surfaces)
       B0OverBBar = dot_product(Boozer_radial_weights, B0_surfaces)
       R0 = dot_product(Boozer_radial_weights, R0_surfaces)
       pPrimeHat = dot_product(Boozer_radial_weights, pPrimeHat_surfaces)

       !Unnnecessary step:
       !BHarmonics_amplitudes = BHarmonics_amplitudes / B0OverBBar

       if (GHat*psiAHat>0) then
          !Note that GHat and psiAHat already have the opposite sign to the corresponding quantities in the .bc file
          !Therefore, the flip is performed if they have the same sign here.
          !print *,"This is a stellarator symmetric file from Joachim Geiger. It will now be turned 180 degrees around a horizontal axis <=> flip the sign of G and I, so that it matches the sign of its total toroidal flux."
          GHat    =-GHat
          GHat_surfaces = -GHat_surfaces
          IHat    =-IHat
          IHat_surfaces = -IHat_surfaces
          !dGdpHat=-dGdpHat
       end if

       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       !(The toroidal direction sign switch psiAHat=psiAHat*(-1) was already made in the initializeGeometry routine)
       GHat          = -GHat
       GHat_surfaces = -GHat_surfaces
       IHat          = -IHat
       IHat_surfaces = -IHat_surfaces
       iota          = -iota
       iota_surfaces = -iota_surfaces
       modes_n = -modes_n
       modes_delta_zeta = -modes_delta_zeta

       ! Compute radial derivatives. 
       diotadpsiHat= (iota_surfaces(1)-iota_surfaces(2))/DeltapsiHat

    case default
       print *,"Error! Invalid geometryScheme"
       stop
    end select


  end subroutine load_B_Fourier_modes_Boozer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine evaluate_B_on_grids_Boozer()

    implicit none

    integer :: level, isurf, iln, itheta, izeta, l, n
    real(prec) :: angle, sinangle, cosangle

    do level = 1,N_levels

       ! Initialize arrays that are set by sums in the loops below:
       levels(level)%BHat = 0
       levels(level)%dBHatdtheta = 0
       levels(level)%dBHatdzeta = 0

       do isurf = 1,2
          do iln = 1,num_modes(isurf)
             l = modes_l(iln,isurf)
             n = modes_n(iln,isurf)
             include_mn = .false.
             if ((abs(n)<=int(Nzeta/2.0)).and.(abs(l)<=int(Ntheta/2.0))) include_mn = .true.
             if ((.not. modes_parity(iln,isurf)) .and. (l==0 .or. real(abs(l))==Ntheta/2.0)) then
                if (l==0 .or. abs(real(n))==Nzeta/2.0 ) include_mn=.false.
             end if
             if (Nzeta==1) include_mn = .true.
             if (level==1 .and. masterProc) then
                if (include_mn) then
                   print *,"Including the mode with l=",l," m=",m" isurf=",isurf
                else
                   print *,"NOT including the mode with l=",l," m=",m" isurf=",isurf
                end if
             end if
             if (include_mn) then
                do itheta = 1,levels(level)%Ntheta
                   do izeta = 1,levels(level)%Nzeta
                      angle = l * levels(level)%theta(itheta) - NPeriods * n * levels(level)%zeta(izeta)
                      sinangle = sin(angle)
                      cosangle = cos(angle)

                      if (modes_parity(iln,isurf)) then
                         ! Stellarator-symmetric modes
                         levels(level)%BHat(itheta,izeta) = levels(level)%BHat(itheta,izeta) &
                              + modes_B(iln,isurf) * cosangle * Boozer_radial_weights(isurf)

                         levels(level)%dBHatdtheta(itheta,izeta) = levels(level)%dBHatdtheta(itheta,izeta) &
                              - l * modes_B(iln,isurf) * sinangle * Boozer_radial_weights(isurf)

                         levels(level)%dBHatdzeta(itheta,izeta) = levels(level)%dBHatdzeta(itheta,izeta) &
                              + n * NPeriods * modes_B(iln,isurf) * sinangle * Boozer_radial_weights(isurf)
                      else
                         ! Stellarator-antisymmetric modes
                         levels(level)%BHat(itheta,izeta) = levels(level)%BHat(itheta,izeta) &
                              + modes_B(iln,isurf) * sinangle * Boozer_radial_weights(isurf)

                         levels(level)%dBHatdtheta(itheta,izeta) = levels(level)%dBHatdtheta(itheta,izeta) &
                              + l * modes_B(iln,isurf) * cosangle * Boozer_radial_weights(isurf)

                         levels(level)%dBHatdzeta(itheta,izeta) = levels(level)%dBHatdzeta(itheta,izeta) &
                              - n * NPeriods * modes_B(iln,isurf) * cosangle * Boozer_radial_weights(isurf)
                      end if
                   end do
                end do
             end if
          end do
       end do

       ! Set the Jacobian and various other components of B:
       levels(level)%sqrt_g = (GHat + iota * IHat) / (levels(level)%BHat ** 2)
       levels(level)%BHat_sup_theta = iota / levels(level)%sqrt_g
       levels(level)%BHat_sup_zeta = 1/level(level)%sqrt_g
       levels(level)%BHat_sub_theta = IHat
       levels(level)%BHat_sub_zeta = GHat

       ! These next arrays turn out to be 0 for Boozer coordinates:
       levels(level)%dBHat_sub_theta_dzeta = 0
       levels(level)%dBHat_sub_zeta_dtheta = 0

       ! These next arrays are generally nonzero, but I haven't gotten around to setting them properly:
       levels(level)%dBHatdpsiHat = 0
       levels(level)%BHat_sub_psi = 0
       levels(level)%dBHat_sub_psi_dtheta = 0
       levels(level)%dBHat_sub_psi_dzeta = 0
       levels(level)%dBHat_sub_theta_dpsiHat = 0
       levels(level)%dBHat_sub_zeta_dpsiHat = 0
    end do

!!$
!!$
!!$    if (nearbyRadiiGiven) then
!!$       dBHat_sup_zeta_dpsiHat = 2.0 * BHat *dBHatdpsiHat / (GHat + iota * IHat) &
!!$            -(dBHat_sub_zeta_dpsiHat + iota*dBHat_sub_theta_dpsiHat + diotadpsiHat*IHat) &
!!$            / (GHat + iota * IHat) / (GHat + iota * IHat)
!!$       dBHat_sup_zeta_dtheta = 2.0 * BHat *dBHatdtheta / (GHat + iota * IHat)
!!$
!!$       dBHat_sup_theta_dpsiHat = iota * dBHat_sup_zeta_dpsiHat + diotadpsiHat / sqrt_g
!!$       dBHat_sup_theta_dzeta = iota * 2.0 * BHat *dBHatdzeta / (GHat + iota * IHat)
!!$    end if
!!$
!!$    !possible double-check
!!$    !dBHatdpsiHat= sqrt(hHat)/2.0*( dBHat_sup_theta_dpsiHat*BHat_sub_theta &
!!$    !                              +dBHat_sub_theta_dpsiHat*BHat_sup_theta &
!!$    !                              +dBHat_sup_zeta_dpsiHat *BHat_sub_zeta &
!!$    !                              +dBHat_sub_zeta_dpsiHat *BHat_sup_zeta)
  end subroutine evaluate_B_on_grids_Boozer

  ! -----------------------------------------------------------------------------------------

  subroutine init_B_Fourier_modes_VMEC

    implicit none

    real(prec), dimension(:), allocatable :: dr2, psiN_full, psiN_half
    integer :: i, j, index
    real(prec) :: min_dr2, b, temp, dphi, dpsi

    ! This subroutine is written so that only psiN_wish is used, not the other *_wish quantities.

    if (masterProc) then
       print *,"Reading VMEC geometry from file ",trim(equilibriumFile)
    end if

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
       if (abs(vmec%phi(j)-vmec%phi(j-1)-dphi) > 1d-11) then
          if (masterProc) then
             print *,"Error! VMEC phi array is not uniformly spaced."
          end if
          stop
       end if
    end do

    ! phips is on the half-mesh, so skip first point.
    do j=2,vmec%ns
       if (abs(vmec%phips(j)+vmec%phi(vmec%ns)/(2*pi)) > 1d-11) then
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
       vmecRadialIndex_half(1) = floor(psiN*(vmec%ns-1) + 0.5d+0)+1
       if (vmecRadialIndex_half(1) < 2) then
          ! This can occur sometimes due to roundoff error.
          vmecRadialIndex_half(1) = 2
       end if
       vmecRadialIndex_half(2) = vmecRadialIndex_half(1) + 1
       vmecRadialWeight_half(1) = vmecRadialIndex_half(1) - psiN*(vmec%ns-one) - (0.5d+0)
    end if
    vmecRadialWeight_half(2) = one-vmecRadialWeight_half(1)

    if (masterProc) then
       if (abs(vmecRadialWeight_half(1)) < 1e-14) then
          print "(a,i3,a,i3,a)"," Using radial index ",vmecRadialIndex_half(2)," of ",vmec%ns," from vmec's half mesh."
       elseif (abs(vmecRadialWeight_half(2)) < 1e-14) then
          print "(a,i3,a,i3,a)"," Using radial index ",vmecRadialIndex_half(1)," of ",vmec%ns," from vmec's half mesh."
       else
          print "(a,i3,a,i3,a,i3,a)", " Interpolating using radial indices ",vmecRadialIndex_half(1)," and ",vmecRadialIndex_half(2),&
               " of ",vmec%ns," from vmec's half mesh."
          print "(a,f17.14,a,f17.14)", " Weights for half mesh = ",vmecRadialWeight_half(1)," and ",vmecRadialWeight_half(2)
          print "(a,i3,a,i3,a,i3,a)", " Interpolating using radial indices ",vmecRadialIndex_full(1)," and ",vmecRadialIndex_full(2),&
               " of ",vmec%ns," from vmec's full mesh."
          print "(a,f17.14,a,f17.14)", " Weights for full mesh = ",vmecRadialWeight_full(1)," and ",vmecRadialWeight_full(2)
       end if
    end if

    iota = vmec%iotas(vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
         + vmec%iotas(vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)

    rN = sqrt(psiN)

    deallocate(psiN_full)
    deallocate(psiN_half)

  end subroutine init_B_Fourier_modes_VMEC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine evaluate_B_on_grids_VMEC()

    implicit none

    integer :: level
    real(prec), dimension(:), allocatable :: vmec_dBHatdpsiHat, vmec_dBHat_sub_theta_dpsiHat, vmec_dBHat_sub_zeta_dpsiHat
    integer :: i, j, isurf, itheta, izeta, m, n
    real(prec) :: angle, sin_angle, cos_angle, b, b00, temp, dphi, dpsi
    integer :: numSymmetricModesIncluded, numAntisymmetricModesIncluded
    real(prec) :: scale_factor

    ! This subroutine is written so that only psiN_wish is used, not the other *_wish quantities.

    allocate(vmec_dBHatdpsiHat(vmec%ns))
    allocate(vmec_dBHat_sub_theta_dpsiHat(vmec%ns))
    allocate(vmec_dBHat_sub_zeta_dpsiHat(vmec%ns))

    do level = 1,N_levels

       ! Initialize arrays to 0:
       levels(level)%BHat = zero
       levels(level)%sqrt_g = zero
       levels(level)%dBHatdtheta = zero
       levels(level)%dBHatdzeta = zero
       levels(level)%dBHatdpsiHat = zero

       levels(level)%BHat_sub_psi = zero
       levels(level)%dBHat_sub_psi_dtheta = zero
       levels(level)%dBHat_sub_psi_dzeta = zero

       levels(level)%BHat_sub_theta = zero
       levels(level)%dBHat_sub_theta_dpsiHat = zero
       levels(level)%dBHat_sub_theta_dzeta = zero

       levels(level)%BHat_sub_zeta = zero
       levels(level)%dBHat_sub_zeta_dpsiHat = zero
       levels(level)%dBHat_sub_zeta_dtheta = zero

       levels(level)%BHat_sup_theta = zero
       levels(level)%BHat_sup_zeta = zero

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

             ! Set scale_factor to rippleScale for non-axisymmetric or non-quasisymmetric modes
             scaleFactor = set_scale_factor(n,m)
             b = b*scale_factor

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

                do itheta = 1,levels(level)%Ntheta
                   do izeta = 1,levels(level)%Nzeta
                      angle = m * levels(level)%theta(itheta) - n * NPeriods * levels(level)%zeta(izeta)
                      cos_angle = cos(angle)
                      sin_angle = sin(angle)

                      levels(level)%BHat(itheta,izeta) = levels(level)%BHat(itheta,izeta) + b * cos_angle

                      levels(level)%dbHatdtheta(itheta,izeta) = levels(level)%dBHatdtheta(itheta,izeta) - m * b * sin_angle

                      levels(level)%dbHatdzeta(itheta,izeta) = levels(level)%dBHatdzeta(itheta,izeta) + n * NPeriods * b * sin_angle

                      do isurf = 1,2
                         ! Handle Jacobian:
                         ! SFINCS's sqrt_g is the Jacobian of the (psiHat, theta, zeta) coordinates.
                         ! VMEC's gmnc and gmns are the Jacobian of the (psiN, theta, zeta) coordinates.
                         ! Because one uses psiHat and the other uses psiN, we need a factor of psiAHat for conversion.
                         temp = vmec%gmnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf) / (psiAHat)
                         temp = temp*scale_factor
                         levels(level)%sqrt_g(itheta,izeta) = levels(level)%sqrt_g(itheta,izeta) + temp * cos_angle

                         ! Handle B sup theta:
                         ! Note that VMEC's bsupumnc and bsupumns are exactly the same as SFINCS's BHat_sup_theta, with no conversion factors of 2pi needed.
                         temp = vmec%bsupumnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         temp = temp*scale_factor
                         levels(level)%BHat_sup_theta(itheta,izeta) = levels(level)%BHat_sup_theta(itheta,izeta) + temp * cos_angle
                         levels(level)%dBHat_sup_theta_dzeta(itheta,izeta) = levels(level)%dBHat_sup_theta_dzeta(itheta,izeta) + n * NPeriods * temp * sin_angle

                         ! Handle B sup zeta:
                         ! Note that VMEC's bsupvmnc and bsupvmns are exactly the same as SFINCS's BHat_sup_zeta, with no conversion factors of 2pi or Nperiods needed.
                         temp = vmec%bsupvmnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         temp = temp*scale_factor
                         levels(level)%BHat_sup_zeta(itheta,izeta) = levels(level)%BHat_sup_zeta(itheta,izeta) + temp * cos_angle
                         levels(level)%dBHat_sup_zeta_dtheta(itheta,izeta) = levels(level)%dBHat_sup_zeta_dtheta(itheta,izeta) - m * temp * sin_angle

                         ! Handle B sub theta:
                         ! Note that VMEC's bsubumnc and bsubumns are exactly the same as SFINCS's BHat_sub_theta, with no conversion factors of 2pi needed.
                         temp = vmec%bsubumnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         temp = temp*scale_factor
                         levels(level)%BHat_sub_theta(itheta,izeta) = levels(level)%BHat_sub_theta(itheta,izeta) + temp * cos_angle
                         levels(level)%dBHat_sub_theta_dzeta(itheta,izeta) = levels(level)%dBHat_sub_theta_dzeta(itheta,izeta) + n * NPeriods * temp * sin_angle

                         ! Handle B sub zeta:
                         ! Note that VMEC's bsubvmnc and bsubvmns are exactly the same as SFINCS's BHat_sub_zeta, with no conversion factors of 2pi needed.
                         temp = vmec%bsubvmnc(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         temp = temp*scale_factor
                         levels(level)%BHat_sub_zeta(itheta,izeta) = levels(level)%BHat_sub_zeta(itheta,izeta) + temp * cos_angle
                         levels(level)%dBHat_sub_zeta_dtheta(itheta,izeta) = levels(level)%dBHat_sub_zeta_dtheta(itheta,izeta) - m * temp * sin_angle

                         ! Handle B sub psi.
                         ! Unlike the other components of B, this one is on the full mesh.
                         ! Notice B_psi = B_s * (d s / d psi), and (d s / d psi) = 1 / psiAHat
                         temp = vmec%bsubsmns(n,m,vmecRadialIndex_full(isurf)) / psiAHat * vmecRadialWeight_full(isurf)
                         temp = temp*scale_factor
                         levels(level)%BHat_sub_psi(itheta,izeta) = levels(level)%BHat_sub_psi(itheta,izeta) + temp * sin_angle
                         levels(level)%dBHat_sub_psi_dtheta(itheta,izeta) = levels(level)%dBHat_sub_psi_dtheta(itheta,izeta) + m * temp * cos_angle
                         levels(level)%dBHat_sub_psi_dzeta(itheta,izeta)  = levels(level)%dBHat_sub_psi_dzeta(itheta,izeta) - n * NPeriods * temp * cos_angle

                         ! Handle dBHatdpsiHat.
                         ! Since bmnc is on the half mesh, its radial derivative is on the full mesh.
                         temp = vmec_dBHatdpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         temp = temp*scale_factor
                         levels(level)%dBHatdpsiHat(itheta,izeta) = levels(level)%dBHatdpsiHat(itheta,izeta) + temp * cos_angle

                         ! Handle dBHat_sub_theta_dpsiHat.
                         ! Since bsubumnc is on the half mesh, its radial derivative is on the full mesh.
                         temp = vmec_dBHat_sub_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         temp = temp*scale_factor
                         levels(level)%dBHat_sub_theta_dpsiHat(itheta,izeta) = levels(level)%dBHat_sub_theta_dpsiHat(itheta,izeta) + temp * cos_angle

                         ! Handle dBHat_sub_zeta_dpsiHat.
                         ! Since bsubvmnc is on the half mesh, its radial derivative is on the full mesh.
                         temp = vmec_dBHat_sub_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         temp = temp*scale_factor
                         levels(level)%dBHat_sub_zeta_dpsiHat(itheta,izeta) = levels(level)%dBHat_sub_zeta_dpsiHat(itheta,izeta) + temp * cos_angle

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
             ! NOTE: This functionality has not been tested as thoroughly as the stellarator-symmetric case.
             ! -----------------------------------------------------

             if (vmec%iasym > 0) then

                b = vmec%bmns(n,m,vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
                     + vmec%bmns(n,m,vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)

                ! Set scale_factor to rippleScale for non-axisymmetric or non-quasisymmetric modes
                scale_factor = set_scale_factor(n,m)
                b = b*scale_factor

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

                   do itheta = 1,levels(level)%Ntheta
                      do izeta = 1,levels(level)%Nzeta
                         angle = m * levels(level)%theta(itheta) - n * NPeriods * levels(level)%zeta(izeta)
                         cos_angle = cos(angle)
                         sin_angle = sin(angle)

                         levels(level)%BHat(itheta,izeta) = levels(level)%BHat(itheta,izeta) + b * sin_angle
                         levels(level)%dbHatdtheta(itheta,izeta) = levels(level)%dBHatdtheta(itheta,izeta) + m * b * cos_angle
                         levels(level)%dbHatdzeta(itheta,izeta) = levels(level)%dBHatdzeta(itheta,izeta) - n * NPeriods * b * cos_angle

                         do isurf = 1,2
                            ! Handle Jacobian:
                            ! SFINCS's sqrt_g is the Jacobian of the (psiHat, theta, zeta) coordinates.
                            ! VMEC's gmnc and gmns are the Jacobian of the (psiN, theta, zeta) coordinates.
                            ! Because one uses psiHat and the other uses psiN, we need a factor of psiAHat for conversion.
                            temp = vmec%gmns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf) / (psiAHat)
                            temp = temp*scale_factor
                            levels(level)%sqrt_g(itheta,izeta) = levels(level)%sqrt_g(itheta,izeta) + temp * sin_angle

                            ! Handle B sup theta:
                            ! Note that VMEC's bsupumnc and bsupumns are exactly the same as SFINCS's BHat_sup_theta, with no conversion factors of 2pi needed.
                            temp = vmec%bsupumns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                            temp = temp*scale_factor	
                            levels(level)%BHat_sup_theta(itheta,izeta) = levels(level)%BHat_sup_theta(itheta,izeta) + temp * sin_angle
                            levels(level)%dBHat_sup_theta_dzeta(itheta,izeta) = levels(level)%dBHat_sup_theta_dzeta(itheta,izeta) - n * NPeriods * temp * cos_angle

                            ! Handle B sup zeta:
                            ! Note that VMEC's bsupvmnc and bsupvmns are exactly the same as SFINCS's BHat_sup_zeta, with no conversion factors of 2pi or Nperiods needed.
                            temp = vmec%bsupvmns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                            temp = temp*scale_factor
                            levels(level)%BHat_sup_zeta(itheta,izeta) = levels(level)%BHat_sup_zeta(itheta,izeta) + temp * sin_angle
                            levels(level)%dBHat_sup_zeta_dtheta(itheta,izeta) = levels(level)%dBHat_sup_zeta_dtheta(itheta,izeta) + m * temp * cos_angle

                            ! Handle B sub theta:
                            ! Note that VMEC's bsubumnc and bsubumns are exactly the same as SFINCS's BHat_sub_theta, with no conversion factors of 2pi needed.
                            temp = vmec%bsubumns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                            temp = temp*scale_factor
                            levels(level)%BHat_sub_theta(itheta,izeta) = levels(level)%BHat_sub_theta(itheta,izeta) + temp * sin_angle
                            levels(level)%dBHat_sub_theta_dzeta(itheta,izeta) = levels(level)%dBHat_sub_theta_dzeta(itheta,izeta) - n * NPeriods * temp * cos_angle

                            ! Handle B sub zeta:
                            ! Note that VMEC's bsubvmnc and bsubvmns are exactly the same as SFINCS's BHat_sub_zeta, with no conversion factors of 2pi needed.
                            temp = vmec%bsubvmns(n,m,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                            temp = temp*scale_factor
                            levels(level)%BHat_sub_zeta(itheta,izeta) = levels(level)%BHat_sub_zeta(itheta,izeta) + temp * sin_angle
                            levels(level)%dBHat_sub_zeta_dtheta(itheta,izeta) = levels(level)%dBHat_sub_zeta_dtheta(itheta,izeta) + m * temp * cos_angle

                            ! Handle B sub psi.
                            ! Unlike the other components of B, this one is on the full mesh.
                            ! Notice B_psi = B_s * (d s / d psi), and (d s / d psi) = 1 / psiAHat
                            temp = vmec%bsubsmnc(n,m,vmecRadialIndex_full(isurf)) / psiAHat * vmecRadialWeight_full(isurf)
                            temp = temp*scale_factor
                            levels(level)%BHat_sub_psi(itheta,izeta) = levels(level)%BHat_sub_psi(itheta,izeta) + temp * cos_angle
                            levels(level)%dBHat_sub_psi_dtheta(itheta,izeta) = levels(level)%dBHat_sub_psi_dtheta(itheta,izeta) - m * temp * sin_angle
                            levels(level)%dBHat_sub_psi_dzeta(itheta,izeta)  = levels(level)%dBHat_sub_psi_dzeta(itheta,izeta) + n * NPeriods * temp * sin_angle

                            ! Handle dBHatdpsiHat.
                            ! Since bmns is on the half mesh, its radial derivative is on the full mesh.
                            temp = vmec_dBHatdpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                            temp = temp*scale_factor
                            levels(level)%dBHatdpsiHat(itheta,izeta) = levels(level)%dBHatdpsiHat(itheta,izeta) + temp * sin_angle

                            ! Handle dBHat_sub_theta_dpsiHat.
                            ! Since bsubumns is on the half mesh, its radial derivative is on the full mesh.
                            temp = vmec_dBHat_sub_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                            temp = temp*scale_factor
                            levels(level)%dBHat_sub_theta_dpsiHat(itheta,izeta) = levels(level)%dBHat_sub_theta_dpsiHat(itheta,izeta) + temp * sin_angle

                            ! Handle dBHat_sub_zeta_dpsiHat.
                            ! Since bsubvmns is on the half mesh, its radial derivative is on the full mesh.
                            temp = vmec_dBHat_sub_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                            temp = temp*scale_factor
                            levels(level)%dBHat_sub_zeta_dpsiHat(itheta,izeta) = levels(level)%dBHat_sub_zeta_dpsiHat(itheta,izeta) + temp * sin_angle

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

       if (masterProc .and. level==1) then
          print "(a,i4,a,i4,a)"," Including ",numSymmetricModesIncluded," of ",(2*vmec%ntor+1)*vmec%mpol," stellarator-symmetric modes from the VMEC file."
          if (vmec%iasym > 0) then
             print "(a,i4,a,i4,a)"," Including ",numAntisymmetricModesIncluded," of ",(2*vmec%ntor+1)*vmec%mpol," stellarator-antisymmetric modes from the VMEC file."
          else
             print *,"Equilibrium is stellarator-symmetric."
          end if
       end if

       ! These next lines should be replaced eventually with a proper calculation.
       ! We set these arrays to 0 since they are not actually used.
       levels(level)%dBHat_sup_theta_dpsiHat = 0
       levels(level)%dBHat_sup_zeta_dpsiHat = 0

       deallocate(vmec_dBHatdpsiHat)
       deallocate(vmec_dBHat_sub_theta_dpsiHat)
       deallocate(vmec_dBHat_sub_zeta_dpsiHat)

    end do

    if (masterProc) print *,"Successfully set geometry using VMEC file ",trim(equilibriumFile)

  end subroutine evaluate_B_on_grids_VMEC

  ! ------------------------------------------------------------------------------------------

  subroutine compute_B_integrals

    implicit none

    integer :: itheta, izeta

    ! This subroutine computes VPrimeHat, FSABHat2, and (if needed) B0OverBBar, GHat, and IHat.

    VPrimeHat = 0
    FSABHat2 = 0
    do itheta=1,Ntheta
       do izeta=1,Nzeta
          VPrimeHat = VPrimeHat + levels(1)%thetaWeights(itheta) * levels(1)%zetaWeights(izeta) * levels(1)%sqrt_g(itheta,izeta)
          FSABHat2 = FSABHat2 + levels(1)%thetaWeights(itheta) * levels(1)%zetaWeights(izeta) &
               * levels(1)%BHat(itheta,izeta) * levels(1)%BHat(itheta,izeta) * levels(1)%sqrt_g(itheta,izeta)
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
             B0OverBBar = B0OverBBar + levels(1)%thetaWeights(itheta) * levels(1)%zetaWeights(izeta) &
                  * (levels(1)%BHat(itheta,izeta) ** 3) * levels(1)%sqrt_g(itheta,izeta)
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
       GHat = dot_product(levels(1)%thetaWeights, matmul(levels(1)%BHat_sub_zeta,  levels(1)%zetaWeights)) / (4*pi*pi)
       IHat = dot_product(levels(1)%thetaWeights, matmul(levels(1)%BHat_sub_theta, levels(1)%zetaWeights)) / (4*pi*pi)

       if (RHSMode==3) then
          ! Monoenergetic coefficient computation.
          ! Overwrite nu_n and dPhiHatd* using nuPrime and EStar.

          ! 20170331 MJL: The absolute value below is needed to ensure the matrix remains diagonally dominant if GHat<0. If we use RHSMode=3 for a non-stellarator-symmetric plasma, the coefficients may not be independent of sgn(nu), so some thought should be given to the signs.
          nu_n = abs(nuPrime * B0OverBBar / (GHat + iota * IHat))
          dPhiHatdpsiHat = 2 / (gamma * Delta) * EStar * iota * B0OverBBar / GHat
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

    end if

  end subroutine compute_B_integrals

  ! Set scale to value of rippleScale for non-axisymmetric or non-quasisymmetric components
  function set_scale_factor(n,m) result(scale)
    integer :: n
    integer :: m
    real(prec) :: scale 
    scale = 1
    if (helicity_n == 0 .and. n /= 0) then
       scale = rippleScale
    else if ((n /= 0) .and. (helicity_l/helicity_n) /= (m/n)) then
       scale = rippleScale
    else if (helicity_n /= 0 .and. n == 0) then
       scale = rippleScale
    end if
  end function set_scale_factor

end module geometry

