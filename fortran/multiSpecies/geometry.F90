module geometry

  use globalVariables
  use petscsysdef

  implicit none

#include <finclude/petscsysdef.h>

contains

  ! -----------------------------------------------------------------------------------

  subroutine setNPeriods()
    ! This subroutine sets NPeriods, which is the number of identical toroidal segments
    ! in the stellarator (e.g. 5 for W7-X, 10 for LHD, 4 for HSX.)

    implicit none

    integer :: fileUnit, didFileAccessWork
    character(len=200) :: lineOfFile
    integer, dimension(4) :: numbers

    select case (geometryScheme)
    case (1)
       NPeriods = max(1, helicity_n)
    case (2,3)
       NPeriods = 10
    case (4)
       NPeriods = 5
    case (10)
       print *,"Error! This geometryScheme has not been implemented yet."

    case (11)
       fileUnit = 11
       open(unit=fileUnit, file=JGboozer_file, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file."
          stop
       else
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             ! Skip lines that begin with "CC":
             if (lineOfFile(1:2) /= "CC") exit
          end do
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) numbers
          if (didFileAccessWork /= 0) then
                print *,"Unable to read number of toroidal periods from the magnetic equilibrium file."
             stop
          else
             NPeriods = numbers(4)
          end if

       end if

       close(unit = fileUnit)
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,&
               "] Successfully opened magnetic equilibrium file ",trim(JGboozer_file),".  Nperiods = ",Nperiods
       end if

    case (12)
       fileUnit = 11
       open(unit=fileUnit, file=JGboozer_file_NonStelSym, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file."
          stop
       else
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             ! Skip lines that begin with "CC":
             if (lineOfFile(1:2) /= "CC") exit
          end do
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) numbers
          if (didFileAccessWork /= 0) then
                print *,"Unable to read number of toroidal periods from the magnetic equilibrium file."
             stop
          else
             NPeriods = numbers(4)
          end if

       end if

       close(unit = fileUnit)
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,&
               "] Successfully opened magnetic equilibrium file ",trim(JGboozer_file_NonStelSym),".  Nperiods = ",Nperiods
       end if

    case default
       print *,"Error! Invalid setting for geometryScheme."
       stop
    end select

  end subroutine setNPeriods

  ! -----------------------------------------------------------------------------------

  subroutine computeBHat()
    ! This subroutine evaluates BHat, dBHatdtheta, and dBHatdzeta on the (theta,zeta) grid.

    ! Note that the BHarmonics_amplitudes are normalized by B0, not by BBar!

    implicit none

    integer :: itheta, izeta, NHarmonics, i, m, n
    integer, dimension(:), allocatable :: BHarmonics_l, BHarmonics_n
    PetscScalar, dimension(:), allocatable :: BHarmonics_amplitudes
    logical, dimension(:), allocatable :: BHarmonics_parity
    PetscScalar, dimension(:,:), allocatable :: hHat, uHat, duHatdtheta, duHatdzeta
    PetscScalar :: a, R0
    
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
    PetscScalar :: normradius_old,  normradius_new, B0_old, B0_new
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
       normradius = -1 !dummy

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
       a = 0.5585d+0 ! (meters)
       GHat = B0OverBBar * R0
       IHat = 0
       psiAHat = B0OverBBar * (a ** 2) / two
       dGdpHat = 0
       normradius = -1 !dummy
                    
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
       a = 0.5400d+0 ! (meters)
       GHat = B0OverBBar * R0
       IHat = 0
       psiAHat = B0OverBBar * (a ** 2) / two
       dGdpHat = 0
       normradius = -1 !dummy

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
       a = 0.5109d+0 ! (meters)
       GHat = -17.885d+0
       IHat = 0
       psiAHat = -0.384935d+0 ! Tesla * meters^2 / radian
       dGdpHat = 0
       normradius = -1 !dummy

    case (11)
       ! Read Boozer coordinate file in .bc format used at IPP Greifswald

       fileUnit = 11
       open(unit=fileUnit, file=JGboozer_file, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file ",JGboozer_file
          stop
       else
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             if (lineOfFile(1:2) /= "CC") exit
          end do

          ! Read header line:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) headerIntegers, headerReals
          if (didFileAccessWork /= 0) then
             print *,"Unable to read header from the magnetic equilibrium file ",JGboozer_file
             stop
          end if

          NPeriods = headerIntegers(4)
          psiAHat  = headerReals(1)/2/pi !Convert the flux from Tm^2 to Tm^2/rad
          a        = headerReals(2)      !minor radius in meters

          end_of_file = .false.

          normradius_old = 0
          no_of_modes_old = 0
          modesm_old = 0
          modesn_old = 0
          modesb_old = 0
          iota_old = 0
          G_old = 0
          I_old = 0
          B0_old = 0
          pPrimeHat_old = 0

          normradius_new = 0
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
             if ((normradius_new .ge. normradius_wish) .or. end_of_file) exit

             normradius_old = normradius_new
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

             normradius_new = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
             iota_new = surfHeader(2)
             G_new = surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
             I_new = surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
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
                print *,"Error: no (0,0) mode found in magnetic equilibrium file ",JGboozer_file
             else if (numB0s > 1) then
                print *,"Error: more than 1 (0,0) mode found in magnetic equilibrium file ",JGboozer_file
             end if
             no_of_modes_new = modeind
          end do

       end if

       close(unit = fileUnit)
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Successfully read magnetic equilibrium from file ",trim(JGboozer_file)
       end if

       if (abs(normradius_old - normradius_wish) < abs(normradius_new - normradius_wish)) then
          iota = iota_old
          GHat = G_old
          IHat = I_old
          normradius = normradius_old
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
          normradius = normradius_new
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
       dGdpHat=(G_new-G_old)/(normradius_new*normradius_new-normradius_old*normradius_old)/pPrimeHat

       BHarmonics_amplitudes = BHarmonics_amplitudes / B0OverBBar

       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] This computation is for the flux surface with minor radius ",normradius*a, &
               " meters, equivalent to r/a = ",normradius
       end if

    case (12)
       ! Read Boozer coordinate file in a generalisation of the .bc format used at IPP Greifswald for non-stellarator symmetric equilibria 

       fileUnit = 11
       open(unit=fileUnit, file=JGboozer_file_NonStelSym, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file ",JGboozer_file_NonStelSym
          stop
       else
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             if (lineOfFile(1:2) /= "CC") exit
          end do

          ! Read header line:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) headerIntegers, headerReals
          if (didFileAccessWork /= 0) then
             print *,"Unable to read header from the magnetic equilibrium file ",JGboozer_file_NonStelSym
             stop
          end if

          NPeriods = headerIntegers(4)
          psiAHat  = headerReals(1)/2/pi !Convert the flux from Tm^2 to Tm^2/rad
          a        = headerReals(2)      !minor radius in meters

          end_of_file = .false.

          normradius_old = 0
          no_of_modes_old = 0
          modesm_old = 0
          modesn_old = 0
          modesb_old = 0
          iota_old = 0
          G_old = 0
          I_old = 0
          B0_old = 0
          pPrimeHat_old = 0

          normradius_new = 0
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
             if ((normradius_new .ge. normradius_wish) .or. end_of_file) exit

             normradius_old = normradius_new
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

             normradius_new = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
             iota_new = surfHeader(2)
             G_new = surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
             I_new = surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
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
                print *,"Error: no (0,0) mode found in magnetic equilibrium file ",JGboozer_file_NonStelSym
             else if (numB0s > 1) then
                print *,"Error: more than 1 (0,0) mode found in magnetic equilibrium file ",JGboozer_file_NonStelSym
             end if
             no_of_modes_new = modeind
          end do

       end if

       close(unit = fileUnit)
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Successfully read magnetic equilibrium from file ",trim(JGboozer_file_NonStelSym)
       end if

       if (abs(normradius_old - normradius_wish) < abs(normradius_new - normradius_wish)) then
          iota = iota_old
          GHat = G_old
          IHat = I_old
          normradius = normradius_old
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
          normradius = normradius_new
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
       dGdpHat=(G_new-G_old)/(normradius_new*normradius_new-normradius_old*normradius_old)/pPrimeHat

       allocate(BHarmonics_parity(NHarmonics))
       do i = 0, NHarmonics/2-1
          BHarmonics_parity(2*i+1)=.true.
          BHarmonics_parity(2*i+2)=.false.
       end do

       BHarmonics_amplitudes = BHarmonics_amplitudes / B0OverBBar

      if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] This computation is for the flux surface with minor radius ",normradius*a, &
               " meters, equivalent to r/a = ",normradius
      end if

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

  end subroutine computeBHat

end module geometry

