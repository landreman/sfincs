#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

subroutine validateInput()

  use kinds
  use globalVariables
  use xGrid, only: xGrid_k
  use petscksp

  implicit none

  character(len=*), parameter :: line="******************************************************************"
  real(prec) :: chargeDensity
  real(prec) :: maxSingleChargeDensity !!Added by AM 2016-03
  integer :: ispecies
  logical :: flag
  real(prec) :: lnLambda, eC, epsilon0, mproton

  ! General namelist

  if (RHSMode<1) then
     if (masterProc) then
        print *,"Error! RHSMode must be at least 1."
     end if
     stop
  end if
  
  if (RHSMode>3) then
     if (masterProc) then
        print *,"Error! RHSMode must be no more than 3."
     end if
     stop
  end if
  
  !!if (RHSMode == 2 .and. nonlinear) then !!Commented by AM 2016-02
  if (RHSMode == 2 .and. includePhi1) then !!Added by AM 2016-02
     if (masterProc) then
        print *,"Error! RHSMode cannot be 2 for a nonlinear calculation."
     end if
     stop
  end if
  
  if (RHSMode == 2 .and. Nspecies>1) then
     if (masterProc) then
        print *,"Error! The transport matrix is presently only available in SFINCS for a 1-species calculation."
     end if
     stop
  end if

  if (RHSMode == 3) then
     ! Computing monoenergetic transport coefficients.
     ! Make sure the code is configured to use the DKES form of the kinetic equation.

     if (Nxi_for_x_option .ne. 0) then 
        if (masterProc) then
           print *,"Setting Nxi_for_x_option=0, since RHSMode=3."
        end if
        Nxi_for_x_option=0
     end if

     !!if (nonlinear) then !!Commented by AM 2016-02
     if (includePhi1) then !!Added by AM 2016-02
        if (masterProc) then
           print *,line
           print *,line
!!           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with nonlinear = .true., which is incompatble." !!Commented by AM 2016-02
!!           print *,"**            Setting nonlinear = .false." !!Commented by AM 2016-02
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with includePhi1 = .true., which is incompatble." !!Added by AM 2016-02
           print *,"**            Setting includePhi1 = .false." !!Added by AM 2016-02
           print *,line
           print *,line
        end if
!!        nonlinear = .false. !!Commented by AM 2016-02
        includePhi1 = .false. !!Added by AM 2016-02
     end if

     if (Nx > 1) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with Nx > 1, which is incompatble."
           print *,"**            Setting Nx = 1."
           print *,line
           print *,line
        end if
        Nx = 1
     end if

     if (ExB_option .ne. 2) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with ExB_option != 2, which is incompatble."
           print *,"**            Setting ExB_option=2."
           print *,line
           print *,line
        end if
        ExB_option = 2
     end if

     if (includeXDotTerm) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with includeXDotTerm = .true., which is incompatble."
           print *,"**            Setting includeXDotTerm = .false."
           print *,line
           print *,line
        end if
        includeXDotTerm = .false.
     end if

     if (includeElectricFieldTermInXiDot) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with includeElectricFieldTermInXiDot = .true., which is incompatble."
           print *,"**            Setting includeElectricFieldTermInXiDot = .false."
           print *,line
           print *,line
        end if
        includeElectricFieldTermInXiDot = .false.
     end if

     if (NSpecies > 1) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with >1 species, which is incompatble."
           print *,"**            Ignoring all species after the first."
           print *,line
           print *,line
        end if
        Nspecies = 1
     end if

     if (includePhi1) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with includePhi1 = .true., which is incompatble."
           print *,"**            Setting includePhi1 = .false."
           print *,line
           print *,line
        end if
        includePhi1 = .false.
     end if

     if (collisionOperator .ne. 1) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with collisionOperator .ne. 1, which is incompatble."
           print *,"**            Setting collisionOperator = 1."
           print *,line
           print *,line
        end if
        collisionOperator = 1
     end if

     if (includeTemperatureEquilibrationTerm) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with includeTemperatureEquilibrationTerm = .true., which is incompatble."
           print *,"**            Setting temperatureEquilibrationTerm = .false."
           print *,line
           print *,line
        end if
        includeTemperatureEquilibrationTerm = .false.
     end if

     mHats = 1
     nHats = 1
     THats = 1
     dnHatdpsiHats = 0
     dTHatdpsiHats = 0
     Zs = 1
     if (masterProc) then
        print *,"Since RHSMode=3, ignoring the requested values of Zs, nHats, THats, nu_n, Er, and dPhiHatd*."
     end if

     if (abs(nuPrime) < 1e-14) then
        if (masterProc) then
           print *,"Error! When running with RHSMode=3, you must set nuPrime to a nonzero value."
        end if
        stop
     end if

  end if

  
  if (saveMatlabOutput .and. Nspecies*Ntheta*Nzeta*Nxi*Nx > 5000 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You asked to save matlab-format ASCII files for a large matrix size."
     print *,"**            This may take a long time and result in large files."
     print *,line
     print *,line
  end if


  ! geometryParameters namelist:

  if (min_Bmn_to_load > 0.01 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: min_Bmn_to_load = ",min_Bmn_to_load
     print *,"              Are you sure you want min_Bmn_to_load to be that large?"
     print *,line
     print *,line
  end if

  ! species namelist:

  flag = .false.
  do ispecies = 1,Nspecies
     if (Zs(ispecies) < 0) then
        if (flag .and. masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: More than 1 species has negative charge, which is unusual."
           print *,line
           print *,line
        end if
        flag = .true.
     end if

     if (abs(Zs(ispecies) - floor(Zs(ispecies))) > 0 .and. masterProc) then
        print *,line
        print *,line
        print *,"**   WARNING: At least one of the charges Zs is not an integer, which is unusual."
        print *,line
        print *,line
     end if

     if (Zs(ispecies) == 0) then
        if (masterProc) then
           print *,"Error! Charges Zs cannot be zero."
        end if
        stop
     end if

     if (mHats(ispecies) .le. 0) then
        if (masterProc) then
           print *,"Error! Masses mHats must be positive."
        end if
        stop
     end if

     if (THats(ispecies) .le. 0) then
        if (masterProc) then
           print *,"Error! Temperatures THats must be positive."
        end if
        stop
     end if

     if (nHats(ispecies) .le. 0) then
        if (masterProc) then
           print *,"Error! Densities nHats must be positive."
        end if
        stop
     end if
  end do


  !!!!!!!!!!!!!!!!!!!!!!!
  !!Added by AM 2015-11!!
  if (withAdiabatic) then
     if (adiabaticZ == 0) then
     	if (masterProc) then
           print *,"Error! Charge adiabaticZ cannot be zero."
        end if
        stop
     end if

     if (adiabaticMHat .le. 0) then
     	if (masterProc) then
   	   print *,"Error! Mass adiabaticMHat must be positive."
	end if
        stop
     end if

     if (adiabaticNHat .le. 0) then
     	if (masterProc) then
     	   print *,"Error! Density adiabaticNHat must be positive."
	end if
	stop
     end if

     if (adiabaticTHat .le. 0) then
     	if (masterProc) then
           print *,"Error! Temperature adiabaticTHat must be positive."
	end if
        stop
     end if

     if (adiabaticZ < 0) then
        if (flag .and. masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: More than 1 species has negative charge, which is unusual."
           print *,line
           print *,line
        end if
        flag = .true.
     end if

     if (abs(adiabaticZ - floor(adiabaticZ)) > 0 .and. masterProc) then
        print *,line
        print *,line
        print *,"**   WARNING: adiabaticZ is not an integer, which is unusual."
        print *,line
        print *,line
     end if

 
  end if
  !!!!!!!!!!!!!!!!!!!!!!!

  ! Ensure charge neutrality.
  chargeDensity = zero
  maxSingleChargeDensity = 1d-11 !!Added by AM 2016-03
  do ispecies = 1,Nspecies
     chargeDensity = chargeDensity + nHats(ispecies)*Zs(ispecies)
     maxSingleChargeDensity = max(abs(maxSingleChargeDensity), abs(nHats(ispecies)*Zs(ispecies))) !!Added by AM 2016-03
     !!Added by AM 2016-02!!
     if (quasineutralityOption == 2 .and. includePhi1) then
        exit !!If running with EUTERPE equations we only use the first kinetic species in quasi-neutrality.
     end if
     !!!!!!!!!!!!!!!!!!!!!!!
  end do
 

  !!!!!!!!!!!!!!!!!!!!!!!
  !!Added by AM 2015-11!!
  if (withAdiabatic) then
     chargeDensity = chargeDensity + adiabaticNHat*adiabaticZ
     maxSingleChargeDensity = max(abs(maxSingleChargeDensity), abs(adiabaticNHat*adiabaticZ))
  end if

!!  if (includePhi1 .and. (abs(chargeDensity) >1d-15)) then
  if (includePhi1 .and. (abs(chargeDensity)/maxSingleChargeDensity >1d-4)) then
     if (masterProc) then
        print *,"Error! When running with includePhi1=.true. you must ensure that"
        print *,"quasi-neutrality is fulfilled for the input species."
        if (withAdiabatic) then
           print *,"Be aware that also the adiabatic species is included in quasi-neutrality."
           if (quasineutralityOption == 2) then
              print *,"If running with EUTERPE quasi-neutrality equations (quasineutralityOption = 2) only the first kinetic species is used in quasi-neutrality."
           end if
        end if
     end if
     stop
  end if 

!!  if (abs(chargeDensity) > 1d-15 .and. (Nspecies > 1 .or. withAdiabatic) .and. masterProc) then
  if (abs(chargeDensity)/maxSingleChargeDensity > 1d-4 .and. (Nspecies > 1 .or. withAdiabatic) .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: Your input does not fulfill quasi-neutrality,"
     print *,"**            which it typically should when using several species."
     print *,line
     print *,line
  end if
  !!!!!!!!!!!!!!!!!!!!!!!


  ! physicsParameters namelist:

  if (nu_n<0) then
     ! Compute nu_n assuming the normalizations are RBar = 1m, mBar = proton mass, nBar = 1e20/m^3, TBar = 1keV.
     if (abs(Delta-(4.5694e-3)) / (4.5694e-3) > 0.1) then
        print *,line
        print *,line
        print *,"**   WARNING: nu_n will be computed assuming RBar = 1m, mBar = proton mass, nBar = 1e20/m^3, TBar = 1keV."
        print *,"**            However it appears you have set Delta using different normalizations."
        print *,line
        print *,line
     end if
     if (abs(Zs(1)+1) > 1e-6) then
        print *,line
        print *,line
        print *,"**   WARNING: lnLambda will be computed using the density and temperature of the first species."
        print *,"**            However, based on Zs, it appears the first species is not electrons."
        print *,line
        print *,line
     end if

     ! Use the same definition of lnLambda as Yuriy Turkin et al,
     ! assuming species(1) is the electrons:
     lnLambda = (25.3d+0) - (1.15d+0)*log10(nHats(1)*(1e14)) + (2.30d+0)*log10(THats(1)*1000)
     !lnLambda = 17

     eC = 1.6022d-19
     epsilon0 = 8.8542d-12
     mproton = 1.6726d-27
     nu_n = sqrt(mproton/(2*1000*eC)) * 4*sqrt(2*pi)*(1e20)*(eC**4)*lnLambda / (3*((4*pi*epsilon0)**2)*sqrt(mproton)*((1000*eC)**(1.5d+0)))

     if (masterProc) then
        print *,"Computing nu_n using ln(Lambda)=",lnLambda
        print *,"New nu_n:",nu_n
     end if
  end if

  if (constraintScheme .ne. -1 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You chose constraintScheme different from -1, which you should not do"
     print *,"**            unless you know what you are doing."
     print *,line
     print *,line
  end if

  if ((abs(gamma-one)>1d-15) .and. (abs(gamma-1000)>1d-15) .and. (abs(gamma-0.001)>1d-15) .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: Usually, either"
     print *,"**            gamma = 1 (if PhiBar = 1 V and TBar = 1 eV, or PhiBar = 1 kV and TBar = 1 keV)"
     print *,"**            or gamma = 1000 (if PhiBar = 1 kV and TBar = 1 eV)"
     print *,"**            or gamma = 0.001 (if PhiBar = 1 V and TBar = 1 keV)."
     print *,"**            Are you sure you want gamma = ",gamma,"?"
     print *,line
     print *,line
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Region commented by AM 2016-02!!
!!$  if (nonlinear .and. (.not. includePhi1)) then
!!$     if (masterProc) then
!!$        print *,"Error! You requested a nonlinear calculation with includePhi1=.false."
!!$        print *,"These options are inconsistent since the nonlinear terms involve Phi1."
!!$     end if
!!$     stop
!!$  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! By AM 2016-03: I decided to remove this check below because includePhi1InKineticEquation does not matter if includePhi1=.false.

!!$!!  if (includeRadialExBDrive .and. (.not. includePhi1)) then !!Commented by AM 2016-03
!!$  if (includePhi1InKineticEquation .and. (.not. includePhi1)) then !!Added by AM 2016-03
!!$     if (masterProc) then
!!$!!        print *,"Error! You requested a calculation including the radial ExB drive term" !!Commented by AM 2016-03
!!$!!        print *,"(includeRadialExBDrive=.true.) but you set includePhi1=.false." !!Commented by AM 2016-03
!!$!!        print *,"These options are inconsistent since the radial ExB drive term involves Phi1." !!Commented by AM 2016-03
!!$        print *,"Error! You requested a calculation including Phi1 in the kinetic equation" !!Added by AM 2016-03
!!$        print *,"(includePhi1InKineticEquation=.true.) but you set includePhi1=.false." !!Added by AM 2016-03
!!$        print *,"These options are inconsistent." !!Added by AM 2016-03
!!$     end if
!!$     stop
!!$  end if

  if (magneticDriftScheme<0) then
     if (masterProc) then
        print *,"Error! magneticDriftScheme must be >= 0."
     end if
     stop
  end if

  if (magneticDriftScheme>5) then
     if (masterProc) then
        print *,"Error! magneticDriftScheme must be <= 5."
     end if
     stop
  end if

  if (magneticDriftScheme>3) then
     if (.not.( (geometryScheme == 11) .or. (geometryScheme == 12) )) then
        if (masterProc) then
           print *,"Error! magneticDriftScheme 4 and 5 have only been implemented for geometryScheme 11 and 12."
        end if
        stop        
     end if
  end if

  if (magneticDriftScheme>0 .and. includePhi1) then
     if (masterProc) then
        print *,"**   ERROR! Some terms involving Phi1 and the magnetic drifts have not yet been implemented."
        print *,"**          Hence magneticDriftScheme>0 is incompatible with includePhi1."
     end if
     stop
  end if

  if (magneticDriftScheme>0) then
     select case (geometryScheme)
        case (5,6,7,11,12)
           ! No problem, magnetic drifts have been implemented for these geometries.
        case default
           if (masterProc) then
              print *,"Error! You requested that poloidal/toroidal magnetic drifts be included (magneticDriftScheme>0)"
              print *,"       but you selected a geometryScheme for which the required components of the magnetic field"
              print *,"       are not available."
           end if
           stop
        end select
  end if

  !!!!!!!!!!!!!!!!!!!!!!!
  !!Added by AM 2016-02!!
  if ((quasineutralityOption < 1 .or. quasineutralityOption > 2) .and. includePhi1) then
      if (masterProc) then
        print *,"Error! quasineutralityOption must be 1 or 2."
      end if
      stop
  end if

  if (withAdiabatic .and. (.not. includePhi1) .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You are running with an adiabatic species,"
     print *,"**            but with includePhi1 = .false. which implies "
     print *,"**            that the adiabatic species has no impact."
     print *,line
     print *,line
  end if

  if ((.not. withAdiabatic) .and. includePhi1 .and. quasineutralityOption == 1 .and. (Nspecies < 2)) then 
      if (masterProc) then
        print *,"Error! In a nonlinear run (includePhi1 = .true.) you must use at least two species, to be able to fulfill quasi-neutrality."
      end if
      stop
  end if 

  if ((.not. withAdiabatic) .and. includePhi1 .and. quasineutralityOption == 2) then 
      if (masterProc) then
        print *,"Error! If running with EUTERPE quasi-neutrality equations (quasineutralityOption = 2) in a nonlinear run (includePhi1 = .true.) you must use an adiabatic species (withAdiabatic = .true.)."
      end if
      stop
  end if 

!!$  if ((Nspecies .ne. 1) .and. includePhi1 .and. quasineutralityOption == 2) then 
!!$      if (masterProc) then
!!$        print *,"Error! If running with EUTERPE quasi-neutrality equations (quasineutralityOption = 2) in a nonlinear run (includePhi1 = .true.) you must use ONE kinetic species only."
!!$      end if
!!$      stop
!!$  end if 

  !!!!!!!!!!!!!!!!!!!!!!!


  ! resolutionParameters namelist:

  if (Ntheta*Nzeta*Nx*Nxi*Nspecies > 1e7 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You have selected large resolution parameters, leading to a matrix size"
     print *,"**            exceeding 10 million x 10 million.  These parameters will be numerically challenging."
     print *,line
     print *,line
  end if

  if (Ntheta<5) then
     if (masterProc) then
        print *,"Error! Ntheta must be at least 5."
     end if
     stop
  end if

  if (((Ntheta > 100 .and. Nzeta > 1) .or. (Ntheta>250)) .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You chose a very large value for Ntheta."
     print *,line
     print *,line
  end if

  if (Nzeta<1) then
     if (masterProc) then
        print *,"Error! Nzeta must be positive."
     end if
     stop
  end if

  if (Nzeta > 200 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You chose a very large value for Nzeta."
     print *,line
     print *,line
  end if

  if (Nxi<1) then
     if (masterProc) then
        print *,"Error! Nxi must be positive."
     end if
     stop
  end if

  if (Nxi < 4 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You probably should have Nxi at least 4."
     print *,line
     print *,line
  end if

  if (Nxi > 200 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You chose a very large value for Nxi."
     print *,line
     print *,line
  end if

  if (Nx<1) then
     if (masterProc) then
        print *,"Error! Nx must be positive."
     end if
     stop
  end if

  if (masterProc .and. (RHSMode .ne. 3)) then
     if (Nx < 4) then
        print *,line
        print *,line
        print *,"**   WARNING: You almost certainly should have Nx at least 4."
        print *,"              (The exception is when RHSMode = 3, in which case Nx = 1.)"
        print *,line
        print *,line
     elseif ((Nx > 20) .and. (xGridScheme==1 .or. xGridScheme==2)) then
        print *,line
        print *,line
        print *,"**   WARNING: You chose a very large value for Nx."
        print *,"**            For xGridMode=1 or 2, typically Nx should be in the range 5-9."
        print *,line
        print *,line
     elseif (((Nx < 5) .or. (Nx > 9)) .and. (xGridScheme==1 .or. xGridScheme==2)) then
        print *,line
        print *,line
        print *,"**   WARNING: For xGridMode=1 or 2, typically Nx should be in the range 5-9."
        print *,line
        print *,line
     end if
  end if

  if (NL<0) then
     if (masterProc) then
        print *,"Error! NL must be at least 0."
     end if
     stop
  end if

  if (masterProc) then
     if (NL < 2) then
        print *,line
        print *,line
        print *,"**   WARNING: You probably should have NL at least 2."
        print *,"**            A value which almost always works is NL = 4."
        print *,line
        print *,line
     elseif (NL > 8) then
        print *,line
        print *,line
        print *,"**   WARNING: You chose a very large value for NL."
        print *,"**            A value which almost always works is NL = 4."
        print *,line
        print *,line
     elseif (NL .ne. 4) then
        print *,line
        print *,line
        print *,"**   WARNING: Although values of NL in the range [2,8] work well, NL = 4 is usually recommended."
        print *,line
        print *,line
     end if
  end if

  if (NxPotentialsPerVth <= 0) then
     if (masterProc) then
        print *,"Error! NxPotentialsPerVth must be positive."
     end if
     stop
  end if

  if (masterProc) then
     if (NxPotentialsPerVth < 10) then
        print *,line
        print *,line
        print *,"**   WARNING: You probably should have NxPotentialsPerVth at least 10."
        print *,"**            A value which almost always works is NxPotentialsPerVth = 40."
        print *,line
        print *,line
     elseif (NxPotentialsPerVth > 100) then
        print *,line
        print *,line
        print *,"**   WARNING: You chose a very large value for NxPotentialsPerVth."
        print *,"**            A value which almost always works is NxPotentialsPerVth = 40."
        print *,line
        print *,line
     elseif (abs(NxPotentialsPerVth - 40) > 1d-15) then
        print *,line
        print *,line
        print *,"**   WARNING: Although values of NxPotentialsPerVth in the range [10,100] work well,"
        print *,"**            NxPotentialsPerVth = 40 is usually recommended."
        print *,line
        print *,line
     end if
  end if

  if (xMax <= 0) then
     if (masterProc) then
        print *,"Error! xMax must be positive."
     end if
     stop
  end if

  if (xMax < 2 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You probably should have xMax at least 2."
     print *,"**            A value which almost always works is xMax = 5."
     print *,line
     print *,line
  end if

  if (xMax > 10 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You chose a very large value for xMax."
     print *,"**            A value which almost always works is xMax = 5."
     print *,line
     print *,line
  end if

  if (solverTolerance < 0) then
     if (masterProc) then
        print *,"Error! solverTolerance must be positive."
     end if
     stop
  end if

  if (masterProc) then
     if (solverTolerance < 1d-10) then
        print *,line
        print *,line
        print *,"**   WARNING: You selected a very small solverTolerance which may be hard for the solver to achieve."
        print *,"**            Good values for solverTolerance are typically 1e-5 to 1e-7."
        print *,line
        print *,line
     elseif (solverTolerance < 1d-7) then
        print *,line
        print *,line
        print *,"**   WARNING: You selected a small solverTolerance which may require more Krylov/KSP iterations than necessary."
        print *,"**            Good values for solverTolerance are typically 1e-5 to 1e-7."
        print *,line
        print *,line
     elseif (solverTolerance > 2e-5) then
        print *,line
        print *,line
        print *,"**   WARNING: You selected a large solverTolerance, which may cause the Krylov/KSP iteration"
        print *,"**            to stop before an accurate solution is obtained."
        print *,"**            Good values for solverTolerance are typically 1e-5 to 1e-7."
        print *,line
        print *,line
     end if
  end if

  ! otherNumericalParameters namelist:
  
  if (xGridScheme<1) then
     if (masterProc) then
        print *,"Error! xGridScheme cannot be less than 1."
     end if
     stop
  end if
  
  if (xGridScheme>8) then
     if (masterProc) then
        print *,"Error! xGridScheme cannot be more than 8."
     end if
     stop
  end if
  
  if (xPotentialsGridScheme<1) then
     if (masterProc) then
        print *,"Error! xPotentialsGridScheme cannot be less than 1."
     end if
     stop
  end if
  
  if (xPotentialsGridScheme>4) then
     if (masterProc) then
        print *,"Error! xPotentialsGridScheme cannot be more than 4."
     end if
     stop
  end if
  
  if ((xPotentialsGridScheme==3 .or. xPotentialsGridScheme==4) .and. (xGridScheme .ne. 3 .and. xGridScheme .ne. 4)) then
     if (masterProc) then
        print *,"Error! When xPotentialsGridScheme is 3 or 4, xGridScheme must be 3 or 4."
     end if
     stop
  end if
  
  if ((xGridScheme==2 .or. xGridScheme==6) .and. (xGrid_k .ne. 0)) then
     if (masterProc) then
        print *,line
        print *,line
        print *,"** WARNING: Overriding your request for xGrid_k, since"
        print *,"** for xGridScheme of 2 or 6, you must set xGrid_k to 0."
        print *,line
        print *,line
     end if
     xGrid_k = 0
  end if

  ! preconditionerOptions namelist:

  if (preconditioner_species<0) then
     if (masterProc) then
        print *,"Error! preconditioner_species should not be less than 0."
     end if
     stop
  end if
  
  if (preconditioner_species>1) then
     if (masterProc) then
        print *,"Error! preconditioner_species should not be more than 1."
     end if
     stop
  end if
  
  if (preconditioner_x<0) then
     if (masterProc) then
        print *,"Error! preconditioner_x should not be less than 0."
     end if
     stop
  end if
  
  if (preconditioner_x>4) then
     if (masterProc) then
        print *,"Error! preconditioner_x should not be more than 4."
     end if
     stop
  end if
  
  if (preconditioner_x .ne. 1 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: preconditioner_x = 1 is usually the best option."
     print *,line
     print *,line
  end if

#ifdef PETSC_HAVE_HYPRE
  if (masterProc) print *,"Hypre detected."
#else
  if (masterProc) print *,"Hypre not detected."
  if (preconditioning_option==3 .or. preconditioning_option==6) then
     preconditioning_option=5
     print *,line
     print *,line
     print *,"** WARNING: Hypre preconditioner was requested (preconditioning_option=3 or 6),"
     print *,"** but Hypre is not available in this PETSc library. Switching to the GAMG"
     print *,"** preconditioner (preconditioning_option=5)."
     print *,line
     print *,line
  end if
#endif
  
#ifdef PETSC_HAVE_ML
  if (masterProc) print *,"ML detected."
#else
  if (masterProc) print *,"ML not detected."
  if (preconditioning_option==8 .or. preconditioning_option==9) then
     preconditioning_option=5
     print *,line
     print *,line
     print *,"** WARNING: ML preconditioner was requested (preconditioning_option=8 or 9),"
     print *,"** but ML is not available in this PETSc library. Switching to the GAMG"
     print *,"** preconditioner (preconditioning_option=5)."
     print *,line
     print *,line
  end if
#endif
  
end subroutine validateInput

