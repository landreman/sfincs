#include <finclude/petscsysdef.h>
subroutine validateInput()

  use globalVariables

  implicit none

  character(len=*), parameter :: line="******************************************************************"
  PetscScalar :: chargeDensity
  integer :: ispecies
  logical :: flag

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
  
  if (RHSMode == 2 .and. nonlinear) then
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

     if (nonlinear) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with nonlinear = .true., which is incompatble."
           print *,"**            Setting nonlinear = .false."
           print *,line
           print *,line
        end if
        nonlinear = .false.
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

     if (.not. useDKESExBDrift) then
        if (masterProc) then
           print *,line
           print *,line
           print *,"**   WARNING: You asked for RHSMode=3 (monoenergetic transport matrix) with useDKESExBDrift = .false., which is incompatble."
           print *,"**            Setting useDKESExBDrift = .true."
           print *,line
           print *,line
        end if
        useDKESExBDrift = .true.
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

  chargeDensity = zero
  do ispecies = 1,Nspecies
     chargeDensity = chargeDensity + nHats(ispecies)*Zs(ispecies)
  end do
  ! More needed here...  Ensure charge neutrality.

  ! physicsParameters namelist:

  if (constraintScheme .ne. -1 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You chose constraintScheme different from -1, which you should not do"
     print *,"**            unless you know what you are doing."
     print *,line
     print *,line
  end if

  if ((abs(alpha-one)>1d-15) .and. (abs(alpha-1000)>1d-15) .and. (abs(alpha-0.001)>1d-15) .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: Usually, either"
     print *,"**            alpha = 1 (if PhiBar = 1 V and TBar = 1 eV, or PhiBar = 1 kV and TBar = 1 keV)"
     print *,"**            or alpha = 1000 (if PhiBar = 1 kV and TBar = 1 eV)"
     print *,"**            or alpha = 0.001 (if PhiBar = 1 V and TBar = 1 keV)."
     print *,"**            Are you sure you want alpha = ",alpha,"?"
     print *,line
     print *,line
  end if

  if (nonlinear .and. (.not. includePhi1)) then
     if (masterProc) then
        print *,"Error! You requested a nonlinear calculation with includePhi1=.false."
        print *,"These options are inconsistent since the nonlinear terms involve Phi1."
     end if
     stop
  end if


  ! resolutionParameters namelist:

  if (Ntheta*Nzeta*Nx*Nxi*Nspecies > 1e7 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: You have selected large resolution parameters, leading to a matrix size"
     print *,"**            exceeding 10 million x 10 million.  SFINCS will almost certainly run out"
     print *,"**            of memory."
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

  if (Ntheta > 200 .and. masterProc) then
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
     elseif (Nx > 20) then
        print *,line
        print *,line
        print *,"**   WARNING: You chose a very large value for Nx."
        print *,"**            Typically Nx should be in the range 5-9."
        print *,line
        print *,line
     elseif ((Nx < 5) .or. (Nx > 9)) then
        print *,line
        print *,line
        print *,"**   WARNING: Typically Nx should be in the range 5-9."
        print *,line
        print *,line
     end if
  end if

  if (NL<1) then
     if (masterProc) then
        print *,"Error! NL must be positive."
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
     elseif (solverTolerance > 1e-5) then
        print *,line
        print *,line
        print *,"**   WARNING: You selected a large solverTolerance, which may cause the Krylov/KSP iteration"
        print *,"**            to stop before an accurate solution is obtained."
        print *,"**            Good values for solverTolerance are typically 1e-5 to 1e-7."
        print *,line
        print *,line
     end if
  end if

  if ((.not.forceOddNthetaAndNzeta ) .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: forceOddNthetaAndNzeta = .true. is strongly recommended."
     print *,line
     print *,line
  end if

  ! otherNumericalParameters namelist:

  if (thetaDerivativeScheme<0) then
     if (masterProc) then
        print *,"Error! thetaDerivativeScheme cannot be less than 0."
     end if
     stop
  end if
  
  if (thetaDerivativeScheme>2) then
     if (masterProc) then
        print *,"Error! thetaDerivativeScheme cannot be more than 2."
     end if
     stop
  end if
  
  if (thetaDerivativeScheme == 0 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: thetaDerivativeScheme=0 leads to very dense matrices,"
     print *,"**            meaning a lot of time and memory is required to solve the system."
     print *,"**            thetaDerivativeScheme = 2 is strongly recommended."
     print *,line
     print *,line
  end if

  if (thetaDerivativeScheme == 1 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: thetaDerivativeScheme=2 is typically preferred to thetaDerivativeScheme=1"
     print *,"**            since accuracy is higher for relatively little additional computational cost."
     print *,line
     print *,line
  end if

  if (zetaDerivativeScheme<0) then
     if (masterProc) then
        print *,"Error! zetaDerivativeScheme cannot be less than 0."
     end if
     stop
  end if
  
  if (zetaDerivativeScheme>2) then
     if (masterProc) then
        print *,"Error! zetaDerivativeScheme cannot be more than 2."
     end if
     stop
  end if
  
  if (zetaDerivativeScheme == 0 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: zetaDerivativeScheme=0 leads to very dense matrices,"
     print *,"**            meaning a lot of time and memory is required to solve the system."
     print *,"**            zetaDerivativeScheme = 2 is strongly recommended."
     print *,line
     print *,line
  end if

  if (zetaDerivativeScheme == 1 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: zetaDerivativeScheme=2 is typically preferred to zetaDerivativeScheme=1"
     print *,"**            since accuracy is higher for relatively little additional computational cost."
     print *,line
     print *,line
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

  if (preconditioner_x_min_L<0) then
     if (masterProc) then
        print *,"Error! preconditioner_x_min_L should not be less than 0."
     end if
     stop
  end if
  
  if (preconditioner_x_min_L > 2 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: preconditioner_x_min_L should typically be 0, 1, or 2."
     print *,line
     print *,line
  end if

  if (preconditioner_theta<0) then
     if (masterProc) then
        print *,"Error! preconditioner_theta cannot be less than 0."
     end if
     stop
  end if
  
  if (preconditioner_theta>1) then
     if (masterProc) then
        print *,"Error! preconditioner_theta cannot be more than 1."
     end if
     stop
  end if
  
  if (preconditioner_theta == 1 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: preconditioner_theta = 1 often does not work well (i.e. GMRES/KSP does not converge rapidly.)"
     print *,"**            preconditioner_theta = 0 is strongly recommended."
     print *,line
     print *,line
  end if

  if (preconditioner_zeta<0) then
     if (masterProc) then
        print *,"Error! preconditioner_zeta cannot be less than 0."
     end if
     stop
  end if
  
  if (preconditioner_zeta>1) then
     if (masterProc) then
        print *,"Error! preconditioner_zeta cannot be more than 1."
     end if
     stop
  end if
  
  if (preconditioner_zeta == 1 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: preconditioner_zeta = 1 often does not work well (i.e. GMRES/KSP does not converge rapidly.)"
     print *,"**            preconditioner_zeta = 0 is strongly recommended."
     print *,line
     print *,line
  end if

  if (preconditioner_xi<0) then
     if (masterProc) then
        print *,"Error! preconditioner_xi cannot be less than 0."
     end if
     stop
  end if
  
  if (preconditioner_xi>1) then
     if (masterProc) then
        print *,"Error! preconditioner_xi cannot be more than 1."
     end if
     stop
  end if
  
end subroutine validateInput

