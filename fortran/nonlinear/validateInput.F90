#include <finclude/petscsysdef.h>
subroutine validateInput()

  use globalVariables

  implicit none

  character(len=*), parameter :: line="******************************************************************"
  PetscScalar :: chargeDensity
  integer :: ispecies

  ! General namelist

  if (RHSMode<1) then
     if (masterProc) then
        print *,"Error! RHSMode must be at least 1."
     end if
     stop
  end if
  
  if (RHSMode>2) then
     if (masterProc) then
        print *,"Error! RHSMode must be no more than 2."
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

  chargeDensity = zero
  do ispecies = 1,Nspecies
     chargeDensity = chargeDensity + nHats(ispecies)*Zs(ispecies)
  end do
  ! More needed here...

  ! physicsParameters namelist:
  if (nonlinear .and. (.not. includePhi1)) then
     if (masterProc) then
        print *,"Error! You requested a nonlinear calculation with includePhi1=.false."
        print *,"These options are inconsistent since the nonlinear terms involve Phi1."
     end if
     stop
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

