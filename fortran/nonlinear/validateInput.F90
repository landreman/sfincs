#include <finclude/petscsysdef.h>
subroutine validateInput()

  use globalVariables

  implicit none

  character(len=*), parameter :: line="******************************************************************"


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

  ! preconditionerOptions namelist:

  if (preconditioner_x .ne. 1 .and. masterProc) then
     print *,line
     print *,line
     print *,"**   WARNING: preconditioner_x should typically be 1."
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

end subroutine validateInput

