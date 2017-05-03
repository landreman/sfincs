#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

subroutine set_grid_resolutions()

  use petscsys
  use variables

  implicit none

  PetscScalar :: temp_float
  PetscInt :: temp_int, j
  PetscScalar, parameter :: small = 1.0d-12
  integer, parameter :: max_N_levels = 20
  PetscScalar :: Ntheta_candidates(max_N_levels), Nzeta_candidates(max_N_levels), Nxi_candidates(max_N_levels)
  PetscScalar :: theta_coupling, zeta_coupling, xi_coupling, factor, new_max_coupling

  N_levels = 1
  Ntheta_min = min(Ntheta_min,Ntheta)
  Nzeta_min  = min(Nzeta_min ,Nzeta)
  Nxi_min    = min(Nxi_min   ,Nxi)
!  if (coarsen_theta) N_levels = max(N_levels,ceiling(log(Ntheta*one/Ntheta_min) / log(two) - small) +1)
!  if (coarsen_zeta ) N_levels = max(N_levels,ceiling(log(Nzeta *one/ Nzeta_min) / log(two) - small) +1)
!  if (coarsen_xi   ) N_levels = max(N_levels,ceiling(log((Nxi-one)/ (Nxi_min-one)) / log(two) - small) +1)
  if (coarsen_theta) N_levels = max(N_levels,nint(log(Ntheta*one/Ntheta_min) / log(two) - small) +1)
  if (coarsen_zeta ) N_levels = max(N_levels,nint(log(Nzeta *one/ Nzeta_min) / log(two) - small) +1)
  if (coarsen_xi   ) N_levels = max(N_levels,nint(log((Nxi-one)/ (Nxi_min-one)) / log(two) - small) +1)

  if (coarsen_option==3) then
     ! Smart coarsening, in theta & zeta but not xi.
     coarsen_theta = .true.
     coarsen_zeta = .true.
     coarsen_xi = .false.
     ! Coupling is like d/dx, where dx~1/N, so coupling is like N.
     theta_coupling = iota * Ntheta
     zeta_coupling = Nzeta * Nperiods
     Ntheta_candidates(1) = Ntheta
     Nzeta_candidates(1) = Nzeta
     do j=2,max_N_levels
        if (masterProc) print *,"theta coupling:",theta_coupling,"zeta coupling:",zeta_coupling
        new_max_coupling = max(theta_coupling,zeta_coupling)/2

        factor = min(1.0d+0,new_max_coupling/theta_coupling)
        Ntheta_candidates(j) = Ntheta_candidates(j-1) * factor
        theta_coupling = theta_coupling * factor

        factor = min(1.0d+0,new_max_coupling/zeta_coupling)
        Nzeta_candidates(j) = Nzeta_candidates(j-1) * factor
        zeta_coupling = zeta_coupling * factor

!!$        if (theta_coupling > 1.5 * zeta_coupling) then
!!$           if (masterProc) print *,"theta_coupling > 1.5 * zeta_coupling"
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)/2
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)
!!$           theta_coupling = theta_coupling / 2
!!$        elseif (zeta_coupling > 1.5 * theta_coupling) then
!!$           if (masterProc) print *,"zeta_coupling > 1.5 * theta_coupling"
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)/2
!!$           zeta_coupling = zeta_coupling / 2
!!$        elseif (theta_coupling > zeta_coupling) then
!!$           if (masterProc) print *,"theta_coupling > zeta_coupling"
!!$           factor = theta_coupling / zeta_coupling ! factor is between 1 and 1.5
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)/2
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)*factor/2
!!$           theta_coupling = theta_coupling / 2
!!$           zeta_coupling = zeta_coupling*factor/2
!!$        elseif (zeta_coupling > theta_coupling) then
!!$           if (masterProc) print *,"zeta_coupling > theta_coupling"
!!$           factor = zeta_coupling / theta_coupling ! factor is between 1 and 1.5
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)/2
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)*factor/2
!!$           zeta_coupling = zeta_coupling / 2
!!$           theta_coupling = theta_coupling*factor/2
!!$        end if
        if (Ntheta_candidates(j)<=Ntheta_min) then
           Ntheta_candidates(j) = Ntheta_min
           exit
        end if
        if (Nzeta_candidates(j)<=Nzeta_min) then
           Nzeta_candidates(j) = Nzeta_min
           exit
        end if
     end do
     N_levels = j
  elseif (coarsen_option==4) then
     print *,"Nxi:",Nxi
     ! Smart coarsening, in all 3 coordinates.
     coarsen_theta = .true.
     coarsen_zeta = .true.
     coarsen_xi = .true.
     ! Coupling is like d/dx, where dx~1/N, so coupling is like N.
     theta_coupling = iota * Ntheta / (2*pi)
     zeta_coupling = Nzeta * Nperiods / (2*pi)
     xi_coupling = nu * Nxi / 2
     Ntheta_candidates(1) = Ntheta
     Nzeta_candidates(1) = Nzeta
     Nxi_candidates(1) = Nxi
     do j=2,max_N_levels
        if (masterProc) print *,"Coupling in theta:",theta_coupling,"zeta:",zeta_coupling,"xi:",xi_coupling
        new_max_coupling = max(theta_coupling,zeta_coupling,xi_coupling)/2

        factor = min(1.0d+0,new_max_coupling/theta_coupling)
        Ntheta_candidates(j) = Ntheta_candidates(j-1) * factor
        theta_coupling = theta_coupling * factor

        factor = min(1.0d+0,new_max_coupling/zeta_coupling)
        Nzeta_candidates(j) = Nzeta_candidates(j-1) * factor
        zeta_coupling = zeta_coupling * factor

        factor = min(1.0d+0,new_max_coupling/xi_coupling)
        Nxi_candidates(j) = Nxi_candidates(j-1) * factor
        xi_coupling = xi_coupling * factor

!!$        if (theta_coupling > 1.5 * zeta_coupling .and. theta_coupling > 1.5 * xi_coupling) then
!!$           if (masterProc) print *,"theta_coupling > 1.5 * zeta,xi_coupling"
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)/2
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)
!!$           Nxi_candidates(j) = Nxi_candidates(j-1)
!!$           theta_coupling = theta_coupling / 2
!!$        elseif (zeta_coupling > 1.5 * theta_coupling .and. zeta_coupling > 1.5 * xi_coupling) then
!!$           if (masterProc) print *,"zeta_coupling > 1.5 * theta,xi_coupling"
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)/2
!!$           Nxi_candidates(j) = Nxi_candidates(j-1)
!!$           zeta_coupling = zeta_coupling / 2
!!$        elseif (xi_coupling > 1.5 * theta_coupling .and. xi_coupling > 1.5 * zeta_coupling) then
!!$           if (masterProc) print *,"xi_coupling > 1.5 * theta,zeta_coupling"
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)
!!$           Nxi_candidates(j) = Nxi_candidates(j-1)/2
!!$           xi_coupling = xi_coupling / 2
!!$        elseif (theta_coupling > zeta_coupling .and. theta_coupling > xi_coupling) then
!!$           if (masterProc) print *,"theta_coupling > zeta,xi_coupling"
!!$           factor = theta_coupling / zeta_coupling ! factor is between 1 and 1.5
!!$           factor2 = theta_coupling / xi_coupling ! factor is between 1 and 1.5
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)/2
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)*factor/2
!!$           Nxi_candidates(j) = Nxi_candidates(j-1)*factor2/2
!!$           theta_coupling = theta_coupling / 2
!!$           zeta_coupling = zeta_coupling*factor/2
!!$           xi_coupling = xi_coupling*factor2/2
!!$        elseif (zeta_coupling > theta_coupling .and. zeta_coupling > xi_coupling) then
!!$           if (masterProc) print *,"zeta_coupling > theta,xi_coupling"
!!$           factor = zeta_coupling / theta_coupling ! factor is between 1 and 1.5
!!$           factor2 = zeta_coupling / xi_coupling ! factor is between 1 and 1.5
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)/2
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)*factor/2
!!$           Nxi_candidates(j) = Nxi_candidates(j-1)*factor2/2
!!$           zeta_coupling = zeta_coupling / 2
!!$           theta_coupling = theta_coupling*factor/2
!!$           xi_coupling = xi_coupling*factor2/2
!!$        elseif (xi_coupling > theta_coupling .and. xi_coupling > zeta_coupling) then
!!$           if (masterProc) print *,"xi_coupling > theta,zeta_coupling"
!!$           factor = xi_coupling / theta_coupling ! factor is between 1 and 1.5
!!$           factor2 = xi_coupling / zeta_coupling ! factor is between 1 and 1.5
!!$           Ntheta_candidates(j) = Ntheta_candidates(j-1)*factor/2
!!$           Nzeta_candidates(j) = Nzeta_candidates(j-1)*factor2/2
!!$           Nxi_candidates(j) = Nxi_candidates(j-1)/2
!!$           theta_coupling = theta_coupling*factor/2
!!$           zeta_coupling = zeta_coupling*factor2/2
!!$           xi_coupling = xi_coupling/2
!!$        end if
        if (Ntheta_candidates(j)<=Ntheta_min) then
           Ntheta_candidates(j) = Ntheta_min
           exit
        end if
        if (Nzeta_candidates(j)<=Nzeta_min) then
           Nzeta_candidates(j) = Nzeta_min
           exit
        end if
        if (Nxi_candidates(j)<=Nxi_min) then
           Nxi_candidates(j) = Nxi_min
           exit
        end if
     end do
     N_levels = j
  end if

  if (N_levels<1) stop "Error! N_levels<1"
  allocate(levels(N_levels))

  if (N_levels==1) then
     levels(1)%Ntheta = Ntheta
     levels(1)%Nzeta = Nzeta
     levels(1)%Nxi = Nxi
  else
     ! Handle theta:
     if (coarsen_theta) then
        do j=1,N_levels
           select case (coarsen_option)
           case (1)
              temp_float = max(one*Ntheta_min, Ntheta * (0.5 ** (j-1)))
           case (2)
              temp_float = exp(log(Ntheta*one) - (j-one)/(N_levels-one)*log(Ntheta*one/Ntheta_min))
           case (3,4)
              temp_float = Ntheta_candidates(j)
           case default
              print *,"Invalid coarsen_option:",coarsen_option
              stop
           end select

           temp_int = nint(temp_float)
           if (mod(temp_int,2)==1) then
              levels(j)%Ntheta = temp_int
           else
              ! Check whether temp_int+1 or temp_int-1 is closer to the ideal value:
              if (abs(temp_int+1 - temp_float) > abs(temp_int-1 - temp_float)) then
                 levels(j)%Ntheta = temp_int-1
              else
                 levels(j)%Ntheta = temp_int+1
              end if
           end if
        end do
     else
        do j=1,N_levels
           levels(j)%Ntheta = Ntheta
        end do
     end if

     ! Handle zeta:
     if (coarsen_zeta) then
        do j=1,N_levels
           select case (coarsen_option)
           case (1)
              temp_float = max(one*Nzeta_min, Nzeta * (0.5 ** (j-1)))
           case (2)
              temp_float = exp(log(Nzeta*one) - (j-one)/(N_levels-one)*log(Nzeta*one/Nzeta_min))
           case (3,4)
              temp_float = Nzeta_candidates(j)
           case default
              print *,"Invalid coarsen_option:",coarsen_option
              stop
           end select

           temp_int = nint(temp_float)
           if (mod(temp_int,2)==1) then
              levels(j)%Nzeta = temp_int
           else
              ! Check whether temp_int+1 or temp_int-1 is closer to the ideal value:
              if (abs(temp_int+1 - temp_float) > abs(temp_int-1 - temp_float)) then
                 levels(j)%Nzeta = temp_int-1
              else
                 levels(j)%Nzeta = temp_int+1
              end if
           end if
        end do
     else
        do j=1,N_levels
           levels(j)%Nzeta = Nzeta
        end do
     end if

     ! Handle xi:
     if (coarsen_xi) then
        select case (coarsen_option)
        case (1)
           do j=1,N_levels
              levels(j)%Nxi = max(Nxi_min, nint((Nxi-1) * (0.5 ** (j-1)) + 1))
           end do
        case (2)
           do j=1,N_levels
              levels(j)%Nxi = nint(exp(log(Nxi-one) - (j-one)/(N_levels-one)*log((Nxi-one)/(Nxi_min-one)))) + 1
           end do
        case (3,4)
           do j=1,N_levels
              levels(j)%Nxi = nint(Nxi_candidates(j))
           end do
        case default
           print *,"Invalid coarsen_option:",coarsen_option
           stop
        end select

     else
        do j=1,N_levels
           levels(j)%Nxi = Nxi
        end do
     end if

  end if

  ! Sometimes, the last 2 levels turn out to have exactly the same parameters, due to rounding. Prevent this:
  if (levels(N_levels)%Ntheta==levels(N_levels-1)%Ntheta .and. levels(N_levels)%Nzeta==levels(N_levels-1)%Nzeta &
       .and. levels(N_levels)%Nxi==levels(N_levels-1)%Nxi) N_levels = N_levels-1

  if (masterProc) then
     print *,"----- Computed parameters for multigrid: ----"
     print *," Level Ntheta  Nzeta    Nxi"
     do j=1,N_levels
        print "(i7,i7,i7,i7)",j,levels(j)%Ntheta,levels(j)%Nzeta,levels(j)%Nxi
     end do
  end if

  if (Ntheta .ne. levels(1)%Ntheta) stop "Error! Ntheta should equal Ntheta_levels(1)"
  if (Nzeta  .ne. levels(1)%Nzeta ) stop "Error! Nzeta should equal Nzeta_levels(1)"
  if (Nxi    .ne. levels(1)%Nxi   ) stop "Error! Nxi should equal Nxi_levels(1)"

end subroutine set_grid_resolutions


