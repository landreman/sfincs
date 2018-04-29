#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

! This subroutine is used for cubic spline interpolation on a periodic interval
! x_coarse and y_coarse are arrays of length nx_coarse (allocated)
! x_fine, y_fine, and y_interp are arrays of length nx_fine
! period is the the value such that y(x=period) = y(x=0)
! y_interp will be returned with the interpolated values of y_coarse

subroutine periodicCubicInterpolation(nx_coarse,nx_fine,period,x_coarse,x_fine,y_coarse,y_interp)

  use globalVariables, only: zero, y2

  implicit none

  integer, intent(in) :: nx_coarse, nx_fine
  PetscScalar, intent(in) :: x_coarse(nx_coarse), x_fine(nx_fine), y_coarse(nx_coarse), period
  PetscScalar, intent(out) :: y_interp(nx_fine)
  PetscScalar, dimension(:), allocatable :: x_appended, u
  PetscScalar :: dx_coarse, yp1, ypn, p, un, A, B, C, D, x0, x1
  integer :: i, k, ix_fine, ix_coarse, index, index1, index2
  logical :: flag

  allocate(x_appended(nx_coarse+1))
  x_appended(1:nx_coarse) = x_coarse
  x_appended(nx_coarse+1) = period

  ! Compute vector of second derivatives
  allocate(y2(nx_coarse))
  allocate(u(nx_coarse))
  y2 = zero
  u = zero
  y_interp = zero

  dx_coarse = x_coarse(2)-x_coarse(1)
  ! Estimate of derivative at index = 1 using centered differencing
  yp1 = (-y_coarse(3)+8*y_coarse(2)-8*y_coarse(nx_coarse)+y_coarse(nx_coarse-1))/(12*dx_coarse)
  ! Estimate of derivative at index = nx_coarse
  ypn = (-y_coarse(2)+8*y_coarse(1)-8*y_coarse(nx_coarse-1)+y_coarse(nx_coarse-2))/(12*dx_coarse)
  y2(1) = -0.5
  u(1) = (3.0/dx_coarse)*((y_coarse(2)-y_coarse(1))/dx_coarse-yp1)
  do i=2,(nx_coarse-1)
    p = 0.5*y2(i-1)+2.0
    y2(i) = -0.5/p
    u(i) = (y_coarse(i+1)-y_coarse(i))/dx_coarse - (y_coarse(i)-y_coarse(i-1))/dx_coarse
    u(i) = (6.0*u(i)/(2*dx_coarse)-0.5*u(i-1))/p
  end do
  un = (3.0/dx_coarse)*(ypn-(y_coarse(nx_coarse)-y_coarse(nx_coarse-1))/dx_coarse)
  y2(nx_coarse) = (un-0.5*u(nx_coarse-1))/(0.5*y2(nx_coarse-1)+1.0)
  do k=(nx_coarse-1),1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  end do

  ! For each element in fine grid, find nearest point in coarse grid to the right
  do ix_fine=1,nx_fine
    flag = .true.
    do ix_coarse=1,nx_coarse+1
      if (flag .and. (x_appended(ix_coarse) >= x_fine(ix_fine))) then
        flag = .false.
        index = ix_coarse
      end if
    end do
    if (flag) then
      print *,"Error! Should not reach this point in periodicCubicInterpolation."
      stop
    end if
    if (index==1) then
      x0 = x_coarse(nx_coarse)
      x1 = x_coarse(1)
      index1 = nx_coarse
      index2 = 1
    elseif (index <= nx_coarse+1) then
      x0 = x_appended(index-1)
      x1 = x_appended(index)
      index1 = index-1
      if (index==nx_coarse+1) then
        index2 = 1
      else
        index2 = index
      end if
    else
      print *,"Unexpected index in periodicCubicInterpolation."
      stop
    end if
    A = (x1-x_fine(ix_fine))/(x1-x0)
    B = 1-A
    C = (1/6.0)*(A**3-A)*(x1-x0)**2
    D = (1/6.0)*(B**3-B)*(x1-x0)**2
    y_interp(ix_fine) = A*y_coarse(index1) + B*y_coarse(index2) + C*y2(index1) + D*y2(index2)
  end do

end subroutine periodicCubicInterpolation
