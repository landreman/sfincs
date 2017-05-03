#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

subroutine nonperiodic_interpolation(N, M, x, y, matrix)
  ! Builds a matrix for interpolating from a non-periodic grid to any
  ! other grid, using a 2-point stencil (linear interpolation).
  !
  ! Inputs:
  ! N = number of points in the uniform grid on which we know a function.
  ! M = number of points on the grid to which we want to interpolate.
  ! x(N) = locations of the (old) uniform grid points.
  ! y(M) = locations of the (new) grid onto which we interpolate.
  !
  ! Outputs:
  ! matrix(M,N) = interpolation matrix.
  !
  ! Outside the range of the input x grid, we use a constant extrapolation.

  implicit none

  integer, intent(in) :: N, M
  PetscScalar, intent(in) :: x(N), y(M)
  PetscScalar, intent(out) :: matrix(M,N)
  integer :: i, j, index=1
  logical flag
  integer :: indicesToUse(2), k
  PetscScalar :: x0, x1, xx

  ! Initialize matrix to 0:
  matrix=0d+0

  if (N<2) then
     print *,"Error! Uniform grid x must contain at least 2 points."
     stop
  end if

  ! Test to make sure the x grid points are sorted in increasing order: 
!!$  interval = x(2)-x(1)
  do j=2,N
     if (x(j-1) >= x(j)) then 
        print *,"Error! x grid points are not sorted in increasing order."
        stop
     end if
!!$     if (abs(x(j)-x(j-1)-interval) > 1d-10) then
!!$        print *,"Error! x grid points must be uniformly spaced."
!!$        print *,"Grid points:",x
!!$        stop
!!$     end if
  end do

  ! Loop over each point in the new grid:
  do i=1,M
     if (y(i) <= x(1)) then
        ! We are extrapolating off the left end of the x grid.
        matrix(i,1) = 1
     elseif (y(i) >= x(N)) then
        ! We are extrapolating off the right end of the x grid.
        matrix(i,N) = 1
     else
        ! We are in the interior of the x grid.

        ! Set 'index' to point to the first element in the x grid that is
        ! >= to the desired location:
        flag=.true.
        do j=1,N
           if (flag .and. (x(j) >= y(i))) then
              flag = .false.
              index = j
           end if
        end do

        if (flag) then
           ! y(i) exceeds all points in the x grid, so we are extrapolating
           ! off the right end of the x grid.
           stop "Error 1 in nonperiodic_interpolation."
        end if

        if (index==1) then
           indicesToUse = [1,2]
        else if (index .le. N) then
           indicesToUse = [index-1, index]
        else
           print *,"Error 2 in nonperiodic_interpolation! Unexpected value for index."
           print *,"index = ",index
           print *,"N = ",N
           stop
        end if
        
        x0 = x(indicesToUse(1))
        x1 = x(indicesToUse(2))
        xx = y(i)
        
        matrix(i,indicesToUse(1))=(xx-x1)/(x0-x1)
        matrix(i,indicesToUse(2))=(xx-x0)/(x1-x0)
     end if
  end do

  ! Sanity test, which can be commented out for speed: row sums should all be 1.
  do j=1,M
     if (abs(sum(matrix(j,:))-1) > 1d-12) then
        print *,"Error in periodic_interpolation! Row sums are not 1"
        print *,"Here comes interpolation matrix:"
        do k=1,M
           print *,matrix(k,:)
        end do
        stop
     end if
  end do

end subroutine nonperiodic_interpolation

