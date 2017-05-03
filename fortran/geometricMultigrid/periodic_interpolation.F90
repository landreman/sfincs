#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

subroutine periodic_interpolation(N, M, period, y, matrix)
  ! Builds a matrix for interpolating from a uniform periodic grid to any
  ! other grid, using a 2-point stencil (linear interpolation). It is assumed
  ! that the original grid is uniform with a grid point at 0 but not at period.
  !
  ! Inputs:
  ! N = number of points in the uniform grid on which we know a function.
  ! M = number of points on the grid to which we want to interpolate.
  ! period = periodicity of the grid.
  ! y(M) = locations of the (new) grid onto which we interpolate.
  !
  ! Outputs:
  ! matrix(M,N) = interpolation matrix.
  !
  ! Outside the range of the input x grid, we use a constant extrapolation.

  implicit none

  integer, intent(in) :: N, M
  PetscScalar, intent(in) :: y(M), period
  PetscScalar, intent(out) :: matrix(M,N)
  integer :: i, j, index=1
  logical flag
  integer :: indicesToUse(2), k
  PetscScalar :: x0, x1, xx
  PetscScalar, allocatable, dimension(:) :: x, y_copy

  ! Initialize matrix to 0:
  matrix=0d+0

  if (N<2) stop "Error! Uniform grid x must contain at least 2 points."
  if (period<=0) stop "Error! Period must be > 0."

  allocate(x(N+1))
  x = [( (period*i)/N, i=0,N )]

  allocate(y_copy(M))
  ! Note that we should use 'modulo' and not 'mod' in the next line, since the latter returns negative results if any y values are <0.
  y_copy = modulo(y,period)

  ! Loop over each point in the new grid:
  do i=1,M
     if (y_copy(i) < x(1)) then
        ! We are extrapolating off the left end of the x grid.
        stop "Error 3! Program should not get here."
     elseif (y_copy(i) > x(N+1)) then
        ! We are extrapolating off the right end of the x grid.
        stop "Error 4! Program should not get here."
     else
        ! We are in the interior of the x grid.

        ! Set 'index' to point to the first element in the x grid that is
        ! >= to the desired location:
        flag=.true.
        do j=1,N+1
           if (flag .and. (x(j) >= y_copy(i))) then
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
           x0 = x(indicesToUse(1))
           x1 = x(indicesToUse(2))
        else if (index .le. N) then
           indicesToUse = [index-1, index]
           x0 = x(indicesToUse(1))
           x1 = x(indicesToUse(2))
        elseif (index == N+1) then
           indicesToUse = [index-1, index]
           x0 = x(indicesToUse(1))
           x1 = x(indicesToUse(2))
           indicesToUse = [index-1, 1]
        else
           print *,"Error 2 in nonperiodic_interpolation! Unexpected value for index."
           print *,"index = ",index
           print *,"N = ",N
           stop
        end if
        
        xx = y_copy(i)
        
        matrix(i,indicesToUse(1))=(xx-x1)/(x0-x1)
        matrix(i,indicesToUse(2))=(xx-x0)/(x1-x0)
     end if
  end do

  deallocate(x, y_copy)

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

end subroutine periodic_interpolation

