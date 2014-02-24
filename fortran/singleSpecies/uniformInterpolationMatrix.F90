#include <finclude/petscsysdef.h>

subroutine interpolationMatrix(N, M, x, y, matrix, scheme, L)
  ! Builds a matrix for interpolating from a uniform grid to any
  ! other grid.
  ! Based on the Matlab function 
  ! m20121127_02_makeHighOrderInterpolationMatrix.m.
  ! 
  ! Created by Matt Landreman, 
  ! Massachusetts Institute of Technology, Plasma Science & Fusion Center, 2012.
  !
  ! Inputs:
  ! N = number of points in the uniform grid on which we know a function.
  ! M = number of points on the grid to which we want to interpolate.
  ! x(N) = locations of the (old) uniform grid points.
  ! y(M) = locations of the (new) grid onto which we interpolate.
  ! scheme = switch that determines the extrapolation method (see below).
  ! L = Legendre index, used for some extrapolation methods.
  !
  ! Outputs:
  ! matrix(M,N) = interpolation matrix.
  !
  ! The x grid points must be sorted in increasing order.
  ! Extrapolation to the right of the x grid is allowed but extrapolation
  ! to the left of the x grid is not allowed, since the former
  ! can occur but the latter should not.
  !
  ! Possible values for scheme:
  ! -1 = Do not allow extrapolation to the right.
  ! 0 = return 0 whenever extrapolating to the right.
  ! 1 = extrapolate as 1/(y^(L+1)) 
  !     (Suitable for the first Rosenbluth potential.)
  ! 2 = extrapolate using equation (27) of
  !     http://arxiv.org/pdf/1210.5289v1.pdf
  !     (suitable for d^2/dv^2 of the 2nd Rosenbluth potential.)
  !     Note: for the pure plasma code, it is not necessary to extrapolate the
  !     potentials off the grid, so the 'L' and 'GOrH' parameters do not
  !     matter.  They do matter in the impure-plasma code however.

  implicit none

  integer, intent(in) :: N, M, scheme, L
  PetscScalar, intent(in) :: x(N), y(M)
  PetscScalar, intent(out) :: matrix(M,N)
  integer :: i, j, index=1
  logical flag
  integer :: indicesToUse(4), k
  PetscScalar :: x0, x1, x2, x3, xx, xi, xj, denomP, denomQ, interval

  ! Initialize matrix to 0:
  matrix=0d+0

  if (N<4) then
     print *,"Error! Uniform grid x must contain at least x points."
     stop
  end if

  ! Test to make sure the x grid points are sorted in increasing order: 
  interval = x(2)-x(1)
  do j=2,N
     if (x(j-1) >= x(j)) then 
        print *,"Error! x grid points are not sorted in increasing order."
        stop
     end if
     if (abs(x(j)-x(j-1)-interval) > 1d-10) then
        print *,"Error! x grid points must be uniformly spaced."
        print *,"Grid points:",x
        stop
     end if
  end do

  ! Loop over each point in the new grid:
  do i=1,M
     if (y(i) < x(1)) then
        ! We are extrapolating off the left end of the x grid.
        ! This should never happen in sfincs.
        print *,"Error! Attempt to extrapolate to the left.  This should not happen."
        stop
     else

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

           select case (scheme)
           case (-1)
              print *,"Error! Attempt to extrapolate to the right."
              print *,"x grid:",x
              print *,"y grid:",y
              stop
           case (0)
              ! Nothing to do here.

           case (1)
              ! Extrapolation for the first Rosenbluth potential.
              ! H and x * dH/dx behave as \propto 1/(x ^ (L+1)).
              matrix(i,N) = (x(N)/y(i))**(L+1)

           case (2)
              ! Extrapolation for d^2/dx^2 of the second Rosenbluth potential.
              xi = x(N)
              xj = x(N-1)
              denomQ = xi*xi-xj*xj
              denomP = 1/(xi*xi) - 1/(xj*xj)

              matrix(i,N) = xi**(L+1)/denomP/y(i)**(L+3) &
                   + xi**(L+3)/denomQ/y(i)**(L+1)

              matrix(i,N-1) = -xj**(L+1)/denomP/y(i)**(L+3) &
                   - xj**(L+3)/denomQ/y(i)**(L+1)

           case default
              print *,"Error! Invalid setting for scheme"
              stop
           end select
        else
           ! We are interpolating rather than extrapolating.
           ! This is the normal condition.

           if (index > N-2) then
              indicesToUse = [( k+N, k=(-3),0 )]
           else if (index < 3) then
              indicesToUse = [( k, k=1,4 )]
           else
              indicesToUse = [( k+index, k=(-2),1 )]
           end if

           x0 = x(indicesToUse(1))
           x1 = x(indicesToUse(2))
           x2 = x(indicesToUse(3))
           x3 = x(indicesToUse(4))
           xx = y(i)

           matrix(i,indicesToUse(1))=(xx-x1)*(xx-x2)*(xx-x3)/((x0-x1)*(x0-x2)*(x0-x3))
           matrix(i,indicesToUse(2))=(xx-x0)*(xx-x2)*(xx-x3)/((x1-x0)*(x1-x2)*(x1-x3))
           matrix(i,indicesToUse(3))=(xx-x0)*(xx-x1)*(xx-x3)/((x2-x0)*(x2-x1)*(x2-x3))
           matrix(i,indicesToUse(4))=(xx-x0)*(xx-x1)*(xx-x2)/((x3-x0)*(x3-x1)*(x3-x2))

        end if
     end if
  end do

end subroutine interpolationMatrix

