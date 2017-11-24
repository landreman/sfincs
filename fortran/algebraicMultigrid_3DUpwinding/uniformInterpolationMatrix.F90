subroutine interpolationMatrix(N, M, x, y, scheme, matrix, extrapMatrix)
  ! Builds a matrix for interpolating from a uniform grid to any
  ! other grid.
  ! Based on the Matlab function 
  ! m20131119_07_makeHighOrderUniformRegriddingMatrix_extrap.m
  ! 
  ! Created by Matt Landreman, 
  ! Massachusetts Institute of Technology, Plasma Science & Fusion Center, 2012.
  !
  ! Inputs:
  ! N = number of points in the uniform grid on which we know a function.
  ! M = number of points on the grid to which we want to interpolate.
  ! x(N) = locations of the (old) uniform grid points.
  ! y(M) = locations of the (new) grid onto which we interpolate.
  ! scheme = method for interpolation. Available options are 1 and 2:
  !          1: Use a 2-point stencil. (Low order)
  !          2: Use a 4-point stencil. (High order)
  !
  ! Outputs:
  ! matrix(M,N) = interpolation matrix.
  ! extrapMatrix(M,N) = a matrix of 1s and 0s which will be used for extrapolating the
  !     Rosenbluth potentials in the collision operator
  !
  ! The x grid points must be sorted in increasing order.
  ! Extrapolation to the right of the x grid is allowed but extrapolation
  ! to the left of the x grid is not allowed, since the former
  ! can occur in the collision operator but the latter should not.

  use kinds

  implicit none

  integer, intent(in) :: N, M, scheme
  real(prec), intent(in) :: x(N), y(M)
  real(prec), intent(out) :: matrix(M,N), extrapMatrix(M,N)
  integer :: i, j, index=1
  logical flag
  integer :: indicesToUse(4), k
  real(prec) :: x0, x1, x2, x3, xx, xi, xj, denomP, denomQ, interval

  ! Initialize matrix to 0:
  matrix=0d+0
  extrapMatrix=0d+0

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
        ! This should never happen in the collision operator.
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

           extrapMatrix(i,N)=1

        else
           ! We are interpolating rather than extrapolating.
           ! This is the normal condition.

           select case (scheme)
           case (1)
              if (index==1) then
                 indicesToUse = [1,2,-1,-1]
                 ! Only the first 2 entries of indicesToUse matter.
              else if (index .le. N) then
                 indicesToUse = [index-1, index,-1,-1]
              else
                 print *,"Error! Unexpected value for index."
                 print *,"index = ",index
                 print *,"N = ",N
                 stop
              end if

              x0 = x(indicesToUse(1))
              x1 = x(indicesToUse(2))
              xx = y(i)

              matrix(i,indicesToUse(1))=(xx-x1)/(x0-x1)
              matrix(i,indicesToUse(2))=(xx-x0)/(x1-x0)

           case (2)
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

           case default
              print *,"Error! Invalid value for scheme"
              stop
           end select
        end if
     end if
  end do

end subroutine interpolationMatrix

