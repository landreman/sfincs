subroutine uniformDiffMatrices(N, xMin, xMax, scheme, x, weights, ddx, d2dx2)
  ! Finite difference and spectral differentiation matrices and integration
  ! weights for a uniform grid.
  !
  ! Created by Matt Landreman, 
  ! Massachusetts Institute of Technology, Plasma Science & Fusion Center, 2012.
  !
  ! Inputs:
  !   N = number of grid points.
  !   xMin = minimum value in the domain.
  !   xMax = maximum value in the domain.
  !   scheme = switch for controlling order of accuracy for differentiation
  !            and handling of endpoints.
  !
  ! Options for scheme:
  ! 0 =  The domain [xMin, xMax] is assumed to be periodic. A 3-point stencil
  !      is used everywhere. A grid point will be placed at xMin but not 
  !      xMax.
  ! 1 =  Same as scheme=0, except a grid point will be placed at xMax but not
  !      xMin.
  ! 2 =  The domain [xMin, xMax] is assumed to be non-periodic. A 3-point 
  !      stencil is used everywhere.  The first and last row of the
  !      differentiation matrices will use one-sided differences, so they
  !      will each have a non-tridiagonal element.
  ! 3 =  The same as scheme=2, except that the first differentiation matrix
  !      will use a 2-point 1-sided stencil for the first and last elements
  !      so the matrix is strictly tri-diagonal.  The 2nd derivative matrix
  !      is the same as for option 2, since it is not possible to compute 
  !      the 2nd derivative with only a 2-point stencil.
  ! 10 = The domain [xMin, xMax] is assumed to be periodic. A 5-point stencil
  !      is used everywhere. A grid point will be placed at xMin but not 
  !      xMax.  This option is like scheme=0 but more accurate.
  ! 11 = Same as scheme=10, except a grid point will be placed at xMax but
  !      not xMin.  This option is like scheme=1 but more accurate.
  ! 12 = The domain [xMin, xMax] is assumed to be non-periodic. A 5-point 
  !      stencil is used everywhere.  The first two and last two rows of 
  !      the differentiation matrices will then each have non-pentadiagonal 
  !      elements.
  ! 13 = The same as option 12, except that 3-point stencils are used for the
  !      first and last rows of the differentiation matrices, and 4-point 
  !      stencils are used for the 2nd and penultimate rows of the 
  !      differentiation matrices.  With this option, both differentiation 
  !      matrices are strictly penta-diagonal.
  ! 20 = The domain [xMin, xMax] is assumed to be periodic. Spectral
  !      differentiation matrices are returned. A grid point will be placed 
  !      at xMin but not xMax.
  ! 21 = Same as scheme=20, except a grid point will be placed at xMax but not
  !      xMin.
  ! 30 = Periodic with a grid point at xMin but not xMax.  Upwinding to the
  !      left. A 2-point stencil is used for the first derivative and a 
  !      3-point stencil is used for the second derivative.
  ! 31 = Periodic with a grid point at xMax but not xMin.  Upwinding to the
  !      left. A 2-point stencil is used for the first derivative and a 
  !      3-point stencil is used for the second derivative.
  ! 32 = Aperiodic.  Upwinding to the left. A 2-point stencil is used for the
  !      first derivative and a 3-point stencil is used for the second 
  !      derivative.  The top row of D and the top two rows of DD are zero.
  ! 40 = Same as 30 but upwinding to the right.
  ! 41 = Same as 31 but upwinding to the right.
  ! 42 = Same as 32 but upwinding to the right.  The bottom row of D and the
  !      bottom two rows of DD are zero.
  ! 50 = Periodic with a grid point at xMin but not xMax.  Upwinding to the
  !      left. A 3-point stencil is used for both derivatives.
  ! 51 = Periodic with a grid point at xMax but not xMin.  Upwinding to the
  !      left. A 3-point stencil is used for both derivatives.
  ! 52 = Aperiodic.  Upwinding to the left. A 3-point stencil is used for 
  !      both derivatives.  The top row of both matrices is all zero. The
  !      second row of D uses a 2-point stencil.
  ! 60 = Same as 50 but upwinding to the right.
  ! 61 = Same as 51 but upwinding to the right.
  ! 62 = Same as 52 but upwinding to the right.  The bottom row of both
  !      derivative matrices is all zero. The penultimate row of D uses a
  !      2-point stencil.
  ! 80 = Periodic with a grid point at xMin but not xMax.  The first derivative is upwinded to the
  !      left. A stencil is used with 1 point on 1 side and 2 points on the
  !      other side. The second derivative is the same as in scheme 0.
  ! 81 = Same as 80 but with a grid point at xMax and not xMin.
  ! 82 = Like 80 but not periodic, with a grid point at both xMin and xMax.
  !      The top row of D is zero.
  ! 90 = Same as 80 but upwinding to the right.
  ! 91 = Same as 90 but with a grid point at xMax and not xMin.
  ! 92 = Like 90 but not periodic, with a grid point at both xMin and xMax.
  !      The bottom row of D is zero.
  ! 100 = Periodic with a grid point at xMin but not xMax.
  !       1st derivative only. Upwinded to the left.
  !       Stencil has 1 point on 1 side, 3 points on the other side
  ! 101 = Same as 100, but with no grid point at xMin and with a grid point at xMax.
  ! 102 = Same as 100, but aperiodic, with grid points at both xMin and xMax.
  ! 110 = Periodic with a grid point at xMin but not xMax.
  !       1st derivative only. Upwinded to the right.
  !       Stencil has 1 point on 1 side, 3 points on the other side
  ! 111 = Same as 110, but with no grid point at xMin and with a grid point at xMax.
  ! 112 = Same as 100, but aperiodic, with grid points at both xMin and xMax.
  ! 120 = Periodic with a grid point at xMin but not xMax.
  !       1st derivative only. Upwinded to the left.
  !       Stencil has 2 points on 1 side, 3 points on the other side
  ! 121 = Same as 120, but with no grid point at xMin and with a grid point at xMax.
  ! 122 = Same as 120, but aperiodic, with grid points at both xMin and xMax.
  ! 130 = Periodic with a grid point at xMin but not xMax.
  !       1st derivative only. Upwinded to the right.
  !       Stencil has 2 points on 1 side, 3 points on the other side
  ! 131 = Same as 130, but with no grid point at xMin and with a grid point at xMax.
  ! 132 = Same as 130, but aperiodic, with grid points at both xMin and xMax.
  
  !
  ! Outputs:
  !   x = column vector with the grid points.
  !   weights = column vector with the weights for integration using the trapezoid rule.
  !   ddx = matrix for differentiation.
  !   d2dx2 = matrix for the 2nd derivative.

  use globalVariables, only: prec

  implicit none

  integer, intent(in) :: N, scheme
  real(prec), intent(in) :: xMin, xMax
  real(prec), intent(out), dimension(N) :: x, weights
  real(prec), intent(out), dimension(N,N) :: ddx, d2dx2
  integer :: i
  real(prec) :: dx, dx2
  real(prec) :: h
  integer :: n1, n2
  real(prec), allocatable :: topc(:), col1(:)
  real(prec), parameter :: pi = 3.1415926535897932384626433d+0

  ! ***************************************************************
  ! Validate input
  ! ***************************************************************

  if (N<2) then
     print *,"Error! N must be at least 2."
     print *,"N = ",N
     stop
  end if
  if (xMin > xMax) then 
     print *,"Error! xMax must be larger than xMin"
     print *,"xMax=",xMax
     print *,"xMin=",xMin
     stop 
  end if
  if (xMin == xMax) then 
     print *,"Error! xMax cannot equal xMin" 
     print *,"xMin=xMax=",xMin
     stop 
  end if

  ! ***************************************************************
  ! Set gridpoints
  ! ***************************************************************

  select case (scheme)
  case (2, 3, 12, 13, 32, 42, 52, 62, 82, 92, 102, 112, 122, 132)
     ! Include points at both xMin and xMax:
     x = [( (xMax-xMin)*i/(N-1)+xMin, i=0,N-1 )]
  case (0,10,20,30,40,50,60,80,90,100,110,120,130)
     ! Include a point at xMin but not xMax:
     x = [( (xMax-xMin)*i/(N)+xMin, i=0,N-1 )]
  case (1,11,21,31,41,51,61,81,91,101,111,121,131)
     ! Include a point at xMax but not xMin:
     x = [( (xMax-xMin)*i/(N)+xMin, i=1,N )]
  case default
     print *,"Error! Invalid value for scheme"
     stop
  end select

  dx=x(2)-x(1)
  dx2=dx*dx

  ! ***************************************************************
  ! Set integration weights
  ! ***************************************************************

  weights=dx
  select case (scheme)
  case (2, 3, 12, 13, 32, 42, 52, 62, 82, 92, 102, 112, 122, 132)
     ! Grid is aperiodic
     weights(1)=weights(1)/2
     weights(N)=weights(N)/2
  end select

  ! ***************************************************************
  ! Fill the interior of the differentiation matrices
  ! ***************************************************************

  ddx=0d+0
  d2dx2=0d+0

  select case (scheme)
  case (0,1,2,3)
     ! 2nd order (3 point stencil):
     if (N<3) then
        print *,"Error! N must be at least 3 for the 3-point stencil methods"
        stop
     end if
     do i=2,N-1
        ddx(i,i+1)=1/(2*dx)
        ddx(i,i-1)=-1/(2*dx)

        d2dx2(i,i+1)=1/(dx2)
        d2dx2(i,i)=-2/(dx2)
        d2dx2(i,i-1)=1/(dx2)
     end do

  case (10,11,12,13)
     ! 4th order (5 point stencil):
     if (N<5) then
        print *,"Error! N must be at least 5 for the 5-point stencil methods"
        stop
     end if
     do i=3,N-2
        ddx(i,i+2)=-1/(6*2*dx)
        ddx(i,i+1)=4/(3*2*dx)
        ddx(i,i-1)=-4/(3*2*dx)
        ddx(i,i-2)=1/(6*2*dx)

        d2dx2(i,i+2)=-1/(12*dx2)
        d2dx2(i,i+1)=4/(3*dx2)
        d2dx2(i,i)=-5/(2*dx2)
        d2dx2(i,i-1)=4/(3*dx2)
        d2dx2(i,i-2)=-1/(12*dx2)
     end do

  case (30,31,32)
     ! 2-point stencil for D and 3-point stencil for DD,
     ! upwinding to the left.
     if (N<3) then
        print *,"Error! N must be at least 3 for this scheme."
        stop
     end if
     do i=2,N
        ddx(i,i) = 1/dx
        ddx(i,i-1) = -1/dx
     end do
     do i=3,N
        d2dx2(i,i) = 1/dx2
        d2dx2(i,i-1) = -2/dx2
        d2dx2(i,i-2) = 1/dx2
     end do

  case (40,41,42)
     ! 2-point stencil for D and 3-point stencil for DD,
     ! upwinding to the right.
     if (N<3) then
        print *,"Error! N must be at least 3 for this scheme."
        stop
     end if
     do i=1,N-1
        ddx(i,i) = -1/dx
        ddx(i,i+1) = 1/dx
     end do
     do i=1,N-2
        d2dx2(i,i) = 1/dx2
        d2dx2(i,i+1) = -2/dx2
        d2dx2(i,i+2) = 1/dx2
     end do

  case (50,51,52)
     ! 3-point stencil for both D and DD,
     ! upwinding to the left.
     if (N<3) then
        print *,"Error! N must be at least 3 for this scheme."
        stop
     end if
     do i=3,N
        ddx(i,i) = (1.5d+0)/dx
        ddx(i,i-1) = -2/dx
        ddx(i,i-2) = 1/(2*dx)

        d2dx2(i,i) = 1/dx2
        d2dx2(i,i-1) = -2/dx2
        d2dx2(i,i-2) = 1/dx2
     end do

  case (60,61,62)
     ! 3-point stencil for both D and DD,
     ! upwinding to the right.
     if (N<3) then
        print *,"Error! N must be at least 3 for this scheme."
        stop
     end if
     do i=1,N-2
        ddx(i,i) = -(1.5d+0)/dx
        ddx(i,i+1) = 2/dx
        ddx(i,i+2) = -1/(2*dx)

        d2dx2(i,i) = 1/dx2
        d2dx2(i,i+1) = -2/dx2
        d2dx2(i,i+2) = 1/dx2
     end do

  case (20,22)
     ! Create spectral differentiation matrices.
     ! Here I've translated the Matlab fourdif.m routine from the
     ! DMSuite package by S.C. Reddy and J.A.C. Weideman, available at
     ! http://www.mathworks.com/matlabcentral/fileexchange/29
     ! or here:
     ! http://dip.sun.ac.za/~weideman/research/differ.html

     h = 2*pi/N;
     n1 = floor((N-1.0)/2)
     n2 = ceiling((N-1.0)/2)
     allocate(topc(n2))
     allocate(col1(N))

     ! Create first derivative matrix:
     col1(1)=0d+0
     if (mod(N,2)==0) then
        topc = [( 0.5d+0/tan(i*h/2), i=1,n2 )]
        col1(2:n2+1) = topc
        col1(n2+2:) = -topc(n1:1:-1)
        do i=2,N,2
           col1(i) = -col1(i)
        end do
        col1 = 2*pi/(xMax-xMin)*col1
        ! Create a toeplitz matrix:
        do i=1,N
           ddx(i,i:) = -col1(1:N+1-i)
           ddx(i,1:i-1) = col1(i:2:-1)
        end do
     else
        topc = [( 0.5d+0 / sin(i*h/2), i=1,n2 )]
        col1(2:n2+1) = topc
        col1(n2+2:) = topc(n1:1:-1)
        do i=2,N,2
           col1(i) = -col1(i)
        end do
        col1 = 2*pi/(xMax-xMin)*col1
        ! Create a toeplitz matrix:
        do i=1,N
           ddx(i,i:) = -col1(1:N+1-i)
           ddx(i,1:i-1) = col1(i:2:-1)
        end do
     end if

     ! Create second derivative matrix:
     if (mod(N,2)==0) then
        col1(1)=-pi*pi/(3*h*h)-1d+0/6
        topc = [( -(0.5d+0)/(sin(i*h/2)**2), i=1,n2 )]
        col1(2:n2+1) = topc
        col1(n2+2:) = topc(n1:1:-1)
        do i=2,N,2
           col1(i) = -col1(i)
        end do
        col1 = (2*pi/(xMax-xMin))**2 *col1
        ! Create a toeplitz matrix:
        do i=1,N
           d2dx2(i,i:) = col1(1:N+1-i)
           d2dx2(i,1:i-1) = col1(i:2:-1)
        end do
     else
        col1(1)=-pi*pi/(3*h*h) + 1d+0/12
        topc = [( -(0.5d+0) / (sin(i*h/2) * tan(i*h/2)), i=1,n2 )]
        col1(2:n2+1) = topc
        col1(n2+2:) = -topc(n1:1:-1)
        do i=2,N,2
           col1(i) = -col1(i)
        end do
        col1 = (2*pi/(xMax-xMin))**2 *col1
        ! Create a toeplitz matrix:
        do i=1,N
           d2dx2(i,i:) = col1(1:N+1-i)
           d2dx2(i,1:i-1) = col1(i:2:-1)
        end do
     end if

     deallocate(topc)

  case (80,81)
     ! 4 point stencil (upwinding, with 1 point on 1 side, and 2 points on the other side.)

     if (N<5) then
        print *,"Error! N must be at least 5 for 4 point stencil"
        stop
     end if
     do i=1,N
        ddx(i,modulo(i,N)+1)   =  1/(3*dx)
        ddx(i,i)               =  1/(2*dx)
        ddx(i,modulo(i-2,N)+1) = -1/(dx)
        ddx(i,modulo(i-3,N)+1) =  1/(6*dx)

        d2dx2(i,modulo(i,N)+1)   =  1/(dx2)
        d2dx2(i,i)               = -2/(dx2)
        d2dx2(i,modulo(i-2,N)+1) =  1/(dx2)
     end do

  case (82)
     ! 4 point stencil (upwinding, with 1 point on 1 side, and 2 points on the other side.)

     if (N<5) then
        print *,"Error! N must be at least 5 for 4 point stencil"
        stop
     end if
     do i=3,N-1
        ddx(i,i+1) =  1/(3*dx)
        ddx(i,i)   =  1/(2*dx)
        ddx(i,i-1) = -1/(dx)
        ddx(i,i-2) =  1/(6*dx)
     end do
     do i=2,N-1
        d2dx2(i,i+1) =  1/(dx2)
        d2dx2(i,i)   = -2/(dx2)
        d2dx2(i,i-1) =  1/(dx2)
     end do

  case (90,91)
     ! 4 point stencil (upwinding, with 1 point on 1 side, and 2 points on the other side.)

     if (N<5) then
        print *,"Error! N must be at least 5 for 4 point stencil"
        stop
     end if
     do i=1,N
        ddx(i,modulo(i-2,N)+1) = -1/(3*dx)
        ddx(i,i)               = -1/(2*dx)
        ddx(i,modulo(i,N)+1)   =  1/(dx)
        ddx(i,modulo(i+1,N)+1) = -1/(6*dx)

        d2dx2(i,modulo(i,N)+1)   =  1/(dx2)
        d2dx2(i,i)               = -2/(dx2)
        d2dx2(i,modulo(i-2,N)+1) =  1/(dx2)
     end do

  case (92)
     ! 4 point stencil (upwinding, with 1 point on 1 side, and 2 points on the other side.)

     if (N<5) then
        print *,"Error! N must be at least 5 for 4 point stencil"
        stop
     end if
     do i=2,N-2
        ddx(i,i-1) = -1/(3*dx)
        ddx(i,i)   = -1/(2*dx)
        ddx(i,i+1) =  1/(dx)
        ddx(i,i+2) = -1/(6*dx)
     end do
     do i=2,N-1
        d2dx2(i,i+1) =  1/(dx2)
        d2dx2(i,i)   = -2/(dx2)
        d2dx2(i,i-1) =  1/(dx2)
     end do

  case (100,101)
     ! upwinding, with 1 point on 1 side, and 3 points on the other side.

     if (N<5) then
        print *,"Error! N must be at least 5 for scheme 100,101"
        stop
     end if
     do i=1,N
        ddx(i,modulo(i+0,N)+1) =  1/(4*dx)
        ddx(i,i)               =  5/(6*dx)
        ddx(i,modulo(i-2,N)+1) = -3/(2*dx)
        ddx(i,modulo(i-3,N)+1) =  1/(2*dx)
        ddx(i,modulo(i-4,N)+1) = -1/(12*dx)
     end do

  case (102)
     ! upwinding, with 1 point on 1 side, and 3 points on the other side.

     if (N<5) then
        print *,"Error! N must be at least 5 for scheme 102"
        stop
     end if
     do i=4,N-1
        ddx(i,i+1) =  1/(4*dx)
        ddx(i,i)   =  5/(6*dx)
        ddx(i,i-1) = -3/(2*dx)
        ddx(i,i-2) =  1/(2*dx)
        ddx(i,i-3) = -1/(12*dx)
     end do

  case (110,111)
     ! upwinding, with 1 point on 1 side, and 2 points on the other side.

     if (N<5) then
        print *,"Error! N must be at least 5 for scheme 110,111"
        stop
     end if
     do i=1,N
        ddx(i,modulo(i-2,N)+1) = -1/(4*dx)
        ddx(i,i)               = -5/(6*dx)
        ddx(i,modulo(i+0,N)+1) =  3/(2*dx)
        ddx(i,modulo(i+1,N)+1) = -1/(2*dx)
        ddx(i,modulo(i+2,N)+1) =  1/(12*dx)
     end do

  case (112)
     ! upwinding, with 1 point on 1 side, and 2 points on the other side.

     if (N<5) then
        print *,"Error! N must be at least 5 for scheme 112"
        stop
     end if
     do i=2,N-3
        ddx(i,i-1) = -1/(4*dx)
        ddx(i,i)   = -5/(6*dx)
        ddx(i,i+1) =  3/(2*dx)
        ddx(i,i+2) = -1/(2*dx)
        ddx(i,i+3) =  1/(12*dx)
     end do

  case (120,121)
     ! upwinding, with 2 points on 1 side, and 3 points on the other side.

     if (N<5) then
        print *,"Error! N must be at least 5 for scheme 120,121"
        stop
     end if
     do i=1,N
        ddx(i,modulo(i+1,N)+1) = -1/(20*dx)
        ddx(i,modulo(i+0,N)+1) =  1/(2*dx)
        ddx(i,i)               =  1/(3*dx)
        ddx(i,modulo(i-2,N)+1) = -1/(dx)
        ddx(i,modulo(i-3,N)+1) =  1/(4*dx)
        ddx(i,modulo(i-4,N)+1) = -1/(30*dx)
     end do

!  case (122)
     ! Not implemented yet!

  case (130,131)
     ! upwinding, with 2 points on 1 side, and 3 points on the other side.

     if (N<5) then
        print *,"Error! N must be at least 5 for scheme 130,131"
        stop
     end if
     do i=1,N
        ddx(i,modulo(i-3,N)+1) =  1/(20*dx)
        ddx(i,modulo(i-2,N)+1) = -1/(2*dx)
        ddx(i,i)               = -1/(3*dx)
        ddx(i,modulo(i+0,N)+1) =  1/(dx)
        ddx(i,modulo(i+1,N)+1) = -1/(4*dx)
        ddx(i,modulo(i+2,N)+1) =  1/(30*dx)

     end do

!  case (132)
     ! Not implemented yet!


  end select

  ! ***************************************************************
  ! Handle endpoints of grid in differentiation matrices
  ! ***************************************************************

  select case (scheme)
  case (0,1)
     ddx(1,N) = -1/(2*dx)
     ddx(1,2) = 1/(2*dx)
     ddx(N,1) = 1/(2*dx)
     ddx(N,N-1) = -1/(2*dx)

     d2dx2(1,1) = -2/dx2
     d2dx2(N,N) = -2/dx2
     d2dx2(1,N) = 1/dx2
     d2dx2(1,2) = 1/dx2
     d2dx2(N,1) = 1/dx2
     d2dx2(N,N-1) = 1/dx2


  case (2)
     ! 3-point stencil, aperiodic
     ddx(1,1)=-1.5/dx
     ddx(1,2)=2/dx
     ddx(1,3)=-0.5/dx

     ddx(N,N)=1.5/dx
     ddx(N,N-1)=-2/dx
     ddx(N,N-2)=0.5/dx

     d2dx2(1,1)=1/dx2
     d2dx2(1,2)=-2/dx2
     d2dx2(1,3)=1/dx2

     d2dx2(N,N)=1/dx2
     d2dx2(N,N-1)=-2/dx2
     d2dx2(N,N-2)=1/dx2

  case (3)
     ! Aperiodic.
     ! 2-point stencil for the first and last rows of the first
     ! differentiation matrix, so the matrix is strictly tri-diagonal.
     ! The 2nd derivative matrix is the same as for scheme=0 (i.e. not
     ! strictly tri-diagonal) since it is not possible to approximate
     ! the 2nd derivative with a 2-point stencil.

     ddx(1,1)=-1/dx
     ddx(1,2)=1/dx

     ddx(N,N)=1/dx
     ddx(N,N-1)=-1/dx

     d2dx2(1,1)=1/dx2
     d2dx2(1,2)=-2/dx2
     d2dx2(1,3)=1/dx2

     d2dx2(N,N)=1/dx2
     d2dx2(N,N-1)=-2/dx2
     d2dx2(N,N-2)=1/dx2

  case (10,11)
     ! 5-point stencil, periodic

     ddx(1, N) = -(4.0d+0/3)/(2*dx)
     ddx(1, N-1) = (1.0d+0/6)/(2*dx)
     ddx(2, N) = (1.0d+0/6)/(2*dx)

     ddx(N, 1) = (4.0d+0/3)/(2*dx)
     ddx(N, 2) = -(1.0d+0/6)/(2*dx)
     ddx(N-1, 1) = -(1.0d+0/6)/(2*dx)

     d2dx2(1, N) = (4d+0/3)/dx2
     d2dx2(1, N-1) = -(1d+0/12)/dx2
     d2dx2(2, N) = -(1d+0/12)/dx2

     d2dx2(N, 1) = (4d+0/3)/dx2
     d2dx2(N, 2) = -(1d+0/12)/dx2
     d2dx2(N-1, 1) = -(1d+0/12)/dx2


     ! i=1
     ddx(1,1+2)=-1/(6*2*dx)
     ddx(1,1+1)=4/(3*2*dx)

     d2dx2(1,1+2)=-1/(12*dx2)
     d2dx2(1,1+1)=4/(3*dx2)
     d2dx2(1,1)=-5/(2*dx2)

     ! i=2
     ddx(2,2+2)=-1/(6*2*dx)
     ddx(2,2+1)=4/(3*2*dx)
     ddx(2,2-1)=-4/(3*2*dx)

     d2dx2(2,2+2)=-1/(12*dx2)
     d2dx2(2,2+1)=4/(3*dx2)
     d2dx2(2,2)=-5/(2*dx2)
     d2dx2(2,2-1)=4/(3*dx2)

     ! i=N
     ddx(N,N-1)=-4/(3*2*dx)
     ddx(N,N-2)=1/(6*2*dx)

     d2dx2(N,N)=-5/(2*dx2)
     d2dx2(N,N-1)=4/(3*dx2)
     d2dx2(N,N-2)=-1/(12*dx2)

     ! i=N-1
     ddx(N-1,N-1+1)=4/(3*2*dx)
     ddx(N-1,N-1-1)=-4/(3*2*dx)
     ddx(N-1,N-1-2)=1/(6*2*dx)

     d2dx2(N-1,N-1+1)=4/(3*dx2)
     d2dx2(N-1,N-1)=-5/(2*dx2)
     d2dx2(N-1,N-1-1)=4/(3*dx2)
     d2dx2(N-1,N-1-2)=-1/(12*dx2)

  case (12)
     ! 5 point stencil, aperiodic:

     ddx(1,1)= -25/(12*dx)
     ddx(1,2)= 4/(dx)
     ddx(1,3)=-3/dx
     ddx(1,4)=4/(3*dx)
     ddx(1,5)=-1/(4*dx)

     ddx(2,1)= -1/(4*dx)
     ddx(2,2)= -5/(6*dx)
     ddx(2,3)=3/(2*dx)
     ddx(2,4)=-1/(2*dx)
     ddx(2,5)=1/(12*dx)

     ddx(N,N)= 25/(12*dx)
     ddx(N,N-1)= -4/(dx)
     ddx(N,N-2)=3/dx
     ddx(N,N-3)=-4/(3*dx)
     ddx(N,N-4)=1/(4*dx)

     ddx(N-1,N)= 1/(4*dx)
     ddx(N-1,N-1)= 5/(6*dx)
     ddx(N-1,N-2)=-3/(2*dx)
     ddx(N-1,N-3)=1/(2*dx)
     ddx(N-1,N-4)=-1/(12*dx)


     d2dx2(1,1)=35/(12*dx2)
     d2dx2(1,2)=-26/(3*dx2)
     d2dx2(1,3)=19/(2*dx2)
     d2dx2(1,4)=-14/(3*dx2)
     d2dx2(1,5)=11/(12*dx2)

     d2dx2(2,1)=11/(12*dx2)
     d2dx2(2,2)=-5/(3*dx2)
     d2dx2(2,3)=1/(2*dx2)
     d2dx2(2,4)=1/(3*dx2)
     d2dx2(2,5)=-1/(12*dx2)

     d2dx2(N,N)=35/(12*dx2)
     d2dx2(N,N-1)=-26/(3*dx2)
     d2dx2(N,N-2)=19/(2*dx2)
     d2dx2(N,N-3)=-14/(3*dx2)
     d2dx2(N,N-4)=11/(12*dx2)

     d2dx2(N-1,N-0)=11/(12*dx2)
     d2dx2(N-1,N-1)=-5/(3*dx2)
     d2dx2(N-1,N-2)=1/(2*dx2)
     d2dx2(N-1,N-3)=1/(3*dx2)
     d2dx2(N-1,N-4)=-1/(12*dx2)

  case (13)
     ! Aperiodic.
     ! 3-point stencil for the first and last rows of the
     ! differentiation matrices, and 4-point stencil for the 2nd and
     ! penultimate rows of the differentiation matrices, so the matrices
     ! are strictly penta-diagonal.

     ddx(1,1)=-1.5/dx
     ddx(1,2)=2/dx
     ddx(1,3)=-0.5/dx

     ddx(N,N)=1.5/dx
     ddx(N,N-1)=-2/dx
     ddx(N,N-2)=0.5/dx

     ddx(2,1)=-1/(3*dx)
     ddx(2,2)=-1/(2*dx)
     ddx(2,3)=1/(dx)
     ddx(2,4)=-1/(6*dx)

     ddx(N-1,N-0)=1/(3*dx)
     ddx(N-1,N-1)=1/(2*dx)
     ddx(N-1,N-2)=-1/(dx)
     ddx(N-1,N-3)=1/(6*dx)

     d2dx2(1,1)=1/dx2
     d2dx2(1,2)=-2/dx2
     d2dx2(1,3)=1/dx2

     d2dx2(N,N)=1/dx2
     d2dx2(N,N-1)=-2/dx2
     d2dx2(N,N-2)=1/dx2

     ! It turns out that the 4-point stencil for the second derivative
     ! has a weight of 0 for the most distant point, making it identical
     ! to the 3-point stencil:

     d2dx2(2,1)=1/(dx2)
     d2dx2(2,2)=-2/(dx2)
     d2dx2(2,3)=1/(dx2)
     d2dx2(2,4)=0

     d2dx2(N-1,N-0)=1/(dx2)
     d2dx2(N-1,N-1)=-2/(dx2)
     d2dx2(N-1,N-2)=1/(dx2)
     d2dx2(N-1,N-3)=0

  case (20,21)
     ! Nothing to be done here

  case (30,31)
     ddx(1,1) = 1/dx
     ddx(1,N) = -1/dx
     d2dx2(1,1) = 1/dx2
     d2dx2(2,2) = 1/dx2
     d2dx2(2,N) = 1/dx2
     d2dx2(1,N-1) = 1/dx2
     d2dx2(1,N) = -2/dx2
     d2dx2(2,1) = -2/dx2

  case (32)
     ! Nothing to be done here

  case (40,41)
     ddx(N,1) = 1/dx
     ddx(N,N) = -1/dx

     d2dx2(N-1,N-1) = 1/dx2
     d2dx2(N,N) = 1/dx2
     d2dx2(N-1,1) = 1/dx2
     d2dx2(N,2) = 1/dx2
     d2dx2(N,1) = -2/dx2
     d2dx2(N-1,N) = -2/dx2

  case (42)
     ! Nothing to be done here

  case (50,51)
     ddx(1,1) = (1.5d+0)/(dx)
     ddx(2,2) = (1.5d+0)/(dx)
     ddx(1,N) = -2/dx
     ddx(2,1) = -2/dx
     ddx(2,N) = 1/(2*dx)
     ddx(1,N-1) = 1/(2*dx)

     d2dx2(1,1) = 1/dx2
     d2dx2(2,2) = 1/dx2
     d2dx2(2,N) = 1/dx2
     d2dx2(1,N-1) = 1/dx2
     d2dx2(1,N) = -2/dx2
     d2dx2(2,1) = -2/dx2

  case (52)
     ddx(2,1) = -1/dx
     ddx(2,2) = 1/dx

  case (60,61)
     ddx(N-1,1) = -1/(2*dx)
     ddx(N,2) = -1/(2*dx)
     ddx(N,1) = 2/dx
     ddx(N-1,N) = 2/dx
     ddx(N-1,N-1) = -(1.5d+0)/dx
     ddx(N,N) = -(1.5d+0)/dx

     d2dx2(N-1,N-1) = 1/dx2
     d2dx2(N,N) = 1/dx2
     d2dx2(N-1,1) = 1/dx2
     d2dx2(N,2) = 1/dx2
     d2dx2(N,1) = -2/dx2
     d2dx2(N-1,N) = -2/dx2

  case (62)
     ddx(N-1,N-1) = -1/dx
     ddx(N-1,N) = 1/dx

  case (80,81,90,91)
     ! Handled previously

  case (82)
     ddx(2,2) =  1/dx
     ddx(2,1) = -1/dx

     ddx(N,N) = (1.5d+0)/dx
     ddx(N,N-1) = -2/dx
     ddx(N,N-2) = 1/(2*dx)

  case (92)
     ddx(N-1,N-1) = -1/dx
     ddx(N-1,N)   =  1/dx

     ddx(1,1) = -(1.5d+0)/dx
     ddx(1,2) = 2/dx
     ddx(1,3) = -1/(2*dx)

  case (100,101,110,111)
     ! Handled previously

  case (102)
     ddx(2,2) =  1/dx
     ddx(2,1) = -1/dx

     ddx(3,4) =  1/(3*dx)
     ddx(3,3)   =  1/(2*dx)
     ddx(3,2) = -1/(dx)
     ddx(3,1) =  1/(6*dx)

!!$     ddx(N,N) = (1.5d+0)/dx
!!$     ddx(N,N-1) = -2/dx
!!$     ddx(N,N-2) = 1/(2*dx)

     ddx(N,N)   =  5/(6*dx)
     ddx(N,N-1) = -3/(2*dx)
     ddx(N,N-2) =  1/(2*dx)
     ddx(N,N-3) = -1/(12*dx)

!!$     do i = 2,N
!!$        ddx(N,i) =  ddx(N-1,i-1)
!!$     end do

  case (112)
     ddx(N-1,N-1) = -1/dx
     ddx(N-1,N)   =  1/dx

     ddx(1,1) = -(1.5d+0)/dx
     ddx(1,2) = 2/dx
     ddx(1,3) = -1/(2*dx)

  case (120,121,130,131)
     ! Handled previously

  case default
     print *,"Error! Invalid value for scheme."
     stop
  end select

end subroutine uniformDiffMatrices

