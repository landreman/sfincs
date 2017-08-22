subroutine ChebyshevGrid(N, xMin, xMax, x, weights, D)
  ! Creates a Chebyshev grid, together with the associated weights
  ! for Clenshaw-Curtis integration, and the spectral
  ! differentiation matrix.
  ! Based on the Matlab function 
  ! m20120124_06_multiChebyshevWeightsAndDifferentiation.m
  ! 
  ! Created by Matt Landreman, 
  ! Massachusetts Institute of Technology, Plasma Science & Fusion Center, 2012.
  !
  ! Inputs:
  ! N = number of grid points desired.
  ! xMin = minimum of interval
  ! xMax = maximum of interval
  !
  ! Outputs:
  ! x(N) = grid points (abscissae)
  ! weights(N) = Clenshaw-Curtis integration weights
  ! D(N,N) = differentiation matrix

  use kinds

  implicit none

  integer, intent(in) :: N
  real(prec), intent(in) :: xMin, xMax
  real(prec), intent(out) :: x(N), weights(N), D(N,N)

  integer :: i, j, N1, M
  real(prec), parameter :: pi = 3.1415926535897932384626433d+0
  real(prec), allocatable :: c(:), bigX(:,:), dX(:,:), sumD(:), cw(:), cc(:), f(:)

  if (xMax <= xMin) then
     print *,"Error! xMax should be larger than xMin."
     stop
  end if

  if (N<1) then
     print *,"Error! N must be at least 1."
     stop
  end if

  ! ****************************************************
  ! First, build abscissae and differentiation matrix
  ! ****************************************************

  N1 = N-1
  x = [( cos(pi*i/N1), i=0,N1 )]

  allocate(c(N))
  allocate(sumD(N))
  allocate(bigX(N,N))
  allocate(dX(N,N))

  c = 1
  c(1) = 2
  c(N) = 2
  do i=2, N, 2
     c(i) = -c(i)
  end do

  do i=1,N
     bigX(i,:)=x(i)
  end do

  dX = bigX - transpose(bigX)

  do i=1,N
     dX(i,i) = 1
  end do

  do i=1,N
     do j=1,N
        D(i,j) = c(i)/c(j)
     end do
  end do
  D = D / dX

  sumD = sum(D,2)
  do i=1,N
     D(i,i) = D(i,i) - sumD(i)
  end do

  D = -D * 2/(xMax-xMin)

  x = (1-x) * (xMax-xMin)/2 + xMin

  ! ****************************************************
  ! Done building absicssae and differentiation matrix.
  ! Next, build the integration weights:
  ! ****************************************************

  allocate(cw(N))
  cw = 0
  cw(1) = 2
  do i=3, N, 2
     cw(i) = 2 / ((1d+0) - (i-1)*(i-1))
  end do

  M = N*2-2
  allocate(cc(M))
  cc(1:N) = cw
  cc((N+1):M) = cw((N-1):2:(-1))

  ! Next, we need to form real(ifft(cc)).
  ! Rather than use a library for this, I'll just compute
  ! it in a direct (but slow) way, since speed
  ! is not essential here.

  allocate(f(M))
  f = 0
  do i=1,M
     do j=1,M
        f(i) = f(i) + cc(j)*cos(2*pi*(i-1)*(j-1)/M)
     end do
  end do
  f = f/M

  weights = f
  weights(1) = f(1)/2
  weights(N) = f(N)/2
  weights = weights * (xMax-xMin)

end subroutine ChebyshevGrid

