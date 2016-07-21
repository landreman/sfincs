subroutine FourierDifferentiationMatrices(NFourier, xm, xn, ddtheta, ddzeta)

  use globalVariables, only: prec

  ! ddtheta and ddzeta should have been allocated with size (2*NFourier-1, 2*NFourier-1)

  implicit none

  integer, intent(in) :: NFourier
  integer, intent(in), dimension(NFourier) :: xm, xn
  real(prec), intent(out), dimension(NFourier*2-1,NFourier*2-1) :: ddtheta, ddzeta

  integer :: j, NFourier2

  NFourier2 = NFourier*2-1

  ! ***************************************************************
  ! Validate input
  ! ***************************************************************
  
  if (size(xm) .ne. NFourier) stop "Size of xm is not correct"
  if (size(xn) .ne. NFourier) stop "Size of xn is not correct"
  if (xm(1) .ne. 0) stop "xm(1) should be 0."
  if (xn(1) .ne. 0) stop "xn(1) should be 0."
  if (size(ddtheta,1) .ne. NFourier*2-1) stop "Size of first dimension of ddtheta is not correct"
  if (size(ddtheta,2) .ne. NFourier*2-1) stop "Size of second dimension of ddtheta is not correct"
  if (size(ddzeta, 1) .ne. NFourier*2-1) stop "Size of first dimension of ddzeta is not correct"
  if (size(ddzeta, 2) .ne. NFourier*2-1) stop "Size of second dimension of ddzeta is not correct"

  ! ***************************************************************
  ! Assemble matrices
  ! ***************************************************************

  ddtheta=0.0
  ddzeta=0.0

  do j=2,NFourier
     ! (d/dtheta) sin(m theta - n zeta) =  m cos(m theta - n zeta):
     ddtheta(j, j+NFourier-1) =  xm(j)
     ! (d/dtheta) cos(m theta - n zeta) = -m sin(m theta - n zeta):
     ddtheta(j+NFourier-1, j) = -xm(j)
     ! (d/dzeta) sin(m theta - n zeta) =  -n cos(m theta - n zeta):
     ddzeta(j, j+NFourier-1) = -xn(j)
     ! (d/dzeta) cos(m theta - n zeta) =   n sin(m theta - n zeta):
     ddzeta(j+NFourier-1, j) =  xn(j)
  end do

!!$  print *,"ddtheta:"
!!$  do j=1,NFourier2
!!$     print *,ddtheta(j,:)
!!$  end do
!!$  print *,"ddzeta:"
!!$  do j=1,NFourier2
!!$     print *,ddzeta(j,:)
!!$  end do

end subroutine FourierDifferentiationMatrices
