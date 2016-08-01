subroutine FourierDifferentiationMatrices(NFourier2, A, B, C)

  ! This subroutine forms C = A * d/dtheta + B * d/dzeta

  use globalVariables, only: prec, NFourier, xm, xn

  ! A, B, and C should have been allocated with size (2*NFourier-1, 2*NFourier-1)

  implicit none

  integer, intent(in) :: NFourier2
  real(prec), intent(in), dimension(NFourier2,NFourier2) :: A, B
  real(prec), intent(out), dimension(NFourier2,NFourier2) :: C

  integer :: i,j

  ! ***************************************************************
  ! Assemble matrices
  ! ***************************************************************

  C(:,1)=0
  
  do j=2,NFourier
     do i=1,NFourier2
        ! (d/dtheta) sin(m theta - n zeta) =  m cos(m theta - n zeta):
        ! (d/dzeta) sin(m theta - n zeta) =  -n cos(m theta - n zeta):
        !ddtheta(j, j+NFourier-1) =  xm(j)
        !ddzeta(j, j+NFourier-1) = -xn(j)
        C(i,j+NFourier-1) = A(i,j)*xm(j) - B(i,j)*xn(j)

        ! (d/dtheta) cos(m theta - n zeta) = -m sin(m theta - n zeta):
        ! (d/dzeta) cos(m theta - n zeta) =   n sin(m theta - n zeta):
        !ddtheta(j+NFourier-1, j) = -xm(j)
        !ddzeta(j+NFourier-1, j) =  xn(j)
        C(i,j) = -A(i,j+NFourier-1)*xm(j) + B(i,j+NFourier-1)*xn(j)
     end do
  end do

!!$  ddtheta=0.0
!!$  ddzeta=0.0
!!$
!!$  do j=2,NFourier
!!$     ! (d/dtheta) sin(m theta - n zeta) =  m cos(m theta - n zeta):
!!$     ddtheta(j, j+NFourier-1) =  xm(j)
!!$     ! (d/dtheta) cos(m theta - n zeta) = -m sin(m theta - n zeta):
!!$     ddtheta(j+NFourier-1, j) = -xm(j)
!!$     ! (d/dzeta) sin(m theta - n zeta) =  -n cos(m theta - n zeta):
!!$     ddzeta(j, j+NFourier-1) = -xn(j)
!!$     ! (d/dzeta) cos(m theta - n zeta) =   n sin(m theta - n zeta):
!!$     ddzeta(j+NFourier-1, j) =  xn(j)
!!$  end do


end subroutine FourierDifferentiationMatrices


subroutine FourierDerivative(NFourier2, U, V, whichCoordinate)

  ! This subroutine forms 
  ! U = (d/dtheta) * V (if whichCoordinate==1)
  ! or
  ! U =  (d/dzeta) * V  (if whichCoordinate==2)


  use globalVariables, only: prec, xm, xn, NFourier

  ! U and V should have been allocated with size (NFourier2)

  implicit none

  integer, intent(in) :: NFourier2, whichCoordinate
  real(prec), intent(in), dimension(NFourier2) :: U
  real(prec), intent(out), dimension(NFourier2) :: V

  integer :: j

  ! ***************************************************************
  ! Assemble matrices
  ! ***************************************************************

  V(1)=0
 
  select case (whichCoordinate)
  case (1)
     do j=2,NFourier
        ! (d/dtheta) sin(m theta - n zeta) =  m cos(m theta - n zeta):
        !ddtheta(j, j+NFourier-1) =  xm(j)
        V(j) = xm(j) * U(j+NFourier-1)

        ! (d/dtheta) cos(m theta - n zeta) = -m sin(m theta - n zeta):
        !ddtheta(j+NFourier-1, j) = -xm(j)
        V(j+NFourier-1) = -xm(j) * U(j)
     end do

  case (2)
     do j=2,NFourier
        ! (d/dzeta) sin(m theta - n zeta) =  -n cos(m theta - n zeta):
        !ddzeta(j, j+NFourier-1) = -xn(j)
        V(j) = -xn(j) * U(j+NFourier-1)

        ! (d/dzeta) cos(m theta - n zeta) =   n sin(m theta - n zeta):
        !ddzeta(j+NFourier-1, j) =  xn(j)
        V(j+NFourier-1) = xn(j) * U(j)
     end do

  case default
     print *,"Invalid whichCoordinate:",whichCoordinate
     stop
  end select

!!$  ddtheta=0.0
!!$  ddzeta=0.0
!!$
!!$  do j=2,NFourier
!!$     ! (d/dtheta) sin(m theta - n zeta) =  m cos(m theta - n zeta):
!!$     ddtheta(j, j+NFourier-1) =  xm(j)
!!$     ! (d/dtheta) cos(m theta - n zeta) = -m sin(m theta - n zeta):
!!$     ddtheta(j+NFourier-1, j) = -xm(j)
!!$     ! (d/dzeta) sin(m theta - n zeta) =  -n cos(m theta - n zeta):
!!$     ddzeta(j, j+NFourier-1) = -xn(j)
!!$     ! (d/dzeta) cos(m theta - n zeta) =   n sin(m theta - n zeta):
!!$     ddzeta(j+NFourier-1, j) =  xn(j)
!!$  end do


end subroutine FourierDerivative
