subroutine legendrePolynomials(P, xi, Nxi, NL)

#include "PETScVersions.F90"

  implicit none

  integer, intent(in) :: NL, Nxi
  PetscScalar, dimension(Nxi,NL), intent(out) :: P
  PetscScalar, intent(in), dimension(Nxi) :: xi
  integer :: L
  PetscScalar :: one = 1.0

  L=0
  P(:,L+1) = one  ! P_0(xi) = 1
  
  if (NL>1) then
     L=1
     P(:,L+1) = xi  ! P_1(xi) = xi
  end if
  
  do L=2,NL-1
     ! Recursion relation for Legendre polynomials:
     P(:,L+1) = ((2*L-1)* xi * P(:,L) &
          - (L-1)*P(:,L-1))/(one*L)
       
  end do
end subroutine legendrePolynomials

