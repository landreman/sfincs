! This module computes the first and second spectral differentiation matrices.
! Any weight and grid points are allowed.  However, the method is poorly
! conditioned unless the grid points correspond to Gaussian integration
! abscissae for the given weight.  No warning is given in the case of
! poor conditioning.  Presently, the weight function is
! set to exp(-x^2), but you are free to change it.  When changing the
! weight, you must also change the functions dweightdxOverWeight
! and d2weightdx2OverWeight appropriately.
!
! This algorithm closely follows poldif.m for MATLAB
! Part of DMSuite
! written by J.A.C. Weideman, S.C. Reddy 1998
! Available here
! http://www.mathworks.com/matlabcentral/fileexchange/29
! or here
! http://dip.sun.ac.za/~weideman/research/differ.html
!
! The type PetscScalar is used, so this module can be used in a PETSc
! application. However, no other PETSc functionality is used, so you can
! replace the type with e.g. real if you want to build a non-PETSc
! application.
!
! Matt Landreman
! Massachusetts Institute of Technology
! Plasma Science & Fusion Center
! November, 2012
!
module polynomialDiffMatrices

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

  public :: makeXPolynomialDiffMatrices
  private :: weight, dweightdxOverWeight, d2weightdx2OverWeight

contains

  ! ------------------------------------------------------------------

  function weight(x)

    use xGrid, only: xGrid_k

    implicit none

    PetscScalar, intent(in) :: x
    PetscScalar :: weight

    weight = exp(-x*x)*(x ** xGrid_k)

  end function weight

  ! ------------------------------------------------------------------

  function dweightdxOverWeight(x)

    use xGrid, only: xGrid_k

    implicit none

    PetscScalar, intent(in) :: x
    PetscScalar :: dweightdxOverWeight

    if (abs(x)<1e-12) then
       ! Handle possible point at x=0, avoiding divide-by-0:
       dweightdxOverWeight = 0.0d+0
    else
       !dweightdxOverWeight = -2*x
       dweightdxOverWeight = xGrid_k / x - 2*x
    end if

  end function dweightdxOverWeight

  ! ------------------------------------------------------------------

  function d2weightdx2OverWeight(x)

    use xGrid, only: xGrid_k

    implicit none

    PetscScalar, intent(in) :: x
    PetscScalar :: d2weightdx2OverWeight

    if (abs(x)<1e-12) then
       ! Handle possible point at x=0, avoiding divide-by-0:
       d2weightdx2OverWeight = -2.0d+0
    else
       !d2weightdx2OverWeight = -2 + 4*x*x
       d2weightdx2OverWeight = xGrid_k*(xGrid_k-1)/(x*x)-2*(2*xGrid_k+1)+4*x*x
    end if

  end function d2weightdx2OverWeight

  ! ------------------------------------------------------------------

  subroutine makeXPolynomialDiffMatrices(x,ddx,d2dx2)
    ! This is the main function of the module.
    !
    ! Inputs:
    !  x = abscissae (grid point locations)
    !
    ! Outputs:
    !   ddx = matrix for the first derivative
    !   d2dx2 = matrix for the second derivative
    !
    ! Both ddx and d2dx should be pre-allocated as 2D arrays of size N x N,
    ! where N is the size of the 1D array x.

    implicit none

    PetscScalar, intent(in), dimension(:) :: x
    PetscScalar, intent(out) :: ddx(:,:), d2dx2(:,:)
    integer :: i, N
    PetscScalar, dimension(:,:), allocatable :: XX, DX, CC, CCC, Z, XXX, Y, D, oldY, repmatDiagD
    PetscScalar, dimension(:), allocatable :: c
    ! XXX in this subroutine corresponds to X in the DMSuite version of poldif

    N = size(x)
    allocate(Z(N,N))
    allocate(XX(N,N))
    allocate(CC(N,N))
    allocate(CCC(N,N))
    allocate(DX(N,N))
    allocate(c(N))
    do i=1,N
       XX(i,:)=x(i)
    end do
    DX = XX - transpose(XX)
    do i=1,N
       DX(i,i)=1d+0
    end do


    c = product(DX,2)


    do i=1,N
       c(i) = c(i)*weight(x(i))
       CC(i,:) = c(i)
    end do
    CCC = CC / transpose(CC)


    Z = 1 / DX
    do i=1,N
       Z(i,i)=0d+0
    end do
    allocate(XXX(N-1,N))
    XXX=0d+0
    do i=1,N
       XXX(i:N-1, i) = Z(i, i+1:N)
       XXX(1:i-1, i) = Z(i, 1:i-1)
    end do

    allocate(Y(N,N))
    allocate(oldY(N,N))
    ! ell = 1 (i.e. first derivative)
    do i=1,N
       Y(1,i) = dweightdxOverWeight(x(i))
    end do
    do i=2,N
       Y(i,:) = Y(i-1,:) + XXX(i-1,:)
    end do

    ddx = Z * CCC
    do i=1,N
       ddx(i,i) = Y(N,i)
    end do

    ! ell = 2 (i.e. second derivative)
    oldY = Y
    do i=1,N
       Y(1,i) = d2weightdx2OverWeight(x(i))
    end do
    do i=2,N
       Y(i,:) = Y(i-1,:) + 2*oldY(i-1,:)*XXX(i-1,:)
    end do
    allocate(repmatDiagD(N,N))
    do i=1,N
       repmatDiagD(i,:) = ddx(i,i)
    end do
    d2dx2 = 2*Z*(CCC*repmatDiagD-ddx)
    do i=1,N
       d2dx2(i,i) = Y(N,i)
    end do

  end subroutine makeXPolynomialDiffMatrices

end module polynomialDiffMatrices
