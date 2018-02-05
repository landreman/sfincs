#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif
  
subroutine polynomialInterpolationMatrix(N, M, xk, x, alpxk, alpx, matrix)
  ! This module returns the matrix for spectral interpolation from any grid to
  ! any other grid, using any weights. However, the method is poorly
  ! conditioned unless the points of the initial grid (on which the function
  ! is already known) correspond to Gaussian integration
  ! abscissae for the given weight.  No warning is given in the case of
  ! poor conditioning.  The grid points onto which you interpolate need not
  ! be nicely distributed.  The weight function does not appear in this module;
  ! rather, the weight evaluated on both grids is specified as an input.
  !
  ! This algorithm closely follows polint.m for MATLAB
  ! Part of DMSuite
  ! written by J.A.C. Weideman, S.C. Reddy 1998
  ! Available here
  ! http://www.mathworks.com/matlabcentral/fileexchange/29
  ! or here
  ! http://dip.sun.ac.za/~weideman/research/differ.html
  !
  ! Note from DMSuite's polint function:
  !     The code implements the barycentric formula; see page 252 in
  !     P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
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

  ! Inputs:
  !   N = number of grid points on which the function is known.
  !   M = number of grid points onto which we would like to interpolate.
  !   xk(N) = grid points on which the function is known.
  !   alpxk(N) = weight function evaluated on the xk grid.
  !   x(M) = grid onto which we would like to interpolate.
  !   alpx(M) = weight function evaluated on the x grid.
  !
  ! Outputs:
  !   matrix(M,N) = interpolation matrix.
  !
  ! xk and alpxk should be 1D arrays with the same size. (Call it N.)
  ! x and alpx should be 1D arrays with the same size. (Call it M.)
  ! matrix should be preallocated with size (M rows) x (N columns).

  implicit none

  integer, intent(in) :: N, M
  PetscScalar, dimension(N), intent(in) :: xk, alpxk
  PetscScalar, dimension(M), intent(in) :: x, alpx
  PetscScalar, intent(out) :: matrix(M,N)
  integer :: i
  PetscScalar, dimension(:,:), allocatable :: xkMatrixified, D, xMatrixified
  PetscScalar, dimension(:), allocatable :: w

  allocate(xkMatrixified(N,N))
  allocate(D(N,N))
  allocate(w(N))
  do i=1,N
     xkMatrixified(i,:)=xk(i)
  end do
  D = xkMatrixified - transpose(xkMatrixified)
  do i=1,N
     D(i,i)=1d+0
  end do

  w = 1/product(D,1)

  allocate(xMatrixified(M,N))
  do i=1,M
     xMatrixified(i,:) = x(i)
  end do
  deallocate(xkMatrixified)
  allocate(xkMatrixified(N,M))
  do i=1,N
     xkMatrixified(i,:) = xk(i)
  end do
  deallocate(D)
  allocate(D(M,N))
  D = xMatrixified - transpose(xkMatrixified)

  where (D==0e+0)
     D=1d-15
  end where
  D=1/D

  matrix = D

  do i=1,M
     matrix(i,:) = matrix(i,:) &
          & * alpx(i)/sum(D(i,:)*w)
  end do

  do i=1,N
     matrix(:,i) = matrix(:,i) &
          & * w(i)/alpxk(i)
  end do

end subroutine polynomialInterpolationMatrix

