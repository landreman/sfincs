module sparsify

  ! This module provides wrappers for Petsc's MatSetValue and MatSetValues subroutines.
  ! The wrappers act exactly the same as the original Petsc subroutines, except that
  ! values are only added when they exceed a small positive threshhold. This threshholding
  ! prevents zero entries from being added to the nonzero structure of the matrix, increasing
  ! sparsity.

  use kinds

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

  PetscScalar :: threshholdForInclusion = 1d-12

contains

  subroutine MatSetValueSparse(myMat, row, col, valueToSet, mode, err)

    implicit none

    Mat :: myMat
    integer :: row, col
    InsertMode :: mode
    real(prec) :: valueToSet
    PetscScalar :: valueToSet_PetscScalar
    PetscErrorCode :: err

    if (abs(valueToSet) > threshholdForInclusion) then
       valueToSet_PetscScalar = valueToSet ! Cast type from double to single if needed. 
       call MatSetValue(myMat, row, col, valueToSet_PetscScalar, mode, err)
    end if

  end subroutine MatSetValueSparse


  subroutine MatSetValuesSparse(myMat, m, idxm, n, idxn, v, mode, err)

    implicit none

    Mat :: myMat
    integer :: m, n, row, col
    integer :: idxm(m), idxn(n)
    real(prec) :: v(n,m)
    PetscScalar :: valueToSet
    InsertMode :: mode
    PetscErrorCode :: err

    do row=1,m
       do col=1,n
          valueToSet = v(col,row) ! I'll use PETSc's ordering instead of Fortran's.

          if (abs(valueToSet) > threshholdForInclusion) then
             call MatSetValue(myMat, idxm(row), idxn(col), valueToSet, mode, err)
          end if
       end do
    end do

  end subroutine MatSetValuesSparse

end module sparsify
