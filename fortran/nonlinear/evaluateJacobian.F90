#include <finclude/petscsnesdef.h>
#include "PETScVersions.F90"


  subroutine evaluateJacobian(mysnes, stateVec, jacobian, jacobianPC, userContext, ierr)

!    use petscmat
    use petscsnes

    implicit none

    SNES :: mysnes
    Vec :: stateVec
    Mat :: jacobian, jacobianPC
    PetscErrorCode :: ierr
    integer :: userContext(*)

!!$    PetscScalar :: x, y
!!$    PetscInt :: rowIndices(1), colIndices(2)
!!$    PetscScalar, pointer, dimension(:) :: stateArray
!!$    PetscScalar :: values(2)
!!$    PetscInt :: matrixSize=2
!!$
!!$    call VecGetArrayF90(stateVec,stateArray,ierr)
!!$
!!$    x = stateArray(1)
!!$    y = stateArray(2)
!!$    print *,"formJacobian called, x=",x,", y=",y
!!$
!!$    colIndices(1) = 0
!!$    colIndices(2) = 1
!!$
!!$    ! Set first row
!!$    rowIndices = 0
!!$    values(1) = y - 2*x
!!$    values(2) = x + 1
!!$    call MatSetValues(jacobian, 1, rowIndices, matrixSize, colIndices, values, INSERT_VALUES, ierr)
!!$
!!$    ! Set second row
!!$    rowIndices = 1
!!$    values(1) = exp(x+y-1)
!!$    values(2) = exp(x+y-1)
!!$    call MatSetValues(jacobian, 1, rowIndices, matrixSize, colIndices, values, INSERT_VALUES, ierr)
!!$
!!$    call VecRestoreArrayF90(stateVec, stateArray, ierr)
!!$
!!$    call MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY, ierr)
!!$    call MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY, ierr)

    call populateMatrix(jacobianPC, 0)
    call populateMatrix(jacobian, 1)

  end subroutine evaluateJacobian
