#include <finclude/petscsnesdef.h>


  subroutine evaluateResidual(mysnes, stateVec, residualVec, userContext, ierr)

!    use petscmat
    use petscsnes

    implicit none

    SNES :: mysnes
    Vec :: stateVec, residualVec
    PetscErrorCode :: ierr
    integer :: userContext(*)
    PetscScalar :: x, y
    PetscScalar, pointer, dimension(:) :: stateArray, residualArray

    call VecGetArrayF90(stateVec,stateArray,ierr)
    call VecGetArrayF90(residualVec,residualArray,ierr)


    x = stateArray(1)
    y = stateArray(2)

    residualArray(1) = x*y - x*x + y + 1
    residualArray(2) = exp(x+y-1) - 1

    print *,"FormFunction called: x=",x,", y=",y, &
         ", residual(1)=",residualArray(1),", residual(2)=",residualArray(2)

    call VecRestoreArrayF90(stateVec,stateArray,ierr)
    call VecRestoreArrayF90(residualVec,residualArray,ierr)

  end subroutine evaluateResidual
