#include <finclude/petscsnesdef.h>
#include "PETScVersions.F90"


  subroutine evaluateJacobian(mysnes, stateVec, jacobian, jacobianPC, userContext, ierr)

    use petscsnes
    use globalVariables, only: masterProc, useIterativeLinearSolver, firstMatrixCreation, reusePreconditioner

    implicit none

    SNES :: mysnes
    Vec :: stateVec
    Mat :: jacobian, jacobianPC
    PetscErrorCode :: ierr
    integer :: userContext(*)

    if (masterProc) then
       print *,"evaluateJacobian called."
    end if

    ! When PETSc assembles a matrix, it reduces the structure of nonzeros to the actual number of nonzeros.
    ! If we try to re-assemble the matrix with additional nonzero entries without first re-allocating space for the nonzeros,
    ! we get the error about 'new nonzero caused a malloc'. Therefore, here we destroy the matrices and reallocate them.

    if (useIterativeLinearSolver) then
       ! If reusePreconditioner = true, then we only need to assemble the preconditioner in the first iteration.
       ! If reusePreconditioner = false, then we need to assemble the preconditioner in every iteration.
       if (firstMatrixCreation .or. .not. reusePreconditioner) then
          call populateMatrix(jacobianPC, 0, stateVec)
       end if
    end if
    call populateMatrix(jacobian, 1, stateVec)

    firstMatrixCreation = .false.

  end subroutine evaluateJacobian
