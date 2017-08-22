#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif

#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
       ! Syntax for PETSc versions up through 3.4
  subroutine evaluateJacobian(mysnes, stateVec, jacobian, jacobianPC, flag, userContext, ierr)
#else
       ! Syntax for PETSc version 3.5 and later
  subroutine evaluateJacobian(mysnes, stateVec, jacobian, jacobianPC, userContext, ierr)
#endif

    use petscsnes
    use globalVariables, only: masterProc, useIterativeLinearSolver, firstMatrixCreation, reusePreconditioner

    implicit none

    SNES :: mysnes
    Vec :: stateVec
    Mat :: jacobian, jacobianPC
    PetscErrorCode :: ierr
    integer :: userContext(*)
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
    ! Syntax for PETSc versions up through 3.4
    MatStructure :: flag
#endif

    if (masterProc) then
       print *,"evaluateJacobian called."
    end if

#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
    ! Syntax for PETSc versions up through 3.4.
    ! For PETSc versions 3.5 and later, the code related to 'reusePreconditioner' is implemented in solver.F90.
    if (reusePreconditioner) then
       flag = SAME_PRECONDITIONER
    else
       flag = DIFFERENT_NONZERO_PATTERN
    end if
#endif

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
