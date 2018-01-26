#include "PETSCVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif
#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif

subroutine testingAdjointOperator(forwardSolution,adjointSolution,RHS,adjointRHS,matrix,adjointMatrix,whichAdjointRHS,whichSpecies)

  use globalVariables
  use petscmat
  use adjointDiagnostics
  use petscsnes

  implicit none

  Vec :: forwardSolution, adjointSolution, RHS, adjointRHS
  integer :: whichAdjointRHS, whichSpecies
  Vec :: dummy, result
  Mat :: adjointMatrix, matrix
  PetscErrorCode :: ierr
  Vec :: forwardSolution1, forwardSolution2, forwardSolution3
  PetscScalar :: norm
  Vec :: result1, result2, result3
  PetscScalar :: innerProduct1, innerProduct2
  Vec :: adjointSolution2, RHS2, RHS1, adjointSolution1, adjointRHS1
  SNES :: mysnes
  integer :: userContext(1)
  Vec :: solutionVec

  if (masterProc) then
    print *,"Testing adjoint operator by ensuring that (forwardSolution,matrix*forwardSolution) = (adjointMatrix*forwardSolution,forwardSolution)"
  end if

  ! Copy forwardSolution
!  call VecDuplicate(forwardSolution,forwardSolution1,ierr)
!  call VecCopy(forwardSolution,forwardSolution1,ierr)
!  call VecDuplicate(forwardSolution,forwardSolution2,ierr)
!  call VecCopy(forwardSolution,forwardSolution2,ierr)
!  call VecDuplicate(forwardSolution,forwardSolution3,ierr)
!  call VecCopy(forwardSolution,forwardSolution3,ierr)

!  ! Populate adjoint operator
  call preallocateMatrix(adjointMatrix, 1)
  call populateMatrix(adjointMatrix, 4, dummy)

!  ! Populate forward operator
  call preallocateMatrix(matrix, 1)
  call populateMatrix(matrix, 1, dummy)

  ! Create result1 and set to zero
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, result1, ierr)
  call VecSet(result1, zero, ierr)

  ! Create result2 and set to zero
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, result2, ierr)
  call VecSet(result2, zero, ierr)

  ! result1 = (matrix)(forwardSolution1)
  call MatMult(matrix,forwardSolution,result1,ierr)

  ! result2 = (adjointMatrix)(forwardSolution2)
  call MatMult(adjointMatrix,forwardSolution,result2,ierr)

  ! Computer inner products
  call innerProduct(forwardSolution,result1,innerProduct1)
  call innerProduct(forwardSolution,result2,innerProduct2)
  if (masterproc) then
    print *,"Inner product 1: ", innerProduct1
    print *,"Inner product 2: ", innerProduct2

    print *,"Testing adjoint operator by ensuring that (matrix*adjointRHS,forwardSolution) = (adjointRHS,adjointMatrix*forwardSolution) "
    !print "(a,i4,a,i4,a)","Using whichAdjointRHS = ", whichAdjointRHS, &
            !    " and whichSpecies = ", whichSpecies," -----------------------------"
  end if

  ! Create result3 and set to zero
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, result3, ierr)
  call VecSet(result3, zero, ierr)

  ! Allocate adjointRHS
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointRHS, ierr)
  call VecSet(adjointRHS, zero, ierr)
  call populateAdjointRHS(adjointRHS, whichAdjointRHS, whichSpecies)

  ! Compute inner product
  call innerProduct(adjointRHS, result2, innerProduct1)

  ! Create adjointRHS1
!  call VecDuplicate(adjointRHS,adjointRHS,ierr)
!  call VecCopy(adjointRHS,adjointRHS,ierr)

  ! result3 = (matrix)(adjointRHS)

!  call VecNorm(adjointRHS,2,norm,ierr)
!  if (masterProc) then
!    print *,"norm before mult: ",norm
!  end if
  call MatMult(matrix,adjointRHS,result3,ierr)
!  call VecNorm(adjointRHS,2,norm,ierr)
!  if (masterProc) then
!    print *,"norm after mult: ",norm
!  end if

  call innerProduct(forwardSolution,result3,innerProduct2)
  if (masterProc) then
    print *,"Inner product 1: ", innerProduct1
    print *,"Inner product 2: ", innerProduct2
  end if

  ! ((matrix)^{-1} RHS, adjointRHS) = (RHS, (adjointMatrix)^{-1} adjointRHS)
  if (masterProc) then
    print *,"Testing adjoint operator by ensuring that (matrix^{-1}RHS, adjointRHS) = (RHS, adjointMatrix^{-1}adjointRHS) "
   ! print "(a,i4,a,i4,a)","Using whichAdjointRHS = ", whichAdjointRHS, &
    !" and whichSpecies = ", whichSpecies," -----------------------------"
  end if

  ! Copy adjointSolution
!  call VecDuplicate(adjointSolution,adjointSolution1,ierr)
!  call VecCopy(adjointSolution,adjointSolution1,ierr)

  ! Allocate and populate RHS
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, RHS, ierr)
  call VecSet(RHS, zero, ierr)
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, solutionVec, ierr)
  call VecSet(solutionVec, zero, ierr)
  call evaluateResidual(mysnes, solutionVec, RHS, userContext, ierr)
  ! Multiply the residual by (-1):
  call VecScale(RHS, -1d+0, ierr)

  ! Copy RHS
!  call VecDuplicate(RHS,RHS1,ierr)
!  call VecCopy(RHS,RHS1,ierr)

  call innerProduct(adjointSolution,RHS, innerProduct1)
  call innerProduct(forwardSolution,adjointRHS, innerProduct2)
  if (masterProc) then
    print *,"Inner product 1: ", innerProduct1
    print *,"Inner product 2: ", innerProduct2
  end if

end subroutine testingAdjointOperator
