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

subroutine testingAdjointOperator(forwardSolution,adjointSolution,whichAdjointRHS,whichSpecies)

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
  Vec :: solutionVec, result4, adjointRHS2
  integer :: iSpecies, iAdjointRHS

  if (masterProc) then
    print *,"ispecies = ", whichSpecies
    print *,"whichAdjointRHS = ", whichAdjointRHS
    print *,"Testing adjoint operator by ensuring that (forwardSolution,matrix*forwardSolution) = (adjointMatrix*forwardSolution,forwardSolution)"
  end if

  ! Populate forward operator
  call preallocateMatrix(matrix, 1)
  call populateMatrix(matrix, 6, dummy)

  ! Populate adjoint operator
  call preallocateMatrix(adjointMatrix, 1)
  if (discreteAdjointOption) then
    call MatTranspose(matrix, MAT_INPLACE_MATRIX,adjointMatrix,ierr)
  else
    call populateMatrix(adjointMatrix, 4, dummy)
  end if

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

  call MatMult(matrix,adjointRHS,result3,ierr)

  call innerProduct(forwardSolution,result3,innerProduct2)
  if (masterProc) then
    print *,"Inner product 1: ", innerProduct1
    print *,"Inner product 2: ", innerProduct2
  end if

  if (masterproc) then
    print *,"Testing adjoint operator by ensuring that (matrix*adjointRHS,adjointRHS2) = (adjointRHS,adjointMatrix*adjointRHS2) "
  end if

  ! Create result3 and set to zero
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, result4, ierr)
  call VecSet(result4, zero, ierr)

  ! Allocate adjointRHS
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointRHS2, ierr)
  do iAdjointRHS = 1,3
    do iSpecies = 0,NSpecies
      call VecSet(adjointRHS2, zero, ierr)
      call populateAdjointRHS(adjointRHS2, iAdjointRHS, iSpecies)
      call VecSet(result4,zero,ierr)
      call MatMult(adjointMatrix,adjointRHS2,result4,ierr)
      call innerProduct(adjointRHS,result4,innerProduct1)
      call innerProduct(result3,adjointRHS2,innerProduct2)
      if (masterProc) then
        print *,"whichAdjointRHS = ", iAdjointRHS
        print *,"whichSpecies = ", iSpecies
        print *,"innerProduct1 = ", innerProduct1
        print *,"innerProduct2 = ", innerProduct2
      end if
    end do
  end do

  ! ((matrix)^{-1} RHS, adjointRHS) = (RHS, (adjointMatrix)^{-1} adjointRHS)
!  if (masterProc) then
!    print *,"Testing adjoint operator by ensuring that (matrix^{-1}RHS, adjointRHS) = (RHS, adjointMatrix^{-1}adjointRHS) "
!   ! print "(a,i4,a,i4,a)","Using whichAdjointRHS = ", whichAdjointRHS, &
!    !" and whichSpecies = ", whichSpecies," -----------------------------"
!  end if

  ! Copy adjointSolution
!  call VecDuplicate(adjointSolution,adjointSolution1,ierr)
!  call VecCopy(adjointSolution,adjointSolution1,ierr)

  ! Allocate and populate RHS 
  ! residualVec may be overwritte by call to evaluateResidual?
!  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, RHS, ierr)
!  call VecSet(RHS, zero, ierr)
!  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, solutionVec, ierr)
!  call VecSet(solutionVec, zero, ierr)
!  call evaluateResidual(mysnes, solutionVec, RHS, userContext, ierr)
!  ! Multiply the residual by (-1):
!  call VecScale(RHS, -1.d+0, ierr)

  ! Copy RHS
!  call VecDuplicate(RHS,RHS1,ierr)
!  call VecCopy(RHS,RHS1,ierr)

!  call innerProduct(adjointSolution,RHS, innerProduct1)
!  call innerProduct(forwardSolution,adjointRHS, innerProduct2)
!  if (masterProc) then
!    print *,"Inner product 1: ", innerProduct1
!    print *,"Inner product 2: ", innerProduct2
!  end if

end subroutine testingAdjointOperator
