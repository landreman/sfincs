
subroutine testingAdjointOperator()

  use globalVariables

  implicit none

  Vec :: dummy
  Mat :: adjointMatrix

  ! Populate adjoint operator
  call preallocateMatrix(adjointMatrix, 1)
  call populateMatrix(adjointMatrix, 4, dummy)

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointRHSVec1, ierr)
  call VecSet(adjointRHSVec1, zero, ierr)

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointRHSVec2, ierr)
  call VecSet(adjointRHSVec2, zero, ierr)

  !> Populate RHS vec
  call populateAdjointRHS(adjointRHSVec1, 1, 1)

  ! Populate dMatrixdLambda
!  call preallocateMatrix(dMatrixdLambda, 1)
!  call populatedMatrixdLambda(dMatrixdLambda,whichLambda,whichMode)
!
!  ! Populate dRHSdLambda
!  call VecCreateMPI(MPIComm,PETSC_DECIDE, matrixSize, dRHSdLambda,ierr)
!  call populatedRHSdLambda(dRHSdLambda, whichLambda, whichMode)

  ! dL/dlambda F - dS/dlambda


end subroutine testingAdjointOperator
