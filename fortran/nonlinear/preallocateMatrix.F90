#include <finclude/petscmatdef.h>
#include "PETScVersions.F90"

subroutine preallocateMatrix(matrix, whichMatrix)

  use petscmat
  use globalVariables, only: Nx, Nxi, Ntheta, Nzeta, Nspecies, matrixSize, &
       constraintScheme, PETSCPreallocationStrategy, MPIComm, numProcs, masterProc

  implicit none

  integer, intent(in) :: whichMatrix
  Mat :: matrix
  integer :: predictedNNZForEachRowOfPreconditioner, predictedNNZForEachRowOfTotalMatrix
  integer, dimension(:), allocatable :: predictedNNZsForEachRow, predictedNNZsForEachRowDiagonal
  PetscErrorCode :: ierr
  integer :: tempInt1, i
  integer :: firstRowThisProcOwns, lastRowThisProcOwns, numLocalRows

  if (masterProc) then
     print *,"Beginning preallocation for whichMatrix = ",whichMatrix
  end if

  !predictedNNZForEachRowOfTotalMatrix = 4*(3*Nx + 5*3 + 5*3 + 5 + Nx)
  !predictedNNZForEachRowOfTotalMatrix = 4*(3*Nx + 5*3 + 5*3 + 5 + Nx + Ntheta*Nzeta)
  tempInt1 = Nspecies*Nx + 5*3 + 5*3 + 5 + 3*Nx + 2 + Nx*Ntheta*Nzeta
  if (tempInt1 > matrixSize) then
     tempInt1 = matrixSize
  end if
  predictedNNZForEachRowOfTotalMatrix = tempInt1

  ! In principle, we do not need to preallocate as many nonzeros in the preconditioner.
  ! But for simplicity, just use the same estimate for now.
  predictedNNZForEachRowOfPreconditioner = predictedNNZForEachRowOfTotalMatrix
  
  allocate(predictedNNZsForEachRow(matrixSize))
  allocate(predictedNNZsForEachRowDiagonal(matrixSize))
  tempInt1 = 3*Nx + (Nspecies-1)*Nx + (5*3-1) + (5*3-1)
  if (tempInt1 > matrixSize) then
     tempInt1 = matrixSize
  end if
  predictedNNZsForEachRow = tempInt1
  
  select case (constraintScheme)
  case (0)
  case (1)
     ! The rows for the constraints have more nonzeros:
     predictedNNZsForEachRow((Nspecies*Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta*Nx
  case (2)
     ! The rows for the constraints have more nonzeros:
     predictedNNZsForEachRow((Nspecies*Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta
  case default
  end select
  predictedNNZsForEachRowDiagonal = predictedNNZsForEachRow
  
  
  ! Allocate the main global matrix:
  select case (PETSCPreallocationStrategy)
  case (0)
     ! 0 = Old method with high estimated number-of-nonzeros.
     ! This method is simpler and works consistently but uses unnecessary memory.
     if (whichMatrix==0) then
        call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
             predictedNNZForEachRowOfPreconditioner, PETSC_NULL_INTEGER, &
             predictedNNZForEachRowOfPreconditioner, PETSC_NULL_INTEGER, &
             matrix, ierr)
     else
        call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
             predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER, &
             predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER, &
             matrix, ierr)
     end if
     
  case (1)
     ! 1 = New method with lower, more precise estimated number-of-nonzeros.
     ! This method is more complicated, but it should use much less memory.
     
     call MatCreate(MPIComm, matrix, ierr)
     !call MatSetType(matrix, MATMPIAIJ, ierr)
     call MatSetType(matrix, MATAIJ, ierr)
     
     numLocalRows = PETSC_DECIDE
     call PetscSplitOwnership(MPIComm, numLocalRows, matrixSize, ierr)
     
     call MatSetSizes(matrix, numLocalRows, numLocalRows, PETSC_DETERMINE, PETSC_DETERMINE, ierr)
     
     ! We first pre-allocate assuming number-of-nonzeros = 0, because due to a quirk in PETSc,
     ! MatGetOwnershipRange only works after MatXXXSetPreallocation is called:
     if (numProcs == 1) then
        call MatSeqAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER, ierr)
     else
        call MatMPIAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, ierr)
     end if
     
     call MatGetOwnershipRange(matrix, firstRowThisProcOwns, lastRowThisProcOwns, ierr)
     !print *,"I am proc ",myRank," and I own rows ",firstRowThisProcOwns," to ",lastRowThisProcOwns-1
     
     ! To avoid a PETSc error message, the predicted nnz for each row of the diagonal blocks must be no greater than the # of columns this proc owns:
     ! But we must not lower the predicted nnz for the off-diagonal blocks, because then the total predicted nnz for the row
     ! would be too low.
     tempInt1 = lastRowThisProcOwns - firstRowThisProcOwns
     do i=firstRowThisProcOwns+1,lastRowThisProcOwns
        if (predictedNNZsForEachRowDiagonal(i) > tempInt1) then
           predictedNNZsForEachRowDiagonal(i) = tempInt1
        end if
     end do
     
     ! Now, set the real estimated number-of-nonzeros:
     if (numProcs == 1) then
        call MatSeqAIJSetPreallocation(matrix, 0, predictedNNZsForEachRow(firstRowThisProcOwns+1:lastRowThisProcOwns), ierr)
     else
        call MatMPIAIJSetPreallocation(matrix, &
             0, predictedNNZsForEachRowDiagonal(firstRowThisProcOwns+1:lastRowThisProcOwns), &
             0, predictedNNZsForEachRow(firstRowThisProcOwns+1:lastRowThisProcOwns), ierr)
     end if
     
  case default
     print *,"Error! Invalid setting for PETSCPreallocationStrategy."
     stop
  end select
  
  ! If any mallocs are required during matrix assembly, do not generate an error:
  !call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
  
end subroutine preallocateMatrix
