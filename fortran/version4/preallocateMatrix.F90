#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

subroutine preallocateMatrix(matrix, whichMatrix)

  use petscmat
  use globalVariables, only: Nx, Nxi, NFourier2, Nspecies, matrixSize, includePhi1, &
       constraintScheme, PETSCPreallocationStrategy, MPIComm, numProcs, masterProc, & 
       includePhi1InKineticEquation, quasineutralityOption
  use indices

  implicit none

  integer, intent(in) :: whichMatrix
  Mat :: matrix
  integer :: predictedNNZForEachRowOfPreconditioner, predictedNNZForEachRowOfTotalMatrix
  integer, dimension(:), allocatable :: predictedNNZsForEachRow, predictedNNZsForEachRowDiagonal
  PetscErrorCode :: ierr
  integer :: tempInt1, i, imn, ispecies, ix, index
  integer :: firstRowThisProcOwns, lastRowThisProcOwns, numLocalRows

  if (masterProc) then
     print *,"Beginning preallocation for whichMatrix = ",whichMatrix
  end if

  tempInt1 = Nspecies*Nx & ! Collision operator.
       + NFourier2*5 & ! Streaming, ExB, and xiDot terms are pentadiagonal in L and are dense in Fourier.
       + NFourier2*3*Nx & ! xDot term is dense in Fourier and x, and occupies (L-2,L,L+2). Subtract 1 to avoid double-counting
       + 2 ! Sources
  if (tempInt1 > matrixSize) then
     tempInt1 = matrixSize
  end if
  predictedNNZForEachRowOfTotalMatrix = tempInt1

  ! In principle, we do not need to preallocate as many nonzeros in the preconditioner.
  ! But for simplicity, just use the same estimate for now.
  predictedNNZForEachRowOfPreconditioner = predictedNNZForEachRowOfTotalMatrix
  
  allocate(predictedNNZsForEachRow(matrixSize))
  allocate(predictedNNZsForEachRowDiagonal(matrixSize))
  ! Set tempInt1 to the expected number of nonzeros in a row of the kinetic equation block:


  tempInt1 = Nspecies*Nx & ! Collision operator
       + NFourier2*5 & ! Streaming, ExB, and xiDot terms are pentadiagonal in L and are dense in Fourier.
       + NFourier2*3*Nx & ! xDot term is dense in Fourier and x, and occupies (L-2,L,L+2).
       + 2 ! Sources


! THIS TERM HAS BEEN REMOVED BY AM 2016-03 !
!!$  if (includePhi1) then
!!$     tempInt1 = tempInt1 &
!!$!!       + 4 &               ! d Phi_1 / d theta term at L=1
!!$!!       + 4                 ! d Phi_1 / d zeta term at L=1, -1 because diagonal was accounted for in the previous line.
!!$  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The next block of code (for Phi1) has not been updated for the new Fourier discretization.
!!  if (includeRadialExBDrive) then !!Commented by AM 2016-03
  if (includePhi1InKineticEquation .and. includePhi1) then !!Added by AM 2016-03
     tempInt1 = tempInt1 &
       + 4 &               ! d Phi_1 / d theta term at L=0
       + 4                 ! d Phi_1 / d zeta term at L=0
  !!end if !!Commented by AM 2016-04
     tempInt1 = tempInt1 + 1 !!Added by AM 2016-04, for row = BLOCK_F, col = BLOCK_QN)
  !!if (nonlinear) then !!Commented by AM 2016-02
  !!if (includePhi1) then !!Added by AM 2016-02
     tempInt1 = tempInt1 + 2*Nx -2 ! Nonlinear term is dense in x with ell = L +/- 1, which we have not yet counted. Subtract 2 for the diagonal we already counted.
  end if
  ! Note: we do not need extra nonzeros for the d/dxi terms or for the diagonal, since these have already been accounted for in the ddx and ddtheta parts.
  if (tempInt1 > matrixSize) then
     tempInt1 = matrixSize
  end if
  predictedNNZsForEachRow = tempInt1

  ! Handle the rows for constraints
  select case (constraintScheme)
  case (0)
  case (1, 3, 4)
     ! The rows for the constraints have a different number of nonzeros:
     do ispecies=1,Nspecies
        index = getIndex(ispecies,1,1,1,BLOCK_DENSITY_CONSTRAINT)
        predictedNNZsForEachRow(index+1) = NFourier2*Nx + 1 ! +1 for diagonal
        index = getIndex(ispecies,1,1,1,BLOCK_PRESSURE_CONSTRAINT)
        predictedNNZsForEachRow(index+1) = NFourier2*Nx + 1 ! +1 for diagonal
     end do
  case (2)
     do ispecies=1,Nspecies
        do ix = 1, Nx
           index = getIndex(ispecies,ix,1,1,BLOCK_F_CONSTRAINT)
           predictedNNZsForEachRow(index+1) = NFourier2 + 1 ! +1 for diagonal
        end do
     end do
  case default
  end select
  
  ! MJL 20160720: I think the code block below is appropriately updated for Fourier, but should be checked.
  if (includePhi1) then
     ! Set rows for the quasineutrality condition:
     do imn = 1,NFourier2
        index = getIndex(1,1,1,imn,BLOCK_QN)

        !!Added by AM 2016-03!!
        if (quasineutralityOption == 1) then
!!!!!!!!!!!!!!!!!!!!!!!
           ! Add 1 because we are indexing a Fortran array instead of a PETSc matrix:
           predictedNNZsForEachRow(index+1) = Nx*Nspecies &  ! Integrals over f to get the density
                + 1 &          ! lambda
                + 1            ! Diagonal entry
              
           !!Added by AM 2016-03!!
        else
           ! Add 1 because we are indexing a Fortran array instead of a PETSc matrix:
           predictedNNZsForEachRow(index+1) = Nx*1 &  ! Integrals over f to get the density (only 1 kinetic species for EUTERPE equations)
                + 1 &          ! lambda
                + 1            ! Diagonal entry
        end if
!!!!!!!!!!!!!!!!!!!!!!!
     end do
     ! Set row for lambda constraint:
     index = getIndex(1,1,1,1,BLOCK_PHI1_CONSTRAINT)
     predictedNNZsForEachRow(index+1) = NFourier2 + 1 ! <Phi_1>, plus 1 for diagonal.
  end if
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
  
!!$  if (whichMatrix==0 .and. masterProc) then
!!$     print *,"Here comes predictedNNZsForEachRow:"
!!$     print *,predictedNNZsForEachRow
!!$     print *,"Here comes predictedNNZsForEachRowDiagonal:"
!!$     print *,predictedNNZsForEachRowDiagonal
!!$  end if

  ! If any mallocs are required during matrix assembly, do not generate an error:
!!$  call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
  
  !if (masterProc) then
  !   print *,"Done with preallocation for whichMatrix = ",whichMatrix
  !end if

end subroutine preallocateMatrix
