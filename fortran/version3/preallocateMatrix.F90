subroutine preallocateMatrix(matrix, whichMatrix)

#include "PETScVersions.F90"

  use globalVariables, only: Nx, Nxi, Ntheta, Nzeta, Nspecies, matrixSize, includePhi1, &
       !!constraintScheme, PETSCPreallocationStrategy, MPIComm, numProcs, masterProc, nonlinear, & !!Commented by AM 2016-02
       constraintScheme, PETSCPreallocationStrategy, MPIComm, numProcs, masterProc, & !!Added by AM 2016-02
       !!thetaDerivativeScheme, zetaDerivativeScheme, includeRadialExBDrive !!Commented by AM 2016-03
       thetaDerivativeScheme, zetaDerivativeScheme, includePhi1InKineticEquation, quasineutralityOption, readExternalPhi1 !!Added by AM 2016-03 and 2018-12
  use indices

  implicit none

  integer, intent(in) :: whichMatrix
  Mat :: matrix
  integer :: predictedNNZForEachRowOfPreconditioner, predictedNNZForEachRowOfTotalMatrix
  integer, dimension(:), allocatable :: predictedNNZsForEachRow, predictedNNZsForEachRowDiagonal
  PetscErrorCode :: ierr
  integer :: tempInt1, i, itheta, izeta, ispecies, ix, index
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
  ! Set tempInt1 to the expected number of nonzeros in a row of the kinetic equation block:

  tempInt1 = 3*Nx &        ! ddx terms on diagonal from collision operator, and ell=L +/- 2 in xDot. (Dense in x, and tridiagonal in L.)
       + (Nspecies-1)*Nx & ! inter-species collisional coupling. (Dense in x and species, -Nx since we already counted the diagonal-in-species block)
       + 2                 ! particle and heat sources

  if (thetaDerivativeScheme==0) then
     tempInt1 = tempInt1 + Ntheta*5-1  ! d/dtheta terms (dense in theta, pentadiagonal in L, -1 since we already counted the diagonal)
  else
     tempInt1 = tempInt1 + 5*5-1       ! d/dtheta terms (pentadiagonal in theta, pentadiagonal in L, -1 since we already counted the diagonal)
  end if

  if (zetaDerivativeScheme==0) then
     tempInt1 = tempInt1 + Nzeta*5-1  ! d/dzeta terms (dense in theta, pentadiagonal in L, -1 since we already counted the diagonal)
  else
     tempInt1 = tempInt1 + 5*5-1      ! d/dzeta terms (pentadiagonal in theta, pentadiagonal in L, -1 since we already counted the diagonal)
  end if

  ! We don't need to separately count the d/dxi terms, since they just add to the diagonals of the d/dtheta and d/dzeta terms we already counted.

! THIS TERM HAS BEEN REMOVED BY AM 2016-03 !
!!$  if (includePhi1) then
!!$     tempInt1 = tempInt1 &
!!$!!       + 4 &               ! d Phi_1 / d theta term at L=1
!!$!!       + 4                 ! d Phi_1 / d zeta term at L=1, -1 because diagonal was accounted for in the previous line.
!!$  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!  if (includeRadialExBDrive) then !!Commented by AM 2016-03
  if (includePhi1InKineticEquation .and. includePhi1) then 
     if (.not. readExternalPhi1) then !!Added by AM 2018-12
        tempInt1 = tempInt1 &
             + 4 &               ! d Phi_1 / d theta term at L=0
             + 4                 ! d Phi_1 / d zeta term at L=0
        !!end if !!Commented by AM 2016-04
        tempInt1 = tempInt1 + 1 !!Added by AM 2016-04, for row = BLOCK_F, col = BLOCK_QN)
     end if !!Added by AM 2018-12
  !!if (nonlinear) then !!Commented by AM 2016-02
  !!if (includePhi1) then !!Added by AM 2016-02
     tempInt1 = tempInt1 + 2*Nx -2 ! Nonlinear term is dense in x with ell = L +/- 1, which we have not yet counted. Subtract 2 for the diagonal we already counted.
  end if
  ! Note: we do not need extra nonzeros for the d/dxi terms or for the diagonal, since these have already been accounted for in the ddx and ddtheta parts.
  if (tempInt1 > matrixSize) then
     tempInt1 = matrixSize
  end if
  predictedNNZsForEachRow = tempInt1
  
  select case (constraintScheme)
  case (0)
  case (1, 3, 4)
     ! The rows for the constraints have more nonzeros:
     !predictedNNZsForEachRow((Nspecies*Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta*Nx + 1
     !predictedNNZsForEachRow((Nspecies*Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta*Nx
     do ispecies=1,Nspecies
        index = getIndex(ispecies,1,1,1,1,BLOCK_DENSITY_CONSTRAINT)
        predictedNNZsForEachRow(index+1) = Ntheta*Nzeta*Nx + 1 ! +1 for diagonal
        index = getIndex(ispecies,1,1,1,1,BLOCK_PRESSURE_CONSTRAINT)
        predictedNNZsForEachRow(index+1) = Ntheta*Nzeta*Nx + 1 ! +1 for diagonal
     end do
  case (2)
     ! The rows for the constraints have more nonzeros:
     !predictedNNZsForEachRow((Nspecies*Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta + 1
     !predictedNNZsForEachRow((Nspecies*Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta
     do ispecies=1,Nspecies
        do ix = 1, Nx
           index = getIndex(ispecies,ix,1,1,1,BLOCK_F_CONSTRAINT)
           predictedNNZsForEachRow(index+1) = Ntheta*Nzeta + 1 ! +1 for diagonal
        end do
     end do
  case default
  end select
  
  !!if (includePhi1) then !!Commented by AM 2018-12
  if (includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
     ! Set rows for the quasineutrality condition:
     do itheta=1,Ntheta
        do izeta=1,Nzeta
           index = getIndex(1,1,1,itheta,izeta,BLOCK_QN)

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
     end do
     ! Set row for lambda constraint:
     index = getIndex(1,1,1,1,1,BLOCK_PHI1_CONSTRAINT)
     predictedNNZsForEachRow(index+1) = Ntheta*Nzeta + 1 ! <Phi_1>, plus 1 for diagonal.
  end if
  predictedNNZsForEachRowDiagonal = predictedNNZsForEachRow
  
  
  ! Allocate the main global matrix:
  select case (PETSCPreallocationStrategy)
  case (0)
     ! 0 = Old method with high estimated number-of-nonzeros.
     ! This method is simpler and works consistently but uses unnecessary memory.
     if (whichMatrix==0) then
        call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
             predictedNNZForEachRowOfPreconditioner, PETSC_NULL_INTEGER_ARRAY, &
             predictedNNZForEachRowOfPreconditioner, PETSC_NULL_INTEGER_ARRAY, &
             matrix, ierr)
     else
        call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
             predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER_ARRAY, &
             predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER_ARRAY, &
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
        call MatSeqAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER_ARRAY, ierr)
     else
        call MatMPIAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER_ARRAY, 0, PETSC_NULL_INTEGER_ARRAY, ierr)
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
