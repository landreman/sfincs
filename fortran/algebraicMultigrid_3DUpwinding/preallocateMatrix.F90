#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

subroutine preallocateMatrix(matrix, whichMatrix)

  use kinds
  use petscmat
  use globalVariables, only: Nx, Nxi, Ntheta, Nzeta, Nspecies, matrixSize, includePhi1, &
       constraintScheme, PETSCPreallocationStrategy, MPIComm, numProcs, masterProc, & 
       includePhi1InKineticEquation, quasineutralityOption, collisionOperator, &
       preconditioner_field_term_xi_option, first_derivative_stencil, first_derivative_stencil_preconditioner, &
       second_derivative_stencil, second_derivative_stencil_preconditioner, stencil_width, upwinding_option
  use indices

  implicit none

  integer, intent(in) :: whichMatrix
  Mat :: matrix
  integer :: predictedNNZForEachRowOfPreconditioner, predictedNNZForEachRowOfTotalMatrix
  integer, dimension(:), allocatable :: predictedNNZsForEachRow, predictedNNZsForEachRowDiagonal
  PetscErrorCode :: ierr
  integer :: tempInt1, i, itheta, izeta, ispecies, ix, index, num_terms
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

!!$  tempInt1 = 3*Nx &        ! ddx terms on diagonal from collision operator, and ell=L +/- 2 in xDot. (Dense in x, and tridiagonal in L.)
!!$       + (Nspecies-1)*Nx & ! inter-species collisional coupling. (Dense in x and species, -Nx since we already counted the diagonal-in-species block)
!!$       + 2                 ! particle and heat sources
!!$
!!$  if (thetaDerivativeScheme==0) then
!!$     tempInt1 = tempInt1 + Ntheta*5-1  ! d/dtheta terms (dense in theta, pentadiagonal in L, -1 since we already counted the diagonal)
!!$  else
!!$     tempInt1 = tempInt1 + 5*5-1       ! d/dtheta terms (pentadiagonal in theta, pentadiagonal in L, -1 since we already counted the diagonal)
!!$  end if
!!$
!!$  if (zetaDerivativeScheme==0) then
!!$     tempInt1 = tempInt1 + Nzeta*5-1  ! d/dzeta terms (dense in theta, pentadiagonal in L, -1 since we already counted the diagonal)
!!$  else
!!$     tempInt1 = tempInt1 + 5*5-1      ! d/dzeta terms (pentadiagonal in theta, pentadiagonal in L, -1 since we already counted the diagonal)
!!$  end if
!!$
!!$  ! We don't need to separately count the d/dxi terms, since they just add to the diagonals of the d/dtheta and d/dzeta terms we already counted.

  tempInt1 = 2 & ! particle and heat sources
       + Nx ! xdot
  ! Below, the 1 is subtracted because we already counted the diagonal, above.
  !tempInt1 = tempInt1 + max(max_nnz_per_row(Ntheta,ddtheta_plus), max_nnz_per_row(Ntheta,ddtheta_plus_preconditioner)) - 1
  !tempInt1 = tempInt1 + max(max_nnz_per_row(Nzeta,ddzeta_plus), max_nnz_per_row(Nzeta,ddzeta_plus_preconditioner)) - 1
  !tempInt1 = tempInt1 + max(max_nnz_per_row(Nxi,ddxi_plus+100*pitch_angle_scattering_operator), &
  !     max_nnz_per_row(Nxi,ddxi_plus_preconditioner+100*pitch_angle_scattering_operator_preconditioner)) - 1
  num_terms = 3
  if (upwinding_option == 3) num_terms = 13
  do i = 1, num_terms ! Add nonzeros for each of the three derivative directions:
     tempInt1 = tempInt1 + max(nnz_in_stencil(first_derivative_stencil), nnz_in_stencil(first_derivative_stencil_preconditioner)) - 1
  end do
  ! Add nonzeros for the 2nd derivative in xi, which might not have been one of the 3 derivative directions above:
  tempInt1 = tempInt1 + max(nnz_in_stencil(second_derivative_stencil), nnz_in_stencil(second_derivative_stencil_preconditioner)) - 1
  ! I need to add the nonzeros for the Fokker-Planck operator.
  if (masterProc) then
     print *,"nnz for first_derivative_stencil:",nnz_in_stencil(first_derivative_stencil)
     print *,"nnz for first_derivative_stencil_preconditioner:",nnz_in_stencil(first_derivative_stencil_preconditioner)
     print *,"nnz for second_derivative_stencil:",nnz_in_stencil(second_derivative_stencil)
     print *,"nnz for second_derivative_stencil_preconditioner:",nnz_in_stencil(second_derivative_stencil_preconditioner)
  end if
  if (collisionOperator==0) then
     ! Eventually, add a test so these terms are only added if preconditioner_x=0.
     select case (preconditioner_field_term_xi_option)
     case (0)
        tempInt1 = tempInt1 + Nspecies*Nx*Nxi - 1
     case (1)
        tempInt1 = tempInt1 + Nspecies*Nx - 1
     case (2)
        tempInt1 = tempInt1 + Nspecies*Nx*3 - 1
     case default
        print *,"Error! Invalid preconditioner_field_term_xi_option:",preconditioner_field_term_xi_option
     end select
  end if

  if (includePhi1InKineticEquation .and. includePhi1) then !!Added by AM 2016-03
     tempInt1 = tempInt1 &
       + 4 &               ! d Phi_1 / d theta term at L=0
       + 4                 ! d Phi_1 / d zeta term at L=0
     tempInt1 = tempInt1 + 1 !!Added by AM 2016-04, for row = BLOCK_F, col = BLOCK_QN)
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
     do ispecies=1,Nspecies
        index = getIndex(ispecies,1,1,1,1,BLOCK_DENSITY_CONSTRAINT)
        predictedNNZsForEachRow(index+1) = Ntheta*Nzeta*Nxi*Nx + 1 ! +1 for diagonal
        index = getIndex(ispecies,1,1,1,1,BLOCK_PRESSURE_CONSTRAINT)
        predictedNNZsForEachRow(index+1) = Ntheta*Nzeta*Nxi*Nx + 1 ! +1 for diagonal
     end do
  case (2)
     ! The rows for the constraints have more nonzeros:
     do ispecies=1,Nspecies
        do ix = 1, Nx
           index = getIndex(ispecies,ix,1,1,1,BLOCK_F_CONSTRAINT)
           predictedNNZsForEachRow(index+1) = Ntheta*Nzeta*Nxi + 1 ! +1 for diagonal
        end do
     end do
  case default
  end select
  
  if (includePhi1) then
     ! Set rows for the quasineutrality condition:
     do itheta=1,Ntheta
        do izeta=1,Nzeta
           index = getIndex(1,1,1,itheta,izeta,BLOCK_QN)

           !!Added by AM 2016-03!!
           if (quasineutralityOption == 1) then
           !!!!!!!!!!!!!!!!!!!!!!!
              ! Add 1 because we are indexing a Fortran array instead of a PETSc matrix:
              predictedNNZsForEachRow(index+1) = Nx*Nxi*Nspecies &  ! Integrals over f to get the density
                   + 1 &          ! lambda
                   + 1            ! Diagonal entry
              
           !!Added by AM 2016-03!!
           else
              ! Add 1 because we are indexing a Fortran array instead of a PETSc matrix:
              predictedNNZsForEachRow(index+1) = Nx*Nxi*1 &  ! Integrals over f to get the density (only 1 kinetic species for EUTERPE equations)
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

contains

  function nnz_in_stencil(stencil)
    ! This function intentionally always includes a nonzero for the middle element, even if it is 0!
    implicit none
    
    integer :: nnz_in_stencil
    real(prec), dimension(-stencil_width:stencil_width) :: stencil
    integer :: index, counter
    
    nnz_in_stencil = 0
    do index = -stencil_width, stencil_width
       if (abs(stencil(index))>1e-12 .or. (index==0)) nnz_in_stencil = nnz_in_stencil + 1
    end do
    
  end function nnz_in_stencil
end subroutine preallocateMatrix

! -----------------------------------------------------

