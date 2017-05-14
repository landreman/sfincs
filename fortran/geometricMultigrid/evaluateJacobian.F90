#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif

  ! Syntax for PETSc version 3.5 and later
  subroutine evaluateJacobian(mysnes, stateVec, jacobian, jacobianPC, userContext, ierr)

    use petscsnes
    use globalVariables, only: masterProc, firstMatrixCreation, reusePreconditioner, &
         Mat_for_Jacobian, Mat_for_preconditioner

    implicit none

    SNES :: mysnes
    Vec :: stateVec
    Mat :: jacobian, jacobianPC
    PetscErrorCode :: ierr
    integer :: userContext(*)

    if (masterProc) print *,"evaluateJacobian called."

    ! When PETSc assembles a matrix, it reduces the structure of nonzeros to the actual number of nonzeros.
    ! If we try to re-assemble the matrix with additional nonzero entries without first re-allocating space for the nonzeros,
    ! we get the error about 'new nonzero caused a malloc'. Therefore, here we destroy the matrices and reallocate them.

    ! If reusePreconditioner = true, then we only need to assemble the preconditioner in the first iteration.
    ! If reusePreconditioner = false, then we need to assemble the preconditioner in every iteration.
    if (firstMatrixCreation .or. .not. reusePreconditioner) then
       !call populateMatrix(jacobianPC, 0, stateVec)
       call populateMatrix(Mat_for_preconditioner, 0, stateVec)
    end if
    if (firstMatrixCreation) then
       call preallocateMatrix(Mat_for_Jacobian, 1)
       call populateMatrix(Mat_for_Jacobian, 1, stateVec)
    end if

  ! Always build the high-order matrix at the finest level:
  call populateMatrix(levels(1)%high_order_matrix,1,stateVec,1)

  ! We need smoothing matrices on every level except the coarsest level, where we do a direct solve.
  do level = 1,N_levels-1
     call populateMatrix(levels(level)%smoothing_matrix,4,stateVec,level)
  end do

  select case (defect_option)
  case (1)
     ! In this case, we need low order matrices on every level
     do level = 1,N_levels
        ! If we re-use the factorization for the coarsest level, then only populate it on the first SNES iteration.
        if (level<N_levels .or. firstMatrixCreation .or. .not. reusePreconditioner) then
           call populateMatrix(levels(level)%low_order_matrix,0,stateVec,level)
        end if
     end do
  case (2)
     ! In this case, we need high order matrices on every level. We already did level 1, so start at level 2.
     do level = 2,N_levels
        ! If we re-use the factorization for the coarsest level, then only populate it on the first SNES iteration.
        if (level<N_levels .or. firstMatrixCreation .or. .not. reusePreconditioner) then
           call populateMatrix(levels(level)%high_order_matrix,1,stateVec,level)
        end if
     end do
  case default
     print *,"Error! Invalid defect_option:",defect_option
     stop
  end select


    firstMatrixCreation = .false.

  end subroutine evaluateJacobian


  subroutine apply_Jacobian(matrix, inputVec, outputVec, ierr)

    use globalVariables, only: masterProc, Mat_for_Jacobian, collisionOperator, preconditioner_field_term_xi_option

    implicit none

    Mat :: matrix
    Vec :: inputVec, outputVec
    PetscErrorCode :: ierr

    if (masterProc) print *,"Applying Jacobian matrix to a vector."

    !if (collisionOperator==0 .and. preconditioner_field_term_xi_option>0) then
    call VecZeroEntries(outputVec,ierr)
    call apply_dense_terms(inputVec, outputVec, 1)
    call MatMultAdd(Mat_for_Jacobian, inputVec, outputVec, outputVec, ierr)

  end subroutine apply_Jacobian
