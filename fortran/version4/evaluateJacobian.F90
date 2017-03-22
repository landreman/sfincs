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
    use globalVariables, only: masterProc, firstMatrixCreation, reusePreconditioner, &
         Mat_for_Jacobian, fieldsplit, Nspecies, Nx, null_space_option, &
         Mat_for_preconditioner, inner_preconditioner, inner_KSP

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
    KSP :: myksp
    !PC :: preconditionerContext
    KSP, dimension(:), allocatable :: sub_ksps
    Mat :: sub_Amat, sub_Pmat
    MatNullSpace :: nullspace
    integer :: j, num_fieldsplits

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

    firstMatrixCreation = .false.

    if (fieldsplit .and. null_space_option>0) then
       if (null_space_option==1 .and. masterProc) print *,"Adding null space"
       if (null_space_option==2 .and. masterProc) print *,"Adding NEAR null space"
!!$       call SNESGetKSP(mysnes, myksp, ierr)
!!$       call KSPSetUp(myksp, ierr)
!!$       call KSPGetPC(myksp, preconditionerContext, ierr)
       call KSPSetUp(inner_KSP, ierr)
       allocate(sub_ksps(Nspecies*Nx+1))
       !call PCFieldSplitGetSubKSP(preconditionerContext, num_fieldsplits, sub_ksps, ierr)
       call PCFieldSplitGetSubKSP(inner_preconditioner, num_fieldsplits, sub_ksps, ierr)
       do j = 1,Nspecies*Nx
          call KSPGetOperators(sub_ksps(j), sub_Amat, sub_Pmat, ierr)
          print *,"Does sub_Amat==sub_Pmat?",sub_Amat==sub_Pmat
          call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
          select case (null_space_option)
          case (1)
             call MatSetNullSpace(sub_Pmat,nullspace,ierr)
          case (2)
             call MatSetNearNullSpace(sub_Pmat,nullspace,ierr)
          case default
             print *,"Invalid null_space_option:",null_space_option
             stop
          end select
          call MatNullSpaceDestroy(nullspace,ierr)
          !print *,"Here comes the Pmat for fieldsplit",j-1
          !call MatView(sub_Pmat,PETSC_VIEWER_STDOUT_WORLD,ierr)
       end do
       ! The next lines are temporary:
       j = Nspecies*Nx+1
       call KSPGetOperators(sub_ksps(j), sub_Amat, sub_Pmat, ierr)
       print *,"Does sub_Amat==sub_Pmat?",sub_Amat==sub_Pmat
       print *,"Here comes the Pmat for fieldsplit",j-1
       call MatView(sub_Pmat,PETSC_VIEWER_STDOUT_WORLD,ierr)
       deallocate(sub_ksps)
    else
       if (masterProc) print *,"NOT adding null space"
    end if

!!$    call SNESGetKSP(mysnes, myksp, ierr)
!!$    call KSPGetOperators(myksp, sub_Amat, sub_Pmat, ierr)
!!$    print *,"Here comes main Pmat:"
!!$    call MatView(sub_Pmat,PETSC_VIEWER_STDOUT_WORLD,ierr)

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
