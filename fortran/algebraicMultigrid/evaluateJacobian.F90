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
         Mat_for_Jacobian, fieldsplit, Nspecies, Nx, null_space_option, &
         Mat_for_preconditioner, inner_preconditioner, inner_KSP, &
         myRank, preconditioning_option, ISs

    implicit none

    SNES :: mysnes
    Vec :: stateVec
    Mat :: jacobian, jacobianPC
    PetscErrorCode :: ierr
    integer :: userContext(*)
    KSP :: myksp, Richardson_ksp
    KSP, dimension(:), allocatable :: sub_ksps
    Mat :: sub_Amat, sub_Pmat
    MatNullSpace :: nullspace
    integer :: j, num_fieldsplits, level, N_levels
    integer :: firstRowThisProcOwns, lastRowThisProcOwns, num_rows, num_columns
    PC :: gamg_pc
    Mat :: Amat, Pmat, restriction_matrix, restriction_matrix_transpose, interpolation_matrix, residual_matrix, last_matrix
    integer :: N_rows, N_cols
    integer :: index_of_min, index_of_max
    Vec :: diagonal_vec
    PetscReal :: diagonal_min, diagonal_max

    if (masterProc) then
       print *,"evaluateJacobian called."
    end if

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

    if (fieldsplit) then
       call KSPSetUp(inner_KSP, ierr) 
       allocate(sub_ksps(Nspecies*Nx+1))
       !call PCFieldSplitGetSubKSP(preconditionerContext, num_fieldsplits, sub_ksps, ierr)
       call PCFieldSplitGetSubKSP(inner_preconditioner, num_fieldsplits, sub_ksps, ierr)
       do j = 1,Nspecies*Nx+1
          call KSPGetOperators(sub_ksps(j), sub_Amat, sub_Pmat, ierr)
          print *,"Does sub_Amat==sub_Pmat?",sub_Amat==sub_Pmat
          call MatGetOwnershipRange(sub_Pmat, firstRowThisProcOwns, lastRowThisProcOwns, ierr)
          call MatGetSize(sub_Pmat, num_rows, num_columns, ierr)
          print "(a,i2,a,i3,a,i5,a,i7,a,i7)","Fieldsplit ",j,": Proc",myRank," owns indices",firstRowThisProcOwns," to",lastRowThisProcOwns-1," of",num_rows
          ! Print info about the max and min values along the diagonal, to verify the matrix scaling makes sense.
          call MatCreateVecs(sub_Pmat, diagonal_vec, PETSC_NULL_OBJECT, ierr)
          call MatGetDiagonal(sub_Pmat, diagonal_vec, ierr)
          call VecMin(diagonal_vec, index_of_min, diagonal_min, ierr)
          call VecMax(diagonal_vec, index_of_max, diagonal_max, ierr)
          if (masterProc) print "(a,i2,a,es11.3,a,i4,a,es11.3,a,i4)"," Fieldsplit ",j,": min=",diagonal_min," at index ",index_of_min,", max=",diagonal_max," at index ",index_of_max
          call VecDestroy(diagonal_vec, ierr)
       end do
    end if

    if (fieldsplit .and. null_space_option>0) then
       if (null_space_option==1 .and. masterProc) print *,"Adding null space"
       if (null_space_option==2 .and. masterProc) print *,"Adding NEAR null space"
!!$       call SNESGetKSP(mysnes, myksp, ierr)
!!$       call KSPSetUp(myksp, ierr)
!!$       call KSPGetPC(myksp, preconditionerContext, ierr)
       call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
       do j = 1,Nspecies*Nx
          select case (null_space_option)
          case (1)
             call MatSetNullSpace(sub_Pmat,nullspace,ierr)
          case (2)
             call MatSetNearNullSpace(sub_Pmat,nullspace,ierr)
          case default
             print *,"Invalid null_space_option:",null_space_option
             stop
          end select
          !print *,"Here comes the Pmat for fieldsplit",j-1
          !call MatView(sub_Pmat,PETSC_VIEWER_STDOUT_WORLD,ierr)
       end do
       call MatNullSpaceDestroy(nullspace,ierr)
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


    if (preconditioning_option==2 .or. preconditioning_option==5) then
       ! For GAMG, we need to change the sub-KSPs from Chebyshev (which doesn't seem to work for sfincs) to Richardson.
       ! We put this loop after the addition of the null space because it seems necessary to call KSPSetUp on the gamg block
       ! before extracting the Chebyshev KSPs inside.
       do j = 1,Nspecies*Nx
          call KSPSetUp(sub_ksps(j), ierr)
          call KSPGetPC(sub_ksps(j), gamg_pc, ierr)
          call PCMGGetLevels(gamg_pc, N_levels, ierr)
          if (masterProc) print *,"N_levels=",N_levels
          do level = 1,N_levels-1 ! Level 0 is the coarse level, and we don't want to change the KSP for that level.
             call PCMGGetSmoother(gamg_pc, level, Richardson_ksp, ierr)
             call KSPSetType(Richardson_ksp, KSPRICHARDSON, ierr)
          end do

          ! The following lines are my method for AMG evaluating the residual using Amat rather than Pmat.
!!$          !call KSPGetOperators(sub_ksps(j), Amat, Pmat, ierr)
!!$          !last_matrix = Pmat
!!$          call MatGetSubMatrix(Mat_for_Jacobian, ISs(j), ISs(j), MAT_INITIAL_MATRIX, last_matrix, ierr)
!!$          !call MatGetSubMatrix(Mat_for_preconditioner, ISs(j), ISs(j), MAT_INITIAL_MATRIX, last_matrix, ierr)
!!$          ! Level 0 is the coarse level
!!$          do level = N_levels-2,0,-1
!!$             if (masterProc) print *,"Setting residual matrix for level",level
!!$             ! The restriction matrix between level i and i-1 is indexed by i.
!!$             call PCMGGetRestriction(gamg_pc, level+1, restriction_matrix, ierr)
!!$             print *,"BBB"
!!$             !call MatTranspose(restriction_matrix,MAT_INITIAL_MATRIX,restriction_matrix,ierr) ! Could also use MatCreateTranspose?
!!$             !call MatCreateTranspose(restriction_matrix, restriction_matrix_transpose, ierr)
!!$             !call MatDuplicate(restriction_matrix, MAT_COPY_VALUES, restriction_matrix_transpose, ierr)
!!$             call MatTranspose(restriction_matrix, MAT_INITIAL_MATRIX, restriction_matrix_transpose, ierr)
!!$             call MatGetSize(restriction_matrix_transpose, N_rows, N_cols, ierr)
!!$             print *,"Restriction_matrix_transpose is ",N_rows,"x",N_cols
!!$             call PCMGGetInterpolation(gamg_pc, level+1, interpolation_matrix, ierr)
!!$             print *,"CCC"
!!$             call MatGetSize(last_matrix, N_rows, N_cols, ierr)
!!$             print *,"last_matrix is ",N_rows,"x",N_cols
!!$             call MatGetSize(interpolation_matrix, N_rows, N_cols, ierr)
!!$             print *,"Interpolation matrix is ",N_rows,"x",N_cols
!!$             !call MatMatMult(Amat,interpolation_matrix,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,temp_matrix,ierr) ! temp_matrix = last_matrix * interpolation_matrix
!!$             print *,"DDD"
!!$             !call MatMatMult(restriction_matrix,temp_matrix,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,Amat, ierr) ! last_matrix = restriction_matrix * temp_matrix
!!$             call MatMatMatMult(restriction_matrix_transpose, last_matrix, interpolation_matrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, residual_matrix, ierr)
!!$             call MatGetSize(residual_matrix, N_rows, N_cols, ierr)
!!$             print *,"Residual matrix is ",N_rows,"x",N_cols
!!$             call PCMGSetResidual(gamg_pc, level, PCMGResidualDefault, residual_matrix, ierr)
!!$             last_matrix = residual_matrix
!!$          end do
       end do


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
