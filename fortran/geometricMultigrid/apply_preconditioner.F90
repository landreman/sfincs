#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

  ! This subroutine could probably be sped up by re-using Vecs instead of creating/destroying them over and over again?

  subroutine apply_preconditioner(outer_preconditioner, input_Vec, output_Vec, ierr)

    use petscksp
    use globalVariables, only: levels, masterProc, numProcs, myRank, inner_KSP, matrixSize, constraintScheme, Nx, Nspecies, one, &
         preconditioning_option, MPIComm
    use indices

    implicit none

    PetscErrorCode :: ierr
    PC :: outer_preconditioner
    Vec :: input_Vec, output_Vec
    KSPConvergedReason :: reason
    integer, dimension(:), allocatable :: IS_array
    integer :: IS_array_index, ix, ispecies, j
    ! For these next variables, the 'save' attribute means these variables can be initialized in the first pass
    ! through this subroutine, and used again on subsequent calls to this subroutine.
    KSP, save :: constraints_times_sources_KSP
    Mat, save :: sources_Mat, constraints_Mat, constraints_times_sources_Mat
    IS, save :: IS_all, IS_source_constraint
    logical, save :: first_call = .true.
    Vec :: temp_Vec_1, temp_Vec_2, ainv_times_stuff_Vec, c_times_r_Vec, y_Vec, s_Vec
    PC :: constraints_times_sources_PC
!    logical :: verbose = .true.
    logical :: verbose = .false.
    integer :: first_row_this_proc_owns, last_row_this_proc_owns
    PetscLogEvent :: event

    call PetscLogEventRegister("apply_preconditi", 0, event, ierr)
    call PetscLogEventBegin(event,ierr)

    select case (preconditioning_option)
       case (1)
          if (masterProc) print *,"apply_preconditioner called. Only applying inner_KSP."
          !print *,"Here comes KSPView on inner_KSP:"
          !call KSPView(inner_KSP, PETSC_VIEWER_STDOUT_WORLD,ierr)
          if (verbose) then
             print *,"000 Here comes input_Vec for the inner_KSP solve:"
             call VecView(input_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
             print *,"0.3 0.3 0.3 Here comes output_Vec before solve"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
             print *,"0.7 0.7 0.7 about to call KSPSolve."
          end if
          call KSPSolve(inner_KSP, input_Vec, output_Vec, ierr)
          if (verbose) then
             print *,"0.9 0.9 0.9 Here comes output_Vec after solve"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
             print *,"111"
          end if
          call KSPGetConvergedReason(inner_KSP, reason, ierr)
          print *,"222"
          if (reason <= 0) then
             print *,"WARNING: inner KSP failed with reason",reason
          end if

       case (2)
          if (masterProc) print *,"apply_preconditioner called. Handling sources/constraints."
          
          ! We will solve the block system
          !
          ! | a  b | |x|   |r|
          ! |      | | | = | |
          ! | c  0 | |y|   |s|
          !
          ! using the properties c a = 0 and a b = 0. The solution turns out to be
          ! x = (I - b (c b)^{-1} c) a^{-1} (I - b (c b)^{-1} c) r + b (c b)^{-1} s,
          ! y = (c b)^{-1} c r.

          ! First, if we haven't already done so, extract the sub-matrices b and c.
          if (first_call) then
             if (masterProc) print *,"Extracting sources and constraints for the shell preconditioner."

             ! Create an index set 'IS_all' which represents all indices of the big matrix and vectors, which each processor
             ! owning the indices it usually owns.
             allocate(IS_array(matrixSize))
             call VecGetOwnershipRange(input_Vec, first_row_this_proc_owns, last_row_this_proc_owns, ierr)
             IS_array(1:last_row_this_proc_owns-first_row_this_proc_owns) = [( j, j = first_row_this_proc_owns, last_row_this_proc_owns-1 )]
             call ISCreateGeneral(PETSC_COMM_WORLD,last_row_this_proc_owns-first_row_this_proc_owns,IS_array,PETSC_COPY_VALUES,IS_all,ierr)

             ! The next 2 lines only work in serial.
             !IS_array = [( j, j=0,matrixSize-1 )]
             !call ISCreateGeneral(PETSC_COMM_WORLD,matrixSize,IS_array,PETSC_COPY_VALUES,IS_all,ierr)
          
             ! Create an index set 'IS_source_constraint' that represents the indices for the sources and constraints.
             ! In this index set, the master processor owns everything, unlike the global matrix and vectors in which
             ! one or more processors at the end of the communicator own the sources/constraints.
             select case (constraintScheme)
             case (0)
                ! Nothing to do
             case (1,3,4)
                !if (masterProc) then
                if (myRank == numProcs-1) then
                   IS_array_index = 1
                   do ispecies = 1,Nspecies
                      IS_array(IS_array_index) = getIndex(1,ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
                      IS_array_index = IS_array_index + 1
                      IS_array(IS_array_index) = getIndex(1,ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
                      IS_array_index = IS_array_index + 1
                   end do
                   call ISCreateGeneral(PETSC_COMM_WORLD,2*Nspecies,IS_array,PETSC_COPY_VALUES,IS_source_constraint,ierr)
                else
                   IS_array=0
                   call ISCreateGeneral(PETSC_COMM_WORLD,0,IS_array,PETSC_COPY_VALUES,IS_source_constraint,ierr)
                end if
             case (2)
                !if (masterProc) then
                if (myRank == numProcs-1) then
                   IS_array_index = 1
                   do ispecies = 1,Nspecies
                      do ix = 1,Nx
                         IS_array(IS_array_index) = getIndex(1,ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT)
                         IS_array_index = IS_array_index + 1
                      end do
                   end do
                   call ISCreateGeneral(PETSC_COMM_WORLD,Nx*Nspecies,IS_array,PETSC_COPY_VALUES,IS_source_constraint,ierr)
                else
                   IS_array=0
                   call ISCreateGeneral(PETSC_COMM_WORLD,0,IS_array,PETSC_COPY_VALUES,IS_source_constraint,ierr)
                end if
             case default
                if (masterProc) print *,"Invalid constraintScheme:",constraintScheme
                stop
             end select

             ! Now 'slice' the big matrix to get the smaller non-square matrices that represent the sources and constraints:
             ! All rows, last columns -> sources_Mat
             call MatGetSubMatrix(levels(1)%high_order_matrix, IS_all, IS_source_constraint, MAT_INITIAL_MATRIX, sources_Mat, ierr)
             ! Last rows, all columns -> constraints_Mat
             call MatGetSubMatrix(levels(1)%high_order_matrix, IS_source_constraint, IS_all, MAT_INITIAL_MATRIX, constraints_Mat, ierr)

             ! The matrix given by the product of the constraints and sources will be used later. Compute it now.
             call MatMatMult(constraints_Mat, sources_Mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, constraints_times_sources_Mat, ierr)
             if (verbose) then
                if (masterProc) print *,"Here comes sources_Mat:"
                call MatView(sources_Mat,PETSC_VIEWER_STDOUT_WORLD,ierr)
                if (masterProc) print *,"Here comes constraints_Mat:"
                call MatView(constraints_Mat,PETSC_VIEWER_STDOUT_WORLD,ierr)
                if (masterProc) print *,"Here comes constraints_times_sources_Mat:"
                call MatView(constraints_times_sources_Mat,PETSC_VIEWER_STDOUT_WORLD,ierr)
             end if

             ! Set up the KSP for applying (c b)^{-1}:
             call KSPCreate(MPIComm, constraints_times_sources_KSP, ierr)
             call KSPSetType(constraints_times_sources_KSP, KSPPREONLY, ierr)
             call KSPGetPC(constraints_times_sources_KSP, constraints_times_sources_PC, ierr)
             call PCSetType(constraints_times_sources_PC, PCLU, ierr)
             call KSPSetOperators(constraints_times_sources_KSP, constraints_times_sources_Mat, constraints_times_sources_Mat, ierr)
          end if

          if (verbose) then
             if (masterProc) print *,"Here is input_Vec:"
             call VecView(input_Vec, PETSC_VIEWER_STDOUT_WORLD, ierr)
          end if

          ! Compute the part of the solution in the DKE rows (x) arising from s:
          ! x = b (c b)^{-1} s
          call VecGetSubVector(input_Vec, IS_source_constraint, s_Vec, ierr) ! Extract s from the global input_Vec.
          if (verbose) then
             if (masterProc) print *,"Here is s_Vec:"
             call VecView(s_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          
          if (masterProc) print *,"Here is temp_Vec_1:"
          call VecView(temp_Vec_1, PETSC_VIEWER_STDOUT_WORLD,ierr)

          call KSPSolve(constraints_times_sources_KSP, s_Vec, temp_Vec_1, ierr) ! Apply (c b)^{-1}
          call KSPGetConvergedReason(constraints_times_sources_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: constraints_times_sources_KSP failed (call 1) with reason",reason
          if (verbose) then
             if (masterProc) print *,"Here is (c b)^{-1} s:"
             call VecView(temp_Vec_1, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if
          if (first_call) then
             if (masterProc) print *,"Here comes constraints_times_sources_KSP:"
             call KSPView(constraints_times_sources_KSP, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if
          call VecRestoreSubVector(input_Vec, IS_source_constraint,s_Vec, ierr)
          call MatMult(sources_Mat, temp_Vec_1, output_Vec, ierr) ! Left-multiply by b, and store the result in output_Vec.
          call VecDestroy(temp_Vec_1, ierr)
          if (verbose) then
             if (masterProc) print *,"Here is output_Vec after the first term has been added:"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if

          ! Add the part of the solution in the DKE rows (x) arising from r:
          ! x += (I - b (c b)^{-1} c) a^{-1} (I - b (c b)^{-1} c) r
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, c_times_r_Vec, ierr) ! Create c_times_r_Vec of the appropriate size.
          call MatMult(constraints_Mat, input_Vec, c_times_r_Vec, ierr) ! Form c*r.
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          call KSPSolve(constraints_times_sources_KSP, c_times_r_Vec, temp_Vec_1, ierr) ! Apply (c b)^{-1}
          call KSPGetConvergedReason(constraints_times_sources_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: constraints_times_sources_KSP failed (call 2) with reason",reason
          if (verbose) then
             if (masterProc) print *,"Here is (c b)^{-1} c r:"
             call VecView(temp_Vec_1, PETSC_VIEWER_STDOUT_WORLD, ierr)
          end if
          call VecDuplicate(input_Vec, temp_Vec_2, ierr) ! Create temp_Vec_2 of the appropriate size.
          call MatMult(sources_Mat, temp_Vec_1, temp_Vec_2, ierr) ! Multiply by b to get b (c b)^{-1} c r. Result is now in temp_Vec_2.
          call VecDestroy(temp_Vec_1, ierr)
          call VecDuplicate(input_Vec, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          call VecWAXPY(temp_Vec_1, -one, temp_Vec_2, input_Vec, ierr) ! Now temp_Vec_1 holds (I - b (c b)^{-1} c) r.
          call VecDestroy(temp_Vec_2, ierr)
          call VecDuplicate(input_Vec, ainv_times_stuff_Vec, ierr) ! Create ainv_times_stuff_Vec of the appropriate size. 
          ! VVVVVVV  Next comes the big solve, involving multigrid. VVVVVVVVV
          call KSPSolve(inner_KSP, temp_Vec_1, ainv_times_stuff_Vec, ierr)
          ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          call KSPGetConvergedReason(inner_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: inner KSP failed with reason",reason
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, temp_Vec_2, ierr) ! Create temp_Vec_2 of the appropriate size.
          call MatMult(constraints_Mat, ainv_times_stuff_Vec, temp_Vec_2, ierr) ! Form c*a^{-1}*....
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          call KSPSolve(constraints_times_sources_KSP, temp_Vec_2, temp_Vec_1, ierr) ! Apply (c b)^{-1}
          call KSPGetConvergedReason(constraints_times_sources_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: constraints_times_sources_KSP failed (call 3) with reason",reason
          call VecDestroy(temp_Vec_2, ierr)
          call VecDuplicate(input_Vec, temp_Vec_2, ierr) ! Create temp_Vec_2 of the appropriate size.
          call MatMult(sources_Mat, temp_Vec_1, temp_Vec_2, ierr) ! Multiply by b to get b (c b)^{-1} c r. Result is now in temp_Vec_2.
          call VecDestroy(temp_Vec_1, ierr)
          call VecDuplicate(input_Vec, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          call VecWAXPY(temp_Vec_1, -one, temp_Vec_2, ainv_times_stuff_Vec, ierr) ! Now temp_Vec_1 holds (I - b (c b)^{-1} c) a^{-1} ....
          call VecDestroy(temp_Vec_2, ierr)
          call VecAXPY(output_Vec, one, temp_Vec_1, ierr) ! Add the result to the term in output_Vec we computed previously.
          call VecDestroy(temp_Vec_1, ierr)

          if (verbose) then
             if (masterProc) print *,"Here is output_Vec after the DKE rows have been finished:"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if

          ! Compute the part of the solution in the constraint rows (y):
          ! y = (c b)^{-1} c r
          ! I think this over-writes what was previously in the constraint rows of output_Vec, but I'm not certain.
          call VecGetSubVector(output_Vec, IS_source_constraint, y_Vec, ierr) ! Extract y from the global output_Vec.
          call KSPSolve(constraints_times_sources_KSP, c_times_r_Vec, y_Vec, ierr) ! Apply (c b)^{-1}
          call KSPGetConvergedReason(constraints_times_sources_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: constraints_times_sources_KSP failed (call 4) with reason",reason
          call VecRestoreSubVector(output_Vec, IS_source_constraint, y_Vec, ierr)
          call VecDestroy(c_times_r_Vec, ierr)

          if (verbose) then
             if (masterProc) print *,"Here is the final output_Vec:"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if

       case default
          print *,"Error! Invalid preconditioning_option:",preconditioning_option
          stop
    end select

    if (first_call) then
       first_call = .false.
       if (masterProc) print *,"Here comes inner_KSP:"
       call KSPView(inner_KSP, PETSC_VIEWER_STDOUT_WORLD,ierr)
    end if

    call PetscLogEventEnd(event,ierr)

  end subroutine apply_preconditioner
