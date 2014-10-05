! Next come some definitions required because the syntax for several PETSc objects
! has changed from version to version.
  
! For PETSc versions prior to 3.3, the MatCreateAIJ subroutine was called MatCreateMPIAIJ.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 3))
#define MatCreateAIJ MatCreateMPIAIJ
#endif
! Hereafter in this code, use MatCreateAIJ.

! For PETSc versions prior to 3.4, the PetscTime subroutine was called PetscGetTime.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 4))
#define PetscTime PetscGetTime
#endif
!Hereafter in this code, use PetscTime.

! For PETSc versions prior to 3.5, PETSC_DEFAULT_DOUBLE_PRECISION was used in place of PETSC_DEFAULT_REAL.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
#define PETSC_DEFAULT_REAL PETSC_DEFAULT_DOUBLE_PRECISION
#endif
!Hereafter in this code, use PETSC_DEFAULT_REAL.

module solver

  use globalVariables
  use petscsysdef

  implicit none

#include <finclude/petscsysdef.h>

  Mat :: matrix, preconditionerMatrix

  subroutine solveSystem

    implicit none

    PetscErrorCode :: ierr
    Vec :: solutionVec, residualVec
    SNES :: mysnes

    external evaluateJacobian, evaluateResidual

    call setMatrixPreallocation()

    call SNESCreate(PETSC_COMM_WORLD, mysnes, ierr)
    call SNESSetFunction(mysnes, residualVec, evaluateResidual, PETSC_NULL_OBJECT, ierr)
    call SNESSetJacobian(mysnes, matrix, preconditionerMatrix, evaluateJacobian, PETSC_NULL_OBJECT, ierr)

    call SNESGetKSP(mysnes, myksp, ierr)


  end subroutine solveSystem


  ! ------------------------------------------------------------------------

  subroutine setMatrixPreallocation()

    implicit none

    !predictedNNZForEachRowOfTotalMatrix = 4*(3*Nx + 5*3 + 5*3 + 5 + Nx)
    !predictedNNZForEachRowOfTotalMatrix = 4*(3*Nx + 5*3 + 5*3 + 5 + Nx + Ntheta*Nzeta)
    tempInt1 = Nspecies*Nx + 5*3 + 5*3 + 5 + 3*Nx + 2 + Nx*Ntheta*Nzeta
    if (tempInt1 > matrixSize) then
       tempInt1 = matrixSize
    end if
    predictedNNZForEachRowOfTotalMatrix = tempInt1

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

    if (useIterativeSolver) then
       whichMatrixMin = 0
    else
       whichMatrixMin = 1
    end if
    do whichMatrix = whichMatrixMin,1
       ! When whichMatrix = 0, build the preconditioner.
       ! When whichMatrix = 1, build the final matrix.

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
          ! This method is les thoroughly tested, but it should use much less memory.

          call MatCreate(MPIComm, matrix, ierr)
          !call MatSetType(matrix, MATMPIAIJ, ierr)
          call MatSetType(matrix, MATAIJ, ierr)

          numLocalRows = PETSC_DECIDE
          call PetscSplitOwnership(MPIComm, numLocalRows, matrixSize, ierr)

          call MatSetSizes(matrix, numLocalRows, numLocalRows, PETSC_DETERMINE, PETSC_DETERMINE, ierr)

          ! We first pre-allocate assuming number-of-nonzeros = 0, because due to a quirk in PETSc,
          ! MatGetOwnershipRange only works after MatXXXSetPreallocation is called:
          if (numProcsInSubComm == 1) then
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
          if (numProcsInSubComm == 1) then
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

       if (whichMatrix==0) then
          preconditionerMatrix = matrix
       end if

    end do

  end subroutine setMatrixPreallocation

  ! ------------------------------------------------------------------------

  subroutine chooseParallelDirectSolver()

    implicit none

    logical :: isAParallelDirectSolverInstalled

    isAParallelDirectSolverInstalled = .false.

    if ((whichParallelSolverToFactorPreconditioner<1) .or. (whichParallelSolverToFactorPreconditioner>2)) then
       print *,"Error! Invalid setting for whichParallelSolverToFactorPreconditioner"
       stop
    end if

#ifdef PETSC_HAVE_MUMPS
    isAParallelDirectSolverInstalled = .true.
    if (masterProc) then
       print *,"mumps detected"
    end if
#else
    whichParallelSolverToFactorPreconditioner = 2
    if (masterProc) then
       print *,"mumps not detected"
    end if
#endif

#ifdef PETSC_HAVE_SUPERLU_DIST
    isAParallelDirectSolverInstalled = .true.
    if (masterProc) then
       print *,"superlu_dist detected"
    end if
#else
    if (masterProc) then
       print *,"superlu_dist not detected"
    end if
    if (whichParallelSolverToFactorPreconditioner==2) then
       whichParallelSolverToFactorPreconditioner = 1
    end if
#endif

    if ((.not. isAParallelDirectSolverInstalled) .and. (numProcs > 1)) then
       if (masterProc) then
          print *,"Error! To run with more than 1 processors, you must have either"
          print *,"mumps no superlu_dist installed."
       end if
       stop
    end if

  end subroutine chooseParallelDirectSolver



end module scan

