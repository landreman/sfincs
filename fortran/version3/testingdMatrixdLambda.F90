#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

!> This subroutine is used to test the derivatives of the
!! matrix computed in the populatedMatrixdLambda subroutine
!! using finite difference derivatives
!! @param whichMode Which mode to differentiate with respect to (index of ms and and ns).
!! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
!! @param deltaLambda Finite difference step size.
subroutine testingdMatrixdLambda(forwardSolution,whichMode, whichLambda)

  use globalVariables
  use geometry
  use petscmat
  use indices

  implicit none

  Vec :: forwardSolution, forwardSolution2, forwardSolution1
  integer :: whichMode, whichLambda
  Mat :: initMatrix, newMatrix, dMatrixdLambda_analytic
  PetscErrorCode :: ierr
  Vec :: result, result_analytic
  PetscScalar, pointer :: resultArray(:), result_analyticArray(:)
  Vec :: resultOnProc0, result_analyticOnProc0
  VecScatter :: VecScatterContext
  Vec :: dummyVec
  PetscScalar, dimension(:), allocatable :: percentError
  integer :: i

  ! Populate matrix
  call preallocateMatrix(initMatrix, 1) ! the whichMatrix argument doesn't matter here
  call populateMatrix(initMatrix,1,dummyVec)

  ! Populate dMatrixdLambda_analytic
  call preallocateMatrix(dMatrixdLambda_analytic, 1)
  call populatedMatrixdLambda(dMatrixdLambda_analytic, whichLambda, whichMode)

  ! Recompute geometry
  call updateVMECGeometry(whichMode, whichLambda,.false.)

  ! Populate new matrix
  call preallocateMatrix(newMatrix, 1) ! the whichMatrix argument doesn't matter here
  call populateMatrix(newMatrix,1,dummyVec)

  ! newMatrix = (newMatrix - initMatrix)
  call MatAXPY(newMatrix, -one, initMatrix, DIFFERENT_NONZERO_PATTERN,ierr)

  ! dMatrixdLambda = (newMatrix - matrix)/deltaLambda
  call MatScale(newMatrix,one/deltaLambda,ierr)

  ! Create result_analytic and set to zero
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, result_analytic, ierr)
  call VecSet(result_analytic, zero, ierr)

  ! Copy forwardSolution - it will be cleared by MatMult
  call VecDuplicate(forwardSolution,forwardSolution2,ierr)
  call VecCopy(forwardSolution,forwardSolution2,ierr)
  call VecDuplicate(forwardSolution,forwardSolution1,ierr)
  call VecCopy(forwardSolution,forwardSolution1,ierr)

  ! (dMatrixdLambda_analytic)(forwardSolution)
  call MatMult(dMatrixdLambda_analytic,forwardSolution1,result_analytic)

  ! Create result and set to zero
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, result, ierr)
  call VecSet(result, zero, ierr)

  ! (dMatrixdLambda)(forwardSolution)
  call MatMult(newMatrix, forwardSolution2, result)

  ! Destroy matrices
  call MatDestroy(newMatrix,ierr)
  call MatDestroy(dMatrixdLambda_analytic,ierr)
  call MatDestroy(initMatrix,ierr)

  ! Scatter result to masterProc
  call VecScatterCreateToZero(result, VecScatterContext, resultOnProc0, ierr)
  call VecScatterBegin(VecScatterContext, result, resultOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, result, resultOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  if (masterProc) then
    ! Convert the PETSc vector into a normal Fortran array
    call VecGetArrayF90(resultOnProc0, resultArray, ierr)
  end if

  ! Scatter result_analytic to masterProc
  call VecScatterCreateToZero(result_analytic, VecScatterContext, result_analyticOnProc0, ierr)
  call VecScatterBegin(VecScatterContext, result_analytic, result_analyticOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, result_analytic, result_analyticOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  if (masterProc) then
    ! Convert the PETSc vector into a normal Fortran array
    call VecGetArrayF90(result_analyticOnProc0, result_analyticArray, ierr)
  end if

  ! Free results
  call VecDestroy(result_analytic,ierr)
  call VecDestroy(result,ierr)

  ! Compute percentError
  allocate(percentError(matrixSize))
  do i=1,matrixSize
    if (abs(resultArray(i)) < 1.d-12) then
      if (abs(result_analyticArray(i)) < 1.d-12) then
        percentError(i) = zero
      else
        percentError(i) = 1.d6
      end if
    else
      percentError(i) = 100.0*abs(result_analyticArray(i) - resultArray(i))/abs(resultArray(i))
    end if
    if (percentError(i) > 1.0) then
      print *,"percentError: ", percentError(i)
      print *,"result: ", resultArray(i)
      print *,"result_analytic: ", result_analyticArray(i)
    end if
  end do

  ! Compute metrics for percentError
  print *,"Max result: ", maxval(resultArray)
  print *,"Max result_analytic: ", maxval(result_analyticArray)
  print *,"Maximum percentError: ", maxval(percentError), "%"
  print *,"resultArray at max: ", resultArray(maxloc(percentError))
  print *,"result_analyticArray at max: ", result_analyticArray(maxloc(percentError))

  ! Recompute geometry
  call updateVMECGeometry(whichMode, whichLambda,.true.)

end subroutine testingdMatrixdLambda
