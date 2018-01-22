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
subroutine testingdMatrixdLambda(whichMode, whichLambda, deltaLambda)

  use globalVariables
  use geometry
  use petscmat

  implicit none

  integer :: whichMode, whichLambda
  PetscScalar :: deltaLambda
  Vec :: dummyVec
  Mat :: initMatrix, newMatrix, dMatrixdLambda

  ! Populate matrix
  call preallocateMatrix(initMatrix, 1) ! the whichMatrix argument doesn't matter here
  call populateMatrix(initMatrix,1,dummyVec)

  ! Populate dMatrixdLambda
  call preallocateMatrix(dMatrixdLambda, 1)
  call populatedMatrixdLambda(dMatrixdLambda, whichLambda, whichMode)

  ! Recompute geometry
  call updateVMECGeometry(whichMode, whichLambda, deltaLambda)

  ! Populate new matrix
  call preallocateMatrix(newMatrix, 1) ! the whichMatrix argument doesn't matter here
  call populateMatrix(newMatrix,1,dummyVec)

end subroutine testingdMatrixdLambda
