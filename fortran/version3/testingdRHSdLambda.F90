#include "../PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscsnesdef.h>
#endif

!> This subroutine is used to test the derivatives of the 
!! RHS computed in the populatedRHSdLambda subroutine
!! using finite difference derivatives
!! @param whichMode Which mode to differentiate with respect to (index of ms and and ns).
!! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
!! @param deltaLambda Finite difference step size.
subroutine testingdRHSdLambda(whichMode, whichLambda, deltaLambda)

  use globalVariables
  use geometry
  use petscsnes
  use petscvec

  implicit none

  integer :: whichMode, whichLambda
  PetscScalar :: deltaLambda
  SNES :: mysnes
  Vec :: stateVec, RHS, dRHSdLambda, dRHSdLambda_analytic, dRHSdLambdaOnProc0, dRHSdLambda_analyticOnProc0
  PetscErrorCode :: ierr
  integer :: userContext(1)
  PetscScalar, pointer :: dRHSdLambdaArray(:), dRHSdLambda_analyticArray(:)
  PetscScalar, dimension(:), allocatable :: percentError
  integer :: i
  VecScatter :: VecScatterContext

  ! Create stateVec and set to zero
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, stateVec, ierr)
  call VecSet(stateVec, zero, ierr)

  ! Create RHS
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, RHS, ierr)

  ! Pass stateVec=0 to evaluateResidual to get -RHS
  call evaluateResidual(mysnes, stateVec, RHS, userContext, ierr)

  ! Multiply the residual by (-1)
  call VecScale(RHS, -1d+0, ierr)

  ! Update geometry
  call updateVMECGeometry(whichMode, whichLambda, deltaLambda)

  ! Compute new RHS
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, dRHSdLambda, ierr)

  ! Pass stateVec=0 to evaluateResidual to get -RHS
  call evaluateResidual(mysnes, stateVec, dRHSdLambda, userContext, ierr)

  ! Multiply the residual by (-1) to get new RHS
  call VecScale(dRHSdLambda, -1d+0, ierr)

  ! Form finite diff derivative
  call VecAXPY(dRHSdLambda,-1d+0,RHS,ierr)
  call VecScale(dRHSdLambda,one/deltaLambda,ierr)

  ! Compute analytic derivatives
  call VecCreateMPI(MPIComm,PETSC_DECIDE, matrixSize, dRHSdLambda_analytic,ierr)
  call populatedRHSdLambda(dRHSdLambda_analytic, whichLambda, whichMode)

  ! Send to masterProc
  !> Create a scattering context for dRHSdLambda
  call VecScatterCreateToZero(dRHSdLambda, VecScatterContext, dRHSdLambdaOnProc0, ierr)
  !> Send dRHSdLambda to master proc
  call VecScatterBegin(VecScatterContext, dRHSdLambda, dRHSdLambdaOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, dRHSdLambda, dRHSdLambdaOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  if (masterProc) then
    ! Convert the PETSc vector into a normal Fortran array
    call VecGetARrayF90(dRHSdLambdaOnProc0, dRHSdLambdaArray, ierr)
  end if

  !> Create a scattering context for dRHSdLambda_analytic
  call VecScatterCreateToZero(dRHSdLambda_analytic, VecScatterContext, dRHSdLambda_analyticOnProc0, ierr)
  !> Send dRHSdLambda to master proc
  call VecScatterBegin(VecScatterContext, dRHSdLambda_analytic, dRHSdLambda_analyticOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, dRHSdLambda_analytic, dRHSdLambda_analyticOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  if (masterProc) then
    ! Convert the PETSc vector into a normal Fortran array
    call VecGetARrayF90(dRHSdLambda_analyticOnProc0, dRHSdLambda_analyticArray, ierr)
  end if

  ! Compute percentError
  allocate(percentError(matrixSize))
  do i=1,matrixSize
    if (abs(dRHSdLambdaArray(i)) < 1.d-16) then
      if (abs(dRHSdLambdaArray(i)) < 1.d-16) then
        percentError(i) = zero
      else
        percentError(i) = 1.d6
      end if
    else
      percentError(i) = 100.0*abs(dRHSdLambda_analyticArray(i) - dRHSdLambdaArray(i))/abs(dRHSdLambdaArray(i))
    end if
  end do

  ! Compute metrics for percentError
  print *,"Maximum percentError: ", maxval(percentError), "%"
  print *,"dRHSdLambdaArray at max: ", dRHSdLambdaArray(maxloc(percentError))
  print *,"dRHSdLambda_analyticArray at max: ", dRHSdLambda_analyticArray(maxloc(percentError))
!  print *,"Minimum percentError: ", minval(percentError), "%"
!  print *,"Maximum dRHSdLambdaArray: ", maxval(dRHSdLambdaArray)
!  print *,"Maximum dRHSdLambda_analyticArray: ", maxval(dRHSdLambda_analyticArray)

end subroutine testingdRHSdLambda
