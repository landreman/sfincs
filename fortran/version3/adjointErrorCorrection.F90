#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#include <finclude/petscsnesdef.h>
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscsnesdef.h>
#endif

subroutine adjointErrorCorrection(forwardSolution, adjointSolution,whichAdjointRHS,whichSpecies)

  use globalVariables
  use petscsnes
  use petscvec
  use adjointDiagnostics
  use indices

  Vec :: forwardSolution, adjointSolution
  integer :: whichAdjointRHS, whichSpecies
  Mat :: matrix_fine
  Vec :: RHS_fine, RHS
  integer :: userContext(1)
  SNES :: mysnes
  Vec :: dummy, residual, adjoint_interp, forward_interp
  PetscScalar :: innerProductResult
  PetscErrorCode :: ierr
  PetscScalar :: norm
  PetscScalar :: fineMeshMoment, coarseMeshMoment
  Vec :: adjointRHSVec, adjointRHSVec_coarse

  ! Prolongation matrix is build in createGrids_fine()
  if (masterProc) then
    print *,"This is adjointErrorCorrection."
    print *,"whichSpecies = ", whichSpecies
    print *,"whichAdjointRHS = ", whichAdjointRHS
  end if

  ! Create forward_interp Vec
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, forward_interp, ierr)
  call VecSet(forward_interp,zero,ierr)
  ! Interpolate forward solution
  call MatMult(adjointECProlongationMatrix, forwardSolution, forward_interp, ierr)

  if (debugAdjointEC) then
    adjoint_interp = adjointSolutino
  else
    ! Interpolate adjoint solution
    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, adjoint_interp,ierr)
    call VecSet(adjoint_interp,zero,ierr)
    call MatMult(adjointECProlongationmatrix, adjointSolution, adjoint_interp, ierr)
  end if

  ! Construct fine matrix and RHS
  ! Create dummy and set to zero
  if (masterProc) then
    print *,"Building fine grid matrix & RHS."
  end if
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, dummy, ierr)
  call VecSet(dummy, zero, ierr)

  ! whichMatrix = 7 -> matrixSize_fine
  call preallocateMatrix(matrix_fine,7) ! matrix_fine is initialized to 0 here
  call populateMatrix(matrix_fine,7,dummy) ! dummy is ignored here

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, RHS_fine, ierr)
  call VecSet(RHS_fine,zero,ierr)
  userContext(1) = 7
  call evaluateResidual(mysnes,dummy,RHS_fine,userContext,ierr)

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, RHS, ierr)
  call VecSet(RHS,zero,ierr)
  userContext(1) = 0
  call evaluateResidual(mysnes,dummy,RHS,userContext,ierr)

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, residual, ierr)
  call VecSet(residual,zero,ierr)
  ! Evaluate residual
  call MatMultAdd(matrix_fine, forward_interp, RHS_fine, residual, ierr)

  ! Evaluate inner product
  call innerProduct(adjoint_interp,residual,innerProductResult,1)

  !> Evaluate moment on fine grid
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, adjointRHSVec, ierr)
  call VecSet(adjointRHSVec, zero, ierr)
  ! 1 indicates fine grid
  call populateAdjointRHS(adjointRHSVec, whichAdjointRHS, whichSpecies, 1)
  call innerProduct(adjointRHSVec,forward_interp,fineMeshMoment,1)

  ! Evaluate moment on coarse grid
  call VecCreateMPI(MPIComm, PETSC_DECIdE, matrixSize, adjointRHSVec_coarse,ierr)
  call VecSet(adjointRHSVec_coarse, zero, ierr)
  call populateAdjointRHS(adjointRHSVec_coarse, whichAdjointRHS, whichSpecies, 0)
  call innerProduct(adjointRHSVec_coarse,forwardSolution,coarseMeshMoment,0)

  if (masterProc) then
    if (whichSpecies == 0) then
      select case (whichAdjointRHS)
        case (1) ! Particle flux
          radialCurrent_corrected = fineMeshMoment-innerProductResult
        case (2) ! Heat Flux
          totalHeatFlux_corrected = fineMeshMoment-innerProductResult
        case (3) ! Bootstrap
          bootstrap_corrected = fineMeshMoment-innerProductResult
      end select
    else
      select case (whichAdjointRHS)
        case (1) ! Particle flux
          print *,"Coarse mesh particleFlux: ", coarseMeshMoment
          particleFlux_corrected(whichSpecies) = fineMeshMoment-innerProductResult
        case (2) ! Heat Flux
          print *,"Coarse mesh heatFlux: ", coarseMeshMoment
          heatFlux_corrected(whichSpecies) = fineMeshMoment-innerProductResult
        case (3) ! Parallel Flow
          print *,"Coarse mesh parallel flow: ", coarseMeshMoment
          parallelFlow_corrected(whichSpecies) = fineMeshMoment-innerProductResult
      end select
    end if
  end if

  call MatDestroy(matrix_fine,ierr)
  call VecDestroy(RHS_fine,ierr)

end subroutine adjointErrorCorrection
