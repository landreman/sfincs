#include "../PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

!> This subroutine is used to test the innerProduct subroutine. The inner product
!! of the forward solution with the corresponding adjoint RHS should yield the same
!! result as heatFlux_vm_rHat, particleFlux_vm_rHat, or radial current using
!! FSABVelocityUsingFSADensityOverRootFSAB2.
!! @param forwardSolution Solution to forward equation
!! @param whichAdjointRHS Indicates which integrated quantity is being differentiated. If = 1 (particle flux), = 2 (heat flux), = 3 (bootstrap current).
!! @param whichSpecies Indicates species used for inner product. If = 0, corresponds to a species-summed quantity. If nonzero, indicates number of species.
subroutine testingInnerProduct(forwardSolution,whichAdjointRHS,whichSpecies)

  use globalVariables
  use petscvec
  use adjointDiagnostics

  implicit none

  integer :: whichSpecies, whichAdjointRHS
  Vec :: forwardSolution
  integer :: ispecies
  PetscErrorCode :: ierr
  PetscScalar :: result
  PetscScalar :: quantityToCompare
  Vec :: adjointRHSVec
  VecScatter :: VecScatterContext
  Vec :: adjointRHSOnProc0, forwardSolutionOnProc0
  PetscScalar, pointer :: adjointRHSArray(:), forwardSolutionArray(:)
  PetscScalar :: innerProductResult, residual

  if (masterProc) then
    print "(a,i4,a,i4,a)","--------- Computing innerProduct test with species ", whichSpecies," and whichAdjointRHS ", whichAdjointRHS,"-----------------------------------"
  end if

  !> Allocate adjointRHSVec
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointRHSVec, ierr)
  call VecSet(adjointRHSVec, zero, ierr)

  !> Populate RHS vec
  if (masterProc) then
    print*, "populating adjointRHS"
  end if
  call populateAdjointRHS(adjointRHSVec, whichAdjointRHS, whichSpecies)

  !> Create a scattering context for adjointRHSVec
  call VecScatterCreateToZero(adjointRHSVec, VecScatterContext, adjointRHSOnProc0, ierr)
  !> Create forwardSolutionOnProc0 of same type
  call VecDuplicate(adjointRHSOnProc0, forwardSolutionOnProc0, ierr)
  !> Send forwardSolution to master proc
  call VecScatterBegin(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  !> Send adjointRHS to master proc
  call VecScatterBegin(VecScatterContext, adjointRHSVec, adjointRHSOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, adjointRHSVec, adjointRHSOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  if (masterProc) then
    ! Convert the PETSc vector into a normal Fortran array
    call VecGetArrayF90(adjointRHSOnProc0, adjointRHSArray, ierr)
    call VecGetARrayF90(forwardSolutionOnProc0, forwardSolutionArray, ierr)
    print*, "calling innerProduct"
    ! Compute inner product
    call innerProduct(forwardSolutionArray,adjointRHSArray,innerProductResult)
  end if

  ! Compare with appropriate quantity
  if (masterProc) then
    print*, "selecting quantityToCompare"
  end if
  quantityToCompare = zero
  select case (whichAdjointRHS)
    case(1) ! Particle Flux
      if (whichSpecies == 0) then ! Radial current
        do ispecies=1,Nspecies
          quantityToCompare = quantityToCompare + particleFlux_vm_rN(ispecies)*Zs(ispecies)
        end do
      else
        quantityToCompare = particleFlux_vm_rN(whichSpecies)
      end if
    case(2) ! Heat Flux
      if (whichSpecies == 0) then ! Species-summed heat flux
        do ispecies=1,Nspecies
          quantityToCompare = quantityToCompare + heatFlux_vm_rN(ispecies)
        end do
      else
        quantityToCompare = heatFlux_vm_rN(whichSpecies)
      end if
    case(3) ! Bootstrap
      if (whichSpecies == 0) then ! Bootstrap
        do ispecies=1,Nspecies
          quantityToCompare = quantityToCompare + FSABVelocityUsingFSADensityOverRootFSAB2(ispecies)*Zs(ispecies)
        end do
      else
         quantityToCompare = FSABVelocityUsingFSADensityOverRootFSAB2(whichSpecies)
      end if
  end select

  residual = abs(quantityToCompare - innerProductResult)
  if (masterProc) then
    print "(a,es14.7,a)","quantityToCompare: ", quantityToCompare," -----------------------------"
    print "(a,es14.7,a)","innerProductResult: ", innerProductResult," -----------------------------"
    print "(a,es14.7,a)","--------- Residual: ",residual," -----------------------------"
  end if

end subroutine testingInnerProduct


