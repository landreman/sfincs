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

  if (masterProc) then
      print*, "calling innerProduct"
  end if
  ! Compute inner product
  call innerProduct(forwardSolution,adjointRHSVec,innerProductResult,0)

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
        quantityToCompare = sum(FSABVelocityUsingFSADensityOverRootFSAB2(1:Nspecies)*Zs(1:Nspecies))
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


