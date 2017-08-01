#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

subroutine populateAdjointRHS(rhs, whichAdjointRHS, whichSpecies)

  use globalVariables
  use indices
  use petscvec

  implicit none

  Vec :: rhs
  integer :: whichAdjointRHS, whichSpecies

  PetscErrorCode :: ierr
  integer :: ix, L, itheta, izeta, ispecies, index, ixMin
  PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, xPartOfRHS, factor, nHat

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, rhs, ierr)
  call VecSet(rhs, zero, ierr)

  if (pointAtX0) then
     ixMin = 2
  else
     ixMin = 1
  end if

  THat = THats(whichSpecies)
  mHat = mHats(whichSpecies)
  nHat = nHats(whichSpecies)
  sqrtTHat = sqrt(THat)
  sqrtMHat = sqrt(mHat)

  select case (whichAdjointRHS)
    case (1) ! particle flux

      do ix=ixMin,Nx
        xPartOfRHS = exp(-x2(ix))*nHat*mHat*sqrtMHat*Delta*x2(ix)/ &
          (pi*sqrtpi*THat*sqrtTHat*Zs(whichSpecies))
        do itheta = ithetaMin,ithetaMax
           do izeta = izetaMin,izetaMax
              factor = xPartOfRHS/(BHat(itheta,izeta)**3)*(BHat_sub_theta(itheta,izeta) &
                *dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)* &
                dBHatdtheta(itheta,izeta))

              L = 0
              index = getIndex(whichSpecies, ix, L+1, itheta, izeta, BLOCK_F)
              call VecSetValue(rhs, index, (4/three)*factor, INSERT_VALUES, ierr)

              L = 2
              index = getIndex(whichSpecies, ix, L+1, itheta, izeta, BLOCK_F)
              call VecSetValue(rhs, index, (two/three)*factor, INSERT_VALUES, ierr)
           enddo
        enddo
      enddo

    case (2) ! heat flux

      do ix=ixMin,Nx
        xPartOfRHS = exp(-x2(ix))*nHat*mHat*sqrtMHat*Delta*x2(ix)*x2(ix)/ &
          (2*pi*sqrtpi*sqrtTHat*Zs(whichSpecies))
        do itheta = ithetaMin,ithetaMax
           do izeta = izetaMin,izetaMax
              factor = xPartOfRHS/(BHat(itheta,izeta)**3)*(BHat_sub_theta(itheta,izeta) &
                *dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)* &
                dBHatdtheta(itheta,izeta))

              L = 0
              index = getIndex(whichSpecies, ix, L+1, itheta, izeta, BLOCK_F)
              call VecSetValue(rhs, index, (4/three)*factor, INSERT_VALUES, ierr)

              L = 2
              index = getIndex(whichSpecies, ix, L+1, itheta, izeta, BLOCK_F)
              call VecSetValue(rhs, index, (two/three)*factor, INSERT_VALUES, ierr)
           enddo
        enddo
      enddo

    case (3) ! radial current of species

      do ix=ixMin,Nx
        xPartOfRHS = exp(-x2(ix))*2*Zs(whichSpecies)*nHat*mHat*sqrtMHat*x(ix)/ &
          (pi*sqrtpi*THat*THat*sqrtTHat*sqrt(FSABHat2))
        do itheta = ithetaMin,ithetaMax
           do izeta = izetaMin,izetaMax
              factor = xPartOfRHS*BHat(itheta,izeta)

              L = 1
              index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
              call VecSetValue(rhs, index, factor, INSERT_VALUES, ierr)
           enddo
        enddo
      enddo

  end select

  ! Done inserting values.
  ! Finally, assemble the RHS vector:
  call VecAssemblyBegin(rhs, ierr)
  call VecAssemblyEnd(rhs, ierr)


end subroutine populateAdjointRHS
