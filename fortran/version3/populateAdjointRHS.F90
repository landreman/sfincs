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

  ! Validate input
  if (whichAdjointRHS<1 .or. whichAdjointRHS>3 .or. whichSpecies<0 .or. whichSpecies>Nspecies) then
    if (masterProc) then
      print *,"Error! Incorrect input to populateAdjointRHS"
    end if
    stop
  end if
  if (whichAdjointRHS==3 .and. whichSpecies/=0) then
    if (masterProc) then
      print *,"Error! Incorrect input to populateAdjointRHS"
    end if
  end if

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, rhs, ierr)
  call VecSet(rhs, zero, ierr)

  if (pointAtX0) then
     ixMin = 2
  else
     ixMin = 1
  end if

  select case (whichAdjointRHS)
    case (1) ! particle flux

    do ispecies=1,Nspecies
      if (whichSpecies == ispecies .or. whichSpecies == 0) then
        THat = THats(ispecies)
        mHat = mHats(ispecies)
        nHat = nHats(ispecies)
        sqrtTHat = sqrt(THat)
        sqrtMHat = sqrt(mHat)
        do ix=ixMin,Nx
          xPartOfRHS = exp(-x2(ix))*nHat*mHat*sqrtMHat*Delta*x2(ix)/ &
            (pi*sqrtpi*THat*sqrtTHat*Zs(ispecies))
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                factor = xPartOfRHS/(BHat(itheta,izeta)**3)*(BHat_sub_theta(itheta,izeta) &
                  *dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)* &
                  dBHatdtheta(itheta,izeta))
                ! For species summed radial current, weighted by Zs
                if (whichSpecies == 0) then
                  factor = factor*Zs(ispecies)
                end if

                L = 0
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, (4/three)*factor, INSERT_VALUES, ierr)

                L = 2
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, (two/three)*factor, INSERT_VALUES, ierr)

             enddo
            enddo
          enddo
        end if
      enddo

    case (2) ! heat flux

    do ispecies=1,Nspecies
      if (whichSpecies == ispecies .or. whichSpecies == 0) then
        THat = THats(ispecies)
        mHat = mHats(ispecies)
        nHat = nHats(ispecies)
        sqrtTHat = sqrt(THat)
        sqrtMHat = sqrt(mHat)
        do ix=ixMin,Nx
          xPartOfRHS = exp(-x2(ix))*nHat*mHat*sqrtMHat*Delta*x2(ix)*x2(ix)/ &
            (2*pi*sqrtpi*sqrtTHat*Zs(ispecies))
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                factor = xPartOfRHS/(BHat(itheta,izeta)**3)*(BHat_sub_theta(itheta,izeta) &
                  *dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)* &
                  dBHatdtheta(itheta,izeta))
                ! factor is the samed for species-summed heat flux

                L = 0
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, (4/three)*factor, INSERT_VALUES, ierr)

                L = 2
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, (two/three)*factor, INSERT_VALUES, ierr)

             end do
            end do
          end do
        end if
      end do

    case (3) ! radial current of species
      ! This should only occur with whichSpecies == 0
      do ispecies=1,Nspecies
        THat = THats(ispecies)
        mHat = mHats(ispecies)
        nHat = nHats(ispecies)
        sqrtTHat = sqrt(THat)
        sqrtMHat = sqrt(mHat)

        do ix=ixMin,Nx
          xPartOfRHS = exp(-x2(ix))*2*Zs(ispecies)*nHat*mHat*sqrtMHat*x(ix)/ &
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
      enddo

  end select

  ! Done inserting values.
  ! Finally, assemble the RHS vector:
  call VecAssemblyBegin(rhs, ierr)
  call VecAssemblyEnd(rhs, ierr)


end subroutine populateAdjointRHS
