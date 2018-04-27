#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#endif

! This subroutine is used for testing 
! build_adjointECProlongationMatrix
subroutine testingInterpolation()

  use globalVariables
  use indices
  use petscvec

  implicit none

  Vec :: coarse_vec, coarse_interp
  integer :: ispecies,ix,itheta,izeta,L,ixMin
  PetscScalar :: value, index, norm
  PetscErrorCode :: ierr

  if (pointAtX0) then
     ixMin = 2
  else
     ixMin = 1
  end if

  ! Allocate coarse_vec
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, coarse_vec, ierr)
  call VecSet(coarse_vec, zero, ierr)

  ! Populate coarse_vec
  do ispecies=1,Nspecies
    do ix=ixMin,Nx
      do itheta= ithetaMin,ithetaMax
        do izeta= izetaMin,izetaMax
          do L=0,Nxi_for_x(ix)-1
            value = sin(theta(itheta))
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)
            call VecSetValue(coarse_vec, index, value, INSERT_VALUES, ierr)
          end do
        end do
      end do
    end do
  end do

  ! Assemble coarse_vec
  call VecAssemblyBegin(coarse_vec, ierr)
  call VecAssemblyEnd(coarse_vec, ierr)

  ! Create coarse_interp Vec
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, coarse_interp, ierr)
  call VecSet(coarse_interp,zero,ierr)
  ! Interpolate coarse_vec
  call MatMult(adjointECProlongationMatrix, coarse_vec, coarse_interp, ierr)

  call VecNorm(coarse_vec,NORM_1,norm,ierr)
  print *,"norm: ", (norm)*thetaWeights(1)
  call VecNorm(coarse_interp,NORM_1,norm,ierr)
  print *,"norm: ", (norm)*thetaWeights_fine(1)

  call VecSet(coarse_vec, zero, ierr)

  ! Populate coarse_vec
  do ispecies=1,Nspecies
    do ix=ixMin,Nx
      do itheta= ithetaMin,ithetaMax
        do izeta= izetaMin,izetaMax
          do L=0,Nxi_for_x(ix)-1
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)
            call VecSetValue(coarse_vec, index, value, INSERT_VALUES, ierr)
          end do
        end do
      end do
    end do
  end do

  ! Assemble coarse_vec
  call VecAssemblyBegin(coarse_vec, ierr)
  call VecAssemblyEnd(coarse_vec, ierr)


end subroutine testingInterpolation
