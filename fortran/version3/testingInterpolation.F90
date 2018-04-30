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
  integer :: ispecies,ix,itheta,izeta,L,ixMin,index
  PetscScalar :: value, norm
  PetscErrorCode :: ierr
  PetscScalar, dimension(:,:), allocatable :: theta_prolongation
  Vec :: coarse_interpOnProc0
  VecScatter :: VecScatterContext
  integer, dimension(:), allocatable :: thetaIndices
  PetscScalar, pointer :: coarse_interpArray(:)

  allocate(sin_interp(Ntheta_fine))

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
  ! Interpolate coarse_vec -> coarse_interp
  if (adjointECInterpOption==1) then
    call MatMult(adjointECProlongationMatrix, coarse_vec, coarse_interp, ierr)
  else
    call prolongateVecCubic(coarse_vec,coarse_interp)
  end if

  ! Scatter coarse_interp to masterProc
  call VecScatterCreateToZero(coarse_interp, VecScatterContext, coarse_interpOnProc0, ierr)
  call VecScatterBegin(VecScatterContext, coarse_interp, coarse_interpOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, coarse_interp, coarse_interpOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  if (masterProc) then
    ! Convert the PETSc vector into a normal Fortran array
    call VecGetArrayF90(coarse_interpOnProc0, coarse_interpArray, ierr)
    ! Extract sin_interp
    ispecies = Nspecies
    ix = Nx
    izeta = Nzeta
    L = 0
    allocate(thetaIndices(Ntheta_fine))
    do itheta=1,Ntheta_fine
      thetaIndices(itheta) = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,7)
    end do
    call VecGetValues(coarse_interpOnProc0,Ntheta_fine,thetaIndices,sin_interp,ierr)
  end if

end subroutine testingInterpolation
