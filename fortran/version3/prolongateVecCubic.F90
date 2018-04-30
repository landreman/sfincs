#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

subroutine prolongateVecCubic(coarseVec,coarse_interpVec)

  use globalVariables
  use indices
  use petscvec

  implicit none

  Vec :: coarseVec, coarse_interpVec
  PetscErrorCode :: ierr
  integer :: itheta,izeta, ix, L,ispecies,index, ixMin, ixMax
  Vec :: coarseVecOnProc0
  VecScatter :: VecScatterContext
  integer, dimension(:), allocatable :: thetaIndices, thetaIndices_fine
  PetscScalar, dimension(:), allocatable :: coarseArray, coarse_interpArray

  allocate(thetaIndices(Ntheta))
  allocate(thetaIndices_fine(Ntheta_fine))
  allocate(coarseArray(Ntheta))
  allocate(coarse_interpArray(Ntheta_fine))

  if (pointAtX0) then
    ixMin = 2
  else
    ixMin = 1
  end if

  ! Scatter coarseVec to masterProc
  call VecScatterCreateToZero(coarseVec, VecScatterContext, coarseVecOnProc0, ierr)
  call VecScatterBegin(VecScatterContext, coarseVec, coarseVecOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, coarseVec, coarseVecOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  if (masterProc) then
    do ispecies=1,Nspecies
      do ix=1,Nx
        do L=1,Nxi_for_x(ix)
          do izeta=1,Nzeta
            do itheta=1,Ntheta
              thetaIndices(itheta) = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)
            end do
            do itheta=1,Ntheta_fine
              thetaIndices_fine(itheta) = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,7)
            end do
            call VecGetValues(coarseVecOnProc0, Ntheta, thetaIndices,coarseArray,ierr)
            call periodicCubicInterpolation(Ntheta,Ntheta_fine,2*pi,theta,theta_fine,coarseArray,coarse_interpArray)
            call VecSetValues(coarse_interpVec,Ntheta_fine,thetaIndices_fine,coarse_interpArray,INSERT_VALUES,ierr)
          end do
        end do
      end do
    end do
  end if

end subroutine prolongateVecCubic
