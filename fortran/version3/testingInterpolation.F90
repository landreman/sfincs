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

!  if (pointAtX0) then
!     ixMin = 2
!  else
!     ixMin = 1
!  end if
!
!  ! Allocate coarse_vec
!  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, coarse_vec, ierr)
!  call VecSet(coarse_vec, zero, ierr)
!
!  ! Populate coarse_vec
!  do ispecies=1,Nspecies
!    do ix=ixMin,Nx
!      do itheta= ithetaMin,ithetaMax
!        do izeta= izetaMin,izetaMax
!          do L=0,Nxi_for_x(ix)-1
!            value = sin(theta(itheta))
!            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)
!            call VecSetValue(coarse_vec, index, value, INSERT_VALUES, ierr)
!          end do
!        end do
!      end do
!    end do
!  end do
!
!  ! Assemble coarse_vec
!  call VecAssemblyBegin(coarse_vec, ierr)
!  call VecAssemblyEnd(coarse_vec, ierr)
!
!  ! Create coarse_interp Vec
!  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, coarse_interp, ierr)
!  call VecSet(coarse_interp,zero,ierr)
!  ! Interpolate coarse_vec -> coarse_interp
!  call periodicCubicInterpolation(n)
!!  call MatMult(adjointECProlongationMatrix, coarse_vec, coarse_interp, ierr)
!
!  call VecSum(coarse_vec,norm,ierr)
!  print *,"sum(coarse_vec)/DKE_size: ", norm/(sum(Nxi_for_x)*Ntheta*Nzeta)
!  call VecSum(coarse_interp,norm,ierr)
!  print *,"sum(coarse_interp)/DKE_size_fine: ", norm/(sum(Nxi_for_x)*Ntheta_fine*Nzeta_fine)
!
!  call VecSet(coarse_vec, zero, ierr)

  allocate(sin_interp(Ntheta_fine))
  if (adjointECInterpOption==1) then ! Linear interpolation
    allocate(theta_prolongation(Ntheta_fine,Ntheta))
    call  periodic_interpolation(Ntheta,Ntheta_fine,2*pi,theta_fine,theta_prolongation)
    sin_interp = matmul(theta_prolongation,sin(theta))
  else ! Cubic spline interpolation
    call periodicCubicInterpolation(Ntheta,Ntheta_fine,2*pi,theta,theta_fine,sin(theta),sin_interp)
  end if
!  allocate(zeta_prolongation(Nzeta_fine,Nzeta))
!  call periodic_interpolation(Nzeta,Nzeta_fine,zetaMax,zeta_fine,zeta_prolongation)

  ! Check geometry arrays
!  print *,"mean(BHat_fine): ", sum(BHat_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(BHat): ", sum(BHat)/(Ntheta*Nzeta)
!  print *,"mean(DHat_fine): ", sum(DHat_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(DHat): ", sum(DHat)/(Ntheta*Nzeta)
!  print *,"mean(dBHatdtheta_fine): ", sum(dBhatdtheta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(dBhatdtheta): ", sum(dBHatdtheta)/(Ntheta*Nzeta)
!  print *,"mean(dBHatdzeta_fine): ", sum(dBhatdzeta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(dBhatdzeta): ", sum(dBHatdzeta)/(Ntheta*Nzeta)
!  print *,"mean(BHat_sub_theta_fine): ", sum(BHat_sub_theta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(BHat_sub_theta): ", sum(BHat_sub_theta)/(Ntheta*Nzeta)
!  print *,"mean(dBHat_sub_theta_dzeta_fine): ", sum(dBHat_sub_theta_dzeta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(dBHat_sub_theta_dzeta): ", sum(dBHat_sub_theta_dzeta)/(Ntheta*Nzeta)
!  print *,"mean(BHat_sub_zeta_fine): ", sum(BHat_sub_zeta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(BHat_sub_zeta): ", sum(BHat_sub_zeta)/(Ntheta*Nzeta)
!  print *,"mean(dBHat_sub_zeta_dtheta_fine): ", sum(dBHat_sub_zeta_dtheta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(dBHat_sub_zeta_dtheta): ", sum(dBHat_sub_zeta_dtheta)/(Ntheta*Nzeta)
!  print *,"mean(BHat_sup_theta_fine): ", sum(BHat_sup_theta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(BHat_sup_theta): ", sum(BHat_sup_theta)/(Ntheta*Nzeta)
!  print *,"mean(dBHat_sup_theta_dzeta_fine): ", sum(dBHat_sup_theta_dzeta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(dBHat_sup_theta_dzeta): ", sum(dBHat_sup_theta_dzeta)/(Ntheta*Nzeta)
!  print *,"mean(BHat_sup_zeta_fine): ", sum(BHat_sup_zeta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(BHat_sup_zeta): ", sum(BHat_sup_zeta)/(Ntheta*Nzeta)
!  print *,"mean(dBHat_sup_zeta_dtheta_fine): ", sum(dBHat_sup_zeta_dtheta_fine)/(Ntheta_fine*Nzeta_fine)
!  print *,"mean(dBHat_sup_zeta_dtheta): ", sum(dBHat_sup_zeta_dtheta)/(Ntheta*Nzeta)
!
!  ! Testing weights
!  print *,"sum(sin(theta)^2)*thetaWeights: ", sum(sin(theta)*sin(theta)*thetaWeights)
!  print *,"sum(sin(theta_fine)^2)*thetaWeights_fine: ", sum(sin(theta_fine)*sin(theta_fine)*thetaWeights_fine)
!  print *,"sum(sin(zeta)^2)*zetaWeights: ", sum(sin(zeta)*sin(zeta)*zetaWeights)
!  print *,"sum(sin(zeta_fine)^2)*zetaWeights_fine: ", sum(sin(zeta_fine)*sin(zeta_fine)*zetaWeights_fine)

end subroutine testingInterpolation
