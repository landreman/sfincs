#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

subroutine build_adjointECProlongationMatrix()

  use petscmat
  use indices
  use globalVariables
  use xGrid

  implicit none

  PetscScalar, allocatable, dimension(:,:) :: theta_prolongation, zeta_prolongation, x_prolongation
  integer :: itheta_fine, itheta_coarse, nnz_per_row
  PetscErrorCode :: ierr
  integer :: ispecies, ix_coarse, L_coarse, izeta_coarse, izeta_fine
  PetscScalar :: theta_value, zeta_value
  !! Testing  
  PetscScalar, allocatable, dimension(:,:) :: DHat_interp
  integer :: itheta, izeta

  if (masterProc) then
    print *,"Assembling prolongation matrix for adjoint EC."
  end if

  ! First, generate prolongation matrices for each of the coordinates individually:
  allocate(theta_prolongation(Ntheta_fine,Ntheta))
  call  periodic_interpolation(Ntheta, Ntheta_fine,2*pi,theta_fine,theta_prolongation)
  allocate(zeta_prolongation(Nzeta_fine,Nzeta))
  call periodic_interpolation(Nzeta,Nzeta_fine,zetaMax,zeta_fine,zeta_prolongation)

  ! Initialize the 3D prolongation matrix
  call MatCreate(PETSC_COMM_WORLD, adjointECProlongationMatrix, ierr)
  call MatSetType(adjointECProlongationMatrix, MATAIJ, ierr)
  call MatSetSizes(adjointECProlongationMatrix, PETSC_DECIDE, PETSC_DECIDE, &
    matrixSize_fine, matrixSize, ierr)

  ! Allocate 8 nonzero entries per row:
  ! This is a very conservative estimate: 2 off-diagonal elements in theta, zeta, x, and xi: 2*2*2=8.
  ! This estimate could surely be made tighter if it mattered significantly, which it probably does not.
  nnz_per_row = 4
  if (numProcs == 1) then
     call MatSeqAIJSetPreallocation(adjointECProlongationMatrix, nnz_per_row, PETSC_NULL_INTEGER, ierr)
  else
     call MatMPIAIJSetPreallocation(adjointECProlongationMatrix, nnz_per_row, PETSC_NULL_INTEGER, nnz_per_row, PETSC_NULL_INTEGER, ierr)
  end if

  ! Transfer the source/constraint elements without change:
  print *,"Transferring source/constraints."
  if (procThatHandlesConstraints) then
    do ispecies = 1,Nspecies
       call MatSetValue(adjointECProlongationMatrix, &
            getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT,7), &
            getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT,0), &
            one, ADD_VALUES, ierr)
       call MatSetValue(adjointECProlongationMatrix, &
            getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT,7), &
            getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT,0), &
            one, ADD_VALUES, ierr)
    end do
  end if

!  print *,"Populating matrix."
  ! Now populate the 3D prolongation matrix.
  ! Parallelize the looping over rows (fine grid) but not the looping over columns (coarse grid)
  ! No prolongation in xi or x as of yet
!  do itheta_coarse = 1,Ntheta
    do itheta_fine = 1,Ntheta_fine
   do itheta_coarse  = ithetaMin,ithetaMax
!     do itheta_fine = ithetaMin_fine,ithetaMax_fine
      theta_value = theta_prolongation(itheta_fine,itheta_coarse)
      if (abs(theta_value)<1e-12) cycle
!      do izeta_coarse = 1,Nzeta
        do izeta_fine = 1,Nzeta_fine
        do izeta_coarse = izetaMin,izetaMax
!           do izeta_fine = izetaMin_fine,izetaMax_fine
          zeta_value = zeta_prolongation(izeta_fine, izeta_coarse)
          if (abs(zeta_value)<1e-12) cycle
          do ix_coarse = 1,Nx
!            do ix_fine = 1,Nx_fine
!              x_value = x_prolongation(ix_fine,ix_coarse)
!              x_value = 1
!              if (abs(x_value)<1e-12) cycle
              do L_coarse=0,(Nxi_for_x(ix_coarse)-1)
!                do L_fine=0,(Nxi_for_x_fine(ix_fine)-1)
                  do ispecies = 1,Nspecies
                    call MatSetValue(adjointECProlongationMatrix, &
                            getIndex(ispecies,ix_coarse,L_coarse+1,itheta_fine,izeta_fine,BLOCK_F,7), &
                            getIndex(ispecies,ix_coarse,L_coarse+1,itheta_coarse,izeta_coarse,BLOCK_F,0), &
                            theta_value*zeta_value, ADD_VALUES, ierr)
                  end do ! ispecies
!                end do ! L_fine
              end do ! L_coarse
!            end do ! ix_fine
          end do ! ix_coarse
        end do ! izeta_fine
      end do ! izeta_coarse
    end do ! itheta_fine
  end do ! itheta_coarse

  call MatAssemblyBegin(adjointECProlongationMatrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(adjointECProlongationMatrix, MAT_FINAL_ASSEMBLY, ierr)

  deallocate(theta_prolongation, zeta_prolongation)

end subroutine build_adjointECProlongationMatrix
