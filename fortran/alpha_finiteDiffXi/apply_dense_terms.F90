#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif


  subroutine apply_dense_terms(inputVec, outputVec, whichMatrix)

    use petscvec
    use globalVariables 
    use indices

    implicit none

    Vec :: inputVec, outputVec
    integer, intent(in) :: whichMatrix

    PetscErrorCode :: ierr
    VecScatter :: VecScatterContext
    Vec :: inputVecLocal
    PetscScalar, pointer :: stateArray(:)
    PetscLogDouble :: time1, time2, time3, time4, time5, time6
    integer :: iSpeciesA, iSpeciesB, L, ialpha, izeta, ixi, ix_row, ix_col, index
    PetscScalar :: factor
    PetscScalar, dimension(:), allocatable :: xi_vector_in, xi_vector_out, species_factor

    ! Only execute the rest of this subroutine if we need to:
    if (collisionOperator .ne. 0 .or. preconditioner_field_term_xi_option==0) return

    if (masterProc) print *,"Applying dense terms for whichMatrix=",whichMatrix
    call PetscTime(time1,ierr)

    ! Send the entire input vector to each proc. I believe this is typically a small amount of memory
    ! compared to the memory used on each proc by the mumps factorization. Perhaps we
    ! could eventually move to a more efficient approach in which each proc does not 'know'
    ! the entire input vector.
    call VecScatterCreateToAll(inputVec, VecScatterContext, inputVecLocal, ierr)
    call VecScatterBegin(VecScatterContext, inputVec, inputVecLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, inputVec, inputVecLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterDestroy(VecScatterContext, ierr)

!    if (iFourierMax>=iFourierMin) then ! For procs that do not own any Fourier indices, do not try allocating arrays of size 0 or calling BLAS with
       ! arrays of size 0. This causes errors.
    
    call VecGetArrayF90(inputVecLocal, stateArray, ierr)

    allocate(xi_vector_in(Nxi))
    allocate(xi_vector_out(Nxi))
    allocate(species_factor(Nspecies))

    do iSpeciesA = 1,Nspecies
       species_factor(iSpeciesA) = sqrt(mHats(ispeciesA)/THats(ispeciesA))
    end do
       
    ! *********************************************************
    ! Add the collision operator
    ! *********************************************************
       
    call PetscTime(time3,ierr)
    do ialpha = ialphaMin,ialphaMax
       do izeta = izetaMinDKE,izetaMaxDKE
          !factor = -nu_n*BHat(ialpha,izeta)*species_factor(ispeciesA)/DHat(ialpha,izeta)
          factor = -nu_n*BHat(ialpha,izeta)/abs(DHat(ialpha,izeta))
          do iSpeciesB = 1,Nspecies
             do ix_col = 1,Nx
                do ixi = 1,Nxi
                   index = getIndex(ispeciesB,ix_col,ixi,ialpha,izeta,BLOCK_F) + 1 ! +1 to go from PETsc to Fortran indices
                   xi_vector_in(ixi) = stateArray(index)
                end do
                do iSpeciesA = 1,Nspecies
                   do ix_row = 1,Nx
                      xi_vector_out = 0
                      do L = 0,NL-1
                         call dgemv('n',Nxi,Nxi,species_factor(iSpeciesA)*factor* RosenbluthPotentialTerms(iSpeciesA,iSpeciesB,L+1,ix_row,ix_col),&
                              Legendre_projection(:,:,L+1),&
                              Nxi,xi_vector_in,1,1.0d+0,xi_vector_out,1)
                      end do
                      do ixi = 1,Nxi
                         index = getIndex(ispeciesA,ix_row,ixi,ialpha,izeta,BLOCK_F)
                         call VecSetValue(outputVec,index,xi_vector_out(ixi), ADD_VALUES, ierr)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    call PetscTime(time4,ierr)
    if (masterProc) print *," main loop:   ",time4-time3,"sec"
    time3=time4
    
    call VecRestoreArrayF90(inputVecLocal, stateArray, ierr)
    deallocate(xi_vector_in, xi_vector_out, species_factor)
    
    !end if  ! End test to see if this proc owns any Fourier indices

    call VecAssemblyBegin(outputVec, ierr)
    call VecAssemblyEnd(outputVec, ierr)

    call VecDestroy(inputVecLocal, ierr)

    call PetscTime(time2,ierr)
    if (masterProc) print *,"Total time for apply_dense_terms:",time2-time1

  end subroutine apply_dense_terms
