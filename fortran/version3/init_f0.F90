  subroutine init_f0()

#include "PETScVersions.F90"
    
    use globalVariables
    use petscmat
    use indices

    implicit none

    integer :: L, ix, itheta, izeta, ispecies, index
    PetscScalar :: factor
    PetscErrorCode :: ierr

    if (masterProc) then
       print *,"Initializing f0"
    end if

    call VecSet(f0, zero, ierr)
    
    L = 0
    do ispecies = 1,Nspecies
       do ix = 1,Nx
          factor = nHats(ispecies)*mHats(ispecies)/(pi*THats(ispecies)) &
               * sqrt(mHats(ispecies)/(pi*THats(ispecies))) * expx2(ix)
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(f0, index, &
                     exp(-Zs(ispecies)*alpha*Phi1Hat(itheta,izeta)/THats(ispecies))*factor, INSERT_VALUES, ierr) !!Added by AM 2016-06
             end do
          end do
       end do
    end do

    call VecAssemblyBegin(f0, ierr)
    call VecAssemblyEnd(f0, ierr)
   

  end subroutine init_f0
