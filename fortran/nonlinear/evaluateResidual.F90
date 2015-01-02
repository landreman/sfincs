#include <finclude/petscsnesdef.h>
#include "PETScVersions.F90"


  subroutine evaluateResidual(mysnes, stateVec, residualVec, userContext, ierr)

    use petscsnes
    use globalVariables
    use indices

    implicit none

    SNES :: mysnes
    Vec :: stateVec, residualVec
    PetscErrorCode :: ierr
    integer :: userContext(*)
    Vec :: rhs
    PetscScalar :: scalar, xPartOfRHS, factor
    integer :: ix, L, itheta, izeta, ispecies, index
    PetscScalar :: THat, mHat, sqrtTHat, sqrtmHat
    Mat :: residualMatrix
    PetscScalar :: EParallelHatToUse, dPhiHatdpsiHatToUse


    if (masterProc) then
       print *,"evaluateResidual called."
    end if

!!$    PetscScalar :: x, y
!!$    PetscScalar, pointer, dimension(:) :: stateArray, residualArray
!!$
!!$    call VecGetArrayF90(stateVec,stateArray,ierr)
!!$    call VecGetArrayF90(residualVec,residualArray,ierr)
!!$
!!$
!!$    x = stateArray(1)
!!$    y = stateArray(2)
!!$
!!$    residualArray(1) = x*y - x*x + y + 1
!!$    residualArray(2) = exp(x+y-1) - 1
!!$
!!$    print *,"FormFunction called: x=",x,", y=",y, &
!!$         ", residual(1)=",residualArray(1),", residual(2)=",residualArray(2)
!!$
!!$    call VecRestoreArrayF90(stateVec,stateArray,ierr)
!!$    call VecRestoreArrayF90(residualVec,residualArray,ierr)

    call preallocateMatrix(residualMatrix, 2)
    call populateMatrix(residualMatrix, 2)
    call MatMult(residualMatrix, stateVec, residualVec, ierr)

    call MatDestroy(residualMatrix, ierr)

    ! Next, evaluate the "right-hand side", and subtract the result from the residual.

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, rhs, ierr)
    call VecSet(rhs, zero, ierr)

    ! If we ever run with multiple RHS's, these next lines might change:
    EParallelHatToUse = EParallelHat
    dPhiHatdpsiHatToUse = dPhiHatdpsiHat

    ! First add the term arising from radial gradients:
    x2 = x*x
    do ispecies = 1,Nspecies
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)
       
       do ix=1,Nx
          !xPartOfRHS = x2(ix)*exp(-x2(ix))*( dnHatdpsiNs(ispecies)/nHats(ispecies) &
          !     + alpha*Zs(ispecies)/THats(ispecies)*dPhiHatdpsiNToUse &
          !     + (x2(ix) - three/two)*dTHatdpsiNs(ispecies)/THats(ispecies))
          xPartOfRHS = x2(ix)*exp(-x2(ix))*( dnHatdpsiHats(ispecies)/nHats(ispecies) &
               + alpha*Zs(ispecies)/THats(ispecies)*dPhiHatdpsiHatToUse &
               + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THats(ispecies))
          do itheta = ithetaMin,ithetaMax
             do izeta = 1,Nzeta
                
                !factor = Delta*nHats(ispecies)*mHat*sqrtMHat &
                !     /(2*pi*sqrtpi*Zs(ispecies)*psiAHat*(BHat(itheta,izeta)**3)*sqrtTHat) &
                !     *(GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))&
                !     *xPartOfRHS
                factor = Delta*nHats(ispecies)*mHat*sqrtMHat &
                     /(2*pi*sqrtpi*Zs(ispecies)*(BHat(itheta,izeta)**3)*sqrtTHat) &
                     *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                     - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                     * DHat(itheta,izeta) * xPartOfRHS
                
                L = 0
                index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                call VecSetValue(rhs, index, (4/three)*factor, INSERT_VALUES, ierr)
                
                L = 2
                index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                call VecSetValue(rhs, index, (two/three)*factor, INSERT_VALUES, ierr)
             end do
          end do
       end do
    end do

    ! Add the inductive electric field term:
    L=1
    do ispecies = 1,Nspecies
       do ix=1,Nx
          factor = alpha*Zs(ispecies)*x(ix)*exp(-x2(ix))*EParallelHatToUse*(GHat+iota*IHat)&
               *nHats(ispecies)*mHats(ispecies)/(pi*sqrtpi*THats(ispecies)*THats(ispecies)*FSABHat2)
          do itheta=ithetaMin,ithetaMax
             do izeta = 1,Nzeta
                index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                call VecSetValue(rhs, index, &
                     factor/BHat(itheta,izeta), INSERT_VALUES, ierr)
             end do
          end do
       end do
    end do

    ! Done inserting values.
    ! Finally, assemble the RHS vector:
    call VecAssemblyBegin(rhs, ierr)
    call VecAssemblyEnd(rhs, ierr)


    ! Subtract the RHS from the residual:
    scalar = -1
    call VecAXPY(residualVec, scalar, rhs, ierr)
    call VecDestroy(rhs, ierr)

  end subroutine evaluateResidual
