#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

! See evaluateResidual.f90 (where rhs is constructed). Much of the code has been copied

! dRHSdLambda is a Vec which should be allocated by the calling subroutine
! whichLambda corresponds to the component of B that is going to vary
! whichMode corresponds to the index imn of ns and ms of the desired Fourier mode

  subroutine populatedRHSdLambda(dRHSdLambda, whichLambda, whichMode)

    use petscvec
    use globalVariables
    use indices

    implicit none

    Vec :: dRHSdLambda
    integer :: whichLambda, whichMode

    integer :: ix, L, itheta, izeta, ispecies, index
    integer :: ixMin
    PetscScalar :: THat, mHat, sqrtTHat, sqrtmHat, xPartOfRHS, dFactordLambda
    PetscErrorCode :: ierr

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, dRHSdLambda, ierr)

    if (pointAtX0) then
      ixMin = 2
    else
      ixMin = 1
    end if

    do ispecies = 1, Nspecies
      THat = THats(ispecies)
      mHat = mHats(ispecies)
      sqrtTHat = sqrt(THat)
      sqrtMHat = sqrt(mHat)

      do ix=ixMin,Nx
      !xPartOfRHS is independent of BHat
        xPartOfRHS = x2(ix)*exp(-x2(ix))*( dnHatdpsiHats(ispecies)/nHats(ispecies) &
          + alpha*Zs(ispecies)/THats(ispecies)*dPhiHatdpsiHat &
          + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THats(ispecies))

        do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
            select case(whichLambda)
              case (0) ! Er
                if (masterProc) then
                  print *,"Error! Er sensitivity not yet implemented."
                end if
                stop
              case (1) ! BHat
                ! dBHatdBmn = cos(m*theta - n N_p zeta)
                dFactordLambda = &
                  ! Term from BHat**(-3)
                  -3*Delta*nHats(ispecies)*mHat*sqrtMHat &
                  /(2*pi*sqrtpi*Zs(ispecies)*(BHat(itheta,izeta)**4)*sqrtTHat) &
                  *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                  - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                  * DHat(itheta,izeta) * xPartOfRHS * dBHatdFourier(itheta,izeta,whichMode) &
                  ! Term from dBHatdtheta
                  + Delta*nHats(ispecies)*mHat*sqrtMHat &
                  /(2*pi*sqrtpi*Zs(ispecies)*(BHat(itheta,izeta)**4)*sqrtTHat) &
                  *BHat_sub_zeta(itheta,izeta)* DHat(itheta,izeta) * xPartOfRHS &
                  * dBHatdthetadFourier(itheta,izeta,whichMode) &
                  ! Term from dBHatdzeta
                  - Delta*nHats(ispecies)*mHat*sqrtMHat &
                  /(2*pi*sqrtpi*Zs(ispecies)*(BHat(itheta,izeta)**4)*sqrtTHat) &
                  *BHat_sub_theta(itheta,izeta)* DHat(itheta,izeta) * xPartOfRHS &
                  * dBHatdzetadFourier(itheta,izeta,whichMode)
              case (2) ! BHat_sup_theta
                dFactordLambda = 0
              case (3) ! BHat_sup_zeta
                dFactordLambda = 0
              case (4) ! BHat_sub_theta
                dFactordLambda = -Delta*nHats(ispecies)*mHat*sqrtMHat &
                  /(2*pi*sqrtpi*Zs(ispecies)*(BHat(itheta,izeta)**3)*sqrtTHat) &
                  * dBHatdzeta(itheta,izeta)*dBHat_sub_thetadFourier(itheta,izeta,whichMode) &
                  * DHat(itheta,izeta)*xPartOfRHS
              case (5) ! BHat_sup_zeta
                dFactordLambda = Delta*nHats(ispecies)*mHat*sqrtMHat &
                  /(2*pi*sqrtpi*Zs(ispecies)*(BHat(itheta,izeta)**3)*sqrtTHat) &
                  * dBHatdtheta(itheta,izeta)*dBHat_sub_zetadFourier(itheta,izeta,whichMode) &
                  * DHat(itheta,izeta)*xPartOfRHS
              case (6) ! DHat
                dFactordLambda = Delta*nHats(ispecies)*mHat*sqrtMHat &
                  /(2*pi*sqrtpi*Zs(ispecies)*(BHat(itheta,izeta)**3)*sqrtTHat) &
                  *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                  - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta)) &
                  * dDHatdFourier(itheta,izeta,whichMode) * xPartOfRHS
            end select

            L = 0
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
            call VecSetValue(dRHSdLambda, index, (4/three)*dFactordLambda, INSERT_VALUES, ierr)

            L = 2
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
            call VecSetValue(dRHSdLambda, index, (two/three)*dFactordLambda, INSERT_VALUES, ierr)

          enddo
        enddo
      enddo
    enddo

!    call VecAssemblyBegin(dRHSdLambda, ierr)
!    call VecAssemblyEnd(dRHSdLambda, ierr)

  end subroutine populatedRHSdLambda
