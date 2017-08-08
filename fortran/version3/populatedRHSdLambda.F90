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
    PetscScalar :: dBHatdLambda, dBHatdthetadLambda, dBHatdzetadLambda, dDHatdLambda
    PetscScalar :: dBHat_sub_thetadLambda, dBHat_sub_zetadLambda

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
        xPartOfRHS = (x2(ix)*exp(-x2(ix))*( dnHatdpsiHats(ispecies)/nHats(ispecies) &
          + alpha*Zs(ispecies)/THats(ispecies)*dPhiHatdpsiHat &
          + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THats(ispecies))) &
          *(Delta*nHats(ispecies)*mHat*sqrtMHat/(2*pi*sqrtpi*Zs(ispecies)*sqrtTHat))

!       factor = 1/(BHat(itheta,izeta)**3) &
!            *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
!            - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
!            * DHat(itheta,izeta) * xPartOfRHS

        do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
            select case(whichLambda)
              case (0) ! Er
                if (masterProc) then
                  print *,"Error! Er sensitivity not yet implemented."
                end if
                stop
              case (1) ! BHat
                dBHatdLambda = cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
                dBHatdthetadLambda = -ms(whichMode)*sin(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
                dBHatdzetadLambda = ns(whichMode)*Nperiods*sin(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
                dFactordLambda = -3/(BHat(itheta,izeta)**4) & ! Term from BHat**(-3)
                  *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                  - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                  * DHat(itheta,izeta)*dBHatdLambda &
                  ! Term from dBHatdtheta
                  + DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                  *(BHat_sub_zeta(itheta,izeta)*dBHatdthetadLambda &
                  ! Term from dBHatdzeta
                  - BHat_sub_theta(itheta,izeta)*dBHatdzetadLambda)
              case (2) ! BHat_sup_theta
                dFactordLambda = 0
              case (3) ! BHat_sup_zeta
                dFactordLambda = 0
              case (4) ! BHat_sub_theta
                dBHat_sub_thetadLambda = cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
                dFactordLambda = -DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                  *dBHatdzeta(itheta,izeta)*dBHat_sub_thetadLambda
              case (5) ! BHat_sub_zeta
                dBHat_sub_zetadLambda = cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
                dFactordLambda = DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                  * dBHatdtheta(itheta,izeta)*dBHat_sub_zetadLambda
              case (6) ! DHat
                dDHatdLambda = -DHat(itheta,izeta)*DHat(itheta,izeta)*cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
                dFactordLambda = 1/(BHat(itheta,izeta)**3) &
                  *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                  - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta)) &
                  * dDHatdLambda
            end select

            L = 0
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
            call VecSetValue(dRHSdLambda, index, (4/three)*dFactordLambda*xPartOfRHS, INSERT_VALUES, ierr)

            L = 2
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
            call VecSetValue(dRHSdLambda, index, (two/three)*dFactordLambda*xPartOfRHS, INSERT_VALUES, ierr)

          enddo
        enddo
      enddo
    enddo

    call VecAssemblyBegin(dRHSdLambda, ierr)
    call VecAssemblyEnd(dRHSdLambda, ierr)

  end subroutine populatedRHSdLambda
