#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

  !> Populates \f$\partial \mathbb{S}/\partial \lambda\f$ to compute explicit dependence of integrated quantiteis on geometry. See evaluateResidual.f90 (where \f$\mathbb{S}\f$ is constructed). Much of the code has been copied.
  !! @param dMatrixdLambda Matrix to be populated. Should be allocated by calling subroutine.
  !! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
  !! @param whichMode Indicates index of ms and ns for derivative.
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
    PetscScalar :: angle, cos_angle, sin_angle, stuffToAdd, factor
    integer :: m,n

    if (whichLambda > 0) then
      m = ms(whichMode)
      n = ns(whichMode)
    end if

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
          *(Delta*nHats(ispecies)*mHat*sqrtMHat/(two*pi*sqrtpi*Zs(ispecies)*sqrtTHat))

       factor = 1/(BHat(itheta,izeta)**3) &
            *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
            - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
            * DHat(itheta,izeta) * xPartOfRHS

        do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            sin_angle = sin(angle)

            select case(whichLambda)
              case (0) ! Er
                stuffToAdd = 1/(BHat(itheta,izeta)**3) &
                  *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                  - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta)) &
                  * DHat(itheta,izeta)*(x2(ix)*exp(-x2(ix))*(alpha*Zs(ispecies)/THats(ispecies))) &
                  *(Delta*nHats(ispecies)*mHat*sqrtMHat/(two*pi*sqrtpi*Zs(ispecies)*sqrtTHat))
              case (1) ! BHat
                dBHatdLambda = cos_angle
                dBHatdthetadLambda = -m*sin_angle
                dBHatdzetadLambda = n*Nperiods*sin_angle
                dFactordLambda = -three/(BHat(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)) & ! Term from BHat**(-3)
                  *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                  - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                  * DHat(itheta,izeta)*dBHatdLambda &
                  ! Term from dBHatdtheta
                  + DHat(itheta,izeta)/(BHat(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                  *(BHat_sub_zeta(itheta,izeta)*dBHatdthetadLambda &
                  ! Term from dBHatdzeta
                  - BHat_sub_theta(itheta,izeta)*dBHatdzetadLambda)
                stuffToAdd = dFactordLambda*xPartOfRHS
                if (geometryScheme /= 5) then ! Boozer
                  dDHatdLambda = 2*DHat(itheta,izeta)*cos_angle/BHat(itheta,izeta)
                  dFactordLambda = dFactordLambda + 1/(BHat(itheta,izeta)**3) &
                    *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                    - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                    * dDHatdLambda
                  stuffToAdd = dFactordLambda*xPartOfRHS
                end if
              case (2) ! BHat_sup_theta
                dFactordLambda = 0
                stuffToAdd = dFactordLambda*xPartOfRHS
              case (3) ! BHat_sup_zeta
                dFactordLambda = 0
                stuffToAdd = dFactordLambda*xPartOfRHS
              case (4) ! BHat_sub_theta
                dBHat_sub_thetadLambda = cos_angle
                dFactordLambda = -DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                  *dBHatdzeta(itheta,izeta)*dBHat_sub_thetadLambda
                stuffToAdd = dFactordLambda*xPartOfRHS
              case (5) ! BHat_sub_zeta
                dBHat_sub_zetadLambda = cos_angle
                dFactordLambda = DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                  * dBHatdtheta(itheta,izeta)*dBHat_sub_zetadLambda
                stuffToAdd = dFactordLambda*xPartOfRHS
              case (6) ! DHat
                dDHatdLambda = -DHat(itheta,izeta)*DHat(itheta,izeta)*cos_angle
                dFactordLambda = one/(BHat(itheta,izeta)**3) &
                  *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                  - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta)) &
                  * dDHatdLambda
                stuffToAdd = dFactordLambda*xPartOfRHS
            end select

            L = 0
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)
            call VecSetValue(dRHSdLambda, index, (four/three)*stuffToAdd, INSERT_VALUES, ierr)

            L = 2
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)
            call VecSetValue(dRHSdLambda, index, (two/three)*stuffToAdd, INSERT_VALUES, ierr)

          enddo
        enddo
      enddo
    enddo

    call VecAssemblyBegin(dRHSdLambda, ierr)
    call VecAssemblyEnd(dRHSdLambda, ierr)

  end subroutine populatedRHSdLambda
