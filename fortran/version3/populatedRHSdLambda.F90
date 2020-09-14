  !> Populates \f$\partial \mathbb{S}/\partial \lambda\f$ to compute explicit dependence of integrated quantiteis on geometry. See evaluateResidual.f90 (where \f$\mathbb{S}\f$ is constructed). Much of the code has been copied.
  !! @param dMatrixdLambda Matrix to be populated. Should be allocated by calling subroutine.
  !! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
  !! @param whichMode Indicates index of ms and ns for derivative.
  subroutine populatedRHSdLambda(dRHSdLambda, whichLambda, whichMode)

#include "PETScVersions.F90"

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
    PetscScalar :: dBHatdLambda, dBHatdthetadLambda, dBHatdzetadLambda
    PetscScalar :: angle, cos_angle, sin_angle, stuffToAdd
    integer :: m,n

    if (whichLambda > 0) then
      m = ms_sensitivity(whichMode)
      n = ns_sensitivity(whichMode)
    end if

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

         ! factor = 1/(BHat(itheta,izeta)**3) &
         !    *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
         !    - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
         !    * DHat(itheta,izeta) * xPartOfRHS

        ! factor = (1/(BHat*(GHat+iota*IHat)))*(GHat*dBHatdtheta-IHat*dBHatdzeta)*xPartOfRHS
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
                dFactordLambda = -dBHatdLambda/(BHat(itheta,izeta)*BHat(itheta,izeta)*(GHat+iota*IHat)) &
                  * (GHat*dBHatdtheta(itheta,izeta)-IHat*dBHatdzeta(itheta,izeta)) &
                  + (GHat*dBHatdthetadLambda-IHat*dBHatdzetadLambda)/(BHat(itheta,izeta)*(GHat+iota*IHat))
                stuffToAdd = dFactordLambda*xPartOfRHS
              case (2) ! IHat
                dFactordLambda = -iota/(BHat(itheta,izeta)*(GHat+iota*IHat)**2) &
                  *(GHat*dBHatdtheta(itheta,izeta)-IHat*dBHatdzeta(itheta,izeta)) &
                  -dBHatdzeta(itheta,izeta)/(BHat(itheta,izeta)*(GHat+iota*IHat))
                stuffToAdd = dFactordLambda*xPartOfRHS
              case (3) ! GHat
                dFactordLambda = -1/(BHat(itheta,izeta)*(GHat+iota*IHat)**2) &
                  *(GHat*dBHatdtheta(itheta,izeta)-IHat*dBHatdzeta(itheta,izeta)) &
                  + dBHatdtheta(itheta,izeta)/(BHat(itheta,izeta)*(GHat+iota*IHat))
                stuffToAdd = dFactordLambda*xPartOfRHS
              case (4) ! iota
                dFactordLambda = -IHat/(BHat(itheta,izeta)*(GHat+iota*IHat)**2) &
                  *(GHat*dBHatdtheta(itheta,izeta)-IHat*dBHatdzeta(itheta,izeta))
                stuffToAdd = dFactordLambda*xPartOfRHS
            end select

            L = 0
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
            call VecSetValue(dRHSdLambda, index, (four/three)*stuffToAdd, INSERT_VALUES, ierr)

            L = 2
            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
            call VecSetValue(dRHSdLambda, index, (two/three)*stuffToAdd, INSERT_VALUES, ierr)

          enddo
        enddo
      enddo
    enddo

    call VecAssemblyBegin(dRHSdLambda, ierr)
    call VecAssemblyEnd(dRHSdLambda, ierr)

  end subroutine populatedRHSdLambda
