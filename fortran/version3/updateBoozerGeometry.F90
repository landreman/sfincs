!> This subroutine is used to update BHat, BHat_sub_theta,
!! BHat_sub_zeta, BHat_sup_theta, BHat_sup_zeta, and DHat
!! for finite difference testing of the adjoint subroutines. 
!! @param whichMode Which mode to differentiate with respect to (index of ms and and ns).
!! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
!! @param deltaLambda Finite difference step size.
subroutine updateBoozerGeometry(whichMode, whichLambda, reset)

#include "PETScVersions.F90"

  use globalVariables
  use geometry
  use petscvec

  implicit none

  integer :: whichMode, whichLambda
  integer :: itheta, izeta, m, n
  PetscScalar :: angle, cos_angle, sin_angle
  logical :: reset

  PetscScalar, dimension(:,:), allocatable :: deltaBHat, deltadBHatdtheta, deltadBHatdzeta, deltaBHat_sup_theta, deltadBHat_sup_theta_dzeta, deltaBHat_sup_zeta, deltadBHat_sup_zeta_dtheta, deltaBHat_sub_theta, deltadBHat_sub_theta_dzeta, deltaBHat_sub_zeta, deltadBHat_sub_zeta_dtheta, deltaDHat

  m = ms_sensitivity(whichMode)
  n = ns_sensitivity(whichMode)

  if (reset) then
    DHat = DHat_init
    BHat = BHat_init
    dBHatdtheta = dBHatdtheta_init
    dBHatdzeta = dBHatdzeta_init
    BHat_sup_theta = BHat_sup_theta_init
    dBHat_sup_theta_dzeta = dBHat_sup_theta_dzeta_init
    BHat_sup_zeta = BHat_sup_zeta_init
    dBHat_sup_zeta_dtheta = dBHat_sup_zeta_dtheta_init
    GHat = GHat_init
    IHat = IHat_init
    iota = iota_init
    BHat_sub_theta = BHat_sub_theta_init
    BHat_sub_zeta = BHat_sub_zeta_init
  else

    allocate(deltaBHat(Ntheta,Nzeta))
    allocate(deltadBHatdtheta(Ntheta,Nzeta))
    allocate(deltadBHatdzeta(Ntheta,Nzeta))
    allocate(deltaBHat_sup_theta(Ntheta,Nzeta))
    allocate(deltadBHat_sup_theta_dzeta(Ntheta,Nzeta))
    allocate(deltaBHat_sup_zeta(Ntheta,Nzeta))
    allocate(deltadBHat_sup_zeta_dtheta(Ntheta,Nzeta))
    allocate(deltaDHat(Ntheta,Nzeta))

    select case (whichLambda)
    case (1) ! BHat
      do itheta = 1,Ntheta
        do izeta = 1,Nzeta
          angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
          cos_angle = cos(angle)
          sin_angle = sin(angle)
          deltaBHat(itheta,izeta) = deltaLambda*bmnc_init(whichMode)*cos_angle
          deltadBHatdtheta(itheta,izeta) = -deltaLambda*bmnc_init(whichMode)*m*sin_angle
          deltadBHatdzeta(itheta,izeta) = deltaLambda*bmnc_init(whichMode)*n*Nperiods*sin_angle
          deltadBHat_sup_theta_dzeta(itheta,izeta) = dBHat_sup_theta_dzeta(itheta,izeta)*(deltaBHat(itheta,izeta)/BHat(itheta,izeta) + deltadBHatdzeta(itheta,izeta)/dBHatdzeta(itheta,izeta))
          deltadBHat_sup_zeta_dtheta(itheta,izeta) = dBHat_sup_zeta_dtheta(itheta,izeta)*(deltaBHat(itheta,izeta)/BHat(itheta,izeta) + deltadBHatdtheta(itheta,izeta)/dBHatdtheta(itheta,izeta))
        end do
      end do
      BHat = BHat + deltaBHat
      dBHatdtheta = dBHatdtheta + deltadBHatdtheta
      dBHatdzeta = dBHatdzeta + deltadBHatdzeta
      DHat = BHat*BHat/(GHat+iota*IHat)
      BHat_sup_theta = iota*DHat
      BHat_sup_zeta = DHat
      dBHat_sup_theta_dzeta = dBHat_sup_theta_dzeta + deltadBHat_sup_theta_dzeta
      dBHat_sup_zeta_dtheta = dBHat_sup_zeta_dtheta + deltadBHat_sup_zeta_dtheta
    case(2) ! IHat
      IHat = IHat*(one+deltaLambda)
      BHat_sub_theta = IHat
      DHat = BHat*BHat/(GHat+iota*IHat)
      BHat_sup_theta = iota*DHat
      BHat_sup_zeta = DHat
      dBHat_sup_theta_dzeta = iota*two*BHat*dBHatdzeta/(GHat+iota*IHat)
      dBHat_sup_zeta_dtheta = two*BHat*dBHatdtheta/(GHat+iota*IHat)
    case(3) ! GHat
      GHat = GHat*(one+deltaLambda)
      BHat_sub_zeta = GHat
      DHat = BHat*BHat/(GHat+iota*IHat)
      BHat_sup_theta = iota*DHat
      BHat_sup_zeta = DHat
      dBHat_sup_theta_dzeta = iota*two*BHat*dBHatdzeta/(GHat+iota*IHat)
      dBHat_sup_zeta_dtheta = two*BHat*dBHatdtheta/(GHat+iota*IHat)
    case(4) ! iota
      iota = iota*(one+deltaLambda)
      DHat = BHat*BHat/(GHat+iota*IHat)
      BHat_sup_theta = iota*DHat
      BHat_sup_zeta = DHat
      dBHat_sup_theta_dzeta = iota*two*BHat*dBHatdzeta/(GHat+iota*IHat)
      dBHat_sup_zeta_dtheta = two*BHat*dBHatdtheta/(GHat+iota*IHat)
    end select

    deallocate(deltaBHat)
    deallocate(deltadBHatdtheta)
    deallocate(deltadBHatdzeta)
    deallocate(deltadBHat_sup_theta_dzeta)
    deallocate(deltadBHat_sup_zeta_dtheta)
  end if

  ! Update B integrals
  call computeBIntegrals()

end subroutine updateBoozerGeometry
