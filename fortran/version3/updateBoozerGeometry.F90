#include "../PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscsnesdef.h>
#endif

!> This subroutine is used to update BHat, BHat_sub_theta, 
!! BHat_sub_zeta, BHat_sup_theta, BHat_sup_zeta, and DHat
!! for finite difference testing of the adjoint subroutines. 
!! @param whichMode Which mode to differentiate with respect to (index of ms and and ns).
!! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
!! @param deltaLambda Finite difference step size.
subroutine updateBoozerGeometry(whichMode, reset)

  use globalVariables
  use geometry
  use petscvec

  implicit none

  integer :: whichMode
  integer :: itheta, izeta, m, n
  PetscScalar :: angle, cos_angle, sin_angle
  logical :: reset

  PetscScalar, dimension(:,:), allocatable :: deltaBHat, deltadBHatdtheta, deltadBHatdzeta, deltaBHat_sup_theta, deltadBHat_sup_theta_dzeta, deltaBHat_sup_zeta, deltadBHat_sup_zeta_dtheta, deltaBHat_sub_theta, deltadBHat_sub_theta_dzeta, deltaBHat_sub_zeta, deltadBHat_sub_zeta_dtheta, deltaDHat

  m = ms(whichMode)
  n = ns(whichMode)

  if (reset) then
    ! BHat_sub_theta, BHat_sub_zeta are fixed
    DHat = DHat_init
    BHat = BHat_init
    dBHatdtheta = dBHatdtheta_init
    dBHatdzeta = dBHatdzeta_init
    BHat_sup_theta = BHat_sup_theta_init
    dBHat_sup_theta_dzeta = dBHat_sup_theta_dzeta_init
    BHat_sup_zeta = BHat_sup_zeta_init
    dBHat_sup_zeta_dtheta = dBHat_sup_zeta_dtheta_init
  else

    allocate(deltaBHat(Ntheta,Nzeta))
    allocate(deltadBHatdtheta(Ntheta,Nzeta))
    allocate(deltadBHatdzeta(Ntheta,Nzeta))
    allocate(deltaBHat_sup_theta(Ntheta,Nzeta))
    allocate(deltadBHat_sup_theta_dzeta(Ntheta,Nzeta))
    allocate(deltaBHat_sup_zeta(Ntheta,Nzeta))
    allocate(deltadBHat_sup_zeta_dtheta(Ntheta,Nzeta))
    allocate(deltaDHat(Ntheta,Nzeta))

    ! Only BHat is changed
    do itheta = 1,Ntheta
      do izeta = 1,Nzeta
        angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
        cos_angle = cos(angle)
        sin_angle = sin(angle)
        deltaBHat(itheta,izeta) = deltaLambda*bmnc(whichMode)*cos_angle
        deltadBHatdtheta(itheta,izeta) = -deltaLambda*bmnc(whichMode)*m*sin_angle
        deltadBHatdzeta(itheta,izeta) = deltaLambda*bmnc(whichMode)*n*Nperiods*sin_angle
        deltadBHat_sup_theta_dzeta(itheta,izeta) = dBHat_sup_theta_dzeta(itheta,izeta)*(deltaBHat(itheta,izeta)/BHat(itheta,izeta) + deltadBHatdzeta(itheta,izeta)/dBHatdzeta(itheta,izeta))
        deltadBHat_sup_zeta_dtheta(itheta,izeta) = dBHat_sup_zeta_dtheta(itheta,izeta)*(deltaBHat(itheta,izeta)/BHat(itheta,izeta) + deltadBHatdtheta(itheta,izeta)/dBHatdtheta(itheta,izeta))
      end do
    end do
    ! Update BHat
    BHat = BHat + deltaBHat
    dBHatdtheta = dBHatdtheta + deltadBHatdtheta
    dBHatdzeta = dBHatdzeta + deltadBHatdzeta
    DHat = BHat*BHat/(GHat+iota*IHat)
    BHat_sup_theta = iota*DHat
    BHat_sup_zeta = DHat
    dBHat_sup_theta_dzeta = dBHat_sup_theta_dzeta + deltadBHat_sup_theta_dzeta
    dBHat_sup_zeta_dtheta = dBHat_sup_zeta_dtheta + deltadBHat_sup_zeta_dtheta

    deallocate(deltaBHat)
    deallocate(deltadBHatdtheta)
    deallocate(deltadBHatdzeta)
    deallocate(deltadBHat_sup_theta_dzeta)
    deallocate(deltadBHat_sup_zeta_dtheta)
  end if

  ! Update B integrals
  call computeBIntegrals()

end subroutine updateBoozerGeometry
