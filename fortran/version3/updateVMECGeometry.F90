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
subroutine updateVMECGeometry(whichMode, whichLambda, deltaLambda)

  use globalVariables
  use geometry
  use petscvec

  implicit none

  integer :: whichMode, whichLambda
  PetscScalar :: deltaLambda
  integer :: itheta, izeta, m, n
  PetscScalar :: angle, cos_angle, sin_angle

  PetscScalar, dimension(:,:), allocatable :: deltaBHat, deltadBHatdtheta, deltadBHatdzeta, deltaBHat_sup_theta, deltadBHat_sup_theta_dzeta, deltaBHat_sup_zeta, deltadBHat_sup_zeta_dtheta, deltaBHat_sub_theta, deltadBHat_sub_theta_dzeta, deltaBHat_sub_zeta, deltadBHat_sub_zeta_dtheta, deltaDHat

  m = ms(whichMode)
  n = ns(whichMode)

  allocate(deltaBHat(Ntheta,Nzeta))
  allocate(deltadBHatdtheta(Ntheta,Nzeta))
  allocate(deltadBHatdzeta(Ntheta,Nzeta))
  allocate(deltaBHat_sup_theta(Ntheta,Nzeta))
  allocate(deltadBHat_sup_theta_dzeta(Ntheta,Nzeta))
  allocate(deltaBHat_sup_zeta(Ntheta,Nzeta))
  allocate(deltadBHat_sup_zeta_dtheta(Ntheta,Nzeta))
  allocate(deltaBHat_sub_theta(Ntheta,Nzeta))
  allocate(deltadBHat_sub_theta_dzeta(Ntheta,Nzeta))
  allocate(deltaBHat_sub_zeta(Ntheta,Nzeta))
  allocate(deltadBHat_sub_zeta_dtheta(Ntheta,Nzeta))
  allocate(deltaDHat(Ntheta,Nzeta))

  select case(whichLambda)
  case(1) ! BHat
    do itheta = 1,Ntheta
      do izeta = 1,Nzeta
        angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
        cos_angle = cos(angle)
        sin_angle = sin(angle)
        deltaBHat(itheta,izeta) = deltaLambda*cos_angle
        deltadBHatdtheta(itheta,izeta) = -deltaLambda*m*sin_angle
        deltadBHatdzeta(itheta,izeta) = deltaLambda*n*Nperiods*sin_angle
      end do
    end do
    BHat = BHat + deltaBHat
    dBHatdtheta = dBHatdtheta + deltadBHatdtheta
    dBHatdzeta = dBHatdzeta + deltadBHatdzeta

  case(2) ! BHat_sup_theta
    do itheta = 1,Ntheta
      do izeta = 1,Nzeta
        angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
        cos_angle = cos(angle)
        sin_angle = sin(angle)
        deltaBHat_sup_theta(itheta,izeta) = deltaLambda*cos_angle
        deltadBHat_sup_theta_dzeta(itheta,izeta) = deltaLambda*n*Nperiods*sin_angle
      end do
    end do
    BHat_sup_theta = BHat_sup_theta + deltaBHat_sup_theta
    dBHat_sup_theta_dzeta = dBHat_sup_theta_dzeta + deltadBHat_sup_theta_dzeta

  case(3) ! BHat_sup_zeta
    do itheta = 1,Ntheta
      do izeta = 1,Nzeta
        angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
        cos_angle = cos(angle)
        sin_angle = sin(angle)
        deltaBHat_sup_zeta(itheta,izeta) = deltaLambda*cos_angle
        deltadBHat_sup_zeta_dtheta(itheta,izeta) = -deltaLambda*m*sin_angle
      end do
    end do
    BHat_sup_zeta = BHat_sup_zeta + deltaBHat_sup_zeta
    dBHat_sup_zeta_dtheta = dBHat_sup_zeta_dtheta + deltadBHat_sup_zeta_dtheta

  case(4) ! BHat_sub_theta
    do itheta = 1,Ntheta
      do izeta = 1,Nzeta
        angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
        cos_angle = cos(angle)
        sin_angle = sin(angle)
        deltaBHat_sub_theta(itheta,izeta) = deltaLambda*cos_angle
        deltadBHat_sub_theta_dzeta(itheta,izeta) = deltaLambda*n*Nperiods*sin_angle
      end do
    end do
    BHat_sub_theta = BHat_sub_theta + deltaBHat_sub_theta
    dBHat_sub_theta_dzeta = dBHat_sub_theta_dzeta + deltadBHat_sub_theta_dzeta

  case (5) ! BHat_sub_zeta
    do itheta = 1,Ntheta
      do izeta = 1,Nzeta
        angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
        cos_angle = cos(angle)
        sin_angle = sin(angle)
        deltaBHat_sub_zeta(itheta,izeta) = deltaLambda*cos_angle
        deltadBHat_sub_zeta_dtheta(itheta,izeta) = -deltaLambda*m*sin_angle
      end do
    end do
    BHat_sub_zeta = BHat_sub_zeta + deltaBHat_sub_zeta
    dBHat_sub_zeta_dtheta = dBHat_sub_zeta_dtheta + deltadBHat_sub_zeta_dtheta

  case (6) ! DHat
    do itheta = 1,Ntheta
      do izeta = 1,Nzeta
        angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
        cos_angle = cos(angle)
        deltaDHat(itheta,izeta) = one/(one/DHat(itheta,izeta) + deltaLambda*cos_angle) - DHat(itheta,izeta)
      end do
    end do
    DHat = DHat + deltaDHat
  end select

  ! Update B integrals
  call computeBIntegrals()

end subroutine updateVMECGeometry
