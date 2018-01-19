#include "../PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

!> This subroutine is used to test the explicit derivatives of
!! the particle flux, heat flux, and flow on geometry with
!! finite difference derivatives. This has been written to test
!! the subroutines heatFluxSensitivity, particleFluxSensitivity,
!! and parallelFlowSensitivity.
!! @param forwardSolution Petsc Vec solution to forward equation. Derivatives
!! are taken with this solution held constant.
!! @whichMode Which mode to differentiate with respect to (index of ms and and ns).
!! @whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
!! @deltaLambda Finite difference step size.
subroutine testingDiagnosticSensitivity(forwardSolution, whichMode, whichLambda, deltaLambda)

  use globalVariables
  use petscvec
  use adjointDiagnostics
  use geometry

  implicit none

  Vec :: forwardSolution, forwardSolutionOnProc0
  integer :: whichMode, whichLambda
  PetscScalar :: deltaLambda
  PetscScalar, dimension(:,:), allocatable :: deltaBHat, deltaBHat_sup_theta, deltaBHat_sup_zeta, deltaBHat_sub_theta, deltaBHat_sub_zeta, deltaDHat
  PetscScalar, dimension(:,:), allocatable :: deltadBHatdtheta, deltadBHatdzeta, deltadBHat_sub_zeta_dtheta, deltadBHat_sub_theta_dzeta, deltadBHat_sup_zeta_dtheta, deltadBHat_sup_theta_dzeta
  integer :: m, n, itheta,izeta, iterationNum, whichSpecies
  PetscScalar :: angle, cos_angle, sin_angle
  PetscScalar, dimension(:), allocatable :: particleFluxInit, heatFluxInit, parallelFlowInit
  PetscScalar :: finiteDiffDerivative
  VecScatter :: VecScatterContext
  PetscScalar, pointer :: forwardSolutionArray(:)
  PetscErrorCode :: ierr
  PetscScalar :: percentError
  PetscScalar :: dParticleFluxdLambda_analytic, dHeatFluxdLambda_analytic, dparallelFlowdLambda_analytic

  allocate(deltaBHat(Ntheta,Nzeta))
  allocate(deltaBHat_sub_theta(Ntheta,Nzeta))
  allocate(deltaBHat_sub_zeta(Ntheta,Nzeta))
  allocate(deltaBHat_sup_theta(Ntheta,Nzeta))
  allocate(deltaBHat_sup_zeta(Ntheta,Nzeta))
  allocate(deltaDHat(Ntheta,Nzeta))
  allocate(deltadBHatdtheta(Ntheta,Nzeta))
  allocate(deltadBHatdzeta(Ntheta,Nzeta))
  allocate(deltadBHat_sub_theta_dzeta(Ntheta,Nzeta))
  allocate(deltadBHat_sub_zeta_dtheta(Ntheta,Nzeta))
  allocate(deltadBHat_sup_theta_dzeta(Ntheta,Nzeta))
  allocate(deltadBHat_sup_zeta_dtheta(Ntheta,Nzeta))

  deltaBHat = zero
  deltaBHat_sub_theta = zero
  deltaBHat_sub_zeta = zero
  deltaBHat_sup_theta = zero
  deltaBHat_sup_zeta = zero
  deltaDHat = zero
  deltadBHatdtheta = zero
  deltadBHatdzeta = zero
  deltadBHat_sub_theta_dzeta = zero
  deltadBHat_sub_zeta_dtheta = zero
  deltadBHat_sup_theta_dzeta = zero
  deltadBHat_sup_zeta_dtheta = zero

  allocate(particleFluxInit(Nspecies))
  allocate(heatFluxInit(Nspecies))
  allocate(parallelFlowInit(Nspecies))

  m = ms(whichMode)
  n = ns(whichMode)

  ! Iteration number for computing diagnostics
  iterationNum = 2

  !> Create a scattering context for forwardSolution
  call VecScatterCreateToZero(forwardSolution, VecScatterContext, forwardSolutionOnProc0, ierr)
  !> Send forwardSolution to master proc
  call VecScatterBegin(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  if (masterProc) then
    ! Convert the PETSc vector into a normal Fortran array
    call VecGetARrayF90(forwardSolutionOnProc0, forwardSolutionArray, ierr)
  end if

  ! Get initial values of fluxes
  particleFluxInit = particleFlux_vm_rN
  heatFluxInit = heatFlux_vm_rN
  parallelFlowInit = FSABVelocityUsingFSADensityOverRootFSAB2

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

  case (5) ! BHat_sub_zeta ! not working
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

  ! Compute diagnostics with new geometry
  call diagnostics(forwardSolution, iterationNum)

  ! Reset geometry to original values
  select case(whichLambda)
  case(1) ! BHat
    BHat = BHat - deltaBHat
    dBHatdtheta = dBHatdtheta - deltadBHatdtheta
    dBHatdzeta = dBHatdzeta - deltadBHatdzeta
  case(2) ! BHat_sup_theta
    BHat_sup_theta = BHat_sup_theta - deltaBHat_sup_theta
    dBHat_sup_theta_dzeta = dBHat_sup_theta_dzeta - deltadBHat_sup_theta_dzeta
  case(3) ! BHat_sup_zeta
    BHat_sup_zeta = BHat_sup_zeta - deltaBHat_sup_zeta
    dBHat_sup_zeta_dtheta = dBHat_sup_zeta_dtheta - deltadBHat_sup_zeta_dtheta
  case(4) ! BHat_sub_theta
    BHat_sub_theta = BHat_sub_theta - deltaBHat_sub_theta
    dBHat_sub_theta_dzeta = dBHat_sub_theta_dzeta - deltadBHat_sub_theta_dzeta
  case (5) ! BHat_sub_zeta
    BHat_sub_zeta = BHat_sub_zeta - deltaBHat_sub_zeta
    dBHat_sub_zeta_dtheta = dBHat_sub_zeta_dtheta - deltadBHat_sub_zeta_dtheta
  case (6) ! DHat
    DHat = DHat - deltaDHat
  end select
  ! Update B integrals
  call computeBIntegrals()

  percentError = zero
  if (masterProc) then
    do whichSpecies = 1,NSpecies
      call parallelFlowSensitivity(dparallelFlowdLambda_analytic, forwardSolutionArray, whichSpecies, whichLambda, whichMode)
      call heatFluxSensitivity(dHeatFluxdLambda_analytic, forwardSolutionArray, whichSpecies, whichLambda, whichMode)
      call particleFluxSensitivity(dParticleFluxdLambda_analytic, forwardSolutionArray, whichSpecies, whichLambda, whichMode)
      print "(a,i4,a)","Benchmarking fluxes for ispecies: ", whichSpecies," -----------------------------"
      finiteDiffDerivative = (particleFlux_vm_rN(whichSpecies)-particleFluxInit(whichSpecies))/deltaLambda
      if (abs(finiteDiffDerivative) < 1e-16 .and. abs(dparticleFluxdLambda_analytic) < 1e-16) then
        percentError = zero
      else if (abs(finiteDiffDerivative) > 1e-16) then
        percentError = 100*abs(dParticleFluxdLambda_analytic-finiteDiffDerivative)/abs(finiteDiffDerivative)
      else
        percentError = 1e6
      end if
      if (percentError > 1.0) then
        print "(a,es14.7,a)","percent error: ", percentError,"%"
        print "(a,es14.7)","dparticleFluxdLambda (finite diff): ", finiteDiffDerivative
        print "(a,es14.7)","dparticleFluxdLambda (analytic): ", dParticleFluxdLambda_analytic
      end if

      finiteDiffDerivative = (heatFlux_vm_rN(whichSpecies)-heatFluxInit(whichSpecies))/deltaLambda
      if (abs(finiteDiffDerivative) < 1e-16 .and. abs(dHeatFluxdLambda_analytic) < 1e-16) then
        percentError = zero
      else if (abs(finiteDiffDerivative) > 1e-16) then
        percentError = 100*abs(dHeatFluxdLambda_analytic-finiteDiffDerivative)/abs(finiteDiffDerivative)
      else
        percentError = 1e6
      end if

      if (percentError > 1.0) then
        print "(a,es14.7,a)","percent error: ", percentError,"%"
        print "(a,es14.7)","dheatFluxdLambda (finite diff): ", finiteDiffDerivative
        print "(a,es14.7)","dheatFluxdLambda (analytic): ", dHeatFluxdLambda_analytic
      end if

      finiteDiffDerivative = (FSABVelocityUsingFSADensityOverRootFSAB2(whichSpecies)-parallelFlowInit(whichSpecies))/deltaLambda
      if (abs(finiteDiffDerivative) < 1e-16 .and. abs(dparallelFlowdLambda_analytic) < 1e-16) then
        percentError = zero
      else if (abs(finiteDiffDerivative) > 1e-16) then
        percentError = 100*abs(dparallelFlowdLambda_analytic-finiteDiffDerivative)/abs(finiteDiffDerivative)
      else
        percentError = 1e6
      end if

      if (percentError > 1.0) then
        print "(a,es14.7,a)","percent error: ", percentError,"%"
        print "(a,es14.7)","dparallelFlowdLambda (finite diff): ", finiteDiffDerivative
        print "(a,es14.7)","dparallelFlowdLambda (analytic): ", dparallelFlowdLambda_analytic
      end if
    end do

    ! Compute diagnostics with old geometry for testing
    call diagnostics(forwardSolution, iterationNum)
!
!    do whichSpecies=1,Nspecies
!      print *,"Delta particleFlux (should be 0): ", particleFlux_vm_rN(whichSpecies)-particleFluxInit(whichSpecies)
!      print *,"Delta heatFlux (should be 0): ",heatFlux_vm_rN(whichSpecies)-heatFluxInit(whichSpecies)
!      print *,"Delta parallelFlow (should be 0):",FSABVelocityUsingFSADensityOverRootFSAB2(whichSpecies)-parallelFlowInit(whichSpecies)
!    end do
  end if

end subroutine testingDiagnosticSensitivity
