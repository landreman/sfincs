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
!! @param whichMode Which mode to differentiate with respect to (index of ms and and ns).
!! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
!! @param deltaLambda Finite difference step size.
subroutine testingAdjointDiagnostics()

  use globalVariables
  use petscvec
  use adjointDiagnostics
  use geometry
  use solver

  implicit none

!  Vec :: forwardSolution
  integer :: whichMode, whichLambda
  PetscScalar :: deltaLambda
  integer :: iterationNum, whichSpecies
  PetscScalar, dimension(:), allocatable :: particleFluxInit, heatFluxInit, parallelFlowInit
  PetscScalar :: finiteDiffDerivative
  VecScatter :: VecScatterContext
  PetscErrorCode :: ierr
  PetscScalar :: percentError
  PetscScalar, dimension(:,:,:), allocatable :: dParticleFluxdLambda_finiteDiff, dHeatFluxdLambda_finiteDiff, dParallelFlowdLambda_finitediff
  PetscScalar, dimension(:,:,:), allocatable :: dParticleFluxdLambda_analytic, dHeatFluxdLambda_analytic, dParallelFlowdLambda_analytic
  integer :: ispecies, whichQuantity
  PetscScalar :: analyticResult, finiteDiffResult
  PetscLogDouble :: startTime, time1

  allocate(particleFluxInit(Nspecies))
  allocate(heatFluxInit(Nspecies))
  allocate(parallelFlowInit(Nspecies))

  allocate(dParticleFluxdLambda_finiteDiff(Nspecies,NLambdas,NModesAdjoint))
  allocate(dHeatFluxdLambda_finiteDiff(Nspecies,NLambdas,NModesAdjoint))
  allocate(dParallelFlowdLambda_finiteDiff(Nspecies,NLambdas,NModesAdjoint))

  allocate(dParticleFluxdLambda_analytic(Nspecies,NLambdas,NModesAdjoint))
  allocate(dHeatFluxdLambda_analytic(Nspecies,NLambdas,NModesAdjoint))
  allocate(dParallelFlowdLambda_analytic(Nspecies,NLambdas,NModesAdjoint))

  deltaLambda = 1.d-4

  ! Change settings so adjoint solve occurs - derivatives of all fluxes computed
  RHSMode = 4
  adjointParticleFluxOption = .true.
  adjointHeatFluxOption = .true.
  adjointParallelFlowOption = .true.

  call PetscTime(time1, ierr)
  startTime = time1
  call mainSolverLoop()
  call PetscTime(time1, ierr)
  if (masterProc) then
    print *,"Time for adjoint solve: ", time1-startTime
  end if

  ! Get initial values of fluxes
  particleFluxInit = particleFlux_vm_rN
  heatFluxInit = heatFlux_vm_rN
  parallelFlowInit = FSABVelocityUsingFSADensityOverRootFSAB2

  ! Get analytic fluxes from adjoint solve
  dParticleFluxdLambda_analytic = dParticleFluxdLambda
  dHeatFluxdLambda_analytic = dHeatFluxdLambda
  dParallelFlowdLambda_analytic = dParallelFlowdLambda

  ! Set RHSMode = 1 so call to solver does not include adjoint solve
  RHSMode = 1
  call PetscTime(time1, ierr)
  startTime = time1
  do whichMode = 1, 1
    do whichLambda = 1, 1
      do ispecies = 1, 1
        ! Update geometry
        call updateVMECGeometry(whichMode, whichLambda, deltaLambda)

        ! Compute solutionVec and diagnostics with new geometry
        call mainSolverLoop()

        ! Compute finite difference derivatives
        dParticleFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (particleFlux_vm_rN(whichSpecies)-particleFluxInit(whichSpecies))/deltaLambda
        dHeatFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (heatFlux_vm_rN(whichSpecies)-heatFluxInit(whichSpecies))/deltaLambda
        dParallelFlowdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (FSABVelocityUsingFSADensityOverRootFSAB2(whichSpecies)-parallelFlowInit(whichSpecies))/deltaLambda

        ! Reset geometry to original values
        call updateVMECGeometry(whichMode, whichLambda, -deltaLambda)
      end do
    end do
  end do
  call PetscTime(time1, ierr)
  if (masterProc) then
    print *,"Time for finite difference: ", time1-startTime
  end if

  percentError = zero

  do whichMode = 1, 1
    do whichLambda = 1, 1
      do ispecies = 1, 1
        if (masterProc) then
          print "(a,i4,a,i4,a,i4,a)","Benchmarking fluxes for ispecies: ", ispecies," whichLambda: ", whichLambda," whichMode: ",whichMode," -----------------------------"
        do whichQuantity = 1,3
          select case(whichQuantity)
            case(1)
            analyticResult = dParticleFluxdLambda_analytic(ispecies,whichLambda,whichMode)
            finiteDiffResult = dParticleFluxdLambda_finitediff(ispecies,whichLambda,whichMode)
            if (masterProc) then
              print *,"particle flux"
            end if
            case(2)
            analyticResult = dHeatFluxdLambda_analytic(ispecies,whichLambda,whichMode)
            finiteDiffResult = dHeatFluxdLambda_finitediff(ispecies,whichLambda,whichMode)
            if (masterProc) then
              print *,"heat flux"
            end if
            case(3)
            analyticResult = dParallelFlowdLambda_analytic(ispecies,whichLambda,whichMode)
            finiteDiffResult = dParallelFlowdLambda_finitediff(ispecies,whichLambda,whichMode)
            if (masterProc) then
              print *,"parallel flow"
            end if
          end select
          if (abs(finiteDiffResult) > 1.d-12) then
            percentError = 100.*abs(finiteDiffResult - analyticResult)/abs(finiteDiffResult)
          else if (abs(analyticResult) < 1.d-12) then
            percentError = zero
          else
            percentError = 1.d-6
          end if
          if (masterProc) then
            print *,"percent error: ", percentError
          end if
          end do
        end if
      end do
    end do
  end do
!
!  do whichSpecies = 1,NSpecies
!    call parallelFlowSensitivity(dparallelFlowdLambda_analytic, forwardSolution, whichSpecies, whichLambda, whichMode)
!    call heatFluxSensitivity(dHeatFluxdLambda_analytic, forwardSolution, whichSpecies, whichLambda, whichMode)
!    call particleFluxSensitivity(dParticleFluxdLambda_analytic, forwardSolution, whichSpecies, whichLambda, whichMode)
!    if (masterProc) then
!      print "(a,i4,a)","Benchmarking fluxes for ispecies: ", whichSpecies," -----------------------------"
!    end if
!    finiteDiffDerivative = (particleFlux_vm_rN(whichSpecies)-particleFluxInit(whichSpecies))/deltaLambda
!    if (abs(finiteDiffDerivative) < 1e-16 .and. abs(dparticleFluxdLambda_analytic) < 1e-16) then
!      percentError = zero
!    else if (abs(finiteDiffDerivative) > 1e-16) then
!      percentError = 100*abs(dParticleFluxdLambda_analytic-finiteDiffDerivative)/abs(finiteDiffDerivative)
!    else
!      percentError = 1e6
!    end if
!    if (percentError > 1.0) then
!      if (masterProc) then
!        print "(a,es14.7,a)","percent error: ", percentError,"%"
!        print "(a,es14.7)","dparticleFluxdLambda (finite diff): ", finiteDiffDerivative
!        print "(a,es14.7)","dparticleFluxdLambda (analytic): ", dParticleFluxdLambda_analytic
!      end if
!    end if
!
!    finiteDiffDerivative = (heatFlux_vm_rN(whichSpecies)-heatFluxInit(whichSpecies))/deltaLambda
!    if (abs(finiteDiffDerivative) < 1e-16 .and. abs(dHeatFluxdLambda_analytic) < 1e-16) then
!      percentError = zero
!    else if (abs(finiteDiffDerivative) > 1e-16) then
!      percentError = 100*abs(dHeatFluxdLambda_analytic-finiteDiffDerivative)/abs(finiteDiffDerivative)
!    else
!      percentError = 1e6
!    end if
!
!    if (percentError > 1.0 .and. masterProc) then
!      print "(a,es14.7,a)","percent error: ", percentError,"%"
!      print "(a,es14.7)","dheatFluxdLambda (finite diff): ", finiteDiffDerivative
!      print "(a,es14.7)","dheatFluxdLambda (analytic): ", dHeatFluxdLambda_analytic
!    end if
!
!    finiteDiffDerivative = (FSABVelocityUsingFSADensityOverRootFSAB2(whichSpecies)-parallelFlowInit(whichSpecies))/deltaLambda
!    if (abs(finiteDiffDerivative) < 1e-16 .and. abs(dparallelFlowdLambda_analytic) < 1e-16) then
!      percentError = zero
!    else if (abs(finiteDiffDerivative) > 1e-16) then
!      percentError = 100*abs(dparallelFlowdLambda_analytic-finiteDiffDerivative)/abs(finiteDiffDerivative)
!    else
!      percentError = 1e6
!    end if
!
!    if (percentError > 1.0 .and. masterProc) then
!      print "(a,es14.7,a)","percent error: ", percentError,"%"
!      print "(a,es14.7)","dparallelFlowdLambda (finite diff): ", finiteDiffDerivative
!      print "(a,es14.7)","dparallelFlowdLambda (analytic): ", dparallelFlowdLambda_analytic
!    end if
!  end do
!
!    ! Compute diagnostics with old geometry for testing
!    call diagnostics(forwardSolution, iterationNum)
!
!    do whichSpecies=1,Nspecies
!      print *,"Delta particleFlux (should be 0): ", particleFlux_vm_rN(whichSpecies)-particleFluxInit(whichSpecies)
!      print *,"Delta heatFlux (should be 0): ",heatFlux_vm_rN(whichSpecies)-heatFluxInit(whichSpecies)
!      print *,"Delta parallelFlow (should be 0):",FSABVelocityUsingFSADensityOverRootFSAB2(whichSpecies)-parallelFlowInit(whichSpecies)
!    end do

end subroutine testingAdjointDiagnostics
