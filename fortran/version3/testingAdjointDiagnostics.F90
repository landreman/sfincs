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
  integer :: iterationNum
  PetscScalar, dimension(:), allocatable :: particleFluxInit, heatFluxInit, parallelFlowInit
  PetscScalar :: finiteDiffDerivative
  VecScatter :: VecScatterContext
  PetscErrorCode :: ierr
  PetscScalar :: percentError
  PetscScalar, dimension(:,:,:), allocatable :: dParticleFluxdLambda_analytic, dHeatFluxdLambda_analytic, dParallelFlowdLambda_analytic
  integer :: ispecies, whichQuantity
  PetscScalar :: analyticResult, finiteDiffResult
  PetscLogDouble :: startTime, time1

  allocate(particleFluxInit(Nspecies))
  allocate(heatFluxInit(Nspecies))
  allocate(parallelFlowInit(Nspecies))

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
      ! Update geometry
      call updateVMECGeometry(whichMode, whichLambda, deltaLambda)

      do ispecies = 1, 2

        ! Compute solutionVec and diagnostics with new geometry
        call mainSolverLoop()

        ! Compute finite difference derivatives
        dParticleFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (particleFlux_vm_rN(iSpecies)-particleFluxInit(iSpecies))/deltaLambda
        dHeatFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (heatFlux_vm_rN(iSpecies)-heatFluxInit(iSpecies))/deltaLambda
        dParallelFlowdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (FSABVelocityUsingFSADensityOverRootFSAB2(iSpecies)-parallelFlowInit(iSpecies))/deltaLambda

      end do
      ! Reset geometry to original values
      call updateVMECGeometry(whichMode, whichLambda, -deltaLambda)
    end do
  end do
  call PetscTime(time1, ierr)
  if (masterProc) then
    print *,"Time for finite difference: ", time1-startTime
  end if

  percentError = zero

  do whichMode = 1, 1
    do whichLambda = 1, 1
      do ispecies = 1, 2
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
            percentError = 1.d6
          end if
          if (masterProc) then
            print *,"percent error: ", percentError
          end if
          end do
        end if
      end do
    end do
  end do

  ! Change RHSMode so adjoint-related quantities are written to output
  RHSMode = 4
  call updateOutputFile(1, .false.)
  call finalizeHDF5()

end subroutine testingAdjointDiagnostics
