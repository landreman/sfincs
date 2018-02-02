#include "../PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

module testingAdjointDiagnostics

  implicit none

  contains

  function percentError(resultAnalytic,resultFiniteDiff)
    implicit none

    PetscScalar :: percentError
    PetscScalar, intent(in) :: resultAnalytic, resultFiniteDiff

    if (abs(resultFiniteDiff) > 1.d-12) then
      percentError = abs(resultFiniteDiff-resultAnalytic)/abs(resultFiniteDiff)
    else if (abs(resultAnalytic) < 1.d-12) then
      percentError = 0
    else
      percentError = 1.d6
    end if

  end function percentError

  !> This subroutine is used to test the explicit derivatives of
  !! the particle flux, heat flux, and flow on geometry with
  !! finite difference derivatives. This has been written to test
  !! the subroutines heatFluxSensitivity, particleFluxSensitivity,
  !! and parallelFlowSensitivity in addition to the adjoint solves 
  !! performed in solver.F90 and innerProduct.
  subroutine compareAdjointDiagnostics()

    use globalVariables
    use petscvec
    use adjointDiagnostics
    use geometry
    use solver

    implicit none

  !  Vec :: forwardSolution
    integer :: whichMode, whichLambda
    integer :: iterationNum
    PetscScalar, dimension(:), allocatable :: particleFluxInit, heatFluxInit, parallelFlowInit
    PetscScalar :: finiteDiffDerivative
    VecScatter :: VecScatterContext
    PetscErrorCode :: ierr
    PetscScalar, dimension(:,:,:), allocatable :: dParticleFluxdLambda_analytic, dHeatFluxdLambda_analytic, dParallelFlowdLambda_analytic
    PetscScalar, dimension(:,:), allocatable :: dRadialCurrentdLambda_analytic, dTotalHeatFluxdLambda_analytic, dBootstrapdLambda_analytic
    integer :: ispecies, whichQuantity
    PetscScalar :: analyticResult, finiteDiffResult
    PetscLogDouble :: startTime, time1

    allocate(particleFluxInit(Nspecies))
    allocate(heatFluxInit(Nspecies))
    allocate(parallelFlowInit(Nspecies))

    allocate(dParticleFluxdLambda_analytic(Nspecies,NLambdas,NModesAdjoint))
    allocate(dHeatFluxdLambda_analytic(Nspecies,NLambdas,NModesAdjoint))
    allocate(dParallelFlowdLambda_analytic(Nspecies,NLambdas,NModesAdjoint))

    allocate(particleFluxPercentError(Nspecies,NLambdas,NModesAdjoint))
    allocate(heatFluxPercentError(Nspecies,NLambdas,NModesAdjoint))
    allocate(parallelFlowPercentError(Nspecies,NLambdas,NModesAdjoint))

    allocate(radialCurrentPercentError(NLambdas,NModesAdjoint))
    allocate(bootstrapPercentError(NLambdas,NModesAdjoint))
    allocate(totalHeatFluxPercentError(NLambdas,NModesAdjoint))

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
    dRadialCurrentdLambda_analytic = dRadialCurrentdLambda
    dTotalHeatFluxdLambda_analytic = dTotalHeatFluxdLambda
    dBootstrapdLambda_analytic = dBootstrapdLambda

    ! Set RHSMode = 1 so call to solver does not include adjoint solve
    RHSMode = 1
    call PetscTime(time1, ierr)
    startTime = time1
    do whichMode = 1, NModesAdjoint
      do whichLambda = 1, NLambdas
        ! Update geometry
        call updateVMECGeometry(whichMode, whichLambda, .false.)

        ! Compute solutionVec and diagnostics with new geometry
        call mainSolverLoop()

        do ispecies = 1, Nspecies

          ! Compute finite difference derivatives
          dParticleFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (particleFlux_vm_rN(iSpecies)-particleFluxInit(iSpecies))/deltaLambda
          dHeatFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (heatFlux_vm_rN(iSpecies)-heatFluxInit(iSpecies))/deltaLambda
          dParallelFlowdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (FSABVelocityUsingFSADensityOverRootFSAB2(iSpecies)-parallelFlowInit(iSpecies))/deltaLambda
        end do
        dTotalHeatFluxdLambda_finiteDiff(whichLambda,whichMode) = (sum(heatFlux_vm_rN)-sum(heatFluxInit))/deltaLambda
        dRadialCurrentdLambda_finiteDiff(whichLambda,whichMode) = (sum(Zs*particleFlux_vm_rN)-sum(Zs*particleFluxInit))/deltaLambda
        dBootstrapdLambda_finiteDiff(whichLambda,whichMode) = (sum(Zs*FSABVelocityUsingFSADensityOverRootFSAB2)-sum(Zs*parallelFlowInit))/deltaLambda
        ! Reset geometry to original values
        call updateVMECGeometry(whichMode, whichLambda, .true.)
      end do
    end do
    call PetscTime(time1, ierr)
    if (masterProc) then
      print *,"Time for finite difference: ", time1-startTime
    end if

    do whichMode = 1, NModesAdjoint
      do whichLambda = 1, NLambdas
        radialCurrentPercentError(whichLambda,whichMode) = percentError(dRadialCurrentdLambda_analytic(whichLambda,whichMode),dRadialCurrentdLambda_finitediff(whichLambda,whichMode))
        totalHeatFluxPercentError(whichLambda,whichMode) = percentError(dTotalHeatFluxdLambda_analytic(whichLambda,whichMode),dTotalHeatFluxdLambda_finitediff(whichLambda,whichMode))
        bootstrapPercentError(whichLambda,whichMode) = percentError(dBootstrapdLambda_analytic(whichLambda, whichMode), dBootstrapdLambda_finitediff(whichLambda,whichMode))

        do ispecies = 1, NSpecies
              particleFluxPercentError(ispecies,whichLambda,whichMode) = percentError(dParticleFluxdLambda_analytic(ispecies,whichLambda,whichMode),dParticleFluxdLambda_finitediff(ispecies,whichLambda,whichMode))

              heatFluxPercentError(ispecies,whichLambda,whichMode) = percentError(dHeatFluxdLambda_analytic(ispecies,whichLambda,whichMode),dHeatFluxdLambda_finitediff(ispecies,whichLambda,whichMode))

              parallelFlowPercentError(ispecies,whichLambda,whichMode) = percentError(dParallelFlowdLambda_analytic(ispecies,whichLambda,whichMode),dParallelFlowdLambda_finitediff(ispecies,whichLambda,whichMode))
        end do ! ispecies
      end do ! whichLambda
    end do ! whichMode

    ! Change RHSMode so adjoint-related quantities are written to output
    RHSMode = 4
    call updateOutputFile(1, .false.)
    call finalizeHDF5()

  end subroutine compareAdjointDiagnostics

end module testingAdjointDiagnostics
