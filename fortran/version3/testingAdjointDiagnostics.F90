#include "../PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

module testingAdjointDiagnostics

  implicit none

  contains

  function percentError(resultAnalytic,resultFiniteDiff,deltaFactor)

    implicit none

    PetscScalar :: percentError
    PetscScalar, intent(in) :: resultAnalytic, resultFiniteDiff, deltaFactor

    if (abs(resultFiniteDiff) > 1.d-12) then
      percentError = 100*abs(resultFiniteDiff-resultAnalytic)/abs(resultFiniteDiff)
    else if (abs(resultAnalytic) < 1.d-12) then
      percentError = 0
    else
      percentError = 1.d6
    end if
    if (abs(deltaFactor) < 1.d-12) then
      percentError = 0
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
    PetscScalar :: analyticResult, finiteDiffResult, deltaFactor
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
    adjointBootstrapOption = .true.
    adjointTotalHeatFluxOption = .true.
    adjointRadialCurrentOption = .true.

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
        select case(whichLambda)
          case(1)
            deltaFactor = bmnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(2)
            deltaFactor = bsupthetamnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(3)
            deltaFactor = bsupzetamnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(4)
            deltaFactor = bsubthetamnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(5)
            deltaFactor = bsubzetamnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(6)
            deltaFactor = gmnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
        end select

        ! Update geometry
        call updateVMECGeometry(whichMode, whichLambda, .false.)

        ! Compute solutionVec and diagnostics with new geometry
        call mainSolverLoop()

        do ispecies = 1, Nspecies
          ! Compute finite difference derivatives
          dParticleFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (particleFlux_vm_rN(iSpecies)-particleFluxInit(iSpecies))/(deltaLambda*deltaFactor)
          dHeatFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (heatFlux_vm_rN(iSpecies)-heatFluxInit(iSpecies))/(deltaLambda*deltaFactor)
          dParallelFlowdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (FSABVelocityUsingFSADensityOverRootFSAB2(iSpecies)-parallelFlowInit(iSpecies))/(deltaLambda*deltaFactor)
        end do ! ispecies
        dTotalHeatFluxdLambda_finiteDiff(whichLambda,whichMode) = (sum(heatFlux_vm_rN)-sum(heatFluxInit))/(deltaLambda*deltaFactor)
        dRadialCurrentdLambda_finiteDiff(whichLambda,whichMode) = (dot_product(Zs(1:Nspecies), particleFlux_vm_rN)-dot_product(Zs(1:Nspecies),particleFluxInit))/(deltaLambda*deltaFactor)
        dBootstrapdLambda_finiteDiff(whichLambda,whichMode) = (dot_product(Zs(1:Nspecies),FSABVelocityUsingFSADensityOverRootFSAB2)-dot_product(Zs(1:Nspecies),parallelFlowInit))/(deltaLambda*deltaFactor)

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
        select case(whichLambda)
          case(1)
            deltaFactor = bmnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(2)
            deltaFactor = bsupthetamnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(3)
            deltaFactor = bsupzetamnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(4)
            deltaFactor = bsubthetamnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(5)
            deltaFactor = bsubzetamnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
          case(6)
            deltaFactor = gmnc(ns(whichMode)+vmec%ntor+1,ms(whichMode)+1)
        end select
        radialCurrentPercentError(whichLambda,whichMode) = percentError(dRadialCurrentdLambda_analytic(whichLambda,whichMode),dRadialCurrentdLambda_finitediff(whichLambda,whichMode),deltaFactor)
        totalHeatFluxPercentError(whichLambda,whichMode) = percentError(dTotalHeatFluxdLambda_analytic(whichLambda,whichMode),dTotalHeatFluxdLambda_finitediff(whichLambda,whichMode),deltaFactor)
        bootstrapPercentError(whichLambda,whichMode) = percentError(dBootstrapdLambda_analytic(whichLambda, whichMode), dBootstrapdLambda_finitediff(whichLambda,whichMode),deltaFactor)

        do ispecies = 1, NSpecies
              particleFluxPercentError(ispecies,whichLambda,whichMode) = percentError(dParticleFluxdLambda_analytic(ispecies,whichLambda,whichMode),dParticleFluxdLambda_finitediff(ispecies,whichLambda,whichMode),deltaFactor)

              heatFluxPercentError(ispecies,whichLambda,whichMode) = percentError(dHeatFluxdLambda_analytic(ispecies,whichLambda,whichMode),dHeatFluxdLambda_finitediff(ispecies,whichLambda,whichMode),deltaFactor)

              parallelFlowPercentError(ispecies,whichLambda,whichMode) = percentError(dParallelFlowdLambda_analytic(ispecies,whichLambda,whichMode),dParallelFlowdLambda_finitediff(ispecies,whichLambda,whichMode),deltaFactor)
        end do ! ispecies
      end do ! whichLambda
    end do ! whichMode

    ! Change RHSMode so adjoint-related quantities are written to output
    RHSMode = 4
    call updateOutputFile(1, .false.)
    call finalizeHDF5()

  end subroutine compareAdjointDiagnostics

end module testingAdjointDiagnostics
