module testingAdjointDiagnostics

#include "PETScVersions.F90"

  use globalVariables

  implicit none

  contains

  function percentError(resultAnalytic,resultFiniteDiff)

    implicit none

    PetscScalar :: percentError
    PetscScalar, intent(in) :: resultAnalytic, resultFiniteDiff

    if (abs(resultFiniteDiff) > 1.d-12) then
      percentError = 100*abs(resultFiniteDiff-resultAnalytic)/abs(resultFiniteDiff)
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
    use ambipolarSolver
		use writeHDF5Output

    implicit none

    integer :: whichMode, whichLambda
    integer :: iterationNum
    PetscScalar, dimension(:), allocatable :: particleFluxInit, heatFluxInit, parallelFlowInit
    PetscScalar :: dPhidPsiInit
    PetscScalar :: finiteDiffDerivative
    VecScatter :: VecScatterContext
    PetscErrorCode :: ierr
    PetscScalar, dimension(:,:,:), allocatable :: dParticleFluxdLambda_analytic, dHeatFluxdLambda_analytic, dParallelFlowdLambda_analytic
    PetscScalar, dimension(:,:), allocatable :: dRadialCurrentdLambda_analytic, dTotalHeatFluxdLambda_analytic, dBootstrapdLambda_analytic, dPhidPsidLambda_analytic
    integer :: ispecies, whichQuantity, this_NmodesAdjoint
    PetscScalar :: analyticResult, finiteDiffResult, deltaFactor
    PetscLogDouble :: time1, time2
    logical :: constantJr

    allocate(particleFluxInit(Nspecies))
    allocate(heatFluxInit(Nspecies))
    allocate(parallelFlowInit(Nspecies))

    allocate(dParticleFluxdLambda_analytic(Nspecies,NLambdas,NModesAdjoint))
    allocate(dHeatFluxdLambda_analytic(Nspecies,NLambdas,NModesAdjoint))
    allocate(dParallelFlowdLambda_analytic(Nspecies,NLambdas,NModesAdjoint))

    if (RHSMode==5) then
      allocate(dPhidPsidLambda_analytic(NLambdas,NModesAdjoint))
      allocate(dPhidPsiPercentError(NLambdas,NModesAdjoint))
      dPhidPsiPercentError = zero
    end if

    ! RHSMode will be overwritten later
    if (RHSMode == 5) then
      constantJr = .true.
    else
      constantJr = .false.
    end if


    ! Change settings so adjoint solve occurs - derivatives of all fluxes computed
    adjointParticleFluxOption = .true.
    adjointHeatFluxOption = .true.
    adjointParallelFlowOption = .true.
    adjointBootstrapOption = .true.
    adjointTotalHeatFluxOption = .true.
    if (constantJr .eqv. .false.) then
      adjointRadialCurrentOption = .true.
    end if

    call PetscTime(time1, ierr)
    if (constantJr) then
			ambipolarSolve = .true.
      call mainAmbipolarSolver()
    else
      call mainSolverLoop()
    end if
    call PetscTime(time2, ierr)
    if (masterProc) then
      print *,"Time for adjoint solve: ", time2-time1
    end if

    ! Get initial values of fluxes
    particleFluxInit = particleFlux_vm_rN
    heatFluxInit = heatFlux_vm_rN
    parallelFlowInit = FSABVelocityUsingFSADensityOverRootFSAB2

    ! Get analytic fluxes from adjoint solve
    dParticleFluxdLambda_analytic = dParticleFluxdLambda
    dHeatFluxdLambda_analytic = dHeatFluxdLambda
    dParallelFlowdLambda_analytic = dParallelFlowdLambda
    dTotalHeatFluxdLambda_analytic = dTotalHeatFluxdLambda
    dBootstrapdLambda_analytic = dBootstrapdLambda
    if (constantJr .eqv. .false.) then
      dRadialCurrentdLambda_analytic = dRadialCurrentdLambda
    end if
    if (constantJr) then
      dPhidPsidLambda_analytic = dPhidPsidLambda
      dPhidPsiInit = dPhiHatdPsiHat
    end if

    ! Set RHSMode = 1 so call to solver does not include adjoint solve
    RHSMode = 1
    call PetscTime(time1, ierr)
		do whichLambda = 1, NLambdas
			if (whichLambda .ne. 1) then
				this_NmodesAdjoint = 1
			else
				this_NmodesAdjoint = NModesAdjoint
			end if
    	do whichMode = 1, this_NModesAdjoint
        select case(whichLambda)
          case(1)
            deltaFactor = bmnc_init(whichMode)
          case(2)
						deltaFactor = IHat_init
          case(3)
						deltaFactor = GHat_init
          case(4)
						deltaFactor = iota_init
        end select

        ! Update geometry
		call updateBoozerGeometry(whichMode, whichLambda, .false.)

        ! Compute solutionVec and diagnostics with new geometry
        if (ambipolarSolve) then
		    ! Need to set this to true so adjoint things are called for ambipolar
		  ambipolarSolve = .true.
          call mainAmbipolarSolver()
        else
          call mainSolverLoop()
        end if

        if (constantJr) then
          dPhidPsidLambda_finitediff(whichLambda,whichMode) = (dPhiHatdPsiHat-dPhidPsiInit)/(deltaLambda*deltaFactor)
        end if

        do ispecies = 1, Nspecies
          ! Compute finite difference derivatives
            dParticleFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (particleFlux_vm_rN(iSpecies)-particleFluxInit(iSpecies))/(deltaLambda*deltaFactor)
			dHeatFluxdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (heatFlux_vm_rN(iSpecies)-heatFluxInit(iSpecies))/(deltaLambda*deltaFactor)
			dParallelFlowdLambda_finiteDiff(ispecies,whichLambda,whichMode) = (FSABVelocityUsingFSADensityOverRootFSAB2(iSpecies)-parallelFlowInit(iSpecies))/(deltaLambda*deltaFactor)
        end do ! ispecies
        dTotalHeatFluxdLambda_finiteDiff(whichLambda,whichMode) = (sum(heatFlux_vm_rN(1:Nspecies))-sum(heatFluxInit(1:Nspecies)))/(deltaLambda*deltaFactor)
        if (constantJr .eqv. .false.) then
            dRadialCurrentdLambda_finiteDiff(whichLambda,whichMode) = (dot_product(Zs(1:Nspecies),particleFlux_vm_rN(1:Nspecies))-dot_product(Zs(1:Nspecies),particleFluxInit(1:Nspecies)))/(deltaLambda*deltaFactor)
        end if
        dBootstrapdLambda_finiteDiff(whichLambda,whichMode) = (dot_product(Zs(1:Nspecies),FSABVelocityUsingFSADensityOverRootFSAB2(1:Nspecies))-dot_product(Zs(1:Nspecies),parallelFlowInit(1:Nspecies)))/(deltaLambda*deltaFactor)
        ! Reset geometry to original values
		call updateBoozerGeometry(whichMode, whichLambda, .true.)
      end do ! whichMode
    end do ! whichLambda
    call PetscTime(time2, ierr)
    if (masterProc) then
      print *,"Time for finite difference: ", time2-time1
    end if

	do whichLambda = 1, NLambdas
		if (whichLambda .ne. 1) then
			this_NmodesAdjoint = 1
		else
			this_NmodesAdjoint = NModesAdjoint
		end if
		do whichMode = 1, this_NModesAdjoint

        select case(whichLambda)
          case(1)
            deltaFactor = bmnc_init(whichMode)
          case(2)
						deltaFactor = IHat_init
          case(3)
						deltaFactor = GHat_init
          case(4)
						deltaFactor = iota_init
        end select
        ! If magnitude of fourier mode is nearly zero, don't use for benchmarking tests
        if (abs(deltaFactor) > 1.d-7) then
          if (constantJr .eqv. .false.) then
            radialCurrentPercentError(whichLambda,whichMode) = percentError(dRadialCurrentdLambda_analytic(whichLambda,whichMode),dRadialCurrentdLambda_finitediff(whichLambda,whichMode))
          end if
          if (constantJr) then
            dPhidPsiPercentError(whichLambda,whichMode) = percentError(dPhidPsidLambda_analytic(whichLambda,whichMode),dPhidPsidLambda_finitediff(whichLambda,whichMode))
          end if
          totalHeatFluxPercentError(whichLambda,whichMode) = percentError(dTotalHeatFluxdLambda_analytic(whichLambda,whichMode),dTotalHeatFluxdLambda_finitediff(whichLambda,whichMode))
          bootstrapPercentError(whichLambda,whichMode) = percentError(dBootstrapdLambda_analytic(whichLambda, whichMode), dBootstrapdLambda_finitediff(whichLambda,whichMode))
          do ispecies = 1, NSpecies
            particleFluxPercentError(ispecies,whichLambda,whichMode) = &
							percentError(dParticleFluxdLambda_analytic(ispecies,whichLambda,whichMode), &
							dParticleFluxdLambda_finitediff(ispecies,whichLambda,whichMode))

            heatFluxPercentError(ispecies,whichLambda,whichMode) = &
							percentError(dHeatFluxdLambda_analytic(ispecies,whichLambda,whichMode), &
							dHeatFluxdLambda_finitediff(ispecies,whichLambda,whichMode))

            parallelFlowPercentError(ispecies,whichLambda,whichMode) = &
							percentError(dParallelFlowdLambda_analytic(ispecies,whichLambda,whichMode), &
							dParallelFlowdLambda_finitediff(ispecies,whichLambda,whichMode))
          end do ! ispecies
        end if
      end do ! whichMode
    end do ! whichLambda

    ! Change RHSMode so adjoint-related quantities are written to output
    if (constantJr) then
      RHSMode = 5
    else
      RHSMode = 4
    end if
    call updateOutputFile(1, .false.)

  end subroutine compareAdjointDiagnostics

end module testingAdjointDiagnostics
