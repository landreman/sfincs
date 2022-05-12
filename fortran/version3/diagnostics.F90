! This next subroutine is called as a "Monitor" of SNES, set in solver.F90 using SNESSetMonitor.
  subroutine diagnosticsMonitor(mysnes, iterationNum, residual, userContext, ierr)

#include "PETScVersions.F90"

    use globalVariables, only: masterProc, iterationForMatrixOutput
    
    implicit none

    SNES :: mysnes
    PetscInt :: iterationNum
    PetscReal :: residual
    integer :: userContext(*)
    PetscErrorCode :: ierr
    KSP :: myKSP

    Vec :: soln

    if (masterProc) then
       if (iterationNum > 0) then
          print "(a,i4,a)",    "--------- Completed iteration ",iterationNum," of SNES -----------------------------------"
       end if
       print "(a,es14.7,a)","--------- Residual function norm: ",residual," -----------------------------"
    end if

    if (iterationNum==0) then
       return
    end if


    call SNESGetKSP(mysnes, myKSP, ierr)
    call checkIfKSPConverged(myKSP)

    iterationForMatrixOutput = iterationNum
    call SNESGetSolution(mysnes, soln, ierr)
    call diagnostics(soln, iterationNum)

  end subroutine diagnosticsMonitor

  ! *******************************************************************************************
  ! *******************************************************************************************

  subroutine checkIfKSPConverged(myKSP)

#include "PETScVersions.F90"

    use globalVariables, only: integerToRepresentTrue, integerToRepresentFalse, useIterativeLinearSolver, masterProc
    implicit none

    KSP :: myKSP
    KSPConvergedReason :: reason
    integer :: ierr
    integer :: didLinearCalculationConverge

    if (useIterativeLinearSolver) then
       call KSPGetConvergedReason(myKSP, reason, ierr)
       if (reason>0) then
          if (masterProc) then
             print *,"Linear iteration (KSP) converged.  KSPConvergedReason = ", reason
             select case (reason)
             case (1)
                print *,"  KSP_CONVERGED_RTOL_NORMAL: "
             case (9)
                print *,"  KSP_CONVERGED_ATOL_NORMAL: "
             case (2)
                print *,"  KSP_CONVERGED_RTOL: Norm decreased by rtol."
             case (3)
                print *,"  KSP_CONVERGED_ATOL: Norm is < abstol."
             case (4)
                print *,"  KSP_CONVERGED_ITS: "
             case (5)
                print *,"  KSP_CONVERGED_CG_NEG_CURVE: "
             case (6)
                print *,"  KSP_CONVERGED_CG_CONSTRAINED: "
             case (7)
                print *,"  KSP_CONVERGED_STEP_LENGTH: "
             case (8)
                print *,"  KSP_CONVERGED_HAPPY_BREAKDOWN: "
             end select
          end if
          didLinearCalculationConverge = integerToRepresentTrue
       else
          if (masterProc) then
             print *,"Linear iteration (KSP) did not converge :(   KSPConvergedReason = ", reason
             select case (reason)
             case (-2)
                print *,"  KSP_DIVERGED_NULL: "
             case (-3)
                print *,"  KSP_DIVERGED_ITS: "
             case (-4)
                print *,"  KSP_DIVERGED_DTOL: "
             case (-5)
                print *,"  KSP_DIVERGED_BREAKDOWN: "
             case (-6)
                print *,"  KSP_DIVERGED_BREAKDOWN_BICG: "
             case (-7)
                print *,"  KSP_DIVERGED_NONSYMMETRIC: "
             case (-8)
                print *,"  KSP_DIVERGED_INDEFINITE_PC: "
             case (-9)
                print *,"  KSP_DIVERGED_NAN: "
             case (-10)
                print *,"  KSP_DIVERGED_INDEFINITE_MAT: "
             case (0)
                print *,"  KSP still iterating ?!?!"
             end select
          end if
          didLinearCalculationConverge = integerToRepresentFalse
       end if
    else
       if (masterProc) then
          print *,"Direct linear solver succeeded."
       end if
       didLinearCalculationConverge = integerToRepresentTrue
    end if

  end subroutine checkIfKSPConverged

  ! *******************************************************************************************
  ! *******************************************************************************************

  subroutine extractPhi1(myVec)

#include "PETScVersions.F90"

    use globalVariables, only: Phi1Hat, dPhi1Hatdtheta, dPhi1Hatdzeta, MPIComm, masterProc, ddtheta, ddzeta, Ntheta, Nzeta
    !!use globalVariables, only: includePhi1, zero !!Commented by AM 2018-12
    use globalVariables, only: includePhi1, zero, readExternalPhi1 !!Added by AM 2018-12
    use indices

    implicit none

    Vec :: myVec
    VecScatter :: VecScatterContext
    Vec :: solnOnProc0
    PetscScalar, pointer :: solnArray(:)
    PetscErrorCode :: ierr

    integer :: itheta, izeta, index

    !!if (includePhi1) then !!Commented by AM 2018-12
    if (includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
       if (masterProc) then
          print *,"Computing Phi1"
       end if

       ! Send the entire solution vector to the master process:
       call VecScatterCreateToZero(myVec, VecScatterContext, solnOnProc0, ierr)
       call VecScatterBegin(VecScatterContext, myVec, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecScatterEnd(VecScatterContext, myVec, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
       
       if (masterProc) then
          ! Convert the PETSc vector into a normal Fortran array:
          call VecGetArrayF90(solnOnProc0, solnArray, ierr)
          
          do itheta = 1,Ntheta
             do izeta = 1,Nzeta
                index = getIndex(1,1,1,itheta,izeta,BLOCK_QN)+1
                ! Add 1 because getIndex returns 0-based PETSc indices, not 1-based fortran indices.
                Phi1Hat(itheta,izeta) = solnarray(index)
             end do
          end do
          
          call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)
       end if

       ! Send Phi1Hat from the masterProc to all procs:
       call MPI_Bcast(Phi1Hat, Ntheta*Nzeta, MPI_DOUBLE_PRECISION, 0, MPIComm, ierr)

       dPhi1Hatdtheta = matmul(ddtheta,Phi1Hat)
       dPhi1Hatdzeta = transpose(matmul(ddzeta,transpose(Phi1Hat)))
    else if (includePhi1 .and. readExternalPhi1) then !!Added by AM 2018-12
       continue !!Do nothing with Phi1 !!Added by AM 2018-12
    else
       ! We are not including Phi_1 in the calculation
       Phi1Hat = zero
       dPhi1Hatdtheta = zero
       dPhi1Hatdzeta = zero
    end if

  end subroutine extractPhi1

  ! *******************************************************************************************
  ! *******************************************************************************************

  subroutine diagnostics(solutionWithDeltaF, iterationNum)

#include "PETScVersions.F90"

    use globalVariables
    use indices
    use writeHDF5Output
    use export_f
    use classicalTransport

    implicit none

    PetscErrorCode :: ierr
    PetscInt :: iterationNum

    VecScatter :: VecScatterContext
    Vec :: solutionWithFullF, solutionWithDeltaF
    Vec :: solutionWithDeltaFOnProc0, solutionWithFullFOnProc0, f0OnProc0
    !!!Vec :: expPhi1 !!Added by AM 2016-06
    PetscScalar, pointer :: solutionWithFullFArray(:), solutionWithDeltaFArray(:), f0Array(:)
    !!PetscScalar, pointer :: expPhi1Array(:) !!Added by AM 2016-06

    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, nHat
    PetscScalar, dimension(:), allocatable :: B2
    integer :: i, j, ix, ispecies, itheta, izeta, L, index
    PetscScalar :: densityFactor, flowFactor, pressureFactor
    PetscScalar :: particleFluxFactor_vm, particleFluxFactor_vE
    PetscScalar :: momentumFluxFactor_vm, momentumFluxFactor_vE
    PetscScalar :: heatFluxFactor_vm, heatFluxFactor_vE
    PetscScalar :: NTVFactor
    PetscScalar, dimension(:), allocatable :: densityIntegralWeights
    PetscScalar, dimension(:), allocatable :: flowIntegralWeights
    PetscScalar, dimension(:), allocatable :: pressureIntegralWeights
    PetscScalar, dimension(:), allocatable :: particleFluxIntegralWeights_vm
    PetscScalar, dimension(:), allocatable :: particleFluxIntegralWeights_vE
    PetscScalar, dimension(:), allocatable :: momentumFluxIntegralWeights_vm
    PetscScalar, dimension(:), allocatable :: momentumFluxIntegralWeights_vE
    PetscScalar, dimension(:), allocatable :: heatFluxIntegralWeights_vm
    PetscScalar, dimension(:), allocatable :: heatFluxIntegralWeights_vE
    PetscScalar, dimension(:), allocatable :: NTVIntegralWeights
    PetscScalar :: factor, factor2, factor_vE, temp1, temp2, temp3
    integer :: itheta1, izeta1, ixi1, ix1
    integer :: itheta2, izeta2, ixi2, ix2
    PetscLogDouble :: time1, time2
    PetscViewer :: viewer
    character(len=200) :: filename

    if (masterProc) then
       print *,"Computing diagnostics."
    end if

    ! Find Phi_1 in the PETSc Vec, and store Phi_1 in a standard Fortran 2D array:
    call extractPhi1(solutionWithDeltaF)

    ! Calculate the classical transport:
    call calculateClassicalFlux(.true.,classicalParticleFlux_psiHat,classicalHeatFlux_psiHat)

    ! The solution vector contains the departure from a Maxwellian, not the "full f" distribution function.
    ! Form the full f:
    call VecDuplicate(solutionWithDeltaF, solutionWithFullF, ierr)
    call VecCopy(solutionWithDeltaF, solutionWithFullF, ierr)

    !!!!!!!!!!!!!!!!!!!!!!!
    !!Added by AM 2016-06!!
    !!!!!!!!!!!!!!!!!!!!!!!
    !!!call VecDuplicate(f0, expPhi1, ierr)
    !!!call VecSet(expPhi1, zero, ierr)
    !!!L = 0
    !!!do ispecies = 1,Nspecies
    !!!   do ix = 1,Nx
    !!!      do itheta = ithetaMin,ithetaMax
    !!!         do izeta = izetaMin,izetaMax
    !!!            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
    !!!            call VecSetValue(expPhi1, index, &
    !!!                 exp(-Zs(ispecies)*alpha*Phi1Hat(itheta,izeta)/THats(ispecies)), INSERT_VALUES, ierr)
    !!!         end do
    !!!      end do
    !!!   end do
    !!!end do

    !!!call VecAssemblyBegin(expPhi1, ierr)
    !!!call VecAssemblyEnd(expPhi1, ierr)

    !!!call VecPointwiseMult(f0, f0, expPhi1, ierr)
    

    call init_f0()
    !!!!!!!!!!!!!!!!!!!!!!!

    call VecAXPY(solutionWithFullF, one, f0, ierr)

    ! Create a "scattering context" for sending vectors to the masterProc, and set up Vecs on the masterProc:
    call VecScatterCreateToZero(solutionWithFullF, VecScatterContext, solutionWithFullFOnProc0, ierr)
    call VecDuplicate(solutionWithFullFOnProc0, solutionWithDeltaFOnProc0, ierr)
    call VecDuplicate(solutionWithFullFOnProc0, f0OnProc0, ierr)
    ! Send the vectors to the master process:
    call VecScatterBegin(VecScatterContext, solutionWithFullF, solutionWithFullFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, solutionWithFullF, solutionWithFullFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterBegin(VecScatterContext, solutionWithDeltaF, solutionWithDeltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, solutionWithDeltaF, solutionWithDeltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterBegin(VecScatterContext, f0, f0OnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, f0, f0OnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)

    if (masterProc) then
       ! All computation of moments of the distribution function is then done on the master process:

       allocate(densityIntegralWeights(Nx))
       allocate(flowIntegralWeights(Nx))
       allocate(pressureIntegralWeights(Nx))
       allocate(particleFluxIntegralWeights_vm(Nx))
       allocate(particleFluxIntegralWeights_vE(Nx))
       allocate(momentumFluxIntegralWeights_vm(Nx))
       allocate(momentumFluxIntegralWeights_vE(Nx))
       allocate(heatFluxIntegralWeights_vm(Nx))
       allocate(heatFluxIntegralWeights_vE(Nx))
       allocate(NTVIntegralWeights(Nx))

       allocate(B2(Ntheta))

       densityPerturbation=0
       flow=0
       pressurePerturbation=0
       pressureAnisotropy=0
       particleFluxBeforeSurfaceIntegral_vm0=0
       particleFluxBeforeSurfaceIntegral_vm=0
       particleFluxBeforeSurfaceIntegral_vE0=0
       particleFluxBeforeSurfaceIntegral_vE=0
       momentumFluxBeforeSurfaceIntegral_vm0=0
       momentumFluxBeforeSurfaceIntegral_vm=0
       momentumFluxBeforeSurfaceIntegral_vE0=0
       momentumFluxBeforeSurfaceIntegral_vE=0
       heatFluxBeforeSurfaceIntegral_vm0=0
       heatFluxBeforeSurfaceIntegral_vm=0
       heatFluxBeforeSurfaceIntegral_vE0=0
       heatFluxBeforeSurfaceIntegral_vE=0
       NTVBeforeSurfaceIntegral=0

       FSADensityPerturbation=0
       FSABFlow=0
       FSAPressurePerturbation=0
       particleFlux_vm0_psiHat=0
       particleFlux_vm_psiHat=0
       particleFlux_vE0_psiHat=0
       particleFlux_vE_psiHat=0
       momentumFlux_vm0_psiHat=0
       momentumFlux_vm_psiHat=0
       momentumFlux_vE0_psiHat=0
       momentumFlux_vE_psiHat=0
       heatFlux_vm0_psiHat=0
       heatFlux_vm_psiHat=0
       heatFlux_vE0_psiHat=0
       heatFlux_vE_psiHat=0
       NTV=0 
       jHat=0

       particleFlux_vm_psiHat_vs_x=0
       heatFlux_vm_psiHat_vs_x=0
       FSABFlow_vs_x=0

       densityIntegralWeights = x*x
       flowIntegralWeights = x*x*x
       pressureIntegralWeights = x*x*x*x
       particleFluxIntegralWeights_vm = x*x*x*x
       particleFluxIntegralWeights_vE = x*x
       momentumFluxIntegralWeights_vm = x*x*x*x*x
       momentumFluxIntegralWeights_vE = x*x*x
       heatFluxIntegralWeights_vm = x*x*x*x*x*x
       heatFluxIntegralWeights_vE = x*x*x*x
       NTVIntegralWeights = x*x*x*x 

       ! Convert the PETSc vectors into normal Fortran arrays:
       call VecGetArrayF90(solutionWithFullFOnProc0, solutionWithFullFArray, ierr)
       call VecGetArrayF90(solutionWithDeltaFOnProc0, solutionWithDeltaFArray, ierr)
       call VecGetArrayF90(f0OnProc0, f0Array, ierr)
       !!call VecGetArrayF90(expPhi1, expPhi1Array, ierr) !!Added by AM 2016-06

!!$    if (whichRHS == numRHSs) then
       select case (constraintScheme)
       case (0)
       case (1,3,4)
          do ispecies = 1,Nspecies
             sources(ispecies,1) = solutionWithDeltaFArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)+1)
             sources(ispecies,2) = solutionWithDeltaFArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)+1)
             ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
          end do
       case (2)
          do ispecies = 1,Nspecies
             do ix=1,Nx
                sources(ispecies,ix) = solutionWithDeltaFArray(getIndex(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT)+1)
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
             end do
          end do
       case default
          print *,"Error! Invalid setting for constraintScheme."
          stop
       end select

       !!if (includePhi1) then !!Commented by AM 2018-12
       if (includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
          lambda = solutionWithDeltaFArray(getIndex(1, 1, 1, 1, 1, BLOCK_PHI1_CONSTRAINT)+1)
       else
          lambda = zero
       end if

       do ispecies = 1,Nspecies
          THat = THats(ispecies)
          mHat = mHats(ispecies)
          sqrtTHat = sqrt(THat)
          sqrtMHat = sqrt(mHat)

          densityFactor = 4*pi*THat*sqrtTHat/(mHat*sqrtMHat)
          flowFactor = 4*pi*THat*THat/(three*mHat*mHat)
          pressureFactor = 8*pi*THat*THat*sqrtTHat/(three*mHat*sqrtMHat)
          !particleFluxFactor = - pi*Delta*THat*THat*sqrtTHat/(Zs(ispecies)*VPrimeHat*mHat*sqrtMHat*(GHat+iota*IHat))
          !momentumFluxFactor = - pi*Delta*THat*THat*THat/(Zs(ispecies)*VPrimeHat*mHat*(GHat+iota*IHat))
          !heatFluxFactor = - pi*Delta*THat*THat*THat*sqrtTHat/(2*Zs(ispecies)*VPrimeHat*mHat*sqrtMHat*(GHat+iota*IHat))

          particleFluxFactor_vm = pi*Delta*THat*THat*sqrtTHat/(Zs(ispecies)*VPrimeHat*mHat*sqrtMHat)
          particleFluxFactor_vE = 2*alpha*pi*Delta*THat*sqrtTHat/(VPrimeHat*mHat*sqrtMHat)
          momentumFluxFactor_vm = pi*Delta*THat*THat*THat/(Zs(ispecies)*VPrimeHat*mHat)
          momentumFluxFactor_vE = 2*alpha*pi*Delta*THat*THat/(VPrimeHat*mHat)
          heatFluxFactor_vm = pi*Delta*THat*THat*THat*sqrtTHat/(2*Zs(ispecies)*VPrimeHat*mHat*sqrtMHat)
          heatFluxFactor_vE = 2*alpha*pi*Delta*THat*THat*sqrtTHat/(2*VPrimeHat*mHat*sqrtMHat)

          ! I haven't looked at how the NTV should be computed in the new units.
          ! Here is the way it was done in the multispecies linear version:
          !NTVFactor = 4*pi*THat*THat*sqrtTHat/(mHat*sqrtMHat*VPrimeHat*(GHat+iota*IHat))
          ! Here is my guess at how it should be done in this version:
          NTVFactor = 4*pi*THat*THat*sqrtTHat/(mHat*sqrtMHat*VPrimeHat)

          do itheta=1,Ntheta
             do izeta=1,Nzeta

                !factor = (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))/(BHat(itheta,izeta) ** 3)

                factor = (BHat_sub_theta(itheta,izeta) * dBHatdzeta(itheta,izeta) &
                     - BHat_sub_zeta(itheta,izeta) * dBHatdtheta(itheta,izeta)) / (BHat(itheta,izeta) ** 3)

                if (force0RadialCurrentInEquilibrium) then
                   factor2 = 0
                else
                   ! Note: this next line has not been tested, since I haven't used eqilibria with a radial current.
                   factor2 = 2 * (dBHat_sub_zeta_dtheta(itheta,izeta) - dBHat_sub_theta_dzeta(itheta,izeta)) &
                        / (BHat(itheta,izeta) ** 2)
                end if

                factor_vE = (BHat_sub_theta(itheta,izeta) * dPhi1Hatdzeta(itheta,izeta) &
                     -BHat_sub_zeta(itheta,izeta) * dPhi1Hatdtheta(itheta,izeta)) &
                     / (BHat(itheta,izeta) ** 2)

                do ix=1,Nx
                   L = 0
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   densityPerturbation(ispecies,itheta,izeta) = densityPerturbation(ispecies,itheta,izeta) &
                        + densityFactor*xWeights(ix)*densityIntegralWeights(ix)*solutionWithDeltaFArray(index)

                   pressurePerturbation(ispecies,itheta,izeta) = pressurePerturbation(ispecies,itheta,izeta) &
                        + pressureFactor*xWeights(ix)*pressureIntegralWeights(ix)*solutionWithDeltaFArray(index)

                   particleFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        + (factor * (8/three) + factor2 * (two/three)) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*f0Array(index)

                   particleFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor * (8/three) + factor2 * (two/three)) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)

                   particleFluxBeforeSurfaceIntegral_vE0(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral_vE0(ispecies,itheta,izeta) &
                        + factor_vE * particleFluxFactor_vE &
                        * xWeights(ix)*particleFluxIntegralWeights_vE(ix)*f0Array(index)

                   particleFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        + factor_vE * particleFluxFactor_vE &
                        * xWeights(ix)*particleFluxIntegralWeights_vE(ix)*solutionWithFullFArray(index)

                   !print *,particleFluxBeforeSurfaceIntegral_vE0(ispecies,itheta,izeta),particleFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta)

                   heatFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        + (factor * (8/three) + factor2 * (two/three)) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*f0Array(index)

                   heatFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor * (8/three) + factor2 * (two/three)) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)

                   heatFluxBeforeSurfaceIntegral_vE0(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral_vE0(ispecies,itheta,izeta) &
                        + factor_vE * heatFluxFactor_vE &
                        * xWeights(ix)*heatFluxIntegralWeights_vE(ix)*f0Array(index)

                   heatFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        + factor_vE * heatFluxFactor_vE &
                        * xWeights(ix)*heatFluxIntegralWeights_vE(ix)*solutionWithFullFArray(index)

                   particleFlux_vm_psiHat_vs_x(ispecies,ix) &
                        = particleFlux_vm_psiHat_vs_x(ispecies,ix) &
                        + (factor * (8/three) + factor2 * (two/three)) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index) &
                        * thetaWeights(itheta) * zetaWeights(izeta)

                   heatFlux_vm_psiHat_vs_x(ispecies,ix) &
                        = heatFlux_vm_psiHat_vs_x(ispecies,ix) &
                        + (factor * (8/three) + factor2 * (two/three)) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index) &
                        * thetaWeights(itheta) * zetaWeights(izeta)



                   L = 1
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   flow(ispecies,itheta,izeta) = flow(ispecies,itheta,izeta) &
                        + flowFactor*xWeights(ix)*flowIntegralWeights(ix)*solutionWithDeltaFArray(index)

                   FSABFlow_vs_x(ispecies,ix) = FSABFlow_vs_x(ispecies,ix) &
                        + flowFactor*xWeights(ix)*flowIntegralWeights(ix)*solutionWithDeltaFArray(index) &
                        * thetaWeights(itheta) * zetaWeights(izeta) * BHat(itheta,izeta) / DHat(itheta,izeta)

                   momentumFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        + (factor * (16d+0/15) + factor2 * (two/5)) * momentumFluxFactor_vm * BHat(itheta,izeta) &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*f0Array(index)

                   momentumFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor * (16d+0/15) + factor2 * (two/5)) * momentumFluxFactor_vm * BHat(itheta,izeta) &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)

                   momentumFluxBeforeSurfaceIntegral_vE0(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral_vE0(ispecies,itheta,izeta) &
                        + factor_vE * (two/3) * momentumFluxFactor_vE * BHat(itheta,izeta) &
                        * xWeights(ix)*momentumFluxIntegralWeights_vE(ix)*f0Array(index)

                   momentumFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        + factor_vE * (two/3) * momentumFluxFactor_vE * BHat(itheta,izeta) &
                        * xWeights(ix)*momentumFluxIntegralWeights_vE(ix)*solutionWithFullFArray(index)

                   L = 2
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   particleFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        + (factor+factor2) * (four/15) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*f0Array(index)

                   particleFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor+factor2) * (four/15) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)

                   heatFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        + (factor+factor2) * (four/15) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*f0Array(index)

                   heatFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor+factor2) * (four/15) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)

                   NTVBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        = NTVBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        + NTVFactor * NTVKernel(itheta,izeta)&
                        * xWeights(ix)*NTVIntegralWeights(ix)*solutionWithDeltaFArray(index)

                   pressureAnisotropy(ispecies,itheta,izeta) = pressureAnisotropy(ispecies,itheta,izeta) &
                        + pressureFactor*(-three/5)* &
                        xWeights(ix)*pressureIntegralWeights(ix)*solutionWithDeltaFArray(index)

                   particleFlux_vm_psiHat_vs_x(ispecies,ix) &
                        = particleFlux_vm_psiHat_vs_x(ispecies,ix) &
                        + (factor+factor2) * (four/15) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index) &
                        * thetaWeights(itheta) * zetaWeights(izeta)
                   
                   heatFlux_vm_psiHat_vs_x(ispecies,ix) &
                        = heatFlux_vm_psiHat_vs_x(ispecies,ix) &
                        + (factor+factor2) * (four/15) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index) &
                        * thetaWeights(itheta) * zetaWeights(izeta)


                   L = 3
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   momentumFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral_vm0(ispecies,itheta,izeta) &
                        + (factor+factor2) * (four/35) * momentumFluxFactor_vm * BHat(itheta,izeta) &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*f0Array(index)

                   momentumFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor+factor2) * (four/35) * momentumFluxFactor_vm * BHat(itheta,izeta) &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)

                end do
             end do
          end do

          do izeta=1,Nzeta
             B2 = BHat(:,izeta)*BHat(:,izeta)

             FSADensityPerturbation(ispecies) = FSADensityPerturbation(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, densityPerturbation(ispecies,:,izeta)/DHat(:,izeta))

             FSABFlow(ispecies) = FSABFlow(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, flow(ispecies,:,izeta)*BHat(:,izeta)/DHat(:,izeta))

             FSAPressurePerturbation(ispecies) = FSAPressurePerturbation(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, pressurePerturbation(ispecies,:,izeta)/DHat(:,izeta))

             particleFlux_vm0_psiHat(ispecies) = particleFlux_vm0_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, particleFluxBeforeSurfaceIntegral_vm0(ispecies,:,izeta))

             particleFlux_vm_psiHat(ispecies) = particleFlux_vm_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, particleFluxBeforeSurfaceIntegral_vm(ispecies,:,izeta))

             particleFlux_vE0_psiHat(ispecies) = particleFlux_vE0_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, particleFluxBeforeSurfaceIntegral_vE0(ispecies,:,izeta))

             particleFlux_vE_psiHat(ispecies) = particleFlux_vE_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, particleFluxBeforeSurfaceIntegral_vE(ispecies,:,izeta))

             momentumFlux_vm0_psiHat(ispecies) = momentumFlux_vm0_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, momentumFluxBeforeSurfaceIntegral_vm0(ispecies,:,izeta))

             momentumFlux_vm_psiHat(ispecies) = momentumFlux_vm_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, momentumFluxBeforeSurfaceIntegral_vm(ispecies,:,izeta))

             momentumFlux_vE0_psiHat(ispecies) = momentumFlux_vE0_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, momentumFluxBeforeSurfaceIntegral_vE0(ispecies,:,izeta))

             momentumFlux_vE_psiHat(ispecies) = momentumFlux_vE_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, momentumFluxBeforeSurfaceIntegral_vE(ispecies,:,izeta))

             heatFlux_vm0_psiHat(ispecies) = heatFlux_vm0_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, heatFluxBeforeSurfaceIntegral_vm0(ispecies,:,izeta))

             heatFlux_vm_psiHat(ispecies) = heatFlux_vm_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, heatFluxBeforeSurfaceIntegral_vm(ispecies,:,izeta))

             heatFlux_vE0_psiHat(ispecies) = heatFlux_vE0_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, heatFluxBeforeSurfaceIntegral_vE0(ispecies,:,izeta))

             heatFlux_vE_psiHat(ispecies) = heatFlux_vE_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, heatFluxBeforeSurfaceIntegral_vE(ispecies,:,izeta))

             NTV(ispecies) = NTV(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, NTVBeforeSurfaceIntegral(ispecies,:,izeta)) 

          end do

          jHat = jHat + Zs(ispecies)*flow(ispecies,:,:)

          !!totalDensity(ispecies,:,:) = nHats(ispecies) + densityPerturbation(ispecies,:,:) !!Commented by AM 2016-06
          !!totalPressure(ispecies,:,:) = nHats(ispecies)*THats(ispecies) + pressurePerturbation(ispecies,:,:) !!Commented by AM 2016-06
          totalDensity(ispecies,:,:) = nHats(ispecies)*exp(-Zs(ispecies)*alpha*Phi1Hat(:,:)/THats(ispecies)) + densityPerturbation(ispecies,:,:) !!Added by AM 2016-06
          totalPressure(ispecies,:,:) = nHats(ispecies)*exp(-Zs(ispecies)*alpha*Phi1Hat(:,:)/THats(ispecies))*THats(ispecies) + pressurePerturbation(ispecies,:,:) !!Added by AM 2016-06
          velocityUsingFSADensity(ispecies,:,:) = flow(ispecies,:,:) / nHats(ispecies)
          velocityUsingTotalDensity(ispecies,:,:) = flow(ispecies,:,:) / totalDensity(ispecies,:,:)
          MachUsingFSAThermalSpeed(ispecies,:,:) = velocityUsingFSADensity(ispecies,:,:) * sqrtMHat/sqrtTHat

       end do


       particleFlux_vd_psiHat = particleFlux_vm_psiHat + particleFlux_vE_psiHat
       momentumFlux_vd_psiHat = momentumFlux_vm_psiHat + momentumFlux_vE_psiHat
       heatFlux_vd_psiHat = heatFlux_vm_psiHat + heatFlux_vE_psiHat

       particleFlux_vd1_psiHat = particleFlux_vm_psiHat + particleFlux_vE0_psiHat
       momentumFlux_vd1_psiHat = momentumFlux_vm_psiHat + momentumFlux_vE0_psiHat
       heatFlux_vd1_psiHat = heatFlux_vm_psiHat + heatFlux_vE0_psiHat

       heatFlux_withoutPhi1_psiHat = heatFlux_vm_psiHat + (5/three)*heatFlux_vE0_psiHat

       particleFlux_vm0_psiN = ddpsiN2ddpsiHat * particleFlux_vm0_psiHat
       particleFlux_vm_psiN = ddpsiN2ddpsiHat * particleFlux_vm_psiHat
       particleFlux_vE0_psiN = ddpsiN2ddpsiHat * particleFlux_vE0_psiHat
       particleFlux_vE_psiN = ddpsiN2ddpsiHat * particleFlux_vE_psiHat
       particleFlux_vd1_psiN = ddpsiN2ddpsiHat * particleFlux_vd1_psiHat
       particleFlux_vd_psiN = ddpsiN2ddpsiHat * particleFlux_vd_psiHat
       momentumFlux_vm0_psiN = ddpsiN2ddpsiHat * momentumFlux_vm0_psiHat
       momentumFlux_vm_psiN = ddpsiN2ddpsiHat * momentumFlux_vm_psiHat
       momentumFlux_vE0_psiN = ddpsiN2ddpsiHat * momentumFlux_vE0_psiHat
       momentumFlux_vE_psiN = ddpsiN2ddpsiHat * momentumFlux_vE_psiHat
       momentumFlux_vd1_psiN = ddpsiN2ddpsiHat * momentumFlux_vd1_psiHat
       momentumFlux_vd_psiN = ddpsiN2ddpsiHat * momentumFlux_vd_psiHat
       heatFlux_vm0_psiN = ddpsiN2ddpsiHat * heatFlux_vm0_psiHat
       heatFlux_vm_psiN = ddpsiN2ddpsiHat * heatFlux_vm_psiHat
       heatFlux_vE0_psiN = ddpsiN2ddpsiHat * heatFlux_vE0_psiHat
       heatFlux_vE_psiN = ddpsiN2ddpsiHat * heatFlux_vE_psiHat
       heatFlux_vd1_psiN = ddpsiN2ddpsiHat * heatFlux_vd1_psiHat
       heatFlux_vd_psiN = ddpsiN2ddpsiHat * heatFlux_vd_psiHat
       heatFlux_withoutPhi1_psiN = ddpsiN2ddpsiHat * heatFlux_withoutPhi1_psiHat
       classicalParticleFlux_psiN = ddpsiN2ddpsiHat * classicalParticleFlux_psiHat
       classicalHeatFlux_psiN = ddpsiN2ddpsiHat * classicalHeatFlux_psiHat

       particleFlux_vm0_rHat = ddrHat2ddpsiHat * particleFlux_vm0_psiHat
       particleFlux_vm_rHat = ddrHat2ddpsiHat * particleFlux_vm_psiHat
       particleFlux_vE0_rHat = ddrHat2ddpsiHat * particleFlux_vE0_psiHat
       particleFlux_vE_rHat = ddrHat2ddpsiHat * particleFlux_vE_psiHat
       particleFlux_vd1_rHat = ddrHat2ddpsiHat * particleFlux_vd1_psiHat
       particleFlux_vd_rHat = ddrHat2ddpsiHat * particleFlux_vd_psiHat
       momentumFlux_vm0_rHat = ddrHat2ddpsiHat * momentumFlux_vm0_psiHat
       momentumFlux_vm_rHat = ddrHat2ddpsiHat * momentumFlux_vm_psiHat
       momentumFlux_vE0_rHat = ddrHat2ddpsiHat * momentumFlux_vE0_psiHat
       momentumFlux_vE_rHat = ddrHat2ddpsiHat * momentumFlux_vE_psiHat
       momentumFlux_vd1_rHat = ddrHat2ddpsiHat * momentumFlux_vd1_psiHat
       momentumFlux_vd_rHat = ddrHat2ddpsiHat * momentumFlux_vd_psiHat
       heatFlux_vm0_rHat = ddrHat2ddpsiHat * heatFlux_vm0_psiHat
       heatFlux_vm_rHat = ddrHat2ddpsiHat * heatFlux_vm_psiHat
       heatFlux_vE0_rHat = ddrHat2ddpsiHat * heatFlux_vE0_psiHat
       heatFlux_vE_rHat = ddrHat2ddpsiHat * heatFlux_vE_psiHat
       heatFlux_vd1_rHat = ddrHat2ddpsiHat * heatFlux_vd1_psiHat
       heatFlux_vd_rHat = ddrHat2ddpsiHat * heatFlux_vd_psiHat
       heatFlux_withoutPhi1_rHat = ddrHat2ddpsiHat * heatFlux_withoutPhi1_psiHat
       classicalParticleFlux_rHat = ddrHat2ddpsiHat * classicalParticleFlux_psiHat
       classicalHeatFlux_rHat = ddrHat2ddpsiHat * classicalHeatFlux_psiHat

       particleFlux_vm0_rN = ddrN2ddpsiHat * particleFlux_vm0_psiHat
       particleFlux_vm_rN = ddrN2ddpsiHat * particleFlux_vm_psiHat
       particleFlux_vE0_rN = ddrN2ddpsiHat * particleFlux_vE0_psiHat
       particleFlux_vE_rN = ddrN2ddpsiHat * particleFlux_vE_psiHat
       particleFlux_vd1_rN = ddrN2ddpsiHat * particleFlux_vd1_psiHat
       particleFlux_vd_rN = ddrN2ddpsiHat * particleFlux_vd_psiHat
       momentumFlux_vm0_rN = ddrN2ddpsiHat * momentumFlux_vm0_psiHat
       momentumFlux_vm_rN = ddrN2ddpsiHat * momentumFlux_vm_psiHat
       momentumFlux_vE0_rN = ddrN2ddpsiHat * momentumFlux_vE0_psiHat
       momentumFlux_vE_rN = ddrN2ddpsiHat * momentumFlux_vE_psiHat
       momentumFlux_vd1_rN = ddrN2ddpsiHat * momentumFlux_vd1_psiHat
       momentumFlux_vd_rN = ddrN2ddpsiHat * momentumFlux_vd_psiHat
       heatFlux_vm0_rN = ddrN2ddpsiHat * heatFlux_vm0_psiHat
       heatFlux_vm_rN = ddrN2ddpsiHat * heatFlux_vm_psiHat
       heatFlux_vE0_rN = ddrN2ddpsiHat * heatFlux_vE0_psiHat
       heatFlux_vE_rN = ddrN2ddpsiHat * heatFlux_vE_psiHat
       heatFlux_vd1_rN = ddrN2ddpsiHat * heatFlux_vd1_psiHat
       heatFlux_vd_rN = ddrN2ddpsiHat * heatFlux_vd_psiHat
       heatFlux_withoutPhi1_rN = ddrN2ddpsiHat * heatFlux_withoutPhi1_psiHat
       classicalParticleFlux_rN = ddrN2ddpsiHat * classicalParticleFlux_psiHat
       classicalHeatFlux_rN = ddrN2ddpsiHat * classicalHeatFlux_psiHat

       FSADensityPerturbation = FSADensityPerturbation / VPrimeHat
       FSABFlow = FSABFlow / VPrimeHat
       FSABFlow_vs_x = FSABFlow_vs_x / VPrimeHat
       FSAPressurePerturbation = FSAPressurePerturbation / VPrimeHat
       FSABjHat = dot_product(Zs(1:Nspecies), FSABFlow)
       FSABjHatOverB0 = FSABjHat / B0OverBBar
       FSABjHatOverRootFSAB2 = FSABjHat / sqrt(FSABHat2)

       do ispecies = 1,Nspecies
          FSABVelocityUsingFSADensity(ispecies) = FSABFlow(ispecies) / nHats(ispecies)
       end do
       FSABVelocityUsingFSADensityOverB0 = FSABVelocityUsingFSADensity / B0OverBBar
       FSABVelocityUsingFSADensityOverRootFSAB2 = FSABVelocityUsingFSADensity / sqrt(FSABHat2)

       if (RHSMode==2) then
          ispecies = 1
          nHat = nHats(ispecies)
          THat = THats(ispecies)
          mHat = mHats(ispecies)
          sqrtTHat = sqrt(THat)
          sqrtMHat = sqrt(mHat)

          select case (whichRHS)
          case (1)
             transportMatrix(1,1) = 4/(Delta*Delta)*sqrtTHat/sqrtmHat*Zs(1)*Zs(1)*(GHat+iota*IHat)&
                  *particleFlux_vm_psiHat(1)*B0OverBBar/(THat*THat*GHat*GHat)
             transportMatrix(2,1) = 8/(Delta*Delta)*sqrtTHat/sqrtmHat*Zs(1)*Zs(1)*(GHat+iota*IHat)&
                  *heatFlux_vm_psiHat(1)*B0OverBBar/(THat*THat*THat*GHat*GHat)
             transportMatrix(3,1) = 2*Zs(1)*FSABFlow(1)/(Delta*GHat*THat)
             !transportMatrix(1,1) = 4*(GHat+iota*IHat)*particleFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*sqrtTHat)*GHat)
             !transportMatrix(2,1) = 8*(GHat+iota*IHat)*heatFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*THat*sqrtTHat)*GHat)
             !transportMatrix(3,1) = 2*nHat*FSABFlow/(GHat*THat)
          case (2)
             transportMatrix(1,2) = 4/(Delta*Delta)*sqrtTHat/sqrtmHat*Zs(1)*Zs(1)*(GHat+iota*IHat)&
                  *particleFlux_vm_psiHat(1)*B0OverBBar/(nHat*THat*GHat*GHat)
             transportMatrix(2,2) = 8/(Delta*Delta)*sqrtTHat/sqrtmHat*Zs(1)*Zs(1)*(GHat+iota*IHat)&
                  *heatFlux_vm_psiHat(1)*B0OverBBar/(nHat*THat*THat*GHat*GHat)
             transportMatrix(3,2) = 2*Zs(1)*FSABFlow(1)/(Delta*GHat*nHat)
             !transportMatrix(1,2) = 4*(GHat+iota*IHat)*particleFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*GHat)
             !transportMatrix(2,2) = 8*(GHat+iota*IHat)*heatFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*THat*GHat)
             !transportMatrix(3,2) = 2*FSABFlow/(GHat)
          case (3)
             transportMatrix(1,3) = particleFlux_vm_psiHat(1)*2*FSABHat2/(nHat*alpha*Delta*GHat)
             transportMatrix(2,3) = heatFlux_vm_psiHat(1)*4*FSABHat2/(nHat*THat*alpha*Delta*GHat)
             transportMatrix(3,3) = FSABFlow(1)*sqrtTHat*sqrtMHat*FSABHat2/((GHat+iota*IHat)*alpha*Zs(1)*nHat*B0OverBBar)
             !transportMatrix(1,3) = particleFlux*Delta*Delta*FSABHat2/(VPrimeHatWithG*GHat*psiAHat*omega)
             !transportMatrix(2,3) = 2*Delta*Delta*heatFlux*FSABHat2/(GHat*VPrimeHatWithG*psiAHat*THat*omega)
             !transportMatrix(3,3) = FSABFlow*Delta*Delta*sqrtTHat*FSABHat2/((GHat+iota*IHat)*2*psiAHat*omega*B0OverBBar)
          end select
       end if

       if (RHSMode==3) then
          ! Monoenergetic transport matrix
          ispecies = 1
          nHat = nHats(ispecies)
          THat = THats(ispecies)
          mHat = mHats(ispecies)
          sqrtTHat = sqrt(THat)
          sqrtMHat = sqrt(mHat)

          ! The factors of THat, mHat, nHat, and Z are unnecessary below (since all are 1).
          select case (whichRHS)
          
          case (1)
             transportMatrix(1,1) = 4/(Delta*Delta)*sqrtTHat/sqrtmHat*Zs(1)*Zs(1)*(GHat+iota*IHat)&
                  *particleFlux_vm_psiHat(1)*B0OverBBar/(THat*THat*GHat*GHat)
             transportMatrix(2,1) = 2*Zs(1)*FSABFlow(1)/(Delta*GHat*THat)
          case (2)
             transportMatrix(1,2) = particleFlux_vm_psiHat(1)*2*FSABHat2/(nHat*alpha*Delta*GHat)
             transportMatrix(2,2) = FSABFlow(1)*sqrtTHat*sqrtMHat*FSABHat2/((GHat+iota*IHat)*alpha*Zs(1)*nHat*B0OverBBar)
          end select
       end if


       ! Interpolate the distribution function from the original grids (used for solving the kinetic equation)
       ! onto whichever grids are requested in the export_f namelist.  I do this here by multiplying by a dense
       ! matrix in each of the 4 coordinates (theta, zeta, x, xi).  This is not the fastest way to do what we want,
       ! but it is relatively simple, and the time required (up to a few seconds) is negligible compared to the time
       ! required for solving the kinetic equation.
       call PetscTime(time1, ierr)
       if (export_full_f) then
          full_f = zero
          do ispecies = 1,Nspecies
             do itheta1 = 1,Ntheta
                do izeta1 = 1,Nzeta
                   do ix1 = 1,Nx
                      do ixi1 = 1,Nxi_for_x(ix1)
                         index = getIndex(ispecies, ix1, ixi1, itheta1, izeta1, BLOCK_F)+1
                         temp1 = solutionWithFullFArray(index)
                         do itheta2 = 1,N_export_f_theta
                            temp2 = temp1 * map_theta_to_export_f_theta(itheta2, itheta1)
                            do izeta2 = 1,N_export_f_zeta
                               temp3 = temp2 * map_zeta_to_export_f_zeta(izeta2, izeta1)
                               do ix2 = 1,N_export_f_x
                                  ! I arbitrarily chose to replace the loop over export_f_xi with ":"
                                  ! We could pick any of the 4 coordinates for this.
                                  full_f(ispecies, itheta2, izeta2, :, ix2) = &
                                       full_f(ispecies, itheta2, izeta2, :, ix2) + temp3 &
                                       * map_x_to_export_f_x(ix2, ix1) &
                                       * map_xi_to_export_f_xi(:, ixi1)
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end if

       if (export_delta_f) then
          delta_f = zero
          do ispecies = 1,Nspecies
             do itheta1 = 1,Ntheta
                do izeta1 = 1,Nzeta
                   do ix1 = 1,Nx
                      do ixi1 = 1,Nxi_for_x(ix1)
                         index = getIndex(ispecies, ix1, ixi1, itheta1, izeta1, BLOCK_F)+1
                         temp1 = solutionWithDeltaFArray(index)
                         do itheta2 = 1,N_export_f_theta
                            temp2 = temp1 * map_theta_to_export_f_theta(itheta2, itheta1)
                            do izeta2 = 1,N_export_f_zeta
                               temp3 = temp2 * map_zeta_to_export_f_zeta(izeta2, izeta1)
                               do ix2 = 1,N_export_f_x
                                  ! I arbitrarily chose to replace the loop over export_f_xi with ":"
                                  ! We could pick any of the 4 coordinates for this.
                                  delta_f(ispecies, itheta2, izeta2, :, ix2) = &
                                       delta_f(ispecies, itheta2, izeta2, :, ix2) + temp3 &
                                       * map_x_to_export_f_x(ix2, ix1) &
                                       * map_xi_to_export_f_xi(:, ixi1)
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end if
       call PetscTime(time2, ierr)
       if (export_delta_f .or. export_full_f) then
          print *,"Time for exporting f: ", time2-time1, " seconds."
       end if

       call VecRestoreArrayF90(solutionWithFullFOnProc0, solutionWithFullFArray, ierr)
       call VecRestoreArrayF90(solutionWithDeltaFOnProc0, solutionWithDeltaFArray, ierr)
       call VecRestoreArrayF90(f0OnProc0, f0Array, ierr)
       !!call VecRestoreArrayF90(expPhi1, expPhi1Array, ierr) !!Added by AM 2016-06

    if (debugAdjoint .eqv. .false.) then
       do ispecies=1,Nspecies
          if (Nspecies>1) then
             print *,"Results for species ",ispecies,":"
          end if
          print *,"   FSADensityPerturbation:  ", FSADensityPerturbation(ispecies)
          print *,"   FSABFlow:                ", FSABFlow(ispecies)
          print *,"   max and min Mach #:      ", maxval(MachUsingFSAThermalSpeed(ispecies,:,:)),&
               minval(MachUsingFSAThermalSpeed(ispecies,:,:))
          print *,"   FSAPressurePerturbation: ", FSAPressurePerturbation(ispecies)
          print *,"   NTV:                     ", NTV(ispecies)
          print *,"   particleFlux_vm0_psiHat  ", particleFlux_vm0_psiHat(ispecies)
          print *,"   particleFlux_vm_psiHat   ", particleFlux_vm_psiHat(ispecies)
          print *,"   classicalParticleFlux    ", classicalParticleFlux_psiHat(ispecies)
          print *,"   classicalHeatFlux        ", classicalHeatFlux_psiHat(ispecies)
          if (includePhi1) then
             print *,"   particleFlux_vE0_psiHat  ", particleFlux_vE0_psiHat(ispecies)
             print *,"   particleFlux_vE_psiHat   ", particleFlux_vE_psiHat(ispecies)
             print *,"   particleFlux_vd1_psiHat  ", particleFlux_vd1_psiHat(ispecies)
             print *,"   particleFlux_vd_psiHat   ", particleFlux_vd_psiHat(ispecies)
          end if
          print *,"   momentumFlux_vm0_psiHat  ", momentumFlux_vm0_psiHat(ispecies)
          print *,"   momentumFlux_vm_psiHat   ", momentumFlux_vm_psiHat(ispecies)
          if (includePhi1) then
             print *,"   momentumFlux_vE0_psiHat  ", momentumFlux_vE0_psiHat(ispecies)
             print *,"   momentumFlux_vE_psiHat   ", momentumFlux_vE_psiHat(ispecies)
             print *,"   momentumFlux_vd1_psiHat  ", momentumFlux_vd1_psiHat(ispecies)
             print *,"   momentumFlux_vd_psiHat   ", momentumFlux_vd_psiHat(ispecies)
          end if
          print *,"   heatFlux_vm0_psiHat      ", heatFlux_vm0_psiHat(ispecies)
          print *,"   heatFlux_vm_psiHat       ", heatFlux_vm_psiHat(ispecies)
          if (includePhi1) then
             print *,"   heatFlux_vE0_psiHat      ", heatFlux_vE0_psiHat(ispecies)
             print *,"   heatFlux_vE_psiHat       ", heatFlux_vE_psiHat(ispecies)
             print *,"   heatFlux_vd1_psiHat      ", heatFlux_vd1_psiHat(ispecies)
             print *,"   heatFlux_vd_psiHat       ", heatFlux_vd_psiHat(ispecies)
             print *,"   heatFlux_withoutPhi1_psiHat ", heatFlux_withoutPhi1_psiHat(ispecies)
          end if
          select case (constraintScheme)
          case (0)
             ! Nothing to print.
          case (1,3,4)
             print *,"   particle source          ", sources(ispecies,1)
             print *,"   heat source              ", sources(ispecies,2)
          case (2)
             print *,"   sources: ", sources(ispecies,:)
          end select
       end do
       print *,"FSABjHat (bootstrap current): ", FSABjHat
       !!if (includePhi1) then !!Commented by AM 2018-12
       if (includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
          print *,"lambda: ", lambda
       end if

       if (rhsMode > 1 .and. RHSMode<4) then
          print *,"Transport matrix:"
          do i=1,transportMatrixSize
             print *,"   ", transportMatrix(i,:)
          end do
       end if
    end if


       deallocate(densityIntegralWeights)
       deallocate(flowIntegralWeights)
       deallocate(pressureIntegralWeights)
       deallocate(particleFluxIntegralWeights_vm)
       deallocate(particleFluxIntegralWeights_vE)
       deallocate(momentumFluxIntegralWeights_vm)
       deallocate(momentumFluxIntegralWeights_vE)
       deallocate(heatFluxIntegralWeights_vm)
       deallocate(heatFluxIntegralWeights_vE)
       deallocate(NTVIntegralWeights)

       deallocate(B2)

    end if

    call VecDestroy(solutionWithFullF, ierr)
    call VecDestroy(solutionWithFullFOnProc0, ierr)
    call VecDestroy(solutionWithDeltaFOnProc0, ierr)
    call VecDestroy(f0OnProc0, ierr)
    !!!call VecDestroy(expPhi1, ierr) !!Added by AM 2016-06

    ! updateOutputFile should be called by all procs since it contains MPI_Barrier
    ! (in order to be sure the HDF5 file is safely closed before moving on to the next computation.)
    ! if (RHSMode >1 .and. whichRHS==transportMatrixSize .and. RHSMode<4) then !This caused hdf5 outputs not to be initiated
    !   call updateOutputFile(iterationNum, .true.)                            !because it is done only at iterationNum=1
    if (RHSMode >1 .and. RHSMode<4) then                                    !HS new version 20220512
       call updateOutputFile(iterationNum, (whichRHS==transportMatrixSize)) !HS new version 20220512
    else if (RHSMode==1 .and. (ambipolarSolve .eqv. .false.) .and. (debugAdjoint .eqv. .false.)) then
       call updateOutputFile(iterationNum, .false.)
    end if

    if (saveMatlabOutput) then
       write (filename,fmt="(a,i3.3,a)") trim(MatlabOutputFilename) // "_iteration_", iterationForStateVectorOutput, &
            "_stateVector.m"
       if (masterProc) then
          print *,"Saving state vector in matlab format: ",trim(filename)
       end if
       call PetscViewerASCIIOpen(MPIComm, trim(filename), viewer, ierr)
       call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB, ierr)

       call PetscObjectSetName(solutionWithDeltaF, "stateVector", ierr)
       call VecView(solutionWithDeltaF, viewer, ierr)

       call PetscViewerDestroy(viewer, ierr)
    end if

    if (saveMatricesAndVectorsInBinary) then
       write (filename,fmt="(a,i3.3,a)") trim(binaryOutputFilename) // "_iteration_", iterationForStateVectorOutput, &
            "_stateVector"
       if (masterProc) then
          print *,"Saving state vector in binary format: ",trim(filename)
       end if
       call PetscViewerBinaryOpen(MPIComm, trim(filename), FILE_MODE_WRITE, viewer, ierr)
       call VecView(solutionWithDeltaF, viewer, ierr)
       call PetscViewerDestroy(viewer, ierr)
    end if

    iterationForStateVectorOutput = iterationForStateVectorOutput + 1

  end subroutine diagnostics

