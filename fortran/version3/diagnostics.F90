#include <finclude/petscsnesdef.h>
#include "PETScVersions.F90"

! This next subroutine is called as a "Monitor" of SNES, set in solver.F90 using SNESSetMonitor.
  subroutine diagnosticsMonitor(mysnes, iterationNum, residual, userContext, ierr)

    use globalVariables, only: masterProc
    use petscsnes

    implicit none

    SNES :: mysnes
    PetscInt :: iterationNum
    PetscReal :: residual
    integer :: userContext(*)
    PetscErrorCode :: ierr

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

    call SNESGetSolution(mysnes, soln, ierr)

    call diagnostics(soln, iterationNum)

  end subroutine diagnosticsMonitor

  ! *******************************************************************************************
  ! *******************************************************************************************

  subroutine extractPhi1(myVec)

    use globalVariables, only: Phi1Hat, dPhi1Hatdtheta, dPhi1Hatdzeta, MPIComm, masterProc, ddtheta, ddzeta, Ntheta, Nzeta
    use globalVariables, only: includePhi1, zero
    use indices
    use petscvec

    implicit none

    Vec :: myVec
    VecScatter :: VecScatterContext
    Vec :: solnOnProc0
    PetscScalar, pointer :: solnArray(:)
    PetscErrorCode :: ierr

    integer :: itheta, izeta, index

    if (includePhi1) then
       ! Send the entire solution vector to the master process:
       call VecScatterCreateToZero(myVec, VecScatterContext, solnOnProc0, ierr)
       call VecScatterBegin(VecScatterContext, myVec, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecScatterEnd(VecScatterContext, myVec, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
       
       if (masterProc) then
          ! Convert the PETSc vector into a normal Fortran array:
          call VecGetArrayF90(solnOnProc0, solnArray, ierr)
          
          do itheta = 1,Ntheta
             do izeta = 1,Nzeta
                index = getIndex(1,1,1,itheta,izeta,BLOCK_QN)
                Phi1Hat(itheta,izeta) = solnarray(index)
             end do
          end do
          
          call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)
       end if

! I need to fix this next line!!!
!    call MPI_Bcast(Phi1Hat, Ntheta*Nzeta, 0, MPI_DOUBLE_PRECISION, 0, MPIComm, ierr)

!!$    do itheta = 1,Ntheta
!!$       dPhi1Hatdzeta(itheta,:) = matmul(ddzeta, Phi1Hat(itheta,:))
!!$    end do
!!$
!!$    do izeta = 1,Nzeta
!!$       dPhi1Hatdtheta(:,izeta) = matmul(ddtheta, Phi1Hat(:,izeta))
!!$    end do

       dPhi1Hatdtheta = matmul(ddtheta,Phi1Hat)
       dPhi1Hatdzeta = transpose(matmul(ddzeta,transpose(Phi1Hat)))

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

    use globalVariables
    use indices
    use writeHDF5Output
    use petscvec

    implicit none

    PetscErrorCode :: ierr
    PetscInt :: iterationNum

    VecScatter :: VecScatterContext
    Vec :: solutionWithFullF, solutionWithDeltaF
    Vec :: solutionWithDeltaFOnProc0, solutionWithFullFOnProc0, f0OnProc0
    PetscScalar, pointer :: solutionWithFullFArray(:), solutionWithDeltaFArray(:), f0Array(:)

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
    PetscScalar :: factor, factor2, factor_vE


    if (masterProc) then
       print *,"Computing diagnostics."
    end if

    ! Find Phi_1 in the PETSc Vec, and store Phi_1 in a standard Fortran 2D array:
    call extractPhi1(solutionWithDeltaF)

    ! The solution vector contains the departure from a Maxwellian, not the "full f" distribution function.
    ! Form the full f:
    call VecDuplicate(solutionWithDeltaF, solutionWithFullF, ierr)
    call VecCopy(solutionWithDeltaF, solutionWithFullF, ierr)
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
       Phi1Hat=0

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

!!$    if (whichRHS == numRHSs) then
       select case (constraintScheme)
       case (0)
       case (1)
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

       if (includePhi1) then
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

                   L = 1
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   flow(ispecies,itheta,izeta) = flow(ispecies,itheta,izeta) &
                        + flowFactor*xWeights(ix)*flowIntegralWeights(ix)*solutionWithDeltaFArray(index)

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
                        = NTVFactor * NTVKernel(itheta,izeta)&
                        * xWeights(ix)*NTVIntegralWeights(ix)*solutionWithDeltaFArray(index) 

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

!!$       print *," "
!!$       print *,"particleFluxBeforeSurfaceIntegral_vE0:",particleFluxBeforeSurfaceIntegral_vE0
!!$       print *," "
!!$       print *,"particleFluxBeforeSurfaceIntegral_vE: ",particleFluxBeforeSurfaceIntegral_vE
!!$       print *," "
!!$       print *,"heatFluxBeforeSurfaceIntegral_vE0:",heatFluxBeforeSurfaceIntegral_vE0
!!$       print *," "
!!$       print *,"heatFluxBeforeSurfaceIntegral_vE: ",heatFluxBeforeSurfaceIntegral_vE
!!$       print *," "

          jHat = jHat + Zs(ispecies)*flow(ispecies,:,:)

          totalDensity(ispecies,:,:) = nHats(ispecies) + densityPerturbation(ispecies,:,:)
          totalPressure(ispecies,:,:) = nHats(ispecies)*THats(ispecies) + pressurePerturbation(ispecies,:,:)
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

       FSADensityPerturbation = FSADensityPerturbation / VPrimeHat
       FSABFlow = FSABFlow / VPrimeHat
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

          ! There are a bunch of overall factors in the matrix elements below which need to be worked out:
          select case (whichRHS)
          
          case (1)
             transportMatrix(1,1) = particleFlux_vm_psiHat(ispecies) ! * ???
             transportMatrix(2,1) = FSABFlow(ispecies) ! * ???
          case (2)
             transportMatrix(1,2) = particleFlux_vm_psiHat(ispecies) ! * ???
             transportMatrix(2,2) = FSABFlow(ispecies) ! * ???
          end select
       end if

       call VecRestoreArrayF90(solutionWithFullFOnProc0, solutionWithFullFArray, ierr)
       call VecRestoreArrayF90(solutionWithDeltaFOnProc0, solutionWithDeltaFArray, ierr)
       call VecRestoreArrayF90(f0OnProc0, f0Array, ierr)

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
          if (constraintScheme==1) then
             print *,"   particle source          ", sources(ispecies,1)
             print *,"   heat source              ", sources(ispecies,2)
          end if
          if (constraintScheme==2) then
             print *,"   sources: ", sources(ispecies,:)
          end if
       end do
       print *,"FSABjHat (bootstrap current): ", FSABjHat
       if (includePhi1) then
          print *,"lambda: ", lambda
       end if

       if (rhsMode > 1) then
          print *,"Transport matrix:"
          do i=1,transportMatrixSize
             print *,"   ", transportMatrix(i,:)
          end do
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

    ! updateOutputFile should be called by all procs since it contains MPI_Barrier
    ! (in order to be sure the HDF5 file is safely closed before moving on to the next computation.)
    if (RHSMode >1 .and. whichRHS==transportMatrixSize) then
       call updateOutputFile(iterationNum, .true.)
    else
       call updateOutputFile(iterationNum, .false.)
    end if

  end subroutine diagnostics

