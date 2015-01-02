#include <finclude/petscsnesdef.h>
#include "PETScVersions.F90"

! This subroutine is called as a "Monitor" of SNES, set in solver.F90 using SNESSetMonitor.
  subroutine diagnostics(mysnes, iterationNum, residual, userContext, ierr)

    use petscsnes
    use globalVariables
    use indices
    use writeHDF5Output

    implicit none

    SNES :: mysnes
    PetscInt :: iterationNum
    PetscReal :: residual
    integer :: userContext(*)
    PetscErrorCode :: ierr

    VecScatter :: VecScatterContext
    Vec :: soln, solnOnProc0
    PetscScalar, pointer :: solnArray(:)

    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat
    PetscScalar, dimension(:), allocatable :: B2
    integer :: i, j, ix, ispecies, itheta, izeta, L, index
    PetscScalar :: densityFactor, flowFactor, pressureFactor, Phi1HatDenominator
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
       if (iterationNum > 0) then
          print "(a,i4,a)",    "--------- Completed iteration ",iterationNum," of SNES (nonlinear solver) ----------------"
       end if
       print "(a,es14.7,a)","--------- Residual function norm: ",residual," -----------------------------"
    end if

    if (iterationNum==0) then
       return
    end if

    if (masterProc) then
       print *,"Computing diagnostics."
    end if

    call SNESGetSolution(mysnes, soln, ierr)

    ! Send the entire solution vector to the master process:
    call VecScatterCreateToZero(soln, VecScatterContext, solnOnProc0, ierr)
    call VecScatterBegin(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)

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
       particleFluxBeforeSurfaceIntegral_vm=0
       particleFluxBeforeSurfaceIntegral_vE=0
       momentumFluxBeforeSurfaceIntegral_vm=0
       momentumFluxBeforeSurfaceIntegral_vE=0
       heatFluxBeforeSurfaceIntegral_vm=0
       heatFluxBeforeSurfaceIntegral_vE=0
       NTVBeforeSurfaceIntegral=0

       FSADensityPerturbation=0
       FSABFlow=0
       FSAPressurePerturbation=0
       particleFlux_vm_psiHat=0
       particleFlux_vE_psiHat=0
       momentumFlux_vm_psiHat=0
       momentumFlux_vE_psiHat=0
       heatFlux_vm_psiHat=0
       heatFlux_vE_psiHat=0
       NTV=0 
       jHat=0
       Phi1Hat=0
       Phi1HatDenominator = 0

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

       ! Convert the PETSc vector into a normal Fortran array:
       call VecGetArrayF90(solnOnProc0, solnArray, ierr)

       ! --------------------------------------
       ! Around this point, I should extract Phi_1 from the solution array
       ! And store it in the (theta,zeta) array Phi1Hat.
       ! ---------------------------------------

       ! Temporary hack:
       Phi1Hat = 0

       dPhi1Hatdtheta = matmul(ddtheta,Phi1Hat)
       dPhi1Hatdzeta = transpose(matmul(ddzeta,transpose(Phi1Hat)))

!!$    if (whichRHS == numRHSs) then
       select case (constraintScheme)
       case (0)
       case (1)
          do ispecies = 1,Nspecies
             sources(ispecies,1) = solnArray(getIndex(ispecies, 1, 1, 1, 1, 1)+1)
             sources(ispecies,2) = solnArray(getIndex(ispecies, 1, 1, 1, 1, 2)+1)
             ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
          end do
       case (2)
          do ispecies = 1,Nspecies
             do ix=1,Nx
                sources(ispecies,ix) = solnArray(getIndex(ispecies, ix, 1, 1, 1, 3)+1)
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
             end do
          end do
       case default
          print *,"Error! Invalid setting for constraintScheme."
          stop
       end select

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
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   densityPerturbation(ispecies,itheta,izeta) = densityPerturbation(ispecies,itheta,izeta) &
                        + densityFactor*xWeights(ix)*densityIntegralWeights(ix)*solnArray(index)

                   pressurePerturbation(ispecies,itheta,izeta) = pressurePerturbation(ispecies,itheta,izeta) &
                        + pressureFactor*xWeights(ix)*pressureIntegralWeights(ix)*solnArray(index)

                   particleFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor * (8/three) + factor2 * (two/three)) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solnArray(index)

                   particleFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        + factor_vE * particleFluxFactor_vE &
                        * xWeights(ix)*particleFluxIntegralWeights_vE(ix)*solnArray(index)

                   heatFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor * (8/three) + factor2 * (two/three)) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solnArray(index)

                   heatFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        + factor_vE * heatFluxFactor_vE &
                        * xWeights(ix)*heatFluxIntegralWeights_vE(ix)*solnArray(index)

                   L = 1
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   flow(ispecies,itheta,izeta) = flow(ispecies,itheta,izeta) &
                        + flowFactor*xWeights(ix)*flowIntegralWeights(ix)*solnArray(index)

                   momentumFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor * (16d+0/15) + factor2 * (two/5)) * momentumFluxFactor_vm * BHat(itheta,izeta) &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*solnArray(index)

                   momentumFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral_vE(ispecies,itheta,izeta) &
                        + factor_vE * (two/3) * momentumFluxFactor_vE * BHat(itheta,izeta) &
                        * xWeights(ix)*momentumFluxIntegralWeights_vE(ix)*solnArray(index)

                   L = 2
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   particleFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor+factor2) * (four/15) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solnArray(index)

                   heatFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor+factor2) * (four/15) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solnArray(index)

                   NTVBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        = NTVFactor * NTVKernel(itheta,izeta)&
                        * xWeights(ix)*NTVIntegralWeights(ix)*solnArray(index) 

                   L = 3
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   momentumFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral_vm(ispecies,itheta,izeta) &
                        + (factor+factor2) * (four/35) * momentumFluxFactor_vm * BHat(itheta,izeta) &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*solnArray(index)

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

             particleFlux_vm_psiHat(ispecies) = particleFlux_vm_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, particleFluxBeforeSurfaceIntegral_vm(ispecies,:,izeta))

             particleFlux_vE_psiHat(ispecies) = particleFlux_vE_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, particleFluxBeforeSurfaceIntegral_vE(ispecies,:,izeta))

             momentumFlux_vm_psiHat(ispecies) = momentumFlux_vm_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, momentumFluxBeforeSurfaceIntegral_vm(ispecies,:,izeta))

             momentumFlux_vE_psiHat(ispecies) = momentumFlux_vE_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, momentumFluxBeforeSurfaceIntegral_vE(ispecies,:,izeta))

             heatFlux_vm_psiHat(ispecies) = heatFlux_vm_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, heatFluxBeforeSurfaceIntegral_vm(ispecies,:,izeta))

             heatFlux_vE_psiHat(ispecies) = heatFlux_vE_psiHat(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, heatFluxBeforeSurfaceIntegral_vE(ispecies,:,izeta))

             NTV(ispecies) = NTV(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, NTVBeforeSurfaceIntegral(ispecies,:,izeta)) 

          end do

          jHat = jHat + Zs(ispecies)*flow(ispecies,:,:)
          Phi1Hat = Phi1Hat + Zs(ispecies)*densityPerturbation(ispecies,:,:)
          Phi1HatDenominator = Phi1HatDenominator + Zs(ispecies)*Zs(ispecies)*nHats(ispecies)/THats(ispecies)

          totalDensity(ispecies,:,:) = nHats(ispecies) + densityPerturbation(ispecies,:,:)
          totalPressure(ispecies,:,:) = nHats(ispecies)*THats(ispecies) + pressurePerturbation(ispecies,:,:)
          velocityUsingFSADensity(ispecies,:,:) = flow(ispecies,:,:) / nHats(ispecies)
          velocityUsingTotalDensity(ispecies,:,:) = flow(ispecies,:,:) / totalDensity(ispecies,:,:)
          MachUsingFSAThermalSpeed(ispecies,:,:) = velocityUsingFSADensity(ispecies,:,:) * sqrtMHat/sqrtTHat

       end do

       particleFlux_vd_psiHat = particleFlux_vm_psiHat + particleFlux_vE_psiHat
       momentumFlux_vd_psiHat = momentumFlux_vm_psiHat + momentumFlux_vE_psiHat
       heatFlux_vd_psiHat = heatFlux_vm_psiHat + heatFlux_vE_psiHat

       particleFlux_vm_psiN = ddpsiN2ddpsiHat * particleFlux_vm_psiHat
       particleFlux_vE_psiN = ddpsiN2ddpsiHat * particleFlux_vE_psiHat
       particleFlux_vd_psiN = ddpsiN2ddpsiHat * particleFlux_vd_psiHat
       momentumFlux_vm_psiN = ddpsiN2ddpsiHat * momentumFlux_vm_psiHat
       momentumFlux_vE_psiN = ddpsiN2ddpsiHat * momentumFlux_vE_psiHat
       momentumFlux_vd_psiN = ddpsiN2ddpsiHat * momentumFlux_vd_psiHat
       heatFlux_vm_psiN = ddpsiN2ddpsiHat * heatFlux_vm_psiHat
       heatFlux_vE_psiN = ddpsiN2ddpsiHat * heatFlux_vE_psiHat
       heatFlux_vd_psiN = ddpsiN2ddpsiHat * heatFlux_vd_psiHat

       particleFlux_vm_rHat = ddrHat2ddpsiHat * particleFlux_vm_psiHat
       particleFlux_vE_rHat = ddrHat2ddpsiHat * particleFlux_vE_psiHat
       particleFlux_vd_rHat = ddrHat2ddpsiHat * particleFlux_vd_psiHat
       momentumFlux_vm_rHat = ddrHat2ddpsiHat * momentumFlux_vm_psiHat
       momentumFlux_vE_rHat = ddrHat2ddpsiHat * momentumFlux_vE_psiHat
       momentumFlux_vd_rHat = ddrHat2ddpsiHat * momentumFlux_vd_psiHat
       heatFlux_vm_rHat = ddrHat2ddpsiHat * heatFlux_vm_psiHat
       heatFlux_vE_rHat = ddrHat2ddpsiHat * heatFlux_vE_psiHat
       heatFlux_vd_rHat = ddrHat2ddpsiHat * heatFlux_vd_psiHat

       particleFlux_vm_rN = ddrN2ddpsiHat * particleFlux_vm_psiHat
       particleFlux_vE_rN = ddrN2ddpsiHat * particleFlux_vE_psiHat
       particleFlux_vd_rN = ddrN2ddpsiHat * particleFlux_vd_psiHat
       momentumFlux_vm_rN = ddrN2ddpsiHat * momentumFlux_vm_psiHat
       momentumFlux_vE_rN = ddrN2ddpsiHat * momentumFlux_vE_psiHat
       momentumFlux_vd_rN = ddrN2ddpsiHat * momentumFlux_vd_psiHat
       heatFlux_vm_rN = ddrN2ddpsiHat * heatFlux_vm_psiHat
       heatFlux_vE_rN = ddrN2ddpsiHat * heatFlux_vE_psiHat
       heatFlux_vd_rN = ddrN2ddpsiHat * heatFlux_vd_psiHat

       Phi1Hat = Phi1Hat / (alpha * Phi1HatDenominator)

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

!!$          if (RHSMode==2) then
!!$             VPrimeHatWithG = VPrimeHat*(GHat+iota*IHat)
!!$             select case (whichRHS)
!!$             case (1)
!!$                transportMatrix(1,1) = 4*(GHat+iota*IHat)*particleFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*sqrtTHat)*GHat)
!!$                transportMatrix(2,1) = 8*(GHat+iota*IHat)*heatFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*THat*sqrtTHat)*GHat)
!!$                transportMatrix(3,1) = 2*nHat*FSABFlow/(GHat*THat)
!!$             case (2)
!!$                transportMatrix(1,2) = 4*(GHat+iota*IHat)*particleFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*GHat)
!!$                transportMatrix(2,2) = 8*(GHat+iota*IHat)*heatFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*THat*GHat)
!!$                transportMatrix(3,2) = 2*FSABFlow/(GHat)
!!$             case (3)
!!$                transportMatrix(1,3) = particleFlux*Delta*Delta*FSABHat2/(VPrimeHatWithG*GHat*psiAHat*omega)
!!$                transportMatrix(2,3) = 2*Delta*Delta*heatFlux*FSABHat2/(GHat*VPrimeHatWithG*psiAHat*THat*omega)
!!$                transportMatrix(3,3) = FSABFlow*Delta*Delta*sqrtTHat*FSABHat2/((GHat+iota*IHat)*2*psiAHat*omega*B0OverBBar)
!!$             end select
!!$          end if

       call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)



       do ispecies=1,Nspecies
          if (Nspecies>1) then
             print *,"Results for species ",ispecies,":"
          end if
          print *,"   FSADensityPerturbation:  ", FSADensityPerturbation(ispecies)
          print *,"   FSABFlow:                ", FSABFlow(ispecies)
          print *,"   FSAPressurePerturbation: ", FSAPressurePerturbation(ispecies)
          print *,"   NTV:                     ", NTV(ispecies)
          print *,"   particleFlux_vd_psiHat   ", particleFlux_vd_psiHat(ispecies)
          print *,"   momentumFlux_vd_psiHat   ", momentumFlux_vd_psiHat(ispecies)
          print *,"   heatFlux_vd_psiHat       ", heatFlux_vd_psiHat(ispecies)
       end do
       print *,"FSABjHat (bootstrap current): ", FSABjHat
!!$    if (rhsMode == 2) then
!!$       print *,"Transport matrix:"
!!$       do i=1,3
!!$          print *,"   ", transportMatrix(i,:)
!!$       end do
!!$    end if


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

    ! updateOutputFile should be called by all procs since it contains MPI_Barrier
    ! (in order to be sure the HDF5 file is safely closed before moving on to the next computation.)
    call updateOutputFile(iterationNum)

  end subroutine diagnostics

