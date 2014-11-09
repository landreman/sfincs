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
    PetscScalar :: particleFluxFactor, momentumFluxFactor, heatFluxFactor, NTVFactor
    PetscScalar, dimension(:), allocatable :: densityIntegralWeights
    PetscScalar, dimension(:), allocatable :: flowIntegralWeights
    PetscScalar, dimension(:), allocatable :: pressureIntegralWeights
    PetscScalar, dimension(:), allocatable :: particleFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: momentumFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: heatFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: NTVIntegralWeights
    PetscScalar :: factor


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
       allocate(particleFluxIntegralWeights(Nx))
       allocate(momentumFluxIntegralWeights(Nx))
       allocate(heatFluxIntegralWeights(Nx))
       allocate(NTVIntegralWeights(Nx))

       allocate(B2(Ntheta))

       densityPerturbation=0
       flow=0
       pressurePerturbation=0
       particleFluxBeforeSurfaceIntegral=0
       momentumFluxBeforeSurfaceIntegral=0
       heatFluxBeforeSurfaceIntegral=0
       NTVBeforeSurfaceIntegral=0

       FSADensityPerturbation=0
       FSABFlow=0
       FSAPressurePerturbation=0
       particleFlux=0
       momentumFlux=0
       heatFlux=0
       NTV=0 
       jHat=0
       Phi1Hat=0
       Phi1HatDenominator = 0

       densityIntegralWeights = x*x
       flowIntegralWeights = x*x*x
       pressureIntegralWeights = x*x*x*x
       particleFluxIntegralWeights = x*x*x*x
       momentumFluxIntegralWeights = x*x*x*x*x
       heatFluxIntegralWeights = x*x*x*x*x*x
       NTVIntegralWeights = x*x*x*x 

       ! Convert the PETSc vector into a normal Fortran array:
       call VecGetArrayF90(solnOnProc0, solnArray, ierr)

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
          particleFluxFactor = - pi*Delta*THat*THat*sqrtTHat/(Zs(ispecies)*VPrimeHat*mHat*sqrtMHat*(GHat+iota*IHat))
          momentumFluxFactor = - pi*Delta*THat*THat*THat/(Zs(ispecies)*VPrimeHat*mHat*(GHat+iota*IHat))
          heatFluxFactor = - pi*Delta*THat*THat*THat*sqrtTHat/(2*Zs(ispecies)*VPrimeHat*mHat*sqrtMHat*(GHat+iota*IHat))
          NTVFactor = 4*pi*THat*THat*sqrtTHat/(mHat*sqrtMHat*VPrimeHat*(GHat+iota*IHat))

          do itheta=1,Ntheta
             do izeta=1,Nzeta

                factor = (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))/(BHat(itheta,izeta) ** 3)

                do ix=1,Nx
                   L = 0
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   densityPerturbation(ispecies,itheta,izeta) = densityPerturbation(ispecies,itheta,izeta) &
                        + densityFactor*xWeights(ix)*densityIntegralWeights(ix)*solnArray(index)

                   pressurePerturbation(ispecies,itheta,izeta) = pressurePerturbation(ispecies,itheta,izeta) &
                        + pressureFactor*xWeights(ix)*pressureIntegralWeights(ix)*solnArray(index)

                   particleFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        + factor * (8/three) * particleFluxFactor &
                        * xWeights(ix)*particleFluxIntegralWeights(ix)*solnArray(index)

                   heatFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        + factor * (8/three) * heatFluxFactor &
                        * xWeights(ix)*heatFluxIntegralWeights(ix)*solnArray(index)

                   L = 1
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   flow(ispecies,itheta,izeta) = flow(ispecies,itheta,izeta) &
                        + flowFactor*xWeights(ix)*flowIntegralWeights(ix)*solnArray(index)

                   momentumFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        + factor * (16d+0/15) * momentumFluxFactor &
                        * xWeights(ix)*momentumFluxIntegralWeights(ix)*solnArray(index)

                   L = 2
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   particleFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        = particleFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        + factor * (four/15) * particleFluxFactor &
                        * xWeights(ix)*particleFluxIntegralWeights(ix)*solnArray(index)

                   heatFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        = heatFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        + factor * (four/15) * heatFluxFactor &
                        * xWeights(ix)*heatFluxIntegralWeights(ix)*solnArray(index)

                   NTVBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        = NTVFactor * NTVKernel(itheta,izeta)&
                        * xWeights(ix)*NTVIntegralWeights(ix)*solnArray(index) 

                   L = 3
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                   momentumFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        = momentumFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                        + factor * (four/35) * momentumFluxFactor &
                        * xWeights(ix)*momentumFluxIntegralWeights(ix)*solnArray(index)

                end do
             end do
          end do

          do izeta=1,Nzeta
             B2 = BHat(:,izeta)*BHat(:,izeta)

             FSADensityPerturbation(ispecies) = FSADensityPerturbation(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, densityPerturbation(ispecies,:,izeta)/B2)

             FSABFlow(ispecies) = FSABFlow(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, flow(ispecies,:,izeta)/BHat(:,izeta))

             FSAPressurePerturbation(ispecies) = FSAPressurePerturbation(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, pressurePerturbation(ispecies,:,izeta)/B2)

             particleFlux(ispecies) = particleFlux(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, particleFluxBeforeSurfaceIntegral(ispecies,:,izeta))

             momentumFlux(ispecies) = momentumFlux(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, momentumFluxBeforeSurfaceIntegral(ispecies,:,izeta))

             heatFlux(ispecies) = heatFlux(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, heatFluxBeforeSurfaceIntegral(ispecies,:,izeta))

             NTV(ispecies) = NTV(ispecies) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, NTVBeforeSurfaceIntegral(ispecies,:,izeta)) 

          end do

          jHat = jHat + Zs(ispecies)*flow(ispecies,:,:)
          Phi1Hat = Phi1Hat + Zs(ispecies)*densityPerturbation(ispecies,:,:)
          Phi1HatDenominator = Phi1HatDenominator + Zs(ispecies)*Zs(ispecies)*nHats(ispecies)/THats(ispecies)
       end do

       Phi1Hat = Phi1Hat / (alpha * Phi1HatDenominator)

       FSADensityPerturbation = FSADensityPerturbation / VPrimeHat
       FSABFlow = FSABFlow / VPrimeHat
       FSAPressurePerturbation = FSAPressurePerturbation / VPrimeHat
       FSABjHat = dot_product(Zs(1:Nspecies), FSABFlow)

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
          print *,"   particleFlux:            ", particleflux(ispecies)
          print *,"   momentumFlux:            ", momentumflux(ispecies)
          print *,"   heatFlux:                ", heatflux(ispecies)
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
       deallocate(particleFluxIntegralWeights)
       deallocate(momentumFluxIntegralWeights)
       deallocate(heatFluxIntegralWeights)
       deallocate(NTVIntegralWeights)

       deallocate(B2)
    end if

    ! updateOutputFile should be called by all procs since it contains MPI_Barrier
    ! (in order to be sure the HDF5 file is safely closed before moving on to the next computation.)
    call updateOutputFile(iterationNum)

  end subroutine diagnostics

