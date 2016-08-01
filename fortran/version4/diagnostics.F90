#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif

! This next subroutine is called as a "Monitor" of SNES, set in solver.F90 using SNESSetMonitor.
  subroutine diagnosticsMonitor(mysnes, iterationNum, residual, userContext, ierr)

    use globalVariables, only: masterProc, iterationForMatrixOutput
    use petscsnes

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
       print "(a,es14.7,a)","--------- Residual function norm: ",dble(residual)," -----------------------------"
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

    use globalVariables, only: integerToRepresentTrue, integerToRepresentFalse, masterProc
    use petscksp

    implicit none

    KSP :: myKSP
    KSPConvergedReason :: reason
    integer :: ierr
    integer :: didLinearCalculationConverge

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

  end subroutine checkIfKSPConverged

  ! *******************************************************************************************
  ! *******************************************************************************************

  subroutine extractPhi1(myVec)

    use globalVariables, only: Phi1Hat_Fourier, dPhi1Hatdtheta_Fourier, dPhi1Hatdzeta_Fourier
    use globalVariables, only: Phi1Hat_realSpace, dPhi1Hatdtheta_realSpace, dPhi1Hatdzeta_realSpace
    use globalVariables, only: includePhi1, zero, MPIComm, masterProc, NFourier2
    use indices
    use FourierTransformMod, only: inverseFourierTransform
    use petscvec

    implicit none

    Vec :: myVec
    VecScatter :: VecScatterContext
    Vec :: solnOnProc0
    PetscScalar, pointer :: solnArray(:)
    PetscErrorCode :: ierr

    integer :: imn, index

    if (includePhi1) then
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
          
          do imn = 1,NFourier2
             index = getIndex(1,1,1,imn,BLOCK_QN)+1
             ! Add 1 because getIndex returns 0-based PETSc indices, not 1-based fortran indices.
             Phi1Hat_Fourier(imn) = solnarray(index) ! Note cast from PetscScalar -> real(prec) may happen here.
          end do
          
          call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)
       end if

       ! Send Phi1Hat from the masterProc to all procs:
       call MPI_Bcast(Phi1Hat_Fourier, NFourier2, MPI_DOUBLE_PRECISION, 0, MPIComm, ierr)

       call FourierDerivative(NFourier2,Phi1Hat_Fourier,dPhi1Hatdtheta_Fourier,1)
       call FourierDerivative(NFourier2,Phi1Hat_Fourier,dPhi1Hatdzeta_Fourier, 2)
       call inverseFourierTransform(Phi1Hat_Fourier,Phi1Hat_realSpace)
       call inverseFourierTransform(dPhi1Hatdtheta_Fourier,dPhi1Hatdtheta_realSpace)
       call inverseFourierTransform(dPhi1Hatdzeta_Fourier, dPhi1Hatdzeta_realSpace)

    end if

  end subroutine extractPhi1

  ! *******************************************************************************************
  ! *******************************************************************************************

  subroutine diagnostics(solutionWithDeltaF, iterationNum)

    use FourierConvolutionMatrixMod
    use FourierTransformMod
    use globalVariables
    use indices
    use writeHDF5Output
    use petscvec
    !use export_f

    implicit none

    PetscErrorCode :: ierr
    PetscInt :: iterationNum

    VecScatter :: VecScatterContext
    Vec :: solutionWithFullF, solutionWithDeltaF
    Vec :: solutionWithDeltaFOnProc0, solutionWithFullFOnProc0, f0OnProc0
    !!!Vec :: expPhi1 !!Added by AM 2016-06
    PetscScalar, pointer :: solutionWithFullFArray(:), solutionWithDeltaFArray(:), f0Array(:)
    !!PetscScalar, pointer :: expPhi1Array(:) !!Added by AM 2016-06

    real(prec) :: THat, mHat, sqrtTHat, sqrtMHat, nHat
    integer :: i, j, ix, ispecies, imn, imn2, L, index
    real(prec) :: densityFactor, flowFactor, pressureFactor
    real(prec) :: particleFluxFactor_vm, particleFluxFactor_vE
    real(prec) :: momentumFluxFactor_vm, momentumFluxFactor_vE
    real(prec) :: heatFluxFactor_vm, heatFluxFactor_vE
    real(prec) :: NTVFactor
    real(prec), dimension(:), allocatable :: densityIntegralWeights
    real(prec), dimension(:), allocatable :: flowIntegralWeights
    real(prec), dimension(:), allocatable :: pressureIntegralWeights
    real(prec), dimension(:), allocatable :: particleFluxIntegralWeights_vm
    real(prec), dimension(:), allocatable :: particleFluxIntegralWeights_vE
    real(prec), dimension(:), allocatable :: momentumFluxIntegralWeights_vm
    real(prec), dimension(:), allocatable :: momentumFluxIntegralWeights_vE
    real(prec), dimension(:), allocatable :: heatFluxIntegralWeights_vm
    real(prec), dimension(:), allocatable :: heatFluxIntegralWeights_vE
    real(prec), dimension(:), allocatable :: NTVIntegralWeights

    real(prec), dimension(:), allocatable :: FourierVector
    real(prec), dimension(:,:), allocatable :: FourierMatrix_DInverse
    real(prec), dimension(:,:), allocatable :: FourierMatrix_BOverD
    real(prec), dimension(:,:), allocatable :: FourierMatrix_particleHeat_vm
    real(prec), dimension(:,:), allocatable :: FourierMatrix_momentum_vm
    real(prec), dimension(:,:), allocatable :: FourierMatrix_particleHeat_vE
    real(prec), dimension(:,:), allocatable :: FourierMatrix_momentum_vE

    real(prec) :: factor, factor2, factor_vE, temp1, temp2, temp3
    integer :: itheta1, izeta1, ixi1, ix1
    integer :: itheta2, izeta2, ixi2, ix2
    PetscLogDouble :: time1, time2, presentTime
    PetscViewer :: viewer
    character(len=200) :: filename
    integer :: whichMatrix

    if (masterProc) then
       print *,"Computing diagnostics."
    end if

    ! Find Phi_1 in the PETSc Vec, and store Phi_1 in a standard Fortran 2D array:
    call extractPhi1(solutionWithDeltaF)

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

       densityNonadiabaticPerturbation_Fourier=0
       flow_Fourier=0
       pressureNonadiabaticPerturbation_Fourier=0
       pressureAnisotropy_Fourier=0
       particleFluxBeforeSurfaceIntegral_vm0_Fourier=0
       particleFluxBeforeSurfaceIntegral_vm_Fourier=0
       particleFluxBeforeSurfaceIntegral_vE0_Fourier=0
       particleFluxBeforeSurfaceIntegral_vE_Fourier=0
       momentumFluxBeforeSurfaceIntegral_vm0_Fourier=0
       momentumFluxBeforeSurfaceIntegral_vm_Fourier=0
       momentumFluxBeforeSurfaceIntegral_vE0_Fourier=0
       momentumFluxBeforeSurfaceIntegral_vE_Fourier=0
       heatFluxBeforeSurfaceIntegral_vm0_Fourier=0
       heatFluxBeforeSurfaceIntegral_vm_Fourier=0
       heatFluxBeforeSurfaceIntegral_vE0_Fourier=0
       heatFluxBeforeSurfaceIntegral_vE_Fourier=0

       FSADensityNonadiabaticPerturbation=0
       FSABFlow=0
       FSAPressureNonadiabaticPerturbation=0
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
       jHat_realSpace=0
       jHat_Fourier=0

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
             sources(ispecies,1) = solutionWithDeltaFArray(getIndex(ispecies, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)+1)
             sources(ispecies,2) = solutionWithDeltaFArray(getIndex(ispecies, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)+1)
             ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
          end do
       case (2)
          do ispecies = 1,Nspecies
             do ix=1,Nx
                sources(ispecies,ix) = solutionWithDeltaFArray(getIndex(ispecies, ix, 1, 1, BLOCK_F_CONSTRAINT)+1)
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
             end do
          end do
       case default
          print *,"Error! Invalid setting for constraintScheme."
          stop
       end select

       if (includePhi1) then
          lambda = solutionWithDeltaFArray(getIndex(1, 1, 1, 1, BLOCK_PHI1_CONSTRAINT)+1)
       else
          lambda = zero
       end if


       allocate(FourierVector(NFourier2))
       allocate(FourierMatrix_DInverse(NFourier2,NFourier2))
       allocate(FourierMatrix_BOverD(NFourier2,NFourier2))
       allocate(FourierMatrix_particleHeat_vm(NFourier2,NFourier2))
       allocate(FourierMatrix_momentum_vm(NFourier2,NFourier2))
       allocate(FourierMatrix_particleHeat_vE(NFourier2,NFourier2))
       allocate(FourierMatrix_momentum_vE(NFourier2,NFourier2))

       whichMatrix = 1 ! This value means the convolution matrices will not be simplified, as they would be for the preconditioner.
       ! Spatial weight for flux-surface-average perturbation to density and pressure, and for flow:
       call FourierTransform(1/DHat,FourierVector)
       call FourierConvolutionMatrix(FourierVector,FourierMatrix_DInverse,whichMatrix)
       ! Spatial weight for FSABFlow:
       call FourierTransform(BHat/DHat,FourierVector)
       call FourierConvolutionMatrix(FourierVector,FourierMatrix_BOverD,whichMatrix)
       ! Spatial weight for particle and heat fluxes from magnetic drifts:
       call FourierTransform((BHat_sub_theta*dBHatdzeta - BHat_sub_zeta*dBHatdtheta)/(BHat*BHat*BHat),FourierVector)
       call FourierConvolutionMatrix(FourierVector,FourierMatrix_particleHeat_vm,whichMatrix)
       ! Spatial weight for momentum flux from magnetic drifts:
       call FourierTransform((BHat_sub_theta*dBHatdzeta - BHat_sub_zeta*dBHatdtheta)/(BHat*BHat),FourierVector)
       call FourierConvolutionMatrix(FourierVector,FourierMatrix_momentum_vE,whichMatrix)
       ! Spatial weight for particle and heat fluxes from ExB drift:
       call FourierTransform((BHat_sub_theta*dPhi1Hatdzeta_realSpace - BHat_sub_zeta*dPhi1Hatdtheta_realSpace) &
            /(BHat*BHat),FourierVector)
       call FourierConvolutionMatrix(FourierVector,FourierMatrix_particleHeat_vE,whichMatrix)
       ! Spatial weight for momentum flux from ExB drift:
       call FourierTransform((BHat_sub_theta*dPhi1Hatdzeta_realSpace - BHat_sub_zeta*dPhi1Hatdtheta_realSpace) &
            /(BHat),FourierVector)
       call FourierConvolutionMatrix(FourierVector,FourierMatrix_momentum_vE,whichMatrix)

       do ispecies = 1,Nspecies
          THat = THats(ispecies)
          mHat = mHats(ispecies)
          sqrtTHat = sqrt(THat)
          sqrtMHat = sqrt(mHat)

          densityFactor = 4*pi*THat*sqrtTHat/(mHat*sqrtMHat)
          flowFactor = 4*pi*THat*THat/(three*mHat*mHat)
          pressureFactor = 8*pi*THat*THat*sqrtTHat/(three*mHat*sqrtMHat)

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


          ! First do the density, pressure, and flow, since these do not require any spatial convolution:
          do imn = 1,NFourier2
             do ix = 1,Nx
                L = 0
                index = getIndex(ispecies, ix, L+1, imn, BLOCK_F)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                
                densityNonadiabaticPerturbation_Fourier(ispecies,imn) = densityNonadiabaticPerturbation_Fourier(ispecies,imn) &
                     + densityFactor*xWeights(ix)*densityIntegralWeights(ix)*solutionWithDeltaFArray(index)
                
                pressureNonadiabaticPerturbation_Fourier(ispecies,imn) = pressureNonadiabaticPerturbation_Fourier(ispecies,imn) &
                     + pressureFactor*xWeights(ix)*pressureIntegralWeights(ix)*solutionWithDeltaFArray(index)
                
                
                L = 1
                index = getIndex(ispecies, ix, L+1, imn, BLOCK_F)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                
                flow_Fourier(ispecies,imn) = flow_Fourier(ispecies,imn) &
                     + flowFactor*xWeights(ix)*flowIntegralWeights(ix)*solutionWithDeltaFArray(index)
                
                
                L = 2
                index = getIndex(ispecies, ix, L+1, imn, BLOCK_F)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                
                pressureAnisotropy_Fourier(ispecies,imn) = pressureAnisotropy_Fourier(ispecies,imn) &
                     + pressureFactor*(-three/5)* &
                     xWeights(ix)*pressureIntegralWeights(ix)*solutionWithDeltaFArray(index)
             end do
          end do

          ! Next do BFlow and the fluxes, which require spatial convolutions:

          factor2 = 0 ! Could eventually replace factor2 with a FourierMatrix associated with radial current in the equilibrium.

          do imn = 1,NFourier2     ! imn indexes the diagnostic output
             do imn2 = 1,NFourier2 ! imn2 is the index that is contracted when the top row of the convolution matrix is multiplied by the solution vector.
                do ix = 1,Nx
                   L = 0
                   index = getIndex(ispecies, ix, L+1, imn2, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                   
                   particleFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        = particleFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        + (FourierMatrix_particleHeat_vm(1,imn2) * (8/three) + factor2 * (two/three)) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*f0Array(index)
                   
                   particleFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        = particleFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        + (FourierMatrix_particleHeat_vm(1,imn2) * (8/three) + factor2 * (two/three)) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)
                   
                   particleFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,imn) &
                        = particleFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,imn) &
                        + FourierMatrix_particleHeat_vE(1,imn2) * particleFluxFactor_vE &
                        * xWeights(ix)*particleFluxIntegralWeights_vE(ix)*f0Array(index)
                   
                   particleFluxBeforeSurfaceIntegral_vE_Fourier(ispecies,imn) &
                        = particleFluxBeforeSurfaceIntegral_vE_Fourier(ispecies,imn) &
                        + FourierMatrix_particleHeat_vE(1,imn2) * particleFluxFactor_vE &
                        * xWeights(ix)*particleFluxIntegralWeights_vE(ix)*solutionWithFullFArray(index)
                   
                   heatFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        = heatFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        + (FourierMatrix_particleHeat_vm(1,imn2) * (8/three) + factor2 * (two/three)) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*f0Array(index)
                   
                   heatFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        = heatFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        + (FourierMatrix_particleHeat_vm(1,imn2) * (8/three) + factor2 * (two/three)) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)
                   
                   heatFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,imn) &
                        = heatFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,imn) &
                        + FourierMatrix_particleHeat_vE(1,imn2) * heatFluxFactor_vE &
                        * xWeights(ix)*heatFluxIntegralWeights_vE(ix)*f0Array(index)
                   
                   heatFluxBeforeSurfaceIntegral_vE_Fourier(ispecies,imn) &
                        = heatFluxBeforeSurfaceIntegral_vE_Fourier(ispecies,imn) &
                        + FourierMatrix_particleHeat_vE(1,imn2) * heatFluxFactor_vE &
                        * xWeights(ix)*heatFluxIntegralWeights_vE(ix)*solutionWithFullFArray(index)
                   
!!$                   particleFlux_vm_psiHat_vs_x_Fourier(ispecies,ix) &
!!$                        = particleFlux_vm_psiHat_vs_x_Fourier(ispecies,ix) &
!!$                        + (FourierMatrix_particleHeat_vm(1,imn2) * (8/three) + factor2 * (two/three)) * particleFluxFactor_vm &
!!$                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index) &
!!$                        * thetaWeights(itheta) * zetaWeights(izeta)
!!$                   
!!$                   heatFlux_vm_psiHat_vs_x_Fourier(ispecies,ix) &
!!$                        = heatFlux_vm_psiHat_vs_x_Fourier(ispecies,ix) &
!!$                        + (FourierMatrix_particleHeat_vm(1,imn2) * (8/three) + factor2 * (two/three)) * heatFluxFactor_vm &
!!$                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index) &
!!$                        * thetaWeights(itheta) * zetaWeights(izeta)
                
                
                   L = 1
                   index = getIndex(ispecies, ix, L+1, imn2, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                   
!!$                   FSABFlow_vs_x_Fourier(ispecies,ix) = FSABFlow_vs_x_Fourier(ispecies,ix) &
!!$                        + flowFactor*xWeights(ix)*flowIntegralWeights(ix)*solutionWithDeltaFArray(index) &
!!$                        * thetaWeights(itheta) * zetaWeights(izeta) * BHat(itheta,izeta) / DHat(itheta,izeta)
                   
                   momentumFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        = momentumFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        + (FourierMatrix_momentum_vm(1,imn2) * (16d+0/15) + factor2 * (two/5)) * momentumFluxFactor_vm &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*f0Array(index)
                   
                   momentumFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        = momentumFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        + (FourierMatrix_momentum_vm(1,imn2) * (16d+0/15) + factor2 * (two/5)) * momentumFluxFactor_vm &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)
                   
                   momentumFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,imn) &
                        = momentumFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,imn) &
                        + FourierMatrix_momentum_vE(1,imn2) * (two/3) * momentumFluxFactor_vE &
                        * xWeights(ix)*momentumFluxIntegralWeights_vE(ix)*f0Array(index)
                   
                   momentumFluxBeforeSurfaceIntegral_vE_Fourier(ispecies,imn) &
                        = momentumFluxBeforeSurfaceIntegral_vE_Fourier(ispecies,imn) &
                        + FourierMatrix_momentum_vE(1,imn2) * (two/3) * momentumFluxFactor_vE &
                        * xWeights(ix)*momentumFluxIntegralWeights_vE(ix)*solutionWithFullFArray(index)
                   
                   L = 2
                   index = getIndex(ispecies, ix, L+1, imn2, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                   
                   particleFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        = particleFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        + (FourierMatrix_particleHeat_vm(1,imn2)+factor2) * (four/15) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*f0Array(index)
                   
                   particleFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        = particleFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        + (FourierMatrix_particleHeat_vm(1,imn2)+factor2) * (four/15) * particleFluxFactor_vm &
                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)
                   
                   heatFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        = heatFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        + (FourierMatrix_particleHeat_vm(1,imn2)+factor2) * (four/15) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*f0Array(index)
                   
                   heatFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        = heatFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        + (FourierMatrix_particleHeat_vm(1,imn2)+factor2) * (four/15) * heatFluxFactor_vm &
                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)
                   
!!$                NTVBeforeSurfaceIntegral_Fourier(ispecies,imn) &
!!$                     = NTVBeforeSurfaceIntegral_Fourier(ispecies,imn) &
!!$                     + NTVFactor * NTVKernel(itheta,izeta)&
!!$                     * xWeights(ix)*NTVIntegralWeights(ix)*solutionWithDeltaFArray(index)
                   
!!$                   particleFlux_vm_psiHat_vs_x_Fourier(ispecies,ix) &
!!$                        = particleFlux_vm_psiHat_vs_x_Fourier(ispecies,ix) &
!!$                        + (FourierMatrix_particleHeat_vm(1,imn2)+factor2) * (four/15) * particleFluxFactor_vm &
!!$                        * xWeights(ix)*particleFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index) &
!!$                        * thetaWeights(itheta) * zetaWeights(izeta)
!!$                   
!!$                   heatFlux_vm_psiHat_vs_x_Fourier(ispecies,ix) &
!!$                        = heatFlux_vm_psiHat_vs_x_Fourier(ispecies,ix) &
!!$                        + (FourierMatrix_particleHeat_vm(1,imn2)+factor2) * (four/15) * heatFluxFactor_vm &
!!$                        * xWeights(ix)*heatFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index) &
!!$                        * thetaWeights(itheta) * zetaWeights(izeta)
                   
                   
                   L = 3
                   index = getIndex(ispecies, ix, L+1, imn2, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                   
                   momentumFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        = momentumFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,imn) &
                        + (FourierMatrix_momentum_vm(1,imn2)+factor2) * (four/35) * momentumFluxFactor_vm &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*f0Array(index)
                   
                   momentumFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        = momentumFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,imn) &
                        + (FourierMatrix_momentum_vm(1,imn2)+factor2) * (four/35) * momentumFluxFactor_vm &
                        * xWeights(ix)*momentumFluxIntegralWeights_vm(ix)*solutionWithFullFArray(index)
                end do
             end do
          end do

          ! For quantities that vary on a flux surface, convert from Fourier to real space:
          call inverseFourierTransform(densityNonadiabaticPerturbation_Fourier(ispecies,:), densityNonadiabaticPerturbation_realSpace(ispecies,:,:))
          call inverseFourierTransform(flow_Fourier(ispecies,:), flow_realSpace(ispecies,:,:))
          call inverseFourierTransform(pressureNonadiabaticPerturbation_Fourier(ispecies,:), pressureNonadiabaticPerturbation_realSpace(ispecies,:,:))
          call inverseFourierTransform(pressureAnisotropy_Fourier(ispecies,:), pressureAnisotropy_realSpace(ispecies,:,:))
          call inverseFourierTransform(particleFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,:), particleFluxBeforeSurfaceIntegral_vm0_realSpace(ispecies,:,:))
          call inverseFourierTransform(particleFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,:), particleFluxBeforeSurfaceIntegral_vm_realSpace(ispecies,:,:))
          call inverseFourierTransform(particleFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,:), particleFluxBeforeSurfaceIntegral_vE0_realSpace(ispecies,:,:))
          call inverseFourierTransform(particleFluxBeforeSurfaceIntegral_vE_Fourier(ispecies,:), particleFluxBeforeSurfaceIntegral_vE_realSpace(ispecies,:,:))
          call inverseFourierTransform(momentumFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,:), momentumFluxBeforeSurfaceIntegral_vm0_realSpace(ispecies,:,:))
          call inverseFourierTransform(momentumFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,:), momentumFluxBeforeSurfaceIntegral_vm_realSpace(ispecies,:,:))
          call inverseFourierTransform(momentumFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,:), momentumFluxBeforeSurfaceIntegral_vE0_realSpace(ispecies,:,:))
          call inverseFourierTransform(momentumFluxBeforeSurfaceIntegral_vE_Fourier(ispecies,:), momentumFluxBeforeSurfaceIntegral_vE_realSpace(ispecies,:,:))
          call inverseFourierTransform(heatFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,:), heatFluxBeforeSurfaceIntegral_vm0_realSpace(ispecies,:,:))
          call inverseFourierTransform(heatFluxBeforeSurfaceIntegral_vm_Fourier(ispecies,:), heatFluxBeforeSurfaceIntegral_vm_realSpace(ispecies,:,:))
          call inverseFourierTransform(heatFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,:), heatFluxBeforeSurfaceIntegral_vE0_realSpace(ispecies,:,:))
          call inverseFourierTransform(heatFluxBeforeSurfaceIntegral_vE_Fourier(ispecies,:), heatFluxBeforeSurfaceIntegral_vE_realSpace(ispecies,:,:))

          ! Compute surface-averaged quantities:
          ! (The dot_product between the top row of the Fourier matrix and the Fourier vector gives the constant term of the Fourier expansion, which integrates to 4*pi*pi
          ! when you integrate over theta and zeta.)
          FSADensityNonadiabaticPerturbation(ispecies) = dot_product(FourierMatrix_DInverse(1,:),densityNonadiabaticPerturbation_Fourier(ispecies,:))*4*pi*pi/VPrimeHat
          FSABFlow(ispecies) = dot_product(FourierMatrix_BOverD(1,:),flow_Fourier(ispecies,:))*4*pi*pi/VPrimeHat
          FSAPressureNonadiabaticPerturbation(ispecies) = dot_product(FourierMatrix_DInverse(1,:),pressureNonadiabaticPerturbation_Fourier(ispecies,:))*4*pi*pi/VPrimeHat

          particleFlux_vm0_psiHat(ispecies) = 4*pi*pi*particleFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,1)
          particleFlux_vm_psiHat( ispecies) = 4*pi*pi*particleFluxBeforeSurfaceIntegral_vm_Fourier( ispecies,1)
          particleFlux_vE0_psiHat(ispecies) = 4*pi*pi*particleFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,1)
          particleFlux_vE_psiHat( ispecies) = 4*pi*pi*particleFluxBeforeSurfaceIntegral_vE_Fourier( ispecies,1)

          momentumFlux_vm0_psiHat(ispecies) = 4*pi*pi*momentumFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,1)
          momentumFlux_vm_psiHat( ispecies) = 4*pi*pi*momentumFluxBeforeSurfaceIntegral_vm_Fourier( ispecies,1)
          momentumFlux_vE0_psiHat(ispecies) = 4*pi*pi*momentumFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,1)
          momentumFlux_vE_psiHat( ispecies) = 4*pi*pi*momentumFluxBeforeSurfaceIntegral_vE_Fourier( ispecies,1)

          heatFlux_vm0_psiHat(ispecies) = 4*pi*pi*heatFluxBeforeSurfaceIntegral_vm0_Fourier(ispecies,1)
          heatFlux_vm_psiHat( ispecies) = 4*pi*pi*heatFluxBeforeSurfaceIntegral_vm_Fourier( ispecies,1)
          heatFlux_vE0_psiHat(ispecies) = 4*pi*pi*heatFluxBeforeSurfaceIntegral_vE0_Fourier(ispecies,1)
          heatFlux_vE_psiHat( ispecies) = 4*pi*pi*heatFluxBeforeSurfaceIntegral_vE_Fourier( ispecies,1)

          jHat_realSpace = jHat_realSpace + Zs(ispecies)*flow_realSpace(ispecies,:,:)
          jHat_Fourier   = jHat_Fourier   + Zs(ispecies)*flow_Fourier(  ispecies,:)

          !!totalDensity(ispecies,:,:) = nHats(ispecies) + densityPerturbation(ispecies,:,:) !!Commented by AM 2016-06
          !!totalPressure(ispecies,:,:) = nHats(ispecies)*THats(ispecies) + pressurePerturbation(ispecies,:,:) !!Commented by AM 2016-06

          totalDensity_realSpace(ispecies,:,:) = nHats(ispecies)*exp(-Zs(ispecies)*alpha*Phi1Hat_realSpace/THats(ispecies)) &
               + densityNonadiabaticPerturbation_realSpace(ispecies,:,:)
          
          totalPressure_realSpace(ispecies,:,:) = nHats(ispecies)*exp(-Zs(ispecies)*alpha*Phi1Hat_realSpace/THats(ispecies))*THats(ispecies) &
               + pressureNonadiabaticPerturbation_realSpace(ispecies,:,:)

          call FourierTransform( totalDensity_realSpace(ispecies,:,:),  totalDensity_Fourier(ispecies,:))
          call FourierTransform(totalPressure_realSpace(ispecies,:,:), totalPressure_Fourier(ispecies,:))

          velocityUsingFSADensity_realSpace(ispecies,:,:) = flow_realSpace(ispecies,:,:) / nHats(ispecies)
          velocityUsingFSADensity_Fourier(ispecies,:) = flow_Fourier(ispecies,:) / nHats(ispecies)

          MachUsingFSAThermalSpeed_realSpace(ispecies,:,:) = velocityUsingFSADensity_realSpace(ispecies,:,:) * sqrtMHat/sqrtTHat
          MachUsingFSAThermalSpeed_Fourier(ispecies,:) = velocityUsingFSADensity_Fourier(ispecies,:) * sqrtMHat/sqrtTHat

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

       ! export_f has not been updated to reflect the new Fourier discretization.
!!$       ! Interpolate the distribution function from the original grids (used for solving the kinetic equation)
!!$       ! onto whichever grids are requested in the export_f namelist.  I do this here by multiplying by a dense
!!$       ! matrix in each of the 4 coordinates (theta, zeta, x, xi).  This is not the fastest way to do what we want,
!!$       ! but it is relatively simple, and the time required (up to a few seconds) is negligible compared to the time
!!$       ! required for solving the kinetic equation.
!!$       call PetscTime(time1, ierr)
!!$       if (export_full_f) then
!!$          full_f = zero
!!$          do ispecies = 1,Nspecies
!!$             do itheta1 = 1,Ntheta
!!$                do izeta1 = 1,Nzeta
!!$                   do ixi1 = 1,Nxi
!!$                      do ix1 = 1,Nx
!!$                         index = getIndex(ispecies, ix1, ixi1, itheta1, izeta1, BLOCK_F)+1
!!$                         temp1 = solutionWithFullFArray(index)
!!$                         do itheta2 = 1,N_export_f_theta
!!$                            temp2 = temp1 * map_theta_to_export_f_theta(itheta2, itheta1)
!!$                            do izeta2 = 1,N_export_f_zeta
!!$                               temp3 = temp2 * map_zeta_to_export_f_zeta(izeta2, izeta1)
!!$                               do ix2 = 1,N_export_f_x
!!$                                  ! I arbitrarily chose to replace the loop over export_f_xi with ":"
!!$                                  ! We could pick any of the 4 coordinates for this.
!!$                                  full_f(ispecies, itheta2, izeta2, :, ix2) = &
!!$                                       full_f(ispecies, itheta2, izeta2, :, ix2) + temp3 &
!!$                                       * map_x_to_export_f_x(ix2, ix1) &
!!$                                       * map_xi_to_export_f_xi(:, ixi1)
!!$                               end do
!!$                            end do
!!$                         end do
!!$                      end do
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$       end if
!!$
!!$       if (export_delta_f) then
!!$          delta_f = zero
!!$          do ispecies = 1,Nspecies
!!$             do itheta1 = 1,Ntheta
!!$                do izeta1 = 1,Nzeta
!!$                   do ixi1 = 1,Nxi
!!$                      do ix1 = 1,Nx
!!$                         index = getIndex(ispecies, ix1, ixi1, itheta1, izeta1, BLOCK_F)+1
!!$                         temp1 = solutionWithDeltaFArray(index)
!!$                         do itheta2 = 1,N_export_f_theta
!!$                            temp2 = temp1 * map_theta_to_export_f_theta(itheta2, itheta1)
!!$                            do izeta2 = 1,N_export_f_zeta
!!$                               temp3 = temp2 * map_zeta_to_export_f_zeta(izeta2, izeta1)
!!$                               do ix2 = 1,N_export_f_x
!!$                                  ! I arbitrarily chose to replace the loop over export_f_xi with ":"
!!$                                  ! We could pick any of the 4 coordinates for this.
!!$                                  delta_f(ispecies, itheta2, izeta2, :, ix2) = &
!!$                                       delta_f(ispecies, itheta2, izeta2, :, ix2) + temp3 &
!!$                                       * map_x_to_export_f_x(ix2, ix1) &
!!$                                       * map_xi_to_export_f_xi(:, ixi1)
!!$                               end do
!!$                            end do
!!$                         end do
!!$                      end do
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$       end if
!!$       call PetscTime(time2, ierr)
!!$       if (export_delta_f .or. export_full_f) then
!!$          print *,"Time for exporting f: ", time2-time1, " seconds."
!!$       end if

       call VecRestoreArrayF90(solutionWithFullFOnProc0, solutionWithFullFArray, ierr)
       call VecRestoreArrayF90(solutionWithDeltaFOnProc0, solutionWithDeltaFArray, ierr)
       call VecRestoreArrayF90(f0OnProc0, f0Array, ierr)
       !!call VecRestoreArrayF90(expPhi1, expPhi1Array, ierr) !!Added by AM 2016-06

       do ispecies=1,Nspecies
          if (Nspecies>1) then
             print *,"Results for species ",ispecies,":"
          end if
          print *,"   FSADensityNonadiabaticPerturbation:  ", FSADensityNonadiabaticPerturbation(ispecies)
          print *,"   FSABFlow:                ", FSABFlow(ispecies)
          print *,"   max and min Mach #:      ", maxval(MachUsingFSAThermalSpeed_realSpace(ispecies,:,:)),&
               minval(MachUsingFSAThermalSpeed_realSpace(ispecies,:,:))
          print *,"   FSAPressureNonadiabaticPerturbation: ", FSAPressureNonadiabaticPerturbation(ispecies)
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

       deallocate(FourierVector)
       deallocate(FourierMatrix_DInverse)
       deallocate(FourierMatrix_BOverD)
       deallocate(FourierMatrix_particleHeat_vm)
       deallocate(FourierMatrix_momentum_vm)
       deallocate(FourierMatrix_particleHeat_vE)
       deallocate(FourierMatrix_momentum_vE)


    end if

    call VecDestroy(solutionWithFullF, ierr)
    call VecDestroy(solutionWithFullFOnProc0, ierr)
    call VecDestroy(solutionWithDeltaFOnProc0, ierr)
    call VecDestroy(f0OnProc0, ierr)
    !!!call VecDestroy(expPhi1, ierr) !!Added by AM 2016-06

    call PetscTime(presentTime, ierr)
    elapsedTime = startTime-presentTime


    ! updateOutputFile should be called by all procs since it contains MPI_Barrier
    ! (in order to be sure the HDF5 file is safely closed before moving on to the next computation.)
    if (RHSMode >1 .and. whichRHS==transportMatrixSize) then
       call updateOutputFile(iterationNum, .true.)
    else
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

