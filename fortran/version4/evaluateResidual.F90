#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif

  subroutine evaluateResidual(mysnes, stateVec, residualVec, userContext, ierr)

    use petscsnes
    use FourierTransformMod
    use FourierConvolutionMatrixMod
    use globalVariables
    use indices

    implicit none

    SNES :: mysnes
    Vec :: stateVec, residualVec
    PetscErrorCode :: ierr
    integer :: userContext(*)
    Vec :: termsInResidual
    real(prec) :: factor, stuffToAdd
    integer :: ix, L, imn, ispecies, index
    real(prec) :: THat, mHat, sqrtTHat, sqrtmHat, nHat, Z
    Mat :: residualMatrix
    real(prec) :: dPhiHatdpsiHatToUseInDriveTerm
    PetscReal :: norm
    integer :: ixMin
    PetscViewer :: viewer
    character(len=200) :: filename
    PetscScalar :: PetscScalarZero, PetscScalarValue
    real(prec), dimension(:), allocatable :: FourierVector, FourierVector2
    real(prec), dimension(:,:), allocatable :: spatialFactor
    real(prec), dimension(:), allocatable :: xPart

    if (masterProc) then
       print *,"evaluateResidual called."
    end if

    ! This next line is used to cast 0 into either double or single precision, whichever is used in PETSc:
    PetscScalarZero = 0

    allocate(FourierVector(NFourier2))
    allocate(FourierVector2(NFourier2))
    allocate(spatialFactor(Ntheta,Nzeta)
    allocate(xPart(Nx))

    ! Often, evaluateResidual is called when the state vector is 0.
    ! In this case, there is no need to build the first matrix.
    call VecNorm(stateVec, NORM_INFINITY, norm, ierr)
    if (norm > 0) then

       ! Some terms in the residual are computed by calling populateMatrix(...,3)
       ! and multiplying the result by the state vector:
       call preallocateMatrix(residualMatrix, 3)
       call populateMatrix(residualMatrix, 3, stateVec)
       call MatMult(residualMatrix, stateVec, residualVec, ierr)
       call MatDestroy(residualMatrix, ierr)

    else
       if (masterProc) then
          print *,"State vector is 0 so I will skip building the first matrix when evaluating the residual."
       end if
       call VecSet(residualVec, PetscScalarZero, ierr)
    end if

    ! The collision term (temperature equilibration) in the residual is computed by calling populateMatrix(...,2)
    ! any multiplying the result by the Vec f0:
    if (includeTemperatureEquilibrationTerm) then
       call preallocateMatrix(residualMatrix, 2)
       call populateMatrix(residualMatrix, 2, stateVec)
       call MatMultAdd(residualMatrix, f0, residualVec, residualVec, ierr)
       call MatDestroy(residualMatrix, ierr)
    end if

    ! --------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------
    ! MODIFIED BY AM 2016-02/03
    ! Next, evaluate the terms (other than C[f_M]) that are independent of both f_1 and Phi_1:
    ! THESE TERMS USED TO BE INDEPENDENT OF PHI1, BUT NOT ANYMORE
    ! --------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, termsInResidual, ierr)
    call VecSet(termsInResidual, PetscScalarZero, ierr)

    if (RHSMode==1) then
       dPhiHatdpsiHatToUseInDriveTerm = dPhiHatdpsiHat
    else
       dPhiHatdpsiHatToUseInDriveTerm = 0
    end if

    if (pointAtX0) then
       ixMin = 2
    else
       ixMin = 1
    end if

    !!Added by AM 2016-03!!
    if (includePhi1) then
       call extractPhi1(stateVec)
    end if
    !!!!!!!!!!!!!!!!!!!!!!!

    ! Do the next part all on the master proc.
    if (masterProc) then
       ! First add the terms arising from \dot{\psi} df_M/d\psi
       x2 = x*x
       xPart = x2*exp(-x2)
       do ispecies = 1,Nspecies
          THat = THats(ispecies)
          mHat = mHats(ispecies)
          nHat = nHats(ispecies)
          Z    =    Zs(ispecies)
          sqrtTHat = sqrt(THat)
          sqrtMHat = sqrt(mHat)
          
          ! Handle spatial part of the residual terms R_m:
          ! This part needs to be inside the species loop because the exp(-Z e Phi1/T) factor depends on species.
          spatialFactor = exp(-Phi1Hat_realSpace*Z*alpha/THat) &
               * DHat*(BHat_sub_theta*dBHatdzeta-BHat_sub_zeta*dBHatdtheta)/(BHat*BHat*BHat)
          call FourierTransform(spatialFactor,FourierVector)
          if (includePhi1 .and. includePhi1InKineticEquation) then
             call FourierTransform((Z*alpha*dTHatdpsiHat(ispecies)/(THat*THat))*Phi1Hat_realSpace * spatialFactor, FourierVector2)
          end if

          factor = Delta*nHat*mHat*sqrtMHat/(2*pi*sqrtpi*Z*sqrtTHat)

          do ix=ixMin,Nx
             do imn = 1,NFourier2
                ! stuffToAdd includes the dependence on species, space, and x, but not the Legendre dependence.
                stuffToAdd = factor*xPart(ix)*FourierVector(imn)*(dNHatdpsiHats(ispecies)/nHat &
                     + Z*alpha*dPhiHatdpsiHatToUseInDriveTerm/THat+(x2(ix)-three/two)*dTHatdpsiHat(ispecies)/THat)
                if (includePhi1 .and. includePhi1InKineticEquation) then
                   stuffToAdd = stuffToAdd + factor*xPart(ix)*FourierVector2(imn)
                end if

                L = 0
                index = getIndex(ispecies, ix, L+1, imn, BLOCK_F)
                PetscScalarValue = (4/three)*stuffToAdd ! Cast real(prec) to PetscScalar
                call VecSetValue(termsInResidual, index, PetscScalarValue, INSERT_VALUES, ierr)
                
                L = 2
                index = getIndex(ispecies, ix, L+1, imn, BLOCK_F)
                PetscScalarValue = (two/three)*stuffToAdd ! Cast real(prec) to PetscScalar
                call VecSetValue(termsInResidual, index, PetscScalarValue, INSERT_VALUES, ierr)

             end do
          end do
       end do

! This next block of code has not yet been updated for the Fourier discretization.       
!!$       ! *******************************************************************************
!!$       ! SECTION ADDED BY AM 2016-02/03
!!$       ! Add part of the residual related to Phi1 in the quasineutrality equation
!!$       ! *******************************************************************************
!!$       
!!$       if (includePhi1 .and. quasineutralityOption == 1) then
!!$          L=0
!!$          do itheta = ithetaMin,ithetaMax
!!$             do izeta = izetaMin,izetaMax
!!$                index = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN)
!!$                
!!$                PetscScalarValue = 0
!!$                do ispecies = 1,Nspecies
!!$                   
!!$!                call VecSetValue(rhs, index, - Zs(ispecies) * NHats(ispecies) &
!!$!                     * exp (- Zs(ispecies)* alpha * Phi1Hat(itheta,izeta) / THats(ispecies)), ADD_VALUES, ierr)
!!$                   
!!$                   PetscScalarValue = PetscScalarValue - Zs(ispecies) * NHats(ispecies)* exp (- Zs(ispecies) * alpha * Phi1Hat(itheta,izeta) / THats(ispecies))
!!$                end do
!!$                
!!$                if (withAdiabatic) then
!!$                   
!!$!                call VecSetValue(rhs, index, - adiabaticZ * adiabaticNHat &
!!$!                     * exp (- adiabaticZ* alpha * Phi1Hat(itheta,izeta) / adiabaticTHat), ADD_VALUES, ierr)
!!$                   
!!$                   PetscScalarValue = PetscScalarValue - adiabaticZ * adiabaticNHat * exp (- adiabaticZ* alpha * Phi1Hat(itheta,izeta) / adiabaticTHat)
!!$                end if
!!$                call VecSetValue(rhs, index, PetscScalarValue, INSERT_VALUES, ierr)
!!$             end do
!!$          end do
!!$       end if
!!$       ! *******************************************************************************
       
       
       ! Add the inductive electric field term:
       call FourierTransform(BHat,FourierVector)
       L=1
       do ispecies = 1,Nspecies
          do ix=ixMin,Nx
             factor = -alpha*Zs(ispecies)*x(ix)*exp(-x2(ix))*EParallelHat &
                  *nHats(ispecies)*mHats(ispecies)/(pi*sqrtpi*THats(ispecies)*THats(ispecies)*FSABHat2)
             do imn = 1,NFourier2
                index = getIndex(ispecies, ix, L+1, imn, BLOCK_F)
                PetscScalarValue = factor * FourierVector(imn)
                call VecSetValue(termsInResidual, index, &
                     PetscScalarValue, INSERT_VALUES, ierr)
             end do
          end do
       end do

    end if ! End of section for only the master proc.

    ! Done inserting values.
    ! Finally, assemble the TERMSINRESIDUAL vector:
    call VecAssemblyBegin(termsInResidual, ierr)
    call VecAssemblyEnd(termsInResidual, ierr)


    ! Add the new terms to the residual:
    PetscScalarValue = 1 ! Cast 1 into type PetscScalar before passing it as a parameter.
    call VecAXPY(residualVec, PetscScalarValue, termsInResidual, ierr)
    call VecDestroy(termsInResidual, ierr)

    if (saveMatlabOutput) then
       write (filename,fmt="(a,i3.3,a)") trim(MatlabOutputFilename) // "_iteration_", iterationForResidualOutput, &
            "_residual.m"
       if (masterProc) then
          print *,"Saving residual in matlab format: ",trim(filename)
       end if
       call PetscViewerASCIIOpen(MPIComm, trim(filename), viewer, ierr)
       call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB, ierr)

       call PetscObjectSetName(residualVec, "residual", ierr)
       call VecView(residualVec, viewer, ierr)

       call PetscViewerDestroy(viewer, ierr)
    end if

    if (saveMatricesAndVectorsInBinary) then
       write (filename,fmt="(a,i3.3,a)") trim(binaryOutputFilename) // "_iteration_", iterationForResidualOutput, &
            "_residual"
       if (masterProc) then
          print *,"Saving residual in binary format: ",trim(filename)
       end if
       call PetscViewerBinaryOpen(MPIComm, trim(filename), FILE_MODE_WRITE, viewer, ierr)
       call VecView(residualVec, viewer, ierr)
       call PetscViewerDestroy(viewer, ierr)
    end if

    iterationForResidualOutput = iterationForResidualOutput + 1
    deallocate(FourierVector,FourierVector2,spatialFactor)

  end subroutine evaluateResidual
