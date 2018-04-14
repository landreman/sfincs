#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif

  subroutine evaluateResidual(mysnes, stateVec, residualVec, userContext, ierr)

    use petscsnes
    use globalVariables
    use indices

    implicit none

    integer :: whichMatrix
    SNES :: mysnes
    Vec :: stateVec, residualVec
    PetscErrorCode :: ierr
    integer :: userContext(*)
    Vec :: rhs
!!    PetscScalar :: scalar, xPartOfRHS, factor !!Commented by AM 2016-03
    PetscScalar :: scalar, xPartOfRHS, factor, xPartOfRHS2 !!Added by AM 2016-03
    integer :: ix, L, itheta, izeta, ispecies, index
    PetscScalar :: THat, mHat, sqrtTHat, sqrtmHat, dfMdx
    Mat :: residualMatrix
    PetscScalar :: dPhiHatdpsiHatToUseInRHS
    PetscReal :: norm
    integer :: ixMin
    PetscViewer :: viewer
    character(len=200) :: filename
    integer :: this_matrixSize, this_ithetaMin, this_ithetaMax, this_izetaMin, this_izetaMax
    PetscScalar, dimension(:,:), allocatable :: this_DHat, this_BHat, this_BHat_sub_zeta, this_dBHatdtheta, this_dBHatdzeta, this_BHat_sub_theta

    whichMatrix = userContext(1)
    print *,"whichMatrix: ", whichMatrix

    if (whichMatrix == 7) then
      allocate(this_DHat(Ntheta_fine,Nzeta_fine))
      allocate(this_BHat(Ntheta_fine,Nzeta_fine))
      allocate(this_BHat_sub_zeta(Ntheta_fine,Nzeta_fine))
      allocate(this_BHat_sub_theta(Ntheta_fine,Nzeta_fine))
      allocate(this_dBHatdtheta(Ntheta_fine,Nzeta_fine))
      allocate(this_dBHatdzeta(Ntheta_fine,Nzeta_fine))

      this_matrixSize = matrixSize_fine
      this_ithetaMin = ithetaMin_fine
      this_ithetaMax = ithetaMax_fine
      this_izetaMin = izetaMin_fine
      this_izetaMax = izetaMax_fine
      this_DHat = DHat_fine
      this_BHat = BHat_fine
      this_BHat_sub_zeta = BHat_sub_zeta_fine
      this_BHat_sub_theta = BHat_sub_theta_fine
      this_dBHatdtheta = dBHatdtheta_fine
      this_dBhatdzeta = dBHatdzeta_fine
    else
      allocate(this_DHat(Ntheta,Nzeta))
      allocate(this_BHat(Ntheta,Nzeta))
      allocate(this_BHat_sub_zeta(Ntheta,Nzeta))
      allocate(this_BHat_sub_theta(Ntheta,Nzeta))
      allocate(this_dBHatdtheta(Ntheta,Nzeta))
      allocate(this_dBHatdzeta(Ntheta,Nzeta))

      this_matrixSize = matrixSize
      this_ithetaMin = ithetaMin
      this_ithetaMax = ithetaMax
      this_izetaMin = izetaMin
      this_izetaMax = izetaMax
      this_DHat = DHat
      this_BHat = BHat
      this_BHat_sub_zeta = BHat_sub_zeta
      this_BHat_sub_theta = BHat_sub_theta
      this_dBHatdtheta = dBHatdtheta
      this_dBhatdzeta = dBHatdzeta
    end if
!    print *,"this_dBHatdzeta: ", this_dBHatdzeta

    if (masterProc) then
       print *,"evaluateResidual called."
    end if

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
       call VecSet(residualVec, zero, ierr)
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

    call VecCreateMPI(MPIComm, PETSC_DECIDE, this_matrixSize, rhs, ierr)
    call VecSet(rhs, zero, ierr)

    if (RHSMode==1 .or. RHSMode>3) then
       dPhiHatdpsiHatToUseInRHS = dPhiHatdpsiHat
    else
       dPhiHatdpsiHatToUseInRHS = 0
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

    ! First add the terms arising from \dot{\psi} df_M/d\psi
    x2 = x*x
    do ispecies = 1,Nspecies
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)
       
       do ix=ixMin,Nx
          ! Old linear version:
          !xPartOfRHS = x2(ix)*exp(-x2(ix))*( dnHatdpsiNs(ispecies)/nHats(ispecies) &
          !     + alpha*Zs(ispecies)/THats(ispecies)*dPhiHatdpsiNToUse &
          !     + (x2(ix) - three/two)*dTHatdpsiNs(ispecies)/THats(ispecies))

          ! Old nonlinear version:
          xPartOfRHS = x2(ix)*exp(-x2(ix))*( dnHatdpsiHats(ispecies)/nHats(ispecies) &
               + alpha*Zs(ispecies)/THats(ispecies)*dPhiHatdpsiHatToUseInRHS &
               + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THats(ispecies))

          !!Added by AM 2016-02!!
          if (includePhi1 .and. includePhi1InKineticEquation) then
             xPartOfRHS2 = x2(ix)*exp(-x2(ix))*dTHatdpsiHats(ispecies)/(THats(ispecies)*THats(ispecies))
          !!else
          !!   xPartOfRHS2 = 0.0
          end if
          !!!!!!!!!!!!!!!!!!!!!!!

          !xPartOfRHS = x2(ix)*exp(-x2(ix))*( dnHatdpsiHats(ispecies)/nHats(ispecies) &
          !     + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THats(ispecies))

          do itheta = this_ithetaMin,this_ithetaMax
             do izeta = this_izetaMin,this_izetaMax
                
                !factor = Delta*nHats(ispecies)*mHat*sqrtMHat &
                !     /(2*pi*sqrtpi*Zs(ispecies)*psiAHat*(BHat(itheta,izeta)**3)*sqrtTHat) &
                !     *(GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))&
                !     *xPartOfRHS
                if (includePhi1 .and. includePhi1InKineticEquation) then !!Added by AM 2016-03
                   factor = Delta*nHats(ispecies)*mHat*sqrtMHat &
                        /(2*pi*sqrtpi*Zs(ispecies)*(BHat(itheta,izeta)**3)*sqrtTHat) &
                        *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                        - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                        !!                     * DHat(itheta,izeta) * xPartOfRHS !!Commented by AM 2016-02
                        * DHat(itheta,izeta) * (xPartOfRHS + xPartOfRHS2*Zs(ispecies)*alpha*Phi1Hat(itheta,izeta))  & !!Added by AM 2016-02
                        * exp(-Zs(ispecies)*alpha*Phi1Hat(itheta,izeta)/THats(ispecies)) !!Added by AM 2016-02
                else !!Added by AM 2016-03
                   factor = Delta*nHats(ispecies)*mHat*sqrtMHat &
                        /(2*pi*sqrtpi*Zs(ispecies)*(this_BHat(itheta,izeta)**3)*sqrtTHat) &
                        *(this_BHat_sub_zeta(itheta,izeta)*this_dBHatdtheta(itheta,izeta) &
                        - this_BHat_sub_theta(itheta,izeta)*this_dBHatdzeta(itheta,izeta))&
                        * this_DHat(itheta,izeta) * xPartOfRHS
                end if !!Added by AM 2016-03
                
                
                L = 0
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,whichMatrix)
                call VecSetValue(rhs, index, (4/three)*factor, INSERT_VALUES, ierr)
                
                L = 2
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,whichMatrix)
                call VecSetValue(rhs, index, (two/three)*factor, INSERT_VALUES, ierr)
             end do
          end do
       end do
    end do
    
    ! *******************************************************************************
    ! SECTION ADDED BY AM 2016-02/03
    ! Add part of the residual related to Phi1 in the quasineutrality equation
    ! *******************************************************************************

    if (includePhi1 .and. quasineutralityOption == 1) then
       L=0
       do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
             index = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN,whichMatrix)

             factor = 0
             do ispecies = 1,Nspecies

!!$                call VecSetValue(rhs, index, - Zs(ispecies) * NHats(ispecies) &
!!$                     * exp (- Zs(ispecies)* alpha * Phi1Hat(itheta,izeta) / THats(ispecies)), ADD_VALUES, ierr)

                factor = factor - Zs(ispecies) * NHats(ispecies)* exp (- Zs(ispecies) * alpha * Phi1Hat(itheta,izeta) / THats(ispecies))
             end do
             
             if (withAdiabatic) then

!!$                call VecSetValue(rhs, index, - adiabaticZ * adiabaticNHat &
!!$                     * exp (- adiabaticZ* alpha * Phi1Hat(itheta,izeta) / adiabaticTHat), ADD_VALUES, ierr)

                factor = factor - adiabaticZ * adiabaticNHat * exp (- adiabaticZ* alpha * Phi1Hat(itheta,izeta) / adiabaticTHat)
             end if
             call VecSetValue(rhs, index, factor, INSERT_VALUES, ierr)
          end do
       end do
    end if
    ! *******************************************************************************


!!$    ! Next add the terms arising from \dot{x} df_M/dx
!!$    do ispecies = 1,Nspecies
!!$       THat = THats(ispecies)
!!$       mHat = mHats(ispecies)
!!$       sqrtTHat = sqrt(THat)
!!$       sqrtMHat = sqrt(mHat)
!!$       
!!$       do ix=1,Nx
!!$          dfMdx = -2*x(ix)*mHat*sqrtmHat*nHat*expx2(ix)/(pi*sqrtpi*THat*sqrtTHat)
!!$
!!$          do itheta = ithetaMin,ithetaMax
!!$             do izeta = 1,Nzeta
!!$                
!!$                factor = alpha*Delta*x(ix)*DHat(itheta,izeta)/(4*(BHat(itheta,izeta)**3)) &
!!$                     * (BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
!!$                     - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta)) &
!!$                     * dPhiHatdpsiHat
!!$
!!$                L = 0
!!$                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
!!$                call VecSetValue(rhs, index, (4/three)*factor, INSERT_VALUES, ierr)
!!$                
!!$                L = 2
!!$                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
!!$                call VecSetValue(rhs, index, (two/three)*factor, INSERT_VALUES, ierr)
!!$             end do
!!$          end do
!!$       end do
!!$    end do

    ! Add the inductive electric field term:
    L=1
    do ispecies = 1,Nspecies
       do ix=ixMin,Nx
          !factor = alpha*Zs(ispecies)*x(ix)*exp(-x2(ix))*EParallelHatToUse*(GHat+iota*IHat)&
          !     *nHats(ispecies)*mHats(ispecies)/(pi*sqrtpi*THats(ispecies)*THats(ispecies)*FSABHat2)
          factor = alpha*Zs(ispecies)*x(ix)*exp(-x2(ix))*(EParallelHat+EParallelHatSpec(ispecies)) & !!EParallelHatSpec added 2017-02/HS
               *nHats(ispecies)*mHats(ispecies)/(pi*sqrtpi*THats(ispecies)*THats(ispecies)*FSABHat2)
          do itheta=this_ithetaMin,this_ithetaMax
             do izeta = this_izetaMin,this_izetaMax
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,whichMatrix)
                call VecSetValue(rhs, index, &
                     factor * this_BHat(itheta,izeta), INSERT_VALUES, ierr)
                     !factor/BHat(itheta,izeta), INSERT_VALUES, ierr)
             end do
          end do
       end do
    end do

    ! Done inserting values.
    ! Finally, assemble the RHS vector:
    call VecAssemblyBegin(rhs, ierr)
    call VecAssemblyEnd(rhs, ierr)


    ! Subtract the RHS from the residual:
    scalar = -1
    call VecAXPY(residualVec, scalar, rhs, ierr)
    call VecDestroy(rhs, ierr)

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

  end subroutine evaluateResidual
