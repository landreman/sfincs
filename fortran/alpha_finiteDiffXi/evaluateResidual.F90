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

    SNES :: mysnes
    Vec :: stateVec, residualVec
    PetscErrorCode :: ierr
    integer :: userContext(*)
    Vec :: inhomogeneous_terms
    PetscScalar :: scalar, x_part, factor, x_part_2
    integer :: ix, ixi, L, ialpha, izeta, ispecies, index
    PetscScalar :: THat, mHat, sqrtTHat, sqrtmHat, dfMdx
    Mat :: residualMatrix
    PetscScalar :: dPhiHatdpsiHatToUseInRHS
    PetscReal :: norm
    integer :: ixMin
    PetscViewer :: viewer
    character(len=200) :: filename

    if (masterProc) then
       print *,"evaluateResidual called."
    end if

    ! Often, evaluateResidual is called when the state vector is 0.
    ! In this case, there is no need to build the first matrix.
    call VecNorm(stateVec, NORM_INFINITY, norm, ierr)
    if (norm > 0) then

       ! Begin with the linear terms that are dense:
       call VecZeroEntries(residualVec,ierr)
       call apply_dense_terms(stateVec, residualVec, 1)

       ! Other terms in the residual are computed by calling populateMatrix(...,3)
       ! and multiplying the result by the state vector:
       call preallocateMatrix(residualMatrix, 3)
       call populateMatrix(residualMatrix, 3, stateVec)
       call MatMultAdd(residualMatrix, stateVec, residualVec, residualVec, ierr)
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

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, inhomogeneous_terms, ierr)
    call VecSet(inhomogeneous_terms, zero, ierr)

    if (RHSMode==1) then
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
          x_part = x2(ix)*exp(-x2(ix))*( dnHatdpsiHats(ispecies)/nHats(ispecies) &
               + gamma*Zs(ispecies)/THats(ispecies)*dPhiHatdpsiHatToUseInRHS &
               + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THats(ispecies))

          if (includePhi1 .and. includePhi1InKineticEquation) then
             x_part_2 = x2(ix)*exp(-x2(ix))*dTHatdpsiHats(ispecies)/(THats(ispecies)*THats(ispecies))
          end if

          do ialpha = ialphaMin,ialphaMax
             do izeta = izetaMinDKE,izetaMaxDKE                
                if (includePhi1 .and. includePhi1InKineticEquation) then !!Added by AM 2016-03
                   stop "This part not yet ready for alpha_finiteDiffXi"
                   factor = -Delta*nHats(ispecies)*mHat*sqrtMHat &
                        /(2*pi*sqrtpi*Zs(ispecies)*(BHat(ialpha,izeta)**3)*sqrtTHat) &
                        *(BHat_sub_zeta(ialpha,izeta)*dBHatdtheta(ialpha,izeta) &
                        - BHat_sub_theta(ialpha,izeta)*dBHatdzeta(ialpha,izeta))&
                        !!                     * DHat(ialpha,izeta) * x_part !!Commented by AM 2016-02
                        * DHat(ialpha,izeta) * (x_part + x_part_2*Zs(ispecies)*gamma*Phi1Hat(ialpha,izeta))  & !!Added by AM 2016-02
                        * exp(-Zs(ispecies)*gamma*Phi1Hat(ialpha,izeta)/THats(ispecies)) !!Added by AM 2016-02
                else
!!$                   factor = Delta*nHats(ispecies)*mHat*sqrtMHat &
!!$                        /(2*pi*sqrtpi*Zs(ispecies)*(BHat(ialpha,izeta)**3)*sqrtTHat) &
!!$                        *(BHat_sub_zeta(ialpha,izeta)*dBHatdtheta(ialpha,izeta) &
!!$                        - BHat_sub_theta(ialpha,izeta)*dBHatdzeta(ialpha,izeta))&
!!$                        * DHat(ialpha,izeta) * x_part
                   factor = -sqrt_g_sign*Delta*sqrtTHat*sqrtMHat &
                        /(2*pi*sqrtpi*Zs(ispecies)*BHat(ialpha,izeta)*BHat(ialpha,izeta)) &
                        *(BHat_sub_zeta(ialpha,izeta)*dBHatdtheta(ialpha,izeta) &
                        - BHat_sub_theta(ialpha,izeta)*dBHatdzeta(ialpha,izeta))&
                        * x_part
                end if 
                
                do ixi = 1,Nxi_for_x(ix)
                   index = getIndex(ispecies, ix, ixi, ialpha, izeta, BLOCK_F)
                   call VecSetValue(inhomogeneous_terms, index, (1+xi(ixi)*xi(ixi))*factor, ADD_VALUES, ierr)
                end do
             end do
          end do
       end do
    end do
    
    ! *******************************************************************************
    ! SECTION ADDED BY AM 2016-02/03
    ! Add part of the residual related to Phi1 in the quasineutrality equation
    ! *******************************************************************************

    if (includePhi1 .and. quasineutralityOption == 1) then
       stop "This section not ready yet for the alpha_finiteDiffXi version"
       L=0
       do ialpha = ialphaMin,ialphaMax
          do izeta = izetaMin,izetaMax
             index = getIndex(1, 1, 1, ialpha, izeta, BLOCK_QN)

             factor = 0
             do ispecies = 1,Nspecies

!!$                call VecSetValue(rhs, index, - Zs(ispecies) * NHats(ispecies) &
!!$                     * exp (- Zs(ispecies)* gamma * Phi1Hat(ialpha,izeta) / THats(ispecies)), ADD_VALUES, ierr)

                factor = factor - Zs(ispecies) * NHats(ispecies)* exp (- Zs(ispecies) * gamma * Phi1Hat(ialpha,izeta) / THats(ispecies))
             end do
             
             if (withAdiabatic) then

!!$                call VecSetValue(rhs, index, - adiabaticZ * adiabaticNHat &
!!$                     * exp (- adiabaticZ* gamma * Phi1Hat(ialpha,izeta) / adiabaticTHat), ADD_VALUES, ierr)

                factor = factor - adiabaticZ * adiabaticNHat * exp (- adiabaticZ* gamma * Phi1Hat(ialpha,izeta) / adiabaticTHat)
             end if
             call VecSetValue(inhomogeneous_terms, index, factor, ADD_VALUES, ierr)
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
!!$          do ialpha = ialphaMin,ialphaMax
!!$             do izeta = 1,Nzeta
!!$                
!!$                factor = gamma*Delta*x(ix)*DHat(ialpha,izeta)/(4*(BHat(ialpha,izeta)**3)) &
!!$                     * (BHat_sub_zeta(ialpha,izeta)*dBHatdtheta(ialpha,izeta) &
!!$                     - BHat_sub_theta(ialpha,izeta)*dBHatdzeta(ialpha,izeta)) &
!!$                     * dPhiHatdpsiHat
!!$
!!$                L = 0
!!$                index = getIndex(ispecies, ix, L+1, ialpha, izeta, BLOCK_F)
!!$                call VecSetValue(rhs, index, (4/three)*factor, ADD_VALUES, ierr)
!!$                
!!$                L = 2
!!$                index = getIndex(ispecies, ix, L+1, ialpha, izeta, BLOCK_F)
!!$                call VecSetValue(rhs, index, (two/three)*factor, ADD_VALUES, ierr)
!!$             end do
!!$          end do
!!$       end do
!!$    end do

    ! Add the inductive electric field term:
    do ispecies = 1,Nspecies
       do ix=ixMin,Nx
          do ialpha=ialphaMin,ialphaMax
             do izeta = izetaMinDKE,izetaMaxDKE
                !factor = gamma*Zs(ispecies)*x(ix)*exp(-x2(ix))*EParallelHat &
                !     *nHats(ispecies)*mHats(ispecies)/(pi*sqrtpi*THats(ispecies)*THats(ispecies)*FSABHat2)
                factor = -gamma*Zs(ispecies)*x(ix)*exp(-x2(ix))*EParallelHat &
                     *BHat(ialpha,izeta)*BHat(ialpha,izeta)/(pi*sqrtpi*THats(ispecies)*abs(DHat(ialpha,izeta))*FSABHat2)
                do ixi = 1,Nxi
                   index = getIndex(ispecies, ix, ixi, ialpha, izeta, BLOCK_F)
                   call VecSetValue(inhomogeneous_terms, index, &
                        factor * xi(ixi), ADD_VALUES, ierr)
                end do
             end do
          end do
       end do
    end do

    ! Done inserting values.
    ! Finally, assemble the RHS vector:
    call VecAssemblyBegin(inhomogeneous_terms, ierr)
    call VecAssemblyEnd(inhomogeneous_terms, ierr)


    ! Add the inhomogeneous terms to the residual:
    scalar = 1
    call VecAXPY(residualVec, scalar, inhomogeneous_terms, ierr)
    call VecDestroy(inhomogeneous_terms, ierr)

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
