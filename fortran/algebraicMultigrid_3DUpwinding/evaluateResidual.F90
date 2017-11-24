#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif

  subroutine evaluateResidual(mysnes, stateVec, residualVec, userContext, ierr)

    use kinds
    use petscsnes
    use globalVariables
    use indices

    implicit none

    SNES :: mysnes
    Vec :: stateVec, residualVec
    PetscErrorCode :: ierr
    integer :: userContext(*)
    Vec :: inhomogeneous_terms
    PetscScalar :: scalar
    real(prec) :: x_part, factor, x_part_2, species_factor
    integer :: ix, ixi, L, itheta, izeta, ispecies, index
    real(prec) :: THat, mHat, nHat, dfMdx
    Mat :: residualMatrix
    real(prec) :: dPhiHatdpsiHatToUseInRHS
    PetscReal :: norm
    integer :: ixMin
    PetscViewer :: viewer
    character(len=200) :: filename

    if (masterProc) then
       print *,"evaluateResidual called."
    end if

    ! In the first iteration of SNES, evaluateResidual is called when the state vector is 0.
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
       call VecZeroEntries(residualVec, ierr)
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
    call VecZeroEntries(inhomogeneous_terms, ierr)

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
       nHat = nHats(ispecies)
       
       do ix=ixMin,Nx
          x_part = x2(ix)*exp(-x2(ix))*( dnHatdpsiHats(ispecies)/nHats(ispecies) &
               + Zs(ispecies)/THats(ispecies)*dPhiHatdpsiHatToUseInRHS &
               + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THats(ispecies))

          if (includePhi1 .and. includePhi1InKineticEquation) then
             x_part_2 = x2(ix)*exp(-x2(ix))*dTHatdpsiHats(ispecies)/(THats(ispecies)*THats(ispecies))
          end if

          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax                
                if (includePhi1 .and. includePhi1InKineticEquation) then !!Added by AM 2016-03
                   stop "This part not yet ready for theta_finiteDiffXi"
!!$                   factor = -spatial_scaling(itheta,izeta)*Delta*nHats(ispecies)*mHat*sqrtMHat &
!!$                        /(2*Zs(ispecies)*(BHat(itheta,izeta)**3)*sqrtTHat) &
!!$                        *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
!!$                        - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
!!$                        !!                     * DHat(itheta,izeta) * x_part !!Commented by AM 2016-02
!!$                        / sqrt_g(itheta,izeta) * (x_part + x_part_2*Zs(ispecies)*gamma*Phi1Hat(itheta,izeta))  & !!Added by AM 2016-02
!!$                        * exp(-Zs(ispecies)*gamma*Phi1Hat(itheta,izeta)/THats(ispecies)) !!Added by AM 2016-02
                else
!!$                   factor = Delta*nHats(ispecies)*mHat*sqrtMHat &
!!$                        /(2*pi*sqrtpi*Zs(ispecies)*(BHat(itheta,izeta)**3)*sqrtTHat) &
!!$                        *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
!!$                        - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
!!$                        * DHat(itheta,izeta) * x_part
                   factor = -nHat*THat &
                        /(pi*sqrtpi*thermal_speeds(ispecies)*thermal_speeds(ispecies)*thermal_speeds(ispecies)*Zs(ispecies) &
                        *BHat(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)*sqrt_g(itheta,izeta)) &
                        *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                        - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                        * x_part * spatial_scaling(itheta,izeta) * x_scaling(ix,ispecies) / f_scaling(ix,ispecies)
                end if 
                
                do ixi = 1,Nxi_for_x(ix)
                   index = getIndex(ispecies, ix, ixi, itheta, izeta, BLOCK_F)
                   scalar = (1+xi(ixi)*xi(ixi)) * factor * xi_scaling(ixi)
                   call VecSetValue(inhomogeneous_terms, index, scalar, ADD_VALUES, ierr)
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
       stop "This section not ready yet for the theta_finiteDiffXi version"
       L=0
       do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
             index = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN)

             factor = 0
             do ispecies = 1,Nspecies

!!$                call VecSetValue(rhs, index, - Zs(ispecies) * NHats(ispecies) &
!!$                     * exp (- Zs(ispecies)* gamma * Phi1Hat(itheta,izeta) / THats(ispecies)), ADD_VALUES, ierr)

                factor = factor - Zs(ispecies) * NHats(ispecies)* exp (- Zs(ispecies) * gamma * Phi1Hat(itheta,izeta) / THats(ispecies))
             end do
             
             if (withAdiabatic) then

!!$                call VecSetValue(rhs, index, - adiabaticZ * adiabaticNHat &
!!$                     * exp (- adiabaticZ* gamma * Phi1Hat(itheta,izeta) / adiabaticTHat), ADD_VALUES, ierr)

                factor = factor - adiabaticZ * adiabaticNHat * exp (- adiabaticZ* gamma * Phi1Hat(itheta,izeta) / adiabaticTHat)
             end if
             scalar = factor
             call VecSetValue(inhomogeneous_terms, index, scalar, ADD_VALUES, ierr)
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
!!$                factor = gamma*Delta*x(ix)*DHat(itheta,izeta)/(4*(BHat(itheta,izeta)**3)) &
!!$                     * (BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
!!$                     - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta)) &
!!$                     * dPhiHatdpsiHat
!!$
!!$                L = 0
!!$                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
!!$                call VecSetValue(rhs, index, (4/three)*factor, ADD_VALUES, ierr)
!!$                
!!$                L = 2
!!$                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
!!$                call VecSetValue(rhs, index, (two/three)*factor, ADD_VALUES, ierr)
!!$             end do
!!$          end do
!!$       end do
!!$    end do

    ! Add the inductive electric field term:
    do ispecies = 1,Nspecies
       do ix=ixMin,Nx
          do itheta=ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                factor = -nHats(ispecies)*Zs(ispecies)*x(ix)*exp(-x2(ix))*EParallelHat * BHat(itheta,izeta) &
                     / (pi*sqrtpi*thermal_speeds(ispecies)*thermal_speeds(ispecies)*THats(ispecies)*FSABHat2) &
                     * spatial_scaling(itheta,izeta) * x_scaling(ix,ispecies) / f_scaling(ix,ispecies)
                do ixi = 1,Nxi
                   index = getIndex(ispecies, ix, ixi, itheta, izeta, BLOCK_F)
                   scalar = factor * xi(ixi) * xi_scaling(ixi)
                   call VecSetValue(inhomogeneous_terms, index, scalar, ADD_VALUES, ierr)
                end do
             end do
          end do
       end do
    end do

    ! Done inserting values.
    ! Finally, assemble the RHS vector:
    call VecAssemblyBegin(inhomogeneous_terms, ierr)
    call VecAssemblyEnd(inhomogeneous_terms, ierr)

!    print *,"Here comes inhomogeneous_terms:"
!    call VecView(inhomogeneous_terms, PETSC_VIEWER_STDOUT_WORLD,ierr)

    ! Add the inhomogeneous terms to the residual:
    scalar = 1
    call VecAXPY(residualVec, scalar, inhomogeneous_terms, ierr)
    call VecDestroy(inhomogeneous_terms, ierr)

!    print *,"Here comes residualVec:"
!    call VecView(residualVec, PETSC_VIEWER_STDOUT_WORLD,ierr)

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
