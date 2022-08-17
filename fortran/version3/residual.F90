module residual

#include "PETScVersions.F90"

  use globalVariables
  use indices
  use DKEMatrix
  use diagnostics
  
  implicit none

  private
  PetscScalar :: dPhiHatdpsiHatToUseInRHS
  logical, parameter :: localUsePhi1 = .true.

  public ::  evaluateResidual
  public ::  DKERHS

contains
  
subroutine evaluateResidual(mysnes, stateVec, residualVec, userContext, ierr)

    SNES :: mysnes
    Vec :: stateVec, residualVec
    PetscErrorCode :: ierr
    integer :: userContext(*)
    Vec :: rhs
    PetscScalar :: scalar, xPartOfRHS, factor, xPartOfRHS2 !!Added by AM 2016-03
    integer :: ix, L, itheta, izeta, ispecies, index
    PetscScalar :: THat, mHat, sqrtTHat, sqrtmHat
    Mat :: residualMatrix
    PetscReal :: norm
    integer :: ixMin, LMax
    PetscViewer :: viewer
    character(len=200) :: filename

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
       call populateMatrix(residualMatrix, 3, stateVec, 0, 0)
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
       call populateMatrix(residualMatrix, 2, stateVec, 0, 0)
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

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, rhs, ierr)
    call VecSet(rhs, zero, ierr)

    if (RHSMode==1 .or. RHSMode==4 .or. RHSMode==5) then
       dPhiHatdpsiHatToUseInRHS = dPhiHatdpsiHat
    else
       dPhiHatdpsiHatToUseInRHS = 0
    end if

    if (pointAtX0) then
       ixMin = 2
    else
       ixMin = 1
    end if
    
    if (includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
       call extractPhi1(stateVec)
    end if

    call DKERHS(rhs, -1, -1)
    
    ! *******************************************************************************
    ! SECTION ADDED BY AM 2016-02/03
    ! Add part of the residual related to Phi1 in the quasineutrality equation
    ! *******************************************************************************

    if (includePhi1 .and. quasineutralityOption == 1 .and. (.not. readExternalPhi1)) then 
       L=0
       do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
             index = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN)

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
             
             ! Added by HS 2017-09
             if (withNBIspec) then
                factor = factor - NBIspecZ * NBIspecNHat * BHat(itheta,izeta)/FSABHat
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
    if (trim(EParallelHatSpec_bcdatFile)=="") then !This is the normal case
       do ispecies = 1,Nspecies
          do ix=ixMin,Nx
             factor = alpha*Zs(ispecies)*x(ix)*exp(-x2(ix))*(EParallelHat+EParallelHatSpec(ispecies)) & !!EParallelHatSpec added 2017-02/HS
                  *nHats(ispecies)*mHats(ispecies)/(pi*sqrtpi*THats(ispecies)*THats(ispecies)*FSABHat2)
             do itheta=ithetaMin,ithetaMax
                do izeta = izetaMin,izetaMax
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                   call VecSetValue(rhs, index, &
                        factor * BHat(itheta,izeta), INSERT_VALUES, ierr)
                   !factor/BHat(itheta,izeta), INSERT_VALUES, ierr)
                end do
             end do
          end do
       end do
    else !This is a very specialized case, added 2017-08/HS
      do ispecies = 1,Nspecies
          do ix=ixMin,Nx
             factor = alpha*Zs(ispecies)*x(ix)*exp(-x2(ix)) & 
                  *nHats(ispecies)*mHats(ispecies)/(pi*sqrtpi*THats(ispecies)*THats(ispecies)*FSABHat2)
             do itheta=ithetaMin,ithetaMax
                do izeta = izetaMin,izetaMax
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                   call VecSetValue(rhs, index, &
                        factor * (EParallelHat+EParallelHatSpec(ispecies)*bcdata(itheta,izeta)) &
                        *BHat(itheta,izeta), INSERT_VALUES, ierr)
                   !factor/BHat(itheta,izeta), INSERT_VALUES, ierr)
                end do
             end do
          end do
       end do 
    end if

    ! Done inserting values.

     ! After inserting, we can add values
    if (readExternalF) then
       do ispecies = 1,Nspecies
          do ix=ixMin,Nx
             LMax = min(Nxi_for_x(ix),externalNL) -1
             do L=0,Lmax
                do itheta=ithetaMin,ithetaMax
                   do izeta = izetaMin,izetaMax
                      index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                      if (isnan(externalRosenPotentialTerms(ispecies,itheta,izeta,L+1,ix))) then
                         print *,"NANANANANA"
                      end if

                      call VecSetValue(rhs, index, &
                          externalRosenPotentialTerms(ispecies,itheta,izeta,L+1,ix), ADD_VALUES, ierr)
                   end do

                end do
             end do
          end do
       end do
    end if
    
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

  ! Add the terms arising from \dot{\psi} df_M/d\psi
  subroutine DKERHS(rhs, whichLambda, whichMode)
    ! The part of the right hand side
    ! from the drift kinetic equation (DKE)
    Vec :: rhs
    integer, intent(in) :: whichLambda, whichMode
    PetscErrorCode :: ierr
    integer :: ix, L, itheta, izeta, ispecies, index
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat
    PetscScalar :: scalar, xPartOfRHS, factor, xPartOfRHS1, xPartOfRHS2, stuffToAdd
    PetscScalar :: factorExternalPhi1, xPartOfRHSExternalPhi1, stuffToAddExternalPhi1
    PetscScalar :: angle, cos_angle, sin_angle, dBHatdLambda, dBHatdthetadLambda, dBHatdzetadLambda
    PetscScalar :: dPhi1HatdLambda, dPhi1HatdthetadLambda, dPhi1HatdzetadLambda
    PetscScalar :: dxPartOfRHSExternalPhi1dLambda, dfactorExternalPhi1dLambda
    integer :: m,n

    if (whichLambda > 0) then
       m = ms_sensitivity(whichMode)
       n = ns_sensitivity(whichMode)
    end if

    !!!!!!!!!!1
    if (RHSMode==1 .or. RHSMode==4 .or. RHSMode==5) then
       dPhiHatdpsiHatToUseInRHS = dPhiHatdpsiHat
    else
       dPhiHatdpsiHatToUseInRHS = 0
    end if

    !!!!!!!!!!!!!!

    ! zero is default value, will be overriden if needed below
    stuffToAddExternalPhi1 = zero

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)
       scalar = Delta*nHat*mHat*sqrtMHat/(two*pi*sqrtpi*Z*sqrtTHat)

       do ix=ixMin,Nx
          ! Old nonlinear version:
          xPartOfRHS1 = x2(ix)*exp(-x2(ix))*( dnHatdpsiHats(ispecies)/nHat &
               + alpha*Z/THat*dPhiHatdpsiHatToUseInRHS &
               + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THat)

          !!Added by AM 2016-02!!
          if (includePhi1 .and. includePhi1InKineticEquation .and. localUsePhi1) then
             xPartOfRHS2 = x2(ix)*exp(-x2(ix))*dTHatdpsiHats(ispecies)/(THat*THat)
          end if
!!!!!!!!!!!!!!!!!!!!!!!

          ! factor = (1/(BHat*(GHat+iota*IHat)))*(GHat*dBHatdtheta-IHat*dBHatdzeta)
          !
          
          
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                if (whichLambda >= 0) then
                   select case(whichLambda)
                   case (0) ! Er
                      factor =  1/(BHat(itheta,izeta)**3) &
                        *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                        - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                        * DHat(itheta,izeta)
                      xPartOfRHS = x2(ix)* exp(-x2(ix))* (alpha*Z/THat)
                      
                      if (readExternalPhi1 .and. includePhi1 .and. includePhi1InKineticEquation  .and. localUsePhi1) then
                         factorExternalPhi1 = (alpha*Z/THat) &
                              *1/(BHat(itheta,izeta)**2) & 
                              *(BHat_sub_zeta(itheta,izeta)*dPhi1Hatdtheta(itheta,izeta) & 
                              - BHat_sub_theta(itheta,izeta)*dPhi1Hatdzeta(itheta,izeta))& 
                              * DHat(itheta,izeta) 
                         xPartOfRHSExternalPhi1 = xPartOfRHS/x2(ix)
                         stuffToAddExternalPhi1 = scalar * factorExternalPhi1 * xPartOfRHSExternalPhi1
                      end if
                      
                         
                   case (1) ! BHat cos
                      angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                      cos_angle = cos(angle)
                      sin_angle = sin(angle)                   
                      dBHatdLambda = cos_angle
                      dBHatdthetadLambda = -m*sin_angle
                      dBHatdzetadLambda = n*Nperiods*sin_angle
                      factor = -dBHatdLambda/(BHat(itheta,izeta)*BHat(itheta,izeta)*(GHat+iota*IHat)) &
                           * (GHat*dBHatdtheta(itheta,izeta)-IHat*dBHatdzeta(itheta,izeta)) &
                           + (GHat*dBHatdthetadLambda-IHat*dBHatdzetadLambda)/(BHat(itheta,izeta)*(GHat+iota*IHat))
                      if (includePhi1 .and. includePhi1InKineticEquation  .and. localUsePhi1) then
                         xPartOfRHS = (xPartOfRHS1 + xPartOfRHS2*Z*alpha*Phi1Hat(itheta,izeta))  & 
                              * exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)
                         if (readExternalPhi1) then
                            factorExternalPhi1 = 0
                            ! everything zero below
                            xPartOfRHSExternalPhi1  = xPartOfRHS/x2(ix)
                            stuffToAddExternalPhi1 = scalar * factorExternalPhi1 * xPartOfRHSExternalPhi1
                         end if   
                      else
                         xPartOfRHS = xPartOfRHS1
                      end if
                                   
                   case (2) ! IHat
                      factor = -iota/(BHat(itheta,izeta)*(GHat+iota*IHat)**2) &
                           *(GHat*dBHatdtheta(itheta,izeta)-IHat*dBHatdzeta(itheta,izeta)) &
                           -dBHatdzeta(itheta,izeta)/(BHat(itheta,izeta)*(GHat+iota*IHat))
                      if (includePhi1 .and. includePhi1InKineticEquation  .and. localUsePhi1) then
                         xPartOfRHS = (xPartOfRHS1 + xPartOfRHS2*Z*alpha*Phi1Hat(itheta,izeta))  & 
                              * exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)
                         if (readExternalPhi1) then
                            factorExternalPhi1 =  (alpha*Z/THat) &
                              * (-dPhi1Hatdzeta(itheta,izeta)/(GHat + iota * IHat) &
                              - (iota/(GHat + iota * IHat)**2) * (GHat*dPhi1Hatdtheta(itheta,izeta) & 
                              - IHat*dPhi1Hatdzeta(itheta,izeta)))
                            xPartOfRHSExternalPhi1  = xPartOfRHS/x2(ix)
                            stuffToAddExternalPhi1 = scalar * factorExternalPhi1 * xPartOfRHSExternalPhi1
                         end if   
                      else
                         xPartOfRHS = xPartOfRHS1
                      end if
                                        
                   case (3) ! GHat
                      factor = -1/(BHat(itheta,izeta)*(GHat+iota*IHat)**2) &
                           *(GHat*dBHatdtheta(itheta,izeta)-IHat*dBHatdzeta(itheta,izeta)) &
                           + dBHatdtheta(itheta,izeta)/(BHat(itheta,izeta)*(GHat+iota*IHat))
                      if (includePhi1 .and. includePhi1InKineticEquation .and. localUsePhi1) then
                         xPartOfRHS = (xPartOfRHS1 + xPartOfRHS2*Z*alpha*Phi1Hat(itheta,izeta))  & 
                              * exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)
                         if (readExternalPhi1) then
                            factorExternalPhi1 =  (alpha*Z/THat) &
                              * (dPhi1Hatdtheta(itheta,izeta)/(GHat + iota * IHat) &
                              - (1/(GHat + iota * IHat)**2) * (GHat*dPhi1Hatdtheta(itheta,izeta) & 
                              - IHat*dPhi1Hatdzeta(itheta,izeta))) 
                            xPartOfRHSExternalPhi1  = xPartOfRHS/x2(ix)
                            stuffToAddExternalPhi1 = scalar * factorExternalPhi1 * xPartOfRHSExternalPhi1
                         end if   
                      else
                         xPartOfRHS = xPartOfRHS1
                      end if
                      
                   case (4) ! iota
                      factor = -IHat/(BHat(itheta,izeta)*(GHat+iota*IHat)**2) &
                           *(GHat*dBHatdtheta(itheta,izeta)-IHat*dBHatdzeta(itheta,izeta))
                      if (includePhi1 .and. includePhi1InKineticEquation .and. localUsePhi1) then
                         xPartOfRHS = (xPartOfRHS1 + xPartOfRHS2*Z*alpha*Phi1Hat(itheta,izeta))  & 
                              * exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)
                         if (readExternalPhi1) then
                            factorExternalPhi1 =  (alpha*Z/THat) &
                              * (- (IHat/(GHat + iota * IHat)**2) * (GHat*dPhi1Hatdtheta(itheta,izeta) & 
                              - IHat*dPhi1Hatdzeta(itheta,izeta)))
                            xPartOfRHSExternalPhi1  = xPartOfRHS/x2(ix)
                            stuffToAddExternalPhi1 = scalar * factorExternalPhi1 * xPartOfRHSExternalPhi1
                         end if   
                      else
                         xPartOfRHS = xPartOfRHS1
                      end if
                   end select
                
                else
                   
                   factor =  1/(BHat(itheta,izeta)**3) &
                        *(BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                        - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                        * DHat(itheta,izeta)
                   if (includePhi1 .and. includePhi1InKineticEquation .and. localUsePhi1) then !!Added by AM 2016-03
                      xPartOfRHS = (xPartOfRHS1 + xPartOfRHS2*Z*alpha*Phi1Hat(itheta,izeta))  & 
                              * exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)
                      
                      if (readExternalPhi1) then !!Added by AM 2018-12
                         factorExternalPhi1 = (alpha*Z/THat) &
                              *1/(BHat(itheta,izeta)**2) & 
                              *(BHat_sub_zeta(itheta,izeta)*dPhi1Hatdtheta(itheta,izeta) & 
                              - BHat_sub_theta(itheta,izeta)*dPhi1Hatdzeta(itheta,izeta))& 
                              * DHat(itheta,izeta) 
                         xPartOfRHSExternalPhi1   = xPartOfRHS/x2(ix)
                         stuffToAddExternalPhi1 = scalar * factorExternalPhi1 * xPartOfRHSExternalPhi1
                      end if 
                   else 
                      xPartOfRHS = xPartOfRHS1
                   end if 
                end if
                stuffToAdd = scalar * factor * xPartOfRHS

                L = 0
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, (four/three)*stuffToAdd + stuffToAddExternalPhi1, INSERT_VALUES, ierr) 

                L = 2
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, (two/three)*stuffToAdd, INSERT_VALUES, ierr)
             end do
          end do
       end do
    end do

    
  end subroutine DKERHS
  
end module residual
