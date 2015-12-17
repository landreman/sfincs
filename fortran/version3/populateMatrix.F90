! For compilers that do not include the error function erf(x), the line
! below should be un-commented:
!#define USE_GSL_ERF
  
#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

  subroutine populateMatrix(matrix, whichMatrix, stateVec)

    use petscmat
    use globalVariables
    use sparsify
    use indices
    use xGrid, only: xGrid_k

    implicit none

    Mat :: matrix
    integer, intent(in) :: whichMatrix
    Vec :: stateVec ! stateVec is ignored when nonlinear=false or when evaluating the residual.

    ! Allowed values for whichMatrix:
    ! 0 = preconditioner for Jacobian
    ! 1 = Jacobian
    ! 2 = matrix which multiplies f0 when evaluating the residual
    ! 3 = matrix which multiplies f1 when evaluating the residual

    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, speciesFactor, speciesFactor2
    PetscScalar :: T32m, factor, LFactor, temp, temp1, temp2, xDotFactor, xDotFactor2, stuffToAdd
    PetscScalar, dimension(:), allocatable :: xb, expxb2
    PetscScalar, dimension(:,:), allocatable :: thetaPartOfTerm, localThetaPartOfTerm
    PetscScalar, dimension(:,:), allocatable :: xPartOfXDot_plus, xPartOfXDot_minus, xPartOfXDot
    PetscScalar, dimension(:,:), allocatable :: ddxToUse_plus, ddxToUse_minus
    integer :: i, j, ix, ispecies, itheta, izeta, L, ixi, index, ix_row, ix_col
    integer :: rowIndex, colIndex
    integer :: ell, iSpeciesA, iSpeciesB, maxL
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar, dimension(:,:), allocatable :: ddxToUse, d2dx2ToUse, zetaPartOfTerm, localZetaPartOfTerm
    PetscScalar, dimension(:,:), allocatable :: fToFInterpolationMatrix
    PetscScalar, dimension(:,:), allocatable :: potentialsToFInterpolationMatrix
    PetscScalar, dimension(:,:,:,:), allocatable :: CECD
    PetscScalar :: xPartOfSource1, xPartOfSource2, geometricFactor1, geometricFactor2, geometricFactor3
    PetscScalar, dimension(:,:), allocatable :: M11, M21, M32, LaplacianTimesX2WithoutL, nuDHat
    PetscScalar, dimension(:), allocatable :: erfs, Psi_Chandra
    PetscScalar, dimension(:,:), allocatable :: CHat, M22, M33, M12, M13
    PetscScalar, dimension(:,:), allocatable :: M22BackslashM21, M33BackslashM32
    PetscScalar, dimension(:,:,:), allocatable :: M22BackslashM21s, M33BackslashM32s
    integer, dimension(:), allocatable :: IPIV  ! Needed by LAPACK
    integer :: LAPACKInfo
    PetscLogDouble :: time1, time2
    PetscScalar, dimension(:,:), allocatable :: ddthetaToUse, ddzetaToUse
    PetscScalar, dimension(:,:), allocatable :: tempMatrix, tempMatrix2, extrapMatrix
    double precision :: myMatInfo(MAT_INFO_SIZE)
    integer :: NNZ, NNZAllocated, NMallocs
    PetscScalar :: CHat_element, dfMdx
    character(len=200) :: whichMatrixName, filename
    PetscViewer :: viewer
    integer :: ithetaRow, ithetaCol, izetaRow, izetaCol, ixMin, ixMinCol
    VecScatter :: vecScatterContext
    Vec :: vecOnEveryProc
    PetscScalar, pointer :: stateArray(:)
    logical :: useStateVec
    PetscScalar, dimension(:,:), allocatable :: nonlinearTerm_Lp1, nonlinearTerm_Lm1
    PetscScalar, dimension(:), allocatable :: tempVector1, tempVector2
    PetscScalar, dimension(:,:), allocatable :: tempExtrapMatrix, fToFInterpolationMatrix_plus1

    ! *******************************************************************************
    ! Do a few sundry initialization tasks:
    ! *******************************************************************************

    ! This next line only matters for nonlinear calculations, in which the Mat objects for the matrix and preconditioner matrix 
    ! are reused at each iteration of SNES. In this case we need to clear all the previous entries, or else we would add the new
    ! values on top of the previous values:
    call MatZeroEntries(matrix,ierr)

    if (masterProc) then
       print *,"Running populateMatrix with whichMatrix = ",whichMatrix
    end if

    select case (whichMatrix)
    case (0)
       ! Preconditioner matrix for the Jacobian matrix
       whichMatrixName = "Jacobian preconditioner"
    case (1)
       ! Jacobian matrix
       whichMatrixName = "Jacobian"
    case (2)
       ! The matrix which is multiplied by f0 when evaluating the residual.
       ! This matrix should only contain the collision operator, and it is used for computing temperature equilibration.
       whichMatrixName = "residual f0"
    case (3)
       ! The matrix which is multiplied by f1 when evaluating the residual.
       ! This matrix is quite similar to the Jacobian matrix (whichMatrix=1), since most terms in the system of equations are linear.
       ! However there are a few differences related to the nonlinear term.
       whichMatrixName = "residual f1"
    case default
       if (masterProc) then
          print *,"Error! whichMatrix must be 0, 1, 2, or 3."
       end if
       stop
    end select

    call PetscTime(time1, ierr)

    ! Sometimes PETSc complains if any of the diagonal elements are not set.
    ! Therefore, set the entire diagonal to 0 to be safe.
    if (masterProc) then
       do i=1,matrixSize
          call MatSetValue(matrix, i-1, i-1, zero, ADD_VALUES, ierr)
       end do
    end if

    ! Since PETSc's direct sparse solver complains if there are any zeros on the diagonal
    ! (unlike mumps or superlu_dist), then if we're using this solver
    ! add some values to the diagonals of the preconditioner.  By trial-and-error, I found it works
    ! best to shift the diagonal of the quasineutrality and constraint blocks but not for the kinetic-equation block.
    if ((.not. isAParallelDirectSolverInstalled) .and. masterProc .and. whichMatrix==0) then
       print *,"Since PETSc's built-in solver is being used instead of superlu_dist or mumps, and this"
       print *,"   fragile solver often gives spurious error messages about zero pivots, the diagonal"
       print *,"   of the preconditioner is being shifted."

       ! Amount of shift:
       temp = 1d+0

       if (constraintScheme==1) then
          do ispecies = 1,Nspecies
             index = getIndex(ispecies,1,1,1,1,BLOCK_DENSITY_CONSTRAINT)
             call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
             index = getIndex(ispecies,1,1,1,1,BLOCK_PRESSURE_CONSTRAINT)
             call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
          end do
       elseif (constraintScheme==2) then
          do ispecies = 1,Nspecies
             do ix = 1,Nx
                index = getIndex(ispecies,ix,1,1,1,BLOCK_F_CONSTRAINT)
                call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
             end do
          end do
       end if
       if (includePhi1) then
          index = getIndex(1,1,1,1,1,BLOCK_PHI1_CONSTRAINT)
          call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
          do itheta = 1,Ntheta
             do izeta = 1,Nzeta
                index = getIndex(1,1,1,itheta,izeta,BLOCK_QN)
                call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
             end do
          end do
       end if
    end if

    useStateVec = (nonlinear .and. (whichMatrix==0 .or. whichMatrix==1))
    if (useStateVec) then
       ! We need delta f to evaluate the Jacobian, so send a copy to every proc:
       call VecScatterCreateToAll(stateVec, vecScatterContext, vecOnEveryProc, ierr)
       call VecScatterBegin(vecScatterContext, stateVec, vecOnEveryProc, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecScatterEnd(vecScatterContext, stateVec, vecOnEveryProc, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecGetArrayF90(vecOnEveryProc, stateArray, ierr)
    end if

    ! In nonlinear runs, the Jacobian and residual require Phi1:
    if (nonlinear .and. (whichMatrix .ne. 2)) then
       call extractPhi1(stateVec)
    end if

    if (pointAtX0) then
       ixMin = 2
    else
       ixMin = 1
    end if

    ! *********************************************************
    ! Allocate small matrices:
    ! *********************************************************

    allocate(ddxToUse(Nx,Nx))
    allocate(ddxToUse_plus(Nx,Nx))
    allocate(ddxToUse_minus(Nx,Nx))
    allocate(d2dx2ToUse(Nx,Nx))
    allocate(ddthetaToUse(Ntheta, Ntheta))
    allocate(ddzetaToUse(Nzeta, Nzeta))

    allocate(xb(Nx))
    allocate(expxb2(Nx))
    allocate(erfs(Nx))
    allocate(Psi_Chandra(Nx))
    allocate(nuDHat(Nspecies, Nx))
    allocate(fToFInterpolationMatrix(Nx,Nx))
    allocate(potentialsToFInterpolationMatrix(Nx, NxPotentials))
    allocate(CECD(Nspecies, Nspecies, Nx, Nx))

    allocate(M21(NxPotentials, Nx))
    allocate(M32(NxPotentials, NxPotentials))
    allocate(M22BackslashM21(NxPotentials, Nx))
    allocate(M33BackslashM32(NxPotentials, NxPotentials))
    allocate(M22BackslashM21s(NL,NxPotentials, Nx))
    allocate(M33BackslashM32s(NL,NxPotentials, NxPotentials))
    allocate(LaplacianTimesX2WithoutL(NxPotentials, NxPotentials))

    allocate(M12(Nx,NxPotentials))
    allocate(M13(Nx,NxPotentials))
    allocate(M22(NxPotentials,NxPotentials))
    allocate(M33(NxPotentials,NxPotentials))

    allocate(M11(Nx,Nx))
    allocate(CHat(Nx,Nx))
    allocate(IPIV(NxPotentials))


    ! ************************************************************
    ! ************************************************************
    ! Begin adding the collisionless terms of the kinetic equation
    ! ************************************************************
    ! ************************************************************

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)

       ! *********************************************************
       ! Add the streaming d/dtheta term:
       ! *********************************************************
       
       if (whichMatrix .ne. 2) then
          allocate(thetaPartOfTerm(Ntheta,Ntheta))
          allocate(localThetaPartOfTerm(Ntheta,localNtheta))
          allocate(rowIndices(localNtheta))
          allocate(colIndices(Ntheta))
          do L=0,(Nxi-1)

             if (whichMatrix>0 .or. L < preconditioner_theta_min_L) then
                ddthetaToUse = ddtheta
             else
                ddthetaToUse = ddtheta_preconditioner
             end if

             do izeta=izetaMin,izetaMax
                do itheta=1,Ntheta
                   thetaPartOfTerm(itheta,:) = BHat_sup_theta(itheta,izeta) &
                        * sqrtTHat/sqrtMHat * ddthetaToUse(itheta,:) &
                        / BHat(itheta,izeta)
                end do
                
                ! PETSc uses the opposite convention to Fortran:
                thetaPartOfTerm = transpose(thetaPartOfTerm)
                localThetaPartOfTerm = thetaPartOfTerm(:,ithetaMin:ithetaMax)
                
                do ix=ixMin,Nx
                   do itheta=1,localNtheta
                      rowIndices(itheta) = getIndex(ispecies, ix, L+1, ithetaMin+itheta-1, izeta, BLOCK_F)
                   end do
                   
                   ! Super-diagonal-in-L term
                   if (L < Nxi-1) then
                      ell = L+1
                      do itheta=1,Ntheta
                         colIndices(itheta) = getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                      end do
                      
                      call MatSetValuesSparse(matrix, localNtheta, rowIndices, Ntheta, colIndices, &
                           (L+1)/(2*L+three)*x(ix)*localThetaPartOfTerm, ADD_VALUES, ierr)
                   end if
                   
                   ! Sub-diagonal-in-L term
                   if (L > 0) then
                      ell = L-1
                      do itheta=1,Ntheta
                         colIndices(itheta) = getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                      end do
                      
                      call MatSetValuesSparse(matrix, localNtheta, rowIndices, Ntheta, colIndices, &
                           L/(2*L-one)*x(ix)*localThetaPartOfTerm, ADD_VALUES, ierr)
                   end if
                   
                   
                end do
             end do
          end do
          deallocate(rowIndices)
          deallocate(colIndices)
          deallocate(thetaPartOfTerm)
          deallocate(localThetaPartOfTerm)
       end if

       ! *********************************************************
       ! Add the streaming d/dzeta term:
       ! *********************************************************

       if (whichMatrix .ne. 2) then
          allocate(zetaPartOfTerm(Nzeta,Nzeta))
          allocate(localZetaPartOfTerm(Nzeta,localNzeta))
          allocate(rowIndices(localNzeta))
          allocate(colIndices(Nzeta))
          do L=0,(Nxi-1)

             if (whichMatrix>0 .or. L < preconditioner_zeta_min_L) then
                ddzetaToUse = ddzeta
             else
                ddzetaToUse = ddzeta_preconditioner
             end if

             do itheta=ithetaMin, ithetaMax
                do izeta=1,Nzeta
                   zetaPartOfTerm(izeta,:) = sqrtTHat/sqrtMHat * BHat_sup_zeta(itheta,izeta) &
                        * ddzetaToUse(izeta,:) / BHat(itheta,izeta)
                end do
                
                ! PETSc uses the opposite convention to Fortran:
                zetaPartOfTerm = transpose(zetaPartOfTerm)
                localZetaPartOfTerm = zetaPartOfTerm(:,izetaMin:izetaMax)
                
                do ix=ixMin,Nx
                   do izeta = 1,localNzeta
                      rowIndices(izeta)=getIndex(ispecies, ix, L+1, itheta, izetaMin+izeta-1, BLOCK_F)
                   end do
                   
                   ! Super-diagonal-in-L term
                   if (L < Nxi-1) then
                      ell = L + 1
                      do izeta = 1,Nzeta
                         colIndices(izeta)=getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                      end do
                      
                      call MatSetValuesSparse(matrix, localNzeta, rowIndices, Nzeta, colIndices, &
                           (L+1)/(2*L+three)*x(ix)*localZetaPartOfTerm, ADD_VALUES, ierr)
                   end if
                   
                   ! Sub-diagonal-in-L term
                   if (L > 0) then
                      ell = L - 1
                      do izeta = 1,Nzeta
                         colIndices(izeta)=getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                      end do
                      
                      call MatSetValuesSparse(matrix, localNzeta, rowIndices, Nzeta, colIndices, &
                           L/(2*L-one)*x(ix)*localZetaPartOfTerm, ADD_VALUES, ierr)
                   end if
                   
                   
                end do
             end do
          end do
          deallocate(rowIndices)
          deallocate(colIndices)
          deallocate(zetaPartOfTerm)
          deallocate(localZetaPartOfTerm)
       end if

       ! *********************************************************
       ! Add the ExB d/dtheta term:
       ! *********************************************************

       if (whichMatrix .ne. 2) then
          factor = alpha*Delta/two*dPhiHatdpsiHat
          allocate(thetaPartOfTerm(Ntheta,Ntheta))
          allocate(localThetaPartOfTerm(Ntheta,localNtheta))
          allocate(rowIndices(localNtheta))
          allocate(colIndices(Ntheta))
          do L=0,(Nxi-1)

             if (ExBDerivativeScheme==0) then
                if (whichMatrix>0 .or. L < preconditioner_theta_min_L) then
                   ddthetaToUse = ddtheta
                else
                   ddthetaToUse = ddtheta_preconditioner
                end if
             else
                ! Assume DHat has the same sign everywhere. (Should be true even for VMEC coordinates.)
                ! Assume BHat_sub_zeta has the same sign everywhere. (True for Boozer but in principle could be violated for VMEC coordinates?)
                if (factor*DHat(1,1)*BHat_sub_zeta(1,1) > 0) then
                   ddthetaToUse = ddtheta_ExB_plus
                else
                   ddthetaToUse = ddtheta_ExB_minus
                end if
             end if

             do izeta=izetaMin,izetaMax
                if (useDKESExBDrift) then
                   do itheta=1,Ntheta
                      thetaPartOfTerm(itheta,:) = ddthetaToUse(itheta,:) / FSABHat2 &
                           * DHat(itheta,izeta) * BHat_sub_zeta(itheta,izeta)
                   end do
                else
                   do itheta=1,Ntheta
                      thetaPartOfTerm(itheta,:) = ddthetaToUse(itheta,:) / (BHat(itheta,izeta) ** 2) &
                           * DHat(itheta,izeta) * BHat_sub_zeta(itheta,izeta)
                   end do
                end if
                
                ! PETSc uses the opposite convention to Fortran:
                thetaPartOfTerm = transpose(thetaPartOfTerm*factor)
                localThetaPartOfTerm = thetaPartOfTerm(:,ithetaMin:ithetaMax)
                
                do ix=ixMin,Nx
                   do itheta=1,localNtheta
                      rowIndices(itheta)=getIndex(ispecies,ix,L+1,itheta+ithetaMin-1,izeta,BLOCK_F)
                   end do
                   do itheta=1,Ntheta
                      colIndices(itheta)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                   end do
                   
                   call MatSetValuesSparse(matrix, localNtheta, rowIndices, Ntheta, colIndices, &
                        localThetaPartOfTerm, ADD_VALUES, ierr)
                end do
             end do
          end do
          deallocate(rowIndices)
          deallocate(colIndices)
          deallocate(thetaPartOfTerm)
          deallocate(localThetaPartOfTerm)
       end if

       ! *********************************************************
       ! Add the ExB d/dzeta term:
       ! *********************************************************

       if (whichMatrix .ne. 2) then
          factor = -alpha*Delta/two*dPhiHatdpsiHat
          allocate(zetaPartOfTerm(Nzeta,Nzeta))
          allocate(localZetaPartOfTerm(Nzeta,localNzeta))
          allocate(rowIndices(localNzeta))
          allocate(colIndices(Nzeta))
          do L=0,(Nxi-1)

             if (ExBDerivativeScheme==0) then
                if (whichMatrix>0 .or. L < preconditioner_zeta_min_L) then
                   ddzetaToUse = ddzeta
                else
                   ddzetaToUse = ddzeta_preconditioner
                end if
             else
                ! Assume DHat has the same sign everywhere. (Should be true even for VMEC coordinates.)
                ! Assume BHat_sub_theta has the same sign everywhere. (True for Boozer but could be violated for VMEC coordinates?)
                if (factor*DHat(1,1)*BHat_sub_theta(1,1) > 0) then
                   ddzetaToUse = ddzeta_ExB_plus
                else
                   ddzetaToUse = ddzeta_ExB_minus
                end if
             end if

             do itheta=ithetaMin, ithetaMax
                if (useDKESExBDrift) then
                   do izeta=1,Nzeta
                      zetaPartOfTerm(izeta,:) = ddzetaToUse(izeta,:) / FSABHat2 &
                           * DHat(itheta,izeta) * BHat_sub_theta(itheta,izeta)
                   end do
                else
                   do izeta=1,Nzeta
                      zetaPartOfTerm(izeta,:) = ddzetaToUse(izeta,:) / (BHat(itheta,izeta) ** 2) &
                           * DHat(itheta,izeta) * BHat_sub_theta(itheta,izeta)
                   end do
                end if
                
                ! PETSc uses the opposite convention to Fortran:
                zetaPartOfTerm = transpose(zetaPartOfTerm*factor)
                localzetaPartOfTerm = zetaPartOfTerm(:,izetaMin:izetaMax)
                
                do ix=ixMin,Nx
                   do izeta=1,localNzeta
                      rowIndices(izeta)=getIndex(ispecies,ix,L+1,itheta,izeta+izetaMin-1,BLOCK_F)
                   end do
                   do izeta=1,Nzeta
                      colIndices(izeta)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                   end do
                   
                   call MatSetValuesSparse(matrix, localNzeta, rowIndices, Nzeta, colIndices, &
                        localZetaPartOfTerm, ADD_VALUES, ierr)
                end do
             end do
          end do
          deallocate(rowIndices)
          deallocate(colIndices)
          deallocate(zetaPartOfTerm)
          deallocate(localZetaPartOfTerm)
       end if

       ! *********************************************************
       ! Add the magnetic drift d/dtheta term:
       ! *********************************************************

       itheta = -1  ! So itheta is not used in place of ithetaRow or ithetaCol by mistake.
       izetaRow = -1 ! So izetaRow is not used in place of izeta by mistake.
       izetaCol = -1 ! So izetaCol is not used in place of izeta by mistake.
       if ((whichMatrix .ne. 2) .and. (magneticDriftScheme>0)) then
          if (whichMatrix==0) then
             maxL = min(preconditioner_magnetic_drifts_max_L,Nxi-1)
          else
             maxL = Nxi-1
          end if
          do L = 0, maxL

             do izeta = izetaMin, izetaMax                
                do ithetaRow = ithetaMin, ithetaMax
                   geometricFactor1 = (BHat_sub_zeta(ithetaRow,izeta)*dBHatdpsiHat(ithetaRow,izeta) &
                        - BHat_sub_psi(ithetaRow,izeta)*dBHatdzeta(ithetaRow,izeta))
                   
                   geometricFactor2 = 2 * BHat(ithetaRow,izeta) &
                        * (dBHat_sub_psi_dzeta(ithetaRow,izeta) - dBHat_sub_zeta_dpsiHat(ithetaRow,izeta))

                   if (magneticDriftScheme==2) then
                      geometricFactor3 = BDotCurlB(ithetaRow,izeta)*BHat_sup_theta(ithetaRow,izeta) &
                           /(BHat(ithetaRow,izeta)*DHat(ithetaRow,izeta))
                   else
                      geometricFactor3 = 0
                   end if

                   if (magneticDriftDerivativeScheme==0) then
                      if (whichMatrix>0 .or. L < preconditioner_theta_min_L) then
                         ddthetaToUse = ddtheta
                      else
                         ddthetaToUse = ddtheta_preconditioner
                      end if
                   else
                      ! Assume DHat has the same sign everywhere. (Should be true even for VMEC coordinates.)
                      ! We assume here that geometricFactor1*factor sets the direction of upwinding to use. This should be correct at beta=0
                      ! since then geometricFactor2=0, but may need modification when beta>0 (in which case geometricFactor2 is nonzero.)
                      if (geometricFactor1*DHat(1,1)/Z > 0) then
                         ddthetaToUse = ddtheta_magneticDrift_plus
                      else
                         ddthetaToUse = ddtheta_magneticDrift_minus
                      end if
                   end if

                   do ix = ixMin, Nx
                      rowIndex = getIndex(ispecies, ix, L+1, ithetaRow, izeta, BLOCK_F)
                      
                      factor = Delta*THat*DHat(ithetaRow,izeta)*x(ix)*x(ix) &
                           / (2*Z*BHat(ithetaRow,izeta)*BHat(ithetaRow,izeta)*BHat(ithetaRow,izeta))
                   
                      do ithetaCol = 1, Ntheta

                         ! Diagonal-in-L term
                         ell = L
                         colIndex = getIndex(ispecies, ix, ell+1, ithetaCol, izeta, BLOCK_F)

                         stuffToAdd = factor * (2*(3*L*L+3*L-2)/((two*L+3)*(2*L-1)) * geometricFactor1 &
                              + (2*L*L+2*L-1)/((two*L+3)*(2*L-1)) * geometricFactor2 &
                              + (-2)*L*(L+1)/((two*L+3)*(2*L-1)) * geometricFactor3)

                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              stuffToAdd * ddthetaToUse(ithetaRow,ithetaCol), &
                              ADD_VALUES, ierr)

                         ! Drop the off-by-2 diagonal terms in L if this is the preconditioner
                         ! and preconditioner_xi = 1:
                         if (whichMatrix .ne. 0 .or. preconditioner_xi==0) then

                            ! Super-super-diagonal-in-L term
                            if (L < Nxi-2) then
                               ell = L+2
                               colIndex = getIndex(ispecies, ix, ell+1, ithetaCol, izeta, BLOCK_F)

                               stuffToAdd = factor*(L+2)*(L+1)/((two*L+5)*(2*L+3)) &
                                    * (geometricFactor1 + geometricFactor2 - 3*geometricFactor3)

                               call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                    stuffToAdd * ddthetaToUse(ithetaRow,ithetaCol), &
                                    ADD_VALUES, ierr)
                            end if

                            ! Sub-sub-diagonal-in-L term
                            if (L > 1) then
                               ell = L-2
                               colIndex = getIndex(ispecies, ix, ell+1, ithetaCol, izeta, BLOCK_F)

                               stuffToAdd = factor*(L-1)*L/((two*L-3)*(2*L-1)) &
                                    * (geometricFactor1 + geometricFactor2 - 3*geometricFactor3)

                               call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                    stuffToAdd * ddthetaToUse(ithetaRow,ithetaCol), &
                                    ADD_VALUES, ierr)
                            end if
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end if

       ! *********************************************************
       ! Add the magnetic drift d/dzeta term:
       ! *********************************************************

       izeta = -1  ! So izeta is not used in place of izetaRow or izetaCol by mistake.
       ithetaRow = -1 ! So ithetaRow is not used in place of itheta by mistake.
       ithetaCol = -1 ! So ithetaCol is not used in place of itheta by mistake.
       if ((whichMatrix .ne. 2) .and. (magneticDriftScheme>0)) then
          if (whichMatrix==0) then
             maxL = min(preconditioner_magnetic_drifts_max_L,Nxi-1)
          else
             maxL = Nxi-1
          end if
          do L = 0, maxL

             do itheta = ithetaMin, ithetaMax                
                do izetaRow = izetaMin, izetaMax
                   geometricFactor1 = (BHat_sub_psi(itheta,izetaRow)*dBHatdtheta(itheta,izetaRow) &
                        - BHat_sub_theta(itheta,izetaRow)*dBHatdpsiHat(itheta,izetaRow))
                   
                   geometricFactor2 = 2 * BHat(itheta,izetaRow) &
                        * (dBHat_sub_theta_dpsiHat(itheta,izetaRow) - dBHat_sub_psi_dtheta(itheta,izetaRow))

                   if (magneticDriftScheme==2) then
                      geometricFactor3 = BDotCurlB(itheta,izetaRow)*BHat_sup_zeta(itheta,izetaRow) &
                           /(BHat(itheta,izetaRow)*DHat(itheta,izetaRow))
                   else
                      geometricFactor3 = 0
                   end if

                   if (magneticDriftDerivativeScheme==0) then
                      if (whichMatrix>0 .or. L < preconditioner_zeta_min_L) then
                         ddzetaToUse = ddzeta
                      else
                         ddzetaToUse = ddzeta_preconditioner
                      end if
                   else
                      ! Assume DHat has the same sign everywhere. (Should be true even for VMEC coordinates.)
                      ! We assume here that geometricFactor1*factor sets the direction of upwinding to use. This should be correct at beta=0
                      ! since then geometricFactor2=0, but may need modification when beta>0 (in which case geometricFactor2 is nonzero.)
                      if (geometricFactor1*DHat(1,1)/Z > 0) then
                         ddzetaToUse = ddzeta_magneticDrift_plus
                      else
                         ddzetaToUse = ddzeta_magneticDrift_minus
                      end if
                   end if

                   do ix = ixMin, Nx
                      rowIndex = getIndex(ispecies, ix, L+1, itheta, izetaRow, BLOCK_F)
                      
                      factor = Delta*THat*DHat(itheta,izetaRow)*x(ix)*x(ix) &
                           / (2*Z*BHat(itheta,izetaRow)*BHat(itheta,izetaRow)*BHat(itheta,izetaRow))
                   
                      do izetaCol = 1, Nzeta

                         ! Diagonal-in-L term
                         ell = L
                         colIndex = getIndex(ispecies, ix, ell+1, itheta, izetaCol, BLOCK_F)

                         stuffToAdd = factor * (2*(3*L*L+3*L-2)/((two*L+3)*(2*L-1)) * geometricFactor1 &
                              + (2*L*L+2*L-1)/((two*L+3)*(2*L-1)) * geometricFactor2 &
                              + (-2)*L*(L+1)/((two*L+3)*(2*L-1)) * geometricFactor3)

                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              stuffToAdd * ddzetaToUse(izetaRow,izetaCol), &
                              ADD_VALUES, ierr)

                         ! Drop the off-by-2 diagonal terms in L if this is the preconditioner
                         ! and preconditioner_xi = 1:
                         if (whichMatrix .ne. 0 .or. preconditioner_xi==0) then

                            ! Super-super-diagonal-in-L term
                            if (L < Nxi-2) then
                               ell = L+2
                               colIndex = getIndex(ispecies, ix, ell+1, itheta, izetaCol, BLOCK_F)

                               stuffToAdd = factor*(L+2)*(L+1)/((two*L+5)*(2*L+3)) &
                                    * (geometricFactor1 + geometricFactor2 - 3*geometricFactor3)

                               call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                    stuffToAdd * ddzetaToUse(izetaRow,izetaCol), &
                                    ADD_VALUES, ierr)
                            end if

                            ! Sub-sub-diagonal-in-L term
                            if (L > 1) then
                               ell = L-2
                               colIndex = getIndex(ispecies, ix, ell+1, itheta, izetaCol, BLOCK_F)

                               stuffToAdd = factor*(L-1)*L/((two*L-3)*(2*L-1)) &
                                    * (geometricFactor1 + geometricFactor2 - 3*geometricFactor3)

                               call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                    stuffToAdd * ddzetaToUse(izetaRow,izetaCol), &
                                    ADD_VALUES, ierr)
                            end if
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end if

       ! *********************************************************
       ! Add the standard mirror term:
       ! *********************************************************

       if (whichMatrix .ne. 2) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax
                factor = -sqrtTHat/(2*sqrtMHat*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                     * (BHat_sup_theta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                     + BHat_sup_zeta(itheta,izeta) * dBHatdzeta(itheta,izeta))
                
                do ix=ixMin,Nx
                   do L=0,(Nxi-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                      
                      if (L<Nxi-1) then
                         ! Super-diagonal-in-L term:
                         ell = L+1
                         colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              (L+1)*(L+2)/(2*L+three)*x(ix)*factor, ADD_VALUES, ierr)
                      end if
                      
                      if (L>0) then
                         ! Sub-diagonal-in-L term:
                         ell = L-1
                         colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              -L*(L-1)/(2*L-one)*x(ix)*factor, ADD_VALUES, ierr)
                      end if
                   end do
                end do
             end do
          end do
       end if

       ! *********************************************************
       ! Add the non-standard d/dxi term associated with E_r:
       ! *********************************************************

       if (includeElectricFieldTermInXiDot .and. (whichMatrix .ne. 2)) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax

                temp = BHat_sub_zeta(itheta,izeta) * dBHatdtheta(itheta,izeta) &
                     - BHat_sub_theta(itheta,izeta) * dBHatdzeta(itheta,izeta)

                if (.not. force0RadialCurrentInEquilibrium) then
                   temp = temp - 2 * BHat(itheta,izeta) &
                        * (dBHat_sub_zeta_dtheta(itheta,izeta) - dBHat_sub_theta_dzeta(itheta,izeta))
                end if

                factor = alpha*Delta*dPhiHatdpsiHat/(4*(BHat(itheta,izeta)**3)) &
                     * DHat(itheta,izeta) * temp
                     
                do ix=ixMin,Nx
                   do L=0,(Nxi-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)

                      ! Diagonal-in-L term
                      call MatSetValueSparse(matrix, rowIndex, rowIndex, &
                           (L+1)*L/((2*L-one)*(2*L+three))*factor, ADD_VALUES, ierr)

                      ! Drop the off-by-2 diagonal terms in L if this is the preconditioner
                      ! and preconditioner_xi = 1:
                      if (whichMatrix .ne. 0 .or. preconditioner_xi==0) then

                         if (L<Nxi-2) then
                            ! Super-super-diagonal-in-L term:
                            ell = L+2
                            colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                 (L+3)*(L+2)*(L+1)/((two*L+5)*(2*L+three))*factor, ADD_VALUES, ierr)
                         end if

                         if (L>1) then
                            ! Sub-sub-diagonal-in-L term:
                            ell = L-2
                            colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                 -L*(L-1)*(L-2)/((2*L-3)*(2*L-one))*factor, ADD_VALUES, ierr)
                         end if
                      end if
                   end do
                end do
             end do
          end do

       end if

       ! ****************************************************************
       ! Add the non-standard d/dxi term associated with magnetic drifts:
       ! ****************************************************************

       if ((magneticDriftScheme>0) .and. (whichMatrix .ne. 2)) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax

                temp = (dBHat_sub_psi_dzeta(itheta,izeta) - dBHat_sub_zeta_dpsiHat(itheta,izeta)) * dBHatdtheta(itheta,izeta) &
                     + (dBHat_sub_theta_dpsiHat(itheta,izeta) - dBHat_sub_psi_dtheta(itheta,izeta)) * dBHatdzeta(itheta,izeta)

                if (.not. force0RadialCurrentInEquilibrium) then
                   temp = temp + (dBHat_sub_zeta_dtheta(itheta,izeta) - dBHat_sub_theta_dzeta(itheta,izeta)) * dBHatdpsiHat(itheta,izeta)
                end if

                do ix=ixMin,Nx
                   factor = -Delta*DHat(itheta,izeta)*THat*x(ix)*x(ix)/(2*Z*(BHat(itheta,izeta)**3)) * temp
                     
                   do L=0,(Nxi-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)

                      ! Diagonal-in-L term
                      call MatSetValueSparse(matrix, rowIndex, rowIndex, &
                           (L+1)*L/((2*L-one)*(2*L+three))*factor, ADD_VALUES, ierr)

                      ! Drop the off-by-2 diagonal terms in L if this is the preconditioner
                      ! and preconditioner_xi = 1:
                      if (whichMatrix .ne. 0 .or. preconditioner_xi==0) then

                         if (L<Nxi-2) then
                            ! Super-super-diagonal-in-L term:
                            ell = L+2
                            colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                 (L+3)*(L+2)*(L+1)/((two*L+5)*(2*L+three))*factor, ADD_VALUES, ierr)
                         end if

                         if (L>1) then
                            ! Sub-sub-diagonal-in-L term:
                            ell = L-2
                            colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                 -L*(L-1)*(L-2)/((2*L-3)*(2*L-one))*factor, ADD_VALUES, ierr)
                         end if
                      end if
                   end do
                end do
             end do
          end do

       end if

       ! *********************************************************
       ! Add the collisionless d/dx term associated with E_r
       ! *********************************************************

       ! This term always operates on f0, so it should always be included when whichMatrix=2 even if includeXDotTerm=.false.!
       ! This term also operates on f1 if includeXDotTerm=.t.
       !if (includeXDotTerm .or. whichMatrix==2) then

       if (includeXDotTerm .and. (whichMatrix .ne. 2)) then

          allocate(xPartOfXDot(Nx,Nx))
          allocate(xPartOfXDot_plus(Nx,Nx))
          allocate(xPartOfXDot_minus(Nx,Nx))
          allocate(rowIndices(Nx))
          allocate(colIndices(Nx))
          !factor = alpha*Delta/(4*psiAHat)*dPhiHatdpsiN
          factor = -alpha*Delta*dPhiHatdpsiHat/4

          do L=0,(Nxi-1)
             if (L>0 .and. pointAtX0) then
                ixMinCol = 2
             else
                ixMinCol = 1
             end if

             ! Upwind in x
             if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                ddxToUse_plus = ddx_xDot_preconditioner_plus
                ddxToUse_minus = ddx_xDot_preconditioner_minus
             else
                ddxToUse_plus = ddx_xDot_plus
                ddxToUse_minus = ddx_xDot_minus
             end if

             do ix=1,Nx
                xPartOfXDot_plus(ix,:) = x(ix) * ddxToUse_plus(ix,:)
                xPartOfXDot_minus(ix,:) = x(ix) * ddxToUse_minus(ix,:)
             end do

             ! Note: in previous versions I take the transpose of xPartOfXDot here,
             ! but since I have switched to using MatSetValueSparse instead of MatSetValuesSparse,
             ! the transpose should no longer be applied here.
             !xPartOfXDot = transpose(xPartOfXDot)  ! PETSc uses the opposite convention of Fortran

             do itheta=ithetaMin,ithetaMax

                do izeta=izetaMin,izetaMax
                   !xDotFactor = factor/(BHat(itheta,izeta)**3) &
                   !     * (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))

                   xDotFactor = factor*DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                        * (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) &
                        - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))

                   if (force0RadialCurrentInEquilibrium) then
                      xDotFactor2 = zero
                   else
                      xDotFactor2 = factor*DHat(itheta,izeta)/(BHat(itheta,izeta)**2) * 2 &
                           * (dBHat_sub_zeta_dtheta(itheta,izeta) - dBHat_sub_theta_dzeta(itheta,izeta))
                   end if

                   if (xDotFactor>0) then  ! This is the correct inequality
                   !if (xDotFactor<0) then  ! Incorrect inequality!
                      xPartOfXDot = xPartOfXDot_plus
                   else
                      xPartOfXDot = xPartOfXDot_minus
                   end if

                   do ix=1,Nx
                      rowIndices(ix)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                   end do

                   ! Term that is diagonal in L:
                   colIndices = rowIndices
                   stuffToAdd = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*xDotFactor &
                        + (2*L*L+2*L-one)/((two*L+3)*(2*L-1))*xDotFactor2
                   do ix_col=ixMinCol,Nx
                      do ix_row=ixMin,Nx
                         call MatSetValueSparse(matrix, rowIndices(ix_row), colIndices(ix_col), &
                              stuffToAdd*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                      end do
                   end do

                   ! Drop the off-by-2 diagonal terms in L if this is the preconditioner
                   ! and preconditioner_xi = 1:
                   if (whichMatrix>0 .or. preconditioner_xi==0) then

                      ! Term that is super-super-diagonal in L:
                      if (L<(Nxi-2)) then
                         ell = L + 2
                         stuffToAdd = (L+1)*(L+2)/((two*L+5)*(2*L+3))*(xDotFactor+xDotFactor2)
                         do ix_col=ixMinCol,Nx
                            colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
                            do ix_row=ixMin,Nx
                               call MatSetValueSparse(matrix, rowIndices(ix_row), colIndex, &
                                    stuffToAdd*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                            end do
                         end do
                      end if

                      ! Term that is sub-sub-diagonal in L:
                      if (L>1) then
                         ell = L - 2
                         stuffToAdd = L*(L-1)/((two*L-3)*(2*L-1))*(xDotFactor+xDotFactor2)
                         do ix_col=ixMinCol,Nx
                            colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
                            do ix_row=ixMin,Nx
                               call MatSetValueSparse(matrix, rowIndices(ix_row), colIndex, &
                                    stuffToAdd*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                            end do
                         end do
                      end if
                   end if

                end do
             end do
          end do
          deallocate(rowIndices)
          deallocate(colIndices)
          deallocate(xPartOfXDot)
          deallocate(xPartOfXDot_plus)
          deallocate(xPartOfXDot_minus)
       end if

       ! *********************************************************
       ! Add the optional term which does not involve derivatives of f,
       ! Sometimes useful for restoring Liouville's theorem:
       ! *********************************************************

!!$       if (include_fDivVE_term) then
!!$          do itheta = ithetaMin,ithetaMax
!!$             do izeta = izetaMin,izetaMax
!!$                do ix = 1,Nx
!!$                   do ixi = 1,Nxi
!!$                      index=getIndex(ispecies,ix,ixi,itheta,izeta,BLOCK_F)
!!$                      call MatSetValueSparse(matrix,index,index, &
!!$                           -dPhiHatdpsiN*Delta*alpha/(psiAHat*(BHat(itheta,izeta)**3)) &
!!$                           *(GHat*dBHatdtheta(itheta,izeta)- IHat*dBHatdzeta(itheta,izeta)), &
!!$                           ADD_VALUES, ierr)
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$       end if


       ! *********************************************************
       ! Add the large Phi_1 terms acting on d f_0 / d x
       ! (These are the terms that give the adiabatic response.)
       ! *********************************************************

       if ((whichMatrix .ne. 2) .and. includePhi1) then
          L=1
          do ix = ixMin,Nx

             dfMdx = -2*x(ix)*expx2(ix)*nHat*mHat*sqrtmHat &
                  /(pi*sqrtpi*THat*sqrtTHat)

             ! d Phi_1 / d theta term:
             itheta = -1 ! Just so a runtime error occurs if I use itheta by mistake
             do ithetaRow = ithetaMin, ithetaMax
                do izeta = izetaMin, izetaMax
                   rowIndex = getIndex(ispecies, ix, L+1, ithetaRow, izeta, BLOCK_F)

                   factor = -alpha*Zs(ispecies)*BHat_sup_theta(ithetaRow,izeta) &
                        /(two*sqrtTHat*sqrtMHat*BHat(ithetaRow,izeta))

                   do ithetaCol = 1,Ntheta
                      colIndex = getIndex(1,1,1,ithetaCol, izeta, BLOCK_QN)
                      call MatSetValueSparse(matrix, rowIndex, colIndex,&
                           dfMdx*factor*ddtheta(ithetaRow, ithetaCol), ADD_VALUES, ierr)
                   end do
                end do
             end do

             ! d Phi_1 / d zeta term:
             izeta = -1 ! Just so a runtime error occurs if I use itheta by mistake
             ithetaRow = -1 ! Just so a runtime error occurs if I use itheta by mistake
             ithetaCol = -1 ! Just so a runtime error occurs if I use itheta by mistake
             do itheta = ithetaMin, ithetaMax
                do izetaRow = izetaMin, izetaMax
                   rowIndex = getIndex(ispecies, ix, L+1, itheta, izetaRow, BLOCK_F)

                   factor = -alpha*Zs(ispecies)*BHat_sup_zeta(itheta,izetaRow) &
                        /(two*sqrtTHat*sqrtMHat*BHat(itheta,izetaRow))

                   do izetaCol = 1,Nzeta
                      colIndex = getIndex(1,1,1,itheta, izetaCol, BLOCK_QN)
                      call MatSetValueSparse(matrix, rowIndex, colIndex,&
                           dfMdx*factor*ddzeta(izetaRow, izetaCol), ADD_VALUES, ierr)
                   end do
                end do
             end do

          end do
       end if

       ! *********************************************************
       ! Add the radial ExB drive term:
       ! (v_E dot grad psi) f_M [(1/n)(dn/dpsi) + (x^2 - 3/2)(1/T)(dT/dpsi)]
       ! This is the new linear term considered by Garcia-Regana et al PPCF (2013).
       ! *********************************************************

       if ((whichMatrix .ne. 2) .and. includeRadialExBDrive) then
          L=0
          do ix = ixMin,Nx

             ! d Phi_1 / d theta term:
             itheta = -1 ! Just so a runtime error occurs if I use itheta by mistake
             do ithetaRow = ithetaMin, ithetaMax
                do izeta = izetaMin, izetaMax
                   rowIndex = getIndex(ispecies, ix, L+1, ithetaRow, izeta, BLOCK_F)

                   factor = -alpha*Delta*DHat(ithetaRow,izeta)*BHat_sub_zeta(ithetaRow,izeta) &
                        /(two*BHat(ithetaRow,izeta)*BHat(ithetaRow,izeta)) &
                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)
                   
                   do ithetaCol = 1,Ntheta
                      colIndex = getIndex(1,1,1,ithetaCol, izeta, BLOCK_QN)
                      call MatSetValueSparse(matrix, rowIndex, colIndex,&
                           factor*ddtheta(ithetaRow, ithetaCol), ADD_VALUES, ierr)
                   end do
                end do
             end do

             ! d Phi_1 / d zeta term:
             izeta = -1 ! Just so a runtime error occurs if I use itheta by mistake
             ithetaRow = -1 ! Just so a runtime error occurs if I use itheta by mistake
             ithetaCol = -1 ! Just so a runtime error occurs if I use itheta by mistake
             do itheta = ithetaMin, ithetaMax
                do izetaRow = izetaMin, izetaMax
                   rowIndex = getIndex(ispecies, ix, L+1, itheta, izetaRow, BLOCK_F)

                   factor = alpha*Delta*DHat(ithetaRow,izeta)*BHat_sub_theta(ithetaRow,izeta) &
                        /(two*BHat(ithetaRow,izeta)*BHat(ithetaRow,izeta)) &
                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)

                   do izetaCol = 1,Nzeta
                      colIndex = getIndex(1,1,1,itheta, izetaCol, BLOCK_QN)
                      call MatSetValueSparse(matrix, rowIndex, colIndex,&
                           factor*ddzeta(izetaRow, izetaCol), ADD_VALUES, ierr)
                   end do
                end do
             end do

          end do
       end if

       ! *********************************************************
       ! Add the nonlinear term in the residual, which also is
       ! the block in the Jacobian from d(kinetic eqn) / d f.
       ! *********************************************************
       ! Note that this term is absent in the first iteration of a nonlinear calculation, because the term is proportional to Phi1.
       ! Therefore, if reusePreconditioner=true, we do not include this term in the preconditioner (which is great because this term introduces a lot of nonzeros.)
       ! If reusePreconditioner=false, then we DO include this term in the preconditioner.
       ! We must use MatSetValue instead of MatSetValueSparse in this term so that we can add some "nonzero" elements whose value is actually 0
       ! in the first iteration (due to Phi1=0). The reason is that PETSc will re-use the pattern of nonzeros from the first iteration in subsequent iterations.
       ! However, we need not add elements which are 0 due to ddtheta=0 as opposed to because Phi1=0, since such elements will remain 0 at every iteration of SNES.

       !if (nonlinear .and. (whichMatrix .ne. 2)) then
       if (nonlinear .and. (whichMatrix == 1 .or. whichMatrix == 3 .or. (whichMatrix==0 .and. .not. reusePreconditioner))) then

          !print *,"@@@@@@ ",myRank," max(abs(dPhi1Hatdtheta)): ",maxval(abs(dPhi1Hatdtheta)),maxval(abs(dPhi1Hatdzeta))

          allocate(nonlinearTerm_Lp1(Nx,Nx))
          allocate(nonlinearTerm_Lm1(Nx,Nx))

          do L=0,(Nxi-1)
             if (L>0 .and. pointAtX0) then
                ixMinCol = 2
             else
                ixMinCol = 1
             end if

             if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                ddxToUse = ddx_preconditioner
             else
                ddxToUse = ddx
             end if

             nonlinearTerm_Lp1 = (L+1)/(two*L+3) * ddxToUse
             nonlinearTerm_Lm1 = L/(two*L-1)     * ddxToUse
             do ix=1,Nx
                nonlinearTerm_Lp1(ix,ix) = nonlinearTerm_Lp1(ix,ix) + (L+1)*(L+2)/(two*L+3)/x(ix)
                nonlinearTerm_Lm1(ix,ix) = nonlinearTerm_Lm1(ix,ix) - (L-1)*L/(two*L-1)/x(ix)
             end do
             do itheta=ithetaMin,ithetaMax

                do izeta=izetaMin,izetaMax
                   factor = -alpha*Zs(ispecies)/(2*BHat(itheta,izeta)*sqrtTHat*sqrtMHat) &
                        * (BHat_sup_theta(itheta,izeta)*dPhi1Hatdtheta(itheta,izeta) &
                        + BHat_sup_zeta(itheta,izeta)*dPhi1Hatdzeta(itheta,izeta))

                   do ix_row=ixMin,Nx
                      rowIndex = getIndex(ispecies,ix_row,L+1,itheta,izeta,BLOCK_F)

                      ! Term that is super-diagonal in L:
                      if (L<(Nxi-1)) then
                         ell = L + 1
                         do ix_col=ixMinCol,Nx
                            if (abs(nonlinearTerm_Lp1(ix_row,ix_col))>threshholdForInclusion) then
                               colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
                               ! We must use MatSetValue instead of MatSetValueSparse here!!
                               call MatSetValue(matrix, rowIndex, colIndex, &
                                    factor*nonlinearTerm_Lp1(ix_row,ix_col), ADD_VALUES, ierr)
                            end if
                         end do
                      end if
                   
                      ! Term that is sub-diagonal in L:
                      if (L>0) then
                         ell = L - 1
                         do ix_col=ixMinCol,Nx
                            if (abs(nonlinearTerm_Lm1(ix_row,ix_col))>threshholdForInclusion) then
                               colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
                               ! We must use MatSetValue instead of MatSetValueSparse here!!
                               call MatSetValue(matrix, rowIndex, colIndex, &
                                    factor*nonlinearTerm_Lm1(ix_row,ix_col), ADD_VALUES, ierr)
                            end if
                         end do
                      end if
                   end do
                end do
             end do
          end do
          deallocate(nonlinearTerm_Lp1)
          deallocate(nonlinearTerm_Lm1)
          
       end if

       ! *********************************************************
       ! Add the block in the Jacobian from d(kinetic eqn) / d Phi1
       ! associated with the nonlinear term.
       ! *********************************************************
       ! Note that this term is absent in the first iteration of a nonlinear calculation, because the term is proportional to delta f.
       ! Therefore, if reusePreconditioner=true, we do not include this term in the preconditioner (which is great because this term introduces a lot of nonzeros.)
       ! If reusePreconditioner=false, then we DO include this term in the preconditioner.
       ! We must use MatSetValue instead of MatSetValueSparse in this term so that we can add some "nonzero" elements whose value is actually 0
       ! in the first iteration (due to delta f = 0). The reason is that PETSc will re-use the pattern of nonzeros from the first iteration in subsequent iterations.
       ! However, we need not add elements which are 0 due to ddtheta=0 as opposed to because delta f = 0, since such elements will remain 0 at every iteration of SNES.

       if (nonlinear .and. (whichMatrix==1 .or. (whichMatrix==0 .and. .not. reusePreconditioner))) then
       !if (nonlinear .and. (whichMatrix==1)) then
       !if (nonlinear .and. (whichMatrix==0 .or. whichMatrix==1)) then
          allocate(tempVector1(Nx))
          allocate(tempVector2(Nx))
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                factor = -alpha*Zs(ispecies)/(2*BHat(itheta,izeta)*sqrtTHat*sqrtMHat)
                do L=0,(Nxi-1)
                   tempVector2=0

                   if (L>0) then
                      ! Add the delta_{L-1, ell} terms:
                      ell = L-1
                      do ix=1,Nx
                         index = getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                         ! Add 1 because we are indexing a Fortran array instead of a PETSc object
                         tempVector1(ix) = stateArray(index+1)
                         tempVector2(ix) = tempVector2(ix) - (L-1)*L/(two*L-1)*stateArray(index+1)/x(ix)
                      end do
                      tempVector2 = tempVector2 + L/(two*L-1) * matmul(ddx,tempVector1)
                   end if

                   if (L<Nxi-1) then
                      ! Add the delta_{L+1, ell} terms:
                      ell = L+1
                      do ix=1,Nx
                         index = getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                         ! Add 1 because we are indexing a Fortran array instead of a PETSc object
                         tempVector1(ix) = stateArray(index+1)
                         tempVector2(ix) = tempVector2(ix) + (L+1)*(L+2)/(two*L+3)*stateArray(index+1)/x(ix)
                      end do
                      tempVector2 = tempVector2 + (L+1)/(two*L+3) * matmul(ddx,tempVector1)
                   end if

                   do ix=ixMin,Nx
                      rowIndex = getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)

                      ! Add the d Phi_1 / d theta term:
                      do j=1,Ntheta
                         if (abs(ddtheta(itheta,j))>threshholdForInclusion) then
                            colIndex = getIndex(1,1,1,j,izeta,BLOCK_QN)
                            ! We must use MatSetValue instead of MatSetValueSparse here!!
                            call MatSetValue(matrix, rowIndex, colIndex, &
                                 factor*BHat_sup_theta(itheta,izeta)*ddtheta(itheta,j)*tempVector2(ix), &
                                 ADD_VALUES, ierr)
                         end if
                      end do

                      ! Add the d Phi_1 / d zeta term:
                      do j=1,Nzeta
                         if (abs(ddzeta(izeta,j))>threshholdForInclusion) then
                            colIndex = getIndex(1,1,1,itheta,j,BLOCK_QN)
                            ! We must use MatSetValue instead of MatSetValueSparse here!!
                            call MatSetValue(matrix, rowIndex, colIndex, &
                                 factor*BHat_sup_zeta(itheta,izeta)*ddzeta(izeta,j)*tempVector2(ix), &
                                 ADD_VALUES, ierr)
                         end if
                      end do
                   end do
                end do
             end do
          end do
          deallocate(tempVector1)
          deallocate(tempVector2)
       end if

    end do ! End of loop over species for the collisionless terms.

    ! *********************************************************
    ! *********************************************************
    !
    ! End of adding the collisionless kinetic terms.
    !
    ! *********************************************************
    ! *********************************************************
    !
    ! Next, we add the collision operator.
    !
    ! *********************************************************
    ! *********************************************************

    ! The collision operator always acts on f1.
    ! The collision operator also acts on f0 if includeTemperatureEquilibrationTerm=.t.

    if (whichMatrix .ne. 2 .or. includeTemperatureEquilibrationTerm) then

       select case (collisionOperator)
          
       case (0)
          ! *********************************************************
          ! Full linearized Fokker-Planck operator
          ! *********************************************************
          
          
          ! *********************************************************
          ! In preparation for adding the collision operator,
          ! create several matrices which will be needed.
          ! *********************************************************
          
          allocate(tempMatrix(Nx, NxPotentials))
          allocate(tempMatrix2(NxPotentials, NxPotentials))
          allocate(extrapMatrix(Nx, NxPotentials))
          
          ! For future possible preconditioners, I might want the change the following 2 lines.
          ddxToUse = ddx
          d2dx2ToUse = d2dx2
          
          ! First assemble rows 2 and 3 of the block linear system, since they
          ! are independent of psi and independent of species.
          
          M32 = zero
          M21 = 4*pi*interpolateXToXPotentials
          do i=2,NxPotentials-1
             M21(i,:) = M21(i,:)*xPotentials(i)*xPotentials(i)
             M32(i,i) = -2*xPotentials(i)*xPotentials(i)
          end do
          M21(1,:)=zero
          M21(NxPotentials,:)=zero
          M32(1,:)=zero
          M32(NxPotentials,:)=zero
          do i=1,NxPotentials
             LaplacianTimesX2WithoutL(i,:) = xPotentials(i)*xPotentials(i)*d2dx2Potentials(i,:) &
                  + 2 * xPotentials(i) * ddxPotentials(i,:)
          end do
          
          do L=0,(NL-1)
             M22 = LaplacianTimesX2WithoutL
             do i=1,NxPotentials
                M22(i,i) = M22(i,i) - L*(L+1)
             end do
             
             ! Add Dirichlet or Neumann boundary condition for potentials at x=0:
             if (L==0) then
                M22(1,:)=ddxPotentials(1,:)
             else
                M22(1,:) = 0
                M22(1,1) = 1
             end if
             M33 = M22;
             
             ! Add Robin boundary condition for potentials at x=xMax:
             M22(NxPotentials,:) = xMaxNotTooSmall*ddxPotentials(NxPotentials,:)
             M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1
             
             ! Boundary condition for G:
             M33(NxPotentials,:) = xMaxNotTooSmall*xMaxNotTooSmall*d2dx2Potentials(NxPotentials,:) &
                  + (2*L+1)*xMaxNotTooSmall*ddxPotentials(NxPotentials,:)
             M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1)
             
             if (L /= 0) then
                M22(NxPotentials,1)=0
                M33(NxPotentials,1)=0
             end if
             
             ! Call LAPACK subroutine DGESV to solve a linear system
             ! Note: this subroutine changes M22 and M33!
             M22BackslashM21 = M21  ! This will be overwritten by LAPACK.
#if defined(PETSC_USE_REAL_SINGLE)
             call SGESV(NxPotentials, Nx, M22, NxPotentials, IPIV, M22BackslashM21, NxPotentials, LAPACKInfo)
#else
             call DGESV(NxPotentials, Nx, M22, NxPotentials, IPIV, M22BackslashM21, NxPotentials, LAPACKInfo)
#endif
             if (LAPACKInfo /= 0) then
                print *, "Error in LAPACK call: info = ", LAPACKInfo
                stop
             end if
             M33BackslashM32 = M32  ! This will be overwritten by LAPACK.
#if defined(PETSC_USE_REAL_SINGLE)
             call SGESV(NxPotentials, NxPotentials, M33, NxPotentials, IPIV, M33BackslashM32, NxPotentials, LAPACKInfo)
#else
             call DGESV(NxPotentials, NxPotentials, M33, NxPotentials, IPIV, M33BackslashM32, NxPotentials, LAPACKInfo)
#endif
             if (LAPACKInfo /= 0) then
                print *, "Error in LAPACK call: info = ", LAPACKInfo
                stop
             end if
             
             M33BackslashM32s(L+1,:,:) = M33BackslashM32
             M22BackslashM21s(L+1,:,:) = M22BackslashM21
          end do
          
          
          nuDHat = zero
          CECD = zero
          ! Before adding the collision operator, we must loop over both species
          ! to build several terms in the operator.
          ! row is species a, column is species b
          do iSpeciesA = 1,Nspecies
             do iSpeciesB = 1,Nspecies
                speciesFactor = sqrt(THats(iSpeciesA)*mHats(iSpeciesB) &
                     / (THats(iSpeciesB) * mHats(iSpeciesA)))
                xb =  x * speciesFactor
                expxb2 = exp(-xb*xb)
                do ix=1,Nx
                   ! erf is vectorized in gfortran but not pathscale
                   temp1 = xb(ix)
#ifdef USE_GSL_ERF
                   call erf(temp1, temp2)
#else
                   temp2 = erf(temp1)
#endif
                   erfs(ix) = temp2
                end do
                Psi_Chandra = (erfs - 2/sqrtpi * xb * expxb2) / (2*xb*xb)
                
                T32m = THats(iSpeciesA) * sqrt(THats(iSpeciesA)*mHats(ispeciesA))
                
                ! Build the pitch-angle scattering frequency:
                nuDHat(iSpeciesA, :) =  nuDHat(iSpeciesA, :) &
                     + (three*sqrtpi/four) / T32m &
                     * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) &
                     * nHats(iSpeciesB)*(erfs - Psi_Chandra)/(x*x*x)
                
                ! Given a vector of function values on the species-B grid, multiply the vector
                ! by this interpolation matrix to obtain its values on the species-A grid:
                if (iSpeciesA /= iSpeciesB) then
                   select case (xGridScheme)
                   case (1,2,5,6)
                      call polynomialInterpolationMatrix(Nx, Nx, x, xb, expx2*(x**xGrid_k), &
                           expxb2*(xb**xGrid_k), fToFInterpolationMatrix)
                   case (3,4)
                      allocate(tempExtrapMatrix(Nx, Nx+1))
                      allocate(fToFInterpolationMatrix_plus1(Nx, Nx+1))
                      call interpolationMatrix(Nx+1, Nx, x_plus1, xb, &
                           xInterpolationScheme, fToFInterpolationMatrix_plus1, tempExtrapMatrix)
                      fToFInterpolationMatrix = fToFInterpolationMatrix_plus1(:,1:Nx)
                      deallocate(tempExtrapMatrix)
                      deallocate(fToFInterpolationMatrix_plus1)
                   case (7)
                      allocate(fToFInterpolationMatrix_plus1(Nx, Nx+1))
                      call ChebyshevInterpolationMatrix(Nx+1, Nx, x_plus1, xb, fToFInterpolationMatrix_plus1)
                      fToFInterpolationMatrix = fToFInterpolationMatrix_plus1(:,1:Nx)
                      deallocate(fToFInterpolationMatrix_plus1)
                   case (8)
                      call ChebyshevInterpolationMatrix(Nx, Nx, x, xb, fToFInterpolationMatrix)
                   case default
                      print *,"Error! Invalid xGridScheme"
                      stop
                   end select
                else
                   fToFInterpolationMatrix = zero
                   do i=1,Nx
                      fToFInterpolationMatrix(i, i) = one
                   end do
                end if
                
!!$                if (masterProc) then
!!$                   print *,"Here comes fToFInterpolationMatrix for ispeciesA=",iSpeciesA,", iSpeciesB=",iSpeciesB
!!$                   do ix=1,Nx
!!$                      print *,fToFInterpolationMatrix(ix,:)
!!$                   end do
!!$                end if

                ! Using the resulting interpolation matrix,
                ! add CD (the part of the field term independent of Rosenbluth potentials.
                ! CD is dense in the species indices.
                
                speciesFactor = 3 * nHats(iSpeciesA)  * mHats(iSpeciesA)/mHats(iSpeciesB) &
                     * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) / T32m
                
                do ix=1,Nx
                   CECD(iSpeciesA, iSpeciesB, ix, :) = CECD(iSpeciesA, iSpeciesB, ix, :) &
                        + speciesFactor * expx2(ix) * fToFInterpolationMatrix(ix, :)
                end do
                
                ! Done adding CD. Now add energy scattering (CE).
                ! Unlike CD, CE is diagonal in the species index.
                
                speciesFactor = 3*sqrtpi/four * nHats(iSpeciesB)  &
                     * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) / T32m
                
                do ix=1,Nx
                   !Now add the d2dx2 and ddx terms in CE:
                   !CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                   CECD(iSpeciesA, iSpeciesA, ix, :) = CECD(iSpeciesA, iSpeciesA, ix, :) &
                        + speciesFactor * (Psi_Chandra(ix)/x(ix)*d2dx2ToUse(ix,:) &
                        + (-2*THats(iSpeciesA)*mHats(iSpeciesB)/(THats(iSpeciesB)*mHats(iSpeciesA)) &
                        * Psi_Chandra(ix)*(1-mHats(iSpeciesA)/mHats(iSpeciesB)) &
                        + (erfs(ix)-Psi_Chandra(ix))/x2(ix)) * ddxToUse(ix,:))
                   
                   ! Lastly, add the part of CE for which f is not differentiated:
                   ! CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                   CECD(iSpeciesA, iSpeciesA, ix, ix) = CECD(iSpeciesA, iSpeciesA, ix, ix) &
                        + speciesFactor *4/sqrtpi*THats(iSpeciesA)/THats(iSpeciesB) &
                        *sqrt(THats(iSpeciesA)*mHats(iSpeciesB)/(THats(iSpeciesB)*mHats(iSpeciesA))) &
                        * expxb2(ix)
                   
                end do
                
             end do
          end do
          
          
          ! *****************************************************************
          ! Now we are ready to add the collision operator to the main matrix.
          ! *****************************************************************
          
          do L=0, Nxi-1
             if (L>0 .and. pointAtX0) then
                ixMinCol = 2
             else
                ixMinCol = 1
             end if

             do iSpeciesB = 1,Nspecies
                do iSpeciesA = 1,Nspecies
                   if (iSpeciesA==iSpeciesB .or. whichMatrix>0 .or. preconditioner_species==0) then
                      
                      speciesFactor = sqrt(THats(iSpeciesA)*mHats(iSpeciesB) &
                           / (THats(iSpeciesB) * mHats(iSpeciesA)))
                      xb =  x * speciesFactor
                      
                      ! Build M11
                      M11 = CECD(iSpeciesA, iSpeciesB,:,:)
                      if (iSpeciesA == iSpeciesB) then
                         do i=1,Nx
                            M11(i,i) = M11(i,i) + (-oneHalf*nuDHat(iSpeciesA,i)*L*(L+1))
                         end do
                      end if
                      
                      !   if (.false.) then
                      if (L < NL) then
                         ! Add Rosenbluth potential terms.

                         if (xGridScheme==5 .or. xGridScheme==6) then
                            ! New scheme for the Rosenbluth potential terms.

                            speciesFactor = 3/(2*pi)*nHats(iSpeciesA) &
                                 * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) &
                                 / (THats(iSpeciesA) * sqrt(THats(iSpeciesA)*mHats(ispeciesA))) &
                                 * THats(iSpeciesB)*mHats(iSpeciesA)/(THats(iSpeciesA)*mHats(iSpeciesB))

                            ! Add terms involving H and d H / d x_b:
                            temp = 1 - mHats(iSpeciesA)/mHats(iSpeciesB)
                            do i=1,Nx
                               M11(i,:) = M11(i,:) - speciesFactor*expx2(i)*( &
                                    Rosenbluth_H(iSpeciesA,iSpeciesB,L+1,i,:) &
                                    + temp * xb(i) * Rosenbluth_dHdxb(iSpeciesA,iSpeciesB,L+1,i,:))
                            end do

                            ! Add term involving d^2 G / d x_b^2:
                            do i=1,Nx
                               M11(i, :) = M11(i,:) + speciesFactor*expx2(i)*x2(i)&
                                    * Rosenbluth_d2Gdxb2(iSpeciesA,iSpeciesB,L+1,i,:)
                            end do

                            CHat = M11

                         else
                            ! Original scheme for the Rosenbluth potential terms.

                            speciesFactor2 = sqrt(THats(iSpeciesA)*mHats(iSpeciesB) &
                                 / (THats(iSpeciesB) * mHats(iSpeciesA)))
                         
                            ! Build M13:
                            call interpolationMatrix(NxPotentials, Nx, xPotentials, x*speciesFactor2, &
                                 xPotentialsInterpolationScheme, potentialsToFInterpolationMatrix, extrapMatrix)
                         
                            speciesFactor = 3/(2*pi)*nHats(iSpeciesA) &
                                 * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) &
                                 / (THats(iSpeciesA) * sqrt(THats(iSpeciesA)*mHats(ispeciesA))) &
                                 * THats(iSpeciesB)*mHats(iSpeciesA)/(THats(iSpeciesA)*mHats(iSpeciesB))
                         
                            tempMatrix = matmul(potentialsToFInterpolationMatrix, d2dx2Potentials)
                            do i=1,Nx
                               !M13(i, :) = speciesFactor*expx2(i)*x2(i)*tempMatrix(i,:)
                               M13(i, :) = speciesFactor*expx2(i) * (x2(i)*tempMatrix(i,:) &
                                    + THats(ispeciesB)*mHats(ispeciesA)/(THats(ispeciesA)*mHats(ispeciesB)) &
                                    *(L+1)*(L+2)*(maxxPotentials ** (L+1)) * (xb(i) ** (-L-1))*extrapMatrix(i,:))
                            end do
                         
                            temp = 1-mHats(iSpeciesA)/mHats(iSpeciesB)
                            do i=1,NxPotentials
                               tempMatrix2(i,:) = temp*xPotentials(i)*ddxPotentials(i,:)
                               tempMatrix2(i,i) = tempMatrix2(i,i) + one
                            end do
                            tempMatrix = matmul(potentialsToFInterpolationMatrix, tempMatrix2)
                            do i=1,Nx
                               !M12(i,:) = -speciesFactor*expx2(i)*tempMatrix(i,:)
                               M12(i,:) = -speciesFactor*expx2(i) * ( tempMatrix(i,:) &
                                    +( -((maxxPotentials/xb(i)) ** (L+1)) &
                                    * ((L+1)*(1-mHats(ispeciesA)/mHats(ispeciesB)) - 1) &
                                    -THats(ispeciesB)*mHats(ispeciesA)/(THats(ispeciesA)*mHats(ispeciesB))&
                                    *((L+1)*(L+2)/(2*L-1) * (maxxPotentials**(L+3))*(xb(i) ** (-L-1)) &
                                    -L*(L-1)/(2*L-1) * (maxxPotentials ** (L+1))*(xb(i)**(-L+1)))) &
                                    *extrapMatrix(i,:))
                            end do
                         
                            ! Possibly add Dirichlet boundary condition for potentials at x=0:
                            if (L /= 0) then
                               M12(:,1) = 0
                               M13(:,1) = 0
                            end if
                         
                            !CHat = M11 -  (M12 - M13 * (M33 \ M32)) * (M22 \ M21);
                            CHat = M11 - matmul(M12 - matmul(M13, M33BackslashM32s(L+1,:,:)),&
                                 M22BackslashM21s(L+1,:,:))
                         end if
                      else
                         CHat = M11;
                      end if
                      
                      if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                         ! We're making the preconditioner, so simplify the x part of the matrix if desired.
                         select case (preconditioner_x)
                         case (0)
                            ! Do nothing.
                         case (1)
                            ! Keep only diagonal in x:
                            do i=1,Nx
                               do j=1,Nx
                                  if (i /= j) then
                                     CHat(i,j) = zero
                                  end if
                               end do
                            end do
                         case (2)
                            ! Keep only upper-triangular part:
                            do i=2,Nx
                               do j=1,(i-1)
                                  CHat(i,j) = zero
                               end do
                            end do
                         case (3,5)
                            ! Keep only tridiagonal part:
                            do i=1,Nx
                               do j=1,Nx
                                  if (abs(i-j)>1) then
                                     CHat(i,j) = zero
                                  end if
                               end do
                            end do
                         case (4)
                            ! Keep only the diagonal and super-diagonal:
                            do i=1,Nx
                               do j=1,Nx
                                  if (i /= j .and. j /= (i+1)) then
                                     CHat(i,j) = zero
                                  end if
                               end do
                            end do
                         case default
                            print *,"Error! Invalid preconditioner_x"
                            stop
                         end select
                         
                      end if
                      
                      ! Note: in previous versions I take the transpose of CHat here,
                      ! but since I have switched to using MatSetValueSparse instead of MatSetValuesSparse,
                      ! the transpose should no longer be applied here.
                      !CHat = transpose(CHat)
                      
                      ! At this point, CHat contains the collision operator normalized by
                      ! \bar{nu}, (the collision frequency at the reference mass, density, and temperature.)
                      
                      do itheta=ithetaMin,ithetaMax
                         do izeta=izetaMin,izetaMax
                            do ix_row=ixMin,Nx
                               rowIndex=getIndex(iSpeciesA,ix_row,L+1,itheta,izeta,BLOCK_F)
                               do ix_col = ixMinCol,Nx
                                  colIndex=getIndex(iSpeciesB,ix_col,L+1,itheta,izeta,BLOCK_F)
                                  call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                       -nu_n*CHat(ix_row,ix_col), ADD_VALUES, ierr)
                               end do
                            end do
                         end do
                      end do
                      
                   end if
                end do
             end do
          end do
          
          deallocate(tempMatrix)
          deallocate(tempMatrix2)
          deallocate(extrapMatrix)
          
          ! *******************************************************************************
          ! *******************************************************************************
          !
          ! Done adding the multi-species Fokker-Planck collision operator.
          !
          ! *******************************************************************************
          ! *******************************************************************************
          
       case (1)
          ! *********************************************************
          ! Pure pitch-angle scattering collision operator
          ! *********************************************************
          
          nuDHat = zero
          ! row is species A, column is species B
          do iSpeciesA = 1,Nspecies
             do iSpeciesB = 1,Nspecies
                speciesFactor = sqrt(THats(iSpeciesA)*mHats(iSpeciesB) &
                     / (THats(iSpeciesB) * mHats(iSpeciesA)))
                xb =  x * speciesFactor
                expxb2 = exp(-xb*xb)
                do ix=1,Nx
                   ! erf is vectorized in gfortran but not pathscale
                   temp1 = xb(ix)
#ifdef USE_GSL_ERF
                   call erf(temp1, temp2)
#else
                   temp2 = erf(temp1)
#endif
                   erfs(ix) = temp2
                end do
                Psi_Chandra = (erfs - 2/sqrtpi * xb * expxb2) / (2*xb*xb)
                
                T32m = THats(iSpeciesA) * sqrt(THats(iSpeciesA)*mHats(ispeciesA))
                
                ! Build the pitch-angle scattering frequency:
                nuDHat(iSpeciesA, :) =  nuDHat(iSpeciesA, :) &
                     + (three*sqrtpi/four) / T32m &
                     * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) &
                     * nHats(iSpeciesB)*(erfs - Psi_Chandra)/(x*x*x)
                
             end do
             
             do L=1, Nxi-1
                do ix=ixMin,Nx
                   CHat_element = -oneHalf*nuDHat(iSpeciesA,ix)*L*(L+1)
                   
                   ! At this point, CHat contains the collision operator normalized by
                   ! \bar{nu}, (the collision frequency at the reference mass, density, and temperature.)
                   
                   do itheta=ithetaMin,ithetaMax
                      do izeta=izetaMin,izetaMax
                         index=getIndex(iSpeciesA,ix,L+1,itheta,izeta,BLOCK_F)
                         call MatSetValueSparse(matrix, index, index, &
                              -nu_n*CHat_element, ADD_VALUES, ierr)
                         !-nu_n*(GHat+iota*IHat)/(BHat(itheta,izeta)*BHat(itheta,izeta))*CHat_element, &
                         !ADD_VALUES, ierr)
                      end do
                   end do
                end do
                
             end do
          end do
          
       case default
          print *,"Error! collisionOperator must be 0 or 1."
          stop
          
       end select
    end if

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Done adding the collision operator.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    ! *******************************************************************************
    ! If there is a grid point at x=0, add the boundary conditions for f at x=0.
    ! *******************************************************************************

    if (pointAtX0 .and. whichMatrix .ne. 2) then
       ! For L > 0 modes, impose f=0 at x=0:
       ix = 1
       do L = 1,(Nxi-1)
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                do ispecies = 1,Nspecies
                   index = getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                   call MatSetValue(matrix, index, index, one, ADD_VALUES, ierr)
                end do
             end do
          end do
       end do

       ! For L=0 mode, impose regularity (df/dx=0) at x=0:
       L=0
       if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
          ddxToUse = ddx_preconditioner
       else
          ddxToUse = ddx
       end if
       ix_row = 1
       do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
             do ispecies = 1,Nspecies
                rowIndex = getIndex(ispecies,ix_row,L+1,itheta,izeta,BLOCK_F)
                do ix_col = 1,Nx
                   colIndex = getIndex(ispecies,ix_col,L+1,itheta,izeta,BLOCK_F)
                   call MatSetValueSparse(matrix, rowIndex, colIndex, ddxToUse(ix_row,ix_col), ADD_VALUES, ierr)
                end do
             end do
          end do
       end do
    end if

    ! *******************************************************************************
    ! Add sources:
    ! *******************************************************************************

    if (whichMatrix .ne. 2) then
       select case (constraintScheme)
       case (0)
          ! Do nothing here.
          
       case (1)
          ! Add a heat source and a particle source.
          
          L=0
          do ix = ixMin,Nx
             xPartOfSource1 = (x2(ix)-5/two)*exp(-x2(ix)) ! Provides particles but no heat
             xPartOfSource2 = (x2(ix)-3/two)*exp(-x2(ix)) ! Provides heat but no particles
             do itheta = ithetaMin,ithetaMax
                do izeta = izetaMin,izetaMax
                   do ispecies = 1,Nspecies
                      rowIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                      
                      colIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource1, ADD_VALUES, ierr)
                      
                      colIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource2, ADD_VALUES, ierr)
                   end do
                end do
             end do
          end do
          
       case (2)
          ! Add a L=0 source (which is constant on the flux surface) at each x.
          L=0
          do ix = ixMin,Nx
             do itheta = ithetaMin,ithetaMax
                do izeta = izetaMin,izetaMax
                   do ispecies = 1,Nspecies
                      rowIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                      colIndex = getIndex(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT)
                      call MatSetValue(matrix, rowIndex, colIndex, one, ADD_VALUES, ierr)
                   end do
                end do
             end do
          end do
          
       case default
          print *,"Error! Invalid constraintScheme."
          stop
       end select
    end if
       
    ! *******************************************************************************
    ! Add the density and pressure constraints:
    ! *******************************************************************************

    if (whichMatrix .ne. 2 .and. procThatHandlesConstraints) then
       select case (constraintScheme)
       case (0)
          ! Do nothing here.

       case (1)
          ! Force the flux-surface-averaged perturbed density and 
          ! flux-surface-averaged perturbed pressure to be zero.

          L=0
          do itheta=1,Ntheta
             do izeta=1,Nzeta
                factor = thetaWeights(itheta)*zetaWeights(izeta)/DHat(itheta,izeta)

                do ix=1,Nx
                   do ispecies=1,Nspecies
                      colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)

                      rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)

                      rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           x2(ix)*x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)
                   end do
                end do
             end do
          end do

       case (2)
          ! Force the flux-surface-averaged distribution function to be zero
          ! at each value of x:

          L=0
          do itheta=1,Ntheta
             do izeta=1,Nzeta
                factor = thetaWeights(itheta)*zetaWeights(izeta)/DHat(itheta,izeta)

                do ix=ixMin,Nx
                   do ispecies = 1,Nspecies
                      colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                      rowIndex = getIndex(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           factor, ADD_VALUES, ierr)
                   end do
                end do
             end do
          end do

          if (pointAtX0) then
             ! If the "normal" <f>=0 constraint is imposed at x=0, it gives a singular matrix.
             ! I think this is because the x=0 constraint can be expressed as a linear combination of the constraints
             ! imposed at other x values. So instead, just set the diagonal to 1 so this row of the matrix does nothing.

             ix = 1
             do ispecies = 1,Nspecies
                rowIndex = getIndex(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT)
                colIndex = rowIndex
                call MatSetValue(matrix, rowIndex, colIndex, &
                     one, ADD_VALUES, ierr)
             end do
          end if

       case default
          print *,"Error! Invalid constraintScheme."
          stop
       end select
    end if

    ! *******************************************************************************
    ! Add the quasineutrality equation
    ! *******************************************************************************

    if (whichMatrix .ne. 2 .and. includePhi1) then
       L=0
       do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
             rowIndex = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN)

             ! Add the charge density of each species:
             do ispecies = 1,Nspecies
                speciesFactor = Zs(ispecies)*THats(ispecies)/mHats(ispecies) &
                     * sqrt(THats(ispecies)/mHats(ispecies))
                do ix = 1,Nx
                   colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                   call MatSetValueSparse(matrix, rowIndex, colIndex, &
                        speciesFactor*x2(ix)*xWeights(ix), ADD_VALUES, ierr)
                end do
             end do

             ! Add the Lagrange multiplier lambda:
             colIndex = getIndex(1, 1, 1, 1, 1, BLOCK_PHI1_CONSTRAINT)
             call MatSetValueSparse(matrix, rowIndex, colIndex, &
                  one, ADD_VALUES, ierr)

          end do
       end do
    end if


    ! *******************************************************************************
    ! Add the constraint < Phi_1 > = 0
    ! *******************************************************************************

    if (whichMatrix .ne. 2 .and. procThatHandlesConstraints .and. includePhi1) then
       allocate(rowIndices(1))
       allocate(colIndices(Nzeta))

       rowIndices(1) = getIndex(1, 1, 1, 1, 1, BLOCK_PHI1_CONSTRAINT)
       do itheta=1,Ntheta
          do izeta=1,Nzeta
             colIndices(izeta) = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN)
          end do

          call MatSetValuesSparse(matrix, 1, rowIndices(1), Nzeta, colIndices, &
               thetaWeights(itheta)*zetaWeights/DHat(itheta,:), ADD_VALUES, ierr)
       end do

       deallocate(rowIndices)
       deallocate(colIndices)
    end if


    ! *******************************************************************************
    ! Done inserting values into the matrices.
    ! Now finalize the matrix
    ! *******************************************************************************

    call PetscTime(time2, ierr)
    if (masterProc) then
       print *,"Time to pre-assemble ",trim(whichMatrixName)," matrix: ", time2-time1, " seconds."
    end if
    call PetscTime(time1, ierr)

    call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

    call PetscTime(time2, ierr)
    if (masterProc) then
       print *,"Time to assemble ",trim(whichMatrixName)," matrix: ", time2-time1, " seconds."
    end if
    call PetscTime(time1, ierr)


    call MatGetInfo(matrix, MAT_GLOBAL_SUM, myMatInfo, ierr)
    NNZ = nint(myMatInfo(MAT_INFO_NZ_USED))
    NNZAllocated = nint(myMatInfo(MAT_INFO_NZ_ALLOCATED))
    NMallocs = nint(myMatInfo(MAT_INFO_MALLOCS))
    if (masterProc) then
       print *,"# of nonzeros in ",trim(whichMatrixName)," matrix:",NNZ, ", allocated:",NNZAllocated, &
            ", mallocs:",NMallocs," (should be 0)"
    end if

    if (saveMatlabOutput) then
       write (filename,fmt="(a,i3.3,a,i1,a)") trim(MatlabOutputFilename) // "_iteration_", iterationForMatrixOutput, &
            "_whichMatrix_",whichMatrix,".m"
       if (masterProc) then
          print *,"Saving matrix in matlab format: ",trim(filename)
       end if
       call PetscViewerASCIIOpen(MPIComm, trim(filename), viewer, ierr)
       call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB, ierr)

       call PetscObjectSetName(matrix, "matrix", ierr)
       call MatView(matrix, viewer, ierr)

       call PetscViewerDestroy(viewer, ierr)
    end if

    if (saveMatricesAndVectorsInBinary) then
       write (filename,fmt="(a,i3.3,a,i1)") trim(binaryOutputFilename) // "_iteration_", iterationForMatrixOutput, &
            "_whichMatrix_",whichMatrix
       if (masterProc) then
          print *,"Saving matrix in binary format: ",trim(filename)
       end if
       call PetscViewerBinaryOpen(MPIComm, trim(filename), FILE_MODE_WRITE, viewer, ierr)
       call MatView(matrix, viewer, ierr)
       call PetscViewerDestroy(viewer, ierr)
    end if

    if (useStateVec) then
       call VecRestoreArrayF90(vecOnEveryProc, stateArray, ierr)
       call VecScatterDestroy(vecScatterContext, ierr)
       call VecDestroy(vecOnEveryProc, ierr)
    end if

  end subroutine populateMatrix

! -----------------------------------------------------------------------
! -----------------------------------------------------------------------


  subroutine init_f0()

    use globalVariables
    use petscmat
    use indices

    implicit none

    integer :: L, ix, itheta, izeta, ispecies, index
    PetscScalar :: factor
    PetscErrorCode :: ierr

    if (masterProc) then
       print *,"Initializing f0"
    end if

    call VecSet(f0, zero, ierr)
    
    L = 0
    do ispecies = 1,Nspecies
       do ix = 1,Nx
          factor = nHats(ispecies)*mHats(ispecies)/(pi*THats(ispecies)) &
               * sqrt(mHats(ispecies)/(pi*THats(ispecies))) * expx2(ix)
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(f0, index, &
                     factor, INSERT_VALUES, ierr)
             end do
          end do
       end do
    end do

    call VecAssemblyBegin(f0, ierr)
    call VecAssemblyEnd(f0, ierr)

  end subroutine init_f0
