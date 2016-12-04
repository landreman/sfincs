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
    Vec :: stateVec ! stateVec is ignored when includePhi1=false or when evaluating the residual.

    ! Allowed values for whichMatrix:
    ! 0 = preconditioner for Jacobian
    ! 1 = Jacobian
    ! 2 = matrix which multiplies f0 when evaluating the residual
    ! 3 = matrix which multiplies f1 when evaluating the residual

    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, speciesFactor, speciesFactor2
    PetscScalar :: T32m, factor, LFactor, temp, temp1, temp2, xDotFactor, xDotFactor2, stuffToAdd
    PetscScalar :: factor1, factor2, factorJ1, factorJ2, factorJ3, factorJ4, factorJ5  !!Added by AM 2016-03
    PetscScalar, dimension(:), allocatable :: xb, expxb2
    PetscScalar, dimension(:,:), allocatable :: xPartOfXDot_plus, xPartOfXDot_minus, xPartOfXDot
    PetscScalar, dimension(:,:), allocatable :: ddx_to_use_plus, ddx_to_use_minus
    integer :: i, j, ix, ispecies, ialpha, izeta, L, ixi, index, ix_row, ix_col, ixi_row, ixi_col
    integer :: rowIndex, colIndex
    integer :: rowIndex2, L2 !!Added by AM 2016-03
    integer :: ell, iSpeciesA, iSpeciesB, maxL
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar, dimension(:,:), allocatable :: ddx_to_use, d2dx2ToUse, zetaPartOfTerm, alphaPartOfTerm
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
    PetscScalar, dimension(:,:), allocatable :: ddalpha_to_use, ddzeta_to_use, ddxi_to_use
    PetscScalar, dimension(:,:), allocatable :: pitch_angle_scattering_operator_to_use
    PetscScalar, dimension(:,:), allocatable :: tempMatrix, tempMatrix2, extrapMatrix
    double precision :: myMatInfo(MAT_INFO_SIZE)
    integer :: NNZ, NNZAllocated, NMallocs
    PetscScalar :: CHat_element, dfMdx
    character(len=200) :: whichMatrixName, filename
    PetscViewer :: viewer
    integer :: ialpha_row, ialpha_col, izeta_row, izeta_col, ixMin, ixMinCol
    VecScatter :: vecScatterContext
    Vec :: vecOnEveryProc
    PetscScalar, pointer :: stateArray(:)
    logical :: useStateVec
    PetscScalar, dimension(:,:), allocatable :: nonlinearTerm_Lp1, nonlinearTerm_Lm1
    PetscScalar, dimension(:), allocatable :: tempVector1, tempVector2
    PetscScalar, dimension(:,:), allocatable :: tempExtrapMatrix, fToFInterpolationMatrix_plus1
    PetscScalar :: dPhiHatdpsiHatToUseInRHS, xPartOfRHS, xPartOfRHS2 !!Added by AM 2016-03
    PetscScalar, dimension(:,:), allocatable :: alpha_interpolation_matrix
    integer :: zeta_shift, zeta_pad_size
    integer :: stencil, sign, izeta_interpolation

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
       ! However there are a few differences related to the nonlinear terms.
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

       select case (constraintScheme)
       case (0)
          ! Nothing to do here.
       case (1,3,4)
          do ispecies = 1,Nspecies
             index = getIndex(ispecies,1,1,1,1,BLOCK_DENSITY_CONSTRAINT)
             call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
             index = getIndex(ispecies,1,1,1,1,BLOCK_PRESSURE_CONSTRAINT)
             call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
          end do
       case (2)
          do ispecies = 1,Nspecies
             do ix = 1,Nx
                index = getIndex(ispecies,ix,1,1,1,BLOCK_F_CONSTRAINT)
                call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
             end do
          end do
       case default
          stop "Invalid constraintScheme!"
       end select
       if (includePhi1) then
          index = getIndex(1,1,1,1,1,BLOCK_PHI1_CONSTRAINT)
          call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
          do ialpha = 1,Nalpha
             do izeta = 1,Nzeta
                index = getIndex(1,1,1,ialpha,izeta,BLOCK_QN)
                call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
             end do
          end do
       end if
    end if

    useStateVec = (includePhi1 .and. (whichMatrix==0 .or. whichMatrix==1)) !!Added by AM 2016-02
    if (useStateVec) then
       ! We need delta f to evaluate the Jacobian, so send a copy to every proc:
       call VecScatterCreateToAll(stateVec, vecScatterContext, vecOnEveryProc, ierr)
       call VecScatterBegin(vecScatterContext, stateVec, vecOnEveryProc, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecScatterEnd(vecScatterContext, stateVec, vecOnEveryProc, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecGetArrayF90(vecOnEveryProc, stateArray, ierr)
    end if

    ! In nonlinear runs, the Jacobian and residual require Phi1:
    if (includePhi1 .and. (whichMatrix .ne. 2)) then
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

    allocate(d2dx2ToUse(Nx,Nx))

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
    ! Add rows for interpolation in alpha
    ! ************************************************************
    ! ************************************************************

    if (masterProc) print *,"Beginning alpha interpolation"
    if (whichMatrix .ne. 2) then
       if (whichMatrix==0) then
          stencil = preconditioner_alpha_interpolation_stencil
       else
          stencil = alpha_interpolation_stencil
       end if

       zeta_shift = Nzeta - 2*buffer_zeta_points_on_each_side
       
       allocate(alpha_interpolation_matrix(Nalpha,Nalpha))
       do izeta = izetaMin,izetaMax
          if (izeta <= buffer_zeta_points_on_each_side) then
             sign = -1
          elseif (izeta > Nzeta - buffer_zeta_points_on_each_side) then
             sign = 1
          else
             !print *,"skipping izeta=",izeta
             cycle
          end if
          
          call periodicInterpolation(Nalpha, alpha_interpolation_matrix, sign*iota*zetaMax, stencil)
          alpha_interpolation_matrix = -alpha_interpolation_matrix

          izeta_interpolation = izeta-sign*zeta_shift
          !print *,"izeta=",izeta,", izeta_interpolation=",izeta_interpolation
          do ix=ixMin,Nx
             do ixi=1,Nxi_for_x(ix)
                do ialpha_row = ialphaMin,ialphaMax
                   do ispecies = 1,Nspecies
                      rowIndex = getIndex(ispecies, ix, ixi, ialpha_row, izeta, BLOCK_F)
                      ! Add identity along the diagonal:
                      call MatSetValue(matrix,rowIndex,rowIndex,one,ADD_VALUES,ierr)
                      
                      do ialpha_col = 1,Nalpha
                         ! Add interpolation:
                         colIndex = getIndex(ispecies, ix, ixi, ialpha_col, izeta_interpolation, BLOCK_F)
                         call MatSetValueSparse(matrix,rowIndex,colIndex, &
                              alpha_interpolation_matrix(ialpha_row,ialpha_col),ADD_VALUES,ierr)
                      end do
                   end do
                end do
             end do
          end do
       end do
       deallocate(alpha_interpolation_matrix)
    end if
    if (masterProc) print *,"Done with alpha interpolation"
   
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
       ! Add the d/dzeta term:
       ! *********************************************************
       
       if (masterProc) print *,"Beginning d/dzeta"
       allocate(ddzeta_to_use(Nzeta,Nzeta))
       if ((whichMatrix .ne. 2) .and. (Nzeta > 1)) then
          do izeta_row = izetaMinDKE,izetaMaxDKE
             do ialpha = ialphaMin,ialphaMax
                do ix = 1,Nx
                   do ixi = 1,Nxi_for_x(ix)
                      factor = x(ix)*xi(ixi)
                      select case (ExB_option)
                      case (1)
                         factor = factor - (gamma*Delta*sqrtMHat/(2*sqrtTHat))*BHat_sub_theta(ialpha,izeta_row)/BHat(ialpha,izeta_row)*dPhiHatdpsiHat
                      case (2)
                         factor = factor - (gamma*Delta*sqrtMHat/(2*sqrtTHat))*BHat_sub_theta(ialpha,izeta_row)*BHat(ialpha,izeta_row)/FSABHat2*dPhiHatdpsiHat
                      end select

                      if (whichMatrix>0) then
                         ! Not preconditioner
                         if (factor>0) then
                            ddzeta_to_use = ddzeta_plus
                         else
                            ddzeta_to_use = ddzeta_minus
                         end if
                      else
                         ! Preconditioner
                         if (factor>0) then
                            ddzeta_to_use = ddzeta_plus_preconditioner
                         else
                            ddzeta_to_use = ddzeta_minus_preconditioner
                         end if
                      end if
                      rowIndex = getIndex(ispecies, ix, ixi, ialpha, izeta_row, BLOCK_F)
                      do izeta_col = 1,Nzeta
                         colIndex = getIndex(ispecies, ix, ixi, ialpha, izeta_col, BLOCK_F)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              factor * ddzeta_to_use(izeta_row,izeta_col), ADD_VALUES, ierr)
                      end do
                   end do
                end do
             end do
          end do
       end if
       ! To generate an error if these variables are used by mistake later:
       izeta_row = -1
       izeta_col = -1
       deallocate(ddzeta_to_use)
       if (masterProc) print *,"Done with d/dzeta"

       
       ! *********************************************************
       ! Add the d/dalpha term:
       ! *********************************************************

       allocate(ddalpha_to_use(Nalpha,Nalpha))
       if (whichMatrix .ne. 2) then
          do izeta = izetaMinDKE,izetaMaxDKE
             do ialpha_row = ialphaMin,ialphaMax
                do ix = 1,Nx
                   do ixi = 1,Nxi_for_x(ix)
                      if (Nzeta > 1) then
                         ! We are in a stellarator
                         select case (ExB_option)
                         case (1)
                            factor = gamma*Delta*sqrtMHat/(2*sqrtTHat)*BHat(ialpha_row,izeta)/DHat(ialpha_row,izeta)*dPhiHatdpsiHat
                         case (2)
                            factor = gamma*Delta*sqrtMHat/(2*sqrtTHat)*(BHat(ialpha_row,izeta) ** 3)/(DHat(ialpha_row,izeta)*FSABHat2)*dPhiHatdpsiHat
                         end select
                      else
                         ! Nzeta==1, so we are in a tokamak.
                         factor = iota * x(ix) * xi(ixi)
                         select case (ExB_option)
                         case (1)
                            factor = factor + gamma*Delta*sqrtMHat*dPhiHatdpsiHat*BHat_sub_zeta(ialpha,izeta)/(2*BHat(ialpha,izeta)*sqrtTHat)
                         case (2)
                            factor = factor + gamma*Delta*sqrtMHat*dPhiHatdpsiHat*BHat_sub_zeta(ialpha,izeta)*BHat(ialpha,izeta)/(2*FSABHat2*sqrtTHat)
                         end select
                      end if

                      if (whichMatrix>0) then
                         ! Not preconditioner
                         if (factor>0) then
                            ddalpha_to_use = ddalpha_plus
                         else
                            ddalpha_to_use = ddalpha_minus
                         end if
                      else
                         ! Preconditioner
                         if (factor>0) then
                            ddalpha_to_use = ddalpha_plus_preconditioner
                         else
                            ddalpha_to_use = ddalpha_minus_preconditioner
                         end if
                      end if

                      rowIndex = getIndex(ispecies, ix, ixi, ialpha_row, izeta, BLOCK_F)
                      do ialpha_col = 1,Nalpha
                         colIndex = getIndex(ispecies, ix, ixi, ialpha_col, izeta, BLOCK_F)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              factor * ddalpha_to_use(ialpha_row,ialpha_col), ADD_VALUES, ierr)
                      end do
                   end do
                end do
             end do
          end do
       end if
       ! To generate an error if these variables are used by mistake later:
       ialpha_row = -1
       ialpha_col = -1
       deallocate(ddalpha_to_use)


       ! *********************************************************
       ! Add the d/dxi term:
       ! *********************************************************

       allocate(ddxi_to_use(Nxi,Nxi))
       if (whichMatrix .ne. 2) then
          do izeta = izetaMinDKE,izetaMaxDKE
             do ialpha = ialphaMin,ialphaMax
                do ix = 1,Nx
                   do ixi_row = 1,Nxi_for_x(ix)
                      ! The next line contains \dot{\hat{\xi}}:
                      factor = - x(ix)*(1-xi(ixi_row)*xi(ixi_row))/(2*BHat(ialpha,izeta)) * (dBHatdzeta(ialpha,izeta)+iota*dBHatdtheta(ialpha,izeta))
                      if (includeElectricFieldTermInXiDot) then
                         factor = factor + xi(ixi_row)*(1-xi(ixi_row)*xi(ixi_row))*gamma*Delta*sqrtMHat/(4*sqrtTHat*BHat(ialpha,izeta)*BHat(ialpha,izeta)) &
                              * (BHat_sub_zeta(ialpha,izeta)*dBHatdtheta(ialpha,izeta) - BHat_sub_theta(ialpha,izeta)*dBHatdzeta(ialpha,izeta))*dPhiHatdpsiHat
                      end if

                      if (whichMatrix>0) then
                         ! Not preconditioner
                         if (factor>0) then
                            ddxi_to_use = ddxi_plus
                         else
                            ddxi_to_use = ddxi_minus
                         end if
                      else
                         ! Preconditioner
                         if (factor>0) then
                            ddxi_to_use = ddxi_plus_preconditioner
                         else
                            ddxi_to_use = ddxi_minus_preconditioner
                         end if
                      end if
                      rowIndex = getIndex(ispecies, ix, ixi_row, ialpha, izeta, BLOCK_F)
                      do ixi_col = 1,Nxi_for_x(ix)
                         colIndex = getIndex(ispecies, ix, ixi_col, ialpha, izeta, BLOCK_F)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              factor * ddxi_to_use(ixi_row,ixi_col), ADD_VALUES, ierr)
                      end do
                   end do
                end do
             end do
          end do
       end if
       ! To generate an error if these variables are used by mistake later:
       ixi_row = -1
       ixi_col = -1
       deallocate(ddxi_to_use)


       ! *********************************************************
       ! Add the d/dx term:
       ! *********************************************************

       allocate(ddx_to_use(Nx,Nx))
       if ((whichMatrix .ne. 2) .and. includeXDotTerm) then
          do izeta = izetaMinDKE,izetaMaxDKE
             do ialpha = ialphaMin,ialphaMax
                do ix_row = 1,Nx
                   do ixi = 1,Nxi
                      ! The next line contains \dot{\hat{x}}:
                      factor = -gamma*Delta*sqrtMHat*dPhiHatdpsiHat/(4*sqrtTHat*BHat(ialpha,izeta)*BHat(ialpha,izeta))
                      if (whichMatrix>0) then
                         ! Not preconditioner
                         ddx_to_use = ddx
                      else
                         ! Preconditioner
                         ddx_to_use = ddx_preconditioner
                      end if
                      rowIndex = getIndex(ispecies, ix_row, ixi, ialpha, izeta, BLOCK_F)
                      do ix_col = 1,Nx
                         colIndex = getIndex(ispecies, ix_col, ixi, ialpha, izeta, BLOCK_F)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              factor * ddx_to_use(ix_row,ix_col), ADD_VALUES, ierr)
                      end do
                   end do
                end do
             end do
          end do
       end if
       ! To generate an error if these variables are used by mistake later:
       ix_row = -1
       ix_col = -1
       deallocate(ddx_to_use)





!!$       ! *********************************************************
!!$       ! SECTION MODIFIED BY AM 2016-02/03
!!$       ! Before this section added the radial ExB drive term:
!!$       ! (v_E dot grad psi) f_M [(1/n)(dn/dpsi) + (x^2 - 3/2)(1/T)(dT/dpsi)],
!!$       ! considered by Garcia-Regana et al PPCF (2013),
!!$       ! but now it is multiplied by exp(-Ze Phi1 / T) according to 
!!$       ! Garcia-Regana et al arxiv:1501.03967 (2015).
!!$       ! We also add the rest of the terms in the drift-kinetic equation 
!!$       ! that contain a linear operator acting on Phi1 (but not f1).
!!$       ! These terms are the same in the residual and Jacobian.
!!$       ! We also add additional terms for the Jacobian.
!!$       ! *********************************************************
!!$
!!$       !!THE includeRadialExBDrive FLAG IS NOT USED ANYMORE
!!$       !!if ((whichMatrix .ne. 2) .and. includeRadialExBDrive) then !!Commented by AM 2016-02
!!$       if ((whichMatrix .ne. 2) .and. includePhi1 .and. includePhi1InKineticEquation) then !!Added by AM 2016-02
!!$          L=0
!!$          L2=2 !!Added by AM 2016-03 to add extra P_2 terms
!!$          
!!$          !!Added by AM 2016-03!!
!!$          !!The naming "RHS" here is a remnant from that term was first introduced in evaluateResidual.F90 in an earlier version of SFINCS
!!$          if (RHSMode==1) then
!!$             dPhiHatdpsiHatToUseInRHS = dPhiHatdpsiHat
!!$          else
!!$             dPhiHatdpsiHatToUseInRHS = 0
!!$          end if
!!$          !!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$          do ix = ixMin,Nx
!!$
!!$             !!Added by AM 2016-03!!
!!$             !!The naming "RHS" here is a remnant from that term was first introduced in evaluateResidual.F90 in an earlier version of SFINCS
!!$             xPartOfRHS = x2(ix)*expx2(ix)*( dnHatdpsiHats(ispecies)/nHat &
!!$                  + gamma*Z/THat*dPhiHatdpsiHatToUseInRHS &
!!$                  + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THat)
!!$
!!$             xPartOfRHS2 = x2(ix)*expx2(ix)*dTHatdpsiHats(ispecies)/(THat*THat)
!!$             !!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$             ! d Phi_1 / d theta term:
!!$             ialpha = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$             do ialpha_row = ialphaMin, ialphaMax
!!$                do izeta = izetaMinDKE, izetaMaxDKE
!!$                   rowIndex = getIndex(ispecies, ix, L+1, ialpha_row, izeta, BLOCK_F)
!!$
!!$! COMMENTED BY AM 2016-02
!!$!                   factor = -gamma*Delta*DHat(ialpha_row,izeta)*BHat_sub_zeta(ialpha_row,izeta) &
!!$!                        /(two*BHat(ialpha_row,izeta)*BHat(ialpha_row,izeta)) &
!!$!                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$!                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)
!!$
!!$                   !!Added by AM 2016-03!!
!!$                   factor1 = -exp(-Z*gamma*Phi1Hat(ialpha_row,izeta)/THat)*gamma*Delta*DHat(ialpha_row,izeta)*BHat_sub_zeta(ialpha_row,izeta) &
!!$                        /(two*BHat(ialpha_row,izeta)*BHat(ialpha_row,izeta)) &
!!$                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)
!!$
!!$                   factor2 = - gamma*gamma*Delta*DHat(ialpha_row,izeta)/(two*pi*sqrtpi*BHat(ialpha_row,izeta)*BHat(ialpha_row,izeta)) &
!!$                        * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
!!$                        * exp(-Z*gamma*Phi1Hat(ialpha_row,izeta)/THat)*expx2(ix)*BHat_sub_zeta(ialpha_row,izeta) &
!!$                        * (dPhiHatdpsiHat + Phi1Hat(ialpha_row,izeta)*dTHatdpsiHats(ispecies)/THat)
!!$
!!$                   do ialpha_col = 1,Nalpha
!!$                      colIndex = getIndex(1,1,1,ialpha_col, izeta, BLOCK_QN)
!!$!!                      call MatSetValueSparse(matrix, rowIndex, colIndex,& !!Commented by AM 2016-03
!!$!!                           factor*ddalpha(ialpha_row, ialpha_col), ADD_VALUES, ierr) !!Commented by AM 2016-03
!!$                      !!! SHOULD THE NEXT LINE BE MatSetValueSparse???
!!$                      !!NEXT CALL TEMPORARILY COMMENTED BY AM!!
!!$!                      call MatSetValue(matrix, rowIndex, colIndex,& !!Added by AM 2016-03
!!$!                           (factor1 + factor2)*ddalpha(ialpha_row, ialpha_col), ADD_VALUES, ierr) !!Added by AM 2016-03
!!$
!!$                      !!TEST BY AM 2016-04-11!!
!!$                      call MatSetValueSparse(matrix, rowIndex, colIndex,& !!Added by AM 2016-03
!!$                           (factor1 + factor2)*ddalpha(ialpha_row, ialpha_col), ADD_VALUES, ierr) !!Added by AM 2016-03
!!$                      !!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$                   end do
!!$                end do
!!$             end do
!!$
!!$             ! d Phi_1 / d zeta term:
!!$             izeta = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$             ialpha_row = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$             ialpha_col = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$             do ialpha = ialphaMin, ialphaMax
!!$                do izeta_row = izetaMinDKE, izetaMaxDKE
!!$                   rowIndex = getIndex(ispecies, ix, L+1, ialpha, izeta_row, BLOCK_F)
!!$
!!$! COMMENTED BY AM 2016-02
!!$!                   factor = gamma*Delta*DHat(ialpha_row,izeta)*BHat_sub_theta(ialpha_row,izeta) &
!!$!                        /(two*BHat(ialpha_row,izeta)*BHat(ialpha_row,izeta)) &
!!$!                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$!                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)
!!$                   
!!$                   !! ADDED BY AM 2016-02!!
!!$                   factor1 = exp(-Z*gamma*Phi1Hat(ialpha,izeta_row)/THat)*gamma*Delta*DHat(ialpha,izeta_row)*BHat_sub_theta(ialpha,izeta_row) &
!!$                        /(two*BHat(ialpha,izeta_row)*BHat(ialpha,izeta_row)) &
!!$                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)
!!$                   
!!$                   factor2 = gamma*gamma*Delta*DHat(ialpha,izeta_row)/(two*pi*sqrtpi*BHat(ialpha,izeta_row)*BHat(ialpha,izeta_row)) &
!!$                        * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
!!$                        * exp(-Z*gamma*Phi1Hat(ialpha,izeta_row)/THat)*expx2(ix)*BHat_sub_theta(ialpha,izeta_row) &
!!$                        * (dPhiHatdpsiHat + Phi1Hat(ialpha,izeta_row)*dTHatdpsiHats(ispecies)/THat)
!!$                   
!!$                   do izeta_col = 1,Nzeta
!!$                      colIndex = getIndex(1,1,1,ialpha, izeta_col, BLOCK_QN)
!!$                      !!                      call MatSetValueSparse(matrix, rowIndex, colIndex,& !!Commented by AM 2016-03
!!$                      !!                           factor*ddzeta(izeta_row, izeta_col), ADD_VALUES, ierr) !!Commented by AM 2016-03
!!$                      !!! SHOULD THE NEXT LINE BE MatSetValueSparse ????
!!$                      !!NEXT CALL TEMPORARILY COMMENTED BY AM!!
!!$!                      call MatSetValue(matrix, rowIndex, colIndex,& !!Added by AM 2016-03
!!$!                           (factor1 + factor2)*ddzeta(izeta_row, izeta_col), ADD_VALUES, ierr) !!Added by AM 2016-03
!!$
!!$                      !!TEST BY AM 2016-04-11!!
!!$                      call MatSetValueSparse(matrix, rowIndex, colIndex,& !!Added by AM 2016-03
!!$                           (factor1 + factor2)*ddzeta(izeta_row, izeta_col), ADD_VALUES, ierr) !!Added by AM 2016-03
!!$                      !!!!!!!!!!!!!!!!!!!!!!!!!
!!$                      
!!$                   end do
!!$                end do
!!$             end do
!!$             !!!!!!!!!!!!!!!!!!!!!!!!
!!$             !!Added by AM 2016-02!!
!!$             !!Add additional terms in the Jacobian
!!$             if (whichMatrix == 0 .or. whichMatrix == 1) then
!!$                ialpha = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$                izeta = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$                ialpha_row = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$                ialpha_col = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$                izeta_row = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$                izeta_col = -1 ! Just so a runtime error occurs if I use ialpha by mistake
!!$                do ialpha = ialphaMin, ialphaMax
!!$                   do izeta = izetaMinDKE, izetaMaxDKE
!!$                      rowIndex = getIndex(ispecies, ix, L+1, ialpha, izeta, BLOCK_F)
!!$                      rowIndex2 = getIndex(ispecies, ix, L2+1, ialpha, izeta, BLOCK_F) !!Added by AM 2016-03 to add extra P_2 terms
!!$                      
!!$                      !!These terms are diagonal in theta and zeta
!!$                      colIndex = getIndex(1,1,1,ialpha, izeta, BLOCK_QN)
!!$                      
!!$                      !!The following factors are only used in the Jacobian. These terms do not contain d/dtheta or d/dzeta
!!$                      !!but d(Phi1)/dtheta, d(Phi1)/dzeta, dB/dtheta, dB/dzeta.
!!$                      !!Therefore they are the same when adding the d/dtheta terms as when adding the d/dzeta terms.
!!$                      
!!$                      !!factorJ1 stems from factor1 but contains both dphi1hatdalpha and dPhi1Hatdzeta parts
!!$                      factorJ1 = exp(-Z*gamma*Phi1Hat(ialpha,izeta)/THat)*gamma*Delta*DHat(ialpha,izeta) &
!!$                           *(- BHat_sub_zeta(ialpha,izeta)* dPhi1Hatdalpha(ialpha, izeta) + BHat_sub_theta(ialpha, izeta)*dPhi1Hatdzeta(ialpha, izeta))&
!!$                           /(two*BHat(ialpha,izeta)*BHat(ialpha,izeta)) &
!!$                           * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$                           * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)  
!!$                      
!!$                      !!factorJ2 stems from factor2 but contains both dphi1hatdalpha and dPhi1Hatdzeta parts
!!$                      factorJ2 = gamma*gamma*Delta*DHat(ialpha,izeta)/(two*pi*sqrtpi*BHat(ialpha,izeta)*BHat(ialpha,izeta)) &
!!$                           * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
!!$                           * exp(-Z*gamma*Phi1Hat(ialpha,izeta)/THat)*expx2(ix) &
!!$                           * (- BHat_sub_zeta(ialpha,izeta)*dphi1hatdalpha(ialpha, izeta) + BHat_sub_theta(ialpha, izeta)*dPhi1Hatdzeta(ialpha, izeta)) &
!!$                           * (dPhiHatdpsiHat + Phi1Hat(ialpha,izeta)*dTHatdpsiHats(ispecies)/THat)
!!$                      
!!$                      !!factorJ3 is only used in the Jacobian, it corresponds to a part \propto exp(-Ze Phi1/T) which is
!!$                      !!implemented in evaluateResidual.F90 for the residual
!!$                      factorJ3 = Delta*nHat*mHat*sqrtMHat &
!!$                           /(2*pi*sqrtpi*Z*(BHat(ialpha,izeta)**3)*sqrtTHat) &
!!$                           *(- BHat_sub_zeta(ialpha,izeta)*dBHatdtheta(ialpha,izeta) + BHat_sub_theta(ialpha,izeta)*dBHatdzeta(ialpha,izeta))&
!!$                           * DHat(ialpha,izeta) * (xPartOfRHS + xPartOfRHS2*Z*gamma*Phi1Hat(ialpha,izeta))  & 
!!$                           * exp(-Z*gamma*Phi1Hat(ialpha,izeta)/THat)
!!$                      
!!$                      !!factorJ4 is only used in the Jacobian, it corresponds to the second term in factorJ2 divided by Phi1
!!$                      factorJ4 =  gamma*gamma*Delta*DHat(ialpha,izeta)/(two*pi*sqrtpi*BHat(ialpha,izeta)*BHat(ialpha,izeta)) &
!!$                           * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
!!$                           * exp(-Z*gamma*Phi1Hat(ialpha,izeta)/THat)*expx2(ix) &
!!$                           *(- BHat_sub_zeta(ialpha,izeta)*dphi1hatdalpha(ialpha,izeta) + BHat_sub_theta(ialpha,izeta)*dPhi1Hatdzeta(ialpha,izeta))&
!!$                           * dTHatdpsiHats(ispecies)/THat
!!$                      
!!$                      !!factorJ5 is only used in the Jacobian, it corresponds to the second term in factorJ3 divided by Phi1
!!$                      factorJ5 = Delta*nHat*mHat*sqrtMHat &
!!$                           /(2*pi*sqrtpi*(BHat(ialpha,izeta)**3)*sqrtTHat) &
!!$                           *(- BHat_sub_zeta(ialpha,izeta)*dBHatdtheta(ialpha,izeta) + BHat_sub_theta(ialpha,izeta)*dBHatdzeta(ialpha,izeta))&
!!$                           * DHat(ialpha,izeta) * xPartOfRHS2*gamma  & 
!!$                           * exp(-Z*gamma*Phi1Hat(ialpha,izeta)/THat)
!!$                      
!!$                      
!!$                      !!SECTION TEMPORARILY COMMENTED BY AM!!
!!$                      !!Add L=0 component
!!$                      call MatSetValue(matrix, rowIndex, colIndex,&
!!$                           -Z*gamma*(factorJ1 + factorJ2)/THat &
!!$                           -Z*gamma* (4/three)*factorJ3/THat &
!!$                           + factorJ4 &
!!$                           + (4/three)*factorJ5, ADD_VALUES, ierr) 
!!$                      
!!$                      !!Add L=2 component
!!$                      call MatSetValue(matrix, rowIndex2, colIndex,&
!!$                           -Z*gamma* (two/three)*factorJ3/THat &
!!$                           + (two/three)*factorJ5, ADD_VALUES, ierr)
!!$                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$                      !!TEST BY AM 2016-04-11!!
!!$!                      call MatSetValueSparse(matrix, rowIndex, colIndex,&
!!$!                           -Z*gamma*(factorJ1 + factorJ2)/THat &
!!$!                           -Z*gamma* (4/three)*factorJ3/THat &
!!$!                           + factorJ4 &
!!$!                           + (4/three)*factorJ5, ADD_VALUES, ierr) 
!!$!                      
!!$!                      !!Add L=2 component
!!$!                      call MatSetValueSparse(matrix, rowIndex2, colIndex,&
!!$!                           -Z*gamma* (two/three)*factorJ3/THat &
!!$!                           + (two/three)*factorJ5, ADD_VALUES, ierr)
!!$                      !!!!!!!!!!!!!!!!!!!!!!!!!
!!$                      
!!$                   end do
!!$                end do
!!$             end if
!!$             !!!!!!!!!!!!!!!!!!!!!!!
!!$          end do
!!$       end if       
!!$
!!$       ! *********************************************************
!!$       ! Add the nonlinear term in the residual, which also is
!!$       ! the block in the Jacobian from d(kinetic eqn) / d f.
!!$       ! *********************************************************
!!$       ! Note that this term is absent in the first iteration of a nonlinear calculation, because the term is proportional to Phi1.
!!$       ! Therefore, if reusePreconditioner=true, we do not include this term in the preconditioner (which is great because this term introduces a lot of nonzeros.)
!!$       ! If reusePreconditioner=false, then we DO include this term in the preconditioner.
!!$       ! We must use MatSetValue instead of MatSetValueSparse in this term so that we can add some "nonzero" elements whose value is actually 0
!!$       ! in the first iteration (due to Phi1=0). The reason is that PETSc will re-use the pattern of nonzeros from the first iteration in subsequent iterations.
!!$       ! However, we need not add elements which are 0 due to ddalpha=0 as opposed to because Phi1=0, since such elements will remain 0 at every iteration of SNES.
!!$
!!$       !if (nonlinear .and. (whichMatrix .ne. 2)) then
!!$       !!if (nonlinear .and. (whichMatrix == 1 .or. whichMatrix == 3 .or. (whichMatrix==0 .and. .not. reusePreconditioner))) then !!Commented by AM 2016-02
!!$       !if (.false.) then
!!$       if (includePhi1 .and. includePhi1InKineticEquation .and. (whichMatrix == 1 .or. whichMatrix == 3 .or. (whichMatrix==0 .and. .not. reusePreconditioner))) then !!Added by AM 2016-02
!!$
!!$          !print *,"@@@@@@ ",myRank," max(abs(dphi1hatdalpha)): ",maxval(abs(dphi1hatdalpha)),maxval(abs(dPhi1Hatdzeta))
!!$
!!$          allocate(nonlinearTerm_Lp1(Nx,Nx))
!!$          allocate(nonlinearTerm_Lm1(Nx,Nx))
!!$
!!$          do L=0,(Nxi-1)
!!$             if (L>0 .and. pointAtX0) then
!!$                ixMinCol = 2
!!$             else
!!$                ixMinCol = 1
!!$             end if
!!$
!!$             if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
!!$                ddx_to_use = ddx_preconditioner
!!$             else
!!$                ddx_to_use = ddx
!!$             end if
!!$
!!$             nonlinearTerm_Lp1 = (L+1)/(two*L+3) * ddx_to_use
!!$             nonlinearTerm_Lm1 = L/(two*L-1)     * ddx_to_use
!!$             do ix=1,Nx
!!$                nonlinearTerm_Lp1(ix,ix) = nonlinearTerm_Lp1(ix,ix) + (L+1)*(L+2)/(two*L+3)/x(ix)
!!$                nonlinearTerm_Lm1(ix,ix) = nonlinearTerm_Lm1(ix,ix) - (L-1)*L/(two*L-1)/x(ix)
!!$             end do
!!$             do ialpha=ialphaMin,ialphaMax
!!$
!!$                do izeta=izetaMinDKE,izetaMaxDKE
!!$                   factor = -gamma*Zs(ispecies)/(2*BHat(ialpha,izeta)*sqrtTHat*sqrtMHat) &
!!$                        * (BHat_sup_theta(ialpha,izeta)*dphi1hatdalpha(ialpha,izeta) &
!!$                        + BHat_sup_zeta(ialpha,izeta)*dPhi1Hatdzeta(ialpha,izeta))
!!$
!!$                   !do ix_row=ixMin,Nx
!!$                   do ix_row=max(ixMin,min_x_for_L(L)),Nx
!!$                      rowIndex = getIndex(ispecies,ix_row,L+1,ialpha,izeta,BLOCK_F)
!!$
!!$                      ! Term that is super-diagonal in L:
!!$                      if (L<Nxi-1) then
!!$                         ell = L + 1
!!$                         !do ix_col=ixMinCol,Nx
!!$                         do ix_col=max(ixMinCol,min_x_for_L(ell)),Nx
!!$                            if (abs(nonlinearTerm_Lp1(ix_row,ix_col))>threshholdForInclusion) then
!!$                               colIndex=getIndex(ispecies,ix_col,ell+1,ialpha,izeta,BLOCK_F)
!!$                               ! We must use MatSetValue instead of MatSetValueSparse here!!
!!$                               call MatSetValue(matrix, rowIndex, colIndex, &
!!$                                    factor*nonlinearTerm_Lp1(ix_row,ix_col), ADD_VALUES, ierr)
!!$                            end if
!!$                         end do
!!$                      end if
!!$                   
!!$                      ! Term that is sub-diagonal in L:
!!$                      if (L>0) then
!!$                         ell = L - 1
!!$                         !do ix_col=ixMinCol,Nx
!!$                         do ix_col=max(ixMinCol,min_x_for_L(ell)),Nx
!!$                            if (abs(nonlinearTerm_Lm1(ix_row,ix_col))>threshholdForInclusion) then
!!$                               colIndex=getIndex(ispecies,ix_col,ell+1,ialpha,izeta,BLOCK_F)
!!$                               ! We must use MatSetValue instead of MatSetValueSparse here!!
!!$                               call MatSetValue(matrix, rowIndex, colIndex, &
!!$                                    factor*nonlinearTerm_Lm1(ix_row,ix_col), ADD_VALUES, ierr)
!!$                            end if
!!$                         end do
!!$                      end if
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$          deallocate(nonlinearTerm_Lp1)
!!$          deallocate(nonlinearTerm_Lm1)
!!$          
!!$       end if
!!$
!!$       ! *********************************************************
!!$       ! Add the block in the Jacobian from d(kinetic eqn) / d Phi1
!!$       ! associated with the nonlinear term.
!!$       ! *********************************************************
!!$       ! Note that this term is absent in the first iteration of a nonlinear calculation, because the term is proportional to delta f.
!!$       ! Therefore, if reusePreconditioner=true, we do not include this term in the preconditioner (which is great because this term introduces a lot of nonzeros.)
!!$       ! If reusePreconditioner=false, then we DO include this term in the preconditioner.
!!$       ! We must use MatSetValue instead of MatSetValueSparse in this term so that we can add some "nonzero" elements whose value is actually 0
!!$       ! in the first iteration (due to delta f = 0). The reason is that PETSc will re-use the pattern of nonzeros from the first iteration in subsequent iterations.
!!$       ! However, we need not add elements which are 0 due to ddalpha=0 as opposed to because delta f = 0, since such elements will remain 0 at every iteration of SNES.
!!$
!!$       !!if (nonlinear .and. (whichMatrix==1 .or. (whichMatrix==0 .and. .not. reusePreconditioner))) then !!Commented by AM 2016-02
!!$       
!!$       ! This next line should be replaced with the line after!!!
!!$       !if (.false.) then
!!$       if (includePhi1 .and. includePhi1InKineticEquation .and. (whichMatrix==1 .or. (whichMatrix==0 .and. .not. reusePreconditioner))) then !!Added by AM 2016-02
!!$       !if (nonlinear .and. (whichMatrix==1)) then
!!$       !if (nonlinear .and. (whichMatrix==0 .or. whichMatrix==1)) then
!!$          allocate(tempVector1(Nx))
!!$          allocate(tempVector2(Nx))
!!$          do ialpha = ialphaMin,ialphaMax
!!$             do izeta = izetaMinDKE,izetaMaxDKE
!!$                factor = -gamma*Zs(ispecies)/(2*BHat(ialpha,izeta)*sqrtTHat*sqrtMHat)
!!$                do L=0,(Nxi-1)
!!$                   tempVector2=0
!!$
!!$                   if (L>0) then
!!$                      ! Add the delta_{L-1, ell} terms:
!!$                      ell = L-1
!!$                      tempVector1=0
!!$                      tempVector2=0
!!$                      !do ix=1,Nx
!!$                      do ix=min_x_for_L(ell),Nx
!!$                         index = getIndex(ispecies,ix,ell+1,ialpha,izeta,BLOCK_F)
!!$                         ! Add 1 because we are indexing a Fortran array instead of a PETSc object
!!$                         tempVector1(ix) = stateArray(index+1)
!!$                         tempVector2(ix) = tempVector2(ix) - (L-1)*L/(two*L-1)*stateArray(index+1)/x(ix)
!!$                      end do
!!$                      tempVector2 = tempVector2 + L/(two*L-1) * matmul(ddx,tempVector1)
!!$                   end if
!!$
!!$                   if (L<Nxi-1) then
!!$                      ! Add the delta_{L+1, ell} terms:
!!$                      ell = L+1
!!$                      tempVector1=0
!!$                      tempVector2=0
!!$                      !do ix=1,Nx
!!$                      do ix=min_x_for_L(ell),Nx
!!$                         index = getIndex(ispecies,ix,ell+1,ialpha,izeta,BLOCK_F)
!!$                         ! Add 1 because we are indexing a Fortran array instead of a PETSc object
!!$                         tempVector1(ix) = stateArray(index+1)
!!$                         tempVector2(ix) = tempVector2(ix) + (L+1)*(L+2)/(two*L+3)*stateArray(index+1)/x(ix)
!!$                      end do
!!$                      tempVector2 = tempVector2 + (L+1)/(two*L+3) * matmul(ddx,tempVector1)
!!$                   end if
!!$
!!$                   !do ix=ixMin,Nx
!!$                   do ix=max(ixMin,min_x_for_L(L)),Nx
!!$                      rowIndex = getIndex(ispecies,ix,L+1,ialpha,izeta,BLOCK_F)
!!$
!!$                      ! Add the d Phi_1 / d theta term:
!!$                      do j=1,Nalpha
!!$                         if (abs(ddalpha(ialpha,j))>threshholdForInclusion) then
!!$                            colIndex = getIndex(1,1,1,j,izeta,BLOCK_QN)
!!$                            ! We must use MatSetValue instead of MatSetValueSparse here!!
!!$                            call MatSetValue(matrix, rowIndex, colIndex, &
!!$                                 factor*BHat_sup_theta(ialpha,izeta)*ddalpha(ialpha,j)*tempVector2(ix), &
!!$                                 ADD_VALUES, ierr)
!!$                         end if
!!$                      end do
!!$
!!$                      ! Add the d Phi_1 / d zeta term:
!!$                      do j=1,Nzeta
!!$                         if (abs(ddzeta(izeta,j))>threshholdForInclusion) then
!!$                            colIndex = getIndex(1,1,1,ialpha,j,BLOCK_QN)
!!$                            ! We must use MatSetValue instead of MatSetValueSparse here!!
!!$                            call MatSetValue(matrix, rowIndex, colIndex, &
!!$                                 factor*BHat_sup_zeta(ialpha,izeta)*ddzeta(izeta,j)*tempVector2(ix), &
!!$                                 ADD_VALUES, ierr)
!!$                         end if
!!$                      end do
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$          deallocate(tempVector1)
!!$          deallocate(tempVector2)
!!$       end if

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
          ddx_to_use = ddx
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
                        + (erfs(ix)-Psi_Chandra(ix))/x2(ix)) * ddx_to_use(ix,:))
                   
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
                            M11 = M11 + RosenbluthPotentialTerms(iSpeciesA,iSpeciesB,L+1,:,:)
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
                      
                      !if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                      if (whichMatrix==0) then
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
                      
                      do ialpha=ialphaMin,ialphaMax
                         do izeta=izetaMinDKE,izetaMaxDKE
                            !do ix_row=ixMin,Nx
                            do ix_row=max(ixMin,min_x_for_L(L)),Nx
                               rowIndex=getIndex(iSpeciesA,ix_row,L+1,ialpha,izeta,BLOCK_F)
                               !do ix_col = ixMinCol,Nx
                               do ix_col = max(ixMinCol,min_x_for_L(L)),Nx
                                  colIndex=getIndex(iSpeciesB,ix_col,L+1,ialpha,izeta,BLOCK_F)
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

          allocate(pitch_angle_scattering_operator_to_use(Nxi,Nxi))
          if (whichMatrix==0) then
             pitch_angle_scattering_operator_to_use = pitch_angle_scattering_operator_preconditioner
          else
             pitch_angle_scattering_operator_to_use = pitch_angle_scattering_operator
          end if

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
             
             do ix=ixMin,Nx
                do ialpha=ialphaMin,ialphaMax
                   do izeta=izetaMinDKE,izetaMaxDKE
                      factor = -nu_n*BHat(ialpha,izeta)*sqrt(mHats(ispeciesA)/THats(ispeciesA))/DHat(ialpha,izeta)*nuDHat(iSpeciesA,ix)
                      do ixi_row = 1,Nxi
                         rowIndex=getIndex(iSpeciesA,ix,ixi_row,ialpha,izeta,BLOCK_F)
                         do ixi_col = 1,Nxi
                            colIndex=getIndex(iSpeciesA,ix,ixi_col,ialpha,izeta,BLOCK_F)
                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                 factor*pitch_angle_scattering_operator_to_use(ixi_row,ixi_col), ADD_VALUES, ierr)
                         end do
                      end do
                   end do
                end do
             end do
          end do
          deallocate(pitch_angle_scattering_operator_to_use)
          
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

!!$    ! *******************************************************************************
!!$    ! If there is a grid point at x=0, add the boundary conditions for f at x=0.
!!$    ! *******************************************************************************
!!$
!!$    if (pointAtX0 .and. whichMatrix .ne. 2) then
!!$       stop "Not set up for pointAtX0"
!!$       ! For L > 0 modes, impose f=0 at x=0:
!!$       ix = 1
!!$       do L = 1,(Nxi_for_x(ix)-1)
!!$          do ialpha = ialphaMin,ialphaMax
!!$             do izeta = izetaMin,izetaMax
!!$                do ispecies = 1,Nspecies
!!$                   index = getIndex(ispecies,ix,L+1,ialpha,izeta,BLOCK_F)
!!$                   call MatSetValue(matrix, index, index, one, ADD_VALUES, ierr)
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$
!!$       ! For L=0 mode, impose regularity (df/dx=0) at x=0:
!!$       L=0
!!$       if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
!!$          ddx_to_use = ddx_preconditioner
!!$       else
!!$          ddx_to_use = ddx
!!$       end if
!!$       ix_row = 1
!!$       do ialpha = ialphaMin,ialphaMax
!!$          do izeta = izetaMin,izetaMax
!!$             do ispecies = 1,Nspecies
!!$                rowIndex = getIndex(ispecies,ix_row,L+1,ialpha,izeta,BLOCK_F)
!!$                do ix_col = 1,Nx
!!$                   colIndex = getIndex(ispecies,ix_col,L+1,ialpha,izeta,BLOCK_F)
!!$                   call MatSetValueSparse(matrix, rowIndex, colIndex, ddx_to_use(ix_row,ix_col), ADD_VALUES, ierr)
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end if

    ! *******************************************************************************
    ! Add sources:
    ! *******************************************************************************

    if (whichMatrix .ne. 2) then
       select case (constraintScheme)
       case (0)
          ! Do nothing here.
          
       case (1,3,4)
          ! Add a heat source and a particle source.
          
          do ix = ixMin,Nx
             ! The particle and energy sources S_p and S_e are normalized so that
             ! \int dx_0^inf (4 pi x^2) S_p = 1
             ! \int dx_0^inf (4 pi x^2) S_e = 0
             ! \int dx_0^inf (4 pi x^2) S_p x^2 = 0
             ! \int dx_0^inf (4 pi x^2) S_e x^2 = 1
             ! See Mathematica notebook 20160107-01 for calculation of the coefficients used here.
             select case (constraintScheme)
             case (1)
                ! Constant and quadratic terms:
                xPartOfSource1 = (         -x2(ix) + 5/two) * exp(-x2(ix)) / (pi*sqrtpi) ! Provides particles but no heat
                xPartOfSource2 = (two/three*x2(ix) -     1) * exp(-x2(ix)) / (pi*sqrtpi) ! Provides heat but no particles
                ! Definition prior to 2016-01-07, which differs just in the overall constant:
                !xPartOfSource1 = (x2(ix)-5/two)*exp(-x2(ix)) ! Provides particles but no heat
                !xPartOfSource2 = (x2(ix)-3/two)*exp(-x2(ix)) ! Provides heat but no particles
             case (3)
                ! Constant and quartic terms:
                xPartOfSource1 = (       -one/5*x2(ix)*x2(ix) + 7/(4.0d-0)) * exp(-x2(ix)) / (pi*sqrtpi) ! Provides particles but no heat
                xPartOfSource2 = (two/(15.0d-0)*x2(ix)*x2(ix) -   (0.5d-0)) * exp(-x2(ix)) / (pi*sqrtpi) ! Provides heat but no particles
             case (4)
                ! Quadratic and quartic terms:
                xPartOfSource1 = ( -two/three*x2(ix)*x2(ix) +   7/(3.0d-0)*x2(ix)) * exp(-x2(ix)) / (pi*sqrtpi) ! Provides particles but no heat
                xPartOfSource2 = (4/(15.0d-0)*x2(ix)*x2(ix) - two/(3.0d-0)*x2(ix)) * exp(-x2(ix)) / (pi*sqrtpi) ! Provides heat but no particles
             case default
                stop "Invalid constraintScheme!"
             end select
             do ialpha = ialphaMin,ialphaMax
                do izeta = izetaMinDKE,izetaMaxDKE
                   factor = (B0OverBBar/GHat)*BHat(ialpha,izeta)/DHat(ialpha,izeta)
                   do ispecies = 1,Nspecies
                      do ixi = 1,Nxi_for_x(ix)
                         rowIndex = getIndex(ispecies, ix, ixi, ialpha, izeta, BLOCK_F)
                         
                         colIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, factor*xPartOfSource1, ADD_VALUES, ierr)
                         
                         colIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, factor*xPartOfSource2, ADD_VALUES, ierr)
                      end do
                   end do
                end do
             end do
          end do
          
       case (2)
          ! Add a source (which is constant on the flux surface and independent of xi) at each x.
          do ialpha = ialphaMin,ialphaMax
             do izeta = izetaMinDKE,izetaMaxDKE
                factor = (B0OverBBar/GHat)*BHat(ialpha,izeta)/DHat(ialpha,izeta)
                do ix = ixMin,Nx
                   do ixi = 1,Nxi_for_x(ix)
                      do ispecies = 1,Nspecies
                         rowIndex = getIndex(ispecies, ix, ixi, ialpha, izeta, BLOCK_F)
                         colIndex = getIndex(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT)
                         call MatSetValue(matrix, rowIndex, colIndex, factor, ADD_VALUES, ierr)
                      end do
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

       case (1,3,4)
          ! Force the flux-surface-averaged perturbed density and 
          ! flux-surface-averaged perturbed pressure to be zero.

          do ialpha=1,Nalpha
             do izeta=1,Nzeta
                do ixi=1,Nxi
                   ! The matrix elements must be proportional to 1/DHat, and the sum should exclude the repeated values of zeta,
                   ! but otherwise the row could probably be scaled any way you like. Here we include a factor 1/VPrimeHat so the 
                   ! row is dimensionless and the row sum is O(1), which seems like a reasonable scaling.
                   !factor = alphaWeights(ialpha)*zetaWeights(izeta)/DHat(ialpha,izeta)
                   factor = xiWeights(ixi)*alphaWeights(ialpha)*zetaWeights(izeta)/(DHat(ialpha,izeta)*VPrimeHat)

                   do ix=1,Nx
                      do ispecies=1,Nspecies
                         colIndex = getIndex(ispecies, ix, ixi, ialpha, izeta, BLOCK_F)
                         
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
          end do

       case (2)
          ! Force the flux-surface-averaged distribution function to be zero
          ! at each value of x:

          do ialpha=1,Nalpha
             do izeta=1,Nzeta
                do ixi=1,Nxi
                   ! The matrix elements must be proportional to 1/DHat, and the sum should exclude the repeated values of zeta,
                   ! but otherwise the row could probably be scaled any way you like. Here we include a factor 1/VPrimeHat so the 
                   ! row is dimensionless and the row sum is O(1), which seems like a reasonable scaling.
                   !factor = alphaWeights(ialpha)*zetaWeights(izeta)/DHat(ialpha,izeta)
                   factor = xiWeights(ixi)*alphaWeights(ialpha)*zetaWeights(izeta)/(DHat(ialpha,izeta)*VPrimeHat)

                   do ix=ixMin,Nx
                      do ispecies = 1,Nspecies
                         colIndex = getIndex(ispecies, ix, ixi, ialpha, izeta, BLOCK_F)
                         rowIndex = getIndex(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              factor, ADD_VALUES, ierr)
                      end do
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
    ! SECTION MODIFIED BY AM 2016-02/03
    ! Add the quasineutrality equation
    ! *******************************************************************************

    if (whichMatrix .ne. 2 .and. includePhi1) then
       L=0
       do ialpha = ialphaMin,ialphaMax
          do izeta = izetaMinDKE,izetaMaxDKE
             rowIndex = getIndex(1, 1, 1, ialpha, izeta, BLOCK_QN)

             ! Add the charge density of the kinetic species (the part of the residual/Jacobian related to f1):
             do ispecies = 1,Nspecies
                !!speciesFactor = Zs(ispecies)*THats(ispecies)/mHats(ispecies) & !!Commented by AM 2016-02
                !!     * sqrt(THats(ispecies)/mHats(ispecies)) !!Commented by AM 2016-02
                speciesFactor = four*pi*Zs(ispecies)*THats(ispecies)/mHats(ispecies) & !!Added by AM 2016-02 !!The factor 4pi comes from velocity integration over xi.
                     * sqrt(THats(ispecies)/mHats(ispecies)) !!Added by AM 2016-02

                do ix = 1,Nx
                   colIndex = getIndex(ispecies, ix, L+1, ialpha, izeta, BLOCK_F)
                   call MatSetValueSparse(matrix, rowIndex, colIndex, &
                        speciesFactor*x2(ix)*xWeights(ix), ADD_VALUES, ierr)
                end do

                !!Added by AM 2016-02!!
                if (quasineutralityOption == 2) then
                   exit !!If running with EUTERPE equations we only use the first kinetic species in quasi-neutrality.
                end if
                !!!!!!!!!!!!!!!!!!!!!!!

             end do
             !!!!!!!!!!!!!!!!!!!!!!!
             !!Added by AM 2016-02!!

             !Add the part of the residual/Jacobian related to Phi1
             colIndex = getIndex(1, 1, 1, ialpha, izeta, BLOCK_QN)
             if (quasineutralityOption == 1 .and. (whichMatrix == 0 .or. whichMatrix == 1)) then !!Only add in Jacobian matrix
                do ispecies = 1,Nspecies
                   call MatSetValueSparse(matrix, rowIndex, colIndex, &
                        - gamma * Zs(ispecies) * Zs(ispecies) * NHats(ispecies) * exp (- Zs(ispecies)* gamma * Phi1Hat(ialpha,izeta) / THats(ispecies)) / THats(ispecies), ADD_VALUES, ierr)
                end do
                if (withAdiabatic) then
                   call MatSetValueSparse(matrix, rowIndex, colIndex, &
                        - gamma * adiabaticZ * adiabaticZ *adiabaticNHat * exp (- adiabaticZ* gamma * Phi1Hat(ialpha,izeta) / adiabaticTHat) / adiabaticTHat, ADD_VALUES, ierr)
                end if
             else if (quasineutralityOption == 2 .and. withAdiabatic .and. Nspecies > 0) then
                call MatSetValueSparse(matrix, rowIndex, colIndex, &
                     - gamma * (Zs(1)*Zs(1)*NHats(1)/THats(1) + adiabaticZ*adiabaticZ*adiabaticNHat/adiabaticTHat), ADD_VALUES, ierr)                
             end if
	     !!!!!!!!!!!!!!!!!!!!!!!

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
       do ialpha=1,Nalpha
          do izeta=1,Nzeta
             colIndices(izeta) = getIndex(1, 1, 1, ialpha, izeta, BLOCK_QN)
          end do

          call MatSetValuesSparse(matrix, 1, rowIndices(1), Nzeta, colIndices, &
               alphaWeights(ialpha)*zetaWeights/DHat(ialpha,:), ADD_VALUES, ierr)
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

    ! *******************************************************************************
    ! If requested, save the matrix.
    ! *******************************************************************************

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

    integer :: ix, ixi, ialpha, izeta, ispecies, index
    PetscScalar :: factor, pi32
    PetscErrorCode :: ierr

    if (masterProc) then
       print *,"Initializing f0"
    end if

    call VecSet(f0, zero, ierr)
    
    pi32 = pi*sqrt(pi)
    do ispecies = 1,Nspecies
       do ix = 1,Nx
          !factor = nHats(ispecies)*mHats(ispecies)/(pi*THats(ispecies)) &
          !     * sqrt(mHats(ispecies)/(pi*THats(ispecies))) * expx2(ix)
          factor = expx2(ix) / (pi32)
          do ialpha = ialphaMin,ialphaMax
             do izeta = izetaMin,izetaMax
                do ixi = 1,Nxi_for_x(ix)
                   index = getIndex(ispecies, ix, ixi, ialpha, izeta, BLOCK_F)
                   call VecSetValue(f0, index, &
                !!     factor, INSERT_VALUES, ierr) !!Commented by AM 2016-06
                        exp(-Zs(ispecies)*gamma*Phi1Hat(ialpha,izeta)/THats(ispecies))*factor, INSERT_VALUES, ierr) !!Added by AM 2016-06
                end do
             end do
          end do
       end do
    end do

    call VecAssemblyBegin(f0, ierr)
    call VecAssemblyEnd(f0, ierr)

  end subroutine init_f0
