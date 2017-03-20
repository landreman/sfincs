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
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, speciesFactor, speciesFactor2, v_s
    PetscScalar :: T32m, factor, spatial_factor, LFactor, temp, temp1, temp2, xDotFactor, xDotFactor2, stuffToAdd
    PetscScalar :: factor1, factor2, factorJ1, factorJ2, factorJ3, factorJ4, factorJ5  !!Added by AM 2016-03
    PetscScalar, dimension(:), allocatable :: xb, expxb2
    PetscScalar, dimension(:,:), allocatable :: xPartOfXDot_plus, xPartOfXDot_minus, xPartOfXDot
    !PetscScalar, dimension(:,:), allocatable :: ddx_to_use_plus, ddx_to_use_minus
    integer :: i, j, ix, ispecies, itheta, izeta, L, ixi, index, ix_row, ix_col, ixi_row, ixi_col
    integer :: rowIndex, colIndex
    integer :: rowIndex1, rowIndex2, L2 !!Added by AM 2016-03
    integer :: ell, iSpeciesA, iSpeciesB, maxL
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar, dimension(:,:), allocatable :: d2dx2_to_use, zetaPartOfTerm, thetaPartOfTerm
    PetscScalar, dimension(:,:), allocatable :: fToFInterpolationMatrix
    PetscScalar, dimension(:,:), allocatable :: potentialsToFInterpolationMatrix
    PetscScalar, dimension(:,:,:,:), allocatable :: CECD
    PetscScalar :: xPartOfSource1, xPartOfSource2, geometricFactor1, geometricFactor2, geometricFactor3
    PetscScalar, dimension(:,:), allocatable ::  nuDHat
    PetscScalar, dimension(:), allocatable :: erfs, Psi_Chandra
    PetscLogDouble :: time1, time2, time3, time4
    PetscScalar, dimension(:,:), pointer :: ddtheta_to_use, ddzeta_to_use, ddxi_to_use, pitch_angle_scattering_operator_to_use, ddx_to_use
    double precision :: myMatInfo(MAT_INFO_SIZE)
    integer :: NNZ, NNZAllocated, NMallocs
    PetscScalar :: dfMdx
    character(len=200) :: whichMatrixName, filename
    PetscViewer :: viewer
    integer :: itheta_row, itheta_col, izeta_row, izeta_col, ixMin, ixMinCol
    VecScatter :: vecScatterContext
    Vec :: vecOnEveryProc
    PetscScalar, pointer :: stateArray(:)
    logical :: useStateVec
    PetscScalar, dimension(:,:), allocatable :: nonlinearTerm_Lp1, nonlinearTerm_Lm1
    PetscScalar, dimension(:), allocatable :: tempVector1, tempVector2
    PetscScalar, dimension(:,:), allocatable :: tempExtrapMatrix, fToFInterpolationMatrix_plus1
    PetscScalar :: dPhiHatdpsiHatToUseInRHS, xPartOfRHS, xPartOfRHS2 !!Added by AM 2016-03
    PetscScalar, dimension(:,:), allocatable :: theta_interpolation_matrix
    integer :: zeta_shift, zeta_pad_size
    integer :: stencil, sign, izeta_interpolation
    integer :: iSpecies_min, iSpecies_max, ix_min, ix_max
    PetscScalar, dimension(:,:), allocatable :: collision_operator_xi_block
    PetscScalar, dimension(:,:,:), allocatable :: Legendre_projection_to_use
    PetscScalar :: Er_factor

    ! *******************************************************************************
    ! Do a few sundry initialization tasks:
    ! *******************************************************************************

    Er_factor = Delta * gamma / 2 ! When I switch to SI units, I can replace Er_factor with 1

    call PetscTime(time3, ierr)

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
          do itheta = 1,Ntheta
             do izeta = 1,Nzeta
                index = getIndex(1,1,1,itheta,izeta,BLOCK_QN)
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

    allocate(d2dx2_to_use(Nx,Nx))

    allocate(xb(Nx))
    allocate(expxb2(Nx))
    allocate(erfs(Nx))
    allocate(Psi_Chandra(Nx))
    allocate(nuDHat(Nspecies, Nx))
    allocate(fToFInterpolationMatrix(Nx,Nx))
    allocate(potentialsToFInterpolationMatrix(Nx, NxPotentials))
    allocate(CECD(Nspecies, Nspecies, Nx, Nx))

    call PetscTime(time4, ierr)
    if (masterProc) print *,"  Time for init:",time4-time3
    call PetscTime(time3, ierr)
    
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
       v_s = sqrt(THat/mHat) ! Using the version3 normalizations, there is no 2 here. But when I move to SI units, add a 2 inside the sqrt.

       ! *********************************************************
       ! Add the d/dzeta term:
       ! *********************************************************
       
       if (masterProc) print *,"Beginning d/dzeta"
       if ((whichMatrix .ne. 2) .and. (Nzeta > 1)) then
          do izeta_row = izetaMin,izetaMax
             do itheta = ithetaMin,ithetaMax
                do ix = 1,Nx
                   do ixi = 1,Nxi_for_x(ix)
                      ! Next line is the parallel streaming term:
                      factor = v_s * x(ix) * xi(ixi) * BHat_sup_zeta(itheta,izeta_row) / BHat(itheta,izeta_row)
                      ! Add the ExB term:
                      select case (ExB_option)
                      case (1)
                         factor = factor - BHat_sub_theta(itheta,izeta_row)/(sqrt_g(itheta,izeta_row)*BHat(itheta,izeta_row)*BHat(itheta,izeta_row))*dPhiHatdpsiHat*Er_factor
                      case (2)
                         factor = factor - BHat_sub_theta(itheta,izeta_row)/(sqrt_g(itheta,izeta_row)*FSABHat2)*dPhiHatdpsiHat*Er_factor
                      end select
                      factor = factor * spatial_scaling(itheta,izeta_row) * x_scaling(ix,ispecies)

                      if (whichMatrix>0) then
                         ! Not preconditioner
                         if (factor>0) then
                            ddzeta_to_use => ddzeta_plus
                         else
                            ddzeta_to_use => ddzeta_minus
                         end if
                      else
                         ! Preconditioner
                         if (factor>0) then
                            ddzeta_to_use => ddzeta_plus_preconditioner
                         else
                            ddzeta_to_use => ddzeta_minus_preconditioner
                         end if
                      end if
                      rowIndex = getIndex(ispecies, ix, ixi, itheta, izeta_row, BLOCK_F)
                      do izeta_col = 1,Nzeta
                         colIndex = getIndex(ispecies, ix, ixi, itheta, izeta_col, BLOCK_F)
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
       if (masterProc) print *,"Done with d/dzeta"

       
       ! *********************************************************
       ! Add the d/dtheta term:
       ! *********************************************************

       if (whichMatrix .ne. 2) then
          do izeta = izetaMin,izetaMax
             do itheta_row = ithetaMin,ithetaMax
                do ix = 1,Nx
                   do ixi = 1,Nxi_for_x(ix)
                      ! Next line is the parallel streaming term:
                      factor = v_s * x(ix) * xi(ixi) * BHat_sup_theta(itheta_row,izeta) / BHat(itheta_row,izeta)
                      ! Add the ExB term:
                      select case (ExB_option)
                      case (1)
                         factor = factor + BHat_sub_zeta(itheta_row,izeta)/(sqrt_g(itheta_row,izeta)*BHat(itheta_row,izeta)*BHat(itheta_row,izeta))*dPhiHatdpsiHat*Er_factor
                      case (2)
                         factor = factor + BHat_sub_zeta(itheta_row,izeta)/(sqrt_g(itheta_row,izeta)*FSABHat2)*dPhiHatdpsiHat*Er_factor
                      end select
                      factor = factor * spatial_scaling(itheta_row,izeta) * x_scaling(ix,ispecies)

                      if (whichMatrix>0) then
                         ! Not preconditioner
                         if (factor>0) then
                            ddtheta_to_use => ddtheta_plus
                         else
                            ddtheta_to_use => ddtheta_minus
                         end if
                      else
                         ! Preconditioner
                         if (factor>0) then
                            ddtheta_to_use => ddtheta_plus_preconditioner
                         else
                            ddtheta_to_use => ddtheta_minus_preconditioner
                         end if
                      end if

                      rowIndex = getIndex(ispecies, ix, ixi, itheta_row, izeta, BLOCK_F)
                      do itheta_col = 1,Ntheta
                         colIndex = getIndex(ispecies, ix, ixi, itheta_col, izeta, BLOCK_F)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              factor * ddtheta_to_use(itheta_row,itheta_col), ADD_VALUES, ierr)
                      end do
                   end do
                end do
             end do
          end do
       end if
       ! To generate an error if these variables are used by mistake later:
       itheta_row = -1
       itheta_col = -1


       ! *********************************************************
       ! Add the d/dxi term:
       ! *********************************************************

       if (whichMatrix .ne. 2) then
          do izeta = izetaMin,izetaMax
             do itheta = ithetaMin,ithetaMax
                do ix = 1,Nx
                   do ixi_row = 1,Nxi_for_x(ix)
                      ! The next line is the standard mirror term:
                      factor = - v_s*x(ix)*(1-xi(ixi_row)*xi(ixi_row))/(2*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                           * (BHat_sup_theta(itheta,izeta)*dBHatdtheta(itheta,izeta) + BHat_sup_zeta(itheta,izeta)*dBHatdzeta(itheta,izeta))
                      if (includeElectricFieldTermInXiDot) then
                         factor = factor + xi(ixi_row)*(1-xi(ixi_row)*xi(ixi_row))/(2*sqrt_g(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                              * (BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))*dPhiHatdpsiHat*Er_factor
                      end if
                      factor = factor * spatial_scaling(itheta,izeta) * x_scaling(ix,ispecies)

                      if (whichMatrix>0) then
                         ! Not preconditioner
                         if (factor>0) then
                            ddxi_to_use => ddxi_plus
                         else
                            ddxi_to_use => ddxi_minus
                         end if
                      else
                         ! Preconditioner
                         if (factor>0) then
                            ddxi_to_use => ddxi_plus_preconditioner
                         else
                            ddxi_to_use => ddxi_minus_preconditioner
                         end if
                      end if
                      rowIndex = getIndex(ispecies, ix, ixi_row, itheta, izeta, BLOCK_F)
                      do ixi_col = 1,Nxi_for_x(ix)
                         colIndex = getIndex(ispecies, ix, ixi_col, itheta, izeta, BLOCK_F)
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


       ! *********************************************************
       ! Add the d/dx term:
       ! *********************************************************

       if ((whichMatrix .ne. 2) .and. includeXDotTerm) then
          do izeta = izetaMin,izetaMax
             do itheta = ithetaMin,ithetaMax
                spatial_factor = spatial_scaling(itheta,izeta)*dPhiHatdpsiHat/(2*sqrt_g(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                     * (BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))*Er_factor
                do ix_row = 1,Nx
                   do ixi = 1,Nxi
                      factor = spatial_factor * x_scaling(ix_row,ispecies) * x(ix_row) * (1 + xi(ixi)*xi(ixi))
                      if (whichMatrix>0) then
                         ! Not preconditioner
                         ddx_to_use => ddx
                      else
                         ! Preconditioner
                         ddx_to_use => ddx_preconditioner
                      end if
                      rowIndex = getIndex(ispecies, ix_row, ixi, itheta, izeta, BLOCK_F)
                      do ix_col = 1,Nx
                         colIndex = getIndex(ispecies, ix_col, ixi, itheta, izeta, BLOCK_F)
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
!!$             itheta = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$             do itheta_row = ithetaMin, ithetaMax
!!$                do izeta = izetaMin, izetaMax
!!$                   rowIndex = getIndex(ispecies, ix, L+1, itheta_row, izeta, BLOCK_F)
!!$
!!$! COMMENTED BY AM 2016-02
!!$!                   factor = -gamma*Delta*DHat(itheta_row,izeta)*BHat_sub_zeta(itheta_row,izeta) &
!!$!                        /(two*BHat(itheta_row,izeta)*BHat(itheta_row,izeta)) &
!!$!                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$!                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)
!!$
!!$                   !!Added by AM 2016-03!!
!!$                   factor1 = -exp(-Z*gamma*Phi1Hat(itheta_row,izeta)/THat)*gamma*Delta*DHat(itheta_row,izeta)*BHat_sub_zeta(itheta_row,izeta) &
!!$                        /(two*BHat(itheta_row,izeta)*BHat(itheta_row,izeta)) &
!!$                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)
!!$
!!$                   factor2 = - gamma*gamma*Delta*DHat(itheta_row,izeta)/(two*pi*sqrtpi*BHat(itheta_row,izeta)*BHat(itheta_row,izeta)) &
!!$                        * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
!!$                        * exp(-Z*gamma*Phi1Hat(itheta_row,izeta)/THat)*expx2(ix)*BHat_sub_zeta(itheta_row,izeta) &
!!$                        * (dPhiHatdpsiHat + Phi1Hat(itheta_row,izeta)*dTHatdpsiHats(ispecies)/THat)
!!$
!!$                   do itheta_col = 1,Ntheta
!!$                      colIndex = getIndex(1,1,1,itheta_col, izeta, BLOCK_QN)
!!$!!                      call MatSetValueSparse(matrix, rowIndex, colIndex,& !!Commented by AM 2016-03
!!$!!                           factor*ddtheta(itheta_row, itheta_col), ADD_VALUES, ierr) !!Commented by AM 2016-03
!!$                      !!! SHOULD THE NEXT LINE BE MatSetValueSparse???
!!$                      !!NEXT CALL TEMPORARILY COMMENTED BY AM!!
!!$!                      call MatSetValue(matrix, rowIndex, colIndex,& !!Added by AM 2016-03
!!$!                           (factor1 + factor2)*ddtheta(itheta_row, itheta_col), ADD_VALUES, ierr) !!Added by AM 2016-03
!!$
!!$                      !!TEST BY AM 2016-04-11!!
!!$                      call MatSetValueSparse(matrix, rowIndex, colIndex,& !!Added by AM 2016-03
!!$                           (factor1 + factor2)*ddtheta(itheta_row, itheta_col), ADD_VALUES, ierr) !!Added by AM 2016-03
!!$                      !!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$                   end do
!!$                end do
!!$             end do
!!$
!!$             ! d Phi_1 / d zeta term:
!!$             izeta = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$             itheta_row = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$             itheta_col = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$             do itheta = ithetaMin, ithetaMax
!!$                do izeta_row = izetaMin, izetaMax
!!$                   rowIndex = getIndex(ispecies, ix, L+1, itheta, izeta_row, BLOCK_F)
!!$
!!$! COMMENTED BY AM 2016-02
!!$!                   factor = gamma*Delta*DHat(itheta_row,izeta)*BHat_sub_theta(itheta_row,izeta) &
!!$!                        /(two*BHat(itheta_row,izeta)*BHat(itheta_row,izeta)) &
!!$!                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$!                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)
!!$                   
!!$                   !! ADDED BY AM 2016-02!!
!!$                   factor1 = exp(-Z*gamma*Phi1Hat(itheta,izeta_row)/THat)*gamma*Delta*DHat(itheta,izeta_row)*BHat_sub_theta(itheta,izeta_row) &
!!$                        /(two*BHat(itheta,izeta_row)*BHat(itheta,izeta_row)) &
!!$                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)
!!$                   
!!$                   factor2 = gamma*gamma*Delta*DHat(itheta,izeta_row)/(two*pi*sqrtpi*BHat(itheta,izeta_row)*BHat(itheta,izeta_row)) &
!!$                        * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
!!$                        * exp(-Z*gamma*Phi1Hat(itheta,izeta_row)/THat)*expx2(ix)*BHat_sub_theta(itheta,izeta_row) &
!!$                        * (dPhiHatdpsiHat + Phi1Hat(itheta,izeta_row)*dTHatdpsiHats(ispecies)/THat)
!!$                   
!!$                   do izeta_col = 1,Nzeta
!!$                      colIndex = getIndex(1,1,1,itheta, izeta_col, BLOCK_QN)
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
!!$                itheta = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$                izeta = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$                itheta_row = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$                itheta_col = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$                izeta_row = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$                izeta_col = -1 ! Just so a runtime error occurs if I use itheta by mistake
!!$                do itheta = ithetaMin, ithetaMax
!!$                   do izeta = izetaMin, izetaMax
!!$                      rowIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
!!$                      rowIndex2 = getIndex(ispecies, ix, L2+1, itheta, izeta, BLOCK_F) !!Added by AM 2016-03 to add extra P_2 terms
!!$                      
!!$                      !!These terms are diagonal in theta and zeta
!!$                      colIndex = getIndex(1,1,1,itheta, izeta, BLOCK_QN)
!!$                      
!!$                      !!The following factors are only used in the Jacobian. These terms do not contain d/dtheta or d/dzeta
!!$                      !!but d(Phi1)/dtheta, d(Phi1)/dzeta, dB/dtheta, dB/dzeta.
!!$                      !!Therefore they are the same when adding the d/dtheta terms as when adding the d/dzeta terms.
!!$                      
!!$                      !!factorJ1 stems from factor1 but contains both dphi1hatdtheta and dPhi1Hatdzeta parts
!!$                      factorJ1 = exp(-Z*gamma*Phi1Hat(itheta,izeta)/THat)*gamma*Delta*DHat(itheta,izeta) &
!!$                           *(- BHat_sub_zeta(itheta,izeta)* dPhi1Hatdtheta(itheta, izeta) + BHat_sub_theta(itheta, izeta)*dPhi1Hatdzeta(itheta, izeta))&
!!$                           /(two*BHat(itheta,izeta)*BHat(itheta,izeta)) &
!!$                           * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
!!$                           * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)  
!!$                      
!!$                      !!factorJ2 stems from factor2 but contains both dphi1hatdtheta and dPhi1Hatdzeta parts
!!$                      factorJ2 = gamma*gamma*Delta*DHat(itheta,izeta)/(two*pi*sqrtpi*BHat(itheta,izeta)*BHat(itheta,izeta)) &
!!$                           * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
!!$                           * exp(-Z*gamma*Phi1Hat(itheta,izeta)/THat)*expx2(ix) &
!!$                           * (- BHat_sub_zeta(itheta,izeta)*dphi1hatdtheta(itheta, izeta) + BHat_sub_theta(itheta, izeta)*dPhi1Hatdzeta(itheta, izeta)) &
!!$                           * (dPhiHatdpsiHat + Phi1Hat(itheta,izeta)*dTHatdpsiHats(ispecies)/THat)
!!$                      
!!$                      !!factorJ3 is only used in the Jacobian, it corresponds to a part \propto exp(-Ze Phi1/T) which is
!!$                      !!implemented in evaluateResidual.F90 for the residual
!!$                      factorJ3 = Delta*nHat*mHat*sqrtMHat &
!!$                           /(2*pi*sqrtpi*Z*(BHat(itheta,izeta)**3)*sqrtTHat) &
!!$                           *(- BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) + BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
!!$                           * DHat(itheta,izeta) * (xPartOfRHS + xPartOfRHS2*Z*gamma*Phi1Hat(itheta,izeta))  & 
!!$                           * exp(-Z*gamma*Phi1Hat(itheta,izeta)/THat)
!!$                      
!!$                      !!factorJ4 is only used in the Jacobian, it corresponds to the second term in factorJ2 divided by Phi1
!!$                      factorJ4 =  gamma*gamma*Delta*DHat(itheta,izeta)/(two*pi*sqrtpi*BHat(itheta,izeta)*BHat(itheta,izeta)) &
!!$                           * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
!!$                           * exp(-Z*gamma*Phi1Hat(itheta,izeta)/THat)*expx2(ix) &
!!$                           *(- BHat_sub_zeta(itheta,izeta)*dphi1hatdtheta(itheta,izeta) + BHat_sub_theta(itheta,izeta)*dPhi1Hatdzeta(itheta,izeta))&
!!$                           * dTHatdpsiHats(ispecies)/THat
!!$                      
!!$                      !!factorJ5 is only used in the Jacobian, it corresponds to the second term in factorJ3 divided by Phi1
!!$                      factorJ5 = Delta*nHat*mHat*sqrtMHat &
!!$                           /(2*pi*sqrtpi*(BHat(itheta,izeta)**3)*sqrtTHat) &
!!$                           *(- BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) + BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
!!$                           * DHat(itheta,izeta) * xPartOfRHS2*gamma  & 
!!$                           * exp(-Z*gamma*Phi1Hat(itheta,izeta)/THat)
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
!!$       ! However, we need not add elements which are 0 due to ddtheta=0 as opposed to because Phi1=0, since such elements will remain 0 at every iteration of SNES.
!!$
!!$       !if (nonlinear .and. (whichMatrix .ne. 2)) then
!!$       !!if (nonlinear .and. (whichMatrix == 1 .or. whichMatrix == 3 .or. (whichMatrix==0 .and. .not. reusePreconditioner))) then !!Commented by AM 2016-02
!!$       !if (.false.) then
!!$       if (includePhi1 .and. includePhi1InKineticEquation .and. (whichMatrix == 1 .or. whichMatrix == 3 .or. (whichMatrix==0 .and. .not. reusePreconditioner))) then !!Added by AM 2016-02
!!$
!!$          !print *,"@@@@@@ ",myRank," max(abs(dphi1hatdtheta)): ",maxval(abs(dphi1hatdtheta)),maxval(abs(dPhi1Hatdzeta))
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
!!$             do itheta=ithetaMin,ithetaMax
!!$
!!$                do izeta=izetaMin,izetaMax
!!$                   factor = -gamma*Zs(ispecies)/(2*BHat(itheta,izeta)*sqrtTHat*sqrtMHat) &
!!$                        * (BHat_sup_theta(itheta,izeta)*dphi1hatdtheta(itheta,izeta) &
!!$                        + BHat_sup_zeta(itheta,izeta)*dPhi1Hatdzeta(itheta,izeta))
!!$
!!$                   !do ix_row=ixMin,Nx
!!$                   do ix_row=max(ixMin,min_x_for_L(L)),Nx
!!$                      rowIndex = getIndex(ispecies,ix_row,L+1,itheta,izeta,BLOCK_F)
!!$
!!$                      ! Term that is super-diagonal in L:
!!$                      if (L<Nxi-1) then
!!$                         ell = L + 1
!!$                         !do ix_col=ixMinCol,Nx
!!$                         do ix_col=max(ixMinCol,min_x_for_L(ell)),Nx
!!$                            if (abs(nonlinearTerm_Lp1(ix_row,ix_col))>threshholdForInclusion) then
!!$                               colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
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
!!$                               colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
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
!!$       ! However, we need not add elements which are 0 due to ddtheta=0 as opposed to because delta f = 0, since such elements will remain 0 at every iteration of SNES.
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
!!$          do itheta = ithetaMin,ithetaMax
!!$             do izeta = izetaMin,izetaMax
!!$                factor = -gamma*Zs(ispecies)/(2*BHat(itheta,izeta)*sqrtTHat*sqrtMHat)
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
!!$                         index = getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
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
!!$                         index = getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
!!$                         ! Add 1 because we are indexing a Fortran array instead of a PETSc object
!!$                         tempVector1(ix) = stateArray(index+1)
!!$                         tempVector2(ix) = tempVector2(ix) + (L+1)*(L+2)/(two*L+3)*stateArray(index+1)/x(ix)
!!$                      end do
!!$                      tempVector2 = tempVector2 + (L+1)/(two*L+3) * matmul(ddx,tempVector1)
!!$                   end if
!!$
!!$                   !do ix=ixMin,Nx
!!$                   do ix=max(ixMin,min_x_for_L(L)),Nx
!!$                      rowIndex = getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
!!$
!!$                      ! Add the d Phi_1 / d theta term:
!!$                      do j=1,Ntheta
!!$                         if (abs(ddtheta(itheta,j))>threshholdForInclusion) then
!!$                            colIndex = getIndex(1,1,1,j,izeta,BLOCK_QN)
!!$                            ! We must use MatSetValue instead of MatSetValueSparse here!!
!!$                            call MatSetValue(matrix, rowIndex, colIndex, &
!!$                                 factor*BHat_sup_theta(itheta,izeta)*ddtheta(itheta,j)*tempVector2(ix), &
!!$                                 ADD_VALUES, ierr)
!!$                         end if
!!$                      end do
!!$
!!$                      ! Add the d Phi_1 / d zeta term:
!!$                      do j=1,Nzeta
!!$                         if (abs(ddzeta(izeta,j))>threshholdForInclusion) then
!!$                            colIndex = getIndex(1,1,1,itheta,j,BLOCK_QN)
!!$                            ! We must use MatSetValue instead of MatSetValueSparse here!!
!!$                            call MatSetValue(matrix, rowIndex, colIndex, &
!!$                                 factor*BHat_sup_zeta(itheta,izeta)*ddzeta(izeta,j)*tempVector2(ix), &
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

    call PetscTime(time4, ierr)
    if (masterProc) print *,"  Time for collisionless terms:",time4-time3
    call PetscTime(time3, ierr)

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
          
          ! For future possible preconditioners, I might want the change the following 2 lines.
          ddx_to_use => ddx
          d2dx2_to_use = d2dx2
          
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
                ! The subtraction in the next line causes a loss of some digits at small x. Is there a better method?
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
                
                ! version 3 normalization:
                !speciesFactor = 3 * nHats(iSpeciesA)  * mHats(iSpeciesA)/mHats(iSpeciesB) &
                !     * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) / T32m
                
                ! alpha_finiteDiffXi normalization
                speciesFactor = 3 * nHats(iSpeciesB)  * sqrt(mHats(iSpeciesB)/THats(iSpeciesB)) &
                     * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) / (THats(iSpeciesB)*mHats(iSpeciesA))
                
                ! WRONG normalization
                !speciesFactor = 3 * nHats(iSpeciesB) * nHats(ispeciesB) * sqrt(mHats(iSpeciesA) * mHats(iSpeciesB)) &
                !     * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) &
                !     / (THats(iSpeciesA)*THats(iSpeciesB)*sqrt(THats(iSpeciesA)*THats(iSpeciesB)))
                
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
                        + speciesFactor * (Psi_Chandra(ix)/x(ix)*d2dx2_to_use(ix,:) &
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

          allocate(collision_operator_xi_block(Nxi,Nxi))
          if (whichMatrix==0) then
             pitch_angle_scattering_operator_to_use => pitch_angle_scattering_operator_preconditioner
          else
             pitch_angle_scattering_operator_to_use => pitch_angle_scattering_operator
          end if

          if (NL>0) then
             allocate(Legendre_projection_to_use(Nxi,Nxi,NL))
             if (whichMatrix==0 .and. preconditioner_field_term_xi_option==1) then
                ! Keep only the part that is diagonal in xi.
                Legendre_projection_to_use = 0
                do ixi = 1,Nxi
                   Legendre_projection_to_use(ixi,ixi,:) = Legendre_projection(ixi,ixi,:)
                end do
             elseif (whichMatrix==0 .and. preconditioner_field_term_xi_option==2) then
                ! Keep only the part that is tridiagonal in xi.
                Legendre_projection_to_use = 0
                do ixi = 1,Nxi ! Handle diagonal
                   Legendre_projection_to_use(ixi,ixi,:) = Legendre_projection(ixi,ixi,:)
                end do
                do ixi = 1,Nxi-1 ! Handle +/- 1 off-diagonal.
                   Legendre_projection_to_use(ixi,ixi+1,:) = Legendre_projection(ixi,ixi+1,:)
                   Legendre_projection_to_use(ixi+1,ixi,:) = Legendre_projection(ixi+1,ixi,:)
                end do
             else
                ! Keep full xi coupling.
                Legendre_projection_to_use = Legendre_projection
             end if
          end if

          do iSpeciesB = 1,Nspecies
             if (whichMatrix==0 .and. preconditioner_species==1) then
                iSpecies_min = iSpeciesB
                iSpecies_max = iSpeciesB
             else
                iSpecies_min = 1
                iSpecies_max = Nspecies
             end if
             do iSpeciesA = iSpecies_min, iSpecies_max
                do ix_row = 1,Nx
                   if (whichMatrix==0.and. preconditioner_x==1) then
                      ix_min = ix_row
                      ix_max = ix_row
                   else
                      ix_min = 1
                      ix_max = Nx
                   end if
                   do ix_col = ix_min, ix_max
                      collision_operator_xi_block = 0
                      do ixi = 1,Nxi
                         ! Add energy scattering, plus the part of the field term that does not depend on the Rosenbluth potentials:
                         ! (These terms are diagonal in xi.)
                         collision_operator_xi_block(ixi,ixi) = CECD(iSpeciesA,iSpeciesB,ix_row,ix_col)
                      end do
                      if (whichMatrix==0 .or. preconditioner_field_term_xi_option==0) then
                         ! Otherwise the field term will be added by apply_dense_terms.F90
                         do L = 0,NL-1
                            ! Add the terms in the field part involving the Rosenbluth potentials:
                            ! (These terms are generally dense in xi.)
                            collision_operator_xi_block = collision_operator_xi_block + Legendre_projection_to_use(:,:,L+1) * RosenbluthPotentialTerms(iSpeciesA,iSpeciesB,L+1,ix_row,ix_col)
                         end do
                      end if
                      if (iSpeciesA==iSpeciesB .and. ix_row==ix_col) then
                         ! Add pitch angle scattering:
                         ! (This operator is usually tri-diagonal or penta-diagonal in xi.)
                         collision_operator_xi_block = collision_operator_xi_block + pitch_angle_scattering_operator_to_use * nuDHat(iSpeciesA,ix_row)
                      end if

                      do itheta = ithetaMin,ithetaMax
                         do izeta = izetaMin,izetaMax
                            !factor = -nu_n*BHat(itheta,izeta)*sqrt(mHats(ispeciesA)/THats(ispeciesA))/abs(DHat(itheta,izeta))
                            factor = -nu_n*spatial_scaling(itheta,izeta)*x_scaling(ix_row,iSpeciesA)
                            do ixi_col = 1,Nxi
                               colIndex = getIndex(iSpeciesB,ix_col,ixi_col,itheta,izeta,BLOCK_F)
                               do ixi_row = 1,Nxi
                                  rowIndex = getIndex(iSpeciesA,ix_row,ixi_row,itheta,izeta,BLOCK_F)
                                  call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                       factor*collision_operator_xi_block(ixi_row,ixi_col), ADD_VALUES, ierr)
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
                      
          if (allocated(Legendre_projection_to_use)) deallocate(Legendre_projection_to_use)
          deallocate(collision_operator_xi_block)
          
          
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

          if (whichMatrix==0) then
             pitch_angle_scattering_operator_to_use => pitch_angle_scattering_operator_preconditioner
          else
             pitch_angle_scattering_operator_to_use => pitch_angle_scattering_operator
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
                ! The subtraction in the next line causes a loss of some digits at small x. Is there a better method?
                Psi_Chandra = (erfs - 2/sqrtpi * xb * expxb2) / (2*xb*xb)
                
                T32m = THats(iSpeciesA) * sqrt(THats(iSpeciesA)*mHats(ispeciesA))
                
                ! Build the pitch-angle scattering frequency:
                nuDHat(iSpeciesA, :) =  nuDHat(iSpeciesA, :) &
                     + (three*sqrtpi/four) / T32m &
                     * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) &
                     * nHats(iSpeciesB)*(erfs - Psi_Chandra)/(x*x*x)
                
             end do

             do ix=ixMin,Nx
                do itheta=ithetaMin,ithetaMax
                   do izeta=izetaMin,izetaMax
                      !factor = -nu_n*BHat(itheta,izeta)*sqrt(mHats(ispeciesA)/THats(ispeciesA))/abs(DHat(itheta,izeta))*nuDHat(iSpeciesA,ix)
                      factor = -nu_n*spatial_scaling(itheta,izeta)*x_scaling(ix,iSpeciesA)*nuDHat(iSpeciesA,ix)
                      do ixi_row = 1,Nxi
                         rowIndex=getIndex(iSpeciesA,ix,ixi_row,itheta,izeta,BLOCK_F)
                         do ixi_col = 1,Nxi
                            colIndex=getIndex(iSpeciesA,ix,ixi_col,itheta,izeta,BLOCK_F)
                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                 factor*pitch_angle_scattering_operator_to_use(ixi_row,ixi_col), ADD_VALUES, ierr)
                         end do
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

    call PetscTime(time4, ierr)
    if (masterProc) print *,"  Time for collision op:",time4-time3
    call PetscTime(time3, ierr)

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
!!$          do itheta = ithetaMin,ithetaMax
!!$             do izeta = izetaMin,izetaMax
!!$                do ispecies = 1,Nspecies
!!$                   index = getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
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
!!$       do itheta = ithetaMin,ithetaMax
!!$          do izeta = izetaMin,izetaMax
!!$             do ispecies = 1,Nspecies
!!$                rowIndex = getIndex(ispecies,ix_row,L+1,itheta,izeta,BLOCK_F)
!!$                do ix_col = 1,Nx
!!$                   colIndex = getIndex(ispecies,ix_col,L+1,itheta,izeta,BLOCK_F)
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
             temp = Ntheta*Nzeta/sum(spatial_scaling)
             do ispecies = 1,Nspecies
                speciesFactor = sqrt(THats(ispecies)/mHats(ispecies)) ! Include 2 when I move to SI units?
                do itheta = ithetaMin,ithetaMax
                   do izeta = izetaMin,izetaMax
                      factor = spatial_scaling(itheta,izeta) * temp * x_scaling(ix,ispecies) * speciesFactor ! This quantity is scaled so it should be O(1)
                      do ixi = 1,Nxi_for_x(ix)
                         rowIndex = getIndex(ispecies, ix, ixi, itheta, izeta, BLOCK_F)
                         
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
          temp = Ntheta*Nzeta/sum(spatial_scaling)
          do ispecies = 1,Nspecies
             speciesFactor = sqrt(THats(ispecies)/mHats(ispecies)) ! Include 2 when I move to SI units?
             do itheta = ithetaMin,ithetaMax
                do izeta = izetaMin,izetaMax
                   do ix = ixMin,Nx
                      factor = spatial_scaling(itheta,izeta) * temp * x_scaling(ix,ispecies) * speciesFactor ! This quantity is scaled so it should be O(1) 
                      do ixi = 1,Nxi_for_x(ix)
                         rowIndex = getIndex(ispecies, ix, ixi, itheta, izeta, BLOCK_F)
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
       
    call PetscTime(time4, ierr)
    if (masterProc) print *,"  Time for sources:",time4-time3
    call PetscTime(time3, ierr)

    ! *******************************************************************************
    ! Add the density and pressure constraints:
    ! *******************************************************************************

    !if (whichMatrix .ne. 2 .and. procThatHandlesConstraints) then
    if (whichMatrix .ne. 2) then
       select case (constraintScheme)
       case (0)
          ! Do nothing here.

       case (1,3,4)
          ! Force the flux-surface-averaged perturbed density and 
          ! flux-surface-averaged perturbed pressure to be zero.

          do ispecies=1,Nspecies
             rowIndex1 = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
             rowIndex2 = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
             !do itheta=1,Ntheta
             do itheta = ithetaMin, ithetaMax
                !do izeta=1,Nzeta
                do izeta = izetaMin, izetaMax
                   do ixi=1,Nxi
                      ! The matrix elements must be proportional to sqrt_g,
                      ! but otherwise the row could probably be scaled any way you like. Here we include a factor 1/VPrimeHat so the 
                      ! row is dimensionless and the row sum is O(1), which seems like a reasonable scaling.
                      factor = xiWeights(ixi)*thetaWeights(itheta)*zetaWeights(izeta)*sqrt_g(itheta,izeta)/VPrimeHat
                      
                      do ix=1,Nx
                         colIndex = getIndex(ispecies, ix, ixi, itheta, izeta, BLOCK_F)
                         
                         !call MatSetValueSparse(matrix, rowIndex1, colIndex, &
                         call MatSetValue(matrix, rowIndex1, colIndex, &
                              x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)
                         
                         !call MatSetValueSparse(matrix, rowIndex2, colIndex, &
                         call MatSetValue(matrix, rowIndex2, colIndex, &
                              x2(ix)*x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)
                      end do
                   end do
                end do
             end do
          end do

       case (2)
          ! Force the flux-surface-averaged distribution function to be zero
          ! at each value of x:

          temp = Ntheta*Nzeta/sum(sqrt_g)
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                do ixi=1,Nxi
                   ! The matrix elements must be proportional to sqrt_g,
                   ! but otherwise the row could probably be scaled any way you like. Here we include a factor 1/VPrimeHat so the 
                   ! row is dimensionless and the row sum is O(1), which seems like a reasonable scaling.
                   select case (constraint_scaling_option)
                   case (1)
                      factor = xiWeights(ixi)*thetaWeights(itheta)*zetaWeights(izeta)*sqrt_g(itheta,izeta)/VPrimeHat
                   case (2)
                      factor = sqrt_g(itheta,izeta)/VPrimeHat
                   case (3)
                      factor = sqrt_g(itheta,izeta) * temp
                   case default
                      print *,"Invalid constraint_scaling_option:",constraint_scaling_option
                   end select

                   do ix=ixMin,Nx
                      do ispecies = 1,Nspecies
                         colIndex = getIndex(ispecies, ix, ixi, itheta, izeta, BLOCK_F)
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

    call PetscTime(time4, ierr)
    if (masterProc) print *,"  Time for constraints:",time4-time3
    call PetscTime(time3, ierr)

    ! *******************************************************************************
    ! SECTION MODIFIED BY AM 2016-02/03
    ! Add the quasineutrality equation
    ! *******************************************************************************

    if (whichMatrix .ne. 2 .and. includePhi1) then
       L=0
       do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
             rowIndex = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN)

             ! Add the charge density of the kinetic species (the part of the residual/Jacobian related to f1):
             do ispecies = 1,Nspecies
                !!speciesFactor = Zs(ispecies)*THats(ispecies)/mHats(ispecies) & !!Commented by AM 2016-02
                !!     * sqrt(THats(ispecies)/mHats(ispecies)) !!Commented by AM 2016-02
                speciesFactor = four*pi*Zs(ispecies)*THats(ispecies)/mHats(ispecies) & !!Added by AM 2016-02 !!The factor 4pi comes from velocity integration over xi.
                     * sqrt(THats(ispecies)/mHats(ispecies)) !!Added by AM 2016-02

                do ix = 1,Nx
                   colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
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
             colIndex = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN)
             if (quasineutralityOption == 1 .and. (whichMatrix == 0 .or. whichMatrix == 1)) then !!Only add in Jacobian matrix
                do ispecies = 1,Nspecies
                   call MatSetValueSparse(matrix, rowIndex, colIndex, &
                        - gamma * Zs(ispecies) * Zs(ispecies) * NHats(ispecies) * exp (- Zs(ispecies)* gamma * Phi1Hat(itheta,izeta) / THats(ispecies)) / THats(ispecies), ADD_VALUES, ierr)
                end do
                if (withAdiabatic) then
                   call MatSetValueSparse(matrix, rowIndex, colIndex, &
                        - gamma * adiabaticZ * adiabaticZ *adiabaticNHat * exp (- adiabaticZ* gamma * Phi1Hat(itheta,izeta) / adiabaticTHat) / adiabaticTHat, ADD_VALUES, ierr)
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
       do itheta=1,Ntheta
          do izeta=1,Nzeta
             colIndices(izeta) = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN)
          end do

          call MatSetValuesSparse(matrix, 1, rowIndices(1), Nzeta, colIndices, &
               thetaWeights(itheta)*zetaWeights*sqrt_g(itheta,:), ADD_VALUES, ierr)
       end do

       deallocate(rowIndices)
       deallocate(colIndices)
    end if

    
    ! *******************************************************************************
    ! If we are going to use a fieldsplit preconditioner, then put 1's on the
    ! diagonal for the source/constraint blocks of the preconditioner.
    ! *******************************************************************************

    if (whichMatrix==0 .and. procThatHandlesConstraints .and. fieldsplit) then
       select case (constraintScheme)
       case (0)
          ! Do nothing here.
       case (1,3,4)
          print *,"Adding 1s on the source/constraint diagonal block for fieldsplit."
          do ispecies=1,Nspecies
             rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
             call MatSetValue(matrix, rowIndex, rowIndex, one, ADD_VALUES, ierr)
             rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
             call MatSetValue(matrix, rowIndex, rowIndex, one, ADD_VALUES, ierr)
          end do
       case (2)
          print *,"Adding 1s on the source/constraint diagonal block for fieldsplit."
          do ix=ixMin,Nx
             do ispecies = 1,Nspecies
                rowIndex = getIndex(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT)
                call MatSetValue(matrix, rowIndex, rowIndex, one, ADD_VALUES, ierr)
             end do
          end do
       end select
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

    call PetscTime(time4, ierr)
    if (masterProc) print *,"  Time for remaining stuff:",time4-time3

  end subroutine populateMatrix

! -----------------------------------------------------------------------
! -----------------------------------------------------------------------


  subroutine init_f0()

    use globalVariables
    use petscmat
    use indices

    implicit none

    integer :: ix, ixi, itheta, izeta, ispecies, index
    PetscScalar :: factor
    PetscErrorCode :: ierr

    if (masterProc) then
       print *,"Initializing f0"
    end if

    call VecSet(f0, zero, ierr)
    
    do ispecies = 1,Nspecies
       do ix = 1,Nx
          factor = expx2(ix)
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                do ixi = 1,Nxi_for_x(ix)
                   index = getIndex(ispecies, ix, ixi, itheta, izeta, BLOCK_F)
                   call VecSetValue(f0, index, &
                        exp(-Zs(ispecies)*gamma*Phi1Hat(itheta,izeta)/THats(ispecies))*factor, INSERT_VALUES, ierr) ! This line needs fixing when Phi1 is included.
                end do
             end do
          end do
       end do
    end do

    call VecAssemblyBegin(f0, ierr)
    call VecAssemblyEnd(f0, ierr)

  end subroutine init_f0
