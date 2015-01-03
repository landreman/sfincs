! For compilers that do not include the error function erf(x), the line
! below should be un-commented:
!#define USE_GSL_ERF
  
#include <finclude/petscmatdef.h>
#include "PETScVersions.F90"


  subroutine populateMatrix(matrix, whichMatrix)

    use petscmat
    use globalVariables
    use sparsify
    use indices

    implicit none

    Mat :: matrix
    integer, intent(in) :: whichMatrix
    ! Allowed values for whichMatrix:
    ! 0 = preconditioner for Jacobian
    ! 1 = Jacobian
    ! 2 = matrix used to evaluate the residual

    PetscErrorCode :: ierr
    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, speciesFactor, speciesFactor2
    PetscScalar :: T32m, factor, LFactor, temp, temp1, temp2, xDotFactor
    PetscScalar, dimension(:), allocatable :: xb, expxb2
    PetscScalar, dimension(:,:), allocatable :: thetaPartOfTerm, localThetaPartOfTerm, xPartOfXDot
    integer :: i, j, ix, ispecies, itheta, izeta, L, ixi, index
    integer :: rowIndex, colIndex
    integer :: ell, iSpeciesA, iSpeciesB
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar, dimension(:,:), allocatable :: ddxToUse, d2dx2ToUse, zetaPartOfTerm, localZetaPartOfTerm
    PetscScalar, dimension(:,:), allocatable :: fToFInterpolationMatrix
    PetscScalar, dimension(:,:), allocatable :: potentialsToFInterpolationMatrix
    PetscScalar, dimension(:,:,:,:), allocatable :: CECD
    PetscScalar :: xPartOfSource1, xPartOfSource2
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
    PetscScalar :: CHat_element
    character(len=50) :: whichMatrixName


    ! *******************************************************************************
    ! Do a few sundry initialization tasks:
    ! *******************************************************************************

    select case (whichMatrix)
    case (0)
       ! Preconditioner for Jacobian
       whichMatrixName = "Jacobian preconditioner"
    case (1)
       ! Jacobian
       whichMatrixName = "Jacobian"
    case (2)
       ! Matrix used to evaluate residual
       whichMatrixName = "residual"
    case default
       if (masterProc) then
          print *,"Error! whichMatrix must be 0, 1, or 2."
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
    ! (unlike mumps or superlu_dist), then if we're running with just 1 proc,
    ! add some values to the diagonals of the preconditioner.
    if (masterProc .and. numProcs==1) then
       ! Value to set:
       temp = 1d+0
       do i=Nspecies*Nx*Nxi*Ntheta*Nzeta, matrixSize-1
          call MatSetValue(matrix, i, i, temp, ADD_VALUES, ierr)
       end do
    end if


    ! *********************************************************
    ! Allocate small matrices:
    ! *********************************************************

    allocate(ddxToUse(Nx,Nx))
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


    ! *********************************************************
    ! Select appropriate differentiation matrices depending on
    ! whether this is the preconditioner or the final matrix:
    ! *********************************************************

    if (whichMatrix>0 .or. preconditioner_theta==0) then
       ddthetaToUse = ddtheta
    else
       ddthetaToUse = ddtheta_preconditioner
    end if

    if (whichMatrix>0 .or. preconditioner_zeta==0) then
       ddzetaToUse = ddzeta
    else
       ddzetaToUse = ddzeta_preconditioner
    end if

    do ispecies = 1,Nspecies
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)

       ! *********************************************************
       ! Add the streaming d/dtheta term:
       ! *********************************************************

       allocate(thetaPartOfTerm(Ntheta,Ntheta))
       allocate(localThetaPartOfTerm(Ntheta,localNtheta))
       allocate(rowIndices(localNtheta))
       allocate(colIndices(Ntheta))
       do izeta=izetaMin,izetaMax
          do itheta=1,Ntheta
             !thetaPartOfTerm(itheta,:) = iota*sqrtTHat/sqrtMHat * ddthetaToUse(itheta,:) &
             thetaPartOfTerm(itheta,:) = BHat_sup_theta(itheta,izeta) &
                  * sqrtTHat/sqrtMHat * ddthetaToUse(itheta,:) &
                  / BHat(itheta,izeta)
          end do

          ! PETSc uses the opposite convention to Fortran:
          thetaPartOfTerm = transpose(thetaPartOfTerm)
          localThetaPartOfTerm = thetaPartOfTerm(:,ithetaMin:ithetaMax)

          do ix=1,Nx
             do L=0,(Nxi-1)
                do itheta=1,localNtheta
                   rowIndices(itheta) = getIndex(ispecies, ix, L+1, ithetaMin+itheta-1, izeta, 0)
                end do

                ! Super-diagonal term
                if (L < Nxi-1) then
                   ell = L+1
                   do itheta=1,Ntheta
                      colIndices(itheta) = getIndex(ispecies, ix, ell+1, itheta, izeta, 0)
                   end do

                   call MatSetValuesSparse(matrix, localNtheta, rowIndices, Ntheta, colIndices, &
                        (L+1)/(2*L+three)*x(ix)*localThetaPartOfTerm, ADD_VALUES, ierr)
                end if

                ! Sub-diagonal term
                if (L > 0) then
                   ell = L-1
                   do itheta=1,Ntheta
                      colIndices(itheta) = getIndex(ispecies, ix, ell+1, itheta, izeta, 0)
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

       ! *********************************************************
       ! Add the streaming d/dzeta term:
       ! *********************************************************

       allocate(zetaPartOfTerm(Nzeta,Nzeta))
       allocate(localZetaPartOfTerm(Nzeta,localNzeta))
       allocate(rowIndices(localNzeta))
       allocate(colIndices(Nzeta))
       do itheta=ithetaMin, ithetaMax
          do izeta=1,Nzeta
             !zetaPartOfTerm(izeta,:) = sqrtTHat/sqrtMHat * ddzetaToUse(izeta,:) / BHat(itheta,izeta)
             zetaPartOfTerm(izeta,:) = sqrtTHat/sqrtMHat * BHat_sup_zeta(itheta,izeta) &
                  * ddzetaToUse(izeta,:) / BHat(itheta,izeta)
          end do

          ! PETSc uses the opposite convention to Fortran:
          zetaPartOfTerm = transpose(zetaPartOfTerm)
          localZetaPartOfTerm = zetaPartOfTerm(:,izetaMin:izetaMax)

          do ix=1,Nx
             do L=0,(Nxi-1)
                do izeta = 1,localNzeta
                   rowIndices(izeta)=getIndex(ispecies, ix, L+1, itheta, izetaMin+izeta-1, 0)
                end do

                ! Super-diagonal term
                if (L < Nxi-1) then
                   ell = L + 1
                   do izeta = 1,Nzeta
                      colIndices(izeta)=getIndex(ispecies, ix, ell+1, itheta, izeta, 0)
                   end do

                   call MatSetValuesSparse(matrix, localNzeta, rowIndices, Nzeta, colIndices, &
                        (L+1)/(2*L+three)*x(ix)*localZetaPartOfTerm, ADD_VALUES, ierr)
                end if

                ! Sub-diagonal term
                if (L > 0) then
                   ell = L - 1
                   do izeta = 1,Nzeta
                      colIndices(izeta)=getIndex(ispecies, ix, ell+1, itheta, izeta, 0)
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

       ! *********************************************************
       ! Add the ExB d/dtheta term:
       ! *********************************************************

       !factor = alpha*Delta*GHat/(2*psiAHat)*dPhiHatdpsiN
       factor = alpha*Delta/two*dPhiHatdpsiHat
       allocate(thetaPartOfTerm(Ntheta,Ntheta))
       allocate(localThetaPartOfTerm(Ntheta,localNtheta))
       allocate(rowIndices(localNtheta))
       allocate(colIndices(Ntheta))
       do izeta=izetaMin,izetaMax
          if (useDKESExBDrift) then
             !thetaPartOfTerm = ddthetaToUse / FSABHat2
             do itheta=1,Ntheta
                thetaPartOfTerm(itheta,:) = ddthetaToUse(itheta,:) / FSABHat2 &
                     * DHat(itheta,izeta) * BHat_sub_zeta(itheta,izeta)
             end do
          else
             do itheta=1,Ntheta
                !thetaPartOfTerm(itheta,:) = ddthetaToUse(itheta,:) / (BHat(itheta,izeta) ** 2)
                thetaPartOfTerm(itheta,:) = ddthetaToUse(itheta,:) / (BHat(itheta,izeta) ** 2) &
                     * DHat(itheta,izeta) * BHat_sub_zeta(itheta,izeta)
             end do
          end if

          ! PETSc uses the opposite convention to Fortran:
          thetaPartOfTerm = transpose(thetaPartOfTerm*factor)
          localThetaPartOfTerm = thetaPartOfTerm(:,ithetaMin:ithetaMax)

          do ix=1,Nx
             do L=0,(Nxi-1)
                do itheta=1,localNtheta
                   rowIndices(itheta)=getIndex(ispecies,ix,L+1,itheta+ithetaMin-1,izeta,0)
                end do
                do itheta=1,Ntheta
                   colIndices(itheta)=getIndex(ispecies,ix,L+1,itheta,izeta,0)
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

       ! *********************************************************
       ! Add the ExB d/dzeta term:
       ! *********************************************************

       !factor = -alpha*Delta*IHat/(2*psiAHat)*dPhiHatdpsiN
       factor = -alpha*Delta/two*dPhiHatdpsiHat
       allocate(zetaPartOfTerm(Nzeta,Nzeta))
       allocate(localZetaPartOfTerm(Nzeta,localNzeta))
       allocate(rowIndices(localNzeta))
       allocate(colIndices(Nzeta))
       do itheta=ithetaMin, ithetaMax
          if (useDKESExBDrift) then
             !zetaPartOfTerm = ddzetaToUse / FSABHat2
             do izeta=1,Nzeta
                zetaPartOfTerm(izeta,:) = ddzetaToUse(izeta,:) / FSABHat2 &
                     * DHat(itheta,izeta) * BHat_sub_theta(itheta,izeta)
             end do
          else
             do izeta=1,Nzeta
                !zetaPartOfTerm(izeta,:) = ddzetaToUse(izeta,:) / (BHat(itheta,izeta) ** 2)
                zetaPartOfTerm(izeta,:) = ddzetaToUse(izeta,:) / (BHat(itheta,izeta) ** 2) &
                     * DHat(itheta,izeta) * BHat_sub_theta(itheta,izeta)
             end do
          end if

          ! PETSc uses the opposite convention to Fortran:
          zetaPartOfTerm = transpose(zetaPartOfTerm*factor)
          localzetaPartOfTerm = zetaPartOfTerm(:,izetaMin:izetaMax)

          do ix=1,Nx
             do L=0,(Nxi-1)
                do izeta=1,localNzeta
                   rowIndices(izeta)=getIndex(ispecies,ix,L+1,itheta,izeta+izetaMin-1,0)
                end do
                do izeta=1,Nzeta
                   colIndices(izeta)=getIndex(ispecies,ix,L+1,itheta,izeta,0)
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

       ! *********************************************************
       ! Add the standard mirror term:
       ! *********************************************************

       do itheta=ithetaMin,ithetaMax
          do izeta=izetaMin,izetaMax
             factor = -sqrtTHat/(2*sqrtMHat*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                  * (BHat_sup_theta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                  + BHat_sup_zeta(itheta,izeta) * dBHatdzeta(itheta,izeta))
                 !* (iota*dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta))

             do ix=1,Nx
                do L=0,(Nxi-1)
                   rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,0)

                   if (L<Nxi-1) then
                      ! Super-diagonal term:
                      ell = L+1
                      colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,0)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           (L+1)*(L+2)/(2*L+three)*x(ix)*factor, ADD_VALUES, ierr)
                   end if

                   if (L>0) then
                      ! Sub-diagonal term:
                      ell = L-1
                      colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,0)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           -L*(L-1)/(2*L-one)*x(ix)*factor, ADD_VALUES, ierr)
                   end if
                end do
             end do
          end do
       end do

       ! *********************************************************
       ! Add the non-standard d/dxi term:
       ! *********************************************************

       if (includeElectricFieldTermInXiDot) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax
                !factor = alpha*Delta*dPhiHatdpsiN/(4*psiAHat*(BHat(itheta,izeta)**3)) &
                !     * (GHat*dBHatdtheta(itheta,izeta) - IHat* dBHatdzeta(itheta,izeta))

                temp = BHat_sub_zeta(itheta,izeta) * dBHatdtheta(itheta,izeta) &
                     - BHat_sub_theta(itheta,izeta) * dBHatdzeta(itheta,izeta)

                if (.not. force0RadialCurrentInEquilibrium) then
                   temp = temp - 2 * BHat(itheta,izeta) &
                        * (dBHat_sub_zeta_dtheta(itheta,izeta) - dBHat_sub_theta_dzeta(itheta,izeta))
                end if

                factor = alpha*Delta*dPhiHatdpsiHat/(4*(BHat(itheta,izeta)**3)) &
                     * DHat(itheta,izeta) * temp
                     

                do ix=1,Nx
                   do L=0,(Nxi-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,0)

                      ! Diagonal term
                      call MatSetValueSparse(matrix, rowIndex, rowIndex, &
                           (L+1)*L/((2*L-one)*(2*L+three))*factor, ADD_VALUES, ierr)

                      if (whichMatrix>0 .or. preconditioner_xi==0) then
                         if (L<Nxi-2) then
                            ! Super-super-diagonal term:
                            ell = L+2
                            colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,0)
                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                 (L+3)*(L+2)*(L+1)/((two*L+5)*(2*L+three))*factor, ADD_VALUES, ierr)
                         end if

                         if (L>1) then
                            ! Sub-sub-diagonal term:
                            ell = L-2
                            colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,0)
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
       ! Add the collisionless d/dx term:
       ! *********************************************************

       if (includeXDotTerm) then

          allocate(xPartOfXDot(Nx,Nx))
          allocate(rowIndices(Nx))
          allocate(colIndices(Nx))
          factor = alpha*Delta/(4*psiAHat)*dPhiHatdpsiN

          do L=0,(Nxi-1)
             if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                ddxToUse = ddx_preconditioner
             else
                ddxToUse = ddx
             end if
             do ix=1,Nx
                xPartOfXDot(ix,:) = x(ix) * ddxToUse(ix,:)
             end do
             xPartOfXDot = transpose(xPartOfXDot)  ! PETSc uses the opposite convention of Fortran

             do itheta=ithetaMin,ithetaMax

                do izeta=izetaMin,izetaMax
                   xDotFactor = factor/(BHat(itheta,izeta)**3) &
                        * (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))

                   do ix=1,Nx
                      rowIndices(ix)=getIndex(ispecies,ix,L+1,itheta,izeta,0)
                   end do

                   ! Term that is diagonal in L:
                   colIndices = rowIndices
                   LFactor = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*xDotFactor
                   call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                        LFactor*xPartOfXDot, ADD_VALUES, ierr)

                   if (whichMatrix>0 .or. preconditioner_xi==0) then
                      ! Term that is super-super-diagonal in L:
                      if (L<(Nxi-2)) then
                         ell = L + 2
                         do ix=1,Nx
                            colIndices(ix)=getIndex(ispecies,ix,ell+1,itheta,izeta,0)
                         end do
                         LFactor = (L+1)*(L+2)/((two*L+5)*(2*L+3))*xDotFactor
                         call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                              LFactor*xPartOfXDot, ADD_VALUES, ierr)
                      end if

                      ! Term that is sub-sub-diagonal in L:
                      if (L>1) then
                         ell = L - 2
                         do ix=1,Nx
                            colIndices(ix)=getIndex(ispecies,ix,ell+1,itheta,izeta,0)
                         end do
                         LFactor = L*(L-1)/((two*L-3)*(2*L-1))*xDotFactor
                         call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                              LFactor*xPartOfXDot, ADD_VALUES, ierr)
                      end if
                   end if

                end do
             end do
          end do
          deallocate(rowIndices)
          deallocate(colIndices)
          deallocate(xPartOfXDot)
       end if

       ! *********************************************************
       ! Add the optional term which does not involve derivatives of f,
       ! Sometimes useful for restoring Liouville's theorem:
       ! *********************************************************

       if (include_fDivVE_term) then
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                do ix = 1,Nx
                   do ixi = 1,Nxi
                      index=getIndex(ispecies,ix,ixi,itheta,izeta,0)
                      call MatSetValueSparse(matrix,index,index, &
                           -dPhiHatdpsiN*Delta*alpha/(psiAHat*(BHat(itheta,izeta)**3)) &
                           *(GHat*dBHatdtheta(itheta,izeta)- IHat*dBHatdzeta(itheta,izeta)), &
                           ADD_VALUES, ierr)
                   end do
                end do
             end do
          end do
       end if


    end do

    ! End of adding the collisionless kinetic terms

    ! *********************************************************
    ! *********************************************************
    !
    ! Next, we add the collision operator.
    !
    ! *********************************************************
    ! *********************************************************


    select case (collisionOperator)

    case (0)
       ! *********************************************************
       ! Full linearized Fokker-Planck operator
       ! *********************************************************


       ! *********************************************************
       ! In preparation for adding the collision operator,
       ! create several matrices which will be needed.
       ! *********************************************************

       allocate(rowIndices(Nx))
       allocate(colIndices(Nx))
       allocate(tempMatrix(Nx, NxPotentials))
       allocate(tempMatrix2(NxPotentials, NxPotentials))
       allocate(extrapMatrix(Nx, NxPotentials))

       ! For future possible preconditioners, I might want the change the following 2 lines.
       ddxToUse = ddx
       d2dx2ToUse = d2dx2

       ! First assemble rows 2 and 3 of the block linear system, since they
       ! are independent of psi and independent of species.

       M32 = zero
       M21 = 4*pi*regridPolynomialToUniform
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
             ! by this regridding matrix to obtain its values on the species-A grid:
             if (iSpeciesA /= iSpeciesB) then
                call polynomialInterpolationMatrix(Nx, Nx, x, xb, expx2, &
                     expxb2, fToFInterpolationMatrix)
             else
                fToFInterpolationMatrix = zero
                do i=1,Nx
                   fToFInterpolationMatrix(i, i) = one
                end do
             end if

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

                   if (L < NL) then
                      !   if (.false.) then
                      ! Add Rosenbluth potential terms.

                      speciesFactor2 = sqrt(THats(iSpeciesA)*mHats(iSpeciesB) &
                           / (THats(iSpeciesB) * mHats(iSpeciesA)))

                      ! Build M13:
                      call interpolationMatrix(NxPotentials, Nx, xPotentials, x*speciesFactor2, &
                           potentialsToFInterpolationMatrix, extrapMatrix)

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

                   ! PETSc and Fortran use row-major vs column-major:
                   CHat = transpose(CHat)

                   ! At this point, CHat contains the collision operator normalized by
                   ! \bar{nu}, (the collision frequency at the reference mass, density, and temperature.)

                   do itheta=ithetaMin,ithetaMax
                      do izeta=izetaMin,izetaMax
                         do ix=1,Nx
                            rowIndices(ix)=getIndex(iSpeciesA,ix,L+1,itheta,izeta,0)
                            colIndices(ix)=getIndex(iSpeciesB,ix,L+1,itheta,izeta,0)
                         end do
                         call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                              -nu_n*CHat, ADD_VALUES, ierr)
                              !-nu_n*(GHat+iota*IHat)/(BHat(itheta,izeta)*BHat(itheta,izeta))*CHat, &
                              ! ADD_VALUES, ierr)
                      end do
                   end do

                end if
             end do
          end do
       end do

       deallocate(rowIndices)
       deallocate(colIndices)
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
             do ix=1,Nx
                CHat_element = -oneHalf*nuDHat(iSpeciesA,ix)*L*(L+1)

                ! At this point, CHat contains the collision operator normalized by
                ! \bar{nu}, (the collision frequency at the reference mass, density, and temperature.)

                do itheta=ithetaMin,ithetaMax
                   do izeta=izetaMin,izetaMax
                      index=getIndex(iSpeciesA,ix,L+1,itheta,izeta,0)
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

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Done adding the collision operator.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    ! *******************************************************************************
    ! Add sources:
    ! *******************************************************************************

    select case (constraintScheme)
    case (0)
       ! Do nothing here.

    case (1)
       ! Add a heat source and a particle source.

       L=0
       do ix=1,Nx
          xPartOfSource1 = (x2(ix)-5/two)*exp(-x2(ix)) ! Provides particles but no heat
          xPartOfSource2 = (x2(ix)-3/two)*exp(-x2(ix)) ! Provides heat but no particles
          do itheta=ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                do ispecies = 1,Nspecies
                   rowIndex = getIndex(ispecies, ix, L+1, itheta, izeta, 0)

                   colIndex = getIndex(ispecies, 1, 1, 1, 1, 1)
                   !call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource1 / (BHat(itheta,izeta) ** 2), &
                   !     ADD_VALUES, ierr)
                   call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource1, ADD_VALUES, ierr)

                   colIndex = getIndex(ispecies, 1, 1, 1, 1, 2)
                   !call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource2 / (BHat(itheta,izeta) ** 2), &
                   !     ADD_VALUES, ierr)
                   call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource2, ADD_VALUES, ierr)
                end do
             end do
          end do
       end do

    case (2)
       ! Add a L=0 source (which is constant on the flux surface) at each x.
       L=0
       do ix=1,Nx
          do itheta=ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                do ispecies = 1,Nspecies
                   rowIndex = getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                   colIndex = getIndex(ispecies, ix, 1, 1, 1, 3)
                   !call MatSetValue(matrix, rowIndex, colIndex, one / (BHat(itheta,izeta) ** 2), &
                   !     ADD_VALUES, ierr)
                   call MatSetValue(matrix, rowIndex, colIndex, one, ADD_VALUES, ierr)
                end do
             end do
          end do
       end do

    case default
       print *,"Error! Invalid constraintScheme."
       stop
    end select


    ! *******************************************************************************
    ! Add constraints:
    ! *******************************************************************************

    if (procThatHandlesConstraints) then
       select case (constraintScheme)
       case (0)
          ! Do nothing here.

       case (1)
          ! Force the flux-surface-averaged perturbed density and 
          ! flux-surface-averaged perturbed pressure to be zero.

          L=0
          do itheta=1,Ntheta
             do izeta=1,Nzeta
                factor = 1/(BHat(itheta,izeta) * BHat(itheta,izeta))

                do ix=1,Nx
                   do ispecies=1,Nspecies
                      colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, 0)

                      rowIndex = getIndex(ispecies, 1, 1, 1, 1, 1)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)

                      rowIndex = getIndex(ispecies, 1, 1, 1, 1, 2)
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
                factor = 1/(BHat(itheta,izeta) * BHat(itheta,izeta))
                do ix=1,Nx
                   do ispecies = 1,Nspecies
                      colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                      rowIndex = getIndex(ispecies, ix, 1, 1, 1, 3)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           factor, ADD_VALUES, ierr)
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


  end subroutine populateMatrix

