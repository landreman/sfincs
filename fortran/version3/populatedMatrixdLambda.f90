! See populateMatrix.f90 (where matrix is constructed). Much of the code has been copied.

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

  subroutine populatdMatrixdLambda(dMatrixdLambda, whichLambda, whichMode)

    use petscmat
    use globalVariables
    use sparsify
    use indices
    use xGrid, only: xGrid_k

    implicit none

    Mat :: dMatrixdLambda
    integer :: whichLambda, whichMode
    integer :: i, ixMin, ixMax, L, ispecies, izeta, itheta, ix, ixMinCol
    PetscScalar, dimension(:,:), allocatable :: dThetaPartOfTermdLambda, localdThetaPartOfTermdLambda, &
      dZetaPartOfTermdLambda, localdZetaPartOfTermdLambda
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor
    PetscScalar, dimension(:,:), allocatable :: ddthetaToUse, ddzetaToUse, ddxToUse, ddxToUse_plus, &
      ddxToUse_minus
    PetscScalar :: factor, temp, xDotFactor, stuffToAdd
    PetscScalar :: xPartOfXDot, xPartOfXDot_plus, xPartOfXDot_minus
    integer :: ix_row, ix_col, rowIndex, colIndex, ell
    PetscLogDouble :: time1, time2
    double precision :: myMatInfo(MAT_INFO_SIZE)
    integer :: NNZ, NNZAllocated, NMallocs

    ! Sometimes PETSc complains if any of the diagonal elements are not set.
    ! Therefore, set the entire diagonal to 0 to be safe.
    if (masterProc) then
      do i=1,matrixSize
        call MatSetValue(dMatrixdLambda, i-1, i-1, zero, ADD_VALUES, ierr)
      end do
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
    allocate(ddthetaToUse(Ntheta, Ntheta))
    allocate(ddzetaToUse(Nzeta, Nzeta))

    do ispecies = 1,Nspecies
      nHat = nHats(ispecies)
      THat = THats(ispecies)
      mHat = mHats(ispecies)
      Z = Zs(ispecies)
      sqrtTHat = sqrt(THat)
      sqrtMHat = sqrt(mHat)

      ! *********************************************************
      ! Add the sensitivity of the streaming d/dtheta term:
      ! *********************************************************

      allocate(dThetaPartOfTermdLambda(Ntheta,Ntheta))
      allocate(dThetaPartOfTermdLambda(Ntheta,localNtheta))
      allocate(rowIndices(localNtheta))
      allocate(colIndices(Ntheta))

      do L=0,(Nxi-1)
        do izeta=izetaMin,izetaMax
          do itheta=1,Ntheta

            !thetaPartOfTerm(itheta,:) = BHat_sup_theta(itheta,izeta) &
            !   * sqrtTHat/sqrtMHat * ddthetaToUse(itheta,:) &
            !   / BHat(itheta,izeta)

            select case(whichLambda)
              case (-1) ! Er
                if (masterProc) then
                  print *,"Error! Er sensitivity not yet implemented."
                end if
                stop
              case (0) ! BHat
                dThetaPartOfTermdLambda(itheta,:) = - BHat_sup_theta(itheta,izeta) &
                  * sqrtTHat/sqrtMHat * ddtheta(itheta,:) * dBHatdFourier(itheta,izeta,whichMode) &
                  / (BHat(itheta,izeta)*BHat(itheta,izeta))
              case (1) ! BHat_sup_theta
                dThetaPartOfTermdLambda(itheta,:) = sqrtTHat/sqrtMHat * ddtheta(itheta,:) &
                  * dBHat_sup_thetadFourier / BHat(itheta,izeta)
              case (2) ! BHat_sup_zeta
                dThetaPartOfTermdLambda(itheta,:) = 0
              case (3) ! BHat_sub_theta
                dThetaPartOfTermdLambda(itheta,:) = 0
              case (4) ! BHat_sub_zeta
                dThetaPartOfTermdLambda(itheta,:) = 0
              case (5) ! DHat
                dThetaPartOfTermdLambda(itheta,:) = 0
            end select

            ! PETSc uses the opposite convention to Fortran:
            dThetaPartOfTermdLambda = transpose(dThetaPartOfTermdLambda)
            localdThetaPartOfTermdLambda = dThetaPartOfTermdLambda(:,ithetaMin:ithetaMax)

            !do ix=ixMin,Nx
            do ix=max(ixMin,min_x_for_L(L)),Nx

              do itheta=1,localNtheta
                rowIndices(itheta) = getIndex(ispecies, ix, L+1, ithetaMin+itheta-1, izeta, BLOCK_F)
              end do

          ! Super-diagonal-in-L term
              if (L < Nxi_for_x(ix)-1) then
                ell = L+1
                do itheta=1,Ntheta
                  colIndices(itheta) = getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                end do

                call MatSetValuesSparse(dMatrixdLambda, localNtheta, rowIndices, Ntheta, colIndices, &
                  (L+1)/(2*L+three)*x(ix)*localdThetaPartOfTermdLambda, ADD_VALUES, ierr)
              end if

              ! Sub-diagonal-in-L term
              if (L > 0) then
                ell = L-1
                do itheta=1,Ntheta
                  colIndices(itheta) = getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                end do

                call MatSetValuesSparse(dMatrixdLambda, localNtheta, rowIndices, Ntheta, colIndices, &
                  L/(2*L-one)*x(ix)*dThetaPartOfTermdLambda, ADD_VALUES, ierr)
              end if

            end do
          end do
        end do
      end do

      deallocate(rowIndices)
      deallocate(colIndices)
      deallocate(dThetaPartOfTermdLambda)
      deallocate(localdThetaPartOfTermdLambda)

      ! *********************************************************
      ! Add the sensitivity of the streaming d/dzeta term:
      ! *********************************************************

      allocate(dZetaPartOfTermdLambda(Nzeta,Nzeta))
      allocate(localdZetaPartOfTermdLambda(Nzeta,localNzeta))
      allocate(rowIndices(localNzeta))
      allocate(colIndices(Nzeta))

      do L=0,(Nxi-1)
        do itheta=ithetaMin, ithetaMax

          do izeta=1,Nzeta
            !zetaPartOfTerm(izeta,:) = sqrtTHat/sqrtMHat * BHat_sup_zeta(itheta,izeta) &
            !  * ddzetaToUse(izeta,:) / BHat(itheta,izeta)

            select case(whichLambda)
              case (-1) ! Er
                if (masterProc) then
                  print *,"Error! Er sensitivity not yet implemented."
                end if
                stop
              case (0) ! BHat
                dZetaPartOfTermdLambda(izeta,:) = - BHat_sup_zeta(itheta,izeta) &
                  * sqrtTHat/sqrtMHat * ddzeta(izeta,:) * dBHatdFourier(itheta,izeta,whichMode) &
                  / (BHat(itheta,izeta)*BHat(itheta,izeta))
              case (1) ! BHat_sup_theta
                dZetaPartOfTermdLambda(izeta,:) = 0
              case (2) ! BHat_sup_zeta
                dZetaPartOfTermdLambda(izeta,:) = sqrtTHat/sqrtMHat * ddzeta(izeta,:) &
                  * dBHat_sup_zetadFourier(itheta,izeta,whichMode)/ BHat(itheta,izeta)
              case (3) ! BHat_sub_theta
                dZetaPartOfTermdLambda(izeta,:) = 0
              case (4) ! BHat_sub_zeta
                dZetaPartOfTermdLambda(izeta,:) = 0
              case (5) ! DHat
                dZetaPartOfTermdLambda(izeta,:) = 0
            end select

          end do

          ! PETSc uses the opposite convention to Fortran:
          dZetaPartOfTermdLambda = transpose(dZetaPartOfTermdLambda)
          localdZetaPartOfTermdLambda = dZetaPartOfTermdLambda(:,izetaMin:izetaMax)

          do ix=max(ixMin,min_x_for_L(L)),Nx
             do izeta = 1,localNzeta
                rowIndices(izeta)=getIndex(ispecies, ix, L+1, itheta, izetaMin+izeta-1, BLOCK_F)
             end do
             
             ! Super-diagonal-in-L term
             if (L < Nxi_for_x(ix)-1) then
                ell = L + 1
                do izeta = 1,Nzeta
                   colIndices(izeta)=getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                end do
                
                call MatSetValuesSparse(dMatrixdLambda, localNzeta, rowIndices, Nzeta, colIndices, &
                     (L+1)/(2*L+three)*x(ix)*localdZetaPartOfTermdLambda, ADD_VALUES, ierr)
             end if
             
             ! Sub-diagonal-in-L term
             if (L > 0) then
                ell = L - 1
                do izeta = 1,Nzeta
                   colIndices(izeta)=getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                end do
                
                call MatSetValuesSparse(matrix, localNzeta, rowIndices, Nzeta, colIndices, &
                     L/(2*L-one)*x(ix)*localdZetaPartOfTermdLambda, ADD_VALUES, ierr)
             end if

            end do
          end do
        end do
        deallocate(rowIndices)
        deallocate(colIndices)
        deallocate(dZetaPartOfTermdLambda)
        deallocate(localdZetaPartOfTermdLambda)

      ! *********************************************************
      ! Add the sensitivity of the ExB d/dtheta term:
      ! *********************************************************

      factor = alpha*Delta/two*dPhiHatdpsiHat
      allocate(dThetaPartOfTerm(Ntheta,Ntheta))
      allocate(localdThetaPartOfTerm(Ntheta,localNtheta))
      allocate(rowIndices(localNtheta))
      allocate(colIndices(Ntheta))

      do L=0,(Nxi-1)

         if (ExBDerivativeSchemeTheta==0) then
            ddthetaToUse = ddtheta
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
           do itheta=1,Ntheta
            !thetaPartOfTerm(itheta,:) = ddthetaToUse(itheta,:) / (BHat(itheta,izeta) ** 2) &
              !* DHat(itheta,izeta) * BHat_sub_zeta(itheta,izeta)
            select case(whichLambda)
              case (-1) ! Er
                if (masterProc) then
                  print *,"Error! Er sensitivity not yet implemented."
                end if
                stop
              case (0) ! BHat
                dThetaPartOfTermdLambda(itheta,:) = - 2 * ddthetaToUse(itheta,:) &
                  * dBHatdFourier(itheta,izeta,whichMode) / (BHat(itheta,izeta) ** 3) &
                  * DHat(itheta,izeta) * BHat_sub_zeta(itheta,izeta)
              case (1) ! BHat_sup_theta
                dThetaPartOfTermdLambda(itheta,:) = 0
              case (2) ! BHat_sup_zeta
                dThetaPartOfTermdLambda(itheta,:) = 0
              case (3) ! BHat_sub_theta
                dThetaPartOfTermdLambda(itheta,:) = 0
              case (4) ! BHat_sub_zeta
                dThetaPartOfTermdLambda(itheta,:) = ddthetaToUse(itheta,:) / (BHat(itheta,izeta) ** 2) &
                  * DHat(itheta,izeta) * dBHat_sub_zetadFourier(itheta,izeta,whichMode)
              case (5) ! DHat
                dThetaPartOfTermdLambda(itheta,:) = ddthetaToUse(itheta,:) / (BHat(itheta,izeta) ** 2) &
                  * dDHatdFourier(itheta,izeta,whichMode) * BHat_sub_zeta(itheta,izeta)
            end select

           end do

            ! PETSc uses the opposite convention to Fortran:
            dThetaPartOfTermdLambda = transpose(dThetaPartOfTermdLambda)
            localdThetaPartOfTermdLambda = dThetaPartOfTermdLambda(:,ithetaMin:ithetaMax)
            
            !do ix=ixMin,Nx
            do ix=max(ixMin,min_x_for_L(L)),Nx
               do itheta=1,localNtheta
                  rowIndices(itheta)=getIndex(ispecies,ix,L+1,itheta+ithetaMin-1,izeta,BLOCK_F)
               end do
               do itheta=1,Ntheta
                  colIndices(itheta)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
               end do
               
               call MatSetValuesSparse(dMatrixdLambda, localNtheta, rowIndices, Ntheta, colIndices, &
                    localdThetaPartOfTermdLambda, ADD_VALUES, ierr)
            end do
         end do
      end do
      deallocate(rowIndices)
      deallocate(colIndices)
      deallocate(dThetaPartOfTermdLambda)
      deallocate(localdThetaPartOfTermdLambda)

     ! *********************************************************
     ! Add the sensitivity of the ExB d/dzeta term:
     ! *********************************************************

      factor = -alpha*Delta/two*dPhiHatdpsiHat
      allocate(dZetaPartOfTermdLambda(Nzeta,Nzeta))
      allocate(localdZetaPartOfTermdLambda(Nzeta,localNzeta))
      allocate(rowIndices(localNzeta))
      allocate(colIndices(Nzeta))
      do L=0,(Nxi-1)

         if (ExBDerivativeSchemeZeta==0) then
            ddzetaToUse = ddzeta
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
           do izeta=1,Nzeta
            !  zetaPartOfTerm(izeta,:) = ddzetaToUse(izeta,:) / (BHat(itheta,izeta) ** 2) &
            !       * DHat(itheta,izeta) * BHat_sub_theta(itheta,izeta)
            select case(whichLambda)
              case (-1) ! Er
                if (masterProc) then
                  print *,"Error! Er sensitivity not yet implemented."
                end if
                stop
              case (0) ! BHat
                dZetaPartOfTermdLambda(izeta,:) = -2 * ddzetaToUse(izeta,:) &
                * dBHatdFourier(itheta,izeta,whichMode) / (BHat(itheta,izeta) ** 3) &
                * DHat(itheta,izeta) * BHat_sub_theta(itheta,izeta)
              case (1) ! BHat_sup_theta
                dZetaPartOfTermdLambda(izeta,:) = 0
              case (2) ! BHat_sup_zeta
                dZetaPartOfTermdLambda(izeta,:) = 0
              case (3) ! BHat_sub_theta
                dZetaPartOfTermdLambda(izeta,:) = ddzetaToUse(izeta,:) / (BHat(itheta,izeta) ** 2) &
                  * DHat(itheta,izeta) * dBHat_sub_thetadFourier(itheta,izeta,whichMode)
              case (4) ! BHat_sub_zeta
                dZetaPartOfTermdLambda(izeta,:) = 0
              case (5) ! DHat
                dZetaPartOfTermdLambda(izeta,:) = ddzetaToUse(izeta,:) / (BHat(itheta,izeta) ** 2) &
                  * dDHatdFourier(itheta,izeta,whichMode) * BHat_sub_theta(itheta,izeta)
            end select

           end do

            ! PETSc uses the opposite convention to Fortran:
            dZetaPartOfTermdLambda = transpose(dZetaPartOfTermdLambda)
            localdZetaPartOfTermdLambda = dZetaPartOfTermdLambda(:,izetaMin:izetaMax)
            
            !do ix=ixMin,Nx
            do ix=max(ixMin,min_x_for_L(L)),Nx
               do izeta=1,localNzeta
                  rowIndices(izeta)=getIndex(ispecies,ix,L+1,itheta,izeta+izetaMin-1,BLOCK_F)
               end do
               do izeta=1,Nzeta
                  colIndices(izeta)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
               end do
               
               call MatSetValuesSparse(dMatrixdLambda, localNzeta, rowIndices, Nzeta, colIndices, &
                    localdZetaPartOfTermdLambda, ADD_VALUES, ierr)
            end do
         end do
      end do
      deallocate(rowIndices)
      deallocate(colIndices)
      deallocate(dZetaPartOfTermdLambda)
      deallocate(localdZetaPartOfTermdLambda)

       ! *********************************************************
       ! Add the sensitivity of the standard mirror term:
       ! *********************************************************

      do itheta=ithetaMin,ithetaMax
         do izeta=izetaMin,izetaMax
            !factor = -sqrtTHat/(2*sqrtMHat*BHat(itheta,izeta)*BHat(itheta,izeta)) &
            !     * (BHat_sup_theta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
            !     + BHat_sup_zeta(itheta,izeta) * dBHatdzeta(itheta,izeta))
            select case(whichLambda)
              case (-1) ! Er
                if (masterProc) then
                  print *,"Error! Er sensitivity not yet implemented."
                end if
                stop
              case (0) ! BHat
                dFactordLambda = &
                  ! Term from 1/(BHat**2)
                  sqrtTHat/(sqrtMHat*(BHat(itheta,izeta)**3)) &
                  * dBHatdFourier(itheta,izeta,whichMode) &
                  * (BHat_sup_theta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                  + BHat_sup_zeta(itheta,izeta) * dBHatdzeta(itheta,izeta)) &
                  ! Term from dBHatdtheta
                  -sqrtTHat/(2*sqrtMHat*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                  * BHat_sup_theta(itheta,izeta)*dBhatdthetadFourier(itheta,izeta,whichMode) &
                  ! Term from dBHatdzeta 
                  -sqrtTHat/(2*sqrtMHat*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                  * BHat_sup_zeta(itheta,izeta)*dBhatdzetadFourier(itheta,izeta,whichMode) &
              case (1) ! BHat_sup_theta
                dFactordLambda = -sqrtTHat/(2*sqrtMHat*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                  * (dBHat_sup_thetadFourier(itheta,izeta,whichMode)*dBHatdtheta(itheta,izeta))
              case (2) ! BHat_sup_zeta
                dFactordLambda = -sqrtTHat/(2*sqrtMHat*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                  * (dBHat_sup_zetadFourier(itheta,izeta,whichMode)*dBHatdzeta(itheta,izeta))
              case (3) ! BHat_sub_theta
                dFactordLambda = 0
              case (4) ! BHat_sub_zeta
                dFactordLambda = 0
              case (5) ! DHat
                dFactordLambda = 0
            end select
            
            do ix=ixMin,Nx
               do L=0,(Nxi_for_x(ix)-1)
                  rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                  
                  if (L<Nxi_for_x(ix)-1) then
                     ! Super-diagonal-in-L term:
                     ell = L+1
                     colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                     call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                          (L+1)*(L+2)/(2*L+three)*x(ix)*dFactordLambda, ADD_VALUES, ierr)
                  end if
                  
                  if (L>0) then
                     ! Sub-diagonal-in-L term:
                     ell = L-1
                     colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                     call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                          -L*(L-1)/(2*L-one)*x(ix)*dFactordLambda, ADD_VALUES, ierr)
                  end if
               end do
            end do
         end do
      end do

       ! *********************************************************
       ! Add the sensitivity of the non-standard d/dxi term associated with E_r:
       ! *********************************************************

       if (includeElectricFieldTermInXiDot) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax

               temp = BHat_sub_zeta(itheta,izeta) * dBHatdtheta(itheta,izeta) &
                     - BHat_sub_theta(itheta,izeta) * dBHatdzeta(itheta,izeta)

               !factor = alpha*Delta*dPhiHatdpsiHat/(4*(BHat(itheta,izeta)**3)) &
               !   * DHat(itheta,izeta) * temp

              select case(whichLambda)
                case (-1) ! Er
                  if (masterProc) then
                    print *,"Error! Er sensitivity not yet implemented."
                  end if
                  stop
                case (0) ! BHat
                  dTempdLambda = BHat_sub_zeta(itheta,izeta) * dBHatdthetadFourier(itheta,izeta,whichMode) &
                    - BHat_sub_theta(itheta,izeta) * dBHatdzetadFourier(itheta,izeta,whichMode)
                  dFactordLambda = &
                    ! Term from 1/(BHat**3)
                    -3*alpha*Delta*dPhiHatdpsiHat/(4*(BHat(itheta,izeta)**4)) &
                    * DHat(itheta,izeta) * temp * dBHatdFourier(itheta,izeta,whichMode) &
                    ! Term from dBHatdtheta and dBHatdzeta
                    + alpha*Delta*dPhiHatdpsiHat/(4*(BHat(itheta,izeta)**3)) &
                    * DHat(itheta,izeta) * dTempdLambda
                case (1) ! BHat_sup_theta
                  dFactordLambda = 0
                case (2) ! BHat_sup_zeta
                  dFactordLambda = 0
                case (3) ! BHat_sub_theta
                  dTempdLambda = - dBHat_sub_thetadLambda(itheta,izeta,whichMode) * dBHatdzeta(itheta,izeta)
                  dFactordLambda = alpha*Delta*dPhiHatdpsiHat/(4*(BHat(itheta,izeta)**3)) &
                    * DHat(itheta,izeta) * dTempdLambda
                case (4) ! BHat_sub_zeta
                  dTempdLambda = dBHat_sub_zetadLambda(itheta,izeta,whichMode) * dBHatdtheta(itheta,izeta)
                  dFactordLambda = alpha*Delta*dPhiHatdpsiHat/(4*(BHat(itheta,izeta)**3)) &
                    * DHat(itheta,izeta) * dTempdLambda
                case (5) ! DHat
                  dFactordLambda = alpha*Delta*dPhiHatdpsiHat/(4*(BHat(itheta,izeta)**3)) &
                  * dDHatdFourier(itheta,izeta,whichMode) * temp
              end select

                do ix=ixMin,Nx
                   do L=0,(Nxi_for_x(ix)-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)

                      ! Diagonal-in-L term
                      call MatSetValueSparse(dMatrixdLambda, rowIndex, rowIndex, &
                           (L+1)*L/((2*L-one)*(2*L+three))*dFactordLambda, ADD_VALUES, ierr)

                       if (L<Nxi_for_x(ix)-2) then
                          ! Super-super-diagonal-in-L term:
                          ell = L+2
                          colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                          call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                               (L+3)*(L+2)*(L+1)/((two*L+5)*(2*L+three))*dFactordLambda, ADD_VALUES, ierr)
                       end if

                       if (L>1) then
                          ! Sub-sub-diagonal-in-L term:
                          ell = L-2
                          colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                          call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                               -L*(L-1)*(L-2)/((2*L-3)*(2*L-one))*dFactordLambda, ADD_VALUES, ierr)
                       end if
                   end do
                end do
             end do
          end do
       end if

    ! *********************************************************
    ! Add the sensitivity of the collisionless d/dx term associated with E_r
    ! *********************************************************

    allocate(xPartOfXDot(Nx,Nx))
    allocate(xPartOfXDot_plus(Nx,Nx))
    allocate(xPartOfXDot_minus(Nx,Nx))
    allocate(rowIndices(Nx))
    allocate(colIndices(Nx))
    factor = -alpha*Delta*dPhiHatdpsiHat/4

    do L=0,(Nxi-1)
       if (L>0 .and. pointAtX0) then
          ixMinCol = 2
       else
          ixMinCol = 1
       end if

       ! Upwind in x
       ddxToUse_plus = ddx_xDot_plus
       ddxToUse_minus = ddx_xDot_minus

       do ix=1,Nx
          xPartOfXDot_plus(ix,:) = x(ix) * ddxToUse_plus(ix,:)
          xPartOfXDot_minus(ix,:) = x(ix) * ddxToUse_minus(ix,:)
       end do

       do itheta=ithetaMin,ithetaMax
          do izeta=izetaMin,izetaMax

             xDotFactor = factor*DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                  * (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) &
                  - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))

             ! Should use same upwinding direction as matrix
             if (xDotFactor>0) then
                xPartOfXDot = xPartOfXDot_plus
             else
                xPartOfXDot = xPartOfXDot_minus
             end if

              select case(whichLambda)
                case (-1) ! Er
                  if (masterProc) then
                    print *,"Error! Er sensitivity not yet implemented."
                  end if
                  stop
                case (0) ! BHat
                  dXDotFactordLambda = &
                    ! Term from 1/(BHat**3)
                    -3*factor*DHat(itheta,izeta)*dBHatdFourier(itheta,izeta,whichMode) &
                    /(BHat(itheta,izeta)**4) &
                    * (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) &
                    - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))
                    ! Term from dBHatdzeta
                    + factor*DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                    * BHat_sub_theta(itheta,izeta)*dBHatdzetadFourier(itheta,izeta,whichMode)
                    ! Term from dBHatdtheta
                    - factor*Dhat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                    * Bhat_sub_zeta(itheta,izeta)*dBHatdthetadFourier(itheta,izeta,whichMode)
                case (1) ! BHat_sup_theta
                  dXDotFactordLambda = 0
                case (2) ! BHat_sup_zeta
                  dXDotFactordLambda = 0
                case (3) ! BHat_sub_theta
                  dXDotFactordLambda = factor*DHat/(BHat(itheta,izeta)**3) &
                    * dBHat_sub_thetadFourier(itheta,izeta,whichMode)*dBHatdzeta(itheta,izeta)
                case (4) ! BHat_sub_zeta
                  dXDotFactordLambda = -factor*DHat/(BHat(itheta,izeta)**3) &
                    * dBHat_sub_zetadFourier(itheta,izeta,whichMode)*dBHatdtheta(itheta,izeta)
                case (5) ! DHat
                  dXDotFactordLambda = factor*dDHatdFourier(itheta,izeta,whichMode)/(BHat(itheta,izeta)**3) &
                    *(BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) &
                    - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))
              end select

             rowIndices = -1
             do ix=min_x_for_L(L),Nx
                rowIndices(ix)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
             end do

             ! Term that is diagonal in L:
             colIndices = rowIndices
             stuffToAdd = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*dXDotFactordLambda
             do ix_col=max(ixMinCol,min_x_for_L(L)),Nx
                do ix_row=max(ixMin,min_x_for_L(L)),Nx
                   call MatSetValueSparse(dMatrixdLambda, rowIndices(ix_row), colIndices(ix_col), &
                        stuffToAdd*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                end do
             end do


              ! Term that is super-super-diagonal in L:
              if (L<(Nxi-2)) then
                 ell = L + 2
                 stuffToAdd = (L+1)*(L+2)/((two*L+5)*(2*L+3))*(dXDotFactordLambda)
                 !do ix_col=ixMinCol,Nx
                 do ix_col=max(ixMinCol,min_x_for_L(ell)),Nx
                    colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
                    !do ix_row=ixMin,Nx
                    do ix_row=max(ixMin,min_x_for_L(L)),Nx
                       call MatSetValueSparse(matrix, rowIndices(ix_row), colIndex, &
                            stuffToAdd*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                    end do
                 end do
              end if

              ! Term that is sub-sub-diagonal in L:
              if (L>1) then
                 ell = L - 2
                 stuffToAdd = L*(L-1)/((two*L-3)*(2*L-1))*(dXDotFactordLambda)
                 do ix_col=max(ixMinCol,min_x_for_L(ell)),Nx
                    colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
                    do ix_row=max(ixMin,min_x_for_L(L)),Nx
                       call MatSetValueSparse(matrix, rowIndices(ix_row), colIndex, &
                            stuffToAdd*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                    end do
                 end do
              end if

          end do
       end do
    end do
    deallocate(rowIndices)
    deallocate(colIndices)
    deallocate(xPartOfXDot)
    deallocate(xPartOfXDot_plus)
    deallocate(xPartOfXDot_minus)

    ! *********************************************************************************
    ! I do not need the collision operator, boundary conditions for f at x=0
    ! or sources as they are independent of geometry
    ! *********************************************************************************

    ! *******************************************************************************
    ! Add the sensitivity of the density and pressure constraints:
    ! *******************************************************************************

    if (procThatHandlesConstraints) then
      L=0
      do itheta=1,Ntheta
         do izeta=1,Nzeta
            !factor = thetaWeights(itheta)*zetaWeights(izeta)/DHat(itheta,izeta)
            dFactordLambda = 0
            if (whichLambda == 5) then
              dFactordLambda = -thetaWeights(itheta)*zetaWeights(izeta) &
                *dDHatdFourier(itheta,izeta,whichMode)/(DHat(itheta,izeta)**2)
            end if

            do ix=1,Nx
               do ispecies=1,Nspecies
                  colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)

                  rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
                  call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                       x2(ix)*xWeights(ix)*dFactordLambda, ADD_VALUES, ierr)

                  rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
                  call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                       x2(ix)*x2(ix)*xWeights(ix)*dFactordLambda, ADD_VALUES, ierr)
               end do
            end do

         end do
      end do
    endif

    ! *******************************************************************************
    ! Done inserting values into the matrices.
    ! Now finalize the matrix
    ! *******************************************************************************

    call PetscTime(time2, ierr)
    if (masterProc) then
       print *,"Time to pre-assemble sensitivity of matrix: ", time2-time1, " seconds."
    end if
    call PetscTime(time1, ierr)

    call MatAssemblyBegin(dMatrixdLambda, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(dMatrixdLambda, MAT_FINAL_ASSEMBLY, ierr)

    call PetscTime(time2, ierr)
    if (masterProc) then
       print *,"Time to assemble sensitivity of matrix: ", time2-time1, " seconds."
    end if
    call PetscTime(time1, ierr)


    call MatGetInfo(matrix, MAT_GLOBAL_SUM, myMatInfo, ierr)
    NNZ = nint(myMatInfo(MAT_INFO_NZ_USED))
    NNZAllocated = nint(myMatInfo(MAT_INFO_NZ_ALLOCATED))
    NMallocs = nint(myMatInfo(MAT_INFO_MALLOCS))
    if (masterProc) then
       print *,"# of nonzeros in sensitivity of matrix:",NNZ, ", allocated:",NNZAllocated, &
            ", mallocs:",NMallocs," (should be 0)"
    end if

  end subroutine populatedMatrixdLambda
