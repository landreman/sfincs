! See populateMatrix.f90 (where matrix is constructed). Much of the code has been copied.

  !> Populates \f$\partial \mathbb{L}/\partial \lambda\f$ to compute explicit dependence of integrated quantiteis on geometry.
  !! @param dMatrixdLambda Matrix to be populated. Should be allocated by calling subroutine.
  !! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
  !! @param whichMode Indicates index of ms and ns for derivative.
  subroutine populatedMatrixdLambda(dMatrixdLambda, whichLambda, whichMode)

#include "PETScVersions.F90"

    use petscmat
    use globalVariables
    use sparsify
    use indices
    use xGrid, only: xGrid_k

    implicit none

    Mat :: dMatrixdLambda
    integer :: whichLambda, whichMode
    integer :: i, ixMin, ixMax, L, ispecies, izeta, itheta, ix, ixMinCol
    PetscScalar, dimension(:,:), allocatable :: thetaPartOfTerm, localThetaPartOfTerm, &
      zetaPartOfTerm, localZetaPartOfTerm
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor
    PetscScalar, dimension(:,:), allocatable :: ddthetaToUse, ddzetaToUse, ddxToUse, ddxToUse_plus, &
      ddxToUse_minus
    PetscScalar :: temp, xDotFactor, stuffToAdd, dTempdLambda
    PetscScalar, dimension(:,:), allocatable :: xPartOfXDot, xPartOfXDot_plus, xPartOfXDot_minus
    integer :: ix_row, ix_col, rowIndex, colIndex, ell
    PetscLogDouble :: time1, time2
    integer :: NNZ, NNZAllocated, NMallocs
    PetscScalar :: dBHatdLambda, dBHat_sub_thetadLambda, dBHat_sub_zetadLambda, dBHat_sup_thetadLambda
    PetscScalar :: dBHat_sup_zetadLambda, dBHatdthetadLambda, dBHatdzetadLambda, dDHatdLambda
    PetscScalar :: geometricFactor, dFSABHat2dBmn, dFSABHat2dDmn, dVPrimeHatdDmn
    integer :: m,n
    PetscScalar :: angle, cos_angle, sin_angle, dVPrimeHatdLambda, dFSABHat2dLambda, dVPrimeHatdBmn

    if (whichLambda > 0) then
      m = ms_sensitivity(whichMode)
      n = ns_sensitivity(whichMode)
    end if

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

    if (useDKESExBDrift) then
      dVPrimeHatdBmn  = zero
      do itheta=1,Ntheta
        do izeta=1,Nzeta
          angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
          cos_angle = cos(angle)
          dVPrimeHatdBmn = dVPrimeHatdBmn - (two*(GHat+iota*IHat)) * thetaWeights(itheta) * zetaWeights(izeta) &
            * BHat(itheta,izeta)**(-3) * cos_angle
        end do
      end do
      dFSABHat2dBmn = -4*pi*pi*(GHat+iota*IHat)*dVPrimeHatdBmn/(VPrimeHat**2)
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

      allocate(thetaPartOfTerm(Ntheta,Ntheta))
      allocate(localThetaPartOfTerm(Ntheta,localNtheta))
      allocate(rowIndices(localNtheta))
      allocate(colIndices(Ntheta))

      ! *********************************************************
      ! Add the sensitivity of the streaming d/dtheta term:
      ! *********************************************************

      do L=0,(Nxi-1)
        ddthetaToUse = ddtheta

        do izeta=izetaMin,izetaMax
          do itheta=1,Ntheta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            sin_angle = sin(angle)

            ! geometricFactor = BHat_sup_theta/BHat = iota*BHat/(GHat+iota*IHat)
            select case(whichLambda)
              case (0) ! Er
                geometricFactor = zero
              case (1) ! BHat
                dBHatdLambda = cos_angle
                geometricFactor = dBHatdLambda*iota/(GHat+iota*IHat)
              case (2) ! IHat
                geometricFactor = -iota*iota*BHat(itheta,izeta)/(GHat+iota*IHat)**2
              case (3) ! GHat
                geometricFactor = -iota*BHat(itheta,izeta)/(GHat+iota*IHat)**2
              case (4) ! iota
                geometricFactor = BHat(itheta,izeta)/(GHat+iota*IHat) &
                  -iota*IHat*BHat(itheta,izeta)/(GHat+iota*IHat)**2
            end select
            thetaPartOfTerm(itheta,:) = sqrtTHat/sqrtMHat * ddthetaToUse(itheta,:) * geometricFactor
          end do ! itheta

            ! PETSc uses the opposite convention to Fortran:
            thetaPartOfTerm = transpose(thetaPartOfTerm)
            localThetaPartOfTerm = thetaPartOfTerm(:,ithetaMin:ithetaMax)

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
                  (L+1)/(two*L+three)*x(ix)*localThetaPartOfTerm, ADD_VALUES, ierr)
              end if

              ! Sub-diagonal-in-L term
              if (L > 0) then
                ell = L-1
                do itheta=1,Ntheta
                  colIndices(itheta) = getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                end do

                call MatSetValuesSparse(dMatrixdLambda, localNtheta, rowIndices, Ntheta, colIndices, &
                  L/(two*L-one)*x(ix)*localThetaPartOfTerm, ADD_VALUES, ierr)
              end if

            end do ! ix
          end do ! izeta
        end do ! L

      deallocate(rowIndices)
      deallocate(colIndices)
      deallocate(thetaPartOfTerm)
      deallocate(localThetaPartOfTerm)

      ! *********************************************************
      ! Add the sensitivity of the streaming d/dzeta term:
      ! *********************************************************

      allocate(zetaPartOfTerm(Nzeta,Nzeta))
      allocate(localZetaPartOfTerm(Nzeta,localNzeta))
      allocate(rowIndices(localNzeta))
      allocate(colIndices(Nzeta))

!      if (.false.) then

      do L=0,(Nxi-1)
        ddzetaToUse = ddzeta
        do itheta=ithetaMin, ithetaMax

          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            sin_angle = sin(angle)

            ! geometricFactor = BHat_sup_zeta/BHat = BHat/(GHat+iota*IHat)
            select case(whichLambda)
              case (0) ! Er
                geometricFactor = zero
              case (1) ! BHat
                dBHatdLambda = cos_angle
                geometricFactor = dBHatdLambda/(GHat+iota*IHat)
              case (2) ! IHat
                geometricFactor = -BHat(itheta,izeta)*iota/(GHat+iota*IHat)**2
              case (3) ! GHat
                geometricFactor = -BHat(itheta,izeta)/(GHat+iota*IHat)**2
              case (4) ! iota
                geometricFactor = -IHat*BHat(itheta,izeta)/(GHat+iota*IHat)
            end select
            zetaPartOfTerm(izeta,:) = (sqrtTHat/sqrtMHat)*ddzetaToUse(izeta,:)*geometricFactor
          end do

          ! PETSc uses the opposite convention to Fortran:
          zetaPartOfTerm = transpose(zetaPartOfTerm)
          localZetaPartOfTerm = zetaPartOfTerm(:,izetaMin:izetaMax)

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
                     (L+1)/(two*L+three)*x(ix)*localZetaPartOfTerm, ADD_VALUES, ierr)
             end if
             
             ! Sub-diagonal-in-L term
             if (L > 0) then
                ell = L - 1
                do izeta = 1,Nzeta
                   colIndices(izeta)=getIndex(ispecies, ix, ell+1, itheta, izeta, BLOCK_F)
                end do
                
                call MatSetValuesSparse(dMatrixdLambda, localNzeta, rowIndices, Nzeta, colIndices, &
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
      ! Add the sensitivity of the ExB d/dtheta term:
      ! *********************************************************


      allocate(thetaPartOfTerm(Ntheta,Ntheta))
      allocate(localThetaPartOfTerm(Ntheta,localNtheta))
      allocate(rowIndices(localNtheta))
      allocate(colIndices(Ntheta))

     if (whichLambda==1) then ! BHat
        dVPrimeHatdLambda = zero
        do itheta=1,Ntheta
          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)

            dVPrimeHatdLambda = dVPrimeHatdLambda - two*thetaWeights(itheta)*zetaWeights(izeta)*cos_angle/(DHat(itheta,izeta)*BHat(itheta,izeta))
          end do
        end do
      end if

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
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            sin_angle = sin(angle)

            ! factor = alpha*Delta/two*dPhiHatdpsiHat
            ! geometricFactor = DHat(itheta,izeta)*BHat_sub_zeta(itheta,izeta)/(BHat*BHat) = GHat/(GHat + iota*IHat)
            ! useDKESExBDrift = .true.
              ! geometricFactor = DHat(itheta,izeta)*BHat_sub_zeta(itheta,izeta)/(FSABHat2) = (GHat/(GHat + iota*IHat))*(BHat*BHat/FSABHat2)
            select case(whichLambda)
              case (0) ! Er
                if (useDKESExBDrift) then
                  geometricFactor = DHat(itheta,izeta)*BHat_sub_zeta(itheta,izeta)/(FSABHat2)
                else
                  geometricFactor = DHat(itheta,izeta)*BHat_sub_zeta(itheta,izeta)/(BHat(itheta,izeta)*BHat(itheta,izeta))
                end if
                factor = alpha*Delta/two
              case (1) ! BHat
                if (useDKESExBDrift) then
                  dBHatdLambda = cos_angle
                  geometricFactor = (GHat/(GHat + iota*IHat))*(two*BHat(itheta,izeta)*dBHatdLambda/FSABHat2 &
                    - BHat(itheta,izeta)*BHat(itheta,izeta)*dFSABHat2dBmn/(FSABHat2**2))
                else
                  geometricFactor = zero
                end if
                factor = alpha*Delta/two*dPhiHatdpsiHat
              case (2) ! IHat
                if (useDKESExBDrift) then
                  geometricFactor = -(GHat*iota/(GHat+iota*IHat)**2)*(BHat(itheta,izeta)*BHat(itheta,izeta)/FSABHat2)
                else
                  geometricFactor = -GHat*iota/(GHat+iota*IHat)**2
                end if
                factor = alpha*Delta/two*dPhiHatdpsiHat
              case (3) ! GHat
                if (useDKESExBDrift) then
                  geometricFactor = (one/(GHat+iota*IHat)-GHat/(GHat+iota*IHat)**2) &
                    * (BHat(itheta,izeta)*BHat(itheta,izeta)/FSABHat2)
                else
                  geometricFactor = (one/(GHat+iota*IHat)-GHat/(GHat+iota*IHat)**2)
                end if
                factor = alpha*Delta/two*dPhiHatdpsiHat
              case (4) ! iota
                if (useDKESExBDrift) then
                  geometricFactor = (-GHat*IHat/(GHat+iota*IHat)**2)*(BHat(itheta,izeta)*BHat(itheta,izeta)/FSABHat2)
                else
                  geometricFactor = (-GHat*IHat/(GHat+iota*IHat)**2)
                end if
                factor = alpha*Delta/two*dPhiHatdpsiHat
            end select
              thetaPartOfTerm(itheta,:) = ddthetaToUse(itheta,:)*geometricFactor*factor

           end do ! itheta

            ! PETSc uses the opposite convention to Fortran:
            thetaPartOfTerm = transpose(thetaPartOfTerm)
            localThetaPartOfTerm = thetaPartOfTerm(:,ithetaMin:ithetaMax)
            
            !do ix=ixMin,Nx
            do ix=max(ixMin,min_x_for_L(L)),Nx
               do itheta=1,localNtheta
                  rowIndices(itheta)=getIndex(ispecies,ix,L+1,itheta+ithetaMin-1,izeta,BLOCK_F)
               end do
               do itheta=1,Ntheta
                  colIndices(itheta)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
               end do
               
               call MatSetValuesSparse(dMatrixdLambda, localNtheta, rowIndices, Ntheta, colIndices, &
                    localThetaPartOfTerm, ADD_VALUES, ierr)
            end do
         end do
      end do
      deallocate(rowIndices)
      deallocate(colIndices)
      deallocate(thetaPartOfTerm)
      deallocate(localThetaPartOfTerm)

     ! *********************************************************
     ! Add the sensitivity of the ExB d/dzeta term:
     ! *********************************************************

      allocate(zetaPartOfTerm(Nzeta,Nzeta))
      allocate(localZetaPartOfTerm(Nzeta,localNzeta))
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
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            sin_angle = sin(angle)

            ! geometryFactor = DHat(itheta,izeta)*BHat_sub_theta(itheta,izeta)/(BHat(itheta,izeta)*BHat(itheta,izeta)) = IHat/(GHat+iota*IHat)
            ! geometryFactor = DHat(itheta,izeta)*BHat_sub_theta(itheta,izeta)/FSABHat2 = (IHat/(GHat+iota*IHat))*(BHat*BHat/FSABHat2)
            select case(whichLambda)
              case (0) ! Er
                if (useDKESExBDrift) then
                  geometricFactor = DHat(itheta,izeta)*BHat_sub_theta(itheta,izeta)/FSABHat2
                else
                  geometricFactor = DHat(itheta,izeta)*BHat_sub_theta(itheta,izeta)/(BHat(itheta,izeta)*BHat(itheta,izeta))
                end if
                factor = -alpha*Delta/two
              case (1) ! BHat
                dBHatdLambda = cos_angle
                if (useDKESExBDrift) then
                  geometricFactor = (2*BHat(itheta,izeta)*dBHatdLambda/FSABHat2)*(IHat/(GHat+iota*IHat))
                else
                  geometricFactor = 0
                end if
                factor = -alpha*Delta/two*dPhiHatdpsiHat
              case (2) ! BHat_sub_theta / IHat
                if (useDKESExBDrift) then
                  geometricFactor = (one/(GHat+iota*IHat)-(IHat*iota/(GHat+iota*IHat)**2)) &
                    *(BHat(itheta,izeta)*BHat(itheta,izeta))/FSABHat2
                else
                  geometricFactor = (one/(GHat+iota*IHat)-(IHat*iota/(GHat+iota*IHat)**2))
                end if
                factor = -alpha*Delta/two*dPhiHatdpsiHat
              case (3) ! BHat_sub_zeta / GHat
                if (useDKESExBDrift) then
                  geometricFactor = -(IHat/(GHat+iota*IHat)**2)*(BHat(itheta,izeta)*BHat(itheta,izeta))/FSABHat2
                else
                  geometricFactor = -(IHat/(GHat+iota*IHat)**2)
                end if
                factor = -alpha*Delta/two*dPhiHatdpsiHat
              case (4) ! BHat_sup_theta / iota
                if (useDKESExBDrift) then
                  geometricFactor = -(IHat*IHat/(GHat+iota*IHat)**2)*(BHat(itheta,izeta)*BHat(itheta,izeta)/FSABHat2)
                else
                  geometricFactor = -(IHat*IHat/(GHat+iota*IHat)**2)
                end if
                factor = -alpha*Delta/two*dPhiHatdpsiHat
            end select
            zetaPartOfTerm(izeta,:) = ddzetaToUse(izeta,:)*geometricFactor*factor
           end do

            ! PETSc uses the opposite convention to Fortran:
            zetaPartOfTerm = transpose(zetaPartOfTerm)
            localZetaPartOfTerm = zetaPartOfTerm(:,izetaMin:izetaMax)
            
            !do ix=ixMin,Nx
            do ix=max(ixMin,min_x_for_L(L)),Nx
               do izeta=1,localNzeta
                  rowIndices(izeta)=getIndex(ispecies,ix,L+1,itheta,izeta+izetaMin-1,BLOCK_F)
               end do
               do izeta=1,Nzeta
                  colIndices(izeta)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
               end do
               
               call MatSetValuesSparse(dMatrixdLambda, localNzeta, rowIndices, Nzeta, colIndices, &
                    localZetaPartOfTerm, ADD_VALUES, ierr)
            end do
         end do
      end do
      deallocate(rowIndices)
      deallocate(colIndices)
      deallocate(zetaPartOfTerm)
      deallocate(localZetaPartOfTerm)

       ! *********************************************************
       ! Add the sensitivity of the standard mirror term:
       ! *********************************************************
       ! (1/(BHat*BHat))*(BHat_sup_theta*dBHatdtheta+BHat_sup_zeta*dBHatdzeta)
          ! = (1/(GHat+iota*IHat))*(iota*dBHatdtheta + dBHatdzeta)

      do itheta=ithetaMin,ithetaMax
         do izeta=izetaMin,izetaMax
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            sin_angle = sin(angle)

            select case(whichLambda)
              case (0) ! Er
                geometricFactor = zero
              case (1) ! BHat
                dBHatdLambda = cos_angle
                dBHatdthetadLambda = -m*sin_angle
                dBHatdzetadLambda = n*Nperiods*sin_angle
          geometricFactor = (1/(GHat+iota*IHat))*(iota*dBHatdthetadLambda + dBHatdzetadLambda)
              case (2) ! BHat_sub_theta / IHat
          geometricFactor = -(iota/(GHat+iota*IHat)**2)*(iota*dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta))
              case (3) ! BHat_sub_zeta / GHat
        geometricFactor = -(one/(GHat+iota*IHat)**2)*(iota*dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta))
              case (4) ! BHat_sup_theta / iota
        geometricFactor = -IHat/(GHat+iota*IHat)**2*(iota*dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta)) &
            + (dBHatdtheta(itheta,izeta)/(GHat+iota*IHat))
            end select
            factor = (-sqrtTHat/(two*sqrtMHat))*geometricFactor

            do ix=ixMin,Nx
               do L=0,(Nxi_for_x(ix)-1)
                  rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                  
                  if (L<Nxi_for_x(ix)-1) then
                     ! Super-diagonal-in-L term:
                     ell = L+1
                     colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                     call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                          (L+1)*(L+2)/(2*L+three)*x(ix)*factor, ADD_VALUES, ierr)
                  end if
                  
                  if (L>0) then
                     ! Sub-diagonal-in-L term:
                     ell = L-1
                     colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                     call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                          -L*(L-1)/(2*L-one)*x(ix)*factor, ADD_VALUES, ierr)
                  end if
               end do
            end do
         end do
      end do

       ! *********************************************************
       ! Add the sensitivity of the non-standard d/dxi term associated with E_r:
       ! *********************************************************

       ! (BHat_sub_zeta(itheta,izeta) * dBHatdtheta(itheta,izeta) &
       !   - BHat_sub_theta(itheta,izeta) * dBHatdzeta(itheta,izeta))*DHat/(BHat(itheta,izeta)**3)
       ! = (GHat*dBHatdtheta-IHat*dBHatdzeta)*(1/(BHat*(GHat+iota*IHat)))

       if (includeElectricFieldTermInXiDot) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax
                angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                cos_angle = cos(angle)
                sin_angle = sin(angle)

                temp = BHat_sub_zeta(itheta,izeta) * dBHatdtheta(itheta,izeta) &
                     - BHat_sub_theta(itheta,izeta) * dBHatdzeta(itheta,izeta)

              select case(whichLambda)
                case (0) ! Er
                  geometricFactor = DHat(itheta,izeta)*temp/(BHat(itheta,izeta)**3)
                  factor = (alpha*Delta/four)*geometricFactor
                case (1) ! BHat
                  dBHatdthetadLambda = -m*sin_angle
                  dBHatdzetadLambda = n*Nperiods*sin_angle
                  dBHatdLambda = cos_angle
          dTempdLambda = BHat_sub_zeta(itheta,izeta) * dBHatdthetadLambda &
            - BHat_sub_theta(itheta,izeta) * dBHatdzetadLambda
!                  factor = (alpha*Delta*dPhiHatdpsiHat/four)*geometricFactor
          geometricFactor = dTempdLambda/ (BHat(itheta,izeta)*(GHat+iota*IHat)) &
            - temp/(BHat(itheta,izeta)*BHat(itheta,izeta)*(GHat+iota*IHat))
          factor = (alpha*Delta*dPhiHatdpsiHat/four)*geometricFactor
                case (2) ! BHat_sub_theta / IHat
          geometricFactor = -dBHatdzeta(itheta,izeta)/(BHat(itheta,izeta)*(GHat+iota*IHat)) &
          - temp*iota/(BHat(itheta,izeta)*(GHat+iota*IHat)**2)
                  factor = (alpha*Delta*dPhiHatdpsiHat/four)*geometricFactor
                case (3) ! BHat_sub_zeta / GHat
          geometricFactor = dBHatdtheta(itheta,izeta)/(BHat(itheta,izeta)*(GHat+iota*IHat)) &
            - temp/(BHat(itheta,izeta)*(GHat+iota*IHat)**2)
                  factor = (alpha*Delta*dPhiHatdpsiHat/four)*geometricFactor
                case (4) ! BHat_sup_theta / iota
          geometricFactor = -IHat*temp/(BHat(itheta,izeta)*(GHat+iota*IHat)**2)
                  factor = (alpha*Delta*dPhiHatdpsiHat/four)*geometricFactor
              end select

                do ix=ixMin,Nx
                   do L=0,(Nxi_for_x(ix)-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)

                      ! Diagonal-in-L term
                      call MatSetValueSparse(dMatrixdLambda, rowIndex, rowIndex, &
                           (L+1)*L/((2*L-one)*(2*L+three))*factor, ADD_VALUES, ierr)

                       if (L<Nxi_for_x(ix)-2) then
                          ! Super-super-diagonal-in-L term:
                          ell = L+2
                          colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                          call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                               (L+3)*(L+2)*(L+1)/((two*L+5)*(2*L+three))*factor, ADD_VALUES, ierr)
                       end if

                       if (L>1) then
                          ! Sub-sub-diagonal-in-L term:
                          ell = L-2
                          colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                          call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                               -L*(L-1)*(L-2)/((2*L-3)*(2*L-one))*factor, ADD_VALUES, ierr)
                       end if
                   end do
                end do
             end do
          end do
       end if

    ! *********************************************************
    ! Add the sensitivity of the collisionless d/dx term associated with E_r
    ! *********************************************************

    ! DHat/BHat**3 * (BHat_sub_theta*dBHatdzeta-BHat_sub_zeta*dBHatdtheta)
    ! = 1/(BHat*(GHat+iota*IHat)) * (IHat*dBHatdzeta-GHat*dBHatdtheta)
    if (includeXDotTerm) then

      allocate(xPartOfXDot(Nx,Nx))
      allocate(xPartOfXDot_plus(Nx,Nx))
      allocate(xPartOfXDot_minus(Nx,Nx))
      allocate(rowIndices(Nx))
      allocate(colIndices(Nx))


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
              angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
              cos_angle = cos(angle)
              sin_angle = sin(angle)
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
                  case (0) ! Er
                    factor = -alpha*Delta/four
                    geometricFactor = DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                      * (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) &
                      - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))
                  case (1) ! BHat
                    dBHatdLambda = cos_angle
          dBHatdzetadLambda = n*Nperiods*sin_angle
          dBHatdthetadLambda = -m*sin_angle
          geometricFactor = -dBHatdLambda/(BHat(itheta,izeta)*BHat(itheta,izeta)*(GHat + iota*IHat)) &
              * (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta)) &
            + 1/(BHat(itheta,izeta)*(GHat+iota*IHat)) * (IHat*dBHatdzetadlambda-GHat*dBHatdthetadLambda)
                    factor = -alpha*Delta*dPhiHatdPsiHat/four
                  case (2) ! IHat
          geometricFactor = -iota/(BHat(itheta,izeta)*(GHat+iota*IHat)**2) &
            * (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta)) &
            + one/(BHat(itheta,izeta)*(GHat+iota*IHat)) * dBHatdzeta(itheta,izeta)
                    factor = -alpha*Delta*dPhiHatdPsiHat/four
                  case (3) ! GHat
          geometricFactor = -one/(BHat(itheta,izeta)*(GHat+iota*IHat)**2) &
            * (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta)) &
              - dBHatdtheta(itheta,izeta)/(BHat(itheta,izeta)*(GHat+iota*IHat))
                    factor = -alpha*Delta*dPhiHatdPsiHat/four
                  case (4) ! iota
          geometricFactor = -IHat/(BHat(itheta,izeta)*(GHat+iota*IHat)**2) &
            * (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta))
                    factor = -alpha*Delta*dPhiHatdPsiHat/four
                end select

               rowIndices = -1
               do ix=min_x_for_L(L),Nx
                  rowIndices(ix)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
               end do

               ! Term that is diagonal in L:
               colIndices = rowIndices
               stuffToAdd = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*factor*geometricFactor
               do ix_col=max(ixMinCol,min_x_for_L(L)),Nx
                  do ix_row=max(ixMin,min_x_for_L(L)),Nx
                     call MatSetValueSparse(dMatrixdLambda, rowIndices(ix_row), colIndices(ix_col), &
                          stuffToAdd*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                  end do
               end do

                ! Term that is super-super-diagonal in L:
                if (L<(Nxi-2)) then
                   ell = L + 2
                   stuffToAdd = (L+1)*(L+2)/((two*L+5)*(2*L+3))*factor*geometricFactor
                   do ix_col=max(ixMinCol,min_x_for_L(ell)),Nx
                      colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
                      do ix_row=max(ixMin,min_x_for_L(L)),Nx
                         call MatSetValueSparse(dMatrixdLambda, rowIndices(ix_row), colIndex, &
                              stuffToAdd*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                      end do
                   end do
                end if

                ! Term that is sub-sub-diagonal in L:
                if (L>1) then
                   ell = L - 2
                   stuffToAdd = L*(L-1)/((two*L-3)*(2*L-1))*factor*geometricFactor
                   do ix_col=max(ixMinCol,min_x_for_L(ell)),Nx
                      colIndex=getIndex(ispecies,ix_col,ell+1,itheta,izeta,BLOCK_F)
                      do ix_row=max(ixMin,min_x_for_L(L)),Nx
                         call MatSetValueSparse(dMatrixdLambda, rowIndices(ix_row), colIndex, &
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
    end if

!    end if ! comment out

    end do ! ispecies

    ! *********************************************************************************
    ! I do not need the collision operator, boundary conditions for f at x=0,
    ! or sources as they are independent of geometry
    ! *********************************************************************************

    ! *******************************************************************************
    ! Add the sensitivity of the density and pressure constraints:
    ! *******************************************************************************

    ! geometricFactor = 1/DHat = (GHat+iota*IHat)/(BHat*BHat)
    if (procThatHandlesConstraints) then
      L=0
      do itheta=1,Ntheta
         do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            sin_angle = sin(angle)
            geometricFactor = zero
            if (whichLambda == 1) then  ! BHat
              dBHatdLambda = cos_angle
              geometricFactor = -two*dBHatdLambda*(GHat+iota*IHat)/(BHat(itheta,izeta)**3)
            else if (whichLambda == 2) then ! IHat
              geometricFactor = iota/(BHat(itheta,izeta)*BHat(itheta,izeta))
            else if (whichLambda == 3) then ! GHat
              geometricFactor = 1/(BHat(itheta,izeta)*BHat(itheta,izeta))
            else if (whichLambda == 4) then ! iota
              geometricFactor = IHat/(BHat(itheta,izeta)*BHat(itheta,izeta))
            end if
            factor = thetaWeights(itheta)*zetaWeights(izeta)*geometricFactor

            do ix=1,Nx
               do ispecies=1,Nspecies
                  colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)

                  rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
                  call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                       x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)

                  rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
                  call MatSetValueSparse(dMatrixdLambda, rowIndex, colIndex, &
                       x2(ix)*x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)
               end do
            end do
         end do
      end do
    endif
    ! *******************************************************************************
    ! Done inserting values into the matrices.
    ! Now finalize the matrix
    ! *******************************************************************************

    call MatAssemblyBegin(dMatrixdLambda, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(dMatrixdLambda, MAT_FINAL_ASSEMBLY, ierr)

  end subroutine populatedMatrixdLambda
