#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif


  subroutine applyDenseTerms(inputVec, outputVec, whichMatrix)

    use petscvec
    use globalVariables 
    use indices

    implicit none

    Vec :: inputVec, outputVec
    integer, intent(in) :: whichMatrix

    PetscErrorCode :: ierr
    VecScatter :: VecScatterContext
    Vec :: inputVecLocal
    PetscScalar, pointer :: stateArray(:)
    PetscLogDouble :: time1, time2, time3, time4, time5, time6
    integer :: ispecies, L, ell, imn, ix, ix_row, ix_col, index
    real(prec) :: factor, LFactor
    real(prec), dimension(:), allocatable :: FourierVector, FourierVector2
    real(prec), dimension(:,:), allocatable :: xPartOfXDot


    if (masterProc) print *,"Applying dense terms for whichMatrix=",whichMatrix
    call PetscTime(time1,ierr)

    ! Send the entire input vector to each proc. I believe this is typically a small amount of memory
    ! compared to the memory used on each proc by the mumps factorization. Perhaps we
    ! could eventually move to a more efficient approach in which each proc does not 'know'
    ! the entire input vector.
    call VecScatterCreateToAll(inputVec, VecScatterContext, inputVecLocal, ierr)
    call VecScatterBegin(VecScatterContext, inputVec, inputVecLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, inputVec, inputVecLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterDestroy(VecScatterContext, ierr)
    
    call VecGetArrayF90(inputVecLocal, stateArray, ierr)

    !if (masterProc) print *," Beginning species",ispecies
    !call PetscTime(time3,ierr)
!!$       nHat = nHats(ispecies)
!!$       THat = THats(ispecies)
!!$       mHat = mHats(ispecies)
!!$       Z = Zs(ispecies)
!!$       sqrtTHat = sqrt(THat)
!!$       sqrtMHat = sqrt(mHat)

    allocate(FourierVector(NFourier2))
    allocate(FourierVector2(NFourier2))

    ! Note that the Fourier convolution matrices
    ! FourierMatrix_streaming
    ! FourierMatrix_ExB
    ! FourierMatrix_mirror
    ! FourierMatrix_xDot
    ! FourierMatrix_xiDot
    ! are generated in createGrids.F90

    ! *********************************************************
    ! Add the streaming d/dtheta and d/dzeta terms:
    ! *********************************************************
       
    call PetscTime(time3,ierr)
    time5=time3
    do ispecies = 1,Nspecies
       factor = sqrt(THats(ispecies)/mHats(ispecies))
       do ell=LMin,LMax
          do ix=1,Nx
             do imn=1,NFourier2
                index = getIndex(ispecies, ix, ell+1, imn, BLOCK_F) + 1 ! +1 to go from PETsc to Fortran indices
                FourierVector(imn) = stateArray(index)
             end do
             !FourierVector2 = factor*matmul(FourierMatrix_streaming,FourierVector) ! The BLAS2 dgemv call in the next line is a faster version of this "matmul" command.
             call dgemv('n',NFourier2,NFourier2,factor,FourierMatrix_streaming, &
                  NFourier2,FourierVector,1,0.0d+0,FourierVector2,1)

             ! Super-diagonal-in-L term
             !if (L < Nxi-1) then
             if (ell > 0) then
                !ell = L+1
                L = ell-1
                LFactor = (L+1)/(2*L+three)*x(ix)
                do imn=1,NFourier2
                   index = getIndex(ispecies, ix, L+1, imn, BLOCK_F)
                   call VecSetValue(outputVec, index,  &
                        LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                end do
             end if

             ! Sub-diagonal-in-L term
             !if (L > 0) then
             if (ell < Nxi-1) then
                !ell = L-1
                L = ell+1
                LFactor = L/(2*L-one)*x(ix)
                do imn = 1,NFourier2
                   index = getIndex(ispecies, ix, L+1, imn, BLOCK_F)
                   call VecSetValue(outputVec, index, &
                        LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                end do
             end if
          end do
       end do
    end do
    call PetscTime(time4,ierr)
    if (masterProc) print *," streaming:",time4-time3,"sec"
    time3=time4

    ! *********************************************************
    ! Add the ExB d/dtheta and d/dzeta terms:
    ! *********************************************************

    factor = alpha*Delta/two*dPhiHatdpsiHat
    do ispecies = 1,Nspecies
       do L=LMin,LMax
          do ix=1,Nx
             do imn = 1,NFourier2
                index = getIndex(ispecies, ix, L+1, imn, BLOCK_F) + 1 ! +1 to go from PETsc to Fortran indices
                FourierVector(imn) = stateArray(index)
             end do
             !FourierVector2 = factor*matmul(FourierMatrix_ExB,FourierVector) ! The BLAS2 dgemv call in the next line is a faster version of this "matmul" command.
             call dgemv('n',NFourier2,NFourier2,factor,FourierMatrix_ExB, &
                  NFourier2,FourierVector,1,0.0d+0,FourierVector2,1)

             do imn = 1,NFourier2
                index = getIndex(ispecies, ix, L+1, imn, BLOCK_F)
                call VecSetValue(outputVec, index, &
                     FourierVector2(imn), ADD_VALUES, ierr)
             end do
          end do
       end do
    end do
    call PetscTime(time4,ierr)
    if (masterProc) print *," ExB:      ",time4-time3,"sec"
    time3=time4

    ! *********************************************************
    ! Add the standard mirror term:
    ! *********************************************************

    do ispecies = 1,Nspecies
       factor = sqrt(THats(ispecies)/mHats(ispecies))
       do ell=LMin,LMax
          do ix=1,Nx
             do imn = 1,NFourier2
                index = getIndex(ispecies,ix,ell+1,imn,BLOCK_F) + 1 ! +1 to go from PETsc to Fortran indices
                FourierVector(imn) = stateArray(index)
             end do
             !FourierVector2 = factor*matmul(FourierMatrix,FourierVector) ! The BLAS2 dgemv call in the next line is a faster version of this "matmul" command.
             call dgemv('n',NFourier2,NFourier2,factor,FourierMatrix_mirror,&
                  NFourier2,FourierVector,1,0.0d+0,FourierVector2,1)

             !if (L<Nxi-1) then
             if (ell>0) then
                ! Super-diagonal-in-L term:    
                !ell = L+1
                L = ell-1
                LFactor = (L+1)*(L+2)/(2*L+three)*x(ix)
                do imn = 1,NFourier2
                   index = getIndex(ispecies,ix,L+1,imn,BLOCK_F)
                   call VecSetValue(outputVec,index,&
                        LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                end do
             end if

             !if (L>0) then
             if (ell<Nxi-1) then
                ! Sub-diagonal-in-L term:
                !ell = L-1
                L = ell + 1
                LFactor = -L*(L-1)/(2*L-one)*x(ix)
                do imn = 1,NFourier2
                   index = getIndex(ispecies,ix,L+1,imn,BLOCK_F)
                   call VecSetValue(outputVec,index,&
                        LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                end do
             end if
          end do
       end do
    end do
    call PetscTime(time4,ierr)
    if (masterProc) print *," mirror:   ",time4-time3,"sec"
    time3=time4

    ! *********************************************************
    ! Add the non-standard d/dxi term associated with E_r:
    ! *********************************************************
    
    if (includeElectricFieldTermInXiDot .and. abs(dPhiHatdpsiHat)>0) then
       factor = alpha*Delta*dPhiHatdpsiHat/4
       do ispecies=1,Nspecies
          do ell=LMin,LMax
             do ix=1,Nx
                do imn=1,NFourier2
                   index = getIndex(ispecies,ix,ell+1,imn,BLOCK_F) + 1 ! +1 to go from PETsc to Fortran indices
                   FourierVector(imn)=stateArray(index)
                end do
                !FourierVector2 = factor*matmul(FourierMatrix_xiDot,FourierVector) ! The BLAS2 dgemv call in the next line is a faster version of this "matmul" command.
                call dgemv('n',NFourier2,NFourier2,factor,FourierMatrix_xiDot, &
                     NFourier2,FourierVector,1,0.0d+0,FourierVector2,1)

                ! Diagonal-in-L term
                L = ell
                LFactor = (L+1)*L/((2*L-one)*(2*L+three))
                do imn=1,NFourier2
                   index = getIndex(ispecies,ix,ell+1,imn,BLOCK_F)
                   call VecSetValue(outputVec, index, &
                        LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                end do

                !if (L<Nxi-2) then
                if (ell>1) then
                   ! Super-super-diagonal-in-L term:
                   !ell = L+2
                   L = ell-2
                   LFactor = (L+3)*(L+2)*(L+1)/((two*L+5)*(2*L+three))
                   do imn=1,NFourier2
                      index=getIndex(ispecies,ix,L+1,imn,BLOCK_F)
                      call VecSetValue(outputVec, index, &
                           LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                   end do
                end if

                !if (L>1) then
                if (ell<Nxi-2) then
                   ! Sub-sub-diagonal-in-L term:
                   !ell = L-2
                   L = ell+2
                   LFactor = -L*(L-1)*(L-2)/((2*L-3)*(2*L-one))
                   do imn=1,NFourier2
                      index=getIndex(ispecies,ix,L+1,imn,BLOCK_F)
                      call VecSetValue(outputVec, index, &
                           LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                   end do
                end if

             end do
          end do
       end do
    end if
    call PetscTime(time4,ierr)
    if (masterProc) print *," Er xiDot: ",time4-time3,"sec"
    time3=time4


!!$       ! ****************************************************************
!!$       ! Add the non-standard d/dxi term associated with magnetic drifts:
!!$       ! ****************************************************************
!!$
!!$       if ((magneticDriftScheme>0) .and. (whichMatrix .ne. 2)) then
!!$          do itheta=ithetaMin,ithetaMax
!!$             do izeta=izetaMin,izetaMax
!!$
!!$                temp = (dBHat_sub_psi_dzeta(itheta,izeta) - dBHat_sub_zeta_dpsiHat(itheta,izeta)) * dBHatdtheta(itheta,izeta) &
!!$                     + (dBHat_sub_theta_dpsiHat(itheta,izeta) - dBHat_sub_psi_dtheta(itheta,izeta)) * dBHatdzeta(itheta,izeta)
!!$
!!$                if (.not. force0RadialCurrentInEquilibrium) then
!!$                   temp = temp + (dBHat_sub_zeta_dtheta(itheta,izeta) - dBHat_sub_theta_dzeta(itheta,izeta)) * dBHatdpsiHat(itheta,izeta)
!!$                end if
!!$
!!$                do ix=ixMin,Nx
!!$                   factor = -Delta*DHat(itheta,izeta)*THat*x(ix)*x(ix)/(2*Z*(BHat(itheta,izeta)**3)) * temp
!!$                     
!!$                   do L=0,(Nxi-1)
!!$                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
!!$
!!$                      ! Diagonal-in-L term
!!$                      call MatSetValueSparse(matrix, rowIndex, rowIndex, &
!!$                           (L+1)*L/((2*L-one)*(2*L+three))*factor, ADD_VALUES, ierr)
!!$
!!$                      ! Drop the off-by-2 diagonal terms in L if this is the preconditioner
!!$                      ! and preconditioner_xi = 1:
!!$                      if (whichMatrix .ne. 0 .or. preconditioner_xi==0) then
!!$
!!$                         if (L<Nxi-2) then
!!$                            ! Super-super-diagonal-in-L term:
!!$                            ell = L+2
!!$                            colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
!!$                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
!!$                                 (L+3)*(L+2)*(L+1)/((two*L+5)*(2*L+three))*factor, ADD_VALUES, ierr)
!!$                         end if
!!$
!!$                         if (L>1) then
!!$                            ! Sub-sub-diagonal-in-L term:
!!$                            ell = L-2
!!$                            colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
!!$                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
!!$                                 -L*(L-1)*(L-2)/((2*L-3)*(2*L-one))*factor, ADD_VALUES, ierr)
!!$                         end if
!!$                      end if
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$
!!$       end if

    ! *********************************************************
    ! Add the collisionless d/dx term associated with E_r
    ! *********************************************************

    if (includeXDotTerm .and. abs(dPhiHatdpsiHat)>0) then

!!$       if (force0RadialCurrentInEquilibrium) then
!!$          FourierMatrix2=0
!!$       else
!!$          call FourierTransform(2*factor*DHat*(dBHat_sub_zeta_dtheta - dBHat_sub_theta_dzeta)/(BHat*BHat), &
!!$               FourierVector)
!!$          call FourierConvolutionMatrix(FourierVector,FourierMatrix2,thresh)
!!$       end if

       allocate(xPartOfXDot(Nx,Nx))
       do ix=1,Nx
          xPartOfXDot(ix,:) = x(ix) * ddx(ix,:)
       end do

       factor = -alpha*Delta*dPhiHatdpsiHat/4
       do ispecies=1,Nspecies
          do ell=LMin,LMax
             do ix_col = 1,Nx
                do imn = 1,NFourier2
                   index = getIndex(ispecies,ix_col,ell+1,imn,BLOCK_F) + 1 ! +1 to go from PETsc to Fortran indices
                   FourierVector(imn) = stateArray(index)
                end do
                !FourierVector2 = factor*matmul(FourierMatrix_xDot,FourierVector) ! The BLAS2 dgemv call in the next line is a faster version of this "matmul" command.
                call dgemv('n',NFourier2,NFourier2,factor,FourierMatrix_xDot, &
                     NFourier2,FourierVector,1,0.0d+0,FourierVector2,1)

                ! Term that is diagonal in L:
                L=ell
                !stuffToAdd = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*FourierMatrix(imn_row,imn_col) &
                !     + (2*L*L+2*L-one)/((two*L+3)*(2*L-1))*FourierMatrix2(imn_row,imn_col)
                do ix_row=1,Nx
                   LFactor = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*xPartOfXDot(ix_row,ix_col)
                   do imn = 1,NFourier2
                      index=getIndex(ispecies,ix_row,L+1,imn,BLOCK_F)
                      call VecSetValue(outputVec,index, &
                           LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                   end do
                end do
                
                ! Term that is super-super-diagonal in L:
                !if (L<(Nxi-2)) then
                if (ell>1) then
                   !ell = L + 2
                   L = ell-2
                   !stuffToAdd = (L+1)*(L+2)/((two*L+5)*(2*L+3)) * &
                   !     (FourierMatrix(imn_row,imn_col)+FourierMatrix2(imn_row,imn_col))
                   do ix_row=1,Nx
                      LFactor = (L+1)*(L+2)/((two*L+5)*(2*L+3)) * xPartOfXDot(ix_row,ix_col)
                      do imn = 1,NFourier2
                         index=getIndex(ispecies,ix_row,L+1,imn,BLOCK_F)
                         call VecSetValue(outputVec,index, &
                              LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                      end do
                   end do
                end if
                
                ! Term that is sub-sub-diagonal in L:
                !if (L>1) then
                if (ell<Nxi-2) then
                   !ell = L - 2
                   L = ell + 2
                   !stuffToAdd = L*(L-1)/((two*L-3)*(2*L-1)) * &
                   !     (FourierMatrix(imn_row,imn_col)+FourierMatrix2(imn_row,imn_col))
                   do ix_row=1,Nx
                      LFactor = L*(L-1)/((two*L-3)*(2*L-1)) * xPartOfXDot(ix_row,ix_col)
                      do imn = 1,NFourier2
                         index=getIndex(ispecies,ix_row,L+1,imn,BLOCK_F)
                         call VecSetValue(outputVec,index, &
                              LFactor*FourierVector2(imn), ADD_VALUES, ierr)
                      end do
                   end do
                end if

             end do
          end do
       end do
       deallocate(xPartOfXDot)
    end if
    call PetscTime(time4,ierr)
    if (masterProc) print *," Er xDot:  ",time4-time3,"sec"
    time3=time4

    call VecRestoreArrayF90(inputVecLocal, stateArray, ierr)
    call VecDestroy(inputVecLocal, ierr)
    deallocate(FourierVector,FourierVector2)

    call VecAssemblyBegin(outputVec, ierr)
    call VecAssemblyEnd(outputVec, ierr)

    call PetscTime(time2,ierr)
    if (masterProc) print *,"Total time for applyDenseTerms:",time2-time1

  end subroutine applyDenseTerms
