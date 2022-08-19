! For compilers that do not include the error function erf(x), the line
! below should be un-commented:
!#define USE_GSL_ERF

module DKEMatrix

#include "PETScVersions.F90"

  use globalVariables
  use sparsify
  use indices
  use xGrid, only: xGrid_k
  use diagnostics
  
  implicit none

  private
  PetscScalar, dimension(:,:), allocatable :: ddthetaToUse, ddzetaToUse, ddxToUse, d2dx2ToUse
  PetscScalar :: adjointFactor
  PetscScalar, pointer :: stateArray(:)
  integer :: m,n ! only used for whichMatrix = 6

  logical, parameter :: localUsePhi1 = .true.

  public ::  populateMatrix
 
contains
  
  subroutine populateMatrix(matrix, whichMatrix, stateVec, whichLambda, whichMode)

    Mat  :: matrix
    integer, intent(in) :: whichMatrix, whichLambda, whichMode ! whichLambda, whichMode are only used for whichMatrix = 6
    Vec :: stateVec ! stateVec is ignored when includePhi1=false or when evaluating the residual.

    ! Allowed values for whichMatrix:
    ! 0 = preconditioner for Jacobian
    ! 1 = Jacobian
    ! 2 = matrix which multiplies f0 when evaluating the residual
    ! 3 = matrix which multiplies f1 when evaluating the residual
    ! 4 = Adjoint Jacobian
    ! 5 = preconditioner for adjoint Jacobian
    ! 6 = d Jacobian / d lambda

    PetscErrorCode :: ierr
     
    PetscLogDouble :: time1, time2
    double precision :: myMatInfo(MAT_INFO_SIZE)
    integer :: NNZ, NNZAllocated, NMallocs
    character(len=200) :: whichMatrixName, filename
    PetscViewer :: viewer
    VecScatter :: vecScatterContext
    Vec :: vecOnEveryProc
    logical :: useStateVec

    !!!!!!!!!!!!! from dMatrixdlambda !!!!!!!!!!!!
    PetscScalar, dimension(:,:), allocatable :: ddthetaToUse, ddzetaToUse, ddxToUse
    !!!!!!!!!!!!! END from dMatrixdlambda END !!!!!!!!!!!!


       
    if (whichMatrix==4 .or. whichMatrix==5) then
       adjointFactor = -one
    else
       adjointFactor = one
    end if

    ! *******************************************************************************
    ! Do a few sundry initialization tasks:
    ! *******************************************************************************

    ! This next line only matters for nonlinear calculations, in which the Mat objects for the matrix and preconditioner matrix 
    ! are reused at each iteration of SNES. In this case we need to clear all the previous entries, or else we would add the new
    ! values on top of the previous values:
    call MatZeroEntries(matrix,ierr)

    if (masterProc) then
       print *,"Running populateMatrix with whichMatrix = ",whichMatrix
       if (whichMatrix == 6) then
          print *,"   whichLambda = ",whichLambda
          if (.true.) then
             print *,"      whichMode = ",whichMode
          end if
       end if
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
    case (4)
       ! The matrix that is used for the adjoint solve for sensitivity computations
       ! This is the same matrix used to evaluate the residual for the adjoint solve
       whichMatrixName = "adjoint Jacobian" 
    case (5)
       ! The preconditioner for the matrix that is used for the adjoint solve for sensitivity computations
       whichMatrixName = "adjoint Jacobian preconditioner"
    case (6)
       ! The matrix used in adjoint diagnostics to calculate sensitivities to the parameter lambda
       whichMatrixName = "d Jacobian / d lambda"
    case default
       if (masterProc) then
          print *,"Error! whichMatrix must be 0, 1, 2, 3, 4, 5, or 6."
       end if
       stop
    end select

    call PetscTime(time1, ierr)

    call setDiagonalsToZero(matrix)

    if (whichMatrix == 6) then
       if (whichLambda == 1) then
          m = ms_sensitivity(whichMode)
          n = ns_sensitivity(whichMode)
       end if
    end if
    


    if (whichMatrix < 6) then
       useStateVec = (includePhi1 .and. (.not. readExternalPhi1) .and. (whichMatrix==0 .or. whichMatrix==1)) !!Added by AM 2018-12

       call shiftDiagonals(matrix, whichMatrix)
       
       if (useStateVec) then
          ! We need delta f to evaluate the Jacobian, so send a copy to every proc:
          call VecScatterCreateToAll(stateVec, vecScatterContext, vecOnEveryProc, ierr)
          call VecScatterBegin(vecScatterContext, stateVec, vecOnEveryProc, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call VecScatterEnd(vecScatterContext, stateVec, vecOnEveryProc, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call VecGetArrayF90(vecOnEveryProc, stateArray, ierr)
       end if

       ! In nonlinear runs, the Jacobian and residual require Phi1:
       if (includePhi1 .and. (.not. readExternalPhi1) .and. (whichMatrix .ne. 2)) then !!Added by AM 2018-12
          call extractPhi1(stateVec)
       end if
    end if

    ! *********************************************************
    ! Allocate module variables
    ! *********************************************************
    allocate(ddthetaToUse(Ntheta, Ntheta))
    allocate(ddzetaToUse(Nzeta, Nzeta))
    allocate(ddxToUse(Nx,Nx))
    if (whichMatrix < 6) then
       allocate(d2dx2ToUse(Nx,Nx))
    end if
    
    ! ************************************************************
    ! ************************************************************
    ! Begin adding the collisionless terms of the kinetic equation
    ! ************************************************************
    ! ************************************************************


    call streamingDdtheta(matrix, whichMatrix, whichLambda)

    call streamingDdzeta(matrix, whichMatrix, whichLambda)

    call ExBDdtheta(matrix, whichMatrix, whichLambda)

    call ExBDdzeta(matrix, whichMatrix, whichLambda)

    ! not implemented for adjoint
    call magneticDriftDdtheta(matrix, whichMatrix)

    ! not implemented for adjoint
    call magneticDriftDdzeta(matrix, whichMatrix)

    call standardMirror(matrix, whichMatrix, whichLambda)

    call nonstandardErDdxi(matrix, whichMatrix, whichLambda)

    ! not implemented for adjoint
    call nonstandardMagneticDriftDdxi(matrix, whichMatrix)

    call collisionlessErDdx(matrix, whichMatrix, whichLambda)

    ! nothing to implement for adjoint diagnostics
    call adjointErFullTrajectories(matrix, whichMatrix)

    ! not implemented for adjoint.
    ! only needed for non-linear
    call radialExBDrive(matrix, whichMatrix)
    
    call parallelE(matrix, whichMatrix, whichLambda)

    ! not implemented for adjoint.
    ! only needed for non-linear
    call nonlinearParallelE(matrix, whichMatrix)

    ! *********************************************************
    ! *********************************************************
    !
    ! End of adding the collisionless kinetic terms.
    !
    ! *********************************************************
    ! *********************************************************

    if (whichMatrix .ne. 6) then
       ! nothing to implement for adjoint diagnostics
       call addCollisionOperator(matrix, whichMatrix)

       ! nothing to implement for adjoint diagnostics
       call boundaryConditionAtX0(matrix, whichMatrix)

       ! nothing to implement for adjoint diagnostics
       call addSources(matrix, whichMatrix)
    end if


    call densityPressureConstraints(matrix, whichMatrix, whichLambda)

    if (whichMatrix .ne. 6) then
       ! not implemented for adjoint.
       ! only needed for non-linear
       call quasineutrality(matrix, whichMatrix)

       ! not implemented for adjoint.
       ! only needed for non-linear
       call phi1Constraints(matrix, whichMatrix)
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

    if (whichMatrix /= 6) then

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

    end if
    ! *********************************************************
    ! Deallocate module variables
    ! *********************************************************
    deallocate(ddthetaToUse)
    deallocate(ddzetaToUse)
    deallocate(ddxToUse)
    if (whichMatrix < 6) then
       deallocate(d2dx2ToUse)
    end if
 
  end subroutine populateMatrix

! -----------------------------------------------------------------------
! -----------------------------------------------------------------------

  subroutine setDiagonalsToZero(matrix)
    ! Sometimes PETSc complains if any of the diagonal elements are not set.
    ! Therefore, set the entire diagonal to 0 to be safe.
    Mat :: matrix
    PetscErrorCode :: ierr
    integer :: i
    if (masterProc) then
       do i=1,matrixSize
          call MatSetValue(matrix, i-1, i-1, zero, ADD_VALUES, ierr)
       end do
    end if
    
  end subroutine setDiagonalsToZero

  subroutine shiftDiagonals(matrix, whichMatrix)
    ! Since PETSc's direct sparse solver complains if there are any zeros on the diagonal
    ! (unlike mumps or superlu_dist), then if we're using this solver
    ! add some values to the diagonals of the preconditioner.  By trial-and-error, I found it works
    ! best to shift the diagonal of the quasineutrality and constraint blocks but not for the kinetic-equation block.
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    PetscScalar :: temp
    integer :: ispecies, ix, itheta, izeta, index

    if ((.not. isAParallelDirectSolverInstalled) .and. masterProc .and. (whichMatrix==0 .or. whichMatrix==5)) then
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
       !!if (includePhi1) then
 !!Commented by AM 2018-12
       if (includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
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
    
  end subroutine shiftDiagonals

  subroutine streamingDdtheta(matrix, whichMatrix, whichLambda)
    ! Add the streaming d/dtheta term
    Mat :: matrix
    integer, intent(in) :: whichMatrix, whichLambda
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat
    integer :: ispecies, ix, itheta, izeta, L, ell 
    PetscScalar, dimension(:,:), allocatable :: thetaPartOfTerm, localThetaPartOfTerm
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar :: geometricFactor, angle, cos_angle, dBHatdLambda

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)
       if (whichMatrix .ne. 2) then
          allocate(thetaPartOfTerm(Ntheta,Ntheta))
          allocate(localThetaPartOfTerm(Ntheta,localNtheta))
          allocate(rowIndices(localNtheta))
          allocate(colIndices(Ntheta))
          do L=0,(Nxi-1)

             if (((whichMatrix>0) .and. (whichMatrix .ne. 5)) .or. L < preconditioner_theta_min_L .or. (whichMatrix == 6)) then
                ddthetaToUse = ddtheta
             else
                ddthetaToUse = ddtheta_preconditioner
             end if

             do izeta=izetaMin,izetaMax
                do itheta=1,Ntheta
                   if (whichMatrix < 6) then
                      ! Normal matrix
                      geometricFactor = BHat_sup_theta(itheta,izeta)/BHat(itheta,izeta)
                   else
                      ! d matrix /d lambda
                      select case(whichLambda)
                      case (0) ! Er
                         geometricFactor = zero
                      case (1) ! BHat
                         angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                         cos_angle = cos(angle)
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
                   end if
                   thetaPartOfTerm(itheta,:) = adjointFactor * sqrtTHat/sqrtMHat &
                        * ddthetaToUse(itheta,:) * geometricFactor
                end do
                
                ! PETSc uses the opposite convention to Fortran:
                thetaPartOfTerm = transpose(thetaPartOfTerm)
                localThetaPartOfTerm = thetaPartOfTerm(:,ithetaMin:ithetaMax)
                
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
    end do
  end subroutine streamingDdtheta
  
  subroutine streamingDdzeta(matrix, whichMatrix, whichLambda)
    ! Add the streaming d/dzeta term
    
    Mat :: matrix
    integer, intent(in) :: whichMatrix, whichLambda
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat
    integer :: ispecies, ix, itheta, izeta, L, ell 
    PetscScalar, dimension(:,:), allocatable :: zetaPartOfTerm, localZetaPartOfTerm
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar :: geometricFactor, angle, cos_angle, dBHatdLambda

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)

       if (whichMatrix .ne. 2) then
          allocate(zetaPartOfTerm(Nzeta,Nzeta))
          allocate(localZetaPartOfTerm(Nzeta,localNzeta))
          allocate(rowIndices(localNzeta))
          allocate(colIndices(Nzeta))
          do L=0,(Nxi-1)

             if (((whichMatrix>0) .and. (whichMatrix .ne. 5)) .or. L < preconditioner_zeta_min_L) then
                ddzetaToUse = ddzeta
             else
                ddzetaToUse = ddzeta_preconditioner
             end if

             do itheta=ithetaMin, ithetaMax
                do izeta=1,Nzeta
                   if (whichMatrix < 6) then
                      ! Normal matrix
                      geometricFactor = BHat_sup_zeta(itheta,izeta)/BHat(itheta,izeta)
                   else
                      ! d matrix /d lambda
                      select case(whichLambda)
                      case (0) ! Er
                         geometricFactor = zero
                      case (1) ! BHat
                         angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                         cos_angle = cos(angle)                         
                         dBHatdLambda = cos_angle
                         geometricFactor = dBHatdLambda/(GHat+iota*IHat)
                      case (2) ! IHat
                         geometricFactor = -iota*BHat(itheta,izeta)/(GHat+iota*IHat)**2
                      case (3) ! GHat
                         geometricFactor = -BHat(itheta,izeta)/(GHat+iota*IHat)**2
                      case (4) ! iota
                         geometricFactor = -IHat*BHat(itheta,izeta)/(GHat+iota*IHat)
                      end select
                   end if
                  
                   zetaPartOfTerm(izeta,:) = adjointFactor * sqrtTHat/sqrtMHat  &
                        * ddzetaToUse(izeta,:) * geometricFactor
                end do
                
                ! PETSc uses the opposite convention to Fortran:
                zetaPartOfTerm = transpose(zetaPartOfTerm)
                localZetaPartOfTerm = zetaPartOfTerm(:,izetaMin:izetaMax)
                
                !do ix=ixMin,Nx
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
    end do
  end subroutine streamingDdzeta


  ! *********************************************************
  ! Add the ExB d/dtheta term:
  ! *********************************************************
  subroutine ExBDdtheta(matrix, whichMatrix, whichLambda)
    Mat :: matrix
    integer, intent(in) :: whichMatrix, whichLambda
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor
    integer :: ispecies, ix, itheta, izeta, L
    PetscScalar, dimension(:,:), allocatable :: thetaPartOfTerm, localThetaPartOfTerm
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar :: geometricFactor, angle, cos_angle, dBHatdLambda, dVPrimeHatdLambda, dVPrimeHatdBmn, dFSABHat2dBmn

    if (whichMatrix == 6 .and. whichLambda==1) then ! BHat sensitivity
       dVPrimeHatdLambda = zero
       dVPrimeHatdBmn  = zero
        do itheta=1,Ntheta
          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            dVPrimeHatdLambda = dVPrimeHatdLambda - two*thetaWeights(itheta)*zetaWeights(izeta)*cos_angle/(DHat(itheta,izeta)*BHat(itheta,izeta))
            if (useDKESExBDrift) then
            dVPrimeHatdBmn = dVPrimeHatdBmn - (two*(GHat+iota*IHat)) * thetaWeights(itheta) * zetaWeights(izeta) &
                 * BHat(itheta,izeta)**(-3) * cos_angle
         end if
         end do
      end do
      if (useDKESExBDrift) then
         dFSABHat2dBmn = -4*pi*pi*(GHat+iota*IHat)*dVPrimeHatdBmn/(VPrimeHat**2)
      end if
   end if

   
  
    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)

       if (whichMatrix .ne. 2) then
          allocate(thetaPartOfTerm(Ntheta,Ntheta))
          allocate(localThetaPartOfTerm(Ntheta,localNtheta))
          allocate(rowIndices(localNtheta))
          allocate(colIndices(Ntheta))
          do L=0,(Nxi-1)

             if (ExBDerivativeSchemeTheta==0) then
                if (((whichMatrix>0) .and. (whichMatrix .ne. 5)) .or. L < preconditioner_theta_min_L) then
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
                      if (whichMatrix < 6) then
                         ! Normal matrix
                         geometricFactor = DHat(itheta,izeta) * BHat_sub_zeta(itheta,izeta)/ FSABHat2
                         factor = alpha*Delta/two * dPhiHatdpsiHat
                      else
                         select case(whichLambda)
                         case (0) ! Er
                            geometricFactor = DHat(itheta,izeta)*BHat_sub_zeta(itheta,izeta)/(FSABHat2)
                            factor = alpha*Delta/two
                         case (1) ! BHat
                            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                            cos_angle = cos(angle)                                                     
                            dBHatdLambda = cos_angle
                            geometricFactor = (GHat/(GHat + iota*IHat))*(two*BHat(itheta,izeta)*dBHatdLambda/FSABHat2 &
                                 - BHat(itheta,izeta)*BHat(itheta,izeta)*dFSABHat2dBmn/(FSABHat2**2))
                            factor = alpha*Delta/two*dPhiHatdpsiHat
                         case (2) ! IHat
                            geometricFactor = -(GHat*iota/(GHat+iota*IHat)**2)*(BHat(itheta,izeta)*BHat(itheta,izeta)/FSABHat2)
                            factor = alpha*Delta/two*dPhiHatdpsiHat
                         case (3) ! GHat
                            geometricFactor = (one/(GHat+iota*IHat)-GHat/(GHat+iota*IHat)**2) &
                                 * (BHat(itheta,izeta)*BHat(itheta,izeta)/FSABHat2)
                            factor = alpha*Delta/two*dPhiHatdpsiHat
                         case (4) ! iota
                            geometricFactor = (-GHat*IHat/(GHat+iota*IHat)**2)*(BHat(itheta,izeta)*BHat(itheta,izeta)/FSABHat2)
                            factor = alpha*Delta/two*dPhiHatdpsiHat
                         end select
                      end if
                   
                      thetaPartOfTerm(itheta,:) = adjointFactor*ddthetaToUse(itheta,:) * geometricFactor * factor
                   end do
                else
                   do itheta=1,Ntheta
                      if (whichMatrix < 6) then
                         ! Normal matrix
                         geometricFactor = DHat(itheta,izeta)*BHat_sub_zeta(itheta,izeta)/(BHat(itheta,izeta)*BHat(itheta,izeta))
                         factor = alpha*Delta/two*dPhiHatdpsiHat
                      else
                         select case(whichLambda)
                         case (0) ! Er
                            geometricFactor = DHat(itheta,izeta)*BHat_sub_zeta(itheta,izeta)/(BHat(itheta,izeta)*BHat(itheta,izeta))
                            factor = alpha*Delta/two
                         case (1) ! BHat
                            geometricFactor = zero
                            factor = alpha*Delta/two*dPhiHatdpsiHat
                         case (2) ! IHat
                            geometricFactor = -GHat*iota/(GHat+iota*IHat)**2
                            factor = alpha*Delta/two*dPhiHatdpsiHat
                         case (3) ! GHat
                            geometricFactor = (one/(GHat+iota*IHat)-GHat/(GHat+iota*IHat)**2)
                            factor = alpha*Delta/two*dPhiHatdpsiHat
                         case (4) ! iota
                            geometricFactor = (-GHat*IHat/(GHat+iota*IHat)**2)
                            factor = alpha*Delta/two*dPhiHatdpsiHat
                         end select
                      end if
                         
                         
                      thetaPartOfTerm(itheta,:) = adjointFactor*ddthetaToUse(itheta,:) * geometricFactor * factor
                   end do
                end if
                
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
    end do
  end subroutine ExBDdtheta
  
  ! *********************************************************
  ! Add the ExB d/dzeta term:
  ! *********************************************************
  subroutine ExBDdzeta(matrix, whichMatrix, whichLambda)
    Mat :: matrix
    integer, intent(in) :: whichMatrix, whichLambda
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor
    integer :: ispecies, ix, itheta, izeta, L
    PetscScalar, dimension(:,:), allocatable :: zetaPartOfTerm, localZetaPartOfTerm
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar :: geometricFactor, dBHatdLambda

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)


       if (whichMatrix .ne. 2) then
          factor = -alpha*Delta/two*dPhiHatdpsiHat
          allocate(zetaPartOfTerm(Nzeta,Nzeta))
          allocate(localZetaPartOfTerm(Nzeta,localNzeta))
          allocate(rowIndices(localNzeta))
          allocate(colIndices(Nzeta))
          do L=0,(Nxi-1)

             if (ExBDerivativeSchemeZeta==0) then
                if (((whichMatrix>0) .and. (whichMatrix .ne. 5)) .or. L < preconditioner_zeta_min_L) then
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
                   ! geometryFactor = DHat(itheta,izeta)*BHat_sub_theta(itheta,izeta)/FSABHat2 = (IHat/(GHat+iota*IHat))*(BHat*BHat/FSABHat2)
                   do izeta=1,Nzeta
                      if (whichMatrix < 6) then
                         ! Normal matrix
                         geometricFactor = DHat(itheta,izeta) * BHat_sub_theta(itheta,izeta) / FSABHat2
                         factor = -alpha*Delta/two * dPhiHatdpsiHat
                      else
                         select case(whichLambda)
                         case (0) ! Er
                            geometricFactor = DHat(itheta,izeta)*BHat_sub_theta(itheta,izeta)/FSABHat2
                            factor = -alpha*Delta/two
                         case (1) ! BHat
                            geometricFactor = (2*BHat(itheta,izeta)*dBHatdLambda/FSABHat2)*(IHat/(GHat+iota*IHat))
                            factor = -alpha*Delta/two*dPhiHatdpsiHat
                         case (2) ! IHat
                            geometricFactor = (one/(GHat+iota*IHat)-(IHat*iota/(GHat+iota*IHat)**2)) &
                                 *(BHat(itheta,izeta)*BHat(itheta,izeta))/FSABHat2
                            factor = -alpha*Delta/two*dPhiHatdpsiHat
                         case (3) ! GHat
                            geometricFactor = -(IHat/(GHat+iota*IHat)**2)*(BHat(itheta,izeta)*BHat(itheta,izeta))/FSABHat2
                            factor = -alpha*Delta/two*dPhiHatdpsiHat
                         case (4) ! iota
                            geometricFactor = -(IHat*IHat/(GHat+iota*IHat)**2)*(BHat(itheta,izeta)*BHat(itheta,izeta)/FSABHat2)
                            factor = -alpha*Delta/two*dPhiHatdpsiHat
                         end select
                         
                      end if
                      zetaPartOfTerm(izeta,:) = adjointFactor*ddzetaToUse(izeta,:)  * geometricFactor * factor
                   end do
                else
                   ! geometryFactor = DHat(itheta,izeta)*BHat_sub_theta(itheta,izeta)/(BHat(itheta,izeta)*BHat(itheta,izeta)) = IHat/(GHat+iota*IHat)
                   do izeta=1,Nzeta
                      if (whichMatrix < 6) then
                         ! Normal matrix
                         geometricFactor = DHat(itheta,izeta) * BHat_sub_theta(itheta,izeta)/ (BHat(itheta,izeta) ** 2)
                         factor = alpha*Delta/two * dPhiHatdpsiHat
                      else
                         select case(whichLambda)
                         case (0) ! Er
                            geometricFactor = DHat(itheta,izeta)*BHat_sub_theta(itheta,izeta)/(BHat(itheta,izeta)*BHat(itheta,izeta))
                            factor = -alpha*Delta/two
                         case (1) ! BHat
                            geometricFactor = 0
                            factor = -alpha*Delta/two*dPhiHatdpsiHat
                         case (2) ! IHat
                            geometricFactor = (one/(GHat+iota*IHat)-(IHat*iota/(GHat+iota*IHat)**2))
                            factor = -alpha*Delta/two*dPhiHatdpsiHat
                         case (3) ! GHat
                            geometricFactor = -(IHat/(GHat+iota*IHat)**2)
                            factor = -alpha*Delta/two*dPhiHatdpsiHat
                         case (4) ! iota
                            geometricFactor = -(IHat*IHat/(GHat+iota*IHat)**2)
                            factor = -alpha*Delta/two*dPhiHatdpsiHat
                         end select
                      end if
                      zetaPartOfTerm(izeta,:) = adjointFactor*ddzetaToUse(izeta,:) * geometricFactor * factor
                   end do
                end if
                
                ! PETSc uses the opposite convention to Fortran:
                zetaPartOfTerm = transpose(zetaPartOfTerm)
                localzetaPartOfTerm = zetaPartOfTerm(:,izetaMin:izetaMax)
                
                !do ix=ixMin,Nx
                do ix=max(ixMin,min_x_for_L(L)),Nx
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
    end do
       
  end subroutine ExBDdzeta
  
  ! *********************************************************
  ! Add the magnetic drift d/dtheta term:
  ! *********************************************************
  subroutine magneticDriftDdtheta(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor, geometricFactor1, geometricFactor2, geometricFactor3, stuffToAdd
    integer :: ispecies, ix, itheta, izeta, L, ell, izetaRow, izetaCol, ithetaRow, ithetaCol, maxL, rowIndex, colIndex
    
    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)

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
                   select case (magneticDriftScheme)
                   case (1,2,7,8)
                      geometricFactor1 = (BHat_sub_zeta(ithetaRow,izeta)*dBHatdpsiHat(ithetaRow,izeta) &
                      - BHat_sub_psi(ithetaRow,izeta)*dBHatdzeta(ithetaRow,izeta))

                      geometricFactor2 = 2.0 * BHat(ithetaRow,izeta) &
                      * (dBHat_sub_psi_dzeta(ithetaRow,izeta) - dBHat_sub_zeta_dpsiHat(ithetaRow,izeta))
                   case (3,4)
                      geometricFactor1 = 0.0
                      geometricFactor2 = 0.0
                   case (5)
                      geometricFactor1 = BHat_sub_zeta(ithetaRow,izeta) * gradpsidotgradB_overgpsipsi(ithetaRow,izeta)
                      geometricFactor2 = - 2.0 * geometricFactor1
                   case (6)
                      geometricFactor1 = BHat_sub_zeta(ithetaRow,izeta) * gradpsidotgradB_overgpsipsi(ithetaRow,izeta)
                      geometricFactor2 = BHat_sub_zeta(ithetaRow,izeta) * 2.0 * pPrimeHat / BHat(ithetaRow,izeta) !unregularized
                   case (9)
                      geometricFactor1 = (BHat_sub_zeta(ithetaRow,izeta)*dBHatdpsiHat(ithetaRow,izeta) &
                      - BHat_sub_psi(ithetaRow,izeta)*dBHatdzeta(ithetaRow,izeta))

                      geometricFactor2 = 2.0 * BHat(ithetaRow,izeta) &
                      * (dBHat_sub_psi_dzeta(ithetaRow,izeta) - dBHat_sub_zeta_dpsiHat(ithetaRow,izeta) &
                      + diotadpsiHat*GHat*BHat_sup_theta(itheta,izetaRow)/(DHat(itheta,izetaRow)*iota*iota) )
                      ! Need to add shear term here!
                   case default
                      stop "Invalid magneticDriftScheme in d/dtheta term"
                   end select

                   if (magneticDriftScheme==2) then
                      geometricFactor3 = BDotCurlB(ithetaRow,izeta)*BHat_sup_theta(ithetaRow,izeta) &
                           /(BHat(ithetaRow,izeta)*DHat(ithetaRow,izeta))
                   else
                      geometricFactor3 = 0.0
                   end if

                   if (magneticDriftDerivativeScheme==0) then
                      if (((whichMatrix>0) .and. (whichMatrix .ne. 5)) .or. L < preconditioner_theta_min_L) then
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

                   !do ix = ixMin, Nx
                   do ix = max(ixMin,min_x_for_L(L)), Nx
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
                            if (L < Nxi_for_x(ix)-2) then
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
    end do
  end subroutine magneticDriftDdtheta

  ! *********************************************************
  ! Add the magnetic drift d/dzeta term:
  ! *********************************************************
  subroutine magneticDriftDdzeta(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor, geometricFactor1, geometricFactor2, geometricFactor3, stuffToAdd
    integer :: ispecies, ix, itheta, izeta, L, ell, izetaRow, izetaCol, ithetaRow, ithetaCol, maxL, rowIndex, colIndex
    
    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)
      
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

                   select case (magneticDriftScheme)
                   case (1,2,7,9)
                      geometricFactor1 = (BHat_sub_psi(itheta,izetaRow)*dBHatdtheta(itheta,izetaRow) &
                           - BHat_sub_theta(itheta,izetaRow)*dBHatdpsiHat(itheta,izetaRow))
                      
                      geometricFactor2 = 2.0 * BHat(itheta,izetaRow) &
                           * (dBHat_sub_theta_dpsiHat(itheta,izetaRow) - dBHat_sub_psi_dtheta(itheta,izetaRow))
                   case (3) 
                      geometricFactor1 = &
                           (BHat_sub_psi(itheta,izetaRow)* &
                           (dBHatdtheta(itheta,izetaRow) + dBHatdzeta(itheta,izetaRow)/iota)&
                           - (BHat_sub_theta(itheta,izetaRow)+BHat_sub_zeta(itheta,izetaRow)/iota) &
                           *dBHatdpsiHat(itheta,izetaRow))
                      
                      geometricFactor2 = 2.0 * BHat(itheta,izetaRow) &
                           * (dBHat_sub_theta_dpsiHat(itheta,izetaRow) + dBHat_sub_zeta_dpsiHat(itheta,izetaRow)/iota &
                           - (dBHat_sub_psi_dtheta(itheta,izetaRow)+dBHat_sub_psi_dzeta(itheta,izetaRow)/iota))
                           !- (dBHat_sub_psi_dtheta(itheta,izetaRow)+dBHat_sub_psi_dtheta(itheta,izetaRow)/iota))
                   case (4)
                      geometricFactor1 = &
                           (BHat_sub_psi(itheta,izetaRow)* &
                           (dBHatdtheta(itheta,izetaRow) + dBHatdzeta(itheta,izetaRow)/iota)&
                           - (BHat_sub_theta(itheta,izetaRow)+BHat_sub_zeta(itheta,izetaRow)/iota) &
                           *dBHatdpsiHat(itheta,izetaRow))
                      
                      geometricFactor2 = 2.0 * BHat(itheta,izetaRow) &
                           * (dBHat_sub_theta_dpsiHat(itheta,izetaRow) + dBHat_sub_zeta_dpsiHat(itheta,izetaRow)/iota  &
                           - diotadpsiHat / (iota*iota) * BHat_sub_zeta(itheta,izetaRow) &
                           - (dBHat_sub_psi_dtheta(itheta,izetaRow)+dBHat_sub_psi_dzeta(itheta,izetaRow)/iota))
                           !- (dBHat_sub_psi_dtheta(itheta,izetaRow)+dBHat_sub_psi_dtheta(itheta,izetaRow)/iota))
                   case (5)
                      geometricFactor1 = -BHat_sub_theta(itheta,izetaRow) * gradpsidotgradB_overgpsipsi(itheta,izetaRow)
                      geometricFactor2 = - 2.0 * geometricFactor1
                   case (6)
                      geometricFactor1 = -BHat_sub_theta(itheta,izetaRow) * gradpsidotgradB_overgpsipsi(itheta,izetaRow)
                      geometricFactor2 = -BHat_sub_theta(itheta,izetaRow) * 2.0 * pPrimeHat / BHat(itheta,izetaRow) !unregularized
                   case (8)
                      geometricFactor1 = (BHat_sub_psi(itheta,izetaRow)*dBHatdtheta(itheta,izetaRow) &
                           - BHat_sub_theta(itheta,izetaRow)*dBHatdpsiHat(itheta,izetaRow))
                      
                      geometricFactor2 = 2.0 * BHat(itheta,izetaRow) &
                           * (dBHat_sub_theta_dpsiHat(itheta,izetaRow) - dBHat_sub_psi_dtheta(itheta,izetaRow) &
                           -diotadpsiHat*GHat*BHat_sup_theta(itheta,izetaRow)/(DHat(itheta,izetaRow)*iota*iota*iota) )
                   case default
                      stop "Invalid magneticDriftScheme in d/dzeta term"
                   end select

                   if (magneticDriftScheme==2) then
                      geometricFactor3 = BDotCurlB(itheta,izetaRow)*BHat_sup_zeta(itheta,izetaRow) &
                           /(BHat(itheta,izetaRow)*DHat(itheta,izetaRow))
                   else
                      geometricFactor3 = 0.0
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

                   !do ix = ixMin, Nx
                   do ix = max(ixMin,min_x_for_L(L)), Nx
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
                            if (L < Nxi_for_x(ix)-2) then
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
    end do
  end subroutine magneticDriftDdzeta

  ! *********************************************************
  ! Add the standard mirror term:
  ! *********************************************************
  subroutine standardMirror(matrix, whichMatrix, whichLambda)
    Mat :: matrix
    integer, intent(in) :: whichMatrix, whichLambda
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor, geometricFactor
    integer :: ispecies, ix, itheta, izeta, L, ell, rowIndex, colIndex
    PetscScalar :: angle, cos_angle, sin_angle, dBHatdLambda, dBHatdthetadLambda, dBHatdzetadLambda

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)

       if (whichMatrix .ne. 2) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax
                if (whichMatrix < 6) then
                   ! Normal matrix
                   geometricFactor = (BHat_sup_theta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                        + BHat_sup_zeta(itheta,izeta) * dBHatdzeta(itheta,izeta)) &
                        /(BHat(itheta,izeta)*BHat(itheta,izeta))
                else
                   ! d matrix /d lambda
                   select case(whichLambda)
                   case (0) ! Er
                      geometricFactor = zero
                   case (1) ! BHat
                      angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                      cos_angle = cos(angle)
                      sin_angle = sin(angle)                   
                      dBHatdLambda = cos_angle
                      dBHatdthetadLambda = -m*sin_angle
                      dBHatdzetadLambda = n*Nperiods*sin_angle
                      geometricFactor = (1/(GHat+iota*IHat))*(iota*dBHatdthetadLambda + dBHatdzetadLambda)
                   case (2) ! IHat
                      geometricFactor = -(iota/(GHat+iota*IHat)**2)*(iota*dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta))
                   case (3) ! GHat
                      geometricFactor = -(one/(GHat+iota*IHat)**2)*(iota*dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta))
                   case (4) ! iota
                      geometricFactor = -IHat/(GHat+iota*IHat)**2*(iota*dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta)) &
                           + (dBHatdtheta(itheta,izeta)/(GHat+iota*IHat))
                   end select
                end if
                   
                factor = -adjointFactor * geometricFactor * sqrtTHat/(two*sqrtMHat)
                
                do ix=ixMin,Nx
                   do L=0,(Nxi_for_x(ix)-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                      
                      if (L<Nxi_for_x(ix)-1) then
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
    end do
  end subroutine standardMirror

  ! *********************************************************
  ! Add the non-standard d/dxi term associated with E_r:
  ! *********************************************************
  subroutine nonstandardErDdxi(matrix, whichMatrix, whichLambda)
    Mat :: matrix
    integer, intent(in) :: whichMatrix, whichLambda
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor, temp, geometricFactor
    integer :: ispecies, ix, itheta, izeta, L, ell, rowIndex, colIndex
    PetscScalar :: angle, cos_angle, sin_angle, dBHatdLambda, dBHatdthetadLambda, dBHatdzetadLambda, dTempdLambda

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)

       if (includeElectricFieldTermInXiDot .and. (whichMatrix .ne. 2)) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax

                temp = adjointFactor*(BHat_sub_zeta(itheta,izeta) * dBHatdtheta(itheta,izeta) &
                     - BHat_sub_theta(itheta,izeta) * dBHatdzeta(itheta,izeta))

                if (.not. force0RadialCurrentInEquilibrium) then
                   ! this should never happen, but add check anyway
                   if (whichMatrix == 6) then
                      if (masterProc) then
                         print *, "force0RadialCurrentInEquilibrium must be true for adjoint diagnostics"
                      end if
                   end if
                   temp = temp - 2 * BHat(itheta,izeta) &
                        * (dBHat_sub_zeta_dtheta(itheta,izeta) - dBHat_sub_theta_dzeta(itheta,izeta))
                end if

                ! geometricFactor = 
                ! (BHat_sub_zeta(itheta,izeta) * dBHatdtheta(itheta,izeta) &
                !   - BHat_sub_theta(itheta,izeta) * dBHatdzeta(itheta,izeta))*DHat/(BHat(itheta,izeta)**3)
                ! = (GHat*dBHatdtheta-IHat*dBHatdzeta)*(1/(BHat*(GHat+iota*IHat)))
                ! factor = (alpha*Delta/four)* dPhiHatdpsiHat *geometricFactor
                if (whichMatrix < 6) then
                   ! Normal matrix
                   geometricFactor = DHat(itheta,izeta)*temp/(BHat(itheta,izeta)**3)
                   factor = (alpha*Delta/four) * dPhiHatdpsiHat * geometricFactor
                else
                   select case(whichLambda)
                   case (0) ! Er
                      geometricFactor = DHat(itheta,izeta)*temp/(BHat(itheta,izeta)**3)
                      factor = (alpha*Delta/four)*geometricFactor
                   case (1) ! BHat
                      angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                      cos_angle = cos(angle)
                      sin_angle = sin(angle)                   
                      dBHatdthetadLambda = -m*sin_angle
                      dBHatdzetadLambda = n*Nperiods*sin_angle
                      dBHatdLambda = cos_angle
                      dTempdLambda = BHat_sub_zeta(itheta,izeta) * dBHatdthetadLambda &
                           - BHat_sub_theta(itheta,izeta) * dBHatdzetadLambda
                      geometricFactor = dTempdLambda/ (BHat(itheta,izeta)*(GHat+iota*IHat)) &
                           - temp/(BHat(itheta,izeta)*BHat(itheta,izeta)*(GHat+iota*IHat))
                      factor = (alpha*Delta*dPhiHatdpsiHat/four)*geometricFactor
                   case (2) ! BHat_sub_theta = IHat
                      geometricFactor = -dBHatdzeta(itheta,izeta)/(BHat(itheta,izeta)*(GHat+iota*IHat)) &
                           - temp*iota/(BHat(itheta,izeta)*(GHat+iota*IHat)**2)
                      factor = (alpha*Delta*dPhiHatdpsiHat/four)*geometricFactor
                   case (3) ! BHat_sub_zeta = GHat
                      geometricFactor = dBHatdtheta(itheta,izeta)/(BHat(itheta,izeta)*(GHat+iota*IHat)) &
                           - temp/(BHat(itheta,izeta)*(GHat+iota*IHat)**2)
                      factor = (alpha*Delta*dPhiHatdpsiHat/four)*geometricFactor
                   case (4) ! BHat_sup_theta = iota
                      geometricFactor = -IHat*temp/(BHat(itheta,izeta)*(GHat+iota*IHat)**2)
                      factor = (alpha*Delta*dPhiHatdpsiHat/four)*geometricFactor
                   end select

                end if
                
                     
                do ix=ixMin,Nx
                   do L=0,(Nxi_for_x(ix)-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)

                      ! Diagonal-in-L term
                      call MatSetValueSparse(matrix, rowIndex, rowIndex, &
                           (L+1)*L/((2*L-one)*(2*L+three))*factor, ADD_VALUES, ierr)

                      ! Drop the off-by-2 diagonal terms in L if this is the preconditioner
                      ! and preconditioner_xi = 1:
                      if (((whichMatrix .ne. 0) .and. (whichMatrix .ne. 5)) .or. preconditioner_xi==0) then

                         if (L<Nxi_for_x(ix)-2) then
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
    end do
  end subroutine nonstandardErDdxi

  ! ****************************************************************
  ! Add the non-standard d/dxi term associated with magnetic drifts:
  ! ****************************************************************
  subroutine nonstandardMagneticDriftDdxi(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor, temp
    integer :: ispecies, ix, itheta, izeta, L, ell, rowIndex, colIndex

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)
     
       if ((magneticDriftScheme>0) .and. (whichMatrix .ne. 2)) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax
                select case (magneticDriftScheme)
                case (1,2,8,9)
                   temp = (dBHat_sub_psi_dzeta(itheta,izeta) - dBHat_sub_zeta_dpsiHat(itheta,izeta)) * dBHatdtheta(itheta,izeta) &
                        + (dBHat_sub_theta_dpsiHat(itheta,izeta) - dBHat_sub_psi_dtheta(itheta,izeta)) * dBHatdzeta(itheta,izeta)
                   
                   if (.not. force0RadialCurrentInEquilibrium) then
                      temp = temp + (dBHat_sub_zeta_dtheta(itheta,izeta) - dBHat_sub_theta_dzeta(itheta,izeta)) * dBHatdpsiHat(itheta,izeta)
                   end if
                case (3,4,7)
                   temp = 0
                case (5,6)
                   temp = - (  BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) &
                             - BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta)) &
                          / BHat(itheta,izeta) * gradpsidotgradB_overgpsipsi(itheta,izeta)
                case default
                   stop "Invalid magneticDriftScheme in d/dxi term"
                end select

                do ix=ixMin,Nx
                   factor = -Delta*DHat(itheta,izeta)*THat*x(ix)*x(ix)/(2*Z*(BHat(itheta,izeta)**3)) * temp
                     
                   do L=0,(Nxi_for_x(ix)-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)

                      ! Diagonal-in-L term
                      call MatSetValueSparse(matrix, rowIndex, rowIndex, &
                           (L+1)*L/((2*L-one)*(2*L+three))*factor, ADD_VALUES, ierr)

                      ! Drop the off-by-2 diagonal terms in L if this is the preconditioner
                      ! and preconditioner_xi = 1:
                      if (whichMatrix .ne. 0 .or. preconditioner_xi==0) then

                         if (L<Nxi_for_x(ix)-2) then
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
    end do
  end subroutine nonstandardMagneticDriftDdxi
  
  ! *********************************************************
  ! Add the collisionless d/dx term associated with E_r
  ! *********************************************************
  ! This term always operates on f0, so it should always be included when whichMatrix=2 even if includeXDotTerm=.false.!
  ! This term also operates on f1 if includeXDotTerm=.t.
  subroutine collisionlessErDdx(matrix, whichMatrix, whichLambda)
    Mat :: matrix
    integer, intent(in) :: whichMatrix, whichLambda
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor, stuffToAdd, xDotFactor, xDotFactor2, geometricFactor
    integer :: ispecies, ix, itheta, izeta, L, ell, colIndex, ix_row, ix_col, ixMinCol
    PetscScalar, dimension(:,:), allocatable :: xPartOfXDot_plus, xPartOfXDot_minus, xPartOfXDot
    PetscScalar, dimension(:,:), allocatable :: ddxToUse_plus, ddxToUse_minus
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar :: angle, cos_angle, sin_angle, dBHatdLambda, dBHatdthetadLambda, dBHatdzetadLambda

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)

       if (includeXDotTerm .and. (whichMatrix .ne. 2)) then
          allocate(ddxToUse_plus(Nx,Nx))
          allocate(ddxToUse_minus(Nx,Nx))
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
             if ((whichMatrix==0 .or. whichMatrix==5) .and. L >= preconditioner_x_min_L) then
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

             factor = -adjointFactor*alpha*Delta*dPhiHatdpsiHat/4
             do itheta=ithetaMin,ithetaMax
                do izeta=izetaMin,izetaMax

                   ! xDotFactor is only used to calculate
                   ! which upwinding matrix to use
                   xDotFactor = factor*DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                        * (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) &
                        - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))

                   if (force0RadialCurrentInEquilibrium) then
                      xDotFactor2 = zero
                   else
                      ! This should not happen
                      if (whichMatrix == 6) then
                         if (masterProc) then
                            print *, "force0RadialCurrentInEquilibrium must be true for adjoint diagnostics"
                         end if
                      end if
                      xDotFactor2 = factor*DHat(itheta,izeta)/(BHat(itheta,izeta)**2) * 2 &
                           * (dBHat_sub_zeta_dtheta(itheta,izeta) - dBHat_sub_theta_dzeta(itheta,izeta))
                   end if

                   if (xDotFactor>0) then  
                      xPartOfXDot = xPartOfXDot_plus
                   else
                      xPartOfXDot = xPartOfXDot_minus
                   end if

                   if (whichMatrix < 6) then
                      ! Normal matrix
                      geometricFactor = DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                      * (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) &
                      - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))
                      factor = -adjointFactor*alpha*Delta*dPhiHatdpsiHat/4
                   else
                      select case(whichLambda)
                      case (0) ! Er
                         geometricFactor = DHat(itheta,izeta)/(BHat(itheta,izeta)**3) &
                              * (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) &
                              - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))
                         factor = -alpha*Delta/four
                      case (1) ! BHat
                         angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                         cos_angle = cos(angle)
                         sin_angle = sin(angle)                      
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

                   end if
                   
                      

                   rowIndices = -1
                   do ix=min_x_for_L(L),Nx
                      rowIndices(ix)=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                   end do

                   ! Term that is diagonal in L:
                   colIndices = rowIndices
                   stuffToAdd = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*factor*geometricFactor &
                        + (2*L*L+2*L-one)/((two*L+3)*(2*L-1))*xDotFactor2
                   do ix_col=max(ixMinCol,min_x_for_L(L)),Nx
                      do ix_row=max(ixMin,min_x_for_L(L)),Nx
                         call MatSetValueSparse(matrix, rowIndices(ix_row), colIndices(ix_col), &
                              stuffToAdd*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                      end do
                   end do

                   ! Drop the off-by-2 diagonal terms in L if this is the preconditioner
                   ! and preconditioner_xi = 1:
                   if ((whichMatrix>0 .and. whichMatrix .ne. 5) .or. preconditioner_xi==0) then

                      ! Term that is super-super-diagonal in L:
                      if (L<(Nxi-2)) then
                         ell = L + 2
                         stuffToAdd = (L+1)*(L+2)/((two*L+5)*(2*L+3))*(factor*geometricFactor+xDotFactor2)
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
                         stuffToAdd = L*(L-1)/((two*L-3)*(2*L-1))*(factor*geometricFactor+xDotFactor2)
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
                   end if

                end do
             end do
          end do
          deallocate(rowIndices)
          deallocate(colIndices)
          deallocate(xPartOfXDot)
          deallocate(xPartOfXDot_plus)
          deallocate(xPartOfXDot_minus)
          deallocate(ddxToUse_plus)
          deallocate(ddxToUse_minus)
       end if
    end do
  end subroutine collisionlessErDdx

  ! *********************************************************
  ! Add the adjoint operator term associated with E_r and
  ! full trajectories
  ! *********************************************************
  subroutine adjointErFullTrajectories(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat, factor
    integer :: ispecies, ix, itheta, izeta, L, ell, colIndex, rowIndex
    
    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)
       
       if (RHSMode>3 .and. (discreteAdjointOption .eqv. .false.) .and. (includeXDotTerm .eqv. .true.) .and. &
          (useDKESExBDrift .eqv. .false.) .and. (includeElectricFieldTermInXiDot .eqv. .true.) .and. ((whichMatrix == 4) .or. (whichMatrix == 5))) then
          do itheta=ithetaMin,ithetaMax
             do izeta=izetaMin,izetaMax
                factor = (alpha*delta/2)*(BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta) &
                    *dBHatdtheta(itheta,izeta))*dPhiHatdPsiHat*DHat(itheta,izeta)/BHat(itheta,izeta)**3
                do ix=ixMin,Nx
                   do L=0,(Nxi_for_x(ix)-1)
                      rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,BLOCK_F)
                      if (L<Nxi_for_x(ix)-2) then
                         ! Super-super-diagonal-in-L term:
                         ell = L+2
                         colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              (L+2)*(L+1)/((2*L+five)*(2*L+three))*x2(ix)*factor, ADD_VALUES, ierr)
                      end if
                      if (L>1) then
                         ! Sub-sub-diagonal-in-L term:
                         ell = L-2
                         colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              (L-one)*L/((two*L-three)*(two*L-one))*x2(ix)*factor, ADD_VALUES, ierr)
                      end if
                      ! diagonal-in-L term
                      ell = L
                      colIndex=getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                        two*(three*L*L+three*L-two)/((two*L+three)*(two*L-one))*x2(ix)*factor, ADD_VALUES, ierr)
                   end do
                end do
             end do
          end do
       end if
    end do
  end subroutine adjointErFullTrajectories

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
  subroutine radialExBDrive(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat
    integer :: ispecies, ix, itheta, izeta, L, colIndex, rowIndex, ithetaRow, ithetaCol, izetaRow, izetaCol
    PetscScalar :: dPhiHatdpsiHatToUseInRHS, factor1, factor2, factorJ1, factorJ2, factorJ3, factorJ4, factorJ5  !!Added by AM 2016-03
    integer :: rowIndex2, L2 !!Added by AM 2016-03
    PetscScalar :: xPartOfRHS, xPartOfRHS2 !!Added by AM 2016-03

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)

       if ((whichMatrix .ne. 2) .and. includePhi1 .and. includePhi1InKineticEquation .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
          L=0
          L2=2 !!Added by AM 2016-03 to add extra P_2 terms

          !!Added by AM 2016-03!!
          !!The naming "RHS" here is a remnant from that term was first introduced in evaluateResidual.F90 in an earlier version of SFINCS
          if (RHSMode==1) then
             dPhiHatdpsiHatToUseInRHS = dPhiHatdpsiHat
          else
             dPhiHatdpsiHatToUseInRHS = 0
          end if
          !!!!!!!!!!!!!!!!!!!!!!!

          do ix = ixMin,Nx

             !!Added by AM 2016-03!!
             !!The naming "RHS" here is a remnant from that term was first introduced in evaluateResidual.F90 in an earlier version of SFINCS
             xPartOfRHS = x2(ix)*expx2(ix)*( dnHatdpsiHats(ispecies)/nHat &
                  + alpha*Z/THat*dPhiHatdpsiHatToUseInRHS &
                  + (x2(ix) - three/two)*dTHatdpsiHats(ispecies)/THat)

             xPartOfRHS2 = x2(ix)*expx2(ix)*dTHatdpsiHats(ispecies)/(THat*THat)
             !!!!!!!!!!!!!!!!!!!!!!!

             ! d Phi_1 / d theta term:
             itheta = -1 ! Just so a runtime error occurs if I use itheta by mistake
             do ithetaRow = ithetaMin, ithetaMax
                do izeta = izetaMin, izetaMax
                   rowIndex = getIndex(ispecies, ix, L+1, ithetaRow, izeta, BLOCK_F)

                   !!Added by AM 2016-03!!
                   factor1 = -exp(-Z*alpha*Phi1Hat(ithetaRow,izeta)/THat)*alpha*Delta*DHat(ithetaRow,izeta)*BHat_sub_zeta(ithetaRow,izeta) &
                        /(two*BHat(ithetaRow,izeta)*BHat(ithetaRow,izeta)) &
                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)

                   factor2 = - alpha*alpha*Delta*DHat(ithetaRow,izeta)/(two*pi*sqrtpi*BHat(ithetaRow,izeta)*BHat(ithetaRow,izeta)) &
                        * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
                        * exp(-Z*alpha*Phi1Hat(ithetaRow,izeta)/THat)*expx2(ix)*BHat_sub_zeta(ithetaRow,izeta) &
                        * (dPhiHatdpsiHat + Phi1Hat(ithetaRow,izeta)*dTHatdpsiHats(ispecies)/THat)

                   do ithetaCol = 1,Ntheta
                      colIndex = getIndex(1,1,1,ithetaCol, izeta, BLOCK_QN)

                      call MatSetValueSparse(matrix, rowIndex, colIndex,& !!Added by AM 2016-03
                           (factor1 + factor2)*ddtheta(ithetaRow, ithetaCol), ADD_VALUES, ierr) !!Added by AM 2016-03
                     !!!!!!!!!!!!!!!!!!!!!!!!!

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

                   !! ADDED BY AM 2016-02!!
                   factor1 = exp(-Z*alpha*Phi1Hat(itheta,izetaRow)/THat)*alpha*Delta*DHat(itheta,izetaRow)*BHat_sub_theta(itheta,izetaRow) &
                        /(two*BHat(itheta,izetaRow)*BHat(itheta,izetaRow)) &
                        * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
                        * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)

                   factor2 = alpha*alpha*Delta*DHat(itheta,izetaRow)/(two*pi*sqrtpi*BHat(itheta,izetaRow)*BHat(itheta,izetaRow)) &
                        * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
                        * exp(-Z*alpha*Phi1Hat(itheta,izetaRow)/THat)*expx2(ix)*BHat_sub_theta(itheta,izetaRow) &
                        * (dPhiHatdpsiHat + Phi1Hat(itheta,izetaRow)*dTHatdpsiHats(ispecies)/THat)

                   do izetaCol = 1,Nzeta
                      colIndex = getIndex(1,1,1,itheta, izetaCol, BLOCK_QN)

                      call MatSetValueSparse(matrix, rowIndex, colIndex,& !!Added by AM 2016-03
                           (factor1 + factor2)*ddzeta(izetaRow, izetaCol), ADD_VALUES, ierr) !!Added by AM 2016-03
                   !!!!!!!!!!!!!!!!!!!!!!!!!

                   end do
                end do
             end do
             !!!!!!!!!!!!!!!!!!!!!!!!
             !!Added by AM 2016-02!!
             !!Add additional terms in the Jacobian
             if (whichMatrix == 0 .or. whichMatrix == 1) then
                itheta = -1 ! Just so a runtime error occurs if I use itheta by mistake
                izeta = -1 ! Just so a runtime error occurs if I use itheta by mistake
                ithetaRow = -1 ! Just so a runtime error occurs if I use itheta by mistake
                ithetaCol = -1 ! Just so a runtime error occurs if I use itheta by mistake
                izetaRow = -1 ! Just so a runtime error occurs if I use itheta by mistake
                izetaCol = -1 ! Just so a runtime error occurs if I use itheta by mistake
                do itheta = ithetaMin, ithetaMax
                   do izeta = izetaMin, izetaMax
                      rowIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                      rowIndex2 = getIndex(ispecies, ix, L2+1, itheta, izeta, BLOCK_F) !!Added by AM 2016-03 to add extra P_2 terms

                      !!These terms are diagonal in theta and zeta
                      colIndex = getIndex(1,1,1,itheta, izeta, BLOCK_QN)

                      !!The following factors are only used in the Jacobian. These terms do not contain d/dtheta or d/dzeta
                      !!but d(Phi1)/dtheta, d(Phi1)/dzeta, dB/dtheta, dB/dzeta.
                      !!Therefore they are the same when adding the d/dtheta terms as when adding the d/dzeta terms.

                      !!factorJ1 stems from factor1 but contains both dPhi1Hatdtheta and dPhi1Hatdzeta parts
                      factorJ1 = exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)*alpha*Delta*DHat(itheta,izeta) &
                           *(- BHat_sub_zeta(itheta,izeta)* dPhi1Hatdtheta(itheta, izeta) + BHat_sub_theta(itheta, izeta)*dPhi1Hatdzeta(itheta, izeta))&
                           /(two*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                           * nHat * (mHat*sqrtMHat)/(THat*sqrtTHat*pi*sqrtpi) * expx2(ix) & ! This line is f_M
                           * (dNHatdpsiHats(ispecies)/nHat + (x2(ix)-3/two)*dTHatdpsiHats(ispecies)/THat)  

                      !!factorJ2 stems from factor2 but contains both dPhi1Hatdtheta and dPhi1Hatdzeta parts
                      factorJ2 = alpha*alpha*Delta*DHat(itheta,izeta)/(two*pi*sqrtpi*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                           * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
                           * exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)*expx2(ix) &
                           * (- BHat_sub_zeta(itheta,izeta)*dPhi1Hatdtheta(itheta, izeta) + BHat_sub_theta(itheta, izeta)*dPhi1Hatdzeta(itheta, izeta)) &
                           * (dPhiHatdpsiHat + Phi1Hat(itheta,izeta)*dTHatdpsiHats(ispecies)/THat)

                      !!factorJ3 is only used in the Jacobian, it corresponds to a part \propto exp(-Ze Phi1/T) which is
                      !!implemented in evaluateResidual.F90 for the residual
                      factorJ3 = Delta*nHat*mHat*sqrtMHat &
                           /(2*pi*sqrtpi*Z*(BHat(itheta,izeta)**3)*sqrtTHat) &
                           *(- BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) + BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                           * DHat(itheta,izeta) * (xPartOfRHS + xPartOfRHS2*Z*alpha*Phi1Hat(itheta,izeta))  & 
                           * exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)

                      !!factorJ4 is only used in the Jacobian, it corresponds to the second term in factorJ2 divided by Phi1
                      factorJ4 =  alpha*alpha*Delta*DHat(itheta,izeta)/(two*pi*sqrtpi*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                           * Z * nHat*mHat*sqrtMHat/(THat*THat*sqrtTHat) &
                           * exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)*expx2(ix) &
                           *(- BHat_sub_zeta(itheta,izeta)*dPhi1Hatdtheta(itheta,izeta) + BHat_sub_theta(itheta,izeta)*dPhi1Hatdzeta(itheta,izeta))&
                           * dTHatdpsiHats(ispecies)/THat

                      !!factorJ5 is only used in the Jacobian, it corresponds to the second term in factorJ3 divided by Phi1
                      factorJ5 = Delta*nHat*mHat*sqrtMHat &
                           /(2*pi*sqrtpi*(BHat(itheta,izeta)**3)*sqrtTHat) &
                           *(- BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta) + BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta))&
                           * DHat(itheta,izeta) * xPartOfRHS2*alpha  & 
                           * exp(-Z*alpha*Phi1Hat(itheta,izeta)/THat)


                      !!SECTION TEMPORARILY COMMENTED BY AM!!
                      !!Add L=0 component
                      call MatSetValue(matrix, rowIndex, colIndex,&
                           -Z*alpha*(factorJ1 + factorJ2)/THat &
                           -Z*alpha* (4/three)*factorJ3/THat &
                           + factorJ4 &
                           + (4/three)*factorJ5, ADD_VALUES, ierr) 

                      !!Add L=2 component
                      call MatSetValue(matrix, rowIndex2, colIndex,&
                           -Z*alpha* (two/three)*factorJ3/THat &
                           + (two/three)*factorJ5, ADD_VALUES, ierr)
                   end do
                end do
             end if
          end do
       end if
    end do
  end subroutine radialExBDrive



  !NOTE BY AM 2018-12: This term should be added even if readExternalPhi1 = .true.
  subroutine parallelE(matrix, whichMatrix, whichLambda)
    Mat :: matrix
    integer, intent(in) :: whichMatrix, whichLambda
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat
    integer :: ispecies, ix, itheta, izeta, L, ell, colIndex, rowIndex, ixMinCol, ix_row, ix_col
    PetscScalar :: factor, geometricFactor, geometricFactor2
    PetscScalar, dimension(:,:), allocatable :: nonlinearTerm_Lp1, nonlinearTerm_Lm1
    PetscScalar :: angle, cos_angle, dBHatdLambda

    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)
         
       if (includePhi1 .and. includePhi1InKineticEquation .and. (whichMatrix == 1 .or. whichMatrix == 3 .or. whichMatrix == 4 .or. (whichMatrix==0 .and. ((.not. reusePreconditioner) .or. readExternalPhi1)) .or. (whichMatrix==5 .and. ((.not. reusePreconditioner) .or. readExternalPhi1)) .or. whichMatrix == 6) .and. localUsePhi1) then 
          allocate(nonlinearTerm_Lp1(Nx,Nx))
          allocate(nonlinearTerm_Lm1(Nx,Nx))
          nonlinearTerm_Lp1 = zero
          nonlinearTerm_Lm1 = zero
          do L=0,(Nxi-1)
             if (L>0 .and. pointAtX0) then
                ixMinCol = 2
             else
                ixMinCol = 1
             end if

             if ((whichMatrix==0 .and. L >= preconditioner_x_min_L) .or. (whichMatrix==5 .and. L >= preconditioner_x_min_L)) then
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
                   ! geometricFactor = (BHat^theta * dPhi1Hatdtheta(itheta,izeta) 
                   !    + BHat^zeta * dPhi1Hatdzeta(itheta,izeta))/BHat
                   ! = BHat/(GHat+iota * IHat) * (iota * dPhi1Hatdtheta
                   !    + dPhi1Hatdzeta)

                   if (whichMatrix < 6) then
                      ! normal matrix
                      geometricFactor = (BHat_sup_theta(itheta,izeta)*dPhi1Hatdtheta(itheta,izeta) &
                           + BHat_sup_zeta(itheta,izeta)*dPhi1Hatdzeta(itheta,izeta))/BHat(itheta,izeta)
                      geometricFactor2 = BHat(itheta,izeta)*(iota * dPhi1Hatdtheta(itheta,izeta) + dPhi1Hatdzeta(itheta,izeta))/(GHat + iota * IHat) 
                   else
                      ! dMatrix/dlambda
                      select case(whichLambda)
                      case (0) ! Er
                         geometricFactor = 0
                      case (1) ! BHat cos
                         angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                         cos_angle = cos(angle)                      
                         dBHatdLambda = cos_angle
                         ! dBHatdzetadLambda = n*Nperiods*sin_angle
                         ! dBHatdthetadLambda = -m*sin_angle
                         geometricFactor = (iota*dPhi1Hatdtheta(itheta,izeta) &
                              + dPhi1Hatdzeta(itheta,izeta))*dBHatdLambda/(GHat + iota * IHat)
                      case (2) ! IHat
                         geometricFactor = -(iota * dPhi1Hatdtheta(itheta,izeta) + dPhi1Hatdzeta(itheta,izeta)) * BHat(itheta,izeta) * iota/((GHat + iota * IHat)**2)
                      case (3) ! GHat
                         geometricFactor = -(iota * dPhi1Hatdtheta(itheta,izeta) + dPhi1Hatdzeta(itheta,izeta)) * BHat(itheta,izeta)/((GHat + iota * IHat)**2)
                      case (4) ! iota
                         geometricFactor = dPhi1Hatdtheta(itheta,izeta) * BHat(itheta,izeta)/(GHat + iota * IHat) &
                              - (iota * dPhi1Hatdtheta(itheta,izeta) + dPhi1Hatdzeta(itheta,izeta)) * BHat(itheta,izeta) * IHat/((GHat + iota * IHat)**2)
                      end select
                   end if
                   
                   factor = -adjointFactor * alpha * Z * geometricFactor/(2 * sqrtTHat * sqrtMHat)
                   ! what is this factor?
                   !factor = -1.5*adjointFactor * alpha * Z * geometricFactor/(2 * sqrtTHat) 
                       

                   !do ix_row=ixMin,Nx
                   do ix_row=max(ixMin,min_x_for_L(L)),Nx
                      rowIndex = getIndex(ispecies,ix_row,L+1,itheta,izeta,BLOCK_F)

                      ! Term that is super-diagonal in L:
                      if (L<Nxi-1) then
                         ell = L + 1
                         !do ix_col=ixMinCol,Nx
                         do ix_col=max(ixMinCol,min_x_for_L(ell)),Nx
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
                         !do ix_col=ixMinCol,Nx
                         do ix_col=max(ixMinCol,min_x_for_L(ell)),Nx
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
    end do
  end subroutine parallelE

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
  
  !NOTE BY AM 2018-12: This term should not be added if readExternalPhi1 = .true.
  subroutine nonlinearParallelE(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    PetscScalar :: Z, nHat, THat, mHat, sqrtTHat, sqrtMHat
    integer :: ispecies, ix, itheta, izeta, L, ell, colIndex, rowIndex, index, j
    PetscScalar :: factor
    PetscScalar, dimension(:), allocatable :: tempVector1, tempVector2


    do ispecies = 1,Nspecies
       nHat = nHats(ispecies)
       THat = THats(ispecies)
       mHat = mHats(ispecies)
       Z = Zs(ispecies)
       sqrtTHat = sqrt(THat)
       sqrtMHat = sqrt(mHat)
         
       if (includePhi1 .and. includePhi1InKineticEquation .and. (.not. readExternalPhi1) .and. (whichMatrix==1 .or. (whichMatrix==0 .and. .not. reusePreconditioner))) then !!Added by AM 2018-12
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
                      tempVector1=0
                      tempVector2=0
                      !do ix=1,Nx
                      do ix=min_x_for_L(ell),Nx
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
                      tempVector1=0
                      tempVector2=0
                      !do ix=1,Nx
                      do ix=min_x_for_L(ell),Nx
                         index = getIndex(ispecies,ix,ell+1,itheta,izeta,BLOCK_F)
                         ! Add 1 because we are indexing a Fortran array instead of a PETSc object
                         tempVector1(ix) = stateArray(index+1)
                         tempVector2(ix) = tempVector2(ix) + (L+1)*(L+2)/(two*L+3)*stateArray(index+1)/x(ix)
                      end do
                      tempVector2 = tempVector2 + (L+1)/(two*L+3) * matmul(ddx,tempVector1)
                   end if

                   !do ix=ixMin,Nx
                   do ix=max(ixMin,min_x_for_L(L)),Nx
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
    end do
  end subroutine nonlinearParallelE


  ! *********************************************************
  !
  ! Next, we add the collision operator.
  !
  ! *********************************************************
  ! The collision operator always acts on f1.
  ! The collision operator also acts on f0 if includeTemperatureEquilibrationTerm=.t.
  subroutine addCollisionOperator(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    integer :: iSpeciesA, iSpeciesB, ix, itheta, izeta, L, colIndex, rowIndex, index, j, i, ix_col, ix_row, ixMinCol
    PetscScalar :: temp, temp1, temp2, speciesFactor, speciesFactor2, T32m
    PetscScalar, dimension(:), allocatable :: xb, expxb2
    PetscScalar, dimension(:,:), allocatable :: tempMatrix, tempMatrix2, extrapMatrix
    PetscScalar, dimension(:), allocatable :: CHatTimesf, fM, f1b ! Added by AI (2017-09) 
    PetscScalar, dimension(:,:), allocatable :: M21, M32, M22, M33 
    PetscScalar, dimension(:,:), allocatable :: LaplacianTimesX2WithoutL
    PetscScalar, dimension(:,:), allocatable :: M22BackslashM21, M33BackslashM32
    PetscScalar, dimension(:,:,:), allocatable :: M22BackslashM21s, M33BackslashM32s
    PetscScalar, dimension(:), allocatable :: erfs, Psi_Chandra
    PetscScalar, dimension(:,:), allocatable :: tempExtrapMatrix, fToFInterpolationMatrix_plus1
    PetscScalar, dimension(:,:), allocatable :: fToFInterpolationMatrix
    PetscScalar, dimension(:,:), allocatable :: M11, nuDHat, M11J !! M11J Added by AI (2017-09)
    PetscScalar, dimension(:,:,:,:), allocatable :: nuDHatpol,nuDHatpolJ ! Added by AI (2017-09)
    PetscScalar, dimension(:,:), allocatable :: CHat, M12, M13, CHatJ !! CHatJ, Added by AI (2017-09)
    PetscScalar, dimension(:,:,:,:), allocatable :: CECD
    PetscScalar, dimension(:,:,:,:,:,:), allocatable :: CECDpol, CECDpolJ !! Added by AI (2017-09)
    integer, dimension(:), allocatable :: IPIV  ! Needed by LAPACK
    PetscScalar, dimension(:,:), allocatable :: potentialsToFInterpolationMatrix
    integer :: LAPACKInfo
    PetscScalar :: preFactor, preFactorJ !! preFactor, preFactorJ Added by AI (2017-09) 


    ! Needed for Pitch-angle scattering operator only
    PetscScalar :: CHat_element, CHat_elementJ


    if (whichMatrix .ne. 2 .or. includeTemperatureEquilibrationTerm) then
       allocate(xb(Nx))
       allocate(expxb2(Nx))
       allocate(erfs(Nx))
       allocate(Psi_Chandra(Nx))
       allocate(nuDHat(Nspecies, Nx))


       
       select case (collisionOperator)
          
       case (0)
          ! *********************************************************
          ! Full linearized Fokker-Planck operator
          ! *********************************************************
          
          
          ! *********************************************************
          ! In preparation for adding the collision operator,
          ! create several matrices which will be needed.
          ! *********************************************************
          allocate(fToFInterpolationMatrix(Nx,Nx))
          allocate(CHat(Nx,Nx))
          allocate(CECD(Nspecies, Nspecies, Nx, Nx))
          allocate(IPIV(NxPotentials))
          allocate(potentialsToFInterpolationMatrix(Nx, NxPotentials))

          
          allocate(M11(Nx,Nx))
          allocate(M12(Nx,NxPotentials))
          allocate(M13(Nx,NxPotentials))
          allocate(M21(NxPotentials, Nx))
          allocate(M22(NxPotentials,NxPotentials))
          allocate(M32(NxPotentials, NxPotentials))
          allocate(M33(NxPotentials,NxPotentials))

          allocate(M22BackslashM21(NxPotentials, Nx))
          allocate(M33BackslashM32(NxPotentials, NxPotentials))
          allocate(M22BackslashM21s(NL,NxPotentials, Nx))
          allocate(M33BackslashM32s(NL,NxPotentials, NxPotentials))

          allocate(LaplacianTimesX2WithoutL(NxPotentials, NxPotentials))
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
          
          if (includePhi1InCollisionOperator .and. includePhi1 .and. includePhi1InKineticEquation) then   !! Added by AI (2017-09), Do with Phi1
             ! if this is false, proceed to original code     
             
             ! ************************************************************************************
             ! This section has been added by AI (2017-09) in order to include Phi1
             ! in the collision operator. Phi1 is included by setting 
             ! includePhi1InCollisionOperator = .true. Note that this section replaces
             ! the original block when includePhi1InCollisionOperator = .true.
             ! See the documentation at
             ! https://github.com/landreman/sfincs/blob/poloidalVariationInCollisionOperator/doc/PoloidalVariationInCollisionOperator_code.pdf
             ! ************************************************************************************
             
             ! Allocate matrices
             allocate(CHatTimesf(Nx)) 
             allocate(fM(Nx)) 
             allocate(f1b(Nx))
             !! Added by AI (2017-09) Required for Phi1 in collision operator
             allocate(CHatJ(Nx,Nx))
             allocate(M11J(Nx,Nx))
             allocate(nuDHatpol(Nspecies, Nx,Ntheta,Nzeta))
             allocate(nuDHatpolJ(Nspecies, Nx,Ntheta,Nzeta))
             allocate(CECDpol(Nspecies, Nspecies, Nx, Nx,Ntheta,Nzeta))
             allocate(CECDpolJ(Nspecies, Nspecies, Nx, Nx,Ntheta,Nzeta)) 
             
             ! Initiate matrices (same as in the original, except that these terms are now
             ! theta, zeta dependent)                             
             nuDHatpol = zero
             CECDpol = zero 
             CECDpolJ = zero
             nuDHatpolJ = zero
             CHatTimesf = zero !!Added by AM 2018-01
             fM = zero !!Added by AM 2018-01
             f1b = zero !!Added by AM 2018-01
             
             
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
                   
                   ! In order to get the correct preFactor and introduce the Phi1Hat dependence
                   ! we need to make a loop over itheta, izeta                                                                
                   do itheta=ithetaMin,ithetaMax
                      do izeta=izetaMin,izetaMax
                         
                         ! Generate preFactor for nHats(iSpeciesA) terms
                         preFactor =  exp(-Zs(iSpeciesA)*alpha*Phi1Hat(itheta,izeta)/Thats(iSpeciesA))
                         
                         speciesFactor = 3 * nHats(iSpeciesA)*preFactor* mHats(iSpeciesA)/mHats(iSpeciesB) &
                         * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) / T32m
                         
                         ! Using the resulting interpolation matrix,
                         ! add CD (the part of the field term independent of Rosenbluth potentials.
                         ! CD is dense in the species indices.
                         
                         do ix=1,Nx
                            CECDpol(iSpeciesA, iSpeciesB, ix, :,itheta,izeta) = CECDpol(iSpeciesA, iSpeciesB, ix, :,itheta,izeta) &
                            + speciesFactor * expx2(ix) * fToFInterpolationMatrix(ix, :)
                         end do ! Should be the same if we put outside loop and use (iSpeciesA, iSpeciesB, ix, :,:,:)
                         
                         ! Generate preFactor for nHats(iSpeciesB) terms
                         preFactor =  exp(-Zs(iSpeciesB)*alpha*Phi1Hat(itheta,izeta)/Thats(iSpeciesB))
                         
                         ! Build the pitch-angle scattering frequency:
                         nuDHatpol(iSpeciesA, :,itheta,izeta) =  nuDHatpol(iSpeciesA, :,itheta,izeta) &
                         + (three*sqrtpi/four) / T32m &
                         * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) &
                         * nHats(iSpeciesB)*preFactor*(erfs - Psi_Chandra)/(x*x*x)
                         
                         speciesFactor = 3*sqrtpi/four * nHats(iSpeciesB)* preFactor  &
                         * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) / T32m
                         
                         do ix=1,Nx
                            !Now add the d2dx2 and ddx terms in CE:
                            !CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                            CECDpol(iSpeciesA, iSpeciesA, ix, :,itheta,izeta) = CECDpol(iSpeciesA, iSpeciesA, ix, :,itheta,izeta) &
                            + speciesFactor * (Psi_Chandra(ix)/x(ix)*d2dx2ToUse(ix,:) &
                            + (-2*THats(iSpeciesA)*mHats(iSpeciesB)/(THats(iSpeciesB)*mHats(iSpeciesA)) &
                            * Psi_Chandra(ix)*(1-mHats(iSpeciesA)/mHats(iSpeciesB)) &
                            + (erfs(ix)-Psi_Chandra(ix))/x2(ix)) * ddxToUse(ix,:))
                            
                            ! Lastly, add the part of CE for which f is not differentiated:
                            ! CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                            CECDpol(iSpeciesA, iSpeciesA, ix, ix,itheta,izeta) = CECDpol(iSpeciesA, iSpeciesA, ix, ix,itheta,izeta) &
                            + speciesFactor *4/sqrtpi*THats(iSpeciesA)/THats(iSpeciesB) &
                            *sqrt(THats(iSpeciesA)*mHats(iSpeciesB)/(THats(iSpeciesB)*mHats(iSpeciesA))) &
                            * expxb2(ix)
                         end do ! ix
                         
                         ! If we are calculating the Jacobian, we need to redo the same, but using a different preFactor                              
                         if (whichMatrix == 1 .or. whichMatrix == 0) then
                            
                            ! Generate preFactorJ for nHats(iSpeciesB) terms              
                            preFactorJ =  (-Zs(iSpeciesA)*alpha/Thats(iSpeciesA)) &
                            *exp(-Zs(iSpeciesA)*alpha*Phi1Hat(itheta,izeta)/Thats(iSpeciesA))
                            
                            speciesFactor = 3 * nHats(iSpeciesA)*preFactorJ* mHats(iSpeciesA)/mHats(iSpeciesB) &
                            * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) / T32m
                            
                            ! Using the resulting interpolation matrix,
                            ! add CD (the part of the field term independent of Rosenbluth potentials.
                            ! CD is dense in the species indices.
                            
                            do ix=1,Nx
                               CECDpolJ(iSpeciesA, iSpeciesB, ix, :,itheta,izeta) = CECDpolJ(iSpeciesA, iSpeciesB, ix, :,itheta,izeta) &
                               + speciesFactor * expx2(ix) * fToFInterpolationMatrix(ix, :)
                            end do
                            
                            ! Generate preFactorJ for nHats(iSpeciesB) terms
                            preFactorJ =  (-Zs(iSpeciesB)*alpha/Thats(iSpeciesB)) &
                            *exp(-Zs(iSpeciesB)*alpha*Phi1Hat(itheta,izeta)/Thats(iSpeciesB))
                            
                            ! Build the pitch-angle scattering frequency:
                            nuDHatpolJ(iSpeciesA, :,itheta,izeta) =  nuDHatpolJ(iSpeciesA, :,itheta,izeta) &
                            + (three*sqrtpi/four) / T32m &
                            * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) &
                            * nHats(iSpeciesB)*preFactorJ*(erfs - Psi_Chandra)/(x*x*x)
                            
                            !! speciesFactor = 3*sqrtpi/four * nHats(iSpeciesB)* preFactor  &
 !!Commented by AM 2018-01
                            speciesFactor = 3*sqrtpi/four * nHats(iSpeciesB)* preFactorJ  & !!Added by AM 2018-01
                            * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) / T32m
                            
                            do ix=1,Nx
                               !Now add the d2dx2 and ddx terms in CE:
                               !CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                               CECDpolJ(iSpeciesA, iSpeciesA, ix, :,itheta,izeta) = CECDpolJ(iSpeciesA, iSpeciesA, ix, :,itheta,izeta) &
                               + speciesFactor * (Psi_Chandra(ix)/x(ix)*d2dx2ToUse(ix,:) &
                               + (-2*THats(iSpeciesA)*mHats(iSpeciesB)/(THats(iSpeciesB)*mHats(iSpeciesA)) &
                               * Psi_Chandra(ix)*(1-mHats(iSpeciesA)/mHats(iSpeciesB)) &
                               + (erfs(ix)-Psi_Chandra(ix))/x2(ix)) * ddxToUse(ix,:))
                               
                               ! Lastly, add the part of CE for which f is not differentiated:
                               ! CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                               CECDpolJ(iSpeciesA, iSpeciesA, ix, ix,itheta,izeta) = CECDpolJ(iSpeciesA, iSpeciesA, ix, ix,itheta,izeta) &
                               + speciesFactor *4/sqrtpi*THats(iSpeciesA)/THats(iSpeciesB) &
                               *sqrt(THats(iSpeciesA)*mHats(iSpeciesB)/(THats(iSpeciesB)*mHats(iSpeciesA))) &
                               * expxb2(ix)
                               
                            end do      ! ix                           
                         end if ! jacobian part
                      end do ! izeta
                   end do ! itheta
                end do ! SpeciesB
             end do ! Species A
             
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
                         
                         ! Because of the new theta,zeta dependence, we need to iterate over itheta and izeta
                         ! already here
                         do itheta=ithetaMin,ithetaMax 
                            do izeta=izetaMin,izetaMax
                               
                               ! Generate preFactor and preFactorJ for nHats(iSpeciesA) terms
                               preFactor =  exp(-Zs(iSpeciesA)*alpha*Phi1Hat(itheta,izeta)/Thats(iSpeciesA))
                               preFactorJ =  (-Zs(iSpeciesA)*alpha/Thats(iSpeciesA)) &
                               *exp(-Zs(iSpeciesA)*alpha*Phi1Hat(itheta,izeta)/Thats(iSpeciesA))
                               
                               speciesFactor = sqrt(THats(iSpeciesA)*mHats(iSpeciesB) &
                               / (THats(iSpeciesB) * mHats(iSpeciesA)))
                               xb =  x * speciesFactor
                               
                               ! Build M11
                               M11 = CECDpol(iSpeciesA, iSpeciesB,:,:,itheta,izeta)
                               M11J = CECDpolJ(iSpeciesA, iSpeciesB,:,:,itheta,izeta)                     
                               if (iSpeciesA == iSpeciesB) then
                                  do i=1,Nx
                                     !!M11(i,i) = M11(i,i) + (-oneHalf*nuDHatpol(iSpeciesA,i,itheta,izeta)*L*(L+1))
 !!Commented by AM 2018-01
                                     M11(i,i) = M11(i,i) + (-oneHalf*nuDHatpol(iSpeciesA,i,itheta,izeta)*L*(L+1)  + Krook*2) !!Added by AM 2018-01
                                     !!M11J(i,i) = M11J(i,i) + (-oneHalf*nuDHatpolJ(iSpeciesA,i,itheta,izeta)*L*(L+1))
 !!Commented by AM 2018-01
                                     M11J(i,i) = M11J(i,i) + (-oneHalf*nuDHatpolJ(iSpeciesA,i,itheta,izeta)*L*(L+1)  + Krook*2) !!Added by AM 2018-01
                                  end do
                               end if
                               
                               !   if (.false.) then
                               if (L < NL) then
                                  ! Add Rosenbluth potential terms.
                                  
                                  if (xGridScheme==5 .or. xGridScheme==6) then
                                     ! New scheme for the Rosenbluth potential terms.
                                     
                                     M11 = M11 + preFactor*RosenbluthPotentialTerms(iSpeciesA,iSpeciesB,L+1,:,:)
                                     M11J = M11J + preFactorJ*RosenbluthPotentialTerms(iSpeciesA,iSpeciesB,L+1,:,:)
                                     
                                     CHat = M11
                                     CHatJ = M11J                           
                                     
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
                                     
                                     CHat = M11 - preFactor*(matmul(M12 - matmul(M13, M33BackslashM32s(L+1,:,:)),&
                                     M22BackslashM21s(L+1,:,:)))
                                     CHatJ = M11J - preFactorJ*(matmul(M12 - matmul(M13, M33BackslashM32s(L+1,:,:)),&
                                     M22BackslashM21s(L+1,:,:)))      
                                  end if
                               else
                                  CHat = M11; ! This is with preFactors for the residual
                                  CHatJ = M11J; ! This is with the preFactors for the Jacobian
                               end if
                               
                               if ((whichMatrix==0 .or. whichMatrix==5).and. L >= preconditioner_x_min_L) then
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
                                              CHatJ(i,j) = zero                                                                                
                                           end if
                                        end do
                                     end do
                                  case (2)
                                     ! Keep only upper-triangular part:
                                     do i=2,Nx
                                        do j=1,(i-1)
                                           CHat(i,j) = zero   
                                           CHatJ(i,j) = zero 
                                        end do
                                     end do
                                  case (3,5)
                                     ! Keep only tridiagonal part:
                                     do i=1,Nx
                                        do j=1,Nx
                                           if (abs(i-j)>1) then
                                              CHat(i,j) = zero   
                                              CHatJ(i,j) = zero                                   
                                           end if
                                        end do
                                     end do
                                  case (4)
                                     ! Keep only the diagonal and super-diagonal:
                                     do i=1,Nx
                                        do j=1,Nx
                                           if (i /= j .and. j /= (i+1)) then
                                              CHat(i,j) = zero   
                                              CHatJ(i,j) = zero                         
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
                               
                               ! Save the residual, and equivalently, the d(collision op.) / d f1 terms
                               
                               do ix_row=max(ixMin,min_x_for_L(L)),Nx 
                                  rowIndex=getIndex(iSpeciesA,ix_row,L+1,itheta,izeta,BLOCK_F)
                                  do ix_col = max(ixMinCol,min_x_for_L(L)),Nx
                                     colIndex=getIndex(iSpeciesB,ix_col,L+1,itheta,izeta,BLOCK_F)
                                     call MatSetValueSparse(matrix, rowIndex, colIndex, &
 
!!                                     call MatSetValue(matrix, rowIndex, colIndex, & 
                                     -nu_n*CHat(ix_row,ix_col), ADD_VALUES, ierr)
                                  end do ! ix_col
                               end do ! ix_row
                               
                               
                               ! The temperature equilibration part is already implemented since f0, which we multiply with in evaluateResidual.F90
                               ! vec, already contains the extra exp(Phi1Hat) factor (see subroutine init_f0)
                               
                               
                               ! ************************************************************************************
                               ! Calculate d(collision op.) / d Phi1 contribution to the Jacobian
                               ! Because of the extra Phi1Hat factors, we need to add more terms when calculating 
                               ! the Jacobian
                               ! ************************************************************************************
                               
                               
                               if (whichMatrix == 1 .or. whichMatrix == 0) then ! Jacobian or Preconditioner
                                  
                                  ! First we generate the distribution function which has to be included
                                  ! in the d(collision op.) / d Phi1 terms
                                  
                                  do ix= max(ixMinCol,min_x_for_L(L)),Nx
                                     ! Generate f1b from state vector
                                     index = getIndex(iSpeciesB,ix,L+1,itheta,izeta,BLOCK_F) 
                                     f1b(ix) = stateArray(index + 1)
                                     
                                     ! If includeTemperatureEquilibrationTerm = .true., also generate the Maxwellian for species B
                                     if (includeTemperatureEquilibrationTerm .and. L==0) then

                                        expxb2 = exp(-xb*xb) !!Added by AM 2018-01  (expxb2 was not defined in this loop)

                                        fM(ix) = sqrt(mhats(iSpeciesB)/Thats(iSpeciesB))*mhats(iSpeciesB)/Thats(iSpeciesB) & 
                                        *nhats(iSpeciesB)/(pi*sqrtpi)*expxb2(ix)
                                     end if
                                  end do
                                  
                                  ! In contrary to the residual, now we need to include the distribution function.
                                  ! Multiply the total collision operator with the distribution function.                       
                                  CHatTimesf = matmul(CHatJ,f1b)
                                  
                                  ! Save into the main matrix
                                  do ix_row=max(ixMin,min_x_for_L(L)),Nx
                                     rowIndex=getIndex(iSpeciesA,ix_row,L+1,itheta,izeta,BLOCK_F)
                                     ! Get column index for the d/dPhi1 terms
                                     colIndex=getIndex(1,1,1,itheta,izeta,BLOCK_QN)
                                     ! Save into the main matrix, note that here we only use ix_row since CHatTimesf is now a vector
                                     call MatSetValue(matrix, rowIndex, colIndex, & 
 
                                     -nu_n*CHatTimesf(ix_row), ADD_VALUES, ierr) 
                                     ! need to use MatSetValue, otherwise petsc gives error
                                  end do ! ix_row
                                  
                                  ! If includeTemperatureEquilibrationTerm = .true., we need to add
                                  ! more terms to the Jacobian 
                                  if (includeTemperatureEquilibrationTerm .and. L == 0) then
                                     
                                     ! Since we now should use f0, which has a Phi1Hat dependence we get
                                     !  d(collision op.*f0) / d Phi1 = CHatJ*f0 + CHat*(d f0 / d Phi1)
                                     
                                     CHatTimesf = matmul((CHatJ + CHat*(-Zs(iSpeciesB)*alpha/Thats(iSpeciesB)))&
                                     *exp(-Zs(iSpeciesB)*alpha*Phi1Hat(itheta,izeta)/Thats(iSpeciesB)),fM)
                                     
                                     ! Save into the main matrix
                                     do ix_row=max(ixMin,min_x_for_L(L)),Nx
                                        rowIndex=getIndex(iSpeciesA,ix_row,L+1,itheta,izeta,BLOCK_F)
                                        colIndex=getIndex(1,1,1,itheta,izeta,BLOCK_QN)                                                                                        
                                        call MatSetValue(matrix, rowIndex, colIndex, &
                                        -nu_n*CHatTimesf(ix_row), ADD_VALUES, ierr)
                                     end do !ix_row
                                  end if ! includeTemperatureEquilibrationTerm
                               end if ! Jacobian part
                            end do ! izeta
                         end do !itheta
                      end if ! whichMatrix > 0, iSpeciesA==iSpeciesB                               
                   end do !iSpeciesA
                end do !iSpeciesB
             end do ! L
             
             
             
             deallocate(CHatTimesf) 
             deallocate(CECDpol) 
             deallocate(CECDpolJ)  
             deallocate(fM)  
             deallocate(f1b)
             
          else ! original code (without Phi1 variations)
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
                   !CECD=0 ! For debugging
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
                      if (iSpeciesA==iSpeciesB .or. (whichMatrix>0 .and. (whichMatrix .ne. 5)) .or. preconditioner_species==0) then
                         
                         speciesFactor = sqrt(THats(iSpeciesA)*mHats(iSpeciesB) &
                         / (THats(iSpeciesB) * mHats(iSpeciesA)))
                         xb =  x * speciesFactor
                         
                         ! Build M11
                         M11 = CECD(iSpeciesA, iSpeciesB,:,:)
                         if (iSpeciesA == iSpeciesB) then
                            do i=1,Nx
                               M11(i,i) = M11(i,i) + (-oneHalf*nuDHat(iSpeciesA,i)*(L*(L+1) + Krook*2))
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
                               
                               CHat = M11 - matmul(M12 - matmul(M13, M33BackslashM32s(L+1,:,:)),&
                               M22BackslashM21s(L+1,:,:))
                            end if
                         else
                            CHat = M11;
                         end if
                         
                         if ((whichMatrix==0 .or. whichMatrix==5) .and. L >= preconditioner_x_min_L) then
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
                               !do ix_row=ixMin,Nx
                               do ix_row=max(ixMin,min_x_for_L(L)),Nx
                                  rowIndex=getIndex(iSpeciesA,ix_row,L+1,itheta,izeta,BLOCK_F)
                                  !do ix_col = ixMinCol,Nx
                                  do ix_col = max(ixMinCol,min_x_for_L(L)),Nx
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
          end if ! Phi1 part   
          
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
          
          ! ************************************************************************************
          ! This section has been modified by AI (2017-09) in order to include Phi1
          ! in the collision operator. Phi1 is included by setting 
          ! includePhi1InCollisionOperator = .true.
          ! See the documentation at
          ! https://github.com/landreman/sfincs/blob/poloidalVariationInCollisionOperator/doc/PoloidalVariationInCollisionOperator_code.pdf
          ! ************************************************************************************

          ! ************************************ !
          ! *** SECTION CHANGED BY AM 2018-01 *** !
          ! ************************************ !


          ! *** WITH PHI1 *** !
          ! With Phi1 the deflection frequency is a function of (species, velocity, theta, zeta)
          if (includePhi1InCollisionOperator .and. includePhi1 .and. includePhi1InKineticEquation) then 
             !nuDHat = zero
 !!Commented by AM 2018-01
             nuDHatpol = zero !!Added by AM 2018-01
             nuDHatpolJ = zero !!Added by AM 2018-01
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
                   
                   ! Build the pitch-angle scattering frequencies:
                   do itheta=ithetaMin,ithetaMax !!Added by AM 2018-01
                      do izeta=izetaMin,izetaMax !!Added by AM 2018-01
                         ! Generate preFactor for nHats(iSpeciesB) terms
                         preFactor = exp(-Zs(iSpeciesB)*alpha*Phi1Hat(itheta,izeta)/Thats(iSpeciesB)) !!Added by AM 2018-01
                         ! Build the pitch-angle scattering frequency:
                         nuDHatpol(iSpeciesA, :,itheta,izeta) = nuDHatpol(iSpeciesA, :,itheta,izeta) & !!Added by AM 2018-01
                              + (three*sqrtpi/four) / T32m & !!Added by AM 2018-01
                              * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) & !!Added by AM 2018-01
                              * nHats(iSpeciesB)*preFactor*(erfs - Psi_Chandra)/(x*x*x) !!Added by AM 2018-01

                         ! Generate preFactorJ for nHats(iSpeciesB) terms
                         preFactorJ = (-Zs(iSpeciesB)*alpha/Thats(iSpeciesB)) & !!Added by AM 2018-01
                              *exp(-Zs(iSpeciesB)*alpha*Phi1Hat(itheta,izeta)/Thats(iSpeciesB)) !!Added by AM 2018-01
                         ! Build the pitch-angle scattering frequency:
                         nuDHatpolJ(iSpeciesA, :,itheta,izeta) = nuDHatpolJ(iSpeciesA, :,itheta,izeta) & !!Added by AM 2018-01
                              + (three*sqrtpi/four) / T32m & !!Added by AM 2018-01
                              * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) & !!Added by AM 2018-01
                              * nHats(iSpeciesB)*preFactorJ*(erfs - Psi_Chandra)/(x*x*x) !!Added by AM 2018-01
                      end do !!Added by AM 2018-01
                   end do !!Added by AM 2018-01

                   
                end do
                
                do ix=ixMin,Nx
                   !do L=1, Nxi-1
                   do L=0, Nxi_for_x(ix)-1                      
                      ! At this point, CHat contains the collision operator normalized by
                      ! \bar{nu}, (the collision frequency at the reference mass, density, and temperature.)
                      
                      do itheta=ithetaMin,ithetaMax
                         do izeta=izetaMin,izetaMax
                           
                            CHat_element = -oneHalf*nuDHatpol(iSpeciesA, ix,itheta,izeta)*(L*(L+1) + Krook*2) !!Added by AM 2018-01

                            index=getIndex(iSpeciesA,ix,L+1,itheta,izeta,BLOCK_F)
                            call MatSetValueSparse(matrix, index, index, &
 

                            -nu_n*CHat_element, ADD_VALUES, ierr)
                            
                            ! ************************************************************************************
                            !  Calculate d(collision op.) / d Phi1 contribution to the Jacobian.
                            !  Added by AI (2017-09) 
                            !  Required since we now have a Phi1Hat dependence in the collision operator
                            ! ************************************************************************************
                            
                            
                            if (whichMatrix == 1 .or. whichMatrix == 0) then !!Added by AM 2018-01
                               
                               ! Generate pre-factor together with f1b from state vector
 
                               CHat_elementJ = (-oneHalf*nuDHatpolJ(iSpeciesA, ix,itheta,izeta)*(L*(L+1) + Krook*2))*stateArray(index + 1) !!Added by AM 2018-01

                               ! Now we need to use a different col index since we save in the phi1 part
                               colIndex=getIndex(1,1,1,itheta,izeta,BLOCK_QN)
                               call MatSetValue(matrix, index, colIndex, & 
                               -nu_n*CHat_elementJ, ADD_VALUES, ierr) !!Added by AM 2018-01
                               ! use MatSetValue, otherwise petc error
                            end if
                         end do
                      end do
                   end do
                end do
             end do
             
             ! *** WITHOUT PHI1 *** !

          else
             ! Original part without Phi1
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
                   !do L=1, Nxi-1
                   do L=0, Nxi_for_x(ix)-1
                      CHat_element = -oneHalf*nuDHat(iSpeciesA,ix)*(L*(L+1) + Krook*2)
                      
                      ! At this point, CHat contains the collision operator normalized by
                      ! \bar{nu}, (the collision frequency at the reference mass, density, and temperature.)
                      
                      do itheta=ithetaMin,ithetaMax
                         do izeta=izetaMin,izetaMax
                            index=getIndex(iSpeciesA,ix,L+1,itheta,izeta,BLOCK_F)
                            call MatSetValueSparse(matrix, index, index, &
                            -nu_n*CHat_element, ADD_VALUES, ierr)
                         end do
                      end do
                   end do
                   
                end do
             end do
             
          end if
          ! *** END OF NEW SECTION ADDED BY AM 2018-01 *** !

          ! **************************************** !
          ! *** END SECTION CHANGED BY AM 2018-01 *** !
          ! **************************************** !
          
       case default
          print *,"Error! collisionOperator must be 0 or 1."
          stop
          
       end select
    end if
  end subroutine addCollisionOperator

  ! *******************************************************************************
  ! If there is a grid point at x=0, add the boundary conditions for f at x=0.
  ! *******************************************************************************
  ! WARNING: NOT COVERED BY ANY TESTS 2020-04
  subroutine boundaryConditionAtX0(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    integer :: iSpecies, ix, itheta, izeta, L, index, ix_row, ix_col, rowIndex, colIndex
    if ((pointAtX0 .and. whichMatrix .ne. 2)) then
       ! For L > 0 modes, impose f=0 at x=0:
       ix = 1
       do L = 1,(Nxi_for_x(ix)-1)
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
       if ((whichMatrix==0 .or. whichMatrix==5) .and. L >= preconditioner_x_min_L) then
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
  end subroutine boundaryConditionAtX0

  ! *******************************************************************************
  ! Add sources:
  ! *******************************************************************************
  subroutine addSources(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    integer :: iSpecies, ix, itheta, izeta, L, rowIndex, colIndex
    PetscScalar :: xPartOfSource1, xPartOfSource2
    
    if (whichMatrix .ne. 2) then
       select case (constraintScheme)
       case (0)
          ! Do nothing here.
          
       case (1,3,4)
          ! Add a heat source and a particle source.
          
          L=0
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
                if (whichMatrix == 4 .or. whichMatrix == 5) then
                  xPartOfSource1 = exp(-x2(ix))/(pi*sqrtpi)
                  xPartOfSource2 = x2(ix)*exp(-x2(ix))/(pi*sqrtpi)
                else
                  xPartOfSource1 = (         -x2(ix) + 5/two) * exp(-x2(ix)) / (pi*sqrtpi) ! Provides particles but no heat
                  xPartOfSource2 = (two/three*x2(ix) -     1) * exp(-x2(ix)) / (pi*sqrtpi) ! Provides heat but no particles
                end if
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
  end subroutine addSources

  ! *******************************************************************************
  ! Add the density and pressure constraints:
  ! *******************************************************************************
  subroutine densityPressureConstraints(matrix, whichMatrix, whichLambda)
    Mat :: matrix
    integer, intent(in) :: whichMatrix, whichLambda
    PetscErrorCode :: ierr
    integer :: iSpecies, ix, itheta, izeta, L, rowIndex, colIndex
    PetscScalar :: factor, geometricFactor
    PetscScalar :: angle, cos_angle, dBHatdLambda
    
    if (whichMatrix .ne. 2 .and. procThatHandlesConstraints) then
       select case (constraintScheme)
       case (0)
          ! Do nothing here.

       case (1,3,4)
          ! Force the flux-surface-averaged perturbed density and 
          ! flux-surface-averaged perturbed pressure to be zero.

          L=0
          ! geometricFactor = 1/DHat = (GHat+iota*IHat)/(BHat*BHat)
          do itheta=1,Ntheta
             do izeta=1,Nzeta
                if (whichMatrix < 6) then
                   ! normal matrix
                   geometricFactor = 1/DHat(itheta,izeta)
                else
                   select case (whichLambda)
                   case (0)
                      geometricFactor = zero
                   case (1)
                      angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                      cos_angle = cos(angle)                   
                      dBHatdLambda = cos_angle
                      geometricFactor = -two*dBHatdLambda*(GHat+iota*IHat)/(BHat(itheta,izeta)**3)
                   case (2)
                      geometricFactor = iota/(BHat(itheta,izeta)*BHat(itheta,izeta))
                   case (3)
                      geometricFactor = 1/(BHat(itheta,izeta)*BHat(itheta,izeta))
                   case (4)
                      geometricFactor = IHat/(BHat(itheta,izeta)*BHat(itheta,izeta))
                   end select
                   
                end if
                factor = thetaWeights(itheta)*zetaWeights(izeta)*geometricFactor

                
                do ix=1,Nx
                   do ispecies=1,Nspecies
                      if ((whichMatrix == 4) .or. (whichMatrix==5)) then
                        colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)

                        rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
                        call MatSetValueSparse(matrix, rowIndex, colIndex, &
                             factor*xWeights(ix)*(x2(ix)*five/two-x2(ix)*x2(ix)), ADD_VALUES, ierr)

                        rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
                        call MatSetValueSparse(matrix, rowIndex, colIndex, &
                             factor*xWeights(ix)*(two/three*x2(ix)*x2(ix)-x2(ix)), ADD_VALUES, ierr)
                      else
                        colIndex = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)

                        rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
                        call MatSetValueSparse(matrix, rowIndex, colIndex, &
                             x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)

                        rowIndex = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
                        call MatSetValueSparse(matrix, rowIndex, colIndex, &
                             x2(ix)*x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)
                    end if
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
  end subroutine densityPressureConstraints
  
  ! *******************************************************************************
  ! SECTION MODIFIED BY AM 2016-02/03
  ! Add the quasineutrality equation
  ! *******************************************************************************
  subroutine quasineutrality(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    integer :: iSpecies, ix, itheta, izeta, L, rowIndex, colIndex
    PetscScalar :: speciesFactor
           
   
    if (whichMatrix .ne. 2 .and. includePhi1  .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
       L=0
       do itheta = ithetaMin,ithetaMax
          do izeta = izetaMin,izetaMax
             rowIndex = getIndex(1, 1, 1, itheta, izeta, BLOCK_QN)

             ! Add the charge density of the kinetic species (the part of the residual/Jacobian related to f1):
             do ispecies = 1,Nspecies
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
                        - alpha * Zs(ispecies) * Zs(ispecies) * NHats(ispecies) * exp (- Zs(ispecies)* alpha * Phi1Hat(itheta,izeta) / THats(ispecies)) / THats(ispecies), ADD_VALUES, ierr)
                end do
                if (withAdiabatic) then
                   call MatSetValueSparse(matrix, rowIndex, colIndex, &
                        - alpha * adiabaticZ * adiabaticZ *adiabaticNHat * exp (- adiabaticZ* alpha * Phi1Hat(itheta,izeta) / adiabaticTHat) / adiabaticTHat, ADD_VALUES, ierr)
                end if
             else if (quasineutralityOption == 2 .and. withAdiabatic .and. Nspecies > 0) then
                call MatSetValueSparse(matrix, rowIndex, colIndex, &
                     - alpha * (Zs(1)*Zs(1)*NHats(1)/THats(1) + adiabaticZ*adiabaticZ*adiabaticNHat/adiabaticTHat), ADD_VALUES, ierr)                
             end if
       !!!!!!!!!!!!!!!!!!!!!!!

             ! Add the Lagrange multiplier lambda:
             colIndex = getIndex(1, 1, 1, 1, 1, BLOCK_PHI1_CONSTRAINT)
             call MatSetValueSparse(matrix, rowIndex, colIndex, &
                  one, ADD_VALUES, ierr)

          end do
       end do
    end if
  end subroutine quasineutrality

  ! *******************************************************************************
  ! Add the constraint < Phi_1 > = 0
  ! *******************************************************************************
  subroutine phi1Constraints(matrix, whichMatrix)
    Mat :: matrix
    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    integer :: itheta, izeta
    integer, dimension(:), allocatable :: rowIndices, colIndices
  
    if (whichMatrix .ne. 2 .and. procThatHandlesConstraints .and. includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
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
  end subroutine phi1Constraints
  
  
end module DKEMatrix

  
