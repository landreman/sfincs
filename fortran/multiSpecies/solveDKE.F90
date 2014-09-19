! For compilers that do not include the error function erf(x), the line
! below should be un-commented:
!#define USE_GSL_ERF
  
#include <finclude/petsckspdef.h>
#include <finclude/petscdmdadef.h>
#include <petscversion.h>

! Next come some definitions required because the syntax for several PETSc objects
! has changed from version to version.
  
! For PETSc versions prior to 3.3, the MatCreateAIJ subroutine was called MatCreateMPIAIJ.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 3))
#define MatCreateAIJ MatCreateMPIAIJ
#endif
! Hereafter in this code, use MatCreateAIJ.

! For PETSc versions prior to 3.4, the PetscTime subroutine was called PetscGetTime.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 4))
#define PetscTime PetscGetTime
#endif
!Hereafter in this code, use PetscTime.

! For PETSc versions prior to 3.5, PETSC_DEFAULT_DOUBLE_PRECISION was used in place of PETSC_DEFAULT_REAL.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
#define PETSC_DEFAULT_REAL PETSC_DEFAULT_DOUBLE_PRECISION
#endif
!Hereafter in this code, use PETSC_DEFAULT_REAL.

! For PETSc versions prior to 3.5, DMDA_BOUNDARY_NONE was used in place of DM_BOUNDARY_NONE.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
#define DM_BOUNDARY_NONE DMDA_BOUNDARY_NONE
#endif
!Hereafter in this code, use DM_BOUNDARY_NONE.

  subroutine solveDKE()

    use globalVariables
    use polynomialDiffMatrices
    use xGrid
    use printToStdout
    use petscksp
    use petscdmda
    use sparsify
    use geometry
    use indices

    implicit none

    PetscErrorCode :: ierr
    Vec :: rhs, soln, solnOnProc0
    integer :: whichRHS
    Mat :: matrix, preconditionerMatrix
    PetscViewer :: MatlabOutput, binaryOutputViewer
    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, xPartOfRHS, speciesFactor, speciesFactor2
    !!Modified by AM 2014-09!!
    !PetscScalar :: dnHatdpsiToUse, dTHatdpsiToUse, EParallelHatToUse, dPhiHatdpsiNToUse, T32m
    PetscScalar :: EParallelHatToUse, dPhiHatdpsiNToUse, T32m
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights, B2, xb, expxb2
    PetscScalar, dimension(:,:), allocatable :: ddtheta, d2dtheta2
    PetscScalar, dimension(:,:), allocatable :: ddzeta, d2dzeta2
    PetscScalar, dimension(:,:), allocatable :: thetaPartOfTerm, localThetaPartOfTerm, xPartOfXDot
    integer :: i, j, ix, ispecies, itheta, izeta, L, NxPotentials, index
    integer :: ithetaRow, ithetaCol, scheme, ell, iSpeciesA, iSpeciesB
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar, dimension(:), allocatable :: x, xWeights, xPotentials, xWeightsPotentials
    PetscScalar, dimension(:), allocatable :: x2
    PetscScalar, dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
    PetscScalar, dimension(:,:), allocatable :: ddxPreconditioner, ddxToUse, d2dx2ToUse, zetaPartOfTerm
    PetscScalar, dimension(:,:), allocatable :: regridPolynomialToUniform, fToFInterpolationMatrix
    PetscScalar, dimension(:,:), allocatable :: potentialsToFInterpolationMatrix
    PetscScalar, dimension(:,:,:,:), allocatable :: CECD
    PetscScalar :: dtheta, xMaxNotTooSmall, BMax, BMin, xPartOfSource1, xPartOfSource2
    PetscScalar, dimension(:), allocatable :: thetaPartOfRHS
    PetscScalar, dimension(:,:), allocatable :: M11, M21, M32, LaplacianTimesX2WithoutL, nuDHat
    PetscScalar, dimension(:,:), allocatable :: xPartOfCECD, M12IncludingX0, M13IncludingX0
    PetscScalar, dimension(:), allocatable :: erfs, x3, expx2, Psi_Chandra, nuD, PsiPrime
    PetscScalar, dimension(:,:), allocatable :: CHat, M22, M33, M12, M13
    PetscScalar, dimension(:), allocatable :: diagonalOfKWithoutThetaPart
    PetscScalar, dimension(:,:), allocatable :: M22BackslashM21, M33BackslashM32
    !PetscScalar, dimension(:,:), allocatable :: fieldTerm
    PetscScalar, dimension(:,:,:), allocatable :: M22BackslashM21s, M33BackslashM32s
    integer, dimension(:), allocatable :: IPIV  ! Needed by LAPACK
    integer :: LAPACKInfo, predictedNNZForEachRowOfPreconditioner, predictedNNZForEachRowOfTotalMatrix
    integer, dimension(:), allocatable :: predictedNNZsForEachRow, predictedNNZsForEachRowDiagonal
    PetscScalar :: collisionTermFactor, xDotFactor, LFactor, temp, temp1, temp2
    integer :: rowIndex, colIndex, ixi
    PetscScalar :: densityFactor, flowFactor, pressureFactor, Phi1HatDenominator
    PetscScalar :: particleFluxFactor, momentumFluxFactor, heatFluxFactor, NTVFactor
    PetscScalar, dimension(:), allocatable :: densityIntegralWeights
    PetscScalar, dimension(:), allocatable :: flowIntegralWeights
    PetscScalar, dimension(:), allocatable :: pressureIntegralWeights
    PetscScalar, dimension(:), allocatable :: particleFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: momentumFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: heatFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: NTVIntegralWeights
    character :: trans='n'
    PetscLogDouble :: time1, time2, startTime
    KSP :: KSPInstance
    PC :: preconditionerContext
    KSPConvergedReason :: reason
    PetscScalar, pointer :: solnArray(:)
    DM :: myDM
    integer :: ithetaMin, ithetaMax, localNtheta
    VecScatter :: VecScatterContext
    logical :: procThatHandlesConstraints
    integer :: whichMatrix, whichMatrixMin, rowIndexArray(1), tempInt1, tempInt2, keepXCoupling
    PetscScalar :: singleValueArray(1), factor, zetaMax, VPrimeHatWithG
    PetscScalar, dimension(:), allocatable :: theta_preconditioner, thetaWeights_preconditioner
    PetscScalar, dimension(:,:), allocatable :: ddtheta_preconditioner, d2dtheta2_preconditioner, ddthetaToUse
    PetscScalar, dimension(:), allocatable :: zeta_preconditioner, zetaWeights_preconditioner
    PetscScalar, dimension(:,:), allocatable :: ddzeta_preconditioner, d2dzeta2_preconditioner, ddzetaToUse
    PetscScalar, dimension(:,:), allocatable :: tempMatrix, tempMatrix2, extrapMatrix
    Mat :: permutationMatrix, tempMat
    Vec :: tempVec
    double precision :: myMatInfo(MAT_INFO_SIZE)
    integer :: NNZMain, NNZPreconditioner, NNZAllocatedMain, NNZAllocatedPreconditioner
    integer :: mallocsMain, mallocsPreconditioner
    integer :: firstRowThisProcOwns, lastRowThisProcOwns, numLocalRows
    PetscScalar :: maxxPotentials, CHat_element


    !!Lines added by AM 2014-09-17
    PetscScalar, dimension(:), allocatable :: dnHatdpsiNsToUse, dTHatdpsiNsToUse
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! *******************************************************************************
    ! Do a few sundry initialization tasks:
    ! *******************************************************************************

    call PetscTime(time1, ierr)
    startTime = time1

    if ((.not. isAParallelDirectSolverInstalled) .and. (numProcsInSubComm > 1)) then
       if (masterProcInSubComm) then
          print *,"Error! Neither mumps nor superlu_dist appears to be installed,"
          print *," yet you have asked for a matrix to be distributed across processsors."
       end if
       stop
    end if

    if (forceOddNthetaAndNzeta) then
       if (mod(Ntheta, 2) == 0) then
          Ntheta = Ntheta + 1
       end if
       if (mod(Nzeta, 2) == 0) then
          Nzeta = Nzeta + 1
       end if
    end if

    call printInputs()

    matrixSize = Ntheta * Nzeta * Nxi * Nx
    select case (constraintScheme)
    case (0)
    case (1)
       matrixSize = matrixSize + 2
    case (2)
       matrixSize = matrixSize + Nx
    case default
       print *,"Error! Invalid constraintScheme"
       stop
    end select
    matrixSize = matrixSize * Nspecies

    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] The matrix is ",matrixSize,"x",matrixSize," elements."
    end if

    call validateInput()

    if (RHSMode==2) then
       print *,"Error! RHSMode 2 is not yet implemented in this version."
       stop
    end if

    transportMatrix = 0

    !!Added by AM 2014-09!!
    ArrayFirstSpeciesParticleFluxCoefficients = 0
    !!!!!!!!!!!!!!!!!!!!!!!

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Create grids, integration weights, and differentiation matrices
    !
    ! *******************************************************************************
    ! *******************************************************************************

    ! Assign a range of theta indices to each processor.
    ! This is done by creating a PETSc DM that is not actually used for anything else.
    call DMDACreate1d(MPIComm, DM_BOUNDARY_NONE, Ntheta, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)
    call DMDAGetCorners(myDM, ithetaMin, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         localNtheta, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
    ! Switch to 1-based indices:
    ithetaMin = ithetaMin + 1

    ithetaMax = ithetaMin+localNtheta-1
    procThatHandlesConstraints = masterProcInSubComm

    CHKERRQ(ierr)
    print *,"[",myCommunicatorIndex,"] Processor ",myRank," owns theta indices ",ithetaMin," to ",ithetaMax

    ! Each processor is responsible for building the rows of the matrix and rhs corresponding
    ! to its ithetaMin:ithetaMax, and each processor is resposible for all columns of the matrix.

    ! *******************************************************************************
    ! Build theta grid, integration weights, and differentiation matrices:
    ! *******************************************************************************

    select case (thetaDerivativeScheme)
    case (0)
       scheme = 20
    case (1)
       scheme = 0
    case (2)
       scheme = 10
    case default
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Error! Invalid setting for thetaDerivativeScheme"
       end if
       stop
    end select

    allocate(theta(Ntheta))
    allocate(thetaWeights(Ntheta))
    allocate(ddtheta(Ntheta,Ntheta))
    allocate(ddthetaToUse(Ntheta,Ntheta))
    allocate(d2dtheta2(Ntheta,Ntheta))
    allocate(theta_preconditioner(Ntheta))
    allocate(thetaWeights_preconditioner(Ntheta))
    allocate(ddtheta_preconditioner(Ntheta,Ntheta))
    allocate(d2dtheta2_preconditioner(Ntheta,Ntheta))

    call uniformDiffMatrices(Ntheta, 0, two*pi, scheme, theta, thetaWeights, ddtheta, d2dtheta2)

    ! If needed, also make a sparser differentiation matrix for the preconditioner:
    if (preconditioner_theta==1) then
       scheme = 0
       call uniformDiffMatrices(Ntheta, 0, two*pi, scheme, theta_preconditioner, &
            thetaWeights_preconditioner, ddtheta_preconditioner, d2dtheta2_preconditioner)

    end if

    ! *******************************************************************************
    ! Build zeta grid, integration weights, and differentiation matrices:
    ! *******************************************************************************

    call setNPeriods()

    zetaMax = 2*pi/NPeriods

    select case (thetaDerivativeScheme)
    case (0)
       scheme = 20
    case (1)
       scheme = 0
    case (2)
       scheme = 10
    case default
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Error! Invalid setting for thetaDerivativeScheme"
       end if
       stop
    end select

    allocate(zeta(Nzeta))
    allocate(zetaWeights(Nzeta))
    allocate(ddzeta(Nzeta,Nzeta))
    allocate(ddzetaToUse(Nzeta,Nzeta))
    allocate(d2dzeta2(Nzeta,Nzeta))
    allocate(zeta_preconditioner(Nzeta))
    allocate(zetaWeights_preconditioner(Nzeta))
    allocate(ddzeta_preconditioner(Nzeta,Nzeta))
    allocate(d2dzeta2_preconditioner(Nzeta,Nzeta))

    if (Nzeta==1) then
       ! Axisymmetry:
       zeta = 0
       zetaWeights = 2*pi
       ddzeta = 0
       d2dzeta2 = 0 ! d2dzeta2 is not actually used.
    else
       call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta, zetaWeights, ddzeta, d2dzeta2)
    end if

    ! If needed, also make a sparser differentiation matrix for the preconditioner:
    if (preconditioner_zeta==1) then
       if (Nzeta==1) then
          zeta_preconditioner = 0
          zetaWeights_preconditioner = 2*pi
          ddzeta_preconditioner = 0
          d2dzeta2_preconditioner = 0
       else
          scheme = 0
          call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta_preconditioner, &
               zetaWeights_preconditioner, ddzeta_preconditioner, d2dzeta2_preconditioner)
       end if
    end if

    zetaWeights = zetaWeights * NPeriods

    ! *******************************************************************************
    ! Build x grids, integration weights, and differentiation matrices.
    ! Also build interpolation matrices to map functions from one x grid to the other.
    ! *******************************************************************************

    allocate(x(Nx))
    allocate(xWeights(Nx))
    call makeXGrid(Nx, x, xWeights)
    xWeights = xWeights / exp(-x*x)
    xMaxNotTooSmall = max(x(Nx), xMax)
    allocate(x2(Nx))
    x2=x*x

    allocate(ddx(Nx,Nx))
    allocate(d2dx2(Nx,Nx))
    allocate(ddxPreconditioner(Nx,Nx))
    allocate(ddxToUse(Nx,Nx))
    allocate(d2dx2ToUse(Nx,Nx))
    call makeXPolynomialDiffMatrices(x,ddx,d2dx2)


    NxPotentials = ceiling(xMaxNotTooSmall*NxPotentialsPerVth)


    allocate(xPotentials(NxPotentials))
    allocate(xWeightsPotentials(NxPotentials))
    allocate(ddxPotentials(NxPotentials, NxPotentials))
    allocate(d2dx2Potentials(NxPotentials, NxPotentials))
    call uniformDiffMatrices(NxPotentials, zero, xMaxNotTooSmall, 12, xPotentials, &
         xWeightsPotentials, ddxPotentials, d2dx2Potentials)
    maxxPotentials = xPotentials(NxPotentials)

    allocate(regridPolynomialToUniform(NxPotentials, Nx))
    call polynomialInterpolationMatrix(Nx, NxPotentials, x, xPotentials, &
         exp(-x*x), exp(-xPotentials*xPotentials), regridPolynomialToUniform)
    !    allocate(regridUniformToPolynomial(Nx,NxPotentials))
    !    call interpolationMatrix(NxPotentials, Nx, xPotentials, x, regridUniformToPolynomial, -1, 0)

    allocate(expx2(Nx))
    expx2 = exp(-x*x)

!!$    ! We need to evaluate the error function on the x grid.
!!$    ! Some Fortran compilers have this function built in.
!!$    ! If not, call the gnu scientific library:
!!$    allocate(erfs(Nx))
!!$    do i=1,Nx
!!$#ifdef USE_GSL_ERF
!!$       call erf(x(i), temp1)
!!$#else
!!$       temp1 = erf(x(i))
!!$#endif
!!$       erfs(i) = temp1
!!$    end do

    ddxPreconditioner = 0
    select case (preconditioner_x)
    case (0)
       ! No simplification in x:
       ddxPreconditioner = ddx
    case (1)
       ! Keep only diagonal terms in x:
       do i=1,Nx
          ddxPreconditioner(i,i) = ddx(i,i)
       end do
    case (2)
       ! Keep only upper-triangular terms in x:
       do i=1,Nx
          do j=i,Nx
             ddxPreconditioner(i,j) = ddx(i,j)
          end do
       end do
    case (3)
       ! Keep only tridiagonal terms in x:
       do i=1,Nx
          do j=1,Nx
             if (abs(i-j) <= 1) then
                ddxPreconditioner(i,j) = ddx(i,j)
             end if
          end do
       end do
    case (4)
       ! Keep only diagonal and super-diagonal in x:
       do i=1,Nx
          ddxPreconditioner(i,i) = ddx(i,i)
       end do
       do i=1,(Nx-1)
          ddxPreconditioner(i,i+1) = ddx(i,i+1)
       end do
    case default
       print *,"Error! Invalid preconditioner_x"
       stop
    end select

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Evaluate the magnetic field on the (theta, zeta) grid.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    allocate(BHat(Ntheta,Nzeta))
    allocate(dBHatdtheta(Ntheta,Nzeta))
    allocate(dBHatdzeta(Ntheta,Nzeta))
    allocate(NTVKernel(Ntheta,Nzeta))

    call computeBHat()

    ! *********************************************************
    ! Compute a few quantities related to the magnetic field:
    ! *********************************************************

    VPrimeHat = 0
    do itheta=1,Ntheta
       do izeta=1,Nzeta
          VPrimeHat = VPrimeHat + thetaWeights(itheta)*zetaWeights(izeta) / (BHat(itheta,izeta) * BHat(itheta,izeta))
       end do
    end do

    FSABHat2 = 4*pi*pi/VPrimeHat

    ! *********************************************************
    ! *********************************************************
    !
    ! Now build the main matrix.
    !
    ! *********************************************************
    ! *********************************************************


    ! *********************************************************
    ! Allocate matrices:
    ! *********************************************************

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

    !predictedNNZForEachRowOfTotalMatrix = 4*(3*Nx + 5*3 + 5*3 + 5 + Nx)
    !predictedNNZForEachRowOfTotalMatrix = 4*(3*Nx + 5*3 + 5*3 + 5 + Nx + Ntheta*Nzeta)
    tempInt1 = Nspecies*Nx + 5*3 + 5*3 + 5 + 3*Nx + 2 + Nx*Ntheta*Nzeta
    if (tempInt1 > matrixSize) then
       tempInt1 = matrixSize
    end if
    predictedNNZForEachRowOfTotalMatrix = tempInt1

    predictedNNZForEachRowOfPreconditioner = predictedNNZForEachRowOfTotalMatrix

    allocate(predictedNNZsForEachRow(matrixSize))
    allocate(predictedNNZsForEachRowDiagonal(matrixSize))
    tempInt1 = 3*Nx + (Nspecies-1)*Nx + (5*3-1) + (5*3-1)
    if (tempInt1 > matrixSize) then
       tempInt1 = matrixSize
    end if
    predictedNNZsForEachRow = tempInt1

    select case (constraintScheme)
    case (0)
    case (1)
       ! The rows for the constraints have more nonzeros:
       predictedNNZsForEachRow((Nspecies*Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta*Nx
    case (2)
       ! The rows for the constraints have more nonzeros:
       predictedNNZsForEachRow((Nspecies*Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta
    case default
    end select
    predictedNNZsForEachRowDiagonal = predictedNNZsForEachRow

    if (useIterativeSolver) then
       whichMatrixMin = 0
    else
       whichMatrixMin = 1
    end if
    do whichMatrix = whichMatrixMin,1
       ! When whichMatrix = 0, build the preconditioner.
       ! When whichMatrix = 1, build the final matrix.

       ! Allocate the main global matrix:
       select case (PETSCPreallocationStrategy)
       case (0)
          ! 0 = Old method with high estimated number-of-nonzeros.
          ! This method is simpler and works consistently but uses unnecessary memory.
          if (whichMatrix==0) then
             call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
                  predictedNNZForEachRowOfPreconditioner, PETSC_NULL_INTEGER, &
                  predictedNNZForEachRowOfPreconditioner, PETSC_NULL_INTEGER, &
                  matrix, ierr)
          else
             call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
                  predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER, &
                  predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER, &
                  matrix, ierr)
          end if

       case (1)
          ! 1 = New method with lower, more precise estimated number-of-nonzeros.
          ! This method is les thoroughly tested, but it should use much less memory.

          call MatCreate(MPIComm, matrix, ierr)
          !call MatSetType(matrix, MATMPIAIJ, ierr)
          call MatSetType(matrix, MATAIJ, ierr)

          numLocalRows = PETSC_DECIDE
          call PetscSplitOwnership(MPIComm, numLocalRows, matrixSize, ierr)

          call MatSetSizes(matrix, numLocalRows, numLocalRows, PETSC_DETERMINE, PETSC_DETERMINE, ierr)

          ! We first pre-allocate assuming number-of-nonzeros = 0, because due to a quirk in PETSc,
          ! MatGetOwnershipRange only works after MatXXXSetPreallocation is called:
          if (numProcsInSubComm == 1) then
             call MatSeqAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER, ierr)
          else
             call MatMPIAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, ierr)
          end if
          
          call MatGetOwnershipRange(matrix, firstRowThisProcOwns, lastRowThisProcOwns, ierr)
          !print *,"I am proc ",myRank," and I own rows ",firstRowThisProcOwns," to ",lastRowThisProcOwns-1
          
          ! To avoid a PETSc error message, the predicted nnz for each row of the diagonal blocks must be no greater than the # of columns this proc owns:
          ! But we must not lower the predicted nnz for the off-diagonal blocks, because then the total predicted nnz for the row
          ! would be too low.
          tempInt1 = lastRowThisProcOwns - firstRowThisProcOwns
          do i=firstRowThisProcOwns+1,lastRowThisProcOwns
             if (predictedNNZsForEachRowDiagonal(i) > tempInt1) then
                predictedNNZsForEachRowDiagonal(i) = tempInt1
             end if
          end do

          ! Now, set the real estimated number-of-nonzeros:
          if (numProcsInSubComm == 1) then
             call MatSeqAIJSetPreallocation(matrix, 0, predictedNNZsForEachRow(firstRowThisProcOwns+1:lastRowThisProcOwns), ierr)
          else
             call MatMPIAIJSetPreallocation(matrix, &
                  0, predictedNNZsForEachRowDiagonal(firstRowThisProcOwns+1:lastRowThisProcOwns), &
                  0, predictedNNZsForEachRow(firstRowThisProcOwns+1:lastRowThisProcOwns), ierr)
          end if

       case default
          print *,"Error! Invalid setting for PETSCPreallocationStrategy."
          stop
       end select

       ! If any mallocs are required during matrix assembly, do not generate an error:
       !call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)

       CHKERRQ(ierr)

       ! *********************************************************
       ! Select appropriate differentiation matrices depending on
       ! whether this is the preconditioner or the final matrix:
       ! *********************************************************

       if (whichMatrix==1 .or. preconditioner_theta==0) then
          ddthetaToUse = ddtheta
       else
          ddthetaToUse = ddtheta_preconditioner
       end if

       if (whichMatrix==1 .or. preconditioner_zeta==0) then
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
          do izeta=1,Nzeta
             do itheta=1,Ntheta
                thetaPartOfTerm(itheta,:) = iota*sqrtTHat/sqrtMHat * ddthetaToUse(itheta,:) &
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
          allocate(rowIndices(Nzeta))
          allocate(colIndices(Nzeta))
          do itheta=ithetaMin, ithetaMax
             do izeta=1,Nzeta
                zetaPartOfTerm(izeta,:) = sqrtTHat/sqrtMHat * ddzetaToUse(izeta,:) / BHat(itheta,izeta)
             end do

             ! PETSc uses the opposite convention to Fortran:
             zetaPartOfTerm = transpose(zetaPartOfTerm)

             do ix=1,Nx
                do L=0,(Nxi-1)
                   do izeta = 1,Nzeta
                      rowIndices(izeta)=getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                   end do

                   ! Super-diagonal term
                   if (L < Nxi-1) then
                      ell = L + 1
                      do izeta = 1,Nzeta
                         colIndices(izeta)=getIndex(ispecies, ix, ell+1, itheta, izeta, 0)
                      end do

                      call MatSetValuesSparse(matrix, Nzeta, rowIndices, Nzeta, colIndices, &
                           (L+1)/(2*L+three)*x(ix)*zetaPartOfTerm, ADD_VALUES, ierr)
                   end if

                   ! Sub-diagonal term
                   if (L > 0) then
                      ell = L - 1
                      do izeta = 1,Nzeta
                         colIndices(izeta)=getIndex(ispecies, ix, ell+1, itheta, izeta, 0)
                      end do

                      call MatSetValuesSparse(matrix, Nzeta, rowIndices, Nzeta, colIndices, &
                           L/(2*L-one)*x(ix)*zetaPartOfTerm, ADD_VALUES, ierr)
                   end if


                end do
             end do
          end do
          deallocate(rowIndices)
          deallocate(colIndices)
          deallocate(zetaPartOfTerm)

          ! *********************************************************
          ! Add the ExB d/dtheta term:
          ! *********************************************************

          factor = alpha*Delta*GHat/(2*psiAHat)*dPhiHatdpsiN
          allocate(thetaPartOfTerm(Ntheta,Ntheta))
          allocate(localThetaPartOfTerm(Ntheta,localNtheta))
          allocate(rowIndices(localNtheta))
          allocate(colIndices(Ntheta))
          do izeta=1,Nzeta
             if (useDKESExBDrift) then
                thetaPartOfTerm = ddthetaToUse / FSABHat2
             else
                do itheta=1,Ntheta
                   thetaPartOfTerm(itheta,:) = ddthetaToUse(itheta,:) / (BHat(itheta,izeta) ** 2)
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

          factor = -alpha*Delta*IHat/(2*psiAHat)*dPhiHatdpsiN
          allocate(zetaPartOfTerm(Nzeta,Nzeta))
          allocate(rowIndices(Nzeta))
          do itheta=ithetaMin, ithetaMax
             if (useDKESExBDrift) then
                zetaPartOfTerm = ddzetaToUse / FSABHat2
             else
                do izeta=1,Nzeta
                   zetaPartOfTerm(izeta,:) = ddzetaToUse(izeta,:) / (BHat(itheta,izeta) ** 2)
                end do
             end if

             ! PETSc uses the opposite convention to Fortran:
             zetaPartOfTerm = transpose(zetaPartOfTerm*factor)

             do ix=1,Nx
                do L=0,(Nxi-1)
                   do izeta=1,Nzeta
                      rowIndices(izeta)=getIndex(ispecies,ix,L+1,itheta,izeta,0)
                   end do

                   call MatSetValuesSparse(matrix, Nzeta, rowIndices, Nzeta, rowIndices, &
                        zetaPartOfTerm, ADD_VALUES, ierr)
                end do
             end do
          end do
          deallocate(rowIndices)
          deallocate(zetaPartOfTerm)

          ! *********************************************************
          ! Add the standard mirror term:
          ! *********************************************************

          do itheta=ithetaMin,ithetaMax
             do izeta=1,Nzeta
                factor = -sqrtTHat/(2*sqrtMHat*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                     * (iota*dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta))

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
                do izeta=1,Nzeta
                   factor = alpha*Delta*dPhiHatdpsiN/(4*psiAHat*(BHat(itheta,izeta)**3)) &
                        * (GHat*dBHatdtheta(itheta,izeta) - IHat* dBHatdzeta(itheta,izeta))

                   do ix=1,Nx
                      do L=0,(Nxi-1)
                         rowIndex=getIndex(ispecies,ix,L+1,itheta,izeta,0)

                         ! Diagonal term
                         call MatSetValueSparse(matrix, rowIndex, rowIndex, &
                              (L+1)*L/((2*L-one)*(2*L+three))*factor, ADD_VALUES, ierr)

                         if (whichMatrix==1 .or. preconditioner_xi==0) then
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
                   ddxToUse = ddxPreconditioner
                else
                   ddxToUse = ddx
                end if
                do ix=1,Nx
                   xPartOfXDot(ix,:) = x(ix) * ddxToUse(ix,:)
                end do
                xPartOfXDot = transpose(xPartOfXDot)  ! PETSc uses the opposite convention of Fortran

                do itheta=ithetaMin,ithetaMax

                   do izeta=1,Nzeta
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

                      if (whichMatrix==1 .or. preconditioner_xi==0) then
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
                do izeta = 1,Nzeta
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
                   if (iSpeciesA==iSpeciesB .or. whichMatrix==1 .or. preconditioner_species==0) then

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
                         do izeta=1,Nzeta
                            do ix=1,Nx
                               rowIndices(ix)=getIndex(iSpeciesA,ix,L+1,itheta,izeta,0)
                               colIndices(ix)=getIndex(iSpeciesB,ix,L+1,itheta,izeta,0)
                            end do
                            call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                                 -nu_n*(GHat+iota*IHat)/(BHat(itheta,izeta)*BHat(itheta,izeta))*CHat, &
                                 ADD_VALUES, ierr)
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
                      do izeta=1,Nzeta
                         index=getIndex(iSpeciesA,ix,L+1,itheta,izeta,0)
                         call MatSetValueSparse(matrix, index, index, &
                              -nu_n*(GHat+iota*IHat)/(BHat(itheta,izeta)*BHat(itheta,izeta))*CHat_element, &
                              ADD_VALUES, ierr)
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
                do izeta = 1,Nzeta
                   do ispecies = 1,Nspecies
                      rowIndex = getIndex(ispecies, ix, L+1, itheta, izeta, 0)

                      colIndex = getIndex(ispecies, 1, 1, 1, 1, 1)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource1 / (BHat(itheta,izeta) ** 2), &
                           ADD_VALUES, ierr)

                      colIndex = getIndex(ispecies, 1, 1, 1, 1, 2)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource2 / (BHat(itheta,izeta) ** 2), &
                           ADD_VALUES, ierr)
                   end do
                end do
             end do
          end do

       case (2)
          ! Add a L=0 source (which is constant on the flux surface) at each x.
          L=0
          do ix=1,Nx
             do itheta=ithetaMin,ithetaMax
                do izeta = 1,Nzeta
                   do ispecies = 1,Nspecies
                      rowIndex = getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                      colIndex = getIndex(ispecies, ix, 1, 1, 1, 3)
                      call MatSetValue(matrix, rowIndex, colIndex, one / (BHat(itheta,izeta) ** 2), &
                           ADD_VALUES, ierr)
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
       ! Now finalize the matrices:
       ! *******************************************************************************

       call PetscTime(time2, ierr)
       if (masterProcInSubComm) then
          if (whichMatrix==0) then
             print *,"[",myCommunicatorIndex,"] Time to pre-assemble preconditioner matrix: ", time2-time1, " seconds."
          else
             print *,"[",myCommunicatorIndex,"] Time to pre-assemble matrix: ", time2-time1, " seconds."
          end if
       end if
       call PetscTime(time1, ierr)

       call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
       call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
       CHKERRQ(ierr)

       if (whichMatrix==0) then
          preconditionerMatrix = matrix
       end if

       call PetscTime(time2, ierr)
       if (masterProcInSubComm) then
          if (whichMatrix==0) then
             print *,"[",myCommunicatorIndex,"] Time to assemble preconditioner matrices: ", time2-time1, " seconds."
          else
             print *,"[",myCommunicatorIndex,"] Time to assemble matrices: ", time2-time1, " seconds."
          end if
       end if
       call PetscTime(time1, ierr)

    end do



    ! ***********************************************************************
    ! ***********************************************************************
    ! 
    !  Initialize the solver:
    !
    ! ***********************************************************************
    ! ***********************************************************************

    call KSPCreate(MPIComm, KSPInstance, ierr)
    CHKERRQ(ierr)

    if (useIterativeSolver) then
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
       ! Syntax for PETSc versions up through 3.4
       call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, SAME_PRECONDITIONER, ierr)
#else
       ! Syntax for PETSc version 3.5 and later
       call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, ierr)
#endif
       CHKERRQ(ierr)
       call KSPGetPC(KSPInstance, preconditionerContext, ierr)
       CHKERRQ(ierr)
       call PCSetType(preconditionerContext, PCLU, ierr)
       CHKERRQ(ierr)
       !call KSPSetType(KSPInstance, KSPBCGSL, ierr)  ! Set the Krylov solver algorithm to BiCGStab(l)
       call KSPSetType(KSPInstance, KSPGMRES, ierr)   ! Set the Krylov solver algorithm to GMRES
       call KSPGMRESSetRestart(KSPInstance, 500, ierr)
       CHKERRQ(ierr)
       call KSPSetTolerances(KSPInstance, solverTolerance, PETSC_DEFAULT_REAL, &
            PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
       CHKERRQ(ierr)
       call KSPSetFromOptions(KSPInstance, ierr)
       CHKERRQ(ierr)
       call KSPMonitorSet(KSPInstance, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
    else
       ! Direct solver:
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
       ! Syntax for PETSc versions up through 3.4
       call KSPSetOperators(KSPInstance, matrix, matrix, SAME_PRECONDITIONER, ierr)
#else
       ! Syntax for PETSc version 3.5 and later
       call KSPSetOperators(KSPInstance, matrix, matrix, ierr)
#endif
       CHKERRQ(ierr)
       call KSPGetPC(KSPInstance, preconditionerContext, ierr)
       CHKERRQ(ierr)
       call PCSetType(preconditionerContext, PCLU, ierr)
       CHKERRQ(ierr)
       call KSPSetType(KSPInstance, KSPPREONLY, ierr)
       CHKERRQ(ierr)
       call KSPSetFromOptions(KSPInstance, ierr)
       CHKERRQ(ierr)
    end if

    if (numProcsInSubComm > 1) then
       select case (whichParallelSolverToFactorPreconditioner)
       case (1)
          call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERMUMPS, ierr)
          if (masterProcInSubComm) then
             print *,"[",myCommunicatorIndex,"] Using mumps to factorize the preconditioner."
          end if
       case (2)
          call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERSUPERLU_DIST, ierr)
          if (masterProcInSubComm) then
             print *,"[",myCommunicatorIndex,"] Using superlu_dist to factorize the preconditioner."
          end if
       case default
          if (masterProcInSubComm) then
             print *,"Error: Invalid setting for whichParallelSolverToFactorPreconditioner"
             stop
          end if
       end select
    else
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Using PETSc's serial sparse direct solver to factorize the preconditioner."
       end if

       ! When using PETSc's built-in solver (which is only done when running with a single processor),
       ! I originally always got an error message that there was a zero pivot. The following line seems to solve
       ! this problem.  The "zero pivot" error seemed to only arise with the "nd" ordering (nested dissection), which is the 
       ! default.  I'm not sure which of the other orderings is most efficient. I picked "rcm" for no particular reason.
       call PCFactorSetMatOrderingType(preconditionerContext, MATORDERINGRCM, ierr)

       ! I'm not sure this next line actually accomplishes anything, since it doesn't solve the "zero pivot" problem.
       ! But it doesn't seem to cost anything, and perhaps it reduces the chance of getting the "zero pivot" error in the future.
       call PCFactorReorderForNonzeroDiagonal(preconditionerContext, 1d-12, ierr) 

    end if

    call MatGetInfo(matrix, MAT_GLOBAL_SUM, myMatInfo, ierr)
    NNZMain = nint(myMatInfo(MAT_INFO_NZ_USED))
    NNZAllocatedMain = nint(myMatInfo(MAT_INFO_NZ_ALLOCATED))
    mallocsMain = nint(myMatInfo(MAT_INFO_MALLOCS))
    if (useIterativeSolver) then
       call MatGetInfo(preconditionerMatrix, MAT_GLOBAL_SUM, myMatInfo, ierr)
       NNZPreconditioner = nint(myMatInfo(MAT_INFO_NZ_USED))
       NNZAllocatedPreconditioner = nint(myMatInfo(MAT_INFO_NZ_ALLOCATED))
       mallocsPreconditioner = nint(myMatInfo(MAT_INFO_MALLOCS))
    end if
    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] # of nonzeros in matrix:",NNZMain, ", allocated:",NNZAllocatedMain, &
            ", mallocs:",mallocsMain," (should be 0)"
       if (useIterativeSolver) then
          print *,"[",myCommunicatorIndex,"] # of nonzeros in preconditioner:", NNZPreconditioner,&
               ", allocated:",NNZAllocatedPreconditioner,", mallocs:",mallocsPreconditioner, " (should be 0)"
          print *,"[",myCommunicatorIndex,"] nnz(preconditioner)/nnz(matrix):",((one*NNZPreconditioner)/NNZMain)
       end if
    end if


    ! ***********************************************************************
    ! Create some objects that will be needed soon for the solve:
    ! ***********************************************************************

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, rhs, ierr)
    CHKERRQ(ierr)
    call VecDuplicate(rhs, soln, ierr)
    CHKERRQ(ierr)

    call VecScatterCreateToZero(soln, VecScatterContext, solnOnProc0, ierr)
    CHKERRQ(ierr)

    allocate(FSADensityPerturbation(Nspecies))
    allocate(FSABFlow(Nspecies))
    allocate(FSAPressurePerturbation(Nspecies))
    allocate(particleFlux(Nspecies))
    allocate(momentumFlux(Nspecies))
    allocate(heatFlux(Nspecies))
    allocate(NTV(Nspecies)) 

    allocate(densityPerturbation(Nspecies,Ntheta,Nzeta))
    allocate(flow(Nspecies,Ntheta,Nzeta))
    allocate(pressurePerturbation(Nspecies,Ntheta,Nzeta))
    allocate(particleFluxBeforeSurfaceIntegral(Nspecies,Ntheta,Nzeta))
    allocate(momentumFluxBeforeSurfaceIntegral(Nspecies,Ntheta,Nzeta))
    allocate(heatFluxBeforeSurfaceIntegral(Nspecies,Ntheta,Nzeta))
    allocate(NTVBeforeSurfaceIntegral(Nspecies,Ntheta,Nzeta)) 

    allocate(jHat(Ntheta,Nzeta))
    allocate(Phi1Hat(Ntheta,Nzeta))

    allocate(densityIntegralWeights(Nx))
    allocate(flowIntegralWeights(Nx))
    allocate(pressureIntegralWeights(Nx))
    allocate(particleFluxIntegralWeights(Nx))
    allocate(momentumFluxIntegralWeights(Nx))
    allocate(heatFluxIntegralWeights(Nx))
    allocate(NTVIntegralWeights(Nx))


    ! ***********************************************************************
    ! ***********************************************************************
    ! 
    !  Beginning of the main solver loop:
    !
    ! ***********************************************************************
    ! ***********************************************************************
    !MODIFIED BY AM FROM HERE
    !Nspecies contains number of species
    allocate(dnHatdpsiNsToUse(Nspecies))
    allocate(dTHatdpsiNsToUse(Nspecies))

    select case (RHSMode)
    case (1)
       numRHSs = 1
    case (2)
       numRHSs = 3
    case default
       print *,"Error! Invalid setting for RHSMode."
       stop
    end select

    do whichRHS = 1,numRHSs
       if (RHSMode /= 1 .and. masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Solving system with right-hand side ",whichRHS," of ",numRHSs
       end if

       ! *******************************************************************************
       ! *******************************************************************************
       !
       ! Create the right-hand side vector
       !
       ! *******************************************************************************
       ! *******************************************************************************

       !dnHatdpsiToUse = dnHatdpsiNs(1)
       !dTHatdpsiToUse = dTHatdpsiNs(1)
       EParallelHatToUse = EParallelHat
       dPhiHatdpsiNToUse = dPhiHatdpsiN

       dnHatdpsiNsToUse(:) = 0
       dTHatdpsiNsToUse(:) = 0

       select case (RHSMode)
       case (1)
          ! Single solve: nothing more to do here.

       case (2)
          !print *,"RHSMode=2 is not yet implemented"
          print *,"Solve  for 3 linearly independent right-hand sides to get L_{11}^{11}, L_{11}^{12} and L_{12}.\n"
          !stop

          ! Solve for 3 linearly independent right-hand sides to get the full 3x3 transport matrix:
          dPhiHatdpsiNToUse = 0
          EParallelHatToUse = 0
          select case (whichRHS)
          case (1)
             !dnHatdpsiToUse = 1
             !dTHatdpsiToUse = 0
             !EParallelHatToUse = 0
             dnHatdpsiNsToUse(1) = 1
          case (2)
             ! The next 2 lines ensure (1/n)*dn/dpsi + (3/2)*dT/dpsi = 0 while dT/dpsi is nonzero.
             !dnHatdpsiToUse = (3/two)*nHats(1)/THats(1)
             !dTHatdpsiToUse = 1
             !EParallelHatToUse = 0

             !!Added by AM 2014-09!!
             if (Nspecies < 2) then !!Can not do solve because there is only 1 species
                print *,"WARNING! Trying to calculate transport coefficients with only 1 species.\n"
                ArrayFirstSpeciesParticleFluxCoefficients(2) = 0
                cycle
             end if
             !!!!!!!!!!!!!!!!!!!!!!!

             dnHatdpsiNsToUse(2) = 1
          case (3)
             !dnHatdpsiToUse = 0
             !dTHatdpsiToUse = 0
             !EParallelHatToUse = 1
             dnHatdpsiNsToUse(1) = (3/two)*nHats(1)/THats(1)
             dTHatdpsiNsToUse(1) = 1

             !!Added by AM 2014-09!!
             if (Nspecies > 1) then
                dnHatdpsiNsToUse(2) = (3/two)*nHats(2)/THats(2)
                dTHatdpsiNsToUse(2) = 1
             end if
             !!!!!!!!!!!!!!!!!!!!!!!

          case default
             print *,"Program should not get here"
             stop
          end select
       end select

       call VecSet(rhs, zero, ierr)

       ! First add the term arising from radial gradients:
       CHKERRQ(ierr)
       x2 = x*x
       do ispecies = 1,Nspecies
          THat = THats(ispecies)
          mHat = mHats(ispecies)
          sqrtTHat = sqrt(THat)
          sqrtMHat = sqrt(mHat)

          do ix=1,Nx
             xPartOfRHS = x2(ix)*exp(-x2(ix))*( dnHatdpsiNsToUse(ispecies)/nHats(ispecies) &
                  + alpha*Zs(ispecies)/THats(ispecies)*dPhiHatdpsiNToUse &
                  + (x2(ix) - three/two)*dTHatdpsiNsToUse(ispecies)/THats(ispecies))
             do itheta = ithetaMin,ithetaMax
                do izeta = 1,Nzeta

                   factor = Delta*nHats(ispecies)*mHat*sqrtMHat &
                        /(2*pi*sqrtpi*Zs(ispecies)*psiAHat*(BHat(itheta,izeta)**3)*sqrtTHat) &
                        *(GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))&
                        *xPartOfRHS

                   L = 0
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                   call VecSetValue(rhs, index, (4/three)*factor, INSERT_VALUES, ierr)

                   L = 2
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                   call VecSetValue(rhs, index, (two/three)*factor, INSERT_VALUES, ierr)
                end do
             end do
          end do
       end do
       CHKERRQ(ierr)

       ! Add the inductive electric field term:
       L=1
       do ispecies = 1,Nspecies
          do ix=1,Nx
             factor = alpha*Zs(ispecies)*x(ix)*exp(-x2(ix))*EParallelHatToUse*(GHat+iota*IHat)&
                  *nHats(ispecies)*mHats(ispecies)/(pi*sqrtpi*THats(ispecies)*THats(ispecies)*FSABHat2)
             do itheta=ithetaMin,ithetaMax
                do izeta = 1,Nzeta
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)
                   call VecSetValue(rhs, index, &
                        factor/BHat(itheta,izeta), INSERT_VALUES, ierr)
                end do
             end do
          end do
       end do

       ! Done inserting values.
       ! Finally, assemble the RHS vector:
       call VecAssemblyBegin(rhs, ierr)
       call VecAssemblyEnd(rhs, ierr)

       ! ***********************************************************************
       ! ***********************************************************************
       ! 
       !  Solve the main linear system:
       !
       ! ***********************************************************************
       ! ***********************************************************************
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Beginning the main solve.  This could take a while ..."
       end if

       if (solveSystem) then
          call KSPSolve(KSPInstance, rhs, soln, ierr)
       end if
       CHKERRQ(ierr)

       call PetscTime(time2, ierr)
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Done with the main solve.  Time to solve: ", time2-time1, " seconds."
       end if
       call PetscTime(time1, ierr)

       if (useIterativeSolver) then
          call KSPGetConvergedReason(KSPInstance, reason, ierr)
          if (reason>0) then
             if (masterProcInSubComm) then
                print *,"[",myCommunicatorIndex,"] Converged!  KSPConvergedReason = ", reason
             end if
             didItConverge = integerToRepresentTrue
          else
             if (masterProcInSubComm) then
                print *,"[",myCommunicatorIndex,"] Did not converge :(   KSPConvergedReason = ", reason
             end if
             didItConverge = integerToRepresentFalse
          end if
       else
          didItConverge = integerToRepresentTrue
       end if


       !**************************************************************************
       !**************************************************************************
       ! 
       !  Calculate moments of the solution:
       !
       !**************************************************************************
       !**************************************************************************


       ! First, send the entire solution vector to the master process:
       call VecScatterBegin(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecScatterEnd(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
       CHKERRQ(ierr)


       if (masterProcInSubComm) then
          ! All computation of moments of the distribution function is then done on the master process:

          densityPerturbation=0
          flow=0
          pressurePerturbation=0
          particleFluxBeforeSurfaceIntegral=0
          momentumFluxBeforeSurfaceIntegral=0
          heatFluxBeforeSurfaceIntegral=0
          NTVBeforeSurfaceIntegral=0

          FSADensityPerturbation=0
          FSABFlow=0
          FSAPressurePerturbation=0
          particleFlux=0
          momentumFlux=0
          heatFlux=0
          NTV=0 
          jHat=0
          Phi1Hat=0
          Phi1HatDenominator = 0

          densityIntegralWeights = x*x
          flowIntegralWeights = x*x*x
          pressureIntegralWeights = x*x*x*x
          particleFluxIntegralWeights = x*x*x*x
          momentumFluxIntegralWeights = x*x*x*x*x
          heatFluxIntegralWeights = x*x*x*x*x*x
          NTVIntegralWeights = x*x*x*x 

          ! Convert the PETSc vector into a normal Fortran array:
          call VecGetArrayF90(solnOnProc0, solnArray, ierr)
          CHKERRQ(ierr)

          if (whichRHS == numRHSs) then
             select case (constraintScheme)
             case (0)
             case (1)
                allocate(sources(Nspecies,2))
                do ispecies = 1,Nspecies
                   sources(ispecies,1) = solnArray(getIndex(ispecies, 1, 1, 1, 1, 1)+1)
                   sources(ispecies,2) = solnArray(getIndex(ispecies, 1, 1, 1, 1, 2)+1)
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                end do
             case (2)
                allocate(sources(Nspecies,Nx))
                do ispecies = 1,Nspecies
                   do ix=1,Nx
                      sources(ispecies,ix) = solnArray(getIndex(ispecies, ix, 1, 1, 1, 3)+1)
                      ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                   end do
                end do
             case default
                print *,"Error! Invalid setting for constraintScheme."
                stop
             end select
          end if

          allocate(B2(Ntheta))
          do ispecies = 1,Nspecies
             THat = THats(ispecies)
             mHat = mHats(ispecies)
             sqrtTHat = sqrt(THat)
             sqrtMHat = sqrt(mHat)

             densityFactor = 4*pi*THat*sqrtTHat/(mHat*sqrtMHat)
             flowFactor = 4*pi*THat*THat/(three*mHat*mHat)
             pressureFactor = 8*pi*THat*THat*sqrtTHat/(three*mHat*sqrtMHat)
             particleFluxFactor = - pi*Delta*THat*THat*sqrtTHat/(Zs(ispecies)*VPrimeHat*mHat*sqrtMHat*(GHat+iota*IHat))
             momentumFluxFactor = - pi*Delta*THat*THat*THat/(Zs(ispecies)*VPrimeHat*mHat*(GHat+iota*IHat))
             heatFluxFactor = - pi*Delta*THat*THat*THat*sqrtTHat/(2*Zs(ispecies)*VPrimeHat*mHat*sqrtMHat*(GHat+iota*IHat))
             NTVFactor = 4*pi*THat*THat*sqrtTHat/(mHat*sqrtMHat*VPrimeHat*(GHat+iota*IHat))

             do itheta=1,Ntheta
                do izeta=1,Nzeta

                   factor = (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))/(BHat(itheta,izeta) ** 3)

                   do ix=1,Nx
                      L = 0
                      index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                      ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                      densityPerturbation(ispecies,itheta,izeta) = densityPerturbation(ispecies,itheta,izeta) &
                           + densityFactor*xWeights(ix)*densityIntegralWeights(ix)*solnArray(index)

                      pressurePerturbation(ispecies,itheta,izeta) = pressurePerturbation(ispecies,itheta,izeta) &
                           + pressureFactor*xWeights(ix)*pressureIntegralWeights(ix)*solnArray(index)

                      particleFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           = particleFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           + factor * (8/three) * particleFluxFactor &
                           * xWeights(ix)*particleFluxIntegralWeights(ix)*solnArray(index)

                      heatFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           = heatFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           + factor * (8/three) * heatFluxFactor &
                           * xWeights(ix)*heatFluxIntegralWeights(ix)*solnArray(index)

                      L = 1
                      index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                      ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                      flow(ispecies,itheta,izeta) = flow(ispecies,itheta,izeta) &
                           + flowFactor*xWeights(ix)*flowIntegralWeights(ix)*solnArray(index)

                      momentumFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           = momentumFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           + factor * (16d+0/15) * momentumFluxFactor &
                           * xWeights(ix)*momentumFluxIntegralWeights(ix)*solnArray(index)

                      L = 2
                      index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                      ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                      particleFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           = particleFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           + factor * (four/15) * particleFluxFactor &
                           * xWeights(ix)*particleFluxIntegralWeights(ix)*solnArray(index)

                      heatFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           = heatFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           + factor * (four/15) * heatFluxFactor &
                           * xWeights(ix)*heatFluxIntegralWeights(ix)*solnArray(index)

                      NTVBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           = NTVFactor * NTVKernel(itheta,izeta)&
                           * xWeights(ix)*NTVIntegralWeights(ix)*solnArray(index) 

                      L = 3
                      index = getIndex(ispecies, ix, L+1, itheta, izeta, 0)+1
                      ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                      momentumFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           = momentumFluxBeforeSurfaceIntegral(ispecies,itheta,izeta) &
                           + factor * (four/35) * momentumFluxFactor &
                           * xWeights(ix)*momentumFluxIntegralWeights(ix)*solnArray(index)

                   end do
                end do
             end do

             do izeta=1,Nzeta
                B2 = BHat(:,izeta)*BHat(:,izeta)

                FSADensityPerturbation(ispecies) = FSADensityPerturbation(ispecies) + zetaWeights(izeta) &
                     * dot_product(thetaWeights, densityPerturbation(ispecies,:,izeta)/B2)

                FSABFlow(ispecies) = FSABFlow(ispecies) + zetaWeights(izeta) &
                     * dot_product(thetaWeights, flow(ispecies,:,izeta)/BHat(:,izeta))

                FSAPressurePerturbation(ispecies) = FSAPressurePerturbation(ispecies) + zetaWeights(izeta) &
                     * dot_product(thetaWeights, pressurePerturbation(ispecies,:,izeta)/B2)

                particleFlux(ispecies) = particleFlux(ispecies) + zetaWeights(izeta) &
                     * dot_product(thetaWeights, particleFluxBeforeSurfaceIntegral(ispecies,:,izeta))

                momentumFlux(ispecies) = momentumFlux(ispecies) + zetaWeights(izeta) &
                     * dot_product(thetaWeights, momentumFluxBeforeSurfaceIntegral(ispecies,:,izeta))

                heatFlux(ispecies) = heatFlux(ispecies) + zetaWeights(izeta) &
                     * dot_product(thetaWeights, heatFluxBeforeSurfaceIntegral(ispecies,:,izeta))

                NTV(ispecies) = NTV(ispecies) + zetaWeights(izeta) &
                     * dot_product(thetaWeights, NTVBeforeSurfaceIntegral(ispecies,:,izeta)) 

             end do

             jHat = jHat + Zs(ispecies)*flow(ispecies,:,:)
             Phi1Hat = Phi1Hat + Zs(ispecies)*densityPerturbation(ispecies,:,:)
             Phi1HatDenominator = Phi1HatDenominator + Zs(ispecies)*Zs(ispecies)*nHats(ispecies)/THats(ispecies)
          end do
          deallocate(B2)

          Phi1Hat = Phi1Hat / (alpha * Phi1HatDenominator)

          FSADensityPerturbation = FSADensityPerturbation / VPrimeHat
          FSABFlow = FSABFlow / VPrimeHat
          FSAPressurePerturbation = FSAPressurePerturbation / VPrimeHat
          FSABjHat = dot_product(Zs(1:Nspecies), FSABFlow)

          !!Section modified by AM 2014-09!!
          if (RHSMode==2) then
             !!VPrimeHatWithG = VPrimeHat*(GHat+iota*IHat)
             select case (whichRHS)
             case (1)
!!$                transportMatrix(1,1) = 4*(GHat+iota*IHat)*particleFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*sqrtTHat)*GHat)
!!$                transportMatrix(2,1) = 8*(GHat+iota*IHat)*heatFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*THat*sqrtTHat)*GHat)
!!$                transportMatrix(3,1) = 2*nHat*FSABFlow/(GHat*THat)
                ArrayFirstSpeciesParticleFluxCoefficients(1) = particleFlux*4*(Zs(1)**2)*psiAHat*(GHat+iota*IHat)*B0OverBBar*sqrt(THats(1)/mHats(1)) / ( (GHat**2)*(THats(1)**2)*(Delta**2) )
             case (2)
!!$                transportMatrix(1,2) = 4*(GHat+iota*IHat)*particleFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*GHat)
!!$                transportMatrix(2,2) = 8*(GHat+iota*IHat)*heatFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*THat*GHat)
!!$                transportMatrix(3,2) = 2*FSABFlow/(GHat)
                ArrayFirstSpeciesParticleFluxCoefficients(2) = particleFlux*(NHats(2)/NHats(1))*4*(Zs(1)**2)*psiAHat*(GHat+iota*IHat)*B0OverBBar*sqrt(THats(1)/mHats(1)) / ( (GHat**2)*(THats(1)**2)*(Delta**2) )
             case (3)
!!$                transportMatrix(1,3) = particleFlux*Delta*Delta*FSABHat2/(VPrimeHatWithG*GHat*psiAHat*omega)
!!$                transportMatrix(2,3) = 2*Delta*Delta*heatFlux*FSABHat2/(GHat*VPrimeHatWithG*psiAHat*THat*omega)
!!$                transportMatrix(3,3) = FSABFlow*Delta*Delta*sqrtTHat*FSABHat2/((GHat+iota*IHat)*2*psiAHat*omega*B0OverBBar)
                ArrayFirstSpeciesParticleFluxCoefficients(3) = particleFlux*4*(Zs(1)**2)*psiAHat*(GHat+iota*IHat)*B0OverBBar*sqrt(THats(1)/mHats(1)) / ( (GHat**2)*(THats(1)*NHats(1))*(Delta**2) )
             end select
          end if
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)
          CHKERRQ(ierr)

       end if

    end do

    !!Added by AM 2014-09
    deallocate(dnHatdpsiNsToUse)
    deallocate(dTHatdpsiNsToUse)
    !!!!!!!!!!!!!!!!!!!!!

    ! ***********************************************************************
    ! ***********************************************************************
    ! 
    !  End of the main solver loop.
    !
    ! ***********************************************************************
    ! ***********************************************************************


    ! *********************************************************
    ! Create a PETSc viewer to record output
    ! *********************************************************

    if (saveMatlabOutput) then
       call PetscViewerASCIIOpen(MPIComm, &
            & MatlabOutputFilename,&
            & MatlabOutput, ierr)
       CHKERRQ(ierr)
       call PetscViewerSetFormat(MatlabOutput, PETSC_VIEWER_ASCII_MATLAB, ierr)
       CHKERRQ(ierr)

       call PetscObjectSetName(rhs, "rhs", ierr)
       CHKERRQ(ierr)
       call VecView(rhs, MatlabOutput, ierr)
       CHKERRQ(ierr)
       call PetscObjectSetName(soln, "soln", ierr)
       CHKERRQ(ierr)
       call VecView(soln, MatlabOutput, ierr)
       CHKERRQ(ierr)

       call PetscObjectSetName(matrix, "matrix", ierr)
       CHKERRQ(ierr)
       call MatView(matrix, MatlabOutput, ierr)
       CHKERRQ(ierr)
       if (useIterativeSolver) then
          call PetscObjectSetName(preconditionerMatrix, "preconditionerMatrix", ierr)
          CHKERRQ(ierr)
          call MatView(preconditionerMatrix, MatlabOutput, ierr)
          CHKERRQ(ierr)
       end if

       call PetscTime(time2, ierr)
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Time to write output: ", time2-time1, " seconds."
       end if
       call PetscTime(time1, ierr)

       call PetscViewerDestroy(MatlabOutput, ierr)
       CHKERRQ(ierr)
    end if

    if (saveMatricesAndVectorsInBinary) then
       call PetscViewerBinaryOpen(MPIComm, &
            & trim(binaryOutputFilename) // '_rhs', &
            & FILE_MODE_WRITE, &
            & binaryOutputViewer, ierr)
       CHKERRQ(ierr)
       call VecView(rhs, binaryOutputViewer, ierr)       
       call PetscViewerDestroy(binaryOutputViewer, ierr)       

       call PetscViewerBinaryOpen(MPIComm, &
            & trim(binaryOutputFilename) // '_matrix', &
            & FILE_MODE_WRITE, &
            & binaryOutputViewer, ierr)
       CHKERRQ(ierr)
       call MatView(matrix, binaryOutputViewer, ierr)       
       call PetscViewerDestroy(binaryOutputViewer, ierr)       

       if (useIterativeSolver) then
          call PetscViewerBinaryOpen(MPIComm, &
               & trim(binaryOutputFilename) // '_pc', &
               & FILE_MODE_WRITE, &
               & binaryOutputViewer, ierr)
          CHKERRQ(ierr)
          call MatView(preconditionerMatrix, binaryOutputViewer, ierr)       
          call PetscViewerDestroy(binaryOutputViewer, ierr)       
       end if

       call PetscViewerBinaryOpen(MPIComm, &
            & trim(binaryOutputFilename) // '_soln', &
            & FILE_MODE_WRITE, &
            & binaryOutputViewer, ierr)
       CHKERRQ(ierr)
       call VecView(soln, binaryOutputViewer, ierr)       
       call PetscViewerDestroy(binaryOutputViewer, ierr)       
    end if

    !    call VecDestroy(rhs, ierr)
    call VecScatterDestroy(VecScatterContext, ierr)
    call VecDestroy(soln, ierr)
    CHKERRQ(ierr)
    call MatDestroy(matrix, ierr)
    if (useIterativeSolver) then
       call MatDestroy(preconditionerMatrix, ierr)
    end if
    call KSPDestroy(KSPInstance, ierr)
    CHKERRQ(ierr)


    call PetscTime(time2, ierr)
    elapsedTime = time2 - startTime

    call printOutputs()

  end subroutine solveDKE

