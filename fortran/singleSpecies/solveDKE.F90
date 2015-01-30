! Some compilers (e.g. gfortran) include a built-in erf() function, but others do not. 
! You should set this variable accordingly:
#define ERF_TO_USE 0
! 0:         Use compiler's built-in erf()
! 1:         Use Gnu Scientific Library's erf()
! otherwise: Use a rough approximation erf(x) \approx x / (exp(-2*x)+x)

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

    implicit none

    PetscErrorCode :: ierr, ierr2
    Vec :: rhs, soln, solnOnProc0
    integer :: whichRHS
    Mat :: matrix, preconditionerMatrix
    PetscViewer :: MatlabOutput, binaryOutputViewer
    PetscScalar :: dnHatdpsiToUse, dTHatdpsiToUse, EHatToUse, dPhiHatdpsiToUse
    PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights, B2
    PetscScalar, dimension(:,:), allocatable :: ddtheta, d2dtheta2
    PetscScalar, dimension(:,:), allocatable :: ddzeta, d2dzeta2
    PetscScalar, dimension(:,:), allocatable :: thetaPartOfTerm, localThetaPartOfTerm, xPartOfXDot
    integer :: i, j, ix, itheta, izeta, L, NxPotentials, matrixSize, index, ixi
    integer :: ithetaRow, ithetaCol, scheme, ell
    PetscScalar, dimension(:), allocatable :: xWeights, xPotentials, xWeightsPotentials
    PetscScalar, dimension(:), allocatable :: x2, xPartOfRHS
    PetscScalar, dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
    PetscScalar, dimension(:,:), allocatable :: ddxPreconditioner, ddxToUse, zetaPartOfTerm
    PetscScalar, dimension(:,:), allocatable :: regridPolynomialToUniform, regridUniformToPolynomial
    PetscScalar :: dtheta, xMaxNotTooSmall, BMax, BMin, xPartOfSource1, xPartOfSource2
    PetscScalar, dimension(:), allocatable :: thetaPartOfRHS
    integer, dimension(:), allocatable :: indices, rowIndices, colIndices
    PetscScalar, dimension(:,:), allocatable :: M11, M21, M32, LaplacianTimesX2WithoutL
    PetscScalar, dimension(:,:), allocatable :: xPartOfCECD, M12IncludingX0, M13IncludingX0
    PetscScalar, dimension(:), allocatable :: erfs, x3, expx2, Psi_Chandra, nuD, PsiPrime
    PetscScalar, dimension(:,:), allocatable :: KWithoutThetaPart, M22, M33, M12, M13
    PetscScalar, dimension(:), allocatable :: diagonalOfKWithoutThetaPart
    PetscScalar, dimension(:,:), allocatable :: M22BackslashM21, M33BackslashM32, fieldTerm
    integer, dimension(:), allocatable :: IPIV  ! Needed by LAPACK
    integer :: LAPACKInfo, predictedNNZForEachRowOfPreconditioner, predictedNNZForEachRowOfTotalMatrix
    integer, dimension(:), allocatable :: predictedNNZsForEachRow, predictedNNZsForEachRowDiagonal
    PetscScalar :: collisionTermFactor, xDotFactor, LFactor, temp1, temp2
    integer :: rowIndex, colIndex
    PetscScalar :: densityFactor, flowFactor, pressureFactor
    PetscScalar :: particleFluxFactor, momentumFluxFactor, heatFluxFactor, NTVFactor, fNormFactor
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
    PetscScalar :: singleValueArray(1), sqrtTHat, factor, zetaMax, VPrimeHatWithG
    PetscScalar, dimension(:), allocatable :: theta_preconditioner, thetaWeights_preconditioner
    PetscScalar, dimension(:,:), allocatable :: ddtheta_preconditioner, d2dtheta2_preconditioner, ddthetaToUse
    PetscScalar, dimension(:), allocatable :: zeta_preconditioner, zetaWeights_preconditioner
    PetscScalar, dimension(:,:), allocatable :: ddzeta_preconditioner, d2dzeta2_preconditioner, ddzetaToUse
    Mat :: permutationMatrix, tempMat
    Vec :: tempVec
    double precision :: myMatInfo(MAT_INFO_SIZE)
    integer :: NNZMain, NNZPreconditioner, NNZAllocatedMain, NNZAllocatedPreconditioner
    integer :: mallocsMain, mallocsPreconditioner
    integer :: firstRowThisProcOwns, lastRowThisProcOwns, numLocalRows

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

    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] The matrix is ",matrixSize,"x",matrixSize," elements."
    end if

    call validateInput()
    sqrtTHat = sqrt(THat)

    transportMatrix = 0
    transportCoeffs = 0
    NTVMatrix = 0

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
    if (Nx==1) then  !Mono-energetic calculation
       xWeights = exp(1.0)
       x=1.0
    else
       call makeXGrid(Nx, x, xWeights)
       xWeights = xWeights / exp(-x*x)
    end if
    xMaxNotTooSmall = max(x(Nx), xMax)
    allocate(x2(Nx))
    x2=x*x

    allocate(ddx(Nx,Nx))
    allocate(d2dx2(Nx,Nx))
    allocate(ddxPreconditioner(Nx,Nx))
    allocate(ddxToUse(Nx,Nx))
    if (Nx==1) then  !Mono-energetic calculation
       ddx = 0
       d2dx2 = 0
       NxPotentials = 1
    else
       call makeXPolynomialDiffMatrices(x,ddx,d2dx2)
       NxPotentials = ceiling(xMaxNotTooSmall*NxPotentialsPerVth)
    end if

 



    allocate(xPotentials(NxPotentials))
    allocate(xWeightsPotentials(NxPotentials))
    allocate(ddxPotentials(NxPotentials, NxPotentials))
    allocate(d2dx2Potentials(NxPotentials, NxPotentials))
    if (Nx==1) then  !Mono-energetic calculation
       xPotentials = 0
       xWeightsPotentials = 0
       ddxPotentials = 0
       d2dx2Potentials = 0
    else
       call uniformDiffMatrices(NxPotentials, zero, xMaxNotTooSmall, 12, xPotentials, &
            xWeightsPotentials, ddxPotentials, d2dx2Potentials)
    end if
    allocate(regridPolynomialToUniform(NxPotentials, Nx))
    allocate(regridUniformToPolynomial(Nx,NxPotentials))
    if (Nx==1) then  !Mono-energetic calculation
       regridPolynomialToUniform = 0
       regridUniformToPolynomial = 0
    else
       call polynomialInterpolationMatrix(Nx, NxPotentials, x, xPotentials, &
            exp(-x*x), exp(-xPotentials*xPotentials), regridPolynomialToUniform)
       call interpolationMatrix(NxPotentials, Nx, xPotentials, x, regridUniformToPolynomial, -1, 0)
    end if

    ! We need to evaluate the error function on the x grid.
    ! Some Fortran compilers have this function built in.
    ! If not, call the gnu scientific library:
    allocate(erfs(Nx))
    do i=1,Nx
#if ERF_TO_USE == 0
       temp1 = erf(x(i))
#elif ERF_TO_USE == 1
       call erf(x(i), temp1)
#else
       temp1 = x(i) / (exp(-2*x(i)*x(i)) + x(i))
#endif
       erfs(i) = temp1
    end do

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
    ! Depending on RHSMode, specify collisionality and E_r by
    ! one of two methods:
    ! *********************************************************

    if (RHSMode==1) then

       ! Single RHS
       nuPrime = nuN / sqrt(THat) / B0OverBBar * (GHat+iota*IHat)
       EStar = dPhiHatdpsi / iota / sqrt(THat) / psiAHat / B0OverBBar * (omega * GHat)

    else

       ! Compute transport matrix
       nuN = nuPrime * sqrt(THat) * B0OverBBar / (GHat+iota*IHat)
       dPhiHatdpsi = EStar * iota * sqrt(THat) * psiAHat * B0OverBBar / (omega * GHat)
    end if

    ! *********************************************************
    ! *********************************************************
    !
    ! Now build the main matrix, as well as local matrices for 
    ! the left and right boundaries.
    !
    ! *********************************************************
    ! *********************************************************


    ! *********************************************************
    ! Allocate matrices:
    ! *********************************************************

    !predictedNNZForEachRowOfTotalMatrix = 4*(3*Nx + 5*3 + 5*3 + 5 + Nx)
    !predictedNNZForEachRowOfTotalMatrix = 4*(3*Nx + 5*3 + 5*3 + 5 + Nx + Ntheta*Nzeta)
    !predictedNNZForEachRowOfTotalMatrix = 2*(3*Nx + 5*3 + 5*3 + 5 + Nx)
    tempInt1 = 3*Nx + 5*3 + 5*3 + 5 + Nx + Ntheta*Nzeta*Nx
    if (tempInt1 > matrixSize) then
       tempInt1 = matrixSize
    end if
    predictedNNZForEachRowOfTotalMatrix = tempInt1

    predictedNNZForEachRowOfPreconditioner = predictedNNZForEachRowOfTotalMatrix

    allocate(predictedNNZsForEachRow(matrixSize))
    allocate(predictedNNZsForEachRowDiagonal(matrixSize))
    tempInt1 = 5*3 + 5*3 + 3*Nx
    if (tempInt1 > matrixSize) then
       tempInt1 = matrixSize
    end if
    predictedNNZsForEachRow = tempInt1

    select case (constraintScheme)
    case (0)
    case (1)
       ! The rows for the constraints have more nonzeros:
       predictedNNZsForEachRow((Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta*Nx + 1
    case (2)
       ! The rows for the constraints have more nonzeros:
       predictedNNZsForEachRow((Nx*Ntheta*Nzeta*Nxi+1):matrixSize) = Ntheta*Nzeta + 1
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
          ! This method is less thoroghly tested, but it should use much less memory.

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

       ! Sometimes PETSc's sparse direct solver, which is used only when running with 1 proc,
       ! fails with a zero-pivot error (which mumps and superlu_dist do not do.)
       ! To avoid this error, shift the diagonal for the constraint rows:
       if (numProcsInSubComm .eq. 1 .and. whichMatrix==0) then
          temp1 = 1d+0
          do i = Ntheta*Nzeta*Nx*Nxi,matrixSize-1
             call MatSetValue(matrix,i,i,temp1,ADD_VALUES,ierr)
          end do
       end if

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

       ! *********************************************************
       ! Add the streaming d/dtheta term:
       ! *********************************************************

       allocate(thetaPartOfTerm(Ntheta,Ntheta))
       allocate(localThetaPartOfTerm(Ntheta,localNtheta))
       allocate(rowIndices(localNtheta))
       allocate(colIndices(Ntheta))
       do izeta=1,Nzeta
          do itheta=1,Ntheta
             thetaPartOfTerm(itheta,:) = iota*sqrtTHat * ddthetaToUse(itheta,:) / BHat(itheta,izeta)
          end do

          ! PETSc uses the opposite convention to Fortran:
          thetaPartOfTerm = transpose(thetaPartOfTerm)
          localThetaPartOfTerm = thetaPartOfTerm(:,ithetaMin:ithetaMax)

          do ix=1,Nx
             do L=0,(Nxi-1)
                rowIndices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta &
                     + [(j-1, j=ithetaMin,ithetaMax)]*Nzeta + izeta - 1

                ! Super-diagonal term
                if (L < Nxi-1) then
                   ell = L+1
                   colIndices = (ix-1)*Nxi*Ntheta*Nzeta + ell*Ntheta*Nzeta &
                        + [(j-1, j=1,Ntheta)]*Nzeta + izeta - 1

                   call MatSetValuesSparse(matrix, localNtheta, rowIndices, Ntheta, colIndices, &
                        (L+1)/(2*L+three)*x(ix)*localThetaPartOfTerm, ADD_VALUES, ierr)
                end if

                ! Sub-diagonal term
                if (L > 0) then
                   ell = L-1
                   colIndices = (ix-1)*Nxi*Ntheta*Nzeta + ell*Ntheta*Nzeta &
                        + [(j-1, j=1,Ntheta)]*Nzeta + izeta - 1

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
             zetaPartOfTerm(izeta,:) = sqrtTHat * ddzetaToUse(izeta,:) / BHat(itheta,izeta)
          end do

          ! PETSc uses the opposite convention to Fortran:
          zetaPartOfTerm = transpose(zetaPartOfTerm)

          do ix=1,Nx
             do L=0,(Nxi-1)
                rowIndices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + [(i, i=1,Nzeta)] - 1

                ! Super-diagonal term
                if (L < Nxi-1) then
                   colIndices = rowIndices + Ntheta*Nzeta

                   call MatSetValuesSparse(matrix, Nzeta, rowIndices, Nzeta, colIndices, &
                        (L+1)/(2*L+three)*x(ix)*zetaPartOfTerm, ADD_VALUES, ierr)
                end if

                ! Sub-diagonal term
                if (L > 0) then
                   colIndices = rowIndices - Ntheta*Nzeta

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

       factor = omega*GHat/psiAHat*dPhiHatdpsi
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
                rowIndices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta &
                     + [(j-1, j=ithetaMin,ithetaMax)]*Nzeta + izeta - 1

                colIndices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta &
                     + [(j-1, j=1,Ntheta)]*Nzeta + izeta - 1

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

       factor = -omega*IHat/psiAHat*dPhiHatdpsi
       allocate(zetaPartOfTerm(Nzeta,Nzeta))
       allocate(indices(Nzeta))
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
                indices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + [(i, i=1,Nzeta)] - 1
                call MatSetValuesSparse(matrix, Nzeta, indices, Nzeta, indices, &
                     zetaPartOfTerm, ADD_VALUES, ierr)
             end do
          end do
       end do
       deallocate(indices)
       deallocate(zetaPartOfTerm)

       ! *********************************************************
       ! Add the standard mirror term:
       ! *********************************************************

       do itheta=ithetaMin,ithetaMax
          do izeta=1,Nzeta
             factor = -sqrtTHat/(2*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                  * (iota*dBHatdtheta(itheta,izeta) + dBHatdzeta(itheta,izeta))

             do ix=1,Nx
                do L=0,(Nxi-1)
                   rowIndex = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta - 1

                   if (L<Nxi-1) then
                      ! Super-diagonal term:
                      colIndex = rowIndex + Ntheta*Nzeta
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           (L+1)*(L+2)/(2*L+three)*x(ix)*factor, ADD_VALUES, ierr)
                   end if

                   if (L>0) then
                      ! Sub-diagonal term:
                      colIndex = rowIndex - Ntheta*Nzeta
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
                factor = omega*dPhiHatdpsi/(2*psiAHat*(BHat(itheta,izeta)**3)) &
                     * (GHat*dBHatdtheta(itheta,izeta) - IHat* dBHatdzeta(itheta,izeta))

                do ix=1,Nx
                   do L=0,(Nxi-1)
                      rowIndex = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta - 1

                      ! Diagonal term
                      call MatSetValueSparse(matrix, rowIndex, rowIndex, &
                           (L+1)*L/((2*L-one)*(2*L+three))*factor, ADD_VALUES, ierr)

                      if (whichMatrix==1 .or. preconditioner_xi==0) then
                         if (L<Nxi-2) then
                            ! Super-super-diagonal term:
                            colIndex = rowIndex + 2*Ntheta*Nzeta
                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                 (L+3)*(L+2)*(L+1)/((two*L+5)*(2*L+three))*factor, ADD_VALUES, ierr)
                         end if

                         if (L>1) then
                            ! Sub-sub-diagonal term:
                            colIndex = rowIndex - 2*Ntheta*Nzeta
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

          allocate(rowIndices(Nx))
          allocate(colIndices(Nx))
          allocate(xPartOfXDot(Nx,Nx))
          factor = omega/(two*psiAHat)*dPhiHatdpsi

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

                   rowIndices = [(ix-1,ix=1,Nx)]*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta - 1

                   ! Term that is diagonal in L:
                   colIndices = rowIndices
                   LFactor = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*xDotFactor
                   call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                        LFactor*xPartOfXDot, ADD_VALUES, ierr)

                   if (whichMatrix==1 .or. preconditioner_xi==0) then
                      ! Term that is super-super-diagonal in L:
                      if (L<(Nxi-2)) then
                         colIndices = rowIndices + 2*Ntheta*Nzeta
                         LFactor = (L+1)*(L+2)/((two*L+5)*(2*L+3))*xDotFactor
                         call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                              LFactor*xPartOfXDot, ADD_VALUES, ierr)
                      end if

                      ! Term that is sub-sub-diagonal in L:
                      if (L>1) then
                         colIndices = rowIndices - 2*Ntheta*Nzeta
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
                      index = (ix-1)*Nxi*Ntheta*Nzeta + (ixi-1)*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta - 1
                      call MatSetValueSparse(matrix,index,index, &
                           -dPhiHatdpsi*2*omega/(psiAHat*(BHat(itheta,izeta)**3)) &
                           *(GHat*dBHatdtheta(itheta,izeta)- IHat*dBHatdzeta(itheta,izeta)), &
                           ADD_VALUES, ierr)
                   end do
                end do
             end do
          end do
       end if



       ! *********************************************************
       ! Add the collision operator:
       ! *********************************************************

       allocate(x3(Nx))
       allocate(expx2(Nx))
       allocate(Psi_Chandra(Nx))
       allocate(nuD(Nx))
       allocate(PsiPrime(Nx))
       x3 = x*x2
       expx2 = exp(-x2)
       Psi_Chandra = (erfs - 2/sqrtpi*x*expx2) / (2*x2)
       nuD = 3*sqrtpi/4*(erfs - Psi_Chandra) / x3
       PsiPrime = (-erfs + 2/sqrtpi*x*(1+x2)*expx2) / x3

       select case (collisionOperator)
       case (0)
          ! Full linearized Fokker-Planck operator:
          ! For a detailed description of the implementation used here,
          ! see Landreman & Ernst, arXiv:1210.5289 (2012)

          allocate(M21(NxPotentials, Nx))
          allocate(M32(NxPotentials, NxPotentials))
          allocate(M22BackslashM21(NxPotentials, Nx))
          allocate(M33BackslashM32(NxPotentials, NxPotentials))
          allocate(LaplacianTimesX2WithoutL(NxPotentials, NxPotentials))
          allocate(diagonalOfKWithoutThetaPart(Nx))
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

          allocate(xPartOfCECD(Nx,Nx))
          allocate(M12IncludingX0(Nx,NxPotentials))
          allocate(M13IncludingX0(Nx,NxPotentials))
          allocate(M12(Nx,NxPotentials))
          allocate(M13(Nx,NxPotentials))
          allocate(M22(NxPotentials,NxPotentials))
          allocate(M33(NxPotentials,NxPotentials))
          M12IncludingX0 = nuN * 3/(2*pi)*regridUniformToPolynomial
          M13IncludingX0 = -nuN * 3/(2*pi) * matmul(regridUniformToPolynomial, d2dx2Potentials)
          do i=1,Nx
             xPartOfCECD(i,:) = (3*sqrtpi/4)*((Psi_Chandra(i)/x(i))*d2dx2(i,:)  &
                  + (PsiPrime(i)*x(i)  + Psi_Chandra(i) + 2*Psi_Chandra(i)*x2(i))/x2(i) * ddx(i,:))
             xPartOfCECD(i,i) = xPartOfCECD(i,i) + (3*sqrtpi/4)*(2*PsiPrime(i) + 4*Psi_Chandra(i)/x(i)) + 3*expx2(i)
             M12IncludingX0(i,:) = M12IncludingX0(i,:) * expx2(i)
             M13IncludingX0(i,:) = M13IncludingX0(i,:) * x2(i) * expx2(i)
          end do


          allocate(M11(Nx,Nx))
          allocate(KWithoutThetaPart(Nx,Nx))
          allocate(IPIV(NxPotentials))
          allocate(indices(Nx))
          do L=0, Nxi-1
             ! Build M11
             do i=1,Nx
                M11(i,:) = -nuN * xPartOfCECD(i,:)
                M11(i,i) = M11(i,i) - nuN * (-oneHalf*nuD(i)*L*(L+1))
             end do

             if (L < NL) then
                ! Add Rosenbluth potential stuff
                M13 = M13IncludingX0
                M12 = M12IncludingX0

                M22 = LaplacianTimesX2WithoutL + zero
                do i=1,NxPotentials
                   M22(i,i) = M22(i,i) - L*(L+1)
                end do

                ! Add Dirichlet or Neumann boundary condition for potentials at x=0:
                if (L==0) then
                   M22(1,:)=ddxPotentials(1,:)
                else
                   M22(1,:) = 0
                   M22(1,1) = 1
                   M12(:,1) = 0
                   M13(:,1) = 0
                end if
                M33 = M22;

                ! Add Robin boundary condition for potentials at x=xMax:
                M22(NxPotentials,:) = xMaxNotTooSmall*ddxPotentials(NxPotentials,:)
                M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1

                ! My original boundary condition for G:
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

                !KWithoutThetaPart = M11 -  (M12 - M13 * (M33 \ M32)) * (M22 \ M21);
                KWithoutThetaPart = M11 - matmul(M12 - matmul(M13, M33BackslashM32), M22BackslashM21)
             else
                KWithoutThetaPart = M11
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
                            KWithoutThetaPart(i,j) = zero
                         end if
                      end do
                   end do
                case (2)
                   ! Keep only upper-triangular part:
                   do i=2,Nx
                      do j=1,(i-1)
                         KWithoutThetaPart(i,j) = zero
                      end do
                   end do
                case (3)
                   ! Keep only tridiagonal part:
                   do i=1,Nx
                      do j=1,Nx
                         if (abs(i-j)>1) then
                            KWithoutThetaPart(i,j) = zero
                         end if
                      end do
                   end do
                case (4)
                   ! Keep only the diagonal and super-diagonal:
                   do i=1,Nx
                      do j=1,Nx
                         if (i /= j .and. j /= (i+1)) then
                            KWithoutThetaPart(i,j) = zero
                         end if
                      end do
                   end do
                case default
                   print *,"Error! Invalid preconditioner_x"
                   stop
                end select

             end if

             ! PETSc and Fortran use row-major vs column-major:
             KWithoutThetaPart = transpose(KWithoutThetaPart)

             ! Finally, insert values into the main matrix:
             do itheta = ithetaMin,ithetaMax
                do izeta=1,Nzeta
                   indices = [(i, i=0,Nx-1)]*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta - 1
                   call MatSetValuesSparse(matrix, Nx, indices, Nx, indices, &
                        (GHat+iota*IHat)/(BHat(itheta,izeta)**2)*KWithoutThetaPart, ADD_VALUES, ierr)
                end do
             end do
          end do

          deallocate(indices)
          deallocate(M11)
          deallocate(KWithoutThetaPart)
          deallocate(IPIV)

          deallocate(M21)
          deallocate(M32)
          deallocate(M22BackslashM21)
          deallocate(M33BackslashM32)
          deallocate(LaplacianTimesX2WithoutL)
          deallocate(diagonalOfKWithoutThetaPart)

          deallocate(xPartOfCECD)
          deallocate(M12IncludingX0)
          deallocate(M13IncludingX0)
          deallocate(M12)
          deallocate(M13)
          deallocate(M22)
          deallocate(M33)

       case (1,2)
          ! Pitch-angle scattering operator

          do itheta=ithetaMin,ithetaMax
             do izeta=1,Nzeta
                factor = -nuN*(GHat+iota*IHat)/(BHat(itheta,izeta)**2)
                do ix=1,Nx
                   do L=1,(Nxi-1)
                      index = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta-1
                      call MatSetValueSparse(matrix, index, index, -nuD(ix)*oneHalf*L*(L+1)*factor, ADD_VALUES, ierr)
                   end do
                end do
             end do
          end do

          if (collisionOperator==2) then
             ! Add a model momentum-conserving term:

             L=1
             allocate(fieldTerm(Nx,Nx))
             do i=1,Nx
                do j=1,Nx
                   fieldTerm(i,j) = nuD(i)*expx2(i)*x(i)*xWeights(j)*nuD(j)*x3(j)
                end do
             end do
             ! The constant below is int_0^{\infty} dx x^4 \nu_D(x) \exp(-x^2)
             fieldTerm = fieldTerm / 0.354162849836926d+0

             if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                ! we're making the preconditioner, so if desired, simplify the matrix:
                select case (preconditioner_x)
                case (0)
                   ! Do nothing here.
                case (1)
                   ! Zero out everything that is off-diagonal:
                   do i=1,Nx
                      do j=1,Nx
                         if (i /= j) then
                            fieldTerm(i,j)=zero
                         end if
                      end do
                   end do
                case (2)
                   ! Keep only the upper-triangular terms:
                   do i=2,Nx
                      do j=1,(i-1)
                         fieldTerm(i,j) = zero
                      end do
                   end do
                case (3)
                   ! Keep only tridiagonal part:
                   do i=1,Nx
                      do j=1,Nx
                         if (abs(i-j)>1) then
                            fieldTerm(i,j) = zero
                         end if
                      end do
                   end do
                case (4)
                   ! Keep only the diagonal and super-diagonal:
                   do i=1,Nx
                      do j=1,Nx
                         if (i /= j .and. j /= (i+1)) then
                            fieldTerm(i,j) = zero
                         end if
                      end do
                   end do
                case default
                   print *,"Invalid preconditioner_x!"
                   stop
                end select
             end if

             fieldTerm = transpose(fieldTerm) ! PETSc and Fortran use opposite array conventions.

             L=1
             allocate(indices(Nx))
             do itheta=ithetaMin,ithetaMax
                do izeta = 1,Nzeta
                   factor = -nuN*(GHat+iota*IHat)/(BHat(itheta,izeta)**2)
                   indices = [(ix-1,ix=1,Nx)]*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta-1
                   call MatSetValuesSparse(matrix, Nx, indices, Nx, indices, factor*fieldTerm, ADD_VALUES, ierr)
                end do
             end do
             deallocate(indices)
             deallocate(fieldTerm)

          end if

       case default
          print *,"Error! Invalid collisionOperator"
          stop
       end select

       deallocate(x3)
       deallocate(expx2)
       deallocate(Psi_Chandra)
       deallocate(nuD)
       deallocate(PsiPrime)

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
                   rowIndex = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta-1

                   colIndex = matrixSize-2
                   call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource1 / (BHat(itheta,izeta) ** 2),&
                        ADD_VALUES, ierr)

                   colIndex = matrixSize-1
                   call MatSetValueSparse(matrix, rowIndex, colIndex, xPartOfSource2 / (BHat(itheta,izeta) ** 2),&
                        ADD_VALUES, ierr)
                end do
             end do
          end do

       case (2)
          ! Add a L=0 source (which is constant on the flux surface) at each x.
          L=0
          do ix=1,Nx
             do itheta=ithetaMin,ithetaMax
                do izeta = 1,Nzeta
                   rowIndex = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta-1
                   colIndex = Nx*Nxi*Ntheta*Nzeta + ix - 1
                   call MatSetValue(matrix, rowIndex, colIndex, one / (BHat(itheta,izeta) ** 2),&
                        ADD_VALUES, ierr)
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
                      colIndex = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta - 1

                      rowIndex = matrixSize-2
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)

                      rowIndex = matrixSize-1
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           x2(ix)*x2(ix)*xWeights(ix)*factor, ADD_VALUES, ierr)
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
                      colIndex = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta - 1
                      rowIndex = Nx*Nxi*Ntheta*Nzeta + ix - 1
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           factor, ADD_VALUES, ierr)
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
       if (ierr/=0) then
         print *,"[",myCommunicatorIndex,"] after MatAssemblyBegin: ierr=",ierr
       end if
       if ((masterProcInSubComm).and.(ierr/=0)) then
          call PetscTime(time2, ierr2)
          if (whichMatrix==0) then
             print *,"[",myCommunicatorIndex,"] Time to fail to MatAssemblyBegin preconditioner matrices: ", time2-time1, " seconds."
          else
             print *,"[",myCommunicatorIndex,"] Time to fail to MatAssemblyBegin matrices: ", time2-time1, " seconds."
          end if
       end if
       CHKERRQ(ierr)

       call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
       if (ierr/=0) then
         print *,"[",myCommunicatorIndex,"] after MatAssemblyEnd: ierr=",ierr
       end if
       if ((masterProcInSubComm).and.(ierr/=0)) then
          call PetscTime(time2, ierr2)
          if (whichMatrix==0) then
             print *,"[",myCommunicatorIndex,"] Time to fail to MatAssemblyEnd preconditioner matrices: ", time2-time1, " seconds."
          else
             print *,"[",myCommunicatorIndex,"] Time to fail to MatAssemblyEnd matrices: ", time2-time1, " seconds."
          end if
       end if
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
    
    allocate(xPartOfRHS(Nx))

    allocate(densityPerturbation(Ntheta,Nzeta))
    allocate(flow(Ntheta,Nzeta))
    allocate(pressurePerturbation(Ntheta,Nzeta))
    allocate(particleFluxBeforeSurfaceIntegral(Ntheta,Nzeta))
    allocate(momentumFluxBeforeSurfaceIntegral(Ntheta,Nzeta))
    allocate(heatFluxBeforeSurfaceIntegral(Ntheta,Nzeta))
    allocate(NTVBeforeSurfaceIntegral(Ntheta,Nzeta))
    allocate(fNormIsotropicBeforeSurfaceIntegral(Ntheta,Nzeta,Nx))

    allocate(densityIntegralWeights(Nx))
    allocate(flowIntegralWeights(Nx))
    allocate(pressureIntegralWeights(Nx))
    allocate(particleFluxIntegralWeights(Nx))
    allocate(momentumFluxIntegralWeights(Nx))
    allocate(heatFluxIntegralWeights(Nx))
    allocate(NTVIntegralWeights(Nx))
    allocate(fNormIsotropic(Nx))


    ! ***********************************************************************
    ! ***********************************************************************
    ! 
    !  Beginning of the main solver loop:
    !
    ! ***********************************************************************
    ! ***********************************************************************

    select case (RHSMode)
    case (1)
       numRHSs = 1
    case (2)
       numRHSs = 3
    case (3)
       numRHSs = 2
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

       dnHatdpsiToUse = dnHatdpsi
       dTHatdpsiToUse = dTHatdpsi
       EHatToUse = EHat
       dPhiHatdpsiToUse = dPhiHatdpsi
       
       select case (RHSMode)
       case (1)
          ! Single solve: nothing more to do here.

       case (2)
          ! Solve for 3 linearly independent right-hand sides to get the full 3x3 transport matrix:
          dPhiHatdpsiToUse = 0
          select case (whichRHS)
          case (1)
             dnHatdpsiToUse = 1
             dTHatdpsiToUse = 0
             EHatToUse = 0
          case (2)
             ! The next 2 lines ensure (1/n)*dn/dpsi + (3/2)*dT/dpsi = 0 while dT/dpsi is nonzero.
             dnHatdpsiToUse = (3/two)*nHat/THat
             dTHatdpsiToUse = 1
             EHatToUse = 0
          case (3)
             dnHatdpsiToUse = 0
             dTHatdpsiToUse = 0
             EHatToUse = 1
          case default
             print *,"Program should not get here"
             stop
          end select
       case (3)
          ! Solve for 2 linearly independent right-hand sides to get the 2x2 mono-energetic transport coefficients:
          dPhiHatdpsiToUse = 0
          select case (whichRHS)
          case (1)
             dnHatdpsiToUse = 1
             dTHatdpsiToUse = 0
             EHatToUse = 0
          case (2)
             dnHatdpsiToUse = 0
             dTHatdpsiToUse = 0
             EHatToUse = 1
          case default
             print *,"Program should not get here"
             stop
          end select
       end select

       call VecSet(rhs, zero, ierr)

       ! First add the term arising from radial gradients:
       do i=1,Nx
          xPartOfRHS(i) = x2(i)*exp(-x2(i))*( dnHatdpsiToUse/nHat + 2*omega/(Delta*THat)*dPhiHatdpsiToUse &
               + (x2(i) - three/two)*dTHatdpsiToUse/THat)
       end do
       allocate(indices(Nzeta))

       CHKERRQ(ierr)
       do itheta = ithetaMin,ithetaMax
          do ix=1,Nx
             L = 0
             indices = (ix-1)*Ntheta*Nzeta*Nxi + L*Ntheta*Nzeta + (itheta-1)*Nzeta + [(j, j=1,Nzeta )] -1
             call VecSetValues(rhs, Nzeta, indices, &
                  (4/three)/(2 * (BHat(itheta,:)**3) * sqrtTHat) &
                  *xPartOfRHS(ix)*(GHat*dBHatdtheta(itheta,:) - IHat*dBHatdzeta(itheta,:)), &
                  INSERT_VALUES, ierr)

             L = 2
             indices = (ix-1)*Ntheta*Nzeta*Nxi + L*Ntheta*Nzeta + (itheta-1)*Nzeta + [(j, j=1,Nzeta )] -1
             call VecSetValues(rhs, Nzeta, indices, &
                  (2/three) / (2 * (BHat(itheta,:)**3) * sqrtTHat) &
                  *xPartOfRHS(ix)*(GHat*dBHatdtheta(itheta,:) - IHat*dBHatdzeta(itheta,:)), &
                  INSERT_VALUES, ierr)
          end do
       end do
       CHKERRQ(ierr)

       ! Add inductive electric field term:
       L=1
       do i=1,Nx
          xPartOfRHS(i) = x(i) * exp(-x2(i))
       end do
       factor = 2*omega*psiAHat*EHatToUse*(GHat+iota*IHat)/(Delta*Delta*THat*THat*FSABHat2)
       do itheta=ithetaMin,ithetaMax
          do ix=1,Nx
             indices = (ix-1)*Ntheta*Nzeta*Nxi + L*Ntheta*Nzeta + (itheta-1)*Nzeta + [(j, j=1,Nzeta )] -1
             call VecSetValues(rhs, Nzeta, indices, &
                  xPartOfRHS(ix)*factor/BHat(itheta,:), INSERT_VALUES, ierr)
          end do
       end do
       deallocate(indices)

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

!!$    if (layout /= 0) then
!!$       call MatMult(permutationMatrix, soln, tempVec, ierr)
!!$       call VecDestroy(soln, ierr)
!!$       soln = tempVec
!!$    end if

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

          densityIntegralWeights = x*x
          flowIntegralWeights = x*x*x
          pressureIntegralWeights = x*x*x*x
          particleFluxIntegralWeights = x*x*x*x
          momentumFluxIntegralWeights = x*x*x*x*x
          heatFluxIntegralWeights = x*x*x*x*x*x
          NTVIntegralWeights = x*x*x*x

          densityFactor = Delta*4*THat*sqrtTHat/(sqrtpi*psiAHat)
          flowFactor = 4/(three*sqrtpi)*THat*THat
          pressureFactor = Delta*8/(three*sqrtpi*psiAHat)*THat*sqrtTHat
          particleFluxFactor = - (THat ** (5/two))/(sqrtpi)
          momentumFluxFactor = - (THat ** 3)/(sqrtpi)
          heatFluxFactor = - (THat ** (7/two))/(2*sqrtpi)
          NTVFactor =  2 / iota * (THat ** (5/two))/(sqrtpi)
          fNormFactor = Delta*THat*sqrtTHat/psiAHat

          ! Convert the PETSc vector into a normal Fortran array:
          call VecGetArrayF90(solnOnProc0, solnArray, ierr)
          CHKERRQ(ierr)

          if (whichRHS == numRHSs) then
             select case (constraintScheme)
             case (0)
             case (1)
                allocate(sources(2))
                sources = solnArray((matrixSize-1):matrixSize)
             case (2)
                allocate(sources(Nx))
                sources = solnArray((Nx*Nxi*Ntheta*Nzeta+1):matrixSize)
             case default
                print *,"Error! Invalid setting for constraintScheme."
                stop
             end select
          end if

          allocate(indices(Nx))

          L = 0
          do itheta=1,Ntheta
             do izeta=1,Nzeta

                indices = [(ix-1, ix=1,Nx)]*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta

                densityPerturbation(itheta,izeta) = dot_product(xWeights, densityIntegralWeights * solnArray(indices)) &
                     * densityFactor

                pressurePerturbation(itheta,izeta) = dot_product(xWeights, pressureIntegralWeights * solnArray(indices)) &
                     * pressureFactor

                factor = (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))/(BHat(itheta,izeta) ** 3)

                particleFluxBeforeSurfaceIntegral(itheta,izeta) = factor * (8/three) * particleFluxFactor &
                     * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices))

                heatFluxBeforeSurfaceIntegral(itheta,izeta) = factor * (8/three) * heatFluxFactor &
                     * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices))

                fNormIsotropicBeforeSurfaceIntegral(itheta,izeta,1:Nx) = fNormFactor * solnArray(indices)

             end do
          end do

          L = 1
          do itheta=1,Ntheta
             do izeta=1,Nzeta
                indices = [(ix-1, ix=1,Nx)]*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta

                flow(itheta,izeta) = dot_product(xWeights, flowIntegralWeights * solnArray(indices)) &
                     * flowFactor

                factor = (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))/(BHat(itheta,izeta) ** 3)
                momentumFluxBeforeSurfaceIntegral(itheta,izeta) = factor * (16d+0/15) * momentumFluxFactor &
                     * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices))

             end do
          end do

          L = 2
          do itheta=1,Ntheta
             do izeta=1,Nzeta
                indices = [(ix-1, ix=1,Nx)]*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta

                factor = (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))/(BHat(itheta,izeta) ** 3)

                particleFluxBeforeSurfaceIntegral(itheta,izeta) = particleFluxBeforeSurfaceIntegral(itheta,izeta) &
                     + factor * (four/15) * particleFluxFactor &
                     * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices))

                heatFluxBeforeSurfaceIntegral(itheta,izeta) = heatFluxBeforeSurfaceIntegral(itheta,izeta) &
                     + factor * (four/15) * heatFluxFactor &
                     * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices))

                NTVBeforeSurfaceIntegral(itheta,izeta) =  &
                      NTVFactor * NTVKernel(itheta,izeta) &
                      * dot_product(xWeights, NTVIntegralWeights * solnArray(indices))
             end do
          end do

          L = 3
          do itheta=1,Ntheta
             do izeta=1,Nzeta
                indices = [(ix-1, ix=1,Nx)]*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta

                factor = (GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))/(BHat(itheta,izeta) ** 3)

                momentumFluxBeforeSurfaceIntegral(itheta,izeta) = momentumFluxBeforeSurfaceIntegral(itheta,izeta) &
                     + factor * (four/35) * momentumFluxFactor &
                     * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices))

             end do
          end do

          deallocate(indices)

          FSADensityPerturbation=0
          FSAFlow=0
          FSAPressurePerturbation=0
          particleFlux=0
          momentumFlux=0
          heatFlux=0
          NTV=0
          fNormIsotropic=0
          allocate(B2(Ntheta))
          do izeta=1,Nzeta
             B2 = BHat(:,izeta)*BHat(:,izeta)

             FSADensityPerturbation = FSADensityPerturbation + zetaWeights(izeta) &
                  * dot_product(thetaWeights, densityPerturbation(:,izeta)/B2)

             FSAFlow = FSAFlow + zetaWeights(izeta) &
                  * dot_product(thetaWeights, flow(:,izeta)/BHat(:,izeta))

             FSAPressurePerturbation = FSAPressurePerturbation + zetaWeights(izeta) &
                  * dot_product(thetaWeights, pressurePerturbation(:,izeta)/B2)

             particleFlux = particleFlux + zetaWeights(izeta) &
                  * dot_product(thetaWeights, particleFluxBeforeSurfaceIntegral(:,izeta))

             momentumFlux = momentumFlux + zetaWeights(izeta) &
                  * dot_product(thetaWeights, momentumFluxBeforeSurfaceIntegral(:,izeta))

             heatFlux = heatFlux + zetaWeights(izeta) &
                  * dot_product(thetaWeights, heatFluxBeforeSurfaceIntegral(:,izeta))

             NTV = NTV + zetaWeights(izeta) &
                  * dot_product(thetaWeights, NTVBeforeSurfaceIntegral(:,izeta))

          end do
          do ix=1,Nx
             do izeta=1,Nzeta
                fNormIsotropic(ix) = fNormIsotropic(ix) + zetaWeights(izeta) &
                  * dot_product(thetaWeights, fNormIsotropicBeforeSurfaceIntegral(:,izeta,ix))
             end do
          end do
          NTVmulti = NTV * iota * nHat * 2*Delta / (psiAHat * (GHat+iota*IHat) * VPrimeHat)
          deallocate(B2)

          FSADensityPerturbation = FSADensityPerturbation / VPrimeHat
          FSAFlow = FSAFlow / VPrimeHat
          FSAPressurePerturbation = FSAPressurePerturbation / VPrimeHat

          if (RHSMode==2) then
             VPrimeHatWithG = VPrimeHat*(GHat+iota*IHat)
             select case (whichRHS)
             case (1)
                transportMatrix(1,1) = 4*(GHat+iota*IHat)*particleFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*sqrtTHat)*GHat)
                transportMatrix(2,1) = 8*(GHat+iota*IHat)*heatFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*THat*sqrtTHat)*GHat)
                transportMatrix(3,1) = 2*nHat*FSAFlow/(GHat*THat)
                NTVMatrix(1)         = 4*(GHat+iota*IHat)*    NTV     *nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*sqrtTHat)*GHat)
             case (2)
                transportMatrix(1,2) = 4*(GHat+iota*IHat)*particleFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*GHat)
                transportMatrix(2,2) = 8*(GHat+iota*IHat)*heatFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*THat*GHat)
                transportMatrix(3,2) = 2*FSAFlow/(GHat)
                NTVMatrix(2)         = 4*(GHat+iota*IHat)*    NTV     *B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*GHat)
             case (3)
                transportMatrix(1,3) = particleFlux*Delta*Delta*FSABHat2/(VPrimeHatWithG*GHat*psiAHat*omega)
                transportMatrix(2,3) = 2*Delta*Delta*heatFlux*FSABHat2/(GHat*VPrimeHatWithG*psiAHat*THat*omega)
                transportMatrix(3,3) = FSAFlow*Delta*Delta*sqrtTHat*FSABHat2/((GHat+iota*IHat)*2*psiAHat*omega*B0OverBBar)
                NTVMatrix(3)         =     NTV     *Delta*Delta*FSABHat2/(VPrimeHatWithG*GHat*psiAHat*omega)
             end select
          elseif (RHSMode==3) then
             VPrimeHatWithG = VPrimeHat*(GHat+iota*IHat)
             select case (whichRHS)
             case (1)
                transportCoeffs(1,1) = 4*(GHat+iota*IHat)*particleFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat*sqrtTHat)*GHat)
                transportCoeffs(2,1) = 2*nHat*FSAFlow/(GHat*THat)
             case (2)
                transportCoeffs(1,2) = particleFlux*Delta*Delta*FSABHat2/(VPrimeHatWithG*GHat*psiAHat*omega)
                transportCoeffs(2,2) = FSAFlow*Delta*Delta*sqrtTHat*FSABHat2/((GHat+iota*IHat)*2*psiAHat*omega*B0OverBBar)
             end select
          end if

          call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)
          CHKERRQ(ierr)
       end if

    end do

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

