#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscdmdadef.h>
#else
#include <petsc/finclude/petscdmdadef.h>
#endif

  subroutine createGrids()

    use globalVariables
    use polynomialDiffMatrices
    use xGrid
    use petscdmda
    use geometry
    use indices
    use export_f

    implicit none

    PetscErrorCode :: ierr
    integer :: i, j, itheta, izeta, scheme
    real(prec), dimension(:,:), allocatable :: d2dtheta2, d2dzeta2
    real(prec), dimension(:), allocatable :: theta_preconditioner, thetaWeights_preconditioner
    real(prec), dimension(:,:), allocatable :: d2dtheta2_preconditioner
    real(prec), dimension(:), allocatable :: zeta_preconditioner, zetaWeights_preconditioner
    real(prec), dimension(:,:), allocatable :: d2dzeta2_preconditioner
    real(prec), dimension(:), allocatable :: xWeightsPotentials

    real(prec), dimension(:), allocatable :: xWeights_plus1
    real(prec), dimension(:,:), allocatable :: ddx_plus1, d2dx2_plus1
    real(prec), dimension(:,:), allocatable :: interpolateXToXPotentials_plus1, extrapMatrix
    real(prec), dimension(:), allocatable :: x_subset, xWeights_subset
    real(prec), dimension(:,:), allocatable :: ddx_subset, d2dx2_subset


    DM :: myDM
    integer, parameter :: bufferLength = 200
    character(len=bufferLength) :: procAssignments

    integer :: tag, dummy(1)
    integer :: status(MPI_STATUS_SIZE)


    ! *******************************************************************************
    ! Do a few sundry initialization tasks:
    ! *******************************************************************************

    if (constraintScheme < 0) then
       if (collisionOperator == 0) then
          constraintScheme = 1
       else
          constraintScheme = 2
       end if
    end if

    if (masterProc) then
       print *,"---- Numerical parameters: ----"
       print *,"NFourier           = ", NFourier
       print *,"mmax               = ", mmax
       print *,"nmax               = ", nmax
       print *,"Nxi                = ", Nxi
       print *,"NL                 = ", NL
       print *,"Nx                 = ", Nx
       if (xGridScheme<5) then
          print *,"NxPotentialsPerVth = ", NxPotentialsPerVth
          print *,"xMax               = ",xMax
       end if
       print *,"solverTolerance    = ",solverTolerance
       if (useIterativeLinearSolver) then
          print *,"For solving large linear systems, an iterative Krylov solver will be used."
       else
          print *,"For solving large linear systems, a direct solver will be used."
       end if
    end if

    call computeMatrixSize()

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Set up ranges of indices owned by each processor.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    ! Each processor is responsible for building the rows of the matrix and rhs corresponding
    ! to its iFourierMin:iFourierMax, and each processor is resposible for all columns of the matrix.

    ! Assign a range of theta indices to each processor.
    ! This is done by creating a PETSc DM that is not actually used for anything else.
    call DMDACreate1d(MPIComm, DM_BOUNDARY_NONE, NFourier2, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)
    
    call DMDAGetCorners(myDM, iFourierMin, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         localNFourier, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
    
    call DMDestroy(myDM, ierr)

    ! Switch to 1-based indices:
    iFourierMin = iFourierMin + 1
    iFourierMax = iFourierMin+localNFourier-1

    procThatHandlesConstraints = masterProc

    write (procAssignments,fmt="(a,i4,a,i3,a,i3,a,i3,a,i3,a)") "Processor ",myRank," owns Fourier indices ",iFourierMin," to ",iFourierMax

!    call PetscSynchronizedPrintf(MPIComm, procAssignments, ierr)
!    call PetscSynchronizedFlush(MPIComm, ierr)

    ! PETSc's synchronized printing functions seem buggy, so here I've implemented my own version:
    dummy = 0
    tag = 0
    if (masterProc) then
       print *,trim(procAssignments)
       do i = 1,numProcs-1
          ! To avoid a disordered flood of messages to the masterProc,
          ! ping each proc 1 at a time by sending a dummy value:
          call MPI_SEND(dummy,1,MPI_INT,i,tag,MPIComm,ierr)
          ! Now receive the message from proc i:
          call MPI_RECV(procAssignments,bufferLength,MPI_CHAR,i,MPI_ANY_TAG,MPIComm,status,ierr)
          print *,trim(procAssignments)
       end do
    else
       ! First, wait for the dummy message from proc 0:
       call MPI_RECV(dummy,1,MPI_INT,0,MPI_ANY_TAG,MPIComm,status,ierr)
       ! Now send the message to proc 0:
       call MPI_SEND(procAssignments,bufferLength,MPI_CHAR,0,tag,MPIComm,ierr)
    end if




    ! *******************************************************************************
    ! Build x grids, integration weights, and differentiation matrices.
    ! Also build interpolation matrices to map functions from one x grid to the other.
    ! *******************************************************************************

    select case (xGridScheme)
    case (1,2,5,6,7,8)
       ! For these values of xGridScheme, xInterpolationScheme does not matter.
       xInterpolationScheme = -1
    case (3)
       xInterpolationScheme = 1
    case (4)
       xInterpolationScheme = 2
    case default
       print *,"Error! Invalid setting for xGridScheme."
       stop
    end select

    select case (xPotentialsGridScheme)
    case (1)
       xPotentialsInterpolationScheme = 1
    case (2)
       xPotentialsInterpolationScheme = 2
    case (3)
       xPotentialsInterpolationScheme = 1
    case (4)
       xPotentialsInterpolationScheme = 2
    case default
       if (masterProc) then
          print *,"Error! Invalid setting for xPotentialsGridScheme."
       end if
       stop
    end select

    allocate(x(Nx))
    allocate(xWeights(Nx))
    ! The next few arrays/matrices are used only when there is a point at x=0.
    allocate(x_plus1(Nx+1))
    allocate(xWeights_plus1(Nx+1))
    allocate(ddx_plus1(Nx+1,Nx+1))
    allocate(d2dx2_plus1(Nx+1,Nx+1))
    x_plus1 = -1 ! so we know its value if it is not set otherwise.

    if (RHSMode .ne. 3) then
       select case (xGridScheme)
       case (1,5)
          pointAtX0 = .false.
          call makeXGrid(Nx, x, xWeights, .false.)
          xWeights = xWeights / (exp(-x*x)*(x**xGrid_k))

       case (2,6)
          pointAtX0 = .true.
          call makeXGrid(Nx, x, xWeights, .true.)
          xWeights = xWeights / (exp(-x*x)*(x**xGrid_k))

       case (3,4)
          pointAtX0 = .true.

          scheme = 12
          call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
          x_plus1(1)=0 ! For some reason it usually comes out to be 2d-314
          x = x_plus1(1:Nx)
          xWeights = xWeights_plus1(1:Nx)

       case (7)
          pointAtX0 = .true.
          call ChebyshevGrid(Nx+1, zero, xMax, x_plus1, xWeights_plus1, ddx_plus1)
          x_plus1(1)=0 ! Make sure this is exact.
          x = x_plus1(1:Nx)
          xWeights = xWeights_plus1(1:Nx)

          d2dx2_plus1 = matmul(ddx_plus1,ddx_plus1)

       case (8)
          pointAtX0 = .true.
          deallocate(ddx_plus1)
          allocate(ddx_plus1(Nx,Nx))
          call ChebyshevGrid(Nx, zero, xMax, x, xWeights, ddx_plus1)
          x(1)=0 ! Make sure this is exact.

       case default
          print *,"Error! Invalid xGridScheme."
          stop
       end select

    else
       ! Monoenergetic transport matrix calculation.
       x = one
       xWeights = exp(one)

       ! For monoenergetic calculations, we do not want to impose any regularity condition at the first (and only) x index:
       pointAtX0 = .false.
    end if

    xMaxNotTooSmall = max(x(Nx), xMax)
    allocate(x2(Nx))
    x2=x*x
    allocate(expx2(Nx))
    expx2 = exp(-x*x)


    allocate(ddx(Nx,Nx))
    allocate(d2dx2(Nx,Nx))
    allocate(ddx_preconditioner(Nx,Nx))
    if (RHSMode .ne. 3) then

       select case (xGridScheme)
       case (1,2,5,6)
          call makeXPolynomialDiffMatrices(x,ddx,d2dx2)

       case (3,4,7)
          ddx = ddx_plus1(1:Nx, 1:Nx)
          d2dx2 = d2dx2_plus1(1:Nx, 1:Nx)

       case (8)
          ddx = ddx_plus1
          d2dx2 = matmul(ddx,ddx)

       end select

       if (xPotentialsGridScheme==3 .or. xPotentialsGridScheme==4) then
          ! The potentials have an explicit grid point at xMax, whereas the distribution function does not (since f=0 there.)
          NxPotentials = Nx+1
       else
          NxPotentials = ceiling(xMaxNotTooSmall*NxPotentialsPerVth)
       end if
    else
       ! Monoenergetic transport matrix calculation.
       ddx = zero
       d2dx2 = zero
       NxPotentials = 1
    end if

    ! To allow for upwinding in the xDot term associated with Er, set up some other differentiation matrices:
    allocate(ddx_xDot_plus(Nx,Nx))
    allocate(ddx_xDot_preconditioner_plus(Nx,Nx))
    allocate(ddx_xDot_minus(Nx,Nx))
    allocate(ddx_xDot_preconditioner_minus(Nx,Nx))

    select case (xDotDerivativeScheme)
    case (-2)
       ddx_xDot_plus = zero
       ddx_xDot_minus = zero
       allocate(x_subset(Nx-1))
       allocate(ddx_subset(Nx-1,Nx-1))
       allocate(d2dx2_subset(Nx-1,Nx-1))

       x_subset = x(1:Nx-1)
       call makeXPolynomialDiffMatrices(x_subset,ddx_subset,d2dx2_subset)
       ddx_xDot_plus(1:Nx-1,1:Nx-1) = ddx_subset

       x_subset = x(2:Nx)
       call makeXPolynomialDiffMatrices(x_subset,ddx_subset,d2dx2_subset)
       ddx_xDot_minus(2:Nx,2:Nx) = ddx_subset

       deallocate(x_subset,ddx_subset,d2dx2_subset)

    case (-1)
       ddx_xDot_plus = zero
       ddx_xDot_minus = zero
       do i=i,Nx
          allocate(x_subset(i))
          allocate(ddx_subset(i,i))
          allocate(d2dx2_subset(i,i))

          x_subset = x(1:i)
          call makeXPolynomialDiffMatrices(x_subset,ddx_subset,d2dx2_subset)
          ddx_xDot_plus(i,1:i) = ddx_subset(i,:)

          x_subset = x(Nx-i+1:Nx)
          call makeXPolynomialDiffMatrices(x_subset,ddx_subset,d2dx2_subset)
          ddx_xDot_minus(Nx-i+1,Nx-i+1:Nx) = ddx_subset(1,:)

          deallocate(x_subset,ddx_subset,d2dx2_subset)
       end do

    case (0)
       ddx_xDot_plus = ddx
       ddx_xDot_minus = ddx

    case (1)
       scheme = 32
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)

       scheme = 42
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)

    case (2)
       scheme = 52
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)

       scheme = 62
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)

    case (3)
       scheme = 52
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)

       scheme = 62
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
       do i = 2,Nx
          ddx_xDot_minus(Nx,i) =  ddx_xDot_minus(Nx-1,i-1)
       end do

    case (4)
       scheme = 82
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)

       scheme = 92
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)

    case (5)
       scheme = 82
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)
       ! I'm not sure whether these next lines are good or not
       do i = 1,Nx
          ddx_xDot_plus(2,i) =  ddx(2,i)
       end do

       scheme = 92
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
       do i = 2,Nx
          ddx_xDot_minus(Nx,i) =  ddx_xDot_minus(Nx-1,i-1)
       end do

    case (6)
       do i=1,Nx
          do j=1,Nx
             ddx_xDot_plus(i,j) = expx2(i) * ddx(i,j) / expx2(j)
             if (i==j) then
                ddx_xDot_plus(i,j) = ddx_xDot_plus(i,j) - 2*x(i)
             end if
             ddx_xDot_minus(i,j) = ddx_xDot_plus(i,j)
          end do
       end do

    case (7)

       scheme = 82
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)
!!$       ! I'm not sure whether these next lines are good or not
!!$       do i = 1,Nx
!!$          ddx_xDot_plus(2,i) =  ddx(2,i)
!!$       end do

       scheme = 92
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
       do i = 2,Nx
          ddx_xDot_minus(Nx,i) =  ddx_xDot_minus(Nx-1,i-1)
       end do

       do i=1,Nx
          do j=1,Nx
             ddx_xDot_plus(i,j) = expx2(i) * ddx_xDot_plus(i,j) / expx2(j)
             ddx_xDot_minus(i,j) = expx2(i) * ddx_xDot_minus(i,j) / expx2(j)
             if (i==j) then
                ddx_xDot_plus(i,j) = ddx_xDot_plus(i,j) - 2*x(i)
                ddx_xDot_minus(i,j) = ddx_xDot_minus(i,j) - 2*x(i)
             end if
          end do
       end do

    case (8)
       scheme = 102
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)

       scheme = 112
       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
       do i = 3,Nx
          ddx_xDot_minus(Nx,i)     =  ddx_xDot_minus(Nx-2,i-2)
          ddx_xDot_minus(Nx-1,i) =  ddx_xDot_minus(Nx-2,i-1)
       end do

    case (9)
       ! Where trajectories are going into the domain (ddx_xDot_minus), use the standard ddx, in which the first ghost point is set to 0.
       ! Where trajectories are leaving the domain (ddx_xDot_plus), use scheme=12 without setting any ghost points to 0.
       ddx_xDot_minus = ddx
       
       allocate(x_subset(Nx))
       allocate(xWeights_subset(Nx))
       allocate(d2dx2_subset(Nx,Nx))

       scheme = 12
       call uniformDiffMatrices(Nx, zero, x(Nx), scheme, x_subset, x_subset, ddx_xDot_plus, d2dx2_subset)

       deallocate(x_subset,xWeights_subset,d2dx2_subset)

    case (10)
       ! Same as case 9, but switching plus and minus. This should be backwards.
       ddx_xDot_plus = ddx
       
       allocate(x_subset(Nx))
       allocate(xWeights_subset(Nx))
       allocate(d2dx2_subset(Nx,Nx))

       scheme = 12
       call uniformDiffMatrices(Nx, zero, x(Nx), scheme, x_subset, x_subset, ddx_xDot_minus, d2dx2_subset)

       deallocate(x_subset,xWeights_subset,d2dx2_subset)

    case default
       print *,"Error!  Invalid xDotDerivativeScheme"
       stop
    end select

    allocate(xPotentials(NxPotentials))
    allocate(xWeightsPotentials(NxPotentials))
    allocate(ddxPotentials(NxPotentials, NxPotentials))
    allocate(d2dx2Potentials(NxPotentials, NxPotentials))
    if (RHSMode .ne. 3) then
       call uniformDiffMatrices(NxPotentials, zero, xMaxNotTooSmall, 12, xPotentials, &
            xWeightsPotentials, ddxPotentials, d2dx2Potentials)
    else
       xPotentials = 0
       xWeightsPotentials = 0
       ddxPotentials = 0
       d2dx2Potentials = 0
    end if
    maxxPotentials = xPotentials(NxPotentials)

    deallocate(xWeightsPotentials)

    ! Create matrix to interpolate from the distribution-function grid to the Rosenbluth-potential grid:
    allocate(interpolateXToXPotentials(NxPotentials, Nx))
    if (RHSMode .ne. 3) then
       select case (xGridScheme)
       case (1,2,5,6)
          call polynomialInterpolationMatrix(Nx, NxPotentials, x, xPotentials, &
               expx2*(x**xGrid_k), exp(-xPotentials*xPotentials)*(xPotentials**xGrid_k), interpolateXToXPotentials)
       case (3,4)
          allocate(extrapMatrix(NxPotentials, Nx+1))
          allocate(interpolateXToXPotentials_plus1(NxPotentials, Nx+1))
          call interpolationMatrix(Nx+1, NxPotentials, x_plus1, xPotentials, &
               xInterpolationScheme, interpolateXToXPotentials_plus1, extrapMatrix)
          interpolateXToXPotentials = interpolateXToXPotentials_plus1(:,1:Nx)
          deallocate(extrapMatrix)
          deallocate(interpolateXToXPotentials_plus1)
       case (7)
          allocate(interpolateXToXPotentials_plus1(NxPotentials, Nx+1))
          call ChebyshevInterpolationMatrix(Nx+1, NxPotentials, x_plus1, xPotentials, interpolateXToXPotentials_plus1)
          interpolateXToXPotentials = interpolateXToXPotentials_plus1(:,1:Nx)
          deallocate(interpolateXToXPotentials_plus1)
       case (8)
          call ChebyshevInterpolationMatrix(Nx, NxPotentials, x, xPotentials, interpolateXToXPotentials)
       end select
    else
       interpolateXToXPotentials = zero
    end if

    ddx_preconditioner = 0
    ddx_xDot_preconditioner_plus = 0
    ddx_xDot_preconditioner_minus = 0
    select case (preconditioner_x)
    case (0)
       ! No simplification in x:
       ddx_preconditioner = ddx
       ddx_xDot_preconditioner_plus = ddx_xDot_plus
       ddx_xDot_preconditioner_minus = ddx_xDot_minus
    case (1)
       ! Keep only diagonal terms in x:
       do i=1,Nx
          ddx_preconditioner(i,i) = ddx(i,i)
          ddx_xDot_preconditioner_plus(i,i) = ddx_xDot_plus(i,i)
          ddx_xDot_preconditioner_minus(i,i) = ddx_xDot_minus(i,i)
       end do
    case (2)
       ! Keep only upper-triangular terms in x:
       do i=1,Nx
          do j=i,Nx
             ddx_preconditioner(i,j) = ddx(i,j)
             ddx_xDot_preconditioner_plus(i,j) = ddx_xDot_plus(i,j)
             ddx_xDot_preconditioner_minus(i,j) = ddx_xDot_minus(i,j)
          end do
       end do
    case (3)
       ! Keep only tridiagonal terms in x:
       do i=1,Nx
          do j=1,Nx
             if (abs(i-j) <= 1) then
                ddx_preconditioner(i,j) = ddx(i,j)
                ddx_xDot_preconditioner_plus(i,j) = ddx_xDot_plus(i,j)
                ddx_xDot_preconditioner_minus(i,j) = ddx_xDot_minus(i,j)
             end if
          end do
       end do
    case (4)
       ! Keep only diagonal and super-diagonal in x:
       do i=1,Nx
          ddx_preconditioner(i,i) = ddx(i,i)
          ddx_xDot_preconditioner_plus(i,i) = ddx_xDot_plus(i,i)
          ddx_xDot_preconditioner_minus(i,i) = ddx_xDot_minus(i,i)
       end do
       do i=1,(Nx-1)
          ddx_preconditioner(i,i+1) = ddx(i,i+1)
          ddx_xDot_preconditioner_plus(i,i+1) = ddx_xDot_plus(i,i+1)
          ddx_xDot_preconditioner_minus(i,i+1) = ddx_xDot_minus(i,i+1)
       end do
    case default
       print *,"Error! Invalid preconditioner_x"
       stop
    end select

    if ((xGridScheme==5 .or. xGridScheme==6) .and. (RHSMode .ne. 3)) then
       allocate(Rosenbluth_H(Nspecies,Nspecies,NL,Nx,Nx))
       allocate(Rosenbluth_dHdxb(Nspecies,Nspecies,NL,Nx,Nx))
       allocate(Rosenbluth_d2Gdxb2(Nspecies,Nspecies,NL,Nx,Nx))
       call computeRosenbluthPotentialResponse(Nx, x, xWeights, Nspecies, mHats, THats, NL, &
         Rosenbluth_H, Rosenbluth_dHdxb, Rosenbluth_d2Gdxb2,.false.)
         !Rosenbluth_H, Rosenbluth_dHdxb, Rosenbluth_d2Gdxb2,masterProc)
    end if

!    if (masterProc) then
    if (.false.) then
       print *,"xGridScheme:",xGridScheme
       print *,"xInterpolationScheme:",xInterpolationScheme
       print *,"xPotentialsGridScheme:",xPotentialsGridScheme
       print *,"xPotentialsInterpolationScheme:",xPotentialsInterpolationScheme
       print *,"NxPotentials:",NxPotentials
       print *,"x:"
       print *,x
       print *,"xWeights:"
       print *,xWeights
       print *,"ddx:"
       do i=1,Nx
          print *,ddx(i,:)
       end do
       print *,"ddx_xDot_plus:"
       do i=1,Nx
          print *,ddx_xDot_plus(i,:)
       end do
       print *,"ddx_xDot_minus:"
       do i=1,Nx
          print *,ddx_xDot_minus(i,:)
       end do
       print *,"ddx_preconditioner:"
       do i=1,Nx
          print *,ddx_preconditioner(i,:)
       end do
       print *,"ddx_xDot_preconditioner_plus:"
       do i=1,Nx
          print *,ddx_xDot_preconditioner_plus(i,:)
       end do
       print *,"ddx_xDot_preconditioner_minus:"
       do i=1,Nx
          print *,ddx_xDot_preconditioner_minus(i,:)
       end do
!!$       print *,"d2dx2:"
!!$       do i=1,Nx
!!$          print *,d2dx2(i,:)
!!$       end do
!!$       print *,"xPotentials:"
!!$       print *,xPotentials
!!$       if (NxPotentials < 20) then
!!$          print *,"ddxPotentials:"
!!$          do i=1,NxPotentials
!!$             print *,ddxPotentials(i,:)
!!$          end do
!!$          print *,"d2dx2Potentials:"
!!$          do i=1,NxPotentials
!!$             print *,d2dx2Potentials(i,:)
!!$          end do
!!$       end if
!!$       print *,"interpolateXToXPotentials:"
!!$       do i=1,NxPotentials
!!$          print *,interpolateXToXPotentials(i,:)
!!$       end do
    end if

    deallocate(xWeights_plus1)
    deallocate(ddx_plus1)
    deallocate(d2dx2_plus1)

    ! *******************************************************************************
    ! Initialize quantities related to the poloidal and toroidal angle coordinates
    ! *******************************************************************************

    ! Set up real-space grids for theta and zeta:
    Ntheta = 2*mmax+3
    Nzeta  = 2*nmax+3
    allocate(theta(Ntheta))
    allocate(zeta(Nzeta))
    theta = [( 2*pi*j/Ntheta, j=0,Ntheta-1 )]
    zeta  = [( 2*pi*j/Nzeta,  j=0,Nzeta-1  )] / Nperiods
    thetaWeight = 2*pi/Ntheta
    zetaWeight  = 2*pi/Nzeta

    ! *******************************************************************************
    ! Evaluate the magnetic field (and its derivatives) on the (theta, zeta) grid.
    ! *******************************************************************************

    allocate(DHat(Ntheta,Nzeta))

    allocate(BHat(Ntheta,Nzeta))
    allocate(BDotCurlB(Ntheta,Nzeta))
    allocate(uHat(Ntheta,Nzeta))
    allocate(dBHatdtheta(Ntheta,Nzeta))
    allocate(dBHatdzeta(Ntheta,Nzeta))
    allocate(dBHatdpsiHat(Ntheta,Nzeta))

    allocate(BHat_sub_psi(Ntheta,Nzeta))
    allocate(dBHat_sub_psi_dtheta(Ntheta,Nzeta))
    allocate(dBHat_sub_psi_dzeta(Ntheta,Nzeta))

    allocate(BHat_sub_theta(Ntheta,Nzeta))
    allocate(dBHat_sub_theta_dpsiHat(Ntheta,Nzeta))
    allocate(dBHat_sub_theta_dzeta(Ntheta,Nzeta))

    allocate(BHat_sub_zeta(Ntheta,Nzeta))
    allocate(dBHat_sub_zeta_dpsiHat(Ntheta,Nzeta))
    allocate(dBHat_sub_zeta_dtheta(Ntheta,Nzeta))

    allocate(BHat_sup_theta(Ntheta,Nzeta))
    allocate(dBHat_sup_theta_dpsiHat(Ntheta,Nzeta))
    allocate(dBHat_sup_theta_dzeta(Ntheta,Nzeta))

    allocate(BHat_sup_zeta(Ntheta,Nzeta))
    allocate(dBHat_sup_zeta_dpsiHat(Ntheta,Nzeta))
    allocate(dBHat_sup_zeta_dtheta(Ntheta,Nzeta))

    allocate(NTVKernel(Ntheta,Nzeta))


    call computeBHat()

    ! *******************************************************************************
    ! More initialization related to the theta and zeta coordinates:
    ! *******************************************************************************

    call chooseFourierModes()

    allocate(ddtheta(NFourier2,NFourier2))
    allocate(ddzeta (NFourier2,NFourier2))
    call FourierDifferentiationMatrices(NFourier, xm, xn, ddtheta, ddzeta)

    ! *********************************************************
    ! Compute a few quantities related to the magnetic field:
    ! *********************************************************

    call computeBIntegrals()

    if (masterProc) then
       print *,"---- Geometry parameters: ----"
       print *,"Geometry scheme = ", geometryScheme
       print *,"psiAHat (Normalized toroidal flux at the last closed flux surface) = ", psiAHat
       print *,"aHat (Radius of the last closed flux surface in units of RHat) = ", aHat
       if (geometryScheme==1) then
          print *,"epsilon_t = ", epsilon_t
          print *,"epsilon_h = ", epsilon_h
          print *,"epsilon_antisymm = ", epsilon_antisymm
       end if
       print *,"GHat (Boozer component multiplying grad zeta) = ", GHat
       print *,"IHat (Boozer component multiplying grad theta) = ", IHat
       print *,"iota (Rotational transform) = ", iota
    end if

    ! *********************************************************
    ! Allocate some arrays that will be used later for output quantities:
    ! *********************************************************

    allocate(FSADensityPerturbation(Nspecies))
    allocate(FSABFlow(Nspecies))
    allocate(FSABVelocityUsingFSADensity(Nspecies))
    allocate(FSABVelocityUsingFSADensityOverB0(Nspecies))
    allocate(FSABVelocityUsingFSADensityOverRootFSAB2(Nspecies))
    allocate(FSAPressurePerturbation(Nspecies))

    allocate(particleFlux_vm0_psiHat(Nspecies))
    allocate(particleFlux_vm_psiHat(Nspecies))
    allocate(particleFlux_vE0_psiHat(Nspecies))
    allocate(particleFlux_vE_psiHat(Nspecies))
    allocate(particleFlux_vd1_psiHat(Nspecies))
    allocate(particleFlux_vd_psiHat(Nspecies))
    allocate(particleFlux_vm0_psiN(Nspecies))
    allocate(particleFlux_vm_psiN(Nspecies))
    allocate(particleFlux_vE0_psiN(Nspecies))
    allocate(particleFlux_vE_psiN(Nspecies))
    allocate(particleFlux_vd1_psiN(Nspecies))
    allocate(particleFlux_vd_psiN(Nspecies))
    allocate(particleFlux_vm0_rHat(Nspecies))
    allocate(particleFlux_vm_rHat(Nspecies))
    allocate(particleFlux_vE0_rHat(Nspecies))
    allocate(particleFlux_vE_rHat(Nspecies))
    allocate(particleFlux_vd1_rHat(Nspecies))
    allocate(particleFlux_vd_rHat(Nspecies))
    allocate(particleFlux_vm0_rN(Nspecies))
    allocate(particleFlux_vm_rN(Nspecies))
    allocate(particleFlux_vE0_rN(Nspecies))
    allocate(particleFlux_vE_rN(Nspecies))
    allocate(particleFlux_vd1_rN(Nspecies))
    allocate(particleFlux_vd_rN(Nspecies))

    allocate(momentumFlux_vm0_psiHat(Nspecies))
    allocate(momentumFlux_vm_psiHat(Nspecies))
    allocate(momentumFlux_vE0_psiHat(Nspecies))
    allocate(momentumFlux_vE_psiHat(Nspecies))
    allocate(momentumFlux_vd1_psiHat(Nspecies))
    allocate(momentumFlux_vd_psiHat(Nspecies))
    allocate(momentumFlux_vm0_psiN(Nspecies))
    allocate(momentumFlux_vm_psiN(Nspecies))
    allocate(momentumFlux_vE0_psiN(Nspecies))
    allocate(momentumFlux_vE_psiN(Nspecies))
    allocate(momentumFlux_vd1_psiN(Nspecies))
    allocate(momentumFlux_vd_psiN(Nspecies))
    allocate(momentumFlux_vm0_rHat(Nspecies))
    allocate(momentumFlux_vm_rHat(Nspecies))
    allocate(momentumFlux_vE0_rHat(Nspecies))
    allocate(momentumFlux_vE_rHat(Nspecies))
    allocate(momentumFlux_vd1_rHat(Nspecies))
    allocate(momentumFlux_vd_rHat(Nspecies))
    allocate(momentumFlux_vm0_rN(Nspecies))
    allocate(momentumFlux_vm_rN(Nspecies))
    allocate(momentumFlux_vE0_rN(Nspecies))
    allocate(momentumFlux_vE_rN(Nspecies))
    allocate(momentumFlux_vd1_rN(Nspecies))
    allocate(momentumFlux_vd_rN(Nspecies))

    allocate(heatFlux_vm0_psiHat(Nspecies))
    allocate(heatFlux_vm_psiHat(Nspecies))
    allocate(heatFlux_vE0_psiHat(Nspecies))
    allocate(heatFlux_vE_psiHat(Nspecies))
    allocate(heatFlux_vd1_psiHat(Nspecies))
    allocate(heatFlux_vd_psiHat(Nspecies))
    allocate(heatFlux_vm0_psiN(Nspecies))
    allocate(heatFlux_vm_psiN(Nspecies))
    allocate(heatFlux_vE0_psiN(Nspecies))
    allocate(heatFlux_vE_psiN(Nspecies))
    allocate(heatFlux_vd1_psiN(Nspecies))
    allocate(heatFlux_vd_psiN(Nspecies))
    allocate(heatFlux_vm0_rHat(Nspecies))
    allocate(heatFlux_vm_rHat(Nspecies))
    allocate(heatFlux_vE0_rHat(Nspecies))
    allocate(heatFlux_vE_rHat(Nspecies))
    allocate(heatFlux_vd1_rHat(Nspecies))
    allocate(heatFlux_vd_rHat(Nspecies))
    allocate(heatFlux_vm0_rN(Nspecies))
    allocate(heatFlux_vm_rN(Nspecies))
    allocate(heatFlux_vE0_rN(Nspecies))
    allocate(heatFlux_vE_rN(Nspecies))
    allocate(heatFlux_vd1_rN(Nspecies))
    allocate(heatFlux_vd_rN(Nspecies))

    allocate(heatFlux_withoutPhi1_psiHat(Nspecies))
    allocate(heatFlux_withoutPhi1_psiN(Nspecies))
    allocate(heatFlux_withoutPhi1_rHat(Nspecies))
    allocate(heatFlux_withoutPhi1_rN(Nspecies))

    allocate(NTV(Nspecies)) 

    allocate(densityPerturbation(Nspecies,Ntheta,Nzeta))
    allocate(totalDensity(Nspecies,Ntheta,Nzeta))
    allocate(flow(Nspecies,Ntheta,Nzeta))
    allocate(velocityUsingFSADensity(Nspecies,Ntheta,Nzeta))
    allocate(velocityUsingTotalDensity(Nspecies,Ntheta,Nzeta))
    allocate(MachUsingFSAThermalSpeed(Nspecies,Ntheta,Nzeta))
    allocate(pressurePerturbation(Nspecies,Ntheta,Nzeta))
    allocate(pressureAnisotropy(Nspecies,Ntheta,Nzeta))
    allocate(totalPressure(Nspecies,Ntheta,Nzeta))

    allocate(particleFluxBeforeSurfaceIntegral_vm0(Nspecies,Ntheta,Nzeta))
    allocate(particleFluxBeforeSurfaceIntegral_vm(Nspecies,Ntheta,Nzeta))
    allocate(particleFluxBeforeSurfaceIntegral_vE0(Nspecies,Ntheta,Nzeta))
    allocate(particleFluxBeforeSurfaceIntegral_vE(Nspecies,Ntheta,Nzeta))
    allocate(momentumFluxBeforeSurfaceIntegral_vm0(Nspecies,Ntheta,Nzeta))
    allocate(momentumFluxBeforeSurfaceIntegral_vm(Nspecies,Ntheta,Nzeta))
    allocate(momentumFluxBeforeSurfaceIntegral_vE0(Nspecies,Ntheta,Nzeta))
    allocate(momentumFluxBeforeSurfaceIntegral_vE(Nspecies,Ntheta,Nzeta))
    allocate(heatFluxBeforeSurfaceIntegral_vm0(Nspecies,Ntheta,Nzeta))
    allocate(heatFluxBeforeSurfaceIntegral_vm(Nspecies,Ntheta,Nzeta))
    allocate(heatFluxBeforeSurfaceIntegral_vE0(Nspecies,Ntheta,Nzeta))
    allocate(heatFluxBeforeSurfaceIntegral_vE(Nspecies,Ntheta,Nzeta))
    allocate(NTVBeforeSurfaceIntegral(Nspecies,Ntheta,Nzeta)) 

    allocate(particleFlux_vm_psiHat_vs_x(Nspecies,Nx))
    allocate(heatFlux_vm_psiHat_vs_x(Nspecies,Nx))
    allocate(FSABFlow_vs_x(Nspecies,Nx))

    allocate(jHat(Ntheta,Nzeta))
    allocate(Phi1Hat(Ntheta,Nzeta))
    allocate(dPhi1Hatdtheta(Ntheta,Nzeta))
    allocate(dPhi1Hatdzeta(Ntheta,Nzeta))
    Phi1Hat = zero
    dPhi1Hatdtheta = zero
    dPhi1Hatdzeta = zero

    select case (constraintScheme)
    case (0)
       ! No allocation needed in this case.
    case (1,3,4)
       allocate(sources(Nspecies,2))
    case (2)
       allocate(sources(Nspecies,Nx))
    case default
       print *,"Error! Invalid setting for constraintScheme."
       stop
    end select

    select case (RHSMode)
    case (2)
       transportMatrixSize = 3
       allocate(transportMatrix(transportMatrixSize, transportMatrixSize))
       transportMatrix = 0
    case (3)
       transportMatrixSize = 2
       allocate(transportMatrix(transportMatrixSize, transportMatrixSize))
       transportMatrix = 0
    end select

    ! *******************************************************************************

    if (export_full_f .or. export_delta_f) then
       call setup_grids_for_export_f()
    end if

    if (masterProc) then
       print *,"------------------------------------------------------"
       print *,"Done creating grids."
    end if

  end subroutine createGrids
