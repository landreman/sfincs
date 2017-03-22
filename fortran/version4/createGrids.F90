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
    !use export_f

    implicit none

    PetscErrorCode :: ierr
    integer :: i, j, k, itheta, izeta, ispecies, scheme, L
    PetscScalar, dimension(:,:), allocatable :: d2dtheta2, d2dzeta2, ddxi, d2dxi2
    PetscScalar, dimension(:), allocatable :: xWeightsPotentials

    PetscScalar, dimension(:), allocatable :: xWeights_plus1
    PetscScalar, dimension(:,:), allocatable :: ddx_plus1, d2dx2_plus1
    PetscScalar, dimension(:,:), allocatable :: interpolateXToXPotentials_plus1, extrapMatrix
    PetscScalar, dimension(:), allocatable :: x_subset, xWeights_subset
    PetscScalar, dimension(:,:), allocatable :: ddx_subset, d2dx2_subset
    PetscScalar :: temp, Delta_zeta, v_s
    PetscScalar, dimension(:), allocatable :: xi_to_Legendre
    PetscScalar, dimension(:,:), allocatable :: d2dy2, ddy, ddy_plus, ddy_minus
    PetscScalar, dimension(:), allocatable :: y, y_dummy, yWeights_dummy, yWeights, dxi_dy, d2xi_dy2

    DM :: myDM
    integer, parameter :: bufferLength = 200
    character(len=bufferLength) :: procAssignments

    integer :: tag, dummy(1)
    integer :: status(MPI_STATUS_SIZE)
    logical :: call_uniform_diff_matrices
    integer :: derivative_option_plus, derivative_option_minus, derivative_option, quadrature_option

    PetscScalar :: nonuniform_xi_a = 0.7, nonuniform_xi_b = 0.3 ! b=1-a

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

    if (mod(Ntheta, 2) == 0) then
       Ntheta = Ntheta + 1
    end if
    if (mod(Nzeta, 2) == 0) then
       Nzeta = Nzeta + 1
    end if

    if (masterProc) then
       print *,"---- Numerical parameters: ----"
       print *,"Ntheta             = ", Ntheta
       print *,"Nzeta              = ", Nzeta
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


    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Set up ranges of indices owned by each processor.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    ! Each processor is responsible for building the rows of the matrix and rhs corresponding
    ! to its ithetaMin:ithetaMax and izetaMin:izetaMax, and each processor is resposible for all columns of the matrix.

    ! In principle we could distribute in both theta and zeta at the same time.
    ! However, this would lead to negligible increase in speed, since the bottleneck is
    ! not matrix construction but rather the solve, which is parallelized in a completely different
    ! way (determined internally by superlu_dist or mumps.)
    if (Ntheta > Nzeta) then
       ! Distribute in theta but not in zeta

       ! Assign a range of theta indices to each processor.
       ! This is done by creating a PETSc DM that is not actually used for anything else.
       call DMDACreate1d(MPIComm, DM_BOUNDARY_NONE, Ntheta, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)

       call DMDAGetCorners(myDM, ithetaMin, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
            localNtheta, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

       call DMDestroy(myDM, ierr)

       izetaMin = 0
       izetaMax = Nzeta-1
       localNzeta = Nzeta
    else
       ! Distribute in zeta but not in theta

       ! Assign a range of zeta indices to each processor.
       ! This is done by creating a PETSc DM that is not actually used for anything else.
       call DMDACreate1d(MPIComm, DM_BOUNDARY_NONE, Nzeta, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)

       call DMDAGetCorners(myDM, izetaMin, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
            localNzeta, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

       call DMDestroy(myDM, ierr)
       ithetaMin = 0
       ithetaMax = Ntheta-1
       localNtheta = Ntheta
    end if

    ! Below is some code that breaks up the theta and zeta ranges at the same time.
    ! I'm commented it out because PETSc kept giving an error when the number of
    ! procs was large compared to Ntheta and Nzeta.
!!$    ! Assign a range of theta and zeta indices to each processor.
!!$    ! This is done by creating a PETSc DM that is not actually used for anything else.
!!$    call DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, &
!!$         Ntheta, Nzeta, PETSC_DECIDE, PETSC_DECIDE, 1, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, myDM, ierr)
!!$
!!$    call DMDAGetCorners(myDM, ithetaMin, izetaMin, PETSC_NULL_INTEGER, &
!!$         localNtheta, localNzeta, PETSC_NULL_INTEGER, ierr)
!!$
!!$    call DMDestroy(myDM, ierr)

    ! Switch to 1-based indices:
    ithetaMin = ithetaMin + 1
    ithetaMax = ithetaMin+localNtheta-1
    izetaMin = izetaMin + 1
    izetaMax = izetaMin+localNzeta-1

    procThatHandlesConstraints = masterProc

    write (procAssignments,fmt="(a,i4,a,i3,a,i3,a,i3,a,i3,a)") "Processor ",myRank," owns theta indices ",ithetaMin," to ",ithetaMax,&
         " and zeta indices ",izetaMin," to ",izetaMax

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
    ! *******************************************************************************
    !
    ! Build theta grid, integration weights, and differentiation matrices:
    !
    ! *******************************************************************************
    ! *******************************************************************************

    allocate(theta(Ntheta))
    allocate(thetaWeights(Ntheta))
    allocate(ddtheta_plus(Ntheta,Ntheta))
    allocate(ddtheta_minus(Ntheta,Ntheta))
    allocate(ddtheta_plus_preconditioner(Ntheta,Ntheta))
    allocate(ddtheta_minus_preconditioner(Ntheta,Ntheta))
    allocate(d2dtheta2(Ntheta,Ntheta))

    ! *******************************************************************************
    ! Handle d/dtheta for the main matrix.
    ! *******************************************************************************

    call_uniform_diff_matrices = .true.
    select case (theta_derivative_option)

    case (1)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using Fourier pseudospectral method."
       end if
       derivative_option_plus = 20
       derivative_option_minus = derivative_option_plus

    case (2)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using centered finite differences:"
          print *,"   1 point on either side."
       end if
       derivative_option_plus = 0
       derivative_option_minus = derivative_option_plus

    case (3)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using centered finite differences:"
          print *,"   2 points on either side."
       end if
       derivative_option_plus = 10
       derivative_option_minus = derivative_option_plus

    case (4)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 1 point on the other side."
       end if
       derivative_option_plus  = 30
       derivative_option_minus = 40

    case (5)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 50
       derivative_option_minus = 60

    case (6)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 80
       derivative_option_minus = 90

    case (7)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 100
       derivative_option_minus = 110

    case (8)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   2 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 120
       derivative_option_minus = 130

    case (10)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using partially upwinded finite differences:"
          print *,"   2 point on one side, 3 points on the other side."
          print *,"   upwinding factor:",theta_upwinding_factor
       end if

       call_uniform_diff_matrices = .false.
       quadrature_option = 0
       derivative_option_plus  = 120
       derivative_option_minus = 130
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, ddtheta_plus_preconditioner,  d2dtheta2)
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_minus, quadrature_option, theta, thetaWeights, ddtheta_minus_preconditioner, d2dtheta2)
       derivative_option_plus  = 10
       derivative_option_minus = 10
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, ddtheta_plus,  d2dtheta2)
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_minus, quadrature_option, theta, thetaWeights, ddtheta_minus, d2dtheta2)
       ddtheta_plus  = (1-theta_upwinding_factor) * ddtheta_plus  + theta_upwinding_factor * ddtheta_plus_preconditioner
       ddtheta_minus = (1-theta_upwinding_factor) * ddtheta_minus + theta_upwinding_factor * ddtheta_minus_preconditioner

    case default
       if (masterProc) then
          print *,"Error! Invalid setting for theta_derivative_option:",theta_derivative_option
       end if
       stop
    end select

    if (call_uniform_diff_matrices) then
       quadrature_option = 0
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, ddtheta_plus,  d2dtheta2)
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_minus, quadrature_option, theta, thetaWeights, ddtheta_minus, d2dtheta2)
    end if

    ! *******************************************************************************
    ! Handle d/dtheta for the preconditioner matrix.
    ! *******************************************************************************

    call_uniform_diff_matrices = .true.
    select case (abs(preconditioner_theta_derivative_option))
    case (0)
       if (masterProc) then
          print *,"d/dtheta terms are dropped in the preconditioner."
       end if
       ddtheta_plus_preconditioner = 0
       ddtheta_minus_preconditioner = 0
       call_uniform_diff_matrices = .false.

    case (100)
       if (masterProc) then
          print *,"d/dtheta matrices are the same in the preconditioner."
       end if
       ddtheta_plus_preconditioner  = ddtheta_plus
       ddtheta_minus_preconditioner = ddtheta_minus
       call_uniform_diff_matrices = .false.

    case (1)
       if (masterProc) then
          print *,"Preconditioner d/dtheta derivatives discretized using Fourier pseudospectral method."
       end if
       derivative_option_plus = 20
       derivative_option_minus = derivative_option_plus

    case (2)
       if (masterProc) then
          print *,"Preconditioner d/dtheta derivatives discretized using centered finite differences:"
          print *,"   1 point on either side."
       end if
       derivative_option_plus = 0
       derivative_option_minus = derivative_option_plus

    case (3)
       if (masterProc) then
          print *,"Preconditioner d/dtheta derivatives discretized using centered finite differences:"
          print *,"   2 points on either side."
       end if
       derivative_option_plus = 10
       derivative_option_minus = derivative_option_plus

    case (4)
       if (masterProc) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 1 point on the other side."
       end if
       derivative_option_plus  = 30
       derivative_option_minus = 40

    case (5)
       if (masterProc) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 50
       derivative_option_minus = 60

    case (6)
       if (masterProc) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 80
       derivative_option_minus = 90

    case (7)
       if (masterProc) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 100
       derivative_option_minus = 110

    case (8)
       if (masterProc) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   2 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 120
       derivative_option_minus = 130

    case default
       if (masterProc) then
          print *,"Error! Invalid setting for preconditioner_theta_derivative_option:",preconditioner_theta_derivative_option
       end if
       stop
    end select

    if (call_uniform_diff_matrices) then
       quadrature_option = 0
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, ddtheta_plus_preconditioner,  d2dtheta2)
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_minus, quadrature_option, theta, thetaWeights, ddtheta_minus_preconditioner, d2dtheta2)
    end if

    if (preconditioner_theta_derivative_option<0) then
       if (masterProc) then
          print *,"   But only the diagonal is kept."
       end if
       do j=1,Ntheta
          do k=1,Ntheta
             if (j .ne. k) then
                ddtheta_plus_preconditioner(j,k) = 0
                ddtheta_minus_preconditioner(j,k) = 0
             end if
          end do
       end do
    end if


    ! The following arrays will not be needed:
    deallocate(d2dtheta2)


    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Build zeta grid, integration weights, and differentiation matrices:
    !
    ! *******************************************************************************
    ! *******************************************************************************

    zetaMax = 2*pi/NPeriods

    allocate(zeta(Nzeta))
    allocate(zetaWeights(Nzeta))
    allocate(ddzeta_plus(Nzeta,Nzeta))
    allocate(ddzeta_minus(Nzeta,Nzeta))
    allocate(ddzeta_plus_preconditioner(Nzeta,Nzeta))
    allocate(ddzeta_minus_preconditioner(Nzeta,Nzeta))
    allocate(d2dzeta2(Nzeta,Nzeta))

    if (Nzeta==1) then

       ! *******************************************************************************
       ! Axisymmetry is a special case:
       ! *******************************************************************************
       zeta = 0
       zetaWeights = 2*pi
       ddzeta_plus = 0
       ddzeta_minus = 0
       ddzeta_plus_preconditioner = 0
       ddzeta_minus_preconditioner = 0
       d2dzeta2 = 0 ! d2dzeta2 is not actually used.

    else

       ! *******************************************************************************
       ! Not axisymmetric.
       ! First, handle d/dzeta for the main matrix:
       ! *******************************************************************************

       call_uniform_diff_matrices = .true.

       select case (zeta_derivative_option)

       case (2)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_plus = 0
          derivative_option_minus = derivative_option_plus
          
       case (3)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus = 10
          derivative_option_minus = derivative_option_plus
          
       case (4)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
          end if
          derivative_option_plus  = 30
          derivative_option_minus = 40
          
       case (5)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 50
          derivative_option_minus = 60
          
       case (6)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 80
          derivative_option_minus = 90
          
       case (7)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 100
          derivative_option_minus = 110
          
       case (8)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 120
          derivative_option_minus = 130
          
       case (10)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using partially upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
             print *,"   Upwinding factor:",zeta_upwinding_factor
          end if
          call_uniform_diff_matrices = .false.
          quadrature_option = 0
          derivative_option_plus  = 120
          derivative_option_minus = 130
          call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_plus,  quadrature_option, zeta, zetaWeights, ddzeta_plus_preconditioner,  d2dzeta2)
          call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_minus, quadrature_option, zeta, zetaWeights, ddzeta_minus_preconditioner, d2dzeta2)
          derivative_option_plus  = 10
          derivative_option_minus = 10
          call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_plus,  quadrature_option, zeta, zetaWeights, ddzeta_plus,  d2dzeta2)
          call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_minus, quadrature_option, zeta, zetaWeights, ddzeta_minus, d2dzeta2)
          ddzeta_plus  = (1-zeta_upwinding_factor) * ddzeta_plus  + zeta_upwinding_factor * ddzeta_plus_preconditioner
          ddzeta_minus = (1-zeta_upwinding_factor) * ddzeta_minus + zeta_upwinding_factor * ddzeta_minus_preconditioner

       case default
          if (masterProc) then
             print *,"Error! Invalid setting for zeta_derivative_option:",zeta_derivative_option
          end if
          stop
       end select

       if (call_uniform_diff_matrices) then
          quadrature_option = 0 
          call uniformDiffMatrices(Nzeta, 0, zetaMax, derivative_option_plus,  quadrature_option, zeta, zetaWeights, ddzeta_plus, d2dzeta2)
          call uniformDiffMatrices(Nzeta, 0, zetaMax, derivative_option_minus, quadrature_option, zeta, zetaWeights, ddzeta_minus, d2dzeta2)
       end if

       ! *******************************************************************************
       ! Handle d/dzeta for the preconditioner matrix.
       ! *******************************************************************************

       call_uniform_diff_matrices = .true.
       select case (abs(preconditioner_zeta_derivative_option))
       case (0)
          if (masterProc) then
             print *,"d/dzeta terms are dropped in the preconditioner."
          end if
          ddzeta_plus_preconditioner = 0
          ddzeta_minus_preconditioner = 0
          call_uniform_diff_matrices = .false.
          
       case (100)
          if (masterProc) then
             print *,"d/dzeta matrices are the same in the preconditioner."
          end if
          ddzeta_plus_preconditioner  = ddzeta_plus
          ddzeta_minus_preconditioner = ddzeta_minus
          call_uniform_diff_matrices = .false.          
          
       case (2)
          if (masterProc) then
             print *,"Preconditioner d/dzeta derivatives discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_plus = 0
          derivative_option_minus = derivative_option_plus
          
       case (3)
          if (masterProc) then
             print *,"Preconditioner d/dzeta derivatives discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus = 10
          derivative_option_minus = derivative_option_plus

       case (4)
          if (masterProc) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
          end if
          derivative_option_plus  = 30
          derivative_option_minus = 40
          
       case (5)
          if (masterProc) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 50
          derivative_option_minus = 60
          
       case (6)
          if (masterProc) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 80
          derivative_option_minus = 90
          
       case (7)
          if (masterProc) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 100
          derivative_option_minus = 110
          
       case (8)
          if (masterProc) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 120
          derivative_option_minus = 130
          
       case default
          if (masterProc) then
             print *,"Error! Invalid setting for zeta_derivative_option:",zeta_derivative_option
          end if
          stop
       end select
       
       if (call_uniform_diff_matrices) then
          quadrature_option = 0
          call uniformDiffMatrices(Nzeta, 0, zetaMax, &
               derivative_option_plus,  quadrature_option, zeta, zetaWeights, ddzeta_plus_preconditioner, d2dzeta2)
          call uniformDiffMatrices(Nzeta, 0, zetaMax, &
               derivative_option_minus, quadrature_option, zeta, zetaWeights, ddzeta_minus_preconditioner, d2dzeta2)
       end if

       if (preconditioner_zeta_derivative_option<0) then
          if (masterProc) then
             print *,"   But only the diagonal is kept."
          end if
          do j=1,Nzeta
             do k=1,Nzeta
                if (j .ne. k) then
                   ddzeta_plus_preconditioner(j,k) = 0
                   ddzeta_minus_preconditioner(j,k) = 0
                end if
             end do
          end do
       end if
       
       zetaWeights = zetaWeights * Nperiods
    end if

    ! The following arrays will not be needed:
    deallocate(d2dzeta2)

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Build xi grids, integration weights, and differentiation matrices.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    allocate(xi(Nxi))
    allocate(xiWeights(Nxi))
    allocate(ddxi_plus(Nxi,Nxi))
    allocate(ddxi_minus(Nxi,Nxi))
    allocate(ddxi_plus_preconditioner(Nxi,Nxi))
    allocate(ddxi_minus_preconditioner(Nxi,Nxi))
    allocate(d2dxi2(Nxi,Nxi))
    allocate(ddxi(Nxi,Nxi))
    allocate(pitch_angle_scattering_operator(Nxi,Nxi))
    allocate(pitch_angle_scattering_operator_preconditioner(Nxi,Nxi))

    allocate(y(Nxi))
    allocate(dxi_dy(Nxi))
    allocate(d2xi_dy2(Nxi))
    allocate(y_dummy(Nxi))
    allocate(yWeights(Nxi))
    allocate(yWeights_dummy(Nxi))
    allocate(ddy(Nxi,Nxi))
    allocate(d2dy2(Nxi,Nxi))
    allocate(ddy_plus(Nxi,Nxi))
    allocate(ddy_minus(Nxi,Nxi))

    if (pitch_angle_scattering_option==1 .or. xi_derivative_option==1) then
       if (masterProc) then
          print *,"Since at least one of pitch_angle_scattering_option or xi_derivative_option is 1,"
          print *,"we will use a non-preconditioned Chebyshev grid in xi for both."
       end if
       pitch_angle_scattering_option = 1
       xi_derivative_option = 1
       preconditioner_pitch_angle_scattering_option = 100
       preconditioner_xi_derivative_option = 100

       call ChebyshevGrid(Nxi, -one, one, xi, xiWeights, ddxi)
       ddxi_plus = ddxi
       ddxi_minus = ddxi
       ddxi_plus_preconditioner = ddxi_plus
       ddxi_minus_preconditioner = ddxi_minus
       d2dxi2 = matmul(ddxi,ddxi)
       do j=1,Nxi
          pitch_angle_scattering_operator(j,:) = (1/two)*(1-xi(j)*xi(j))*d2dxi2(j,:) - xi(j)*ddxi(j,:)
       end do
       pitch_angle_scattering_operator_preconditioner = pitch_angle_scattering_operator

    else

       ! *******************************************************************************
       ! Set up uniform y grid and the associated nonuniform xi grid
       ! *******************************************************************************
       
       derivative_option = 2
       call uniformDiffMatrices(Nxi, -one, one, derivative_option, xi_quadrature_option, y, yWeights, ddy, d2dy2)
       do j=1,Nxi
          xi(j) = compute_xi_from_y(y(j))
          dxi_dy(j) = compute_dxi_dy(y(j))
          d2xi_dy2(j) = compute_d2xi_dy2(y(j))
       end do
       xiWeights = dxi_dy * yWeights
       ! Do some validation:
       if (abs(sum(xiWeights)-2) > 1.0e-2) then
          if (masterProc) then
             print *,"Error! xiWeights do not sum to 2!"
             print *,"xiWeights:",xiWeights
             print *,"Sum is",sum(xiWeights)
          end if
          stop
       end if
       if (abs(sum(xiWeights*xi)-0) > 1.0e-12) then
          if (masterProc) then
             print *,"Error! xiWeights*xi does not sum to 0!"
             print *,"xi:",xi
             print *,"xiWeights:",xiWeights
             print *,"Sum is",sum(xiWeights*xi)
          end if
          stop
       end if
       if (abs(sum(xiWeights*xi*xi)-2.0/3) > 1.0e-2) then
          if (masterProc) then
             print *,"Error! xiWeights*xi*xi does not sum to 2/3!"
             print *,"xi:",xi
             print *,"xiWeights:",xiWeights
             print *,"Sum is",sum(xiWeights*xi*xi)
          end if
          stop
       end if

       ! *******************************************************************************
       ! Handle d/dxi for the pitch angle scattering operator in the main matrix.
       ! *******************************************************************************
       
       
       select case (pitch_angle_scattering_option)
          
       case (2)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option = 2
          
       case (3)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option = 12
          
       case default
          if (masterProc) then
             print *,"Error! Invalid setting for pitch_angle_scattering_option:",pitch_angle_scattering_option
          end if
          stop
       end select
       
       call uniformDiffMatrices(Nxi, -one, one, derivative_option, xi_quadrature_option, y_dummy, yWeights_dummy, ddy, d2dy2)
       do j=1,Nxi
          !pitch_angle_scattering_operator(j,:) = (1/two)*(1-xi(j)*xi(j))*d2dxi2(j,:) - xi(j)*ddxi(j,:)
          pitch_angle_scattering_operator(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j)) &
               *(d2dy2(j,:) - d2xi_dy2(j)/dxi_dy(j)*ddy(j,:)) &
               - xi(j)/dxi_dy(j)*ddy(j,:)
       end do

       ! *******************************************************************************
       ! Handle d/dxi for the pitch angle scattering operator in the preconditioner matrix.
       ! *******************************************************************************
       
       call_uniform_diff_matrices = .true.
       select case (abs(preconditioner_pitch_angle_scattering_option))
       case (0)
          if (masterProc) then
             print *,"Pitch angle scattering operator is dropped in the preconditioner."
          end if
          ddy = 0
          d2dy2 = 0
          call_uniform_diff_matrices = .false.
          
       case (100)
          if (masterProc) then
             print *,"Pitch angle scattering operator is the same in the preconditioner."
          end if
          ! ddy and d2dy2 will be carried over from the previous section then.
          call_uniform_diff_matrices = .false.
          
       case (2)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator is discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option = 2
          
       case (3)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator is discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option = 12
          
       case default
          if (masterProc) then
             print *,"Error! Invalid setting for preconditioner_pitch_angle_scattering_option:",preconditioner_pitch_angle_scattering_option
          end if
          stop
       end select

       if (call_uniform_diff_matrices) then
          call uniformDiffMatrices(Nxi, -one, one, derivative_option, xi_quadrature_option, y_dummy, yWeights_dummy, ddy, d2dy2)
       end if
       do j=1,Nxi
          !pitch_angle_scattering_operator_preconditioner(j,:) = (1/two)*(1-xi(j)*xi(j))*d2dxi2(j,:) - xi(j)*ddxi(j,:)
          pitch_angle_scattering_operator_preconditioner(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j)) &
               *(d2dy2(j,:) - d2xi_dy2(j)/dxi_dy(j)*ddy(j,:)) &
               - xi(j)/dxi_dy(j)*ddy(j,:)
       end do
       
       if (preconditioner_pitch_angle_scattering_option<0) then
          if (masterProc) then
             print *,"   But only the diagonal is kept."
          end if
          do j=1,Nxi
             do k=1,Nxi
                if (j .ne. k) then
                   pitch_angle_scattering_operator_preconditioner(j,k) = 0
                end if
             end do
          end do
       end if

       ! *******************************************************************************
       ! Handle d/dxi for the mirror term in the main matrix.
       ! *******************************************************************************
       
       select case (xi_derivative_option)
          
       case (2)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_plus = 2
          derivative_option_minus = derivative_option_plus
          
       case (3)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus  = 12
          derivative_option_minus = 12
          
       case (4)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
          end if
          derivative_option_plus  = 32
          derivative_option_minus = 42
          
       case (5)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 52
          derivative_option_minus = 62
          
       case (6)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 82
          derivative_option_minus = 92
          
       case (7)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 102
          derivative_option_minus = 112
          
       case (8)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
             print *,"   High accuracy at the domain ends, though upwinding breaks down there."
          end if
          derivative_option_plus  = 122
          derivative_option_minus = 132
          
       case (9)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
             print *,"   Upwinding all the way to the domain ends, meaning lower accuracy there."
          end if
          derivative_option_plus  = 123
          derivative_option_minus = 133
          
       case default
          if (masterProc) then
             print *,"Error! Invalid setting for xi_derivative_option:",xi_derivative_option
          end if
          stop
       end select

       call uniformDiffMatrices(Nxi, -one, one, derivative_option_plus,  xi_quadrature_option, y_dummy, yWeights_dummy, ddy_plus,  d2dy2)
       call uniformDiffMatrices(Nxi, -one, one, derivative_option_minus, xi_quadrature_option, y_dummy, yWeights_dummy, ddy_minus, d2dy2)
       do j=1,Nxi
          ddxi_plus(j,:)  =  ddy_plus(j,:) / dxi_dy(j)
          ddxi_minus(j,:) = ddy_minus(j,:) / dxi_dy(j)
       end do

       ! *******************************************************************************
       ! Handle d/dxi for the mirror term in the preconditioner matrix.
       ! *******************************************************************************

       call_uniform_diff_matrices = .true.
       select case (abs(preconditioner_xi_derivative_option))
       case (0)
          if (masterProc) then
             print *,"d/dxi terms are dropped in the preconditioner."
          end if
          ddxi_plus_preconditioner = 0
          ddxi_minus_preconditioner = 0
          call_uniform_diff_matrices = .false.

       case (100)
          if (masterProc) then
             print *,"d/dxi matrices are the same in the preconditioner."
          end if
          ddxi_plus_preconditioner  = ddxi_plus
          ddxi_minus_preconditioner = ddxi_minus
          call_uniform_diff_matrices = .false.
          
       case (2)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_plus = 2
          derivative_option_minus = derivative_option_plus
          
       case (3)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus  = 12
          derivative_option_minus = 12
          
       case (4)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
          end if
          derivative_option_plus  = 32
          derivative_option_minus = 42
          
       case (5)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 52
          derivative_option_minus = 62
          
       case (6)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 82
          derivative_option_minus = 92
          
       case (7)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 102
          derivative_option_minus = 112
          
       case (8)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
             print *,"   High accuracy at the domain ends, though upwinding breaks down there."
          end if
          derivative_option_plus  = 122
          derivative_option_minus = 132
          
       case (9)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
             print *,"   Upwinding all the way to the domain ends, meaning lower accuracy there."
          end if
          derivative_option_plus  = 123
          derivative_option_minus = 133
          
       case default
          if (masterProc) then
             print *,"Error! Invalid setting for xi_derivative_option:",xi_derivative_option
          end if
          stop
       end select
       
       if (call_uniform_diff_matrices) then
          call uniformDiffMatrices(Nxi, -one, one, derivative_option_plus,  xi_quadrature_option, y_dummy, yWeights_dummy, ddy_plus,  d2dy2)
          call uniformDiffMatrices(Nxi, -one, one, derivative_option_minus, xi_quadrature_option, y_dummy, yWeights_dummy, ddy_minus, d2dy2)
          do j=1,Nxi
             ddxi_plus_preconditioner( j,:) = ddy_plus( j,:) / dxi_dy(j)
             ddxi_minus_preconditioner(j,:) = ddy_minus(j,:) / dxi_dy(j)
          end do
       end if
       
       if (preconditioner_xi_derivative_option<0) then
          if (masterProc) then
             print *,"   But only the diagonal is kept."
          end if
          do j=1,Nxi
             do k=1,Nxi
                if (j .ne. k) then
                   ddxi_plus_preconditioner(j,k) = 0
                   ddxi_minus_preconditioner(j,k) = 0
                end if
             end do
          end do
       end if
    end if

    ! The following arrays will not be needed:
    deallocate(d2dxi2,ddxi)
    deallocate(d2dy2,ddy,y_dummy,yWeights_dummy,yWeights,ddy_plus,ddy_minus,y)

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

!!$    ! To allow for upwinding in the xDot term associated with Er, set up some other differentiation matrices:
!!$    allocate(ddx_xDot_plus(Nx,Nx))
!!$    allocate(ddx_xDot_preconditioner_plus(Nx,Nx))
!!$    allocate(ddx_xDot_minus(Nx,Nx))
!!$    allocate(ddx_xDot_preconditioner_minus(Nx,Nx))
!!$
!!$    select case (xDotDerivativeScheme)
!!$    case (-2)
!!$       ddx_xDot_plus = zero
!!$       ddx_xDot_minus = zero
!!$       allocate(x_subset(Nx-1))
!!$       allocate(ddx_subset(Nx-1,Nx-1))
!!$       allocate(d2dx2_subset(Nx-1,Nx-1))
!!$
!!$       x_subset = x(1:Nx-1)
!!$       call makeXPolynomialDiffMatrices(x_subset,ddx_subset,d2dx2_subset)
!!$       ddx_xDot_plus(1:Nx-1,1:Nx-1) = ddx_subset
!!$
!!$       x_subset = x(2:Nx)
!!$       call makeXPolynomialDiffMatrices(x_subset,ddx_subset,d2dx2_subset)
!!$       ddx_xDot_minus(2:Nx,2:Nx) = ddx_subset
!!$
!!$       deallocate(x_subset,ddx_subset,d2dx2_subset)
!!$
!!$    case (-1)
!!$       ddx_xDot_plus = zero
!!$       ddx_xDot_minus = zero
!!$       do i=i,Nx
!!$          allocate(x_subset(i))
!!$          allocate(ddx_subset(i,i))
!!$          allocate(d2dx2_subset(i,i))
!!$
!!$          x_subset = x(1:i)
!!$          call makeXPolynomialDiffMatrices(x_subset,ddx_subset,d2dx2_subset)
!!$          ddx_xDot_plus(i,1:i) = ddx_subset(i,:)
!!$
!!$          x_subset = x(Nx-i+1:Nx)
!!$          call makeXPolynomialDiffMatrices(x_subset,ddx_subset,d2dx2_subset)
!!$          ddx_xDot_minus(Nx-i+1,Nx-i+1:Nx) = ddx_subset(1,:)
!!$
!!$          deallocate(x_subset,ddx_subset,d2dx2_subset)
!!$       end do
!!$
!!$    case (0)
!!$       ddx_xDot_plus = ddx
!!$       ddx_xDot_minus = ddx
!!$
!!$    case (1)
!!$       scheme = 32
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)
!!$
!!$       scheme = 42
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
!!$
!!$    case (2)
!!$       scheme = 52
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)
!!$
!!$       scheme = 62
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
!!$
!!$    case (3)
!!$       scheme = 52
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)
!!$
!!$       scheme = 62
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
!!$       do i = 2,Nx
!!$          ddx_xDot_minus(Nx,i) =  ddx_xDot_minus(Nx-1,i-1)
!!$       end do
!!$
!!$    case (4)
!!$       scheme = 82
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)
!!$
!!$       scheme = 92
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
!!$
!!$    case (5)
!!$       scheme = 82
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)
!!$       ! I'm not sure whether these next lines are good or not
!!$       do i = 1,Nx
!!$          ddx_xDot_plus(2,i) =  ddx(2,i)
!!$       end do
!!$
!!$       scheme = 92
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
!!$       do i = 2,Nx
!!$          ddx_xDot_minus(Nx,i) =  ddx_xDot_minus(Nx-1,i-1)
!!$       end do
!!$
!!$    case (6)
!!$       do i=1,Nx
!!$          do j=1,Nx
!!$             ddx_xDot_plus(i,j) = expx2(i) * ddx(i,j) / expx2(j)
!!$             if (i==j) then
!!$                ddx_xDot_plus(i,j) = ddx_xDot_plus(i,j) - 2*x(i)
!!$             end if
!!$             ddx_xDot_minus(i,j) = ddx_xDot_plus(i,j)
!!$          end do
!!$       end do
!!$
!!$    case (7)
!!$
!!$       scheme = 82
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)
!!$!       ! I'm not sure whether these next lines are good or not
!!$!       do i = 1,Nx
!!$!          ddx_xDot_plus(2,i) =  ddx(2,i)
!!$!       end do
!!$
!!$       scheme = 92
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
!!$       do i = 2,Nx
!!$          ddx_xDot_minus(Nx,i) =  ddx_xDot_minus(Nx-1,i-1)
!!$       end do
!!$
!!$       do i=1,Nx
!!$          do j=1,Nx
!!$             ddx_xDot_plus(i,j) = expx2(i) * ddx_xDot_plus(i,j) / expx2(j)
!!$             ddx_xDot_minus(i,j) = expx2(i) * ddx_xDot_minus(i,j) / expx2(j)
!!$             if (i==j) then
!!$                ddx_xDot_plus(i,j) = ddx_xDot_plus(i,j) - 2*x(i)
!!$                ddx_xDot_minus(i,j) = ddx_xDot_minus(i,j) - 2*x(i)
!!$             end if
!!$          end do
!!$       end do
!!$
!!$    case (8)
!!$       scheme = 102
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_plus = ddx_plus1(1:Nx,1:Nx)
!!$
!!$       scheme = 112
!!$       call uniformDiffMatrices(Nx+1, zero, xMax, scheme, x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1)
!!$       ddx_xDot_minus = ddx_plus1(1:Nx,1:Nx)
!!$       do i = 3,Nx
!!$          ddx_xDot_minus(Nx,i)     =  ddx_xDot_minus(Nx-2,i-2)
!!$          ddx_xDot_minus(Nx-1,i) =  ddx_xDot_minus(Nx-2,i-1)
!!$       end do
!!$
!!$    case (9)
!!$       ! Where trajectories are going into the domain (ddx_xDot_minus), use the standard ddx, in which the first ghost point is set to 0.
!!$       ! Where trajectories are leaving the domain (ddx_xDot_plus), use scheme=12 without setting any ghost points to 0.
!!$       ddx_xDot_minus = ddx
!!$       
!!$       allocate(x_subset(Nx))
!!$       allocate(xWeights_subset(Nx))
!!$       allocate(d2dx2_subset(Nx,Nx))
!!$
!!$       scheme = 12
!!$       call uniformDiffMatrices(Nx, zero, x(Nx), scheme, x_subset, x_subset, ddx_xDot_plus, d2dx2_subset)
!!$
!!$       deallocate(x_subset,xWeights_subset,d2dx2_subset)
!!$
!!$    case (10)
!!$       ! Same as case 9, but switching plus and minus. This should be backwards.
!!$       ddx_xDot_plus = ddx
!!$       
!!$       allocate(x_subset(Nx))
!!$       allocate(xWeights_subset(Nx))
!!$       allocate(d2dx2_subset(Nx,Nx))
!!$
!!$       scheme = 12
!!$       call uniformDiffMatrices(Nx, zero, x(Nx), scheme, x_subset, x_subset, ddx_xDot_minus, d2dx2_subset)
!!$
!!$       deallocate(x_subset,xWeights_subset,d2dx2_subset)
!!$
!!$    case default
!!$       print *,"Error!  Invalid xDotDerivativeScheme"
!!$       stop
!!$    end select

    allocate(xPotentials(NxPotentials))
    allocate(xWeightsPotentials(NxPotentials))
    allocate(ddxPotentials(NxPotentials, NxPotentials))
    allocate(d2dx2Potentials(NxPotentials, NxPotentials))
    if (RHSMode .ne. 3) then
       quadrature_option = 0
       call uniformDiffMatrices(NxPotentials, zero, xMaxNotTooSmall, 12, quadrature_option, xPotentials, &
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
    !ddx_xDot_preconditioner_plus = 0
    !ddx_xDot_preconditioner_minus = 0
    select case (preconditioner_x)
    case (0)
       ! No simplification in x:
       ddx_preconditioner = ddx
       !ddx_xDot_preconditioner_plus = ddx_xDot_plus
       !ddx_xDot_preconditioner_minus = ddx_xDot_minus
    case (1)
       ! Keep only diagonal terms in x:
       do i=1,Nx
          ddx_preconditioner(i,i) = ddx(i,i)
          !ddx_xDot_preconditioner_plus(i,i) = ddx_xDot_plus(i,i)
          !ddx_xDot_preconditioner_minus(i,i) = ddx_xDot_minus(i,i)
       end do
    case (2)
       ! Keep only upper-triangular terms in x:
       do i=1,Nx
          do j=i,Nx
             ddx_preconditioner(i,j) = ddx(i,j)
             !ddx_xDot_preconditioner_plus(i,j) = ddx_xDot_plus(i,j)
             !ddx_xDot_preconditioner_minus(i,j) = ddx_xDot_minus(i,j)
          end do
       end do
    case (3)
       ! Keep only tridiagonal terms in x:
       do i=1,Nx
          do j=1,Nx
             if (abs(i-j) <= 1) then
                ddx_preconditioner(i,j) = ddx(i,j)
                !ddx_xDot_preconditioner_plus(i,j) = ddx_xDot_plus(i,j)
                !ddx_xDot_preconditioner_minus(i,j) = ddx_xDot_minus(i,j)
             end if
          end do
       end do
    case (4)
       ! Keep only diagonal and super-diagonal in x:
       do i=1,Nx
          ddx_preconditioner(i,i) = ddx(i,i)
          !ddx_xDot_preconditioner_plus(i,i) = ddx_xDot_plus(i,i)
          !ddx_xDot_preconditioner_minus(i,i) = ddx_xDot_minus(i,i)
       end do
       do i=1,(Nx-1)
          ddx_preconditioner(i,i+1) = ddx(i,i+1)
          !ddx_xDot_preconditioner_plus(i,i+1) = ddx_xDot_plus(i,i+1)
          !ddx_xDot_preconditioner_minus(i,i+1) = ddx_xDot_minus(i,i+1)
       end do
    case default
       print *,"Error! Invalid preconditioner_x"
       stop
    end select

    if ((xGridScheme==5 .or. xGridScheme==6) .and. (RHSMode .ne. 3)) then
       allocate(RosenbluthPotentialTerms(Nspecies,Nspecies,NL,Nx,Nx))
       call computeRosenbluthPotentialResponse(Nx, x, xWeights, Nspecies, mHats, THats, nHats, Zs, NL, &
         RosenbluthPotentialTerms,.false.)
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
!!$       print *,"ddx_xDot_plus:"
!!$       do i=1,Nx
!!$          print *,ddx_xDot_plus(i,:)
!!$       end do
!!$       print *,"ddx_xDot_minus:"
!!$       do i=1,Nx
!!$          print *,ddx_xDot_minus(i,:)
!!$       end do
       print *,"ddx_preconditioner:"
       do i=1,Nx
          print *,ddx_preconditioner(i,:)
       end do
!!$       print *,"ddx_xDot_preconditioner_plus:"
!!$       do i=1,Nx
!!$          print *,ddx_xDot_preconditioner_plus(i,:)
!!$       end do
!!$       print *,"ddx_xDot_preconditioner_minus:"
!!$       do i=1,Nx
!!$          print *,ddx_xDot_preconditioner_minus(i,:)
!!$       end do
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
    ! Compute the Legendre polynomials recursively
    ! *******************************************************************************
    
    if (NL>0) then
       allocate(Legendre_polynomials(Nxi,NL))
       Legendre_polynomials = 1 ! This line takes care of the L=0 polynomial.
       if (NL>1) Legendre_polynomials(:,2) = xi
       do L = 1, NL-2
          Legendre_polynomials(:,L+1+1) = ((2*L+1)*xi*Legendre_polynomials(:,L+1) - L*Legendre_polynomials(:,L-1+1)) / (L+one)
       end do

!!$       print *,'Here come Legendre polynomials:'
!!$       do k=1,Nxi
!!$          print *,Legendre_polynomials(k,:)
!!$       end do
!!$       print *,"xiWeights:"
!!$       print *,xiWeights

       allocate(Legendre_projection(Nxi,Nxi,NL))
       allocate(xi_to_Legendre(Nxi))
       do j=1,NL
          L = j-1
          xi_to_Legendre = (2*L+one)/2*Legendre_polynomials(:,j) * xiWeights
          !print *,"For L=",L,", here is xi_to_Legendre:"
          !print *,xi_to_Legendre
          do k=1,Nxi
             Legendre_projection(k,:,j) = xi_to_Legendre * Legendre_polynomials(k,j)
          end do
          !print *,"For L=",L,", here is the Legendre projection matrix:"
          !do k=1,Nxi
          !   print *,Legendre_projection(k,:,j)
          !end do
       end do
       deallocate(xi_to_Legendre)
    end if

    ! *******************************************************************************
    ! Set the number of Legendre modes used for each value of x
    ! *******************************************************************************
    
    allocate(Nxi_for_x(Nx))

    if (masterProc) print *,"Nxi_for_x_option:",Nxi_for_x_option
    select case (Nxi_for_x_option)
    case (0)
       Nxi_for_x = Nxi
    case (1)
       do j=1,Nx
          ! Linear ramp from 0.1*Nxi to Nxi as x increases from 0 to 2:
          temp = Nxi*(0.1 + 0.9*x(j)/2)
          ! Always keep at least 3 Legendre modes, for the sake of diagnostics.
          ! Always keep at least NL Legendre modes, to simplify the collision operator loops.
          ! Above the threshold value of x, keep exactly Nxi Legendre modes.
          Nxi_for_x(j) = max(3,NL,min(int(temp),Nxi))
       end do
    case (2)
       do j=1,Nx
          ! Quadratic ramp from 0.1*Nxi to Nxi as x increases from 0 to 2:
          temp = Nxi*(0.1 + 0.9*( (x(j)/2)**2) )
          ! Always keep at least 3 Legendre modes, for the sake of diagnostics.
          ! Always keep at least NL Legendre modes, to simplify the collision operator loops.
          ! Above the threshold value of x, keep exactly Nxi Legendre modes.
          Nxi_for_x(j) = max(3,NL,min(int(temp),Nxi))
       end do
    case default
       if (masterProc) print *,"Error! Invalid Nxi_for_x_option"
       stop
    end select

    allocate(min_x_for_L(0:(Nxi-1)))
    min_x_for_L=1
    do j=1,Nx
       min_x_for_L(Nxi_for_x(j):) = j+1
    end do

    if (masterProc) then
       print *,"x:",x
       print *,"Nxi for each x:",Nxi_for_x
       print *,"min_x_for_L:",min_x_for_L
    end if

    call computeMatrixSize()

    ! *******************************************************************************

!!$    if (export_full_f .or. export_delta_f) then
!!$       call setup_grids_for_export_f()
!!$    end if

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Evaluate the magnetic field (and its derivatives) on the (theta, zeta) grid.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    allocate(sqrt_g(Ntheta,Nzeta))

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

    allocate(gradpsidotgradB_overgpsipsi(Ntheta,Nzeta))
    
    allocate(NTVKernel(Ntheta,Nzeta))

    call computeBHat()

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

    allocate(spatial_scaling(Ntheta,Nzeta))
    select case (spatial_scaling_option)
    case (1)
       if (Nzeta==1) then
          spatial_scaling = abs(BHat / BHat_sup_theta)
       else
          spatial_scaling = abs(BHat / BHat_sup_zeta)
       end if
    case (2)
       if (Nzeta==1) then
          spatial_scaling = abs( (theta(2)-theta(1)) * BHat / BHat_sup_theta )
       else
          spatial_scaling = abs( (zeta(2) - zeta(1)) * BHat / BHat_sub_zeta  )
       end if
    case (3)
       if (Nzeta==1) then
          spatial_scaling = abs(BHat / BHat_sup_theta)
       else
          spatial_scaling = abs(BHat / BHat_sup_zeta)
       end if
       spatial_scaling = sum(spatial_scaling)/(Ntheta*Nzeta)
    case (4)
       if (Nzeta==1) then
          spatial_scaling = abs( (theta(2)-theta(1)) * BHat / BHat_sup_theta )
       else
          spatial_scaling = abs( (zeta(2) - zeta(1)) * BHat / BHat_sub_zeta  )
       end if
       spatial_scaling = sum(spatial_scaling)/(Ntheta*Nzeta)
    case default
       if (masterProc) print *,"Error! Invalid spatial_scaling_option:",spatial_scaling_option
       stop
    end select

    if (masterProc) then
       print *,"Here comes spatial_scaling:"
       print *,spatial_scaling
    end if

    allocate(x_scaling(Nx,Nspecies))
    do ispecies = 1,Nspecies
       !v_s = sqrt(2*THats(ispecies)/mHats(ispecies)) ! Once I switch to SI units, include the 2 here.
       v_s = sqrt(THats(ispecies)/mHats(ispecies))    ! But while using the old units, v_s is measured in units of vBar, so there is no 2.

       select case (x_scaling_option)
       case (1)
          x_scaling(:,ispecies) = 1 / (x * v_s)
       case (2)
          x_scaling(:,ispecies) = 1 / v_s
       case default
          if (masterProc) print *,"Error! Invalid x_scaling_option:",x_scaling_option
          stop
       end select
    end do

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

    if (masterProc) then
       print *,"------------------------------------------------------"
       print *,"Done creating grids."
    end if

  contains

    function compute_xi_from_y(yy)
      implicit none
      PetscScalar, intent(in) :: yy
      PetscScalar :: compute_xi_from_y
      
      compute_xi_from_y = nonuniform_xi_a * yy + nonuniform_xi_b * (yy ** 5)
    end function compute_xi_from_y
    
    function compute_dxi_dy(yy)
      implicit none
      PetscScalar, intent(in) :: yy
      PetscScalar :: compute_dxi_dy
      
      compute_dxi_dy = nonuniform_xi_a + 5*nonuniform_xi_b * (yy ** 4)
    end function compute_dxi_dy
    
    function compute_d2xi_dy2(yy)
      implicit none
      PetscScalar, intent(in) :: yy
      PetscScalar :: compute_d2xi_dy2
      
      compute_d2xi_dy2 = 20*nonuniform_xi_b * (yy ** 3)
    end function compute_d2xi_dy2
    
  end subroutine createGrids
