#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscdmdadef.h>
#else
#include <petsc/finclude/petscdmdadef.h>
#endif

  subroutine createGrids()

    use kinds
    use globalVariables
    use polynomialDiffMatrices
    use xGrid
    use petscdmda
    use geometry
    use indices
    !use export_f

    implicit none

    PetscErrorCode :: ierr
    integer :: i, j, k, itheta, izeta, ispecies, ispecies_a, ispecies_b, ix, ixi, scheme, L
    real(prec), dimension(:,:), allocatable :: d2dtheta2, d2dzeta2, ddxi, d2dxi2
    real(prec), dimension(:), allocatable :: xWeightsPotentials

    real(prec), dimension(:), allocatable :: xWeights_plus1
    real(prec), dimension(:,:), allocatable :: ddx_plus1, d2dx2_plus1
    real(prec), dimension(:,:), allocatable :: interpolateXToXPotentials_plus1, extrapMatrix
    real(prec), dimension(:), allocatable :: x_subset, xWeights_subset
    real(prec), dimension(:,:), allocatable :: ddx_subset, d2dx2_subset
    real(prec) :: temp, Delta_zeta, v_s
    real(prec), dimension(:), allocatable :: xi_to_Legendre, temp_array
    real(prec), dimension(:,:), allocatable :: d2dy2, ddy, ddy_plus, ddy_minus, d2dy2_dummy, ddy_dummy, temp_matrix
    real(prec), dimension(:), allocatable :: y, y_dummy, yWeights_dummy, yWeights, dxi_dy, d2xi_dy2
    real(prec) :: temp_plus, temp_minus, factor

    DM :: myDM
    integer, parameter :: bufferLength = 200
    character(len=bufferLength) :: procAssignments

    integer :: tag, dummy(1)
    integer :: status(MPI_STATUS_SIZE)
    logical :: call_uniform_diff_matrices
    integer :: derivative_option_plus, derivative_option_minus, derivative_option, derivative_option_diffusion, quadrature_option
    real(prec) :: dalpha
    real(prec), dimension(:), allocatable :: alpha, alpha_dummy, alphaWeights_dummy
    real(prec), dimension(:,:), allocatable :: ddalpha, ddalpha_extended, d2dalpha2, d2dalpha2_extended, ddalpha_dummy, d2dalpha2_dummy

    real(prec) :: nonuniform_xi_a = 0.7, nonuniform_xi_b = 0.3 ! b=1-a

    ! Variables needed by LAPACK:                                                                                  
    character :: JOBZ
    integer :: INFO, LDA, LDU, LDVT, LWORK, M, N, num_Legendres_to_orthogonalize, iflag
    real(prec), dimension(:,:), allocatable :: U, VT, A
    real(prec), dimension(:), allocatable :: WORK, singular_values
    integer, dimension(:), allocatable :: IWORK

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

    allocate(thermal_speeds(Nspecies))
    do ispecies=1,Nspecies
       thermal_speeds(ispecies) = sqrt(2*THats(ispecies)*electron_charge/(proton_mass*mHats(ispecies)))
    end do

    allocate(collision_frequencies(Nspecies,Nspecies))
    do ispecies_a = 1,Nspecies
       do ispecies_b = 1,Nspecies
          collision_frequencies(ispecies_a,ispecies_b) = 4 * sqrt(2*pi) / 3 &
               * nHats(ispecies_b) * ((Zs(ispecies_a)**2) * (Zs(ispecies_b)**2) * (electron_charge ** 4) * ln_Lambda &
               / (((4*pi*epsilon_0)**2) * ((electron_charge * THats(ispecies_a))**(1.5d+0)) * sqrt(proton_mass*mHats(ispecies_a))))
       end do
    end do

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
    allocate(temp_matrix(Ntheta,Ntheta))
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
          print *,"   2 points on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 120
       derivative_option_minus = 130

    case (11)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using Fromm scheme (upwinding):"
          print *,"   1 point on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 140
       derivative_option_minus = 150

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

    case (12)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using 3rd order upwinded finite differences:"
          print *,"   2 points on either side."
          print *,"   upwinding factor:",upwinding_factor
       end if

       call_uniform_diff_matrices = .false.
       quadrature_option = 0
       derivative_option_plus  = 160
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, temp_matrix,  d2dtheta2)
       ddtheta_plus  = temp_matrix + upwinding_factor * d2dtheta2
       ddtheta_minus = temp_matrix - upwinding_factor * d2dtheta2

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

    case (11)
       if (masterProc) then
          print *,"Preconditioner d/dtheta derivatives discretized using Fromm scheme (upwinding):"
          print *,"   1 point on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 140
       derivative_option_minus = 150

    case (12)
       if (masterProc) then
          print *,"d/dtheta derivatives discretized using 3rd order upwinded finite differences:"
          print *,"   2 points on either side."
          print *,"   upwinding factor:",upwinding_factor
       end if

       call_uniform_diff_matrices = .false.
       quadrature_option = 0
       derivative_option_plus  = 160
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, temp_matrix,  d2dtheta2)
       ddtheta_plus_preconditioner  = temp_matrix + upwinding_factor * d2dtheta2
       ddtheta_minus_preconditioner = temp_matrix - upwinding_factor * d2dtheta2

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
       !temp_plus = maxval(abs(ddtheta_plus_preconditioner))
       !temp_minus = maxval(abs(ddtheta_minus_preconditioner))
       if (masterProc) then
          !print *,"   But only the diagonal is kept."
          !print *,"   But the diagonal is shifted to maintain diagonal dominance."
          print *,"    But blending with the main operator: factor=",preconditioner_theta_blend
       end if
       ddtheta_plus_preconditioner  = (1-preconditioner_theta_blend)*ddtheta_plus_preconditioner  + preconditioner_theta_blend*ddtheta_plus
       ddtheta_minus_preconditioner = (1-preconditioner_theta_blend)*ddtheta_minus_preconditioner + preconditioner_theta_blend*ddtheta_minus
!!$       do j=1,Ntheta
!          do k=1,Ntheta
!             if (j .ne. k) then
!                ddtheta_plus_preconditioner(j,k) = 0
!                ddtheta_minus_preconditioner(j,k) = 0
!             end if
!          end do
!!$          ddtheta_plus_preconditioner(j,j) = temp_plus
!!$          ddtheta_minus_preconditioner(j,j) = -temp_minus
!!$       end do
    end if


    ! The following arrays will not be needed:
    deallocate(d2dtheta2, temp_matrix)


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
    allocate(temp_matrix(Nzeta,Nzeta))
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
          
       case (11)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using Fromm scheme (upwinding):"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 140
          derivative_option_minus = 150

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

       case (12)
          if (masterProc) then
             print *,"d/dzeta derivatives discretized using 3rd order upwinded finite differences:"
             print *,"   2 points on either side."
             print *,"   upwinding factor:",upwinding_factor
          end if
          
          call_uniform_diff_matrices = .false.
          quadrature_option = 0
          derivative_option_plus  = 160
          call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_plus,  quadrature_option, zeta, zetaWeights, temp_matrix,  d2dzeta2)
          ddzeta_plus  = temp_matrix + upwinding_factor * d2dzeta2
          ddzeta_minus = temp_matrix - upwinding_factor * d2dzeta2

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
          
       case (11)
          if (masterProc) then
             print *,"Preconditioner d/dzeta derivatives discretized using Fromm scheme (upwinding):"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 140
          derivative_option_minus = 150

       case (12)
          if (masterProc) then
             print *,"Preconditioner d/dzeta derivatives discretized using 3rd order upwinded finite differences:"
             print *,"   2 points on either side."
             print *,"   upwinding factor:",upwinding_factor
          end if
          
          call_uniform_diff_matrices = .false.
          quadrature_option = 0
          derivative_option_plus  = 160
          call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_plus,  quadrature_option, zeta, zetaWeights, temp_matrix,  d2dzeta2)
          ddzeta_plus_preconditioner  = temp_matrix + upwinding_factor * d2dzeta2
          ddzeta_minus_preconditioner = temp_matrix - upwinding_factor * d2dzeta2

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
          !temp_plus = maxval(abs(ddzeta_plus_preconditioner))
          !temp_minus = maxval(abs(ddzeta_minus_preconditioner))
          if (masterProc) then
             !print *,"   But only the diagonal is kept."
             !print *,"   But the diagonal is shifted to maintain diagonal dominance."
             print *,"   But blending with the main matrix: factor=",preconditioner_zeta_blend
          end if
          ddzeta_plus_preconditioner  = (1-preconditioner_zeta_blend)*ddzeta_plus_preconditioner  + preconditioner_zeta_blend*ddzeta_plus
          ddzeta_minus_preconditioner = (1-preconditioner_zeta_blend)*ddzeta_minus_preconditioner + preconditioner_zeta_blend*ddzeta_minus
!!$          do j=1,Nzeta
!!$! !!$             do k=1,Nzeta
!!$! !!$                if (j .ne. k) then
!!$! !!$                   ddzeta_plus_preconditioner(j,k) = 0
!!$! !!$                   ddzeta_minus_preconditioner(j,k) = 0
!!$! !!$                end if
!!$! !!$             end do
!!$             ddzeta_plus_preconditioner(j,j) = temp_plus
!!$             ddzeta_minus_preconditioner(j,j) = -temp_minus
!!$          end do
       end if
       
       zetaWeights = zetaWeights * Nperiods
    end if

    ! The following arrays will not be needed:
    deallocate(d2dzeta2, temp_matrix)

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Build xi grids, integration weights, and differentiation matrices.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    allocate(xi(Nxi))
    allocate(xiWeights(Nxi))
    allocate(temp_matrix(Nxi,Nxi))
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
    allocate(d2dy2_dummy(Nxi,Nxi))
    allocate(ddy_plus(Nxi,Nxi))
    allocate(ddy_minus(Nxi,Nxi))
    allocate(ddy_dummy(Nxi,Nxi))

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


!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
! Beginning of new xi grid
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************

    elseif (xi_derivative_option >= 1000 .and. xi_derivative_option < 2000) then
       ! Uniform grid in alpha, where xi = -cos(alpha), and f is treated as an odd periodic function in alpha with period 2pi.
       ! Grid points are spaced on the 'half mesh' so there is no point at xi = +/- 1.

       if (preconditioner_xi_derivative_option < 1000 .or. preconditioner_xi_derivative_option >= 2000) &
            stop "Incompatible xi_derivative_option and xi_derivative_option_preconditioner"
       if (pitch_angle_scattering_option < 1000 .or. pitch_angle_scattering_option >= 2000) &
            stop "Incompatible xi_derivative_option and pitch_angle_scattering_option"
       if (preconditioner_pitch_angle_scattering_option < 1000 .or. preconditioner_pitch_angle_scattering_option >= 2000) &
            stop "Incompatible xi_derivative_option and preconditioner_pitch_angle_scattering_option"

       ! Initialize the alpha and xi grids:
       allocate(alpha(Nxi*2))
       alpha = [( (ixi - one)*pi/Nxi, ixi = 1,2*Nxi )]
       ! Shift by half a grid point:
       dalpha = alpha(2) - alpha(1)
       alpha = alpha + dalpha/2
       xi = -cos(alpha(1:Nxi))
       if (masterProc) print *,"alpha:",alpha
       if (masterProc) print *,"xi:",xi

       ! Initialize shifted Clenshaw-Curtis weights:
       xiWeights = dalpha*2/pi ! a_0 term
       do k = 1,(Nxi/2) ! Integer division rounds down.
          xiWeights = xiWeights + dalpha*2/pi*2/(1-4*k*k)*cos(2*k*alpha(1:Nxi)) ! a_{2k} term
       end do
       if (masterProc) print *,"xiWeights:",xiWeights

       ! Do some validation:
       if (masterProc) print *,"sum(xiWeights):",sum(xiWeights)
       if (abs(sum(xiWeights)-2) > 1.0e-12) then
          if (masterProc) then
             print *,"Error! xiWeights do not sum to 2!"
             print *,"xiWeights:",xiWeights
             print *,"Sum is",sum(xiWeights)
          end if
          stop
       end if
       if (masterProc) print *,"sum(xiWeights*xi):",sum(xiWeights*xi)
       if (abs(sum(xiWeights*xi)-0) > 1.0e-12) then
          if (masterProc) then
             print *,"Error! xiWeights*xi does not sum to 0!"
             print *,"xi:",xi
             print *,"xiWeights:",xiWeights
             print *,"Sum is",sum(xiWeights*xi)
          end if
          stop
       end if
       if (masterProc) print *,"sum(xiWeights*xi*xi):",sum(xiWeights*xi*xi)
       if (abs(sum(xiWeights*xi*xi)-two/3) > 1.0e-12) then
          if (masterProc) then
             print *,"Error! xiWeights*xi*xi do not sum to 2/3!"
             print *,"xiWeights:",xiWeights
             print *,"Sum is",sum(xiWeights)
          end if
          stop
       end if

       allocate(ddalpha(Nxi,Nxi))
       allocate(ddalpha_extended(Nxi*2, Nxi*2))
       allocate(d2dalpha2(Nxi,Nxi))
       allocate(d2dalpha2_extended(Nxi*2, Nxi*2))
       allocate(alpha_dummy(Nxi*2))
       allocate(alphaWeights_dummy(Nxi*2))
       allocate(ddalpha_dummy(Nxi*2, Nxi*2))
       allocate(d2dalpha2_dummy(Nxi*2, Nxi*2))

       ! *******************************************************************************
       ! Handle d/dxi for the mirror term and pitch angle scattering in the main matrix.
       ! *******************************************************************************

       call_uniform_diff_matrices = .true.
       select case (xi_derivative_option)
          
       case (1002)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_plus = 0
          derivative_option_minus = derivative_option_plus
          
       case (1003)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus  = 10
          derivative_option_minus = derivative_option_plus
          
       case (1004)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
          end if
          derivative_option_plus  = 30
          derivative_option_minus = 40
          
       case (1005)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 50
          derivative_option_minus = 60
          
       case (1006)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 80
          derivative_option_minus = 90
          
       case (1007)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 100
          derivative_option_minus = 110
          
       case (1008)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 120
          derivative_option_minus = 130
          
       case (1011)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using Fromm scheme (upwinding):"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 140
          derivative_option_minus = 150

       case (1012)
          if (masterProc) then
             print *,"d/dxi derivatives discretized using 3rd order upwinded finite differences:"
             print *,"   2 points on either side."
             print *,"   upwinding factor:",upwinding_factor
          end if
          
          call_uniform_diff_matrices = .false.
          quadrature_option = 0
          derivative_option_plus  = 160
          call uniformDiffMatrices(Nxi*2, zero, two*pi, derivative_option_plus,  quadrature_option, alpha_dummy, alphaWeights_dummy, ddalpha_dummy,  d2dalpha2_dummy)
          ddalpha_extended = ddalpha_dummy + upwinding_factor * d2dalpha2_dummy
          ddalpha = ddalpha_extended(1:Nxi, 1:Nxi) + ddalpha_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))
          do ixi = 1,Nxi
             ddxi_plus(ixi,:) = ddalpha(ixi,:) / sin(alpha(ixi))
          end do

          ddalpha_extended = ddalpha_dummy - upwinding_factor * d2dalpha2_dummy
          ddalpha = ddalpha_extended(1:Nxi, 1:Nxi) + ddalpha_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))
          do ixi = 1,Nxi
             ddxi_minus(ixi,:) = ddalpha(ixi,:) / sin(alpha(ixi))
          end do

       case default
          if (masterProc) then
             print *,"Error! Invalid setting for xi_derivative_option:",xi_derivative_option
          end if
          stop
       end select

       if (masterProc) print *,"   Uniform grid in alpha, using the HALF mesh."

       if (call_uniform_diff_matrices) then
          call uniformDiffMatrices(Nxi*2, zero, 2*pi, derivative_option_plus,  xi_quadrature_option, alpha_dummy, alphaWeights_dummy, ddalpha_extended, d2dalpha2_dummy)
          ! Evaluate boundary terms by extending the alpha domain beyond [0,pi], using f(-alpha) = f(alpha) and f(pi + alpha0) = f(pi - alpha0):
          ddalpha = ddalpha_extended(1:Nxi, 1:Nxi) + ddalpha_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))
          do ixi = 1,Nxi
             ddxi_plus(ixi,:) = ddalpha(ixi,:) / sin(alpha(ixi))
          end do
          
          call uniformDiffMatrices(Nxi*2, zero, 2*pi, derivative_option_minus, xi_quadrature_option, alpha_dummy, alphaWeights_dummy, ddalpha_extended, d2dalpha2_dummy)
          ! Evaluate boundary terms by extending the alpha domain beyond [0,pi], using f(-alpha) = f(alpha) and f(pi + alpha0) = f(pi - alpha0):
          ddalpha = ddalpha_extended(1:Nxi, 1:Nxi) + ddalpha_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))
          do ixi = 1,Nxi
             ddxi_minus(ixi,:) = ddalpha(ixi,:) / sin(alpha(ixi))
          end do
       end if

       ! *******************************************************************************
       ! Handle d/dxi for the mirror term and pitch angle scattering in the preconditioner matrix.
       ! *******************************************************************************

       call_uniform_diff_matrices = .true.
       select case (abs(preconditioner_xi_derivative_option))

       case (1100)
          if (masterProc) then
             print *,"d/dxi matrices are the same in the preconditioner."
          end if
          ddxi_plus_preconditioner  = ddxi_plus
          ddxi_minus_preconditioner = ddxi_minus
          call_uniform_diff_matrices = .false.
          
       case (1002)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_plus = 0
          derivative_option_minus = derivative_option_plus
          
       case (1003)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus  = 10
          derivative_option_minus = derivative_option_plus
          
       case (1004)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
          end if
          derivative_option_plus  = 30
          derivative_option_minus = 40
          
       case (1005)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 50
          derivative_option_minus = 60
          
       case (1006)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 80
          derivative_option_minus = 90
          
       case (1007)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 100
          derivative_option_minus = 110
          
       case (1008)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 120
          derivative_option_minus = 130
          
       case (1011)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using Fromm scheme (upwinding):"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 140
          derivative_option_minus = 150

       case (1012)
          if (masterProc) then
             print *,"Preconditioner d/dxi derivatives discretized using 3rd order upwinded finite differences:"
             print *,"   2 points on either side."
             print *,"   upwinding factor:",upwinding_factor
          end if
          
          call_uniform_diff_matrices = .false.
          quadrature_option = 0
          derivative_option_plus  = 160
          call uniformDiffMatrices(Nxi*2, zero, two*pi, derivative_option_plus,  quadrature_option, alpha_dummy, alphaWeights_dummy, ddalpha_dummy,  d2dalpha2_dummy)
          ddalpha_extended = ddalpha_dummy + upwinding_factor * d2dalpha2_dummy
          ddalpha = ddalpha_extended(1:Nxi, 1:Nxi) + ddalpha_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))
          do ixi = 1,Nxi
             ddxi_plus_preconditioner(ixi,:) = ddalpha(ixi,:) / sin(alpha(ixi))
          end do

          ddalpha_extended = ddalpha_dummy - upwinding_factor * d2dalpha2_dummy
          ddalpha = ddalpha_extended(1:Nxi, 1:Nxi) + ddalpha_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))
          do ixi = 1,Nxi
             ddxi_minus_preconditioner(ixi,:) = ddalpha(ixi,:) / sin(alpha(ixi))
          end do

       case default
          if (masterProc) then
             print *,"Error! Invalid setting for preconditioner_xi_derivative_option:",preconditioner_xi_derivative_option
          end if
          stop
       end select
       
       if (masterProc) print *,"   Uniform grid in alpha, using the HALF mesh."

       if (call_uniform_diff_matrices) then
          call uniformDiffMatrices(Nxi*2, zero, 2*pi, derivative_option_plus,  xi_quadrature_option, alpha_dummy, alphaWeights_dummy, ddalpha_extended, d2dalpha2_dummy)
          ! Evaluate boundary terms by extending the alpha domain beyond [0,pi], using f(-alpha) = f(alpha) and f(pi + alpha0) = f(pi - alpha0):
          ddalpha = ddalpha_extended(1:Nxi, 1:Nxi) + ddalpha_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))
          do ixi = 1,Nxi
             ddxi_plus_preconditioner(ixi,:) = ddalpha(ixi,:) / sin(alpha(ixi))
          end do

          call uniformDiffMatrices(Nxi*2, zero, 2*pi, derivative_option_minus, xi_quadrature_option, alpha_dummy, alphaWeights_dummy, ddalpha_extended, d2dalpha2_dummy)
          ! Evaluate boundary terms by extending the alpha domain beyond [0,pi], using f(-alpha) = f(alpha) and f(pi + alpha0) = f(pi - alpha0):
          ddalpha = ddalpha_extended(1:Nxi, 1:Nxi) + ddalpha_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))
          do ixi = 1,Nxi
             ddxi_minus_preconditioner(ixi,:) = ddalpha(ixi,:) / sin(alpha(ixi))
          end do
       end if


       ! *******************************************************************************
       ! Handle d^2/dxi^2 for the pitch angle scattering operator in the main matrix.
       ! *******************************************************************************

       select case (pitch_angle_scattering_option)
          
       case (1002)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_diffusion = 0
          
       case (1003)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_diffusion = 10
          
       case default
          if (masterProc) then
             print *,"Error! Invalid setting for pitch_angle_scattering_option:",pitch_angle_scattering_option
          end if
          stop
       end select
       
       if (masterProc) print *,"   Uniform grid in alpha, using the HALF mesh."

       call uniformDiffMatrices(Nxi*2, zero, 2*pi, derivative_option_diffusion,  xi_quadrature_option, alpha_dummy, alphaWeights_dummy, ddalpha_dummy, d2dalpha2_extended)
       ! Evaluate boundary terms by extending the alpha domain beyond [0,pi], using f(-alpha) = f(alpha) and f(pi + alpha0) = f(pi - alpha0):
       d2dalpha2 = d2dalpha2_extended(1:Nxi, 1:Nxi) + d2dalpha2_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))

       ! Now use the fact that pitch_angle_scattering = (1/2)*(d2f/dalpha2) - (xi/2)*(df/dxi)
       do j = 1,Nxi
          if (xi(j)>0) then
             pitch_angle_scattering_operator(j,:) = d2dalpha2(j,:)/two - (xi(j)/2)*ddxi_plus(j,:)
          else
             pitch_angle_scattering_operator(j,:) = d2dalpha2(j,:)/two - (xi(j)/2)*ddxi_minus(j,:)
          end if
       end do

       ! *******************************************************************************
       ! Handle d^2/dxi^2 for the pitch angle scattering operator in the preconditioner matrix.
       ! *******************************************************************************
       
       call_uniform_diff_matrices = .true.
       select case (abs(preconditioner_pitch_angle_scattering_option))

       case (1100)
          if (masterProc) then
             print *,"Pitch angle scattering operator is the same in the preconditioner."
          end if
          ! d2dalpha2 will be carried over from the previous section then.
          call_uniform_diff_matrices = .false.
          
       case (1002)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator is discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_diffusion = 0
          
       case (1003)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator is discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_diffusion = 10
          
       case default
          if (masterProc) then
             print *,"Error! Invalid setting for preconditioner_pitch_angle_scattering_option:",preconditioner_pitch_angle_scattering_option
          end if
          stop
       end select

       if (masterProc) print *,"   Uniform grid in alpha, using the HALF mesh."

       if (call_uniform_diff_matrices) then
          call uniformDiffMatrices(Nxi*2, zero, 2*pi, derivative_option_diffusion,  xi_quadrature_option, alpha_dummy, alphaWeights_dummy, ddalpha_dummy, d2dalpha2_extended)
          ! Evaluate boundary terms by extending the alpha domain beyond [0,pi], using f(-alpha) = f(alpha) and f(pi + alpha0) = f(pi - alpha0):
          d2dalpha2 = d2dalpha2_extended(1:Nxi, 1:Nxi) + d2dalpha2_extended(1:Nxi, (2*Nxi):(Nxi+1):(-1))
       end if
       ! Now use the fact that pitch_angle_scattering = (1/2)*(d2f/dalpha2) - (xi/2)*(df/dxi)
       do j = 1,Nxi
          if (xi(j)>0) then
             pitch_angle_scattering_operator_preconditioner(j,:) = d2dalpha2(j,:)/two - (xi(j)/2)*ddxi_plus_preconditioner(j,:)
          else
             pitch_angle_scattering_operator_preconditioner(j,:) = d2dalpha2(j,:)/two - (xi(j)/2)*ddxi_minus_preconditioner(j,:)
          end if
       end do
 

       deallocate(alpha_dummy, alphaWeights_dummy, ddalpha_dummy, d2dalpha2_dummy)
       deallocate(alpha, ddalpha, ddalpha_extended, d2dalpha2, d2dalpha2_extended)

!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
! End of new xi grid
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************

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
!!$       if (abs(sum(xiWeights*xi*xi)-2.0/3) > 1.0e-2) then
!!$          if (masterProc) then
!!$             print *,"Error! xiWeights*xi*xi does not sum to 2/3!"
!!$             print *,"xi:",xi
!!$             print *,"xiWeights:",xiWeights
!!$             print *,"Sum is",sum(xiWeights*xi*xi)
!!$          end if
!!$          stop
!!$       end if

       ! *******************************************************************************
       ! Handle d/dxi for the pitch angle scattering operator in the main matrix.
       ! *******************************************************************************
       
       select case (pitch_angle_scattering_option)
          
       case (2)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_plus = 2
          derivative_option_minus = 2
          derivative_option_diffusion = 2
          
       case (3)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus = 12
          derivative_option_minus = 12
          derivative_option_diffusion = 12
          
       case (4)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
             print *,"   Diffusion term uses centered differences with 1 point on either side."
          end if
          derivative_option_plus  = 32
          derivative_option_minus = 42
          derivative_option_diffusion = 2
          
       case (5)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 52
          derivative_option_minus = 62
          derivative_option_diffusion = 12
          
       case (6)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 82
          derivative_option_minus = 92
          derivative_option_diffusion = 12
          
       case (7)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 102
          derivative_option_minus = 112
          derivative_option_diffusion = 12
          
       case (8)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
             print *,"   High accuracy at the domain ends, though upwinding breaks down there."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 122
          derivative_option_minus = 132
          derivative_option_diffusion = 12
          
       case (9)
          if (masterProc) then
             print *,"Pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
             print *,"   Upwinding all the way to the domain ends, meaning lower accuracy there."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 123
          derivative_option_minus = 133
          derivative_option_diffusion = 12
          
       case default
          if (masterProc) then
             print *,"Error! Invalid setting for pitch_angle_scattering_option:",pitch_angle_scattering_option
          end if
          stop
       end select
       
       call uniformDiffMatrices(Nxi, -one, one, derivative_option_plus,      xi_quadrature_option, y_dummy, yWeights_dummy, ddy_plus, d2dy2_dummy)
       call uniformDiffMatrices(Nxi, -one, one, derivative_option_minus,     xi_quadrature_option, y_dummy, yWeights_dummy, ddy_minus, d2dy2_dummy)
       call uniformDiffMatrices(Nxi, -one, one, derivative_option_diffusion, xi_quadrature_option, y_dummy, yWeights_dummy, ddy_dummy, d2dy2)
       do j=1,Nxi
          !pitch_angle_scattering_operator(j,:) = (1/two)*(1-xi(j)*xi(j))*d2dxi2(j,:) - xi(j)*ddxi(j,:)
          !pitch_angle_scattering_operator(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j)) &
          !     *(d2dy2(j,:) - d2xi_dy2(j)/dxi_dy(j)*ddy(j,:)) &
          !     - xi(j)/dxi_dy(j)*ddy(j,:)
          factor = -(1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j))*d2xi_dy2(j)/dxi_dy(j) - xi(j)/dxi_dy(j)
          if (factor<0) then
             pitch_angle_scattering_operator(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j))*d2dy2(j,:) + factor*ddy_plus(j,:)
          else
             pitch_angle_scattering_operator(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j))*d2dy2(j,:) + factor*ddy_minus(j,:)
          end if
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
          ddy_plus = 0
          ddy_minus = 0
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
          derivative_option_plus = 2
          derivative_option_minus = 2
          derivative_option_diffusion = 2
          
       case (3)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator is discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus = 12
          derivative_option_minus = 12
          derivative_option_diffusion = 12
          
       case (4)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
             print *,"   Diffusion term uses centered differences with 1 point on either side."
          end if
          derivative_option_plus  = 32
          derivative_option_minus = 42
          derivative_option_diffusion = 2
          
       case (5)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 52
          derivative_option_minus = 62
          derivative_option_diffusion = 12
          
       case (6)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 82
          derivative_option_minus = 92
          derivative_option_diffusion = 12
          
       case (7)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 102
          derivative_option_minus = 112
          derivative_option_diffusion = 12
          
       case (8)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
             print *,"   High accuracy at the domain ends, though upwinding breaks down there."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 122
          derivative_option_minus = 132
          derivative_option_diffusion = 12
          
       case (9)
          if (masterProc) then
             print *,"Preconditioner pitch angle scattering operator discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
             print *,"   Upwinding all the way to the domain ends, meaning lower accuracy there."
             print *,"   Diffusion term uses centered differences with 2 points on either side."
          end if
          derivative_option_plus  = 123
          derivative_option_minus = 133
          derivative_option_diffusion = 12

       case default
          if (masterProc) then
             print *,"Error! Invalid setting for preconditioner_pitch_angle_scattering_option:",preconditioner_pitch_angle_scattering_option
          end if
          stop
       end select

       if (call_uniform_diff_matrices) then
          !call uniformDiffMatrices(Nxi, -one, one, derivative_option, xi_quadrature_option, y_dummy, yWeights_dummy, ddy, d2dy2)
          call uniformDiffMatrices(Nxi, -one, one, derivative_option_plus,      xi_quadrature_option, y_dummy, yWeights_dummy, ddy_plus, d2dy2_dummy)
          call uniformDiffMatrices(Nxi, -one, one, derivative_option_minus,     xi_quadrature_option, y_dummy, yWeights_dummy, ddy_minus, d2dy2_dummy)
          call uniformDiffMatrices(Nxi, -one, one, derivative_option_diffusion, xi_quadrature_option, y_dummy, yWeights_dummy, ddy_dummy, d2dy2)
       end if
       do j=1,Nxi
          !pitch_angle_scattering_operator_preconditioner(j,:) = (1/two)*(1-xi(j)*xi(j))*d2dxi2(j,:) - xi(j)*ddxi(j,:)
          !pitch_angle_scattering_operator_preconditioner(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j)) &
          !     *(d2dy2(j,:) - d2xi_dy2(j)/dxi_dy(j)*ddy(j,:)) &
          !     - xi(j)/dxi_dy(j)*ddy(j,:)
          factor = -(1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j))*d2xi_dy2(j)/dxi_dy(j) - xi(j)/dxi_dy(j)
          if (factor<0) then
             pitch_angle_scattering_operator_preconditioner(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j))*d2dy2(j,:) + factor*ddy_plus(j,:)
          else
             pitch_angle_scattering_operator_preconditioner(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j))*d2dy2(j,:) + factor*ddy_minus(j,:)
          end if
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
             print *,"Error! Invalid setting for preconditioner_xi_derivative_option:",preconditioner_xi_derivative_option
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
             !print *,"   But only the diagonal is kept."
             !print *,"   But the diagonal is shifted to maintain diagonal dominance."
             print *,"   But blending with the main matrix: factor=",preconditioner_xi_blend
          end if
          ddxi_plus_preconditioner  = (1-preconditioner_xi_blend)*ddxi_plus_preconditioner  + preconditioner_xi_blend*ddxi_plus
          ddxi_minus_preconditioner = (1-preconditioner_xi_blend)*ddxi_minus_preconditioner + preconditioner_xi_blend*ddxi_minus
!!$          do j=1,Nxi
!!$             do k=1,Nxi
!!$                if (j .ne. k) then
!!$                   ddxi_plus_preconditioner(j,k) = 0
!!$                   ddxi_minus_preconditioner(j,k) = 0
!!$                end if
!!$             end do
!!$          end do

!!$          allocate(temp_array(Nxi))
!!$          do j=1,Nxi
!!$             temp_array(j) = max(maxval(abs(ddxi_plus_preconditioner(j,:))), maxval(abs(ddxi_plus_preconditioner(:,j))))
!!$          end do
!!$          do j=1,Nxi
!!$             ddxi_plus_preconditioner(j,j) = temp_array(j)
!!$             ddxi_minus_preconditioner(j,j) = -temp_array(j)
!!$          end do
!!$          deallocate(temp_array)
       end if
    end if

    ! The following arrays will not be needed:
    deallocate(d2dxi2,ddxi,temp_matrix)
    deallocate(d2dy2,d2dy2_dummy,ddy,y_dummy,yWeights_dummy,yWeights,ddy_plus,ddy_minus,ddy_dummy,y)

    if (masterProc) then
       print *,"pitch_angle_scattering_operator:"
       do j=1,Nxi
          print "(*(f7.2))",pitch_angle_scattering_operator(j,:)
       end do
       print *,"pitch_angle_scattering_operator_preconditioner:"
       do j=1,Nxi
          print "(*(f7.2))",pitch_angle_scattering_operator_preconditioner(j,:)
       end do
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

       if (xi_quadrature_option<0) then
          num_Legendres_to_orthogonalize = NL
          if (num_Legendres_to_orthogonalize>Nxi) stop "Error! num_Legendres_to_orthogonalize>Nxi"
          JOBZ='S' ! Compute only the first min(M,N) singular vectors
          M = Nxi
          N = num_Legendres_to_orthogonalize
          LDA = M
          LDU = M
          LDVT = N
          ! This next formula comes from the LAPACK documentation at the end of the file.
          LWORK = max( 3*min(M,N) + max(max(M,N),7*min(M,N)), &
               3*min(M,N) + max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)), &
               min(M,N)*(6+4*min(M,N))+max(M,N))
          allocate(WORK(LWORK),stat=iflag)
          allocate(IWORK(8*min(M,N)),stat=iflag)
          allocate(singular_values(num_Legendres_to_orthogonalize))
          ! Matrix is destroyed by LAPACK, so make a copy:
          allocate(A(M,N),stat=iflag)
          A = Legendre_polynomials(:,1:num_Legendres_to_orthogonalize)
          allocate(U(M,M),stat=iflag)
          allocate(VT(N,N),stat=iflag)
          ! Call LAPACK to do the SVD:
          call DGESDD(JOBZ, M, N, A, LDA, singular_values, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO)
          
          if (INFO==0) then
             print *,"SVD (DGESDD) successful."
             print *,"Singular values:",singular_values
          else if (INFO>0) then
             print *,"Error in SVD (DGESDD): Did not converge."
             stop
          else
             print *,"Error in SVD (DGESDD): Argument",INFO," was invalid."
             stop
          end if
          do ixi = 1,Nxi
             ! U^T(1:N, 1:Nxi) = U(1:Nxi, 1:N)
             ! V(1, 1:N) = VT(1:N, 1)
             xiWeights(ixi) = 2*dot_product(VT(1:num_Legendres_to_orthogonalize, 1) / singular_values, U(ixi,1:num_Legendres_to_orthogonalize)) 
          end do
          print *,"xiWeights:"
          print *,xiWeights
          print *,"sum(xiWeights):",sum(xiWeights)
          print *,"sum(xiWeights*xi*xi) (should be 2/3):",sum(xi*xi*xiWeights)
          print *,"sum(xiWeights*xi*xi*xi*xi) (should be 2/5):",sum(xi*xi*xi*xi*xiWeights)

       end if

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

    allocate(f_scaling(Nx,Nspecies))
    do ispecies = 1,Nspecies
       select case (f_scaling_option)
       case (1)
          ! Expected magnitude of the leading-order Maxwellian, without the exponential.
          f_scaling(:,ispecies) = nHats(ispecies) / (pi*sqrtpi*(thermal_speeds(ispecies)**3))
       case (2)
          ! rho* times the leading-order Maxwellian, without the exponential.
          f_scaling(:,ispecies) = sqrt(2 * mHats(ispecies) * proton_mass * THats(ispecies) * electron_charge) &
               / (Zs(ispecies) * electron_charge * sqrt(FSABHat2) * aHat) &
               * nHats(ispecies) / (pi*sqrtpi*(thermal_speeds(ispecies)**3))
       case (3)
          ! Expected magnitude of the leading-order Maxwellian, with the exponential.
          f_scaling(:,ispecies) = expx2 * nHats(ispecies) / (pi*sqrtpi*(thermal_speeds(ispecies)**3))
       case (4)
          ! rho* times the leading-order Maxwellian, with the exponential.
          f_scaling(:,ispecies) = expx2 * sqrt(2 * mHats(ispecies) * proton_mass * THats(ispecies) * electron_charge) &
               / (Zs(ispecies) * electron_charge * sqrt(FSABHat2) * aHat) &
               * nHats(ispecies) / (pi*sqrtpi*(thermal_speeds(ispecies)**3))
       case default
          if (masterProc) print *,"Error! Invalid f_scaling_option:",f_scaling_option
          stop
       end select
    end do
    if (masterProc) then
       print "(a,i2,a)","f_scaling_option =",f_scaling_option,". Here comes f_scaling:"
       do ix=1,Nx
          print "(*(es10.2))",f_scaling(ix,:)
       end do
    end if

    allocate(x_scaling(Nx,Nspecies))
    do ispecies = 1,Nspecies
       !v_s = sqrt(2*THats(ispecies)/mHats(ispecies)) ! Once I switch to SI units, include the 2 here.
       !v_s = sqrt(THats(ispecies)/mHats(ispecies))    ! But while using the old units, v_s is measured in units of vBar, so there is no 2.

       select case (x_scaling_option)
       case (1)
          ! without the exponential, without the x.
          x_scaling(:,ispecies) = 1 / thermal_speeds(ispecies)
       case (2)
          ! without the exponential, with the x.
          x_scaling(:,ispecies) = 1 / (thermal_speeds(ispecies) * x)
       case (3)
          ! with the exponential, without the x.
          x_scaling(:,ispecies) = 1 / (thermal_speeds(ispecies) * expx2)
       case (4)
          ! with the exponential, with the x.
          x_scaling(:,ispecies) = 1 / (thermal_speeds(ispecies) * expx2 * x)
       case default
          if (masterProc) print *,"Error! Invalid x_scaling_option:",x_scaling_option
          stop
       end select

!!$       select case (x_scaling_option)
!!$       case (1)
!!$          x_scaling(:,ispecies) = 1 / (x * v_s)
!!$       case (2)
!!$          x_scaling(:,ispecies) = 1 / v_s
!!$       case (3)
!!$          x_scaling(:,ispecies) = 1 / (x * v_s * expx2)
!!$       case (4)
!!$          x_scaling(:,ispecies) = 1 / (v_s * expx2)
!!$       case default
!!$          if (masterProc) print *,"Error! Invalid x_scaling_option:",x_scaling_option
!!$          stop
!!$       end select
    end do
    if (masterProc) then
       print "(a,i2,a)","x_scaling_option =",x_scaling_option,". Here comes x_scaling:"
       do ix=1,Nx
          print "(*(es10.2))",x_scaling(ix,:)
       end do
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
       print "(a,i2,a)","spatial_scaling_option =",spatial_scaling_option,". Here comes spatial_scaling:"
       do itheta=1,Ntheta
          print "(*(f5.2))",spatial_scaling(itheta,:)
       end do
    end if

    allocate(xi_scaling(Nxi))
    select case (xi_scaling_option)
    case (0)
       xi_scaling = 1
    case (1)
       xi_scaling = sqrt(1 - xi*xi)
       ! Set end points exactly in case the line above yields a tiny nonzero value due to roundoff error.
       !xi_scaling(1) = 0
       !xi_scaling(Nxi) = 0
    case default
       if (masterProc) print *,"Error! Invalid setting for xi_scaling_option:",xi_scaling_option
       stop
    end select
    if (masterProc) print *,"xi_scaling:",xi_scaling

    do ixi = 1,Nxi
       pitch_angle_scattering_operator(ixi,:) = xi_scaling(ixi) * pitch_angle_scattering_operator(ixi,:)
       pitch_angle_scattering_operator_preconditioner(ixi,:) = xi_scaling(ixi) * pitch_angle_scattering_operator_preconditioner(ixi,:)
       do j = 1,NL
          Legendre_projection(ixi,:,j) = xi_scaling(ixi) * Legendre_projection(ixi,:,j)
       end do
    end do
!!$    if (xi_scaling_option==1) then
!!$       ! Fix up the first and last rows of pitch_angle_scattering_operator.
!!$       ! Recall pitch angle scattering is defined as xi_scaling * (1/2) d/dxi [(1-xi^2) df/dxi] = xi_scaling * (1/2)(1/sin(alpha)) d/dalpha [sin(alpha) df/dalpha]
!!$       ! so at the boundaries, if xi_scaling = sin(alpha), then pitch angle scattering = (1/2) cos(alpha) df/dalpha.
!!$       pitch_angle_scattering_operator(  1,:) = -(1/two) * ddalpha_plus(   1,:)
!!$       pitch_angle_scattering_operator(Nxi,:) =  (1/two) * ddalpha_minus(Nxi,:)
!!$    end if

    ! *********************************************************
    ! Compute Rosenbluth potential response matrices
    ! *********************************************************

    if ((xGridScheme==5 .or. xGridScheme==6) .and. (RHSMode .ne. 3)) then
       allocate(RosenbluthPotentialTerms(Nspecies,Nspecies,NL,Nx,Nx))
       !call computeRosenbluthPotentialResponse(Nx, x, xWeights, Nspecies, mHats, THats, nHats, Zs, NL, f_scaling, RosenbluthPotentialTerms,.false.)
       call computeRosenbluthPotentialResponse()
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

    if (masterProc) then
       print *,"Here comes ddtheta_plus:"
       do j=1,Ntheta
          print "(*(f6.1))",ddtheta_plus(j,:)
       end do
       print *,"Here comes ddtheta_minus:"
       do j=1,Ntheta
          print "(*(f6.1))",ddtheta_minus(j,:)
       end do
       print *,"Here comes ddtheta_plus_preconditioner:"
       do j=1,Ntheta
          print "(*(f6.1))",ddtheta_plus_preconditioner(j,:)
       end do
       print *,"Here comes ddtheta_minus_preconditioner:"
       do j=1,Ntheta
          print "(*(f6.1))",ddtheta_minus_preconditioner(j,:)
       end do

       print *,"Here comes ddzeta_plus:"
       do j=1,Nzeta
          print "(*(f6.1))",ddzeta_plus(j,:)
       end do
       print *,"Here comes ddzeta_minus:"
       do j=1,Nzeta
          print "(*(f6.1))",ddzeta_minus(j,:)
       end do
       print *,"Here comes ddzeta_plus_preconditioner:"
       do j=1,Nzeta
          print "(*(f6.1))",ddzeta_plus_preconditioner(j,:)
       end do
       print *,"Here comes ddzeta_minus_preconditioner:"
       do j=1,Nzeta
          print "(*(f6.1))",ddzeta_minus_preconditioner(j,:)
       end do

       print *,"Here comes ddxi_plus:"
       do j=1,Nxi
          print "(*(f6.1))",ddxi_plus(j,:)
       end do
       print *,"Here comes ddxi_minus:"
       do j=1,Nxi
          print "(*(f6.1))",ddxi_minus(j,:)
       end do
       print *,"Here comes ddxi_plus_preconditioner:"
       do j=1,Nxi
          print "(*(f6.1))",ddxi_plus_preconditioner(j,:)
       end do
       print *,"Here comes ddxi_minus_preconditioner:"
       do j=1,Nxi
          print "(*(f6.1))",ddxi_minus_preconditioner(j,:)
       end do

       print *,"Here comes pitch_angle_scattering_operator:"
       do j=1,Nxi
          print "(*(f6.1))",pitch_angle_scattering_operator(j,:)
       end do
       print *,"Here comes pitch_angle_scattering_operator_preconditioner:"
       do j=1,Nxi
          print "(*(f6.1))",pitch_angle_scattering_operator_preconditioner(j,:)
       end do
    end if

    if (masterProc) then
       print *,"------------------------------------------------------"
       print *,"Done creating grids."
    end if

  contains

    function compute_xi_from_y(yy)
      implicit none
      real(prec), intent(in) :: yy
      real(prec) :: compute_xi_from_y
      
      compute_xi_from_y = nonuniform_xi_a * yy + nonuniform_xi_b * (yy ** 5)
    end function compute_xi_from_y
    
    function compute_dxi_dy(yy)
      implicit none
      real(prec), intent(in) :: yy
      real(prec) :: compute_dxi_dy
      
      compute_dxi_dy = nonuniform_xi_a + 5*nonuniform_xi_b * (yy ** 4)
    end function compute_dxi_dy
    
    function compute_d2xi_dy2(yy)
      implicit none
      real(prec), intent(in) :: yy
      real(prec) :: compute_d2xi_dy2
      
      compute_d2xi_dy2 = 20*nonuniform_xi_b * (yy ** 3)
    end function compute_d2xi_dy2
    
  end subroutine createGrids
