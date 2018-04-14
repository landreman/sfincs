#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscdmdadef.h>
#else
#include <petsc/finclude/petscdmdadef.h>
#endif

  ! magnetic drift derivative and preconditioner matrices not needed
  subroutine createGrids_fine()

    use globalVariables
    use petscdmda

    PetscErrorCode :: ierr
    integer :: i, j, k
    PetscScalar, dimension(:,:), allocatable :: d2dtheta2, d2dzeta2, temp_matrix
    PetscScalar :: temp
    logical :: describe_stencils
    DM :: myDM
    integer, parameter :: bufferLength = 200
    character(len=bufferLength) :: procAssignments
    PetscScalar :: zetaMax_fine

    integer :: tag, dummy(1), scheme
    integer :: status(MPI_STATUS_SIZE)
    logical :: call_uniform_diff_matrices
    integer :: derivative_option_plus, derivative_option_minus, derivative_option, quadrature_option, derivative_option_diffusion

    ! We don't actually need to keep these
    PetscScalar, dimension(:), allocatable :: theta_preconditioner, thetaWeights_preconditioner
    PetscScalar, dimension(:,:), allocatable :: d2dtheta2_preconditioner
    PetscScalar, dimension(:), allocatable :: zeta_preconditioner, zetaWeights_preconditioner
    PetscScalar, dimension(:,:), allocatable :: d2dzeta2_preconditioner

    ! Compute matrixSize_fine 
    DKE_size_fine = sum(Nxi_for_x)*Ntheta_fine*Nzeta_fine
    matrixSize_fine = Nspecies * DKE_size_fine
    ! Constraints
    matrixSize_fine = matrixSize_fine + 2 * Nspecies

    if (masterProc) print *, "---- Initializing grids for adjoint error correction"

    ! *******************************************************************************
    ! Do a few sundry initialization tasks:
    ! *******************************************************************************

    if (forceOddNthetaAndNzeta) then
       if (mod(Ntheta_fine, 2) == 0) then
          Ntheta_fine = Ntheta_fine + 1
       end if
       if (mod(Nzeta_fine, 2) == 0) then
          Nzeta_fine = Nzeta_fine + 1
       end if
    end if

    if (Ntheta_fine > Nzeta_fine) then
         ! Distribute in theta but not in zeta

         ! Assign a range of theta indices to each processor.
         ! This is done by creating a PETSc DM that is not actually used for anything else.
         call DMDACreate1d(MPIComm, DM_BOUNDARY_NONE, Ntheta_fine, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)

         call DMDAGetCorners(myDM, ithetaMin_fine, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
              localNtheta_fine, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

         call DMDestroy(myDM, ierr)

         izetaMin_fine = 0
         izetaMax_fine = Nzeta_fine-1
         localNzeta_fine = Nzeta_fine
      else
         ! Distribute in zeta but not in theta

         ! Assign a range of zeta indices to each processor.
         ! This is done by creating a PETSc DM that is not actually used for anything else.
         call DMDACreate1d(MPIComm, DM_BOUNDARY_NONE, Nzeta_fine, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)

         call DMDAGetCorners(myDM, izetaMin_fine, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
              localNzeta_fine, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

         call DMDestroy(myDM, ierr)
         ithetaMin_fine = 0
         ithetaMax_fine = Ntheta_fine-1
         localNtheta_fine = Ntheta_fine
      end if

    ! Switch to 1-based indices:
    ithetaMin_fine = ithetaMin_fine + 1
    ithetaMax_fine = ithetaMin_fine+localNtheta_fine-1
    izetaMin_fine = izetaMin_fine + 1
    izetaMax_fine = izetaMin_fine+localNzeta_fine-1

    procThatHandlesConstraints = masterProc

    write (procAssignments,fmt="(a,i4,a,i3,a,i3,a,i3,a,i3,a)") "Processor ",myRank," owns theta indices ",ithetaMin_fine," to ",ithetaMax_fine,&
         " and zeta indices ",izetaMin_fine," to ",izetaMax_fine

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

    ! global vairables
    allocate(theta_fine(Ntheta_fine))
    allocate(thetaWeights_fine(Ntheta_fine))
    allocate(ddtheta_ExB_plus_fine(Ntheta_fine,Ntheta_fine))
    allocate(ddtheta_ExB_minus_fine(Ntheta_fine,Ntheta_fine))
    allocate(ddtheta_fine(Ntheta_fine,Ntheta_fine))

    ! local variables
    allocate(theta_preconditioner(Ntheta_fine))
    allocate(thetaWeights_preconditioner(Ntheta_fine))
    allocate(d2dtheta2_preconditioner(Ntheta_fine,Ntheta_fine))
    allocate(temp_matrix(Ntheta_fine,Ntheta_fine))
    allocate(d2dtheta2(Ntheta_fine,Ntheta_fine))

    ! *******************************************************************************
    ! Handle d/dtheta for the main matrix.
    ! *******************************************************************************

    select case (thetaDerivativeScheme)
    case (0)
       scheme = 20
    case (1)
       scheme = 0
    case (2)
       scheme = 10
    case default
       if (masterProc) then
          print *,"Error! Invalid setting for thetaDerivativeScheme"
       end if
       stop
    end select

    call uniformDiffMatrices(Ntheta_fine, zero, two*pi, scheme, theta_fine, thetaWeights_fine, ddtheta_fine, d2dtheta2)

    ! Create upwinded matrices for ExB terms:
    !print *,"Creating upwinded matrices for ExB terms, theta"
    select case (ExBDerivativeSchemeTheta)
    case (0)
       ! It should not matter what ddtheta_ExB_plus and ddtheta_ExB_minus are in this case.
       ddtheta_ExB_plus_fine = ddtheta_fine
       ddtheta_ExB_minus_fine = ddtheta_fine
    case (1)
       scheme = 80
       call uniformDiffMatrices(Ntheta_fine, zero, two*pi, scheme, theta_preconditioner, &
            thetaWeights_preconditioner, ddtheta_ExB_plus_fine, d2dtheta2_preconditioner)
       scheme = 90
       call uniformDiffMatrices(Ntheta_fine, zero, two*pi, scheme, theta_preconditioner_fine, &
            thetaWeights_preconditioner, ddtheta_ExB_minus_fine, d2dtheta2_preconditioner)
    case (2)
       scheme = 100
       call uniformDiffMatrices(Ntheta_fine, zero, two*pi, scheme, theta_preconditioner_fine, &
            thetaWeights_preconditioner, ddtheta_ExB_plus_fine, d2dtheta2_preconditioner)
       scheme = 110
       call uniformDiffMatrices(Ntheta_fine, zero, two*pi, scheme, theta_preconditioner, &
            thetaWeights_preconditioner, ddtheta_ExB_minus_fine, d2dtheta2_preconditioner)
    case (3)
       scheme = 120
       call uniformDiffMatrices(Ntheta_fine, zero, two*pi, scheme, theta_preconditioner, &
            thetaWeights_preconditioner, ddtheta_ExB_plus_fine, d2dtheta2_preconditioner)
       scheme = 130
       call uniformDiffMatrices(Ntheta_fine, zero, two*pi, scheme, theta_preconditioner, &
            thetaWeights_preconditioner, ddtheta_ExB_minus_fine, d2dtheta2_preconditioner)
    case default
       print *,"Error! Invalid ExBDerivativeSchemeTheta:",ExBDerivativeSchemeTheta
       stop
    end select

    ! preconditioner and magnetic drift matrices not needed

    ! The following arrays will not be needed:
    deallocate(d2dtheta2)
    deallocate(theta_preconditioner)
    deallocate(thetaWeights_preconditioner)
    deallocate(d2dtheta2_preconditioner)

    ! *******************************************************************************
    ! Build zeta grid, integration weights, and differentiation matrices:
    ! *******************************************************************************

    zetaMax_fine = 2*pi/NPeriods

    allocate(zeta_fine(Nzeta_fine))
    allocate(zetaWeights_fine(Nzeta_fine))
    allocate(ddzeta_fine(Nzeta_fine,Nzeta_fine))
    allocate(ddzeta_ExB_plus_fine(Nzeta_fine,Nzeta_fine))
    allocate(ddzeta_ExB_minus_fine(Nzeta_fine,Nzeta_fine))

    allocate(d2dzeta2(Nzeta_fine,Nzeta_fine))
    allocate(zeta_preconditioner(Nzeta_fine))
    allocate(zetaWeights_preconditioner(Nzeta_fine))
    allocate(d2dzeta2_preconditioner(Nzeta_fine,Nzeta_fine))

    select case (zetaDerivativeScheme)
    case (0)
       scheme = 20
    case (1)
       scheme = 0
    case (2)
       scheme = 10
    case default
       if (masterProc) then
          print *,"Error! Invalid setting for zetaDerivativeScheme"
       end if
       stop
    end select

    if (Nzeta_fine==1) then
       ! Axisymmetry:
       zeta_fine = 0
       zetaWeights_fine = 2*pi
       ddzeta_fine = 0
       d2dzeta2 = 0 ! d2dzeta2 is not actually used.
       ddzeta_ExB_plus_fine = 0
       ddzeta_ExB_minus_fine = 0
    else
       call uniformDiffMatrices(Nzeta_fine, zero, zetaMax_fine, scheme, zeta_fine, zetaWeights_fine, ddzeta_fine, d2dzeta2)
      
       ! Create upwinded matrices for ExB terms:
       !print *,"Creating upwinded matrices for ExB terms, zeta"
       select case (ExBDerivativeSchemeZeta)
       case (0)
          ! It should not matter what ddzeta_ExB_plus and ddzeta_ExB_minus are in this case.
          ddzeta_ExB_plus_fine = ddzeta_fine
          ddzeta_ExB_minus_fine = ddzeta_fine
       case (1)
          scheme = 80
          call uniformDiffMatrices(Nzeta_fine, zero, zetaMax_fine, scheme, zeta_preconditioner, &
               zetaWeights_preconditioner, ddzeta_ExB_plus_fine, d2dzeta2_preconditioner)
          scheme = 90
          call uniformDiffMatrices(Nzeta_fine, zero, zetaMax_fine, scheme, zeta_preconditioner, &
               zetaWeights_preconditioner, ddzeta_ExB_minus_fine, d2dzeta2_preconditioner)
       case (2)
          scheme = 100
          call uniformDiffMatrices(Nzeta_fine, zero, zetaMax_fine, scheme, zeta_preconditioner, &
               zetaWeights_preconditioner, ddzeta_ExB_plus_fine, d2dzeta2_preconditioner)
          scheme = 110
          call uniformDiffMatrices(Nzeta_fine, zero, zetaMax_fine, scheme, zeta_preconditioner, &
               zetaWeights_preconditioner, ddzeta_ExB_minus_fine, d2dzeta2_preconditioner)
       case (3)
          scheme = 120
          call uniformDiffMatrices(Nzeta_fine, zero, zetaMax_fine, scheme, zeta_preconditioner, &
               zetaWeights_preconditioner, ddzeta_ExB_plus_fine, d2dzeta2_preconditioner)
          scheme = 130
          call uniformDiffMatrices(Nzeta_fine, zero, zetaMax_fine, scheme, zeta_preconditioner, &
               zetaWeights_preconditioner, ddzeta_ExB_minus_fine, d2dzeta2_preconditioner)
       case default
          print *,"Error! Invalid ExBDerivativeSchemeZeta:",ExBDerivativeSchemeZeta
          stop
       end select

    end if

    zetaWeights_fine = zetaWeights_fine * NPeriods

    ! The following arrays will not be needed:
    deallocate(d2dzeta2)
    deallocate(zeta_preconditioner)
    deallocate(zetaWeights_preconditioner)
    deallocate(d2dzeta2_preconditioner)

    ! Allocate geometry arrays
    allocate(BHat_fine(Ntheta_fine,Nzeta_fine))
    allocate(DHat_fine(Ntheta_fine,Nzeta_fine))
    allocate(dBHatdtheta_fine(Ntheta_fine,Nzeta_fine))
    allocate(dBHatdzeta_fine(Ntheta_fine,Nzeta_fine))
    allocate(BHat_sub_theta_fine(Ntheta_fine,Nzeta_fine))
    allocate(dBHat_sub_theta_dzeta_fine(Ntheta_fine,Nzeta_fine))
    allocate(BHat_sub_zeta_fine(Ntheta_fine,Nzeta_fine))
    allocate(dBHat_sub_zeta_dtheta_fine(Ntheta_fine,Nzeta_fine))
    allocate(BHat_sup_theta_fine(Ntheta_fine,Nzeta_fine))
    allocate(dBHat_sup_theta_dzeta_fine(Ntheta_fine,Nzeta_fine))
    allocate(BHat_sup_zeta_fine(Ntheta_fine,Nzeta_fine))
    allocate(dBHat_sup_zeta_dtheta_fine(Ntheta_fine,Nzeta_fine))

    ! Buil interpolation matrix using a 4 point stencil
    call build_adjointECProlongationMatrix()

 end subroutine createGrids_fine
