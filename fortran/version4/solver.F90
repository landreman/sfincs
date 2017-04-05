#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif

module solver

  use petscsnesdef
  use petsckspdef
  use globalVariables
  use indices

  implicit none

  contains

  subroutine mainSolverLoop

    implicit none

    PetscErrorCode :: ierr
    Vec :: solutionVec, residualVec
    Mat :: matrix  !, preconditionerMatrix
    SNES :: mysnes
    KSP :: outer_KSP
    PC :: outer_preconditioner, sub_preconditioner
    integer :: numRHSs
    SNESConvergedReason :: reason
    KSPConvergedReason :: KSPReason
    PetscLogDouble :: time1, time2
    integer :: userContext(1) = 0
    Vec :: dummyVec
    Mat :: factorMat
    PetscInt :: mumps_which_cntl
    PetscReal :: mumps_value
    PetscReal :: atol, rtol, stol
    integer :: maxit, maxf
    PetscInt :: VecLocalSize
    IS, dimension(:), allocatable :: ISs
    integer :: IS_index, IS_array_index, ix, ispecies, ixi, itheta, izeta, num_fieldsplits, j
    integer, dimension(:), allocatable :: IS_array
    KSP, dimension(:), allocatable :: sub_ksps
    Mat :: sub_Amat, sub_Pmat
    MatNullSpace :: nullspace
    Vec :: temp_Vec
    integer :: fieldsplit_index_min, fieldsplit_index_max

    external apply_Jacobian
    external apply_preconditioner

!!Following three lines added by AM 2016-07-06
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
    PetscViewerAndFormat vf !!Added by AM 2016-07-06
#endif

    external evaluateJacobian, evaluateResidual, diagnosticsMonitor

    if (masterProc) then
       print *,"Setting up PETSc solver."
    end if
    iterationForMatrixOutput = 0

    call PetscOptionsInsertString(PETSC_NULL_OBJECT,"-options_left -snes_view", ierr)
!    call PetscOptionsInsertString(PETSC_NULL_OBJECT,"-options_left -snes_view -log_summary", ierr)

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, solutionVec, ierr)
    call VecDuplicate(solutionVec, residualVec, ierr)

    call VecDuplicate(solutionVec, f0, ierr)
    call init_f0()

    call SNESCreate(MPIComm, mysnes, ierr)
    call SNESAppendOptionsPrefix(mysnes, 'outer_', ierr)
    call SNESSetFunction(mysnes, residualVec, evaluateResidual, PETSC_NULL_OBJECT, ierr)

    firstMatrixCreation = .true.

    ! Create the Mat object for the Jacobian
    call VecGetLocalSize(solutionVec,VecLocalSize,ierr)
    call MatCreateShell(PETSC_COMM_WORLD,VecLocalSize,VecLocalSize,matrixSize,matrixSize,&
         PETSC_NULL_OBJECT,matrix,ierr)
    call MatShellSetOperation(matrix,MATOP_MULT,apply_Jacobian,ierr)

    call preallocateMatrix(Mat_for_preconditioner, 0)
    call SNESSetJacobian(mysnes, matrix, Mat_for_preconditioner, evaluateJacobian, PETSC_NULL_OBJECT, ierr)

    call SNESGetKSP(mysnes, outer_KSP, ierr)
    call KSPGetPC(outer_KSP, outer_preconditioner, ierr)
    call PCSetType(outer_preconditioner, PCSHELL, ierr)
    call PCShellSetApply(outer_preconditioner, apply_preconditioner, ierr)

    call KSPCreate(MPIComm, inner_ksp, ierr)
    call KSPAppendOptionsPrefix(inner_ksp, 'inner_', ierr)
    call KSPSetType(inner_ksp, KSPPREONLY, ierr)
    call KSPGetPC(inner_ksp, inner_preconditioner, ierr)
    call KSPSetOperators(inner_ksp, Mat_for_preconditioner, Mat_for_preconditioner, ierr)

    select case (preconditioning_option)
    case (1)
       ! Direct solver
       fieldsplit = .false.
       call PCSetType(inner_preconditioner, PCLU, ierr)
       call KSPSetType(outer_KSP, KSPGMRES, ierr)   ! Set the Krylov solver algorithm to GMRES
       call KSPGMRESSetRestart(outer_KSP, 2000, ierr)

       if (isAParallelDirectSolverInstalled) then
          select case (whichParallelSolverToFactorPreconditioner)
          case (1)
             call PCFactorSetMatSolverPackage(inner_preconditioner, MATSOLVERMUMPS, ierr)
             if (masterProc) then
                print *,"We will use mumps to factorize the preconditioner."
             end if
          case (2)
             call PCFactorSetMatSolverPackage(inner_preconditioner, MATSOLVERSUPERLU_DIST, ierr)
             if (masterProc) then
                print *,"We will use superlu_dist to factorize the preconditioner."
             end if
             ! Turn on superlu_dist diagnostic output:
             call PetscOptionsInsertString("-mat_superlu_dist_statprint", ierr)
          case default
             if (masterProc) then
                print *,"Error: Invalid setting for whichParallelSolverToFactorPreconditioner"
                stop
             end if
          end select
       else
          if (masterProc) then
             print *,"We will use PETSc's serial sparse direct solver to factorize the preconditioner."
          end if

          ! When using PETSc's built-in solver (which is only done when running with a single processor),
          ! I often get 2 kinds of errors. One kind of error is a "zero pivot" error.  Other times, there is
          ! no explicit error message, but KSP converges really slowly, and the answers come out wrong.
          ! The next few lines seem to help minimize this problem:

          ! The available orderings are natural, nd, 1wd, rcm, and qmd.
          ! natural (i.e. no reordering) is horrible since the LU factors are very dense.
          ! nd, rcm, and qmd all seem to work with some examples but fail with others. You should try one of these 3 options.
          ! 1wd seems to use more memory.
          call PCFactorSetMatOrderingType(inner_preconditioner, MATORDERINGRCM, ierr)
          call PCFactorReorderForNonzeroDiagonal(inner_preconditioner, 1d-12, ierr) 
          call PCFactorSetZeroPivot(inner_preconditioner, 1d-200, ierr) 
       end if

       ! End of steps used for direct LU factorization of the preconditioner.

    case (2,3,4,5,6,7)
       ! Fieldsplit
       fieldsplit = .true.
       if (masterProc) then
          print *,"Setting up PETSc FIELDSPLIT preconditioning"
          select case (preconditioning_option)
          case (2,5)
             print *,"  Using PETSc's GAMG (algebraic multigrid) on the blocks for the kinetic equation."
          case (3,6)
             print *,"  Using Hypre BoomerAMG (algebraic multigrid) on the blocks for the kinetic equation."
          case (4,7)
             print *,"  Using LU factorization on the blocks for the kinetic equation."
          case default
             print *,"Error! Invalid preconditioning_option:",preconditioning_option
             stop
          end select
       end if
       call PCSetType(inner_preconditioner, PCFIELDSPLIT, ierr)
       call PCFIELDSPLITSetType(inner_preconditioner, PC_COMPOSITE_ADDITIVE, ierr)

       allocate(ISs(Nx*Nspecies+1))
       allocate(is_array(matrixSize))
       IS_index = 1
       ! For splitting the fields among procs, copy the division PETSc uses for a Vec:
       call VecCreateMPI(MPIComm, PETSC_DECIDE, Ntheta*Nzeta*Nxi, temp_Vec, ierr)
       call VecGetOwnershipRange(temp_Vec, fieldsplit_index_min, fieldsplit_index_max, ierr)
       call VecDestroy(temp_Vec, ierr)
       ! The order of these loops must be consistent with getIndices(), or else PETSc gives an error about the IS values not being sorted.
       do ispecies = 1,Nspecies
          do ix = 1,Nx
             IS_array_index = 1
             do ixi = 1,Nxi
                do itheta = 1,Ntheta
                   do izeta = 1,Nzeta
                      IS_array(IS_array_index) = getIndex(ispecies, ix, ixi, itheta, izeta, BLOCK_F)
                      IS_array_index = IS_array_index + 1
                   end do
                end do
             end do
             select case (fieldsplit_parallelization_option)
             case (1)
                call ISCreateGeneral(PETSC_COMM_WORLD,fieldsplit_index_max-fieldsplit_index_min,IS_array(fieldsplit_index_min+1:fieldsplit_index_max),PETSC_COPY_VALUES,ISs(IS_index),ierr)
             case default
                print *,"Error! Invalid fieldsplit_parallelization_option:",fieldsplit_parallelization_option
                stop
             end select
             !if (masterProc) print *,"Here comes the IS for fieldsplit index",IS_index
             !call ISView(ISs(IS_index),PETSC_VIEWER_STDOUT_WORLD,ierr)
             call PCFieldSplitSetIS(inner_preconditioner,'f',ISs(IS_index),ierr)
             IS_index = IS_index + 1
          end do
       end do

       ! Set up an index set (IS) for the sources and constraints
       select case (constraintScheme)
       case (0)
          ! Nothing to do
       case (1,3,4)
          IS_array_index = 1
          do ispecies = 1,Nspecies
             IS_array(IS_array_index) = getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)
             IS_array_index = IS_array_index + 1
             IS_array(IS_array_index) = getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)
             IS_array_index = IS_array_index + 1
          end do
          if (myRank == numProcs-1) then
             ! This rank handles the whole field: 
             call ISCreateGeneral(PETSC_COMM_WORLD,2*Nspecies,IS_array,PETSC_COPY_VALUES,ISs(IS_index),ierr)
          else
             ! Other ranks don't handle any of this field:
             call ISCreateGeneral(PETSC_COMM_WORLD,0,IS_array,PETSC_COPY_VALUES,ISs(IS_index),ierr)
          end if
          call ISView(ISs(IS_index),PETSC_VIEWER_STDOUT_WORLD,ierr)
          call PCFieldSplitSetIS(inner_preconditioner,'constraints',ISs(IS_index),ierr)
       case (2)
          IS_array_index = 1
          do ispecies = 1,Nspecies
             do ix = 1,Nx
                IS_array(IS_array_index) = getIndex(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT)
                IS_array_index = IS_array_index + 1
             end do
          end do
          if (myRank == numProcs-1) then
             ! This rank handles the whole field:
             call ISCreateGeneral(PETSC_COMM_WORLD,Nx*Nspecies,IS_array,PETSC_COPY_VALUES,ISs(IS_index),ierr)
          else
             ! Other ranks don't handle any of this field:
             call ISCreateGeneral(PETSC_COMM_WORLD,0,IS_array,PETSC_COPY_VALUES,ISs(IS_index),ierr)
          end if
          call ISView(ISs(IS_index),PETSC_VIEWER_STDOUT_WORLD,ierr)
          call PCFieldSplitSetIS(inner_preconditioner,'constraints',ISs(IS_index),ierr)
       case default
          if (masterProc) print *,"Invalid constraintScheme:",constraintScheme
          stop
       end select

       call KSPSetType(outer_KSP, KSPFGMRES, ierr)   ! Set the Krylov solver algorithm to Flexible GMRES
       call KSPGMRESSetRestart(outer_KSP, 200, ierr)

       allocate(sub_ksps(Nspecies*Nx+1))
       call PCFieldSplitGetSubKSP(inner_preconditioner, num_fieldsplits, sub_ksps, ierr)
       if (constraintScheme .ne. 0) then
          ! Set the source/constraint diagonal block to use Jacobi
          j = Nspecies*Nx + 1
          call KSPSetType(sub_ksps(j),KSPPREONLY,ierr)
          call KSPGetPC(sub_ksps(j), sub_preconditioner, ierr)
          call PCSetType(sub_preconditioner, PCJACOBI, ierr)
       end if

       do j = 1,Nspecies*Nx
          ! For each of the diagonal blocks corresponding to the kinetic equation, set the appropriate preconditioner:
          call KSPSetType(sub_ksps(j),KSPPREONLY,ierr)
          !call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
          !call KSPSetNullSpace(sub_ksps(j),nullspace,ierr)
          !call MatNullSpaceDestroy(nullspace,ierr)

          !call KSPGetOperators(sub_ksps(j), sub_Amat, sub_Pmat, ierr)
          !print *,"Does sub_Amat==sub_Pmat?",sub_Amat==sub_Pmat
          !call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
          !call MatSetNullSpace(sub_Pmat,nullspace,ierr)
          !call MatNullSpaceDestroy(nullspace,ierr)

          call KSPGetPC(sub_ksps(j), sub_preconditioner, ierr)
          select case (preconditioning_option)
          case (2,5)
             call PCSetType(sub_preconditioner, PCGAMG, ierr)
          case (3,6)
             call PCSetType(sub_preconditioner, PCHYPRE, ierr)
          case (4,7)
             call PCSetType(sub_preconditioner, PCLU, ierr)
             if (isAParallelDirectSolverInstalled) then
                select case (whichParallelSolverToFactorPreconditioner)
                case (1)
                   call PCFactorSetMatSolverPackage(sub_preconditioner, MATSOLVERMUMPS, ierr)
                case (2)
                   call PCFactorSetMatSolverPackage(sub_preconditioner, MATSOLVERSUPERLU_DIST, ierr)
                end select
             end if
          case default
             stop "Should not get here!"
          end select
       end do

    case default
       if (masterProc) print *,"Error! Invalid preconditioning_option:",preconditioning_option
       stop
    end select

    call KSPSetTolerances(outer_KSP, solverTolerance, PETSC_DEFAULT_REAL, &
         PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)

    ! Allow options to be controlled using command-line flags:
    call KSPSetFromOptions(outer_KSP, ierr)

#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 7))
    call KSPMonitorSet(outer_KSP, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
#else
    call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr) !!Added by AM 2016-07-06
    !call KSPMonitorSet(outer_KSP, KSPMonitorDefault, vf, PetscViewerAndFormatDestroy, ierr) !!Added by AM 2016-07-06 
    call KSPMonitorSet(outer_KSP, KSPMonitorTrueResidualNorm, vf, PetscViewerAndFormatDestroy, ierr) !!Added by AM 2016-07-06 
#endif
    
    ! Tell PETSc to call the diagnostics subroutine at each iteration of SNES:
    call SNESMonitorSet(mysnes, diagnosticsMonitor, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)

    if (reusePreconditioner) then
       ! Syntax for PETSc version 3.5 and later
       ! Do I actually need all 4 of these commands?
       call KSPSetReusePreconditioner(outer_KSP, PETSC_TRUE, ierr)
       call KSPSetReusePreconditioner(inner_KSP, PETSC_TRUE, ierr)
       call PCSetReusePreconditioner(outer_preconditioner, PETSC_TRUE, ierr)
       call PCSetReusePreconditioner(inner_preconditioner, PETSC_TRUE, ierr)
    end if

    ! In older versions of PETSC (either <3.5 or <3.4, I'm not certain)
    ! the monitor is never called when snes type = SNESKSPONLY.
    ! Therefore, it is preferable to always have snes type = SNESNEWTONLS but set the 
    ! number of iterations to 1 for a linear run.
    call SNESSetType(mysnes, SNESNEWTONLS, ierr)
    !!if (nonlinear) then !!Commented by AM 2016-02
    if (includePhi1) then !!Added by AM 2016-02
       if (masterProc) then
          print *,"Since this is a nonlinear run, we will use Newton's method."
       end if
    else
       call SNESGetTolerances(mysnes, atol, rtol, stol, maxit, maxf, ierr)
       maxit = 1
       call SNESSetTolerances(mysnes, atol, rtol, stol, maxit, maxf, ierr)
       if (masterProc) then
          print *,"Since this is a linear run, we will only take a single step, and not iterate Newton's method."
       end if
    end if
    ! Below is the old way of setting linear vs nonlinear, which does not work for linear runs in petsc 3.3:
!!$    ! Set the algorithm to use for the nonlinear solver:
!!$    if (nonlinear) then
!!$       ! SNESNEWTONLS = Newton's method with an optional line search.
!!$       ! As of PETSc version 3.5, this is the default algorithm, but I'll set it manually here anyway to be safe.
!!$       call SNESSetType(mysnes, SNESNEWTONLS, ierr)
!!$       if (masterProc) then
!!$          print *,"Since this is a nonlinear run, we will use Newton's method."
!!$       end if
!!$    else
!!$       ! SNESKSPONLY = Only do 1 linear step.
!!$       call SNESSetType(mysnes, SNESKSPONLY, ierr)
!!$       if (masterProc) then
!!$          print *,"Since this is a linear run, we will only take a single step, and not iterate Newton's method."
!!$       end if
!!$    end if

    call SNESSetFromOptions(mysnes, ierr)

    if (preconditioning_option==1 .and. isAParallelDirectSolverInstalled .and. (whichParallelSolverToFactorPreconditioner==1)) then
       ! When mumps is the solver, it is very handy to set mumps's control parameter CNTL(1) to a number like 1e-6.
       ! CNTL(1) is a threshhold for pivoting. For the default value of 0.01, there is a lot of pivoting.
       ! This causes memory demands to increase beyond mumps's initial estimate, causing errors (INFO(1)=-9).
       ! If we set CNTL(1) all the way to 0, there is sometimes an error about a zero pivot.
       ! The setting CNTL(1)=1e-6 seems robust.
       ! I notice that if CNTL(1) is smaller, like 1e-15, the Fokker-Planck calculations work find but pitch-angle-scattering calculations
       ! fail. There might be value in using different values of CNTL(1) for FP vs PAS.

#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
       ! Syntax for PETSc versions up through 3.4
#else
       ! These commands must be AFTER SNESSetFromOptions or else there is a seg fault.
       call PCFactorSetUpMatSolverPackage(inner_preconditioner,ierr)
       call PCFactorGetMatrix(inner_preconditioner,factorMat,ierr)
       mumps_which_cntl = 1
       mumps_value = 1.e-6
       call MatMumpsSetCntl(factorMat,mumps_which_cntl,mumps_value,ierr)

       ! Turn on mumps diagnostic output
       mumps_which_cntl = 4
       call MatMumpsSetIcntl(factorMat,mumps_which_cntl,2,ierr)

       ! Increase amount by which the mumps work array can expand due to near-0 pivots.
       ! (The default value of icntl(14) is 25.)
       ! For many cases it is not necessary to increase icntl(14), but an increase sometimes helps
       ! to prevent the mumps error with info(1)=-9.  There appears to be basically
       ! no significant cost in memory or time to increase this parameter.
       mumps_which_cntl = 14
       call MatMumpsSetIcntl(factorMat,mumps_which_cntl,50,ierr)
#endif
    end if

    call KSPSetFromOptions(inner_ksp, ierr)

!    if (fieldsplit) then
!       print *,"AAAAAAAAAAAAA"
!       call SNESSetUp(mysnes, ierr)
!       call KSPSetUp(outer_KSP, ierr)
!!$       do j = 1,Nspecies*Nx
!!$          call KSPGetOperators(sub_ksps(j), sub_Amat, sub_Pmat, ierr)
!!$          print *,"Does sub_Amat==sub_Pmat?",sub_Amat==sub_Pmat
!!$          call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
!!$          call MatSetNullSpace(sub_Pmat,nullspace,ierr)
!!$          call MatNullSpaceDestroy(nullspace,ierr)
!!$       end do
!    end if

    ! ***********************************************************************
    ! ***********************************************************************
    ! 
    !  Main solver call:
    !
    ! ***********************************************************************
    ! ***********************************************************************

    select case (RHSMode)
    case (1)
       ! Single solve, either linear or nonlinear.

       !  Set initial guess for the solution:
       call VecSet(solutionVec, zero, ierr)

       if (masterProc) then
          print *,"------------------------------------------------------"
          print *,"Finished initialization."
          print *,"Beginning the main solve.  This could take a while ..."
       end if

       call PetscTime(time1, ierr)
       if (solveSystem) then
          ! All the magic happens in this next line!
          call SNESSolve(mysnes,PETSC_NULL_OBJECT, solutionVec, ierr)
       end if

       call PetscTime(time2, ierr)
       if (masterProc) then
          print *,"Done with the main solve.  Time to solve: ", time2-time1, " seconds."
       end if
       call PetscTime(time1, ierr)

       !!if (nonlinear) then !!Commented by AM 2016-02
       if (includePhi1) then !!Added by AM 2016-02
          call SNESGetConvergedReason(mysnes, reason, ierr)
          if (reason>0) then
             if (masterProc) then
                print *,"Nonlinear iteration (SNES) converged!  SNESConvergedReason = ", reason
                select case (reason)
                case (2)
                   print *,"  SNES_CONVERGED_FNORM_ABS: ||F|| < atol"
                case (3)
                   print *,"  SNES_CONVERGED_FNORM_RELATIVE: ||F|| < rtol*||F_initial||"
                case (4)
                   print *,"  SNES_CONVERGED_SNORM_RELATIVE: Newton computed step size small; || delta x || < stol || x ||"
                case (5)
                   print *,"  SNES_CONVERGED_ITS: Maximum iterations reached"
                case (7)
                   print *,"  SNES_CONVERGED_TR_DELTA:"
                end select
             end if
             didNonlinearCalculationConverge = integerToRepresentTrue
          else
             if (masterProc) then
                print *,"Nonlinear iteration (SNES) did not converge :(   SNESConvergedReason = ", reason
                select case (reason)
                case (-1)
                   print *,"  SNES_DIVERGED_FUNCTION_DOMAIN: The new x location passed the function is not in the domain of F"
                case (-2)
                   print *,"  SNES_DIVERGED_FUNCTION_COUNT:"
                case (-3)
                   print *,"  SNES_DIVERGED_LINEAR_SOLVE: The linear solve failed"
                case (-4)
                   print *,"  SNES_DIVERGED_FNORM_NAN"
                case (-5)
                   print *,"  SNES_DIVERGED_MAX_IT"
                case (-6)
                   print *,"  SNES_DIVERGED_LINE_SEARCH: The line search failed"
                case (-7)
                   print *,"  SNES_DIVERGED_INNER: Inner solve failed"
                case (-8)
                   print *,"  SNES_DIVERGED_LOCAL_MIN: || J^T b || is small, implies converged to local minimum of F()"
                end select
             end if
             didNonlinearCalculationConverge = integerToRepresentFalse
          end if
       else
          didNonlinearCalculationConverge = integerToRepresentTrue
       end if


       ! End of RHSMode=1 case, which handles a single linear or nonlinear solve.

    case (2,3)
       ! RHSMode = 2 or 3:
       ! Do a linear solve for multiple right-hand sides to get the transport matrix.

       !  Set f=0:
       call VecSet(solutionVec, zero, ierr)

       call KSPSetOperators(outer_KSP, matrix, Mat_for_preconditioner, ierr)

       ! Call evaluateJacobian, which has the effect of populating the main matrix and preconditioner matrix.
       call evaluateJacobian(mysnes, solutionVec, matrix, Mat_for_preconditioner, userContext, ierr) ! The Mat arguments are not actually used.

!!$       ! Build the main linear system matrix:
!!$       !call populateMatrix(matrix,1,dummyVec)
!!$       call preallocateMatrix(Mat_for_Jacobian, 1)
!!$       call populateMatrix(Mat_for_Jacobian, 1, solutionVec)
!!$
!!$       ! Build the preconditioner:
!!$       call populateMatrix(Mat_for_preconditioner,0,dummyVec)

       ! Syntax for PETSc version 3.5 and later
       !call KSPSetOperators(outer_KSP, matrix, Mat_for_preconditioner, ierr)
       !call KSPSetReusePreconditioner(outer_KSP, PETSC_TRUE, ierr)
       !call PCSetReusePreconditioner(inner_preconditioner, PETSC_TRUE, ierr)

       select case (RHSMode)
       case (2)
          numRHSs = 3
       case (3)
          numRHSs = 2
       end select

       do whichRHS = 1,numRHSs
          if (masterProc) then
             print *,"################################################################"
             print "(a,i1,a,i1)"," Solving system with right-hand side ",whichRHS," of ",numRHSs
             print *,"################################################################"
          end if

          call PetscTime(iteration_start_time, ierr)

          ! To get a transport matrix, change the equilibrium gradients here.
          select case (RHSMode)
          case (2)
             ! Energy-integrated transport matrix
             select case (whichRHS)
             case (1)
                dnHatdpsiHats = 1
                dTHatdpsiHats = 0
                EParallelHat = 0
             case (2)
                ! The next 2 lines ensure (1/n)*dn/dpsi + (3/2)*dT/dpsi = 0 while dT/dpsi is nonzero.
                dnHatdpsiHats = (3/two)*nHats(1)*THats(1)
                dTHatdpsiHats = 1
                EParallelHat = 0
             case (3)
                dnHatdpsiHats = 0
                dTHatdpsiHats = 0
                EParallelHat = 1
             end select

          case (3)
             ! Monoenergetic transport matrix
             select case (whichRHS)
             case (1)
                dnHatdpsiHats = 1
                dTHatdpsiHats = 0
                EParallelHat = 0
             case (2)
                dnHatdpsiHats = 0
                dTHatdpsiHats = 0
                EParallelHat = 1
             end select
          end select

          ! To compute the right-hand side, we note the following:
          ! The right-hand side for a linear problem 
          ! is the same as 
          ! (-1) * the residual when f=0.

          !  Set f=0:
          call VecSet(solutionVec, zero, ierr)

          call evaluateResidual(mysnes, solutionVec, residualVec, userContext, ierr)
          ! Multiply the residual by (-1):
          call VecScale(residualVec, -1d+0, ierr)

          if (masterProc) then
             print *,"Beginning the main solve.  This could take a while ..."
          end if

          call PetscTime(time1, ierr)
          if (solveSystem) then
             ! All the magic happens in this next line!
             call KSPSolve(outer_KSP,residualVec, solutionVec, ierr)
          end if

!!$          if (whichRHS==1) then
!!$             print *,"Here comes Mat_for_Jacobian:"
!!$             call MatView(Mat_for_Jacobian, PETSC_VIEWER_STDOUT_WORLD,ierr)
!!$             print *,"Here comes Mat_for_preconditioner:"
!!$             call MatView(Mat_for_preconditioner, PETSC_VIEWER_STDOUT_WORLD,ierr)
!!$          end if

          call PetscTime(time2, ierr)
          if (masterProc) then
             print *,"Done with the main solve.  Time to solve: ", time2-time1, " seconds."
          end if
          call PetscTime(time1, ierr)

          call checkIfKSPConverged(outer_KSP)

          ! Compute flows, fluxes, etc.:
          call diagnostics(outer_KSP, solutionVec, whichRHS)

       end do

    case default
       print *,"Error! Invalid setting for RHSMode."
       stop
    end select



    ! ***********************************************************************
    ! ***********************************************************************
    ! 
    !  End of the main solver code.
    !
    ! ***********************************************************************
    ! ***********************************************************************


    call VecDestroy(solutionVec, ierr)
    call VecDestroy(residualVec, ierr)
    ! A seg fault occurs in nonlinear runs if we call MatDestroy here, since the matrix was already destroyed (and a different Mat created) in evaluateJacobian().
!!$    call MatDestroy(matrix, ierr)
!!$    if (useIterativeLinearSolver) then
!!$       call MatDestroy(Mat_for_preconditioner, ierr)
!!$    end if
    call SNESDestroy(mysnes,ierr)


  end subroutine mainSolverLoop


  ! ------------------------------------------------------------------------


  subroutine chooseParallelDirectSolver()

    implicit none

    character(len=*), parameter :: line="******************************************************************"
    isAParallelDirectSolverInstalled = .false.

    if ((whichParallelSolverToFactorPreconditioner<1) .or. (whichParallelSolverToFactorPreconditioner>2)) then
       print *,"Error! Invalid setting for whichParallelSolverToFactorPreconditioner"
       stop
    end if

#ifdef PETSC_HAVE_MUMPS
    isAParallelDirectSolverInstalled = .true.
    if (masterProc) then
       print *,"mumps detected"
    end if
#else
    whichParallelSolverToFactorPreconditioner = 2
    if (masterProc) then
       print *,"mumps not detected"
    end if
#endif

#ifdef PETSC_HAVE_SUPERLU_DIST
    isAParallelDirectSolverInstalled = .true.
    if (masterProc) then
       print *,"superlu_dist detected"
    end if
#else
    if (masterProc) then
       print *,"superlu_dist not detected"
    end if
    if (whichParallelSolverToFactorPreconditioner==2) then
       whichParallelSolverToFactorPreconditioner = 1
    end if
#endif

    if ((.not. isAParallelDirectSolverInstalled) .and. (numProcs > 1)) then
       if (masterProc) then
          print *,"Error! To run with more than 1 processors, you must have either"
          print *,"mumps or superlu_dist installed."
       end if
       stop
    end if

    if (.not. isAParallelDirectSolverInstalled) then
       if (masterProc) then
          print *,line
          print *,line
          print *,"**   WARNING: As neither mumps nor superlu_dist seems to be available, SFINCS will have to use PETSc's built-in sparse direct solver,"
          print *,"**            which tends to be fragile.  You may get a zero-pivot error, or sometimes this solver simply gives the wrong"
          print *,"**            factorization of the preconditioner matrix, leading to slow KSP convergence and/or wrong outputs of sfincs."
          print *,"**            It is highly recommended that you use instead use a build of PETSc that recognizes mumps or superlu_dist."
          print *,line
          print *,line
       end if
    end if

  end subroutine chooseParallelDirectSolver



end module solver

