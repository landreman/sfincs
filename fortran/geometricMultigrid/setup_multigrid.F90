#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

subroutine setup_multigrid()

  use petscmat
  use petscksp
  use globalVariables

  implicit none

  integer :: j, k, level
  Vec :: temp_vec
  PetscErrorCode :: ierr
  character(len=100) :: filename
  PetscViewer :: viewer
  KSP :: smoother_ksp, ksp_on_coarsest_level
  PC :: smoother_pc, pc_on_coarsest_level
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
  PetscViewerAndFormat vf 
#endif

  external apply_Jacobian
  external apply_preconditioner
  external evaluateJacobian, evaluateResidual, diagnosticsMonitor

  if (masterProc) print *,"Entering setup_multigrid"

  iterationForMatrixOutput = 0

  call PetscOptionsInsertString(PETSC_NULL_OBJECT,"-options_left", ierr) ! This helps for spotting misspelled options!

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
  call KSPSetType(outer_KSP, KSPGMRES, ierr)   ! 20170503 Should this be FGMRES?
  call KSPGMRESSetRestart(outer_KSP, 2000, ierr)
  
  call KSPCreate(MPIComm, inner_ksp, ierr)
  call KSPAppendOptionsPrefix(inner_ksp, 'inner_', ierr)
  call KSPSetType(inner_ksp, KSPPREONLY, ierr)
  call KSPGetPC(inner_ksp, inner_preconditioner, ierr)
  call KSPSetOperators(inner_ksp, Mat_for_preconditioner, Mat_for_preconditioner, ierr)
  call PCSetType(inner_preconditioner, PCMG, ierr)

  call set_grid_resolutions()
  call PCMGSetLevels(inner_preconditioner,N_levels, PETSC_NULL_OBJECT,ierr)

  do j = 1,N_levels
     call create_multigrid_grids(j)
  end do

  call computeB()

  ! Always build the high-order matrix at the finest level:
  call preallocateMatrix(levels(1)%high_order_matrix,1,1)
  call populateMatrix(levels(1)%high_order_matrix,1,1)

  if (defect_option==4) then
     call preallocateMatrix(levels(N_levels)%mixed_order_matrix,10,N_levels)
     call populateMatrix(levels(N_levels)%mixed_order_matrix,10,N_levels)
  end if

  do j = 1,N_levels
    
     if (defect_option > 1 .and. j>1) then
        ! Build the high-order matrix at this level:
        call preallocateMatrix(levels(j)%high_order_matrix,1,j)
        call populateMatrix(levels(j)%high_order_matrix,1,j)
     end if

     ! Build the low-order matrix at this level:
     if (j==N_levels) then
        ! DO include constraint and source at the coarse level, where we do a direct solve.
        call preallocateMatrix(levels(j)%low_order_matrix,0,j)
        call populateMatrix(levels(j)%low_order_matrix,0,j)
     else
        ! DON'T include constraint and source at finer levels, where we do smoothing.
        call preallocateMatrix(levels(j)%low_order_matrix,4,j)
        call populateMatrix(levels(j)%low_order_matrix,4,j)
     end if

  end do

  call KSPSetOperators(main_ksp, levels(1)%high_order_matrix, levels(1)%low_order_matrix, ierr)
  do j = 1,N_levels
     if (defect_option==1) then
        call PCMGSetResidual(inner_preconditioner, N_levels-j, PCMGResidualDefault, levels(j)%low_order_matrix, ierr)
     else
        call PCMGSetResidual(inner_preconditioner, N_levels-j, PCMGResidualDefault, levels(j)%high_order_matrix, ierr)
     end if
  end do

  matrixSize = levels(1)%matrixSize

  allocate(multigrid_prolongation_matrices(N_levels-1))
  allocate(multigrid_restriction_matrices(N_levels-1))

  do j=1,N_levels-1
     call restriction_prolongation_matrices(j)
     ! My level 1 is the finest. PETSc's level 0 is the coarsest.
     call PCMGSetRestriction(  inner_preconditioner,N_levels-j, multigrid_restriction_matrices(j),ierr)
     call PCMGSetInterpolation(inner_preconditioner,N_levels-j,multigrid_prolongation_matrices(j),ierr)
  end do

  ! *****************************************************
  ! Build matrices and vectors needed for smoothing
  ! *****************************************************

  select case (smoothing_option)
  case (0)
     ! Try setting nothing except the KSP operators.
     do level = 1,N_levels-1  ! We don't need to smooth on the coarsest level, where we do a direct solve.
        call PCMGGetSmoother(inner_preconditioner,N_levels-level,smoother_ksp,ierr)
        !call KSPSetOperators(smoother_ksp, levels(level)%smoothing_matrix, levels(level)%smoothing_matrix, ierr)
        call KSPSetOperators(smoother_ksp, levels(level)%low_order_matrix, levels(level)%low_order_matrix, ierr)
        call KSPSetType(smoother_KSP, KSPRICHARDSON, ierr)
     end do
  case (1)
     do level = 1,N_levels-1  ! We don't need to smooth on the coarsest level, where we do a direct solve.
        if (constraint_option==1) then
           call preallocateMatrix(levels(level)%smoothing_matrix,4,level)
           call populateMatrix(levels(level)%smoothing_matrix,4,level)
        else
           levels(level)%smoothing_matrix = levels(level)%low_order_matrix
        end if

        call PCMGGetSmoother(inner_preconditioner,N_levels-level,smoother_ksp,ierr)
        call KSPSetOperators(smoother_ksp, levels(level)%smoothing_matrix, levels(level)%smoothing_matrix, ierr)

        call KSPSetType(smoother_KSP, KSPRICHARDSON, ierr)
        call KSPSetNormType(smoother_KSP, KSP_NORM_NONE, ierr)
        call KSPSetTolerances(smoother_KSP, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, 2, ierr)
        call KSPGetPC(smoother_ksp, smoother_pc, ierr)
        call PCSetType(smoother_pc, PCJACOBI, ierr)

        ! I don't think I need the next line?
        !call KSPSetInitialGuessNonzero(levels(level)%smoothing_ksp, PETSC_TRUE, ierr)

     end do

  case (3)
     do level = 1,N_levels-1  ! We don't need to smooth on the coarsest level, where we do a direct solve.
        if (constraint_option==1) then
           call preallocateMatrix(levels(level)%smoothing_matrix,4,level)
           call populateMatrix(levels(level)%smoothing_matrix,4,level)
        else
           levels(level)%smoothing_matrix = levels(level)%low_order_matrix
        end if

        call PCMGGetSmoother(inner_preconditioner,N_levels-level,smoother_ksp,ierr)
        call KSPSetOperators(smoother_ksp, levels(level)%smoothing_matrix, levels(level)%smoothing_matrix, ierr)

        call KSPSetType(smoother_KSP, KSPRICHARDSON, ierr)
        call KSPSetNormType(smoother_KSP, KSP_NORM_NONE, ierr)
        call KSPSetTolerances(smoother_KSP, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, 2, ierr)
        call KSPGetPC(smoother_ksp, smoother_pc, ierr)
        call PCSetType(smoother_pc, PCSOR, ierr)

        ! I don't think I need the next line?
        !call KSPSetInitialGuessNonzero(levels(level)%smoothing_ksp, PETSC_TRUE, ierr)
     end do
  case default
     print *,"Invalid smoothing_option:",smoothing_option
  end select


  ! Set up the direct solver for the coarsest level
  call PCMGGetCoarseSolve(inner_preconditioner,ksp_on_coarsest_level,ierr)
  if (defect_option==2) then
     if (masterProc) print *,"Using HIGH order matrix for the direct solve on the coarest level."
     call KSPSetOperators(ksp_on_coarsest_level, levels(N_levels)%high_order_matrix, levels(N_levels)%high_order_matrix, ierr)
  elseif (defect_option==4) then
     if (masterProc) print *,"Using MIXED discretization matrix for the direct solve on the coarest level."
     call KSPSetOperators(ksp_on_coarsest_level, levels(N_levels)%mixed_order_matrix, levels(N_levels)%mixed_order_matrix, ierr)
  else
     if (masterProc) print *,"Using LOW order matrix for the direct solve on the coarest level."
     call KSPSetOperators(ksp_on_coarsest_level, levels(N_levels)%low_order_matrix, levels(N_levels)%low_order_matrix, ierr)
  end if
  call KSPGetPC(ksp_on_coarsest_level,pc_on_coarsest_level,ierr)
  call PCSetType(pc_on_coarsest_level,PCLU,ierr)
  call KSPSetType(ksp_on_coarsest_level, KSPPREONLY, ierr)
  call PCFactorSetMatSolverPackage(pc_on_coarsest_level, MATSOLVERMUMPS, ierr)



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

    scalar = solverTolerance
    call KSPSetTolerances(outer_KSP, scalar, PETSC_DEFAULT_REAL, &
         PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)

#if defined(PETSC_USE_REAL_SINGLE)
    ! When using single precision, it sometimes helps to use the gmres setting below that is less susceptible to roundoff error:
    call KSPGMRESSetCGSRefinementType(outer_KSP, KSP_GMRES_CGS_REFINE_IFNEEDED, ierr)
#endif

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

end subroutine setup_multigrid


