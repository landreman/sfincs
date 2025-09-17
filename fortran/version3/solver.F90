module solver

  use globalVariables
  use adjointDiagnostics
  use writeHDF5Output

#include "PETScVersions.F90"

  implicit none

contains

  subroutine mainSolverLoop

    implicit none

    PetscErrorCode :: ierr
    Vec :: solutionVec, residualVec
    Mat :: matrix, preconditionerMatrix
    SNES :: mysnes
    KSP :: KSPInstance, KSPInstance_adjoint
    PC :: preconditionerContext, preconditionerContext_adjoint
    integer :: numRHSs
    SNESConvergedReason :: reason
    integer :: reason_int
    double precision :: time1, time2
    integer :: userContext(1)
    Vec :: dummyVec
    Mat :: factorMat, factorMat_adjoint
    PetscInt :: mumps_which_cntl, mumps_icntl_14 = 50
    PetscReal :: mumps_value
    PetscReal :: atol, rtol, stol, norm
    integer :: maxit, maxf, factor_err

    ! Related to adjoint solve
    integer :: whichAdjointRHS, iSpecies, whichLambda, whichMode, whichSpecies
    Vec :: adjointSolutionVec, summedSolutionVec, adjointRHSVec, adjointSolutionJr
    Mat :: adjointMatrix, adjointPreconditionerMatrix
    logical :: useSummedSolutionVec
    logical :: summedSolution_particleFlux = .false.
    logical :: summedSolution_heatFlux = .false.
    logical :: summedSolution_parallelFlow = .false.
    logical :: solve_particleFlux = .false.
    logical :: solve_heatFlux = .false.
    logical :: solve_parallelFlow = .false.

    MatSolverType :: actualSolverType

    PetscBool :: is_icntl_14_set
    logical :: doAnotherSolve
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
    PetscViewerAndFormat vf
    PetscViewerAndFormat vf_adjoint
#endif

    external evaluateJacobian, evaluateResidual, diagnosticsMonitor

    ! If mumps ICNTL(14) was set via the command line, we need to clear it, or else
    ! the auto-increase feature will not work.
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 8))
    ! Version 3.8
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,'-mat_mumps_icntl_14',mumps_icntl_14,is_icntl_14_set,ierr)
#elif (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 7)
    ! Version 3.7
    call PetscOptionsGetInt(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER,'-mat_mumps_icntl_14',mumps_icntl_14,is_icntl_14_set,ierr)
#else
    ! Version <= 3.6
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-mat_mumps_icntl_14',mumps_icntl_14,is_icntl_14_set,ierr)
#endif
    if (.not. is_icntl_14_set) then
       ! Default:
       mumps_icntl_14 = 50
    else
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 8))
       ! Verion 3.8
       call PetscOptionsClearValue(PETSC_NULL_OPTIONS,'-mat_mumps_icntl_14',ierr)
#elif (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 7)
       ! Version 3.7
       call PetscOptionsClearValue(PETSC_NULL_OBJECT,'-mat_mumps_icntl_14',ierr)
#else
       ! Versions <= 3.6
       call PetscOptionsClearValue('-mat_mumps_icntl_14',ierr)
#endif
    end if

    doAnotherSolve = .true.
    do while (doAnotherSolve)
       doAnotherSolve = .false.

       if (masterProc) then
          print *,"Entering main solver loop."
       end if
       iterationForMatrixOutput = 0

       call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, solutionVec, ierr)
       call VecDuplicate(solutionVec, residualVec, ierr)

       call VecDuplicate(solutionVec, f0, ierr)
       call init_f0()

       call SNESCreate(MPIComm, mysnes, ierr)
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 8))
#define PETSC_NULL_OBJECT_OR_INTEGER 0
#else
#define PETSC_NULL_OBJECT_OR_INTEGER PETSC_NULL_OBJECT
#endif

       call SNESSetFunction(mysnes, residualVec, evaluateResidual, PETSC_NULL_OBJECT_OR_INTEGER, ierr)

       firstMatrixCreation = .true.
       call preallocateMatrix(matrix, 1)
       if (useIterativeLinearSolver) then
          call preallocateMatrix(preconditionerMatrix, 0)
          call SNESSetJacobian(mysnes, matrix, preconditionerMatrix, evaluateJacobian, PETSC_NULL_OBJECT_OR_INTEGER, ierr)
       else
          call SNESSetJacobian(mysnes, matrix, matrix, evaluateJacobian, PETSC_NULL_OBJECT_OR_INTEGER, ierr)
       end if

       call SNESGetKSP(mysnes, KSPInstance, ierr)
       if (discreteAdjointOption .eqv. .false.) then
          call KSPCreate(MPIComm,KSPInstance_adjoint, ierr)
          call preallocateMatrix(adjointMatrix,1)
          call populateMatrix(adjointMatrix,4,dummyVec)
          if (useIterativeLinearSolver) then
            call preallocateMatrix(adjointPreconditionerMatrix, 1)
            call populateMatrix(adjointPreconditionerMatrix,5,dummyVec)
          end if
       end if

       if (useIterativeLinearSolver) then
!!$#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
!!$       ! Syntax for PETSc versions up through 3.4
!!$       call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, SAME_PRECONDITIONER, ierr)
!!$#else
!!$       ! Syntax for PETSc version 3.5 and later
!!$       call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, ierr)
!!$#endif
          if (discreteAdjointOption .eqv. .false.) then
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
         ! Syntax for PETSc versions up through 3.4
             call KSPSetOperators(KSPInstance_adjoint, adjointMatrix, adjointPreconditionerMatrix, SAME_PRECONDITIONER, ierr)
#else
        ! Syntax for PETSc version 3.5 and later
             call KSPSetOperators(KSPInstance_adjoint, adjointMatrix, adjointPreconditionerMatrix, ierr)
             call KSPSetReusePreconditioner(KSPInstance_adjoint, PETSC_TRUE, ierr)
#endif
          end if
          call KSPGetPC(KSPInstance, preconditionerContext, ierr)
          call PCSetType(preconditionerContext, PCLU, ierr)
          call KSPSetType(KSPInstance, KSPGMRES, ierr)   ! Set the Krylov solver algorithm to GMRES
          call KSPGMRESSetRestart(KSPInstance, 2000, ierr)
          call KSPSetTolerances(KSPInstance, solverTolerance, PETSC_DEFAULT_REAL, &
               PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)

          ! Allow options to be controlled using command-line flags:
          call KSPSetFromOptions(KSPInstance, ierr)
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 7))
          call KSPMonitorSet(KSPInstance, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
#else
          call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr)
#if ((PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 14))
          call KSPMonitorSet(KSPInstance, KSPMonitorResidual, vf, PetscViewerAndFormatDestroy, ierr)
#else
          call KSPMonitorSet(KSPInstance, KSPMonitorDefault, vf, PetscViewerAndFormatDestroy, ierr)
#endif
#endif

        if (discreteAdjointOption .eqv. .false.) then
            call KSPGetPC(KSPInstance_adjoint, preconditionerContext_adjoint, ierr)
            call PCSetType(preconditionerContext_adjoint, PCLU, ierr)
            call KSPSetType(KSPInstance_adjoint, KSPGMRES, ierr)   ! Set the Krylov solver algorithm to GMRES
            call KSPGMRESSetRestart(KSPInstance_adjoint, 2000, ierr)
            call KSPSetTolerances(KSPInstance_adjoint, solverTolerance, PETSC_DEFAULT_REAL, &
                 PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)

            ! Allow options to be controlled using command-line flags:
            call KSPSetFromOptions(KSPInstance_adjoint, ierr)
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 7))
            call KSPMonitorSet(KSPInstance_adjoint, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
#else
            call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf_adjoint, ierr)
#if ((PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 14))
            call KSPMonitorSet(KSPInstance_adjoint, KSPMonitorResidual, vf_adjoint, PetscViewerAndFormatDestroy, ierr)
#else
            call KSPMonitorSet(KSPInstance_adjoint, KSPMonitorDefault, vf_adjoint, PetscViewerAndFormatDestroy, ierr)
#endif
            call PCSetReusePreconditioner(preconditionerContext_adjoint, PETSC_TRUE, ierr)
#endif
          end if

       else
          ! Direct solver:
!!$#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
!!$       ! Syntax for PETSc versions up through 3.4
!!$       call KSPSetOperators(KSPInstance, matrix, matrix, SAME_PRECONDITIONER, ierr)
!!$#else
!!$       ! Syntax for PETSc version 3.5 and later
!!$       call KSPSetOperators(KSPInstance, matrix, matrix, ierr)
!!$#endif
        if (discreteAdjointOption .eqv. .false.) then
          ! Direct solver:
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
       ! Syntax for PETSc versions up through 3.4
       call KSPSetOperators(KSPInstance_adjoint, adjointMatrix, adjointMatrix, SAME_PRECONDITIONER, ierr)
#else
       ! Syntax for PETSc version 3.5 and later
       call KSPSetOperators(KSPInstance_adjoint, adjointMatrix, adjointMatrix, ierr)
#endif
       end if
          call KSPGetPC(KSPInstance, preconditionerContext, ierr)
          call PCSetType(preconditionerContext, PCLU, ierr)
          call KSPSetType(KSPInstance, KSPPREONLY, ierr)
          ! Allow options to be controlled using command-line flags:
          call KSPSetFromOptions(KSPInstance, ierr)
          if (discreteAdjointOption .eqv. .false.) then
            call KSPGetPC(KSPInstance_adjoint, preconditionerContext_adjoint, ierr)
            call PCSetType(preconditionerContext_adjoint, PCLU, ierr)
            call KSPSetType(KSPInstance_adjoint, KSPPREONLY, ierr)
!             Allow options to be controlled using command-line flags:
            call KSPSetFromOptions(KSPInstance_adjoint, ierr)
          end if

       end if

       if (isAParallelDirectSolverInstalled) then
          select case (whichParallelSolverToFactorPreconditioner)
          case (1)
             call PCFactorSetMatSolverType(preconditionerContext, MATSOLVERMUMPS, ierr)
           if (discreteAdjointOption .eqv. .false.) then
              call PCFactorSetMatSolverType(preconditionerContext_adjoint, MATSOLVERMUMPS, ierr)
           end if
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
             ! The functions MatMumpsSetICNTL were introduced in PETSc 3.5.
             ! For earlier versions, we can achieve a similar result with the following hack:
             call PetscOptionsInsertString("-mat_mumps_cntl_1 1e-6 -mat_mumps_icntl_4 2", ierr)
#endif
          case (2)
             call PCFactorSetMatSolverType(preconditionerContext, MATSOLVERSUPERLU_DIST, ierr)
             if (discreteAdjointOption .eqv. .false.) then
                call PCFactorSetMatSolverType(preconditionerContext_adjoint, MATSOLVERSUPERLU_DIST, ierr)
             end if
             ! Turn on superlu_dist diagnostic output:
             !call PetscOptionsInsertString("-mat_superlu_dist_statprint", ierr)
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
          call PCFactorSetMatOrderingType(preconditionerContext, MATORDERINGRCM, ierr)

          call PCFactorReorderForNonzeroDiagonal(preconditionerContext, 1d-12, ierr) 

          call PCFactorSetZeroPivot(preconditionerContext, 1d-200, ierr)
          if (discreteAdjointOption .eqv. .false.) then
            call PCFactorSetMatOrderingType(preconditionerContext_adjoint, MATORDERINGRCM, ierr)

            call PCFactorReorderForNonzeroDiagonal(preconditionerContext_adjoint, 1d-12, ierr)

            call PCFactorSetZeroPivot(preconditionerContext_adjoint, 1d-200, ierr)
          end if
       end if


       ! Tell PETSc to call the diagnostics subroutine at each iteration of SNES:
       ! 20180525 uncomment the next line!
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 8))
       call SNESMonitorSet(mysnes, diagnosticsMonitor, 0, PETSC_NULL_FUNCTION, ierr)
#else
       call SNESMonitorSet(mysnes, diagnosticsMonitor, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
#endif

       if (reusePreconditioner .or. (discreteAdjointOption)) then
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
          ! Syntax for PETSc versions up through 3.4
          ! In this case the associated code appears in evaluateJacobian.F90
#else
          ! Syntax for PETSc version 3.5 and later
          call KSPSetReusePreconditioner(KSPInstance, PETSC_TRUE, ierr)
          call PCSetReusePreconditioner(preconditionerContext, PETSC_TRUE, ierr)
#endif
       end if
      if (discreteAdjointOption .eqv. .false.) then
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
! Syntax for PETSc versions up through 3.4
! In this case the associated code appears in evaluateJacobian.F90
#else
      ! Syntax for PETSc version 3.5 and later
      call KSPSetReusePreconditioner(KSPInstance_adjoint, PETSC_TRUE, ierr)
      call PCSetReusePreconditioner(preconditionerContext_adjoint, PETSC_TRUE, ierr)
#endif
      end if

       ! In older versions of PETSC (either <3.5 or <3.4, I'm not certain)
       ! the monitor is never called when snes type = SNESKSPONLY.
       ! Therefore, it is preferable to always have snes type = SNESNEWTONLS but set the 
       ! number of iterations to 1 for a linear run.
       call SNESSetType(mysnes, SNESNEWTONLS, ierr)
       !!if (nonlinear) then !!Commented by AM 2016-02
       !!if (includePhi1) then !!Added by AM 2016-02 !!Commented by AM 2018-12
       if (includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
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
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 9))
       ! Version 3.9
       call PCFactorGetMatSolverType(preconditionerContext, actualSolverType, ierr)
#else
       call PCFactorGetMatSolverType(preconditionerContext, actualSolverType, ierr)
#endif
       if (masterProc) then
          print *,"Solver package which will be used: ",actualSolverType
       end if
      if (discreteAdjointOption .eqv. .false.) then
          call PCFactorGetMatSolverType(preconditionerContext_adjoint, actualSolverType, ierr)
          if (masterProc) then
            print *,"Solver package which will be used for adjoint: ",actualSolverType
          end if
       end if


       !if (isAParallelDirectSolverInstalled .and. (whichParallelSolverToFactorPreconditioner==1)) then
       if (actualSolverType == MATSOLVERMUMPS) then
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
          call PCFactorSetUpMatSolverType(preconditionerContext,ierr)

          call PCFactorGetMatrix(preconditionerContext,factorMat,ierr)

          if (discreteAdjointOption .eqv. .false.) then
            call PCFactorSetUpMatSolverType(preconditionerContext_adjoint,ierr)
            call PCFactorGetMatrix(preconditionerContext_adjoint,factorMat_adjoint,ierr)
          end if

          ! All options set below can be over-ridden by command-line arguments, even though
          ! SNESSetFromOptions was called above rather than below.

          mumps_which_cntl = 1
          mumps_value = 1.e-6
          call MatMumpsSetCntl(factorMat,mumps_which_cntl,mumps_value,ierr)
          if (discreteAdjointOption .eqv. .false.) then
            call MatMumpsSetCntl(factorMat_adjoint,mumps_which_cntl,mumps_value,ierr)
          end if

          ! Turn on mumps diagnostic output:
          mumps_which_cntl = 4
          call MatMumpsSetIcntl(factorMat,mumps_which_cntl,2,ierr)
          if (discreteAdjointOption .eqv. .false.) then
            call MatMumpsSetIcntl(factorMat_adjoint,mumps_which_cntl,2,ierr)
          end if

!!$       if (numProcs>1) then
!!$          ! Turn on parallel ordering by default.
!!$          ! Only do this if >1 procs, to avoid confusing warning message otherwise.
!!$          mumps_which_cntl = 28
!!$          call MatMumpsSetIcntl(factorMat,mumps_which_cntl,2,ierr)
!!$       end if

          ! Increase amount by which the mumps work array can expand due to near-0 pivots.
          ! (The default value of icntl(14) is 25.)
          ! For many cases it is not necessary to increase icntl(14), but an increase sometimes helps
          ! to prevent the mumps error with info(1)=-9.  There appears to be basically
          ! no significant cost in memory or time to increase this parameter.
          mumps_which_cntl = 14
          call MatMumpsSetIcntl(factorMat,mumps_which_cntl,mumps_icntl_14,ierr)
          if (discreteAdjointOption .eqv. .false.) then
          call MatMumpsSetIcntl(factorMat_adjoint,mumps_which_cntl,mumps_icntl_14,ierr)
          end if
#endif
       end if

       ! ***********************************************************************
       ! ***********************************************************************
       ! 
       !  Main solver call:
       !
       ! ***********************************************************************
       ! ***********************************************************************

       select case (RHSMode)
       case (1,4,5)
          ! Single solve, either linear or nonlinear.
          ! In case of adjoint solve, perform forward solve first.


          !  Set initial guess for the solution:
          call VecSet(solutionVec, zero, ierr)

          if (masterProc) then
             print *,"------------------------------------------------------"
             print *,"Finished initialization."
             print *,"Beginning the main solve.  This could take a while ..."
          end if

          time1 = MPI_Wtime()
          if (solveSystem) then
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 8))
             call SNESSolve(mysnes,PETSC_NULL_VEC, solutionVec, ierr)
#else
             call SNESSolve(mysnes,PETSC_NULL_OBJECT, solutionVec, ierr)
#endif
             call MatMumpsGetInfog(factorMat, 1, factor_err, ierr)

#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
             ! No way to automatically control MUMPS for old versions of PETSC.
#else
             if (actualSolverType == MATSOLVERMUMPS .and. factor_err .ne. 0 .and. mumps_icntl_14<1024) then
                ! For now, assume all failures are due to MUMPS.  This might not always be true....
                ! Try increasing the amount by which the mumps work array can expand due to near-0 pivots.
                if (masterProc) then
                   print *,"Mumps INFOG(1) = ",factor_err
                   print *,"Solve failed, so doubling MUMPS icntl(14), which had been ",mumps_icntl_14
                end if
                mumps_icntl_14 = mumps_icntl_14 * 2
                ierr = 0
                doAnotherSolve = .true.
             end if
#endif
          end if

          time2 = MPI_Wtime()
          if (masterProc) then
             print *,"Done with the main solve.  Time to solve: ", time2-time1, " seconds."
          end if
          time1 = MPI_Wtime()

          !!if (includePhi1) then !!Commented by AM 2018-12
          if (includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
             call SNESGetConvergedReason(mysnes, reason, ierr)
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 23))
             reason_int = reason
#else
             reason_int = reason%v
#endif   
             if (reason_int > 0) then
                if (masterProc) then
                   print *,"Nonlinear iteration (SNES) converged!  SNESConvergedReason = ", reason_int
                   select case (reason_int)
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
                   print *,"Nonlinear iteration (SNES) did not converge :(   SNESConvergedReason = ", reason_int
                   select case (reason_int)
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

          ! Build the main linear system matrix:
          call populateMatrix(matrix,1,dummyVec)

          if (useIterativeLinearSolver) then

             ! Build the preconditioner:
             call populateMatrix(preconditionerMatrix,0,dummyVec)

#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
             ! Syntax for PETSc versions up through 3.4
             call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, SAME_PRECONDITIONER, ierr)
#else
             ! Syntax for PETSc version 3.5 and later
             call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, ierr)
             call KSPSetReusePreconditioner(KSPInstance, PETSC_TRUE, ierr)
             call PCSetReusePreconditioner(preconditionerContext, PETSC_TRUE, ierr)
#endif

          else
             ! Direct solver, so no preconditioner needed.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
             ! Syntax for PETSc versions up through 3.4
             call KSPSetOperators(KSPInstance, matrix, matrix, SAME_PRECONDITIONER, ierr)
#else
             ! Syntax for PETSc version 3.5 and later
             call KSPSetOperators(KSPInstance, matrix, matrix, ierr)
             call KSPSetReusePreconditioner(KSPInstance, PETSC_TRUE, ierr)
             call PCSetReusePreconditioner(preconditionerContext, PETSC_TRUE, ierr)
#endif
          end if

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

             time1 = MPI_Wtime()
             if (solveSystem) then
                ! All the magic happens in this next line!
                call KSPSolve(KSPInstance,residualVec, solutionVec, ierr)
             end if

             time2 = MPI_Wtime()
             if (masterProc) then
                print *,"Done with the main solve.  Time to solve: ", time2-time1, " seconds."
             end if
             time1 = MPI_Wtime()

             call checkIfKSPConverged(KSPInstance)

             ! Compute flows, fluxes, etc.:
             call diagnostics(solutionVec, whichRHS)

          end do

       case default
          print *,"Error! Invalid setting for RHSMode."
          stop
       end select

! Initialize things needed for adjoint solve
    if (RHSMode>3 .or. (ambipolarSolve .and. (ambipolarSolveOption==1 .or. ambipolarSolveOption==3))) then

      ! Forward matrix & preconditioner no longer needed
      if (discreteAdjointOption .eqv. .false.) then
        call VecDestroy(residualVec,ierr)
        if (useIterativeLinearSolver) then
          call MatDestroy(preconditionerMatrix,ierr)
        end if
        call MatDestroy(matrix,ierr)
        call SNESDestroy(mysnes,ierr)
      end if

      !> Allocate adjointRHSVec
      call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointRHSVec, ierr)
      call VecSet(adjointRHSVec, zero, ierr)

      if (masterProc) then
         print *,"Finished allocating adjoint Vecs and matrices."
      end if

    end if

    ! Adjoint solve for Newton method
    if (ambipolarSolve .and. (ambipolarSolveOption == 1 .or. ambipolarSolveOption==3) .or. RHSMode==5) then
      whichAdjointRHS = 1
      ispecies = 0
      if (masterProc) then
        print *,"################################################################"
        print "(a,i1,a,i1,a)"," Solving adjoint system with adjoint RHS ",whichAdjointRHS," and species ",ispecies," for ambipolar solve."
        print *,"################################################################"
      end if

      !> Construct RHS vec
      call VecSet(adjointRHSVec, zero, ierr)
      call populateAdjointRHS(adjointRHSVec, whichAdjointRHS, ispecies, 0)

      !> Construct solution vec
      call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointSolutionJr, ierr)
      call VecSet(adjointSolutionJr, zero, ierr)

      if (masterProc) then
         print *,"Beginning the adjoint solve.  This could take a while ..."
      end if

      time1 = MPI_Wtime()
      if (solveSystem) then
        if (discreteAdjointOption .eqv. .false.) then
          call KSPSolve(KSPInstance_adjoint,adjointRHSVec,adjointSolutionJr, ierr)
          call checkIfKSPConverged(KSPInstance_adjoint)
        else
          call KSPSolveTranspose(KSPInstance,adjointRHSVec,adjointSolutionJr, ierr)
          call checkIfKSPConverged(KSPInstance)
        end if
      end if

      time2 = MPI_Wtime()
      if (masterProc) then
         print *,"Done with the adjoint solve.  Time to solve: ", time2-time1, " seconds."
      end if
      time1 = MPI_Wtime()

      if (ambipolarSolve .and. (ambipolarSolveOption == 1 .or. ambipolarSolveOption == 3) .and. RHSMode < 3) then
        call computedRadialCurrentdEr(solutionVec,adjointSolutionJr)
        call VecDestroy(adjointRHSVec, ierr)
        if (discreteAdjointOption .eqv. .false.) then
          call MatDestroy(adjointMatrix, ierr)
        end if
        call VecDestroy(adjointSolutionJr,ierr)
      end if
    end if

    !> This is where all the adjoint solves happen
    !> If currently looking for ambipolar solution, don't solve adjoint
    if (RHSMode>3 .and. (ambipolarSolve .eqv. .false.)) then

      !> Allocate adjointSolutionVec
      call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointSolutionVec, ierr)
      call VecSet(adjointSolutionVec, zero, ierr)

        if (all(adjointHeatFluxOption) .and. adjointTotalHeatFluxOption) then
          summedSolution_heatFlux = .true.
        end if
        if (all(adjointParticleFluxOption) .and. adjointRadialCurrentOption) then
          summedSolution_particleFlux = .true.
        end if
        if (all(adjointParallelFlowOption) .and. adjointBootstrapOption) then
          summedSolution_parallelFlow = .true.
        end if
        if (any(adjointHeatFluxOption)) then
          solve_heatFlux = .true.
        end if
        if (any(adjointParticleFluxOption)) then
          solve_particleFlux = .true.
        end if
        if (any(adjointParallelFlowOption)) then
          solve_parallelFlow = .true.
        end if
      if (summedSolution_heatFlux .or. summedSolution_particleFlux .or. summedSolution_parallelFlow) then
        call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, summedSolutionVec, ierr)
        call VecSet(summedSolutionVec, zero, ierr)
      end if

      !> First,we'll compute the species-specific fluxes if needed
      do whichAdjointRHS=1,3
        useSummedSolutionVec = .false.

        !> useSummedSolutionVec = .true. if all fluxes computed for all species and sensitivity of total flux also computed
        select case (whichAdjointRHS)
          case (1) ! particle flux
            if (summedSolution_particleFlux) then
              useSummedSolutionVec = .true.
            end if
            if (.not. solve_particleFlux) then
              cycle
            end if
          case (2) ! heat flux
            if (summedSolution_heatFlux) then
              useSummedSolutionVec = .true.
            end if
            if (.not. solve_heatFlux) then
              cycle
            end if
          case (3) ! parallel flow
            if (summedSolution_parallelFlow) then
              useSummedSolutionVec = .true.
            end if
            if (.not. solve_parallelFlow) then
              cycle
            end if
        end select
        do ispecies=1,Nspecies

          if (masterProc) then
            print *,"ispecies = ", ispecies, ", whichAdjointRHS = ", whichAdjointRHS
          end if

          !> Check if adjoint equation needs to be solved for this species
          select case (whichAdjointRHS)
            case (1) ! particle flux
                if (.not. adjointParticleFluxOption(ispecies)) then
                  cycle
                end if
            case (2) ! heat flux
                if (.not. adjointHeatFluxOption(ispecies)) then
                  cycle
                end if
            case (3) ! parallel flow
                if (.not. adjointParallelFlowOption(ispecies)) then
                  cycle
                end if
          end select

          if (masterProc) then
            print *,"################################################################"
            print "(a,i1,a,i1)"," Solving adjoint system with adjoint RHS ",whichAdjointRHS," and species ",ispecies
            print *,"################################################################"
          end if

          !> Construct RHS vec
          call VecSet(adjointRHSVec, zero, ierr)
          call populateAdjointRHS(adjointRHSVec, whichAdjointRHS, ispecies, 0)

          if (masterProc) then
             print *,"Beginning the adjoint solve.  This could take a while ..."
          end if

          time1 = MPI_Wtime()
          if (solveSystem) then
            if (discreteAdjointOption .eqv. .false.) then
              call KSPSolve(KSPInstance_adjoint,adjointRHSVec,adjointSolutionVec, ierr)
              call checkIfKSPConverged(KSPInstance_adjoint)
            else
              call KSPSolveTranspose(KSPInstance,adjointRHSVec,adjointSolutionVec, ierr)
              call checkIfKSPConverged(KSPInstance)
            end if
          end if

          time2 = MPI_Wtime()
          if (masterProc) then
             print *,"Done with the adjoint solve.  Time to solve: ", time2-time1, " seconds."
          end if

          !> Compute diagnostics for species-specific fluxes
          time1 = MPI_Wtime()
          call evaluateDiagnostics(solutionVec, adjointSolutionVec,adjointSolutionJr,whichAdjointRHS,ispecies)
          time2 = MPI_Wtime()
          if (masterProc) then
            print *,"Done with the adjoint diagnostics.  Time: ", time2-time1, " seconds."
          end if

          select case (whichAdjointRHS)
            case (1) ! particle flux
              if (useSummedSolutionVec) then
                call VecAXPY(summedSolutionVec,Zs(ispecies),adjointSolutionVec,ierr)
              end if
            case (2) ! heat flux
              if (useSummedSolutionVec) then
                call VecAXPY(summedSolutionVec,one,adjointSolutionVec,ierr)
              end if
            case (3) ! parallel flow
              if (useSummedSolutionVec) then
                call VecAXPY(summedSolutionVec,Zs(ispecies),adjointSolutionVec,ierr)
              end if
          end select

          !> Done with required adjoint solve and diagnostics. Now clear solutionVec
          call VecSet(adjointSolutionVec, zero, ierr)
        end do ! ispecies

        ! Now compute diagnostics for species-summed quantities
        if (useSummedSolutionVec) then
            time1 = MPI_Wtime()
            call evaluateDiagnostics(solutionVec,summedSolutionVec,adjointSolutionJr,whichAdjointRHS,0)
            time2 = MPI_Wtime()
            if (masterProc) then
              print *,"Done with the adjoint diagnostics.  Time: ", time2-time1, " seconds."
            end if
          call VecSet(summedSolutionVec,zero,ierr)
        end if

      end do ! whichAdjointRHS

      !> Now, we'll compute the sensitivity of the species-summed fluxes if not already computed
      if (adjointTotalHeatFluxOption .or. adjointBootstrapOption .or. adjointRadialCurrentOption) then
        ispecies = 0 ! 0 denotes species-summed RHS

        do whichAdjointRHS=1,3
          select case(whichAdjointRHS)
            case (1) ! particle flux
              if (summedSolution_particleFlux .or. (RHSMode==5)) then
                ! species-summed flux computed already
                ! For RHSMode = 5, derivatives computed at constant radial current, so
                ! it's sensitivity is 0
                cycle
              end if
                if (adjointRadialCurrentOption .eqv. .false.) then
                  cycle
                end if
            case (2) ! heat flux
              if (summedSolution_heatFlux) then
                ! species-summed flux computed already
                cycle
              end if
              if (adjointTotalHeatFluxOption .eqv. .false.) then
                cycle
              end if
            case (3) ! bootstrap
              if (summedSolution_parallelFlow) then
                cycle
              end if
              if (adjointBootstrapOption .eqv. .false.) then
                cycle
              end if
          end select

          if (masterProc) then
            print *,"################################################################"
            print "(a,i1,a,i1)"," Solving adjoint system with adjoint RHS ",whichAdjointRHS," and species ",ispecies
            print *,"################################################################"
          end if

          ! Construct RHS vec
          call populateAdjointRHS(adjointRHSVec, whichAdjointRHS, ispecies, 0)

          if (masterProc) then
             print *,"Beginning the adjoint solve.  This could take a while ..."
          end if

          time1 = MPI_Wtime()
          if (solveSystem) then
             if (discreteAdjointOption .eqv. .false.) then
             ! All the magic happens in this next line!
                call KSPSolve(KSPInstance_adjoint,adjointRHSVec,adjointSolutionVec, ierr)
                call checkIfKSPConverged(KSPInstance_adjoint)
             else
                call KSPSolveTranspose(KSPInstance,adjointRHSVec,adjointSolutionVec, ierr)
                call checkIfKSPConverged(KSPInstance)
             end if
          end if
          time2 = MPI_Wtime()
          if (masterProc) then
             print *,"Done with the adjoint solve.  Time to solve: ", time2-time1, " seconds."
          end if

          time1 = MPI_Wtime()
          call evaluateDiagnostics(solutionVec, adjointSolutionVec, adjointSolutionJr, whichAdjointRHS, ispecies)
          time2 = MPI_Wtime()
          if (masterProc) then
            print *,"Done with the adjoint diagnostics.  Time: ", time2-time1, " seconds."
          end if

          ! Done with required adjoint solve and diagnostics. Now clear solutionVec
          call VecSet(adjointSolutionVec, zero, ierr)
        end do ! whichAdjointRHS
      end if

      ! Now deallocate things
      call VecDestroy(adjointSolutionVec, ierr)
      if (summedSolution_particleFlux .or. summedSolution_heatFlux .or. summedSolution_parallelFlow) then
        call VecDestroy(summedSolutionVec, ierr)
      end if
      call VecDestroy(adjointRHSVec, ierr)
      if (discreteAdjointOption .eqv. .false.) then
        call MatDestroy(adjointMatrix, ierr)
      end if
      if (RHSMode==5) then
        call VecDestroy(adjointSolutionJr,ierr)
      end if

      ! Update HDF5 - this was not done in diagnostics()
      ! In debug mode, updateOutputFile is called from testingAdjointDiagnostics
      if (debugAdjoint .eqv. .false.) then
        call updateOutputFile(1, .false.)
      end if

    end if ! RHSMode < 4

       ! ***********************************************************************
       ! ***********************************************************************
       ! 
       !  End of the main solver code.
       !
       ! ***********************************************************************
       ! ***********************************************************************


       call VecDestroy(solutionVec, ierr)
       call VecDestroy(residualVec, ierr)
       if (RHSMode>3) then
         call VecDestroy(f0,ierr)
       end if

       ! A seg fault occurs in nonlinear runs if we call MatDestroy or KSPDestroy here, since the matrix was already destroyed (and a different Mat created) in evaluateJacobian().
!!$    call MatDestroy(matrix, ierr)
!!$    if (useIterativeLinearSolver) then
!!$       call MatDestroy(preconditionerMatrix, ierr)
!!$    end if
    if ((discreteAdjointOption .eqv. .false.) .and. RHSMode>3) then
        call KSPDestroy(KSPInstance_adjoint,ierr)
    end if

    if ((includePhi1 .eqv. .false.) .and. ((discreteAdjointOption.eqv. .false.) .and. RHSMode>3) .eqv. .false.) then
        call MatDestroy(matrix, ierr)
        if (useIterativeLinearSolver) then
         call MatDestroy(preconditionerMatrix, ierr)
        end if
    end if

    call SNESDestroy(mysnes,ierr)

    end do

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


