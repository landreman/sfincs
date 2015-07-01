#include <finclude/petscsnesdef.h>
#include "PETScVersions.F90"

module solver

  use petscsnesdef
  use globalVariables

  implicit none

  contains

  subroutine mainSolverLoop

    implicit none

    PetscErrorCode :: ierr
    Vec :: solutionVec, residualVec
    Mat :: matrix, preconditionerMatrix
    SNES :: mysnes
    KSP :: KSPInstance
    PC :: preconditionerContext
    integer :: numRHSs
    SNESConvergedReason :: reason
    KSPConvergedReason :: KSPReason
    PetscLogDouble :: time1, time2
    integer :: userContext(1)
    Vec :: dummyVec
    Mat :: factorMat
    PetscInt :: mumps_which_cntl
    PetscReal :: mumps_value
    PetscReal :: atol, rtol, stol
    integer :: maxit, maxf

    external evaluateJacobian, evaluateResidual, diagnosticsMonitor

    if (masterProc) then
       print *,"Entering main solver loop."
    end if
    iterationForMatrixOutput = 0

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, solutionVec, ierr)
    call VecDuplicate(solutionVec, residualVec, ierr)

    call VecDuplicate(solutionVec, f0, ierr)
    call init_f0()

    call SNESCreate(MPIComm, mysnes, ierr)
    call SNESSetFunction(mysnes, residualVec, evaluateResidual, PETSC_NULL_OBJECT, ierr)

    firstMatrixCreation = .true.
    call preallocateMatrix(matrix, 1)
    if (useIterativeLinearSolver) then
       call preallocateMatrix(preconditionerMatrix, 0)
       call SNESSetJacobian(mysnes, matrix, preconditionerMatrix, evaluateJacobian, PETSC_NULL_OBJECT, ierr)
    else
       call SNESSetJacobian(mysnes, matrix, matrix, evaluateJacobian, PETSC_NULL_OBJECT, ierr)
    end if

    call SNESGetKSP(mysnes, KSPInstance, ierr)

    if (useIterativeLinearSolver) then
!!$#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
!!$       ! Syntax for PETSc versions up through 3.4
!!$       call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, SAME_PRECONDITIONER, ierr)
!!$#else
!!$       ! Syntax for PETSc version 3.5 and later
!!$       call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, ierr)
!!$#endif
       call KSPGetPC(KSPInstance, preconditionerContext, ierr)
       call PCSetType(preconditionerContext, PCLU, ierr)
       call KSPSetType(KSPInstance, KSPGMRES, ierr)   ! Set the Krylov solver algorithm to GMRES
       call KSPGMRESSetRestart(KSPInstance, 2000, ierr)
       call KSPSetTolerances(KSPInstance, solverTolerance, PETSC_DEFAULT_REAL, &
            PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)

       ! Allow options to be controlled using command-line flags:
       call KSPSetFromOptions(KSPInstance, ierr)
       call KSPMonitorSet(KSPInstance, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)

    else
       ! Direct solver:
!!$#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
!!$       ! Syntax for PETSc versions up through 3.4
!!$       call KSPSetOperators(KSPInstance, matrix, matrix, SAME_PRECONDITIONER, ierr)
!!$#else
!!$       ! Syntax for PETSc version 3.5 and later
!!$       call KSPSetOperators(KSPInstance, matrix, matrix, ierr)
!!$#endif
       call KSPGetPC(KSPInstance, preconditionerContext, ierr)
       call PCSetType(preconditionerContext, PCLU, ierr)
       call KSPSetType(KSPInstance, KSPPREONLY, ierr)
       ! Allow options to be controlled using command-line flags:
       call KSPSetFromOptions(KSPInstance, ierr)
    end if

    if (isAParallelDirectSolverInstalled) then
       select case (whichParallelSolverToFactorPreconditioner)
       case (1)
          call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERMUMPS, ierr)
          if (masterProc) then
             print *,"We will use mumps to factorize the preconditioner."
          end if
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
          ! The functions MatMumpsSetICNTL were introduced in PETSc 3.5.
          ! For earlier versions, we can achieve a similar result with the following hack:
          call PetscOptionsInsertString("-mat_mumps_cntl_1 1e-6 -mat_mumps_icntl_4 2", ierr)
#endif
       case (2)
          call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERSUPERLU_DIST, ierr)
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
       call PCFactorSetMatOrderingType(preconditionerContext, MATORDERINGRCM, ierr)

       call PCFactorReorderForNonzeroDiagonal(preconditionerContext, 1d-12, ierr) 

       call PCFactorSetZeroPivot(preconditionerContext, 1d-200, ierr) 
    end if

    
    ! Tell PETSc to call the diagnostics subroutine at each iteration of SNES:
    call SNESMonitorSet(mysnes, diagnosticsMonitor, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)

    if (reusePreconditioner) then
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
       ! Syntax for PETSc versions up through 3.4
       ! In this case the associated code appears in evaluateJacobian.F90
#else
       ! Syntax for PETSc version 3.5 and later
       call KSPSetReusePreconditioner(KSPInstance, PETSC_TRUE, ierr)
       call PCSetReusePreconditioner(preconditionerContext, PETSC_TRUE, ierr)
#endif
    end if

    ! In older versions of PETSC (either <3.5 or <3.4, I'm not certain)
    ! the monitor is never called when snes type = SNESKSPONLY.
    ! Therefore, it is preferable to always have snes type = SNESNEWTONLS but set the 
    ! number of iterations to 1 for a linear run.
    call SNESSetType(mysnes, SNESNEWTONLS, ierr)
    if (nonlinear) then
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

    if (isAParallelDirectSolverInstalled .and. (whichParallelSolverToFactorPreconditioner==1)) then
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
       call PCFactorSetUpMatSolverPackage(preconditionerContext,ierr)
       call PCFactorGetMatrix(preconditionerContext,factorMat,ierr)
       mumps_which_cntl = 1
       mumps_value = 1.e-6
       call MatMumpsSetCntl(factorMat,mumps_which_cntl,mumps_value,ierr)

       ! Turn on mumps diagnostic output
       mumps_which_cntl = 4
       call MatMumpsSetIcntl(factorMat,mumps_which_cntl,2,ierr)
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

       if (useIterativeLinearSolver) then
          call SNESGetConvergedReason(mysnes, reason, ierr)
          if (reason>0) then
             if (masterProc) then
                print *,"Converged!  SNESConvergedReason = ", reason
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
                print *,"Did not converge :(   SNESConvergedReason = ", reason
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

          call PetscTime(time1, ierr)
          if (solveSystem) then
             ! All the magic happens in this next line!
             call KSPSolve(KSPInstance,residualVec, solutionVec, ierr)
          end if

          call PetscTime(time2, ierr)
          if (masterProc) then
             print *,"Done with the main solve.  Time to solve: ", time2-time1, " seconds."
          end if
          call PetscTime(time1, ierr)

          if (useIterativeLinearSolver) then
             call KSPGetConvergedReason(KSPInstance, KSPReason, ierr)
             if (KSPReason>0) then
                if (masterProc) then
                   print *,"Converged!  KSPConvergedReason = ", KSPReason
                end if
                didNonlinearCalculationConverge = integerToRepresentTrue
             else
                if (masterProc) then
                   print *,"Did not converge :(   KSPConvergedReason = ", KSPReason
                end if
                didNonlinearCalculationConverge = integerToRepresentFalse
             end if
          else
             didNonlinearCalculationConverge = integerToRepresentTrue
          end if

          ! Compute flows, fluxes, etc.:
          call diagnostics(solutionVec, whichRHS)

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
!!$       call MatDestroy(preconditionerMatrix, ierr)
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

