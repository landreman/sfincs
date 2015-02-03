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

    external evaluateJacobian, evaluateResidual, diagnosticsMonitor

    if (masterProc) then
       print *,"Entering main solver loop."
    end if

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, solutionVec, ierr)
    call VecDuplicate(solutionVec, residualVec, ierr)

    call VecDuplicate(solutionVec, f0, ierr)
    call init_f0()
    !call VecDuplicate(solutionVec, temperatureEquilibrationTerm, ierr)
    !call initTemperatureEquilibrationTerm()

    call SNESCreate(MPIComm, mysnes, ierr)
    call SNESSetFunction(mysnes, residualVec, evaluateResidual, PETSC_NULL_OBJECT, ierr)

    call preallocateMatrix(matrix, 1)
    if (useIterativeSolver) then
       call preallocateMatrix(preconditionerMatrix, 0)
       call SNESSetJacobian(mysnes, matrix, preconditionerMatrix, evaluateJacobian, PETSC_NULL_OBJECT, ierr)
    else
       call SNESSetJacobian(mysnes, matrix, matrix, evaluateJacobian, PETSC_NULL_OBJECT, ierr)
    end if

    call SNESGetKSP(mysnes, KSPInstance, ierr)

    if (useIterativeSolver) then
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
       call KSPSetFromOptions(KSPInstance, ierr)
    end if

    if (numProcs > 1) then
       select case (whichParallelSolverToFactorPreconditioner)
       case (1)
          call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERMUMPS, ierr)
          if (masterProc) then
             print *,"We will use mumps to factorize the preconditioner."
          end if
       case (2)
          call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERSUPERLU_DIST, ierr)
          if (masterProc) then
             print *,"We will use superlu_dist to factorize the preconditioner."
          end if
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
       ! I originally always got an error message that there was a zero pivot. The following line seems to solve
       ! this problem.  The "zero pivot" error seemed to only arise with the "nd" ordering (nested dissection), which is the 
       ! default.  I'm not sure which of the other orderings is most efficient. I picked "rcm" for no particular reason.
       call PCFactorSetMatOrderingType(preconditionerContext, MATORDERINGRCM, ierr)

       ! I'm not sure this next line actually accomplishes anything, since it doesn't solve the "zero pivot" problem.
       ! But it doesn't seem to cost anything, and perhaps it reduces the chance of getting the "zero pivot" error in the future.
       call PCFactorReorderForNonzeroDiagonal(preconditionerContext, 1d-12, ierr) 

    end if

    ! Tell PETSc to call the diagnostics subroutine at each iteration of SNES:
    call SNESMonitorSet(mysnes, diagnosticsMonitor, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
    !call SNESMonitorSet(mysnes, SNESMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)

    ! Set the algorithm to use for the nonlinear solver:
    if (nonlinear) then
       ! SNESNEWTONLS = Newton's method with an optional line search.
       ! As of PETSc version 3.5, this is the default algorithm, but I'll set it manually here anyway to be safe.
       call SNESSetType(mysnes, SNESNEWTONLS, ierr)
       if (masterProc) then
          print *,"Since this is a nonlinear run, we will use Newton's method."
       end if
    else
       ! SNESKSPONLY = Only do 1 linear step.
       call SNESSetType(mysnes, SNESKSPONLY, ierr)
       if (masterProc) then
          print *,"Since this is a linear run, we will only take a single step, and not iterate Newton's method."
       end if
    end if

    call SNESSetFromOptions(mysnes, ierr)



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

       if (useIterativeSolver) then
          call SNESGetConvergedReason(mysnes, reason, ierr)
          if (reason>0) then
             if (masterProc) then
                print *,"Converged!  SNESConvergedReason = ", reason
             end if
             didNonlinearCalculationConverge = integerToRepresentTrue
          else
             if (masterProc) then
                print *,"Did not converge :(   SNESConvergedReason = ", reason
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
       call populateMatrix(matrix,1)


       if (useIterativeSolver) then

          ! Build the preconditioner:
          call populateMatrix(preconditionerMatrix,0)

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

          if (useIterativeSolver) then
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
    call MatDestroy(matrix, ierr)
    if (useIterativeSolver) then
       call MatDestroy(preconditionerMatrix, ierr)
    end if
    call SNESDestroy(mysnes,ierr)


  end subroutine mainSolverLoop


  ! ------------------------------------------------------------------------


  subroutine chooseParallelDirectSolver()

    implicit none

    logical :: isAParallelDirectSolverInstalled

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

  end subroutine chooseParallelDirectSolver



end module solver

