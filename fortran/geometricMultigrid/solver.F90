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
    PetscScalar :: scalar
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
    KSP, dimension(:), allocatable :: sub_ksps
    Mat :: sub_Amat, sub_Pmat
    MatNullSpace :: nullspace
    Vec :: temp_Vec

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
       call VecZeroEntries(solutionVec, ierr)

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
       call VecZeroEntries(solutionVec, ierr)

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
          call VecZeroEntries(solutionVec, ierr)

          call evaluateResidual(mysnes, solutionVec, residualVec, userContext, ierr)
          ! Multiply the residual by (-1):
          scalar = -1
          call VecScale(residualVec, scalar, ierr)

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
       print *,"mumps detected."
    end if
#else
    whichParallelSolverToFactorPreconditioner = 2
    if (masterProc) then
       print *,"mumps not detected."
    end if
#endif

#ifdef PETSC_HAVE_SUPERLU_DIST
    isAParallelDirectSolverInstalled = .true.
    if (masterProc) then
       print *,"superlu_dist detected."
    end if
#else
    if (masterProc) then
       print *,"superlu_dist not detected."
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

