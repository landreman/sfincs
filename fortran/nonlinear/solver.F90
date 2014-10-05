module solver

  use globalVariables
  use petscsysdef

  implicit none

#include <finclude/petscsysdef.h>


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
          print *,"mumps no superlu_dist installed."
       end if
       stop
    end if

  end subroutine chooseParallelDirectSolver



end module scan

