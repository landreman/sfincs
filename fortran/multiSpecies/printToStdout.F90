module printToStdout

  use globalVariables
  use petscsysdef

  implicit none

#include <finclude/petscsysdef.h>

contains

  ! -----------------------------------------------------------------------------------

  subroutine printGreeting()

    implicit none

    if (masterProc) then
       print *,"****************************************************************************"
       print *,"SFINCS: Stellarator Fokker-Plank Iterative Neoclassical Conservative Solver"
       print *,"Grid in theta and zeta, Legendre modal in xi, polynomial spectral collocation in x."
       print *,"Multi-species version."
#if defined(PETSC_USE_REAL_SINGLE)
       print *,"Using single precision."
#else
       print *,"Using double precision."
#endif
       if (numProcs==1) then
          print *,"Serial job (1 process) detected."
       else
          print "(a, i4, a)", " Parallel job (",numProcs," processes) detected."
       end if
    end if
  end subroutine printGreeting

  ! -----------------------------------------------------------------------------------

  subroutine printInputs()

    implicit none

    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] ---- Physics parameters: ----"
       print *,"[",myCommunicatorIndex,"] Number of species = ", Nspecies
       print *,"[",myCommunicatorIndex,"] Delta     = ", Delta
       print *,"[",myCommunicatorIndex,"] alpha     = ", alpha
       print *,"[",myCommunicatorIndex,"] nu_n      = ", nu_n
       print *,"[",myCommunicatorIndex,"] dPhiHatdpsiN = ", dPhiHatdpsiN
       print *,"[",myCommunicatorIndex,"] ---- Geometry parameters: ----"
       print *,"[",myCommunicatorIndex,"] Geometry scheme = ", geometryScheme
       print *,"[",myCommunicatorIndex,"] iota      = ", iota
       print *,"[",myCommunicatorIndex,"] GHat      = ", GHat
       print *,"[",myCommunicatorIndex,"] IHat      = ", IHat
       print *,"[",myCommunicatorIndex,"] psiAHat   = ", psiAHat
       if (geometryScheme==1) then
          print *,"[",myCommunicatorIndex,"] epsilon_t = ", epsilon_t
          print *,"[",myCommunicatorIndex,"] epsilon_h = ", epsilon_h
          print *,"[",myCommunicatorIndex,"] epsilon_antisymm = ", epsilon_antisymm
       end if
       print *,"[",myCommunicatorIndex,"] ---- Numerical parameters: ----"
       print *,"[",myCommunicatorIndex,"] Nzeta              = ", Nzeta
       print *,"[",myCommunicatorIndex,"] Ntheta             = ", Ntheta
       print *,"[",myCommunicatorIndex,"] Nxi                = ", Nxi
       print *,"[",myCommunicatorIndex,"] NL                 = ", NL
       print *,"[",myCommunicatorIndex,"] Nx                 = ", Nx
       print *,"[",myCommunicatorIndex,"] NxPotentialsPerVth = ", NxPotentialsPerVth
       print *,"[",myCommunicatorIndex,"] xMax               = ",xMax
       print *,"[",myCommunicatorIndex,"] solverTolerance    = ",solverTolerance
       select case (thetaDerivativeScheme)
       case (0)
          print *,"[",myCommunicatorIndex,"] Theta derivative: spectral collocation"
       case (1)
          print *,"[",myCommunicatorIndex,"] Theta derivative: centered finite differences, 3-point stencil"
       case (2)
          print *,"[",myCommunicatorIndex,"] Theta derivative: centered finite differences, 5-point stencil"
       case default
          print *,"[",myCommunicatorIndex,"] Error! Invalid setting for thetaDerivativeScheme"
          stop
       end select
       if (useIterativeSolver) then
          print *,"[",myCommunicatorIndex,"] Using iterative solver"
       else
          print *,"[",myCommunicatorIndex,"] Using direct solver"
       end if
    end if
  end subroutine printInputs

  ! -----------------------------------------------------------------------------------

  subroutine printOutputs()

    implicit none

    integer :: i, ispecies

    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] Total elapsed time: ", elapsedTime, " seconds."
       do ispecies=1,Nspecies
          if (Nspecies>1) then
             print *,"[",myCommunicatorIndex,"] Results for species ",ispecies,":"
          end if
          print *,"[",myCommunicatorIndex,"]   FSADensityPerturbation:  ", FSADensityPerturbation(ispecies)
          print *,"[",myCommunicatorIndex,"]   FSABFlow:                ", FSABFlow(ispecies)
          print *,"[",myCommunicatorIndex,"]   FSAPressurePerturbation: ", FSAPressurePerturbation(ispecies)
          print *,"[",myCommunicatorIndex,"]   NTV:                     ", NTV(ispecies)
          print *,"[",myCommunicatorIndex,"]   particleFlux:            ", particleflux(ispecies)
          print *,"[",myCommunicatorIndex,"]   momentumFlux:            ", momentumflux(ispecies)
          print *,"[",myCommunicatorIndex,"]   heatFlux:                ", heatflux(ispecies)
       end do
       print *,"[",myCommunicatorIndex,"] FSABjHat (bootstrap current): ", FSABjHat
       if (rhsMode == 2) then
          print *,"[",myCommunicatorIndex,"] Transport matrix:"
          do i=1,3
             print *,"[",myCommunicatorIndex,"]   ", transportMatrix(i,:)
          end do
       end if

    end if

  end subroutine printOutputs

end module printToStdout

