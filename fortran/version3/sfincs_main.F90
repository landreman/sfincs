module sfincs_main

#include "PETScVersions.F90"

  implicit none

  contains

  subroutine sfincs_init(MPI_comm_to_use)
    
    use globalVariables
    use readInput
    use petscsysdef
    
    implicit none
    
    !MPI_COMM :: MPI_comm_to_use
    integer :: MPI_comm_to_use
    PetscErrorCode ierr
    double precision :: startTime, time1

    PETSC_COMM_WORLD = MPI_comm_to_use
    MPIComm = PETSC_COMM_WORLD
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
    call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)
    masterProc = (myRank==0)
    
    
    if (masterProc) then
       print *,"****************************************************************************"
       print *,"SFINCS: Stellarator Fokker-Plank Iterative Neoclassical Conservative Solver"
       print *,"Version 3"
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
    
    time1 = MPI_Wtime()
    startTime = time1
    
    call readNamelistInput()

  end subroutine sfincs_init

  ! -----------------------------------------------------------------------------------

  subroutine sfincs_prepare

    use globalVariables
    use writeHDF5Output
    use solver
    use geometry
    use radialCoordinates
    use readHDF5Input !!Added by AM 2018-12
    
    implicit none
    
    PetscErrorCode ierr

    call validateInput()
    
    ! If running with >1 proc,
    ! make sure either superlu_dist or mumps is installed, and pick which one
    ! of these packages to use:
    call chooseParallelDirectSolver()
    
    if (masterProc) then
       print *,"---- Physics parameters: ----"
       print *,"Number of particle species = ", Nspecies
       print *,"Delta (rho* at reference parameters)          = ", Delta
       print *,"alpha (e Phi / T at reference parameters)     = ", alpha
       print *,"nu_n (collisionality at reference parameters) = ", nu_n
       !!if (includePhi1) then !!Commented by AM 2018-12
       if (includePhi1 .and. (.not. readExternalPhi1)) then !!Added by AM 2018-12
          print *,"Nonlinear run"
          if (includePhi1InKineticEquation) then
             print *,"with Phi1 included in the kinetic equation"
          else
             print *,"but with Phi1 excluded from the kinetic equation"
          end if
          
          if (quasineutralityOption == 1) then
             print *,"Using full quasi-neutrality equation"
          else
             print *,"Using EUTERPE quasi-neutrality equation"
          end if
       else if (includePhi1  .and. includePhi1InKineticEquation .and. readExternalPhi1) then !!Added by AM 2018-12
          print *,"Linear run but with Phi1 read from external file and included in the kinetic equation" !!Added by AM 2018-12
       else
          print *,"Linear run"
       end if
       
       if (withAdiabatic) then
          print *,"Run with adiabatic species"
       end if
       if (withNBIspec) then
          print *,"Run with NBI species"
       end if
       
    end if
    
    ! Initialize NPeriods, psiAHat, and aHat.  We need to know NPeriods before
    ! we can initialize the zeta grid.
    call initializeGeometry()
    
    ! Do various calculations that will not need to be repeated at each
    ! iteration, such as setting up the coordinate grids and evaluating
    ! the magnetic field and its derivatives on the spatial grid.
    call createGrids()
    
    if (includePhi1 .and. readExternalPhi1) then !!Added by AM 2018-12
       call setPhi1() !!Added by AM 2018-12
    end if !!Added by AM 2018-12

  end subroutine sfincs_prepare

  ! -----------------------------------------------------------------------------------

  subroutine sfincs_run

    use globalVariables
    use geometry
    use writeHDF5Output
    use solver
    use classicalTransport
    use ambipolarSolver
    use testingAdjointDiagnostics
    
    implicit none
    
    PetscErrorCode ierr

    ! Compute a few quantities related to the magnetic field:
! This is already called from create_grids()
    !call computeBIntegrals()

    if (RHSMode==3) then
       ! Monoenergetic coefficient computation.
       ! Overwrite nu_n and dPhiHatd* using nuPrime and EStar.
       
       nu_n = nuPrime * B0OverBBar / (GHat + iota * IHat)
       dPhiHatdpsiHat = 2 / (alpha * Delta) * EStar * iota * B0OverBBar / GHat
    end if

    ! For input quantities that depend on the radial coordinate, pick out the values for the selected
    ! radial coordinate, and use these values to over-write values for the other radial coordinates.
    call setInputRadialCoordinate()

    ! calculate the classical transport without Phi1
    ! this can be done without solving the system, and so is done here
    call calculateClassicalFlux(.false.,classicalParticleFluxNoPhi1_psiHat,classicalHeatFluxNoPhi1_psiHat)
    classicalParticleFluxNoPhi1_psiN = ddpsiN2ddpsiHat * classicalParticleFluxNoPhi1_psiHat
    classicalParticleFluxNoPhi1_rHat = ddrHat2ddpsiHat * classicalParticleFluxNoPhi1_psiHat
    classicalParticleFluxNoPhi1_rN = ddrN2ddpsiHat * classicalParticleFluxNoPhi1_psiHat
    classicalHeatFluxNoPhi1_psiN = ddpsiN2ddpsiHat * classicalHeatFluxNoPhi1_psiHat
    classicalHeatFluxNoPhi1_rHat = ddrHat2ddpsiHat * classicalHeatFluxNoPhi1_psiHat
    classicalHeatFluxNoPhi1_rN = ddrN2ddpsiHat * classicalHeatFluxNoPhi1_psiHat

    ! Create HDF5 data structures, and save the quantities that will not change
    ! at each iteration of the solver (i.e. save all quantities except diagnostics.)
    call initializeOutputFile()

    if (debugAdjoint) then
      call compareAdjointDiagnostics()
    end if

    ! Solve the main system, either linear or nonlinear.
    ! This step takes more time than everything else combined.
	if (debugAdjoint .eqv. .false.) then
	    if (ambipolarSolve) then
		    call mainAmbipolarSolver()
		else
			call mainSolverLoop()
		end if
	end if
    
    call finalizeHDF5()
    call PetscFinalize(ierr)
    
    if (masterProc) then
       print *,"Goodbye!"
    end if
    

  end subroutine sfincs_run

end module sfincs_main
