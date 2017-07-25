
#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

module sfincs_main

  implicit none

  contains

  subroutine init_sfincs

    use globalVariables
    use readInput
    use petscsysdef

    implicit none

    PetscErrorCode ierr
    PetscLogDouble :: startTime, time1

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
    call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)
    masterProc = (myRank==0)

    ! In the future, if we want to divide the processors into sub-communicators, this next line would change:
    MPIComm = PETSC_COMM_WORLD

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

    call PetscTime(time1, ierr)
    startTime = time1

    if(.not.in_GS2) call readNamelistInput()

  end subroutine init_sfincs

  ! -----------------------------------------------------------------------------------

  subroutine prepare_sfincs

    use globalVariables
    use writeHDF5Output
    use solver
    use geometry
    use radialCoordinates
    use export_f

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
       if (includePhi1) then
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
       else
          print *,"Linear run"
       end if

       if (withAdiabatic) then
          print *,"Run with adiabatic species"
       end if

    end if

    ! Initialize NPeriods, psiAHat, and aHat.  We need to know NPeriods before
    ! we can initialize the zeta grid.
    call initializeGeometry()

    ! Do various calculations that will not need to be repeated at each
    ! iteration, such as setting up the coordinate grids and evaluating
    ! the magnetic field and its derivatives on the spatial grid.
    call createGrids()

  end subroutine prepare_sfincs

  ! -----------------------------------------------------------------------------------

  subroutine run_sfincs

    use globalVariables
    use geometry
    use writeHDF5Output
    use solver

    implicit none

    PetscErrorCode ierr

    ! Compute a few quantities related to the magnetic field:
    call computeBIntegrals()

    if (RHSMode==3) then
       ! Monoenergetic coefficient computation.
       ! Overwrite nu_n and dPhiHatd* using nuPrime and EStar.

       nu_n = nuPrime * B0OverBBar / (GHat + iota * IHat)
       dPhiHatdpsiHat = 2 / (alpha * Delta) * EStar * iota * B0OverBBar / GHat
    end if

    ! For input quantities that depend on the radial coordinate, pick out the values for the selected
    ! radial coordinate, and use these values to over-write values for the other radial coordinates.
    call setInputRadialCoordinate()

    ! Create HDF5 data structures, and save the quantities that will not change
    ! at each iteration of the solver (i.e. save all quantities except diagnostics.)
!    call initializeOutputFile()

    ! Solve the main system, either linear or nonlinear.
    ! This step takes more time than everything else combined.
    call mainSolverLoop()

    if(.not.in_GS2) call finalizeHDF5()
    call PetscFinalize(ierr)

    if (masterProc) then
       print *,"Goodbye!"
    end if

  end subroutine run_sfincs

  subroutine finalize_sfincs()
     use globalVariables
     use xGrid
     use indices
     use export_f

     implicit none

  !globalVariables array
  if(allocated(Nxi_for_x)) deallocate( Nxi_for_x, min_x_for_L)
  if(allocated(theta)) deallocate( theta, zeta, x, x_plus1)
  if(allocated(thetaWeights)) deallocate( thetaWeights, zetaWeights)
  if(allocated(ddtheta)) deallocate( ddtheta, ddzeta)
  if(allocated(ddtheta_ExB_plus)) deallocate( ddtheta_ExB_plus, ddtheta_ExB_minus)
  if(allocated(ddzeta_ExB_plus)) deallocate( ddzeta_ExB_plus, ddzeta_ExB_minus)
  if(allocated(ddtheta_magneticDrift_plus)) deallocate( ddtheta_magneticDrift_plus, ddtheta_magneticDrift_minus)
  if(allocated(ddzeta_magneticDrift_plus)) deallocate( ddzeta_magneticDrift_plus, ddzeta_magneticDrift_minus)
  if(allocated(xWeights)) deallocate( xWeights, xPotentials)
  if(allocated(x2)) deallocate( x2, expx2)
  if(allocated(ddx)) deallocate( ddx, d2dx2, ddxPotentials, d2dx2Potentials)
  if(allocated(ddx_xDot_plus)) deallocate( ddx_xDot_plus, ddx_xDot_preconditioner_plus)
  if(allocated(ddx_xDot_minus)) deallocate( ddx_xDot_minus, ddx_xDot_preconditioner_minus)
  if(allocated(ddx_preconditioner)) deallocate( ddx_preconditioner)
  if(allocated(ddtheta_preconditioner)) deallocate( ddtheta_preconditioner)
  if(allocated(ddzeta_preconditioner)) deallocate( ddzeta_preconditioner)
  if(allocated(interpolateXToXPotentials)) deallocate( interpolateXToXPotentials)
  if(allocated(RosenbluthPotentialTerms)) deallocate( RosenbluthPotentialTerms)
  if(allocated(BHat)) deallocate( BHat, dBHatdtheta, dBHatdzeta, dBHatdpsiHat, DHat)
  if(allocated(BHat_sub_psi)) deallocate( BHat_sub_psi, dBHat_sub_psi_dtheta, dBHat_sub_psi_dzeta)
  if(allocated(BHat_sub_theta)) deallocate( BHat_sub_theta, dBHat_sub_theta_dzeta, dBHat_sub_theta_dpsiHat)
  if(allocated(BHat_sub_zeta)) deallocate( BHat_sub_zeta, dBHat_sub_zeta_dtheta, dBHat_sub_zeta_dpsiHat)
  if(allocated(BHat_sup_theta)) deallocate( BHat_sup_theta, dBHat_sup_theta_dzeta, dBHat_sup_theta_dpsiHat)
  if(allocated(BHat_sup_zeta)) deallocate( BHat_sup_zeta, dBHat_sup_zeta_dtheta, dBHat_sup_zeta_dpsiHat)
  if(allocated(BDotCurlB)) deallocate( BDotCurlB, uHat, gradpsidotgradB_overgpsipsi)
  if(allocated(sources)) deallocate( sources, jHat, Phi1Hat, dPhi1Hatdtheta, dPhi1Hatdzeta)
  if(allocated(densityPerturbation)) deallocate( densityPerturbation, totalDensity)
  if(allocated(pressurePerturbation)) deallocate( pressurePerturbation, totalPressure, pressureAnisotropy)
  if(allocated(flow)) deallocate( flow, velocityUsingFSADensity, velocityUsingTotalDensity)
  if(allocated(MachUsingFSAThermalSpeed)) deallocate( MachUsingFSAThermalSpeed)
  if(allocated(particleFluxBeforeSurfaceIntegral_vm0)) deallocate( particleFluxBeforeSurfaceIntegral_vm0)
  if(allocated(particleFluxBeforeSurfaceIntegral_vm)) deallocate( particleFluxBeforeSurfaceIntegral_vm)
  if(allocated(particleFluxBeforeSurfaceIntegral_vE0)) deallocate( particleFluxBeforeSurfaceIntegral_vE0)
  if(allocated(particleFluxBeforeSurfaceIntegral_vE)) deallocate( particleFluxBeforeSurfaceIntegral_vE)
  if(allocated(momentumFluxBeforeSurfaceIntegral_vm0)) deallocate( momentumFluxBeforeSurfaceIntegral_vm0)
  if(allocated(momentumFluxBeforeSurfaceIntegral_vm)) deallocate( momentumFluxBeforeSurfaceIntegral_vm)
  if(allocated(momentumFluxBeforeSurfaceIntegral_vE0)) deallocate( momentumFluxBeforeSurfaceIntegral_vE0)
  if(allocated(momentumFluxBeforeSurfaceIntegral_vE)) deallocate( momentumFluxBeforeSurfaceIntegral_vE)
  if(allocated(heatFluxBeforeSurfaceIntegral_vm0)) deallocate( heatFluxBeforeSurfaceIntegral_vm0)
  if(allocated(heatFluxBeforeSurfaceIntegral_vm)) deallocate( heatFluxBeforeSurfaceIntegral_vm)
  if(allocated(heatFluxBeforeSurfaceIntegral_vE0)) deallocate( heatFluxBeforeSurfaceIntegral_vE0)
  if(allocated(heatFluxBeforeSurfaceIntegral_vE)) deallocate( heatFluxBeforeSurfaceIntegral_vE)
  if(allocated(NTVBeforeSurfaceIntegral)) deallocate( NTVBeforeSurfaceIntegral)
  if(allocated(NTVKernel)) deallocate( NTVKernel)
  if(allocated(FSADensityPerturbation)) deallocate( FSADensityPerturbation, FSAPressurePerturbation)
  if(allocated(FSABFlow)) deallocate( FSABFlow, FSABVelocityUsingFSADensity)
  if(allocated(FSABVelocityUsingFSADensityOverB0)) deallocate( FSABVelocityUsingFSADensityOverB0)
  if(allocated(FSABVelocityUsingFSADensityOverRootFSAB2)) deallocate( FSABVelocityUsingFSADensityOverRootFSAB2)
  if(allocated(particleFlux_vm0_psiHat)) deallocate( particleFlux_vm0_psiHat)
  if(allocated(particleFlux_vm0_psiN)) deallocate( particleFlux_vm0_psiN)
  if(allocated(particleFlux_vm0_rHat)) deallocate( particleFlux_vm0_rHat)
  if(allocated(particleFlux_vm0_rN)) deallocate( particleFlux_vm0_rN)
  if(allocated(particleFlux_vm_psiHat)) deallocate( particleFlux_vm_psiHat)
  if(allocated(particleFlux_vm_psiN)) deallocate( particleFlux_vm_psiN)
  if(allocated(particleFlux_vm_rHat)) deallocate( particleFlux_vm_rHat)
  if(allocated(particleFlux_vm_rN)) deallocate( particleFlux_vm_rN)
  if(allocated(particleFlux_vE0_psiHat)) deallocate( particleFlux_vE0_psiHat)
  if(allocated(particleFlux_vE0_psiN)) deallocate( particleFlux_vE0_psiN)
  if(allocated(particleFlux_vE0_rHat)) deallocate( particleFlux_vE0_rHat)
  if(allocated(particleFlux_vE0_rN)) deallocate( particleFlux_vE0_rN)
  if(allocated(particleFlux_vE_psiHat)) deallocate( particleFlux_vE_psiHat)
  if(allocated(particleFlux_vE_psiN)) deallocate( particleFlux_vE_psiN)
  if(allocated(particleFlux_vE_rHat)) deallocate( particleFlux_vE_rHat)
  if(allocated(particleFlux_vE_rN)) deallocate( particleFlux_vE_rN)
  if(allocated(particleFlux_vd1_psiHat)) deallocate( particleFlux_vd1_psiHat)
  if(allocated(particleFlux_vd1_psiN)) deallocate( particleFlux_vd1_psiN)
  if(allocated(particleFlux_vd1_rHat)) deallocate( particleFlux_vd1_rHat)
  if(allocated(particleFlux_vd1_rN)) deallocate( particleFlux_vd1_rN)
  if(allocated(particleFlux_vd_psiHat)) deallocate( particleFlux_vd_psiHat)
  if(allocated(particleFlux_vd_psiN)) deallocate( particleFlux_vd_psiN)
  if(allocated(particleFlux_vd_rHat)) deallocate( particleFlux_vd_rHat)
  if(allocated(particleFlux_vd_rN)) deallocate( particleFlux_vd_rN)
  if(allocated(momentumFlux_vm0_psiHat)) deallocate( momentumFlux_vm0_psiHat)
  if(allocated(momentumFlux_vm0_psiN)) deallocate( momentumFlux_vm0_psiN)
  if(allocated(momentumFlux_vm0_rHat)) deallocate( momentumFlux_vm0_rHat)
  if(allocated(momentumFlux_vm0_rN)) deallocate( momentumFlux_vm0_rN)
  if(allocated(momentumFlux_vm_psiHat)) deallocate( momentumFlux_vm_psiHat)
  if(allocated(momentumFlux_vm_psiN)) deallocate( momentumFlux_vm_psiN)
  if(allocated(momentumFlux_vm_rHat)) deallocate( momentumFlux_vm_rHat)
  if(allocated(momentumFlux_vm_rN)) deallocate( momentumFlux_vm_rN)
  if(allocated(momentumFlux_vE0_psiHat)) deallocate( momentumFlux_vE0_psiHat)
  if(allocated(momentumFlux_vE0_psiN)) deallocate( momentumFlux_vE0_psiN)
  if(allocated(momentumFlux_vE0_rHat)) deallocate( momentumFlux_vE0_rHat)
  if(allocated(momentumFlux_vE0_rN)) deallocate( momentumFlux_vE0_rN)
  if(allocated(momentumFlux_vE_psiHat)) deallocate( momentumFlux_vE_psiHat)
  if(allocated(momentumFlux_vE_psiN)) deallocate( momentumFlux_vE_psiN)
  if(allocated(momentumFlux_vE_rHat)) deallocate( momentumFlux_vE_rHat)
  if(allocated(momentumFlux_vE_rN)) deallocate( momentumFlux_vE_rN)
  if(allocated(momentumFlux_vd1_psiHat)) deallocate( momentumFlux_vd1_psiHat)
  if(allocated(momentumFlux_vd1_psiN)) deallocate( momentumFlux_vd1_psiN)
  if(allocated(momentumFlux_vd1_rHat)) deallocate( momentumFlux_vd1_rHat)
  if(allocated(momentumFlux_vd1_rN)) deallocate( momentumFlux_vd1_rN)
  if(allocated(momentumFlux_vd_psiHat)) deallocate( momentumFlux_vd_psiHat)
  if(allocated(momentumFlux_vd_psiN)) deallocate( momentumFlux_vd_psiN)
  if(allocated(momentumFlux_vd_rHat)) deallocate( momentumFlux_vd_rHat)
  if(allocated(momentumFlux_vd_rN)) deallocate( momentumFlux_vd_rN)
  if(allocated(heatFlux_vm0_psiHat)) deallocate( heatFlux_vm0_psiHat)
  if(allocated(heatFlux_vm0_psiN)) deallocate( heatFlux_vm0_psiN)
  if(allocated(heatFlux_vm0_rHat)) deallocate( heatFlux_vm0_rHat)
  if(allocated(heatFlux_vm0_rN)) deallocate( heatFlux_vm0_rN)
  if(allocated(heatFlux_vm_psiHat)) deallocate( heatFlux_vm_psiHat)
  if(allocated(heatFlux_vm_psiN)) deallocate( heatFlux_vm_psiN)
  if(allocated(heatFlux_vm_rHat)) deallocate( heatFlux_vm_rHat)
  if(allocated(heatFlux_vm_rN)) deallocate( heatFlux_vm_rN)
  if(allocated(heatFlux_vE0_psiHat)) deallocate( heatFlux_vE0_psiHat)
  if(allocated(heatFlux_vE0_rHat)) deallocate( heatFlux_vE0_rHat)
  if(allocated(heatFlux_vE0_rN)) deallocate( heatFlux_vE0_rN)
  if(allocated(heatFlux_vE0_psiN)) deallocate( heatFlux_vE0_psiN)
  if(allocated(heatFlux_vE_psiHat)) deallocate( heatFlux_vE_psiHat)
  if(allocated(heatFlux_vE_rHat)) deallocate( heatFlux_vE_rHat)
  if(allocated(heatFlux_vE_rN)) deallocate( heatFlux_vE_rN)
  if(allocated(heatFlux_vE_psiN)) deallocate( heatFlux_vE_psiN)
  if(allocated(heatFlux_vd1_psiHat)) deallocate( heatFlux_vd1_psiHat)
  if(allocated(heatFlux_vd1_psiN)) deallocate( heatFlux_vd1_psiN)
  if(allocated(heatFlux_vd1_rHat)) deallocate( heatFlux_vd1_rHat)
  if(allocated(heatFlux_vd1_rN)) deallocate( heatFlux_vd1_rN)
  if(allocated(heatFlux_vd_psiHat)) deallocate( heatFlux_vd_psiHat)
  if(allocated(heatFlux_vd_psiN)) deallocate( heatFlux_vd_psiN)
  if(allocated(heatFlux_vd_rHat)) deallocate( heatFlux_vd_rHat)
  if(allocated(heatFlux_vd_rN)) deallocate( heatFlux_vd_rN)
  if(allocated(heatFlux_withoutPhi1_psiHat)) deallocate( heatFlux_withoutPhi1_psiHat)
  if(allocated(heatFlux_withoutPhi1_psiN)) deallocate( heatFlux_withoutPhi1_psiN)
  if(allocated(heatFlux_withoutPhi1_rHat)) deallocate( heatFlux_withoutPhi1_rHat)
  if(allocated(heatFlux_withoutPhi1_rN)) deallocate( heatFlux_withoutPhi1_rN)
  if(allocated(particleFlux_vm_psiHat_vs_x)) deallocate( particleFlux_vm_psiHat_vs_x)
  if(allocated(heatFlux_vm_psiHat_vs_x)) deallocate( heatFlux_vm_psiHat_vs_x)
  if(allocated(FSABFlow_vs_x)) deallocate( FSABFlow_vs_x)
  if(allocated(NTV)) deallocate( NTV)
  if(allocated(transportMatrix)) deallocate( transportMatrix)

  if(allocated(a)) deallocate(a,b,c)

  if(allocated(first_index_for_x)) deallocate(first_index_for_x)

  if(allocated(map_theta_to_export_f_theta))deallocate(map_theta_to_export_f_theta)
  if(allocated(map_zeta_to_export_f_zeta))deallocate(map_zeta_to_export_f_zeta)
  if(allocated(map_x_to_export_f_x))deallocate(map_x_to_export_f_x)
  if(allocated(map_xi_to_export_f_xi))deallocate(map_xi_to_export_f_xi)
  if(allocated(delta_f))deallocate(delta_f)
  if(allocated(full_f))deallocate(full_f)

end subroutine finalize_sfincs



end module sfincs_main
