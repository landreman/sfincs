#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

subroutine testingNewtonSolve

  use globalVariables
  use solver
  use radialCoordinates

  PetscScalar :: deltaEr, radialCurrent_init, radialCurrent, dradialCurrentdEr_finitediff, dradialCurrentdEr_analytic

  deltaEr = 1.d-1

  ! Change settings so dRadialCurrentdEr is computed
  ambipolarSolve = .true.
  ambipolarSolveOption = 1
  call mainSolverLoop()
  radialCurrent_init = sum(particleFlux_vm_rN(1:Nspecies)*Zs(1:Nspecies))
  dRadialCurrentdEr_analytic = dRadialCurrentdEr
  if (masterProc) then
    print *,"radialCurrent_init: ", radialCurrent_init
    print *,"dPhiHatdpsiHat: ", dPhiHatdpsiHat
    print *,"particleFlux_vm_psiHat: ", particleFlux_vm_psiHat
    print *,"Zs: ", Zs
  end if

  ! Reset settings
  ambipolarSolve = .false.

  ! Update Er
  Er = Er + deltaEr
  !dPhiHatdpsiHat = dPhiHatdpsiHat + deltaEr
  !print *,"dPhiHatdpsiHat: ", dPhiHatdpsiHat
  call setInputRadialCoordinate()

  call mainSolverLoop()
  radialCurrent = sum(particleFlux_vm_rN(1:Nspecies)*Zs(1:Nspecies))

  if (masterProc) then
    print *,"radialCurrent: ", radialCurrent
  end if
  dradialCurrentdEr_finitediff = (radialCurrent-radialCurrent_init)/deltaEr

  if (masterProc) then
    print *,"dradialCurrentdEr_finitediff: ", dRadialCurrentdEr_finitediff
    print *,"dradialCurrentdEr_analytic: ", dRadialCurrentdEr_analytic
  end if
  stop

end subroutine testingNewtonSolve
