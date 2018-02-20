#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

!! This subroutine
subroutine ambipolarSolver

  use globalVariables
  use radialCoordinates
  use solver

  implicit none

  integer :: stage, iEr, next_stage
  logical :: last_above_target, initial_above_target, positive_above_target, negative_above_target
  PetscScalar, dimension(:), allocatable :: Er_search, radialCurrent
  PetscScalar :: Brendt_e, Brendt_tol1, Brendt_fa, Brendt_fb, Brendt_s, Brendt_a, Brendt_c, Brendt_p, Brendt_q, Brendt_xm, Brendt_r, Brendt_b, Brendt_d, Brendt_eps, Brendt_fc, Er_init
  integer :: exit_code

  allocate(Er_search(nEr))
  allocate(radialCurrent(nEr))
  Er_search = zero
  Er_init = Er
  radialCurrent = zero

  stage = 4
  exit_code = -1

  do iEr = 1,nEr
    select case(stage)
    case(1)
      ! Initial guess for Er
      Er_search(iEr) = Er_init
      next_stage = 2
    case(2)
      ! Vary Er by 10 until we bracket radialCurrent = 0
      next_stage = 2

      ! radialCurrent(Er_init)>0 and radialCurrent(+100)>0
      ! or radialCurrent(Er_init)<0 and radialCurrent(+100)<0
      if (initial_above_target .eqv. positive_above_target) then
        Er_search(iEr) = Er_search(iEr-1) - 10
      else
        Er_search(iEr) = Er_search(iEr-1) + 10
      end if
    case(3)
      next_stage = 3

      if (abs(Brendt_e) >= Brendt_tol1 .and. abs(Brendt_fa)>abs(Brendt_fb)) then
           ! Attempt inverse quadratic interpolation
           Brendt_s = Brendt_fb / Brendt_fa
           if (Brendt_a == Brendt_c) then
              Brendt_p = 2.0*Brendt_xm*Brendt_s
              Brendt_q = 1.0-Brendt_s
           else
              Brendt_q = Brendt_fa / Brendt_fc
              Brendt_r = Brendt_fb / Brendt_fc
              Brendt_p = Brendt_s*(2.0*Brendt_xm*Brendt_q*(Brendt_q-Brendt_r)-(Brendt_b-Brendt_a)*(Brendt_r-1.0))
              Brendt_q = (Brendt_q-1.0)*(Brendt_r-1.0)*(Brendt_s-1.0)
           end if
           if (Brendt_p > 0) Brendt_q = -Brendt_q ! Check whether in bounds
           Brendt_p = abs(Brendt_p)
           if (2.0*Brendt_p < min(3.0*Brendt_xm*Brendt_q-abs(Brendt_tol1*Brendt_q),abs(Brendt_e*Brendt_q))) then
             ! Accept interpolation
             Brendt_e = Brendt_d
             Brendt_d = Brendt_p / Brendt_q
           else
             ! Interpolation failed, so use bisection
             Brendt_d = Brendt_xm
             Brendt_e = Brendt_d
           end if
        else
           ! Bounds are decreasing too slowly, so use bisection
           Brendt_d = Brendt_xm
           Brendt_e = Brendt_d
        end if
        ! Move last best guess to a.
        Brendt_a = Brendt_b
        Brendt_fa = Brendt_fb
        if (abs(Brendt_d) > Brendt_tol1) then
           ! Evaluate new trial root
           Brendt_b = Brendt_b + Brendt_d
        else
           Brendt_b = Brendt_b + sign(Brendt_tol1,Brendt_xm)
        end if

        Er_search(iEr) = Brendt_b
      case(4) ! Er = 100
        Er_search(iEr) = 100
        next_stage = 5
      case(5) ! Er = -100
        Er_search(iEr) = -100
        next_stage = 1
    end select

    if (masterProc) then
      print "(a,es10.3,a,i3,a,i3,a)"," Solving system for Er=",Er_search(iEr)," (",iEr," of at most ",nEr,")"
    end if

    Er = Er_search(iEr)
    ! Update dPhiHatdpsiHat, dPhiHatdpsiN, dPhiHatdrHat, dPhiHatdrN
    dPhiHatdpsiHat = ddrHat2ddpsiHat * (-Er)
    dPhiHatdpsiN = ddpsiHat2ddpsiN * dPhiHatdpsiHat
    dPhiHatdrHat = ddpsiHat2ddrHat * dPhiHatdpsiHat
    dPhiHatdrN   = ddpsiHat2ddrN   * dPhiHatdpsiHat

!    call setInputRadialCoordinate()
    call mainSolverLoop()
    radialCurrent(iEr) = sum(particleFlux_vm_psiHat(1:Nspecies)*Zs(1:Nspecies))
    if (masterProc) then
      print "(a,es10.3)","   radial current:",radialCurrent(iEr)
    end if
    last_above_target = (radialCurrent(iEr) > zero)
    if (stage == 4) positive_above_target = last_above_target
    if (stage == 5) negative_above_target = last_above_target
    if (stage == 1) initial_above_target = last_above_target
    if (stage == 2 .and. (last_above_target .neqv. initial_above_target)) then
        ! If we've bracketed the target, move on to stage 3.
        next_stage = 3
        if (masterProc) then
          print *,"Zero radial current has been bracketed."
        end if
    end if
    ! Initialize Brendt's algorithm for root-finding
    if (stage==2 .and. next_stage == 3) then
        Brendt_a = Er_search(iEr-1)
        Brendt_b = Er_search(iEr)
        Brendt_fa = radialCurrent(iEr-1)
        Brendt_fb = radialCurrent(iEr)
        Brendt_c = Brendt_b
        Brendt_fc = Brendt_fb
        Brendt_d = Brendt_b - Brendt_a
        Brendt_e = Brendt_d
    end if ! stage=2 and next_stage = 3
    if (stage==3) Brendt_fb = radialCurrent(iEr)
    if (next_stage==3) then
      if ((Brendt_fb > 0 .and. Brendt_fc > 0) .or. (Brendt_fb < 0 .and. Brendt_fc < 0)) then
           Brendt_c = Brendt_a
           Brendt_fc = Brendt_fa
           Brendt_d = Brendt_b - Brendt_a
           Brendt_e = Brendt_d
        end if
        if (abs(Brendt_fc) < abs(Brendt_fb)) then
           Brendt_a = Brendt_b
           Brendt_b = Brendt_c
           Brendt_c = Brendt_a
           Brendt_fa = Brendt_fb
           Brendt_fb = Brendt_fc
           Brendt_fc = Brendt_fa
        end if
        Brendt_EPS = 1.d-15
        Brendt_tol1 = 2.0*Brendt_EPS*abs(Brendt_b) + 0.5*Er_search_tolerance
        Brendt_xm = 0.5*(Brendt_c - Brendt_b)
        if (abs(Brendt_xm) <= Brendt_tol1 .or. (Brendt_fb==0)) then
           ! We met the requested tolerance
           if (masterProc) then
            print *,"Requested tolerance has been met."
           end if
           exit_code=0
           NEr = iEr
           exit
        end if
    end if ! next_stage=3

    stage = next_stage
  end do

  if (exit_code == -1) then
     if (masterProc) then
       print *,"*******************************************************************************"
       print *,"*******************************************************************************"
       print *,"The Er search did not converge within NEr iterations!"
       print *,"*******************************************************************************"
       print *,"*******************************************************************************"
     end if
      print *,"stage: ", stage
      print *,"Here are the Ers we used: "
      print *,Er_search
      print *,"Here are the radial currents: "
      print *,radialCurrent
     stop
  end if
  if (masterProc) then
    print *,"stage: ", stage
    print *,"Here are the Ers we used: "
    print *,Er_search
    print *,"Here are the radial currents: "
    print *,radialCurrent
  end if

  ! Output was not written previously
  if (debugAdjoint .eqv. .false.) then
    call updateOutputFile(1, .false.)
  end if

  ! Deallocate stuff
  deallocate(Er_search)
  deallocate(radialCurrent)

  ! If adjoint solve required, call solver again with Er solution
  ambipolarSolve = .false.
  call mainSolverLoop()

end subroutine ambipolarSolver
