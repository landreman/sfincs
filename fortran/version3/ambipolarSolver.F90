#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

module ambipolarSolver

  contains

  !! This subroutine begins the search for ambipolarEr by bracketing radialCurrent = zero
  subroutine mainAmbipolarSolver

    use globalVariables
    use solver
    use writeHDF5Output

    implicit none

    integer :: stage, next_stage, iEr, exit_code
    PetscScalar :: Er_init, this_Er, this_radialCurrent, ambipolarEr
    logical :: initial_above_target, last_above_target, positive_above_target, negative_above_target
    PetscScalar, dimension(:), allocatable :: Er_search, radialCurrent
    PetscScalar :: time1, time2, dRadialCurrentdEr_init
    PetscErrorCode :: ierr

    allocate(Er_search(NEr_ambipolarSolve))
    allocate(radialCurrent(NEr_ambipolarSolve))
    Er_search = zero
    radialCurrent = zero

    Er_init = Er
    stage = 1
    exit_code = -1

    do iEr = 1,NEr_ambipolarSolve
      select case(stage)
        case(1) ! Er = 100
          Er_search(iEr) = 100
          next_stage = 2
        case(2) ! Er = -100
          Er_search(iEr) = -100
          next_stage = 3
        case(3) ! Initial guess
          Er_search(iEr) = Er_init
          if (ambipolarSolveOption==1) then
            dRadialCurrentdEr_init = dRadialCurrentdEr
          end if
          next_stage = 4

        case(4)
          ! Vary Er by 10 until we bracket radialCurrent = 0
          next_stage = 4

          if (initial_above_target .eqv. positive_above_target) then
            Er_search(iEr) = Er_search(iEr-1) - 10
          else
            Er_search(iEr) = Er_search(iEr-1) + 10
          end if

      end select

      this_Er = Er_search(iEr)
      call updateEr(this_Er,this_radialCurrent)
      radialCurrent(iEr) = this_radialCurrent
      if (masterProc) then
        print *,"stage = ", stage
        print *,"Er = ", this_Er
        print *,"radial current = ", this_radialCurrent
      end if

      last_above_target = (radialCurrent(iEr) > zero)
      if (stage == 1) positive_above_target = last_above_target
      if (stage == 2) then
        negative_above_target = last_above_target
        if (negative_above_target .eqv. positive_above_target) then
          if (masterProc) then
            print *,"Error in ambipolarSolver! Unexpected behavior at Er = +-100."
            print *,"Here are the Ers we used: "
            print *,Er_search
            print *,"Here are the radial currents: "
            print *,radialCurrent
          end if
          stop
        end if
      end if
      if (stage == 3) initial_above_target = last_above_target
      if (stage == 4 .and. (last_above_target .neqv. initial_above_target)) then
          ! We've bracketed zero radial current
          exit_code = 0
          exit
      end if
      stage = next_stage

  end do ! iEr

  if (exit_code == 0) then
    print *,"Successful bracketing in ambipolar solve."
    print *,"Here are the Ers we used: "
    print *,Er_search
    print *,"Here are the radial currents: "
    print *,radialCurrent

    call PetscTime(time1, ierr)
    if (ambipolarSolveOption==1) then
      call ambipolarSolverNewton(Er_search(iEr),Er_search(3),radialCurrent(iEr),radialCurrent(3),dRadialCurrentdEr_init,ambipolarEr)
    else if (ambipolarSolveOption==2) then
      call ambipolarSolverBrent(Er_search(iEr),Er_search(3),radialCurrent(iEr),radialCurrent(3),ambipolarEr)
    end if
    call PetscTime(time2,ierr)
    if (masterProc) then
      print *,"Time for ambipolar solve: ",time2-time1, " seconds."
    end if
  end if

  if (exit_code == -1) then
    if (masterProc) then
      print *,"Error! Ambipolar solver was not able to bracket zero radial current."
      print *,"Here are the Ers we used: "
      print *,Er_search
      print *,"Here are the radial currents: "
      print *,radialCurrent
    end if
    stop
  end if

    ! Output was not written previously
    if (debugAdjoint .eqv. .false.) then
      call updateOutputFile(1, .false.)
    end if

    ! Deallocate stuff
    deallocate(Er_search)
    deallocate(radialCurrent)

    ! If adjoint solve required, call solver again with Er solution
    if (RHSMode>3) then
      ambipolarSolve = .false.
      call mainSolverLoop()
      ambipolarSolve = .true.
    end if

  end subroutine mainAmbipolarSolver

  !> Ambipolar root is known to lie between x1 and x2 with radial current values
  !> f1 and f2
  subroutine ambipolarSolverBrent(Er1,Er2,f1,f2,ambipolarEr)

    use globalVariables

    implicit none

    PetscScalar :: Er1, Er2, f1, f2, ambipolarEr
    integer :: iEr, exit_code
    PetscScalar :: thisRadialCurrent, Brendt_EPS
    PetscScalar, dimension(:), allocatable :: Er_search, radialCurrent
    PetscScalar :: Brendt_a, Brendt_b, Brendt_fa, Brendt_fb, Brendt_c, Brendt_fc, Brendt_d, Brendt_e, Brendt_tol1, Brendt_p, Brendt_q, Brendt_r, Brendt_s, Brendt_xm
    PetscErrorCode :: ierr
    PetscScalar :: time1, time2


    allocate(Er_search(NEr_ambipolarSolve))
    allocate(radialCurrent(NEr_ambipolarSolve))
    Er_search = zero
    radialCurrent = zero
    exit_code = -1

    Brendt_EPS = 1.d-15

    ! Initialize Brent algorithm
    Brendt_a = Er1
    Brendt_b = Er2
    Brendt_fa = f1
    Brendt_fb = f2
    Brendt_c = Brendt_b

    if ((Brendt_fa > zero .and. Brendt_fb > zero) .or. (Brendt_fa < zero .and. Brendt_fb < zero)) then
      if (masterProc) then
        print *,"Root must be bracketed in Brent solve!"
      end if
      stop
    end if
    Brendt_fc = Brendt_fb

    do iEr = 1, NEr_ambipolarSolve
      call PetscTime(time1, ierr)
      if ((Brendt_fb > zero .and. Brendt_fc > zero) .or. (Brendt_fb < zero .and. Brendt_fc < zero)) then
        Brendt_c = Brendt_a
        Brendt_fc = Brendt_fa
        Brendt_e = Brendt_b - Brendt_a
        Brendt_d = Brendt_b - Brendt_a
      end if
      if (abs(Brendt_fc) < abs(Brendt_fb)) then
        Brendt_a = Brendt_b
        Brendt_b = Brendt_c
        Brendt_c = Brendt_a
        Brendt_fa = Brendt_fb
        Brendt_fb = Brendt_fc
        Brendt_fc = Brendt_fa
      end if
      ! Convergence check
      Brendt_tol1 = two*Brendt_EPS*abs(Brendt_b)+0.5*Er_search_tolerance
      Brendt_xm = 0.5*(Brendt_c-Brendt_b)
      if (abs(Brendt_xm) <= Brendt_tol1 .or. Brendt_fb == zero) then
        ambipolarEr = Brendt_b
        exit_code = 0
        exit
      end if

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
      call updateEr(Brendt_b, Brendt_fb)
      radialCurrent(iEr) = Brendt_fb
      Er_search(iEr) = Brendt_b
      call PetscTime(time2,ierr)
      if (masterProc) then
        print *,"One Brent loop: ", time2-time1, " sec."
      end if
    end do ! iEr

    if (exit_code == -1) then
      if (masterProc) then
         print *,"*******************************************************************************"
         print *,"*******************************************************************************"
         print *,"The Er search did not converge within NEr iterations!"
         print *,"*******************************************************************************"
         print *,"*******************************************************************************"
         print *,"Here are the Ers we used: "
         print *,Er_search
         print *,"Here are the radial currents: "
         print *,radialCurrent
      end if
      stop
  end if
  if (exit_code == 0) then
    if (masterProc) then
      print *,"Brent algorithm successful."
      print *,"Here are the Ers we used: "
      print *,Er_search
      print *,"Here are the radial currents: "
      print *,radialCurrent
    end if
  end if

  end subroutine ambipolarSolverBrent

  ! This subroutine should be called after ambipolar Er has been bracketed. This Newton method uses dRadialCurrentdEr computed
  ! using the adjoint method.
  subroutine ambipolarSolverNewton(Er1, Er2, fl, fh, df, ambipolarEr)

    use globalVariables
    use radialCoordinates
    use solver

    implicit none

    PetscScalar :: Er1, Er2, fl, fh, ambipolarEr, xl, xh
    PetscScalar :: rts, dxold, dx, f, df, temp
    integer :: j, exit_code
    PetscScalar, dimension(:), allocatable :: Er_search, radialCurrent
    PetscErrorCode :: ierr
    PetscScalar :: time1, time2

    allocate(Er_search(NEr_ambipolarSolve))
    allocate(radialCurrent(NEr_ambipolarSolve))
    Er_search = zero
    radialCurrent = zero
!    call updateEr(Er1,fl)
!    call updateEr(Er2,fh)

    if ((fl > zero  .and. fh > zero) .or. (fl < zero .and. fh < zero)) then
      if (masterProc) then
        print *,"Error! Root must be bracketed in newton solve."
      end if
      stop
    end if
    ! Orient the search so that f(x1) < 0
    if (fl < zero) then
      xl = Er1
      xh = Er2
    else
      xh = Er1
      xl = Er2
    end if

    ! Initialize guess for root, "stepsize before last", and last step
    ! Initial guess it taken to be R1, and df is input parameter
!    rts = 0.5*(Er1+Er2)
    dxold = abs(Er2-Er1)
    dx = dxold
   ! call updateEr(rts,f)
!    df = dRadialCurrentdEr
    rts = Er1
    f = fl

    exit_code = -1
    ! Loop over allowed iterations
    do j = 1, NEr_ambipolarSolve
      call PetscTime(time1,ierr)
      ! Bisect if Newton out of range or not decreasing fast enough
      if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) .or. (abs(two*f) > abs(dxold*df))) then
        dxold = dx
        dx = 0.5*(xh-xl)
        rts = xl+dx
        if (xl == rts) then
          ambipolarEr = rts
          exit_code = 0
          exit
        end if
      ! Newton step acceptable. Take it.
      else
        dxold = dx
        dx = f/df
        temp = rts
        rts = rts - dx
        if (temp == rts) then
          ambipolarEr = rts
          exit_code = 0
          exit
        end if
      end if
      ! Convergence criterion
      if (abs(dx) < Er_search_tolerance) then
        ambipolarEr = rts
        exit_code = 0
        exit
      end if
      ! The one new function evaluation per iteration
      call updateEr(rts,f)
      Er_search(j) = rts
      radialCurrent(j) = f
      df = dRadialCurrentdEr
      ! Maintain the bracket on the root
      if (f < zero) then
        xl = rts
      else
        xh = rts
      end if
      call PetscTime(time2,ierr)
      if (masterProc) then
        print *,"Time for one Newton loop: ", time2-time1," sec."
      end if
    end do ! iEr

    if (exit_code == 0) then
      if (masterProc) then
        print *,"Newton ambipolar solve was successful."
        print *,"Here are the Ers we used: "
        print *,Er_search
        print *,"Here are the radial currents: "
        print *,radialCurrent
      end if
    end if

    if (exit_code == -1) then
      if (masterProc) then
         print *,"*******************************************************************************"
         print *,"*******************************************************************************"
         print *,"The Er search did not converge within NEr iterations!"
         print *,"*******************************************************************************"
         print *,"*******************************************************************************"
         print *,"Here are the Ers we used: "
         print *,Er_search
         print *,"Here are the radial currents: "
         print *,radialCurrent
      end if
      stop
    end if

  end subroutine ambipolarSolverNewton

  subroutine updateEr(thisEr,radialCurrent)

    use globalVariables
    use solver

    implicit none

    PetscScalar :: thisEr, radialCurrent
    PetscErrorCode :: ierr

    Er = thisEr

    ! Update dPhiHatdpsiHat, dPhiHatdpsiN, dPhiHatdrHat, dPhiHatdrN
    dPhiHatdpsiHat = ddrHat2ddpsiHat * (-Er)
    dPhiHatdpsiN = ddpsiHat2ddpsiN * dPhiHatdpsiHat
    dPhiHatdrHat = ddpsiHat2ddrHat * dPhiHatdpsiHat
    dPhiHatdrN   = ddpsiHat2ddrN   * dPhiHatdpsiHat

    call mainSolverLoop()

    radialCurrent = sum(particleFlux_vm_rN(1:Nspecies)*Zs(1:Nspecies))

  end subroutine

end module ambipolarSolver
