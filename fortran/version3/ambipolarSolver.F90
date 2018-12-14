#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

module ambipolarSolver

	! Initial guess for Er
	PetscScalar :: Er_init

  contains

  !! This subroutine begins the search for ambipolarEr by bracketing radialCurrent = zero
  subroutine mainAmbipolarSolver

    use globalVariables
    use solver
    use writeHDF5Output

    implicit none

    integer :: exit_code
    PetscScalar :: ambipolarEr, radialCurrent_min, radialCurrent_max
    PetscScalar, dimension(:), allocatable :: Er_search, radialCurrent
    PetscScalar :: dRadialCurrentdEr_init, time1, time2
    PetscErrorCode :: ierr
		integer :: RHSMode_init

    allocate(Er_search(2))
    allocate(radialCurrent(2))
    Er_search = zero
    radialCurrent = zero

		! If RHSMode>3, perform ambipolarSolve with adjoint diagnostics
		if (RHSMode>3) then
			RHSMode_init = RHSmode
			RHSMode = 1
		end if

		Er_init = Er

		exit_code = 0
		! Check that Er_min and Er_max bracket zero for these ambipolarSolveOption's
		if (ambipolarSolveOption==1 .or. ambipolarSolveOption==2) then
			exit_code = -1
			call updateEr(Er_min,radialCurrent_min)
			if (ambipolarSolveOption==1) then
				dRadialCurrentdEr_init = dRadialCurrentdEr
			end if
			call updateEr(Er_max,radialCurrent_max)

			Er_search(1) = Er_min
			Er_search(2) = Er_max
			radialCurrent(1) = radialCurrent_min
			radialCurrent(2) = radialCurrent_max

			if ((radialCurrent_min<0) .neqv. (radialCurrent_max<0)) then
				exit_code = 0
			end if
		end if

  if (exit_code == 0) then
    print *,"Successful bracketing in ambipolar solve."
    print *,"Here are the Ers we used: "
    print *,Er_search
    print *,"Here are the radial currents: "
    print *,radialCurrent

    call PetscTime(time1, ierr)
    if (ambipolarSolveOption==1) then
      call ambipolarSolverNewton_Bisection(Er_search(1),Er_search(2),radialCurrent(1),radialCurrent(2),dRadialCurrentdEr_init,ambipolarEr)
    else if (ambipolarSolveOption==2) then
      call ambipolarSolverBrent(Er_search(1),Er_search(2),radialCurrent(1),radialCurrent(2),ambipolarEr)
		else if (ambipolarSolveOption==3) then
			call ambipolarSolverNewton(ambipolarEr)
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
    if ((debugAdjoint .eqv. .false.) .and. (RHSMode_init<4)) then
      call updateOutputFile(1, .false.)
    end if

    ! Deallocate stuff
    deallocate(Er_search)
    deallocate(radialCurrent)

    ! If adjoint solve required, call solver again with Er solution
    if (RHSMode_init>3) then
			RHSMode = RHSMode_init
      ambipolarSolve = .false.
      call mainSolverLoop()
      ambipolarSolve = .true.
    end if

  end subroutine mainAmbipolarSolver

  !> Ambipolar root is known to lie between x1 and x2 with radial current values
  !> f1 and f2
	! Based on zbrent subroutine in Numerical Recipes, Chapber 9.3
  subroutine ambipolarSolverBrent(Er1,Er2,f1,f2,ambipolarEr)

    use globalVariables

    implicit none

    PetscScalar :: Er1, Er2, f1, f2, ambipolarEr, radialCurrent_init
    integer :: iEr, exit_code
    PetscScalar :: thisRadialCurrent, Brent_EPS
    PetscScalar, dimension(:), allocatable :: Er_search, radialCurrent
    PetscScalar :: Brent_a, Brent_b, Brent_fa, Brent_fb, Brent_c, Brent_fc, Brent_d, Brent_e, Brent_tol1, Brent_p, Brent_q, Brent_r, Brent_s, Brent_xm
    PetscErrorCode :: ierr
    PetscScalar :: time1, time2


    allocate(Er_search(NEr_ambipolarSolve))
    allocate(radialCurrent(NEr_ambipolarSolve))
    Er_search = zero
    radialCurrent = zero
    exit_code = -1

    Brent_EPS = 1.d-15

		Er_search(1) = Er1
		radialCurrent(1) = f1
		Er_search(2) = Er2
		radialCurrent(2) = f2
		! First function eval at initial guess
		call updateEr(Er_init,radialCurrent_init)
		Er_search(3) = Er_init
		radialCurrent(3) = radialCurrent_init

    ! Initialize Brent algorithm
    Brent_a = Er1
		Brent_b = Er2
    Brent_fa = f1
		Brent_fb = f2
		Brent_c = Er_init
		Brent_fc = radialCurrent_init

    if ((Brent_fa > zero .and. Brent_fb > zero) .or. (Brent_fa < zero .and. Brent_fb < zero)) then
      if (masterProc) then
        print *,"Root must be bracketed in Brent solve!"
      end if
      stop
    end if

    do iEr = 1, (NEr_ambipolarSolve-3)
      call PetscTime(time1, ierr)
      if ((Brent_fb > zero .and. Brent_fc > zero) .or. (Brent_fb < zero .and. Brent_fc < zero)) then
        Brent_c = Brent_a
        Brent_fc = Brent_fa
        Brent_e = Brent_b - Brent_a
        Brent_d = Brent_b - Brent_a
      end if
      if (abs(Brent_fc) < abs(Brent_fb)) then
        Brent_a = Brent_b
        Brent_b = Brent_c
        Brent_c = Brent_a
        Brent_fa = Brent_fb
        Brent_fb = Brent_fc
        Brent_fc = Brent_fa
      end if
      ! Convergence check
      Brent_tol1 = two*Brent_EPS*abs(Brent_b)+0.5*Er_search_tolerance_f
      Brent_xm = 0.5*(Brent_c-Brent_b)

			if (abs(Brent_xm) <= Brent_tol1 .or. (abs(Brent_fb) < Er_search_tolerance_f)) then
        ambipolarEr = Brent_b
        exit_code = 0
        exit
      end if

      if (abs(Brent_e) >= Brent_tol1 .and. abs(Brent_fa)>abs(Brent_fb)) then
         ! Attempt inverse quadratic interpolation
         Brent_s = Brent_fb / Brent_fa
         if (Brent_a == Brent_c) then
            Brent_p = 2.0*Brent_xm*Brent_s
            Brent_q = 1.0-Brent_s
         else
            Brent_q = Brent_fa / Brent_fc
            Brent_r = Brent_fb / Brent_fc
            Brent_p = Brent_s*(2.0*Brent_xm*Brent_q*(Brent_q-Brent_r)-(Brent_b-Brent_a)*(Brent_r-1.0))
            Brent_q = (Brent_q-1.0)*(Brent_r-1.0)*(Brent_s-1.0)
         end if
         if (Brent_p > 0) Brent_q = -Brent_q ! Check whether in bounds
         Brent_p = abs(Brent_p)
         if (2.0*Brent_p < min(3.0*Brent_xm*Brent_q-abs(Brent_tol1*Brent_q),abs(Brent_e*Brent_q))) then
           ! Accept interpolation
           Brent_e = Brent_d
           Brent_d = Brent_p / Brent_q
         else
           ! Interpolation failed, so use bisection
           Brent_d = Brent_xm
           Brent_e = Brent_d
         end if
      else
         ! Bounds are decreasing too slowly, so use bisection
         Brent_d = Brent_xm
         Brent_e = Brent_d
      end if
      ! Move last best guess to a.
      Brent_a = Brent_b
      Brent_fa = Brent_fb
      if (abs(Brent_d) > Brent_tol1) then
         ! Evaluate new trial root
         Brent_b = Brent_b + Brent_d
      else
         Brent_b = Brent_b + sign(Brent_tol1,Brent_xm)
      end if
      call updateEr(Brent_b, Brent_fb)
      radialCurrent(iEr+3) = Brent_fb
      Er_search(iEr+3) = Brent_b
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
  subroutine ambipolarSolverNewton_Bisection(Er1, Er2, fl, fh, df, ambipolarEr)

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

		print *, "Er1: ", Er1
		print *, "Er2: ", Er2

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
		rts = Er_init
		call updateEr(rts,f)
		df = dRadialCurrentdEr
    dxold = abs(Er2-Er1)
    dx = dxold
		Er_search(1) = rts
		radialCurrent(1) = df

    exit_code = -1
    ! Loop over allowed iterations
    do j = 2, NEr_ambipolarSolve
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
      if (abs(dx) < Er_search_tolerance_dx .or. (abs(f) < Er_search_tolerance_f)) then
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
        print *,"Time for one Newton bisection loop: ", time2-time1," sec."
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

  end subroutine ambipolarSolverNewton_Bisection

	subroutine ambipolarSolverNewton(ambipolarEr)

		use globalVariables
		use radialCoordinates
		use solver

		implicit none

		PetscScalar, intent(out) :: ambipolarEr
		PetscScalar, dimension(:), allocatable :: Er_search, radialCurrent
		integer :: j, exit_code
		PetscErrorCode :: ierr
		PetscScalar :: time1, time2, this_Er, this_RadialCurrent, dx

		allocate(Er_search(NEr_ambipolarSolve))
		allocate(radialCurrent(NEr_ambipolarSolve))
		Er_search = zero
		radialCurrent = zero

		exit_code = -1
		! Initial guess
		this_Er = Er_init

		do j = 1, NEr_ambipolarSolve
			call PetscTime(time1,ierr)
			! The one new function evaluation per iteration
			call updateEr(this_Er,this_RadialCurrent)
			Er_search(j) = this_Er
			radialCurrent(j) = this_RadialCurrent
			dx = this_RadialCurrent/dRadialCurrentdEr ! this comes from globalVariables
			this_Er = this_Er - dx
			if (this_Er < Er_min .or. this_Er > Er_max) then
				print *,"Newton method iterate exceeded Er bounds!"
				exit_code = -2
				exit
			end if
			if ((abs(dx) < Er_search_tolerance_dx) .or. (abs(this_radialCurrent)< Er_search_tolerance_f)) then
				ambipolarEr = this_Er
				exit_code = 0
				exit
			end if
			call PetscTime(time2,ierr)
			if (masterProc) then
				print *,"Time for one Newton loop: ", time2-time1," sec."
			end if
		end do ! j

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

    if (exit_code == -2) then
      if (masterProc) then
         print *,"*******************************************************************************"
         print *,"*******************************************************************************"
         print *,"The Er search exceeded bounds!"
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
    PetscScalar, dimension(:), allocatable :: array
    integer :: rank, master_rank

    allocate(array(1))

    Er = thisEr

    ! Update dPhiHatdpsiHat, dPhiHatdpsiN, dPhiHatdrHat, dPhiHatdrN
    dPhiHatdpsiHat = ddrHat2ddpsiHat * (-Er)
    dPhiHatdpsiN = ddpsiHat2ddpsiN * dPhiHatdpsiHat
    dPhiHatdrHat = ddpsiHat2ddrHat * dPhiHatdpsiHat
    dPhiHatdrN   = ddpsiHat2ddrN   * dPhiHatdpsiHat

    call mainSolverLoop()

    if (masterProc) then
      array(1) = sum(particleFlux_vm_rN(1:Nspecies)*Zs(1:Nspecies))
    end if
    ! Broadcast value of radialCurrent to all procs
    call MPI_BCAST(array,1,MPI_DOUBLE,0,MPIComm,ierr)
    radialCurrent = array(1)

  end subroutine

end module ambipolarSolver
