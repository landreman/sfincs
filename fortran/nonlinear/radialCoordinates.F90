module radialCoordinates

  use globalVariables
  use petscsysdef

  implicit none

#include <finclude/petscsysdef.h>

contains

  ! -----------------------------------------------------------------------------------

  subroutine setInputRadialCoordinateWish

    ! This subroutine takes the selected "wish" radial coordinate, and uses it to over-write
    ! the other "wish" radial coordinates.

    implicit none

    ! First, set psiHat_wish itself:

    select case (inputRadialCoordinate)
    case (0)
       ! Selected input radial coordinate is psiHat.
       ! Nothing to do here.

       if (masterProc) then
          print *,"Selecting the flux surface to use based on psiHat_wish = ",psiHat_wish
       end if

    case (1)
       ! Selected input radial coordinate is psiN.
       psiHat_wish = psiN_wish * psiAHat

       if (masterProc) then
          print *,"Selecting the flux surface to use based on psiN_wish = ",psiN_wish
       end if

    case (2)
       ! Selected input radial coordinate is rHat.
       psiHat_wish = psiAHat * rHat_wish * rHat_wish / (aHat * aHat)

       if (masterProc) then
          print *,"Selecting the flux surface to use based on rHat_wish = ",rHat_wish
       end if

    case (3)
       ! Selected input radial coordinate is rN.
       psiHat_wish = rN_wish * rN_wish * psiAHat

       if (masterProc) then
          print *,"Selecting the flux surface to use based on rN_wish = ",rN_wish
       end if

    case default
       print *,"Error! Invalid inputRadialCoordinate."
       stop
    end select

    ! Now, use psiHat_wish to set the others:

    psiN_wish = psiHat_wish / psiAHat
    rHat_wish = sqrt(aHat * aHat * psiHat_wish / psiAHat)
    rN_wish = sqrt(psiHat_wish / psiAHat)

    ! Validate input:

    if (rN_wish<0) then
       print *,"Error! Requested flux surface has a negative normalized radius.  (rN < 0)"
       stop
    end if

    if (rN_wish>1) then
       print *,"Error! Requested flux surface is outside the last closed flux surface.  (rN > 1)"
       stop
    end if

  end subroutine setInputRadialCoordinateWish

  ! -----------------------------------------------------------------------------------

  subroutine setInputRadialCoordinate

    ! For input gradients, pick out the values for the selected
    ! radial coordinate, and use these values to over-write values for the other radial coordinates.

    implicit none

    integer :: i

    ! At this point, rN should have been set by computeBHat!

    if (abs(rN+9999) < 1e-6) then
       print *,"Error! It appears that rN was not set by computeBHat() in the geometry module."
       stop
    end if

    if (rN < 0) then
       print *,"Error! rN is negative. Probably an error occured while processing the magnetic geometry."
       stop
    end if

    if (rN > 1) then
       print *,"Error! Selected flux surface is outside the last closed flux surface.  (rN > 1.)"
       print *,"Probably an error occured while processing the magnetic geometry."
       stop
    end if

    ! Now set the other flux surface label coordinates using rN:
    psiHat = psiAHat * rN * rN
    psiN = rN * rN
    rHat = aHat * rN

    if (masterProc) then
       print *,"Requested/actual flux surface for this calculation, in various radial coordinates:"
       print *,"  psiHat = ",psiHat_wish,'/',psiHat
       print *,"  psiN   = ",psiN_wish,  '/',psiN
       print *,"  rHat   = ",rHat_wish,  '/',rHat
       print *,"  rN     = ",rN_wish,    '/',rN
    end if


    ! Next, set the d/dpsiHat quantities:

    select case (inputRadialCoordinateForGradients)
    case (0)
       ! Selected input radial coordinate is psiHat.
       ! Nothing to do here.

       if (masterProc) then
          print *,"Selecting the input gradients (of n, T, & Phi) from the specified ddpsiHat values."
       end if

    case (1)
       ! Selected input radial coordinate is psiN.
       dPhiHatdpsiHat = ddpsiN2ddpsiHat(dPhiHatdpsiN)
       do i=1,Nspecies
          dnHatdpsiHats(i) = ddpsiN2ddpsiHat(dnHatdpsiNs(i))
          dTHatdpsiHats(i) = ddpsiN2ddpsiHat(dTHatdpsiNs(i))
       end do

       if (masterProc) then
          print *,"Selecting the input gradients (of n, T, & Phi) from the specified ddpsiN values."
       end if

    case (2)
       ! Selected input radial coordinate is rHat.
       dPhiHatdpsiHat = ddrN2ddpsiHat(dPhiHatdrHat)
       do i=1,Nspecies
          dnHatdpsiHats(i) = ddrHat2ddpsiHat(dnHatdrHats(i))
          dTHatdpsiHats(i) = ddrHat2ddpsiHat(dTHatdrHats(i))
       end do

       if (masterProc) then
          print *,"Selecting the input gradients (of n, T, & Phi) from the specified ddrHat values."
       end if

    case (3)
       ! Selected input radial coordinate is rN.
       dPhiHatdpsiHat = ddrN2ddpsiHat(dPhiHatdrN)
       do i=1,Nspecies
          dnHatdpsiHats(i) = ddrN2ddpsiHat(dnHatdrNs(i))
          dTHatdpsiHats(i) = ddrN2ddpsiHat(dTHatdrNs(i))
       end do

       if (masterProc) then
          print *,"Selecting the input gradients (of n, T, & Phi) from the specified ddrN values."
       end if

    case default
       print *,"Error! Invalid inputRadialCoordinateForGradients."
       stop
    end select


    ! Finally, convert the input gradients from psiHat to all the other radial coordinates:

    dPhiHatdpsiN = ddpsiHat2ddpsiN(dPhiHatdpsiHat)
    dPhiHatdrHat = ddpsiHat2ddrHat(dPhiHatdpsiHat)
    dPhiHatdrN   = ddpsiHat2ddrN(  dPhiHatdpsiHat)
    
    do i=1,Nspecies   
       dnHatdpsiNs(i)  = ddpsiHat2ddpsiN(dnHatdpsiHats(i))
       dnHatdrHats(i)  = ddpsiHat2ddrHat(dnHatdpsiHats(i))
       dnHatdrNs(i)    = ddpsiHat2ddrN(  dnHatdpsiHats(i))

       dTHatdpsiNs(i)  = ddpsiHat2ddpsiN(dTHatdpsiHats(i))
       dTHatdrHats(i)  = ddpsiHat2ddrHat(dTHatdpsiHats(i))
       dTHatdrNs(i)    = ddpsiHat2ddrN(  dTHatdpsiHats(i))
    end do

  end subroutine setInputRadialCoordinate

  ! -----------------------------------------------------------------------------------------

  PetscScalar function ddpsiN2ddpsiHat(x)

    implicit none

    PetscScalar, intent(in) :: x

    ddpsiN2ddpsiHat = x / psiAHat

  end function ddpsiN2ddpsiHat

  ! -----------------------------------------------------------------------------------------

  PetscScalar function ddrHat2ddpsiHat(x)

    implicit none

    PetscScalar, intent(in) :: x

    ddrHat2ddpsiHat = x / (2 * psiAHat * sqrt(psiN))

  end function ddrHat2ddpsiHat

  ! -----------------------------------------------------------------------------------------

  PetscScalar function ddrN2ddpsiHat(x)

    implicit none

    PetscScalar, intent(in) :: x

    ddrN2ddpsiHat = x * aHat / (2 * psiAHat * sqrt(psiN))

  end function ddrN2ddpsiHat

  ! -----------------------------------------------------------------------------------------

  PetscScalar function ddpsiHat2ddpsiN(x)

    implicit none

    PetscScalar, intent(in) :: x

    ddpsiHat2ddpsiN = x * psiAHat

  end function ddpsiHat2ddpsiN

  ! -----------------------------------------------------------------------------------------

  PetscScalar function ddpsiHat2ddrHat(x)

    implicit none

    PetscScalar, intent(in) :: x

    ddpsiHat2ddrHat = x * (2 * psiAHat * sqrt(psiN))

  end function ddpsiHat2ddrHat

  ! -----------------------------------------------------------------------------------------

  PetscScalar function ddpsiHat2ddrN(x)

    implicit none

    PetscScalar, intent(in) :: x

    ddpsiHat2ddrN = x / aHat * (2 * psiAHat * sqrt(psiN))

  end function ddpsiHat2ddrN

end module radialCoordinates

