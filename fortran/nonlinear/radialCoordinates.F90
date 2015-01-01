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

    case (1)
       ! Selected input radial coordinate is psiN.
       psiHat_wish = psiN_wish * psiAHat

    case (2)
       ! Selected input radial coordinate is rHat.
       psiHat_wish = psiAHat * rHat_wish * rHat_wish / (aHat * aHat)

    case (3)
       ! Selected input radial coordinate is rN.
       psiHat_wish = rN_wish * rN_wish * psiAHat

    case default
       print *,"Error! Invalid inputRadialCoordinate."
       stop
    end select

    ! Now, use psiHat_wish to set the others:
    psiN_wish = psiHat_wish / psiAHat
    rHat_wish = sqrt(aHat * aHat * psiHat_wish / psiAHat)
    rN_wish = sqrt(psiHat_wish / psiAHat)

  end subroutine setInputRadialCoordinateWish

  ! -----------------------------------------------------------------------------------

  subroutine setInputRadialCoordinate

    ! For input quantities that depend on the radial coordinate, pick out the values for the selected
    ! radial coordinate, and use these values to over-write values for the other radial coordinates.

    implicit none

    integer :: i

    ! At this point, rN should have been set by computeBHat!


    ! Next, set the d/dpsiHat quantities:

    select case (inputRadialCoordinate)
    case (0)
       ! Selected input radial coordinate is psiHat.
       ! Nothing to do here.

    case (1)
       ! Selected input radial coordinate is psiN.
       dPhiHatdpsiHat = ddpsiN2ddpsiHat(dPhiHatdpsiN)
       do i=1,Nspecies
          dnHatdpsiHats(i) = ddpsiN2ddpsiHat(dnHatdpsiNs(i))
          dTHatdpsiHats(i) = ddpsiN2ddpsiHat(dTHatdpsiNs(i))
       end do

    case (2)
       ! Selected input radial coordinate is rHat.
       dPhiHatdpsiHat = ddrN2ddpsiHat(dPhiHatdrHat)
       do i=1,Nspecies
          dnHatdpsiHats(i) = ddrHat2ddpsiHat(dnHatdrHats(i))
          dTHatdpsiHats(i) = ddrHat2ddpsiHat(dTHatdrHats(i))
       end do

    case (3)
       ! Selected input radial coordinate is rN.
       dPhiHatdpsiHat = ddrN2ddpsiHat(dPhiHatdrN)
       do i=1,Nspecies
          dnHatdpsiHats(i) = ddrN2ddpsiHat(dnHatdrNs(i))
          dTHatdpsiHats(i) = ddrN2ddpsiHat(dTHatdrNs(i))
       end do

    case default
       print *,"Error! Invalid inputRadialCoordinate."
       stop
    end select


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

end module geometry

