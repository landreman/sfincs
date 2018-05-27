module radialCoordinates

  use globalVariables

  implicit none

contains

  ! -----------------------------------------------------------------------------------

  subroutine setInputRadialCoordinateWish

    ! This subroutine takes the selected "wish" radial coordinate, and uses it to over-write
    ! the other "wish" radial coordinates.

    implicit none

    if (aHat<0) then
       print *,"Error! aHat cannot be <0."
       stop
    end if

    ! First, set psiHat_wish itself:

    select case (inputRadialCoordinate)
    case (0)
       ! Selected input radial coordinate is psiHat.
       ! Nothing to do here.

       if (masterProc) then
          print *,"Selecting the flux surface to use based on psiHat_wish = ",psiHat_wish
       end if
       if (psiHat_wish / psiAHat < 0) then
          print *,"Error! psiHat_wish/psiAHat cannot be <0."
          stop
       end if
       if (psiHat_wish / psiAHat > 1) then
          print *,"Error! psiHat_wish/psiAHat cannot be >1."
          stop
       end if

    case (1)
       ! Selected input radial coordinate is psiN.
       psiHat_wish = psiN_wish * psiAHat

       if (masterProc) then
          print *,"Selecting the flux surface to use based on psiN_wish = ",psiN_wish
       end if
       if (psiN_wish<0) then
          print *,"Error! psiN_wish cannot be <0."
          stop
       end if
       if (psiN_wish>1) then
          print *,"Error! psiN_wish cannot be >1."
          stop
       end if

    case (2)
       ! Selected input radial coordinate is rHat.
       psiHat_wish = psiAHat * rHat_wish * rHat_wish / (aHat * aHat)

       if (masterProc) then
          print *,"Selecting the flux surface to use based on rHat_wish = ",rHat_wish
       end if
       if (rHat_wish<0) then
          print *,"Error! rHat_wish cannot be <0."
          stop
       end if
       if (rHat_wish>aHat) then
          print *,"Error! rHat_wish cannot be > aHat."
          stop
       end if

    case (3)
       ! Selected input radial coordinate is rN.
       psiHat_wish = rN_wish * rN_wish * psiAHat

       if (masterProc) then
          print *,"Selecting the flux surface to use based on rN_wish = ",rN_wish
       end if
       if (rN_wish<0) then
          print *,"Error! rN_wish cannot be <0."
          stop
       end if
       if (rN_wish>1) then
          print *,"Error! rN_wish cannot be >1."
          stop
       end if

    case default
       print *,"Error! Invalid inputRadialCoordinate."
       stop
    end select

    ! Now, use psiHat_wish to set the others:
    
    psiN_wish = psiHat_wish / psiAHat
    rHat_wish = sqrt(aHat * aHat * psiHat_wish / psiAHat)
    rN_wish = sqrt(psiHat_wish / psiAHat)

    ! Validate input again, just to be safe:

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


    ! Conversion factors for converting TO d/dpsiHat:
    ddpsiN2ddpsiHat =  one / psiAHat
    ddrHat2ddpsiHat = aHat / (two * psiAHat * sqrt(psiN))
    ddrN2ddpsiHat   =  one / (two * psiAHat * sqrt(psiN))

    ! Conversion factors for converting FROM d/dpsiHat:
    ddpsiHat2ddpsiN = psiAHat
    ddpsiHat2ddrHat = (two * psiAHat * sqrt(psiN)) / aHat
    ddpsiHat2ddrN   = (two * psiAHat * sqrt(psiN))

    ! Next, set the d/dpsiHat quantities related to input gradients:
    if (RHSMode==3) then
       ! Monoenergetic coefficient computation.
       ! Nothing to do here.
    else
       select case (inputRadialCoordinateForGradients)
       case (0)
          ! Selected input radial coordinate is psiHat.
          ! Nothing to do here.
          
          if (masterProc) then
             print *,"Selecting the input gradients (of n, T, & Phi) from the specified ddpsiHat values."
          end if
          
       case (1)
          ! Selected input radial coordinate is psiN.
          dPhiHatdpsiHat = ddpsiN2ddpsiHat * dPhiHatdpsiN
          dnHatdpsiHats  = ddpsiN2ddpsiHat * dnHatdpsiNs
          dTHatdpsiHats  = ddpsiN2ddpsiHat * dTHatdpsiNs
          
          if (masterProc) then
             print *,"Selecting the input gradients (of n, T, & Phi) from the specified ddpsiN values."
          end if
          
       case (2)
          ! Selected input radial coordinate is rHat.
          dPhiHatdpsiHat = ddrHat2ddpsiHat * dPhiHatdrHat
          dnHatdpsiHats  = ddrHat2ddpsiHat * dnHatdrHats
          dTHatdpsiHats  = ddrHat2ddpsiHat * dTHatdrHats
          
          if (masterProc) then
             print *,"Selecting the input gradients (of n, T, & Phi) from the specified ddrHat values."
          end if
          
       case (3)
          ! Selected input radial coordinate is rN.
          dPhiHatdpsiHat = ddrN2ddpsiHat * dPhiHatdrN
          dnHatdpsiHats  = ddrN2ddpsiHat * dnHatdrNs
          dTHatdpsiHats  = ddrN2ddpsiHat * dTHatdrNs
          
          if (masterProc) then
             print *,"Selecting the input gradients (of n, T, & Phi) from the specified ddrN values."
          end if
          
       case (4)
          ! Selected input radial coordinate is rHat, and use Er instead of dPhiHatdrHat:
          dPhiHatdpsiHat = ddrHat2ddpsiHat * (-Er)
          dnHatdpsiHats  = ddrHat2ddpsiHat * dnHatdrHats
          dTHatdpsiHats  = ddrHat2ddpsiHat * dTHatdrHats
          
          if (masterProc) then
             print *,"Selecting the input gradients of n & T from the specified ddrHat values."
             print *,"Selecting the input gradient of Phi from the specified Er."
          end if
          
       case default
          print *,"Error! Invalid inputRadialCoordinateForGradients."
          stop
       end select
    end if

    ! Finally, convert the input gradients from psiHat to all the other radial coordinates:

    dPhiHatdpsiN = ddpsiHat2ddpsiN * dPhiHatdpsiHat
    dPhiHatdrHat = ddpsiHat2ddrHat * dPhiHatdpsiHat
    dPhiHatdrN   = ddpsiHat2ddrN   * dPhiHatdpsiHat
    Er = - dPhiHatdrHat
    
    dnHatdpsiNs  = ddpsiHat2ddpsiN * dnHatdpsiHats
    dnHatdrHats  = ddpsiHat2ddrHat * dnHatdpsiHats
    dnHatdrNs    = ddpsiHat2ddrN   * dnHatdpsiHats

    dTHatdpsiNs  = ddpsiHat2ddpsiN * dTHatdpsiHats
    dTHatdrHats  = ddpsiHat2ddrHat * dTHatdpsiHats
    dTHatdrNs    = ddpsiHat2ddrN   * dTHatdpsiHats

  end subroutine setInputRadialCoordinate


end module radialCoordinates

