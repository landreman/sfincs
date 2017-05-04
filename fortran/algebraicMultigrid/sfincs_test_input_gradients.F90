program sfincs

  use sfincs_main

  ! If you want to alter input parameters like Ntheta, you can use a line like the one commented out here:
  !use globalVariables, only: Ntheta

  implicit none

  call init_sfincs()

  ! To override any parameters in the input namelist, do it here.

  Ntheta=21
  Nzeta=1
  Nxi=10
  Nx=12

  aHat=0
  psiAHat=0

  inputRadialCoordinateForGradients = 0 ! This line says input gradients are set by d*dpsiHat
  Nspecies = 2

  Zs(1) = 1.6
  Zs(2) = 2.9
  
  mHats(1) = 1.9
  mHats(2) = 2.4

  nHats(1) = 3.1
  nHats(2) = 1.8

  THats(1) = 0.5
  THats(2) = 0.3

  dnHatdpsiHats(1) = -22.2222
  dnHatdpsiHsts(2) = -15.5556

  dTHatdpsiHats(1) = -11.1111
  dTHatdpsiHats(2) = -8.14815

  dPhiHatdpsiHat = 0


  call prepare_sfincs()

  ! To alter sfincs internal arrays like BHat and JHat, do it here.

  call run_sfincs()

end program sfincs
