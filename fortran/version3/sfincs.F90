program sfincs

  use globalVariables
  use sfincs_main

  ! If you want to alter input parameters like Ntheta, you can use a line like the one commented out here:
  !use globalVariables, only: Ntheta

  implicit none

  call init_sfincs()

  ! To override any parameters in the input namelist, do it here.
  ! Here is an example:
  !Ntheta=17
  call prepare_sfincs()

  ! To alter sfincs internal arrays like BHat and JHat, do it here.

  call run_sfincs()

  call finalize_sfincs()



end program sfincs
