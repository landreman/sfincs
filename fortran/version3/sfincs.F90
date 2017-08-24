program sfincs

  use sfincs_main
  use mpi

  ! If you want to alter input parameters like Ntheta, you can use a line like the one commented out here:
  !use globalVariables, only: Ntheta

  implicit none

  integer :: ierr

  call MPI_INIT(ierr)

  call init_sfincs(MPI_COMM_WORLD)

  ! To override any parameters in the input namelist, do it here.
  ! Here is an example:
  !Ntheta=17

  call prepare_sfincs()

  ! To alter sfincs internal arrays like BHat and JHat, do it here.

  call run_sfincs()

  call MPI_FINALIZE(ierr)

end program sfincs
