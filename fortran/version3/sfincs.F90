program sfincs

  use sfincs_main
  use mpi

  ! If you want to alter input parameters like Ntheta, you can use a line like the one commented out here:
  !use globalVariables, only: Ntheta

  implicit none

  integer :: ierr

  call MPI_INIT(ierr)

  call sfincs_init(MPI_COMM_WORLD)

  ! To override any parameters in the input namelist, do it here.
  ! Here is an example:
  !Ntheta=17

  call sfincs_prepare()

  ! To alter sfincs internal arrays like BHat and JHat, do it here.

  call sfincs_run()

  ! This next subroutine deallocates all the arrays.
  ! It is not necessary unless you call sfincs multiple times within a single executable.
  call sfincs_finalize()

  call MPI_FINALIZE(ierr)

end program sfincs
