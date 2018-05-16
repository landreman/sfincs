program sfincs_testAdjointDiagnostics

    use sfincs_main
    use testingAdjointDiagnostics
    use mpi

    implicit none

    integer :: ierr

    call MPI_INIT(ierr)

    call sfincs_init(MPI_COMM_WORLD)

    call sfincs_prepare()

    call compareAdjointDiagnostics()

    call sfincs_finalize()

    call MPI_FINALIZE(ierr)

end program sfincs_testAdjointDiagnostics
