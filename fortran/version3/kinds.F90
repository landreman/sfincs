module kinds

  implicit none

  ! Double precision, for calculations that are always done in double even when PETSc is compiled with single-precision
  integer, parameter :: prec = SELECTED_REAL_KIND(12,100)

end module kinds

