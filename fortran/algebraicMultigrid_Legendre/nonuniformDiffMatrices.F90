subroutine nonuniformDiffMatrices(N, x, option, ddx, d2dx2)
  ! Finite difference matrices for a nonuniform grid.
  ! The coefficients here are derived in the note 20170716-02.
  !
  ! Inputs:
  !   integer :: N = number of grid points.
  !   real(prec) :: x(N) = grid point locations.
  !   option = switch for controlling order of accuracy for differentiation
  !            and handling of endpoints.
  !
  ! Presently, the x grid is always assumed to be non-periodic.
  !
  ! Options for option:
  ! 2 =  The domain [xMin, xMax] is assumed to be non-periodic. A 3-point 
  !      stencil is used everywhere.  The first and last row of the
  !      differentiation matrices will use one-sided differences, so they
  !      will each have a non-tridiagonal element.
  ! 3 =  The same as option=2, except that the first differentiation matrix
  !      will use a 2-point 1-sided stencil for the first and last elements
  !      so the matrix is strictly tri-diagonal.  The 2nd derivative matrix
  !      is the same as for option 2, since it is not possible to compute 
  !      the 2nd derivative with only a 2-point stencil.
  ! 12 = The domain [xMin, xMax] is assumed to be non-periodic. A 5-point 
  !      stencil is used everywhere.  The first two and last two rows of 
  !      the differentiation matrices will then each have non-pentadiagonal 
  !      elements.
  ! 13 = The same as option 12, except that 3-point stencils are used for the
  !      first and last rows of the differentiation matrices, and 4-point 
  !      stencils are used for the 2nd and penultimate rows of the 
  !      differentiation matrices.  With this option, both differentiation 
  !      matrices are strictly penta-diagonal.
  ! 32 = Aperiodic.  Upwinding to the left. A 2-point stencil is used for the
  !      first derivative and a 3-point stencil is used for the second 
  !      derivative.  The top row of D and the top two rows of DD are zero.
  ! 42 = Same as 32 but upwinding to the right.  The bottom row of D and the
  !      bottom two rows of DD are zero.
  ! 50 = Periodic with a grid point at xMin but not xMax.  Upwinding to the
  !      left. A 3-point stencil is used for both derivatives.
  ! 51 = Periodic with a grid point at xMax but not xMin.  Upwinding to the
  !      left. A 3-point stencil is used for both derivatives.
  ! 52 = Aperiodic.  Upwinding to the left. A 3-point stencil is used for 
  !      both derivatives.  The top row of both matrices is all zero. The
  !      second row of D uses a 2-point stencil.
  ! 60 = Same as 50 but upwinding to the right.
  ! 61 = Same as 51 but upwinding to the right.
  ! 62 = Same as 52 but upwinding to the right.  The bottom row of both
  !      derivative matrices is all zero. The penultimate row of D uses a
  !      2-point stencil.
  ! 80 = Periodic with a grid point at xMin but not xMax.  The first derivative is upwinded to the
  !      left. A stencil is used with 1 point on 1 side and 2 points on the
  !      other side. The second derivative is the same as in option 0.
  ! 81 = Same as 80 but with a grid point at xMax and not xMin.
  ! 82 = Like 80 but not periodic, with a grid point at both xMin and xMax.
  !      The top row of D is zero.
  ! 90 = Same as 80 but upwinding to the right.
  ! 91 = Same as 90 but with a grid point at xMax and not xMin.
  ! 92 = Like 90 but not periodic, with a grid point at both xMin and xMax.
  !      The bottom row of D is zero.
  ! 100 = Periodic with a grid point at xMin but not xMax.
  !       1st derivative only. Upwinded to the left.
  !       Stencil has 1 point on 1 side, 3 points on the other side
  ! 101 = Same as 100, but with no grid point at xMin and with a grid point at xMax.
  ! 102 = Same as 100, but aperiodic, with grid points at both xMin and xMax.
  ! 110 = Periodic with a grid point at xMin but not xMax.
  !       1st derivative only. Upwinded to the right.
  !       Stencil has 1 point on 1 side, 3 points on the other side
  ! 111 = Same as 110, but with no grid point at xMin and with a grid point at xMax.
  ! 112 = Same as 100, but aperiodic, with grid points at both xMin and xMax.
  ! 120 = Periodic with a grid point at xMin but not xMax.
  !       1st derivative only. Upwinded to the left.
  !       Stencil has 2 points on 1 side, 3 points on the other side
  ! 121 = Same as 120, but with no grid point at xMin and with a grid point at xMax.
  ! 122 = Same as 120, but aperiodic, with grid points at both xMin and xMax.
  !       The first and last rows use at least a 5 point stencil, so there is
  !       no upwinding or backwards upwinding for 3 rows.
  ! 123 = Same as 122, but the first and last rows all are strictly upwinded
  !       so the diagonal is everywhere positive, except for the first row.
  ! 130 = Periodic with a grid point at xMin but not xMax.
  !       1st derivative only. Upwinded to the right.
  !       Stencil has 2 points on 1 side, 3 points on the other side
  ! 131 = Same as 130, but with no grid point at xMin and with a grid point at xMax.
  ! 132 = Same as 130, but aperiodic, with grid points at both xMin and xMax.
  ! 133 = Same as 132, but the first and last rows all are strictly upwinded
  !       so the diagonal is everywhere negative, except for the last row.
  !
  ! Outputs:
  !   ddx = matrix for differentiation.
  !   d2dx2 = matrix for the 2nd derivative.

  use kinds

  implicit none

  integer, intent(in) :: N, option
  real(prec), intent(in), dimension(N) :: x
  real(prec), intent(out), dimension(N,N) :: ddx, d2dx2

  integer :: j
  real(prec) :: x0, x1, x2, x3, x4, x5

  ! ***************************************************************
  ! Validate input
  ! ***************************************************************

  if (N<2) then
     print *,"Error! N must be at least 2."
     print *,"N = ",N
     stop
  end if

  ! ***************************************************************
  ! Fill the interior of the differentiation matrices
  ! ***************************************************************

  ddx = 0d+0
  d2dx2 = 0d+0

  select case (option)
  case (2,3)
     do j = 2, N-1
        x0 = x(j)
        x1 = x(j-1)
        x2 = x(j+1)

        ! First derivative
        a0 = -(-2*x0 + x1 + x2) / ((x0 - x1) * (x0 - x2))
        a1 = -(x0 - x2) / ((x0 - x1) * (x1 - x2))
        a2 = -(x0 - x1) / ((x0 - x2) * (-x1 + x2))

        ddx(j,j+0) = a0
        ddx(j,j-1) = a1
        ddx(j,j+1) = a2

        ! Second derivative
        a0 =  2 / ((x0 - x1) * (x0 - x2))
        a1 = -2 / ((x0 - x1) * (x1 - x2))
        a2 = -2 / ((x0 - x2) * (-x1 + x2))

        d2dx2(j,j+0) = a0
        d2dx2(j,j-1) = a1
        d2dx2(j,j+1) = a2
     end do

  case default
     print *,"Error 1 in nonuniformDiffMatrices. Invalid option:",option
     stop
  end select

  ! ***************************************************************
  ! Handle endpoints of grid in differentiation matrices
  ! ***************************************************************

  select case (option)
  case (2)
     ! ***************************
     ! First derivative
     ! ***************************
     ! First row
     x0 = x(1)
     x1 = x(2)
     x2 = x(3)

     a0 = -(-2*x0 + x1 + x2) / ((x0 - x1) * (x0 - x2))
     a1 = -(x0 - x2) / ((x0 - x1) * (x1 - x2))
     a2 = -(x0 - x1) / ((x0 - x2) * (-x1 + x2))

     ddx(1,1) = a0
     ddx(1,2) = a1
     ddx(1,3) = a2

     ! Last row
     x0 = x(N)
     x1 = x(N-1)
     x2 = x(N-2)

     a0 = -(-2*x0 + x1 + x2) / ((x0 - x1) * (x0 - x2))
     a1 = -(x0 - x2) / ((x0 - x1) * (x1 - x2))
     a2 = -(x0 - x1) / ((x0 - x2) * (-x1 + x2))

     ddx(N,N-1) = a0
     ddx(N,N-2) = a1
     ddx(N,N-3) = a2

     ! ***************************
     ! Second derivative
     ! ***************************
     ! First row
     x0 = x(1)
     x1 = x(2)
     x2 = x(3)

     a0 =  2 / ((x0 - x1) * (x0 - x2))
     a1 = -2 / ((x0 - x1) * (x1 - x2))
     a2 = -2 / ((x0 - x2) * (-x1 + x2))

     d2dx2(1,1) = a0
     d2dx2(1,2) = a1
     d2dx2(1,3) = a2

     ! Last row
     x0 = x(N)
     x1 = x(N-1)
     x2 = x(N-2)

     a0 =  2 / ((x0 - x1) * (x0 - x2))
     a1 = -2 / ((x0 - x1) * (x1 - x2))
     a2 = -2 / ((x0 - x2) * (-x1 + x2))

     d2dx2(N,N-0) = a0
     d2dx2(N,N-1) = a1
     d2dx2(N,N-2) = a2

  case (3)
     ! ***************************
     ! First derivative
     ! ***************************
     ! First row
     ddx(1,1) = -1 / (x(2) - x(1))
     ddx(1,2) =  1 / (x(2) - x(1))

     ! Last row
     ddx(N,N  ) =  1 / (x(N) - x(N-1))
     ddx(N,N-1) = -1 / (x(N) - x(N-1))

     ! ***************************
     ! Second derivative
     ! ***************************
     ! First row
     x0 = x(1)
     x1 = x(2)
     x2 = x(3)

     a0 =  2 / ((x0 - x1) * (x0 - x2))
     a1 = -2 / ((x0 - x1) * (x1 - x2))
     a2 = -2 / ((x0 - x2) * (-x1 + x2))

     d2dx2(1,1) = a0
     d2dx2(1,2) = a1
     d2dx2(1,3) = a2

     ! Last row
     x0 = x(N)
     x1 = x(N-1)
     x2 = x(N-2)

     a0 =  2 / ((x0 - x1) * (x0 - x2))
     a1 = -2 / ((x0 - x1) * (x1 - x2))
     a2 = -2 / ((x0 - x2) * (-x1 + x2))

     d2dx2(N,N-0) = a0
     d2dx2(N,N-1) = a1
     d2dx2(N,N-2) = a2

  case default
     print *,"Error 2 in nonuniformDiffMatrices. Invalid option:", option
     stop
  end select

end subroutine nonuniformDiffMatrices

