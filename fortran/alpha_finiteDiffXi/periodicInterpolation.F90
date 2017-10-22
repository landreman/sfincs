#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

subroutine periodicInterpolation(N, matrix, shift, stencil)
  ! This subroutine constructs a matrix for interpolating a periodic function
  ! known on N uniformly spaced grid points on the domain [0,2*pi).
  ! 'matrix' should have been allocated with size (N,N).

  use rank_mod

  implicit none

  integer, intent(in) :: N, stencil
  PetscScalar, intent(out), dimension(N,N) :: matrix
  PetscScalar, intent(in) :: shift

  PetscScalar :: shift2, value
  PetscScalar, parameter :: pi=3.14159265358979d+0
  integer :: nearestIndex, indexBelow, j, k, firstIndex, index, N_possible_shifts, shift_to_try
  PetscScalar, dimension(:), allocatable :: pattern, x, possible_shifts, errors
  integer, dimension(:), allocatable :: indices
  logical :: found_a_solution

  ! Validate:
  if (stencil <= 0) then
     print *,"Error! Stencil must be positive, but stencil=",stencil
     stop
  end if
  if (stencil > N) then
     print *,"Error! Stencil cannot exceed N. stencil=",stencil,", N=",N
     stop
  end if


  ! Ensure 'shift' is between [0,2*pi).
  ! Note that we must use 'modulo' and not 'mod', since the latter returns negative results if shift<0.
  shift2 = modulo(shift,2*pi)

  ! Only used if stencil is even. Note that 1 <= indexBelow <= N.
  indexBelow = floor(N*shift2/(2*pi)) + 1
  ! Only used if stencil is odd. Note that 1 <= nearestIndex <= N+1.
  nearestIndex = nint(N*shift2/(2*pi)) + 1
  !print *,"shift2=",shift2,", nearestIndex=",nearestIndex

  allocate(pattern(N))
  pattern=0
  ! 'pattern' is the first _row_ of the interpolation matrix.

  if (stencil==1) then
     ! Here we alter nearestIndex to ensure that the amount of shift (nearestIndex-1)
     ! is relatively prime with N, since otherwise subsets of the alpha grid decouple from each other,
     ! yielding a singular matrix. (It is like having a rational iota.)

     if (gcd(N,nearestIndex-1) .ne. 1) then
        print *,"*******************************************************************"
        print *,"*******************************************************************"
        print *,"Warning: The nearest-integer shift in alpha at the ends of the zeta"
        print *,"domain for this value of iota is not relatively prime with Nalpha."
        print *,"This will cause KSP/GMRES to require more iterations. It is"
        print *,"recommended that you change to either a higher or lower Nalpha."
        print *,"*******************************************************************"
        print *,"*******************************************************************"
     end if

     print *,"Original alpha shift:",nearestIndex-1

     N_possible_shifts = 3*N  ! Almost certainly more than necessary...
     allocate(possible_shifts(N_possible_shifts))
     allocate(errors(N_possible_shifts))
     allocate(indices(N_possible_shifts))

     possible_shifts = [( j-N, j=1,N_possible_shifts )]
     errors = [( abs(N*shift2/(2*pi) - possible_shifts(j)), j=1,N_possible_shifts )]
     call rank(errors,indices)
     print *,"possible_shifts:",possible_shifts
     print *,"errors:",errors
     print *,"indices:",indices
     ! The 'rank' subroutine returns results in descending order.
     found_a_solution = .false.
     do j=1,N_possible_shifts
        shift_to_try = possible_shifts(indices(j))
        print *,"gcd(",N,",",shift_to_try,")=",gcd(N,shift_to_try)
        if (gcd(N,shift_to_try)==1) then
           ! If the greatest common divisor==1, then this shift is relatively prime to N,
           ! so accept this shift.
           found_a_solution = .true.
           nearestIndex = shift_to_try + 1
           exit
        end if
     end do
     if (.not. found_a_solution) then
        print *,"Error in periodicInterpolation for stencil=1."
     end if

     deallocate(possible_shifts,errors,indices)
     print *,"Adjusted alpha shift:",nearestIndex-1

     ! Ensure nearestIndex is in the range [1,N]:
     index = modulo(nearestIndex-1,N)+1
     pattern(index)=1
  else

     if (mod(stencil,2)==0) then
        ! Stencil is even.
        firstIndex = indexBelow - (stencil/2) + 1
     else
        ! Stencil is odd.
        firstIndex = nearestIndex - ((stencil-1)/2)
     end if

     allocate(x(stencil))
     x = [( ((firstIndex+j-2)*2*pi)/N, j=1,stencil )]

     do j=1,stencil
        value = 1
        do k=1,stencil
           if (k .ne. j) then
              value = value * (shift2 - x(k)) / (x(j) - x(k))
           end if
        end do
        index = modulo(firstIndex+j-2,N)+1
        pattern(index) = value
     end do

     deallocate(x)
  end if

  do j=1,N
     do k=1,N
        index = modulo(k-j,N)+1
        matrix(j,k) = pattern(index)
     end do
  end do

  deallocate(pattern)

contains

  ! Greatest common divisor
  function gcd(v, t)
    integer :: gcd
    integer, intent(in) :: v, t
    integer :: c, b, a
    
    b = t
    a = v
    do
       c = mod(a, b)
       if ( c == 0) exit
       a = b
       b = c
    end do
    gcd = b ! abs(b)
    
  end function gcd

end subroutine periodicInterpolation
