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

  implicit none

  integer, intent(in) :: N, stencil
  PetscScalar, intent(out), dimension(N,N) :: matrix
  PetscScalar, intent(in) :: shift

  PetscScalar :: shift2, value
  PetscScalar, parameter :: pi=3.14159265358979d+0
  integer :: nearestIndex, indexBelow, j, k, firstIndex, index
  PetscScalar, dimension(:), allocatable :: pattern, x

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
  print *,"shift2=",shift2,", nearestIndex=",nearestIndex

  allocate(pattern(N))
  pattern=0
  ! 'pattern' is the first _row_ of the interpolation matrix.

  if (stencil==1) then
     ! If nearestIndex = N+1, change it to 1:
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

end subroutine periodicInterpolation

