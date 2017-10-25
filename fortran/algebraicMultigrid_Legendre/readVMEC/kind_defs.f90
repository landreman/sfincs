!DK kind_defs
!***********************************************************************
!! kind definitions commonly used
      module kind_defs
      implicit none
      integer, parameter :: i4b = selected_int_kind(9)
      integer, parameter :: i2b = selected_int_kind(4)
      integer, parameter :: i1b = selected_int_kind(2)
!     integer, parameter :: srp = kind(1.0)
      integer, parameter :: dp = kind(1.0d0)
!     integer, parameter :: srpc = kind((1.0,1.0))
      integer, parameter :: dpc = kind((1.0d0,1.0d0))
!     integer, parameter :: lgt = kind(.true.)
!     REAL(SRP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_srp
!     REAL(SRP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_srp
!     REAL(SRP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_srp
!     REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
!     REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
!     REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp

      end module kind_defs

