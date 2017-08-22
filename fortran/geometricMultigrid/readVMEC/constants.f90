!DK constants
!***********************************************************************
!! contains commonly used constants
      module constants
!***********************************************************************
      use kind_defs

      real(DP), parameter ::      pi = 3.14159265358979323846264_dp
      real(DP), parameter ::   twopi = 6.28318530717958647692528_dp
      real(DP), parameter ::     muo = 1.25663706143591729538506d-06
      real(DP), parameter :: sqrtmuo = 1.1209982432795857d-03
      real(DP), parameter ::     cp0 = 0.0_dp  !for constant parameter
      real(DP), parameter ::     cp1 = 1.0_dp
      real(DP), parameter ::     cp2 = 2.0_dp
      real(DP), parameter ::    cpp5 = 0.5_dp

      save

      end module constants

