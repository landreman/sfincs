module export_f

  ! This module contains subroutines and variables related to exporting the distribution function f.
  ! The exported f is available on various grids that can differ from the grid used for solving the kinetic equation.

  use globalVariables

  implicit none

#include <finclude/petscsysdef.h>

  integer, parameter :: max_N_export_f = 500 ! Max grid size in each one of the 4 coordinates (theta, zeta, xi, x)
  PetscScalar, dimension(max_N_export_f) :: export_f_theta
  PetscScalar, dimension(max_N_export_f) :: export_f_zeta
  PetscScalar, dimension(max_N_export_f) :: export_f_xi
  PetscScalar, dimension(max_N_export_f) :: export_f_x
  integer :: N_export_f_theta, N_export_f_zeta, N_export_f_xi, N_export_f_x
  PetscScalar, dimension(:,:), allocatable :: map_theta_to_export_f_theta
  PetscScalar, dimension(:,:), allocatable :: map_zeta_to_export_f_zeta
  PetscScalar, dimension(:,:), allocatable :: map_x_to_export_f_x
  PetscScalar, dimension(:,:), allocatable :: map_xi_to_export_f_xi
  integer :: export_f_theta_option=0, export_f_zeta_option=0, export_f_xi_option=0, export_f_x_option=0
  logical :: export_full_f = .false., export_delta_f = .false.

  PetscScalar, dimension(:,:,:,:,:), allocatable :: delta_f, full_f

  contains

    subroutine setup_grids_for_export_f

      ! This subroutine does some preparatory work for exporting the distribution function. 
      ! It is called from createGrids.F90 after the theta/zeta/x grids have been assembled.

      implicit none

      integer :: i, j, L, indexOfLeastError, index, index1, index2
      logical, dimension(:), allocatable :: includeThisTheta
      logical, dimension(:), allocatable :: includeThisZeta
      logical, dimension(:), allocatable :: includeThisX
      PetscScalar :: error, leastError, weight1, weight2
      PetscScalar, dimension(:,:), allocatable :: extrapMatrix, regridPolynomialToUniform_plus1
      PetscScalar, dimension(:), allocatable :: x_plus1

      ! --------------------------------------------------------
      ! Handle theta coordinate
      ! --------------------------------------------------------

      select case (export_f_theta_option)
      case (0)
         ! 0: Report the distribution function on the original theta grid used internally by SFINCS (with Ntheta elements).
         N_export_f_theta = Ntheta

         if (Ntheta>max_N_export_f) then
            ! Very unlikely, but handle it just in case, to prevent seg fault in a moment when setting export_f_theta.
            if (masterProc) then
               print *,"Error! When export_f_theta_option=0, max_N_export_f must be >= Ntheta."
            end if
            stop
         end if
         export_f_theta(1:Ntheta) = theta

         ! Build the Ntheta * Ntheta identity matrix:
         allocate(map_theta_to_export_f_theta(Ntheta,Ntheta))
         map_theta_to_export_f_theta = zero
         do j=1,Ntheta
            map_theta_to_export_f_theta(j,j)=1
         end do

      case (1)
         ! 1: Interpolate to a different grid, specified by export_f_theta. 

         ! Validation:
         if (N_export_f_theta < 1) then
            if (masterProc) then
               print *,"Error! When export_f_theta_option=1, you must specify at least 1 value for export_f_theta."
            end if
            stop
         end if

         ! Force all requested theta values to the interval [0,2pi):
         do j=1,N_export_f_theta
            export_f_theta(j) = modulo(export_f_theta(j),2*pi)
         end do

         allocate(map_theta_to_export_f_theta(N_export_f_theta, Ntheta))
         map_theta_to_export_f_theta = zero

         do j=1,N_export_f_theta
            ! In the following calculation, I assume theta is a uniform mesh beginning with 0.
            index1 = floor(export_f_theta(j)*Ntheta/(2*pi))+1
            if (index1 < 1) then
               print *,"Something went wrong!  index1=",index1,", export_f_theta(j)=",export_f_theta(j)
               stop
            elseif (index1 == Ntheta+1) then
               ! This could happen if user requests exactly 2pi
               index1 = Ntheta
               index2 = 1
            elseif (index1 > Ntheta+1) then
               print *,"Something went wrong!  index1=",index1,", export_f_theta(j)=",export_f_theta(j)
               stop
            elseif (index1 == Ntheta) then
               ! Fine, but we are interpolating close to 2pi, so must interpolate between the 1st and last points.
               index2 = 1
            else
               ! All is good.
               index2 = index1 + 1
            end if
            weight1 = index1 - export_f_theta(j)*Ntheta/(2*pi)
            weight2 = one - weight1

            map_theta_to_export_f_theta(j, index1) = weight1
            map_theta_to_export_f_theta(j, index2) = weight2
         end do

      case (2)
         ! 2: Do not interpolate. Use the values of the internal theta grid that are closest to the values requested in export_f_theta.
         ! In this case, we will over-write export_f_theta and N_export_f_theta.

         if (N_export_f_theta < 1) then
            if (masterProc) then
               print *,"Error! When export_f_theta_option=2, you must specify at least 1 value for export_f_theta."
            end if
            stop
         end if

         ! Force all requested theta values to the interval [0,2pi):
         do j=1,N_export_f_theta
            export_f_theta(j) = modulo(export_f_theta(j),2*pi)
         end do

         allocate(includeThisTheta(Ntheta))
         includeThisTheta = .false.

         do i=1,N_export_f_theta
            ! Determine the closest point in the theta grid to this requested theta
            indexOfLeastError = 1
            ! When computing the error, also consider the possibility that the points could straddle 2pi in either direction:
            leastError = min(  (export_f_theta(i)-theta(1)) ** 2, &
                 (export_f_theta(i)-theta(1)-2*pi) ** 2, &
                 (export_f_theta(i)-theta(1)+2*pi) ** 2)

            do j=2,Ntheta
               error = min(  (export_f_theta(i)-theta(j)) ** 2, &
                    (export_f_theta(i)-theta(j)-2*pi) ** 2, &
                    (export_f_theta(i)-theta(j)+2*pi) ** 2)

               if (error < leastError) then
                  indexOfLeastError = j
                  leastError = error
               end if
            end do
            includeThisTheta(indexOfLeastError) = .true.
         end do

         ! Over-write export_f_theta and N_export_f_theta:
         index = 1
         do j=1,Ntheta
            if (includeThisTheta(j)) then
               export_f_theta(index) = theta(j)
               index = index + 1
            end if
         end do
         N_export_f_theta = index-1

         allocate(map_theta_to_export_f_theta(N_export_f_theta,Ntheta))
         map_theta_to_export_f_theta = zero
         index = 1
         do j=1,Ntheta
            if (includeThisTheta(j)) then
               map_theta_to_export_f_theta(index, j) = one
               index = index + 1
            end if
         end do

         deallocate(includeThisTheta)

      case default
         if (masterProc) then
            print *,"Error! Invalid setting for export_f_theta_option"
         end if
         stop
      end select



      ! --------------------------------------------------------
      ! Handle zeta coordinate
      ! --------------------------------------------------------

      if (Nzeta==1) then
         ! Ignore export_f_zeta_option and export_f_zeta.
         N_export_f_zeta = 1
         export_f_zeta = 0
         allocate(map_zeta_to_export_f_zeta(1,1))
         map_zeta_to_export_f_zeta = one

      else
         select case (export_f_zeta_option)
         case (0)
            ! 0: Report the distribution function on the original zeta grid used internally by SFINCS (with Nzeta elements).
            N_export_f_zeta = Nzeta

            if (Nzeta>max_N_export_f) then
               ! Very unlikely, but handle it just in case, to prevent seg fault in a moment when setting export_f_zeta.
               if (masterProc) then
                  print *,"Error! When export_f_zeta_option=0, max_N_export_f must be >= Nzeta."
               end if
               stop
            end if
            export_f_zeta(1:Nzeta) = zeta

            ! Build the Nzeta * Nzeta identity matrix:
            allocate(map_zeta_to_export_f_zeta(Nzeta,Nzeta))
            map_zeta_to_export_f_zeta = zero
            do j=1,Nzeta
               map_zeta_to_export_f_zeta(j,j)=1
            end do

         case (1)
            ! 1: Interpolate to a different grid, specified by export_f_zeta. 

            ! Validation:
            if (N_export_f_zeta < 1) then
               if (masterProc) then
                  print *,"Error! When export_f_zeta_option=1, you must specify at least 1 value for export_f_zeta."
               end if
               stop
            end if

            ! Force all requested zeta values to the interval [0,2pi/Nperiods):
            do j=1,N_export_f_zeta
               export_f_zeta(j) = modulo(export_f_zeta(j),2*pi/Nperiods)
            end do

            allocate(map_zeta_to_export_f_zeta(N_export_f_zeta, Nzeta))
            map_zeta_to_export_f_zeta = zero

            do j=1,N_export_f_zeta
               ! In the following calculation, I assume zeta is a uniform mesh beginning with 0.
               index1 = floor(export_f_zeta(j)*Nzeta/(2*pi/Nperiods))+1 ! Note the factor of Nperiods here is different from the theta variable!
               if (index1 < 1) then
                  print *,"Something went wrong!  index1=",index1,", export_f_zeta(j)=",export_f_zeta(j)
                  stop
               elseif (index1 == Nzeta+1) then
                  ! This could happen if user requests exactly 2pi/Nperiods
                  index1 = Nzeta
                  index2 = 1
               elseif (index1 > Nzeta+1) then
                  print *,"Something went wrong!  index1=",index1,", export_f_zeta(j)=",export_f_zeta(j)
                  stop
               elseif (index1 == Nzeta) then
                  ! Fine, but we are interpolating close to 2pi/Nperiods, so must interpolate between the 1st and last points.
                  index2 = 1
               else
                  ! All is good.
                  index2 = index1 + 1
               end if
               weight1 = index1 - export_f_zeta(j)*Nzeta/(2*pi/Nperiods) ! Note the factor of Nperiods here is different from the theta variable!
               weight2 = one - weight1

               map_zeta_to_export_f_zeta(j, index1) = weight1
               map_zeta_to_export_f_zeta(j, index2) = weight2
            end do

         case (2)
            ! 2: Do not interpolate. Use the values of the internal zeta grid that are closest to the values requested in export_f_zeta.
            ! In this case, we will over-write export_f_zeta and N_export_f_zeta.

            if (N_export_f_zeta < 1) then
               if (masterProc) then
                  print *,"Error! When export_f_zeta_option=2, you must specify at least 1 value for export_f_zeta."
               end if
               stop
            end if

            ! Force all requested zeta values to the interval [0,2pi/Nperiods):
            do j=1,N_export_f_zeta
               export_f_zeta(j) = modulo(export_f_zeta(j),2*pi/Nperiods)
            end do

            allocate(includeThisZeta(Nzeta))
            includeThisZeta = .false.

            do i=1,N_export_f_zeta
               ! Determine the closest point in the zeta grid to this requested zeta
               indexOfLeastError = 1
               ! When computing the error, also consider the possibility that the points could straddle 2pi/Nperiods in either direction:
               leastError = min(  (export_f_zeta(i)-zeta(1)) ** 2, &
                    (export_f_zeta(i)-zeta(1)-2*pi/Nperiods) ** 2, &
                    (export_f_zeta(i)-zeta(1)+2*pi/Nperiods) ** 2)

               do j=2,Nzeta
                  error = min(  (export_f_zeta(i)-zeta(j)) ** 2, &
                       (export_f_zeta(i)-zeta(j)-2*pi/Nperiods) ** 2, &
                       (export_f_zeta(i)-zeta(j)+2*pi/Nperiods) ** 2)

                  if (error < leastError) then
                     indexOfLeastError = j
                     leastError = error
                  end if
               end do
               includeThisZeta(indexOfLeastError) = .true.
            end do

            ! Over-write export_f_zeta and N_export_f_zeta:
            index = 1
            do j=1,Nzeta
               if (includeThisZeta(j)) then
                  export_f_zeta(index) = zeta(j)
                  index = index + 1
               end if
            end do
            N_export_f_zeta = index-1

            allocate(map_zeta_to_export_f_zeta(N_export_f_zeta,Nzeta))
            map_zeta_to_export_f_zeta = zero
            index = 1
            do j=1,Nzeta
               if (includeThisZeta(j)) then
                  map_zeta_to_export_f_zeta(index, j) = one
                  index = index + 1
               end if
            end do

            deallocate(includeThisZeta)

         case default
            if (masterProc) then
               print *,"Error! Invalid setting for export_f_zeta_option"
            end if
            stop
         end select
      end if




      ! --------------------------------------------------------
      ! Handle x coordinate
      ! --------------------------------------------------------


      select case (export_f_x_option)
      case (0)
         ! 0: Report the distribution function on the original x grid used internally by SFINCS (with Nx elements).
         N_export_f_x = Nx

         if (Nx>max_N_export_f) then
            ! Very unlikely, but handle it just in case, to prevent seg fault in a moment when setting export_f_x.
            if (masterProc) then
               print *,"Error! When export_f_x_option=0, max_N_export_f must be >= Nx."
            end if
            stop
         end if
         export_f_x(1:Nx) = x

         ! Build the Nx * Nx identity matrix:
         allocate(map_x_to_export_f_x(Nx,Nx))
         map_x_to_export_f_x = zero
         do j=1,Nx
            map_x_to_export_f_x(j,j)=1
         end do

      case (1)
         ! 1: Interpolate to a different grid, specified by export_f_x. 

         ! Validation:
         if (N_export_f_x < 1) then
            if (masterProc) then
               print *,"Error! When export_f_x_option=1, you must specify at least 1 value for export_f_x."
            end if
            stop
         end if
         do i=1,N_export_f_x
            if (export_f_x(i)<0) then
               if (masterProc) then
                  print *,"Error! Requested values of x (normalized speed) in export_f_x must be >= 0."
                  print *,"Here are the values of export_f_x you requested:"
                  print *,export_f_x(1:N_export_f_x)
               end if
               stop
            end if
         end do

         allocate(map_x_to_export_f_x(N_export_f_x, Nx))
         select case (xGridScheme)
         case (1,2)
            call polynomialInterpolationMatrix(Nx, N_export_f_x, x, export_f_x, &
                 expx2, exp(-export_f_x*export_f_x), map_x_to_export_f_x)
         case (3)
            allocate(extrapMatrix(N_export_f_x, Nx+1))
            allocate(regridPolynomialToUniform_plus1(N_export_f_x, Nx+1))
            allocate(x_plus1(Nx+1))
            x_plus1(1:Nx) = x
            x_plus1(Nx+1) = x(Nx)*2-x(Nx-1)
            call interpolationMatrix(Nx+1, N_export_f_x, x_plus1, export_f_x, regridPolynomialToUniform_plus1, extrapMatrix)
            map_x_to_export_f_x = regridPolynomialToUniform_plus1(:,1:Nx)
            deallocate(extrapMatrix)
            deallocate(regridPolynomialToUniform_plus1)
            deallocate(x_plus1)
         case default
            if (masterProc) then
               print *,"Error! Invalid xGridScheme."
            end if
            stop
         end select

      case (2)
         ! 2: Do not interpolate. Use the values of the internal x grid that are closest to the values requested in export_f_x.
         ! In this case, we will over-write export_f_x and N_export_f_x.

         if (N_export_f_x < 1) then
            if (masterProc) then
               print *,"Error! When export_f_x_option=2, you must specify at least 1 value for export_f_x."
            end if
            stop
         end if

         allocate(includeThisX(Nx))
         includeThisX = .false.

         do i=1,N_export_f_x
            ! Determine the closest point in the spectral x grid to this requested x
            indexOfLeastError = 1
            leastError = (export_f_x(i)-x(1)) ** 2
            do j=2,Nx
               error = (export_f_x(i)-x(j)) ** 2
               if (error < leastError) then
                  indexOfLeastError = j
                  leastError = error
               end if
            end do
            includeThisX(indexOfLeastError) = .true.
         end do


         ! Over-write export_f_x and N_export_f_x:
         index = 1
         do j=1,Nx
            if (includeThisX(j)) then
               export_f_x(index) = x(j)
               index = index + 1
            end if
         end do
         N_export_f_x = index-1

         allocate(map_x_to_export_f_x(N_export_f_x,Nx))
         map_x_to_export_f_x = zero
         index = 1
         do j=1,Nx
            if (includeThisX(j)) then
               map_x_to_export_f_x(index, j) = one
               index = index + 1
            end if
         end do

         deallocate(includeThisX)

      case default
         if (masterProc) then
            print *,"Error! Invalid setting for export_f_x_option"
         end if
         stop
      end select


      ! --------------------------------------------------------
      ! Handle xi coordinate
      ! --------------------------------------------------------


      select case (export_f_xi_option)
      case (0)
         ! 0: Report the distribution function as amplitudes of Nxi Legendre polynomials, as used internally by SFINCS for solving the kinetic equation.
         N_export_f_xi = Nxi
         ! export_f_xi should not appear in the output file, so we do not need to set it.

         ! Construct the Nxi * Nxi identity matrix:
         allocate(map_xi_to_export_f_xi(N_export_f_xi, Nxi))
         map_xi_to_export_f_xi = zero
         do j=1,Nxi
            map_xi_to_export_f_xi(j,j) = one
         end do

      case (1)
         ! 1: Report the distribution function on the values of xi specified by export_f_xi

         ! Validation:
         do i=1,N_export_f_xi
            if (export_f_xi(i)<-1 .or. export_f_xi(i) > 1) then
               if (masterProc) then
                  print *,"Error! Requested values of xi in export_f_xi must lie in the interval [-1,1]."
                  print *,"Here are the values of export_f_xi you requested:"
                  print *,export_f_xi(1:N_export_f_xi)
               end if
               stop
            end if
         end do

         ! Form the matrix which transforms from amplitudes of Legendre polynomials to values on the requested grid.
         ! I.e., evaluate the Legendre polynomials on this grid.
         allocate(map_xi_to_export_f_xi(N_export_f_xi, Nxi))
         L = 0
         map_xi_to_export_f_xi(:,L+1) = one ! P_0(xi) = 1
         if (Nxi>1) then
            ! Hard to believe Nxi would ever be 1, but just in case...
            L = 1
            map_xi_to_export_f_xi(:,L+1) = export_f_xi(1:N_export_f_xi) ! P_1(xi) = xi
         end if
         do L=2,Nxi-1
            ! Recursion relation for Legendre polynomials:
            map_xi_to_export_f_xi(:,L+1) = ((2*L-1)*export_f_xi(1:N_export_f_xi) * map_xi_to_export_f_xi(:,L) &
                 - (L-1)*map_xi_to_export_f_xi(:,L-1))/(one*L)
         end do

      case default
         if (masterProc) then
            print *,"Error! Invalid setting for export_f_xi_option"
         end if
         stop
      end select




      if (masterProc) then
         if (export_full_f) then
            allocate(full_f(Nspecies, N_export_f_theta, N_export_f_zeta, N_export_f_xi, N_export_f_x))
         end if

         if (export_delta_f) then
            allocate(delta_f(Nspecies, N_export_f_theta, N_export_f_zeta, N_export_f_xi, N_export_f_x))
         end if
      end if

    end subroutine setup_grids_for_export_f

end module export_f
