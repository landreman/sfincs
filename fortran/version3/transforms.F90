module transforms

#include "PETScVersions.F90"

  implicit none

  public :: forward_b, inverse_b
  
contains

  subroutine forward_b(real_space_matrix, fourier_space_vec, xn_vec, xm_vec, nharmonics, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)
    ! THIS TRANSFORM IS ONLY CORRECT FOR CASE 11
    use globalVariables, only: Ntheta, Nzeta, Nperiods, theta, zeta
    implicit none
    integer :: nharmonics
    integer :: i, j, k, itheta, izeta, ntheta_fft, nzeta_fft
    logical :: include_mn, sine_term
    integer, dimension(:) :: xn_vec, xm_vec
    PetscScalar :: dzeta_fft, dtheta_fft
    PetscScalar, dimension(:) :: theta_fft
    PetscScalar, dimension(:) :: zeta_fft
    PetscScalar, dimension(:), intent(out) :: fourier_space_vec
    PetscScalar, dimension(:,:), intent(in) :: real_space_matrix
    ! The nharmonics input does not include the (0,0) mode, so we start the loop at 2 and add in the (0,0) contribution at the end
!    print *,"theta = ",theta
!    print *,"zeta = ",zeta
    do i = 2, nharmonics+1
       include_mn = .false.
       if ((abs(xn_vec(i-1))<=int(Nzeta/2.0)).and.(xm_vec(i-1)<=int(Nzeta/2.0))) then
          include_mn = .true.
       end if
       if (Nzeta==1) then
          include_mn = .true.
       end if
       if (include_mn) then
          do itheta = 1,ntheta_fft
             do izeta = 1, nzeta_fft
                fourier_space_vec(i) = fourier_space_vec(i) + 2*real_space_matrix(itheta,izeta) &
                     *cos(xm_vec(i-1)*theta_fft(itheta) - Nperiods*xn_vec(i-1)*zeta_fft(izeta))
                if (sine_term) then
                   fourier_space_vec(i+nharmonics) = fourier_space_vec(i+nharmonics) + 2*real_space_matrix(itheta,izeta)*sin(xm_vec(i-1)*theta_fft(itheta) - Nperiods*xn_vec(i-1)*zeta_fft(izeta))             
                end if

             end do
          end do
       end if
!       print *,"theta = ",theta(itheta)
!       if (sine_term) then
!          print *,"n = ",xn_vec(i-1)
!          print *,"m = ",xm_vec(i-1)
!          print *,"fourier_space_vec(i) = ",fourier_space_vec(i)
!       end if
    end do
    
    fourier_space_vec(1) = sum(real_space_matrix)
!    if (sine_term) then
!       print *,"nothing"
!    else
!       print *,"fourier_space_vec(1) = ",fourier_space_vec(1)
!    end if
    fourier_space_vec = fourier_space_vec / (ntheta_fft*nzeta_fft)
!    print *,"fourier_space_vec = ",fourier_space_vec

  end subroutine forward_b

  ! -----------------------------------------------------------------------------------------

  subroutine inverse_b(fourier_space_vec, real_space_matrix, xn_vec, xm_vec, nharmonics, include_zero, sine_term, ntheta_fft, nzeta_fft, theta_fft, zeta_fft, dzeta_fft, dtheta_fft)
    ! THIS TRANSFORM IS ONLY CORRECT FOR CASE 11
    use globalVariables, only: Ntheta, Nzeta, Nperiods, theta, zeta
    implicit none
    integer :: nharmonics
    integer :: i, itheta, izeta, ntheta_fft, nzeta_fft
    logical :: include_mn, include_zero, sine_term
    integer, dimension(:) :: xn_vec, xm_vec
    PetscScalar :: dzeta_fft, dtheta_fft
    PetscScalar, dimension(:) :: theta_fft
    PetscScalar, dimension(:) :: zeta_fft
    PetscScalar, dimension(:), intent(in) :: fourier_space_vec
    PetscScalar, dimension(:,:), intent(out) :: real_space_matrix

    if (include_zero) then ! True if the (0,0) component is included in this transform
       do i = 2, nharmonics+1
!          print *,"fourier_space_vec = ",fourier_space_vec(i+nharmonics)
!          print *,"xn_vec(i-1) = ",xn_vec(i-1)
!          print *,"xm_vec(i-1) = ",xm_vec(i-1)
          include_mn = .false.
          if ((abs(xn_vec(i-1))<=int(Nzeta/2.0)).and.(xm_vec(i-1)<=int(Nzeta/2.0))) then
             include_mn = .true.
          end if
          if (Nzeta==1) then
             include_mn = .true.
          end if
          if (include_mn) then
             do itheta = 1,ntheta_fft
                do izeta = 1,nzeta_fft
                   real_space_matrix(itheta,izeta) = real_space_matrix(itheta,izeta) + fourier_space_vec(i)*cos(xm_vec(i-1)*theta_fft(itheta) - Nperiods*xn_vec(i-1)*zeta_fft(izeta))
                   if (sine_term) then
                      real_space_matrix(itheta,izeta) = real_space_matrix(itheta,izeta) + fourier_space_vec(i+nharmonics)*sin(xm_vec(i-1)*theta_fft(itheta) - Nperiods*xn_vec(i-1)*zeta_fft(izeta))
                   end if
                end do
             end do
          end if
       end do
       real_space_matrix = real_space_matrix + fourier_space_vec(1)
    else
       do i = 2, nharmonics+1
          include_mn = .false.
          if ((abs(xn_vec(i-1))<=int(Nzeta/2.0)).and.(xm_vec(i-1)<=int(Nzeta/2.0))) then
             include_mn = .true.
          end if
          if (Nzeta==1) then
             include_mn = .true.
          end if
          if (include_mn) then
             do itheta = 1,ntheta_fft
                do izeta = 1,nzeta_fft
                   real_space_matrix(itheta,izeta) = real_space_matrix(itheta,izeta) + fourier_space_vec(i-1)*cos(xm_vec(i-1)*theta_fft(itheta) - Nperiods*xn_vec(i-1)*zeta_fft(izeta))
                   if (sine_term) then
                      real_space_matrix(itheta,izeta) = real_space_matrix(itheta,izeta) + fourier_space_vec(i+nharmonics-1)*sin(xm_vec(i-1)*theta_fft(itheta) - Nperiods*xn_vec(i-1)*zeta_fft(izeta))
                   end if
                end do
             end do
          end if
       end do
    end if

  end subroutine inverse_b  

end module transforms
