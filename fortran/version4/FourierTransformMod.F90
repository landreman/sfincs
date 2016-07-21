! We put this subroutine in a module because otherwise there is a segmentation fault when passing allocatable arrays as parameters.
module FourierTransformMod

  use globalVariables, only: dp

  implicit none

  private
  real(dp), dimension(:,:), allocatable :: sinmtheta, cosmtheta, sinnzeta, cosnzeta

  public :: FourierTransform, initFourierTrig, inverseFourierTransform

contains

  subroutine FourierTransform(realSpaceMatrix, FourierVector)

    use globalVariables, only: Ntheta, Nzeta, theta, zeta, xm, xn, NFourier, NPeriods
    
    implicit none
    
    real(dp), dimension(:,:), intent(in) :: realSpaceMatrix
    real(dp), dimension(:), intent(out) :: FourierVector
    
    real(dp) :: angle
    integer :: itheta, izeta, factor
    integer :: imn, m, n
    integer :: tic, toc, countrate

    ! ***************************************************************
    ! Validate input
    ! ***************************************************************
   
    if (xm(1) .ne. 0) stop "xm(1) should be 0."
    if (xn(1) .ne. 0) stop "xn(1) should be 0."
    if (size(FourierVector) .ne. NFourier*2-1) stop "Size of FourierVector is not correct"
    if (size(realSpaceMatrix, 1) .ne. Ntheta) stop "Size of first dimension of realSpaceMatrix is not correct"
    if (size(realSpaceMatrix, 2) .ne. Nzeta)  stop "Size of second dimension of realSpaceMatrix is not correct"
    
    ! ***************************************************************
    ! Main calculation
    ! ***************************************************************
    
    FourierVector = 0
    call system_clock(tic,countrate)

    ! Rather than use an FFT, let's just use a 'slow' but transparent algorithm. Speed is not crucial for this subroutine.
    do imn = 2,NFourier ! Start at 2 rather than 1. We will handle the imn=1 (m=n=0) case later.
       m = xm(imn)
       n = xn(imn)/NPeriods
       do itheta = 1,Ntheta
          do izeta = 1,Nzeta
             !angle = m*theta(itheta) - n*NPeriods*zeta(izeta)
             
             !! cos(m theta - n zeta) = cos(m theta) cos(n zeta) + sin(m theta) sin(n zeta)
             !FourierVector(imn) = FourierVector(imn) + 2*cos(angle)*realSpaceMatrix(itheta,izeta)
             FourierVector(imn) = FourierVector(imn) + 2*realSpaceMatrix(itheta,izeta) &
                  *(cosmtheta(itheta,m)*cosnzeta(izeta,n) + sinmtheta(itheta,m)*sinnzeta(izeta,n))

             !! sin(m theta - n zeta) = sin(m theta) cos(n zeta) - cos(m theta) sin(n zeta)
             !FourierVector(imn+NFourier-1) = FourierVector(imn+NFourier-1) + 2*sin(angle)*realSpaceMatrix(itheta,izeta)
             FourierVector(imn+NFourier-1) = FourierVector(imn+NFourier-1) + 2*realSpaceMatrix(itheta,izeta) &
                  *(sinmtheta(itheta,m)*cosnzeta(izeta,n) - cosmtheta(itheta,m)*sinnzeta(izeta,n))
          end do
       end do
    end do
    ! Handle the m=n=0 case:
    FourierVector(1) = sum(realSpaceMatrix)
    
    FourierVector = FourierVector / (Ntheta*Nzeta)

    call system_clock(toc)
    print *,"  Fourier transform:",real(toc-tic)/countrate,"sec."

  end subroutine FourierTransform

  ! -------------------------------------------------

  subroutine inverseFourierTransform(FourierVector, realSpaceMatrix)

    use globalVariables, only: Ntheta, Nzeta, theta, zeta, xm, xn, NFourier, NFourier2
    
    implicit none
    
    real(dp), dimension(:), intent(in) :: FourierVector
    real(dp), dimension(:,:), intent(out) :: realSpaceMatrix
    
    integer :: itheta, izeta
    integer :: imn, m, n
    integer :: tic, toc, countrate

    ! ***************************************************************
    ! Validate input
    ! ***************************************************************
   
    if (xm(1) .ne. 0) stop "xm(1) should be 0."
    if (xn(1) .ne. 0) stop "xn(1) should be 0."
    if (size(FourierVector) .ne. NFourier2) stop "Size of FourierVector is not correct"
    if (size(realSpaceMatrix, 1) .ne. Ntheta) stop "Size of first dimension of realSpaceMatrix is not correct"
    if (size(realSpaceMatrix, 2) .ne. Nzeta)  stop "Size of second dimension of realSpaceMatrix is not correct"
    
    ! ***************************************************************
    ! Main calculation
    ! ***************************************************************
    
    call system_clock(tic,countrate)

    ! Rather than use an FFT, let's just use a 'slow' but transparent algorithm. Speed is not crucial for this subroutine.
    realSpaceMatrix = 0
    do imn = 2,NFourier ! Start at 2 rather than 1. We will handle the imn=1 (m=n=0) case later.
       m = xm(imn)
       n = xn(imn)
       do izeta = 1,Nzeta
          do itheta = 1,Ntheta
             !angle = m*theta(itheta) - n*zeta(izeta)
             
             !! cos(m theta - n zeta) = cos(m theta) cos(n zeta) + sin(m theta) sin(n zeta)
             realSpaceMatrix(itheta,izeta) = realSpaceMatrix(itheta,izeta) + FourierVector(imn) &
                  *(cosmtheta(itheta,m)*cosnzeta(izeta,n) + sinmtheta(itheta,m)*sinnzeta(izeta,n))

             !! sin(m theta - n zeta) = sin(m theta) cos(n zeta) - cos(m theta) sin(n zeta)
             realSpaceMatrix(itheta,izeta) = realSpaceMatrix(itheta,izeta) + FourierVector(imn+NFourier-1) &
                  *(sinmtheta(itheta,m)*cosnzeta(izeta,n) - cosmtheta(itheta,m)*sinnzeta(izeta,n))
          end do
       end do
    end do
    ! Handle the m=n=0 case:
    realSpaceMatrix = realSpaceMatrix + FourierVector(1)
    
    call system_clock(toc)
    print *,"  inverse Fourier transform:",real(toc-tic)/countrate,"sec."

  end subroutine inverseFourierTransform

  ! -------------------------------------------------

  subroutine initFourierTrig

    use globalVariables, only: Ntheta, Nzeta, mmax, nmax, theta, zeta, NPeriods

    implicit none

    integer :: m, n, itheta, izeta
    integer :: tic, toc, countrate

!!$    if (allocated(sinmtheta)) deallocate(sinmtheta)
!!$    if (allocated(cosmtheta)) deallocate(cosmtheta)
!!$    if (allocated(sinnzeta))  deallocate(sinnzeta)
!!$    if (allocated(cosnzeta))  deallocate(cosnzeta)

    call system_clock(tic, countrate)
    allocate(sinmtheta(Ntheta,0:mmax))
    allocate(cosmtheta(Ntheta,0:mmax))
    allocate(sinnzeta (Nzeta, -nmax:nmax))
    allocate(cosnzeta (Nzeta, -nmax:nmax))

    do itheta=1,Ntheta
       do m=0,mmax
          sinmtheta(itheta,m) = sin(m*theta(itheta))
          cosmtheta(itheta,m) = cos(m*theta(itheta))
       end do
    end do
    do izeta=1,Nzeta
       do n=-nmax,nmax
          sinnzeta(izeta,n) = sin(n*zeta(izeta)*NPeriods)
          cosnzeta(izeta,n) = cos(n*zeta(izeta)*NPeriods)
       end do
    end do

    call system_clock(toc)
    print *,"  init Fourier trig:",real(toc-tic)/countrate,"sec."
  end subroutine initFourierTrig

end module FourierTransformMod
