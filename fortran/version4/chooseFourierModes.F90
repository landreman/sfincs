subroutine chooseFourierModes()

  use rankMod
  use petscsys
  use FourierTransformMod
  use globalVariables, only: NFourier, NFourier2, xm, xn, mmax, nmax, NPeriods, prec

  implicit none

  integer :: NFourier_original
  integer :: m, n, index, imn
  real(prec), dimension(:,:), allocatable :: patternWithExpectedSpectrum
  real(prec), dimension(:), allocatable :: FourierVector, amplitudes
  integer, dimension(:), allocatable :: xn_sorted, xm_sorted, permutation

  allocate(xn_sorted(NFourier))
  allocate(xm_sorted(NFourier))

  ! Temporarily set up (xm,xn) to include all possible modes.
  NFourier_original = NFourier
  NFourier = mmax*(2*nmax+1) + nmax + 1
  NFourier2 = NFourier*2-1
  allocate(xm(NFourier))
  allocate(xn(NFourier))

  xm=0
  xn=0
  ! Handle the m=0 modes
  do n=1,nmax
     xn(n+1)=n
  end do
  ! Handle the m>0 modes
  index = nmax+1
  do m=1,mmax
     do n=-nmax,nmax
        index = index + 1
        xn(index) = n
        xm(index) = m
     end do
  end do
  xn = xn*NPeriods
!!$  print *,"Initialized xm,xn to all non-constant modes."
!!$  print *,"xm=",xm
!!$  print *,"xn=",xn
!!$
!!$  print *,"B:"
!!$  do n=1,Ntheta
!!$     print *,B(n,:)
!!$  end do

  ! Guess a function that we expect will have a similar Fourier spectrum to the distribution function:
  allocate(patternWithExpectedSpectrum(Ntheta,Nzeta))
  patternWithExpectedSpectrum = (B ** 6) + 1/(B*B)

  allocate(FourierVector(NFourier2))
  call initFourierTrig()
  call FourierTransform(patternWithExpectedSpectrum, FourierVector)
  deallocate(patternWithExpectedSpectrum)

  allocate(amplitudes(NFourier))
  amplitudes=0
  do imn=2,NFourier
     amplitudes(imn) = sqrt(FourierVector(imn)**2 + FourierVector(imn+NFourier-1)**2)
  end do

  ! Add an offset proportional to sqrt(m^2 + n^2) so once the 'real' amplitude gets down to machine precision,
  ! additional modes get added first at low |m| and |n| rather than randomly.
  amplitudes = amplitudes - (1.0d-12)*sqrt(real(xm*xm + ((real(xn)/NPeriods)**2)))
  ! For plotting, it is convenient if the amplitudes are all >0:
  amplitudes = amplitudes - minval(amplitudes) + (1d-16)
  ! Make sure the m=n=0 mode comes first:
  amplitudes(1) = maxval(amplitudes) + 1

  ! Sort amplitudes
  allocate(permutation(NFourier))
  call rank(amplitudes,permutation)
  xm_sorted = xm(permutation(NFourier:NFourier-NFourier_original+1:-1))
  xn_sorted = xn(permutation(NFourier:NFourier-NFourier_original+1:-1))
!!$  print *,"amplitudes:",amplitudes
!!$  print *,"permutation:",permutation
!!$  print *,"xm_sorted:",xm_sorted
!!$  print *,"xn_sorted:",xn_sorted

  deallocate(xm,xn)
  NFourier = NFourier_original
  NFourier2 = NFourier*2-1
  allocate(xm(NFourier))
  allocate(xn(NFourier))
  xm = xm_sorted
  xn = xn_sorted
  if ((xm(1) .ne. 0) .or. (xn(1) .ne. 0)) then
     print *,"Error! initial (m,n) is not (0,0)."
     print *,"xm:",xm
     print *,"xn:",xn
     stop
  end if
  if (masterProc) then
     print *,"xm:",xm
     print *,"xn:",xn
  end if

  deallocate(FourierVector)
  deallocate(amplitudes)
  deallocate(xm_sorted,xn_sorted)

end subroutine chooseFourierModes


