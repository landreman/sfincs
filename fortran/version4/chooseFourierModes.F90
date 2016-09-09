subroutine chooseFourierModes()

  use rankMod
  use petscsys
  use FourierTransformMod
  use globalVariables, only: NFourier, NFourier2, xm, xn, mmax, nmax, NPeriods, prec, BHat, masterProc, &
       Ntheta, Nzeta, predictedAmplitudes, FourierOption, FourierFactor

  implicit none

  integer :: NFourier_original
  integer :: m, n, index, imn
  real(prec), dimension(:,:), allocatable :: patternWithExpectedSpectrum1
  real(prec), dimension(:,:), allocatable :: patternWithExpectedSpectrum2
  real(prec), dimension(:,:), allocatable :: patternWithExpectedSpectrum3
  real(prec), dimension(:,:), allocatable :: BScaledToNearly1
  real(prec), dimension(:), allocatable :: FourierVector1
  real(prec), dimension(:), allocatable :: FourierVector2
  real(prec), dimension(:), allocatable :: FourierVector3
  real(prec), dimension(:), allocatable :: amplitudes
  integer, dimension(:), allocatable :: xn_sorted, xm_sorted, permutation

  call initFourierTrig()

  allocate(predictedAmplitudes(mmax+1,2*nmax+1))
  predictedAmplitudes = 0

  if (nmax==0) then
     ! Axisymmetry
     
     allocate(xm(NFourier))
     allocate(xn(NFourier))

     xn=0
     do m=0,NFourier-1
        xm(m+1)=m
     end do

     if (masterProc) then
        print *,"Axisymmetry is being enforced."
        print *,"xm:",xm
        print *,"xn:",xn
     end if

  else

     select case (FourierOption)
     case (0)
        ! Include all (m,n) such that m<=mmax and |n|<=nmax.
        ! NFourier should already have been set to mmax*(2*nmax+1)+nmax+1 in validateInput.
        ! But let's do it again here to be safe:
        NFourier = mmax*(2*nmax+1)+nmax+1
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

     case (1,2)
        ! Pick certain (m,n) pairs to include
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
        
        allocate(amplitudes(NFourier))
        amplitudes=0

        select case (FourierOption)
        case (1)
           amplitudes = - (xm*xm + ((real(xn)/NPeriods/FourierFactor)**2))
           if (masterProc) then
              print *,"Nonaxisymmetric geometry, keeping all (m,n) modes up to mmax and nmax."
           end if

        case (2)
           ! Guess a function that we expect will have a similar Fourier spectrum to the distribution function:
           allocate(patternWithExpectedSpectrum1(Ntheta,Nzeta))
           allocate(patternWithExpectedSpectrum2(Ntheta,Nzeta))
           allocate(patternWithExpectedSpectrum3(Ntheta,Nzeta))
           allocate(BScaledToNearly1(Ntheta,Nzeta))
           BScaledToNearly1 = BHat / (sum(BHat) / (Ntheta*Nzeta))
           if (masterProc) then
              print *,"max & min of BScaledToNearly1:",maxval(BScaledToNearly1),minval(BScaledToNearly1)
           end if
           patternWithExpectedSpectrum1 = (BScaledToNearly1 ** 6)
           patternWithExpectedSpectrum2 = 1/(BScaledToNearly1 ** (1.3))
           patternWithExpectedSpectrum3 = 1/(BScaledToNearly1 ** 2)
           
           allocate(FourierVector1(NFourier2))
           allocate(FourierVector2(NFourier2))
           allocate(FourierVector3(NFourier2))
           call FourierTransform(patternWithExpectedSpectrum1, FourierVector1)
           call FourierTransform(patternWithExpectedSpectrum2, FourierVector2)
           call FourierTransform(patternWithExpectedSpectrum3, FourierVector3)
           deallocate(BScaledToNearly1)
           deallocate(patternWithExpectedSpectrum1)
           deallocate(patternWithExpectedSpectrum2)
           deallocate(patternWithExpectedSpectrum3)
           
           do imn=2,NFourier
              amplitudes(imn) = max( &
                   sqrt(FourierVector1(imn)**2 + FourierVector1(imn+NFourier-1)**2), &
                   sqrt(FourierVector2(imn)**2 + FourierVector2(imn+NFourier-1)**2), &
                   sqrt(FourierVector3(imn)**2 + FourierVector3(imn+NFourier-1)**2) )
           end do
           deallocate(FourierVector1)
           deallocate(FourierVector2)
           deallocate(FourierVector3)
        
           ! Add an offset proportional to sqrt(m^2 + n^2) so once the 'real' amplitude gets down to machine precision,
           ! additional modes get added first at low |m| and |n| rather than randomly.
           amplitudes = amplitudes - (1.0d-13)*sqrt(real(xm*xm + ((real(xn)/NPeriods)**2)))
        
        end select

        ! For plotting, it is convenient if the amplitudes are all >0:
        amplitudes = amplitudes - minval(amplitudes) + (1d-16)
        ! Make sure the m=n=0 mode comes first:
        amplitudes(1) = maxval(amplitudes) + 1

        ! Copy 1D array to the 2D array that will be saved in the HDF5 file:
        do imn=1,NFourier
           predictedAmplitudes(xm(imn)+1, xn(imn)/NPeriods+nmax+1) = amplitudes(imn)
        end do
        
        ! Sort amplitudes
        allocate(permutation(NFourier))
        call rank(amplitudes,permutation)
        if (NFourier-NFourier_original+1<1) then
           if (masterProc) then
              print *,"ERROR! NFourier exceeds mmax*(nmax*2+1)+mmax+1"
           end if
           stop
        end if
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
           print *,"Nonaxisymmetric geometry, keeping selected (m,n) Fourier modes."
           print *,"xm:",xm
           print *,"xn:",xn
        end if
        
        if (masterProc .and. maxval(xm)==mmax) then
           print *,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
           print *,"WARNING: A Fourier mode with m=mmax is expected to have significant amplitude."
           print *,"You should increase mmax to ensure the spectrum is not being cut off."
           print *,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
        if (masterProc .and. maxval(abs(xn))==nmax .and. nmax>0) then
           ! If nmax=0, don't worry, we are just doing an axisymmetric calculation.
           print *,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
           print *,"WARNING: A Fourier mode with n=nmax is expected to have significant amplitude."
           print *,"You should increase nmax to ensure the spectrum is not being cut off."
           print *,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
        
        deallocate(amplitudes)
        deallocate(xm_sorted,xn_sorted)

     case default
        if (masterProc) print *,"Error! Invalid setting for FourierOption:",FourierOption
        stop

     end select
  end if

end subroutine chooseFourierModes


