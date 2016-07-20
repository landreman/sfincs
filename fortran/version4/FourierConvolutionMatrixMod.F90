! We put this subroutine in a module because otherwise there is a segmentation fault when passing allocatable arrays as parameters.
module FourierConvolutionMatrixMod

  implicit none

contains

  subroutine FourierConvolutionMatrix(vector, matrix)

    ! Suppose we have Fourier expansions for two functions, f(theta,zeta) and B(theta,zeta):
    ! B(theta,zeta) = \sum_{m1,n1} [ B^c_{m1,n1} cos(m1 theta - n1 zeta) + B^s_{m1,n1} sin(m1 theta - n1 zeta)
    ! and
    ! f(theta,zeta) = \sum_{m2,n2} [ f^c_{m2,n2} cos(m2 theta - n2 zeta) + f^2_{m2,n2} sin(m2 theta - n2 zeta).
    ! Then using trig identities it can be shown that the product B*f can be written as follows:
    ! B*f = (1/2) \sum_{m1,n1,m2,n2} *
    !       [ (B^c_{m1,n1} f^c_{m2,n2} - B^s_{m1,n1} f^s_{m2,n2}) cos((m1+m2) theta - (n1+n2) zeta)
    !        +(B^s_{m1,n1} f^c_{m2,n2} + B^c_{m1,n1} f^s_{m2,n2}) sin((m1+m2) theta - (n1+n2) zeta)
    !        +(B^c_{m1,n1} f^c_{m2,n2} + B^s_{m1,n1} f^s_{m2,n2}) cos((m1-m2) theta - (n1-n2) zeta)
    !        +(B^s_{m1,n1} f^c_{m2,n2} - B^c_{m1,n1} f^s_{m2,n2}) sin((m1-m2) theta - (n1-n2) zeta)].
    !
    ! If we think of [f^c; f^2] as a vector, then this formula above indicates that multiplying B*f amounts to
    ! operating on the f vector by a matrix. This subroutine constructs this matrix.
    

    ! vector should have been allocated with size  2*NFourier-1.
    ! matrix should have been allocated with size (2*NFourier-1, 2*NFourier-1).

    use globalVariables, only: NFourier, NFourier2, xm, xn, dp

    implicit none

    real(dp), intent(in), dimension(:) :: vector
    real(dp), intent(out), dimension(:,:) :: matrix

    integer :: j, m, n, imn1, imn2, sign, numMatches, imn_new
    real(dp), dimension(:), allocatable :: halfVector
    integer :: tic, toc, countrate

    ! ***************************************************************
    ! Validate input
    ! ***************************************************************

    if (size(xm) .ne. NFourier) stop "Size of xm is not correct"
    if (size(xn) .ne. NFourier) stop "Size of xn is not correct"
    if (xm(1) .ne. 0) stop "xm(1) should be 0."
    if (xn(1) .ne. 0) stop "xn(1) should be 0."
    if (size(vector) .ne. NFourier*2-1) stop "Size of vector is not correct"
    if (size(matrix, 1) .ne. NFourier*2-1) stop "Size of first dimension of matrix is not correct"
    if (size(matrix, 2) .ne. NFourier*2-1) stop "Size of second dimension of matrix is not correct"

    ! ***************************************************************
    ! Assemble matrix
    ! ***************************************************************

    call system_clock(tic,countrate)
    allocate(halfVector(NFourier2))
    halfVector = 0.5*vector

    matrix = 0.0
    ! imn1 indexes the vector that is provided as input to this subroutine.
    ! imn2 indexes the vector that this matrix will be applied to.
    
    do imn1 = 1,NFourier
       do imn2 = 1,NFourier
          
          ! Handle the terms for which m = m1 + m2, n = n1 + n2
          m = xm(imn1) + xm(imn2)
          n = xn(imn1) + xn(imn2)
          sign = 1
          if (m==0 .and. n<0) then
             n = -n
             sign = -1
          end if
          ! Find the index of (m,n) in the xm and xn arrays, if it is present.
          numMatches=0
          do j=1,NFourier
             if (m == xm(j) .and. n == xn(j)) then
                numMatches = numMatches + 1
                imn_new = j
             end if
          end do
          if (numMatches>1) then
             print *,"Found multiple instances of m=",m," n=",n," in the xm and xn arrays."
             print *,"xm=",xm
             print *,"xn=",xn
          elseif (numMatches==1) then
             ! Terms that result in cos(m theta - n zeta):
             matrix(imn_new, imn2)            = matrix(imn_new, imn2)            + halfVector(imn1)
             matrix(imn_new, imn2+NFourier-1) = matrix(imn_new, imn2+NFourier-1) - halfVector(imn1+NFourier-1)
             if ((n .ne. 0) .or. (m .ne. 0)) then
                ! Terms that result in sin(m theta - n zeta):
                matrix(imn_new+NFourier-1, imn2)            = matrix(imn_new+NFourier-1, imn2)            + sign*halfVector(imn1+NFourier-1)
                matrix(imn_new+NFourier-1, imn2+NFourier-1) = matrix(imn_new+NFourier-1, imn2+NFourier-1) + sign*halfVector(imn1)
             end if
          else
             ! The new (m,n) is outside the range we are keeping, so do nothing.
          end if
          
          ! ----------------------------------------------------------------------------------------
          
          ! Handle the terms for which m = m1 - m2, n = n1 - n2
          m = xm(imn1) - xm(imn2)
          n = xn(imn1) - xn(imn2)
          sign = 1
          if ((m==0 .and. n<0) .or. (m<0)) then
             m = -m
             n = -n
             sign = -1
          end if
          ! Find the index of (m,n) in the xm and xn arrays, if it is present.
          numMatches=0
          do j=1,NFourier
             if (m == xm(j) .and. n == xn(j)) then
                numMatches = numMatches + 1
                imn_new = j
             end if
          end do
          if (numMatches>1) then
             print *,"Found multiple instances of m=",m," n=",n," in the xm and xn arrays."
             print *,"xm=",xm
             print *,"xn=",xn
          elseif (numMatches==1) then
             ! Terms that result in cos(m theta - n zeta):
             matrix(imn_new, imn2)            = matrix(imn_new, imn2)            + halfVector(imn1)
             matrix(imn_new, imn2+NFourier-1) = matrix(imn_new, imn2+NFourier-1) + halfVector(imn1+NFourier-1)
             if ((n .ne. 0) .or. (m .ne. 0)) then
                ! Terms that result in sin(m theta - n zeta):
                matrix(imn_new+NFourier-1, imn2)            = matrix(imn_new+NFourier-1, imn2)            + sign*halfVector(imn1+NFourier-1)
                matrix(imn_new+NFourier-1, imn2+NFourier-1) = matrix(imn_new+NFourier-1, imn2+NFourier-1) - sign*halfVector(imn1)
             end if
          else
             ! The new (m,n) is outside the range we are keeping, so do nothing.
          end if
          
       end do
    end do
    
    deallocate(halfVector)
    call system_clock(toc)
    print *,"  convolution matrix: ",real(toc-tic)/countrate,"sec"

  end subroutine FourierConvolutionMatrix
end module FourierConvolutionMatrixMod
