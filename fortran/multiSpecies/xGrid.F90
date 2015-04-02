  ! This module computes the nodes and weights for Gaussian quadrature for any
  ! weight function on any interval.  The algorithm used is the Stieltjes
  ! procedure described in section 4.6.3 of Numerical Recipes, 3rd edition,
  ! together with the Jacobi matrix procedure of section 4.6.2.  See also
  ! www.nr.com/webnotes?3
  ! The module is presently configured to use the weight exp(-x^2) and to
  ! use the interval [0, Inf), which yields the grid and integration weights
  ! dicussed in the following paper:
  ! Landreman & Ernst, Journal of Computational Physics 243, 130 (2013).
  ! This grid is useful for velocity space.
  ! You can change the "weight" function in the module and the integration
  ! domain in the makeXGrid subroutine if you wish.
  !
  ! The type PetscScalar is used, so this module can be used in a PETSc
  ! application. However, no other PETSc functionality is used, so you can
  ! replace the type with e.g. real if you want to build a non-PETSc
  ! application.
  !
  ! The makeXGrid subroutine calls LAPACK, which is
  ! automatically included in any program that is linked to the PETSc 
  ! libraries.
  !
  ! The makeXGrid subroutine also calls QUADPACK for numerical
  ! integration.  You can download QUADPACK here:
  ! http://www.netlib.org/quadpack/
  !
  ! Matt Landreman
  ! Massachusetts Institute of Technology
  ! Plasma Science & Fusion Center
  ! November, 2012
  !
#define workspaceSize 500  
  
  module xGrid

    implicit none

#include <finclude/petscsysdef.h>

    PetscScalar, private, allocatable :: a(:), b(:)
    integer, private :: j
    private weight

  contains

    subroutine makeXGrid(N, abscissae, weights, includePointAtX0)
      ! This is the main function of the module.
      !
      ! Inputs:
      !  N = number of grid points to generate.
      !  includePointAtX0 = Should a point be included at x=0?
      !
      ! Outputs:
      !   abscissae = grid points
      !   weights = Gaussian integration weights
      !
      ! Both abscissae and weights should be allocated to have size N
      ! before makeXGrid is called.

      implicit none

      integer, intent(in) :: N
      PetscScalar, intent(out) :: abscissae(:), weights(:)
      logical, intent(in) :: includePointAtX0

      integer :: key = 6
      integer :: i, info, N_copy
      PetscScalar :: oldc = 1d+0, amountToAdd
      PetscScalar :: absTol = 1d-5, relTol = 1d-13
      PetscScalar, allocatable :: c(:), d(:), eigenvectors(:,:), LAPACKWorkspace(:)

      PetscScalar :: abserr,finiteBound,EPSABS,EPSREL,WORK,LAPACK_abstol
      integer :: ier, inf, last, limit, neval
      PetscScalar, dimension(workspaceSize) :: alist, blist, rlist, elist
      integer, dimension(workspaceSize) :: iord
      PetscScalar :: X0, lastPolynomialAtX0, penultimatePolynomialAtX0
      integer, allocatable, dimension(:) :: LAPACK_ISUPPZ, LAPACK_iwork
      integer :: workspace_size, iwork_size

      X0 = 0.0d+0  ! Special point to include among the abscissae, if requested.
      lastPolynomialAtX0 = 0.0d+0
      penultimatePolynomialAtX0 = 0.0d+0

      allocate(a(N))
      allocate(b(N+1))
      allocate(c(N))
      allocate(d(N))
      allocate(eigenvectors(N,N))
      ! For the workspace, LAPACK's DPTEQR requires 4*N, while SSTEGR requires 18*N:
      workspace_size = 18*N
      allocate(LAPACKWorkspace(workspace_size))
      iwork_size = 10*N
      allocate(LAPACK_iwork(iwork_size))

      allocate(LAPACK_ISUPPZ(2*N))

      finiteBound = 0.0E0

      inf = 1
      !     inf = 1 corresponds to  (bound,+infinity),
      !     inf = -1            to  (-infinity,bound),
      !     inf = 2             to (-infinity,+infinity).

      EPSABS = 0.0E0
      EPSREL = 1d-3
      LIMIT = workspaceSize

      finiteBound = 5d+0
      ! Call QUADPACK routines to evaluate integrals:
      ! qage = single-precision quadrature on finite domain
      ! dqage = double-precision quadrature on finite domain
      ! qagie = single-precision quadrature on (semi-)infinite domain
      ! dqagie = double-precision quadrature on (semi-)infinite domain
      do j=1,N
#if defined(PETSC_USE_REAL_SINGLE)
         call qagie(integrandWithoutX,finiteBound,inf,epsabs,epsrel,limit,c(j),abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)
#else
         call dqagie(integrandWithoutX,finiteBound,inf,epsabs,epsrel,limit,c(j),abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)
#endif
         if (ier /= 0) then
            print *,"Quadrature error A: ier = ",ier
            stop
         end if
#if defined(PETSC_USE_REAL_SINGLE)
         call qage(integrandWithoutX,0.,finiteBound,epsabs,epsrel,key,limit,amountToAdd,abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)
#else
         call dqage(integrandWithoutX,0d+0,finiteBound,epsabs,epsrel,key,limit,amountToAdd,abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)
#endif
         if (ier /= 0) then
            print *,"Quadrature error B: ier = ",ier
            stop
         end if
         c(j) = c(j) + amountToAdd

#if defined(PETSC_USE_REAL_SINGLE)
         call qagie(integrandWithX,finiteBound,inf,epsabs,epsrel,limit,d(j),abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)
#else
         call dqagie(integrandWithX,finiteBound,inf,epsabs,epsrel,limit,d(j),abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)
#endif
         if (ier /= 0) then
            print *,"Quadrature error C: ier = ",ier
            stop
         end if
#if defined(PETSC_USE_REAL_SINGLE)
         call qage(integrandWithX,0.,finiteBound,epsabs,epsrel,key,limit,amountToAdd,abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)
#else
         call dqage(integrandWithX,0d+0,finiteBound,epsabs,epsrel,key,limit,amountToAdd,abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)
#endif
         if (ier /= 0) then
            print *,"Quadrature error D: ier = ",ier
            stop
         end if
         d(j) = d(j) + amountToAdd

         b(j) = c(j) / oldc
         a(j) = d(j) / c(j)
         oldc = c(j)

         penultimatePolynomialAtX0 = lastPolynomialAtX0
         lastPolynomialAtX0 = evaluatePolynomial(X0)
      end do

      if (includePointAtX0) then
         a(N) = X0 - b(N) * penultimatePolynomialAtX0 / lastPolynomialAtX0
      end if

      b = sqrt(b)

      ! Call LAPACK to get eigenvectors & eigenvalues:
      ! We use the fact that the Jacobi matrix is symmetric and tridiagonal.
      if (includePointAtX0) then
         !LAPACK_abstol = DLAMCH('S')
         LAPACK_abstol = 1.0d-15
         call DSTEGR('V','A',N,a,b(2:N+1),0.0d+0,0.0d+0,0,0,LAPACK_abstol,N_copy,abscissae,eigenvectors,&
              N,LAPACK_ISUPPZ, LAPACKWorkspace, workspace_size, LAPACK_iwork, iwork_size, info)
         if (info /= 0) then
            print *,"Error in LAPACK's DSTEGR routine for finding eigenvalues. info = ",info
            stop
         end if
         if (N_copy < N) then
            print *,"Error! LAPACK subroutine DSTEGR only found ",N_copy," of ",N," requested eigenvalues"
            stop
         end if
         ! The first grid point should be nearly 0, but it often comes out to be ~ 1e-17, and if it ends up being negative,
         ! we get a problem elsewhere with interpolation. Therefore, set the grid point to be exactly 0:
         abscissae(1) = 0.0d+0
         weights = c(1) * eigenvectors(1, :) * eigenvectors(1, :)
      else
         ! In this case, the Jacobi matrix is also positive-definite, so use a LAPACK subroutine which exploits this property for extra accuracy.
#if defined(PETSC_USE_REAL_SINGLE)
         call SPTEQR('I', N, a, b(2:N), eigenvectors, N, LAPACKWorkspace, info)
#else
         call DPTEQR('I', N, a, b(2:N), eigenvectors, N, LAPACKWorkspace, info)
#endif
         if (info /= 0) then
            print *,"Error in LAPACK's DPTEQR routine for finding eigenvalues. info = ",info
            stop
         end if
         ! Reverse order
         do i=1,N
            abscissae(i) = a(N+1-i)
            weights(i) = c(1) * eigenvectors(1, N+1-i) * eigenvectors(1, N+1-i)
         end do
      end if

!!$      deallocate(a)
!!$      deallocate(b)
!!$      deallocate(c)
!!$      deallocate(d)
!!$      deallocate(eigenvectors)
!!$      deallocate(LAPACKWorkspace)
!!$      deallocate(LAPACK_iwork)
!!$      deallocate(LAPACK_ISUPPZ)
    end subroutine makeXGrid

    function weight(x)
      PetscScalar, intent(in) :: x
      PetscScalar :: weight
      weight = exp(-x*x)
    end function weight

    function evaluatePolynomial(x)
      PetscScalar, intent(in) :: x
      PetscScalar :: evaluatePolynomial
      PetscScalar :: pjMinus1, pj, y
      integer :: ii

      y = 0d+0
      if (j == 1) then
         evaluatePolynomial = 1d+0
      else
         pjMinus1 = 0d+0
         pj = 1d+0
         do ii=1,j-1
            y = (x-a(ii)) * pj - b(ii) * pjMinus1
            pjMinus1 = pj
            pj = y
         end do
         evaluatePolynomial = y
      end if

    end function evaluatePolynomial

    function integrandWithoutX(x)
      PetscScalar, intent(in) :: x
      PetscScalar :: integrandWithoutX, p

      p = evaluatePolynomial(x)
      integrandWithoutX = p*weight(x)*p
    end function integrandWithoutX

    function integrandWithX(x)
      PetscScalar, intent(in) :: x
      PetscScalar :: integrandWithX, p

      p = evaluatePolynomial(x)
      integrandWithX = x*p*weight(x)*p
    end function integrandWithX

  end module xGrid

