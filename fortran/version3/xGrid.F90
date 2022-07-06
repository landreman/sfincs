  ! This module computes the nodes and weights for Gaussian quadrature for any
  ! weight function on any interval.  The algorithm used is the Stieltjes
  ! procedure described in section 4.6.3 of Numerical Recipes, 3rd edition,
  ! together with the Jacobi matrix procedure of section 4.6.2.  See also
  ! www.nr.com/webnotes?3
  ! The module is presently configured to use the weight exp(-x^2)*(x^k) and to
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

#define workspaceSize 5000

  ! Call QUADPACK routines to evaluate integrals:
  ! qage = single-precision quadrature on finite domain
  ! dqage = double-precision quadrature on finite domain
  ! qagie = single-precision quadrature on (semi-)infinite domain
  ! dqagie = double-precision quadrature on (semi-)infinite domain

#if defined(PETSC_USE_REAL_SINGLE)
#define integrate sfincs_qage
#define integrate_semiinf sfincs_qagie
#define zero 0.
#else
#define integrate sfincs_dqage
#define integrate_semiinf sfincs_dqagie
#define zero 0.0d+0
#endif

! QUADPACK gives ier.ne.0 when it cannot achieve the requested tolerance.
! When this occurs, it probably just means it is not feasbible to get as many digits
! of precision as we requested, so we can safely ignore the warnings.
!#define showQuadpackWarnings .true.
#define showQuadpackWarnings .false.


#define line "******************************************************************"

  module xGrid

#include "PETScVersions.F90"

    implicit none

    private

    PetscScalar, allocatable :: a(:), b(:), c(:)
    integer :: j, integrationPower
    PetscScalar, public :: xGrid_k = 0
    public :: makeXGrid, computeRosenbluthPotentialResponse, xGrid_finalize

  contains

    subroutine makeXGrid(N, abscissae, weights, includePointAtX0)
      ! This is one of 2 main subroutines of the module.
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
      PetscScalar, allocatable :: a_copy(:), sqrtb(:), d(:), eigenvectors(:,:), LAPACKWorkspace(:)

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
      allocate(a_copy(N))
      allocate(sqrtb(N+1))
      allocate(c(N))
      allocate(d(N))
      allocate(eigenvectors(N,N))
      ! For the workspace, LAPACK's DPTEQR requires 4*N, while SSTEGR requires 18*N:
      workspace_size = 18*N
      allocate(LAPACKWorkspace(workspace_size))
      iwork_size = 10*N
      allocate(LAPACK_iwork(iwork_size))

      allocate(LAPACK_ISUPPZ(2*N))

      inf = 1
      !     inf = 1 corresponds to  (bound,+infinity),
      !     inf = -1            to  (-infinity,bound),
      !     inf = 2             to (-infinity,+infinity).

      ! If i denotes the true value of the integral, QUADPACK will try to satisfy
      ! abs(i-result).le.max(epsabs,epsrel*abs(i))
      EPSABS = 0.0D0
      EPSREL = 1d-13
      !EPSREL = 1d-3

      LIMIT = workspaceSize

      ! For some reason, quadpack returns an error unless I split the integration region into 2 pieces.
      ! Make the split at the following value of x:
      finiteBound = 10d+0

      do j=1,N

         call integrate_semiinf(integrandWithoutX,finiteBound,inf,epsabs,epsrel,limit,c(j),abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)
         if (ier /= 0 .and. showQuadpackWarnings) then
            print *,line
            print *,line
            print *,"** WARNING:  Quadrature error A: ier = ",ier
            print *,line
            print *,line
         end if

         call integrate(integrandWithoutX,zero,finiteBound,epsabs,epsrel,key,limit,amountToAdd,abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)

         if (ier /= 0 .and. showQuadpackWarnings) then
            print *,line
            print *,line
            print *,"** WARNING:  Quadrature error B: ier = ",ier
            print *,line
            print *,line
         end if
         c(j) = c(j) + amountToAdd

         call integrate_semiinf(integrandWithX,finiteBound,inf,epsabs,epsrel,limit,d(j),abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)

         if (ier /= 0 .and. showQuadpackWarnings) then
            print *,line
            print *,line
            print *,"** WARNING:  Quadrature error C: ier = ",ier
            print *,line
            print *,line
         end if

         call integrate(integrandWithX,zero,finiteBound,epsabs,epsrel,key,limit,amountToAdd,abserr, &
              neval,ier,alist,blist,rlist,elist,iord,last)

         if (ier /= 0 .and. showQuadpackWarnings) then
            print *,line
            print *,line
            print *,"** WARNING:  Quadrature error D: ier = ",ier
            print *,line
            print *,line
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

      a_copy = a
      sqrtb = sqrt(b)

      ! Call LAPACK to get eigenvectors & eigenvalues:
      ! We use the fact that the Jacobi matrix is symmetric and tridiagonal.
      if (includePointAtX0) then
         !LAPACK_abstol = DLAMCH('S')
         LAPACK_abstol = 1.0d-15
#if defined(PETSC_USE_REAL_SINGLE)
         call SSTEGR('V','A',N,a_copy,sqrtb(2:N+1),zero,zero,0,0,LAPACK_abstol,N_copy,abscissae,eigenvectors,&
              N,LAPACK_ISUPPZ, LAPACKWorkspace, workspace_size, LAPACK_iwork, iwork_size, info)
#else
         call DSTEGR('V','A',N,a_copy,sqrtb(2:N+1),zero,zero,0,0,LAPACK_abstol,N_copy,abscissae,eigenvectors,&
              N,LAPACK_ISUPPZ, LAPACKWorkspace, workspace_size, LAPACK_iwork, iwork_size, info)
#endif
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
         call SPTEQR('I', N, a_copy, sqrtb(2:N), eigenvectors, N, LAPACKWorkspace, info)
#else
         call DPTEQR('I', N, a_copy, sqrtb(2:N), eigenvectors, N, LAPACKWorkspace, info)
#endif
         if (info /= 0) then
            print *,"Error in LAPACK's DPTEQR routine for finding eigenvalues. info = ",info
            stop
         end if
         ! Reverse order
         do i=1,N
            abscissae(i) = a_copy(N+1-i)
            weights(i) = c(1) * eigenvectors(1, N+1-i) * eigenvectors(1, N+1-i)
         end do
      end if

    end subroutine makeXGrid

    ! ---------------------------------------------------------------

    subroutine computeRosenbluthPotentialResponse(Nx, x, xWeights, Nspecies, mHats, THats, nHats, Zs, NL, &
         RosenbluthPotentialTerms,verbose)

      ! This subroutine is the other major routine in this module.
      ! Here, we compute matrices which, when multiplied by a vector of distribution-function
      ! values on the x grid, yields the Rosenbluth potentials (or derivatives thereof)
      ! on the x grid. The source and destination x grids may "belong" to different species.
      
      logical, intent(in) :: verbose
      integer, intent(in) :: Nx, Nspecies, NL
      PetscScalar, dimension(:), intent(in) :: x, xWeights, mHats, THats, nHats, Zs
      PetscScalar, dimension(:,:,:,:,:), intent(out) :: RosenbluthPotentialTerms
      ! Order of indicies in the Rosenbluth response matrices:
      ! (species_row, species_col, L, x_row, x_col)

      PetscScalar, dimension(:,:), allocatable :: collocation2modal
      PetscScalar, dimension(:,:), allocatable :: tempMatrix_H
      PetscScalar, dimension(:,:), allocatable :: tempMatrix_dHdxb
      PetscScalar, dimension(:,:), allocatable :: tempMatrix_d2Gdxb2
      PetscScalar, dimension(:,:), allocatable :: tempMatrix_combined
      PetscScalar, dimension(:), allocatable :: expx2
      integer :: i, L, iSpeciesA, iSpeciesB, ix, imode
      PetscScalar :: alpha, speciesFactor, speciesFactor2, xb
      PetscScalar :: I_2pL, I_4pL, I_1mL, I_3mL
      PetscScalar, parameter :: pi = 3.14159265358979d+0
      PetscScalar, parameter :: one = 1., two = 2.

      ! Variables needed by quadpack:
      integer :: key = 6
      PetscScalar, dimension(workspaceSize) :: alist, blist, rlist, elist
      integer :: ier, inf, last, limit, neval
      integer, dimension(workspaceSize) :: iord
      PetscScalar :: abserr,EPSABS,EPSREL, amountToAdd, partition



      ! Build a Nx * Nx matrix which, when multiplying a vector of function values on the x grid,
      ! yields the amplitudes multiplying the orthogonal polynomial modes:
      allocate(collocation2modal(Nx,Nx))
      do j=1,Nx
         do i=1,Nx
            collocation2modal(j,i) = xWeights(i) * (x(i)**xGrid_k) * evaluatePolynomial(x(i)) / c(j)
         end do
      end do
      
      if (verbose) then
         print *,"collocation2modal:"
         do ix=1,Nx
            print *,collocation2modal(ix,:)
         end do
      end if

      allocate(tempMatrix_H(Nx,Nx))
      allocate(tempMatrix_dHdxb(Nx,Nx))
      allocate(tempMatrix_d2Gdxb2(Nx,Nx))
      allocate(tempMatrix_combined(Nx,Nx))

      allocate(expx2(Nx))
      expx2=exp(-x*x)

      inf = 1
      !     inf = 1 corresponds to  (bound,+infinity),
      !     inf = -1            to  (-infinity,bound),
      !     inf = 2             to (-infinity,+infinity).

      ! If i denotes the true value of the integral, QUADPACK will try to satisfy
      ! abs(i-result).le.max(epsabs,epsrel*abs(i))

      ! Note: If epsabs is set .le. 1d-13 or so, then sometimes there are AA errors (ier=2) indicating pollution by roundoff.
      EPSABS = 1d-13
      EPSREL = 1d-13
      LIMIT = workspaceSize

      ! This next part could be optimized in several ways: the response for self-collisions
      ! only needs to be computed once, and several of the integrals are computed multiple times.
      ! However, this code takes a negligible amount of time compared to the main solve, so there
      ! is no point in optimizing here.
      do L=0,(NL-1)
         if (verbose) then
            print *,"Computing Rosenbluth potential response matrices for L = ",L
         end if

         alpha = -(2*L-1)/(two*L+3)
         !partition = maxval(x) * 2

         do iSpeciesA = 1,Nspecies
            do iSpeciesB = 1,Nspecies
               speciesFactor = sqrt(THats(iSpeciesA)*mHats(iSpeciesB) &
                    / (THats(iSpeciesB) * mHats(iSpeciesA)))

               speciesFactor2 = 3/(2*pi)*nHats(iSpeciesA) &
                    * Zs(iSpeciesA)*Zs(iSpeciesA)*Zs(iSpeciesB)*Zs(iSpeciesB) &
                    / (THats(iSpeciesA) * sqrt(THats(iSpeciesA)*mHats(ispeciesA))) &
                    * THats(iSpeciesB)*mHats(iSpeciesA)/(THats(iSpeciesA)*mHats(iSpeciesB))

               do ix = 1,Nx
                  xb =  x(ix) * speciesFactor
                  ! j denotes which polynomial, beginning with 1 rather than 0.
                  do j = 1,Nx

                     ! Call QUADPACK routines to evaluate integrals:

                     integrationPower = L+2

                     call integrate(integrandWithPower,zero,xb,epsabs,epsrel,key,limit,I_2pL,abserr, &
                          neval,ier,alist,blist,rlist,elist,iord,last)

                     if (ier /= 0 .and. showQuadpackWarnings) then
                        print *,line
                        print *,line
                        print *,"** WARNING:  Quadrature error AA: ier = ",ier
                        print *,"xb=",xb,", j=",j,", ix=",ix,", iSpeciesA=",iSpeciesA,", iSpeciesB=",iSpeciesB
                        print *,line
                        print *,line
                     end if

                     ! --------------------------------------
                     integrationPower = L+4

                     call integrate(integrandWithPower,zero,xb,epsabs,epsrel,key,limit,I_4pL,abserr, &
                          neval,ier,alist,blist,rlist,elist,iord,last)

                     if (ier /= 0 .and. showQuadpackWarnings) then
                        print *,line
                        print *,line
                        print *,"** WARNING:  Quadrature error BB: ier = ",ier
                        print *,"xb=",xb,", j=",j,", ix=",ix,", iSpeciesA=",iSpeciesA,", iSpeciesB=",iSpeciesB
                        print *,line
                        print *,line
                     end if

                     ! --------------------------------------
                     integrationPower = 1-L
                     partition = max(10.0,2*xb)

                     ! [xb, partition)
                     call integrate(integrandWithPower,xb,partition,epsabs,epsrel,key,limit,I_1mL,abserr, &
                          neval,ier,alist,blist,rlist,elist,iord,last)

                     if (ier /= 0 .and. showQuadpackWarnings) then
                        print *,line
                        print *,line
                        print *,"** WARNING:  Quadrature error CC: ier = ",ier
                        print *,"xb=",xb,", j=",j,", ix=",ix,", iSpeciesA=",iSpeciesA,", iSpeciesB=",iSpeciesB
                        print *,line
                        print *,line
                     end if

                     ! [partition, Inf)
                     call integrate_semiinf(integrandWithPower,partition,inf,epsabs,epsrel,limit,amountToAdd,abserr, &
                          neval,ier,alist,blist,rlist,elist,iord,last)

                     I_1mL = I_1mL + amountToAdd

                     if (ier /= 0 .and. showQuadpackWarnings) then
                        print *,line
                        print *,line
                        print *,"** WARNING:  Quadrature error DD: ier = ",ier
                        print *,"xb=",xb,", j=",j,", ix=",ix,", iSpeciesA=",iSpeciesA,", iSpeciesB=",iSpeciesB
                        print *,line
                        print *,line
                     end if

                     ! --------------------------------------
                     integrationPower = 3-L

                     ! [xb, partition]
                     call integrate(integrandWithPower,xb,partition,epsabs,epsrel,key,limit,I_3mL,abserr, &
                          neval,ier,alist,blist,rlist,elist,iord,last)

                     if (ier /= 0 .and. showQuadpackWarnings) then
                        print *,line
                        print *,line
                        print *,"** WARNING:  Quadrature error EE: ier = ",ier
                        print *,"xb=",xb,", j=",j,", ix=",ix,", iSpeciesA=",iSpeciesA,", iSpeciesB=",iSpeciesB
                        print *,line
                        print *,line
                     end if

                     ! [partition, Inf)
                     call integrate_semiinf(integrandWithPower,partition,inf,epsabs,epsrel,limit,amountToAdd,abserr, &
                          neval,ier,alist,blist,rlist,elist,iord,last)

                     I_3mL = I_3mL + amountToAdd

                     if (ier /= 0 .and. showQuadpackWarnings) then
                        print *,line
                        print *,line
                        print *,"** WARNING:  Quadrature error FF: ier = ",ier
                        print *,"xb=",xb,", j=",j,", ix=",ix,", iSpeciesA=",iSpeciesA,", iSpeciesB=",iSpeciesB
                        print *,line
                        print *,line
                     end if


                     ! The next equations can be found in my notes 20150330-03:

                     tempMatrix_H(ix,j) = 4*pi/(two*L+1) * (I_2pL/(xb**(L+1)) + (xb**L)*I_1mL)

                     tempMatrix_dHdxb(ix,j) = 4*pi/(two*L+1) &
                          * (-(L+1)*I_2pL/(xb**(L+2)) + L*(xb**(L-1))*I_1mL)

                     tempMatrix_d2Gdxb2(ix,j) = -4*pi/(4*L*L-one) * ( &
                          L*(L-1)*(xb**(L-2))*I_3mL &
                          + alpha*(L+1)*(L+2)*(xb**L)*I_1mL &
                          + alpha*(L+1)*(L+2)/(xb**(L+3))*I_4pL &
                          + L*(L-1)/(xb**(L+1))*I_2pL)

                  end do
               end do

               tempMatrix_combined=0
               do i=1,Nx
                  tempMatrix_combined(i,:) = speciesFactor2*expx2(i)*( &
                       - tempMatrix_H(i,:) &
                       - (1 - mHats(iSpeciesA)/mHats(iSpeciesB)) * (x(i) * speciesFactor) * tempMatrix_dHdxb(i,:) &
                       + x(i)*x(i) * tempMatrix_d2Gdxb2(i,:))
               end do

               RosenbluthPotentialTerms(iSpeciesA, iSpeciesB, L+1, :, :) = matmul(tempMatrix_combined, collocation2modal)

               if (verbose) then
                  print *,"For iSpeciesA=",iSpeciesA,",iSpeciesB=",iSpeciesB,",L=",L
                  print *,"RosenbluthPotentialTerms="
                  do ix=1,NX
                     print *,RosenbluthPotentialTerms(iSpeciesA,iSpeciesB,L+1,ix,:)
                  end do
               end if
            end do
         end do
      end do

      deallocate(tempMatrix_H)
      deallocate(tempMatrix_dHdxb)
      deallocate(tempMatrix_d2Gdxb2)
      deallocate(tempMatrix_combined)
      deallocate(collocation2modal)

    end subroutine computeRosenbluthPotentialResponse

    ! ---------------------------------------------------------------

    function weight(x)
      PetscScalar, intent(in) :: x
      PetscScalar :: weight
      weight = exp(-x*x)*(x ** xGrid_k)
    end function weight

    ! ---------------------------------------------------------------

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

    ! ---------------------------------------------------------------

    function integrandWithoutX(x)
      PetscScalar, intent(in) :: x
      PetscScalar :: integrandWithoutX, p

      p = evaluatePolynomial(x)
      integrandWithoutX = p*weight(x)*p
    end function integrandWithoutX

    ! ---------------------------------------------------------------

    function integrandWithX(x)
      PetscScalar, intent(in) :: x
      PetscScalar :: integrandWithX, p

      p = evaluatePolynomial(x)
      integrandWithX = x*p*weight(x)*p
    end function integrandWithX

    ! ---------------------------------------------------------------

    function integrandWithPower(x)
      PetscScalar, intent(in) :: x
      PetscScalar :: integrandWithPower, p

      p = evaluatePolynomial(x)
      ! Note that x**xGrid_k should not be included in the next line!
      integrandWithPower = (x**integrationPower)*p*exp(-x*x)
    end function integrandWithPower

    ! ---------------------------------------------------------------

    subroutine xGrid_finalize()

      if(allocated(a)) deallocate(a,b,c)

    end subroutine xGrid_finalize

  end module xGrid

