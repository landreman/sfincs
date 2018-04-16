#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

module adjointDiagnostics

  contains

    subroutine computedRadialCurrentdEr(solutionVec,adjointSolutionJr)

      use globalVariables
      use petscvec

      Vec :: solutionVec, adjointSolutionJr, adjointResidualEr
      PetscErrorCode :: ierr
      PetscScalar :: innerProductResult

      call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointResidualEr, ierr)
      call VecSet(adjointResidualEr,zero,ierr)
      ! compute (dLdEr f - dSdEr)
      call evaluateAdjointInnerProductFactor(solutionVec,adjointResidualEr,0,0)
      call innerProduct(adjointSolutionJr,adjointResidualEr,innerProductResult,0)

      ! dPhiHatdpsiHat = ddrHat2ddpsiHat * (-Er)
      ! dPhiHatdpsiHatdEr = -ddrHat2ddpsiHat
      dRadialCurrentdEr = innerProductResult*ddrHat2ddpsiHat

      call VecDestroy(adjointResidualEr,ierr)

    end subroutine computedRadialCurrentdEr

    
    !> Compuates \f$ \partial \mathbb{L}/\partial \lambda \hat{F} - \partial \mathbb{S}/\partial \lambda \f$. To compute the derivatives of the moments, this is the quantity that is used in the inner product with the adjoint solution.
    !! @param forwardSolution Solution to forward equation.
    !! @param adjointResidual Result of residual.
    !! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
    !! @param whichMode Indicates index of ms and ns for derivative.
    subroutine evaluateAdjointInnerProductFactor(forwardSolution,adjointResidual,whichLambda,whichMode)

      use globalVariables
      use petscmat
      use petscvec

      implicit none

      Vec :: forwardSolution, adjointResidual
      integer :: whichLambda, whichMode
      Mat :: dMatrixdLambda
      Vec :: dRHSdLambda
      PetscErrorCode :: ierr

      ! Allocate and populate dMatrixdLambda
      call preallocateMatrix(dMatrixdLambda, 3)
      call populatedMatrixdLambda(dMatrixdLambda, whichLambda, whichMode)

      ! Allocate and populate dRHSdLambda
      call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, dRHSdLambda, ierr)
      call VecSet(dRHSdLambda,zero,ierr)
      call populatedRHSdLambda(dRHSdLambda, whichLambda, whichMode)

      ! Multiply dRHSdLambda by -1
      call VecScale(dRHSdLambda, -1.0d+0, ierr)
      ! Performs adjointResidual = dMatrixdLambda*forwardSolution - dRHSdLambda
      call MatMultAdd(dMatrixdLambda, forwardSolution, dRHSdLambda, adjointResidual, ierr)

      call VecDestroy(dRHSdLambda, ierr)
      call MatDestroy(dMatrixdLambda, ierr)

    end subroutine

    !> Computes inner product associated with free energy norm for deltaF and deltaG
    !! @param deltaF First quantity in inner product.
    !! @param deltaG Second quantity in inner product.
    !! @param result Result of inner product.
    !! @param fineGrid Integer denoting which grid to computer inner product on. fineGrid = 1 denotes fine, otherwise coarse.
    subroutine innerProduct(deltaF, deltaG, result, fineGrid)

      use globalVariables
      use indices
      use petscvec

      implicit none

      integer :: fineGrid
      Vec :: deltaF, deltaG
      PetscScalar :: result
      Vec :: deltaFOnProc0, deltaGOnProc0
      PetscScalar, pointer :: deltaFArray(:), deltaGArray(:)
      PetscScalar, dimension(:,:), allocatable :: sourcesF, sourcesG
      integer :: index, ispecies, itheta, izeta, ix, L
      PetscScalar :: THat, mHat, nHat
      PetscScalar, dimension(:), allocatable :: xIntegralFactor
      PetscScalar, dimension(:,:), allocatable :: thetazetaIntegralFactor
      VecScatter :: VecScatterContext
      PetscErrorCode :: ierr
      PetscScalar :: speciesResult
      integer :: whichMatrix, this_Ntheta, this_Nzeta
      PetscScalar, dimension(:), pointer :: this_thetaWeights, this_zetaWeights
      PetscScalar, dimension(:,:), pointer :: this_BHat, this_DHat

      if (fineGrid==1) then
        whichMatrix = 7
        this_Ntheta = Ntheta_fine
        this_Nzeta = Nzeta_fine

        this_DHat => DHat_fine
        this_thetaWeights => thetaWeights_fine
        this_zetaWeights => zetaWeights_fine
      else
        whichMatrix = 0
        this_Ntheta = Ntheta
        this_Nzeta = Nzeta

        this_thetaWeights => thetaWeights
        this_zetaWeights => zetaWeights
        this_DHat => DHat
      end if

      select case (discreteAdjointOption)
        case (0) ! continuous
          ! Scatter deltaF to master proc
          call VecScatterCreateToZero(deltaF, VecScatterContext, deltaFOnProc0, ierr)
          call VecScatterBegin(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call VecScatterEnd(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
          if (masterProc) then
            ! Convert the PETSc vector into a normal Fortran array
            call VecGetArrayF90(deltaFOnProc0, deltaFArray, ierr)
          end if

          ! Scatter deltaG to master proc
          call VecScatterCreateToZero(deltaG, VecScatterContext, deltaGOnProc0, ierr)
          call VecScatterBegin(VecScatterContext, deltaG, deltaGOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call VecScatterEnd(VecScatterContext, deltaG, deltaGOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
          if (masterProc) then
            ! Convert the PETSc vector into a normal Fortran array
            call VecGetArrayF90(deltaGOnProc0, deltaGArray, ierr)
          end if

          if (masterProc) then
            allocate(sourcesF(Nspecies,2))
            allocate(sourcesG(Nspecies,2))
            allocate(xIntegralFactor(Nx))
            allocate(thetazetaIntegralFactor(this_Ntheta,this_Nzeta))

            ! Get sources associated with deltaF and deltaG
            if (constraintScheme == 1) then
              do ispecies = 1,Nspecies
                sourcesF(ispecies,1) = deltaFArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT,0)+1)
                sourcesF(ispecies,2) = deltaFArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT,0)+1)
                sourcesG(ispecies,1) = deltaGArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT,0)+1)
                sourcesG(ispecies,2) = deltaGArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT,0)+1)
              end do
            end if

            !> Compute first term in inner product
            result = zero

            do ispecies = 1,Nspecies
              speciesResult = zero
              THat = THats(ispecies)
              mHat = mHats(ispecies)
              nHat = nHats(ispecies)

              xIntegralFactor = ((pi*pi*sqrtpi*THat*THat*THat*THat)/(mHat*mHat*mHat*nHat*VprimeHat))*x*x*exp(x*x)
              thetazetaIntegralFactor = one/this_DHat
              do itheta=1,this_Ntheta
                do izeta=1,this_Nzeta
                  do ix=1,Nx
                    ! The integral over xi turns into a sum over L (with a factor of 2/(2L+1))
                    do L=0,(Nxi_for_x(ix)-1)
                      index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,whichMatrix)+1

                      speciesResult = speciesResult + thetazetaIntegralFactor(itheta,izeta)*xIntegralFactor(ix) &
                        *(two/(two*real(L)+one))*this_thetaWeights(itheta)*this_zetaWeights(izeta) &
                        *deltaFArray(index)*deltaGArray(index)*xWeights(ix)

                    end do
                  end do
                end do
              end do
              result = result + speciesResult
              ! Now add terms which sum over sources
              result = result + (two*pi*(THat**4)/((mHat**3)*nHat*VPrimeHat))*(sourcesF(ispecies,1)*sourcesG(ispecies,1) + sourcesF(ispecies,2)*sourcesG(ispecies,2))
            end do ! ispecies
          end if ! masterProc

        case (1) ! discrete
          call VecDot(deltaF, deltaG, result, ierr)
        case (2) ! rescaled
          call VecDot(deltaF, deltaG, result, ierr)
      end select
    end subroutine

    !> Evaluates the term in the sensitivity derivative of the heat flux
    !! which arise due to the sensitivity of the inner product and the integrating factor,
    !! not f1 itself, hence this is not called with the adjoint variable
    !! @param result Output of sensitivity calculation.
    !! @param deltaF \f$\hat{F}\f$ solution from forward equation.
    !! @param whichSpecies If = 0, summed over species. If nonzero, indicates species number. In this case should always be 0.
    !! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
    !! @param whichMode Indicates index of ms and ns for derivative.
    subroutine heatFluxSensitivity(result, forwardSolution, whichSpecies, whichLambda, whichMode)

    use globalVariables
    use indices
    use petscvec

    implicit none

    PetscScalar :: result
    Vec :: forwardSolution, forwardSolutionOnProc0
    PetscScalar, pointer :: forwardSolutionArray(:)
    integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, dBHatdThetadLambda, dBHatdZetadLambda
    PetscScalar :: dBHat_sub_thetadLambda, dBHat_sub_zetadLambda, dinvDHatdLambda, factor
    PetscScalar, dimension(:), allocatable :: xIntegralFactor
    PetscScalar :: dBHatdLambda, dVPrimeHatdLambda
    integer :: m, n
    PetscScalar :: cos_angle, sin_angle, angle
    VecScatter :: VecScatterContext
    PetscErrorCode :: ierr

    m = ms(whichMode)
    n = ns(whichMode)

    result = zero

    ! Scatter deltaF to master proc
    call VecScatterCreateToZero(forwardSolution, VecScatterContext, forwardSolutionOnProc0, ierr)
    call VecScatterBegin(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    if (masterProc) then
      ! Convert the PETSc vector into a normal Fortran array
      call VecGetArrayF90(forwardSolutionOnProc0, forwardSolutionArray, ierr)
    end if

    if (whichLambda==1 .and. geometryScheme /= 5) then
      dVPrimeHatdLambda = zero
      do itheta=1,Ntheta
        do izeta=1,Nzeta
          angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
          cos_angle = cos(angle)

          dVPrimeHatdLambda = dVPrimeHatdLambda - two*thetaWeights(itheta)*zetaWeights(izeta)*cos_angle/(DHat(itheta,izeta)*BHat(itheta,izeta))
        end do
      end do
    end if

    if (masterProc) then
      allocate(xIntegralFactor(Nx))

      if (whichSpecies == 0) then
        minSpecies = 1
        maxSpecies = NSpecies
      else
        minSpecies = whichSpecies
        maxSpecies = whichSpecies
      end if

      do ispecies = minSpecies,maxSpecies
        THat = THats(ispecies)
        mHat = mHats(ispecies)
        sqrtTHat = sqrt(THats(ispecies))
        sqrtmHat = sqrt(mHats(ispecies))

        ! This is everything independent of geometry
        xIntegralFactor = x*x*x*x*x*x*pi*THat*THat*THat*sqrtTHat*Delta/(2*mHat*sqrtmHat*Zs(ispecies))*ddrN2ddpsiHat

        do itheta=1,Ntheta
          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            sin_angle = sin(angle)
            select case (whichLambda)
              case (1) ! BHat
                dBHatdThetadLambda = -m*sin_angle
                dBHatdZetadLambda = n*Nperiods*sin_angle
                dBHatdLambda = cos_angle
                factor = (BHat_sub_theta(itheta,izeta)*dBHatdZetadLambda - BHat_sub_zeta(itheta,izeta)*dBHatdThetadLambda)/(VPrimeHat*BHat(itheta,izeta)**3) - 3*(BHat_sub_theta(itheta,izeta)*dBHatdZeta(itheta,izeta) - &
                  BHat_sub_zeta(itheta,izeta)*dBHatdTheta(itheta,izeta))*dBHatdLambda/(VPrimeHat*BHat(itheta,izeta)**4)
                if (geometryScheme /= 5) then ! Boozer
                  factor = factor - (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda/ &
                    (VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
                end if
              case (2) ! BHat_sup_theta
                factor = 0
              case (3) ! BHat_sup_zeta
                factor = 0
              case (4) ! BHat_sub_theta
                dBHat_sub_thetadLambda = cos_angle
                factor = dBHatdzeta(itheta,izeta)*dBHat_sub_thetadLambda/(VPrimeHat*BHat(itheta,izeta)**3)
              case (5) ! BHat_sub_zeta
                dBHat_sub_zetadLambda = cos_angle
                factor = -dBHatdtheta(itheta,izeta)*dBHat_sub_zetadLambda/(VPrimeHat*BHat(itheta,izeta)**3)
              case (6) ! DHat
                if (ns(whichMode)==0 .and. ms(whichMode)==0) then
                  dVPrimeHatdLambda = four*pi*pi
                  factor = - (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda/ &
                    (VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
                else
                  factor = zero
                end if
            end select

            do ix=1,Nx

                L = 0
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                result = result + &
                  (8/three)*factor*xWeights(ix)*forwardSolutionArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

                L = 2
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                result = result + &
                  (four/15)*factor*xWeights(ix)*forwardSolutionArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

              end do
            end do
          end do
        end do
      end if !masterProc

    end subroutine

    !> Evaluates the term in the sensitivity derivative of the particle fluxes
    !! which arise due to the sensitivity of the inner product and the integrating factor,
    !! not f1 itself, hence this is not called with the adjoint variable
    !! @param result Output of sensitivity calculation.
    !! @param deltaF \f$\hat{F}\f$ solution from forward equation.
    !! @param whichSpecies If = 0, summed over species. If nonzero, indicates species number. In this case should always be 0.
    !! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
    !! @param whichMode Indicates index of ms and ns for derivative.
    subroutine particleFluxSensitivity(result, deltaF, whichSpecies, whichLambda, whichMode)

    use globalVariables
    use indices
    use petscvec

    implicit none

    PetscScalar :: result
    Vec :: deltaF, deltaFOnProc0
    PetscScalar, pointer :: deltaFArray(:)
    integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, dBHatdThetadLambda, dBHatdZetadLambda
    PetscScalar :: dBHat_sub_thetadLambda, dBHat_sub_zetadLambda, dinvDHatdLambda, factor
    PetscScalar, dimension(:), allocatable :: xIntegralFactor
    PetscScalar :: dBHatdLambda, dVPrimeHatdLambda
    integer :: m, n
    PetscScalar :: angle, cos_angle, sin_angle
    VecScatter :: VecScatterContext
    PetscErrorCode :: ierr

    result = zero

    m = ms(whichMode)
    n = ns(whichMode)

    ! Scatter deltaF to master proc
    call VecScatterCreateToZero(deltaF, VecScatterContext, deltaFOnProc0, ierr)
    call VecScatterBegin(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    if (masterProc) then
      ! Convert the PETSc vector into a normal Fortran array
      call VecGetArrayF90(deltaFOnProc0, deltaFArray, ierr)
    end if

    if (whichLambda==1 .and. geometryScheme /= 5) then
      dVPrimeHatdLambda = zero
      ! For Boozer coordinats, integral needed to compute dVPrimeHatdLambda
      do itheta=1,Ntheta
        do izeta=1,Nzeta
          angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
          cos_angle = cos(angle)

          dVPrimeHatdLambda = dVPrimeHatdLambda - two*thetaWeights(itheta)*zetaWeights(izeta)*cos_angle/(DHat(itheta,izeta)*BHat(itheta,izeta))
        end do
      end do
    end if

    if (masterProc) then
      allocate(xIntegralFactor(Nx))

      if (whichSpecies == 0) then
        minSpecies = 1
        maxSpecies = NSpecies
      else
        minSpecies = whichSpecies
        maxSpecies = whichSpecies
      end if

      do ispecies = minSpecies,maxSpecies
        THat = THats(ispecies)
        mHat = mHats(ispecies)
        sqrtTHat = sqrt(THats(ispecies))
        sqrtmHat = sqrt(mHats(ispecies))

        xIntegralFactor = x*x*x*x*THat*THat*sqrtTHat*pi*Delta/(mHat*sqrtmHat*Zs(ispecies))*ddrN2ddpsiHat

        do itheta=1,Ntheta
          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            sin_angle = sin(angle)
            select case (whichLambda)
              case (1) ! BHat
                dBHatdThetadLambda = -m*sin_angle
                dBHatdZetadLambda = n*Nperiods*sin_angle
                dBHatdLambda = cos_angle
                factor = (BHat_sub_theta(itheta,izeta)*dBHatdZetadLambda - BHat_sub_zeta(itheta,izeta)*dBHatdThetadLambda)/(VPrimeHat*BHat(itheta,izeta)**3) - 3*(BHat_sub_theta(itheta,izeta)*dBHatdZeta(itheta,izeta) &
                  - BHat_sub_zeta(itheta,izeta)*dBHatdTheta(itheta,izeta))*dBHatdLambda/(VPrimeHat*BHat(itheta,izeta)**4)
                if (geometryScheme /= 5) then !Boozer
                    factor = factor - (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                      /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
                end if
              case (2) ! BHat_sup_theta
                factor = zero
              case (3) ! BHat_sup_zeta
                factor = zero
              case (4) ! BHat_sub_theta
                dBHat_sub_thetadLambda = cos_angle
                factor = dBHatdzeta(itheta,izeta)*dBHat_sub_thetadLambda/(VPrimeHat*BHat(itheta,izeta)**3)
              case (5) ! BHat_sub_zeta
                dBHat_sub_zetadLambda = cos_angle
                factor = -dBHatdtheta(itheta,izeta)*dBHat_sub_zetadLambda/(VPrimeHat*BHat(itheta,izeta)**3)
              case (6) ! DHat
                if (m==0 .and. n==0) then
                  dVPrimeHatdLambda = four*pi*pi
                  factor = - (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                      /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta))
                else
                  factor = zero
                end if
            end select

            ! Summed quantity is weighted by charge
            if (whichSpecies == 0) then
              factor = factor*Zs(ispecies)
            end if

            do ix=1,Nx

                L = 0
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                result = result + &
                  (8/three)*factor*xWeights(ix)*deltaFArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

                L = 2
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                result = result + &
                  (four/15)*factor*xWeights(ix)*deltaFArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

              end do
            end do
          end do
        end do
      end if !masterProc
    end subroutine

    !> Evaluates the term in the sensitivity derivative of the parallel flow
    !! which arise due to the sensitivity of the inner product and the integrating factor,
    !! not \f$\hat{F}\f$ itself, hence this is not called with the adjoint variable
    !! @param result Output of sensitivity calculation.
    !! @param deltaF \f$\hat{F}\f$ solution from forward equation.
    !! @param whichSpecies If = 0, summed over species. If nonzero, indicates species number. In this case should always be 0.
    !! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
    !! @param whichMode Indicates index of ms and ns for derivative.
    subroutine parallelFlowSensitivity(result, deltaF, whichSpecies, whichLambda, whichMode)

    use globalVariables
    use indices
    use petscvec

    implicit none

    PetscScalar :: result
    Vec :: deltaF
    integer :: whichSpecies, whichLambda, whichMode
    Vec :: deltaFOnProc0
    PetscScalar, pointer :: deltaFArray(:)
    integer :: minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, dBHatdThetadLambda, dBHatdZetadLambda
    PetscScalar :: dBHat_sub_thetadLambda, dBHat_sub_zetadLambda, dinvDHatdLambda, factor, dVPrimeHatdLambda
    PetscScalar, dimension(:), allocatable :: xIntegralFactor
    PetscScalar :: dBHatdLambda, dFSAB2dLambda, nHat, sqrtFSAB2
    PetscScalar :: angle, cos_angle, sin_angle
    integer :: m, n
    VecScatter :: VecScatterContext
    PetscErrorCode :: ierr
    PetscScalar :: dDHatdLambda

    sqrtFSAB2 = sqrt(FSABHat2)
    m = ms(whichMode)
    n = ns(whichMode)

    ! Scatter deltaF to master proc
    call VecScatterCreateToZero(deltaF, VecScatterContext, deltaFOnProc0, ierr)
    call VecScatterBegin(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    if (masterProc) then
      ! Convert the PETSc vector into a normal Fortran array
      call VecGetArrayF90(deltaFOnProc0, deltaFArray, ierr)
    end if

    if (whichSpecies == 0) then
      minSpecies = 1
      maxSpecies = NSpecies
    else
      minSpecies = whichSpecies
      maxSpecies = whichSpecies
    end if

    if (masterProc) then
      result = zero
      dFSAB2dLambda = zero
      ! For Boozer geometry, dFSAB2dLambda = 0
      if (whichLambda==1 .and. geometryScheme == 5) then
        do itheta=1,Ntheta
          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)

            dFSAB2dLambda = dFSAB2dLambda + two*(thetaWeights(itheta)*zetaWeights(izeta)*BHat(itheta,izeta)*cos_angle)/(DHat(itheta,izeta)*VPrimeHat)
          end do
        end do
      end if

      if (whichLambda==1 .and. geometryScheme /= 5) then
        dVPrimeHatdLambda = zero
        do itheta=1,Ntheta
          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)

            dVPrimeHatdLambda = dVPrimeHatdLambda - two*thetaWeights(itheta)*zetaWeights(izeta)*cos_angle/(DHat(itheta,izeta)*BHat(itheta,izeta))
          end do
        end do
      end if

      if (whichLambda==6) then ! DHat
        do itheta=1,Ntheta
          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)

            dFSAB2dLambda = dFSAB2dLambda + thetaWeights(itheta)*zetaWeights(izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)*cos_angle/(VPrimeHat)
          end do
        end do
        if (m == zero .and. n == zero) then
          dVPrimeHatdLambda = (four*pi*pi)
          dFSAB2dLambda = dFSAB2dLambda - (VPrimeHat*FSABHat2)*dVPrimeHatdLambda/(VPrimeHat*VPrimeHat)
        end if
      end if

      allocate(xIntegralFactor(Nx))

      do ispecies = minSpecies,maxSpecies
        THat = THats(ispecies)
        mHat = mHats(ispecies)
        sqrtTHat = sqrt(THats(ispecies))
        sqrtmHat = sqrt(mHats(ispecies))
        nHat = nHats(ispecies)

        xIntegralFactor = x*x*x*THat*THat*pi/(mHat*mHat*nHat)
        if (whichSpecies == 0) then
          xIntegralFactor = xIntegralFactor*Zs(ispecies)
        end if

        do itheta=1,Ntheta
          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)

            select case (whichLambda)
              case (1) ! BHat
                dBHatdLambda = cos_angle
                if (geometryScheme /= 5) then ! Boozer
                  dDHatdLambda = two*DHat(itheta,izeta)*dBHatdLambda/BHat(itheta,izeta)
                  dinvDHatdLambda = -dDHatdLambda/(DHat(itheta,izeta)**2)
                  factor = dBHatdLambda/(DHat(itheta,izeta)*VPrimeHat*sqrtFSAB2) + BHat(itheta,izeta)*dinvDHatdLambda/(VPrimeHat*sqrtFSAB2) &
                    -0.5*dVPrimeHatdLambda*BHat(itheta,izeta)/(DHat(itheta,izeta)*VPrimeHat*VPrimeHat*sqrtFSAB2)
                else
                  factor = dBHatdLambda/(DHat(itheta,izeta)*VPrimeHat*sqrtFSAB2) - 0.5*BHat(itheta,izeta)*dFSAB2dLambda/(DHat(itheta,izeta)*sqrtFSAB2*FSABHat2*VPrimeHat)
                end if
              case (2) ! BHat_sup_theta
                factor = 0
              case (3) ! BHat_sup_zeta
                factor = 0
              case (4) ! BHat_sub_theta
                factor = 0
              case (5) ! BHat_sub_zeta
                factor = 0
              case (6) ! DHat
                dinvDHatdLambda = cos_angle
                factor = BHat(itheta,izeta)*dinvDHatdLambda/(VPrimeHat*sqrtFSAB2) &
                  - 0.5*BHat(itheta,izeta)*dFSAB2dLambda/(DHat(itheta,izeta)*VPrimeHat*FSABHat2*sqrtFSAB2)
                if (m==0 .and. n==0) then
                  dVPrimeHatdLambda = (four*pi*pi)
                  factor = factor - BHat(itheta,izeta)*dVPrimeHatdLambda/(DHat(itheta,izeta)*VPrimeHat*VPrimeHat*sqrtFSAB2)
                end if
            end select

            do ix=1,Nx

                L = 1
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,0)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                result = result + &
                  (four/three)*factor*xWeights(ix)*deltaFArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

              end do
            end do
          end do
        end do
      end if !masterProc
    end subroutine

    !> Evaluates derivatives of inner products using the forwardSolution and adjointSolution.
    !! @param forwardSolution Solution to forward equation.
    !! @param adjointSolution Solution to adjoint equation.
    !! @param whichAdjointRHS Indicates which integrated quantity is being differentiated. If = 1 (particle flux), = 2 (heat flux), = 3 (bootstrap current).
    !! @param whichSpecies Indicates species used for inner product. If = 0, corresponds to a species-summed quantity. If nonzero, indicates number of species.
    subroutine evaluateDiagnostics(forwardSolution, adjointSolution, adjointSolutionJr, whichAdjointRHS, whichSpecies)

      use globalVariables
      use indices
      use writeHDF5Output
      use petscvec

      implicit none

      PetscErrorCode :: ierr
      Vec :: forwardSolution, adjointSolution, adjointResidual, adjointRHSVec
      integer :: whichAdjointRHS, whichSpecies
      integer :: whichLambda, whichMode
      PetscScalar :: innerProductResult, sensitivityResult, ErTermToAdd
      Vec :: adjointResidualEr, adjointSolutionJr
      PetscScalar :: innerProductResultEr1, innerProductResultEr2, innerProductResultEr3
      PetscScalar :: radialCurrentSensitivity

      if (masterProc) then
        print *,"evaluateDiagnostics was called."
        print *,"Computing adjoint diagnostics for RHS ", whichAdjointRHS, " and species ", whichSpecies
      end if

      ! Allocate adjointResidual
      call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointResidual, ierr)
      call VecSet(adjointResidual, zero, ierr)

      if (RHSMode == 5) then
        call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, adjointResidualEr, ierr)
        call VecSet(adjointResidualEr,zero,ierr)
        ! compute (dLdEr f - dSdEr)
        call evaluateAdjointInnerProductFactor(forwardSolution,adjointResidualEr,0,0)
        call innerProduct(adjointSolution,adjointResidualEr,innerProductResultEr1,0)
        call innerProduct(adjointSolutionJr,adjointResidualEr,innerProductResultEr2,0)
      end if

      ! Same inner product is formed regardless of whichAdjointMatrix and whichSpecies
      ! Loop over lambda's and perform inner product
      ! rethink this for Er
      do whichLambda=1,NLambdas
        do whichMode=1,NModesAdjoint
          call VecSet(adjointResidual,zero,ierr)
          ! Call function to perform (dLdlambdaf - dSdlambda), which calls populatedMatrixdLambda and populatedRHSdLambda
          call evaluateAdjointInnerProductFactor(forwardSolution,adjointResidual,whichLambda,whichMode)

          if (RHSMode == 5) then
            call innerProduct(adjointSolutionJr,adjointResidual,innerProductResultEr3,0)
            call particleFluxSensitivity(radialCurrentSensitivity,forwardSolution,0,whichLambda,whichMode)
            ErTermToAdd = -(innerProductResultEr1/innerProductResultEr2)*(radialCurrentSensitivity-innerProductResultEr3)
            dPhidPsidLambda(whichLambda,whichMode) = (1/innerProductResultEr2)*(radialCurrentSensitivity-innerProductResultEr3)
          end if

            call innerProduct(adjointSolution,adjointResidual,innerProductResult,0)

            ! Save to appropriate sensitivity array
            if (whichSpecies == 0) then
              select case (whichAdjointRHS)
                case (1) ! Particle flux
                  call particleFluxSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                  dRadialCurrentdLambda(whichLambda,whichMode) = -innerProductResult + sensitivityResult

                case (2) ! Heat Flux
                  call heatFluxSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                  dTotalHeatFluxdLambda(whichLambda,whichMode) = -innerProductResult + sensitivityResult
                  if (RHSMode == 5) then
                    dTotalHeatFluxdLambda(whichLambda,whichMode) = dTotalHeatFluxdLambda(whichLambda,whichMode) + ErTermToAdd
                  end if
                case (3) ! Bootstrap
                  call parallelFlowSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                  dBootstrapdLambda(whichLambda,whichMode) = -innerProductResult + sensitivityResult
                  if (RHSMode == 5) then
                    dBootstrapdLambda(whichLambda,whichMode) = dBootstrapdLambda(whichLambda,whichMode) + ErTermToAdd
                  end if
              end select
            else
              select case (whichAdjointRHS)
                case (1) ! Particle flux
                  call particleFluxSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                  dParticleFluxdLambda(whichSpecies,whichLambda,whichMode) = sensitivityResult - innerProductResult
                  if (RHSMode == 5) then
                    dParticleFluxdLambda(whichSpecies,whichLambda,whichMode) = dParticleFluxdLambda(whichSpecies,whichLambda,whichMode) + ErTermToAdd
                  end if
                case (2) ! Heat Flux
                  call heatFluxSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                  dHeatFluxdLambda(whichSpecies,whichLambda,whichMode) = sensitivityResult - innerProductResult
                  if (RHSMode == 5) then
                    dHeatFluxdLambda(whichSpecies,whichLambda,whichMode) = dHeatFluxdLambda(whichSpecies,whichLambda,whichMode) + ErTermToAdd
                  end if
                case (3) ! Parallel Flow
                  call parallelFlowSensitivity(sensitivityResult, forwardSolution,whichSpecies, whichLambda, whichMode)
                  dParallelFlowdLambda(whichSpecies,whichLambda,whichMode) = sensitivityResult - innerProductResult
                  if (RHSMode == 5) then
                    dParallelFlowdLambda(whichSpecies,whichLambda,whichMode) = dParallelFlowdLambda(whichSpecies,whichLambda,whichMode) + ErTermToAdd
                  end if
              end select
            end if

        end do
      end do

      call VecDestroy(adjointResidual,ierr)
      if (RHSMode==5) then
        call VecDestroy(adjointResidualEr,ierr)
      end if

    end subroutine

end module

