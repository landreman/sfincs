module adjointDiagnostics

#include "PETScVersions.F90"

implicit none

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
      use DKEMatrix, only : populateMatrix
      use residual, only : DKERHS

      implicit none

      Vec :: forwardSolution, adjointResidual
      integer :: whichLambda, whichMode
      Mat :: dMatrixdLambda
      Vec :: dRHSdLambda
      PetscErrorCode :: ierr
      Vec :: dummyVec

      ! Allocate and populate dMatrixdLambda
      call preallocateMatrix(dMatrixdLambda, 6)
      call populateMatrix(dMatrixdLambda, 6, dummyVec, whichLambda, whichMode)
      ! call populatedMatrixdLambda(dMatrixdLambda, whichLambda, whichMode)
      if (masterProc) then
         print *,"Exited populateMatrix"
      end if
    

      ! Allocate and populate dRHSdLambda
      call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, dRHSdLambda, ierr)
      call VecSet(dRHSdLambda,zero,ierr)

      !!!!!!!!!!!!!!!!!!!
      call DKERHS(dRHSdLambda, whichLambda, whichMode)
      call VecAssemblyBegin(dRHSdLambda, ierr)
      call VecAssemblyEnd(dRHSdLambda, ierr)
      
      ! Multiply dRHSdLambda by -1
      call VecScale(dRHSdLambda, -1.0d+0, ierr)
      ! Performs adjointResidual = dMatrixdLambda*forwardSolution - dRHSdLambda
      call MatMultAdd(dMatrixdLambda, forwardSolution, dRHSdLambda, adjointResidual, ierr)

      call VecDestroy(dRHSdLambda, ierr)
      call MatDestroy(dMatrixdLambda, ierr)
      !!!!!!!!!!!!!!!!!!!

      if (masterProc) then
         print *,"set values, will exit evaluateAdjointInnerProductFactor"
      end if
    


    end subroutine

    !> Computes inner product associated with free energy norm for deltaF and deltaG
    !! @param deltaF First quantity in inner product.
    !! @param deltaG Second quantity in inner product.
    !! @param result Result of inner product.
    !! @param fineGrid Integer denoting which grid to computer inner product on. fineGrid = 1 denotes fine, otherwise coarse.
    !! fineGrid does not control anything at this point, but is included for future incorporation of adjointEC
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

      if (discreteAdjointOption .eqv. .false.) then
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
            allocate(thetazetaIntegralFactor(Ntheta,Nzeta))

            ! Get sources associated with deltaF and deltaG
            if (constraintScheme == 1) then
              do ispecies = 1,Nspecies
                sourcesF(ispecies,1) = deltaFArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)+1)
                sourcesF(ispecies,2) = deltaFArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)+1)
                sourcesG(ispecies,1) = deltaGArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)+1)
                sourcesG(ispecies,2) = deltaGArray(getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)+1)
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
              thetazetaIntegralFactor = one/DHat
              do itheta=1,Ntheta
                do izeta=1,Nzeta
                  do ix=1,Nx
                    ! The integral over xi turns into a sum over L (with a factor of 2/(2L+1))
                    do L=0,(Nxi_for_x(ix)-1)
                      index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1

                      speciesResult = speciesResult + thetazetaIntegralFactor(itheta,izeta)*xIntegralFactor(ix) &
                        *(two/(two*real(L)+one))*thetaWeights(itheta)*zetaWeights(izeta) &
                        *deltaFArray(index)*deltaGArray(index)*xWeights(ix)

                    end do
                  end do
                end do
              end do
              result = result + speciesResult
              ! Now add terms which sum over sources
              result = result + (two*pi*(THat**4)/((mHat**3)*nHat*VPrimeHat))*(sourcesF(ispecies,1)*sourcesG(ispecies,1) + sourcesF(ispecies,2)*sourcesG(ispecies,2))
            end do ! ispecies
            call VecDestroy(deltaFOnProc0,ierr)
            call VecDestroy(deltaGOnProc0,ierr)
          end if ! masterProc
       else ! discrete
          call VecDot(deltaF, deltaG, result, ierr)
       end if

      end subroutine innerProduct
      

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

    m = ms_sensitivity(whichMode)
    n = ns_sensitivity(whichMode)

    result = zero

    ! Scatter deltaF to master proc
    call VecScatterCreateToZero(forwardSolution, VecScatterContext, forwardSolutionOnProc0, ierr)
    call VecScatterBegin(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    if (masterProc) then
      ! Convert the PETSc vector into a normal Fortran array
      call VecGetArrayF90(forwardSolutionOnProc0, forwardSolutionArray, ierr)
    end if

        if (whichLambda==1) then ! BHat
                dVPrimeHatdLambda = zero
                do itheta=1,Ntheta
                        do izeta=1,Nzeta
                                angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                                cos_angle = cos(angle)

                                dVPrimeHatdLambda = dVPrimeHatdLambda - two*thetaWeights(itheta)*zetaWeights(izeta)*cos_angle/(DHat(itheta,izeta)*BHat(itheta,izeta))
                        end do
                end do
        else if (whichLambda==2) then ! IHat
                dVPrimeHatdLambda = iota*VprimeHat/(GHat+iota*IHat)
        else if (whichLambda==3) then ! GHat
                dVPrimeHatdLambda = VPrimeHat/(GHat+iota*IHat)
        else if (whichLambda==4) then ! iota
                dVPrimeHatdLambda = IHat*VPrimeHat/(GHat+iota*IHat)
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

        ! factor = (BHat_sub_theta*dBHatdzeta-BHat_sub_zeta*dBHatdtheta)/(VPrimeHat*BHat**3)
        !                 = (IHat*dBHatdzeta-GHat*dBHatdtheta)/(VPrimeHat*BHat**3)
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
                factor = (BHat_sub_theta(itheta,izeta)*dBHatdZetadLambda &
                    - BHat_sub_zeta(itheta,izeta)*dBHatdThetadLambda)/(VPrimeHat*BHat(itheta,izeta)**3) &
                    - 3*(BHat_sub_theta(itheta,izeta)*dBHatdZeta(itheta,izeta) &
                    - BHat_sub_zeta(itheta,izeta)*dBHatdTheta(itheta,izeta))*dBHatdLambda/(VPrimeHat*BHat(itheta,izeta)**4)
                factor = factor - (IHat*dBHatdzeta(itheta,izeta) - GHat*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda/ &
                    (VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
              case (2) ! IHat
                factor = dBHatdzeta(itheta,izeta)/(VPrimeHat*BHat(itheta,izeta)**3) &
                    - dVPrimeHatdLambda*(IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta)) &
                    /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
              case (3) ! GHat
                factor = -dBHatdtheta(itheta,izeta)/(VPrimeHat*BHat(itheta,izeta)**3) &
                    - dVPrimeHatdLambda*(IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta)) &
                    /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
              case (4) ! iota
                factor = -dVPrimeHatdLambda*(IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta)) &
                    /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
            end select

            if (abs(factor)>0) then
                do ix=1,Nx

                    L = 0
                    index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                    ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                    result = result + &
                    (8/three)*factor*xWeights(ix)*forwardSolutionArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

                    L = 2
                    index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                    ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                    result = result + &
                    (four/15)*factor*xWeights(ix)*forwardSolutionArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

                end do
            end if
        end do
    end do
end do
        call VecDestroy(forwardSolutionOnProc0,ierr)
      end if !masterProc
        call VecScatterDestroy(VecScatterContext,ierr)

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

    m = ms_sensitivity(whichMode)
    n = ns_sensitivity(whichMode)

    ! Scatter deltaF to master proc
    call VecScatterCreateToZero(deltaF, VecScatterContext, deltaFOnProc0, ierr)
    call VecScatterBegin(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    if (masterProc) then
      ! Convert the PETSc vector into a normal Fortran array
      call VecGetArrayF90(deltaFOnProc0, deltaFArray, ierr)
    end if

        if (whichLambda==1) then ! BHat
            dVPrimeHatdLambda = zero
            do itheta=1,Ntheta
                do izeta=1,Nzeta
                    angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                    cos_angle = cos(angle)

                    dVPrimeHatdLambda = dVPrimeHatdLambda - two*thetaWeights(itheta)*zetaWeights(izeta)*cos_angle/(DHat(itheta,izeta)*BHat(itheta,izeta))
                end do
            end do
        else if (whichLambda==2) then ! IHat
            dVPrimeHatdLambda = iota*VprimeHat/(GHat+iota*IHat)
        else if (whichLambda==3) then ! GHat
            dVPrimeHatdLambda = VPrimeHat/(GHat+iota*IHat)
        else if (whichLambda==4) then ! iota
            dVPrimeHatdLambda = IHat*VPrimeHat/(GHat+iota*IHat)
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

            ! factor = (BHat_sub_theta*dBHatdzeta-BHat_sub_zeta*dBHatdtheta)/(VPrimeHat*BHat**3)
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
                factor = (IHat*dBHatdzetadLambda-GHat*dBHatdthetadLambda)/(VPrimeHat*BHat(itheta,izeta)**3) &
                    - 3*(IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta))/(VprimeHat*BHat(itheta,izeta)**4) &
                    * dBHatdLambda - (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                    /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
              case (2) ! IHat
                factor = dBHatdzeta(itheta,izeta)/(VPrimeHat*BHat(itheta,izeta)**3) &
                    - (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta)) &
                    * dVPrimeHatdLambda/(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
              case (3) ! GHat
                factor = -dBHatdtheta(itheta,izeta)/(VPrimeHat*BHat(itheta,izeta)**3) &
                    - (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                    / (VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
              case (4) ! iota
                factor = - (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                    / (VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
            end select

            ! Summed quantity is weighted by charge
            if (whichSpecies == 0) then
              factor = factor*Zs(ispecies)
            end if

            do ix=1,Nx

                L = 0
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                result = result + &
                  (8/three)*factor*xWeights(ix)*deltaFArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

                L = 2
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                result = result + &
                  (four/15)*factor*xWeights(ix)*deltaFArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

              end do ! ix
            end do ! izeta
          end do ! itheta
        end do ! ispecies
        call VecDestroy(deltaFOnProc0,ierr)
      end if !masterProc
      call VecScatterDestroy(VecScatterContext,ierr)
    end subroutine particleFluxSensitivity

    !> Evaluates the term in the sensitivity derivative of the particle fluxes
    !! which arise due to the sensitivity of the inner product and the integrating factor,
    !! not f1 itself, hence this is not called with the adjoint variable
    !! @param result Output of sensitivity calculation.
    !! @param deltaF \f$\hat{F}\f$ solution from forward equation.
    !! @param whichSpecies If = 0, summed over species. If nonzero, indicates species number. In this case should always be 0.
    !! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
    !! @param whichMode Indicates index of ms and ns for derivative.
    subroutine particleFlux_vd_Sensitivity(result, deltaF, whichSpecies, whichLambda, whichMode)

      use globalVariables
      use indices
      use petscvec

      implicit none

      PetscScalar :: result
      Vec :: deltaF, deltaFOnProc0, f0OnProc0
      PetscScalar, pointer :: deltaFArray(:), f0Array(:)
      integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
      PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, Z
      PetscScalar :: dBHat_sub_thetadLambda, dBHat_sub_zetadLambda, dinvDHatdLambda, factor, factorExB
      PetscScalar :: dfactordLambda, dfactorExBdLambda, df0dLambda
      PetscScalar, dimension(:), allocatable :: xIntegralFactor,  xIntegralFactorExB
      PetscScalar :: dBHatdLambda, dVPrimeHatdLambda, dBHatdThetadLambda, dBHatdZetadLambda
      PetscScalar :: dPhi1HatdLambda, dPhi1HatdThetadLambda, dPhi1HatdZetadLambda
      integer :: m, n
      PetscScalar :: angle, cos_angle, sin_angle
      VecScatter :: VecScatterContext, VecScatterContext0 !MAYBE REMOVE SECOND
      PetscErrorCode :: ierr

      PetscScalar :: result3
      PetscScalar, dimension(:,:), allocatable :: results2

      result = zero
      allocate(results2(Ntheta,Nzeta))
      results2 = zero
      result3 = zero

      m = ms_sensitivity(whichMode)
      n = ns_sensitivity(whichMode)

      ! Scatter deltaF to master proc
      call VecScatterCreateToZero(deltaF, VecScatterContext, deltaFOnProc0, ierr)
      call VecScatterBegin(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
      if (masterProc) then
         ! Convert the PETSc vector into a normal Fortran array
         call VecGetArrayF90(deltaFOnProc0, deltaFArray, ierr)
      end if

      if (includePhi1) then
         ! f0 is defined in global
         call VecScatterCreateToZero(f0, VecScatterContext0, f0OnProc0, ierr)
         call VecScatterBegin(VecScatterContext0, f0, f0OnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
         call VecScatterEnd(VecScatterContext0, f0, f0OnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
         if (masterProc) then
            call VecGetArrayF90(f0OnProc0, f0Array, ierr)
         end if
      end if


      if (whichLambda==1) then ! BHat
         dVPrimeHatdLambda = zero
         do itheta=1,Ntheta
            do izeta=1,Nzeta
               angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
               cos_angle = cos(angle)

               dVPrimeHatdLambda = dVPrimeHatdLambda - two*thetaWeights(itheta)*zetaWeights(izeta)*cos_angle/(DHat(itheta,izeta)*BHat(itheta,izeta))
            end do
         end do
      else if (whichLambda==2) then ! IHat
         dVPrimeHatdLambda = iota*VprimeHat/(GHat+iota*IHat)
      else if (whichLambda==3) then ! GHat
         dVPrimeHatdLambda = VPrimeHat/(GHat+iota*IHat)
      else if (whichLambda==4) then ! iota
         dVPrimeHatdLambda = IHat*VPrimeHat/(GHat+iota*IHat)
      end if

      if (masterProc) then
         allocate(xIntegralFactor(Nx))
         allocate(xIntegralFactorExB(Nx))

         if (whichSpecies == 0) then
            minSpecies = 1
            maxSpecies = NSpecies
         else
            minSpecies = whichSpecies
            maxSpecies = whichSpecies
         end if

         ! factor = (BHat_sub_theta*dBHatdzeta-BHat_sub_zeta*dBHatdtheta)/(VPrimeHat*BHat**3)
         ! factorExB = (BHat_sub_theta*dPhi1Hatdzeta-BHat_sub_zeta*dPhi1Hatdtheta)/(VPrimeHat*BHat**2)
         do ispecies = minSpecies,maxSpecies
            THat = THats(ispecies)
            mHat = mHats(ispecies)
            Z = Zs(ispecies)
            sqrtTHat = sqrt(THats(ispecies))
            sqrtmHat = sqrt(mHats(ispecies))

            xIntegralFactor = x*x*x*x*THat*THat*sqrtTHat*pi*Delta/(mHat*sqrtmHat*Z)*ddrN2ddpsiHat
            xIntegralFactorExB = x*x*THat*sqrtTHat*2*pi*Delta*alpha/(mHat*sqrtmHat)*ddrN2ddpsiHat

            do itheta=1,Ntheta
               do izeta=1,Nzeta
                  angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                  cos_angle = cos(angle)
                  sin_angle = sin(angle)
                  select case (whichLambda)
                  case (1) ! BHat cos
                     dBHatdThetadLambda = -m*sin_angle
                     dBHatdZetadLambda = n*Nperiods*sin_angle
                     dBHatdLambda = cos_angle
                     dfactordLambda = (IHat*dBHatdzetadLambda-GHat*dBHatdthetadLambda)/(VPrimeHat*BHat(itheta,izeta)**3) &
                          - 3*(IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta))/(VprimeHat*BHat(itheta,izeta)**4) &
                          * dBHatdLambda - (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                          /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
                     dfactorExBdLambda = - 2*(IHat*dPhi1Hatdzeta(itheta,izeta)-GHat*dPhi1Hatdtheta(itheta,izeta))/(VprimeHat*BHat(itheta,izeta)**3) &
                          * dBHatdLambda - (IHat*dPhi1Hatdzeta(itheta,izeta)-GHat*dPhi1Hatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                          /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**2)
                  case (2) ! IHat
                     dfactordLambda = dBHatdzeta(itheta,izeta)/(VPrimeHat*BHat(itheta,izeta)**3) &
                          - (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta)) &
                          * dVPrimeHatdLambda/(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
                     dfactorExBdLambda = dPhi1Hatdzeta(itheta,izeta)/(VPrimeHat*BHat(itheta,izeta)**2) &
                          - (IHat*dPhi1Hatdzeta(itheta,izeta)-GHat*dPhi1Hatdtheta(itheta,izeta)) &
                          * dVPrimeHatdLambda/(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**2)
                  case (3) ! GHat
                     dfactordLambda = -dBHatdtheta(itheta,izeta)/(VPrimeHat*BHat(itheta,izeta)**3) &
                          - (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                          / (VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
                     dfactorExBdLambda = -dPhi1Hatdtheta(itheta,izeta)/(VPrimeHat*BHat(itheta,izeta)**2) &
                          - (IHat*dPhi1Hatdzeta(itheta,izeta)-GHat*dPhi1Hatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                          /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**2)
                  case (4) ! iota
                     dfactordLambda = - (IHat*dBHatdzeta(itheta,izeta)-GHat*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                          / (VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
                     dfactorExBdLambda = - (IHat*dPhi1Hatdzeta(itheta,izeta)-GHat*dPhi1Hatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                          / (VPrimeHat*VPrimeHat*BHat(itheta,izeta)**2)
                  case (5) ! Phi1Hat cos
                     dPhi1HatdLambda = cos_angle
                     dPhi1HatdThetadLambda = -m*sin_angle
                     dPhi1HatdZetadLambda = n*Nperiods*sin_angle
                     factor = (IHat * dBHatdzeta(itheta,izeta) - GHat * dBHatdtheta(itheta,izeta))/(VPrimeHat*BHat(itheta,izeta)**3)
                     dfactordLambda = 0
                     factorExB = (IHat*dPhi1HatdZeta(itheta,izeta)-GHat*dPhi1HatdTheta(itheta,izeta))/(VPrimeHat*BHat(itheta,izeta)**2)
                     dfactorExBdLambda = (IHat*dPhi1HatdZetadLambda-GHat*dPhi1HatdThetadLambda)/(VPrimeHat * BHat(itheta,izeta)**2)
                  case (6) ! Phi1Hat sin
                     dPhi1HatdLambda = sin_angle
                     dPhi1HatdThetadLambda = m*cos_angle
                     dPhi1HatdZetadLambda = -n*Nperiods * cos_angle
                     factor = (IHat * dBHatdzeta(itheta,izeta) - GHat * dBHatdtheta(itheta,izeta))/(VPrimeHat*BHat(itheta,izeta)**3)
                     dfactordLambda = 0
                     factorExB = (IHat*dPhi1HatdZeta(itheta,izeta)-GHat*dPhi1HatdTheta(itheta,izeta))/(VPrimeHat*BHat(itheta,izeta)**2)
                     dfactorExBdLambda = (IHat*dPhi1HatdZetadLambda-GHat*dPhi1HatdThetadLambda)/(VPrimeHat*BHat(itheta,izeta)**2)
                  end select

                  ! Summed quantity is weighted by charge
                  if (whichSpecies == 0) then
                     factor = factor * Z
                     dfactordLambda = dfactordLambda * Z
                     factorExB = factorExB * Z
                     dfactorExBdLambda = dfactorExBdLambda * Z
                  end if


                  do ix=1,Nx

                     L = 0
                     index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                     ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                     
                     ! Add magnetic drift contribution from deltaF
                     if (.false.) then
                        result = result + &
                             (8/three)*dfactordLambda*xWeights(ix)*deltaFArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)
                     end if
                  

                     ! Add ExB drift contribution from deltaF
                     if (.false.) then
                        result = result + &
                             dfactorExBdLambda*xWeights(ix)*deltaFArray(index)*xIntegralFactorExB(ix)*thetaWeights(itheta)*zetaWeights(izeta) ! 
                     end if

                     if (.true.) then
                     ! Add magnetic drift contribution from f0
                     if((whichLambda .eq. 5) .or. (whichLambda .eq. 6)) then
                        df0dLambda = -(Z*alpha/THat) *  f0Array(index) * dPhi1HatdLambda                 
                        result = result + &
                             (8/three)*factor*xWeights(ix)*df0dLambda*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)
                        results2(itheta,izeta) =  results2(itheta,izeta) + &
                             (8/three)*factor*xWeights(ix)*df0dLambda*xIntegralFactor(ix)
                     else
                        result = result + &
                             (8/three)*dfactordLambda*xWeights(ix)*f0Array(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)
                     end if
                     end if
                     
                     
                     ! Add ExB drift contribution from f0
                     if (.false.) then
                     
                     if((whichLambda .eq. 5) .or. (whichLambda .eq. 6)) then
                        result = result + &
                             xWeights(ix)*xIntegralFactorExB(ix)*thetaWeights(itheta)*zetaWeights(izeta)* (dfactorExBdLambda * f0Array(index) + factorExB * df0dLambda)
                     else
                        result = result + &
                             dfactorExBdLambda*xWeights(ix)*f0Array(index)*xIntegralFactorExB(ix)*thetaWeights(itheta)*zetaWeights(izeta)  ! *2
                     end if
                     end if
                  
                     
                     L = 2
                     index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                     ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                     ! Add the magnetic drift contribution from DeltaF
                     if (.false.) then
                        result = result + &
                             (four/15)*dfactordLambda*xWeights(ix)*deltaFArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)
                     end if
                  
                     ! Add the magnetic drift contribution from f0 for L=2
                     ! This is zero if f0 is Maxwellian
                     ! (which is always the case 2020-06)
                     ! but it is added for feature parity with 'diagnotics'
                     if (.true.) then
                        if((whichLambda .eq. 5) .or. (whichLambda .eq. 6)) then
                           df0dLambda = -(Z*alpha/THat) *  f0Array(index) * dPhi1HatdLambda
                           result = result + &
                                (four/15)*factor*xWeights(ix)*df0dLambda*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)
                        else
                           result = result + &
                                (four/15)*dfactordLambda*xWeights(ix)*f0Array(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)
                        end if
                     
                     end if
                  

                  end do ! ix
               end do ! izeta
            end do ! itheta
            
            do izeta=1,Nzeta
               result3 = result3 + zetaWeights(izeta) &
                  * dot_product(thetaWeights, results2(:,izeta))
            end do
            
            if (masterProc) then
               print *, "!!!!!!!!!!!!!!!!!"
               print *, result
               print *, result3
               print *, "!!!!!!!!!!!!!!!!!"
            end if
            

         end do ! ispecies
         call VecDestroy(deltaFOnProc0,ierr)
         call VecDestroy(f0OnProc0,ierr)
         ! TODO: maybe destroy f0onProc0??
      end if !masterProc
      call VecScatterDestroy(VecScatterContext,ierr)
      call VecScatterDestroy(VecScatterContext0,ierr)
    end subroutine particleFlux_vd_Sensitivity


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
    PetscScalar :: dBHatdLambda, dFSABHat2dLambda, nHat, sqrtFSAB2
    PetscScalar :: angle, cos_angle, sin_angle
    integer :: m, n
    VecScatter :: VecScatterContext
    PetscErrorCode :: ierr
    PetscScalar :: dDHatdLambda

    sqrtFSAB2 = sqrt(FSABHat2)
    m = ms_sensitivity(whichMode)
    n = ns_sensitivity(whichMode)

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
      dFSABHat2dLambda = zero
            dVPrimeHatdLambda = zero
            if (whichLambda==1) then ! BHat
                do itheta=1,Ntheta
                    do izeta=1,Nzeta
                        angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
                        cos_angle = cos(angle)
                        dVPrimeHatdLambda = dVPrimeHatdLambda - (two*(GHat+iota*IHat)) * thetaWeights(itheta) * zetaWeights(izeta) &
                        * BHat(itheta,izeta)**(-3) * cos_angle
                    end do
                end do
                dFSABHat2dLambda = -4*pi*pi*(GHat+iota*IHat)*dVPrimeHatdLambda/(VPrimeHat**2)
            else if (whichLambda==2) then ! IHat
                dVPrimeHatdLambda = VPrimeHat*iota/(GHat+iota*IHat)
            else if (whichLambda==3) then ! GHat
                dVPrimeHatdLambda = VPrimeHat/(GHat+iota*IHat)
            else if (whichLambda==4) then ! iota
                dVPrimeHatdLambda = VPrimeHat*IHat/(GHat+iota*IHat)
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

                !factor = BHat/(DHat*VPrimeHat*sqrtFSAB2) = (GHat+iota*IHat)/(BHat*VPrimeHat*sqrtFSAB2)

        do itheta=1,Ntheta
          do izeta=1,Nzeta
            angle = m * theta(itheta) - n * NPeriods * zeta(izeta)
            cos_angle = cos(angle)
            select case (whichLambda)
              case (1) ! BHat
                dBHatdLambda = cos_angle
                factor = -(GHat+iota*IHat)*dBHatdLambda/(BHat(itheta,izeta)*BHat(itheta,izeta)*VPrimeHat*sqrtFSAB2) &
                    - (GHat+iota*IHat)*dVPrimeHatdLambda/(BHat(itheta,izeta)*VPrimeHat*VPrimeHat*sqrtFSAB2) &
                    - 0.5*(GHat+iota*IHat)*dFSABHat2dLambda/(BHat(itheta,izeta)*VPrimeHat*sqrtFSAB2*FSABHat2)
              case (2) ! IHat
                factor = iota/(BHat(itheta,izeta)*VPrimeHat*sqrtFSAB2) &
                    - (GHat+iota*IHat)*dVPrimeHatdLambda/ &
                    (BHat(itheta,izeta)*sqrtFSAB2*VPrimeHat*VPrimeHat)
              case (3) ! GHat
                factor = one/(BHat(itheta,izeta)*VPrimeHat*sqrtFSAB2) &
                    - (GHat+iota*IHat)*dVPrimeHatdLambda/(BHat(itheta,izeta)*sqrtFSAB2*VPrimeHat*VPrimeHat)
              case (4) ! iota
                factor = IHat/(BHat(itheta,izeta)*VPrimeHat*sqrtFSAB2) &
                    - (GHat+iota*IHat)*dVPrimeHatdLambda/(BHat(itheta,izeta)*sqrtFSAB2*VPrimeHat*VPrimeHat)
            end select

            do ix=1,Nx

                L = 1
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                result = result + &
                  (four/three)*factor*xWeights(ix)*deltaFArray(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

              end do
            end do
          end do
        end do
        call VecDestroy(deltaFOnProc0,ierr)
      end if !masterProc
      call VecScatterDestroy(VecScatterContext,ierr)

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
      integer :: whichLambda, whichMode, this_NModesAdjoint
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
         if ((whichLambda .eq. 1) .or. (whichLambda .eq. 5) .or. (whichLambda .eq. 6)) then
            this_NmodesAdjoint = NModesAdjoint
         else
            this_NmodesAdjoint = 1
         end if
         do whichMode=1,this_NModesAdjoint
            call VecSet(adjointResidual,zero,ierr)
            ! Call function to perform (dLdlambdaf - dSdlambda), which calls populatedMatrixdLambda and populatedRHSdLambda

            call evaluateAdjointInnerProductFactor(forwardSolution,adjointResidual,whichLambda,whichMode)

            if (RHSMode == 5) then
               call innerProduct(adjointSolutionJr,adjointResidual,innerProductResultEr3,0)
               if (.not. includePhi1) then
                  call particleFluxSensitivity(radialCurrentSensitivity,forwardSolution,0,whichLambda,whichMode)
               else
                  call particleFlux_vd_Sensitivity(radialCurrentSensitivity,forwardSolution,0,whichLambda,whichMode)
               end if
               ErTermToAdd = -(innerProductResultEr1/innerProductResultEr2)*(radialCurrentSensitivity-innerProductResultEr3)
               dPhidPsidLambda(whichLambda,whichMode) = (1/innerProductResultEr2)*(radialCurrentSensitivity-innerProductResultEr3)
            end if

          call innerProduct(adjointSolution,adjointResidual,innerProductResult,0)
          ! Save to appropriate sensitivity array
          if (whichSpecies == 0) then
             select case (whichAdjointRHS)
             case (1) ! Particle flux
                if (.not. includePhi1) then
                   call particleFluxSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                else
                   call particleFlux_vd_Sensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                end if
                
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
                if (.not. includePhi1) then
                   call particleFluxSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                else
                   call particleFlux_vd_Sensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                   
                end if
                dParticleFluxdLambda(whichSpecies,whichLambda,whichMode) = sensitivityResult - innerProductResult
                if ((whichLambda == 1) .and. (whichMode == 1) .and. masterProc) then
                   print *,"sensitivity and inner !!!"
                   print *, sensitivityResult
                   print *, innerProductResult
                end if

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

      if (masterProc) then
         print *,dParticleFluxdLambda(1,1,1)
         print *,dParticleFluxdLambda(2,1,1)
         print *,dParticleFluxdLambda(1,2,1)
         print *,dParticleFluxdLambda(2,2,1)
         print *,dParticleFluxdLambda(1,3,1)
         print *,dParticleFluxdLambda(2,3,1)
         print *,dParticleFluxdLambda(1,4,1)
         print *,dParticleFluxdLambda(2,4,1)
      end if
          
      call VecDestroy(adjointResidual,ierr)
      if (RHSMode==5) then
        call VecDestroy(adjointResidualEr,ierr)
      end if

    end subroutine

end module

