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
    use diagnostics

    implicit none

    PetscScalar :: result
    Vec :: forwardSolution, forwardSolutionOnProc0
    PetscScalar, pointer :: forwardSolutionArray(:)
    integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
    VecScatter :: VecScatterContext
    PetscErrorCode :: ierr

    PetscScalar, dimension(:), allocatable :: R

    result = zero
    
    ! Scatter deltaF to master proc
    call VecScatterCreateToZero(forwardSolution, VecScatterContext, forwardSolutionOnProc0, ierr)
    call VecScatterBegin(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    if (masterProc) then
      ! Convert the PETSc vector into a normal Fortran array
      call VecGetArrayF90(forwardSolutionOnProc0, forwardSolutionArray, ierr)
    end if

    if (masterProc) then
       allocate(R(matrixSize))

       if (whichSpecies == 0) then
          minSpecies = 1
          maxSpecies = NSpecies
       else
          minSpecies = whichSpecies
          maxSpecies = whichSpecies
       end if

       do ispecies = minSpecies,maxSpecies
          call heatFlux_vm(R, 1, ispecies, whichLambda, whichMode)
          R = R * ddrN2ddpsiHat
          do itheta=1,Ntheta
             do izeta=1,Nzeta
                do ix=1,Nx

                   L = 0
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                   result = result + R(index)*forwardSolutionArray(index)

                   L = 2
                   index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                   ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                   result = result + R(index)*forwardSolutionArray(index)

                end do
             end do
          end do
       end do
       call VecDestroy(forwardSolutionOnProc0,ierr)
       deallocate(R)
    end if !masterProc
    call VecScatterDestroy(VecScatterContext,ierr)

    end subroutine

    subroutine heatFlux_vE_Sensitivity(result, deltaF, whichSpecies, whichLambda, whichMode)

      use globalVariables
      use indices
      use petscvec
      use diagnostics
      
      implicit none

      PetscScalar :: result
      Vec :: deltaF, deltaFOnProc0, f0OnProc0
      PetscScalar, pointer :: deltaFArray(:), f0Array(:)
      integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
      VecScatter :: VecScatterContext, VecScatterContext0 !MAYBE REMOVE SECOND
      PetscErrorCode :: ierr

      PetscScalar, dimension(:), allocatable :: R
      
      result = zero

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

      if (masterProc) then
         allocate(R(matrixSize))
    
         
         if (whichSpecies == 0) then
            minSpecies = 1
            maxSpecies = NSpecies
         else
            minSpecies = whichSpecies
            maxSpecies = whichSpecies
         end if

         do ispecies = minSpecies,maxSpecies

            ! most of the work is now done here
            call heatFlux_vE(R, 1, ispecies, whichLambda, whichMode)
            R = R * ddrN2ddpsiHat

            ! Summed quantity is weighted by charge
            if (whichSpecies == 0) then
               R = R * Zs(ispecies)
            end if
            
            do itheta=1,Ntheta
               do izeta=1,Nzeta
                  do ix=1,Nx

                     L = 0
                     index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                     ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                     result = result + R(index)*deltaFArray(index)

                  end do ! ix
               end do ! izeta
            end do ! itheta
         end do ! ispecies
         call VecDestroy(deltaFOnProc0,ierr)
         call VecDestroy(f0OnProc0,ierr)
         deallocate(R)
      end if !masterProc
      call VecScatterDestroy(VecScatterContext,ierr)
      call VecScatterDestroy(VecScatterContext0,ierr)
    end subroutine heatFlux_vE_Sensitivity

    
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
    use diagnostics

    implicit none

    PetscScalar :: result
    Vec :: deltaF, deltaFOnProc0
    PetscScalar, pointer :: deltaFArray(:)
    integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
    VecScatter :: VecScatterContext
    PetscErrorCode :: ierr

    PetscScalar, dimension(:), allocatable :: R
    
    result = zero

    ! Scatter deltaF to master proc
    call VecScatterCreateToZero(deltaF, VecScatterContext, deltaFOnProc0, ierr)
    call VecScatterBegin(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, deltaF, deltaFOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    if (masterProc) then
      ! Convert the PETSc vector into a normal Fortran array
      call VecGetArrayF90(deltaFOnProc0, deltaFArray, ierr)
    end if

    if (masterProc) then
      allocate(R(matrixSize))
      
      if (whichSpecies == 0) then
        minSpecies = 1
        maxSpecies = NSpecies
      else
        minSpecies = whichSpecies
        maxSpecies = whichSpecies
      end if

      do ispecies = minSpecies,maxSpecies

         ! most of the work is now done here
         call particleFlux_vm(R, 1, ispecies, whichLambda, whichMode)
         R = R * ddrN2ddpsiHat
         
         ! Summed quantity is weighted by charge
         if (whichSpecies == 0) then
            R = R * Zs(ispecies)
         end if

         do itheta=1,Ntheta
            do izeta=1,Nzeta
               do ix=1,Nx
                  
                  L = 0
                  index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                  ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                  result = result + R(index)*deltaFArray(index)

                  L = 2
                  index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                  ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                  result = result + R(index)*deltaFArray(index)
                  
              end do ! ix
            end do ! izeta
          end do ! itheta
        end do ! ispecies
        call VecDestroy(deltaFOnProc0,ierr)
        deallocate(R)

      end if !masterProc
      call VecScatterDestroy(VecScatterContext,ierr)
      
    end subroutine particleFluxSensitivity

    !> Evaluates the term in the sensitivity derivative of the particle fluxes with Phi1
    !! which arise due to the sensitivity of the inner product and the integrating factor,
    !! not f1 itself, hence this is not called with the adjoint variable
    !! NOTE: not thoroughly tested,
    !! do not use in production without comparing to finite difference sensitivty.
    !! @param result Output of sensitivity calculation.
    !! @param deltaF \f$\hat{F}\f$ solution from forward equation.
    !! @param whichSpecies If = 0, summed over species. If nonzero, indicates species number. In this case should always be 0.
    !! @param whichLambda Indicates which component of magnetic field derivative is respect to. If = 0 \f$E_r\f$, = 1 \f$\hat{B}\f$, = 2 \f$\hat{B}^{\theta}\f$, = 3 \f$\hat{B}^{\zeta}\f$, = 4 \f$\hat{B}_{\theta}\f$, = 5 \f$\hat{B}_{\zeta}\f$, = 6 \f$\hat{D}\f$
    !! @param whichMode Indicates index of ms and ns for derivative.
    subroutine particleFlux_vE_Sensitivity(result, deltaF, whichSpecies, whichLambda, whichMode)

      use globalVariables
      use indices
      use petscvec
      use diagnostics
      
      implicit none

      PetscScalar :: result
      Vec :: deltaF, deltaFOnProc0, f0OnProc0
      PetscScalar, pointer :: deltaFArray(:), f0Array(:)
      integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
      VecScatter :: VecScatterContext, VecScatterContext0 !MAYBE REMOVE SECOND
      PetscErrorCode :: ierr

      PetscScalar, dimension(:), allocatable :: R
      
      result = zero
      
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

      if (masterProc) then
         allocate(R(matrixSize))
    
         
         if (whichSpecies == 0) then
            minSpecies = 1
            maxSpecies = NSpecies
         else
            minSpecies = whichSpecies
            maxSpecies = whichSpecies
         end if

         do ispecies = minSpecies,maxSpecies

            ! most of the work is now done here
            call particleFlux_vE(R, 1, ispecies, whichLambda, whichMode)
            R = R * ddrN2ddpsiHat

            ! Summed quantity is weighted by charge
            if (whichSpecies == 0) then
               R = R * Zs(ispecies)
            end if
            
            do itheta=1,Ntheta
               do izeta=1,Nzeta
                  do ix=1,Nx

                     L = 0
                     index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                     ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                     result = result + R(index)*deltaFArray(index)

                  end do ! ix
               end do ! izeta
            end do ! itheta
         end do ! ispecies
         call VecDestroy(deltaFOnProc0,ierr)
         call VecDestroy(f0OnProc0,ierr)
         deallocate(R)
      end if !masterProc
      call VecScatterDestroy(VecScatterContext,ierr)
      call VecScatterDestroy(VecScatterContext0,ierr)
    end subroutine particleFlux_vE_Sensitivity


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
    use diagnostics

    implicit none

    PetscScalar :: result
    Vec :: deltaF
    integer :: whichSpecies, whichLambda, whichMode
    Vec :: deltaFOnProc0
    PetscScalar, pointer :: deltaFArray(:)
    integer :: minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
    integer :: m, n
    VecScatter :: VecScatterContext
    PetscErrorCode :: ierr
 
    PetscScalar, dimension(:), allocatable :: R

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

    if (masterProc) then
       allocate(R(matrixSize))
       
       if (whichSpecies == 0) then
          minSpecies = 1
          maxSpecies = NSpecies
       else
          minSpecies = whichSpecies
          maxSpecies = whichSpecies
       end if


      do ispecies = minSpecies,maxSpecies
         ! most of the work is now done here
         call parallelFlow(R, 1, ispecies, whichLambda, whichMode)

         ! Summed quantity is weighted by charge
         if (whichSpecies == 0) then
            R = R * Zs(ispecies)
         end if
         
         do itheta=1,Ntheta
            do izeta=1,Nzeta
               do ix=1,Nx
                  L = 1
                  index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                  ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
                  result = result + R(index)*deltaFArray(index)
              end do ! ix
            end do ! izeta
          end do ! itheta
        end do ! ispecies
        call VecDestroy(deltaFOnProc0,ierr)
        deallocate(R)
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
      Vec :: forwardSolution, adjointSolution, adjointResidual
      integer :: whichAdjointRHS, whichSpecies
      integer :: whichLambda, whichMode, this_NModesAdjoint
      PetscScalar :: innerProductResult, sensitivityResult, sensitivityResult2, ErTermToAdd
      Vec :: adjointResidualEr, adjointSolutionJr
      PetscScalar :: innerProductResultEr1, innerProductResultEr2, innerProductResultEr3
      PetscScalar :: radialCurrentSensitivity, radialCurrentSensitivity2

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
         if ((whichLambda .eq. 1)) then
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
               call particleFluxSensitivity(radialCurrentSensitivity,forwardSolution,0,whichLambda,whichMode)
               if (includePhi1) then
                  call particleFlux_vE_Sensitivity(radialCurrentSensitivity2,forwardSolution,0,whichLambda,whichMode)
                  radialCurrentSensitivity = radialCurrentSensitivity + radialCurrentSensitivity2
               end if
               ErTermToAdd = -(innerProductResultEr1/innerProductResultEr2)*(radialCurrentSensitivity-innerProductResultEr3)
               dPhidPsidLambda(whichLambda,whichMode) = (1/innerProductResultEr2)*(radialCurrentSensitivity-innerProductResultEr3)
            end if
         
            call innerProduct(adjointSolution,adjointResidual,innerProductResult,0)
      
          ! Save to appropriate sensitivity array
          if (whichSpecies == 0) then
             select case (whichAdjointRHS)
             case (1) ! Particle flux
                call particleFluxSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                if (includePhi1) then
                   call particleFlux_vE_Sensitivity(sensitivityResult2, forwardSolution, whichSpecies, whichLambda, whichMode)
                   sensitivityResult = sensitivityResult + sensitivityResult2
                end if
                
                dRadialCurrentdLambda(whichLambda,whichMode) = -innerProductResult + sensitivityResult
             case (2) ! Heat Flux
                call heatFluxSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                if (includePhi1) then
                   call heatFlux_vE_Sensitivity(sensitivityResult2, forwardSolution, whichSpecies, whichLambda, whichMode)
                   sensitivityResult = sensitivityResult + sensitivityResult2
                end if
                
                
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
                if (includePhi1) then
                   call particleFlux_vE_Sensitivity(sensitivityResult2, forwardSolution, whichSpecies, whichLambda, whichMode)
                   sensitivityResult = sensitivityResult + sensitivityResult2
                end if
                dParticleFluxdLambda(whichSpecies,whichLambda,whichMode) = sensitivityResult - innerProductResult
                
                if (RHSMode == 5) then
                   dParticleFluxdLambda(whichSpecies,whichLambda,whichMode) = dParticleFluxdLambda(whichSpecies,whichLambda,whichMode) + ErTermToAdd
                end if
             case (2) ! Heat Flux
                call heatFluxSensitivity(sensitivityResult, forwardSolution, whichSpecies, whichLambda, whichMode)
                if (includePhi1) then
                   call heatFlux_vE_Sensitivity(sensitivityResult2, forwardSolution, whichSpecies, whichLambda, whichMode)
                   sensitivityResult = sensitivityResult + sensitivityResult2
                end if
                
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

