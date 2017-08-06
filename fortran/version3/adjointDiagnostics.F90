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

    subroutine evaluateAdjointResidual(forwardSolution,adjointResidual,whichLambda,whichMode)

      use globalVariables
      use petscmat

      implicit none

      Vec :: forwardSolution, adjointResidual
      integer :: whichLambda, whichMode
      Mat :: dMatrixdLambda
      Vec :: dRHSdLambda
      PetscErrorCode :: ierr

      !external populatedRHSdLambda, populatedMatrixdLambda

      ! Allocate and populate dMatrixdLambda
      call preallocateMatrix(dMatrixdLambda, 3)
      !call populatedMatrixdLambda(dMatrixdLambda, whichLambda, whichMode)

      ! Allocate and populate dRHSdLambda
      call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, dRHSdLambda, ierr)
      !call populatedRHSdLambda(dRHSdLambda, whichLambda, whichMode)

      ! Multiply dRHSdLambda by -1
      call VecScale(dRHSdLambda, -1.0d+0, ierr)
      ! Performs adjointResidual = dMatrixdLambda*forwardSolution - dRHSdLambda
      call MatMultAdd(dMatrixdLambda, forwardSolution, dRHSdLambda, adjointResidual, ierr)

    end subroutine

    ! Computes inner product associated with free energy norm for deltaF and deltaG
    subroutine innerProduct(deltaF, deltaG, result)

      use globalVariables
      use indices

      implicit none

      PetscScalar, dimension(:) :: deltaF, deltaG
      PetscScalar :: result
      PetscScalar :: innerProductFactor
      PetscScalar, dimension(:,:), allocatable :: sourcesF, sourcesG
      integer :: index, ispecies, itheta, izeta, ix, L
      PetscScalar :: THat, mHat, nHat
      PetscScalar, dimension(:), allocatable :: xIntegralFactor
      PetscScalar, dimension(:,:), allocatable :: thetazetaIntegralFactor

      allocate(sourcesF(Nspecies,2))
      allocate(sourcesG(Nspecies,2))
      allocate(xIntegralFactor(Nx))
      allocate(thetazetaIntegralFactor(Ntheta,Nzeta))

      ! Get sources associated with deltaF and deltaG
      do ispecies = 1,Nspecies
        sourcesF(ispecies,1) = deltaF(getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)+1)
        sourcesF(ispecies,2) = deltaF(getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)+1)
        sourcesG(ispecies,1) = deltaG(getIndex(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT)+1)
        sourcesG(ispecies,2) = deltaG(getIndex(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT)+1)
      end do

      result = 0
      do ispecies = 1,Nspecies
        THat = THats(ispecies)
        mHat = mHats(ispecies)
        nHat = nHats(ispecies)

        xIntegralFactor = x*x*THat*nHat*exp(-x*x)/(mHat*sqrtpi)
        thetazetaIntegralFactor = 1/DHat
        do itheta=1,Ntheta
          do izeta=1,Nzeta
            do ix=1,Nx
              ! The integral over xi turns into a sum over L (with a factor of (2L+1)/2
              do L=0,(Nxi-1)
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
                ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.

                result = result + thetazetaIntegralFactor(itheta,izeta)*xIntegralFactor(ix)*(2/(2*L+1)) &
                      *deltaF(index)*deltaG(index)/VPrimeHat

              end do
            end do
          end do
        end do

        ! Now add terms which sum over sources
        result = result + sourcesF(ispecies,1)*sourcesG(ispecies,1) + sourcesF(ispecies,2)*sourcesG(ispecies,2)

      end do

    end subroutine

    ! This subroutine evaluates the term in the sensitivity derivative of the fluxes 
    ! which arise due to the sensitivity of the inner product and the integrating factor,
    ! not f1 itself, hence this is not called with the adjoint variable
    subroutine heatFluxSensitivityInnerProduct(result, deltaF, whichSpecies, whichLambda, whichMode)

    use globalVariables
    use indices

    implicit none

    PetscScalar :: result
    PetscScalar, dimension(:) :: deltaF
    integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, dBHatdThetadLambda, dBHatdZetadLambda
    PetscScalar :: dBHat_sub_thetadLambda, dBHat_sub_zetadLambda, dinvDHatdLambda, factor
    PetscScalar, dimension(:), allocatable :: xIntegralFactor
    PetscScalar :: dBHatdLambda, dVPrimeHatdLambda

    result = 0

    allocate(xIntegralFactor(Nx))

    if (whichSpecies == 0) then
      minSpecies = 0
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
      xIntegralFactor = x*x*x*x*x*x*pi*THat*THat*THat*sqrtTHat*Delta/(2*mHat*sqrtmHat*Zs(ispecies))

      do itheta=1,Ntheta
        do izeta=1,Nzeta

          select case (whichLambda)
            case (0) ! Er
              if (masterProc) then
                print *,"Error! Er sensitivity not yet implemented."
              end if
              stop
            case (1) ! BHat
              dBHatdThetadLambda = -ms(whichMode)*sin(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              dBHatdZetadLambda = ns(whichMode)*sin(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              dBHatdLambda = cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              factor = (BHat_sub_theta(itheta,izeta)*dBHatdZetadLambda - BHat_sub_zeta(itheta,izeta)*dBHatdThetadLambda)/(VPrimeHat*BHat(itheta,izeta)**3) - 3*(BHat_sub_theta(itheta,izeta)*dBHatdZeta(itheta,izeta) - &
                BHat_sub_zeta(itheta,izeta)*dBHatdTheta(itheta,izeta))*dBHatdLambda/(BHat(itheta,izeta)**4)
            case (2) ! BHat_sup_theta
              factor = 0
            case (3) ! BHat_sup_zeta
              factor = 0
            case (4) ! BHat_sub_theta
              dBHat_sub_thetadLambda = cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              factor = dBHatdzeta(itheta,izeta)*dBHat_sub_thetadLambda/(VPrimeHat*BHat(itheta,izeta)**3)
            case (5) ! BHat_sub_zeta
              dBHat_sub_zetadLambda = cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              factor = -dBHatdtheta(itheta,izeta)*dBHat_sub_zetadLambda/(VPrimeHat*BHat(itheta,izeta)**3)
            case (6) ! DHat
              if (ns(whichMode)==0 .and. ns(whichMode)==0) then
                dVPrimeHatdLambda = 4*pi*pi
                factor = - (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda/ &
                  (VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
              end if
          end select

          do ix=1,Nx

              L = 0
              index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
              ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
              result = result + &
                (8/3)*factor*xWeights(ix)*deltaF(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)*ddrN2ddpsiHat

              L = 2
              index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
              ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
              result = result + &
                (4/15)*factor*xWeights(ix)*deltaF(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)*ddrN2ddpsiHat

            end do
          end do
        end do
      end do

    end subroutine

    ! This subroutine evaluates the term in the sensitivity derivative of the fluxes 
    ! which arise due to the sensitivity of the inner product and the integrating factor,
    ! not f1 itself, hence this is not called with the adjoint variable
    subroutine particleFluxSensitivityInnerProduct(result, deltaF, whichSpecies, whichLambda, whichMode)

    use globalVariables
    use indices

    implicit none

    PetscScalar :: result
    PetscScalar, dimension(:) :: deltaF
    integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, dBHatdThetadLambda, dBHatdZetadLambda
    PetscScalar :: dBHat_sub_thetadLambda, dBHat_sub_zetadLambda, dinvDHatdLambda, factor, dRootFSAB2dLambda
    PetscScalar, dimension(:), allocatable :: xIntegralFactor
    PetscScalar :: dBHatdLambda, dVPrimeHatdLambda

    result = 0

    allocate(xIntegralFactor(Nx))

    if (whichSpecies == 0) then
      minSpecies = 0
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

      xIntegralFactor = x*x*x*x*THat*THat*sqrtTHat*pi*Delta/(mHat*sqrtmHat*Zs(ispecies))

      do itheta=1,Ntheta
        do izeta=1,Nzeta

          select case (whichLambda)
            case (0) ! Er
              if (masterProc) then
                print *,"Error! Er sensitivity not yet implemented."
              end if
              stop
            case (1) ! BHat
              dBHatdThetadLambda = -ms(whichMode)*sin(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              dBHatdZetadLambda = ns(whichMode)*sin(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              dBHatdLambda = cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              factor = (BHat_sub_theta(itheta,izeta)*dBHatdZetadLambda - BHat_sub_zeta(itheta,izeta)*dBHatdThetadLambda)/(VPrimeHat*BHat(itheta,izeta)**3) - 3*(BHat_sub_theta(itheta,izeta)*dBHatdZeta(itheta,izeta) &
                - BHat_sub_zeta(itheta,izeta)*dBHatdTheta(itheta,izeta))*dBHatdLambda/(VPrimeHat*BHat(itheta,izeta)**3)
            case (2) ! BHat_sup_theta
              factor = 0
            case (3) ! BHat_sup_zeta
              factor = 0
            case (4) ! BHat_sub_theta
              dBHat_sub_thetadLambda = cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              factor = dBHatdzeta(itheta,izeta)*dBHat_sub_thetadLambda/(VPrimeHat*BHat(itheta,izeta)**3)
            case (5) ! BHat_sub_zeta
              dBHat_sub_zetadLambda = cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))
              factor = -dBHatdtheta(itheta,izeta)*dBHat_sub_zetadLambda/(VPrimeHat*BHat(itheta,izeta)**3)
            case (6) ! DHat
              if (ms(whichMode)==0 .and. ns(whichMode)== 0) then
                dVPrimeHatdLambda = 4*pi*pi
                factor = - (BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))*dVPrimeHatdLambda &
                    /(VPrimeHat*VPrimeHat*BHat(itheta,izeta)**3)
              end if
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
                (8/3)*factor*xWeights(ix)*deltaF(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)*ddrN2ddpsiHat

              L = 2
              index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
              ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
              result = result + &
                (4/15)*factor*xWeights(ix)*deltaF(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)*ddrN2ddpsiHat

            end do
          end do
        end do
      end do
    end subroutine

    ! This subroutine evaluates the term in the sensitivity derivative of the fluxes 
    ! which arise due to the sensitivity of the inner product and the integrating factor,
    ! not f1 itself, hence this is not called with the adjoint variable
    subroutine bootstrapSensitivityInnerProduct(result, deltaF, whichSpecies, whichLambda, whichMode)

    use globalVariables
    use indices

    implicit none

    PetscScalar :: result
    PetscScalar, dimension(:) :: deltaF
    integer :: whichSpecies, whichLambda, whichMode, minSpecies, maxSpecies, itheta, izeta, L, ix, index, ispecies
    PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, dBHatdThetadLambda, dBHatdZetadLambda
    PetscScalar :: dBHat_sub_thetadLambda, dBHat_sub_zetadLambda, dinvDHatdLambda, factor, dVPrimeHatdLambda
    PetscScalar, dimension(:), allocatable :: xIntegralFactor
    PetscScalar :: dBHatdLambda, dRootFSAB2dLambda, nHat, rootFSAB2

    rootFSAB2 = sqrt(FSABHat2)

    result = 0

   if (whichLambda==1) then
      do itheta=1,Ntheta
        do izeta=1,Nzeta
          dRootFSAB2dLambda = dRootFSAB2dLambda + (thetaWeights(itheta)*zetaWeights(izeta)*BHat(itheta,izeta)*cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta)))/(DHat(itheta,izeta)*VPrimeHat*RootFSAB2)
        end do
      end do
    end if

    if (whichLambda==6) then
      do itheta=1,Ntheta
        do izeta=1,Nzeta
          dRootFSAB2dLambda = dRootFSAB2dLambda + 0.5*thetaWeights(itheta)*zetaWeights(izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)*cos(ms(whichMode)*theta(itheta)-ns(whichMode)*Nperiods*zeta(izeta))/(VPrimeHat*RootFSAB2)
        end do
      end do
    end if

    allocate(xIntegralFactor(Nx))

    do ispecies = 0,Nspecies
      THat = THats(ispecies)
      mHat = mHats(ispecies)
      sqrtTHat = sqrt(THats(ispecies))
      sqrtmHat = sqrt(mHats(ispecies))
      nHat = nHats(ispecies)

      xIntegralFactor = x*x*x*THat*THat*2*pi*Zs(ispecies)/(mHat*mHat*nHat)

      do itheta=1,Ntheta
        do izeta=1,Nzeta

          select case (whichLambda)
            case (0) ! Er
              if (masterProc) then
                print *,"Error! Er sensitivity not yet implemented."
              end if
              stop
            case (1) ! BHat
              dBHatdLambda = cos(ms(whichMode)*theta(itheta)- ns(whichMode)*Nperiods*zeta(izeta))
              factor = dBHatdLambda/(DHat(itheta,izeta)*VPrimeHat*rootFSAB2) - 0.5*BHat(itheta,izeta)*dRootFSAB2dLambda/(DHat(itheta,izeta)*rootFSAB2*FSABHat2*VPrimeHat)
            case (2) ! BHat_sup_theta
              factor = 0
            case (3) ! BHat_sup_zeta
              factor = 0
            case (4) ! BHat_sub_theta
              factor = 0
            case (5) ! BHat_sub_zeta
              factor = 0
            case (6) ! DHat
              dinvDHatdLambda = cos(ms(whichMode)*theta(itheta)- ns(whichMode)*Nperiods*zeta(izeta))
              factor = BHat(itheta,izeta)*dinvDHatdLambda/(VPrimeHat*RootFSAB2) &
                - BHat(itheta,izeta)*dRootFSAB2dLambda/(DHat(itheta,izeta)*VPrimeHat*FSABHat2)
              if (ms(whichMode)==0 .and. ns(whichMode)==0) then
                dVPrimeHatdLambda = (4*pi*pi)
                factor = factor - BHat(itheta,izeta)*dVPrimeHatdLambda/(VPrimeHat*VPrimeHat*RootFSAB2)
              end if
          end select

          do ix=1,Nx

              L = 1
              index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)+1
              ! Add 1 to index to convert from PETSc 0-based index to fortran 1-based index.
              result = result + &
                (4/3)*factor*xWeights(ix)*deltaF(index)*xIntegralFactor(ix)*thetaWeights(itheta)*zetaWeights(izeta)

            end do
          end do
        end do
      end do
    end subroutine

    subroutine evaluateDiagnostics(forwardSolution, adjointSolution, whichAdjointRHS, whichSpecies)

      use globalVariables
      use indices
      use writeHDF5Output
      use petscvec

      implicit none

      PetscErrorCode :: ierr
      Vec :: forwardSolution, adjointSolution, adjointResidual, adjointRHSVec
      VecScatter :: VecScatterContext
      integer :: whichAdjointRHS, whichSpecies
      Vec :: adjointSolutionOnProc0, adjointResidualOnProc0, forwardSolutionOnProc0
      PetscScalar, pointer :: adjointResidualArray(:), adjointSolutionArray(:), forwardSolutionArray(:)
      integer :: whichLambda, whichMode
      PetscScalar :: innerProductResult, sensitivityResult

      ! Allocate appropriate sensitivity arrays and populateAdjointRHS for inner product
      if (whichSpecies == 0) then
        select case (whichAdjointRHS)
          case (1) ! Particle flux
            allocate(dRadialCurrentdLambda(NLambdas,NModesAdjoint))
          case (2) ! Heat Flux
            allocate(dTotalHeatFluxdLambda(NLambdas,NModesAdjoint))
          case (3) ! Bootstrap
            allocate(dBootstrapdLambda(NLambdas,NModesAdjoint))
        end select
      else
        select case (whichAdjointRHS)
          case (1) ! Particle flux
            allocate(dParticleFluxdLambda(NSpecies,NLambdas,NModesAdjoint))
          case (2) ! Heat Flux
            allocate(dHeatFluxdLambda(NSpecies,NLambdas,NModesAdjoint))
        end select
      end if

      if (masterProc) then
        print *,"Computing adjoint diagnostics for RHS ", whichAdjointRHS, " and species ", whichSpecies
      end if

      ! Create a scattering context for adjointResidual and adjointSolution
      call VecScatterCreateToZero(adjointSolution, VecScatterContext, adjointSolutionOnProc0, ierr)
      ! Create adjointResidualOnProc0 of same type
      call VecDuplicate(adjointSolutionOnProc0, adjointResidualOnProc0, ierr)
      ! Create forwardSolution of same type
      call VecDuplicate(adjointSolutionOnProc0, forwardSolutionOnProc0, ierr)
      ! Send the adjointSolution to the master process:
      call VecScatterBegin(VecScatterContext, adjointSolution, adjointSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(VecScatterContext, adjointSolution, adjointSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
      ! Send forwardSolution to master proc
      call VecScatterBegin(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(VecScatterContext, forwardSolution, forwardSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
      if (masterProc) then
        ! Convert the PETSc vector into a normal Fortran array
        call VecGetArrayF90(adjointSolutionOnProc0, adjointSolutionArray, ierr)
        call VecGetARrayF90(forwardSolutionOnProc0, forwardSolutionArray, ierr)
      end if

      ! Same inner product is formed regardless of whichAdjointMatrix and whichSpecies
      ! Loop over lambda's and perform inner product
      ! rethink this for Er
      do whichLambda=1,6
        do whichMode=1,NModesAdjoint

          ! Call function to perform (dLdlambdaf - dSdlambda), which calls populatedMatrixdLambda and populatedRHSdLambda
          call evaluateAdjointResidual(forwardSolution,adjointResidual,whichLambda,whichMode)

          ! Send the adjointResidual to the master proc
          call VecScatterBegin(VecScatterContext, adjointResidual, adjointResidualOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call VecScatterEnd(VecScatterContext, adjointResidual, adjointResidualOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)

          if (masterProc) then
            ! Convert the PETSc vector into a normal Fortran array
            call VecGetArrayF90(adjointResidualOnProc0, adjointResidualArray, ierr)

            call innerProduct(adjointSolutionArray,adjointResidualArray,innerProductResult)

            ! Save to appropriate sensitivity array
            if (whichSpecies == 0) then
              select case (whichAdjointRHS)
                case (1) ! Particle flux
                  call particleFluxSensitivityInnerProduct(sensitivityResult, forwardSolutionArray, whichSpecies, whichLambda, whichMode)
                  dRadialCurrentdLambda(whichLambda,whichMode) = innerProductResult + sensitivityResult
                case (2) ! Heat Flux
                  call heatFluxSensitivityInnerProduct(sensitivityResult, forwardSolutionArray, whichSpecies, whichLambda, whichMode)
                  dTotalHeatFluxdLambda(whichLambda,whichMode) = innerProductResult + sensitivityResult
                case (3) ! Bootstrap
                  call bootstrapSensitivityInnerProduct(sensitivityResult, forwardSolutionArray, whichSpecies, whichLambda, whichMode)
                  dBootstrapdLambda(whichLambda,whichMode) = innerProductResult + sensitivityResult
              end select
            else
              select case (whichAdjointRHS)
                case (1) ! Particle flux
                  call particleFluxSensitivityInnerProduct(sensitivityResult, forwardSolutionArray, whichSpecies, whichLambda, whichMode)
                  dParticleFluxdLambda(whichSpecies,whichLambda,whichMode) = innerProductResult
                case (2) ! Heat Flux
                  call heatFluxSensitivityInnerProduct(sensitivityResult, forwardSolutionArray, whichSpecies, whichLambda, whichMode)
                  dHeatFluxdLambda(whichSpecies,whichLambda,whichMode) = innerProductResult
              end select
            end if

          end if

        end do
      end do

    end subroutine

end module

