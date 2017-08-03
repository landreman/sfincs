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

      !PetscScalar, pointer :: deltaF(:), deltaG(:)
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

                result = result + thetazetaIntegralFactor(itheta,izeta)*xIntegralFactor(ix)*((2*L+1)/2) &
                      *deltaF(index)*deltaG(index)/VPrimeHat

              end do
            end do
          end do
        end do

        ! Now add terms which sum over sources
        result = result + sourcesF(ispecies,1)*sourcesG(ispecies,1) + sourcesF(ispecies,2)*sourcesG(ispecies,2)

      end do

    end subroutine


    subroutine evaluateDiagnostics(forwardSolution, adjointSolution, whichAdjointRHS, whichSpecies)

      use globalVariables
      use indices
      use writeHDF5Output
      use petscvec

      implicit none

      PetscErrorCode :: ierr
      Vec :: forwardSolution, adjointSolution, adjointResidual
      VecScatter :: VecScatterContext
      integer :: whichAdjointRHS, whichSpecies
      Vec :: adjointSolutionOnProc0, adjointResidualOnProc0
      PetscScalar, pointer :: adjointResidualArray(:), adjointSolutionArray(:)
      integer :: whichLambda, whichMode
      PetscScalar :: result

      ! Allocate appropriate sensitivity array
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
      ! Send the adjointSolution to the master process:
      call VecScatterBegin(VecScatterContext, adjointSolution, adjointSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(VecScatterContext, adjointSolution, adjointSolutionOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
      if (masterProc) then
        ! Convert the PETSc vector into a normal Fortran array
        call VecGetArrayF90(adjointSolutionOnProc0, adjointSolutionArray, ierr)
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

            call innerProduct(adjointSolutionArray,adjointResidualArray,result)

            ! Save to appropriate sensitivity array
            if (whichSpecies == 0) then
              select case (whichAdjointRHS)
                case (1) ! Particle flux
                  dRadialCurrentdLambda(whichLambda,whichMode) = result
                case (2) ! Heat Flux
                  dTotalHeatFluxdLambda(whichLambda,whichMode) = result
                case (3) ! Bootstrap
                  dBootstrapdLambda(whichLambda,whichMode) = result
              end select
            else
              select case (whichAdjointRHS)
                case (1) ! Particle flux
                  dParticleFluxdLambda(whichSpecies,whichLambda,whichMode) = result
                case (2) ! Heat Flux
                  dHeatFluxdLambda(whichSpecies,whichLambda,whichMode) = result
              end select
            end if

          end if

        end do
      end do

    end subroutine

end module

