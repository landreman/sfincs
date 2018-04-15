#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#include <finclude/petscsnesdef.h>
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petscsnesdef.h>
#endif

subroutine adjointErrorCorrection(forwardSolution, adjointSolution,whichAdjointRHS,whichSpecies)

  use globalVariables
  use petscsnes
  use petscvec
  use adjointDiagnostics
  use indices

  Vec :: forwardSolution, adjointSolution
  integer :: whichAdjointRHS, whichSpecies
  Mat :: matrix_fine
  Vec :: RHS_fine, RHS
  integer :: userContext(1)
  SNES :: mysnes
  Vec :: dummy, residual, adjoint_interp, forward_interp
  PetscScalar :: innerProductResult
  PetscErrorCode :: ierr
  integer :: size, size1, size2
  PetscScalar :: norm
  ! Related to density calculaiton
  integer :: L, ispecies, izeta, itheta, ix, index
  PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, densityFactor, density_interp
  VecScatter :: VecScatterContext
  Vec :: forward_interpOnProc0
  PetscScalar, pointer :: forward_interpArray(:)
  PetscScalar :: VPrimeHat_fine, integral
  PetscScalar, dimension(:,:), allocatable :: zeta_prolongation, DHat_interp

  ! Prolongation matrix is build in createGrids_fine()
  print *,"This is adjointErrorCorrection."
  print *,"whichSpecies = ", whichSpecies
  print *,"whichAdjointRHS = ", whichAdjointRHS

  ! Create forward_interp Vec
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, forward_interp, ierr)
  call VecSet(forward_interp,zero,ierr)
  ! Interpolate forward solution
  call MatMult(adjointECProlongationMatrix, forwardSolution, forward_interp, ierr)
  call VecNorm(forwardSolution,NORM_2,norm,ierr)
  print *,"Norm of forwardSolution: ", norm
  call VecNorm(forward_interp,NORM_2,norm,ierr)
  print *,"Norm of forward_interp: ",norm


!  print *,"DHat_fine: ", DHat_fine

!  allocate(zeta_prolongation(Nzeta_fine,Nzeta))
!  call  periodic_interpolation(Nzeta, Nzeta_fine,2*pi,zeta_fine,zeta_prolongation)
  !! Testing
!  allocate(DHat_interp(Ntheta_fine,Nzeta_fine))
!  DHat_interp = zero
!  !! Compute DHat_fine
!  do itheta = 1,Ntheta_fine
!    do izeta = 1,Nzeta
!      do izeta_fine = 1,Nzeta_fine
!        DHat_interp(itheta,izeta_fine) = DHat_interp(itheta,izeta_fine) + DHat(itheta,izeta)*zeta_prolongation(izeta_fine,izeta)
!      end do
!    end do
!  end do
!  print *,"DHat_interp: ", DHat_interp
!

  !!!!!! Testing interpolation - should provide same density perbation
  ! Scatter deltaF to master proc
!  call VecScatterCreateToZero(forward_interp, VecScatterContext, forward_interpOnProc0, ierr)
!  call VecScatterBegin(VecScatterContext, forward_interp, forward_interpOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
!  call VecScatterEnd(VecScatterContext, forward_interp, forward_interpOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
!  if (masterProc) then
!    ! Convert the PETSc vector into a normal Fortran array
!    call VecGetArrayF90(forward_interpOnProc0, forward_interpArray, ierr)
!
!    do ispecies = 1,Nspecies
!      density_interp = zero
!      THat = THats(ispecies)
!      mHat = mHats(ispecies)
!      sqrtTHat = sqrt(THat)
!      sqrtMHat = sqrt(mHat)
!      densityFactor = 4*pi*THat*sqrtTHat/(mHat*sqrtMHat)
!      do itheta=1, Ntheta_fine
!        do izeta=1, Nzeta_fine
!          do ix=1, Nx
!            L = 0
!            index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F,7)+1
!            density_interp = density_interp &
!            + densityFactor*xWeights(ix)*x2(ix)*forward_interpArray(index)*thetaWeights_fine(itheta)*zetaWeights_fine(izeta)/DHat_fine(itheta,izeta)
!          end do ! ix
!        end do ! izeta
!      end do ! itheta
!      density_interp = density_interp / VPrimeHat
!      print *,"ispecies = ", ispecies
!      print *,"densityPerturbation = ", density_interp
!    end do
!  end if

!  VPrimeHat_fine = zero
!  integral = zero
!  do itheta = 1,Ntheta_fine
!    do izeta = 1,Nzeta_fine
!      VPrimeHat_fine = VPrimeHat_fine + thetaWeights_fine(itheta)*zetaWeights_fine(izeta)/DHat_fine(itheta,izeta)
!      integral = integral + thetaWeights_fine(itheta)*zetaWeights_fine(izeta)*BHat_fine(itheta,izeta)*BHat_fine(itheta,izeta)/DHat_fine(itheta,izeta)
!    end do
!  end do
!  integral = integral/VPrimeHat_fine
!  print *,"VPrimeHat_fine: ", VPrimeHat_fine
!  print *,"VPrimeHat: ", VPrimeHat
!  print *,"FSABHat2_fine: ", integral
!  print *,"FSABHat2: ", FSABHat2

  ! Interpolate adjoint solution
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, adjoint_interp,ierr)
  call VecSet(adjoint_interp,zero,ierr)
  call MatMult(adjointECProlongationmatrix, adjointSolution, adjoint_interp, ierr)
  call VecNorm(adjoint_interp,NORM_2,norm,ierr)
  print *,"Norm of adjoint_interp: ",norm
  call VecNorm(adjointSolution,NORM_2,norm,err)
  print *,"Norm of adjoint: ", norm

  ! Construct fine matrix and RHS
  ! Create dummy and set to zero
  print *,"Building fine grid matrix & RHS."
  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, dummy, ierr)
  call VecSet(dummy, zero, ierr)

  ! whichMatrix = 7 -> matrixSize_fine
  call preallocateMatrix(matrix_fine,7) ! matrix_fine is initialized to 0 here
  call populateMatrix(matrix_fine,7,dummy) ! dummy is ignored here

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, RHS_fine, ierr)
  call VecSet(RHS_fine,zero,ierr)
  userContext(1) = 7
  call evaluateResidual(mysnes,dummy,RHS_fine,userContext,ierr)
  call VecNorm(RHS_fine,NORM_2,norm,ierr)
  print *,"Norm of RHS_fine: ", norm

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, RHS, ierr)
  call VecSet(RHS,zero,ierr)
  userContext(1) = 0
  call evaluateResidual(mysnes,dummy,RHS,userContext,ierr)
  call VecNorm(RHS,NORM_2,norm,ierr)
  print *,"Norm of RHS: ", norm

  call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize_fine, residual, ierr)
  call VecSet(residual,zero,ierr)
  ! Evaluate residual
  call MatMultAdd(matrix_fine, forward_interp, RHS_fine, residual, ierr)
  call VecNorm(residual,NORM_2,norm,ierr)
  print *,"Norm of residual: ",norm

  ! Evaluate inner product
  call innerProduct(adjoint_interp,residual,innerProductResult)
  if (masterProc) then
    if (whichSpecies == 0) then
      select case (whichAdjointRHS)
        case (1) ! Particle flux
          radialCurrent_corrected = dot_product(Zs(1:Nspecies), particleFlux_vm_rN)-innerProductResult
        case (2) ! Heat Flux
          totalHeatFlux_corrected = sum(heatFlux_vm_rN)-innerProductResult
        case (3) ! Bootstrap
          bootstrap_corrected = dot_product(Zs(1:Nspecies),FSABVelocityUsingFSADensityOverRootFSAB2)-innerProductResult
      end select
    else
      select case (whichAdjointRHS)
        case (1) ! Particle flux
          print *,"particleFlux_vm_rN: ", particleFlux_vm_rN(ispecies)
          print *,"innerProductResult: ", innerProductResult
          particleFlux_corrected(whichSpecies) = particleFlux_vm_rN(whichSpecies)-innerProductResult
        case (2) ! Heat Flux
          print *,"heatFlux_vm_rN: ", heatFlux_vm_rN(ispecies)
          print *,"innerProductResult: ", innerProductResult
          heatFlux_corrected(whichSpecies) = heatFlux_vm_rN(whichSpecies)-innerProductResult
        case (3) ! Parallel Flow
          print *,"particleFlow_vm_rN: ", particleFlux_vm_rN(ispecies)
          print *,"innerProductResult: ", innerProductResult
          parallelFlow_corrected(whichSpecies) = FSABVelocityUsingFSADensityOverRootFSAB2(whichSpecies)-innerProductResult
      end select
    end if
  end if

  call MatDestroy(matrix_fine,ierr)
  call VecDestroy(RHS_fine,ierr)

end subroutine adjointErrorCorrection
