#include <finclude/petscdmdadef.h>
#include "PETScVersions.F90"

  subroutine createGrids()

    use globalVariables
    use polynomialDiffMatrices
    use xGrid
    use petscdmda
    use geometry
    use indices

    implicit none

    PetscErrorCode :: ierr
    integer :: i, j, itheta, izeta, scheme

    PetscScalar, dimension(:,:), allocatable :: d2dtheta2, d2dzeta2
    PetscScalar, dimension(:), allocatable :: theta_preconditioner, thetaWeights_preconditioner
    PetscScalar, dimension(:,:), allocatable :: d2dtheta2_preconditioner
    PetscScalar, dimension(:), allocatable :: zeta_preconditioner, zetaWeights_preconditioner
    PetscScalar, dimension(:,:), allocatable :: d2dzeta2_preconditioner
    PetscScalar, dimension(:), allocatable :: xWeightsPotentials

    DM :: myDM
    character(len=200) :: procAssignments

    ! *******************************************************************************
    ! Do a few sundry initialization tasks:
    ! *******************************************************************************

    if (constraintScheme < 0) then
       if (collisionOperator == 0) then
          constraintScheme = 1
       else
          constraintScheme = 2
       end if
    end if

    if (forceOddNthetaAndNzeta) then
       if (mod(Ntheta, 2) == 0) then
          Ntheta = Ntheta + 1
       end if
       if (mod(Nzeta, 2) == 0) then
          Nzeta = Nzeta + 1
       end if
    end if

    if (masterProc) then
       print *,"---- Numerical parameters: ----"
       print *,"Ntheta             = ", Ntheta
       print *,"Nzeta              = ", Nzeta
       print *,"Nxi                = ", Nxi
       print *,"NL                 = ", NL
       print *,"Nx                 = ", Nx
       print *,"NxPotentialsPerVth = ", NxPotentialsPerVth
       print *,"xMax               = ",xMax
       print *,"solverTolerance    = ",solverTolerance
       select case (thetaDerivativeScheme)
       case (0)
          print *,"Theta derivative: spectral collocation"
       case (1)
          print *,"Theta derivative: centered finite differences, 3-point stencil"
       case (2)
          print *,"Theta derivative: centered finite differences, 5-point stencil"
       case default
          print *,"Error! Invalid setting for thetaDerivativeScheme"
          stop
       end select
       select case (zetaDerivativeScheme)
       case (0)
          print *,"Zeta derivative: spectral collocation"
       case (1)
          print *,"Zeta derivative: centered finite differences, 3-point stencil"
       case (2)
          print *,"Zeta derivative: centered finite differences, 5-point stencil"
       case default
          print *,"Error! Invalid setting for zetaDerivativeScheme"
          stop
       end select
       if (useIterativeSolver) then
          print *,"Using iterative solver"
       else
          print *,"Using direct solver"
       end if
    end if

    matrixSize = Nspecies * Ntheta * Nzeta * Nxi * Nx
    !matrixSize = Nspecies * Ntheta * Nzeta * Nxi * Nx + Ntheta * Nzeta + 1 ! Includes QN
    select case (constraintScheme)
    case (0)
    case (1)
       matrixSize = matrixSize + 2 * Nspecies
    case (2)
       matrixSize = matrixSize + Nx * Nspecies
    case default
       print *,"Error! Invalid constraintScheme"
       stop
    end select

    if (masterProc) then
       print *,"The matrix is ",matrixSize,"x",matrixSize," elements."
    end if

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Set up ranges of indices owned by each processor.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    ! Each processor is responsible for building the rows of the matrix and rhs corresponding
    ! to its ithetaMin:ithetaMax and izetaMin:izetaMax, and each processor is resposible for all columns of the matrix.

    ! Assign a range of zeta indices to each processor.
    ! This is done by creating a PETSc DM that is not actually used for anything else.
    call DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, &
         Ntheta, Nzeta, PETSC_DECIDE, PETSC_DECIDE, 1, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, myDM, ierr)

    call DMDAGetCorners(myDM, ithetaMin, izetaMin, PETSC_NULL_INTEGER, &
         localNtheta, localNzeta, PETSC_NULL_INTEGER, ierr)

    call DMDestroy(myDM, ierr)

    ! Switch to 1-based indices:
    ithetaMin = ithetaMin + 1
    ithetaMax = ithetaMin+localNtheta-1
    izetaMin = izetaMin + 1
    izetaMax = izetaMin+localNzeta-1

    procThatHandlesConstraints = masterProc

    write (procAssignments,fmt="(a,i4,a,i3,a,i3,a,i3,a,i3,a)") " Processor ",myRank," owns theta indices ",ithetaMin," to ",ithetaMax,&
         " and zeta indices ",izetaMin," to ",izetaMax,"\n"
    call PetscSynchronizedPrintf(MPIComm, procAssignments, ierr)
    call PetscSynchronizedFlush(MPIComm, ierr)



    ! *******************************************************************************
    ! Build theta grid, integration weights, and differentiation matrices:
    ! *******************************************************************************

    allocate(theta(Ntheta))
    allocate(thetaWeights(Ntheta))
    allocate(ddtheta(Ntheta,Ntheta))
    allocate(d2dtheta2(Ntheta,Ntheta))
    allocate(theta_preconditioner(Ntheta))
    allocate(thetaWeights_preconditioner(Ntheta))
    allocate(ddtheta_preconditioner(Ntheta,Ntheta))
    allocate(d2dtheta2_preconditioner(Ntheta,Ntheta))

    select case (thetaDerivativeScheme)
    case (0)
       scheme = 20
    case (1)
       scheme = 0
    case (2)
       scheme = 10
    case default
       if (masterProc) then
          print *,"Error! Invalid setting for thetaDerivativeScheme"
       end if
       stop
    end select

    call uniformDiffMatrices(Ntheta, 0, two*pi, scheme, theta, thetaWeights, ddtheta, d2dtheta2)

    ! If needed, also make a sparser differentiation matrix for the preconditioner:
    if (preconditioner_theta==1) then
       scheme = 0
       call uniformDiffMatrices(Ntheta, 0, two*pi, scheme, theta_preconditioner, &
            thetaWeights_preconditioner, ddtheta_preconditioner, d2dtheta2_preconditioner)
    end if

    ! The following arrays will not be needed:
    deallocate(d2dtheta2)
    deallocate(theta_preconditioner)
    deallocate(thetaWeights_preconditioner)
    deallocate(d2dtheta2_preconditioner)

    ! *******************************************************************************
    ! Build zeta grid, integration weights, and differentiation matrices:
    ! *******************************************************************************

    call setNPeriods()

    zetaMax = 2*pi/NPeriods

    allocate(zeta(Nzeta))
    allocate(zetaWeights(Nzeta))
    allocate(ddzeta(Nzeta,Nzeta))
    allocate(d2dzeta2(Nzeta,Nzeta))
    allocate(zeta_preconditioner(Nzeta))
    allocate(zetaWeights_preconditioner(Nzeta))
    allocate(ddzeta_preconditioner(Nzeta,Nzeta))
    allocate(d2dzeta2_preconditioner(Nzeta,Nzeta))

    select case (zetaDerivativeScheme)
    case (0)
       scheme = 20
    case (1)
       scheme = 0
    case (2)
       scheme = 10
    case default
       if (masterProc) then
          print *,"Error! Invalid setting for zetaDerivativeScheme"
       end if
       stop
    end select

    if (Nzeta==1) then
       ! Axisymmetry:
       zeta = 0
       zetaWeights = 2*pi
       ddzeta = 0
       d2dzeta2 = 0 ! d2dzeta2 is not actually used.
    else
       call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta, zetaWeights, ddzeta, d2dzeta2)
    end if

    ! If needed, also make a sparser differentiation matrix for the preconditioner:
    if (preconditioner_zeta==1) then
       if (Nzeta==1) then
          zeta_preconditioner = 0
          zetaWeights_preconditioner = 2*pi
          ddzeta_preconditioner = 0
          d2dzeta2_preconditioner = 0
       else
          scheme = 0
          call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta_preconditioner, &
               zetaWeights_preconditioner, ddzeta_preconditioner, d2dzeta2_preconditioner)
       end if
    end if

    zetaWeights = zetaWeights * NPeriods

    ! The following arrays will not be needed:
    deallocate(d2dzeta2)
    deallocate(zeta_preconditioner)
    deallocate(zetaWeights_preconditioner)
    deallocate(d2dzeta2_preconditioner)

    ! *******************************************************************************
    ! Build x grids, integration weights, and differentiation matrices.
    ! Also build interpolation matrices to map functions from one x grid to the other.
    ! *******************************************************************************

    allocate(x(Nx))
    allocate(xWeights(Nx))
    call makeXGrid(Nx, x, xWeights)
    xWeights = xWeights / exp(-x*x)
    xMaxNotTooSmall = max(x(Nx), xMax)
    allocate(x2(Nx))
    x2=x*x

    allocate(ddx(Nx,Nx))
    allocate(d2dx2(Nx,Nx))
    allocate(ddx_preconditioner(Nx,Nx))
    call makeXPolynomialDiffMatrices(x,ddx,d2dx2)

    NxPotentials = ceiling(xMaxNotTooSmall*NxPotentialsPerVth)


    allocate(xPotentials(NxPotentials))
    allocate(xWeightsPotentials(NxPotentials))
    allocate(ddxPotentials(NxPotentials, NxPotentials))
    allocate(d2dx2Potentials(NxPotentials, NxPotentials))
    call uniformDiffMatrices(NxPotentials, zero, xMaxNotTooSmall, 12, xPotentials, &
         xWeightsPotentials, ddxPotentials, d2dx2Potentials)
    maxxPotentials = xPotentials(NxPotentials)

    deallocate(xWeightsPotentials)

    allocate(expx2(Nx))
    expx2 = exp(-x*x)

    allocate(regridPolynomialToUniform(NxPotentials, Nx))
    call polynomialInterpolationMatrix(Nx, NxPotentials, x, xPotentials, &
         expx2, exp(-xPotentials*xPotentials), regridPolynomialToUniform)

    ddx_preconditioner = 0
    select case (preconditioner_x)
    case (0)
       ! No simplification in x:
       ddx_preconditioner = ddx
    case (1)
       ! Keep only diagonal terms in x:
       do i=1,Nx
          ddx_preconditioner(i,i) = ddx(i,i)
       end do
    case (2)
       ! Keep only upper-triangular terms in x:
       do i=1,Nx
          do j=i,Nx
             ddx_preconditioner(i,j) = ddx(i,j)
          end do
       end do
    case (3)
       ! Keep only tridiagonal terms in x:
       do i=1,Nx
          do j=1,Nx
             if (abs(i-j) <= 1) then
                ddx_preconditioner(i,j) = ddx(i,j)
             end if
          end do
       end do
    case (4)
       ! Keep only diagonal and super-diagonal in x:
       do i=1,Nx
          ddx_preconditioner(i,i) = ddx(i,i)
       end do
       do i=1,(Nx-1)
          ddx_preconditioner(i,i+1) = ddx(i,i+1)
       end do
    case default
       print *,"Error! Invalid preconditioner_x"
       stop
    end select

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Evaluate the magnetic field (and its derivatives) on the (theta, zeta) grid.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    allocate(DHat(Ntheta,Nzeta))

    allocate(BHat(Ntheta,Nzeta))
    allocate(dBHatdtheta(Ntheta,Nzeta))
    allocate(dBHatdzeta(Ntheta,Nzeta))
    allocate(dBHatdrho(Ntheta,Nzeta))

    allocate(BHat_sub_psi(Ntheta,Nzeta))
    allocate(dBHat_sub_psi_dtheta(Ntheta,Nzeta))
    allocate(dBHat_sub_psi_dzeta(Ntheta,Nzeta))

    allocate(BHat_sub_theta(Ntheta,Nzeta))
    allocate(dBHat_sub_theta_dpsiHat(Ntheta,Nzeta))
    allocate(dBHat_sub_theta_dzeta(Ntheta,Nzeta))

    allocate(BHat_sub_zeta(Ntheta,Nzeta))
    allocate(dBHat_sub_zeta_dpsiHat(Ntheta,Nzeta))
    allocate(dBHat_sub_zeta_dtheta(Ntheta,Nzeta))

    allocate(BHat_sup_theta(Ntheta,Nzeta))
    allocate(dBHat_sup_theta_dpsiHat(Ntheta,Nzeta))
    allocate(dBHat_sup_theta_dzeta(Ntheta,Nzeta))

    allocate(BHat_sup_zeta(Ntheta,Nzeta))
    allocate(dBHat_sup_zeta_dpsiHat(Ntheta,Nzeta))
    allocate(dBHat_sup_zeta_dtheta(Ntheta,Nzeta))

    allocate(NTVKernel(Ntheta,Nzeta))

    call computeBHat()

    if (masterProc) then
       print *,"---- Geometry parameters: ----"
       print *,"Geometry scheme = ", geometryScheme
       if (coordinateSystem==0) then
          print *,"iota (rotational transform) = ", iota
          print *,"GHat (Toroidal field Boozer component) = ", GHat
          print *,"IHat (Poloidal field Boozer component) = ", IHat
       end if
       print *,"psiAHat (Normalized toroidal flux at the last closed flux surface) = ", psiAHat
       if (geometryScheme==1) then
          print *,"epsilon_t = ", epsilon_t
          print *,"epsilon_h = ", epsilon_h
          print *,"epsilon_antisymm = ", epsilon_antisymm
       end if
    end if

    ! *********************************************************
    ! Compute a few quantities related to the magnetic field:
    ! *********************************************************

    VPrimeHat = 0
    do itheta=1,Ntheta
       do izeta=1,Nzeta
          VPrimeHat = VPrimeHat + thetaWeights(itheta)*zetaWeights(izeta) / (BHat(itheta,izeta) * BHat(itheta,izeta))
       end do
    end do

    FSABHat2 = 4*pi*pi/VPrimeHat

    ! *********************************************************
    ! Allocate some arrays that will be used later for output quantities:
    ! *********************************************************

    allocate(FSADensityPerturbation(Nspecies))
    allocate(FSABFlow(Nspecies))
    allocate(FSAPressurePerturbation(Nspecies))
    allocate(particleFlux(Nspecies))
    allocate(momentumFlux(Nspecies))
    allocate(heatFlux(Nspecies))
    allocate(NTV(Nspecies)) 

    allocate(densityPerturbation(Nspecies,Ntheta,Nzeta))
    allocate(flow(Nspecies,Ntheta,Nzeta))
    allocate(pressurePerturbation(Nspecies,Ntheta,Nzeta))
    allocate(particleFluxBeforeSurfaceIntegral(Nspecies,Ntheta,Nzeta))
    allocate(momentumFluxBeforeSurfaceIntegral(Nspecies,Ntheta,Nzeta))
    allocate(heatFluxBeforeSurfaceIntegral(Nspecies,Ntheta,Nzeta))
    allocate(NTVBeforeSurfaceIntegral(Nspecies,Ntheta,Nzeta)) 

    allocate(jHat(Ntheta,Nzeta))
    allocate(Phi1Hat(Ntheta,Nzeta))

    select case (constraintScheme)
    case (0)
    case (1)
       allocate(sources(Nspecies,2))
    case (2)
       allocate(sources(Nspecies,Nx))
    case default
       print *,"Error! Invalid setting for constraintScheme."
       stop
    end select


    if (masterProc) then
       print *,"------------------------------------------------------"
       print *,"Done creating grids."
    end if

  end subroutine createGrids
