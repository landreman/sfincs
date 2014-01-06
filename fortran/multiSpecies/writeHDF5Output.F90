#define HAVE_PARALLEL_HDF5

module writeHDF5Output

  use globalVariables
  use scan
  use petscsysdef
  use HDF5

  implicit none

#include <finclude/petscsysdef.h>

  integer, private :: HDF5Error
  integer(HID_T), private :: HDF5FileID, parallelID, dspaceIDForScalar
  integer(HID_T), private :: dspaceIDForSpecies
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForZeta
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForTheta
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForThetaZeta
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForSpeciesThetaZeta
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForSources
  integer(HID_T), private :: dspaceIDForTransportMatrix
  integer(HID_T), dimension(:), allocatable, private :: groupIDs

  integer(HID_T), private :: dsetID_programMode

  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_NSpecies
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Ntheta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Nzeta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Nxi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_NL
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Nx
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_NxPotentialsPerVth
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_xMax
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_solverTolerance
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_theta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_zeta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_thetaDerivativeScheme
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_species
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_x
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_x_min_L
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_xi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_theta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_zeta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_constraintScheme
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_BHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dBHatdzeta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dBHatdtheta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_B0OverBBar
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_GHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_IHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_iota
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_epsilon_t
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_epsilon_h
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_epsilon_antisymm
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_NPeriods
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_helicity_l
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_helicity_n
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_helicity_antisymm_l
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_helicity_antisymm_n
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_speciesMode
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Delta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_alpha
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_psiAHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_nu_n
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Zs
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_mHats
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_THats
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_nHats
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dTHatdpsiNs
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dnHatdpsiNs
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_EParallelHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dPhiHatdpsiN
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_collisionOperator
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_include_fDivVE_term
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_includeXDotTerm
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_includeElectricFieldTermInXiDot
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_useDKESExBDrift
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_sources
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_densityPerturbation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_flow
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_pressurePerturbation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_particleFluxBeforeSurfaceIntegral
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_momentumFluxBeforeSurfaceIntegral
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_heatFluxBeforeSurfaceIntegral
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSADensityPerturbation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSABFlow
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSAPressurePerturbation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_particleFlux
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_momentumFlux
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_heatFlux
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_jHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSABjHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_elapsedTime
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_didItConverge
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_integerToRepresentTrue
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_integerToRepresentFalse
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_VPrimeHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSABHat2
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_useIterativeSolver
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_transportMatrix


  integer(HSIZE_T), dimension(1), parameter, private :: dimForScalar = 1
  integer(HSIZE_T), dimension(1), private :: dimForSpecies
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForZeta
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForTheta
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForThetaZeta
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForSpeciesThetaZeta
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForSources
  integer(HSIZE_T), dimension(2), parameter, private :: dimForTransportMatrix = 3

contains

  ! -----------------------------------------------------------------------------------

  subroutine openOutputFile()

    implicit none

#ifdef HAVE_PARALLEL_HDF5
    if (outputScheme > 0) then
#else
    if (outputScheme > 0 .and. masterProcInSubComm) then
#endif
       call h5open_f(HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error initializing HDF5."
          stop
       end if

#ifdef HAVE_PARALLEL_HDF5
       ! Initialize some stuff related to parallel file access
       call h5pcreate_f(H5P_FILE_ACCESS_F, parallelID, HDF5Error)
       call h5pset_fapl_mpio_f(parallelID, MPI_COMM_WORLD, MPI_INFO_NULL, HDF5Error)

       call h5fcreate_f(outputFilename, H5F_ACC_TRUNC_F, HDF5FileID, HDF5Error, access_prp=parallelID)
#else
       call h5fcreate_f(outputFilename, H5F_ACC_TRUNC_F, HDF5FileID, HDF5Error)
#endif
       if (HDF5Error < 0) then
          print *,"Error opening HDF5 output file."
          stop
       end if

    end if

  end subroutine openOutputFile

  ! -----------------------------------------------------------------------------------

  subroutine createHDF5Structures()

    implicit none

    integer :: i, rank
    character(20) :: groupName

#ifdef HAVE_PARALLEL_HDF5
    if (outputScheme > 0) then
#else
    if (outputScheme > 0 .and. masterProcInSubComm) then
#endif
       allocate(groupIDs(numRunsInScan))

       allocate(dsetIDs_Nspecies(numRunsInScan))
       allocate(dsetIDs_Ntheta(numRunsInScan))
       allocate(dsetIDs_Nzeta(numRunsInScan))
       allocate(dsetIDs_Nxi(numRunsInScan))
       allocate(dsetIDs_NL(numRunsInScan))
       allocate(dsetIDs_Nx(numRunsInScan))
       allocate(dsetIDs_NxPotentialsPerVth(numRunsInScan))
       allocate(dsetIDs_xMax(numRunsInScan))
       allocate(dsetIDs_solverTolerance(numRunsInScan))
       allocate(dsetIDs_theta(numRunsInScan))
       allocate(dsetIDs_zeta(numRunsInScan))
       allocate(dsetIDs_thetaDerivativeScheme(numRunsInScan))
       allocate(dsetIDs_preconditioner_species(numRunsInScan))
       allocate(dsetIDs_preconditioner_x(numRunsInScan))
       allocate(dsetIDs_preconditioner_x_min_L(numRunsInScan))
       allocate(dsetIDs_preconditioner_xi(numRunsInScan))
       allocate(dsetIDs_preconditioner_theta(numRunsInScan))
       allocate(dsetIDs_preconditioner_zeta(numRunsInScan))
       allocate(dsetIDs_constraintScheme(numRunsInScan))
       allocate(dsetIDs_BHat(numRunsInScan))
       allocate(dsetIDs_dBHatdtheta(numRunsInScan))
       allocate(dsetIDs_dBHatdzeta(numRunsInScan))
       allocate(dsetIDs_B0OverBBar(numRunsInScan))
       allocate(dsetIDs_GHat(numRunsInScan))
       allocate(dsetIDs_IHat(numRunsInScan))
       allocate(dsetIDs_iota(numRunsInScan))
       allocate(dsetIDs_epsilon_t(numRunsInScan))
       allocate(dsetIDs_epsilon_h(numRunsInScan))
       allocate(dsetIDs_epsilon_antisymm(numRunsInScan))
       allocate(dsetIDs_NPeriods(numRunsInScan))
       allocate(dsetIDs_helicity_l(numRunsInScan))
       allocate(dsetIDs_helicity_n(numRunsInScan))
       allocate(dsetIDs_helicity_antisymm_l(numRunsInScan))
       allocate(dsetIDs_helicity_antisymm_n(numRunsInScan))
       allocate(dsetIDs_speciesMode(numRunsInScan))
       allocate(dsetIDs_Delta(numRunsInScan))
       allocate(dsetIDs_alpha(numRunsInScan))
       allocate(dsetIDs_psiAHat(numRunsInScan))
       allocate(dsetIDs_nu_n(numRunsInScan))
       allocate(dsetIDs_Zs(numRunsInScan))
       allocate(dsetIDs_mHats(numRunsInScan))
       allocate(dsetIDs_THats(numRunsInScan))
       allocate(dsetIDs_nHats(numRunsInScan))
       allocate(dsetIDs_dTHatdpsiNs(numRunsInScan))
       allocate(dsetIDs_dnHatdpsiNs(numRunsInScan))
       allocate(dsetIDs_EParallelHat(numRunsInScan))
       allocate(dsetIDs_dPhiHatdpsiN(numRunsInScan))
       allocate(dsetIDs_collisionOperator(numRunsInScan))
       allocate(dsetIDs_include_fDivVE_Term(numRunsInScan))
       allocate(dsetIDs_includeXDotTerm(numRunsInScan))
       allocate(dsetIDs_includeElectricFieldTermInXiDot(numRunsInScan))
       allocate(dsetIDs_useDKESExBDrift(numRunsInScan))
       allocate(dsetIDs_sources(numRunsInScan))
       allocate(dsetIDs_densityPerturbation(numRunsInScan))
       allocate(dsetIDs_flow(numRunsInScan))
       allocate(dsetIDs_pressurePerturbation(numRunsInScan))
       allocate(dsetIDs_particleFluxBeforeSurfaceIntegral(numRunsInScan))
       allocate(dsetIDs_momentumFluxBeforeSurfaceIntegral(numRunsInScan))
       allocate(dsetIDs_heatFluxBeforeSurfaceIntegral(numRunsInScan))
       allocate(dsetIDs_FSADensityPerturbation(numRunsInScan))
       allocate(dsetIDs_FSABFlow(numRunsInScan))
       allocate(dsetIDs_FSAPressurePerturbation(numRunsInScan))
       allocate(dsetIDs_particleFlux(numRunsInScan))
       allocate(dsetIDs_momentumFlux(numRunsInScan))
       allocate(dsetIDs_heatFlux(numRunsInScan))
       allocate(dsetIDs_jHat(numRunsInScan))
       allocate(dsetIDs_FSABjHat(numRunsInScan))
       allocate(dsetIDs_elapsedTime(numRunsInScan))
       allocate(dsetIDs_didItConverge(numRunsInScan))
       allocate(dsetIDs_integerToRepresentTrue(numRunsInScan))
       allocate(dsetIDs_integerToRepresentFalse(numRunsInScan))
       allocate(dsetIDs_VPrimeHat(numRunsInScan))
       allocate(dsetIDs_FSABHat2(numRunsInScan))
       allocate(dsetIDs_useIterativeSolver(numRunsInScan))
       allocate(dsetIDs_transportMatrix(numRunsInScan))

       allocate(dspaceIDForZeta(numRunsInScan))
       allocate(dspaceIDForTheta(numRunsInScan))
       allocate(dspaceIDForThetaZeta(numRunsInScan))
       allocate(dspaceIDForSpeciesThetaZeta(numRunsInScan))
       allocate(dspaceIDForSources(numRunsInScan))

       allocate(dimForZeta(numRunsInScan,1))
       allocate(dimForTheta(numRunsInScan,1))
       allocate(dimForThetaZeta(numRunsInScan,2))
       allocate(dimForSpeciesThetaZeta(numRunsInScan,3))
       allocate(dimForSources(numRunsInScan,2))

       ! Create a dataspace for storing single numbers:
       rank = 0
       call h5screate_simple_f(rank, dimForScalar, dspaceIDForScalar, HDF5Error)

       rank = 1
       dimForSpecies = Nspecies
       call h5screate_simple_f(rank, dimForSpecies, dspaceIDForSpecies, HDF5Error)

       rank = 2
       call h5screate_simple_f(rank, dimForTransportMatrix, dspaceIDForTransportMatrix, HDF5Error)

       ! Save programMode in the file:
       call h5dcreate_f(HDF5FileID, "programMode", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID_programMode, HDF5Error)
       if (masterProc) then
          call h5dwrite_f(dsetID_programMode, H5T_NATIVE_INTEGER, programMode, dimForScalar, HDF5Error)
       end if
       call h5dclose_f(dsetID_programMode, HDF5Error)

       do i=1,numRunsInScan
          ! Create a group to hold all data for the run:
          write (groupName, "(a, i3)") "run",i
          call h5gcreate_f(HDF5FileID, groupName, groupIDs(i), HDF5Error)

          ! Create dataspaces that depend on resolution parameters:
          rank = 1
          dimForZeta(i,1)=NzetasForScan(i)
          call h5screate_simple_f(rank, dimForZeta(i,:), dspaceIDForZeta(i), HDF5Error)

          dimForTheta(i,1)=NthetasForScan(i)
          call h5screate_simple_f(rank, dimForTheta(i,:), dspaceIDForTheta(i), HDF5Error)

          rank = 2
          dimForThetaZeta(i,1)=NthetasForScan(i)
          dimForThetaZeta(i,2)=NzetasForScan(i)
          call h5screate_simple_f(rank, dimForThetaZeta(i,:), dspaceIDForThetaZeta(i), HDF5Error)

          rank = 3
          dimForSpeciesThetaZeta(i,1)=Nspecies
          dimForSpeciesThetaZeta(i,2)=NthetasForScan(i)
          dimForSpeciesThetaZeta(i,3)=NzetasForScan(i)
          call h5screate_simple_f(rank, dimForSpeciesThetaZeta(i,:), dspaceIDForSpeciesThetaZeta(i), HDF5Error)

          dimForSources(i,1) = Nspecies
          select case (constraintSchemesForScan(i))
          case (0)
             dimForSources(i,2) = 1
          case (1)
             dimForSources(i,2) = 2
          case (2)
             dimForSources(i,2) = NxsForScan(i)
          case default
             print *,"Error in writeHDF5Output! Invalid setting for constraintScheme."
          end select
          rank = 2
          call h5screate_simple_f(rank, dimForSources(i,:), dspaceIDForSources(i), HDF5Error)

          ! Create datasets for each quantity in each run:
          call h5dcreate_f(groupIDs(i), "Nspecies", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_Nspecies(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Ntheta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_Ntheta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Nzeta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_Nzeta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Nxi", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_Nxi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "NL", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_NL(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Nx", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_Nx(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "NxPotentialsPerVth", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_NxPotentialsPerVth(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "xMax", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_xMax(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "solverTolerance", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_solverTolerance(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "theta", H5T_NATIVE_DOUBLE, dspaceIDForTheta(i), &
               dsetIDs_theta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "zeta", H5T_NATIVE_DOUBLE, dspaceIDForZeta(i), &
               dsetIDs_zeta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "thetaDerivativeScheme", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_thetaDerivativeScheme(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_species", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_species(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_x", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_x(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_x_min_L", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_x_min_L(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_xi", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_xi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_theta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_theta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_zeta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_zeta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "constraintScheme", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_constraintScheme(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "BHat", H5T_NATIVE_DOUBLE, dspaceIDForThetaZeta(i), &
               dsetIDs_BHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(BHat)d(theta)", H5T_NATIVE_DOUBLE, dspaceIDForThetaZeta(i), &
               dsetIDs_dBHatdtheta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(BHat)d(zeta)", H5T_NATIVE_DOUBLE, dspaceIDForThetaZeta(i), &
               dsetIDs_dBHatdzeta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "B0OverBBar", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_B0OverBBar(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "GHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_GHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "IHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_IHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "iota", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_iota(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "epsilon_t", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_epsilon_t(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "epsilon_h", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_epsilon_h(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "epsilon_antisymm", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_epsilon_antisymm(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "NPeriods", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_NPeriods(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "helicity_l", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_helicity_l(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "helicity_n", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_helicity_n(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "helicity_antisymm_l", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_helicity_antisymm_l(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "helicity_antisymm_n", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_helicity_antisymm_n(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "speciesMode", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_speciesMode(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Delta", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_Delta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "alpha", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_alpha(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "psiAHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_psiAHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "nu_n", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_nu_n(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Zs", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_Zs(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "mHats", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_mHats(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "THats", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_THats(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "nHats", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_nHats(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(THat)d(psi_N)", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_dTHatdpsiNs(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(nHat)d(psi_N)", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_dnHatdpsiNs(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "EParallelHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_EParallelHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(PhiHat)d(psi_N)", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_dPhiHatdpsiN(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "collisionOperator", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_collisionOperator(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "include_fDivVE_Term", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_include_fDivVE_Term(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "includeXDotTerm", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_includeXDotTerm(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "includeElectricFieldTermInXiDot", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_includeElectricFieldTermInXiDot(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "useDKESExBDrift", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_useDKESExBDrift(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "sources", H5T_NATIVE_DOUBLE, dspaceIDForSources(i), &
               dsetIDs_sources(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "densityPerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta(i), &
               dsetIDs_densityPerturbation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "flow", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta(i), &
               dsetIDs_flow(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "pressurePerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta(i), &
               dsetIDs_pressurePerturbation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "particleFluxBeforeSurfaceIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta(i), &
               dsetIDs_particleFluxBeforeSurfaceIntegral(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "momentumFluxBeforeSurfaceIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta(i), &
               dsetIDs_momentumFluxBeforeSurfaceIntegral(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "heatFluxBeforeSurfaceIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta(i), &
               dsetIDs_heatFluxBeforeSurfaceIntegral(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSADensityPerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_FSADensityPerturbation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSABFlow", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_FSABFlow(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSAPressurePerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_FSAPressurePerturbation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "particleFlux", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_particleFlux(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "momentumFlux", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_momentumFlux(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "heatFlux", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
               dsetIDs_heatFlux(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "jHat", H5T_NATIVE_DOUBLE, dspaceIDForThetaZeta(i), &
               dsetIDs_jHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSABjHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_FSABjHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "elapsed time (s)", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_elapsedTime(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "didItConverge", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_didItConverge(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "integerToRepresentTrue", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_integerToRepresentTrue(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "integerToRepresentFalse", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_integerToRepresentFalse(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "VPrimeHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_VPrimeHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSABHat2", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_FSABHat2(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "useIterativeSolver", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_useIterativeSolver(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "transportMatrix", H5T_NATIVE_DOUBLE, dspaceIDForTransportMatrix, &
               dsetIDs_transportMatrix(i), HDF5Error)

       end do
    end if

  end subroutine createHDF5Structures

  ! -----------------------------------------------------------------------------------

  subroutine writeRunToOutputFile(runNum)

    implicit none

    integer, intent(in) :: runNum
    integer :: temp

    if (outputScheme > 0 .and. masterProcInSubComm) then

       call h5dwrite_f(dsetIDs_Nspecies(runNum), H5T_NATIVE_INTEGER, &
            Nspecies, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Ntheta(runNum), H5T_NATIVE_INTEGER, &
            Ntheta, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Nzeta(runNum), H5T_NATIVE_INTEGER, &
            Nzeta, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Nxi(runNum), H5T_NATIVE_INTEGER, &
            Nxi, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_NL(runNum), H5T_NATIVE_INTEGER, &
            NL, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Nx(runNum), H5T_NATIVE_INTEGER, &
            Nx, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_NxPotentialsPerVth(runNum), H5T_NATIVE_DOUBLE, &
            NxPotentialsPerVth, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_xMax(runNum), H5T_NATIVE_DOUBLE, &
            xMax, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_solverTolerance(runNum), H5T_NATIVE_DOUBLE, &
            solverTolerance, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_theta(runNum), H5T_NATIVE_DOUBLE, &
            theta, dimForTheta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_zeta(runNum), H5T_NATIVE_DOUBLE, &
            zeta, dimForZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_thetaDerivativeScheme(runNum), H5T_NATIVE_INTEGER, &
            thetaDerivativeScheme, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_species(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_species, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_x(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_x, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_x_min_L(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_x_min_L, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_xi(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_xi, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_theta(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_theta, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_zeta(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_zeta, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_constraintScheme(runNum), H5T_NATIVE_INTEGER, &
            constraintScheme, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_BHat(runNum), H5T_NATIVE_DOUBLE, &
            BHat, dimForThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_dBHatdtheta(runNum), H5T_NATIVE_DOUBLE, &
            dBHatdtheta, dimForThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_dBHatdzeta(runNum), H5T_NATIVE_DOUBLE, &
            dBHatdzeta, dimForThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_B0OverBBar(runNum), H5T_NATIVE_DOUBLE, &
            B0OverBBar, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_GHat(runNum), H5T_NATIVE_DOUBLE, &
            GHat, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_IHat(runNum), H5T_NATIVE_DOUBLE, &
            IHat, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_iota(runNum), H5T_NATIVE_DOUBLE, &
            iota, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_epsilon_t(runNum), H5T_NATIVE_DOUBLE, &
            epsilon_t, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_epsilon_h(runNum), H5T_NATIVE_DOUBLE, &
            epsilon_h, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_epsilon_antisymm(runNum), H5T_NATIVE_DOUBLE, &
            epsilon_antisymm, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_NPeriods(runNum), H5T_NATIVE_INTEGER, &
            NPeriods, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_helicity_l(runNum), H5T_NATIVE_INTEGER, &
            helicity_l, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_helicity_n(runNum), H5T_NATIVE_INTEGER, &
            helicity_n, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_helicity_antisymm_l(runNum), H5T_NATIVE_INTEGER, &
            helicity_antisymm_l, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_helicity_antisymm_n(runNum), H5T_NATIVE_INTEGER, &
            helicity_antisymm_n, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_speciesMode(runNum), H5T_NATIVE_INTEGER, &
            speciesMode, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Delta(runNum), H5T_NATIVE_DOUBLE, &
            Delta, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_alpha(runNum), H5T_NATIVE_DOUBLE, &
            alpha, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_psiAHat(runNum), H5T_NATIVE_DOUBLE, &
            psiAHat, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_nu_n(runNum), H5T_NATIVE_DOUBLE, &
            nu_n, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Zs(runNum), H5T_NATIVE_DOUBLE, &
            Zs, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_mHats(runNum), H5T_NATIVE_DOUBLE, &
            mHats, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_THats(runNum), H5T_NATIVE_DOUBLE, &
            THats, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_nHats(runNum), H5T_NATIVE_DOUBLE, &
            nHats, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_dTHatdpsiNs(runNum), H5T_NATIVE_DOUBLE, &
            dTHatdpsiNs, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_dnHatdpsiNs(runNum), H5T_NATIVE_DOUBLE, &
            dnHatdpsiNs, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_EParallelHat(runNum), H5T_NATIVE_DOUBLE, &
            EParallelHat, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_dPhiHatdpsiN(runNum), H5T_NATIVE_DOUBLE, &
            dPhiHatdpsiN, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_collisionOperator(runNum), H5T_NATIVE_INTEGER, &
            collisionOperator, dimForScalar, HDF5Error)

       if (include_fDivVE_term) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetIDs_include_fDivVE_term(runNum), H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       if (includeXDotTerm) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetIDs_includeXDotTerm(runNum), H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       if (includeElectricFieldTermInXiDot) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetIDs_includeElectricFieldTermInXiDot(runNum), H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       if (useDKESExBDrift) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetIDs_useDKESExBDrift(runNum), H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       if (constraintSchemesForScan(runNum)==0) then
          allocate(sources(1,1))
          sources=zero
          call h5dwrite_f(dsetIDs_sources(runNum), H5T_NATIVE_DOUBLE, &
               sources, dimForSources(runNum,:), HDF5Error)
          deallocate(sources)
       else
          call h5dwrite_f(dsetIDs_sources(runNum), H5T_NATIVE_DOUBLE, &
               sources, dimForSources(runNum,:), HDF5Error)
       end if

       call h5dwrite_f(dsetIDs_densityPerturbation(runNum), H5T_NATIVE_DOUBLE, &
            densityPerturbation, dimForSpeciesThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_flow(runNum), H5T_NATIVE_DOUBLE, &
            flow, dimForSpeciesThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_pressurePerturbation(runNum), H5T_NATIVE_DOUBLE, &
            pressurePerturbation, dimForSpeciesThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_particleFluxBeforeSurfaceIntegral(runNum), H5T_NATIVE_DOUBLE, &
            particleFluxBeforeSurfaceIntegral, dimForSpeciesThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_momentumFluxBeforeSurfaceIntegral(runNum), H5T_NATIVE_DOUBLE, &
            momentumFluxBeforeSurfaceIntegral, dimForSpeciesThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_heatFluxBeforeSurfaceIntegral(runNum), H5T_NATIVE_DOUBLE, &
            heatFluxBeforeSurfaceIntegral, dimForSpeciesThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_FSADensityPerturbation(runNum), H5T_NATIVE_DOUBLE, &
            FSADensityPerturbation, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_FSABFlow(runNum), H5T_NATIVE_DOUBLE, &
            FSABFlow, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_FSAPressurePerturbation(runNum), H5T_NATIVE_DOUBLE, &
            FSAPressurePerturbation, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_particleFlux(runNum), H5T_NATIVE_DOUBLE, &
            particleFlux, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_momentumFlux(runNum), H5T_NATIVE_DOUBLE, &
            momentumFlux, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_heatFlux(runNum), H5T_NATIVE_DOUBLE, &
            heatFlux, dimForSpecies, HDF5Error)

       call h5dwrite_f(dsetIDs_jHat(runNum), H5T_NATIVE_DOUBLE, &
            jHat, dimForThetaZeta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_FSABjHat(runNum), H5T_NATIVE_DOUBLE, &
            FSABjHat, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_elapsedTime(runNum), H5T_NATIVE_DOUBLE, &
            elapsedTime, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_didItConverge(runNum), H5T_NATIVE_INTEGER, &
            didItConverge, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_integerToRepresentTrue(runNum), H5T_NATIVE_INTEGER, &
            integerToRepresentTrue, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_integerToRepresentFalse(runNum), H5T_NATIVE_INTEGER, &
            integerToRepresentFalse, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_VPrimeHat(runNum), H5T_NATIVE_DOUBLE, &
            VPrimeHat, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_FSABHat2(runNum), H5T_NATIVE_DOUBLE, &
            FSABHat2, dimForScalar, HDF5Error)

       if (useIterativeSolver) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetIDs_useIterativeSolver(runNum), H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_transportMatrix(runNum), H5T_NATIVE_DOUBLE, &
            transportMatrix, dimForTransportMatrix, HDF5Error)

    end if

  end subroutine writeRunToOutputFile

  ! -----------------------------------------------------------------------------------

  subroutine closeOutputFile()

    implicit none

    integer :: i

#ifdef HAVE_PARALLEL_HDF5
    if (outputScheme > 0) then
#else
       if (outputScheme > 0 .and. masterProcInSubComm) then
#endif

       do i=1,numRunsInScan

          call h5dclose_f(dsetIDs_Nspecies(i), HDF5Error)
          call h5dclose_f(dsetIDs_Ntheta(i), HDF5Error)
          call h5dclose_f(dsetIDs_Nzeta(i), HDF5Error)
          call h5dclose_f(dsetIDs_Nxi(i), HDF5Error)
          call h5dclose_f(dsetIDs_NL(i), HDF5Error)
          call h5dclose_f(dsetIDs_Nx(i), HDF5Error)
          call h5dclose_f(dsetIDs_NxPotentialsPerVth(i), HDF5Error)
          call h5dclose_f(dsetIDs_xMax(i), HDF5Error)
          call h5dclose_f(dsetIDs_solverTolerance(i), HDF5Error)
          call h5dclose_f(dsetIDs_theta(i), HDF5Error)
          call h5dclose_f(dsetIDs_zeta(i), HDF5Error)
          call h5dclose_f(dsetIDs_thetaDerivativeScheme(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_species(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_x(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_x_min_L(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_xi(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_theta(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_zeta(i), HDF5Error)
          call h5dclose_f(dsetIDs_constraintScheme(i), HDF5Error)
          call h5dclose_f(dsetIDs_BHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_dBHatdtheta(i), HDF5Error)
          call h5dclose_f(dsetIDs_dBHatdzeta(i), HDF5Error)
          call h5dclose_f(dsetIDs_B0OverBBar(i), HDF5Error)
          call h5dclose_f(dsetIDs_GHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_IHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_iota(i), HDF5Error)
          call h5dclose_f(dsetIDs_epsilon_t(i), HDF5Error)
          call h5dclose_f(dsetIDs_epsilon_h(i), HDF5Error)
          call h5dclose_f(dsetIDs_epsilon_antisymm(i), HDF5Error)
          call h5dclose_f(dsetIDs_NPeriods(i), HDF5Error)
          call h5dclose_f(dsetIDs_helicity_l(i), HDF5Error)
          call h5dclose_f(dsetIDs_helicity_n(i), HDF5Error)
          call h5dclose_f(dsetIDs_helicity_antisymm_l(i), HDF5Error)
          call h5dclose_f(dsetIDs_helicity_antisymm_n(i), HDF5Error)
          call h5dclose_f(dsetIDs_speciesMode(i), HDF5Error)
          call h5dclose_f(dsetIDs_Delta(i), HDF5Error)
          call h5dclose_f(dsetIDs_alpha(i), HDF5Error)
          call h5dclose_f(dsetIDs_psiAHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_nu_n(i), HDF5Error)
          call h5dclose_f(dsetIDs_Zs(i), HDF5Error)
          call h5dclose_f(dsetIDs_mHats(i), HDF5Error)
          call h5dclose_f(dsetIDs_THats(i), HDF5Error)
          call h5dclose_f(dsetIDs_nHats(i), HDF5Error)
          call h5dclose_f(dsetIDs_dTHatdpsiNs(i), HDF5Error)
          call h5dclose_f(dsetIDs_dnHatdpsiNs(i), HDF5Error)
          call h5dclose_f(dsetIDs_EParallelHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_dPhiHatdpsiN(i), HDF5Error)
          call h5dclose_f(dsetIDs_collisionOperator(i), HDF5Error)
          call h5dclose_f(dsetIDs_include_fDivVE_term(i), HDF5Error)
          call h5dclose_f(dsetIDs_includeXDotTerm(i), HDF5Error)
          call h5dclose_f(dsetIDs_includeElectricFieldTermInXiDot(i), HDF5Error)
          call h5dclose_f(dsetIDs_useDKESExBDrift(i), HDF5Error)
          call h5dclose_f(dsetIDs_sources(i), HDF5Error)
          call h5dclose_f(dsetIDs_densityPerturbation(i), HDF5Error)
          call h5dclose_f(dsetIDs_flow(i), HDF5Error)
          call h5dclose_f(dsetIDs_pressurePerturbation(i), HDF5Error)
          call h5dclose_f(dsetIDs_particleFluxBeforeSurfaceIntegral(i), HDF5Error)
          call h5dclose_f(dsetIDs_momentumFluxBeforeSurfaceIntegral(i), HDF5Error)
          call h5dclose_f(dsetIDs_heatFluxBeforeSurfaceIntegral(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSADensityPerturbation(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSABFlow(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSAPressurePerturbation(i), HDF5Error)
          call h5dclose_f(dsetIDs_particleFlux(i), HDF5Error)
          call h5dclose_f(dsetIDs_momentumFlux(i), HDF5Error)
          call h5dclose_f(dsetIDs_heatFlux(i), HDF5Error)
          call h5dclose_f(dsetIDs_jHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSABjHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_elapsedTime(i), HDF5Error)
          call h5dclose_f(dsetIDs_didItConverge(i), HDF5Error)
          call h5dclose_f(dsetIDs_integerToRepresentTrue(i), HDF5Error)
          call h5dclose_f(dsetIDs_integerToRepresentFalse(i), HDF5Error)
          call h5dclose_f(dsetIDs_VPrimeHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSABHat2(i), HDF5Error)
          call h5dclose_f(dsetIDs_useIterativeSolver(i), HDF5Error)
          call h5dclose_f(dsetIDs_transportMatrix(i), HDF5Error)


          call h5gclose_f(groupIDs(i), HDF5Error)

          call h5sclose_f(dspaceIDForTheta(i), HDF5Error)
          call h5sclose_f(dspaceIDForZeta(i), HDF5Error)
          call h5sclose_f(dspaceIDForThetaZeta(i), HDF5Error)
          call h5sclose_f(dspaceIDForSpeciesThetaZeta(i), HDF5Error)
          call h5sclose_f(dspaceIDForSources(i), HDF5Error)
       end do

       call h5sclose_f(dspaceIDForScalar, HDF5Error)
       call h5sclose_f(dspaceIDForSpecies, HDF5Error)
       call h5sclose_f(dspaceIDForTransportMatrix, HDF5Error)
       call h5pclose_f(parallelID, HDF5Error)
       call h5fclose_f(HDF5FileID, HDF5Error)
       call h5close_f(HDF5Error)
    end if

  end subroutine closeOutputFile

end module writeHDF5Output

