module writeHDF5Output

  use globalVariables
  use petscsysdef
  use HDF5

  implicit none

#include <finclude/petscsysdef.h>

  integer, private :: HDF5Error
  integer(HID_T), private :: HDF5FileID

  integer(HSIZE_T), dimension(1), parameter, private :: dimForScalar = 1
  integer(HSIZE_T), dimension(1), private :: dimForSpecies
  integer(HSIZE_T), dimension(1), private :: dimForTheta
  integer(HSIZE_T), dimension(1), private :: dimForZeta
  integer(HSIZE_T), dimension(1), private :: dimForx
  integer(HSIZE_T), dimension(2), private :: dimForThetaZeta
  integer(HSIZE_T), dimension(3), private :: dimForSpeciesThetaZeta
  integer(HSIZE_T), dimension(2), private :: dimForSources

  integer(HID_T), private :: dspaceIDForScalar
  integer(HID_T), private :: dspaceIDForSpecies
  integer(HID_T), private :: dspaceIDForTheta
  integer(HID_T), private :: dspaceIDForZeta
  integer(HID_T), private :: dspaceIDForx
  integer(HID_T), private :: dspaceIDForThetaZeta
  integer(HID_T), private :: dspaceIDForSpeciesThetaZeta
  integer(HID_T), private :: dspaceIDForSources

contains

  ! -----------------------------------------------------------------------------------

  subroutine openOutputFile()

    implicit none

    if (masterProc) then
       call h5open_f(HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error initializing HDF5."
          stop
       end if

       call h5fcreate_f(outputFilename, H5F_ACC_TRUNC_F, HDF5FileID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening HDF5 output file."
          stop
       end if

    end if

  end subroutine openOutputFile

  ! -----------------------------------------------------------------------------------

  subroutine createHDF5Structures()

    ! This subroutine creates the dataspaces.

    implicit none

    integer :: i, rank

    if (masterProc) then

       ! Create a dataspace for storing single numbers:
       rank = 0
       call h5screate_simple_f(rank, dimForScalar, dspaceIDForScalar, HDF5Error)

       ! Create dataspaces that depend on resolution parameters:
       rank = 1
       dimForSpecies = Nspecies
       call h5screate_simple_f(rank, dimForSpecies, dspaceIDForSpecies, HDF5Error)

       rank = 1
       dimForZeta = Nzeta
       call h5screate_simple_f(rank, dimForZeta, dspaceIDForZeta, HDF5Error)

       dimForTheta = Ntheta
       call h5screate_simple_f(rank, dimForTheta, dspaceIDForTheta, HDF5Error)

       dimForx = Nx
       call h5screate_simple_f(rank, dimForx, dspaceIDForx, HDF5Error)

       rank = 2
       dimForThetaZeta(1) = Ntheta
       dimForThetaZeta(2) = Nzetas
       call h5screate_simple_f(rank, dimForThetaZeta, dspaceIDForThetaZeta, HDF5Error)

       rank = 3
       dimForSpeciesThetaZeta(1) = Nspecies
       dimForSpeciesThetaZeta(2) = Ntheta
       dimForSpeciesThetaZeta(3) = Nzeta
       call h5screate_simple_f(rank, dimForSpeciesThetaZeta, dspaceIDForSpeciesThetaZeta, HDF5Error)

       dimForSources(1) = Nspecies
       select case (constraintScheme)
       case (0)
          dimForSources(2) = 1
       case (1)
          dimForSources(2) = 2
       case (2)
          dimForSources(2) = Nx
       case default
          print *,"Error in writeHDF5Output! Invalid setting for constraintScheme."
       end select
       rank = 2
       call h5screate_simple_f(rank, dimForSources, dspaceIDForSources, HDF5Error)

    end if

  end subroutine createHDF5Structures

  ! -----------------------------------------------------------------------------------

  subroutine writeOutputFile()

    implicit none

    integer(HID_T), private :: dsetID
    integer :: temp

    call createHDF5Structures()

    if (masterProc) then

       call saveInputFileToHDF5()


       ! For each variable we want in the HDF5 file, we must do 3 steps:
       ! 1. Create a dataset id,
       ! 2. Write the dataset to the file,
       ! 3. Close the dataset.

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "Nspecies", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            Nspecies, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "Ntheta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            Ntheta, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "Nzeta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            Nzeta, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "Nxi", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            Nxi, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "NL", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            NL, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "Nx", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            Nx, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "NxPotentialsPerVth", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            NxPotentialsPerVth, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "xMax", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            xMax, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "solverTolerance", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            solverTolerance, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "theta", H5T_NATIVE_DOUBLE, dspaceIDForTheta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            theta, dimForTheta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "zeta", H5T_NATIVE_DOUBLE, dspaceIDForZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            zeta, dimForZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "x", H5T_NATIVE_DOUBLE, dspaceIDForx, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            x, dimForx, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "thetaDerivativeScheme", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            thetaDerivativeScheme, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "preconditioner_species", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            preconditioner_species, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "preconditioner_x", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            preconditioner_x, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "preconditioner_x_min_L", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            preconditioner_x_min_L, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "preconditioner_xi", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            preconditioner_xi, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "preconditioner_theta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            preconditioner_theta, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "preconditioner_zeta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            preconditioner_zeta, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "constraintScheme", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            constraintScheme, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "BHat", H5T_NATIVE_DOUBLE, dspaceIDForThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            BHat, dimForThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "d(BHat)d(theta)", H5T_NATIVE_DOUBLE, dspaceIDForThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            dBHatdtheta, dimForThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "d(BHat)d(zeta)", H5T_NATIVE_DOUBLE, dspaceIDForThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            dBHatdzeta, dimForThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "B0OverBBar", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            B0OverBBar, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "GHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            GHat, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "IHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            IHat, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "iota", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            iota, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "epsilon_t", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            epsilon_t, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "epsilon_h", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            epsilon_h, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "epsilon_antisymm", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            epsilon_antisymm, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "NPeriods", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            NPeriods, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "helicity_l", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            helicity_l, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "helicity_n", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            helicity_n, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "helicity_antisymm_l", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            helicity_antisymm_l, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "helicity_antisymm_n", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            helicity_antisymm_n, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "Delta", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            Delta, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "alpha", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            alpha, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "psiAHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            psiAHat, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "nu_n", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            nu_n, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "Zs", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            Zs, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "mHats", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            mHats, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "THats", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            THats, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "nHats", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            nHats, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "d(THat)d(psi_N)", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            dTHatdpsiNs, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "d(nHat)d(psi_N)", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            dnHatdpsiNs, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "EParallelHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            EParallelHat, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "d(PhiHat)d(psi_N)", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            dPhiHatdpsiN, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "collisionOperator", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            collisionOperator, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "include_fDivVE_Term", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       if (include_fDivVE_term) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "includeXDotTerm", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       if (includeXDotTerm) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "includeElectricFieldTermInXiDot", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       if (includeElectricFieldTermInXiDot) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "useDKESExBDrift", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       if (useDKESExBDrift) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "sources", H5T_NATIVE_DOUBLE, dspaceIDForSources, &
            dsetID, HDF5Error)

       if (constraintScheme==0) then
          allocate(sources(1,1))
          sources=zero
          call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
               sources, dimForSources, HDF5Error)
          deallocate(sources)
       else
          call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
               sources, dimForSources, HDF5Error)
       end if

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "densityPerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            densityPerturbation, dimForSpeciesThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "flow", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            flow, dimForSpeciesThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "pressurePerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            pressurePerturbation, dimForSpeciesThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "particleFluxBeforeSurfaceIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            particleFluxBeforeSurfaceIntegral, dimForSpeciesThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "momentumFluxBeforeSurfaceIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            momentumFluxBeforeSurfaceIntegral, dimForSpeciesThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "heatFluxBeforeSurfaceIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            heatFluxBeforeSurfaceIntegral, dimForSpeciesThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "NTVBeforeSurfaceIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaZeta, &
            dsetID, HDF5Error) 

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            NTVBeforeSurfaceIntegral, dimForSpeciesThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "FSADensityPerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            FSADensityPerturbation, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "FSABFlow", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            FSABFlow, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "FSAPressurePerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            FSAPressurePerturbation, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "particleFlux", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            particleFlux, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "momentumFlux", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            momentumFlux, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "heatFlux", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            heatFlux, dimForSpecies, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "NTV", H5T_NATIVE_DOUBLE, dspaceIDForSpecies, &
            dsetID, HDF5Error) 

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            NTV, dimForSpecies, HDF5Error) 

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "jHat", H5T_NATIVE_DOUBLE, dspaceIDForThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            jHat, dimForThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "FSABjHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            FSABjHat, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "Phi1Hat", H5T_NATIVE_DOUBLE, dspaceIDForThetaZeta, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            Phi1Hat, dimForThetaZeta, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "elapsed time (s)", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            elapsedTime, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "didItConverge", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            didItConverge, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "integerToRepresentTrue", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            integerToRepresentTrue, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "integerToRepresentFalse", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            integerToRepresentFalse, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "VPrimeHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            VPrimeHat, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "FSABHat2", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
            dsetID, HDF5Error)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            FSABHat2, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

       call h5dcreate_f(HDF5FileID, "useIterativeSolver", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID, HDF5Error)

       if (useIterativeSolver) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       call h5dclose_f(dsetID, HDF5Error)

       ! ----------------------------------

    end if

    call closeOutputFile()

  end subroutine writeOutputFile

  ! -----------------------------------------------------------------------------------

  subroutine closeOutputFile()

    implicit none

    if (masterProc) then

       call h5sclose_f(dspaceIDForTheta, HDF5Error)
       call h5sclose_f(dspaceIDForZeta, HDF5Error)
       call h5sclose_f(dspaceIDForThetaZeta, HDF5Error)
       call h5sclose_f(dspaceIDForSpeciesThetaZeta, HDF5Error)
       call h5sclose_f(dspaceIDForSources, HDF5Error)
       call h5sclose_f(dspaceIDForScalar, HDF5Error)
       call h5sclose_f(dspaceIDForSpecies, HDF5Error)

       call h5fclose_f(HDF5FileID, HDF5Error)
       call h5close_f(HDF5Error)

    end if

  end subroutine closeOutputFile

  ! -----------------------------------------------------------------------------------

  subroutine saveInputFileToHDF5

    implicit none

! If the file size is larger than this, it will be truncated in the HDF5 output file:
#define maxInputFileSize 99999
    character(maxInputFileSize) :: fileContents

    character(100) :: filename
    character(1) :: oneCharacter
    integer :: fileunit, didFileAccessWork, fileSize, numRecords, ios
    integer :: numBytesRead, filePosition, iFileLine, rank
    PetscErrorCode :: ierr
    integer(SIZE_T) :: fileSizeCopy
    integer(HID_T) :: dspaceIDForInputNamelist
    integer(HID_T) :: dtypeID_inputNamelist
    integer(HID_T) :: dsetID



    filename = inputFilename

    ! Note: this subroutine is only run by the master proc.

    ! Read input file into a character array.
    ! This requires several steps.
    fileUnit=11
    open(unit=fileUnit, file=filename,    action="read", status="old", iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
       print *, "Error opening input file ", trim(filename)
       stop
    end if
    
    ! Determine how large the input.namelist file is:
    fileSize = 0
    numRecords = 0
    ! A fortran "record" is one line of the file.
    do
       read (unit=fileUnit, fmt="(a)", advance="no", &
            iostat=ios) oneCharacter
       if (is_iostat_eor(ios)) then
          numRecords = numRecords + 1
          cycle
       else if (is_iostat_end(ios)) then
          exit
       else
          fileSize = fileSize + 1
       end if
    end do
    
    ! For each record, we add a newline:
    fileSize = fileSize + numRecords
    
    if (fileSize > maxInputFileSize) then
       print *,"WARNING: Input file is very large, so only the beginning of it will be stored in the HDF5 output file."
       fileSize = maxInputFileSize
    end if
    
    rewind(unit=fileUnit)
    filePosition = 1
    do iFileLine = 1,numRecords
       read (unit=fileUnit,fmt="(a)",advance="no",iostat=ios, size=numBytesRead) fileContents(filePosition:fileSize)
       filePosition = filePosition + numBytesRead + 1
       ! Insert newline between records:
       fileContents(filePosition-1:filePosition-1) = achar(10)
       if (filePosition>fileSize) then
          exit
       end if
    end do
    
    close(unit = fileUnit)

    ! Done reading the file. Now begin the HDF5 commands.

    ! Create a HDF5 type corresponding to a string of the appropriate length:
    call h5tcopy_f(H5T_FORTRAN_S1, dtypeID_inputNamelist, HDF5Error)
    fileSizeCopy = fileSize  ! This line explicitly converts to the necessary type.
    call h5tset_size_f(dtypeID_inputNamelist, fileSizeCopy, HDF5Error)
    
    rank = 1
    call h5screate_simple_f(rank, dimForScalar, dspaceIDForInputNamelist, HDF5Error)

    call h5dcreate_f(HDF5FileID, "input.namelist", dtypeID_inputNamelist, dspaceIDForInputnamelist, &
         dsetID, HDF5Error)

    call h5dwrite_f(dsetID, dtypeID_inputNamelist, fileContents(1:fileSize), dimForScalar, HDF5Error)

    ! Destroy HDF5 objects
    call h5dclose_f(dsetID, HDF5Error)
    call h5tclose_f(dtypeID_inputNamelist, HDF5Error)
    call h5sclose_f(dspaceIDForInputNamelist, HDF5Error)
  end subroutine saveInputFileToHDF5

end module writeHDF5Output

