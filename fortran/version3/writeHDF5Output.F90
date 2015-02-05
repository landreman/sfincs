
module writeHDF5Output

  use globalVariables
  use petscsysdef
  use HDF5

  implicit none

#include <finclude/petscsysdef.h>

  private
  public :: initializeOutputFile, updateOutputFile, finalizeHDF5

  integer :: HDF5Error
  integer(HID_T) :: HDF5FileID

  integer(HSIZE_T), dimension(1), parameter :: dimForScalar = 1
  integer(HSIZE_T), dimension(1) :: dimForSpecies
  integer(HSIZE_T), dimension(1) :: dimForTheta
  integer(HSIZE_T), dimension(1) :: dimForZeta
  integer(HSIZE_T), dimension(1) :: dimForx
  integer(HSIZE_T), dimension(2) :: dimForThetaZeta
  integer(HSIZE_T), dimension(2) :: dimForTransportMatrix

  integer(HID_T) :: dspaceIDForScalar
  integer(HID_T) :: dspaceIDForSpecies
  integer(HID_T) :: dspaceIDForTheta
  integer(HID_T) :: dspaceIDForZeta
  integer(HID_T) :: dspaceIDForx
  integer(HID_T) :: dspaceIDForThetaZeta
  integer(HID_T) :: dspaceIDForTransportMatrix

  ! Dimension arrays related to arrays that expand with each iteration:
  integer(HSIZE_T), dimension(1) :: dimForIteration
  integer(HSIZE_T), dimension(1) :: maxDimForIteration
  integer(HSIZE_T), dimension(1) :: dimForIterationChunk
  integer(HID_T) :: pForIteration
  integer(HID_T) :: dspaceIDForIteration

  integer(HSIZE_T), dimension(2) :: dimForIterationSpecies
  integer(HSIZE_T), dimension(2) :: maxDimForIterationSpecies
  integer(HSIZE_T), dimension(2) :: dimForIterationSpeciesChunk
  integer(HID_T) :: pForIterationSpecies
  integer(HID_T) :: dspaceIDForIterationSpecies

  integer(HSIZE_T), dimension(3) :: dimForIterationSpeciesSources
  integer(HSIZE_T), dimension(3) :: maxDimForIterationSpeciesSources
  integer(HSIZE_T), dimension(3) :: dimForIterationSpeciesSourcesChunk
  integer(HID_T) :: pForIterationSpeciesSources
  integer(HID_T) :: dspaceIDForIterationSpeciesSources

  integer(HSIZE_T), dimension(3) :: dimForIterationThetaZeta
  integer(HSIZE_T), dimension(3) :: maxDimForIterationThetaZeta
  integer(HSIZE_T), dimension(3) :: dimForIterationThetaZetaChunk
  integer(HID_T) :: pForIterationThetaZeta
  integer(HID_T) :: dspaceIDForIterationThetaZeta

  integer(HSIZE_T), dimension(4) :: dimForIterationSpeciesThetaZeta
  integer(HSIZE_T), dimension(4) :: maxDimForIterationSpeciesThetaZeta
  integer(HSIZE_T), dimension(4) :: dimForIterationSpeciesThetaZetaChunk
  integer(HID_T) :: pForIterationSpeciesThetaZeta
  integer(HID_T) :: dspaceIDForIterationSpeciesThetaZeta

  interface writeHDF5Field
     module procedure writeHDF5Integer
     module procedure writeHDF5Integers
     module procedure writeHDF5Double
     module procedure writeHDF5Doubles
     module procedure writeHDF5Doubles2
     module procedure writeHDF5Boolean
  end interface writeHDF5Field

  interface writeHDF5ExtensibleField
     module procedure writeHDF5ExtensibleField1
     module procedure writeHDF5ExtensibleField2
     module procedure writeHDF5ExtensibleField3
     module procedure writeHDF5ExtensibleField4
  end interface writeHDF5ExtensibleField

  integer, parameter :: ARRAY_ITERATION = 100
  integer, parameter :: ARRAY_ITERATION_SPECIES = 101
  integer, parameter :: ARRAY_ITERATION_THETA_ZETA = 102
  integer, parameter :: ARRAY_ITERATION_SPECIES_THETA_ZETA = 103
  integer, parameter :: ARRAY_ITERATION_SPECIES_SOURCES = 104

contains


  ! -----------------------------------------------------------------------------------

  subroutine initializeOutputFile()

    ! This subroutine does several things:
    ! 1. The output .h5 file is created,
    ! 2. The dataspace arrays are created,
    ! 3. Variables that do not change with each iteration of SNES are saved to the .h5 file.

    implicit none

    integer(HID_T) :: dsetID
    integer :: temp
    PetscErrorCode :: ierr

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

       call createHDF5Structures()

       call saveInputFileToHDF5()

       ! Because writeHDF5Field is an interface, you can call writeHDF5Field with both single values
       ! and arrays. For arrays, you must provide a dspaceID and dimensions; for scalars you do not
       ! provide these 2 parameters.

       call writeHDF5Field("RHSMode", RHSMode)
       call writeHDF5Field("Nspecies", Nspecies)
       call writeHDF5Field("Ntheta", Ntheta)
       call writeHDF5Field("Nzeta", Nzeta)
       call writeHDF5Field("Nxi", Nxi)
       call writeHDF5Field("NL", NL)
       call writeHDF5Field("Nx", Nx)
       call writeHDF5Field("NxPotentialsPerVth", NxPotentialsPerVth)
       call writeHDF5Field("xMax", xMax)
       call writeHDF5Field("solverTolerance", solverTolerance)
       call writeHDF5Field("theta", theta, dspaceIDForTheta, dimForTheta)
       call writeHDF5Field("zeta", zeta, dspaceIDForZeta, dimForZeta)
       call writeHDF5Field("x", x, dspaceIDForX, dimForX)
       call writeHDF5Field("thetaDerivativeScheme", thetaDerivativeScheme)
       call writeHDF5Field("zetaDerivativeScheme", zetaDerivativeScheme)
       call writeHDF5Field("preconditioner_species", preconditioner_species)
       call writeHDF5Field("preconditioner_x", preconditioner_x)
       call writeHDF5Field("preconditioner_x_min_L", preconditioner_x_min_L)
       call writeHDF5Field("preconditioner_xi", preconditioner_xi)
       call writeHDF5Field("preconditioner_theta", preconditioner_theta)
       call writeHDF5Field("preconditioner_zeta", preconditioner_zeta)
       call writeHDF5Field("constraintScheme", constraintScheme)

       call writeHDF5Field("DHat", DHat, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("BHat", BHat, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHatdpsiHat", dBHatdpsiHat, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHatdtheta", dBHatdtheta, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHatdzeta", dBHatdzeta, dspaceIDForThetaZeta, dimForThetaZeta)

       call writeHDF5Field("BHat_sub_psi", BHat_sub_psi, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sub_psi_dtheta", dBHat_sub_psi_dtheta, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sub_psi_dzeta", dBHat_sub_psi_dzeta, dspaceIDForThetaZeta, dimForThetaZeta)

       call writeHDF5Field("BHat_sub_theta", BHat_sub_theta, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sub_theta_dpsiHat", dBHat_sub_theta_dpsiHat, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sub_theta_dzeta", dBHat_sub_theta_dzeta, dspaceIDForThetaZeta, dimForThetaZeta)

       call writeHDF5Field("BHat_sub_zeta", BHat_sub_zeta, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sub_zeta_dpsiHat", dBHat_sub_zeta_dpsiHat, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sub_zeta_dtheta", dBHat_sub_zeta_dtheta, dspaceIDForThetaZeta, dimForThetaZeta)

       call writeHDF5Field("BHat_sup_theta", BHat_sup_theta, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sup_theta_dpsiHat", dBHat_sup_theta_dpsiHat, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sup_theta_dzeta", dBHat_sup_theta_dzeta, dspaceIDForThetaZeta, dimForThetaZeta)

       call writeHDF5Field("BHat_sup_zeta", BHat_sup_zeta, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sup_zeta_dpsiHat", dBHat_sup_zeta_dpsiHat, dspaceIDForThetaZeta, dimForThetaZeta)
       call writeHDF5Field("dBHat_sup_zeta_dtheta", dBHat_sup_zeta_dtheta, dspaceIDForThetaZeta, dimForThetaZeta)

       call writeHDF5Field("B0OverBBar", B0OverBBar)
       call writeHDF5Field("inputRadialCoordinate", inputRadialCoordinate)
       call writeHDF5Field("inputRadialCoordinateForGradients", inputRadialCoordinateForGradients)

       call writeHDF5Field("psiHat", psiHat)
       call writeHDF5Field("psiN", psiN)
       call writeHDF5Field("rHat", rHat)
       call writeHDF5Field("rN", rN)

       call writeHDF5Field("aHat", aHat)
       call writeHDF5Field("psiAHat", psiAHat)
       call writeHDF5Field("GHat", GHat)
       call writeHDF5Field("IHat", IHat)
       call writeHDF5Field("iota", iota)
       call writeHDF5Field("coordinateSystem", coordinateSystem)

       if (geometryScheme==1) then
          call writeHDF5Field("epsilon_t", epsilon_t)
          call writeHDF5Field("epsilon_h", epsilon_h)
          call writeHDF5Field("epsilon_antisymm", epsilon_antisymm)
          call writeHDF5Field("helicity_n", helicity_n)
          call writeHDF5Field("helicity_l", helicity_l)
          call writeHDF5Field("helicity_antisymm_n", helicity_antisymm_n)
          call writeHDF5Field("helicity_antisymm_l", helicity_antisymm_l)
       end if
       call writeHDF5Field("NPeriods", NPeriods)
       call writeHDF5Field("Delta", Delta)
       call writeHDF5Field("alpha", alpha)
       call writeHDF5Field("nu_n", nu_n)
       call writeHDF5Field("EParallelHat", EParallelHat)

       call writeHDF5Field("collisionOperator", collisionOperator)
       call writeHDF5Field("Zs", Zs, dspaceIDForSpecies, dimForSpecies)
       call writeHDF5Field("mHats", mHats, dspaceIDForSpecies, dimForSpecies)
       call writeHDF5Field("THats", THats, dspaceIDForSpecies, dimForSpecies)
       call writeHDF5Field("nHats", nHats, dspaceIDForSpecies, dimForSpecies)

       call writeHDF5Field("dPhiHatdpsiHat", dPhiHatdpsiHat)
       call writeHDF5Field("dPhiHatdpsiN", dPhiHatdpsiN)
       call writeHDF5Field("dPhiHatdrHat", dPhiHatdrHat)
       call writeHDF5Field("dPhiHatdrN", dPhiHatdrN)

       call writeHDF5Field("dTHatdpsiHat", dTHatdpsiHats, dspaceIDForSpecies, dimForSpecies)
       call writeHDF5Field("dTHatdpsiN", dTHatdpsiNs, dspaceIDForSpecies, dimForSpecies)
       call writeHDF5Field("dTHatdrHat", dTHatdrHats, dspaceIDForSpecies, dimForSpecies)
       call writeHDF5Field("dTHatdrN", dTHatdrNs, dspaceIDForSpecies, dimForSpecies)

       call writeHDF5Field("dnHatdpsiHat", dnHatdpsiHats, dspaceIDForSpecies, dimForSpecies)
       call writeHDF5Field("dnHatdpsiN", dnHatdpsiNs, dspaceIDForSpecies, dimForSpecies)
       call writeHDF5Field("dnHatdrHat", dnHatdrHats, dspaceIDForSpecies, dimForSpecies)
       call writeHDF5Field("dnHatdrN", dnHatdrNs, dspaceIDForSpecies, dimForSpecies)

       call writeHDF5Field("includeTemperatureEquilibrationTerm", includeTemperatureEquilibrationTerm)
       call writeHDF5Field("include_fDivVE_Term", include_fDivVE_Term)
       call writeHDF5Field("includeXDotTerm", includeXDotTerm)
       call writeHDF5Field("includeElectricFieldTermInXiDot", includeElectricFieldTermInXiDot)
       call writeHDF5Field("useDKESExBDrift", useDKESExBDrift)
       call writeHDF5Field("includePhi1", includePhi1)
       call writeHDF5Field("integerToRepresentTrue", integerToRepresentTrue)
       call writeHDF5Field("integerToRepresentFalse", integerToRepresentFalse)
       call writeHDF5Field("VPrimeHat", VPrimeHat)
       call writeHDF5Field("FSABHat2", FSABHat2)
       call writeHDF5Field("useIterativeSolver", useIterativeSolver)
       call writeHDF5Field("NIterations", 0)

       ! ----------------------------------------------------------------------
       ! ----------------------------------------------------------------------

       call h5sclose_f(dspaceIDForTheta, HDF5Error)
       call h5sclose_f(dspaceIDForZeta, HDF5Error)
       call h5sclose_f(dspaceIDForThetaZeta, HDF5Error)
       call h5sclose_f(dspaceIDForScalar, HDF5Error)
       call h5sclose_f(dspaceIDForSpecies, HDF5Error)

       call h5fclose_f(HDF5FileID, HDF5Error)

       ! ----------------------------------

    end if

    ! The next line is not strictly needed, but I included it to be safe.
    ! This way we are sure the HDF5 file is safely closed before the other procs
    ! move on to the next part of the calculation.
    call MPI_Barrier(MPIComm, ierr)

  end subroutine initializeOutputFile

  ! -----------------------------------------------------------------------------------

  subroutine updateOutputFile(iterationNum, writeTransportMatrix)

    ! For an example of how to create an extendible array in HDF5, see
    ! http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/fortran/examples/h5_extend.f90

    implicit none

    integer, intent(in) :: iterationNum
    logical, intent(in) :: writeTransportMatrix
    PetscErrorCode :: ierr
    integer(HID_T) :: dsetID
    integer :: rank

    if (masterProc) then

       print *,"Saving diagnostics to h5 file for iteration ",iterationNum

       call h5fopen_f(outputFilename, H5F_ACC_RDWR_F, HDF5FileID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening HDF5 output file."
          stop
       end if

       ! Over-write the previous value for NIterations:
       call h5dopen_f(HDF5FileID, "NIterations", dsetID, HDF5Error)
       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, iterationNum, dimForScalar, HDF5Error)
       call h5dclose_f(dsetID, HDF5Error)

       dimForIteration(1) = iterationNum
       dimForIterationSpecies(1) = iterationNum
       dimForIterationThetaZeta(1) = iterationNum
       dimForIterationSpeciesThetaZeta(1) = iterationNum
       dimForIterationSpeciesSources(1) = iterationNum

       call writeHDF5ExtensibleField(iterationNum, "densityPerturbation", densityPerturbation, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "totalDensity", totalDensity, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "flow", flow, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "velocityUsingFSADensity", velocityUsingFSADensity, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "velocityUsingTotalDensity", velocityUsingTotalDensity, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "MachUsingFSAThermalSpeed", MachUsingFSAThermalSpeed, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "pressurePerturbation", pressurePerturbation, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "totalPressure", totalPressure, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "particleFluxBeforeSurfaceIntegral_vm0", particleFluxBeforeSurfaceIntegral_vm0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "particleFluxBeforeSurfaceIntegral_vm", particleFluxBeforeSurfaceIntegral_vm, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "particleFluxBeforeSurfaceIntegral_vE0", particleFluxBeforeSurfaceIntegral_vE0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "particleFluxBeforeSurfaceIntegral_vE", particleFluxBeforeSurfaceIntegral_vE, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "momentumFluxBeforeSurfaceIntegral_vm0", momentumFluxBeforeSurfaceIntegral_vm0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "momentumFluxBeforeSurfaceIntegral_vm", momentumFluxBeforeSurfaceIntegral_vm, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "momentumFluxBeforeSurfaceIntegral_vE0", momentumFluxBeforeSurfaceIntegral_vE0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "momentumFluxBeforeSurfaceIntegral_vE", momentumFluxBeforeSurfaceIntegral_vE, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "heatFluxBeforeSurfaceIntegral_vm0", heatFluxBeforeSurfaceIntegral_vm0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "heatFluxBeforeSurfaceIntegral_vm", heatFluxBeforeSurfaceIntegral_vm, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "heatFluxBeforeSurfaceIntegral_vE0", heatFluxBeforeSurfaceIntegral_vE0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "heatFluxBeforeSurfaceIntegral_vE", heatFluxBeforeSurfaceIntegral_vE, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "NTVBeforeSurfaceIntegral", NTVBeforeSurfaceIntegral, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA)

       call writeHDF5ExtensibleField(iterationNum, "FSADensityPerturbation", FSADensityPerturbation, ARRAY_ITERATION_SPECIES)

       call writeHDF5ExtensibleField(iterationNum, "FSABFlow", FSABFlow, ARRAY_ITERATION_SPECIES)

       call writeHDF5ExtensibleField(iterationNum, "FSABVelocityUsingFSADensity", &
            FSABVelocityUsingFSADensity, ARRAY_ITERATION_SPECIES)

       call writeHDF5ExtensibleField(iterationNum, "FSABVelocityUsingFSADensityOverB0", &
            FSABVelocityUsingFSADensityOverB0, ARRAY_ITERATION_SPECIES)

       call writeHDF5ExtensibleField(iterationNum, "FSABVelocityUsingFSADensityOverRootFSAB2", &
            FSABVelocityUsingFSADensityOverRootFSAB2, ARRAY_ITERATION_SPECIES)

       call writeHDF5ExtensibleField(iterationNum, "FSAPressurePerturbation", FSAPressurePerturbation, ARRAY_ITERATION_SPECIES)

       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm0_psiHat", particleFlux_vm0_psiHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm0_psiN", particleFlux_vm0_psiN, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm0_rHat", particleFlux_vm0_rHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm0_rN", particleFlux_vm0_rN, ARRAY_ITERATION_SPECIES)

       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm_psiHat", particleFlux_vm_psiHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm_psiN", particleFlux_vm_psiN, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm_rHat", particleFlux_vm_rHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm_rN", particleFlux_vm_rN, ARRAY_ITERATION_SPECIES)

       if (includePhi1) then
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE0_psiHat", particleFlux_vE0_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE0_psiN", particleFlux_vE0_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE0_rHat", particleFlux_vE0_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE0_rN", particleFlux_vE0_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE_psiHat", particleFlux_vE_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE_psiN", particleFlux_vE_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE_rHat", particleFlux_vE_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE_rN", particleFlux_vE_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd1_psiHat", particleFlux_vd1_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd1_psiN", particleFlux_vd1_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd1_rHat", particleFlux_vd1_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd1_rN", particleFlux_vd1_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd_psiHat", particleFlux_vd_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd_psiN", particleFlux_vd_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd_rHat", particleFlux_vd_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd_rN", particleFlux_vd_rN, ARRAY_ITERATION_SPECIES)
       end if

       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm0_psiHat", momentumFlux_vm0_psiHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm0_psiN", momentumFlux_vm0_psiN, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm0_rHat", momentumFlux_vm0_rHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm0_rN", momentumFlux_vm0_rN, ARRAY_ITERATION_SPECIES)

       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm_psiHat", momentumFlux_vm_psiHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm_psiN", momentumFlux_vm_psiN, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm_rHat", momentumFlux_vm_rHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm_rN", momentumFlux_vm_rN, ARRAY_ITERATION_SPECIES)

       if (includePhi1) then
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE0_psiHat", momentumFlux_vE0_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE0_psiN", momentumFlux_vE0_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE0_rHat", momentumFlux_vE0_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE0_rN", momentumFlux_vE0_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE_psiHat", momentumFlux_vE_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE_psiN", momentumFlux_vE_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE_rHat", momentumFlux_vE_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE_rN", momentumFlux_vE_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd1_psiHat", momentumFlux_vd1_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd1_psiN", momentumFlux_vd1_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd1_rHat", momentumFlux_vd1_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd1_rN", momentumFlux_vd1_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd_psiHat", momentumFlux_vd_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd_psiN", momentumFlux_vd_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd_rHat", momentumFlux_vd_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd_rN", momentumFlux_vd_rN, ARRAY_ITERATION_SPECIES)
       end if

       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm0_psiHat", heatFlux_vm0_psiHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm0_psiN", heatFlux_vm0_psiN, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm0_rHat", heatFlux_vm0_rHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm0_rN", heatFlux_vm0_rN, ARRAY_ITERATION_SPECIES)

       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm_psiHat", heatFlux_vm_psiHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm_psiN", heatFlux_vm_psiN, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm_rHat", heatFlux_vm_rHat, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm_rN", heatFlux_vm_rN, ARRAY_ITERATION_SPECIES)

       if (includePhi1) then
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE0_psiHat", heatFlux_vE0_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE0_psiN", heatFlux_vE0_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE0_rHat", heatFlux_vE0_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE0_rN", heatFlux_vE0_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE_psiHat", heatFlux_vE_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE_psiN", heatFlux_vE_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE_rHat", heatFlux_vE_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE_rN", heatFlux_vE_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd1_psiHat", heatFlux_vd1_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd1_psiN", heatFlux_vd1_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd1_rHat", heatFlux_vd1_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd1_rN", heatFlux_vd1_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd_psiHat", heatFlux_vd_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd_psiN", heatFlux_vd_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd_rHat", heatFlux_vd_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd_rN", heatFlux_vd_rN, ARRAY_ITERATION_SPECIES)

          call writeHDF5ExtensibleField(iterationNum, "heatFlux_withoutPhi1_psiHat", heatFlux_withoutPhi1_psiHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_withoutPhi1_psiN", heatFlux_withoutPhi1_psiN, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_withoutPhi1_rHat", heatFlux_withoutPhi1_rHat, ARRAY_ITERATION_SPECIES)
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_withoutPhi1_rN", heatFlux_withoutPhi1_rN, ARRAY_ITERATION_SPECIES)
       end if

       call writeHDF5ExtensibleField(iterationNum, "NTV", NTV, ARRAY_ITERATION_SPECIES)
       call writeHDF5ExtensibleField(iterationNum, "jHat", jHat, ARRAY_ITERATION_THETA_ZETA)
       call writeHDF5ExtensibleField(iterationNum, "FSABjHat", FSABjHat, ARRAY_ITERATION)
       call writeHDF5ExtensibleField(iterationNum, "FSABjHatOverB0", FSABjHatOverB0, ARRAY_ITERATION)
       call writeHDF5ExtensibleField(iterationNum, "FSABjHatOverRootFSAB2", FSABjHatOverRootFSAB2, ARRAY_ITERATION)

       if (includePhi1) then
          call writeHDF5ExtensibleField(iterationNum, "Phi1Hat", Phi1Hat, ARRAY_ITERATION_THETA_ZETA)
          call writeHDF5ExtensibleField(iterationNum, "d(Phi1Hat)d(theta)", dPhi1Hatdtheta, ARRAY_ITERATION_THETA_ZETA)
          call writeHDF5ExtensibleField(iterationNum, "d(Phi1Hat)d(zeta)", dPhi1Hatdzeta, ARRAY_ITERATION_THETA_ZETA)
          call writeHDF5ExtensibleField(iterationNum, "lambda", lambda, ARRAY_ITERATION)
       end if

       call writeHDF5ExtensibleField(iterationNum, "elapsed time (s)", elapsedTime, ARRAY_ITERATION)

       if (constraintScheme .ne. 0) then
          call writeHDF5ExtensibleField(iterationNum, "sources", sources, ARRAY_ITERATION_SPECIES_SOURCES)
       end if

!!$       ! ----------------------------------
!!$
!!$       call h5dcreate_f(HDF5FileID, "didNonlinearCalculationConverge", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
!!$            dsetID, HDF5Error)
!!$
!!$       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, &
!!$            didNonlinearCalculationConverge, dimForScalar, HDF5Error)
!!$


       ! Save transport matrix if needed
       if (writeTransportMatrix) then
          rank = 2
          dimForTransportMatrix(1) = transportMatrixSize
          dimForTransportMatrix(2) = transportMatrixSize
          call h5screate_simple_f(rank, dimForTransportMatrix, dspaceIDForTransportMatrix, HDF5Error)

          call writeHDF5Field("transportMatrix", transportMatrix, dspaceIDForTransportMatrix, dimForTransportMatrix)

          call h5sclose_f(dspaceIDForTransportMatrix, HDF5Error)
       end if

       call h5fclose_f(HDF5FileID, HDF5Error)
    end if

    ! The next line is not strictly needed, but I included it to be safe.
    ! This way we are sure the HDF5 file is safely closed before the other procs
    ! move on to the next part of the calculation.
    call MPI_Barrier(MPIComm, ierr)

  end subroutine updateOutputFile

  ! -----------------------------------------------------------------------------------

  subroutine createHDF5Structures()

    ! This subroutine creates the dataspaces.

    implicit none

    integer :: i, rank

    ! Create a dataspace for storing single numbers:
    rank = 0
    call h5screate_simple_f(rank, dimForScalar, dspaceIDForScalar, HDF5Error)

    ! Create dataspaces that depend on resolution parameters
    ! but which do not expand with the number of SNES iterations:
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
    dimForThetaZeta(2) = Nzeta
    call h5screate_simple_f(rank, dimForThetaZeta, dspaceIDForThetaZeta, HDF5Error)

    ! ------------------------------------------------------------------
    ! Next come the arrays that expand with each iteration of SNES.
    ! These are more complicated.
    ! ------------------------------------------------------------------

    rank = 1
    dimForIteration(1)      = 1
    maxDimForIteration(1)   = H5S_UNLIMITED_F
    dimForIterationChunk(1) = 1

    call h5screate_simple_f(rank, dimForIteration, dspaceIDForIteration, &
         HDF5Error, maxDimForIteration)
    call h5pcreate_f(H5P_DATASET_CREATE_F, pForIteration, HDF5Error)
    call h5pset_chunk_f(pForIteration, rank, dimForIterationChunk, HDF5Error)

    ! -------------------------------------

    rank = 2
    dimForIterationSpecies(1)      = 1
    maxDimForIterationSpecies(1)   = H5S_UNLIMITED_F
    dimForIterationSpeciesChunk(1) = 1

    dimForIterationSpecies(2)      = Nspecies
    maxDimForIterationSpecies(2)   = Nspecies
    dimForIterationSpeciesChunk(2) = Nspecies

    call h5screate_simple_f(rank, dimForIterationSpecies, dspaceIDForIterationSpecies, &
         HDF5Error, maxDimForIterationSpecies)
    call h5pcreate_f(H5P_DATASET_CREATE_F, pForIterationSpecies, HDF5Error)
    call h5pset_chunk_f(pForIterationSpecies, rank, dimForIterationSpeciesChunk, HDF5Error)

    ! -------------------------------------

    rank = 3
    dimForIterationThetaZeta(1)      = 1
    maxDimForIterationThetaZeta(1)   = H5S_UNLIMITED_F
    dimForIterationThetaZetaChunk(1) = 1

    dimForIterationThetaZeta(2)      = Ntheta
    maxDimForIterationThetaZeta(2)   = Ntheta
    dimForIterationThetaZetaChunk(2) = Ntheta

    dimForIterationThetaZeta(3)      = Nzeta
    maxDimForIterationThetaZeta(3)   = Nzeta
    dimForIterationThetaZetaChunk(3) = Nzeta

    call h5screate_simple_f(rank, dimForIterationThetaZeta, dspaceIDForIterationThetaZeta, &
         HDF5Error, maxDimForIterationThetaZeta)
    call h5pcreate_f(H5P_DATASET_CREATE_F, pForIterationThetaZeta, HDF5Error)
    call h5pset_chunk_f(pForIterationThetaZeta, rank, dimForIterationThetaZetaChunk, HDF5Error)

    ! -------------------------------------

    rank = 3
    dimForIterationSpeciesSources(1)      = 1
    maxDimForIterationSpeciesSources(1)   = H5S_UNLIMITED_F
    dimForIterationSpeciesSourcesChunk(1) = 1

    dimForIterationSpeciesSources(2)      = Nspecies
    maxDimForIterationSpeciesSources(2)   = Nspecies
    dimForIterationSpeciesSourcesChunk(2) = Nspecies

    if (constraintScheme==1) then
       dimForIterationSpeciesSources(3)      = 2
       maxDimForIterationSpeciesSources(3)   = 2
       dimForIterationSpeciesSourcesChunk(3) = 2
    else
       ! This block of code will also be run when constraintScheme==0,
       ! but it does not matter since we will not write sources to the file
       ! when constraintScheme==0
       dimForIterationSpeciesSources(3)      = Nx
       maxDimForIterationSpeciesSources(3)   = Nx
       dimForIterationSpeciesSourcesChunk(3) = Nx
    end if

    call h5screate_simple_f(rank, dimForIterationSpeciesSources, dspaceIDForIterationSpeciesSources, &
         HDF5Error, maxDimForIterationSpeciesSources)
    call h5pcreate_f(H5P_DATASET_CREATE_F, pForIterationSpeciesSources, HDF5Error)
    call h5pset_chunk_f(pForIterationSpeciesSources, rank, dimForIterationSpeciesSourcesChunk, HDF5Error)

    ! -------------------------------------

    rank = 4
    dimForIterationSpeciesThetaZeta(1)      = 1
    maxDimForIterationSpeciesThetaZeta(1)   = H5S_UNLIMITED_F
    dimForIterationSpeciesThetaZetaChunk(1) = 1

    dimForIterationSpeciesThetaZeta(2)      = Nspecies
    maxDimForIterationSpeciesThetaZeta(2)   = Nspecies
    dimForIterationSpeciesThetaZetaChunk(2) = Nspecies

    dimForIterationSpeciesThetaZeta(3)      = Ntheta
    maxDimForIterationSpeciesThetaZeta(3)   = Ntheta
    dimForIterationSpeciesThetaZetaChunk(3) = Ntheta

    dimForIterationSpeciesThetaZeta(4)      = Nzeta
    maxDimForIterationSpeciesThetaZeta(4)   = Nzeta
    dimForIterationSpeciesThetaZetaChunk(4) = Nzeta

    call h5screate_simple_f(rank, dimForIterationSpeciesThetaZeta, dspaceIDForIterationSpeciesThetaZeta, &
         HDF5Error, maxDimForIterationSpeciesThetaZeta)
    call h5pcreate_f(H5P_DATASET_CREATE_F, pForIterationSpeciesThetaZeta, HDF5Error)
    call h5pset_chunk_f(pForIterationSpeciesThetaZeta, rank, dimForIterationSpeciesThetaZetaChunk, HDF5Error)

  end subroutine createHDF5Structures

  ! -----------------------------------------------------------------------------------

  subroutine finalizeHDF5()

    implicit none

    integer :: rank

    if (masterProc) then

       ! Re-open the file briefly to add the "finished" variable.
       call h5fopen_f(outputFilename, H5F_ACC_RDWR_F, HDF5FileID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening HDF5 output file."
          stop
       end if

       ! Create a dataspace for storing single numbers:
       rank = 0
       call h5screate_simple_f(rank, dimForScalar, dspaceIDForScalar, HDF5Error)

       call writeHDF5Field("finished", integerToRepresentTrue)

       call h5sclose_f(dspaceIDForScalar, HDF5Error)
       call h5fclose_f(HDF5FileID, HDF5Error)

       ! Done adding the "finished" variable.

       call h5pclose_f(pForIteration, HDF5Error)
       call h5pclose_f(pForIterationSpecies, HDF5Error)
       call h5pclose_f(pForIterationSpeciesThetaZeta, HDF5Error)
       call h5pclose_f(pForIterationThetaZeta, HDF5Error)
       call h5pclose_f(pForIterationSpeciesSources, HDF5Error)

       call h5close_f(HDF5Error)

    end if

  end subroutine finalizeHDF5

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

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Integer(arrayName, data)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer :: data

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_INTEGER, dspaceIDForScalar, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, data, dimForScalar, HDF5Error)
    
    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Integer

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Integers(arrayName, data, dspaceID, dims)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID
    integer(HSIZE_T), dimension(*) :: dims
    integer, dimension(*) :: data

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, data, dims, HDF5Error)
    
    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Integers

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Double(arrayName, data)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    PetscScalar :: data

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, dspaceIDForScalar, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dimForScalar, HDF5Error)
    
    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Double

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Doubles(arrayName, data, dspaceID, dims)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID
    integer(HSIZE_T), dimension(*) :: dims
    PetscScalar, dimension(*) :: data

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dims, HDF5Error)
    
    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Doubles

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Doubles2(arrayName, data, dspaceID, dims)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID
    integer(HSIZE_T), dimension(*) :: dims
    PetscScalar, dimension(:,:) :: data

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dims, HDF5Error)
    
    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Doubles2

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Boolean(arrayName, data)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    logical :: data
    integer :: temp

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_INTEGER, dspaceIDForScalar, dsetID, HDF5Error)
    
    if (data) then
       temp = integerToRepresentTrue
    else
       temp = integerToRepresentFalse
    end if
    call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, temp, dimForScalar, HDF5Error)
    
    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Boolean

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5ExtensibleField1(iterationNum, arrayName, data, arrayType)

    implicit none

    integer, parameter :: rank = 1
    integer, intent(in) :: iterationNum
    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID, memspaceID, originalDspaceID
    integer :: temp, arrayType
    PetscScalar :: data
    integer(HSIZE_T) :: offset(rank)
    integer(HSIZE_T), dimension(rank) :: dim, dimForChunk
    integer(HID_T) :: chunkProperties

    offset = (/ iterationNum-1 /)

    select case (arrayType)
    case (ARRAY_ITERATION)
       originalDspaceID = dspaceIDForIteration
       dim = dimForIteration
       dimForChunk = dimForIterationChunk
       chunkProperties = pForIteration
    case default
       print *,"This is writeHDF5ExtensibleField1"
       print *,"Error! Invalid arrayType:",arrayType
       stop
    end select

    if (iterationNum==1) then
       call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, originalDspaceID, &
            dsetID, HDF5Error, chunkProperties)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            data, dimForSpecies, HDF5Error)
    else
       ! Extend an existing array in the .h5 file:
       call h5dopen_f(HDF5FileID, arrayName, dsetID, HDF5Error)
       call h5dset_extent_f(dsetID, dim, HDF5Error)
       call h5dget_space_f(dsetID, dspaceID, HDF5Error)
       call h5sselect_hyperslab_f(dspaceID, H5S_SELECT_SET_F, offset, &
            dimForChunk, HDF5Error)
       call h5screate_simple_f(rank, dimForChunk, memspaceID, HDF5Error)
       call H5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dimForChunk, HDF5Error, &
            memspaceID, dspaceID)
       call h5sclose_f(dspaceID, HDF5Error)
       call h5sclose_f(memspaceID, HDF5Error)
    end if

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5ExtensibleField1

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5ExtensibleField2(iterationNum, arrayName, data, arrayType)

    implicit none

    integer, parameter :: rank = 2
    integer, intent(in) :: iterationNum
    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID, memspaceID, originalDspaceID
    integer :: temp, arrayType
    PetscScalar, dimension(:) :: data
    integer(HSIZE_T) :: offset(rank)
    integer(HSIZE_T), dimension(rank) :: dim, dimForChunk
    integer(HID_T) :: chunkProperties

    offset = (/ iterationNum-1, 0 /)

    select case (arrayType)
    case (ARRAY_ITERATION_SPECIES)
       originalDspaceID = dspaceIDForIterationSpecies
       dim = dimForIterationSpecies
       dimForChunk = dimForIterationSpeciesChunk
       chunkProperties = pForIterationSpecies
    case default
       print *,"This is writeHDF5ExtensibleField2"
       print *,"Error! Invalid arrayType:",arrayType
       stop
    end select

    if (iterationNum==1) then
       call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, originalDspaceID, &
            dsetID, HDF5Error, chunkProperties)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            data, dimForSpecies, HDF5Error)
    else
       ! Extend an existing array in the .h5 file:
       call h5dopen_f(HDF5FileID, arrayName, dsetID, HDF5Error)
       call h5dset_extent_f(dsetID, dim, HDF5Error)
       call h5dget_space_f(dsetID, dspaceID, HDF5Error)
       call h5sselect_hyperslab_f(dspaceID, H5S_SELECT_SET_F, offset, &
            dimForChunk, HDF5Error)
       call h5screate_simple_f(rank, dimForChunk, memspaceID, HDF5Error)
       call H5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dimForChunk, HDF5Error, &
            memspaceID, dspaceID)
       call h5sclose_f(dspaceID, HDF5Error)
       call h5sclose_f(memspaceID, HDF5Error)
    end if

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5ExtensibleField2

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5ExtensibleField3(iterationNum, arrayName, data, arrayType)

    implicit none

    integer, parameter :: rank = 3
    integer, intent(in) :: iterationNum
    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID, memspaceID, originalDspaceID
    integer :: temp, arrayType
    PetscScalar, dimension(:,:) :: data
    integer(HSIZE_T) :: offset(rank)
    integer(HSIZE_T), dimension(rank) :: dim, dimForChunk
    integer(HID_T) :: chunkProperties

    offset = (/ iterationNum-1, 0, 0 /)

    select case (arrayType)
    case (ARRAY_ITERATION_THETA_ZETA)
       originalDspaceID = dspaceIDForIterationThetaZeta
       dim = dimForIterationThetaZeta
       dimForChunk = dimForIterationThetaZetaChunk
       chunkProperties = pForIterationThetaZeta
    case (ARRAY_ITERATION_SPECIES_SOURCES)
       originalDspaceID = dspaceIDForIterationSpeciesSources
       dim = dimForIterationSpeciesSources
       dimForChunk = dimForIterationSpeciesSourcesChunk
       chunkProperties = pForIterationSpeciesSources
    case default
       print *,"This is writeHDF5ExtensibleField3"
       print *,"Error! Invalid arrayType:",arrayType
       stop
    end select

    if (iterationNum==1) then
       call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, originalDspaceID, &
            dsetID, HDF5Error, chunkProperties)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            data, dimForSpecies, HDF5Error)
    else
       ! Extend an existing array in the .h5 file:
       call h5dopen_f(HDF5FileID, arrayName, dsetID, HDF5Error)
       call h5dset_extent_f(dsetID, dim, HDF5Error)
       call h5dget_space_f(dsetID, dspaceID, HDF5Error)
       call h5sselect_hyperslab_f(dspaceID, H5S_SELECT_SET_F, offset, &
            dimForChunk, HDF5Error)
       call h5screate_simple_f(rank, dimForChunk, memspaceID, HDF5Error)
       call H5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dimForChunk, HDF5Error, &
            memspaceID, dspaceID)
       call h5sclose_f(dspaceID, HDF5Error)
       call h5sclose_f(memspaceID, HDF5Error)
    end if

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5ExtensibleField3

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5ExtensibleField4(iterationNum, arrayName, data, arrayType)

    implicit none

    integer, parameter :: rank = 4
    integer, intent(in) :: iterationNum
    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID, memspaceID, originalDspaceID
    integer :: temp, arrayType
    PetscScalar, dimension(:,:,:) :: data
    integer(HSIZE_T) :: offset(rank)
    integer(HSIZE_T), dimension(rank) :: dim, dimForChunk
    integer(HID_T) :: chunkProperties

    offset = (/ iterationNum-1, 0, 0, 0 /)

    select case (arrayType)
    case (ARRAY_ITERATION_SPECIES_THETA_ZETA)
       originalDspaceID = dspaceIDForIterationSpeciesThetaZeta
       dim = dimForIterationSpeciesThetaZeta
       dimForChunk = dimForIterationSpeciesThetaZetaChunk
       chunkProperties = pForIterationSpeciesThetaZeta
    case default
       print *,"This is writeHDF5ExtensibleField4"
       print *,"Error! Invalid arrayType:",arrayType
       stop
    end select

    if (iterationNum==1) then
       call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, originalDspaceID, &
            dsetID, HDF5Error, chunkProperties)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            data, dimForSpecies, HDF5Error)
    else
       ! Extend an existing array in the .h5 file:
       call h5dopen_f(HDF5FileID, arrayName, dsetID, HDF5Error)
       call h5dset_extent_f(dsetID, dim, HDF5Error)
       call h5dget_space_f(dsetID, dspaceID, HDF5Error)
       call h5sselect_hyperslab_f(dspaceID, H5S_SELECT_SET_F, offset, &
            dimForChunk, HDF5Error)
       call h5screate_simple_f(rank, dimForChunk, memspaceID, HDF5Error)
       call H5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dimForChunk, HDF5Error, &
            memspaceID, dspaceID)
       call h5sclose_f(dspaceID, HDF5Error)
       call h5sclose_f(memspaceID, HDF5Error)
    end if

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5ExtensibleField4

end module writeHDF5Output

