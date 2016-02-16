
#define attribute_name "Description"
#define boolDescription "(Interpret this Boolean data by comparing to integerToRepresentTrue/False)"

module writeHDF5Output

  use globalVariables
  use petscsysdef
  use HDF5
  use H5DS
  use H5LT

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

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
  integer(HSIZE_T), dimension(1) :: dimForExport_f_theta
  integer(HSIZE_T), dimension(1) :: dimForExport_f_zeta
  integer(HSIZE_T), dimension(1) :: dimForExport_f_xi
  integer(HSIZE_T), dimension(1) :: dimForExport_f_x

  integer(HID_T) :: dspaceIDForScalar
  integer(HID_T) :: dspaceIDForSpecies
  integer(HID_T) :: dspaceIDForTheta
  integer(HID_T) :: dspaceIDForZeta
  integer(HID_T) :: dspaceIDForx
  integer(HID_T) :: dspaceIDForThetaZeta
  integer(HID_T) :: dspaceIDForTransportMatrix
  integer(HID_T) :: dspaceIDForExport_f_theta
  integer(HID_T) :: dspaceIDForExport_f_zeta
  integer(HID_T) :: dspaceIDForExport_f_xi
  integer(HID_T) :: dspaceIDForExport_f_x

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

  integer(HSIZE_T), dimension(3) :: dimForIterationSpeciesX
  integer(HSIZE_T), dimension(3) :: maxDimForIterationSpeciesX
  integer(HSIZE_T), dimension(3) :: dimForIterationSpeciesXChunk
  integer(HID_T) :: pForIterationSpeciesX
  integer(HID_T) :: dspaceIDForIterationSpeciesX

  integer(HSIZE_T), dimension(4) :: dimForIterationSpeciesThetaZeta
  integer(HSIZE_T), dimension(4) :: maxDimForIterationSpeciesThetaZeta
  integer(HSIZE_T), dimension(4) :: dimForIterationSpeciesThetaZetaChunk
  integer(HID_T) :: pForIterationSpeciesThetaZeta
  integer(HID_T) :: dspaceIDForIterationSpeciesThetaZeta

  integer(HSIZE_T), dimension(6) :: dimForExport_f
  integer(HSIZE_T), dimension(6) :: maxDimForExport_f
  integer(HSIZE_T), dimension(6) :: dimForExport_fChunk
  integer(HID_T) :: pForExport_f
  integer(HID_T) :: dspaceIDForExport_f

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
     module procedure writeHDF5ExtensibleField6
  end interface writeHDF5ExtensibleField

  integer, parameter :: ARRAY_ITERATION = 100
  integer, parameter :: ARRAY_ITERATION_SPECIES = 101
  integer, parameter :: ARRAY_ITERATION_THETA_ZETA = 102
  integer, parameter :: ARRAY_ITERATION_SPECIES_THETA_ZETA = 103
  integer, parameter :: ARRAY_ITERATION_SPECIES_SOURCES = 104
  integer, parameter :: ARRAY_EXPORT_F = 105
  integer, parameter :: ARRAY_ITERATION_SPECIES_X = 106

contains


  ! -----------------------------------------------------------------------------------

  subroutine initializeOutputFile()

    ! This subroutine does several things:
    ! 1. The output .h5 file is created,
    ! 2. The dataspace arrays are created,
    ! 3. Variables that do not change with each iteration of SNES are saved to the .h5 file.

    use export_f
    use xGrid, only: xGrid_k
    
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

       call writeHDF5Field("RHSMode", RHSMode, &
            "Switch that controls how many times the kinetic equation is solved with different right-hand-side drive terms.")
       call writeHDF5Field("Nspecies", Nspecies, "Number of particle species")
       call writeHDF5Field("Ntheta", Ntheta, "Number of grid points in the poloidal angle theta")
       call writeHDF5Field("Nzeta", Nzeta, "Number of grid points in the toroidal angle zeta")
       call writeHDF5Field("Nxi", Nxi, &
            "Number of Legendre polynomial modes P(v_parallel / v) for representing the distribution functions.")
       call writeHDF5Field("NL", NL, &
            "Number of Legendre polynomial modes P(v_parallel / v) for representing the Rosenbluth potentials in the Fokker-Planck operator.")
       call writeHDF5Field("Nx", Nx, "Number of grid points in speed for representing the distribution functions.")
       call writeHDF5Field("NxPotentialsPerVth", NxPotentialsPerVth, "")
       call writeHDF5Field("xMax", xMax, "")
       call writeHDF5Field("solverTolerance", solverTolerance, "")
       call writeHDF5Field("theta", theta, dspaceIDForTheta, dimForTheta, &
            "Grid points in the poloidal angle, which runs from 0 to 2pi")
       call writeHDF5Field("zeta", zeta, dspaceIDForZeta, dimForZeta, &
            "Grid points in the toroidal angle, which runs from 0 to 2pi/Nperiods")
       call writeHDF5Field("x", x, dspaceIDForX, dimForX, &
            "Grid points in normalized speed, x_s = v / sqrt{2 T_s / m_s}, the same for each species s.")
       call writeHDF5Field("geometryScheme", geometryScheme, "")
       call writeHDF5Field("thetaDerivativeScheme", thetaDerivativeScheme, "")
       call writeHDF5Field("zetaDerivativeScheme", zetaDerivativeScheme, "")
       call writeHDF5Field("ExBDerivativeScheme", ExBDerivativeScheme, "")
       call writeHDF5Field("magneticDriftDerivativeScheme", magneticDriftDerivativeScheme, "")
       call writeHDF5Field("xGridScheme", xGridScheme, "")
       call writeHDF5Field("xGrid_k", xGrid_k, "Exponent of x in the orthogonality relation for the speed polynomials")
       call writeHDF5Field("xPotentialsGridScheme", xPotentialsGridScheme, "")
       call writeHDF5Field("pointAtX0", pointAtX0, "Does the x grid include a point at x=0? " // boolDescription)
       call writeHDF5Field("preconditioner_species", preconditioner_species, "")
       call writeHDF5Field("preconditioner_x", preconditioner_x, "")
       call writeHDF5Field("preconditioner_x_min_L", preconditioner_x_min_L, "")
       call writeHDF5Field("preconditioner_xi", preconditioner_xi, "")
       call writeHDF5Field("preconditioner_theta", preconditioner_theta, "")
       call writeHDF5Field("preconditioner_zeta", preconditioner_zeta, "")
       call writeHDF5Field("preconditioner_magnetic_drifts_max_L", preconditioner_magnetic_drifts_max_L, "")
       call writeHDF5Field("reusePreconditioner", reusePreconditioner, "Use the same preconditioner matrix at each iteration of the Newton solver? " &
            // boolDescription)
       call writeHDF5Field("constraintScheme", constraintScheme, "")

       call writeHDF5Field("DHat", DHat, dspaceIDForThetaZeta, dimForThetaZeta, &
            "Inverse Jacobian (grad psi dot grad theta cross grad zeta, where 2*pi*psi is the toroidal flux), in units of BBar / RBar")
       call writeHDF5Field("BHat", BHat, dspaceIDForThetaZeta, dimForThetaZeta, &
            "Magnitude of the magnetic field, in units of BBar")
       call writeHDF5Field("dBHatdpsiHat", dBHatdpsiHat, dspaceIDForThetaZeta, dimForThetaZeta, "")
       call writeHDF5Field("dBHatdtheta", dBHatdtheta, dspaceIDForThetaZeta, dimForThetaZeta, "")
       call writeHDF5Field("dBHatdzeta", dBHatdzeta, dspaceIDForThetaZeta, dimForThetaZeta, "")
       call writeHDF5Field("BDotCurlB", BDotCurlB, dspaceIDForThetaZeta, dimForThetaZeta, &
            "\vect{B}\cdot\nabla\times\vect{B}, in units of BBar^2 / RBar")
       call writeHDF5Field("uHat", uHat, dspaceIDForThetaZeta, dimForThetaZeta, &
            "\nabla_\parallel u = (2/B^4)\nabla B\times\vector{B}\cdot\iota\nabla\psi")

       call writeHDF5Field("BHat_sub_psi", BHat_sub_psi, dspaceIDForThetaZeta, dimForThetaZeta, &
            "Covariant component of B, the component that multiplies grad psi, in units of 1 / RBar")
       call writeHDF5Field("dBHat_sub_psi_dtheta", dBHat_sub_psi_dtheta, dspaceIDForThetaZeta, dimForThetaZeta, "")
       call writeHDF5Field("dBHat_sub_psi_dzeta", dBHat_sub_psi_dzeta, dspaceIDForThetaZeta, dimForThetaZeta, "")

       call writeHDF5Field("BHat_sub_theta", BHat_sub_theta, dspaceIDForThetaZeta, dimForThetaZeta, &
            "Covariant component of B, the component that multiplies grad theta, in units of BBar * RBar")
       call writeHDF5Field("dBHat_sub_theta_dpsiHat", dBHat_sub_theta_dpsiHat, dspaceIDForThetaZeta, dimForThetaZeta, "")
       call writeHDF5Field("dBHat_sub_theta_dzeta", dBHat_sub_theta_dzeta, dspaceIDForThetaZeta, dimForThetaZeta, "")

       call writeHDF5Field("BHat_sub_zeta", BHat_sub_zeta, dspaceIDForThetaZeta, dimForThetaZeta, &
            "Covariant component of B, the component that multiplies grad zeta, in units of BBar * RBar")
       call writeHDF5Field("dBHat_sub_zeta_dpsiHat", dBHat_sub_zeta_dpsiHat, dspaceIDForThetaZeta, dimForThetaZeta, "")
       call writeHDF5Field("dBHat_sub_zeta_dtheta", dBHat_sub_zeta_dtheta, dspaceIDForThetaZeta, dimForThetaZeta, "")

       call writeHDF5Field("BHat_sup_theta", BHat_sup_theta, dspaceIDForThetaZeta, dimForThetaZeta, &
            "Contravariant component of B, the component that multiplies d r / d theta, in units of BBar / RBar")
       call writeHDF5Field("dBHat_sup_theta_dpsiHat", dBHat_sup_theta_dpsiHat, dspaceIDForThetaZeta, dimForThetaZeta, "")
       call writeHDF5Field("dBHat_sup_theta_dzeta", dBHat_sup_theta_dzeta, dspaceIDForThetaZeta, dimForThetaZeta, "")

       call writeHDF5Field("BHat_sup_zeta", BHat_sup_zeta, dspaceIDForThetaZeta, dimForThetaZeta, &
            "Contravariant component of B, the component that multiplies d r / d zeta, in units of BBar / RBar")
       call writeHDF5Field("dBHat_sup_zeta_dpsiHat", dBHat_sup_zeta_dpsiHat, dspaceIDForThetaZeta, dimForThetaZeta, "")
       call writeHDF5Field("dBHat_sup_zeta_dtheta", dBHat_sup_zeta_dtheta, dspaceIDForThetaZeta, dimForThetaZeta, "")

       call writeHDF5Field("B0OverBBar", B0OverBBar, &
            "m=0, n=0 harmonic of |B| in Boozer coordinates, equivalent to <B^3>/<B^2>, in units of BBar")
       call writeHDF5Field("inputRadialCoordinate", inputRadialCoordinate, "")
       call writeHDF5Field("inputRadialCoordinateForGradients", inputRadialCoordinateForGradients, "")

       call writeHDF5Field("psiHat", psiHat, "")
       call writeHDF5Field("psiN", psiN, "")
       call writeHDF5Field("rHat", rHat, "")
       call writeHDF5Field("rN", rN, "")

       call writeHDF5Field("aHat", aHat, "")
       call writeHDF5Field("psiAHat", psiAHat, "")
       call writeHDF5Field("GHat", GHat, "")
       call writeHDF5Field("IHat", IHat, "")
       call writeHDF5Field("iota", iota, "(Rationalized) rotational transform = 1 / (safety factor q)")
       call writeHDF5Field("coordinateSystem", coordinateSystem, "")
       call writeHDF5Field("magneticDriftScheme", magneticDriftScheme, "Which version of the poloidal and toroidal magnetic drifts to use.")
       call writeHDF5Field("force0RadialCurrentInEquilibrium", force0RadialCurrentInEquilibrium,&
            "If true, assume dBHat_sub_zeta_dtheta - dBHat_sub_theta_dzeta = 0, since this relation is implied by the MHD equilibrium relation "//&
            "curl(B) dot grad psi = 0." //&
            "If false, allow dBHat_sub_zeta_dtheta - dBHat_sub_theta_dzeta to be nonzero. " // boolDescription)

       if (geometryScheme==1) then
          call writeHDF5Field("epsilon_t", epsilon_t, "")
          call writeHDF5Field("epsilon_h", epsilon_h, "")
          call writeHDF5Field("epsilon_antisymm", epsilon_antisymm, "")
          call writeHDF5Field("helicity_n", helicity_n, "")
          call writeHDF5Field("helicity_l", helicity_l, "")
          call writeHDF5Field("helicity_antisymm_n", helicity_antisymm_n, "")
          call writeHDF5Field("helicity_antisymm_l", helicity_antisymm_l, "")
       end if
       call writeHDF5Field("NPeriods", NPeriods, "Number of identical toroidal periods (e.g. 5 for W7-X, 10 for LHD, 4 for HSX)")
       call writeHDF5Field("Delta", Delta, &
            "Dimensionless combination of the normalization constants, resembling rho_*: Delta = mBar * vBar / (e * BBar * RBar) (SI units) " // &
            "or c * mBar * vBar / (e * BBar * RBar) (Gaussian units)")
       call writeHDF5Field("alpha", alpha, "Dimensionless combination of the normalization constants: alpha = e * PhiBar / TBar.")
       call writeHDF5Field("nu_n", nu_n, "")
       call writeHDF5Field("EParallelHat", EParallelHat, "")
       if (RHSMode==3) then
          call writeHDF5Field("nuPrime", nuPrime, "")
          call writeHDF5Field("EStar", EStar, "")
       end if
       call writeHDF5Field("collisionOperator", collisionOperator, "")
       call writeHDF5Field("Zs", Zs, dspaceIDForSpecies, dimForSpecies, "Charge of each species, in units of the unit charge e (which is usually the proton charge.)")
       call writeHDF5Field("mHats", mHats, dspaceIDForSpecies, dimForSpecies, "Mass of each species, in units of mBar.")
       call writeHDF5Field("THats", THats, dspaceIDForSpecies, dimForSpecies, "Average temperature of each species, in units of TBar.")
       call writeHDF5Field("nHats", nHats, dspaceIDForSpecies, dimForSpecies, "Flux surface averaged density of each species, in units of nBar.")


       !!Added by AM 2016-01!!
       call writeHDF5Field("withAdiabatic", withAdiabatic, "")
       if (withAdiabatic) then
       	  call writeHDF5Field("adiabaticZ", adiabaticZ, "Charge of adiabatic species, in units of the unit charge e (which is usually the proton charge.)")
	  call writeHDF5Field("adiabaticMHat", adiabaticMHat, "Mass of adiabatic species, in units of mBar.")
	  call writeHDF5Field("adiabaticNHat", adiabaticNHat, "Flux surface averaged density of adiabatic species, in units of nBar.")
	  call writeHDF5Field("adiabaticTHat", adiabaticTHat, "Average temperature of adiabatic species, in units of TBar.")
       end if
       !!!!!!!!!!!!!!!!!!!!!!!

       call writeHDF5Field("dPhiHatdpsiHat", dPhiHatdpsiHat, "")
       call writeHDF5Field("dPhiHatdpsiN", dPhiHatdpsiN, "")
       call writeHDF5Field("dPhiHatdrHat", dPhiHatdrHat, "")
       call writeHDF5Field("dPhiHatdrN", dPhiHatdrN, "")

       call writeHDF5Field("dTHatdpsiHat", dTHatdpsiHats, dspaceIDForSpecies, dimForSpecies, "")
       call writeHDF5Field("dTHatdpsiN", dTHatdpsiNs, dspaceIDForSpecies, dimForSpecies, "")
       call writeHDF5Field("dTHatdrHat", dTHatdrHats, dspaceIDForSpecies, dimForSpecies, "")
       call writeHDF5Field("dTHatdrN", dTHatdrNs, dspaceIDForSpecies, dimForSpecies, "")

       call writeHDF5Field("dnHatdpsiHat", dnHatdpsiHats, dspaceIDForSpecies, dimForSpecies, "")
       call writeHDF5Field("dnHatdpsiN", dnHatdpsiNs, dspaceIDForSpecies, dimForSpecies, "")
       call writeHDF5Field("dnHatdrHat", dnHatdrHats, dspaceIDForSpecies, dimForSpecies, "")
       call writeHDF5Field("dnHatdrN", dnHatdrNs, dspaceIDForSpecies, dimForSpecies, "")

       call writeHDF5Field("includeTemperatureEquilibrationTerm", includeTemperatureEquilibrationTerm, &
            "Include the inhomogeneous term associated with the collision operator acting on the Maxwellians C[f_M, f_M]? " //&
            "This term is nonzero only when the Fokker-Planck operator is used with unequal temperatures. "// boolDescription)
       call writeHDF5Field("include_fDivVE_Term", include_fDivVE_Term, "")
       call writeHDF5Field("includeXDotTerm", includeXDotTerm, "")
       call writeHDF5Field("includeElectricFieldTermInXiDot", includeElectricFieldTermInXiDot, "")
       call writeHDF5Field("useDKESExBDrift", useDKESExBDrift, "")
       call writeHDF5Field("includePhi1", includePhi1, &
            "Include a quasineutrality equation, and include variation of the electrostatic potential on a flux surface? " // boolDescription)
       call writeHDF5Field("includeRadialExBDrive", includeRadialExBDrive, &
            "Include term $(\vect{v}_{E} \cdot\nabla\psi)f_{Ms} [(1/n_s)(dn_s/d\psi) + (x_s^2-3/2)(1/T_s)(dT_s/d\psi)]$ term? " // boolDescription)
       call writeHDF5Field("integerToRepresentTrue", integerToRepresentTrue, &
            "Since HDF5 does not have a Boolean datatype, this integer value is used in this file for Boolean quantities.")
       call writeHDF5Field("integerToRepresentFalse", integerToRepresentFalse, &
            "Since HDF5 does not have a Boolean datatype, this integer value is used in this file for Boolean quantities.")
       call writeHDF5Field("VPrimeHat", VPrimeHat, "")
       call writeHDF5Field("FSABHat2", FSABHat2, &
            "< B^2 >, the flux-surface-averaged squared magnitude of the magnetic field, in units of BBar^2")
       call writeHDF5Field("useIterativeLinearSolver", useIterativeLinearSolver, "")
       call writeHDF5Field("NIterations", 0, "")

       call writeHDF5Field("export_full_f",  export_full_f,  "Save the full f distribution function in this file? " // boolDescription)
       call writeHDF5Field("export_delta_f", export_delta_f, "Save the delta f distribution function in this file? " // boolDescription)
       if (export_full_f .or. export_delta_f) then
          call writeHDF5Field("export_f_theta_option", export_f_theta_option, &
               "Which theta grid to use for exporting the distribution function.")
          call writeHDF5Field("export_f_zeta_option", export_f_zeta_option, &
               "Which zeta grid to use for exporting the distribution function.")
          call writeHDF5Field("export_f_xi_option", export_f_xi_option, &
               "Which xi discretization and grid to use for exporting the distribution function.")
          call writeHDF5Field("export_f_x_option", export_f_x_option, &
               "Which x grid to use for exporting the distribution function.")

          call writeHDF5Field("N_export_f_theta", N_export_f_theta, &
               "Size of export_f_theta, i.e. the number of theta values on which the distribution function is saved")
          call writeHDF5Field("N_export_f_zeta", N_export_f_zeta, &
               "Size of export_f_zeta, i.e. the number of zeta values on which the distribution function is saved")
          if (export_f_xi_option >0) then
             call writeHDF5Field("N_export_f_xi", N_export_f_xi, &
                  "Size of export_f_xi, i.e. the number of xi values on which the distribution function is saved")
          end if
          call writeHDF5Field("N_export_f_x", N_export_f_x, &
               "Size of export_f_x, i.e. the number of x values on which the distribution function is saved")

          call writeHDF5Field("export_f_theta", export_f_theta, dspaceIDForExport_f_theta, dimForExport_f_theta, &
               "Values of the poloidal angle theta for which the distribution functions delta f or full f are saved")
          call writeHDF5Field("export_f_zeta", export_f_zeta, dspaceIDForExport_f_zeta, dimForExport_f_zeta, &
               "Values of the toroidal angle zeta for which the distribution functions delta f or full f are saved")
          if (export_f_xi_option > 0) then
             call writeHDF5Field("export_f_xi", export_f_xi, dspaceIDForExport_f_xi, dimForExport_f_xi, &
                  "Values of cos(pitch angle) xi for which the distribution functions delta f or full f are saved")
          end if
          call writeHDF5Field("export_f_x", export_f_x, dspaceIDForExport_f_x, dimForExport_f_x, &
               "Values of normalized speed x for which the distribution functions delta f or full f are saved")
       end if

       ! ----------------------------------------------------------------------
       ! ----------------------------------------------------------------------

       call h5sclose_f(dspaceIDForTheta, HDF5Error)
       call h5sclose_f(dspaceIDForZeta, HDF5Error)
       call h5sclose_f(dspaceIDForThetaZeta, HDF5Error)
       call h5sclose_f(dspaceIDForScalar, HDF5Error)
       call h5sclose_f(dspaceIDForSpecies, HDF5Error)
       call h5sclose_f(dspaceIDForx, HDF5Error)
       call h5sclose_f(dspaceIDForExport_f_theta, HDF5Error)
       call h5sclose_f(dspaceIDForExport_f_zeta, HDF5Error)
       call h5sclose_f(dspaceIDForExport_f_xi, HDF5Error)
       call h5sclose_f(dspaceIDForExport_f_x, HDF5Error)

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

    use export_f, only : full_f, delta_f, export_full_f, export_delta_f

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
       dimForIterationSpeciesX(1) = iterationNum
       dimForExport_f(1) = iterationNum

       call writeHDF5ExtensibleField(iterationNum, "densityPerturbation", densityPerturbation, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, &
            "Variation of the density over the flux surface, subtracting the flux surface average, in units of nBar.")

       call writeHDF5ExtensibleField(iterationNum, "totalDensity", totalDensity, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, &
            "Density of each species, including both the average and the variation on a flux surface, in units of nBar.")

       call writeHDF5ExtensibleField(iterationNum, "flow", flow, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "velocityUsingFSADensity", velocityUsingFSADensity, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "velocityUsingTotalDensity", velocityUsingTotalDensity, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "MachUsingFSAThermalSpeed", MachUsingFSAThermalSpeed, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "pressurePerturbation", pressurePerturbation, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "pressureAnisotropy", pressureAnisotropy, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "totalPressure", totalPressure, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "particleFluxBeforeSurfaceIntegral_vm0", particleFluxBeforeSurfaceIntegral_vm0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "particleFluxBeforeSurfaceIntegral_vm", particleFluxBeforeSurfaceIntegral_vm, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "particleFluxBeforeSurfaceIntegral_vE0", particleFluxBeforeSurfaceIntegral_vE0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "particleFluxBeforeSurfaceIntegral_vE", particleFluxBeforeSurfaceIntegral_vE, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "momentumFluxBeforeSurfaceIntegral_vm0", momentumFluxBeforeSurfaceIntegral_vm0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "momentumFluxBeforeSurfaceIntegral_vm", momentumFluxBeforeSurfaceIntegral_vm, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "momentumFluxBeforeSurfaceIntegral_vE0", momentumFluxBeforeSurfaceIntegral_vE0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "momentumFluxBeforeSurfaceIntegral_vE", momentumFluxBeforeSurfaceIntegral_vE, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "heatFluxBeforeSurfaceIntegral_vm0", heatFluxBeforeSurfaceIntegral_vm0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "heatFluxBeforeSurfaceIntegral_vm", heatFluxBeforeSurfaceIntegral_vm, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "heatFluxBeforeSurfaceIntegral_vE0", heatFluxBeforeSurfaceIntegral_vE0, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "heatFluxBeforeSurfaceIntegral_vE", heatFluxBeforeSurfaceIntegral_vE, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "NTVBeforeSurfaceIntegral", NTVBeforeSurfaceIntegral, &
            ARRAY_ITERATION_SPECIES_THETA_ZETA, "")

       call writeHDF5ExtensibleField(iterationNum, "FSADensityPerturbation", FSADensityPerturbation, ARRAY_ITERATION_SPECIES, &
            "Flux-surface-averaged density, minus the requested average density, for each species. Should be within machine precision of 0.")

       call writeHDF5ExtensibleField(iterationNum, "FSABFlow", FSABFlow, ARRAY_ITERATION_SPECIES, "")

       call writeHDF5ExtensibleField(iterationNum, "FSABFlow_vs_x", FSABFlow_vs_x, ARRAY_ITERATION_SPECIES_X, "")

       call writeHDF5ExtensibleField(iterationNum, "FSABVelocityUsingFSADensity", &
            FSABVelocityUsingFSADensity, ARRAY_ITERATION_SPECIES, "")

       call writeHDF5ExtensibleField(iterationNum, "FSABVelocityUsingFSADensityOverB0", &
            FSABVelocityUsingFSADensityOverB0, ARRAY_ITERATION_SPECIES, "")

       call writeHDF5ExtensibleField(iterationNum, "FSABVelocityUsingFSADensityOverRootFSAB2", &
            FSABVelocityUsingFSADensityOverRootFSAB2, ARRAY_ITERATION_SPECIES, "")

       call writeHDF5ExtensibleField(iterationNum, "FSAPressurePerturbation", FSAPressurePerturbation, ARRAY_ITERATION_SPECIES, &
            "Flux-surface-averaged pressure, minus the requested average pressure, for each species. Should be within machine precision of 0.")

       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm0_psiHat", particleFlux_vm0_psiHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm0_psiN", particleFlux_vm0_psiN, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm0_rHat", particleFlux_vm0_rHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm0_rN", particleFlux_vm0_rN, ARRAY_ITERATION_SPECIES, "")

       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm_psiHat", particleFlux_vm_psiHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm_psiN", particleFlux_vm_psiN, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm_rHat", particleFlux_vm_rHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm_rN", particleFlux_vm_rN, ARRAY_ITERATION_SPECIES, "")

       call writeHDF5ExtensibleField(iterationNum, "particleFlux_vm_psiHat_vs_x", particleFlux_vm_psiHat_vs_x, ARRAY_ITERATION_SPECIES_X, "")

       if (includePhi1) then
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE0_psiHat", particleFlux_vE0_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE0_psiN", particleFlux_vE0_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE0_rHat", particleFlux_vE0_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE0_rN", particleFlux_vE0_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE_psiHat", particleFlux_vE_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE_psiN", particleFlux_vE_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE_rHat", particleFlux_vE_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vE_rN", particleFlux_vE_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd1_psiHat", particleFlux_vd1_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd1_psiN", particleFlux_vd1_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd1_rHat", particleFlux_vd1_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd1_rN", particleFlux_vd1_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd_psiHat", particleFlux_vd_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd_psiN", particleFlux_vd_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd_rHat", particleFlux_vd_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "particleFlux_vd_rN", particleFlux_vd_rN, ARRAY_ITERATION_SPECIES, "")
       end if

       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm0_psiHat", momentumFlux_vm0_psiHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm0_psiN", momentumFlux_vm0_psiN, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm0_rHat", momentumFlux_vm0_rHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm0_rN", momentumFlux_vm0_rN, ARRAY_ITERATION_SPECIES, "")

       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm_psiHat", momentumFlux_vm_psiHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm_psiN", momentumFlux_vm_psiN, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm_rHat", momentumFlux_vm_rHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vm_rN", momentumFlux_vm_rN, ARRAY_ITERATION_SPECIES, "")

       if (includePhi1) then
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE0_psiHat", momentumFlux_vE0_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE0_psiN", momentumFlux_vE0_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE0_rHat", momentumFlux_vE0_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE0_rN", momentumFlux_vE0_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE_psiHat", momentumFlux_vE_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE_psiN", momentumFlux_vE_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE_rHat", momentumFlux_vE_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vE_rN", momentumFlux_vE_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd1_psiHat", momentumFlux_vd1_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd1_psiN", momentumFlux_vd1_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd1_rHat", momentumFlux_vd1_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd1_rN", momentumFlux_vd1_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd_psiHat", momentumFlux_vd_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd_psiN", momentumFlux_vd_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd_rHat", momentumFlux_vd_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "momentumFlux_vd_rN", momentumFlux_vd_rN, ARRAY_ITERATION_SPECIES, "")
       end if

       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm0_psiHat", heatFlux_vm0_psiHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm0_psiN", heatFlux_vm0_psiN, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm0_rHat", heatFlux_vm0_rHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm0_rN", heatFlux_vm0_rN, ARRAY_ITERATION_SPECIES, "")

       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm_psiHat", heatFlux_vm_psiHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm_psiN", heatFlux_vm_psiN, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm_rHat", heatFlux_vm_rHat, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm_rN", heatFlux_vm_rN, ARRAY_ITERATION_SPECIES, "")

       call writeHDF5ExtensibleField(iterationNum, "heatFlux_vm_psiHat_vs_x", heatFlux_vm_psiHat_vs_x, ARRAY_ITERATION_SPECIES_X, "")

       if (includePhi1) then
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE0_psiHat", heatFlux_vE0_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE0_psiN", heatFlux_vE0_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE0_rHat", heatFlux_vE0_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE0_rN", heatFlux_vE0_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE_psiHat", heatFlux_vE_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE_psiN", heatFlux_vE_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE_rHat", heatFlux_vE_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vE_rN", heatFlux_vE_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd1_psiHat", heatFlux_vd1_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd1_psiN", heatFlux_vd1_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd1_rHat", heatFlux_vd1_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd1_rN", heatFlux_vd1_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd_psiHat", heatFlux_vd_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd_psiN", heatFlux_vd_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd_rHat", heatFlux_vd_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_vd_rN", heatFlux_vd_rN, ARRAY_ITERATION_SPECIES, "")

          call writeHDF5ExtensibleField(iterationNum, "heatFlux_withoutPhi1_psiHat", heatFlux_withoutPhi1_psiHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_withoutPhi1_psiN", heatFlux_withoutPhi1_psiN, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_withoutPhi1_rHat", heatFlux_withoutPhi1_rHat, ARRAY_ITERATION_SPECIES, "")
          call writeHDF5ExtensibleField(iterationNum, "heatFlux_withoutPhi1_rN", heatFlux_withoutPhi1_rN, ARRAY_ITERATION_SPECIES, "")
       end if

       call writeHDF5ExtensibleField(iterationNum, "NTV", NTV, ARRAY_ITERATION_SPECIES, "")
       call writeHDF5ExtensibleField(iterationNum, "jHat", jHat, ARRAY_ITERATION_THETA_ZETA, &
            "Parallel current j dot B / |B|, in units e * nBar * vBar")
       call writeHDF5ExtensibleField(iterationNum, "FSABjHat", FSABjHat, ARRAY_ITERATION, &
            "Flux surface averaged parallel current, <j dot B>, in units e * nBar * vBar * BBar")
       call writeHDF5ExtensibleField(iterationNum, "FSABjHatOverB0", FSABjHatOverB0, ARRAY_ITERATION, "")
       call writeHDF5ExtensibleField(iterationNum, "FSABjHatOverRootFSAB2", FSABjHatOverRootFSAB2, ARRAY_ITERATION, "")

       if (includePhi1) then
          call writeHDF5ExtensibleField(iterationNum, "Phi1Hat", Phi1Hat, ARRAY_ITERATION_THETA_ZETA, &
               "Electrostatic potential Phi minus its flux-surface-average, in units of PhiBar")
          call writeHDF5ExtensibleField(iterationNum, "dPhi1Hatdtheta", dPhi1Hatdtheta, ARRAY_ITERATION_THETA_ZETA, &
               "Derivative of Phi_1 with respect to theta. Phi_1 = Electrostatic potential Phi minus its flux-surface-average, in units of PhiBar." //&
               "theta = poloidal angle")
          call writeHDF5ExtensibleField(iterationNum, "dPhi1Hatdzeta", dPhi1Hatdzeta, ARRAY_ITERATION_THETA_ZETA,  &
               "Derivative of Phi_1 with respect to zeta. Phi_1 = Electrostatic potential Phi minus its flux-surface-average, in units of PhiBar." //&
               "zeta = toroidal angle")
          call writeHDF5ExtensibleField(iterationNum, "lambda", lambda, ARRAY_ITERATION, &
               "Lagrange multiplier associated with the constraint that <Phi_1>=0. Should be within machine precision of 0.")
       end if

       call writeHDF5ExtensibleField(iterationNum, "elapsed time (s)", elapsedTime, ARRAY_ITERATION, "")

       if (export_full_f) then
          call writeHDF5ExtensibleField(iterationNum,"full_f", full_f, ARRAY_EXPORT_F, &
               "Full distribution function for each species, normalized by nBar / (vBar ^3)")
       end if
       if (export_delta_f) then
          call writeHDF5ExtensibleField(iterationNum,"delta_f", delta_f, ARRAY_EXPORT_F, &
               "Distribution function for each species, normalized by nBar / (vBar ^3), subtracting the leading-order Maxwellian")
       end if

       if (constraintScheme .ne. 0) then
          call writeHDF5ExtensibleField(iterationNum, "sources", sources, ARRAY_ITERATION_SPECIES_SOURCES, "")
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

          call writeHDF5Field("transportMatrix", transportMatrix, dspaceIDForTransportMatrix, dimForTransportMatrix, &
               "")

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

    use export_f, only: N_export_f_theta, N_export_f_zeta, N_export_f_xi, N_export_f_x

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

    dimForExport_f_theta = N_export_f_theta
    call h5screate_simple_f(rank, dimForExport_f_theta, dspaceIDForExport_f_theta, HDF5Error)

    dimForExport_f_zeta = N_export_f_zeta
    call h5screate_simple_f(rank, dimForExport_f_zeta, dspaceIDForExport_f_zeta, HDF5Error)

    dimForExport_f_xi = N_export_f_xi
    call h5screate_simple_f(rank, dimForExport_f_xi, dspaceIDForExport_f_xi, HDF5Error)

    dimForExport_f_x = N_export_f_x
    call h5screate_simple_f(rank, dimForExport_f_x, dspaceIDForExport_f_x, HDF5Error)

    rank = 2
    dimForThetaZeta(1) = Ntheta
    dimForThetaZeta(2) = Nzeta
    call h5screate_simple_f(rank, dimForThetaZeta, dspaceIDForThetaZeta, HDF5Error)

    rank = 5
    dimForExport_f(1) = Nspecies
    dimForExport_f(2) = N_export_f_theta
    dimForExport_f(3) = N_export_f_zeta
    dimForExport_f(4) = N_export_f_xi
    dimForExport_f(5) = N_export_f_x
    call h5screate_simple_f(rank, dimForExport_f, dspaceIDForExport_f, HDF5Error)

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

    select case (constraintScheme)
    case (0)
       ! No action needed since we will not write sources to the file
       ! when constraintScheme==0
    case (1,3,4)
       dimForIterationSpeciesSources(3)      = 2
       maxDimForIterationSpeciesSources(3)   = 2
       dimForIterationSpeciesSourcesChunk(3) = 2
    case (2)
       dimForIterationSpeciesSources(3)      = Nx
       maxDimForIterationSpeciesSources(3)   = Nx
       dimForIterationSpeciesSourcesChunk(3) = Nx
    case default
       stop "Invalid constraintScheme!"
    end select

    call h5screate_simple_f(rank, dimForIterationSpeciesSources, dspaceIDForIterationSpeciesSources, &
         HDF5Error, maxDimForIterationSpeciesSources)
    call h5pcreate_f(H5P_DATASET_CREATE_F, pForIterationSpeciesSources, HDF5Error)
    call h5pset_chunk_f(pForIterationSpeciesSources, rank, dimForIterationSpeciesSourcesChunk, HDF5Error)

    ! -------------------------------------

    rank = 3
    dimForIterationSpeciesX(1)      = 1
    maxDimForIterationSpeciesX(1)   = H5S_UNLIMITED_F
    dimForIterationSpeciesXChunk(1) = 1

    dimForIterationSpeciesX(2)      = Nspecies
    maxDimForIterationSpeciesX(2)   = Nspecies
    dimForIterationSpeciesXChunk(2) = Nspecies

    dimForIterationSpeciesX(3)      = Nx
    maxDimForIterationSpeciesX(3)   = Nx
    dimForIterationSpeciesXChunk(3) = Nx

    call h5screate_simple_f(rank, dimForIterationSpeciesX, dspaceIDForIterationSpeciesX, &
         HDF5Error, maxDimForIterationSpeciesX)
    call h5pcreate_f(H5P_DATASET_CREATE_F, pForIterationSpeciesX, HDF5Error)
    call h5pset_chunk_f(pForIterationSpeciesX, rank, dimForIterationSpeciesXChunk, HDF5Error)

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

    ! -------------------------------------

    rank = 6
    dimForExport_f(1)      = 1
    maxDimForExport_f(1)   = H5S_UNLIMITED_F
    dimForExport_fChunk(1) = 1

    dimForExport_f(2)      = Nspecies
    maxDimForExport_f(2)   = Nspecies
    dimForExport_fChunk(2) = Nspecies

    dimForExport_f(3)      = N_export_f_theta
    maxDimForExport_f(3)   = N_export_f_theta
    dimForExport_fChunk(3) = N_export_f_theta

    dimForExport_f(4)      = N_export_f_zeta
    maxDimForExport_f(4)   = N_export_f_zeta
    dimForExport_fChunk(4) = N_export_f_zeta

    dimForExport_f(5)      = N_export_f_xi
    maxDimForExport_f(5)   = N_export_f_xi
    dimForExport_fChunk(5) = N_export_f_xi

    dimForExport_f(6)      = N_export_f_x
    maxDimForExport_f(6)   = N_export_f_x
    dimForExport_fChunk(6) = N_export_f_x

    call h5screate_simple_f(rank, dimForExport_f, dspaceIDForExport_f, &
         HDF5Error, maxDimForExport_f)
    call h5pcreate_f(H5P_DATASET_CREATE_F, pForExport_f, HDF5Error)
    call h5pset_chunk_f(pForExport_f, rank, dimForExport_fChunk, HDF5Error)

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

       call writeHDF5Field("finished", integerToRepresentTrue, &
            "If this variable exists in sfincsOutput.h5, then SFINCS reached the end of all requested computations and exited gracefully.")

       call h5sclose_f(dspaceIDForScalar, HDF5Error)
       call h5fclose_f(HDF5FileID, HDF5Error)

       ! Done adding the "finished" variable.

       call h5pclose_f(pForIteration, HDF5Error)
       call h5pclose_f(pForIterationSpecies, HDF5Error)
       call h5pclose_f(pForIterationSpeciesThetaZeta, HDF5Error)
       call h5pclose_f(pForIterationThetaZeta, HDF5Error)
       call h5pclose_f(pForIterationSpeciesSources, HDF5Error)
       call h5pclose_f(pForIterationSpeciesX, HDF5Error)
       call h5pclose_f(pForExport_f, HDF5Error)

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

  subroutine writeHDF5Integer(arrayName, data, description)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer :: data
    character(len=*) :: description

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_INTEGER, dspaceIDForScalar, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, data, dimForScalar, HDF5Error)
    
    call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Integer

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Integers(arrayName, data, dspaceID, dims, description)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID
    integer(HSIZE_T), dimension(*) :: dims
    integer, dimension(*) :: data
    character(len=*) :: description

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, data, dims, HDF5Error)
    
    call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Integers

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Double(arrayName, data, description)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    PetscScalar :: data
    character(len=*) :: description

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, dspaceIDForScalar, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dimForScalar, HDF5Error)
    
    call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Double

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Doubles(arrayName, data, dspaceID, dims, description)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID
    integer(HSIZE_T), dimension(*) :: dims
    character(len=*) :: description
    PetscScalar, dimension(*) :: data

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dims, HDF5Error)
    
    if (dspaceID == dspaceIDForSpecies) then
       call h5dsset_label_f(dsetID, 1, "species", HDF5Error)
    elseif (dspaceID == dspaceIDForTheta) then
       ! No labels applied in this case.
    elseif (dspaceID == dspaceIDForZeta) then
       ! No labels applied in this case.
    elseif (dspaceID == dspaceIDForx) then
       ! No labels applied in this case.
    elseif (dspaceID == dspaceIDForExport_f_theta) then
       ! No labels applied in this case.
    elseif (dspaceID == dspaceIDForExport_f_zeta) then
       ! No labels applied in this case.
    elseif (dspaceID == dspaceIDForExport_f_xi) then
       ! No labels applied in this case.
    elseif (dspaceID == dspaceIDForExport_f_x) then
       ! No labels applied in this case.
    else
       print *,"WARNING: PROGRAM SHOULD NOT GET HERE. (writeHDF5Doubles)"
    end if

    call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Doubles

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Doubles2(arrayName, data, dspaceID, dims, description)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID
    integer(HSIZE_T), dimension(*) :: dims
    PetscScalar, dimension(:,:) :: data
    character(len=*) :: description
    character(len=100) :: label

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dims, HDF5Error)

    if (dspaceID == dspaceIDForThetaZeta) then
       call h5dsset_label_f(dsetID, 1, "zeta", HDF5Error)
       call h5dsset_label_f(dsetID, 2, "theta", HDF5Error)
    elseif (dspaceID == dspaceIDForTransportMatrix) then
       ! No labels applied in this case.
    else
       print *,"WARNING: PROGRAM SHOULD NOT GET HERE. (writeHDF5Doubles2)"
    end if
    
    call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Doubles2

  ! -----------------------------------------------------------------------------------

!!$  subroutine writeHDF5Doubles5(arrayName, data, dspaceID, dims, description, iteration)
!!$
!!$    use export_f, only: export_f_xi_option
!!$
!!$    implicit none
!!$
!!$    character(len=*) :: arrayName
!!$    integer(HID_T) :: dsetID
!!$    integer(HID_T) :: dspaceID
!!$    integer(HSIZE_T), dimension(*) :: dims
!!$    PetscScalar, dimension(:,:,:,:,:) :: data
!!$    character(len=*) :: description
!!$    character(len=100) :: label
!!$    integer :: iteration
!!$
!!$    if (iteration==1) then
!!$       call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
!!$    else
!!$       call h5dopen_f(HDF5FileID, arrayName, dsetID, HDF5Error)
!!$    end if
!!$
!!$    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, data, dims, HDF5Error)
!!$
!!$    if (dspaceID == dspaceIDForExport_f) then
!!$       call h5dsset_label_f(dsetID, 1, "export_f_x", HDF5Error)
!!$       if (export_f_xi_option==0) then
!!$          call h5dsset_label_f(dsetID, 2, "Legendre mode number in xi, i.e. n in P_n(xi). First index is P_0(xi)=1.", HDF5Error)
!!$       else
!!$          call h5dsset_label_f(dsetID, 2, "export_f_xi", HDF5Error)
!!$       end if
!!$       call h5dsset_label_f(dsetID, 3, "export_f_zeta", HDF5Error)
!!$       call h5dsset_label_f(dsetID, 4, "export_f_theta", HDF5Error)
!!$       call h5dsset_label_f(dsetID, 5, "species", HDF5Error)
!!$    else
!!$       print *,"WARNING: PROGRAM SHOULD NOT GET HERE. (writeHDF5Doubles5)"
!!$    end if
!!$    
!!$    call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)
!!$
!!$    call h5dclose_f(dsetID, HDF5Error)
!!$
!!$  end subroutine writeHDF5Doubles5

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5Boolean(arrayName, data, description)

    implicit none

    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    logical :: data
    character(len=*) :: description
    integer :: temp

    call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_INTEGER, dspaceIDForScalar, dsetID, HDF5Error)
    
    if (data) then
       temp = integerToRepresentTrue
    else
       temp = integerToRepresentFalse
    end if
    call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, temp, dimForScalar, HDF5Error)
    
    call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)

    call h5dclose_f(dsetID, HDF5Error)

  end subroutine writeHDF5Boolean

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5ExtensibleField1(iterationNum, arrayName, data, arrayType, description)

    implicit none

    integer, parameter :: rank = 1
    integer, intent(in) :: iterationNum
    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID, memspaceID, originalDspaceID
    integer :: temp, arrayType
    character(len=*) :: description
    PetscScalar :: data
    integer(HSIZE_T) :: offset(rank)
    integer(HSIZE_T), dimension(rank) :: dim, dimForChunk
    integer(HID_T) :: chunkProperties
    character(len=100) :: label1

    offset = (/ iterationNum-1 /)

    select case (arrayType)
    case (ARRAY_ITERATION)
       originalDspaceID = dspaceIDForIteration
       dim = dimForIteration
       dimForChunk = dimForIterationChunk
       chunkProperties = pForIteration
       label1 = "iteration"
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

       call h5dsset_label_f(dsetID, 1, trim(label1), HDF5Error)
       call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)
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

  subroutine writeHDF5ExtensibleField2(iterationNum, arrayName, data, arrayType, description)

    implicit none

    integer, parameter :: rank = 2
    integer, intent(in) :: iterationNum
    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID, memspaceID, originalDspaceID
    integer :: temp, arrayType
    character(len=*) :: description
    PetscScalar, dimension(:) :: data
    integer(HSIZE_T) :: offset(rank)
    integer(HSIZE_T), dimension(rank) :: dim, dimForChunk
    integer(HID_T) :: chunkProperties
    character(len=100) :: label1, label2

    offset = (/ iterationNum-1, 0 /)

    select case (arrayType)
    case (ARRAY_ITERATION_SPECIES)
       originalDspaceID = dspaceIDForIterationSpecies
       dim = dimForIterationSpecies
       dimForChunk = dimForIterationSpeciesChunk
       chunkProperties = pForIterationSpecies
       label1 = "species"
       label2 = "iteration"
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

       call h5dsset_label_f(dsetID, 1, trim(label1), HDF5Error)
       call h5dsset_label_f(dsetID, 2, trim(label2), HDF5Error)

       call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)
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

  subroutine writeHDF5ExtensibleField3(iterationNum, arrayName, data, arrayType, description)

    implicit none

    integer, parameter :: rank = 3
    integer, intent(in) :: iterationNum
    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID, memspaceID, originalDspaceID
    integer :: temp, arrayType
    character(len=*) :: description
    PetscScalar, dimension(:,:) :: data
    integer(HSIZE_T) :: offset(rank)
    integer(HSIZE_T), dimension(rank) :: dim, dimForChunk
    integer(HID_T) :: chunkProperties
    character(len=100) :: label1, label2, label3

    offset = (/ iterationNum-1, 0, 0 /)

    select case (arrayType)
    case (ARRAY_ITERATION_THETA_ZETA)
       originalDspaceID = dspaceIDForIterationThetaZeta
       dim = dimForIterationThetaZeta
       dimForChunk = dimForIterationThetaZetaChunk
       chunkProperties = pForIterationThetaZeta
       label1 = "zeta"
       label2 = "theta"
       label3 = "iteration"
    case (ARRAY_ITERATION_SPECIES_SOURCES)
       originalDspaceID = dspaceIDForIterationSpeciesSources
       dim = dimForIterationSpeciesSources
       dimForChunk = dimForIterationSpeciesSourcesChunk
       chunkProperties = pForIterationSpeciesSources
       label1 = "sources"
       label2 = "species"
       label3 = "iteration"
    case (ARRAY_ITERATION_SPECIES_X)
       originalDspaceID = dspaceIDForIterationSpeciesX
       dim = dimForIterationSpeciesX
       dimForChunk = dimForIterationSpeciesXChunk
       chunkProperties = pForIterationSpeciesX
       label1 = "x"
       label2 = "species"
       label3 = "iteration"
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

       call h5dsset_label_f(dsetID, 1, trim(label1), HDF5Error)
       call h5dsset_label_f(dsetID, 2, trim(label2), HDF5Error)
       call h5dsset_label_f(dsetID, 3, trim(label3), HDF5Error)

       call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)
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

  subroutine writeHDF5ExtensibleField4(iterationNum, arrayName, data, arrayType, description)

    implicit none

    integer, parameter :: rank = 4
    integer, intent(in) :: iterationNum
    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID, memspaceID, originalDspaceID
    integer :: temp, arrayType
    character(len=*) :: description
    PetscScalar, dimension(:,:,:) :: data
    integer(HSIZE_T) :: offset(rank)
    integer(HSIZE_T), dimension(rank) :: dim, dimForChunk
    integer(HID_T) :: chunkProperties
    character(len=100) :: label1, label2, label3, label4

    offset = (/ iterationNum-1, 0, 0, 0 /)

    select case (arrayType)
    case (ARRAY_ITERATION_SPECIES_THETA_ZETA)
       originalDspaceID = dspaceIDForIterationSpeciesThetaZeta
       dim = dimForIterationSpeciesThetaZeta
       dimForChunk = dimForIterationSpeciesThetaZetaChunk
       chunkProperties = pForIterationSpeciesThetaZeta

       label1 = "zeta"
       label2 = "theta"
       label3 = "species"
       label4 = "iteration"

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

       call h5dsset_label_f(dsetID, 1, trim(label1), HDF5Error)
       call h5dsset_label_f(dsetID, 2, trim(label2), HDF5Error)
       call h5dsset_label_f(dsetID, 3, trim(label3), HDF5Error)
       call h5dsset_label_f(dsetID, 4, trim(label4), HDF5Error)

       call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)
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

  ! -----------------------------------------------------------------------------------

  subroutine writeHDF5ExtensibleField6(iterationNum, arrayName, data, arrayType, description)

    use export_f, only: export_f_xi_option

    implicit none

    integer, parameter :: rank = 6
    integer, intent(in) :: iterationNum
    character(len=*) :: arrayName
    integer(HID_T) :: dsetID
    integer(HID_T) :: dspaceID, memspaceID, originalDspaceID
    integer :: temp, arrayType
    character(len=*) :: description
    PetscScalar, dimension(:,:,:,:,:) :: data
    integer(HSIZE_T) :: offset(rank)
    integer(HSIZE_T), dimension(rank) :: dim, dimForChunk
    integer(HID_T) :: chunkProperties
    character(len=100) :: label1, label2, label3, label4, label5, label6

    offset = (/ iterationNum-1, 0, 0, 0, 0, 0 /)

    select case (arrayType)
    case (ARRAY_EXPORT_F)
       originalDspaceID = dspaceIDForExport_f
       dim = dimForExport_f
       dimForChunk = dimForExport_fChunk
       chunkProperties = pForExport_f

       label1 = "export_f_x"
       if (export_f_xi_option==0) then
          label2 = "Legendre mode number in xi, i.e. n in P_n(xi). First index is P_0(xi)=1."
       else
          label2 = "export_f_xi"
       end if
       label3 = "export_f_zeta"
       label4 = "export_f_theta"
       label5 = "species"
       label6 = "iteration"

    case default
       print *,"This is writeHDF5ExtensibleField6"
       print *,"Error! Invalid arrayType:",arrayType
       stop
    end select

    if (iterationNum==1) then
       call h5dcreate_f(HDF5FileID, arrayName, H5T_NATIVE_DOUBLE, originalDspaceID, &
            dsetID, HDF5Error, chunkProperties)

       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, &
            data, dimForSpecies, HDF5Error)

       call h5dsset_label_f(dsetID, 1, trim(label1), HDF5Error)
       call h5dsset_label_f(dsetID, 2, trim(label2), HDF5Error)
       call h5dsset_label_f(dsetID, 3, trim(label3), HDF5Error)
       call h5dsset_label_f(dsetID, 4, trim(label4), HDF5Error)
       call h5dsset_label_f(dsetID, 5, trim(label5), HDF5Error)
       call h5dsset_label_f(dsetID, 6, trim(label6), HDF5Error)

       call h5ltset_attribute_string_f(HDF5FileID, arrayName, attribute_name, description, HDF5Error)
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

  end subroutine writeHDF5ExtensibleField6

end module writeHDF5Output

