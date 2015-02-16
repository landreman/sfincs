module readInput

  use globalVariables
  use export_f

  implicit none

contains

  subroutine readNamelistInput()

    implicit none

    character(len=100) :: filename
    integer :: fileUnit, didFileAccessWork, i
    integer :: NMHats, NZs, NNHats, NTHats
    integer :: NDNHatdpsiHats, NDTHatdpsiHats
    integer :: NDNHatdpsiNs,   NDTHatdpsiNs
    integer :: NDNHatdrHats,   NDTHatdrHats
    integer :: NDNHatdrNs,     NDTHatdrNs

    namelist / general / saveMatlabOutput, MatlabOutputFilename, &
         outputFilename, solveSystem, RHSMode, &
         saveMatricesAndVectorsInBinary, binaryOutputFilename

    namelist / geometryParameters / GHat, IHat, iota, epsilon_t, epsilon_h, &
         helicity_l, helicity_n, B0OverBBar, geometryScheme, &
         epsilon_antisymm, helicity_antisymm_l, helicity_antisymm_n, &
         equilibriumFile, min_Bmn_to_load, &
         psiAHat, inputRadialCoordinate, inputRadialCoordinateForGradients, aHat, &
         psiHat_wish, psiN_wish, rHat_wish, rN_wish, &
         VMECRadialOption

    namelist / speciesParameters / mHats, Zs, nHats, THats, &
         dNHatdpsiHats, dTHatdpsiHats, &
         dNHatdpsiNs,   dTHatdpsiNs, &
         dNHatdrHats,   dTHatdrHats, &
         dNHatdrNs,     dTHatdrNs

    namelist / physicsParameters / Delta, alpha, nu_n, EParallelHat, &
         collisionOperator, constraintScheme, includeXDotTerm, &
         includeElectricFieldTermInXiDot, useDKESExBDrift, include_fDivVE_term, nonlinear, &
         dPhiHatdpsiHat, dPhiHatdpsiN, dPhiHatdrHat, dPhiHatdrN, &
         includeTemperatureEquilibrationTerm, includePhi1, &
         nuPrime, EStar

    namelist / resolutionParameters / forceOddNthetaAndNzeta, &
         Ntheta, &
         Nzeta, &
         Nxi, &
         NL, &
         Nx, &
         xMax, &
         solverTolerance, &
         NxPotentialsPerVth

    namelist / otherNumericalParameters /  &
         useIterativeLinearSolver, thetaDerivativeScheme, zetaDerivativeScheme, &
         whichParallelSolverToFactorPreconditioner, &
         PETSCPreallocationStrategy

    namelist / preconditionerOptions / preconditioner_x, preconditioner_x_min_L, preconditioner_zeta, &
         preconditioner_theta, preconditioner_xi, preconditioner_species, reusePreconditioner, &
         preconditioner_theta_min_L, preconditioner_zeta_min_L

    namelist / export_f / export_full_f, export_delta_f, export_f_theta, export_f_zeta, export_f_x, export_f_xi, &
         export_f_theta_option, export_f_zeta_option, export_f_xi_option, export_f_x_option

    Zs = speciesNotInitialized
    mHats = speciesNotInitialized
    nHats = speciesNotInitialized
    THats = speciesNotInitialized
    dNHatdpsiHats = speciesNotInitialized
    dTHatdpsiHats = speciesNotInitialized
    dNHatdpsiNs = speciesNotInitialized
    dTHatdpsiNs = speciesNotInitialized
    dNHatdrHats = speciesNotInitialized
    dTHatdrHats = speciesNotInitialized
    dNHatdrNs = speciesNotInitialized
    dTHatdrNs = speciesNotInitialized

    export_f_theta = speciesNotInitialized
    export_f_zeta = speciesNotInitialized
    export_f_x = speciesNotInitialized
    export_f_xi = speciesNotInitialized

    filename = trim(inputFilename)

    fileUnit=11
    open(unit=fileUnit, file=filename,    action="read", status="old", iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
       print *,"Proc ",myRank,": Error opening ", trim(filename)
       stop
    else
       read(fileUnit, nml=general, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
               " but not read data from the general namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from general namelist in ", trim(filename), "."
       end if

       read(fileUnit, nml=geometryParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
               " but not read data from the geometryParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from geometryParameters namelist in ", trim(filename), "."
       end if

       read(fileUnit, nml=speciesParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
               " but not read data from the speciesParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from speciesParameters namelist in ", trim(filename), "."
       end if

       read(fileUnit, nml=physicsParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
               " but not read data from the physicsParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from physicsParameters namelist in ", trim(filename), "."
       end if

       read(fileUnit, nml=resolutionParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
               " but not read data from the resolutionParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from resolutionParameters namelist in ", trim(filename), "."
       end if

       read(fileUnit, nml=otherNumericalParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
               " but not read data from the otherNumericalParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from otherNumericalParameters namelist in ", trim(filename), "."
       end if

       read(fileUnit, nml=preconditionerOptions, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
               " but not read data from the preconditionerOptions namelist in it."
          print *,"Make sure there is a carriage return at the end of the file."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from preconditionerOptions namelist in ", trim(filename), "."
       end if

       read(fileUnit, nml=export_f, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
               " but not read data from the export_f namelist in it."
          print *,"Make sure there is a carriage return at the end of the file."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from export_f namelist in ", trim(filename), "."
       end if
    end if

    close(unit = fileUnit)

    ! Validate species parameters                                                                                                                                                                            

    NZs = maxNumSpecies
    NMHats = maxNumSpecies
    NNhats = maxNumSpecies
    NTHats = maxNumSpecies
    NDNhatdpsiHats = maxNumSpecies
    NDTHatdpsiHats = maxNumSpecies
    NDNhatdpsiNs = maxNumSpecies
    NDTHatdpsiNs = maxNumSpecies
    NDNhatdrHats = maxNumSpecies
    NDTHatdrHats = maxNumSpecies
    NDNhatdrNs = maxNumSpecies
    NDTHatdrNs = maxNumSpecies

    do i=1,maxNumSpecies
       if (Zs(i) == speciesNotInitialized) then
          NZs = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
       if (MHats(i) == speciesNotInitialized) then
          NMHats = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
       if (nHats(i) == speciesNotInitialized) then
          NNHats = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
       if (THats(i) == speciesNotInitialized) then
          NTHats = i-1
          exit
       end if
    end do

    ! --------------

    do i=1,maxNumSpecies
       if (dNHatdpsiHats(i) == speciesNotInitialized) then
          NDNHatdpsiHats = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
       if (dTHatdpsiHats(i) == speciesNotInitialized) then
          NDTHatdpsiHats = i-1
          exit
       end if
    end do

    ! --------------

    do i=1,maxNumSpecies
       if (dNHatdpsiNs(i) == speciesNotInitialized) then
          NDNHatdpsiNs = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
       if (dTHatdpsiNs(i) == speciesNotInitialized) then
          NDTHatdpsiNs = i-1
          exit
       end if
    end do

    ! --------------

    do i=1,maxNumSpecies
       if (dNHatdrHats(i) == speciesNotInitialized) then
          NDNHatdrHats = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
       if (dTHatdrHats(i) == speciesNotInitialized) then
          NDTHatdrHats = i-1
          exit
       end if
    end do

    ! --------------

    do i=1,maxNumSpecies
       if (dNHatdrNs(i) == speciesNotInitialized) then
          NDNHatdrNs = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
       if (dTHatdrNs(i) == speciesNotInitialized) then
          NDTHatdrNs = i-1
          exit
       end if
    end do

    ! --------------

    if (NZs /= NMHats) then
       print *,"Error: number of species charges (Zs) differs from the number of species masses (mHats)."
       print *,"NZs: ",NZs
       print *,"NMHats: ",NMHats
       stop
    end if

    if (NZs /= NNHats) then
       print *,"Error: number of species charges (Zs) differs from the number of species densities (nHats)."
       print *,"NZs: ",NZs
       print *,"NNHats: ",NNHats
       stop
    end if

    if (NZs /= NTHats) then
       print *,"Error: number of species charges (Zs) differs from the number of species temperatures (THats)."
       print *,"NZs: ",NZs
       print *,"NTHats: ",NTHats
       stop
    end if

    select case (inputRadialCoordinateForGradients)
    case (0)
       ! Input radial coordinate is psiHat:

       if (NZs /= NDNHatdpsiHats) then
          print *,"Error: number of species charges (Zs) differs from the number of species density gradients for radial coordinate psiHat."
          print *,"NZs: ",NZs
          print *,"NDNHatdpsiHats: ",NDNHatdpsiHats
          stop
       end if

       if (NZs /= NDTHatdpsiHats) then
          print *,"Error: number of species charges (Zs) differs from the number of species temperature gradients for radial coordinate psiHat."
          print *,"NZs: ",NZs
          print *,"NDTHatdpsiHats: ",NDTHatdpsiHats
          stop
       end if

    case (1)
       ! Input radial coordinate is psiN:

       if (NZs /= NDNHatdpsiNs) then
          print *,"Error: number of species charges (Zs) differs from the number of species density gradients for radial coordinate psiN."
          print *,"NZs: ",NZs
          print *,"NDNHatdpsiNs: ",NDNHatdpsiNs
          stop
       end if

       if (NZs /= NDTHatdpsiNs) then
          print *,"Error: number of species charges (Zs) differs from the number of species temperature gradients for radial coordinate psiN."
          print *,"NZs: ",NZs
          print *,"NDTHatdpsiNs: ",NDTHatdpsiNs
          stop
       end if

    case (2)
       ! Input radial coordinate is rHat:

       if (NZs /= NDNHatdrHats) then
          print *,"Error: number of species charges (Zs) differs from the number of species density gradients for radial coordinate rHat."
          print *,"NZs: ",NZs
          print *,"NDNHatdrHats: ",NDNHatdrHats
          stop
       end if

       if (NZs /= NDTHatdrHats) then
          print *,"Error: number of species charges (Zs) differs from the number of species temperature gradients for radial coordinate rHat."
          print *,"NZs: ",NZs
          print *,"NDTHatdrHats: ",NDTHatdrHats
          stop
       end if

    case (3)
       ! Input radial coordinate is rN:

       if (NZs /= NDNHatdrNs) then
          print *,"Error: number of species charges (Zs) differs from the number of species density gradients for radial coordinate rN."
          print *,"NZs: ",NZs
          print *,"NDNHatdrNs: ",NDNHatdrNs
          stop
       end if

       if (NZs /= NDTHatdrNs) then
          print *,"Error: number of species charges (Zs) differs from the number of species temperature gradients for radial coordinate rN."
          print *,"NZs: ",NZs
          print *,"NDTHatdrNs: ",NDTHatdrNs
          stop
       end if


    case default
       print *,"Error! Invalid inputRadialCoordinateForGradients."
       stop
    end select

    Nspecies = NZs


    ! Validate some other input quantities:

    if (preconditioner_theta < 0 .or. preconditioner_theta > 2) then
       print *,"Error! preconditioner_theta must be 0 or 1."
       stop
    end if
    if (preconditioner_zeta < 0 .or. preconditioner_zeta > 2) then
       print *,"Error! preconditioner_zeta must be 0 or 1."
       stop
    end if
    if (preconditioner_xi < 0 .or. preconditioner_xi > 1) then
       print *,"Error! preconditioner_xi must be 0 or 1."
       stop
    end if
    if (preconditioner_x < 0 .or. preconditioner_x > 4) then
       print *,"Error! preconditioner_x must be in the range [0, 4]."
       stop
    end if

    if (collisionOperator < 0 .or. collisionOperator > 2) then
       print *,"Error! collisionOperator must be 0, 1, or 2."
       stop
    end if

    if (constraintScheme < -1 .or. constraintScheme > 2) then
       print *,"Error! constraintScheme must be -1, 0, 1, or 2."
       stop
    end if

    ! ----------------------------------------------------------------

    ! Determine the number of points requested for each coordinate related to export_f:

    N_export_f_theta = max_N_export_f
    N_export_f_zeta  = max_N_export_f
    N_export_f_x     = max_N_export_f
    N_export_f_xi    = max_N_export_f

    do i=1,max_N_export_f
       if (export_f_theta(i) == speciesNotInitialized) then
          N_export_f_theta = i-1
          exit
       end if
    end do

    do i=1,max_N_export_f
       if (export_f_zeta(i) == speciesNotInitialized) then
          N_export_f_zeta = i-1
          exit
       end if
    end do

    do i=1,max_N_export_f
       if (export_f_xi(i) == speciesNotInitialized) then
          N_export_f_xi = i-1
          exit
       end if
    end do

    do i=1,max_N_export_f
       if (export_f_x(i) == speciesNotInitialized) then
          N_export_f_x = i-1
          exit
       end if
    end do


  end subroutine readNamelistInput

end module readInput

