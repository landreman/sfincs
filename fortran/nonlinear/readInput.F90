module readInput

  use globalVariables

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

NDNHatdrhos, NDTHatdrhos

    namelist / flowControl / saveMatlabOutput, MatlabOutputFilename, &
         outputFilename, solveSystem, RHSMode, &
         saveMatricesAndVectorsInBinary, binaryOutputFilename

    namelist / geometryParameters / GHat, IHat, iota, epsilon_t, epsilon_h, &
         helicity_l, helicity_n, B0OverBBar, geometryScheme, &
         epsilon_antisymm, helicity_antisymm_l, helicity_antisymm_n, &
         JGboozer_file, JGboozer_file_NonStelSym, min_Bmn_to_load, &
         psiAHat, inputRadialCoordinate, aHat, &
         psiHat_wish, psiN_wish, rHat_wish, rN_wish

    namelist / speciesParameters / mHats, Zs, nHats, THats, &
         dNHatdpsiHats, dTHatdpsiHats, &
         dNHatdpsiNs,   dTHatdpsiNs, &
         dNHatdrHats,   dTHatdrHats, &
         dNHatdrNs,     dTHatdrNs

    namelist / physicsParameters / Delta, alpha, nu_n, EParallelHat, &
         dPhiHatdpsiN, collisionOperator, constraintScheme, includeXDotTerm, &
         includeElectricFieldTermInXiDot, useDKESExBDrift, include_fDivVE_term, nonlinear

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
         useIterativeSolver, thetaDerivativeScheme, zetaDerivativeScheme, &
         whichParallelSolverToFactorPreconditioner, &
         PETSCPreallocationStrategy

    namelist / preconditionerOptions / preconditioner_x, preconditioner_x_min_L, preconditioner_zeta, &
         preconditioner_theta, preconditioner_xi, preconditioner_species

    Zs = speciesNotInitialized
    mHats = speciesNotInitialized
    nHats = speciesNotInitialized
    dNHatdpsiNs = speciesNotInitialized
    THats = speciesNotInitialized
    dTHatdpsiNs = speciesNotInitialized

    filename = trim(inputFilename)

    fileUnit=11
    open(unit=fileUnit, file=filename,    action="read", status="old", iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
       print *,"Proc ",myRank,": Error opening ", trim(filename)
       stop
    else
       read(fileUnit, nml=flowControl, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
               " but not read data from the flowControl namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from flowControl namelist in ", trim(filename), "."
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
    end if
    close(unit = fileUnit)

    ! Validate species parameters                                                                                                                                                                            

    NZs = maxNumSpecies
    NMHats = maxNumSpecies
    NNhats = maxNumSpecies
    NTHats = maxNumSpecies
    NDNhatdpsiNs = maxNumSpecies
    NDTHatdpsiNs = maxNumSpecies
    NDNhatdrhos = maxNumSpecies
    NDTHatdrhos = maxNumSpecies

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
       stop
    end if

    if (NZs /= NNHats) then
       print *,"Error: number of species charges (Zs) differs from the number of species densities (nHats)."
       stop
    end if

    if (NZs /= NTHats) then
       print *,"Error: number of species charges (Zs) differs from the number of species temperatures (THats)."
       stop
    end if

    select case (inputRadialCoordinate)
    case (0)
       ! Input radial coordinate is psiHat:

       if (NZs /= NDNHatdpsiHats) then
          print *,"Error: number of species charges (Zs) differs from the number of species density gradients for radial coordinate psiHat."
          stop
       end if

       if (NZs /= NDTHatdpsiHats) then
          print *,"Error: number of species charges (Zs) differs from the number of species temperature gradients for radial coordinate psiHat."
          stop
       end if

    case (1)
       ! Input radial coordinate is psiN:

       if (NZs /= NDNHatdpsiNs) then
          print *,"Error: number of species charges (Zs) differs from the number of species density gradients for radial coordinate psiN."
          stop
       end if

       if (NZs /= NDTHatdpsiNs) then
          print *,"Error: number of species charges (Zs) differs from the number of species temperature gradients for radial coordinate psiN."
          stop
       end if

    case (2)
       ! Input radial coordinate is rHat:

       if (NZs /= NDNHatdrHats) then
          print *,"Error: number of species charges (Zs) differs from the number of species density gradients for radial coordinate rHat."
          stop
       end if

       if (NZs /= NDTHatdrHats) then
          print *,"Error: number of species charges (Zs) differs from the number of species temperature gradients for radial coordinate rHat."
          stop
       end if

    case (3)
       ! Input radial coordinate is rN:

       if (NZs /= NDNHatdrNs) then
          print *,"Error: number of species charges (Zs) differs from the number of species density gradients for radial coordinate rN."
          stop
       end if

       if (NZs /= NDTHatdrNs) then
          print *,"Error: number of species charges (Zs) differs from the number of species temperature gradients for radial coordinate rN."
          stop
       end if


    case default
       print *,"Error! Invalid gradientInputScheme."
       stop
    end select

    Nspecies = NZs

    ! Using the selected radial coordinate, set input quantities for the other radial coordinates:
    call setInputRadialCoordinateWish()
    ! Note that this call only sets the "wish" radial coordinates, not the final radial coordinates
    ! or the input gradients. These quantities will be set later after we have loaded the magnetic
    ! geometry, in case the final radial coordinate is different from the "wish" values.

    ! Validate some other input quantities:

    if (preconditioner_theta < 0 .or. preconditioner_theta > 1) then
       print *,"Error! preconditioner_theta must be 0 or 1."
       stop
    end if
    if (preconditioner_zeta < 0 .or. preconditioner_zeta > 1) then
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

  end subroutine readNamelistInput

end module readInput

