module readInput

  use globalVariables
  use xGrid, only: xGrid_k

  implicit none

contains

  subroutine readNamelistInput()

    implicit none

    character(len=100) :: filename
    integer :: fileUnit, didFileAccessWork
    integer :: NMHats, NZs, NNHats, NTHats, NDNHatdpsiNs, NDTHatdpsiNs, i

    namelist / flowControl / programMode, saveMatlabOutput, outputScheme, MatlabOutputFilename, &
         outputFilename, parallelizeOverScan, solveSystem, RHSMode, &
         saveMatricesAndVectorsInBinary, binaryOutputFilename

    namelist / geometryParameters / GHat, IHat, iota, epsilon_t, epsilon_h, &
         helicity_l, helicity_n, B0OverBBar, geometryScheme, &
         epsilon_antisymm, helicity_antisymm_l, helicity_antisymm_n, &
         fort996boozer_file, JGboozer_file, JGboozer_file_NonStelSym, normradius_wish, min_Bmn_to_load

    namelist / speciesParameters / mHats, Zs, nHats, THats, dNHatdpsiNs, dTHatdpsiNs

    namelist / physicsParameters / speciesMode, Delta, alpha, psiAHat, nu_n, EParallelHat, &
         dPhiHatdpsiN, collisionOperator, constraintScheme, includeXDotTerm, &
         includeElectricFieldTermInXiDot, useDKESExBDrift, include_fDivVE_term, &
         dPhiHatdpsiN_min, dPhiHatdpsiN_max, NErs

    namelist / resolutionParameters / forceOddNthetaAndNzeta, &
         Ntheta, NthetaMaxFactor, NthetaMinFactor, NthetaNumRuns, &
         Nzeta, NzetaMaxFactor, NzetaMinFactor, NzetaNumRuns, &
         Nxi, NxiMaxFactor, NxiMinFactor, NxiNumRuns, &
         NL, NLMaxFactor, NLMinFactor, NLNumRuns, &
         Nx, NxMaxFactor, NxMinFactor, NxNumRuns, &
         xMax, xMaxMaxFactor, xMaxMinFactor, xMaxNumRuns, &
         solverTolerance, solverToleranceMinFactor, solverToleranceMaxFactor, solverToleranceNumRuns, &
         NxPotentialsPerVth, NxPotentialsPerVthMaxFactor, NxPotentialsPerVthMinFactor, NxPotentialsPerVthNumRuns

    namelist / otherNumericalParameters /  &
         useIterativeSolver, thetaDerivativeScheme, whichParallelSolverToFactorPreconditioner, &
         PETSCPreallocationStrategy, xGridScheme, xPotentialsGridScheme, xGrid_k

    namelist / preconditionerOptions / preconditioner_x, preconditioner_x_min_L, preconditioner_zeta, &
         preconditioner_theta, preconditioner_xi, preconditioner_species, &
         preconditioner_theta_min_L, preconditioner_zeta_min_L

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

    if (NZs /= NDNHatdpsiNs) then
       print *,"Error: number of species charges (Zs) differs from the number of species density gradients (dNHatdpsiNs)."
       stop
    end if

    if (NZs /= NDTHatdpsiNs) then
       print *,"Error: number of species charges (Zs) differs from the number of species temperature gradients (dTHatdpsiNs)."
       stop
    end if

    Nspecies = NZs

  end subroutine readNamelistInput

end module readInput

