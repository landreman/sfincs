module readInput

  use globalVariables

  implicit none

contains

  subroutine readNamelistInput(filename)

    implicit none

    character(len=*), intent(in) :: filename
    integer :: fileUnit, didFileAccessWork

    namelist / flowControl / programMode, saveMatlabOutput, outputScheme, MatlabOutputFilename, &
         outputFilename, parallelizeOverScan, solveSystem, RHSMode, &
         saveMatricesAndVectorsInBinary, binaryOutputFilename

    namelist / geometryParameters / GHat, IHat, iota, epsilon_t, epsilon_h, &
         helicity_l, helicity_n, B0OverBBar, geometryScheme, &
         epsilon_antisymm, helicity_antisymm_l, helicity_antisymm_n, &
         fort996boozer_file, JGboozer_file, JGboozer_file_NonStelSym, normradius_wish, min_Bmn_to_load

    namelist / physicsParameters / speciesMode, Delta, omega, psiAHat, nuN, nuPrime, THat, nHat, EHat, &
         dnHatdpsi, dTHatdpsi, dPhiHatdpsi, EStar, collisionOperator, constraintScheme, includeXDotTerm, &
         includeElectricFieldTermInXiDot, useDKESExBDrift, include_fDivVE_term, EStarMin, EStarMax, NEStar

    namelist / resolutionParameters / forceOddNthetaAndNzeta, &
         Ntheta, NthetaMaxFactor, NthetaMinFactor, NthetaNumRuns, &
         Nzeta, NzetaMaxFactor, NzetaMinFactor, NzetaNumRuns, &
         Nxi, NxiMaxFactor, NxiMinFactor, NxiNumRuns, &
         NL, NLMaxFactor, NLMinFactor, NLNumRuns, &
         Nx, NxMaxFactor, NxMinFactor, NxNumRuns, &
         xMax, xMaxMaxFactor, xMaxMinFactor, xMaxNumRuns, &
         solverTolerance, solverToleranceMinFactor, solverToleranceMaxFactor, solverToleranceNumRuns, &
         NxPotentialsPerVth, NxPotentialsPerVthMaxFactor, NxPotentialsPerVthMinFactor, NxPotentialsPerVthNumRuns

    namelist / otherNumericalParameters / &
         useIterativeSolver, thetaDerivativeScheme, whichParallelSolverToFactorPreconditioner, &
         PETSCPreallocationStrategy

    namelist / preconditionerOptions / preconditioner_x, preconditioner_x_min_L, preconditioner_zeta, &
         preconditioner_theta, preconditioner_xi

    fileUnit=11
    open(unit=fileUnit, file=filename,    action="read", status="old", iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
       print *,"Proc ",myRank,": Error opening ", filename
       stop
    else
       read(fileUnit, nml=flowControl, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the flowControl namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from flowControl namelist in ", filename, "."
       end if

       read(fileUnit, nml=geometryParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the geometryParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from geometryParameters namelist in ", filename, "."
       end if

       read(fileUnit, nml=physicsParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the physicsParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from physicsParameters namelist in ", filename, "."
       end if

       read(fileUnit, nml=resolutionParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the resolutionParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from resolutionParameters namelist in ", filename, "."
       end if

       read(fileUnit, nml=otherNumericalParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the otherNumericalParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from otherNumericalParameters namelist in ", filename, "."
       end if

       read(fileUnit, nml=preconditionerOptions, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the preconditionerOptions namelist in it."
          print *,"Make sure there is a carriage return at the end of the file."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from preconditionerOptions namelist in ", filename, "."
       end if
    end if
    close(unit = fileUnit)
  end subroutine readNamelistInput

end module readInput

