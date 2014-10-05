
    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] Total elapsed time: ", elapsedTime, " seconds."
       do ispecies=1,Nspecies
          if (Nspecies>1) then
             print *,"[",myCommunicatorIndex,"] Results for species ",ispecies,":"
          end if
          print *,"[",myCommunicatorIndex,"]   FSADensityPerturbation:  ", FSADensityPerturbation(ispecies)
          print *,"[",myCommunicatorIndex,"]   FSABFlow:                ", FSABFlow(ispecies)
          print *,"[",myCommunicatorIndex,"]   FSAPressurePerturbation: ", FSAPressurePerturbation(ispecies)
          print *,"[",myCommunicatorIndex,"]   NTV:                     ", NTV(ispecies)
          print *,"[",myCommunicatorIndex,"]   particleFlux:            ", particleflux(ispecies)
          print *,"[",myCommunicatorIndex,"]   momentumFlux:            ", momentumflux(ispecies)
          print *,"[",myCommunicatorIndex,"]   heatFlux:                ", heatflux(ispecies)
       end do
       print *,"[",myCommunicatorIndex,"] FSABjHat (bootstrap current): ", FSABjHat
       if (rhsMode == 2) then
          print *,"[",myCommunicatorIndex,"] Transport matrix:"
          do i=1,3
             print *,"[",myCommunicatorIndex,"]   ", transportMatrix(i,:)
          end do
       end if

    end if


