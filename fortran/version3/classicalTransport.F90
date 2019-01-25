module classicalTransport

#include "PETScVersions.F90"

  implicit none

contains
 
  subroutine calculateClassicalParticleFlux(classicalPF, classicalHF)
    use globalVariables, only: Nspecies, Ntheta, Nzeta, Delta, nu_n, mHats, THats, Zs, thetaWeights, zetaWeights, DHat, gpsipsi, BHat,alpha, Phi1Hat, VPrimeHat, NHats, dNHatdpsiHats, dTHatdpsiHats
    implicit None
    integer :: itheta, izeta, iSpeciesA, iSpeciesB
    PetscScalar :: geometry1, geometry2, xab2, Mab00, Mab01, Mab11, Nab11
    PetscScalar, dimension(Nspecies), intent(out) :: classicalPF, classicalHF

    do ispeciesA=1,Nspecies
       
       classicalPF(ispeciesA) = 0
       classicalHF(ispeciesA) = 0

       do ispeciesB=1,Nspecies
          
          xab2 = mHats(ispeciesA)*THats(ispeciesB)/(mHats(ispeciesB)*THats(ispeciesA))

          ! Braginskii matrix elements
          Mab00 = - (1+mHats(ispeciesA)/mHats(ispeciesB)) * (1 + xab2)
          Mab01 = - 1.5 * (1+mHats(ispeciesA)/mHats(ispeciesB))
          Mab11 = -(13. + 16. * xab2 + 30. * xab2**2)/4.
          Nab11 = 27. * mHats(ispeciesA)/(4.*mHats(ispeciesB))
          Mab00 = Mab00/((1+xab2)**(2.5))
          Mab01 = Mab01/((1+xab2)**(2.5))
          Mab11 = Mab11/((1+xab2)**(2.5))
          Nab11 = Nab11/((1+xab2)**(2.5))

          ! Benchmarked against python impelementation 2019-01
          !print *,Mab00,Mab01,Mab11,Nab11
          
          ! geometry1 = <|\nabla \psi|^2 n_a n_b/(B^2)>
          geometry1 = 0
          ! geometry2 = <|\nabla \psi|^2 n_a n_b Phi1/(B^2)>
          geometry2 = 0
          do itheta=1,Ntheta
             do izeta=1,Nzeta
                geometry1 = geometry1 &
                     + thetaWeights(itheta) * zetaWeights(izeta)/DHat(itheta,izeta) &
                     * (gpsipsi(itheta,izeta)/BHat(itheta,izeta)**2) &
                     * exp(-alpha *( Zs(iSpeciesA)/THats(iSpeciesA) + Zs(iSpeciesB)/THats(iSpeciesB)) * Phi1Hat(itheta,izeta))
                geometry2 = geometry2 &
                     + thetaWeights(itheta) * zetaWeights(izeta)/DHat(itheta,izeta) &
                     * (gpsipsi(itheta,izeta)/BHat(itheta,izeta)**2) * Phi1Hat(itheta,izeta) &
                     * exp(-alpha *( Zs(iSpeciesA)/THats(iSpeciesA) + Zs(iSpeciesB)/THats(iSpeciesB)) * Phi1Hat(itheta,izeta))
             end do
          end do   
          geometry1 = NHats(ispeciesA) * NHats(ispeciesB) * geometry1/VPrimeHat    
          geometry2 = NHats(ispeciesA) * NHats(ispeciesB) * geometry2/VPrimeHat
          ! Benchmarked against python impelementation 2019-01
          ! There is disagreement of about order 10^{-3}, which may be due to have
          ! .bc file is read in the external vs internal implementation.
          ! print *, geometry1, geometry2

          classicalPF(ispeciesA) = classicalPF(ispeciesA) &
               + Zs(ispeciesB)**2 * ( &
               geometry1 * Mab00 * (THats(ispeciesA) * dNHatdpsiHats(ispeciesA)/(NHats(ispeciesA) * Zs(ispeciesA)) - THats(ispeciesB) * dNHatdpsiHats(ispeciesB)/(NHats(ispeciesB) * Zs(ispeciesB))) &
               + geometry2 * alpha * Mab00 * (dTHatdpsiHats(ispeciesA)/THats(ispeciesA) - dTHatdpsiHats(ispeciesB)/THats(ispeciesB)) &
               + geometry1 * ((Mab00 - Mab01) * dTHatdpsiHats(ispeciesA)/Zs(ispeciesA) - (Mab00 - xab2*Mab01) * dTHatdpsiHats(ispeciesB)/Zs(ispeciesB)) &
          )

          classicalHF(ispeciesA) = classicalHF(ispeciesA) &
               + Zs(ispeciesB)**2  *( &
               geometry1 * Mab01 * (THats(ispeciesA) * dNHatdpsiHats(ispeciesA)/(NHats(ispeciesA) * Zs(ispeciesA)) - THats(ispeciesB) * dNHatdpsiHats(ispeciesB)/(NHats(ispeciesB) * Zs(ispeciesB))) &
               + geometry2 * alpha * Mab01 * (dTHatdpsiHats(ispeciesA)/THats(ispeciesA) - dTHatdpsiHats(ispeciesB)/THats(ispeciesB)) &
               + geometry1 * ((Mab01 - Mab11) * dTHatdpsiHats(ispeciesA)/Zs(ispeciesA) - (Mab01 + Nab11) * dTHatdpsiHats(ispeciesB)/Zs(ispeciesB)) &
               )

          ! Benchmarked against python impelementation 2019-01
          ! print *,classicalPF(ispeciesA),classicalHF(ispeciesA)
        
       end do       
       classicalPF(ispeciesA) = Zs(ispeciesA) * Delta**2 * nu_n * sqrt(mHats(ispeciesA)) * classicalPF(ispeciesA)/(2*THats(ispeciesA)**1.5) 
       classicalHF(ispeciesA) = -Zs(ispeciesA)  * Delta**2 * nu_n * sqrt(mHats(ispeciesA)) * classicalHF(ispeciesA)/(2*sqrt(THats(ispeciesA)))
       ! The total heat flux: Q_a = q_a + 2.5 T_a Gamma_a
       classicalHF(ispeciesA) = classicalHF(ispeciesA) + 2.5 * THats(ispeciesA) * classicalPF(ispeciesA)


       
    end do

  end subroutine calculateClassicalParticleFlux
      
end module classicalTransport
