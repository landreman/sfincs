module classicalTransport

#include "PETScVersions.F90"

  implicit none

contains
 
  subroutine calculateClassicalParticleFlux(classicalPF, classicalHF)
    use globalVariables, only: Nspecies, Ntheta, Nzeta, Delta, nu_n, mHats, THats, Zs, thetaWeights, zetaWeights, DHat, gpsipsi, BHat,alpha, Phi1Hat, VPrimeHat, NHats, dNHatdpsiHats, dTHatdpsiHats
    implicit None
    integer :: itheta, izeta, iSpeciesA, iSpeciesB
    PetscScalar :: geometry1, geometry2, xab2, Mab00, Mab01, Mab11, Nab11 !, y, classicalPF_OLD
    PetscScalar, dimension(Nspecies), intent(out) :: classicalPF, classicalHF

    do ispeciesA=1,Nspecies
       
       classicalPF(ispeciesA) = 0
       classicalHF(ispeciesA) = 0
       ! classicalPF_OLD = 0

       do ispeciesB=1,Nspecies
          
          xab2 = mHats(ispeciesA)*THats(ispeciesB)/(mHats(ispeciesB)*THats(ispeciesA))

          ! Braginskii matrix elements
          Mab00 = - (1+mHats(ispeciesA)/mHats(ispeciesB)) * (1 + xab2)
          Mab01 = - 1.5 * (1+mHats(ispeciesA)/mHats(ispeciesB))
          Mab11 = -(13. + 16. * xab2 + 30. * xab2**2)/4.
          Nab11 = 27. * mHats(ispeciesA)/(4.*mHats(ispeciesA))
          Mab00 = Mab00/((1+xab2)**(2.5))
          Mab01 = Mab01/((1+xab2)**(2.5))
          Mab11 = Mab11/((1+xab2)**(2.5))
          Nab11 = Nab11/((1+xab2)**(2.5))
          
          ! geometry1 = <|\nabla \psi|^2 n_a n_b/(B^2 N_a N_b)>
          geometry1 = 0
          ! geometry2 = <|\nabla \psi|^2 n_a n_b Phi1/(B^2 N_a N_b)>
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
          
          geometry1 = geometry1/VPrimeHat    
          geometry2 = geometry2/VPrimeHat

          classicalPF(ispeciesA) = classicalPF(ispeciesA) &
               + Zs(ispeciesB)**2 * NHats(ispeciesB) *( &
               geometry1 * Mab00 * (THats(ispeciesA) * dNHatdpsiHats(ispeciesA)/(NHats(ispeciesA) * Zs(ispeciesA)) - THats(ispeciesB) * dNHatdpsiHats(ispeciesB)/(NHats(ispeciesB) * Zs(ispeciesB))) &
               + geometry2 * alpha * Mab00 * (dTHatdpsiHats(ispeciesA)/THats(ispeciesA) - dTHatdpsiHats(ispeciesB)/THats(ispeciesB)) &
               + geometry1 * ((Mab00 - Mab01) * dTHatdpsiHats(ispeciesA)/Zs(ispeciesA) - (Mab00 - xab2*Mab01) * dTHatdpsiHats(ispeciesB)/Zs(ispeciesB)) &
          )

          classicalHF(ispeciesA) = classicalHF(ispeciesA) &
               + Zs(ispeciesB)**2 * NHats(ispeciesB) *( &
               geometry1 * Mab01 * (THats(ispeciesA) * dNHatdpsiHats(ispeciesA)/(NHats(ispeciesA) * Zs(ispeciesA)) - THats(ispeciesB) * dNHatdpsiHats(ispeciesB)/(NHats(ispeciesB) * Zs(ispeciesB))) &
               + geometry2 * alpha * Mab01 * (dTHatdpsiHats(ispeciesA)/THats(ispeciesA) - dTHatdpsiHats(ispeciesB)/THats(ispeciesB)) &
               + geometry1 * ((Mab01 - Mab11) * dTHatdpsiHats(ispeciesA)/Zs(ispeciesA) - (Mab01 + Nab11) * dTHatdpsiHats(ispeciesB)/Zs(ispeciesB)) &
               )

          ! y = THats(ispeciesA) * mHats(ispeciesB)/(THats(ispeciesB) * mHats(ispeciesA))
          ! classicalPF_OLD = classicalPF_OLD &
       !         + (geometry1 * ( &
       !         + Zs(ispeciesA) * (dNHatdpsiHats(ispeciesB)/NHats(ispeciesB) - 0.5 * dTHatdpsiHats(ispeciesB)/THats(ispeciesB)) &
       !         + Zs(ispeciesA) * 1.5 * (THats(ispeciesA) * mHats(ispeciesB)/(THats(ispeciesB) * mHats(ispeciesA) + THats(ispeciesA) * mHats(ispeciesB))) * dTHatdpsiHats(ispeciesB)/THats(ispeciesB) &
       !         - Zs(ispeciesB) * (THats(ispeciesA)/THats(ispeciesB)) * (dNHatdpsiHats(ispeciesA)/NHats(ispeciesA) - 0.5 * dTHatdpsiHats(ispeciesA)/THats(ispeciesA)) &
       !         -  Zs(ispeciesB) * 1.5 * THats(ispeciesA) * (mHats(ispeciesA)/(THats(ispeciesB) * mHats(ispeciesA) + THats(ispeciesA) * mHats(ispeciesB))) * dTHatdpsiHats(ispeciesA)/THats(ispeciesA)) &
       !         + geometry2 * Zs(ispeciesA) * Zs(ispeciesB) * alpha/THats(ispeciesB) * (dTHatdpsiHats(ispeciesB)/THats(ispeciesB) - dTHatdpsiHats(ispeciesA)/THats(ispeciesA))) &
       !         * NHats(ispeciesB) *Zs(ispeciesB) *sqrt(mHats(ispeciesB)/THats(ispeciesB)) *(1 + mHats(ispeciesB)/mHats(ispeciesA)) &
       !         /(4*(1 + y)**(1.5))          
       end do       
       classicalPF(ispeciesA) = Zs(ispeciesA) * NHats(ispeciesA) * Delta**2 * nu_n * sqrt(mHats(ispeciesA)) * classicalPF(ispeciesA)/(2*THats(ispeciesA)**1.5) 
       classicalHF(ispeciesA) = -Zs(ispeciesA) * NHats(ispeciesA) * Delta**2 * nu_n * sqrt(mHats(ispeciesA)) * classicalHF(ispeciesA)/(2*sqrt(THats(ispeciesA)))
       ! The total heat flux: Q_a = q_a + 2.5 T_a Gamma_a
       classicalHF(ispeciesA) = classicalHF(ispeciesA) + 2.5 * Thats(ispeciesA) * classicalPF(ispeciesA)
       !classicalPF_OLD = NHats(ispeciesA) * 2 * Delta**2 * nu_n * classicalPF_OLD


       
    end do

  end subroutine calculateClassicalParticleFlux
      
end module classicalTransport
