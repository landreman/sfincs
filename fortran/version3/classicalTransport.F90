module classicalTransport

#include "PETScVersions.F90"

  implicit none

contains
 
  subroutine calculateClassicalParticleFlux(classicalPF)
    use globalVariables, only: Nspecies, Ntheta, Nzeta, Delta, nu_n, mHats, THats, Zs, thetaWeights, zetaWeights, DHat, gpsipsi, BHat,alpha, Phi1Hat, VPrimeHat, NHats, dNHatdpsiHats, dTHatdpsiHats
    implicit None
    integer :: itheta, izeta, iSpeciesA, iSpeciesB
    PetscScalar :: geometry1, geometry2, y
    PetscScalar, dimension(Nspecies), intent(out) :: classicalPF

    do ispeciesA=1,Nspecies
       classicalPF(ispeciesA) = 0

       do ispeciesB=1,Nspecies
          if (ispeciesA == ispeciesB) then
             ! term is exactly zero
             cycle
          end if
          
          
          y = THats(ispeciesA) * mHats(ispeciesB)/(THats(ispeciesB) * mHats(ispeciesA))
          
          ! geometry1 = <|\nabla \psi|^2 n_a n_b/B^2>
          geometry1 = 0
          ! geometry2 = <|\nabla \psi|^2 n_a n_b Phi1/B^2>
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
          
          ! include other prefactors in geometry
          !geometry1 = geometry1 *Zs(ispeciesB) *sqrt(mHats(ispeciesB)/THats(ispeciesB)) *(1 + mHats(ispeciesB)/mHats(ispeciesA))/(4*(1 + y)**(1.5))
          !geometry2 = geometry2 
          
          classicalPF(ispeciesA) = classicalPF(ispeciesA) &
               + (geometry1 * ( &
               + Zs(ispeciesA) * (dNHatdpsiHats(ispeciesB)/NHats(ispeciesB) - 0.5 * dTHatdpsiHats(ispeciesB)/THats(ispeciesB)) &
               + Zs(ispeciesA) * 1.5 * (THats(ispeciesA) * mHats(ispeciesB)/(THats(ispeciesB) * mHats(ispeciesA) + THats(ispeciesA) * mHats(ispeciesB))) * dTHatdpsiHats(ispeciesB)/THats(ispeciesB) &
               - Zs(ispeciesB) * (THats(ispeciesA)/THats(ispeciesB)) * (dNHatdpsiHats(ispeciesA)/NHats(ispeciesA) - 0.5 * dTHatdpsiHats(ispeciesA)/THats(ispeciesA)) &
               -  Zs(ispeciesB) * 1.5 * THats(ispeciesA) * (mHats(ispeciesA)/(THats(ispeciesB) * mHats(ispeciesA) + THats(ispeciesA) * mHats(ispeciesB))) * dTHatdpsiHats(ispeciesA)/THats(ispeciesA)) &
               + geometry2 * Zs(ispeciesA) * Zs(ispeciesB) * alpha/THats(ispeciesB) * (dTHatdpsiHats(ispeciesB)/THats(ispeciesB) - dTHatdpsiHats(ispeciesA)/THats(ispeciesA))) &
               * NHats(ispeciesB) *Zs(ispeciesB) *sqrt(mHats(ispeciesB)/THats(ispeciesB)) *(1 + mHats(ispeciesB)/mHats(ispeciesA)) &
               /(4*(1 + y)**(1.5))          
       end do
       classicalPF(ispeciesA) = NHats(ispeciesA) * 2 * Delta**2 * nu_n * classicalPF(ispeciesA)
    end do

  end subroutine calculateClassicalParticleFlux
      
end module classicalTransport
