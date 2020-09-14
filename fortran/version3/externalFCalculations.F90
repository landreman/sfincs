module externalFCalculations

#include "PETScVersions.F90"

  implicit none

  PetscScalar, dimension(:,:), allocatable :: externalFL
  PetscScalar, dimension(:,:), allocatable :: P
  PetscScalar, dimension(:,:), allocatable :: IL2, IL4, I1mL, I3mL

  
contains
 
  subroutine calculateExternalN()
    use globalVariables, only: pi, Ntheta, Nzeta, externalNspecies, externalNE, externalNxi, externalXi, externalE, externalF, externalMasses, externalN
    integer :: ispecies, itheta, izeta, ixi
    PetscScalar :: dE, dxi
    PetscScalar, dimension(:), allocatable :: velocityJacobian

    ! Global variable
    if (.not. allocated(externalN)) then
       allocate(externalN(externalNspecies, Ntheta, Nzeta))
    end if
    
    allocate(velocityJacobian(externalNE))
    
    

    dE = externalE(2) - externalE(1)
    dxi = externalXi(2) - externalXi(1)
    
    do ispecies=1,externalNspecies
       velocityJacobian = pi * sqrt(externalE) / &
            (externalMasses(ispecies) * sqrt(externalMasses(ispecies)))
       do itheta = 1,Ntheta
          do izeta = 1,Nzeta
             externalN(ispecies,itheta,izeta) = 0
             do ixi = 1,externalNxi
                externalN(ispecies,itheta,izeta) = externalN(ispecies,itheta,izeta) &
                     + sum(velocityJacobian * externalF(ispecies, itheta, izeta, ixi, :))
             end do
          end do
       end do
    end do
    externalN = externalN * dE * dxi

    deallocate(velocityJacobian)
    
  end subroutine calculateExternalN

  subroutine computeExternalRosenbluthPotentialResponse()
    use globalVariables, only: externalNspecies, externalNL, externalNE, externalNxi, externalXi, externalE, externalF, externalMasses, externalCharges, externalRosenPotentialTerms, nu_n, NHats, mHats, THats, Zs, x2, x, expx2, Ntheta, Nzeta, Nspecies, Nx, pi, zero
    integer :: ispeciesA, ispeciesB, itheta, izeta, L,  ix
    PetscScalar :: externalZ, externalMHat
    PetscScalar, dimension(:), allocatable :: xFactor
    PetscScalar, dimension(:,:), allocatable :: FL, d2GLdx2, HL, dHLdx

    ! global variable
    if (.not. allocated(externalRosenPotentialTerms)) then
       allocate(externalRosenPotentialTerms(Nspecies,Ntheta,Nzeta,externalNL,Nx))
    end if
    externalRosenPotentialTerms = zero 

    ! module globals
    if (.not. allocated(P)) then
       allocate(P(externalNxi,externalNL))
    end if
    if (.not. allocated(externalFL)) then
       allocate(externalFL(externalNL,externalNE))
    end if
    if (.not. allocated(IL2)) then
       allocate(IL2(externalNL,Nx))
    end if
    if (.not. allocated(IL4)) then
       allocate(IL4(externalNL,Nx))
    end if
    if (.not. allocated(I1mL)) then
       allocate(I1mL(externalNL,Nx))
    end if
    if (.not. allocated(I3mL)) then
       allocate(I3mL(externalNL,Nx))
    end if
    
    
    allocate(xFactor(Nx))
    allocate(FL(externalNL,Nx))
    allocate(d2GLdx2(externalNL,Nx))
    allocate(HL(externalNL,Nx))
    allocate(dHLdx(externalNL,Nx))

    ! calculate the  Legendre polynomial on the
    ! external xi grid
    call legendrePolynomials(P, externalXi, externalNxi, externalNL)
    ! P is verified to work here
    
    do ispeciesA=1,Nspecies
       xFactor = 3 * nu_n * nHats(ispeciesA) * Zs(ispeciesA)**2 * expx2/(2*pi*THats(ispeciesA) * sqrt(THats(ispeciesA) * mHats(ispeciesA)))
       do itheta=1,Ntheta
          do izeta=1,Nzeta
             do ispeciesB =1,externalNspecies
                externalZ = externalCharges(ispeciesB)
                externalMHat = externalMasses(ispeciesB)
                
                call calculateFL(FL,ispeciesA,ispeciesB,itheta,izeta) ! this must be called first
                call calculateGL(d2GLdx2,ispeciesA,ispeciesB,itheta,izeta) ! this must be called second
                call calculateHL(HL,dHLdx,ispeciesA,ispeciesB,itheta,izeta) ! this must be called third
                   
                do L=0,externalNL-1
                   do ix=1,Nx 
                      externalRosenPotentialTerms(ispeciesA,itheta,izeta,L+1,ix) &
                           = externalRosenPotentialTerms(ispeciesA,itheta,izeta,L+1,ix) &
                           + externalZ**2 * 2 * pi  * (mHats(ispeciesA)/externalMHat) * FL(L+1,ix) &
                           + externalZ**2 * x2(ix) * d2GLdx2(L+1,ix) &
                           - externalZ**2 * ((1 - mHats(ispeciesA)/externalMHat) * x(ix) * dHLdx(L+1,ix) + HL(L+1,ix))

                   end do
                   externalRosenPotentialTerms(ispeciesA,itheta,izeta,L+1,:) = externalRosenPotentialTerms(ispeciesA,itheta,izeta,L+1,:) * xFactor
                end do
             end do
          end do
       end do
    end do

    ! locals
    deallocate(xFactor)
    deallocate(FL)
    deallocate(d2GLdx2)
    deallocate(HL)
    deallocate(dHLdx)    
    
    ! module globals
    deallocate(P)
    deallocate(externalFL)
    deallocate(IL2)
    deallocate(IL4)
    deallocate(I1mL)
    deallocate(I3mL)
    
    
    ! These global variables are not needed any more
    if(allocated(externalF))deallocate(externalF)
    if(allocated(externalXi))deallocate(externalXi)
    if(allocated(externalE))deallocate(externalE)
    if(allocated(externalMasses))deallocate(externalMasses)
        
    
  end subroutine computeExternalRosenbluthPotentialResponse
  
  subroutine calculateFL(FL,ispeciesA,ispeciesB,itheta,izeta)
    use globalVariables, only: externalNE, externalNxi, externalNL, externalXi, externalE,  externalF, externalMasses, x2, mHats, THats, Nx, extrapolateExternalF
    integer, intent(in) :: ispeciesA,ispeciesB,itheta,izeta
    PetscScalar, dimension(:,:), intent(out) :: FL
    ! PetscScalar, dimension(:,:), allocatable :: FL2
    
    PetscScalar :: mHatB, mHatA
    PetscScalar, dimension(:), allocatable :: Es, interpolatedFLx
    PetscScalar, dimension(:,:), allocatable :: interpolatedFxix
    integer :: indexLower, indexUpper, iE, ix, ixi, L
    PetscScalar :: deltaELower, deltaEUpper, dxi
    
    allocate(Es(Nx))
    allocate(interpolatedFxix(externalNxi,Nx))
    allocate(interpolatedFLx(Nx))

    ! Used for comparing 2 methods of calc
    ! allocate(FL2(externalNL,Nx))
    ! project on externalF on Legendre polynomials
    dxi = externalXi(2) - externalXi(1)
    do L=0,externalNL-1
       externalFL(L+1,:) = 0.0
       do ixi = 1,externalNxi
          externalFL(L+1,:) = externalFL(L+1,:) + &
               P(ixi,L+1) * externalF(ispeciesB,itheta,izeta,ixi,:)
       end do
       externalFL(L+1,:) = externalFL(L+1,:) * (L+0.5) * dxi
    end do

    mHatB = externalMasses(ispeciesB)
    mHatA = mHats(ispeciesA)
    Es = THats(ispeciesA) * x2 * mHatB/mHatA
        
    ! Find the first index that is higher than Es
    do ix = 1,Nx
       indexUpper = 0
       do iE = 1,externalNE
          if ((externalE(iE)) >= Es(ix)) then
             indexUpper = iE
             exit
          end if
       end do

       if (indexUpper == 1) then
          ! Es(ix) is smaller than all externalE
          ! this will likely happen for ions colliding with
          ! fast external ions
          if (extrapolateExternalF) then
             indexLower = 1
             indexUpper = 2
             deltaELower = (Es(ix) - externalE(indexLower)) / (externalE(indexUpper) - externalE(indexLower))
             ! interpolatedFxix(:,ix) = externalF(ispeciesB,itheta,izeta,:,indexLower) + &
             !      (externalF(ispeciesB,itheta,izeta,:,indexUpper) - externalF(ispeciesB,itheta,izeta,:,indexLower)) * deltaELower

             do L=0,externalNL-1
                FL(L+1,ix) = externalFL(L+1,indexLower) + &
                     (externalFL(L+1,indexUpper) - externalFL(L+1,indexLower)) * deltaELower
             end do
          else
             ! assume that F is zero for values below
             ! its energy grid. This appears to be
             ! a suitable appromixation for my ASCOT data.
             do L=0,externalNL-1
                FL(L+1,ix) = 0
             end do
          end if
          
       else if (indexUpper == 0) then
          ! Es(ix) is larger than all externalE
          ! this will likely happen for electrons colliding with
          ! the external ions
          if (extrapolateExternalF) then
             ! extrapolate to get F
             indexLower = externalNE - 1
             indexUpper = externalNE
             deltaEUpper = (externalE(indexUpper) - Es(ix)) / (externalE(indexUpper) - externalE(indexLower))
             ! interpolatedFxix(:,ix) = externalF(ispeciesB,itheta,izeta,:,indexUpper) + &
             !      (externalF(ispeciesB,itheta,izeta,:,indexUpper) - externalF(ispeciesB,itheta,izeta,:,indexLower)) * deltaEUpper

             do L=0,externalNL-1
                FL(L+1,ix) = externalFL(L+1,indexUpper) + &
                     (externalFL(L+1,indexUpper) + externalFL(L+1,indexLower)) * deltaEUpper
             end do
          else
             ! assume that F is zero for values above
             ! its energy grid. This appears to be
             ! a suitable appromixation for my ASCOT data.   
             do L=0,externalNL-1
                FL(L+1,ix) = 0
             end do
          end if
          
          
       else
          ! interpolation
          ! SHOULD PROBABLY USE interpolationMatrix()
          ! Below is copy-pasted from
          ! 'readHDF5Input.F90 :: interpolateToThetaZetaGrid'
          indexLower = indexUpper - 1
          deltaEUpper = (externalE(indexUpper) - Es(ix)) / (externalE(indexUpper) - externalE(indexLower))
          deltaELower = (Es(ix)-externalE(indexLower)) / (externalE(indexUpper) - externalE(indexLower))

          ! interpolatedFxix(:,ix) = externalF(ispeciesB,itheta,izeta,:,indexLower) * deltaELower + externalF(ispeciesB,itheta,izeta,:,indexUpper) * deltaEUpper

          do L=0,externalNL-1
             FL(L+1,ix) = externalFL(L+1,indexUpper) * deltaELower + externalFL(L+1,indexLower) * deltaEUpper 
          end do
       end if
       
       ! do L=0,externalNL-1
       !    ! project on the Legendre polynomial
       !    FL2(L+1,ix) = 0.0
       !    do ixi = 1,externalNxi
       !       FL2(L+1,ix) = FL2(L+1,ix) + &
       !            P(ixi,L+1) * interpolatedFxix(ixi,ix)
       !    end do
       !    FL2(L+1,ix) = FL2(L+1,ix) * (L+0.5)  * dxi
       ! end do
    end do

    deallocate(Es)
    deallocate(interpolatedFxix)
    deallocate(interpolatedFLx)
    
    
  end subroutine calculateFL

  subroutine calculateGL(d2GLdx2,ispeciesA,ispeciesB,itheta,izeta)
    use globalVariables, only: externalNE, externalNL, externalE,  externalMasses, x, x2, mHats, THats, Nx, pi, extrapolateExternalF
    integer, intent(in) :: ispeciesA,ispeciesB,itheta,izeta
    PetscScalar, dimension(:,:), intent(out) :: d2GLdx2
    PetscScalar :: beta, mHatA, mHatB, THatA
    integer :: indexLower, indexUpper, iE, ix, L, extIntCase
    PetscScalar :: deltaELower, deltaEUpper
    PetscScalar :: IL2Lower, IL4Lower, I1mLLower, I3mLLower
    PetscScalar :: IL2Upper, IL4Upper, I1mLUpper, I3mLUpper
    PetscScalar :: dE, ELower, EUpper
    PetscScalar, dimension(:), allocatable :: Es

    allocate(Es(Nx))
    
    mHatB = externalMasses(ispeciesB)
    mHatA = mHats(ispeciesA)
    THatA = THats(ispeciesA)
    Es = THatA * x2 * mHatB/mHatA
    dE = externalE(2) - externalE(1)
    
    ! calculate the relevant integrals interpolated/extrapolated
    ! to the x grid
    ! First: Find the indices we need to calculate the integral for
    ! NOTE: this is different from finding the
    ! uppper indices in calculateFL
    ! since "integration grid" is offset by dE/2
    ! in our 1D cell-centered finite volume
    do ix=1,Nx
       indexUpper = 0
       do iE = 1,externalNE
          if (externalE(iE) + dE/2>= Es(ix)) then
             indexUpper = iE
             exit
          end if
       end do

       if (indexUpper == 1) then
          indexLower = 1
          indexUpper = 2
          extIntCase = 1
       else if (indexUpper == 0) then
          indexLower = externalNE - 1
          indexUpper = externalNE
          extIntCase = 2
       else
          ! interpolation
          indexLower = indexUpper - 1
          extIntCase = 3
       end if

       ELower = externalE(indexLower) + dE/2
       EUpper = externalE(indexUpper)+ dE/2
       
       
       deltaEUpper = (EUpper - Es(ix)) / (EUpper - ELower)
       deltaELower = (Es(ix) - ELower) / (EUpper - ELower)
       
       ! THERE WILL BE A LOT OF REDUNDANT CALCULATIONS HERE
       ! SINCE IL2 AND IL4 HAS A LOT OF OVERLAPPING INTEGRALS
       ! LIKEWISE FOR I1mL AND I3mL BELOW

       do L=0,externalNL-1          
          IL2Lower = 0.0
          IL4Lower = 0.0
          IL2Upper = 0.0
          IL4Upper = 0.0

          do iE = 1,indexLower
             IL2Lower = IL2Lower + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (L+2))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
             IL4Lower = IL4Lower + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (L+4))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
          end do

          do iE = 1,indexUpper
             IL2Upper = IL2Upper + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (L+2))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
             IL4Upper = IL4Upper + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (L+4))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
          end do
          
          IL4Lower = IL4Lower * dE
          IL2Lower = IL2Lower * dE
          IL4Upper = IL4Upper * dE
          IL2Upper = IL2Upper * dE
          
          I1mLLower = 0.0
          I3mLLower = 0.0
          I1mLUpper = 0.0
          I3mLUpper = 0.0
          
          do iE = indexLower+1,externalNE
             I1mLLower = I1mLLower + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (1-L))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
             I3mLLower = I3mLLower + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (3-L))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
          end do
                    
          do iE = indexUpper+1,externalNE
             I1mLUpper = I1mLUpper + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (1-L))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
             I3mLUpper = I3mLUpper + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (3-L))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
          end do
          I1mLLower = I1mLLower * dE
          I3mLLower = I3mLLower * dE

          I1mLUpper = I1mLUpper * dE
          I3mLUpper = I3mLUpper * dE

          if (extIntCase == 1) then
             if (extrapolateExternalF) then
                ! extrapolation down
                IL2(L+1,ix) = IL2Lower - (IL2Upper - IL2Lower) * deltaELower
                IL4(L+1,ix) = IL4Lower - (IL4Upper - IL4Lower) * deltaELower
                I1mL(L+1,ix) = I1mLLower - (I1mLUpper - I1mLLower) * deltaELower
                I3mL(L+1,ix) = I3mLLower - (I3mLUpper - I3mLLower) * deltaELower
             else
                ! assume the integrand is zero below grid
                IL2(L+1,ix) = IL2Lower
                IL4(L+1,ix) = IL4Lower
                I1mL(L+1,ix) = I1mLLower
                I3mL(L+1,ix) = I3mLLower
             end if
             
          else if (extIntCase == 2) then
             if (extrapolateExternalF) then
                ! extrapolation up
                IL2(L+1,ix) = IL2Upper + (IL2Upper - IL2Lower) * deltaEUpper
                IL4(L+1,ix) = IL4Upper + (IL4Upper - IL4Lower) * deltaEUpper
                I1mL(L+1,ix) = I1mLUpper + (I1mLUpper - I1mLLower) * deltaEUpper
                I3mL(L+1,ix) = I3mLUpper + (I3mLUpper - I3mLLower) * deltaEUpper
             else
                ! assume the integrand is zero below grid
                IL2(L+1,ix) = IL2Upper
                IL4(L+1,ix) = IL4Upper
                I1mL(L+1,ix) = I1mLUpper
                I3mL(L+1,ix) = I3mLUpper
             end if
             
          else if (extIntCase == 3) then
             ! interpolation
             IL2(L+1,ix) = IL2Upper * deltaELower + IL2Lower * deltaEUpper
             IL4(L+1,ix) = IL4Upper * deltaELower + IL4Lower * deltaEUpper
             I1mL(L+1,ix) = I1mLUpper * deltaELower + I1mLLower * deltaEUpper
             I3mL(L+1,ix) = I3mLUpper * deltaELower + I3mLLower * deltaEUpper
          end if
          
          beta = (L+1) * (L+2) * (2*L-1)
          beta = beta/(2*L+3)
          d2GLdx2(L+1,ix) &
            = L*(L-1) * x(ix)**(L-2) * I3mL(L+1,ix) &
            + beta * x(ix)**L * I1mL(L+1,ix) &
            + beta * x(ix)**(-L-3) * IL4(L+1,ix) &
            + L*(L-1) * x(ix)**(-L-1) * IL2(L+1,ix)
          d2GLdx2(L+1,ix) = d2GLdx2(L+1,ix) * 4*pi/(1-4*L**2)
       end do
    end do

    
    if ((itheta==1) .and. (izeta==1)) then
       print *,"!!!!!!!!!!!!!"
       print *,x(3)
       print *,ispeciesB, ispeciesA
       print *, "IL2, IL4, I1mL, I3mL"
       print *,IL2(1,3)
       print *,IL4(1,3)
       print *,I1mL(1,3)
       print *,I3mL(1,3)
       print *,"!!!!!!!!!!!!"
    end if

    deallocate(Es)
    
  end subroutine calculateGL

  
  subroutine calculateHL(HL, dHLdx,ispeciesA,ispeciesB,itheta,izeta)
    use globalVariables, only: externalNL, x, pi
    integer, intent(in) :: ispeciesA,ispeciesB,itheta,izeta ! useless once the IL(x) have been calculated
    PetscScalar, dimension(:,:), intent(out) :: HL, dHLdx
    integer :: L

    do L=0,externalNL-1
       HL(L+1,:) &
            = x**L *I1mL(L+1,:) &
            + x**(-L-1) * IL2(L+1,:)
       HL(L+1,:) = HL(L+1,:) * 4*pi/(2*L+1)

       dHLdx(L+1,:) &
            = -(L+1) * x**(-L-2) * IL2(L+1,:) &
            + L * x**(L-1) * I1mL(L+1,:)
       dHLdx(L+1,:) = dHLdx(L+1,:) * 4*pi/(2*L+1)
    end do
    
  end subroutine calculateHL

  
end module externalFCalculations
