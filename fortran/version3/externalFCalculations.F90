module externalFCalculations

#include "PETScVersions.F90"

  implicit none

  PetscScalar, dimension(:,:), allocatable :: externalFL
  PetscScalar, dimension(:,:), allocatable :: P
  PetscScalar, dimension(:,:), allocatable :: IL2, IL4, I1mL, I3mL
  integer, parameter :: debugtheta = 80
  integer, parameter :: debugzeta = 65
  integer, parameter :: debugix = 4
  integer, parameter :: debugL = 2
  logical, parameter :: printDebug = .false.
  
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
    use globalVariables, only: externalNspecies, externalNL, externalNE, externalNxi, externalXi, externalE, externalF, externalMasses, externalCharges, externalRosenPotentialTerms, nu_n, NHats, mHats, THats, Zs, x2, x, expx2, Ntheta, Nzeta, Nspecies, Nx, pi, zero, masterproc    
    integer :: ispeciesA, ispeciesB, itheta, izeta, L,  ix
    PetscScalar :: externalZ, externalMHat
    PetscScalar, dimension(:), allocatable :: xFactor
    PetscScalar, dimension(:,:), allocatable :: FL, d2GLdx2, HL, dHLdx
    PetscScalar :: debug1, debug2, debug3

    ! global variable
    if (.not. allocated(externalRosenPotentialTerms)) then
       allocate(externalRosenPotentialTerms(Nspecies,Ntheta,Nzeta,externalNL,Nx))
    end if
    externalRosenPotentialTerms = zero 

    if (externalNL > 0) then
       
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

                   L = debugL
                   ix = debugix
                   if ((itheta==debugtheta) .and. (izeta==debugzeta) .and. masterproc .and. printDebug) then
                      print *,"!!!!!!!!!!!!!"
                      print *,L
                      print *,"!!!!!!!!!!!!!"
                      print *,ispeciesB, ispeciesA

                      print *, "d2GLdx2, dHLdx, HL"
                      print *,d2GLdx2(L+1,ix)
                      print *,dHLdx(L+1,ix)
                      print *,HL(L+1,ix)
                      print *,"!!!!!!!!!!!!"
                      debug1 = externalZ**2 * x2(ix) * d2GLdx2(L+1,ix) 
                      debug2 = - externalZ**2 * ((1 - mHats(ispeciesA)/externalMHat) * x(ix) * dHLdx(L+1,ix) + HL(L+1,ix))
                      debug3 = externalZ**2 * 2 * pi  * (mHats(ispeciesA)/externalMHat) * FL(L+1,ix)
                      debug1 = debug1 * xFactor(ix)
                      debug2 = debug2 * xFactor(ix)
                      debug3 = debug3 * xFactor(ix)
                      print *,"col"
                      print *,debug1
                      print *,debug2
                      print *,debug3
                   end if

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

    end if
 
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

    ! Used for comparing 2 methods of calc (NOT USED)
    ! allocate(FL2(externalNL,Nx))

    ! print externalF
    ! if ((itheta == debugtheta) .and. (izeta == debugzeta)) then
    !    do ixi = 1,externalNxi
    !       print *, externalF(ispeciesB,itheta,izeta,ixi,:)
    !    end do
    ! end if
    
    
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
    use globalVariables, only: externalNE, externalNL, externalE,  externalMasses, x, x2, mHats, THats, Nx, pi, masterproc!, extrapolateExternalF
    integer, intent(in) :: ispeciesA,ispeciesB,itheta,izeta
    PetscScalar, dimension(:,:), intent(out) :: d2GLdx2
    PetscScalar :: beta, mHatA, mHatB, THatA
    integer :: Nbins, iE, ix, L
    PetscScalar :: remainder
    PetscScalar :: dE, E0
    PetscScalar, dimension(:), allocatable :: Es

    allocate(Es(Nx))
    
    mHatB = externalMasses(ispeciesB)
    mHatA = mHats(ispeciesA)
    THatA = THats(ispeciesA)
    Es = THatA * x2 * mHatB/mHatA
    dE = externalE(2) - externalE(1)
    E0 = externalE(1) - dE/2
    
    ! calculate the relevant integrals interpolated/extrapolated
    ! to the x grid
    
    ! NOTE: 
    ! "integration limits grid" is offset by dE/2
    ! in our 1D cell-centered finite volume
    do ix=1,Nx
       Nbins = floor((Es(ix)-E0)/dE)
       remainder = (Es(ix)- E0)/dE - Nbins
       if ((itheta==debugtheta) .and. (izeta==debugzeta) .and. masterProc .and. printDebug) then
          print *,"Nbins,Es,E0"
          print *,Nbins,Es(ix),E0
       end if
       
       if (Nbins > externalNE) then
          Nbins = externalNE
          remainder = 0.0
       else if (Nbins < 0) then
          Nbins = 0
          remainder = 0.0
       end if
       
          
       do L=0,externalNL-1          
          IL2(L+1,ix) = 0.0
          IL4(L+1,ix) = 0.0
          do iE = 1,Nbins
             IL2(L+1,ix) = IL2(L+1,ix) + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (L+2))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
             IL4(L+1,ix) = IL4(L+1,ix) + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (L+4))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
          end do
          
          I1mL(L+1,ix) = 0.0
          I3mL(L+1,ix) = 0.0
          !print *,Nbins
          !print *,"!!!!!!!!!!!!!"
          do iE = Nbins+1,externalNE
             !print *,iE
             I1mL(L+1,ix) = I1mL(L+1,ix) + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (1-L))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
             I3mL(L+1,ix) = I3mL(L+1,ix) + (externalFL(L+1,iE) * sqrt(externalE(iE) *mHatA/(mHatB * THatA)) ** (3-L))/(2 * sqrt(externalE(iE) * THatA * mHatB/mHatA))
          end do

          
          if (remainder > 0) then
             IL2(L+1,ix) = IL2(L+1,ix) + remainder * (externalFL(L+1,Nbins + 1) * sqrt(externalE(Nbins+1) *mHatA/(mHatB * THatA)) ** (L+2))/(2 * sqrt(externalE(Nbins+1) * THatA * mHatB/mHatA)) 
             IL4(L+1,ix) = IL4(L+1,ix) + remainder * (externalFL(L+1,Nbins + 1) * sqrt(externalE(Nbins+1) *mHatA/(mHatB * THatA)) ** (L+4))/(2 * sqrt(externalE(Nbins+1) * THatA * mHatB/mHatA))
             
             I1mL(L+1,ix) = I1mL(L+1,ix) - remainder * (externalFL(L+1,Nbins+1) * sqrt(externalE(Nbins+1) *mHatA/(mHatB * THatA)) ** (1-L))/(2 * sqrt(externalE(Nbins+1) * THatA * mHatB/mHatA))
             I3mL(L+1,ix) = I3mL(L+1,ix) - remainder * (externalFL(L+1,Nbins+1) * sqrt(externalE(Nbins+1) *mHatA/(mHatB * THatA)) ** (3-L))/(2 * sqrt(externalE(Nbins+1) * THatA * mHatB/mHatA))
          end if
             
          IL4(L+1,ix) = IL4(L+1,ix) * dE
          IL2(L+1,ix) = IL2(L+1,ix) * dE
          
          I1mL(L+1,ix) = I1mL(L+1,ix) * dE
          I3mL(L+1,ix) = I3mL(L+1,ix) * dE

          
          
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

    L = debugL
    ix = debugix
    if ((itheta==debugtheta) .and. (izeta==debugzeta) .and. masterproc .and. printDebug) then
       print *, "!!!!!!!!!!!!!"
       print *, L
       print *, "!!!!!!!!!!!!!"
       print *, ispeciesB, ispeciesA
       print *, "IL2, IL4, I1mL, I3mL"
       print *, IL2(L+1,ix),IL4(L+1,ix),I1mL(L+1,ix),I3mL(L+1,ix)
       print *, externalFL(L+1,:)
       ! recalculate for this ix
       Nbins = floor((Es(ix)-E0)/dE)
       remainder = (Es(ix)- E0)/dE - Nbins
       if (Nbins > externalNE) then
          Nbins = externalNE
          remainder = 0.0
       else if (Nbins < 0) then
          Nbins = 0
          remainder = 0.0
       end if
       print *,Es(ix), E0 + Nbins * dE + remainder * dE
        
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
