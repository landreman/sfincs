!> Populates right hand side of the adjoint equation.
!! @param rhs Vector to be populated.
!! @param whichAdjointRHS Indicates which integrated quantity is being differentiated. If = 1 (particle flux), = 2 (heat flux), = 3 (bootstrap current).
!! @param whichSpecies whichSpecies Indicates species used for inner product. If = 0, corresponds to a species-summed quantity. If nonzero, indicates number of species.
!! @param fineGrid Indicates whether the RHS is evaluated on the fine grid. If 1, corresponds to the fine grid. Otherwise referes to the coarse grid.
!! fineGrid does not control anything at this point, but is included for future incorporation of adjointEC
subroutine populateAdjointRHS(rhs, whichAdjointRHS, whichSpecies, fineGrid)

#include "PETScVersions.F90"


  use globalVariables
  use indices
  use petscvec

  implicit none

  Vec :: rhs
  integer :: whichAdjointRHS, whichSpecies, fineGrid

  PetscErrorCode :: ierr
  integer :: ix, L, itheta, izeta, ispecies, index
  PetscScalar :: THat, mHat, sqrtTHat, sqrtMHat, xPartOfRHS, factor, nHat, Z, sqrtFSAB2
  PetscScalar :: xPartOfRHSExB, factorExB
  PetscReal :: norm, factor1, factor2
  logical, parameter :: localUsePhi1 = .true.


  ! Validate input
  if (whichAdjointRHS<1 .or. whichAdjointRHS>3 .or. whichSpecies<0 .or. whichSpecies>Nspecies) then
    if (masterProc) then
      print *,"whichAdjointRHS: ", whichAdjointRHS
      print *,"whichSpecies: ", whichSpecies
      print *,"Error! Incorrect input to populateAdjointRHS"
    end if
    stop
  end if

  call VecSet(rhs, zero, ierr)


  select case (whichAdjointRHS)

    case (1) ! particle flux/radial current

    do ispecies=1,Nspecies
      ! For whichSpecies == 0, RHS corresponds to radial current - species summed weighted by Zs
      if (whichSpecies == ispecies .or. whichSpecies == 0) then
        THat = THats(ispecies)
        mHat = mHats(ispecies)
        nHat = nHats(ispecies)
        sqrtTHat = sqrt(THat)
        sqrtMHat = sqrt(mHat)
        Z = Zs(ispecies)

        do ix=ixMin,Nx
            if (discreteAdjointOption .eqv. .false.) then
                xPartOfRHS = exp(-x2(ix))*nHat*mHat*sqrtMHat*Delta*x2(ix)/ &
                    (pi*sqrtpi*THat*sqrtTHat*Z)*ddrN2ddpsiHat
            else
               xPartOfRHS = pi*Delta*THat*THat*sqrtTHat*x2(ix)*x2(ix)*ddrN2ddpsiHat*xweights(ix)/(Z*VPrimeHat*mHat*sqrtMHat)
               if (includePhi1 .and. localUsePhi1) then
                  xPartOfRHSExB = 2*pi*Delta*alpha*THat*sqrtTHat*x2(ix)*ddrN2ddpsiHat*xweights(ix)/(VPrimeHat*mHat*sqrtMHat)
               else
                  xPartOfRHSExB = 0
               end if
            end if
            if (whichSpecies == 0) then
               xPartOfRHS = xPartOfRHS*Z
               xPartOfRHSExB = xPartOfRHSExB*Z
            end if

            do itheta = ithetaMin,ithetaMax
                do izeta = izetaMin,izetaMax
                   if (.false.) then
                    factor = xPartOfRHS/(BHat(itheta,izeta)*BHat(itheta,izeta)*BHat(itheta,izeta)) &
                        *(BHat_sub_theta(itheta,izeta)*dBHatdzeta(itheta,izeta) &
                        - BHat_sub_zeta(itheta,izeta)*dBHatdtheta(itheta,izeta))
                 else
                    factor = 0
                 end if
                 
                    if (includePhi1 .and. localUsePhi1) then
                       factorExB = xPartOfRHSExB/(BHat(itheta,izeta)*BHat(itheta,izeta)) &
                            *(BHat_sub_theta(itheta,izeta)*dPhi1Hatdzeta(itheta,izeta) &
                            - BHat_sub_zeta(itheta,izeta)*dPhi1Hatdtheta(itheta,izeta))
                    else
                       factorExB = 0
                    end if
                    if (discreteAdjointOption .eqv. .false.) then
                        factor = factor*DHat(itheta,izeta)
                        factor1 = four/three
                        factor2 = two/three
                     else
                        factor = factor*thetaWeights(itheta)*zetaWeights(izeta)
                        factorExB = factorExB*thetaWeights(itheta)*zetaWeights(izeta)
                        factor1 = 8/three
                        factor2 = four/15
                end if

                L = 0

                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, factor1*factor + factorExB, INSERT_VALUES, ierr) ! *2

                L = 2

                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, factor2*factor, INSERT_VALUES, ierr)

             enddo
            enddo
          enddo
        end if
      enddo

    case (2) ! heat flux

    do ispecies=1,Nspecies
      ! For whichSpecies==0, RHS corresponds to total heat flux - species summed
      if (whichSpecies == ispecies .or. whichSpecies == 0) then
        THat = THats(ispecies)
        mHat = mHats(ispecies)
        nHat = nHats(ispecies)
        sqrtTHat = sqrt(THat)
        sqrtMHat = sqrt(mHat)
        Z = Zs(ispecies)

        do ix=ixMin,Nx
          if (discreteAdjointOption) then
            xPartOfRHS = pi*Delta*x2(ix)*x2(ix)*x2(ix)*xweights(ix)*THat*THat*THat*sqrtTHat*ddrN2ddpsiHat/(2*Z*VPrimeHat*mHat*sqrtMHat)
          else
            xPartOfRHS = exp(-x2(ix))*nHat*mHat*sqrtMHat*Delta*x2(ix)*x2(ix)/ &
              (2*pi*sqrtpi*sqrtTHat*Z)*ddrN2ddpsiHat
          end if
          do itheta = ithetaMin,ithetaMax
             do izeta = izetaMin,izetaMax
                if (discreteAdjointOption .eqv. .false.) then
                    factor = xPartOfRHS*DHat(itheta,izeta)/(BHat(itheta,izeta)**3)*(BHat_sub_theta(itheta,izeta) &
                      *dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)* &
                    dBHatdtheta(itheta,izeta))
                    factor1 = four/three
                    factor2 = two/three
                else
                    factor = xPartOfRHS*thetaWeights(itheta)*zetaWeights(izeta)*(BHat_sub_theta(itheta,izeta) &
                      *dBHatdzeta(itheta,izeta) - BHat_sub_zeta(itheta,izeta)* &
                    dBHatdtheta(itheta,izeta))/(BHat(itheta,izeta)**3)
                    factor1 = 8/three
                    factor2 = four/15
                end if
                ! factor is the samed for species-summed heat flux

                L = 0
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, factor1*factor, INSERT_VALUES, ierr)

                L = 2
                index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                call VecSetValue(rhs, index, factor2*factor, INSERT_VALUES, ierr)

             end do
            end do
          end do
        end if
      end do

    case (3) ! Bootstrap Current/Flow

      sqrtFSAB2 = sqrt(FSABHat2)
      do ispecies=1,Nspecies
        if (whichSpecies==ispecies .or. whichSpecies == 0) then
          THat = THats(ispecies)
          mHat = mHats(ispecies)
          nHat = nHats(ispecies)
          Z = Zs(ispecies)

          do ix=ixMin,Nx
             if (discreteAdjointOption .eqv. .false.) then
                  xPartOfRHS = exp(-x2(ix))*x(ix)*2*mHat/(THat*THat*pi*sqrtPi*sqrtFSAB2)
             else
                  xPartOfRHS = 4*pi*THat*THat*x(ix)*x(ix)*x(ix)*xWeights(ix)/(three*mHat*mHat*sqrtFSAB2*VPrimeHat*nHat)
             end if

            if (whichSpecies == 0) then
              xPartOfRHS = xPartOfRHS*Z
            end if
            do itheta = ithetaMin,ithetaMax
               do izeta = izetaMin,izetaMax
                  if (discreteAdjointOption .eqv. .false.) then
                      factor = xPartOfRHS*BHat(itheta,izeta)
                  else
                      factor = xPartOfRHS*thetaWeights(itheta)*zetaWeights(izeta)*BHat(itheta,izeta)/DHat(itheta,izeta)
                  end if 

                  L = 1
                  index = getIndex(ispecies, ix, L+1, itheta, izeta, BLOCK_F)
                  call VecSetValue(rhs, index, factor, INSERT_VALUES, ierr)

               enddo
            enddo
          enddo
        end if
      enddo

  end select

  ! Done inserting values.
  ! Finally, assemble the RHS vector:
  call VecAssemblyBegin(rhs, ierr)
  call VecAssemblyEnd(rhs, ierr)

end subroutine populateAdjointRHS
