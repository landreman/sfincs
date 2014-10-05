  ! *********************************************************
  ! For constraintScheme==0,
  ! *********************************************************
  
  ! Order of the rows of the matrix and of the RHS:
  ! --------------------------------
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     for L = 0:(NL-1)
  !       for itheta = 1:Ntheta
  !         for izeta = 1:Nzeta
  !           Enforce the drift-kinetic equation
  
  ! Order of the vector of unknowns & of columns in the matrix:
  ! --------------------------------
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     for L = 0:(NL-1)
  !       for itheta = 1:Ntheta
  !         for izeta = 1:Nzeta
  !           Enforce the drift-kinetic equation
  
  ! *********************************************************
  ! For constraintScheme==1,
  ! *********************************************************
  
  ! Order of the rows of the matrix and of the RHS:
  ! --------------------------------
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     for L = 0:(NL-1)
  !       for itheta = 1:Ntheta
  !         for izeta = 1:Nzeta
  !           Enforce the drift-kinetic equation
  ! for iSpecies = 1:Nspecies
  !   Force <n_1> = 0
  !   Force <p_1> = 0
  
  
  ! Order of the vector of unknowns & of columns in the matrix:
  ! --------------------------------
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     for L = 0:(NL-1)
  !       for itheta = 1:Ntheta
  !         for izeta = 1:Nzeta
  !           Enforce the drift-kinetic equation
  ! for iSpecies = 1:Nspecies
  !   particle source
  !   energy source

  ! *********************************************************
  ! For constraintScheme==2,
  ! *********************************************************

  ! Order of the rows of the matrix and of the RHS:
  ! --------------------------------
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     for L = 0:(NL-1)
  !       for itheta = 1:Ntheta
  !         for izeta = 1:Nzeta
  !           Enforce the drift-kinetic equation
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     Force <f_1> = 0 at that x


  ! Order of the vector of unknowns & of columns in the matrix:
  ! --------------------------------
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     for L = 0:(NL-1)
  !       for itheta = 1:Ntheta
  !         for izeta = 1:Nzeta
  !           Enforce the drift-kinetic equation
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     source at that x


  module indices

    implicit none

  contains

    integer function getIndex(i_species, i_x, i_xi, i_theta, i_zeta, f_or_sources)

      ! This function takes as inputs "local" indices in the species, x, xi, theta, and zeta grids,
      ! and returns the "global" index, i.e. the row or column of the master matrix.
      
      ! Allowed values for f_or_sources:
      ! 0: distribution function or DKE
      ! 1: particle source or <n_1> = 0
      ! 2: heat source or <p_1> = 0
      ! 3: <f_1> = 0 at each x
      
      ! The input "local" indices (i_x, i_xi, etc) are 1-based, since small matrices are stored in Fortran.
      ! the output "global" index is 0-based, since PETSc uses 0-based indexing.

      use globalVariables

      implicit none

      integer, intent(in) :: i_species, i_x, i_xi, i_theta, i_zeta, f_or_sources

      ! Validate inputs:

      if (i_species < 1) then
         print *,"Error: i_species < 1"
         stop
      end if

      if (i_species > Nspecies) then
         print *,"Error: i_species > Nspecies"
         stop
      end if

      if (i_x < 1) then
         print *,"Error: i_x < 1"
         stop
      end if

      if (i_x > Nx) then
         print *,"Error: i_x > Nx"
         stop
      end if

      if (i_xi < 1) then
         print *,"Error: i_xi < 1"
         stop
      end if

      if (i_xi > Nxi) then
         print *,"Error: i_xi > Nxi"
         stop
      end if

      if (i_theta < 1) then
         print *,"Error: i_theta < 1"
         stop
      end if

      if (i_theta > Ntheta) then
         print *,"Error: i_theta > Ntheta"
         stop
      end if

      if (i_zeta < 1) then
         print *,"Error: i_zeta < 1"
         stop
      end if

      if (i_zeta > Nzeta) then
         print *,"Error: i_zeta > Nzeta"
         stop
      end if

      ! Done with validation.

      select case (f_or_sources)
      case (0)
         getIndex = (i_species-1)*Nx*Nxi*Ntheta*Nzeta &
              +(i_x-1)*Nxi*Ntheta*Nzeta &
              +(i_xi-1)*Ntheta*Nzeta &
              +(i_theta-1)*Nzeta &
              +i_zeta -1

      case (1)
         if (constraintScheme .ne. 1) then
            print *,"Error! f_or_sources=1 requires constraintScheme=1"
            print *,"Now, constraintScheme = ",constraintScheme
            stop
         end if
         getIndex = Nspecies*Nx*Nxi*Ntheta*Nzeta &
              + (i_species-1)*2

      case (2)
         if (constraintScheme .ne. 1) then
            print *,"Error! f_or_sources=2 requires constraintScheme=1"
            print *,"Now, constraintScheme = ",constraintScheme
            stop
         end if
         getIndex = Nspecies*Nx*Nxi*Ntheta*Nzeta &
              + (i_species-1)*2 + 1

      case (3)
         if (constraintScheme .ne. 2) then
            print *,"Error! f_or_sources=3 requires constraintScheme=2"
            print *,"Now, constraintScheme = ",constraintScheme
            stop
         end if
         getIndex = Nspecies*Nx*Nxi*Ntheta*Nzeta &
              + (i_species-1)*Nx + i_x - 1

      case default
         print *,"Error: f_or_sources must be 0, 1, 2, or 3."
         stop
      end select

      ! One last sanity check:

      if (getIndex < 0) then
         print *,"Error! Something went wrong, and the index came out less than 1."
         stop
      end if

      if (getIndex >= matrixSize) then
         print *,"Error! Something went wrong, and the index came out larger than the matrix size."
         stop
      end if

    end function getIndex

!!$    subroutine getIndicesThetaRange(ispecies, ix, ixi, izeta, outputIndices)
!!$      integer, intent(in) :: ispecies, ix, ixi, izeta
!!$      integer, dimension(:), intent(out) :: outputIndices
!!$      integer :: itheta
!!$
!!$      do itheta=1,Ntheta
!!$         outputIndices(itheta) = getIndices(ispecies, ix, ixi, itheta, izeta, 0)
!!$      end do
!!$    end subroutine getIndicesThetaRange

  end module indices
