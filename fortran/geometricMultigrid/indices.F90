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
  ! for itheta = 1:Ntheta
  !   for izeta = 1:Nzeta
  !     Enforce quasineutrality
  ! Force <\tilde{\phi}> = 0
  
  ! Order of the vector of unknowns & of columns in the matrix:
  ! --------------------------------
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     for L = 0:(NL-1)
  !       for itheta = 1:Ntheta
  !         for izeta = 1:Nzeta
  !           Distribution function
  ! for itheta = 1:Ntheta
  !   for izeta = 1:Nzeta
  !     \tilde{\phi}
  ! lambda


  
  ! *********************************************************
  ! For constraintScheme==1, 3, or 4:
  ! *********************************************************
  
  ! Order of the rows of the matrix and of the RHS:
  ! --------------------------------
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     for L = 0:(NL-1)
  !       for itheta = 1:Ntheta
  !         for izeta = 1:Nzeta
  !           Enforce the drift-kinetic equation
  ! for itheta = 1:Ntheta
  !   for izeta = 1:Nzeta
  !     Enforce quasineutrality
  ! Force <\tilde{\phi}> = 0
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
  ! for itheta = 1:Ntheta
  !   for izeta = 1:Nzeta
  !     \tilde{\phi}
  ! lambda
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
  ! for itheta = 1:Ntheta
  !   for izeta = 1:Nzeta
  !     Enforce quasineutrality
  ! Force <\tilde{\phi}> = 0
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
  ! for itheta = 1:Ntheta
  !   for izeta = 1:Nzeta
  !     \tilde{\phi}
  ! lambda
  ! for iSpecies = 1:Nspecies
  !   for ix = 1:Nx
  !     source at that x


  module indices

    implicit none

    ! Allowed values for "whichBlock":
    ! --------------------------------------------------

    ! Distribution function, or DKE
    integer, parameter :: BLOCK_F                   = 9990

    ! Quasineutrality condition, or \tilde{\phi}
    integer, parameter :: BLOCK_QN                  = 9991

    ! Enforce <phi_1>=0, or lambda
    integer, parameter :: BLOCK_PHI1_CONSTRAINT     = 9992

    ! Enforce <n> = <n>_desired, or particle source
    integer, parameter :: BLOCK_DENSITY_CONSTRAINT  = 9993

    ! Enforce <p> = <p>_desired, or heat source
    integer, parameter :: BLOCK_PRESSURE_CONSTRAINT = 9994

    ! Enforce <f_1> = 0 at each x
    integer, parameter :: BLOCK_F_CONSTRAINT        = 9995

  contains

    integer function getIndex(level, i_species, i_x, i_xi, i_theta, i_zeta, whichBlock)

      ! This function takes as inputs "local" indices in the species, x, xi, theta, and zeta grids,
      ! and returns the "global" index, i.e. the row or column of the master matrix.
      
      ! The input "local" indices (i_x, i_xi, etc) are 1-based, since small matrices are stored in Fortran.
      ! the output "global" index is 0-based, since PETSc uses 0-based indexing.

      use globalVariables, only: levels, Nspecies, Nx, constraintScheme, includePhi1, masterProc

      implicit none

      integer, intent(in) :: level, i_species, i_x, i_xi, i_theta, i_zeta, whichBlock
      integer :: Ntheta, Nzeta, Nxi, matrixSize

      Ntheta = levels(level)%Ntheta
      Nzeta  = levels(level)%Nzeta
      Nxi    = levels(level)%Nxi
      matrixSize = levels(level)%matrixSize

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

      if (i_xi > levels(level)%Nxi_for_x(i_x)) then
         print *,"Error: i_xi > Nxi_for_x(i_x)"
         print *,"i_xi:",i_xi
         print *,"i_x:",i_x
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
      
      if (whichBlock==BLOCK_QN .and. (.not. includePhi1)) then
         if (masterProc) then
            print *,"Error! whichBlock==BLOCK_QN but includePhi1=.false."
         end if
         stop
      end if

      if (whichBlock==BLOCK_PHI1_CONSTRAINT .and. (.not. includePhi1)) then
         if (masterProc) then
            print *,"Error! whichBlock==BLOCK_PHI1_CONSTRAINT but includePhi1=.false."
         end if
         stop
      end if

      ! Done with validation.

      select case (whichBlock)
      case (BLOCK_F)
!!$         getIndex = (i_species-1)*Nx*Nxi*Ntheta*Nzeta &
!!$              +(i_x-1)*Nxi*Ntheta*Nzeta &
!!$              +(i_xi-1)*Ntheta*Nzeta &
!!$              +(i_theta-1)*Nzeta &
!!$              +i_zeta -1
         getIndex = (i_species-1)*levels(level)%DKE_size &
              +levels(level)%first_index_for_x(i_x)*Ntheta*Nzeta &
              +(i_xi-1)*Ntheta*Nzeta &
              +(i_theta-1)*Nzeta &
              +i_zeta -1

      case (BLOCK_QN)
         !getIndex = Nspecies*Nx*Nxi*Ntheta*Nzeta &
         getIndex = Nspecies*levels(level)%DKE_size &
              +(i_theta-1)*Nzeta &
              +i_zeta-1

      case (BLOCK_PHI1_CONSTRAINT)
         !getIndex = Nspecies*Nx*Nxi*Ntheta*Nzeta &
         getIndex = Nspecies*levels(level)%DKE_size &
              +Ntheta*Nzeta

      case (BLOCK_DENSITY_CONSTRAINT)
         if ((constraintScheme .ne. 1) .and. (constraintScheme .ne. 3) .and. (constraintScheme .ne. 4)) then
            print *,"Error! whichBlock=BLOCK_DENSITY_CONSTRAINT requires constraintScheme = 1, 3, or 4."
            print *,"Now, constraintScheme = ",constraintScheme
            stop
         end if
         if (includePhi1) then
            getIndex = Nspecies*levels(level)%DKE_size &
                 + Ntheta*Nzeta + 1 &
                 + (i_species-1)*2
         else
            getIndex = Nspecies*levels(level)%DKE_size &
                 + (i_species-1)*2
         end if

      case (BLOCK_PRESSURE_CONSTRAINT)
         if ((constraintScheme .ne. 1) .and. (constraintScheme .ne. 3) .and. (constraintScheme .ne. 4)) then
            print *,"Error! whichBlock=BLOCK_PRESSURE_CONSTRAINT requires constraintScheme = 1, 3, or 4."
            print *,"Now, constraintScheme = ",constraintScheme
            stop
         end if
         if (includePhi1) then
            getIndex = Nspecies*levels(level)%DKE_size &
                 + Ntheta*Nzeta + 1 &
                 + (i_species-1)*2 + 1
         else
            getIndex = Nspecies*levels(level)%DKE_size &
                 + (i_species-1)*2 + 1
         end if

      case (BLOCK_F_CONSTRAINT)
         if (constraintScheme .ne. 2) then
            print *,"Error! whichBlock=BLOCK_F_CONSTRAINT requires constraintScheme=2"
            print *,"Now, constraintScheme = ",constraintScheme
            stop
         end if
         if (includePhi1) then
            getIndex = Nspecies*levels(level)%DKE_size &
                 + Ntheta*Nzeta + 1 &
                 + (i_species-1)*Nx + i_x - 1
         else
            getIndex = Nspecies*levels(level)%DKE_size &
                 + (i_species-1)*Nx + i_x - 1
         end if

      case default
         print *,"Error: Invalid value for whichBlock:", whichBlock
         stop
      end select

      ! One last sanity check:

      if (getIndex < 0) then
         print *,"Error! Something went wrong, and the index came out less than 0."
         stop
      end if

      if (getIndex >= matrixSize) then
         print *,"Error! Something went wrong, and the index came out larger than the matrix size."
         stop
      end if

    end function getIndex

    ! ---------------------------------------------------------------
    ! ---------------------------------------------------------------

    subroutine computeMatrixSize(level)

      use globalVariables, only: levels, Nspecies, Nx, includePhi1, constraintScheme, masterProc, matrixSize

      implicit none

      integer, intent(in) :: level
      integer :: ix

      Ntheta = levels(level)%Ntheta
      Nzeta  = levels(level)%Nzeta
      Nxi    = levels(level)%Nxi

      levels(level)%DKE_size = sum(levels(level)%Nxi_for_x)*Ntheta*Nzeta

      levels(level)%matrixSize = Nspecies * levels(level)%DKE_size
      !matrixSize = Nspecies * Ntheta * Nzeta * Nxi * Nx

      allocate(levels(level)%first_index_for_x(Nx))
      levels(level)%first_index_for_x(1) = 0
      do ix=2,Nx
         levels(level)%first_index_for_x(ix) = levels(level)%first_index_for_x(ix-1) + levels(level)%Nxi_for_x(ix-1)
      end do


      if (includePhi1) then
         ! Include quasineutrality at Ntheta*Nzeta grid points, plus the <Phi_1>=0 constraint:
         levels(level)%matrixSize = levels(level)%matrixSize + Ntheta*Nzeta + 1
      end if

      select case (constraintScheme)
      case (0)
      case (1,3,4)
         levels(level)%matrixSize = levels(level)%matrixSize + 2 * Nspecies
      case (2)
         levels(level)%matrixSize = levels(level)%matrixSize + Nx * Nspecies
      case default
         print *,"Error! Invalid constraintScheme"
         stop
      end select

      if (level==1) matrixSize = levels(level)%matrixSize

      if (masterProc) print "(a,i3,a,i10,a,i10)"," The matrix for level",level," is",levels(level)%matrixSize," x",levels(level)%matrixSize," elements."

    end subroutine computeMatrixSize

  end module indices
