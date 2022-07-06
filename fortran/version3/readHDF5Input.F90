module readHDF5Input

  use HDF5

#include "PETScVersions.F90"

  implicit none

  interface readVariable
    module procedure readVariable_integer
    module procedure readVariable_scalar
    module procedure readVariable_1d
    module procedure readVariable_2d
    module procedure readVariable_3d
  end interface readVariable

  integer(HID_T), private :: fileID, groupID

contains

  subroutine openInputFile(fileName, groupName)

    character(len=*), intent(in) :: fileName, groupName
    integer :: HDF5Error

    call h5open_f(HDF5Error) 

    ! Open input file
    call h5fopen_f(trim(fileName), H5F_ACC_RDONLY_F, fileID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening file: ", fileName
      stop
    end if
       
    ! Open group
    call h5gopen_f(fileID, trim(groupName), groupID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening group: ",groupName
      stop
    end if

  end subroutine openInputFile

  subroutine closeInputFile()

    integer :: HDF5Error

    ! Close the group
    call h5gclose_f(groupID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing input group"
      stop
    end if

    ! Close the file
    call h5fclose_f(fileID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing input file"
      stop
    end if

  end subroutine closeInputFile

  subroutine readVariable_integer(variable,varname)

    integer, intent(inout) :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(0) :: dimensions
    integer :: HDF5Error

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable: ",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_INTEGER, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable: ",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable: ",varname
      stop
    end if

  end subroutine readVariable_integer

  subroutine readVariable_scalar(variable,varname)

    PetscScalar, intent(inout) :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(0) :: dimensions
    integer :: HDF5Error

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable: ",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable: ",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable: ",varname
      stop
    end if

  end subroutine readVariable_scalar

  subroutine readVariable_1d(variable,varname)

    PetscScalar, intent(inout), dimension(:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(1) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array: ",varname
      stop
    end if

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable: ",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable: ",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable: ",varname
      stop
    end if

  end subroutine readVariable_1d

  subroutine readVariable_2d(variable,varname)

    PetscScalar, intent(inout), dimension(:,:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(2) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array: ",varname
      stop
    end if

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable: ",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable: ",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable: ",varname
      stop
    end if

  end subroutine readVariable_2d

  subroutine readVariable_3d(variable,varname)

    PetscScalar, intent(inout), dimension(:,:,:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(3) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array: ",varname
      stop
    end if

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable: ",varname
      stop
    end if

    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable: ",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable: ",varname
      stop
    end if

  end subroutine readVariable_3d

  subroutine setPhi1()
    !!Reads Phi1 from an external file and sets it globally

    use globalVariables, only: Phi1Hat, dPhi1Hatdtheta, dPhi1Hatdzeta, externalPhi1Filename, theta, zeta, ddtheta, ddzeta, Ntheta, Nzeta, NPeriods, pi, masterProc!, zetaMax

    integer :: externalNzeta, externalNtheta, externalNIterations, externalNPeriods
    integer :: i, j
    character(len=*), parameter :: line="******************************************************************"

    integer :: itheta, izeta, indexZetaUpper, indexZetaLower, indexThetaUpper, indexThetaLower

    PetscScalar, dimension(:), allocatable :: externalTheta, externalZeta

    PetscScalar, dimension(:,:,:), allocatable :: externalPhi1Hat

    PetscScalar :: intervalZeta, intervalTheta, deltaZetaUpper, deltaZetaLower, deltaThetaUpper, deltaThetaLower, intervalRelTol=1d-4, gridAbsTol=1d-8

    if (masterProc) then
       print *,""
       print *,"-------------------------------------------------------"
       print *,"Reading Phi1Hat from ",externalPhi1Filename
    end if

    call openInputFile(externalPhi1Filename, "/")

    call readVariable(externalNzeta, "Nzeta")
    call readVariable(externalNtheta, "Ntheta")
    call readVariable(externalNIterations, "NIterations")
    call readVariable(externalNPeriods, "NPeriods")

    allocate(externalZeta(externalNzeta))
    allocate(externalTheta(externalNtheta))

    call readVariable(externalZeta, "zeta")
    call readVariable(externalTheta, "theta")


    ! *******************************************************************************
    ! Validate the input data from the external file
    ! ******************************************************************************* 

    if (externalNzeta<1) then !!We allow Nzeta to be 1 in case of toroidal symmetry
       if (masterProc) then
          print *,"Error! Nzeta must be at least 1 in ",externalPhi1Filename
       end if
       stop
    end if

    if (externalNtheta<2) then !!There must be at least two grid points in theta
       if (masterProc) then
          print *,"Error! Ntheta must be at least 2 in ",externalPhi1Filename
       end if
       stop
    end if

    if (externalNIterations<1) then
       if (masterProc) then
          print *,"Error! NIterations must be positive in ",externalPhi1Filename
       end if
       stop
    end if
    
    if (externalNPeriods<1) then
       if (masterProc) then
          print *,"Error! NPeriods must be positive in ",externalPhi1Filename
       end if
       stop
    end if

    if ((size(externalZeta) .ne. externalNzeta) .or. (size(externalTheta) .ne. externalNtheta)) then
       if (masterProc) then
          print *,"Error! The size of the zeta, theta arrays do not agree with the dimension Nzeta,Ntheta in ",externalPhi1Filename
       end if
       stop
    end if

    if ((NPeriods .ne. externalNPeriods) .and. masterProc) then
       print *,line
       print *,line
       print *,"**   WARNING: You are reading an external Phi1Hat "
       print *,"**            from a file with a different number of toroidal periods."
       print *,"**            Is the geometry correct for this Phi1?"
       print *,line
       print *,line
    end if

    if ((externalZeta(1) < 0) .or. (externalZeta(externalNzeta) > 2*pi/externalNPeriods)) then 
       if (masterProc) then
          print *,"Error! All elements of the zeta array must be between 0 and 2*pi/NPeriods in ",externalPhi1Filename
       end if
       stop
    end if

    if ((externalTheta(1) < 0) .or. (externalTheta(externalNtheta) > 2*pi)) then
       if (masterProc) then
          print *,"Error! All elements of the theta array must be between 0 and 2*pi in ",externalPhi1Filename
       end if
       stop
    end if

    if (externalNzeta == 1) then
       intervalZeta = 2*pi/externalNPeriods
    else
       intervalZeta = externalZeta(2) - externalZeta(1)
    end if

    intervalTheta = externalTheta(2) - externalTheta(1)

    if (intervalZeta == 0 .or. intervalTheta == 0) then
       if (masterProc) then
          print *,"Error! The zeta, theta arrays should not contain duplicate numbers in ",externalPhi1Filename
       end if
       stop
    end if

    do j=2,externalNzeta
       if (externalZeta(j-1) >= externalZeta(j)) then
          if (masterProc) then
             print *,"Error! zeta grid points are not sorted in increasing order in ",externalPhi1Filename
          end if
          stop
       end if
       if (abs( (externalZeta(j)-externalZeta(j-1)-intervalZeta)/intervalZeta ) > intervalRelTol) then
          if (masterProc) then
             print *,"Error! zeta grid points must be uniformly spaced in ",externalPhi1Filename
             print *,"Grid points:",externalZeta
          end if
          stop
       end if
    end do

    do i=2,externalNtheta
       if (externalTheta(i-1) >= externalTheta(i)) then
          if (masterProc) then
             print *,"Error! theta grid points are not sorted in increasing order in ",externalPhi1Filename
          end if
          stop
       end if
       if (abs( (externalTheta(i)-externalTheta(i-1)-intervalTheta)/intervalTheta ) > intervalRelTol) then
          if (masterProc) then
             print *,"Error! theta grid points must be uniformly spaced in ",externalPhi1Filename
             print *,"Grid points:",externalTheta
          end if
          stop
       end if
    end do

    ! *******************************************************************************
    ! Read Phi1 from the external file and interpolate to the current theta, zeta grid
    ! *******************************************************************************

    !!FOR SOME REASON READING THE HDF5 DATA HAS TO BE MADE IN REVERSE ORDER
    !allocate(externalPhi1Hat(externalNzeta, externalNtheta, externalNIterations))
    allocate(externalPhi1Hat(externalNIterations, externalNtheta, externalNzeta))

    call readVariable(externalPhi1Hat, "Phi1Hat")
    
    if ((externalNzeta == Nzeta) .and. (externalNtheta == Ntheta) .and. ( abs(externalZeta(1) - zeta(1)) < gridAbsTol/externalNPeriods) .and. ( abs(externalTheta(1) - theta(1)) < gridAbsTol) .and. &
         ( abs(externalZeta(externalNzeta) - zeta(Nzeta)) < gridAbsTol/externalNPeriods) .and. ( abs(externalTheta(externalNtheta) - theta(Ntheta)) < gridAbsTol)) then !The grid is the same so we don't have to interpolate
       !Phi1Hat = transpose(externalPhi1Hat(externalNIterations, :, :))
       Phi1Hat = externalPhi1Hat(externalNIterations, :, :)

    else !!We have to interpolate using 2D linear interpolation
       !!Phi1Hat(Ntheta,Nzeta)

       do izeta = 1,Nzeta
          !We know that the zeta,theta arrays are sorted in increasing order
          ! Set 'index' to point to the first element in the x grid that is
          ! >= to the desired location:

          indexZetaUpper = 0
          indexZetaLower = 0

          do j=1,externalNzeta 
             if (externalZeta(j) >= zeta(izeta)) then
                indexZetaUpper = j
                exit
             end if
          end do

          if (indexZetaUpper == 0) then !This means that zeta(izeta) is larger than all externalZeta, so we must extrapolate periodically
             indexZetaUpper = 1 !We have to shift periodically and read from the lowest point
             indexZetaLower = externalNzeta
             deltaZetaUpper = ( 2*pi/NPeriods + externalZeta(1) - zeta(izeta) ) / (2*pi/NPeriods + externalZeta(1) - externalZeta(externalNzeta))
             deltaZetaLower = ( zeta(izeta) - externalZeta(externalNzeta) ) / (2*pi/NPeriods + externalZeta(1) - externalZeta(externalNzeta))

          else if (indexZetaUpper == 1) then !This means that zeta(izeta) is smaller than (or equal to) all externalZeta, so we must extrapolate periodically
             indexZetaLower = externalNzeta !We have to shift periodically and read from the uppermost point
             deltaZetaUpper = ( externalZeta(1) - zeta(izeta) ) / (2*pi/NPeriods + externalZeta(1) - externalZeta(externalNzeta))
             deltaZetaLower = ( 2*pi/NPeriods + zeta(izeta) - externalZeta(externalNzeta) ) / (2*pi/NPeriods + externalZeta(1) - externalZeta(externalNzeta)) 

          else !We found an upper and a lower interpolation point given by indexZetaUpper and indexZetaUpper-1
             indexZetaLower = indexZetaUpper - 1
             deltaZetaUpper = ( externalZeta(indexZetaUpper) - zeta(izeta) ) / (externalZeta(indexZetaUpper) - externalZeta(indexZetaLower))
             deltaZetaLower = ( zeta(izeta) - externalZeta(indexZetaLower) ) / (externalZeta(indexZetaUpper) - externalZeta(indexZetaLower))

          end if

          if (deltaZetaUpper < 0 .or. deltaZetaLower < 0) then !This is just a check that the coding is correct, the deltas should come out positive
             if (masterProc) then
                print *,"Error! Something is wrong in setPhi1() in readHDF5Input.F90, deltaZetaUpper and deltaZetaLower should be positive."
             end if
             stop
          end if

          do itheta = 1,Ntheta

             indexThetaUpper = 0
             indexThetaLower = 0

             do i=1,externalNtheta
                if (externalTheta(i) >= theta(itheta)) then
                   indexThetaUpper = i
                   exit
                end if
             end do

             if (indexThetaUpper == 0) then !This means that theta(itheta) is larger than all externalTheta, so we must extrapolate periodically
                indexThetaUpper = 1 !We have to shift periodically and read from the lowest point
                indexThetaLower = externalNtheta
                deltaThetaUpper = ( 2*pi + externalTheta(1) - theta(itheta) ) / (2*pi + externalTheta(1) - externalTheta(externalNtheta)) 
                deltaThetaLower = ( theta(itheta) - externalTheta(externalNtheta) ) / (2*pi + externalTheta(1) - externalTheta(externalNtheta))

             else if (indexThetaUpper == 1) then !This means that theta(itheta) is smaller than (or equal to) all externalTheta, so we must extrapolate periodically
                indexThetaLower = externalNtheta !We have to shift periodically and read from the uppermost point
                deltaThetaUpper = ( externalTheta(1) - theta(itheta) ) / (2*pi + externalTheta(1) - externalTheta(externalNtheta))
                deltaThetaLower = ( 2*pi + theta(itheta) - externalTheta(externalNtheta) ) / (2*pi + externalTheta(1) - externalTheta(externalNtheta))

             else !We found an upper and a lower interpolation point given by indexThetaUpper and indexThetaUpper-1
                indexThetaLower = indexThetaUpper - 1
                deltaThetaUpper = ( externalTheta(indexThetaUpper) - theta(itheta) ) / (externalTheta(indexThetaUpper) - externalTheta(indexThetaLower)) 
                deltaThetaLower = ( theta(itheta) - externalTheta(indexThetaLower) ) / (externalTheta(indexThetaUpper) - externalTheta(indexThetaLower))

             end if

             if (deltaThetaUpper < 0 .or. deltaThetaLower < 0) then !This is just a check that the coding is correct, the deltas should come out positive
                if (masterProc) then
                   print *,"Error! Something is wrong in setPhi1() in readHDF5Input.F90, deltaThetaUpper and deltaThetaLower should be positive."
                end if
                stop
             end if

             !Note that deltaThetaUpper goes with indexThetaLower and deltaThetaLower goes with indexThetaUpper in the linear interpolation
             !The corresponding is also true for zeta

             if (externalNzeta == 1) then !There is only one point in zeta, so we can only interpolate in theta
                Phi1Hat(itheta,izeta) = deltaThetaLower*externalPhi1Hat(externalNIterations, indexThetaUpper, 1) + deltaThetaUpper*externalPhi1Hat(externalNIterations, indexThetaLower, 1)

             else !2D interpolation
                Phi1Hat(itheta,izeta) = deltaZetaLower * deltaThetaLower * externalPhi1Hat(externalNIterations, indexThetaUpper, indexZetaUpper) &
                     + deltaZetaUpper * deltaThetaLower * externalPhi1Hat(externalNIterations, indexThetaUpper, indexZetaLower) &
                     + deltaZetaLower * deltaThetaUpper * externalPhi1Hat(externalNIterations, indexThetaLower, indexZetaUpper) &
                     + deltaZetaUpper * deltaThetaUpper * externalPhi1Hat(externalNIterations, indexThetaLower, indexZetaLower)
             end if
             !Phi1Hat(itheta,izeta) = 0
          end do
       end do
    end if

    dPhi1Hatdtheta = matmul(ddtheta,Phi1Hat)
    dPhi1Hatdzeta = transpose(matmul(ddzeta,transpose(Phi1Hat)))

    call closeInputFile()

    if(allocated(externalTheta))deallocate(externalTheta)
    if(allocated(externalZeta))deallocate(externalZeta)
    if(allocated(externalPhi1Hat))deallocate(externalPhi1Hat)

    if (masterProc) then
       !print *,""
       print *,"-------------------------------------------------------" 
       !print *,""
    end if

  end subroutine setPhi1

end module readHDF5Input

