!> Module io_module provides utilities and definitions for input and
!! output operations.
!!
!! @author Joachim Geiger
!! @version 1.1
module pisa_io
implicit none

integer, parameter :: eof=-1      !< end of file reached
integer, parameter :: eol=-2      !< end of line reached
logical, parameter :: l_show=.false.  !< test flag

private :: l_show

contains

!----------------------------------------------------------------------
!>find a unit number available for i/o action
!!taken from the book "Object oriented programming via Fortran 90/95"
!! by Ed Akin, p. 305
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!! History:
!! - 1.0: loop through all allowed units was only done if last_unit > 0
!!        and not if the unit was already open.
!! - 1.1: added loop through all allowed units in case the unit was open!
!! - 1.2: logical variable open changed to is_open (keyword).
!!        change integer variable count to icount (keyword).
!!
function get_next_io_unit () result (next)
  implicit none
  integer :: next   !< the next available unit number
  integer, parameter :: min_unit = 10, max_unit = 999
  integer, save      :: last_unit = 0   !< initalize
  integer            :: icount          !< number of failures
  logical            :: is_open         !< file status

  icount = 0 ; next = min_unit - 1
  if ( last_unit > 0 ) then ! check next in line
    next = last_unit + 1
    inquire (unit=next, opened=is_open)
    if ( .not. is_open )then
      last_unit = next ! found it
      return
    end if ! found it
  end if ! last_unit
! next unit was not open => loop over units
  do ! forever
    next = next + 1
    inquire (unit=next, opened=is_open)
    if ( .not. is_open ) then
      last_unit = next     ! found it
      exit ! the unit loop
    end if
    if ( next == max_unit ) then ! attempt reset 3 times
      last_unit = 0
      icount     = icount + 1
      if ( icount <= 3 ) next = min_unit - 1
    end if ! reset try
    if ( next > max_unit ) then ! abort
      print *,'ERROR: max unit exceeded in get_next_io_unit'
      stop    'ERROR: max unit exceeded in get_next_io_unit'
    end if ! abort
  end do ! over unit numbers
end function get_next_io_unit

!----------------------------------------------------------------------
!> gets dimension for matrix import of a file
subroutine check_file(filename, maxcol, maxrow,ierr)
implicit none

character(len=*), intent(in) :: filename
integer, intent(out) :: maxcol   !< number of columns on output
integer, intent(out) :: maxrow   !< number of rows on output
integer, intent(out) :: ierr     !< error flag
character(len=1) :: char   !< readin character
integer :: readstat !< status of read operation
integer :: m,n      !< counter
integer :: inunit   !< unit from which is read

maxcol=0  ! Added by MJL 20160206
m=0; n=0
inunit = get_next_io_unit()
if(l_show) write(6,*)"open unit ",inunit
open(unit=inunit,file=filename,action='read',status='old',iostat=ierr)
if(ierr == 0)then
 do
  read(inunit,'(a)',iostat=readstat,advance='no')char
  m=m+1

  if(readstat == eol)then
   n=n+1
   if(m >= maxcol) maxcol=m
   m=0
  elseif(readstat == eof)then
   maxrow=n
   exit
  elseif(readstat /= 0)then
   exit
  endif
 enddo
 if(l_show) write(6,*)maxrow,' lines in file with a line length <= ', maxcol
 maxcol=maxcol+1  !add end of line character.
 close(unit=inunit)
endif
end subroutine check_file

!----------------------------------------------------------------------
!> imports a file as a string array
subroutine import_file_mat(filename, filestring, maxcol, maxrow,ierr)
implicit none

character(len=*), intent(in) :: filename  !< file to be read
integer, intent(in)  :: maxcol  !< contains maximum number of columns to read
integer, intent(in)  :: maxrow  !< contains number of lines to read
integer, intent(out) :: ierr    !< error flag
character(len=1), dimension(maxcol,maxrow), intent(out) :: filestring  !< 2-D character matrix to hold file content
integer :: readstat !< status of read operation
integer :: m,n      !< counter
integer :: inunit   !< input unit

inunit = get_next_io_unit()
open(unit=inunit,file=filename,action='read',status='old',iostat=ierr)
 n=1
 m=1
 filestring = ' '
 do n=1,maxrow
  do m=1,maxcol
   read(inunit,'(a)',iostat=readstat,advance='no')filestring(n,m)
   if(readstat == eol) exit
   if(readstat == eof)then
    write(6,*) "Premature exit from import file!"
    exit
   endif
  enddo
 enddo
 close(unit=inunit)
end subroutine import_file_mat

!----------------------------------------------------------------------
!> Imports the file into an array of size maxcol with strings of length maxrow.
subroutine import_file_ar(filename, filestring, maxcol, maxrow,ierr)
implicit none

character(len=*), intent(in) :: filename  !< file to be read
integer, intent(in)  :: maxcol  !< contains maximum number of columns to read
integer, intent(in)  :: maxrow  !< contains number of lines to read
integer, intent(out) :: ierr    !< error flag
character(len=*), dimension(maxrow), intent(out) :: filestring  !< string array
integer :: readstat !< status of read operation
integer :: m,n      !< counter
integer :: inunit   !< input unit

inunit = get_next_io_unit()
open(unit=inunit,file=filename,action='read',status='old',iostat=ierr)
 n=1
 m=1
 filestring = ' '
 do n=1,maxrow
  do m=1,maxcol
   read(inunit,'(a)',iostat=readstat,advance='no')filestring(n)(m:m)
   if(readstat == eol) exit
   if(readstat == eof)then
    write(6,*) "Premature exit from import file!"
    exit
   endif
  enddo
 enddo
 close(unit=inunit)
end subroutine import_file_ar

!----------------------------------------------------------------------
!> Counts the number of characters in the current line.
!! Resets the position in file, so that the line can be read again
!! by the calling program, if go_back=.T. (backspace).
function number_of_char_in_line(iunit,go_back) result(n)
implicit none
integer, intent(in) :: iunit   !< input unit from which to read
logical, intent(in), optional :: go_back !< set file position to read line again
integer :: n        !< number of characters in the line
integer :: readstat !< status of read operation
character :: c      !< dummy character

n=0
do
  read(iunit, '(a)', iostat=readstat, advance='no')c
  if(readstat == eol) exit
  n=n+1
enddo
if(.not.(present(go_back)).or. go_back) backspace(iunit)

end function number_of_char_in_line

!----------------------------------------------------------------------
!> Read an double array of size ns from an array of strings
subroutine pisa_read_dar_from_stringar(array,ns,filestring,l,maxrow)
use kind_defs
implicit none
integer, intent(in) :: ns, maxrow   !< array dimensions
integer, intent(inout) :: l   !< on entry : start string array index
                              !! on return: last read string array index
real(DP), dimension(ns), intent(out) :: array !< array to be filled
character(len=*), dimension(maxrow), intent(in) :: filestring  !< string array from which is read
integer :: lb,ub,ierr
integer :: i

lb=1 ; ub=1; ierr=0
do 
 read(filestring(l)(1:len_trim(filestring(l))),*,iostat=ierr) (array(i),i=lb,ub)
 if(ierr == eof) then
  l=l+1   !end of line reached next filestring-line
  lb=ub
  ub=lb
 else
  ub= ub+1
 endif
 if(ub == ns+1) exit ! next array to be read
enddo
end subroutine pisa_read_dar_from_stringar

end module pisa_io
