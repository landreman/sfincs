MODULE ezcdf_opncls
 
  interface cdfOpn
     module procedure ezcdf_open
  end interface

  interface cdf_open
     module procedure ezcdf_open
  end interface

  interface cdfCls
     module procedure ezcdf_close
  end interface

  interface cdf_close
     module procedure ezcdf_close
  end interface

CONTAINS
 
  subroutine ezcdf_open(ncid,filename,opt,ier)
    ! Create/open cdf file
    ! 03/09/99 C.Ludescher
    !
    !  opt values:
    !    "R" -- readonly, existing file
    !    "W" -- write new file
    !    "M" -- "Modify" -- read/write data in existing file
    !                       (but structure of file does not change)
    !    "A" -- "Append" -- add new structure (define new data items)
    !                       in existing file.  Existing items can also
    !                       be read/written, but only after *new* items
    !                       are first defined and then written.
    !    "V" -- "Verify/Append/Write" -- first try opening in "A" mode; if
    !                       that fails open in "W" mode.
    !
    !  for both "W" and "A" modes, the file is opened in "define data mode".
    !
    include "netcdf.inc"
    INTEGER,       intent(out) :: ncid
    character*(*), intent(in) :: filename
    character*1,   intent(in) :: opt
    integer, optional, intent(out) :: ier
    integer :: status
 
    if (opt .eq. 'w' .or. opt .eq. 'W') then
       ! New file... start in "define mode".
       status = nf_create(filename,IOR(NF_CLOBBER,NF_64BIT_OFFSET),ncid)
       call ezcdf_error(status,filename,'cdfcrt','nf_create')
    else if (opt .eq. 'm' .or. opt .eq. 'M') then
       ! Open existing file for read/write modifications...
       status = nf_open(filename, nf_write, ncid)
       call ezcdf_error(status,filename,'cdfopn','nf_open')
    else if (opt .eq. 'a' .or. opt .eq. 'A') then
       ! Open existing file for read/write modifications...
       status = nf_open(filename, nf_write, ncid)
       call ezcdf_error(status,filename,'cdfopn','nf_open')
       if( status .eq. NF_NOERR ) then
          status = nf_redef(ncid)  ! start in "define mode".
          call ezcdf_error(status,filename,'cdfopn','nf_redef')
       endif
    else if (opt .eq. 'v' .or. opt .eq. 'V') then
       ! Open existing file for read/write modifications...
       status = nf_open(filename, nf_write, ncid)
       if( status .eq. NF_NOERR ) then
          status = nf_redef(ncid)  ! start in "define mode".
       else
          status = nf_create(filename,IOR(NF_CLOBBER,NF_64BIT_OFFSET),ncid)  ! start new file.
       endif
       call ezcdf_error(status,filename,'cdfopn','nf_open')
    else
       ! Open for readonly
       status = nf_open(filename, nf_nowrite, ncid)
       call ezcdf_error(status,filename,'cdfopn','nf_open')
    end if
    if (PRESENT (ier)) then
       if (status .ne. NF_NOERR) then
          ier = 1
       else
          ier = 0
       endif
    endif
    return
  end subroutine ezcdf_open
 

  subroutine ezcdf_close(ncid, ier)
    include "netcdf.inc"
    INTEGER, INTENT(in) ::  ncid
    integer, optional,         intent(out) :: ier
    INTEGER status
    status = nf_close(ncid)
    call ezcdf_error(status,' ','cdfcls','nf_close')
    if(PRESENT (ier)) ier = status
  end subroutine ezcdf_close
 
END MODULE ezcdf_opncls
