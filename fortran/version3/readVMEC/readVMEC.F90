module readVMEC

  use pisa_vmec_module

  implicit none

  !type(vmec_fort8_data) :: vmec
  type(vmec_eqdata) :: vmec

  contains

  subroutine read_VMEC(filename)

    implicit none

    character(200), intent(in) :: filename

    type(vmec_wout_data) :: the_txt_data
    type(vmec_wout_cdf_data) :: the_cdf_data
    type(file_type) :: my_filetype

    call check_wout_file_type(filename, my_filetype)
    if( my_filetype%typo == "txt")then
       !print *,"VMEC file is in ASCII format."
       call read_wout_txt_data(filename, the_txt_data)
       !vmec = extract_fort8_data_from_wout_data(the_txt_data)
       vmec = extract_vmec_eqdata_from_wout_data(the_txt_data)
    elseif( my_filetype%typo == "cdf")then
       !print *,"VMEC file is in netCDF format."
#ifdef USE_NETCDF
       call read_wout_cdf_data(filename, the_cdf_data)
       !vmec = extract_fort8_data_from_wout_cdf_data(the_cdf_data)
       vmec = extract_vmec_eqdata_from_wout_cdf_data(the_cdf_data)
#else
       print *,"Error! You are attempting to read a NetCDF file, but you did not set USE_NETCDF in the makefile."
       stop
#endif
    endif

!!$    print *,"VMEC data:"
!!$    print *,"nfp = ",vmec%nfp
!!$    print *,"mnmax = ",vmec%mnmax
!!$    print *,"VMEC dimensions:"
!!$    print *,"radial   : ns   = ",vmec%ns
!!$    print *,"poloidal : mpol = ",vmec%mpol
!!$    print *,"toroidal : ntor = ",vmec%ntor

  end subroutine read_VMEC
end module readVMEC

