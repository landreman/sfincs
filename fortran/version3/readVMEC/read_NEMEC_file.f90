! Changes by MJL from Erika Strumberger's original readnemec.f90
! to match Joachim Geiger's routines:

! hiota -> iotas
! hphi -> phi  (Despite the name hphi, it appeared to be on the full grid.)

! frmnc -> rmnc
! frmns -> rmns
! fzmnc -> zmnc
! fzmns -> zmns
! flmnc -> lmnc
! flmns -> lmns

! hbumnc_up -> bsupumnc
! hbumns_up -> bsupumns
! hbvmnc_up -> bsupvmnc
! hbvmns_up -> bsupvmns
! hbumnc_dw -> bsubumnc
! hbumns_dw -> bsubumns
! hbvmnc_dw -> bsubvmnc
! hbvmns_dw -> bsubvmns
! hbsmnc_dw -> bsubsmnc
! hbsmns_dw -> bsubsmns

! For 3D array indexing, (m,n,j) -> (n,m,j)



!!$      module mod_nemec
!!$!-----------------------------------------------------------------------------!
!!$!  NEMEC parameters and fields                                                !
!!$!-----------------------------------------------------------------------------!
!!$      implicit none
!!$
!!$      integer :: nfp,nrho,mpnt,nsin
!!$      integer :: mpol,ntor,mpol1,iasym
!!$
!!$      real(DP)    :: gam,phiedge
!!$
!!$      real(DP), allocatable :: hiota(:),hpres(:),hbuco(:),hbvco(:)      
!!$      real(DP), allocatable :: hmass(:),hphip(:),hphi(:),hvp(:)
!!$      real(DP), allocatable :: hoverr(:),fjcuru(:),fjcurv(:),hspecw(:)
!!$
!!$      real(DP), allocatable :: frmnc(:,:,:),frmns(:,:,:)
!!$      real(DP), allocatable :: fzmnc(:,:,:),fzmns(:,:,:)
!!$      real(DP), allocatable :: hbsmnc_dw(:,:,:),hbsmns_dw(:,:,:)
!!$      real(DP), allocatable :: hbumnc_dw(:,:,:),hbumns_dw(:,:,:)
!!$      real(DP), allocatable :: hbvmnc_dw(:,:,:),hbvmns_dw(:,:,:)
!!$      real(DP), allocatable :: hbumnc_up(:,:,:),hbumns_up(:,:,:)
!!$      real(DP), allocatable :: hbvmnc_up(:,:,:),hbvmns_up(:,:,:)
!!$      real(DP), allocatable :: flmnc(:,:,:),flmns(:,:,:)
!!$      
!!$      real(DP),allocatable :: fsve(:),hsve(:)    ! radial mesh
!!$
!!$      end module mod_nemec

!-----------------------------------------------------------------------------!

      subroutine read_nemec_file(vmec,in_equilibrium,format_type,ok)
!-----------------------------------------------------------------------------!
! purpose: reads NEMEC output file wout....                                   !
!                                                                             !
! NOTE:                                                                       !
! The values of hbumnc_up(0:mpol1,-ntor:ntor,0)                               !
!                    :                                                        !
!               hbvmns_dw(0:mpol1,-ntor:ntor,0)                               !
! are dummies. If they are needed, than the corresponding radial functions    !
!             hbumnc_up(:,:,1:nsin), ..., hbvmns_dw(:,:,1:nsin)               !
! have to be extrapolated to the magnetic axis.                               !
!-----------------------------------------------------------------------------!
      !use mod_nemec

      use kind_defs
      use constants
      use pisa_vmec_module

      implicit none
      
! --- input parameters
      type(vmec_eqdata) :: vmec
      character(200), intent(in) :: in_equilibrium
      character*25,  intent(in) :: format_type
      integer, intent(out)      :: ok
      
! --- local parameters
      integer,parameter :: inre3=30, outp6=6
      !integer,parameter :: inre3=8, outp6=6
      integer           :: ierr,itype,j,m,n,nmin0
      real(DP)          :: enfp,enrho,empol,entor,empnt,eiasym
      real(DP)          :: ds
      logical           :: ex 


      integer :: mpol1, mpnt, nsin, nn
      real(DP), allocatable :: hpres(:),hbuco(:),hbvco(:)      
      real(DP), allocatable :: hmass(:),hphip(:),hphi(:),hvp(:)
      real(DP), allocatable :: hoverr(:),fjcuru(:),fjcurv(:),hspecw(:)
      real(DP) :: phiedge
      real(DP), dimension(:,:,:), allocatable :: hbsubsmnc, hbsubsmns

      ok = 0

      ! test input
      if(format_type == 'unformatted') then
        itype=0
      elseif (format_type == 'formatted') then
        itype=1
      else
        write(outp6,1) trim(format_type)
    1   format('********** USER error **********',/,      &
               'format_type:        ',1x,a,/,             &
               '****** no valid key-word: stop *****',/)
        ok=-1
        return
      endif
      
      ! open equilibrium file
      open(unit=inre3,file=trim(in_equilibrium),form=trim(format_type), &
           status="old",iostat=ierr)
      if(ierr.ne.0) then
        write(outp6,2) ierr,trim(in_equilibrium)
    2   format('********** USER error **********',/, &
               'ierr = ',i3,/,                       &
               'could not open file:',/,             &
               3x,A120,/,                            &
               'STOP')
        ok = -1
        return
      endif

! --- read dimensions
      if(itype.eq.0) then
        read(inre3) vmec%gam,enfp,enrho,empol,entor,empnt,eiasym,phiedge
      else
        read(inre3,*) vmec%gam,enfp,enrho,empol,entor,empnt,eiasym,phiedge
      endif
      vmec%nfp    = nint(enfp)
      !nrho   = nint(enrho)
      vmec%ns   = nint(enrho)
      vmec%mpol   = nint(empol)
      vmec%ntor   = nint(entor)
      mpnt   = nint(empnt)
      vmec%iasym  = nint(eiasym)
 
      !nsin   = nrho-1
      nsin   = vmec%ns-1
      mpol1  = vmec%mpol-1
      ds     = 1./float(nsin)
      
      allocate(vmec%rmnc(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%rmns(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%zmnc(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%zmns(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      
      allocate(vmec%bsupumnc(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%bsupumns(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%bsupvmnc(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%bsupvmns(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      
      allocate(hbsubsmnc(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(hbsubsmns(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%bsubsmnc(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%bsubsmns(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%bsubumnc(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%bsubumns(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%bsubvmnc(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%bsubvmns(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      
      allocate(vmec%lmnc(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))
      allocate(vmec%lmns(-vmec%ntor:vmec%ntor,0:mpol1,vmec%ns))

      allocate(vmec%iotas(vmec%ns),hpres(vmec%ns),hbuco(vmec%ns),hbvco(vmec%ns))
      allocate(hmass(vmec%ns),hphip(vmec%ns),vmec%phi(vmec%ns),hvp(vmec%ns))
      allocate(hoverr(vmec%ns),fjcuru(vmec%ns),fjcurv(vmec%ns))         
      allocate(hspecw(vmec%ns))
      allocate(vmec%phips(vmec%ns))

      !allocate(fsve(0:nsin),hsve(1:nsin))

      vmec%rmnc=0
      vmec%rmns=0
      vmec%zmnc=0
      vmec%zmns=0
      vmec%lmnc=0
      vmec%lmns=0
      vmec%bsupumnc=0
      vmec%bsupumns=0
      vmec%bsupvmnc=0
      vmec%bsupvmns=0
      vmec%bsubumnc=0
      vmec%bsubumns=0
      vmec%bsubvmnc=0
      vmec%bsubvmns=0
      vmec%bsubsmnc=0
      vmec%bsubsmns=0
      hbsubsmnc=0
      hbsubsmns=0
      vmec%iotas=0

! --- read VMEC output
      if(itype.eq.0) then
       
        !do j=0,nsin
        do j=1,vmec%ns
          do m=0,mpol1
            nmin0=-vmec%ntor
            if(m.eq.0) nmin0=0
            do n=nmin0,vmec%ntor
               ! Sign flips for consistency with Hakan's conventions:
               nn=-n
! --- full mesh
              read(inre3) vmec%rmnc(nn,m,j),vmec%zmns(nn,m,j),         &  
                          vmec%rmns(nn,m,j),vmec%zmnc(nn,m,j),         &  
! --- half mesh
                          vmec%bsupumnc(nn,m,j),vmec%bsupvmnc(nn,m,j), &  
                          vmec%bsupumns(nn,m,j),vmec%bsupvmns(nn,m,j), &  
! --- full mesh
                          vmec%lmns(nn,m,j),vmec%lmnc(nn,m,j),         &  
! --- half mesh
                          vmec%bsubumnc(nn,m,j),vmec%bsubvmnc(nn,m,j), &  
                          hbsubsmns(nn,m,j),                  &  
                          vmec%bsubumns(nn,m,j),vmec%bsubvmns(nn,m,j), &  
                          hbsubsmnc(nn,m,j)                  
            enddo
          enddo
        enddo

! --- half mesh
        read(inre3) (vmec%iotas(j),hmass(j),hpres(j),hphip(j),hbuco(j), &  
                    hbvco(j),vmec%phi(j),hvp(j),hoverr(j),fjcuru(j),   &
                    fjcurv(j),hspecw(j),j=2,vmec%ns)
                    !fjcurv(j),hspecw(j),j=1,nsin)
      else
       
        !do j=0,nsin
        do j=1,vmec%ns
          do m=0,mpol1
            nmin0=-vmec%ntor
            if(m.eq.0) nmin0=0
            do n=nmin0,vmec%ntor
               ! Sign flips for consistency with Hakan's conventions:
               nn=-n
! --- full mesh
              read(inre3,*) vmec%rmnc(nn,m,j),vmec%zmns(nn,m,j),       &  
                          vmec%rmns(nn,m,j),vmec%zmnc(nn,m,j),         &  
! --- half mesh
                          vmec%bsupumnc(nn,m,j),vmec%bsupvmnc(nn,m,j), &  
                          vmec%bsupumns(nn,m,j),vmec%bsupvmns(nn,m,j), &  
! --- full mesh
                          vmec%lmns(nn,m,j),vmec%lmnc(nn,m,j),         &  
! --- half mesh
                          vmec%bsubumnc(nn,m,j),vmec%bsubvmnc(nn,m,j), &  
                          hbsubsmns(nn,m,j),                  &  
                          vmec%bsubumns(nn,m,j),vmec%bsubvmns(nn,m,j), &  
                          hbsubsmnc(nn,m,j)                  
            enddo
          enddo
        enddo

! --- half mesh
        read(inre3,*) (vmec%iotas(j),hmass(j),hpres(j),hphip(j),hbuco(j),   &  
                      hbvco(j),vmec%phi(j),hvp(j),hoverr(j),fjcuru(j),     &
                      fjcurv(j),hspecw(j),j=2,vmec%ns)
!                      fjcurv(j),hspecw(j),j=1,nsin)
      endif
      
      close(inre3)

      vmec%phips = -phiedge/twopi
      vmec%phips(1) = 0

      !print *,"phiedge:",phiedge
      !print *,"hphip:",hphip
      !print *,"vmec%phi:",vmec%phi
      print *,"Data from NEMEC:"
      print *,"  ns:",vmec%ns
      print *,"  ntor:",vmec%ntor
      print *,"  mpol:",vmec%mpol
      print *,"  nfp:",vmec%nfp


      ! Sign flips for consistency with Hakan's conventions:
      vmec%phips = -vmec%phips
      vmec%phi = -vmec%phi
      vmec%iotas = -vmec%iotas
      vmec%bsupvmnc = -vmec%bsupvmnc
      vmec%bsupvmns = -vmec%bsupvmns
      vmec%bsubvmnc = -vmec%bsubvmnc
      vmec%bsubvmns = -vmec%bsubvmns
      

!!$      print *,"rmnc at outermost radius: ",vmec%rmnc(:,:,vmec%ns)
!!$      print *,"rmns at outermost radius: ",vmec%rmns(:,:,vmec%ns)
!!$      print *,"zmnc at outermost radius: ",vmec%zmnc(:,:,vmec%ns)
!!$      print *,"zmns at outermost radius: ",vmec%zmns(:,:,vmec%ns)

!!$! --- grid in s
!!$      fsve(0)     = 0.
!!$      fsve(nsin)  = 1.
!!$      do j=1,nsin-1
!!$! --- half mesh
!!$       hsve(j)  = (j-0.5)*ds
!!$! --- full mesh
!!$        fsve(j)   = j*ds
!!$      enddo
!!$      hsve(nsin) = (nsin-0.5)*ds

!!$      !test output
!!$      write(6,5)     
!!$    5 format(6x,'hsve',10x'hiota',9x,'hmass',9x,'hpres',9x,'hphip', &
!!$             9x,'hbuco',9x,'hbvco',10x,'hphi',11x,'hvp',9x,'hoverr',&
!!$             8x,'hspecw')
!!$      do j=1,nsin
!!$        write(6,3) hsve(j),hiota(j),hmass(j),hpres(j),hphip(j), &
!!$                   hbuco(j),hbvco(j),hphi(j),hvp(j),hoverr(j),hspecw(j)
!!$      enddo
!!$    3 format(13(2x,e12.4))
!!$      write(6,6)
!!$    6 format(6x,'fsve',10x,'fjcuru',8x,'fjcurv')
!!$      do j=1,nsin
!!$        write(6,4) fsve(j),fjcuru(j),fjcurv(j)
!!$      enddo
!!$    4 format(3(2x,e12.4))


      ! Modern VMEC saves bsubs on the full mesh, so this is what
      ! SFINCS uses, whereas
      ! Strumberger's format has bsubs on the half mesh.  Conversion is done here:
      vmec%bsubsmnc(:,:,2:vmec%ns-1) = (hbsubsmnc(:,:,2:vmec%ns-1) + hbsubsmnc(:,:,3:vmec%ns)) / 2
      vmec%bsubsmns(:,:,2:vmec%ns-1) = (hbsubsmns(:,:,2:vmec%ns-1) + hbsubsmns(:,:,3:vmec%ns)) / 2
      vmec%bsubsmnc(:,:,1) = hbsubsmnc(:,:,2)*1.5_DP - hbsubsmnc(:,:,3)*0.5_DP
      vmec%bsubsmns(:,:,1) = hbsubsmns(:,:,2)*1.5_DP - hbsubsmns(:,:,3)*0.5_DP
      vmec%bsubsmnc(:,:,vmec%ns) = hbsubsmnc(:,:,vmec%ns)*1.5_DP - hbsubsmnc(:,:,vmec%ns-1)*0.5_DP
      vmec%bsubsmns(:,:,vmec%ns) = hbsubsmns(:,:,vmec%ns)*1.5_DP - hbsubsmns(:,:,vmec%ns-1)*0.5_DP
      
    end subroutine read_nemec_file

!-----------------------------------------------------------------------------!
      
!!$      program readnemec
!!$!-----------------------------------------------------------------------------!
!!$! purpose: call of read_nemec.f90                                             !
!!$!-----------------------------------------------------------------------------!
!!$      implicit none
!!$      
!!$      character*250 :: in_equilibrium    ! filename of NEMEC output 
!!$      character*25  :: format_type       ! file format (formatted/unformatted)
!!$      integer       :: ok
!!$      
!!$      in_equilibrium = 'wout.test'
!!$      format_type    = 'formatted'
!!$      
!!$      call read_nemec(in_equilibrium,format_type,ok)
!!$      
!!$      end program readnemec
