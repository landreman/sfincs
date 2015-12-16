subroutine NEMEC_compute_missing_fields(vmec)

  ! This subroutine computes |B| and the Jacobian.

  use kind_defs
  use constants
  use pisa_vmec_module

  implicit none

  type(vmec_eqdata) :: vmec

  integer :: m, n, js, nmin0, nu, nv, j, iu, iv
  real(DP), dimension(:), allocatable :: u, v
  real(DP), dimension(:,:), allocatable :: hBx, hBy, hBz, hB, hg
  real(DP), dimension(:,:,:), allocatable :: fx, fy, fz
  real(DP), dimension(:,:,:), allocatable :: fdxdu, fdydu, fdzdu
  real(DP), dimension(:,:,:), allocatable :: fdxdv, fdydv, fdzdv
  real(DP), dimension(:,:), allocatable :: hdxdu, hdydu, hdzdu
  real(DP), dimension(:,:), allocatable :: hdxdv, hdydv, hdzdv
  real(DP), dimension(:,:), allocatable :: hdxds, hdyds, hdzds
  real(DP), dimension(:,:,:), allocatable :: hbsupu, hbsupv
  real(DP) :: ds, angle, cosangle, sinangle, cosv, sinv
  real(DP) :: bctemp, bstemp, gctemp, gstemp, dnorm
  integer :: tic, toc, countrate

  print *,"Computing B and g from NEMEC data..."
  call system_clock(tic,countrate)

  ! I'm not certain these values for nu and nv in the next 2 lines are the best choice,
  ! but they should be roughly ok.
  nu = 2*vmec%mpol
  nv = 2*(2*vmec%ntor+1)
  allocate(u(nu))
  allocate(v(nv))
  do j=1,nu
     u(j) = twopi*(j-1.0_DP)/nu
  end do
  do j=1,nv
     v(j) = twopi*(j-1.0_DP)/(vmec%nfp * nv)
  end do
  ds = 1.0_DP / (vmec%ns-1)
  !print *,"u: ",u
  !print *,"v: ",v

  if (.not. associated(vmec%bmnc)) allocate(vmec%bmnc(-vmec%ntor:vmec%ntor,0:vmec%mpol-1,vmec%ns))
  if (.not. associated(vmec%bmns)) allocate(vmec%bmns(-vmec%ntor:vmec%ntor,0:vmec%mpol-1,vmec%ns))
  if (.not. associated(vmec%gmnc)) allocate(vmec%gmnc(-vmec%ntor:vmec%ntor,0:vmec%mpol-1,vmec%ns))
  if (.not. associated(vmec%gmns)) allocate(vmec%gmns(-vmec%ntor:vmec%ntor,0:vmec%mpol-1,vmec%ns))

  vmec%bmnc=0
  vmec%bmns=0
  vmec%gmnc=0
  vmec%gmns=0

  allocate(fx(nu,nv,vmec%ns))
  allocate(fy(nu,nv,vmec%ns))
  allocate(fz(nu,nv,vmec%ns))
  allocate(fdxdu(nu,nv,vmec%ns))
  allocate(fdydu(nu,nv,vmec%ns))
  allocate(fdzdu(nu,nv,vmec%ns))
  allocate(fdxdv(nu,nv,vmec%ns))
  allocate(fdydv(nu,nv,vmec%ns))
  allocate(fdzdv(nu,nv,vmec%ns))

  allocate(hbsupu(nu,nv,vmec%ns))
  allocate(hbsupv(nu,nv,vmec%ns))

  allocate(hdxdu(nu,nv))
  allocate(hdydu(nu,nv))
  allocate(hdzdu(nu,nv))
  allocate(hdxdv(nu,nv))
  allocate(hdydv(nu,nv))
  allocate(hdzdv(nu,nv))
  allocate(hdxds(nu,nv))
  allocate(hdyds(nu,nv))
  allocate(hdzds(nu,nv))

  allocate(hBx(nu,nv))
  allocate(hBy(nu,nv))
  allocate(hBz(nu,nv))
  allocate(hB(nu,nv))
  allocate(hg(nu,nv))

  fx=0
  fy=0
  fz=0
  fdxdu=0
  fdydu=0
  fdzdu=0
  fdxdv=0
  fdydv=0
  fdzdv=0

  hdxdu=0
  hdydu=0
  hdzdu=0
  hdxdv=0
  hdydv=0
  hdzdv=0
  hdxds=0
  hdyds=0
  hdzds=0

  hbsupu=0
  hbsupv=0

  hBx=0
  hBy=0
  hBz=0
  hB=0
  hg=0

  ! In case a stellarator-symmetric equilibrium was loaded, allocate the non-symmetric arrays and set them to 0:
  if (vmec%iasym .ne. 1) then
     allocate(vmec%rmns(-vmec%ntor:vmec%ntor,0:vmec%mpol-1,vmec%ns))
     allocate(vmec%zmnc(-vmec%ntor:vmec%ntor,0:vmec%mpol-1,vmec%ns))
     allocate(vmec%bsupumns(-vmec%ntor:vmec%ntor,0:vmec%mpol-1,vmec%ns))
     allocate(vmec%bsupvmns(-vmec%ntor:vmec%ntor,0:vmec%mpol-1,vmec%ns))
     vmec%rmns=0
     vmec%zmnc=0
     vmec%bsupumns=0
     vmec%bsupvmns=0
  end if

  ! rmnc, rmns, zmnc, and zmns are on the FULL grid.
  ! bsupu and bsupv are on the HALF grid.

  ! We want bmnc, bmns, gmnc, and gmns on the HALF grid.

  ! x = [rmnc*cos(m*u-n*v*nfp)+rmns*sin(m*u-n*v*nfp)]*cos(v)
  ! y = [rmnc*cos(m*u-n*v*nfp)+rmns*sin(m*u-n*v*nfp)]*sin(v)
  ! z = [zmnc*cos(m*u-n*v*nfp)+zmns*sin(m*u-n*v*nfp)]
  do js = 1,vmec%ns
     do iv = 1,nv
        cosv = cos(v(iv))
        sinv = sin(v(iv))
        do iu = 1,nu
           do m=0,vmec%mpol-1
              nmin0=-vmec%ntor
              if(m.eq.0) nmin0=0
              do n=nmin0,vmec%ntor
                 angle = m*u(iu)-n*vmec%nfp*v(iv)
                 cosangle = cos(angle)
                 sinangle = sin(angle)

                 fx(iu,iv,js) = fx(iu,iv,js) + (vmec%rmnc(n,m,js)*cosangle + vmec%rmns(n,m,js)*sinangle)*cosv
                 fdxdu(iu,iv,js) = fdxdu(iu,iv,js) + (-m*vmec%rmnc(n,m,js)*sinangle + m*vmec%rmns(n,m,js)*cosangle)*cosv
                 fdxdv(iu,iv,js) = fdxdv(iu,iv,js) + (n*vmec%nfp*vmec%rmnc(n,m,js)*sinangle - n*vmec%nfp*vmec%rmns(n,m,js)*cosangle)*cosv &
                      - (vmec%rmnc(n,m,js)*cosangle + vmec%rmns(n,m,js)*sinangle)*sinv

                 fy(iu,iv,js) = fy(iu,iv,js) + (vmec%rmnc(n,m,js)*cosangle + vmec%rmns(n,m,js)*sinangle)*sinv
                 fdydu(iu,iv,js) = fdydu(iu,iv,js) + (-m*vmec%rmnc(n,m,js)*sinangle + m*vmec%rmns(n,m,js)*cosangle)*sinv
                 fdydv(iu,iv,js) = fdydv(iu,iv,js) + (n*vmec%nfp*vmec%rmnc(n,m,js)*sinangle - n*vmec%nfp*vmec%rmns(n,m,js)*cosangle)*sinv &
                      + (vmec%rmnc(n,m,js)*cosangle + vmec%rmns(n,m,js)*sinangle)*cosv

                 fz(iu,iv,js) = fz(iu,iv,js) + vmec%zmnc(n,m,js)*cosangle + vmec%zmns(n,m,js)*sinangle
                 fdzdu(iu,iv,js) = fdzdu(iu,iv,js) - m*vmec%zmnc(n,m,js)*sinangle + m*vmec%zmns(n,m,js)*cosangle
                 fdzdv(iu,iv,js) = fdzdv(iu,iv,js) + n*vmec%nfp*vmec%zmnc(n,m,js)*sinangle - n*vmec%nfp*vmec%zmns(n,m,js)*cosangle

                 hbsupu(iu,iv,js) = hbsupu(iu,iv,js) + vmec%bsupumnc(n,m,js)*cosangle + vmec%bsupumns(n,m,js)*sinangle
                 hbsupv(iu,iv,js) = hbsupv(iu,iv,js) + vmec%bsupvmnc(n,m,js)*cosangle + vmec%bsupvmns(n,m,js)*sinangle
              end do
           end do
        end do
     end do
  end do
  call system_clock(toc)
  print *,"  Fourier -> grid: ",real(toc-tic)/countrate,"sec."
  call system_clock(tic)

  do js = 1,vmec%ns-1
     ! Evaluate quantities we need on the half mesh:
     hdxdu = 0.5*(fdxdu(:,:,js) + fdxdu(:,:,js+1))
     hdxdv = 0.5*(fdxdv(:,:,js) + fdxdv(:,:,js+1))
     hdydu = 0.5*(fdydu(:,:,js) + fdydu(:,:,js+1))
     hdydv = 0.5*(fdydv(:,:,js) + fdydv(:,:,js+1))
     hdzdu = 0.5*(fdzdu(:,:,js) + fdzdu(:,:,js+1))
     hdzdv = 0.5*(fdzdv(:,:,js) + fdzdv(:,:,js+1))

     hdxds = (fx(:,:,js+1) - fx(:,:,js)) / ds
     hdyds = (fy(:,:,js+1) - fy(:,:,js)) / ds
     hdzds = (fz(:,:,js+1) - fz(:,:,js)) / ds

     ! Compute Jacobian:
     hg = hdxds*hdydu*hdzdv + hdyds*hdzdu*hdxdv + hdzds*hdxdu*hdydv - hdzds*hdydu*hdxdv - hdxds*hdzdu*hdydv - hdyds*hdxdu*hdzdv

     ! Compute |B|:
     hBx = hbsupu(:,:,js+1)*hdxdu + hbsupv(:,:,js+1)*hdxdv
     hBy = hbsupu(:,:,js+1)*hdydu + hbsupv(:,:,js+1)*hdydv
     hBz = hbsupu(:,:,js+1)*hdzdu + hbsupv(:,:,js+1)*hdzdv
     hB = sqrt(hBx*hBx + hBy*hBy + hBz*hBz)

     !print *,"Here comes |B| at js=",js
     !print *,hB

     ! Map g & B from the grid to Fourier amplitudes:
     do m = 0, vmec%mpol-1
        nloop: do n=-vmec%ntor,vmec%ntor
           if (m.eq.0 .and. n.gt.0) cycle nloop
           dnorm = (1.0_dp)/(nu*nv)
           if (m.ne.0 .or. n.ne.0) dnorm = 2*dnorm
           bctemp = 0
           bstemp = 0
           gctemp = 0
           gstemp = 0
           do iu = 1, nu
              do iv = 1, nv
                 cosangle = cos(m*u(iu)-n*vmec%nfp*v(iv))*dnorm
                 sinangle = sin(m*u(iu)-n*vmec%nfp*v(iv))*dnorm
                 bctemp  = bctemp  + hB(iu,iv) * cosangle
                 bstemp  = bstemp  + hB(iu,iv) * sinangle
                 gctemp  = gctemp  + hg(iu,iv) * cosangle
                 gstemp  = gstemp  + hg(iu,iv) * sinangle
              end do
           end do
           vmec%bmnc(n,m,js+1) = bctemp
           vmec%bmns(n,m,js+1) = bstemp
           vmec%gmnc(n,m,js+1) = gctemp
           vmec%gmns(n,m,js+1) = gstemp
        end do nloop
     end do
  end do
  call system_clock(toc)
  print *,"  grid -> Fourier: ",real(toc-tic)/countrate,"sec."


  deallocate(fx,fy,fz)
  deallocate(fdxdu,fdxdv,fdydu,fdydv,fdzdu,fdzdv)
  deallocate(hdxds,hdxdu,hdxdv)
  deallocate(hdyds,hdydu,hdydv)
  deallocate(hdzds,hdzdu,hdzdv)
  deallocate(hg,hBx,hBy,hBz,hB)
  deallocate(hbsupu,hbsupv)

end subroutine NEMEC_compute_missing_fields
