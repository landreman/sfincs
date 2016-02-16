!> $Id: pisa_vmec_module.f90 6192 2010-08-11 14:55:36Z geiger $
!
module pisa_vmec_module
use kind_defs
implicit none

!> Identifier which the vmec-output format has been detected.
type file_type
  character(len=3) :: typo = '   '
end type file_type

type vmec_minimum_eqdata
 integer :: ns    !< number of radial surfaces
 integer :: mpol  !< number of poloidal fourier modes
 integer :: ntor  !< upper bound toroidal fourier modes: -ntor <= n <= ntor
 integer :: nfp   !< number of field periods
 integer :: iasym !<  defines stellarator symmetry for iasym=0, otherwise =1
 real(DP) :: phiedge !<  total toroidal flux within the boundary (s=1)
 real(DP), dimension(:), pointer :: iota => NULL()     !< iota-profile
 real(DP), dimension(:,:,:), pointer :: rmnc => NULL() !< Fourier coefficients for R with cosine
 real(DP), dimension(:,:,:), pointer :: rmns => NULL() !< Fourier coefficients for R with sine
 real(DP), dimension(:,:,:), pointer :: zmns => NULL() !< Fourier coefficients for z with sine
 real(DP), dimension(:,:,:), pointer :: zmnc => NULL() !< Fourier coefficients for z with cosine
end type vmec_minimum_eqdata

type vmec_wout_data
! data object to hold the information of a vmec_wout-file. The implementation
! ends before reading fsqt and wdot. Implementation of these may come later.
 character(len=10) :: version   !! holds information about what vmec wout-version is stored.
 real(DP) :: wb           !<  magnetic energy
 real(DP) :: wp           !<  kinetic energy
 real(DP) :: gam          !<  ratio of specific heats used in vmec for energy principle
 real(DP) :: pfac         !<  in reconstruction: ACTUAL PRESSURE = PRESPEAK * PFAC * P(INTERNAL)
 real(DP) :: rmax_surf    !<  maximum value of R of flux surfaces
 real(DP) :: rmin_surf    !<  minimum value of R of flux surfaces
 real(DP) :: zmax_surf    !<  maximum value of z of flux surfaces
 real(DP) :: aspect       !<  aspect ratio
 real(DP) :: betatot      !<  total beta
 real(DP) :: betapol      !<  poloidal beta
 real(DP) :: betator      !<  toroidal beta
 real(DP) :: betaxis      !<  beta-value on axis
 real(DP) :: b0           !<  Btor on axis (seems to be similar to rbtor0)
 real(DP) :: IonLarmor    !<  Ion-Larmor-radius
 real(DP) :: VolAvgB      !<  volume averaged B
 real(DP) :: rbtor0       !<  mean-Btor on axis
 real(DP) :: rbtor        !<  mean-Btor at boundary (s=1)
 real(DP) :: ctor         !<  toroidal current
 real(DP) :: Aminor_p     !<  minor plasma radius
 real(DP) :: Rmajor_p     !<  major plasma radius
 real(DP) :: volume_p     !<  plasma volume
 integer :: nfp           !<  number of field periods
 integer :: ns            !<  number of radial grid points (flux surfaces)
 integer :: mpol          !<  number of poloidal harmonics: m = 0,...,mpol-1
 integer :: ntor          !<  number of toroidal harmonics: n = -ntor, ntor
 integer :: mnmax         !<  number of fourier coefficients per flux surface\n
                          !!    mnmax = ntor + (mpol-1)*(1+2*ntor)
 integer :: mnmax_nyq0    !<  number of fourier coefficients per flux surface\n
                          !!    mnmax = ntor + (mpol-1)*(1+2*ntor)
 integer :: itfsq         !<  number of fsq-values for plotting
 integer :: niter         !<  maximum number of main iterations in vmec
 integer :: iasym         !<  defines whether stellarator symmetry was used:\n
                          !!    iasym = 0 <==> stellarator symmetry\n
                          !!          = 1 <==> no symmetry
 integer :: ireconstruct  !<  reconstruction flag 0=no, 1=yes
 integer :: ier_flag      !<  error flag (vmec_params.f)\n
                          !!  0    o.k.\n
                          !!  1    bad initial jacobian\n
                          !!  4    more iterations to reach force limit\n
                          !!  6    bad jacobian\n
 integer :: imse2_minus_1 !<  mse-information for reconstruction mode
 integer :: itse          !<  thomson data points used in reconstruction mode
 integer :: nbsets        !<  number of B-coil sets defined in mgrid file   (readin.f)
 integer :: nobd          !<  number of connected flux loop measurements
 integer :: nextcur       !<  number of external currents
 integer :: nstore_seq    !<  number of stored sequences
 integer :: signgs        !<  sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)
 character(len=200) :: mgrid_file       !<  name of mgrid_file used in calculation
 character(len=132) :: input_extension  !<  suffix of the vmec-input file: input.<input_extension>
 character(len=30), dimension(:), pointer :: curlabel  !<  strings with names for the external coil currents
 real(DP), dimension(:), pointer :: extcur   =>NULL()  !<  coil currents used.
 real(DP), dimension(:), pointer :: iotaf    =>NULL()  !<  iota         : full mesh
 real(DP), dimension(:), pointer :: presf    =>NULL()  !<  pressure     : full mesh
 real(DP), dimension(:), pointer :: phipf    =>NULL()  !<  dPhi/ds      : full mesh
 real(DP), dimension(:), pointer :: phi      =>NULL()  !<  Phi          : full mesh
 real(DP), dimension(:), pointer :: jcuru    =>NULL()  !<  jpol(dens)   : full mesh
 real(DP), dimension(:), pointer :: jcurv    =>NULL()  !<  jtor(dens)   : full mesh
 real(DP), dimension(:), pointer :: iotas    =>NULL()  !<  iota         : half mesh
 real(DP), dimension(:), pointer :: mass     =>NULL()  !<  mass         : half mesh
 real(DP), dimension(:), pointer :: pres     =>NULL()  !<  pressure     : half mesh
 real(DP), dimension(:), pointer :: beta_vol =>NULL()  !<  beta_vol     : half mesh, beta-profile but not vol. av.!
 real(DP), dimension(:), pointer :: phip     =>NULL()  !<  2*PI*dPhi/ds : half mesh, in vmec2000 called phips
 real(DP), dimension(:), pointer :: buco     =>NULL()  !<  pol. Curr.   : half mesh
 real(DP), dimension(:), pointer :: bvco     =>NULL()  !<  tor. Curr.   : half mesh
 real(DP), dimension(:), pointer :: vp       =>NULL()  !<  dV/ds        : half mesh
 real(DP), dimension(:), pointer :: overr    =>NULL()  !<  (1/V')int (1/R) dV : half mesh
 real(DP), dimension(:), pointer :: specw    =>NULL()  !<  spectr. width : half mesh
 real(DP), dimension(:), pointer :: Dmerc    =>NULL()  !<  Mercier-criterion        : full mesh
 real(DP), dimension(:), pointer :: Dshear   =>NULL()  !<  Shear part of Mercier    : full mesh
 real(DP), dimension(:), pointer :: Dwell    =>NULL()  !<  Well part of Mercier     : full mesh
 real(DP), dimension(:), pointer :: Dcurr    =>NULL()  !<  Current part of Mercier  : full mesh
 real(DP), dimension(:), pointer :: Dgeod    =>NULL()  !<  geodesic part of Mercier : full mesh
 real(DP), dimension(:), pointer :: equif    =>NULL()  !<  force balance mismatch   : full mesh
 real(DP), dimension(:), pointer :: xm       =>NULL()  !<  poloidal fourier numbers
 real(DP), dimension(:), pointer :: xn       =>NULL()  !<  toroidal fourier numbers (number of field periods can be derived from this!)
 real(DP), dimension(:), pointer :: xm_nyq   =>NULL()  !<  poloidal fourier numbers (Nyquist)
 real(DP), dimension(:), pointer :: xn_nyq   =>NULL()  !<  toroidal fourier numbers (Nyquist)
 real(DP), dimension(:,:,:), pointer :: rmnc =>NULL()    !<  full mesh: Rmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: rmns =>NULL()    !<  full mesh: Rmns sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: zmns =>NULL()    !<  full mesh: zmns sin(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: zmnc =>NULL()    !<  full mesh: zmnc cos(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: lmns =>NULL()    !<  half mesh: lmns sin(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: lmnc =>NULL()    !<  half mesh: lmnc cos(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: bmnc =>NULL()    !<  half mesh: bmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bmns =>NULL()    !<  half mesh: bmns sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: gmnc =>NULL()    !<  half mesh: gmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: gmns =>NULL()    !<  half mesh: gmns sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsubumnc =>NULL()  !<  half mesh: B_u(mn)c cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsubumns =>NULL()  !<  half mesh: B_u(mn)s sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsubvmnc =>NULL()  !<  half mesh: B_v(mn)c cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsubvmns =>NULL()  !<  half mesh: B_v(mn)s sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsubsmnc =>NULL()  !<  half mesh: B_s(mn)c cos(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsubsmns =>NULL()  !<  half mesh: B_s(mn)s sin(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsupumnc =>NULL()  !<  half mesh: B^u(mn)c cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsupumns =>NULL()  !<  half mesh: B^u(mn)s sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsupvmnc =>NULL()  !<  half mesh: B^v(mn)c cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsupvmns =>NULL()  !<  half mesh: B^v(mn)s sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:,:), pointer :: currvmnc =>NULL()  ! full mesh: curv(mn)c cos(mu+nv) stellarator symmetric part
end type vmec_wout_data

type vmec_wout_cdf_data
! data object to hold the information of a vmec_wout-netcdf-file. The implementation
! does not read everything but only the data needed up to now. The rest may come later.
 character(len=10) :: version   !! holds information about what vmec wout-version is stored.
 real(DP) :: wb           !<  magnetic energy
 real(DP) :: wp           !<  kinetic energy
 real(DP) :: gam          !<  ratio of specific heats used in vmec for energy principle
 real(DP) :: pfac         !<  in reconstruction: ACTUAL PRESSURE = PRESPEAK * PFAC * P(INTERNAL)
 real(DP) :: rmax_surf    !<  maximum value of R of flux surfaces
 real(DP) :: rmin_surf    !<  minimum value of R of flux surfaces
 real(DP) :: zmax_surf    !<  maximum value of z of flux surfaces
 real(DP) :: aspect       !<  aspect ratio
 real(DP) :: betatot      !<  total beta
 real(DP) :: betapol      !<  poloidal beta
 real(DP) :: betator      !<  toroidal beta
 real(DP) :: betaxis      !<  beta-value on axis
 real(DP) :: b0           !<  Btor on axis (seems to be similar to rbtor0)
 real(DP) :: IonLarmor    !<  Ion-Larmor-radius
 real(DP) :: VolAvgB      !<  volume averaged B
 real(DP) :: rbtor0       !<  mean-Btor on axis
 real(DP) :: rbtor        !<  mean-Btor at boundary (s=1)
 real(DP) :: ctor         !<  toroidal current
 real(DP) :: Aminor_p     !<  minor plasma radius
 real(DP) :: Rmajor_p     !<  major plasma radius
 real(DP) :: volume_p     !<  plasma volume
 integer :: nfp           !<  number of field periods
 integer :: ns            !<  number of radial grid points (flux surfaces)
 integer :: mpol          !<  number of poloidal harmonics: m = 0,...,mpol-1
 integer :: ntor          !<  number of toroidal harmonics: n = -ntor, ntor
 integer :: mnmax         !<  number of fourier coefficients per flux surface\n
                          !!    mnmax = ntor + (mpol-1)*(1+2*ntor)
 integer :: mnmax_nyq0    !<  number of fourier coefficients per flux surface\n
                          !!    mnmax = ntor + (mpol-1)*(1+2*ntor)
 integer :: itfsq         !<  number of fsq-values for plotting
 integer :: niter         !<  maximum number of main iterations in vmec
 integer :: iasym         !<  defines whether stellarator symmetry was used:\n
                          !!    iasym = 0 <==> stellarator symmetry\n
                          !!          = 1 <==> no symmetry
 integer :: ireconstruct  !<  reconstruction flag 0=no, 1=yes
 integer :: ifreeb        !<  free boundary flag 0=no, 1=yes
 integer :: ier_flag      !<  error flag (vmec_params.f)\n
                          !!  0    o.k.\n
                          !!  1    bad initial jacobian\n
                          !!  4    more iterations to reach force limit\n
                          !!  6    bad jacobian\n
 integer :: imse2_minus_1 !<  mse-information for reconstruction mode
 integer :: itse          !<  thomson data points used in reconstruction mode
 integer :: nbsets        !<  number of B-coil sets defined in mgrid file   (readin.f)
 integer :: nobd          !<  number of connected flux loop measurements
 integer :: nextcur       !<  number of external currents
 integer :: nstore_seq    !<  number of stored sequences
 integer :: signgs        !<  sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)
 character(len=200) :: mgrid_file       !<  name of mgrid_file used in calculation
 character(len=132) :: input_extension  !<  suffix of the vmec-input file: input.<input_extension>
 character(len=30), dimension(:), pointer :: curlabel  !<  strings with names for the external coil currents
 real(DP), dimension(:), pointer :: extcur   =>NULL()  !<  coil currents used.
 real(DP), dimension(:), pointer :: iotaf    =>NULL()  !<  iota         : full mesh
 real(DP), dimension(:), pointer :: presf    =>NULL()  !<  pressure     : full mesh
 real(DP), dimension(:), pointer :: phipf    =>NULL()  !<  2*PI*dPhi/ds : full mesh
 real(DP), dimension(:), pointer :: phi      =>NULL()  !<  Phi          : full mesh
 real(DP), dimension(:), pointer :: jcuru    =>NULL()  !<  jpol(dens)   : full mesh
 real(DP), dimension(:), pointer :: jcurv    =>NULL()  !<  jtor(dens)   : full mesh
 real(DP), dimension(:), pointer :: iotas    =>NULL()  !<  iota         : half mesh
 real(DP), dimension(:), pointer :: mass     =>NULL()  !<  mass         : half mesh
 real(DP), dimension(:), pointer :: pres     =>NULL()  !<  pressure     : half mesh
 real(DP), dimension(:), pointer :: beta_vol =>NULL()  !<  beta_vol     : half mesh, beta-profile but not vol. av.!
 real(DP), dimension(:), pointer :: phips    =>NULL()  !<  dPhi/ds      : half mesh, in vmec2000 called phips
 real(DP), dimension(:), pointer :: buco     =>NULL()  !<  pol. Curr.   : half mesh
 real(DP), dimension(:), pointer :: bvco     =>NULL()  !<  tor. Curr.   : half mesh
 real(DP), dimension(:), pointer :: vp       =>NULL()  !<  dV/ds        : half mesh
 real(DP), dimension(:), pointer :: overr    =>NULL()  !<  (1/V')int (1/R) dV : half mesh
 real(DP), dimension(:), pointer :: specw    =>NULL()  !<  spectr. width : half mesh
 real(DP), dimension(:), pointer :: Dmerc    =>NULL()  !<  Mercier-criterion        : full mesh
 real(DP), dimension(:), pointer :: Dshear   =>NULL()  !<  Shear part of Mercier    : full mesh
 real(DP), dimension(:), pointer :: Dwell    =>NULL()  !<  Well part of Mercier     : full mesh
 real(DP), dimension(:), pointer :: Dcurr    =>NULL()  !<  Current part of Mercier  : full mesh
 real(DP), dimension(:), pointer :: Dgeod    =>NULL()  !<  geodesic part of Mercier : full mesh
 real(DP), dimension(:), pointer :: equif    =>NULL()  !<  force balance mismatch   : full mesh
 real(DP), dimension(:), pointer :: xm       =>NULL()  !<  poloidal mode numbers
 real(DP), dimension(:), pointer :: xn       =>NULL()  !<  toroidal mode numbers (number of field periods can be derived from this!)
 real(DP), dimension(:), pointer :: xm_nyq   =>NULL()  !<  poloidal mode numbers (Nyquist)
 real(DP), dimension(:), pointer :: xn_nyq   =>NULL()  !<  toroidal mode numbers (Nyquist)
 real(DP), dimension(:,:), pointer :: rmnc =>NULL()    !<  full mesh: Rmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: rmns =>NULL()    !<  full mesh: Rmns sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:), pointer :: zmns =>NULL()    !<  full mesh: zmns sin(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: zmnc =>NULL()    !<  full mesh: zmnc cos(mu+nv)            asymmetric part
 real(DP), dimension(:,:), pointer :: lmns =>NULL()    !<  half mesh: lmns sin(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: lmnc =>NULL()    !<  half mesh: lmnc cos(mu+nv)            asymmetric part
 real(DP), dimension(:,:), pointer :: bmnc =>NULL()    !<  half mesh: bmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: bmns =>NULL()    !<  half mesh: bmns sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:), pointer :: gmnc =>NULL()    !<  half mesh: gmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: gmns =>NULL()    !<  half mesh: gmns sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:), pointer :: bsubumnc =>NULL()  !<  half mesh: B_u(mn)c cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: bsubumns =>NULL()  !<  half mesh: B_u(mn)s sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:), pointer :: bsubvmnc =>NULL()  !<  half mesh: B_v(mn)c cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: bsubvmns =>NULL()  !<  half mesh: B_v(mn)s sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:), pointer :: bsubsmnc =>NULL()  !<  half mesh: B_s(mn)c cos(mu+nv)            asymmetric part
 real(DP), dimension(:,:), pointer :: bsubsmns =>NULL()  !<  half mesh: B_s(mn)s sin(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: bsupumnc =>NULL()  !<  half mesh: B^u(mn)c cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: bsupumns =>NULL()  !<  half mesh: B^u(mn)s sin(mu+nv)            asymmetric part
 real(DP), dimension(:,:), pointer :: bsupvmnc =>NULL()  !<  half mesh: B^v(mn)c cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:), pointer :: bsupvmns =>NULL()  !<  half mesh: B^v(mn)s sin(mu+nv)            asymmetric part
end type vmec_wout_cdf_data

type vmec_fort8_data
! data object to hold the information of a vmec_wout-file
 integer :: mpol      !<  number of poloidal harmonics: m = 0,...,mpol-1
 integer :: ntor      !<  number of toroidal harmonics: n = -ntor, ntor
 integer :: ns        !<  number of radial grid points (flux surfaces)
 integer :: nfp       !<  number of field periods
 integer :: mnmax     !<  number of fourier coefficients per flux surface\n
                      !!    mnmax = ntor + (mpol-1)*(1+2*ntor)
 integer :: itfsq     !<  itfsq = number of points for plotting fsq, however not used.
 integer :: nit       !<  niter/100 + 1  = for plotting of convergence, however not used.
 real(DP) :: gam      !<  ratio of specific heats used in vmec for energy principle
 real(DP), dimension(:), pointer :: iotas =>NULL()  !<  half mesh: iota
 real(DP), dimension(:), pointer :: mass =>NULL()   !<  half mesh: mass-profile
 real(DP), dimension(:), pointer :: pres =>NULL()   !<  half mesh: pressure
 real(DP), dimension(:), pointer :: phips =>NULL()  !<  half mesh: dPhi/ds
 real(DP), dimension(:), pointer :: bpco =>NULL()   !<  half mesh: toroidal current
 real(DP), dimension(:), pointer :: bzco =>NULL()   !<  half mesh: poloidal current
 real(DP), dimension(:), pointer :: phi =>NULL()    !<  half mesh: toroidal flux
 real(DP), dimension(:), pointer :: vp =>NULL()     !<  half mesh: dV/ds
 real(DP), dimension(:), pointer :: jtheta =>NULL() !<  half mesh: poloidal current density
 real(DP), dimension(:), pointer :: jzeta =>NULL()  !<  half mesh: toroidal current density
 real(DP), dimension(:), pointer :: specw =>NULL()  !<  half mesh: spectral width of fourier representation
 real(DP), dimension(:,:,:), pointer :: rmnc =>NULL()    !<  full mesh: Rmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: zmns =>NULL()    !<  full mesh: zmns sin(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: lmns =>NULL()    !<  full mesh: lmns sin(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bmod =>NULL()    !<  half mesh: Bmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: gmod =>NULL()    !<  half mesh: gmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsubu =>NULL()   !<  half mesh: B_umnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsubv =>NULL()   !<  half mesh: B_vmnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsubs =>NULL()   !<  half mesh: B_smns sin(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsupu =>NULL()   !<  half mesh: B^umnc cos(mu+nv) stellarator symmetric part
 real(DP), dimension(:,:,:), pointer :: bsupv =>NULL()   !<  half mesh: B^vmnc cos(mu+nv) stellarator symmetric part
end type vmec_fort8_data

type vmec_eqdata
! data object to hold the information of a vmec_wout-file for use in
! post-processing codes like prout, jmc etc..
! Tries to decouple the data needed in the post-processing codes from
! the vmec-output-file format used.
! Allows treatment of non-stellarator-symmetric configurations.
 integer :: mpol      !<  number of poloidal harmonics: m = 0,...,mpol-1
 integer :: ntor      !<  number of toroidal harmonics: n = -ntor, ntor
 integer :: ns        !<  number of radial grid points (flux surfaces)
 integer :: nfp       !<  number of field periods
 integer :: mnmax     !<  number of fourier coefficients per flux surface\n
                      !!    mnmax = ntor + (mpol-1)*(1+2*ntor)
 integer :: itfsq     !<  itfsq = number of points for plotting fsq, however not used.
 integer :: nit       !<  niter/100 + 1  = for plotting of convergence, however not used.
 integer :: iasym     !<  defines whether stellarator symmetry was used:\n
                      !!    iasym = 0 <==> stellarator symmetry\n
                      !!          = 1 <==> no symmetry
 real(DP) :: gam      !<  ratio of specific heats used in vmec for energy principle

 ! MJL 20150129 The original code from Joachim Geiger did not include these next fields in vmec_eqdata:
 real(DP) :: Aminor_p     !<  minor plasma radius
 real(DP), dimension(:), pointer :: iotaf =>NULL()  !<  full mesh: iota

 real(DP), dimension(:), pointer :: iotas =>NULL()  !<  half mesh: iota
 real(DP), dimension(:), pointer :: mass =>NULL()   !<  half mesh: mass-profile
 real(DP), dimension(:), pointer :: pres =>NULL()   !<  half mesh: pressure
 real(DP), dimension(:), pointer :: phips =>NULL()  !<  half mesh: dPhi/ds
 real(DP), dimension(:), pointer :: buco =>NULL()   !<  half mesh: toroidal current
 real(DP), dimension(:), pointer :: bvco =>NULL()   !<  half mesh: poloidal current

! MJL 20150128 The original code from Joachim Geiger converted phi from VMEC's full mesh to the half mesh.
! I've modified the code so phi remains in the full mesh here:
 real(DP), dimension(:), pointer :: phi =>NULL()    !<  FULL mesh: toroidal flux

 real(DP), dimension(:), pointer :: vp =>NULL()     !<  half mesh: dV/ds
 real(DP), dimension(:), pointer :: jcuru =>NULL()  !<  half mesh: jpol(dens)
 real(DP), dimension(:), pointer :: jcurv =>NULL()  !<  half mesh: jtor(dens)
 real(DP), dimension(:), pointer :: specw =>NULL()  !<  half mesh: spectral width of fourier representation
 real(DP), dimension(:,:,:), pointer :: rmnc =>NULL()    !<  full mesh: Rmnc cos(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: rmns =>NULL()    !<  full mesh: Rmns sin(mu+nv) stellarator asymmetric part
 real(DP), dimension(:,:,:), pointer :: zmns =>NULL()    !<  full mesh: zmns sin(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: zmnc =>NULL()    !<  full mesh: zmnc cos(mu+nv) stellarator asymmetric part
 real(DP), dimension(:,:,:), pointer :: lmns =>NULL()    !<  full mesh: lmns sin(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: lmnc =>NULL()    !<  full mesh: lmnc cos(mu+nv) stellarator asymmetric part
 real(DP), dimension(:,:,:), pointer :: bmnc =>NULL()    !<  half mesh: Bmnc cos(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: bmns =>NULL()    !<  half mesh: Bmns sin(mu+nv) stellarator asymmetric part
 real(DP), dimension(:,:,:), pointer :: gmnc =>NULL()    !<  half mesh: gmnc cos(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: gmns =>NULL()    !<  half mesh: gmns sin(mu+nv) stellarator asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsubumnc =>NULL()   !<  half mesh: B_umnc cos(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: bsubumns =>NULL()   !<  half mesh: B_umns sin(mu+nv) stellarator asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsubvmnc =>NULL()   !<  half mesh: B_vmnc cos(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: bsubvmns =>NULL()   !<  half mesh: B_vmns sin(mu+nv) stellarator asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsubsmnc =>NULL()   !<  half mesh: B_smns cos(mu+nv) stellarator asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsubsmns =>NULL()   !<  half mesh: B_smns sin(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: bsupumnc =>NULL()   !<  half mesh: B^umnc cos(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: bsupumns =>NULL()   !<  half mesh: B^umns sin(mu+nv) stellarator asymmetric part
 real(DP), dimension(:,:,:), pointer :: bsupvmnc =>NULL()   !<  half mesh: B^vmnc cos(mu+nv) stellarator  symmetric part
 real(DP), dimension(:,:,:), pointer :: bsupvmns =>NULL()   !<  half mesh: B^vmns sin(mu+nv) stellarator asymmetric part
end type vmec_eqdata

!-------- parameters used in type vmec_indata
 integer, parameter :: mpold = 61     !< maximum number of poloidal harmonics (in r,z,lam fourier series)
 integer, parameter :: ntord = 61     !< maximum number of toroidal harmonics
 integer, parameter :: ndatafmax = 101
!
!     DERIVED (FROM FUNDAMENTAL) PARAMETERS FOR VMEC CODE
!
 integer, parameter :: mpol1d = mpold - 1
 integer, parameter :: nmse = 100        !< number of mse measurements
 integer, parameter :: ntse = 100        !< number of thompson scattering measurements
 integer, parameter :: nfloops = 100     !< number of external poloidal flux loops
 integer, parameter :: nbsetsp = 5       !< number of external b-field loop sets allowed
 integer, parameter :: nbcoilsp = 100    !< number of external b-field coils per set
 integer, parameter :: nigroup = 100     !< number of external current groups
 integer, parameter :: nsdim = 100     !< number of external current groups

type vmec_indata  !<  data structure containing the namelist of the vmec input data
 integer :: nfp, ncurr, nsin, niter, nstep, nvacskip, mpol, ntor, &
            ntheta, nzeta &
          , ipmass, ipiota, ipcurr  !< added for profile control, J.Geiger
 integer, dimension(nsdim) :: ns_array
 integer :: imse, isnodes, itse, ipnodes, iopt_raxis, &
            imatch_phiedge, nflxs
 integer, dimension(nbsetsp) :: nbfld
 integer, dimension(nfloops) :: indxflx
 integer, dimension(nbcoilsp,nbsetsp) :: indxbfld
 real(DP), dimension(-ntord:ntord,0:mpol1d) :: rbs, zbc, rbc, zbs
 real(DP) :: time_slice, curtor, delt, ftol, tcon0, &
                gamma, phiedge, phidiam, sigma_current, sigma_delphid, tensi, &
                tensp, tensi2, fpolyi, presfac, mseangle_offset, pres_offset, &
                mseangle_offsetm, spres_ped, bloat, pres_scale
 real(DP), dimension(0:12) :: am, ai, ac, aphi  !< introduction of more coefficients, J.Geiger
 real(DP), dimension(0:ntord) :: raxis, zaxis                !< Backwards compatibility: Obsolete
 real(DP), dimension(0:ntord) :: raxis_cc, raxis_cs, &
                                    zaxis_cc, zaxis_cs
 real(DP), dimension(nsdim) :: ftol_array
 real(DP), dimension(nigroup) :: extcur
 real(DP), dimension(nmse) :: mseprof
 real(DP), dimension(ntse) :: rthom, datathom, sigma_thom
 real(DP), dimension(nmse) :: rstark, datastark, sigma_stark
 real(DP), dimension(nfloops) :: dsiobt, sigma_flux
 real(DP), dimension(nbcoilsp,nbsetsp) :: bbc, sigma_b
 real(DP), dimension(ndatafmax) :: psa, pfa, isa, ifa
 logical :: lpofr, lmac, lfreeb, lrecon, loldout, ledge_dump, lasym
 logical :: laddout, ldiagno, lmoreiter    !< added by J.Geiger
 logical :: lspectrum_dump, loptim           !< Obsolete
 character*(200) :: mgrid_file
 character*(120) :: arg1
 character*(132) :: input_extension
end type vmec_indata

!type four_array
!end type four_array

contains
!-------------------------------------------------------------------------------

subroutine read_fort8_data(filename, d)
! the subroutine reads the formatted fort.8 file of vmec/nemec and stores the
! information in a data structure of type vmec_fort8_data.
use pisa_io
implicit none
character(len=*), intent(in) :: filename
type(vmec_fort8_data), intent(inout) :: d
integer :: iunit
integer :: mpol,ntor,ns,nfp,mnmax,itfsq,nit
real(DP) :: gam
character :: c_dummy
integer :: ios
real(DP), dimension(:), pointer :: p_to_iotas
real(DP), dimension(:), pointer :: p_to_mass
real(DP), dimension(:), pointer :: p_to_pres
real(DP), dimension(:), pointer :: p_to_phips
real(DP), dimension(:), pointer :: p_to_bpco
real(DP), dimension(:), pointer :: p_to_bzco
real(DP), dimension(:), pointer :: p_to_phi
real(DP), dimension(:), pointer :: p_to_vp
real(DP), dimension(:), pointer :: p_to_jtheta
real(DP), dimension(:), pointer :: p_to_jzeta
real(DP), dimension(:), pointer :: p_to_specw
real(DP), dimension(:,:,:), pointer :: p_to_rmnc
real(DP), dimension(:,:,:), pointer :: p_to_zmns
real(DP), dimension(:,:,:), pointer :: p_to_lmns
real(DP), dimension(:,:,:), pointer :: p_to_bmod
real(DP), dimension(:,:,:), pointer :: p_to_gmod
real(DP), dimension(:,:,:), pointer :: p_to_bsubu
real(DP), dimension(:,:,:), pointer :: p_to_bsubv
real(DP), dimension(:,:,:), pointer :: p_to_bsubs
real(DP), dimension(:,:,:), pointer :: p_to_bsupu
real(DP), dimension(:,:,:), pointer :: p_to_bsupv
integer :: xm, xn
integer :: j,m,n  !< loop variables
integer :: nmin, icount

iunit = get_next_io_unit()
open(unit=iunit, file=filename, action="read")
icount=0
do
  read(iunit,*,iostat=ios)gam,nfp,ns,mpol,ntor,mnmax,itfsq,nit
  if(ios /= 0)exit
  icount=icount+1
enddo
rewind(iunit)
do j=1,icount
  read(iunit,*,iostat=ios)gam,nfp,ns,mpol,ntor,mnmax,itfsq,nit
enddo
!print *, "mpol = ",mpol
!print *, "ntor = ",ntor
!print *, "ns = ",ns
! allocate and initialize
allocate(p_to_rmnc(-ntor:ntor,0:mpol-1,ns),p_to_zmns(-ntor:ntor,0:mpol-1,ns), &
         p_to_lmns(-ntor:ntor,0:mpol-1,ns),&
         p_to_bmod(-ntor:ntor,0:mpol-1,ns),p_to_gmod(-ntor:ntor,0:mpol-1,ns),&
         p_to_bsubu(-ntor:ntor,0:mpol-1,ns),p_to_bsubv(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsubs(-ntor:ntor,0:mpol-1,ns),&
         p_to_bsupu(-ntor:ntor,0:mpol-1,ns),p_to_bsupv(-ntor:ntor,0:mpol-1,ns))
p_to_rmnc  = 0.0_dp ; p_to_zmns  = 0.0_dp ; p_to_lmns  = 0.0_dp
p_to_bmod  = 0.0_dp ; p_to_gmod  = 0.0_dp
p_to_bsubu = 0.0_dp ; p_to_bsubv = 0.0_dp ; p_to_bsubs = 0.0_dp
p_to_bsupu = 0.0_dp ; p_to_bsupv = 0.0_dp

! allocate and initialize
allocate(p_to_iotas(ns),p_to_mass(ns),p_to_pres(ns), &
         p_to_phips(ns),p_to_bpco(ns),p_to_bzco(ns),&
         p_to_phi(ns),p_to_vp(ns),p_to_jtheta(ns), &
         p_to_jzeta(ns),p_to_specw(ns))
p_to_iotas = 0.0_dp ; p_to_mass = 0.0_dp ; p_to_pres = 0.0_dp
p_to_phips = 0.0_dp ; p_to_bpco = 0.0_dp ; p_to_bzco = 0.0_dp
p_to_phi   = 0.0_dp ; p_to_vp   = 0.0_dp
p_to_jtheta= 0.0_dp ; p_to_jzeta= 0.0_dp
p_to_specw = 0.0_dp

do j=1,ns
  do m=0,mpol-1
    nmin=-ntor
    if(m == 0)nmin=0
    do n=nmin, ntor
      if(j == 1) read(iunit,'(2i10)')xm , xn
      read(iunit,'(5e20.13)') p_to_rmnc(n,m,j), p_to_zmns(n,m,j), p_to_lmns(n,m,j), &
                              p_to_bmod(n,m,j), p_to_gmod(n,m,j), &
                             p_to_bsubu(n,m,j),p_to_bsubv(n,m,j),p_to_bsubs(n,m,j), &
                             p_to_bsupu(n,m,j),p_to_bsupv(n,m,j)
    enddo
  enddo
enddo
read(iunit,'(5e20.13)')(p_to_iotas(j),p_to_mass(j),p_to_pres(j),p_to_phips(j), &
                        p_to_bpco(j),p_to_bzco(j),p_to_phi(j),p_to_vp(j), &
                        p_to_jtheta(j),p_to_jzeta(j),p_to_specw(j),j=2,ns)

close(iunit)
! put data into structure
d=vmec_fort8_data(mpol,ntor,ns,nfp,mnmax,itfsq,nit,gam,p_to_iotas,p_to_mass,p_to_pres, &
                  p_to_phips,p_to_bpco,p_to_bzco,p_to_phi,p_to_vp,p_to_jtheta, &
                  p_to_jzeta,p_to_specw, &
                  p_to_rmnc,p_to_zmns,p_to_lmns,p_to_bmod,p_to_gmod, &
                  p_to_bsubu,p_to_bsubv,p_to_bsubs,p_to_bsupu,p_to_bsupv)

end subroutine read_fort8_data

!-------------------------------------------------------------------------------
!> the subroutine reads the ascii wout-output file of vmec2000 and stores
!! the data in a data structure of type vmec_wout_data together with the version
!! number.
subroutine read_wout_txt_data(filename, wout_data)
use pisa_io
implicit none
character(len=*), intent(in) :: filename
type(vmec_wout_data), intent(inout) :: wout_data
integer :: iunit
character :: c_dummy
character(len=15) :: c_vmec_version_string
character(len=10) :: c_version_number
real(DP) :: r_version_number
integer :: ios
real(DP) :: wb
real(DP) :: wp
real(DP) :: gam
real(DP) :: pfac
real(DP) :: rmax_surf
real(DP) :: rmin_surf
real(DP) :: zmax_surf
real(DP) :: aspect
real(DP) :: betatot
real(DP) :: betapol
real(DP) :: betator
real(DP) :: betaxis
real(DP) :: b0
real(DP) :: IonLarmor
real(DP) :: VolAvgB
real(DP) :: rbtor0
real(DP) :: rbtor
real(DP) :: ctor
real(DP) :: Aminor_p
real(DP) :: Rmajor_p
real(DP) :: volume_p
integer :: nfp
integer :: ns
integer :: mpol
integer :: ntor
integer :: mnmax
integer :: mnmax_nyq0
integer :: itfsq
integer :: niter
integer :: iasym
integer :: ireconstruct
integer :: ier_flag
integer :: imse2_minus_1
integer :: itse
integer :: nbsets
integer :: nobd
integer :: nextcur
integer :: nstore_seq
integer :: signgs
character(len=200) :: mgrid_file
character(len=132) :: input_extension
character(len=30), dimension(:), pointer :: curlabel
real(DP), dimension(:), pointer :: extcur   =>NULL()
real(DP), dimension(:), pointer :: xm       =>NULL()
real(DP), dimension(:), pointer :: xn       =>NULL()
real(DP), dimension(:), pointer :: xm_nyq   =>NULL()
real(DP), dimension(:), pointer :: xn_nyq   =>NULL()
real(DP), dimension(:), pointer :: iotaf    =>NULL()
real(DP), dimension(:), pointer :: presf    =>NULL()
real(DP), dimension(:), pointer :: phipf    =>NULL()
real(DP), dimension(:), pointer :: phi      =>NULL()     !<  full mesh
real(DP), dimension(:), pointer :: jcuru    =>NULL()     !<  full mesh
real(DP), dimension(:), pointer :: jcurv    =>NULL()     !<  full mesh
real(DP), dimension(:), pointer :: iotas    =>NULL()
real(DP), dimension(:), pointer :: mass     =>NULL()
real(DP), dimension(:), pointer :: pres     =>NULL()
real(DP), dimension(:), pointer :: beta_vol =>NULL()
real(DP), dimension(:), pointer :: phip     =>NULL()
real(DP), dimension(:), pointer :: buco     =>NULL()
real(DP), dimension(:), pointer :: bvco     =>NULL()
real(DP), dimension(:), pointer :: vp       =>NULL()
real(DP), dimension(:), pointer :: overr    =>NULL()
real(DP), dimension(:), pointer :: specw    =>NULL()
real(DP), dimension(:), pointer :: Dmerc    =>NULL()
real(DP), dimension(:), pointer :: Dshear   =>NULL()
real(DP), dimension(:), pointer :: Dwell    =>NULL()
real(DP), dimension(:), pointer :: Dcurr    =>NULL()
real(DP), dimension(:), pointer :: Dgeod    =>NULL()
real(DP), dimension(:), pointer :: equif    =>NULL()
real(DP), dimension(:,:,:), pointer :: rmnc =>NULL()
real(DP), dimension(:,:,:), pointer :: rmns =>NULL()
real(DP), dimension(:,:,:), pointer :: zmns =>NULL()
real(DP), dimension(:,:,:), pointer :: zmnc =>NULL()
real(DP), dimension(:,:,:), pointer :: lmns =>NULL()  !<  lmns on half mesh
real(DP), dimension(:,:,:), pointer :: lmnc =>NULL()  !<  lmnc on half mesh
real(DP), dimension(:,:,:), pointer :: bmnc =>NULL()
real(DP), dimension(:,:,:), pointer :: bmns =>NULL()
real(DP), dimension(:,:,:), pointer :: gmnc =>NULL()
real(DP), dimension(:,:,:), pointer :: gmns =>NULL()
real(DP), dimension(:,:,:), pointer :: bsubumnc =>NULL()
real(DP), dimension(:,:,:), pointer :: bsubumns =>NULL()
real(DP), dimension(:,:,:), pointer :: bsubvmnc =>NULL()
real(DP), dimension(:,:,:), pointer :: bsubvmns =>NULL()
real(DP), dimension(:,:,:), pointer :: bsubsmnc =>NULL()
real(DP), dimension(:,:,:), pointer :: bsubsmns =>NULL()
real(DP), dimension(:,:,:), pointer :: bsupumnc =>NULL()
real(DP), dimension(:,:,:), pointer :: bsupumns =>NULL()
real(DP), dimension(:,:,:), pointer :: bsupvmnc =>NULL()
real(DP), dimension(:,:,:), pointer :: bsupvmns =>NULL()
real(DP), dimension(:,:,:), pointer :: currvmnc =>NULL()
integer :: nbfld_dummy, i, md,nd
integer :: j, m, n, nmin
real(DP) :: eps_w = 1.0e-4
character(80),parameter :: format730="(5e20.13)"

signgs=1 ! Added 20160206 by MJL
iunit = get_next_io_unit()
open(unit=iunit, file=filename, action="read")
read(iunit,'(a15,(a))',iostat=ios)c_vmec_version_string,c_version_number
! MJL 20150202 Commenting out some write statements here.
!write(6,*) "ios = ",ios
!write(6,*) c_vmec_version_string,"|",trim(c_version_number),"|"
read(c_version_number(1:len(trim(c_version_number))),*)r_version_number
!write(6,*)len(trim(c_version_number))
!write(6,*)"Version im real-format: ",r_version_number

read(iunit,*)wb, wp, gam, pfac, rmax_surf, rmin_surf, zmax_surf
if(r_version_number > 8.00 ) then
  read(iunit,*)nfp, ns, mpol, ntor, mnmax, mnmax_nyq0, &
               itfsq, niter, iasym, ireconstruct, ier_flag
else
  read(iunit,*)nfp, ns, mpol, ntor, mnmax, &
               itfsq, niter, iasym, ireconstruct, ier_flag
  mnmax_nyq0=mnmax ! both are the same
endif
read(iunit,*)imse2_minus_1,itse, nbsets, nobd, nextcur, nstore_seq
if(nbsets .gt. 0) then   !these values are not stored in the data structure.
 read(iunit,*)(nbfld_dummy,i=1,nbsets)
endif
read(iunit,'(a)') mgrid_file
! allocate arrays and initialize them
allocate(xm(mnmax),xn(mnmax),xm_nyq(mnmax_nyq0),xn_nyq(mnmax_nyq0))
allocate(rmnc(-ntor:ntor,0:mpol-1,ns),zmns(-ntor:ntor,0:mpol-1,ns), &
         lmns(-ntor:ntor,0:mpol-1,ns),bmnc(-ntor:ntor,0:mpol-1,ns), &
         gmnc(-ntor:ntor,0:mpol-1,ns),bsubumnc(-ntor:ntor,0:mpol-1,ns), &
         bsubvmnc(-ntor:ntor,0:mpol-1,ns),bsubsmns(-ntor:ntor,0:mpol-1,ns), &
         bsupumnc(-ntor:ntor,0:mpol-1,ns),bsupvmnc(-ntor:ntor,0:mpol-1,ns), &
         currvmnc(-ntor:ntor,0:mpol-1,ns))
rmnc = 0.0_dp ; zmns =0.0_dp ; lmns =0.0_dp
bmnc = 0.0_dp ; gmnc =0.0_dp
bsubumnc = 0.0_dp ; bsubvmnc =0.0_dp ; bsubsmns =0.0_dp
bsupumnc = 0.0_dp ; bsupvmnc =0.0_dp
currvmnc =0.0_dp
if(iasym == 1) then
  allocate(rmns(-ntor:ntor,0:mpol-1,ns),zmnc(-ntor:ntor,0:mpol-1,ns), &
           lmnc(-ntor:ntor,0:mpol-1,ns),bmns(-ntor:ntor,0:mpol-1,ns), &
           gmns(-ntor:ntor,0:mpol-1,ns),bsubumns(-ntor:ntor,0:mpol-1,ns), &
           bsubvmns(-ntor:ntor,0:mpol-1,ns),bsubsmnc(-ntor:ntor,0:mpol-1,ns), &
           bsupumns(-ntor:ntor,0:mpol-1,ns),bsupvmns(-ntor:ntor,0:mpol-1,ns))
  rmns = 0.0_dp ; zmnc =0.0_dp ; lmnc =0.0_dp
  bmns = 0.0_dp ; gmns =0.0_dp
  bsubumns = 0.0_dp ; bsubvmns =0.0_dp ; bsubsmnc =0.0_dp
  bsupumns = 0.0_dp ; bsupvmns =0.0_dp
endif
! read fourier coefficients
if(r_version_number > 8.00 ) then
  do j=1,ns
    do m=0,mpol-1
      nmin=-ntor
      if(m==0)nmin=0
      do n=nmin,ntor
        if(j == 1)then
          read(iunit,*)md,nd
        endif
        read(iunit,*)rmnc(n,m,j),zmns(n,m,j),lmns(n,m,j)
        if(iasym == 1)then
          read(iunit,*)rmns(n,m,j),zmnc(n,m,j),lmnc(n,m,j)
        endif
      enddo
    enddo
    do m=0,mpol-1     ! could be a problem if vmec2000 was run with lnyquist==.true.
      nmin=-ntor
      if(m==0)nmin=0
      do n=nmin,ntor
        if(j == 1)then
          read(iunit,*)md,nd
        endif
        read(iunit,*)bmnc(n,m,j),gmnc(n,m,j), &
                     bsubumnc(n,m,j),bsubvmnc(n,m,j),bsubsmns(n,m,j), &
                     bsupumnc(n,m,j),bsupvmnc(n,m,j)
        if(iasym == 1)then
          read(iunit,*)bmns(n,m,j),gmns(n,m,j), &
                       bsubumns(n,m,j),bsubvmns(n,m,j),bsubsmnc(n,m,j), &
                       bsupumns(n,m,j),bsupvmns(n,m,j)
        endif
      enddo
    enddo
  enddo
else   !for version <= 8.00
  do j=1,ns
    do m=0,mpol-1
      nmin=-ntor
      if(m==0)nmin=0
      do n=nmin,ntor
        if(j == 1)then
          read(iunit,*)md,nd
        endif
        read(iunit,*)rmnc(n,m,j),zmns(n,m,j),lmns(n,m,j), &
                     bmnc(n,m,j),gmnc(n,m,j), &
                     bsubumnc(n,m,j),bsubvmnc(n,m,j),bsubsmns(n,m,j), &
                     bsupumnc(n,m,j),bsupvmnc(n,m,j),currvmnc(n,m,j)
        if(iasym == 1)then
          read(iunit,*)rmns(n,m,j),zmnc(n,m,j),lmnc(n,m,j), &
                       bmns(n,m,j),gmns(n,m,j), &
                       bsubumns(n,m,j),bsubvmns(n,m,j),bsubsmnc(n,m,j), &
                       bsupumns(n,m,j),bsupvmns(n,m,j)
        endif
      enddo
    enddo
  enddo
endif

! allocate arrays and initialize them
  allocate(iotaf(ns),presf(ns),phipf(ns),phi(ns),jcuru(ns),jcurv(ns))
  allocate(iotas(ns),mass(ns),pres(ns),beta_vol(ns),phip(ns),buco(ns),bvco(ns), &
           vp(ns),overr(ns),specw(ns))
  iotaf    = 0.0_dp ; presf = 0.0_dp ; phipf = 0.0_dp
  phi      = 0.0_dp ; jcuru = 0.0_dp ; jcurv = 0.0_dp
  iotas    = 0.0_dp ; mass  = 0.0_dp ; pres  = 0.0_dp
  beta_vol = 0.0_dp ; phip  = 0.0_dp
  buco     = 0.0_dp ; bvco  = 0.0_dp ; vp    = 0.0_dp
  overr    = 0.0_dp ; specw = 0.0_dp

! reading of profiles.
! This part has been copied from the read_wout_mod.f-file of LIBSTELL.
if(real(r_version_number) <= 6.05+eps_w ) then
  read(iunit,format730,iostat=ios)(iotas(j),mass(j),pres(j),phip(j),buco(j),bvco(j),phi(j), &
              vp(j),overr(j),jcuru(j),jcurv(j),specw(j),j=2,ns)
  read(iunit,format730,iostat=ios)aspect, betatot, betapol, betator, betaxis, b0
else if(real(r_version_number) <= 6.20+eps_w ) then
  read(iunit,format730,iostat=ios)(iotas(j),mass(j),pres(j),beta_vol(j),phip(j),buco(j),bvco(j),phi(j), &
              vp(j),overr(j),jcuru(j),jcurv(j),specw(j),j=2,ns)
  read(iunit,format730,iostat=ios)aspect, betatot, betapol, betator, betaxis, b0
else if(real(r_version_number) <= 6.95+eps_w ) then
  read(iunit,*,iostat=ios)(iotas(j),mass(j),pres(j),beta_vol(j),phip(j),buco(j),bvco(j),phi(j), &
              vp(j),overr(j),jcuru(j),jcurv(j),specw(j),j=2,ns)
  read(iunit,*,iostat=ios)aspect, betatot, betapol, betator, betaxis, b0
else
  read(iunit,*,iostat=ios)(iotaf(j),presf(j),phipf(j),phi(j),jcuru(j),jcurv(j),j=1,ns)
  read(iunit,*,iostat=ios)(iotas(j),mass(j),pres(j),beta_vol(j),phip(j),buco(j),bvco(j), &
              vp(j),overr(j),specw(j),j=2,ns)
  read(iunit,*,iostat=ios)aspect, betatot, betapol, betator, betaxis, b0
endif

! from read_wout_mod.f (LIBSTELL)
if(real(r_version_number) > 6.10+eps_w) then
  !write(6,*) "Read on !"
  read(iunit,*,iostat=ios) i  ! nint(signgs) not stored.
  read(iunit, '(a)',iostat=ios) input_extension
  !write(6,*) input_extension
  read(iunit,*) IonLarmor, VolAvgB, rbtor0, rbtor, ctor, &
              Aminor_p, Rmajor_p, volume_p
endif

! Mercier criterion
allocate(Dmerc(ns),Dshear(ns),Dwell(ns),Dcurr(ns),Dgeod(ns),equif(ns))
Dmerc = 0.0_dp ; Dshear = 0.0_dp ; Dwell = 0.0_dp
Dcurr = 0.0_dp ; Dgeod  = 0.0_dp ; equif = 0.0_dp
if(real(r_version_number) > 5.10+eps_w .and. real(r_version_number) < 6.20-eps_w) then
  read(iunit, format730) (Dmerc(j), Dshear(j), Dwell(j), Dcurr(j), &
                  Dgeod(j), equif(j), j=2,ns-1)
else if(real(r_version_number) >= 6.20-eps_w) then
  read(iunit, *) (Dmerc(j), Dshear(j), Dwell(j), Dcurr(j), &
                  Dgeod(j), equif(j), j=2,ns-1)
endif

! in case of free-boundary
if(nextcur > 0) then
  allocate(extcur(nextcur),curlabel(nextcur))
  extcur = 0.0_dp
  read(iunit,*)(extcur(i),i=1,nextcur)
  read(iunit,*)(curlabel(i),i=1,nextcur)
endif
close(iunit)

wout_data = vmec_wout_data(c_version_number, wb, wp, gam, pfac, rmax_surf, &
                           rmin_surf, zmax_surf, aspect, betatot, &
                           betapol, betator, betaxis, b0, &
                           IonLarmor, VolAvgB, rbtor0, rbtor, ctor, &
                           Aminor_p, Rmajor_p, volume_p, &
                           nfp, ns, mpol, ntor, mnmax, mnmax_nyq0, itfsq, niter, &
                           iasym, ireconstruct, ier_flag, imse2_minus_1, itse, &
                           nbsets, nobd, nextcur, nstore_seq, signgs, &
                           mgrid_file, input_extension, &
                           curlabel, extcur, &
                           iotaf, presf, phipf, phi, jcuru, jcurv, &
                           iotas, mass, pres, beta_vol, phip, buco, &
                           bvco, vp, overr, specw, &
                           Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, &
                           xm, xn, xm_nyq, xn_nyq, &
                           rmnc, rmns, zmns, zmnc, lmns, lmnc, bmnc, bmns, gmnc, gmns, &
                           bsubumnc, bsubumns, bsubvmnc, bsubvmns, bsubsmnc, bsubsmns, &
                           bsupumnc, bsupumns, bsupvmnc, bsupvmns, currvmnc)

end subroutine read_wout_txt_data

!-------------------------------------------------------------------------------
!> the subroutine reads the netcdf wout-output file of vmec2000 and stores
!! the data in a data structure of type vmec_wout_cdf_data together with the version
!! number.
#ifdef USE_NETCDF
subroutine read_wout_cdf_data(filename, wout_data)
use pisa_io
use ezcdf
implicit none
character(len=*), intent(in) :: filename
type(vmec_wout_cdf_data), intent(inout) :: wout_data
integer :: iunit
character :: c_dummy
character(len=15) :: c_vmec_version_string
character(len=10) :: c_version_number
real(DP) :: r_version_number
integer :: ios
real(DP) :: wb
real(DP) :: wp
real(DP) :: gam
real(DP) :: pfac
real(DP) :: rmax_surf
real(DP) :: rmin_surf
real(DP) :: zmax_surf
real(DP) :: aspect
real(DP) :: betatot
real(DP) :: betapol
real(DP) :: betator
real(DP) :: betaxis
real(DP) :: b0
real(DP) :: IonLarmor
real(DP) :: VolAvgB
real(DP) :: rbtor0
real(DP) :: rbtor
real(DP) :: ctor
real(DP) :: Aminor_p
real(DP) :: Rmajor_p
real(DP) :: volume_p
integer :: nfp
integer :: ns
integer :: mpol
integer :: ntor
integer :: mnmax
integer :: mnmax_nyq0
integer :: itfsq
integer :: niter
integer :: iasym
integer :: ireconstruct
integer :: ifreeb
integer :: ier_flag
integer :: imse2_minus_1
integer :: itse
integer :: nbsets
integer :: nobd
integer :: nextcur
integer :: nstore_seq
integer :: signgs
character(len=200) :: mgrid_file
character(len=132) :: input_extension
character(len=30), dimension(:), pointer :: curlabel
real(DP), dimension(:), pointer :: extcur   =>NULL()
real(DP), dimension(:), pointer :: xm       =>NULL()
real(DP), dimension(:), pointer :: xn       =>NULL()
real(DP), dimension(:), pointer :: xm_nyq   =>NULL()
real(DP), dimension(:), pointer :: xn_nyq   =>NULL()
real(DP), dimension(:), pointer :: iotaf    =>NULL()
real(DP), dimension(:), pointer :: presf    =>NULL()
real(DP), dimension(:), pointer :: phipf    =>NULL()
real(DP), dimension(:), pointer :: phi      =>NULL()     !<  full mesh
real(DP), dimension(:), pointer :: jcuru    =>NULL()     !<  full mesh
real(DP), dimension(:), pointer :: jcurv    =>NULL()     !<  full mesh
real(DP), dimension(:), pointer :: iotas    =>NULL()
real(DP), dimension(:), pointer :: mass     =>NULL()
real(DP), dimension(:), pointer :: pres     =>NULL()
real(DP), dimension(:), pointer :: beta_vol =>NULL()
real(DP), dimension(:), pointer :: phips    =>NULL()
real(DP), dimension(:), pointer :: buco     =>NULL()
real(DP), dimension(:), pointer :: bvco     =>NULL()
real(DP), dimension(:), pointer :: vp       =>NULL()
real(DP), dimension(:), pointer :: overr    =>NULL()
real(DP), dimension(:), pointer :: specw    =>NULL()
real(DP), dimension(:), pointer :: Dmerc    =>NULL()
real(DP), dimension(:), pointer :: Dshear   =>NULL()
real(DP), dimension(:), pointer :: Dwell    =>NULL()
real(DP), dimension(:), pointer :: Dcurr    =>NULL()
real(DP), dimension(:), pointer :: Dgeod    =>NULL()
real(DP), dimension(:), pointer :: equif    =>NULL()
real(DP), dimension(:,:), pointer :: rmnc =>NULL()
real(DP), dimension(:,:), pointer :: rmns =>NULL()
real(DP), dimension(:,:), pointer :: zmns =>NULL()
real(DP), dimension(:,:), pointer :: zmnc =>NULL()
real(DP), dimension(:,:), pointer :: lmns =>NULL()  !<  lmns on half mesh
real(DP), dimension(:,:), pointer :: lmnc =>NULL()  !<  lmnc on half mesh
real(DP), dimension(:,:), pointer :: bmnc =>NULL()
real(DP), dimension(:,:), pointer :: bmns =>NULL()
real(DP), dimension(:,:), pointer :: gmnc =>NULL()
real(DP), dimension(:,:), pointer :: gmns =>NULL()
real(DP), dimension(:,:), pointer :: bsubumnc =>NULL()
real(DP), dimension(:,:), pointer :: bsubumns =>NULL()
real(DP), dimension(:,:), pointer :: bsubvmnc =>NULL()
real(DP), dimension(:,:), pointer :: bsubvmns =>NULL()
real(DP), dimension(:,:), pointer :: bsubsmnc =>NULL()
real(DP), dimension(:,:), pointer :: bsubsmns =>NULL()
real(DP), dimension(:,:), pointer :: bsupumnc =>NULL()
real(DP), dimension(:,:), pointer :: bsupumns =>NULL()
real(DP), dimension(:,:), pointer :: bsupvmnc =>NULL()
real(DP), dimension(:,:), pointer :: bsupvmns =>NULL()
integer :: nbfld_dummy, i, md,nd
integer :: j, m, n, nmin

! Begin additions by MJL 20160206
c_version_number = ''
itse = 0
imse2_minus_1 = 0
nobd = 0
nbsets = 0
nstore_seq = 0
! End MJL additions from 20160206
call cdf_open(iunit,filename,'r')

call cdf_read(iunit,'version_',r_version_number)
call cdf_read(iunit,'mgrid_file',mgrid_file)
call cdf_read(iunit,'input_extension',input_extension)
call cdf_read(iunit,'wb',wb)
call cdf_read(iunit,'wp',wp)
call cdf_read(iunit,'gamma',gam)
pfac = 0      ! not in netcdf file, may be because I did not use reconstruction mode yet.
call cdf_read(iunit,'rmax_surf',rmax_surf)
call cdf_read(iunit,'rmin_surf',rmin_surf)
call cdf_read(iunit,'zmax_surf',zmax_surf)
call cdf_read(iunit,'nfp',nfp)
call cdf_read(iunit,'ns',ns)
call cdf_read(iunit,'mpol',mpol)
call cdf_read(iunit,'ntor',ntor)
call cdf_read(iunit,'mnmax',mnmax)
call cdf_read(iunit,'mnmax_nyq',mnmax_nyq0)
call cdf_read(iunit,'niter',niter)
call cdf_read(iunit,'itfsq',itfsq)
call cdf_read(iunit,'lasym__logical__',iasym)
call cdf_read(iunit,'lrecon__logical__',ireconstruct)
call cdf_read(iunit,'lfreeb__logical__',ifreeb)
call cdf_read(iunit,'ier_flag',ier_flag)
call cdf_read(iunit,'aspect',aspect)
call cdf_read(iunit,'betatotal',betatot)
call cdf_read(iunit,'betapol',betapol)
call cdf_read(iunit,'betator',betator)
call cdf_read(iunit,'betaxis',betaxis)
call cdf_read(iunit,'b0',b0)
call cdf_read(iunit,'rbtor0',rbtor0)
call cdf_read(iunit,'rbtor',rbtor)
call cdf_read(iunit,'signgs',signgs)
call cdf_read(iunit,'IonLarmor',IonLarmor)
call cdf_read(iunit,'volavgB',volavgB)
call cdf_read(iunit,'ctor',ctor)
call cdf_read(iunit,'Aminor_p',Aminor_p)
call cdf_read(iunit,'Rmajor_p',Rmajor_p)
call cdf_read(iunit,'volume_p',volume_p)
! allocate arrays and initialize them
allocate(xm(mnmax),xn(mnmax), &
         xm_nyq(mnmax_nyq0),xn_nyq(mnmax_nyq0))
xm=0.0_dp ; xn=0.0_dp ; xm_nyq=0.0_dp ; xn_nyq=0.0_dp

allocate(rmnc(mnmax,ns),zmns(mnmax,ns), &
         lmns(mnmax,ns),bmnc(mnmax_nyq0,ns), &
         gmnc(mnmax_nyq0,ns),bsubumnc(mnmax_nyq0,ns), &
         bsubvmnc(mnmax_nyq0,ns),bsubsmns(mnmax_nyq0,ns), &
         bsupumnc(mnmax_nyq0,ns),bsupvmnc(mnmax_nyq0,ns))
rmnc = 0.0_dp ; zmns =0.0_dp ; lmns =0.0_dp
bmnc = 0.0_dp ; bmnc =0.0_dp
bsubumnc = 0.0_dp ; bsubvmnc =0.0_dp ; bsubsmns =0.0_dp
bsupumnc = 0.0_dp ; bsupvmnc =0.0_dp
if(iasym == 1) then
  allocate(rmns(mnmax,ns),zmnc(mnmax,ns), &
           lmnc(mnmax,ns),bmns(mnmax_nyq0,ns), &
           gmns(mnmax_nyq0,ns),bsubumns(mnmax_nyq0,ns), &
           bsubvmns(mnmax_nyq0,ns),bsubsmnc(mnmax_nyq0,ns), &
           bsupumns(mnmax_nyq0,ns),bsupvmns(mnmax_nyq0,ns))
  rmns = 0.0_dp ; zmnc =0.0_dp ; lmnc =0.0_dp
  bmns = 0.0_dp ; bmns =0.0_dp
  bsubumns = 0.0_dp ; bsubvmns =0.0_dp ; bsubsmnc =0.0_dp
  bsupumns = 0.0_dp ; bsupvmns =0.0_dp
endif
call cdf_read(iunit,'xm',xm)
call cdf_read(iunit,'xn',xn)
call cdf_read(iunit,'xm_nyq',xm_nyq)
call cdf_read(iunit,'xn_nyq',xn_nyq)
call cdf_read(iunit,'rmnc',rmnc)
call cdf_read(iunit,'zmns',zmns)
call cdf_read(iunit,'lmns',lmns)
call cdf_read(iunit,'gmnc',gmnc)
call cdf_read(iunit,'bmnc',bmnc)
call cdf_read(iunit,'bsubumnc',bsubumnc)
call cdf_read(iunit,'bsubvmnc',bsubvmnc)
call cdf_read(iunit,'bsubsmns',bsubsmns)
call cdf_read(iunit,'bsupumnc',bsupumnc)
call cdf_read(iunit,'bsupvmnc',bsupvmnc)
if(iasym == 1) then
  call cdf_read(iunit,'rmns',rmns)
  call cdf_read(iunit,'zmnc',zmnc)
  call cdf_read(iunit,'lmnc',lmnc)
  call cdf_read(iunit,'gmns',gmns)
  call cdf_read(iunit,'bmns',bmns)
  call cdf_read(iunit,'bsubumns',bsubumns)
  call cdf_read(iunit,'bsubvmns',bsubvmns)
  call cdf_read(iunit,'bsubsmnc',bsubsmnc)
  call cdf_read(iunit,'bsupumns',bsupumns)
  call cdf_read(iunit,'bsupvmns',bsupvmns)
endif

! allocate profile arrays and initialize them
  allocate(iotaf(ns),presf(ns),phipf(ns),phi(ns),jcuru(ns),jcurv(ns))
  allocate(iotas(ns),mass(ns),pres(ns),beta_vol(ns),phips(ns),buco(ns),bvco(ns), &
           vp(ns),overr(ns),specw(ns))
  iotaf    = 0.0_dp ; presf = 0.0_dp ; phipf = 0.0_dp
  phi      = 0.0_dp ; jcuru = 0.0_dp ; jcurv = 0.0_dp
  iotas    = 0.0_dp ; mass  = 0.0_dp ; pres  = 0.0_dp
  beta_vol = 0.0_dp ; phips = 0.0_dp
  buco     = 0.0_dp ; bvco  = 0.0_dp ; vp    = 0.0_dp
  overr    = 0.0_dp ; specw = 0.0_dp

!! read profiles
call cdf_read(iunit,'iotaf',iotaf)
call cdf_read(iunit,'presf',presf)
call cdf_read(iunit,'phi',phi)
call cdf_read(iunit,'phipf',phipf)
call cdf_read(iunit,'jcuru',jcuru)
call cdf_read(iunit,'jcurv',jcurv)
call cdf_read(iunit,'iotas',iotas)
call cdf_read(iunit,'mass',mass)
call cdf_read(iunit,'pres',pres)
call cdf_read(iunit,'beta_vol',beta_vol)
call cdf_read(iunit,'buco',buco)
call cdf_read(iunit,'bvco',bvco)
call cdf_read(iunit,'phips',phips)
call cdf_read(iunit,'vp',vp)
call cdf_read(iunit,'specw',specw)
call cdf_read(iunit,'over_r',overr)

! allocate stability profile arrays and initialize them
  allocate(Dmerc(ns),Dshear(ns),Dwell(ns),Dcurr(ns),Dgeod(ns),equif(ns))
  Dmerc = 0.0_dp ; Dshear = 0.0_dp ; Dwell = 0.0_dp
  Dcurr = 0.0_dp ; Dgeod  = 0.0_dp ; equif = 0.0_dp
! read stability profiles
call cdf_read(iunit,'DMerc',Dmerc)
call cdf_read(iunit,'DShear',Dshear)
call cdf_read(iunit,'DWell',Dwell)
call cdf_read(iunit,'DCurr',Dcurr)
call cdf_read(iunit,'DGeod',Dgeod)
call cdf_read(iunit,'equif',equif)

! read and allocate external coil currents
call cdf_read(iunit,'nextcur',nextcur)
if(nextcur > 0) then
  allocate(extcur(nextcur),curlabel(nextcur))
  extcur = 0.0_dp
  call cdf_read(iunit,'extcur',extcur)
  call cdf_read(iunit,'curlabel',curlabel)
endif
call cdf_close(iunit)

wout_data = vmec_wout_cdf_data( c_version_number, wb, wp, gam, pfac, rmax_surf, &
                    rmin_surf, zmax_surf, aspect, betatot, &
                    betapol, betator, betaxis, b0, &
                    IonLarmor, VolAvgB, rbtor0, rbtor, ctor, &
                    Aminor_p, Rmajor_p, volume_p, &
                    nfp, ns, mpol, ntor, mnmax, mnmax_nyq0, itfsq, niter, &
                    iasym, ireconstruct, ifreeb, ier_flag, imse2_minus_1, itse, &
                    nbsets, nobd, nextcur, nstore_seq, signgs, &
                    mgrid_file, input_extension, &
                    curlabel, extcur, &
                    iotaf, presf, phipf, phi, jcuru, jcurv, &
                    iotas, mass, pres, beta_vol, phips, buco, &
                    bvco, vp, overr, specw, &
                    Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, &
                    xm, xn, xm_nyq, xn_nyq, &
                    rmnc, rmns, zmns, zmnc, lmns, lmnc, bmnc, bmns, gmnc, gmns, &
                    bsubumnc, bsubumns, bsubvmnc, bsubvmns, bsubsmnc, bsubsmns, &
                    bsupumnc, bsupumns, bsupvmnc, bsupvmns)

end subroutine read_wout_cdf_data
#endif
!-------------------------------------------------------------------------------

subroutine check_wout_file_type(filename, this_type)
!> the subroutine tries to determine automatically which wout-filetype
!! is available. Currently to following types are recognised:\n
!!   - cdf  : the netcdf-type of vmec2000\n
!!   - txt  : the acsii-type of vmec2000\n
!!   - fo8  : the fort.8-type of vmec2000, identical to the ascii-output of nemec93\n
!!   - xdr  : the xdr-type of nemec93\n
!! The old fort.8-files from nemec93 which are in binary format can not yet be read.\n
!! The kind-value of the integers is 8, whereas the kind-value of the reals is not
!! yet known. The value 10 gives the correct number of bytes for a real, but the read
!! values are meaningless. The endianess is BIG!
use pisa_io
implicit none
type(file_type), intent(out) :: this_type
character(len=*), intent(in) :: filename
integer :: iunit
character(len=4) :: c_dummy
!---------------------------------------------------------------
! see http://www.unidata.ucar.edu/software/netcdf/docs/faq.html
! identification string for classic netcdf
character(len=4), parameter :: c_cdf001="CDF"//achar(iachar(' ')-31)
! identification string for classic netcdf 64bit offset format
character(len=4), parameter :: c_cdf002="CDF"//achar(iachar(' ')-30)
! identification string for HDF5 netcdf format
character(len=4), parameter :: c_hdf5=achar(iachar('z')+15)//"HDF"
!---------------------------------------------------------------
integer :: ios
integer, parameter :: magic1=-83891667  !<  xdr-file integer values at beginning
integer, parameter :: magic2= 33554432  !<  xdr-file integer values at beginning
integer  :: fo8_i
integer  :: m1,m2
integer(kind=8) :: i1,i2,i3,i4,i5,i6
real(DP) :: fo8_r
real(kind=8) :: r1,r2
!real(kind=8), dimension(10) :: r10_array
integer, parameter :: s12 = selected_real_kind(p=12)
integer, parameter :: s13 = selected_real_kind(p=13)
integer, parameter :: s14 = selected_real_kind(p=14)
integer, parameter :: s15 = selected_real_kind(p=15)
integer, parameter :: s16 = selected_real_kind(p=16)
integer, parameter :: s17 = selected_real_kind(p=17)
integer, parameter :: s18 = selected_real_kind(p=18)
real(kind=8), dimension(10) :: r10_array

iunit = get_next_io_unit()
open(unit=iunit, file=filename, action="read", iostat=ios)
if(ios .ne. 0) then
   stop "Error! Unable to open the requested equilibriumFile." !MJL 20150204
end if
read(iunit,'(a4)',iostat=ios)c_dummy
!write(6,*) "ios = ",ios
close(iunit)
if(ios < 0) then
  !if(ios == -1) stop "file not existing!"
  if(ios == -1) stop "Error! Unable to read the requested equilibriumFile."  !MJL 20150204
  write(6,*)"Input error:  ios = ",ios
  stop "Program ends with Error!!!"
  write(6,'(a)') "c_dummy = ",c_dummy
endif
if(c_dummy(1:1) .eq. "-") c_dummy(2:3)="  "  !there are otherwise unprintable characters
select case (c_dummy)
 case (c_cdf001,c_cdf002,c_hdf5)
  this_type%typo = "cdf"
 case ("VMEC")
  this_type%typo = "txt"
 case ("   ")
  open(unit=iunit, file=filename, action="read")
  read(iunit,*,iostat=ios)fo8_r,fo8_i
  close(iunit)
  if(ios == 0) then
    if(fo8_r >= 0.0_dp .and. fo8_r <= 2.0_dp .and. fo8_i >= 0) this_type%typo = "fo8"
  else
    stop "Unrecognizable file type! Tried type fort.8!"
  endif
 case ("-  ")
  open(unit=iunit, file=filename, action="read",access="stream")
  read(iunit,iostat=ios)m1,m2
  close(iunit)
  if(ios == 0) then
    if(m1 == magic1 .and. m2 == magic2) this_type%typo = "xdr"
  else
    stop "Unrecognizable file type! Tried type xdr!"
  endif
 case default
!--- binary files from old vmec-calculations are not readable. The real
!--- values do have a yet unknown kind-value. kind=10 seems to give at
!--- least the correct number of bytes, but the values are not meaningful.
  open(unit=iunit, file=filename, action="read",access="stream")
  read(iunit,iostat=ios)r1,i1,i2,i3,i4,i5,i6,r2
  if(ios == 0) then
    if( i4+i3*(1+2*(i4-1)) == i5) this_type%typo = "b93"
    write(6,*)r1,i1,i2,i3,i4,i5,i6,r2
    write(6,*) "Type would be : ", this_type%typo
    i1=-4711 ; i2=-4711
    read(iunit,iostat=ios)i1,i2
    write(6,*)"i1/i2 : ",i1,i2
    r10_array = 0.0
    read(iunit,iostat=ios)r10_array
    write(6,*)"r10_array: ",r10_array
    i1=-4711 ; i2=-4711
    read(iunit,iostat=ios)i1,i2
    write(6,*)"i1/i2 : ",i1,i2
    r10_array = 0.0
    read(iunit,iostat=ios)r10_array
    write(6,*)"r10_array: ",r10_array
    i1=-4711 ; i2=-4711
    read(iunit,iostat=ios)i1,i2
    write(6,*)"i1/i2 : ",i1,i2
    r10_array = 0.0
    read(iunit,iostat=ios)r10_array
    write(6,*)"r10_array: ",r10_array
    i1=-4711 ; i2=-4711
    read(iunit,iostat=ios)i1,i2
    write(6,*)"i1/i2 : ",i1,i2
    r10_array = 0.0
    read(iunit,iostat=ios)r10_array
    write(6,*)"r10_array: ",r10_array
    i1=-4711 ; i2=-4711
    read(iunit,iostat=ios)i1,i2
    write(6,*)"i1/i2 : ",i1,i2
    r10_array = 0.0
    read(iunit,iostat=ios)r10_array
    write(6,*)"r10_array: ",r10_array
    stop "Sorry, this type can not be read yet!"
  endif
  close(iunit)
! write(6,*)"c_dummy = ",c_dummy
  stop "Unknown type of file!"
end select
!write(6,*) "Detected type of input file: ",this_type%typo

end subroutine check_wout_file_type
!-------------------------------------------------------------------------------
function extract_minimum_eqdata_from_fort8_data(indata) result(outdata)
implicit none
type(vmec_fort8_data), intent(in) :: indata
type(vmec_minimum_eqdata) :: outdata
integer :: ns, mpol, ntor, nfp, iasym
real(DP) :: phiedge
real(DP), dimension(:), pointer :: p_to_iota => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmnc => NULL()

ns    = indata%ns
mpol  = indata%mpol
ntor  = indata%ntor
nfp   = indata%nfp

iasym = 0   ! default currently fort.8 only with stellarator symmetry!

phiedge = indata%phi(ns)
allocate(p_to_iota(ns), p_to_rmnc(-ntor:ntor,0:mpol-1,ns), p_to_zmns(-ntor:ntor,0:mpol-1,ns))
p_to_iota = 0.0_DP
p_to_rmnc = 0.0_DP
p_to_zmns = 0.0_DP
p_to_iota = indata%iotas
p_to_rmnc = indata%rmnc
p_to_zmns = indata%zmns

outdata = vmec_minimum_eqdata(ns,mpol,ntor,nfp,iasym,phiedge, &
                      p_to_iota,p_to_rmnc,p_to_rmns,p_to_zmns,p_to_zmnc)

end function extract_minimum_eqdata_from_fort8_data
!-------------------------------------------------------------------------------
function extract_minimum_eqdata_from_wout_data(indata) result(outdata)
implicit none
type(vmec_wout_data), intent(in) :: indata
type(vmec_minimum_eqdata) :: outdata
integer :: ns, mpol, ntor, nfp, iasym
real(DP) :: phiedge
real(DP), dimension(:), pointer :: p_to_iota => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmnc => NULL()

ns    = indata%ns
mpol  = indata%mpol
ntor  = indata%ntor
nfp   = indata%nfp

iasym = 0   ! default currently fort.8 only with stellarator symmetry!

phiedge = indata%phi(ns)

allocate(p_to_iota(ns), p_to_rmnc(-ntor:ntor,0:mpol-1,ns), p_to_zmns(-ntor:ntor,0:mpol-1,ns))
if(iasym == 1) allocate(p_to_rmns(-ntor:ntor,0:mpol-1,ns), p_to_zmnc(-ntor:ntor,0:mpol-1,ns))

p_to_iota = 0.0_DP
p_to_rmnc = 0.0_DP
p_to_zmns = 0.0_DP
p_to_iota = indata%iotas
p_to_rmnc = indata%rmnc
p_to_zmns = indata%zmns
if(iasym == 1) then
  p_to_rmns = 0.0_DP
  p_to_zmnc = 0.0_DP
  p_to_rmns = indata%rmns
  p_to_zmnc = indata%zmnc
endif

outdata = vmec_minimum_eqdata(ns,mpol,ntor,nfp,iasym,phiedge, &
                      p_to_iota,p_to_rmnc,p_to_rmns,p_to_zmns,p_to_zmnc)


end function extract_minimum_eqdata_from_wout_data
!-------------------------------------------------------------------------------
function extract_fort8_data_from_wout_data(indata) result(outdata)
implicit none
type(vmec_wout_data), intent(in) :: indata
type(vmec_fort8_data) :: outdata
integer :: ns, mpol, ntor, nfp, mnmax, itfsq, nit
real(DP) :: gam
real(DP), dimension(:), pointer :: p_to_iotas => NULL()
real(DP), dimension(:), pointer :: p_to_mass => NULL()
real(DP), dimension(:), pointer :: p_to_pres => NULL()
real(DP), dimension(:), pointer :: p_to_phips => NULL()
real(DP), dimension(:), pointer :: p_to_bpco => NULL()
real(DP), dimension(:), pointer :: p_to_bzco => NULL()
real(DP), dimension(:), pointer :: p_to_phi => NULL()
real(DP), dimension(:), pointer :: p_to_vp => NULL()
real(DP), dimension(:), pointer :: p_to_jtheta => NULL()
real(DP), dimension(:), pointer :: p_to_jzeta => NULL()
real(DP), dimension(:), pointer :: p_to_specw => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_lmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bmod => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_gmod => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubu => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubv => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubs => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupu => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupv => NULL()
integer :: j

if(indata%iasym == 1) then
 stop "No stellarator symmetry! fort8-data can not hold this information!"
endif

mpol  = indata%mpol
ntor  = indata%ntor
ns    = indata%ns
nfp   = indata%nfp
mnmax = indata%mnmax
itfsq = indata%itfsq
nit   = indata%niter/100 + 1
gam   = indata%gam

allocate(p_to_iotas(ns), p_to_mass(ns), p_to_pres(ns), p_to_phips(ns), &
         p_to_bpco(ns), p_to_bzco(ns), p_to_phi(ns), p_to_vp(ns), &
         p_to_jtheta(ns), p_to_jzeta(ns), p_to_specw(ns))
allocate(p_to_rmnc (-ntor:ntor,0:mpol-1,ns), p_to_zmns (-ntor:ntor,0:mpol-1,ns), &
         p_to_lmns (-ntor:ntor,0:mpol-1,ns), p_to_bmod (-ntor:ntor,0:mpol-1,ns), &
         p_to_gmod (-ntor:ntor,0:mpol-1,ns), p_to_bsubu(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsubv(-ntor:ntor,0:mpol-1,ns), p_to_bsubs(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsupu(-ntor:ntor,0:mpol-1,ns), p_to_bsupv(-ntor:ntor,0:mpol-1,ns))

! Initialize the data
p_to_iotas  = 0.0_DP
p_to_mass   = 0.0_DP
p_to_pres   = 0.0_DP
p_to_phips  = 0.0_DP
p_to_bpco   = 0.0_DP
p_to_bzco   = 0.0_DP
p_to_phi    = 0.0_DP
p_to_vp     = 0.0_DP
p_to_jtheta = 0.0_DP
p_to_jzeta  = 0.0_DP
p_to_specw  = 0.0_DP

p_to_rmnc = 0.0_DP
p_to_zmns = 0.0_DP
p_to_lmns = 0.0_DP
p_to_bmod = 0.0_DP
p_to_gmod = 0.0_DP
p_to_bsubu = 0.0_DP
p_to_bsubv = 0.0_DP
p_to_bsubs = 0.0_DP
p_to_bsupu = 0.0_DP
p_to_bsupv = 0.0_DP

! Transfer the data
p_to_iotas  = indata%iotas
p_to_mass   = indata%mass
p_to_pres   = indata%pres
p_to_phips  = indata%phip
p_to_bpco   = indata%buco
p_to_bzco   = indata%bvco
do j=2,ns   ! construct half-mesh phi
  p_to_phi(j) = (j-1.5_dp)/(ns-1.0_dp)
enddo
p_to_phi = p_to_phi*indata%phi(ns)
p_to_vp     = indata%vp
! fort.8 has current densities on half grid
p_to_jtheta(2:ns) = 0.5*(indata%jcuru(2:ns)+indata%jcuru(1:ns-1))
p_to_jzeta(2:ns)  = 0.5*(indata%jcurv(2:ns)+indata%jcurv(1:ns-1))
p_to_specw  = indata%specw

p_to_rmnc  = indata%rmnc
p_to_zmns  = indata%zmns
p_to_lmns  = indata%lmns
p_to_bmod  = indata%bmnc
p_to_gmod  = indata%gmnc
p_to_bsubu = indata%bsubumnc
p_to_bsubv = indata%bsubvmnc
p_to_bsubs = indata%bsubsmns
p_to_bsupu = indata%bsupumnc
p_to_bsupv = indata%bsupvmnc

outdata = vmec_fort8_data(mpol,ntor,ns,nfp,mnmax,itfsq,nit,gam, &
                   p_to_iotas,p_to_mass,p_to_pres, &
                   p_to_phips,p_to_bpco,p_to_bzco,p_to_phi,p_to_vp,p_to_jtheta, &
                   p_to_jzeta,p_to_specw, &
                   p_to_rmnc,p_to_zmns,p_to_lmns,p_to_bmod,p_to_gmod, &
                   p_to_bsubu,p_to_bsubv,p_to_bsubs,p_to_bsupu,p_to_bsupv)


end function extract_fort8_data_from_wout_data
!-------------------------------------------------------------------------------
function extract_fort8_data_from_wout_cdf_data(indata) result(outdata)
implicit none
type(vmec_wout_cdf_data), intent(in) :: indata
type(vmec_fort8_data) :: outdata
integer :: ns, mpol, ntor, nfp, mnmax, itfsq, nit
real(DP) :: gam
real(DP), dimension(:), pointer :: p_to_iotas => NULL()
real(DP), dimension(:), pointer :: p_to_mass => NULL()
real(DP), dimension(:), pointer :: p_to_pres => NULL()
real(DP), dimension(:), pointer :: p_to_phips => NULL()
real(DP), dimension(:), pointer :: p_to_bpco => NULL()
real(DP), dimension(:), pointer :: p_to_bzco => NULL()
real(DP), dimension(:), pointer :: p_to_phi => NULL()
real(DP), dimension(:), pointer :: p_to_vp => NULL()
real(DP), dimension(:), pointer :: p_to_jtheta => NULL()
real(DP), dimension(:), pointer :: p_to_jzeta => NULL()
real(DP), dimension(:), pointer :: p_to_specw => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_lmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bmod => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_gmod => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubu => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubv => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubs => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupu => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupv => NULL()
integer :: j, m, n, mn

if(indata%iasym == 1) then
 stop "No stellarator symmetry! fort8-data can not hold this information!"
endif

mpol  = indata%mpol
ntor  = indata%ntor
ns    = indata%ns
nfp   = indata%nfp
mnmax = indata%mnmax
itfsq = indata%itfsq
nit   = indata%niter/100 + 1
gam   = indata%gam

allocate(p_to_iotas(ns), p_to_mass(ns), p_to_pres(ns), p_to_phips(ns), &
         p_to_bpco(ns), p_to_bzco(ns), p_to_phi(ns), p_to_vp(ns), &
         p_to_jtheta(ns), p_to_jzeta(ns), p_to_specw(ns))
allocate(p_to_rmnc (-ntor:ntor,0:mpol-1,ns), p_to_zmns (-ntor:ntor,0:mpol-1,ns), &
         p_to_lmns (-ntor:ntor,0:mpol-1,ns), p_to_bmod (-ntor:ntor,0:mpol-1,ns), &
         p_to_gmod (-ntor:ntor,0:mpol-1,ns), p_to_bsubu(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsubv(-ntor:ntor,0:mpol-1,ns), p_to_bsubs(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsupu(-ntor:ntor,0:mpol-1,ns), p_to_bsupv(-ntor:ntor,0:mpol-1,ns))

! Initialize the data
p_to_iotas  = 0.0_DP
p_to_mass   = 0.0_DP
p_to_pres   = 0.0_DP
p_to_phips  = 0.0_DP
p_to_bpco   = 0.0_DP
p_to_bzco   = 0.0_DP
p_to_phi    = 0.0_DP
p_to_vp     = 0.0_DP
p_to_jtheta = 0.0_DP
p_to_jzeta  = 0.0_DP
p_to_specw  = 0.0_DP

p_to_rmnc = 0.0_DP
p_to_zmns = 0.0_DP
p_to_lmns = 0.0_DP
p_to_bmod = 0.0_DP
p_to_gmod = 0.0_DP
p_to_bsubu = 0.0_DP
p_to_bsubv = 0.0_DP
p_to_bsubs = 0.0_DP
p_to_bsupu = 0.0_DP
p_to_bsupv = 0.0_DP

! Transfer the data
p_to_iotas  = indata%iotas
p_to_mass   = indata%mass
p_to_pres   = indata%pres
p_to_phips  = indata%phips
p_to_bpco   = indata%buco
p_to_bzco   = indata%bvco
do j=2,ns   ! construct half-mesh phi
  p_to_phi(j) = (j-1.5_dp)/(ns-1.0_dp)
enddo
p_to_phi = p_to_phi*indata%phi(ns)
p_to_vp     = indata%vp
p_to_jtheta(2:ns) = 0.5*(indata%jcuru(2:ns)+indata%jcuru(1:ns-1))
p_to_jzeta(2:ns)  = 0.5*(indata%jcurv(2:ns)+indata%jcurv(1:ns-1))
p_to_specw  = indata%specw

! Mapping for rmnc, zmns and lmns straight forward.
do j=1,ns
  do mn=1,mnmax
    m = int(indata%xm(mn))
    n = int(indata%xn(mn))/nfp
    p_to_rmnc(n,m,j) = indata%rmnc(mn,j)
    p_to_zmns(n,m,j) = indata%zmns(mn,j)
    p_to_lmns(n,m,j) = indata%lmns(mn,j)
  enddo
! Rest of 2d-arrays has mnmax_nyq0 as first dimension.
! Cut out the parts used in prout: m < mpol and -ntor <= n <= ntor
!                                                 abs(n) <= ntor !
  do mn=1,indata%mnmax_nyq0
    m = int(indata%xm_nyq(mn))
    n = int(indata%xn_nyq(mn))/nfp
    if(m.lt.indata%mpol .and. abs(n).le.indata%ntor) then
      p_to_bmod(n,m,j)  = indata%bmnc(mn,j)
      p_to_bmod(n,m,j)  = indata%bmnc(mn,j)
      p_to_gmod(n,m,j)  = indata%gmnc(mn,j)
      p_to_bsubu(n,m,j) = indata%bsubumnc(mn,j)
      p_to_bsubv(n,m,j) = indata%bsubvmnc(mn,j)
      p_to_bsubs(n,m,j) = indata%bsubsmns(mn,j)
      p_to_bsupu(n,m,j) = indata%bsupumnc(mn,j)
      p_to_bsupv(n,m,j) = indata%bsupvmnc(mn,j)
    endif
  enddo
enddo

outdata = vmec_fort8_data(mpol,ntor,ns,nfp,mnmax,itfsq,nit,gam, &
                   p_to_iotas,p_to_mass,p_to_pres, &
                   p_to_phips,p_to_bpco,p_to_bzco,p_to_phi,p_to_vp,p_to_jtheta, &
                   p_to_jzeta,p_to_specw, &
                   p_to_rmnc,p_to_zmns,p_to_lmns,p_to_bmod,p_to_gmod, &
                   p_to_bsubu,p_to_bsubv,p_to_bsubs,p_to_bsupu,p_to_bsupv)

end function extract_fort8_data_from_wout_cdf_data

!-------------------------------------------------------------------------------
!   Interfaces to vmec_eqdata
!-------------------------------------------------------------------------------
function extract_vmec_eqdata_from_fort8_data(indata) result(outdata)
implicit none
type(vmec_fort8_data), intent(in) :: indata
type(vmec_eqdata) :: outdata
integer :: ns, mpol, ntor, nfp, mnmax, itfsq, nit, iasym
real(DP) :: gam
real(DP), dimension(:), pointer :: p_to_iotaf => NULL() !MJL 20150129
real(DP), dimension(:), pointer :: p_to_iotas => NULL()
real(DP), dimension(:), pointer :: p_to_mass => NULL()
real(DP), dimension(:), pointer :: p_to_pres => NULL()
real(DP), dimension(:), pointer :: p_to_phips => NULL()
real(DP), dimension(:), pointer :: p_to_buco => NULL()
real(DP), dimension(:), pointer :: p_to_bvco => NULL()
real(DP), dimension(:), pointer :: p_to_phi => NULL()
real(DP), dimension(:), pointer :: p_to_vp => NULL()
real(DP), dimension(:), pointer :: p_to_jcuru => NULL()
real(DP), dimension(:), pointer :: p_to_jcurv => NULL()
real(DP), dimension(:), pointer :: p_to_specw => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_lmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_lmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_gmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_gmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubumnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubumns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubvmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubvmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubsmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubsmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupumnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupumns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupvmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupvmns => NULL()
integer :: j

mpol  = indata%mpol
ntor  = indata%ntor
ns    = indata%ns
nfp   = indata%nfp
mnmax = indata%mnmax
itfsq = indata%itfsq
nit   = indata%nit/100 + 1
gam   = indata%gam
iasym = 0          ! old fort.8 data are always stellarator symmetric

allocate(p_to_iotas(ns), p_to_mass(ns), p_to_pres(ns), p_to_phips(ns), &
         p_to_buco(ns), p_to_bvco(ns), p_to_phi(ns), p_to_vp(ns), &
         p_to_jcuru(ns), p_to_jcurv(ns), p_to_specw(ns))
allocate(p_to_rmnc (-ntor:ntor,0:mpol-1,ns)   , p_to_zmns (-ntor:ntor,0:mpol-1,ns), &
         p_to_lmns (-ntor:ntor,0:mpol-1,ns)   , p_to_bmnc (-ntor:ntor,0:mpol-1,ns), &
         p_to_gmnc (-ntor:ntor,0:mpol-1,ns)   , p_to_bsubumnc(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsubvmnc(-ntor:ntor,0:mpol-1,ns), p_to_bsubsmns(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsupumnc(-ntor:ntor,0:mpol-1,ns), p_to_bsupvmnc(-ntor:ntor,0:mpol-1,ns))

! Initialize the data
p_to_iotas  = 0.0_DP
p_to_mass   = 0.0_DP
p_to_pres   = 0.0_DP
p_to_phips  = 0.0_DP
p_to_buco   = 0.0_DP
p_to_bvco   = 0.0_DP
p_to_phi    = 0.0_DP
p_to_vp     = 0.0_DP
p_to_jcuru  = 0.0_DP
p_to_jcurv  = 0.0_DP
p_to_specw  = 0.0_DP

p_to_rmnc = 0.0_DP
p_to_zmns = 0.0_DP
p_to_lmns = 0.0_DP
p_to_bmnc = 0.0_DP
p_to_gmnc = 0.0_DP
p_to_bsubumnc = 0.0_DP
p_to_bsubvmnc = 0.0_DP
p_to_bsubsmns = 0.0_DP
p_to_bsupumnc = 0.0_DP
p_to_bsupvmnc = 0.0_DP

! Transfer the data
p_to_iotas  = indata%iotas
p_to_mass   = indata%mass
p_to_pres   = indata%pres
p_to_phips  = indata%phips
p_to_buco   = indata%bpco
p_to_bvco   = indata%bzco
do j=2,ns   ! construct half-mesh phi
  p_to_phi(j) = (j-1.5_dp)/(ns-1.0_dp)
enddo
p_to_phi    = p_to_phi*indata%phi(ns)
p_to_vp     = indata%vp
p_to_jcuru  = indata%jtheta
p_to_jcurv  = indata%jzeta
p_to_specw  = indata%specw

p_to_rmnc     = indata%rmnc
p_to_zmns     = indata%zmns
p_to_lmns     = indata%lmns
p_to_bmnc     = indata%bmod
p_to_gmnc     = indata%gmod
p_to_bsubumnc = indata%bsubu
p_to_bsubvmnc = indata%bsubv
p_to_bsubsmns = indata%bsubs
p_to_bsupumnc = indata%bsupu
p_to_bsupvmnc = indata%bsupv

outdata = vmec_eqdata(mpol,ntor,ns,nfp,mnmax,itfsq,nit,iasym,gam, &
                   0_dp, p_to_iotaf, &  !MJL 20150129    fort8_data does not have Aminor_p
                   p_to_iotas,p_to_mass,p_to_pres, &
                   p_to_phips,p_to_buco,p_to_bvco,p_to_phi,p_to_vp, &
                   p_to_jcuru, p_to_jcurv, p_to_specw, &
                   p_to_rmnc,p_to_rmns,p_to_zmns,p_to_zmnc,p_to_lmns,p_to_lmnc, &
                   p_to_bmnc,p_to_bmns,p_to_gmnc,p_to_gmns, &
                   p_to_bsubumnc,p_to_bsubumns,p_to_bsubvmnc,p_to_bsubvmns, &
                   p_to_bsubsmnc,p_to_bsubsmns, &  ! Order of c and s corrected. MJL 20150505
                   p_to_bsupumnc,p_to_bsupumns,p_to_bsupvmnc,p_to_bsupvmns)

end function extract_vmec_eqdata_from_fort8_data
!-------------------------------------------------------------------------------
function extract_vmec_eqdata_from_wout_data(indata) result(outdata)
implicit none
type(vmec_wout_data), intent(in) :: indata
type(vmec_eqdata) :: outdata
integer :: ns, mpol, ntor, nfp, mnmax, itfsq, nit, iasym
real(DP) :: gam
real(DP) :: Aminor_p !MJL 20150129
real(DP), dimension(:), pointer :: p_to_iotaf => NULL() !MJL 20150192
real(DP), dimension(:), pointer :: p_to_iotas => NULL()
real(DP), dimension(:), pointer :: p_to_mass => NULL()
real(DP), dimension(:), pointer :: p_to_pres => NULL()
real(DP), dimension(:), pointer :: p_to_phips => NULL()
real(DP), dimension(:), pointer :: p_to_buco => NULL()
real(DP), dimension(:), pointer :: p_to_bvco => NULL()
real(DP), dimension(:), pointer :: p_to_phi => NULL()
real(DP), dimension(:), pointer :: p_to_vp => NULL()
real(DP), dimension(:), pointer :: p_to_jcuru => NULL()
real(DP), dimension(:), pointer :: p_to_jcurv => NULL()
real(DP), dimension(:), pointer :: p_to_specw => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_lmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_lmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_gmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_gmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubumnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubumns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubvmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubvmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubsmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubsmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupumnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupumns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupvmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupvmns => NULL()
integer :: j

mpol  = indata%mpol
ntor  = indata%ntor
ns    = indata%ns
nfp   = indata%nfp
mnmax = indata%mnmax
itfsq = indata%itfsq
nit   = indata%niter/100 + 1
gam   = indata%gam
Aminor_p   = indata%Aminor_p !MJL 20150129
iasym = indata%iasym

! Allocate and initialize the profile data
allocate(p_to_iotas(ns), p_to_mass(ns), p_to_pres(ns), p_to_phips(ns), &
         p_to_iotaf(ns), & !MJL 20150129
         p_to_buco(ns), p_to_bvco(ns), p_to_phi(ns), p_to_vp(ns), &
         p_to_jcuru(ns), p_to_jcurv(ns), p_to_specw(ns))
p_to_iotaf  = 0.0_DP !MJL 20150129
p_to_iotas  = 0.0_DP
p_to_mass   = 0.0_DP
p_to_pres   = 0.0_DP
p_to_phips  = 0.0_DP
p_to_buco   = 0.0_DP
p_to_bvco   = 0.0_DP
p_to_phi    = 0.0_DP
p_to_vp     = 0.0_DP
p_to_jcuru  = 0.0_DP
p_to_jcurv  = 0.0_DP
p_to_specw  = 0.0_DP

! Allocate and initialize the array data
allocate(p_to_rmnc (-ntor:ntor,0:mpol-1,ns)   , p_to_zmns (-ntor:ntor,0:mpol-1,ns), &
         p_to_lmns (-ntor:ntor,0:mpol-1,ns)   , p_to_bmnc (-ntor:ntor,0:mpol-1,ns), &
         p_to_gmnc (-ntor:ntor,0:mpol-1,ns)   , p_to_bsubumnc(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsubvmnc(-ntor:ntor,0:mpol-1,ns), p_to_bsubsmns(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsupumnc(-ntor:ntor,0:mpol-1,ns), p_to_bsupvmnc(-ntor:ntor,0:mpol-1,ns))
p_to_rmnc = 0.0_DP
p_to_zmns = 0.0_DP
p_to_lmns = 0.0_DP
p_to_bmnc = 0.0_DP
p_to_gmnc = 0.0_DP
p_to_bsubumnc = 0.0_DP
p_to_bsubvmnc = 0.0_DP
p_to_bsubsmns = 0.0_DP
p_to_bsupumnc = 0.0_DP
p_to_bsupvmnc = 0.0_DP
if (iasym .eq. 1)then
  allocate(p_to_rmns (-ntor:ntor,0:mpol-1,ns)   , p_to_zmnc (-ntor:ntor,0:mpol-1,ns), &
           p_to_lmnc (-ntor:ntor,0:mpol-1,ns)   , p_to_bmns (-ntor:ntor,0:mpol-1,ns), &
           p_to_gmns (-ntor:ntor,0:mpol-1,ns)   , p_to_bsubumns(-ntor:ntor,0:mpol-1,ns), &
           p_to_bsubvmns(-ntor:ntor,0:mpol-1,ns), p_to_bsubsmnc(-ntor:ntor,0:mpol-1,ns), &
           p_to_bsupumns(-ntor:ntor,0:mpol-1,ns), p_to_bsupvmns(-ntor:ntor,0:mpol-1,ns))
  p_to_rmns = 0.0_DP
  p_to_zmnc = 0.0_DP
  p_to_lmnc = 0.0_DP
  p_to_bmns = 0.0_DP
  p_to_gmns = 0.0_DP
  p_to_bsubumns = 0.0_DP
  p_to_bsubvmns = 0.0_DP
  p_to_bsubsmnc = 0.0_DP
  p_to_bsupumns = 0.0_DP
  p_to_bsupvmns = 0.0_DP
endif

! Transfer the data from half mesh values
p_to_iotaf  = indata%iotaf !MJL 20150129
p_to_iotas  = indata%iotas
p_to_mass   = indata%mass
p_to_pres   = indata%pres
p_to_phips  = indata%phip
p_to_buco   = indata%buco
p_to_bvco   = indata%bvco

! MJL 20150128  In the original code from Joachim Geiger, phi is converted from the full to half grid here.
! Let's instead keep phi on the full grid, as in VMEC.
!do j=2,ns   ! construct half-mesh phi
!  p_to_phi(j) = (j-1.5_dp)/(ns-1.0_dp)
!enddo
!p_to_phi = p_to_phi*indata%phi(ns)
p_to_phi   = indata%phi

p_to_vp     = indata%vp
! currently fort.8 structure implemented -> jcur.. on half grid
p_to_jcuru(2:ns) = 0.5*(indata%jcuru(2:ns)+indata%jcuru(1:ns-1))
p_to_jcurv(2:ns) = 0.5*(indata%jcurv(2:ns)+indata%jcurv(1:ns-1))
p_to_specw  = indata%specw

p_to_rmnc  = indata%rmnc
p_to_zmns  = indata%zmns
p_to_lmns  = indata%lmns
p_to_bmnc  = indata%bmnc
p_to_gmnc  = indata%gmnc
p_to_bsubumnc = indata%bsubumnc
p_to_bsubvmnc = indata%bsubvmnc
p_to_bsubsmns = indata%bsubsmns
p_to_bsupumnc = indata%bsupumnc
p_to_bsupvmnc = indata%bsupvmnc
if (iasym .eq. 1)then
  p_to_rmns  = indata%rmns
  p_to_zmnc  = indata%zmnc
  p_to_lmnc  = indata%lmnc
  p_to_bmns  = indata%bmns
  p_to_gmns  = indata%gmns
  p_to_bsubumns = indata%bsubumns
  p_to_bsubvmns = indata%bsubvmns
  p_to_bsubsmnc = indata%bsubsmnc
  p_to_bsupumns = indata%bsupumns
  p_to_bsupvmns = indata%bsupvmns
endif

outdata = vmec_eqdata(mpol,ntor,ns,nfp,mnmax,itfsq,nit,iasym,gam, &
                   Aminor_p, p_to_iotaf, &  !MJL 20150129
                   p_to_iotas,p_to_mass,p_to_pres, &
                   p_to_phips,p_to_buco,p_to_bvco,p_to_phi,p_to_vp, &
                   p_to_jcuru, p_to_jcurv, p_to_specw, &
                   p_to_rmnc,p_to_rmns,p_to_zmns,p_to_zmnc,p_to_lmns,p_to_lmnc, &
                   p_to_bmnc,p_to_bmns,p_to_gmnc,p_to_gmns, &
                   p_to_bsubumnc,p_to_bsubumns,p_to_bsubvmnc,p_to_bsubvmns, &
                   p_to_bsubsmnc,p_to_bsubsmns, &  ! Order of s and c corrected. MJL 20150505
                   p_to_bsupumnc,p_to_bsupumns,p_to_bsupvmnc,p_to_bsupvmns)

end function extract_vmec_eqdata_from_wout_data
!-------------------------------------------------------------------------------
function extract_vmec_eqdata_from_wout_cdf_data(indata) result(outdata)
implicit none
type(vmec_wout_cdf_data), intent(in) :: indata
type(vmec_eqdata) :: outdata
integer :: ns, mpol, ntor, nfp, mnmax, itfsq, nit, iasym
real(DP) :: gam
real(DP) :: Aminor_p  !MJL 20150129
real(DP), dimension(:), pointer :: p_to_iotaf => NULL() !MJL 20150129
real(DP), dimension(:), pointer :: p_to_iotas => NULL()
real(DP), dimension(:), pointer :: p_to_mass => NULL()
real(DP), dimension(:), pointer :: p_to_pres => NULL()
real(DP), dimension(:), pointer :: p_to_phips => NULL()
real(DP), dimension(:), pointer :: p_to_buco => NULL()
real(DP), dimension(:), pointer :: p_to_bvco => NULL()
real(DP), dimension(:), pointer :: p_to_phi => NULL()
real(DP), dimension(:), pointer :: p_to_vp => NULL()
real(DP), dimension(:), pointer :: p_to_jcuru => NULL()
real(DP), dimension(:), pointer :: p_to_jcurv => NULL()
real(DP), dimension(:), pointer :: p_to_specw => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_rmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_zmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_lmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_lmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_gmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_gmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubumnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubumns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubvmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubvmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubsmns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsubsmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupumnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupumns => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupvmnc => NULL()
real(DP), dimension(:,:,:), pointer :: p_to_bsupvmns => NULL()
integer :: j, m, n, mn
integer :: ialloc

mpol  = indata%mpol
ntor  = indata%ntor
ns    = indata%ns
nfp   = indata%nfp
mnmax = indata%mnmax
itfsq = indata%itfsq
nit   = indata%niter/100 + 1
iasym = indata%iasym
gam   = indata%gam
Aminor_p = indata%Aminor_p  !MJL 20150129

! Allocate and initialize the profile data
allocate(p_to_iotas(ns), p_to_mass(ns), p_to_pres(ns), p_to_phips(ns), &
         p_to_iotaf(ns), &  !MJL 20150129
         p_to_buco(ns), p_to_bvco(ns), p_to_phi(ns), p_to_vp(ns), &
         p_to_jcuru(ns), p_to_jcurv(ns), p_to_specw(ns),stat=ialloc)
if(ialloc .ne. 0) stop "allocation failure p_to_iotas in extract"
p_to_iotaf  = 0.0_DP !MJL 20150129
p_to_iotas  = 0.0_DP
p_to_mass   = 0.0_DP
p_to_pres   = 0.0_DP
p_to_phips  = 0.0_DP
p_to_buco   = 0.0_DP
p_to_bvco   = 0.0_DP
p_to_phi    = 0.0_DP
p_to_vp     = 0.0_DP
p_to_jcuru  = 0.0_DP
p_to_jcurv  = 0.0_DP
p_to_specw  = 0.0_DP

! Allocate and initialize the array data
allocate(p_to_rmnc (-ntor:ntor,0:mpol-1,ns)   , p_to_zmns (-ntor:ntor,0:mpol-1,ns), &
         p_to_lmns (-ntor:ntor,0:mpol-1,ns)   , p_to_bmnc (-ntor:ntor,0:mpol-1,ns), &
         p_to_gmnc (-ntor:ntor,0:mpol-1,ns)   , p_to_bsubumnc(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsubvmnc(-ntor:ntor,0:mpol-1,ns), p_to_bsubsmns(-ntor:ntor,0:mpol-1,ns), &
         p_to_bsupumnc(-ntor:ntor,0:mpol-1,ns), p_to_bsupvmnc(-ntor:ntor,0:mpol-1,ns), &
         stat=ialloc)
if(ialloc .ne. 0) stop "allocation failure p_to_rmnc in extract"
p_to_rmnc = 0.0_DP
p_to_zmns = 0.0_DP
p_to_lmns = 0.0_DP
p_to_bmnc = 0.0_DP
p_to_gmnc = 0.0_DP
p_to_bsubumnc = 0.0_DP
p_to_bsubvmnc = 0.0_DP
p_to_bsubsmns = 0.0_DP
p_to_bsupumnc = 0.0_DP
p_to_bsupvmnc = 0.0_DP

if (iasym .eq. 1)then
  allocate(p_to_rmns (-ntor:ntor,0:mpol-1,ns)   , p_to_zmnc (-ntor:ntor,0:mpol-1,ns), &
           p_to_lmnc (-ntor:ntor,0:mpol-1,ns)   , p_to_bmns (-ntor:ntor,0:mpol-1,ns), &
           p_to_gmns (-ntor:ntor,0:mpol-1,ns)   , p_to_bsubumns(-ntor:ntor,0:mpol-1,ns), &
           p_to_bsubvmns(-ntor:ntor,0:mpol-1,ns), p_to_bsubsmnc(-ntor:ntor,0:mpol-1,ns), &
           p_to_bsupumns(-ntor:ntor,0:mpol-1,ns), p_to_bsupvmns(-ntor:ntor,0:mpol-1,ns), &
         stat=ialloc)
if(ialloc .ne. 0) stop "allocation failure p_to_rmns in extract"
  p_to_rmns = 0.0_DP
  p_to_zmnc = 0.0_DP
  p_to_lmnc = 0.0_DP
  p_to_bmns = 0.0_DP
  p_to_gmns = 0.0_DP
  p_to_bsubumns = 0.0_DP
  p_to_bsubvmns = 0.0_DP
  p_to_bsubsmnc = 0.0_DP
  p_to_bsupumns = 0.0_DP
  p_to_bsupvmns = 0.0_DP
endif

! Transfer the data
p_to_iotaf  = indata%iotaf !MJL 20150129
p_to_iotas  = indata%iotas
p_to_mass   = indata%mass
p_to_pres   = indata%pres
p_to_phips  = indata%phips
p_to_buco   = indata%buco
p_to_bvco   = indata%bvco

! MJL 20150128  In the original code from Joachim Geiger, phi is converted from the full to half grid here.
! Let's instead keep phi on the full grid, as in VMEC.
!do j=2,ns   ! construct half-mesh phi
!  p_to_phi(j) = (j-1.5_dp)/(ns-1.0_dp)
!enddo
!p_to_phi = p_to_phi*indata%phi(ns)
p_to_phi   = indata%phi

p_to_vp     = indata%vp
! currently fort.8 structure implemented -> jcur.. on half grid
p_to_jcuru(2:ns) = 0.5*(indata%jcuru(2:ns)+indata%jcuru(1:ns-1))
p_to_jcurv(2:ns) = 0.5*(indata%jcurv(2:ns)+indata%jcurv(1:ns-1))
p_to_specw  = indata%specw

! Mapping for rmnc, zmns and lmns straight forward.
do j=1,ns
  do mn=1,mnmax
    m = int(indata%xm(mn))
    n = int(indata%xn(mn))/nfp
    p_to_rmnc(n,m,j) = indata%rmnc(mn,j)
    p_to_zmns(n,m,j) = indata%zmns(mn,j)
    p_to_lmns(n,m,j) = indata%lmns(mn,j)
  enddo
! Rest of 2d-arrays has mnmax_nyq0 as first dimension.
! Cut out the parts used in prout: m < mpol and -ntor <= n <= ntor
!                                                 abs(n) <= ntor !
  do mn=1,indata%mnmax_nyq0
    m = int(indata%xm_nyq(mn))
    n = int(indata%xn_nyq(mn))/nfp
    if(m.lt.indata%mpol .and. abs(n).le.indata%ntor) then
      p_to_bmnc(n,m,j)  = indata%bmnc(mn,j)
      p_to_gmnc(n,m,j)  = indata%gmnc(mn,j)
      p_to_bsubumnc(n,m,j) = indata%bsubumnc(mn,j)
      p_to_bsubvmnc(n,m,j) = indata%bsubvmnc(mn,j)
      p_to_bsubsmns(n,m,j) = indata%bsubsmns(mn,j)
      p_to_bsupumnc(n,m,j) = indata%bsupumnc(mn,j)
      p_to_bsupvmnc(n,m,j) = indata%bsupvmnc(mn,j)
    endif
  enddo
enddo
! Mapping for rmnc, zmns and lmns straight forward.
if(iasym .eq. 1)then
  do j=1,ns
    do mn=1,mnmax
      m = int(indata%xm(mn))
      n = int(indata%xn(mn))/nfp
      p_to_rmnc(n,m,j) = indata%rmnc(mn,j)
      p_to_zmns(n,m,j) = indata%zmns(mn,j)
      p_to_lmns(n,m,j) = indata%lmns(mn,j)
    enddo
! Rest of 2d-arrays has mnmax_nyq0 as first dimension.
! Cut out the parts used in prout: m < mpol and -ntor <= n <= ntor
!                                                 abs(n) <= ntor !
    do mn=1,indata%mnmax_nyq0
      m = int(indata%xm_nyq(mn))
      n = int(indata%xn_nyq(mn))/nfp
      if(m.lt.indata%mpol .and. abs(n).le.indata%ntor) then
        p_to_bmns(n,m,j)  = indata%bmns(mn,j)
        p_to_gmns(n,m,j)  = indata%gmns(mn,j)
        p_to_bsubumns(n,m,j) = indata%bsubumns(mn,j)
        p_to_bsubvmns(n,m,j) = indata%bsubvmns(mn,j)
        p_to_bsubsmnc(n,m,j) = indata%bsubsmnc(mn,j)
        p_to_bsupumns(n,m,j) = indata%bsupumns(mn,j)
        p_to_bsupvmns(n,m,j) = indata%bsupvmns(mn,j)
      endif
    enddo
  enddo
endif

outdata = vmec_eqdata(mpol,ntor,ns,nfp,mnmax,itfsq,nit,iasym,gam, &
                   Aminor_p, p_to_iotaf, &  !MJL 20150129
                   p_to_iotas,p_to_mass,p_to_pres, &
                   p_to_phips,p_to_buco,p_to_bvco,p_to_phi,p_to_vp, &
                   p_to_jcuru, p_to_jcurv, p_to_specw, &
                   p_to_rmnc,p_to_rmns,p_to_zmns,p_to_zmnc,p_to_lmns,p_to_lmnc, &
                   p_to_bmnc,p_to_bmns,p_to_gmnc,p_to_gmns, &
                   p_to_bsubumnc,p_to_bsubumns,p_to_bsubvmnc,p_to_bsubvmns, &
                   p_to_bsubsmnc,p_to_bsubsmns, & !Order of s and c corrected. MJL 20150505
                   p_to_bsupumnc,p_to_bsupumns,p_to_bsupvmnc,p_to_bsupvmns)

end function extract_vmec_eqdata_from_wout_cdf_data

!-------------------------------------------------------------------------------
!> the subroutine writes the formatted fort.8 file of vmec/nemec to the
!! file with the name <filename> with the data contained in the structure d.
subroutine write_fort8_data(filename, d)
use pisa_io
implicit none
character(len=*), intent(in) :: filename
type(vmec_fort8_data), intent(inout) :: d
integer :: iunit
character :: c_dummy
integer :: ios
integer :: j,m,n  !< loop variables
integer :: nmin, icount

iunit = get_next_io_unit()
open(unit=iunit, file=filename, action="write",status="new")

write(iunit,'(f10.3,7i6)',iostat=ios) &
             d%gam,d%nfp,d%ns,d%mpol,d%ntor,d%mnmax,d%itfsq,d%nit

do j=1,d%ns
  do m=0,d%mpol-1
    nmin=-d%ntor
    if(m == 0)nmin=0
    do n=nmin, d%ntor
      if(j == 1) write(iunit,'(2i10)')m , d%nfp*n
      write(iunit,'(5e20.13)') d%rmnc(n,m,j), d%zmns(n,m,j), d%lmns(n,m,j), &
                               d%bmod(n,m,j), d%gmod(n,m,j), &
                               d%bsubu(n,m,j),d%bsubv(n,m,j),d%bsubs(n,m,j), &
                               d%bsupu(n,m,j),d%bsupv(n,m,j)
    enddo
  enddo
enddo
write(iunit,'(5e20.13)')(d%iotas(j),d%mass(j),d%pres(j),d%phips(j), &
                         d%bpco(j),d%bzco(j),d%phi(j),d%vp(j), &
                         d%jtheta(j),d%jzeta(j),d%specw(j),j=2,d%ns)

close(iunit)

end subroutine write_fort8_data

!-------------------------------------------------------------------------------
subroutine read_vmec_indata_namelist (input_extension,input_data)
use pisa_io

type(vmec_indata) :: input_data
real(DP), parameter :: smallno = -HUGE(1._dp)
integer :: iunit, istat
integer :: nfp, ncurr, nsin, niter, nstep, nvacskip, mpol, ntor, &
           ntheta, nzeta &
         , ipmass, ipiota, ipcurr  !< added for profile control, J.Geiger
integer, dimension(nsdim) :: ns_array
integer :: imse, isnodes, itse, ipnodes, iopt_raxis, &
           imatch_phiedge, nflxs
integer, dimension(nbsetsp) :: nbfld
integer, dimension(nfloops) :: indxflx
integer, dimension(nbcoilsp,nbsetsp) :: indxbfld
real(DP), dimension(-ntord:ntord,0:mpol1d) :: rbs, zbc, rbc, zbs
real(DP) :: time_slice, curtor, delt, ftol, tcon0, &
               gamma, phiedge, phidiam, sigma_current, sigma_delphid, tensi, &
               tensp, tensi2, fpolyi, presfac, mseangle_offset, pres_offset, &
               mseangle_offsetm, spres_ped, bloat, pres_scale
real(DP), dimension(0:12) :: am, ai, ac, aphi  !< introduction of more coefficients, J.Geiger
real(DP), dimension(0:ntord) :: raxis, zaxis                !< !Backwards compatibility: Obsolete
real(DP), dimension(0:ntord) :: raxis_cc, raxis_cs, &
                                   zaxis_cc, zaxis_cs
real(DP), dimension(nsdim) :: ftol_array
real(DP), dimension(nigroup) :: extcur
real(DP), dimension(nmse) :: mseprof
real(DP), dimension(ntse) :: rthom, datathom, sigma_thom
real(DP), dimension(nmse) :: rstark, datastark, sigma_stark
real(DP), dimension(nfloops) :: dsiobt, sigma_flux
real(DP), dimension(nbcoilsp,nbsetsp) :: bbc, sigma_b
real(DP), dimension(ndatafmax) :: psa, pfa, isa, ifa
logical :: lpofr, lmac, lfreeb, lrecon, loldout, ledge_dump, lasym
logical :: laddout, ldiagno, lmoreiter    !< added by J.Geiger
logical :: lspectrum_dump, loptim           !< !Obsolete
character*(200) :: mgrid_file
character*(120) :: arg1
character*(132) :: input_extension
character*(120) :: filename

integer, parameter :: mpol_default = 6
integer, parameter :: ntor_default = 0
real(DP), parameter :: cbig   = 0.9e30_dp, zero = 0.0_DP, one = 1.0_DP

NAMELIST /indata/ mgrid_file, time_slice, nfp, ncurr, nsin, &
         niter, nstep, nvacskip, delt, ftol, gamma, am, ai, ac, aphi, &
         rbc, zbs, rbs, zbc, spres_ped, pres_scale, raxis_cc, zaxis_cs,  &
         raxis_cs, zaxis_cc, mpol, ntor, ntheta, nzeta,  &
         ns_array, ftol_array, tcon0, curtor, sigma_current, extcur, &
         phiedge, psa, pfa, isa, ifa, imatch_phiedge, iopt_raxis,  &
         tensi, tensp, mseangle_offset, mseangle_offsetm, imse,  &
         isnodes, rstark, datastark, sigma_stark, itse, ipnodes,  &
         presfac, pres_offset, rthom, datathom, sigma_thom, phidiam,  &
         sigma_delphid, tensi2, fpolyi, nflxs, indxflx, dsiobt,  &
         sigma_flux, nbfld, indxbfld, &
         bbc, sigma_b, lpofr, lfreeb, lrecon, lmac, loldout, lasym, &
         ipmass, ipiota, ipcurr, laddout, ldiagno, lmoreiter, & !< added by J.Geiger
         ledge_dump, lspectrum_dump, loptim, bloat, raxis, zaxis 

! Use initializations of vmec2000 - see vsetup-routine
lasym = .false.
mpol = mpol_default
ntor = ntor_default
ntheta = 0
nzeta  = 0
lfreeb = .true.              !reset later after input file READ
lrecon = .true.
loldout = .false.
laddout = .false.    !added to get txt and fort.8 in presence of netcdf
ldiagno = .false.    !added for diagno, J.Geiger
lmoreiter = .true.   !added to iteration control, J.Geiger
lmac = .false.
ledge_dump = .false.
ipmass = 0     ! profile selection 0=original polynomial, J.Geiger
ipiota = 0     ! profile selection 0=original polynomial, J.Geiger
ipcurr = 0     ! profile selection 0=original polynomial, J.Geiger

mgrid_file = 'NONE'
extcur = zero
rbc = zero ; rbs = zero ; zbc = zero ; zbs = zero

am = zero ; ai = zero ; ac = zero ; aphi = zero
bloat = one
pres_scale = one
ns_array = 0
ftol_array = zero
raxis_cc = zero;  raxis_cs = zero
zaxis_cc = zero;  zaxis_cs = zero
gamma = zero
spres_ped = one

delt = 1.1
tcon0 = one
curtor = 1.e30_dp
time_slice = zero
lpofr = .true.
imse = -1
itse = 0
isnodes = 0
ipnodes = 0
iopt_raxis = 1
imatch_phiedge = 1
nflxs = 0
nbfld = 0
mseangle_offset = zero
mseangle_offsetm = zero
pres_offset = zero
sigma_current = 1.e30_dp
sigma_delphid = 1.e30_dp
tensi = one
tensp = one
tensi2 = zero
fpolyi = one
presfac = one
phidiam = 1.e30_dp
mseprof = one
indxflx = 0
indxbfld = 0

sigma_stark = 1.1*cbig
sigma_thom = 1.1*cbig
sigma_flux = 1.1*cbig
sigma_b = 1.1*cbig
!
!     BACKWARDS COMPATIBILITY
!
raxis = smallno
zaxis = smallno
      
filename = 'input.'//trim(input_extension)
iunit=get_next_io_unit()
open(unit=iunit, file=filename, action='read')
read (iunit, nml=indata, iostat=istat)
close(unit=iunit)

where (raxis .gt. smallno) raxis_cc = raxis
where (zaxis .gt. smallno) zaxis_cs = zaxis
where (raxis .le. smallno) raxis = zero !don't keep this after using it
where (zaxis .le. smallno) zaxis = zero !don't keep this after using it

IF (ntheta .le. 0) ntheta = 2*mpol + 6
IF (nzeta .le. 0) nzeta = 2*ntor + 4

! Store everything in data structure
 input_data%nfp = nfp
 input_data%ncurr = ncurr
 input_data%nsin = nsin
 input_data%niter = niter
 input_data%nstep = nstep
 input_data%nvacskip = nvacskip
 input_data%mpol = mpol
 input_data%ntor = ntor
 input_data%ntheta = ntheta
 input_data%nzeta = nzeta
 input_data%ipmass = ipmass
 input_data%ipiota = ipiota
 input_data%ipcurr = ipcurr
 input_data%ns_array = ns_array
 input_data%imse = imse
 input_data%isnodes = isnodes
 input_data%itse = itse
 input_data%ipnodes = ipnodes
 input_data%iopt_raxis = iopt_raxis
 input_data%imatch_phiedge = imatch_phiedge
 input_data%nflxs = nflxs
 input_data%nbfld = nbfld
 input_data%indxflx = indxflx
 input_data%indxbfld = indxbfld
 input_data%rbs = rbs
 input_data%zbc = zbc
 input_data%rbc = rbc
 input_data%zbs = zbs
 input_data%time_slice = time_slice
 input_data%delt = delt
 input_data%ftol = ftol
 input_data%tcon0 = tcon0
 input_data%gamma = gamma
 input_data%phiedge = phiedge
 input_data%phidiam = phidiam
 input_data%sigma_current = sigma_current
 input_data%sigma_delphid = sigma_delphid
 input_data%tensi = tensi
 input_data%tensp = tensp
 input_data%tensi2 = tensi2
 input_data%fpolyi = fpolyi
 input_data%presfac = presfac
 input_data%mseangle_offset = mseangle_offset
 input_data%pres_offset = pres_offset
 input_data%mseangle_offsetm = mseangle_offsetm
 input_data%spres_ped = spres_ped
 input_data%bloat = bloat
 input_data%pres_scale = pres_scale
 input_data%am = am
 input_data%ai = ai
 input_data%ac = ac
 input_data%raxis = raxis
 input_data%zaxis = zaxis
 input_data%raxis_cc = raxis_cc
 input_data%raxis_cs = raxis_cs
 input_data%zaxis_cc = zaxis_cc
 input_data%zaxis_cs = zaxis_cs
 input_data%ftol_array = ftol_array
 input_data%extcur = extcur
 input_data%mseprof = mseprof
 input_data%rthom = rthom
 input_data%datathom = datathom
 input_data%sigma_thom = sigma_thom
 input_data%rstark = rstark
 input_data%datastark = datastark
 input_data%sigma_stark = sigma_stark
 input_data%dsiobt = dsiobt
 input_data%sigma_flux = sigma_flux
 input_data%bbc = bbc
 input_data%sigma_b = sigma_b
 input_data%psa = psa
 input_data%pfa = pfa
 input_data%isa = isa
 input_data%ifa = ifa
 input_data%lpofr = lpofr
 input_data%lmac = lmac
 input_data%lfreeb = lfreeb
 input_data%lrecon = lrecon
 input_data%loldout = loldout
 input_data%ledge_dump = ledge_dump
 input_data%lasym = lasym
 input_data%laddout = laddout
 input_data%ldiagno = ldiagno
 input_data%lmoreiter = lmoreiter
 input_data%lspectrum_dump = lspectrum_dump
 input_data%loptim = loptim
 input_data%mgrid_file = mgrid_file
 input_data%arg1 = arg1
 input_data%input_extension = input_extension

end subroutine read_vmec_indata_namelist

!-------------------------------------------------------------------------------
subroutine write_vmec_input_mcdb_tables(input_data)
use pisa_io
implicit none

type(vmec_indata) :: input_data
integer :: outunit
character(len=256) :: the_file
character(len=256) :: vmec2000_inputdata_database
character(len=256) :: vmec2000_input_am_database
character(len=256) :: vmec2000_input_ac_database
character(len=256) :: vmec2000_input_ai_database
character(len=256) :: vmec2000_input_extcur_database
character(len=256) :: vmec2000_input_ftol_database
character(len=256) :: vmec2000_input_nsarray_database
character(len=256) :: vmec2000_input_raxis_cc_database
character(len=256) :: vmec2000_input_raxis_cs_database
character(len=256) :: vmec2000_input_zaxis_cc_database
character(len=256) :: vmec2000_input_zaxis_cs_database
character(len=256) :: vmec2000_input_rbc_database
character(len=256) :: vmec2000_input_rbs_database
character(len=256) :: vmec2000_input_zbc_database
character(len=256) :: vmec2000_input_zbs_database
character(len=23) :: suffix1='vmec2000_input_data.csv'
character(len=21) :: suffix2='vmec2000_input_am.csv'
character(len=21) :: suffix3='vmec2000_input_ac.csv'
character(len=21) :: suffix4='vmec2000_input_ai.csv'
character(len=25) :: suffix5='vmec2000_input_extcur.csv'
character(len=23) :: suffix6='vmec2000_input_ftol.csv'
character(len=26) :: suffix7='vmec2000_input_nsarray.csv'
character(len=24) :: suffix8='vmec2000_input_raxis.csv'
character(len=24) :: suffix9='vmec2000_input_zaxis.csv'
character(len=27) :: suffix8a='vmec2000_input_raxis_cc.csv'
character(len=27) :: suffix9a='vmec2000_input_raxis_cs.csv'
character(len=27) :: suffix10a='vmec2000_input_zaxis_cc.csv'
character(len=27) :: suffix11a='vmec2000_input_zaxis_cs.csv'
character(len=22) :: suffix12='vmec2000_input_rbc.csv'
character(len=22) :: suffix13='vmec2000_input_rbs.csv'
character(len=22) :: suffix14='vmec2000_input_zbc.csv'
character(len=22) :: suffix15='vmec2000_input_zbs.csv'
integer :: i, nup, m, n

vmec2000_inputdata_database = trim(input_data%input_extension)//'_'//suffix1
the_file = trim(input_data%input_extension)//'_'//suffix1
outunit=get_next_io_unit()
open(unit=outunit, file=the_file, action='write', recl=1000)
!write(outunit,*)trim(input_data%input_extension),"," &
write(outunit,'(1p,(a),a,(a),a,l1,a,l1,a,l1,a,l1,a,2(e22.14,a) &
           &    ,6(i3,a),i8,a,i4,a,i2,5(a,e22.14))') &
                 trim(input_data%input_extension),"," &
                ,trim(input_data%mgrid_file),"," &
                ,input_data%lfreeb,"," &
                ,input_data%loldout,"," &
                ,input_data%laddout,"," &
                ,input_data%ldiagno,"," &
                ,input_data%delt,"," &
                ,input_data%tcon0,"," &
                ,input_data%nfp,"," &
                ,input_data%ncurr,"," &
                ,input_data%mpol,"," &
                ,input_data%ntor,"," &
                ,input_data%nzeta,"," &
                ,input_data%ntheta,"," &
                ,input_data%niter,"," &
                ,input_data%nstep,"," &
                ,input_data%nvacskip,"," &
                ,input_data%gamma,"," &
                ,input_data%phiedge,"," &
                ,input_data%bloat,"," &
                ,input_data%curtor,"," &
                ,input_data%spres_ped
close(unit=outunit)
the_file = trim(input_data%input_extension)//'_'//suffix2
open(unit=outunit, file=the_file, action='write', recl=1000)
write(outunit,'(1p,(a),a,i3,13(a,e22.14))',advance='no')trim(input_data%input_extension),"," &
               ,input_data%ipmass,(",",(input_data%am(i)),i=0,12)
write(outunit,'(a)')""
close(unit=outunit)

the_file = trim(input_data%input_extension)//'_'//suffix3
open(unit=outunit, file=the_file, action='write', recl=1000)
write(outunit,'(1p,(a),a,i3,13(a,e22.14))',advance='no')trim(input_data%input_extension),"," &
               ,input_data%ipcurr,(",",(input_data%ac(i)),i=0,12)
write(outunit,'(a)')""
close(unit=outunit)

the_file = trim(input_data%input_extension)//'_'//suffix4
open(unit=outunit, file=the_file, action='write', recl=1000)
write(outunit,'(1p,(a),a,i3,13(a,e22.14))',advance='no')trim(input_data%input_extension),"," &
               ,input_data%ipiota,(",",(input_data%ai(i)),i=0,12)
write(outunit,'(a)')""
close(unit=outunit)

nup=nigroup
do i=nigroup,1,-1
 if(input_data%extcur(i) == 0.0_DP) nup=i-1
enddo
if(input_data%input_extension(1:3) == "w7x") nup = max(nup,7)
the_file = trim(input_data%input_extension)//'_'//suffix5
open(unit=outunit, file=the_file, action='write', recl=1000)
!write(outunit,'(1p,(a),((a,e22.14)))',advance='no')trim(input_data%input_extension) &
write(outunit,*)trim(input_data%input_extension) &
               ,(",",(input_data%extcur(i)),i=1,nup)
close(unit=outunit)

nup=nsdim
do i=nsdim,1,-1
 if(input_data%ns_array(i) == 0.0_DP) nup=i-1
enddo
the_file = trim(input_data%input_extension)//'_'//suffix6
open(unit=outunit, file=the_file, action='write', recl=1000)
!write(outunit,'(1p,(a),a,(a,e22.14))',advance='no')trim(input_data%input_extension) &
!write(outunit,'(1p,(a),(a,e22.14))',advance="no")trim(input_data%input_extension) &
write(outunit,*)trim(input_data%input_extension) &
               ,(",",(input_data%ftol_array(i)),i=1,nup)
close(unit=outunit)

the_file = trim(input_data%input_extension)//'_'//suffix7
open(unit=outunit, file=the_file, action='write', recl=1000)
write(outunit,*)trim(input_data%input_extension) &
               ,(",",(input_data%ns_array(i)),i=1,nup)
close(unit=outunit)

nup=ntord
do i=ntord,0,-1
 if(input_data%raxis(i) == 0.0_DP) nup=i-1
enddo
the_file = trim(input_data%input_extension)//'_'//suffix8
open(unit=outunit, file=the_file, action='write', recl=1000)
write(outunit,*)trim(input_data%input_extension) &
               ,(",",(input_data%raxis(i)),i=0,nup)
close(unit=outunit)

nup=ntord
do i=ntord,0,-1
 if(input_data%zaxis(i) == 0.0_DP) nup=i-1
enddo
the_file = trim(input_data%input_extension)//'_'//suffix9
open(unit=outunit, file=the_file, action='write', recl=1000)
write(outunit,*)trim(input_data%input_extension) &
               ,(",",(input_data%zaxis(i)),i=0,nup)
close(unit=outunit)

!nup=ntord
!do i=ntord,0,-1
! if(input_data%raxis_cc(i) == 0.0_DP) nup=i-1
!enddo
!the_file = trim(input_data%input_extension)//'_'//suffix8a
!open(unit=outunit, file=the_file, action='write', recl=1000)
!write(outunit,*)trim(input_data%input_extension) &
!               ,(",",(input_data%raxis_cc(i)),i=0,nup)
!close(unit=outunit)
!
!nup=ntord
!do i=ntord,0,-1
! if(input_data%raxis_cs(i) == 0.0_DP) nup=i-1
!enddo
!the_file = trim(input_data%input_extension)//'_'//suffix9a
!open(unit=outunit, file=the_file, action='write', recl=1000)
!write(outunit,*)trim(input_data%input_extension) &
!               ,(",",(input_data%raxis_cs(i)),i=0,nup)
!close(unit=outunit)
!
!nup=ntord
!do i=ntord,0,-1
! if(input_data%zaxis_cc(i) == 0.0_DP) nup=i-1
!enddo
!the_file = trim(input_data%input_extension)//'_'//suffix10a
!open(unit=outunit, file=the_file, action='write', recl=1000)
!write(outunit,*)trim(input_data%input_extension) &
!               ,(",",(input_data%zaxis_cc(i)),i=0,nup)
!close(unit=outunit)
!
!nup=ntord
!do i=ntord,0,-1
! if(input_data%zaxis_cs(i) == 0.0_DP) nup=i-1
!enddo
!the_file = trim(input_data%input_extension)//'_'//suffix11a
!open(unit=outunit, file=the_file, action='write', recl=1000)
!write(outunit,*)trim(input_data%input_extension) &
!               ,(",",(input_data%zaxis_cs(i)),i=0,nup)
!close(unit=outunit)

the_file = trim(input_data%input_extension)//'_'//suffix12
open(unit=outunit, file=the_file, action='write')
write(outunit,'((a))',advance='no')trim(input_data%input_extension)
do m=0,mpol1d
 do n=-ntord,ntord
  if(input_data%rbc(n,m) /= 0.0_DP) then
    write(outunit,'(a,i3,a,i3,a,1pe22.14)',advance='no')"," &
               ,n,",",m,",",input_data%rbc(n,m)
  endif
 enddo
enddo
write(outunit,'(a)')""
close(unit=outunit)

the_file = trim(input_data%input_extension)//'_'//suffix13
open(unit=outunit, file=the_file, action='write')
write(outunit,'((a))',advance='no')trim(input_data%input_extension)
do m=0,mpol1d
 do n=-ntord,ntord
  if(input_data%rbs(n,m) /= 0.0_DP) then
    write(outunit,'(a,i3,a,i3,a,1pe22.14)',advance='no')"," &
               ,n,",",m,",",input_data%rbs(n,m)
  endif
 enddo
enddo
write(outunit,'(a)')""
close(unit=outunit)

the_file = trim(input_data%input_extension)//'_'//suffix14
open(unit=outunit, file=the_file, action='write')
write(outunit,'((a))',advance='no')trim(input_data%input_extension)
do m=0,mpol1d
 do n=-ntord,ntord
  if(input_data%zbc(n,m) /= 0.0_DP) then
    write(outunit,'(a,i3,a,i3,a,1pe22.14)',advance='no')"," &
               ,n,",",m,",",input_data%zbc(n,m)
  endif
 enddo
enddo
write(outunit,'(a)')""
close(unit=outunit)

the_file = trim(input_data%input_extension)//'_'//suffix15
open(unit=outunit, file=the_file, action='write')
write(outunit,'((a))',advance='no')trim(input_data%input_extension)
do m=0,mpol1d
 do n=-ntord,ntord
  if(input_data%zbs(n,m) /= 0.0_DP) then
    write(outunit,'(a,i3,a,i3,a,1pe22.14)',advance='no')"," &
               ,n,",",m,",",input_data%zbs(n,m)
  endif
 enddo
enddo
write(outunit,'(a)')""
close(unit=outunit)

end subroutine write_vmec_input_mcdb_tables

!-------------------------------------------------------------------------------
subroutine scale_vmec_wout_data(lin_scale, b_scale, wout_data)
implicit none
type(vmec_wout_data), intent(inout) :: wout_data
real(DP), intent(in) :: lin_scale  !> linear scaling of geometry
real(DP), intent(in) :: b_scale    !> scaling of field
real(DP) :: b2_scale
real(DP) :: lin2_scale
real(DP) :: lin3_scale

b2_scale = b_scale*b_scale
lin2_scale = lin_scale*lin_scale
lin3_scale = lin2_scale*lin_scale

wout_data%wb        = wout_data%wb*b2_scale
wout_data%wp        = wout_data%wp*b2_scale
wout_data%rmax_surf = wout_data%rmax_surf*lin_scale
wout_data%rmin_surf = wout_data%rmin_surf*lin_scale
wout_data%zmax_surf = wout_data%zmax_surf*lin_scale
wout_data%zmax_surf = wout_data%zmax_surf*lin_scale
wout_data%b0        = wout_data%b0*b2_scale
wout_data%IonLarmor = wout_data%IonLarmor/b_scale
wout_data%VolAvgB   = wout_data%VolAvgB*b_scale/lin3_scale
wout_data%rbtor     = wout_data%rbtor*b_scale
wout_data%ctor      = wout_data%ctor*b_scale
wout_data%Aminor_p  = wout_data%Aminor_p*lin_scale
wout_data%Rmajor_p  = wout_data%Rmajor_p*lin_scale
wout_data%volume_p  = wout_data%volume_p*lin3_scale

wout_data%extcur  = wout_data%extcur*lin_scale*b_scale
wout_data%presf  = wout_data%presf*b2_scale
wout_data%phipf  = wout_data%phipf*b_scale*lin2_scale
wout_data%phi  = wout_data%phi*b_scale*lin2_scale
wout_data%jcuru  = wout_data%jcuru*b_scale*lin_scale
wout_data%jcurv  = wout_data%jcurv*b_scale*lin_scale
!wout_data%mass  = wout_data%mass*b2_scale
wout_data%pres  = wout_data%pres*b2_scale
wout_data%phip  = wout_data%phip*b_scale*lin2_scale
wout_data%buco  = wout_data%buco*b_scale*lin_scale
wout_data%bvco  = wout_data%bvco*b_scale*lin_scale
wout_data%vp  = wout_data%vp*lin3_scale
wout_data%overr  = wout_data%overr/lin_scale

wout_data%rmnc  = wout_data%rmnc*lin_scale
wout_data%rmns  = wout_data%rmns*lin_scale
wout_data%zmns  = wout_data%zmns*lin_scale
wout_data%zmnc  = wout_data%zmnc*lin_scale
wout_data%bmnc  = wout_data%bmnc*b_scale
wout_data%bmns  = wout_data%bmns*b_scale
wout_data%gmnc  = wout_data%gmnc*lin3_scale
wout_data%gmns  = wout_data%gmns*lin3_scale
wout_data%bsubumnc  = wout_data%bsubumnc*b_scale
wout_data%bsubumns  = wout_data%bsubumns*b_scale
wout_data%bsubvmnc  = wout_data%bsubvmnc*b_scale
wout_data%bsubvmns  = wout_data%bsubvmns*b_scale
wout_data%bsubsmnc  = wout_data%bsubsmnc*b_scale
wout_data%bsubsmns  = wout_data%bsubsmns*b_scale
wout_data%bsupumnc  = wout_data%bsupumnc*b_scale
wout_data%bsupumns  = wout_data%bsupumns*b_scale
wout_data%bsupvmnc  = wout_data%bsupvmnc*b_scale
wout_data%bsupvmns  = wout_data%bsupvmns*b_scale
wout_data%currvmnc  = wout_data%currvmnc*b_scale*lin_scale

end subroutine scale_vmec_wout_data

end module pisa_vmec_module
