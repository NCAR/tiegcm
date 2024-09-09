!
      module oplus_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlat,nlonp4,dz,nlon,dlev,nlevp1,rp
      use magfield_module,only: bx,by,bz,bmod2 ! (nlonp4,-1:nlat+2)
      use addfld_module,only: addfld
      use diags_module,only: mkdiag_BXYZ,mkdiag_BMAG
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
      implicit none
!
      contains
!-----------------------------------------------------------------------
      subroutine oplus(tn,te,ti,o2,o1,he,n2,n2d,ne,u,v,w,ui,vi,wi,
     |  xnmbar,scht,op,optm1,opout,optm1out,xiop2p,xiop2d,Fe,Fn,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Update O+ ion at 3d task subdomain.
! Outputs are opout, optm1out, xiop2p, and xiop2d, all other args
!   are input.
! There are 4 latitude scans, with 3d mpi calls in between the loops.
!   (see also 3d gather/scatter calls in sub filter_op).
!
      use params_module,only: zpmid,spval  ! for O+ minimum
      use cons_module,only: rmass_op,gask,grav,re,cs,dphi,dlamda,
     |  shapiro,dtx2inv,rtd,rmassinv_o2,rmassinv_o1,
     |  rmassinv_n2,rmassinv_n2d,dtsmooth,dtsmooth_div2
      use qrj_module,only: 
     |  qop2p, ! O+(2p) ionization from qrj, used in xiop2p
     |  qop2d, ! O+(2d) ionization from qrj, used in xiop2d
     |  qop    ! O+ ionization from qrj
      use chemrates_module,only: ! needed chemical reaction rates
     |  rk1 ,rk2 ,rk10,rk16,rk17,rk18,rk19,rk20,
     |  rk21,rk22,rk23,rk24,rk25,rk26,rk27
      use input_module,only: nstep_sub,opfloor,oprate,oplev,oplatwidth
      use magfield_module,only: dipmag,sndec,csdec,rlatm
      use mpi_module,only: mp_bndlons_f3d, mp_periodic_f3d,
     |  mp_polelats_f3d,mp_geo_halos_f3d,mp_bndlats_f3d
!
! Args:
      integer,intent(in) :: 
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn, te, ti, ! neutral, electron, and ion temperatures (deg K)
     |  o2, o1,     ! o2, o mass mixing ratios
     |  he,         ! he mass mixing ratio
     |  n2d,        ! n2d
     |  ne,         ! electron density
     |  u,v,w,      ! neutral wind velocities (zonal, meridional, omega)
     |  optm1,      ! O+ at time n-1
     |  op,         ! O+ ion
     |  ui,vi,wi,   ! zonal, meridional, and vertical ion velocities
     |  xnmbar,     ! p0*e(-z)*mbar/kT
     |  scht        ! scale height
!
! N2 is intent(inout) for setting halo points (mp_geo_halos_f3d call below)
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(inout) :: n2
!
! Output fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  opout,    ! O+ output for next timestep
     |  optm1out, ! O+ output for time n-1
     |  xiop2p,xiop2d,Fe,Fn
!
! Local:
      integer :: k,i,lonbeg,lonend,lat,ier,nlevs
      integer :: jm2,jm1,jp1,jp2 ! lat-2, lat-1, lat+1, lat+2 
      real(rp),dimension(lon0:lon1,lat0:lat1) ::
     |  opflux,    ! upward number flux of O+ (returned by sub oplus_flux) (t7)
     |  dvb        ! output of sub divb
      real(rp),dimension(lon0:lon1) ::
     |  ubca, ubcb ! O+ upper boundary condition (were t2,t3)
      real(rp) :: explic = 1._rp, gmr
      real(rp),dimension(lev0:lev1,lon0:lon1) ::
     |  bdzdvb_op,        ! was s7
     |  explicit,         ! was s4
     |  hdz,              ! was s15
     |  tphdz1,tphdz0,    ! were s13,s12 (using gmr)
     |  djint,            ! was s11
     |  divbz,            ! was s7 (DIV(B)+(DH*D*BZ)/(D*BZ)
     |  hdzmbz,hdzpbz,    ! were s10,s9
     |  p_coeff,q_coeff,r_coeff, ! coefficients for tridiagonal solver (s1,s2,s3)
     |  tp1,              ! 0.5*(te+ti)
     |  wd,bdotdh_djbz,op_prod,op_loss_out,dopdt,
     |  diffsum,driftsum,windsum
      real(rp),dimension(lev0:lev1,lon0:lon1,lat0-1:lat1+1) :: hj ! (s10-s12)
!
! For O+ "smooth floor". The floor is a
! spatial gaussian based on the O+ minimum (opmin below).
!
! 1/21/16 btf: O+ minimum increased from 1000 to 3000.
! (This alleviates numerical instability that can develop
!  near the equator at the top of the model in some cases,
!  e.g., it allows the 2.5-deg June solstice solar max 
!  benchmark to run with step=30 rather than step=20, and
!  may improve stability in other cases as well.)
!
      real :: opmin
!
! Local fields at 3d subdomain (must be 3d to bridge latitude scans):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  bdotu,       ! was s7 (B.U)
     |  bvel,
     |  diffj,       ! (D/(H*DZ)*2.*TP+M*G/R)*N(O+) (s7,s8,s9)
     |  tp,          ! Plasma temperature (te+ti)
     |  tr,          ! Reduced temperature 0.5*(tn+ti)
     |  bdotdh_op,   ! (b(h)*del(h))*phi
     |  bdotdh_opj,  ! (b(h)*del(h))*phi
     |  bdotdh_diff, ! (b(h)*del(h))*phi
     |  dj,          ! diffusion coefficients (s13,s14,s15)
     |  optm1_smooth,! op at time n-1, with shapiro smoother (was s1)
     |  optm1_smooth2,
     |  op_loss,     ! was s13
     |  vni,djbz,vdotn_h,bdotdh_bvel,diffexp,
     |  diffp,diffq,diffr,driftp,driftq,driftr,windp,windq,windr
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,5) :: f5
      logical,parameter :: debug=.false. ! if set write print statements to stdout
!
! External:
      real,external :: fslt
!
      if (debug) write(6,"('Enter oplus.')")
#ifdef VT
!     code = 113 ; state = 'oplus' ; activity='ModelCode'
      call vtbegin(113,ier)
#endif
!
! Number of pressure levels (this will equal nlevp1):
      nlevs = lev1-lev0+1 ! for bndlons calls

!
! Make latitudes j=-1,0 from j=2,1, and j=nlat+1,nlat+2 from nlat-1,nlat for n2, 
! then define 2d halos points.  Note that mp_geo_halos_f3d and mp_polelats_f3d 
! do not define "outside corner" halo points, but this should be ok since 
! cross-derivatives are not performed.
!
      call mp_polelats_f3d(n2(:,lon0:lon1,:),   ! 3rd dim is lat0-2:lat1+2
     |  lev0,lev1,lon0,lon1,lat0,lat1,1,(/1._rp/)) ! last arg means no change in polesign
      call mp_geo_halos_f3d(n2,lev0,lev1,lon0,lon1,lat0,lat1,1)      
!
! Save inputs to secondary history file:
! 
!     do lat=lat0,lat1
!       call addfld('TE_OP',' ',' ',te(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('TI_OP',' ',' ',ti(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('N2D_OP',' ',' ',n2d(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('NE_OP',' ',' ',ne(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OPTM1',' ',' ',optm1(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OP_OPLUS',' ',' ',op(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('UI_OP',' ',' ',ui(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VI_OP',' ',' ',vi(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('WI_OP',' ',' ',wi(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! Sub oplus_flux returns upward number flux of O+ in 
!
      call oplus_flux(opflux,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('oplus after oplus_flux.')")
!     call addfld('OPFLUX',' ',' ',opflux,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Divergence is returned in dvb(lon0:lon1,lat0:lat1) by sub divb:
      call divb(dvb,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('oplus after divb.')")
!     call addfld('OP_DIVB',' ',' ',dvb,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!
      call mkdiag_BXYZ('BX',bx(lon0:lon1,lat0:lat1),lon0,lon1,lat0,lat1)
      call mkdiag_BXYZ('BY',by(lon0:lon1,lat0:lat1),lon0,lon1,lat0,lat1)
      call mkdiag_BXYZ('BZ',bz(lon0:lon1,lat0:lat1),lon0,lon1,lat0,lat1)
      call mkdiag_BMAG('BMAG',bmod2(lon0:lon1,lat0:lat1),
     |  lon0,lon1,lat0,lat1)
!
!----------------------- Begin first latitude scan ---------------------
      do lat=lat0,lat1
        if (debug) write(6,"('oplus begin first lat scan: lat=',i3)")lat
        jm2 = lat-2
        jm1 = lat-1
        jp1 = lat+1
        jp2 = lat+2
!
! Set reduced temperature (average of tn and ti)
!
! 1/2/16 btf: 
!   This was previously set incorrectly as tp = 0.5*(te+ti)) 
!   (see also plasma temp tp). Correcting this resulted in a
!   40% change in ambipolar diffusion DJ, but only 10% or so
!   change in O+ and Ne, and very small differences in NMF2,HMF2,TEC.
!   Wenbin feels there should be larger changes in these fields
!   as a consequence of fixing this bug. 
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          tr(k,i,jm1) = 0.5_rp*(tn(k,i,jm1)+ti(k,i,jm1))
          tr(k,i,lat) = 0.5_rp*(tn(k,i,lat)+ti(k,i,lat))
          tr(k,i,jp1) = 0.5_rp*(tn(k,i,jp1)+ti(k,i,jp1))
        enddo
      enddo
      tr(lev1,:,:) = 0._rp
!     call addfld('TR','TR: reduced temperature (0.5*(tn+ti))',' ',
!    |  tr(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! rrk returns djm1,dj,djp1:
!
      call rrk(xnmbar(:,lon0:lon1,jm1),
     |  o2(:,lon0:lon1,jm1),o1(:,lon0:lon1,jm1),
     |  he(:,lon0:lon1,jm1),n2(:,lon0:lon1,jm1),
     |  tr(:,lon0:lon1,jm1),dj(:,lon0:lon1,jm1),
     |  vni(:,lon0:lon1,jm1),lon0,lon1,lev0,lev1)

      call rrk(xnmbar(:,lon0:lon1,lat),
     |  o2(:,lon0:lon1,lat),o1(:,lon0:lon1,lat),
     |  he(:,lon0:lon1,lat),n2(:,lon0:lon1,lat),
     |  tr(:,lon0:lon1,lat),dj(:,lon0:lon1,lat),
     |  vni(:,lon0:lon1,lat),lon0,lon1,lev0,lev1)

      call rrk(xnmbar(:,lon0:lon1,jp1),
     |  o2(:,lon0:lon1,jp1),o1(:,lon0:lon1,jp1),
     |  he(:,lon0:lon1,jp1),n2(:,lon0:lon1,jp1),
     |  tr(:,lon0:lon1,jp1),dj(:,lon0:lon1,jp1),
     |  vni(:,lon0:lon1,jp1),lon0,lon1,lev0,lev1)

      if (debug) write(6,"('oplus after rrk: lat=',i3)") lat
!     call addfld('DJ','DJ: Ambipolar Diffusion of O+',' ',
!    |    dj(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Plasma temperature:
      do i=lon0,lon1
        do k=lev0,lev1-1
          tp(k,i,jm1) = te(k,i,jm1)+ti(k,i,jm1)
          tp(k,i,lat) = te(k,i,lat)+ti(k,i,lat)
          tp(k,i,jp1) = te(k,i,jp1)+ti(k,i,jp1)
        enddo
      enddo
      if (debug) write(6,"('oplus after tpj: lat=',i3)") lat
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          hj(k,i,jm1) = scht(k,i,jm1)
          hj(k,i,lat) = scht(k,i,lat)
          hj(k,i,jp1) = scht(k,i,jp1)
        enddo
      enddo
      if (debug) write(6,"('oplus after hj: lat=',i3)") lat
!     call addfld('HJ  ',' ',' ',hj(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! bdotu = B.U (s7)
      do i=lon0,lon1
        do k=lev0,lev1-1
          bdotu(k,i,jm1) = bx(i,jm1)*u(k,i,jm1)+by(i,jm1)*v(k,i,jm1)+
     |      hj(k,i,jm1)*bz(i,jm1)*0.5*(w(k,i,jm1)+w(k+1,i,jm1))
          bdotu(k,i,lat) = bx(i,lat)*u(k,i,lat)+by(i,lat)*v(k,i,lat)+
     |      hj(k,i,lat)*bz(i,lat)*0.5*(w(k,i,lat)+w(k+1,i,lat))
          bdotu(k,i,jp1) = bx(i,jp1)*u(k,i,jp1)+by(i,jp1)*v(k,i,jp1)+
     |      hj(k,i,jp1)*bz(i,jp1)*0.5*(w(k,i,jp1)+w(k+1,i,jp1))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('BDOTU' ,' ',' ',bdotu(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! bvel @ jm1 = (B.U)*N(O+)    (J-1) (was s6)
! bvel @ j   = (B.U)*N(O+)      (J) (was s7)
! bvel @ jp1 = (B.U)*N(O+)    (J+1) (was s8)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          bvel(k,i,jm1) = bdotu(k,i,jm1)*op(k,i,jm1)
          bvel(k,i,lat) = bdotu(k,i,lat)*op(k,i,lat)
          bvel(k,i,jp1) = bdotu(k,i,jp1)*op(k,i,jp1)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
      if (debug) write(6,"('oplus after bvel: lat=',i3)") lat
!     call addfld('BVEL_J'  ,' ',' ',bvel(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      tp(lev1,:,jm1:jp1) = 0._rp
      call diffus(
     |  tp(:,lon0:lon1,jm1),op(:,lon0:lon1,jm1),hj(:,:,jm1),
     |  diffj(:,lon0:lon1,jm1),lon0,lon1,lev0,lev1)
      call diffus(
     |  tp(:,lon0:lon1,lat),op(:,lon0:lon1,lat),hj(:,:,lat),
     |  diffj(:,lon0:lon1,lat),lon0,lon1,lev0,lev1)
      call diffus(
     |  tp(:,lon0:lon1,jp1),op(:,lon0:lon1,jp1),hj(:,:,jp1),
     |  diffj(:,lon0:lon1,jp1),lon0,lon1,lev0,lev1)

      if (debug) write(6,"('oplus after diffus: lat=',i3)") lat
!     call addfld('DIFFJ','DIFFJ after diffus',' ',
!    |  diffj(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Plasma temperature times O+
      do i=lon0,lon1
        do k=lev0,lev1-1
          tp(k,i,jm2) = op(k,i,jm2)*(te(k,i,jm2)+ti(k,i,jm2))
          tp(k,i,jm1) = tp(k,i,jm1)*op(k,i,jm1)
          tp(k,i,lat) = tp(k,i,lat)*op(k,i,lat)
          tp(k,i,jp1) = tp(k,i,jp1)*op(k,i,jp1)
          tp(k,i,jp2) = op(k,i,jp2)*(te(k,i,jp2)+ti(k,i,jp2))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
      if (debug) write(6,"('oplus after tpj2: lat=',i3)") lat
!     call addfld('TPJ','TP after times OP',' ',
!    |  tp(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Shapiro smoother: optm1 is O+ at time n-1 (optm1_smooth was s1)
! optm1_smooth will be used in explicit terms below.
      do i=lon0,lon1
        do k=lev0,lev1-1
          optm1_smooth(k,i,lat) = optm1(k,i,lat)-shapiro/nstep_sub*
     |      (optm1(k,i,jp2)+optm1(k,i,jm2)-4._rp*
     |      (optm1(k,i,jp1)+optm1(k,i,jm1))+6._rp*
     |       optm1(k,i,lat))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('OPTM1_SM0' ,' ',' ',
!    |  optm1_smooth(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
        if (debug) write(6,"('oplus end first lat scan: lat=',i3)") lat
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            wd(k,i) = vni(k,i,lat)*diffj(k,i,lat)*dj(k,i,lat)*
     |        sin(dipmag(i,lat))*cos(dipmag(i,lat))
            Fe(k,i,lat) = wd(k,i)*sndec(i,lat)
            Fn(k,i,lat) = wd(k,i)*csdec(i,lat)
          enddo
        enddo
      enddo

      do lat=lat0-1,lat1+1
        do i=lon0,lon1
          do k=lev0,lev1-1
            djbz(k,i,lat) = dj(k,i,lat)*bz(i,lat)
          enddo
        enddo
      enddo
!------------------------- End first latitude scan ---------------------
!
! Boundary longitudes:
      f5(:,:,:,1) = djbz(:,:,:)
      f5(:,:,:,2) = bvel(:,:,:)
      f5(:,:,:,3) = diffj(:,:,:)
      f5(:,:,:,4) = tp(:,:,:)
      f5(:,:,:,5) = optm1_smooth(:,:,:)

      call mp_bndlons_f3d(f5,nlevs,lon0,lon1,lat0,lat1,5,0)

      djbz(:,:,:)         = f5(:,:,:,1)
      bvel(:,:,:)         = f5(:,:,:,2)
      diffj(:,:,:)        = f5(:,:,:,3)
      tp(:,:,:)          = f5(:,:,:,4)
      optm1_smooth(:,:,:) = f5(:,:,:,5)
!
!----------------------- Begin second latitude scan --------------------
      do lat=lat0,lat1
        if (debug) 
     |    write(6,"('oplus begin second lat scan: lat=',i3)") lat
        jm2 = lat-2
        jm1 = lat-1
        jp1 = lat+1
        jp2 = lat+2
!
! bdotdh_op = (B(H).DEL(H))*(D/(H*DZ)*TP+M*G/R)*N(O+)
! then bdotdh_op = d*bz*bdotdh_op
!
        call bdotdh(
     |    diffj(:,lon0:lon1,jm1),
     |    diffj(:,:,lat),
     |    diffj(:,lon0:lon1,jp1),
     |    bdotdh_op(:,lon0:lon1,lat),lon0,lon1,lev0,lev1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          bdotdh_op(k,i,lat) = djbz(k,i,lat)*bdotdh_op(k,i,lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('BDOTDH_1',' ',' ',
!    |  bdotdh_op(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! bdotdh_opjm1 = (B(H).DEL(H))*2.*TP*N(O+)    (J-1)
! bdotdh_opj   = (B(H).DEL(H))*2.*TP*N(O+)      (J)
! bdotdh_opjp1 = (B(H).DEL(H))*2.*TP*N(O+)    (J+1)
!
      call bdotdh(
     |  tp(:,lon0:lon1,jm2),tp(:,:,jm1),tp(:,lon0:lon1,lat),
     |  bdotdh_opj(:,lon0:lon1,jm1),lon0,lon1,lev0,lev1,jm1)
      call bdotdh(
     |  tp(:,lon0:lon1,jm1),tp(:,:,lat),tp(:,lon0:lon1,jp1),
     |  bdotdh_opj(:,lon0:lon1,lat),lon0,lon1,lev0,lev1,lat)
      call bdotdh(
     |  tp(:,lon0:lon1,lat),tp(:,:,jp1),tp(:,lon0:lon1,jp2),
     |  bdotdh_opj(:,lon0:lon1,jp1),lon0,lon1,lev0,lev1,jp1)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          bdotdh_opj(k,i,jm1) = bdotdh_opj(k,i,jm1)*dj(k,i,jm1)
          bdotdh_opj(k,i,lat) = bdotdh_opj(k,i,lat)*dj(k,i,lat)
          bdotdh_opj(k,i,jp1) = bdotdh_opj(k,i,jp1)*dj(k,i,jp1)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('BDOTDH_2',' ',' ',
!    |  bdotdh_op(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
        if (debug) 
     |    write(6,"('oplus end second lat scan: lat=',i3)") lat
      enddo ! lat=lat0,lat1
!------------------------ End second latitude scan ---------------------
!
! Periodic points for bdotdh_opj (output from bdotdh above):
      call mp_periodic_f3d(bdotdh_opj(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
!
! Boundary longitudes for bdotdh_opj (input to below call to bdotdh):
      call mp_bndlons_f3d(bdotdh_opj,nlevs,lon0,lon1,lat0,lat1,1,0)
!
!----------------------- Begin third latitude scan ---------------------
      do lat=lat0,lat1
        if (debug) 
     |    write(6,"('oplus begin third lat scan: lat=',i3)") lat
        jm2 = lat-2
        jm1 = lat-1
        jp1 = lat+1
        jp2 = lat+2
!
! bdotdh_opj = (B(H).DEL(H))*D*(B(H).DEL(H))*2.*TP*N(O+)   (J)
! Note bdotdh_opj longitude dimension is lon-2:lon+2. bdotdh_diff 
!   is returned. (periodic points apparently not necessary for
!   bdotdh_diff)
!
      call bdotdh(
     |  bdotdh_opj(:,lon0:lon1,jm1),
     |  bdotdh_opj(:,:,lat),
     |  bdotdh_opj(:,lon0:lon1,jp1),
     |  bdotdh_diff(:,lon0:lon1,lat),lon0,lon1,lev0,lev1,lat)

!     call addfld('BDOT_J'  ,' ',' ',bdotdh_opj(lev0:lev1-1,lon0:lon1,
!    |  lat),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('BDOT_DIF',' ',' ',bdotdh_diff(lev0:lev1-1,lon0:lon1,
!    |  lat),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! bdzdvb_op = (BZ*D/(H*DZ)+DIV(*B))*S2
! bdzdvb returns bdzdvb_op.
!
      call bdzdvb(bdotdh_opj(:,lon0:lon1,lat),dvb(:,lat),hj(:,:,lat),
     |  bdzdvb_op,lev0,lev1,lon0,lon1,lat)

!     call addfld('BDZDVB',' ',' ',
!    |  bdzdvb_op(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Collect explicit terms:
      do i=lon0,lon1
        do k=lev0,lev1-1
          diffexp(k,i,lat) =
     |      bdzdvb_op(k,i)+bdotdh_diff(k,i,lat)+bdotdh_op(k,i,lat)
          explicit(k,i) = -explic*diffexp(k,i,lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('EXPLIC0',' ',' ',explicit(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('BX_OP',' ',' ',bx(lon0:lon1,:),
!    |  'lon',lon0,lon1,'lat',lat,lat,0)
!     call addfld('BY_OP',' ',' ',by(lon0:lon1,:),
!    |  'lon',lon0,lon1,'lat',lat,lat,0)
!     call addfld('UI_VEL',' ',' ',ui(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('VI_VEL',' ',' ',vi(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
 
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2
!
      call bdotdh(
     |  bvel(:,lon0:lon1,jm1),bvel(:,:,lat),bvel(:,lon0:lon1,jp1),
     |  bdotdh_bvel(:,lon0:lon1,lat),lon0,lon1,lev0,lev1,lat)
!
! Note if input flag DYNAMO<=0, then ui,vi,wi velocities will be zero.
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          vdotn_h(k,i,lat) = 1._rp/(2._rp*re)*
     |      (1._rp/(cs(lat)*dlamda)*
     |      (0.5_rp*(ui(k,i,lat)+ui(k+1,i,lat))*bmod2(i,lat)**2*
     |      (op(k,i+1,lat)/bmod2(i+1,lat)**2-
     |       op(k,i-1,lat)/bmod2(i-1,lat)**2))+
!
     |      1._rp/dphi*
     |      (0.5_rp*(vi(k,i,lat)+vi(k+1,i,lat))*bmod2(i,lat)**2*
     |      (op(k,i,jp1)/bmod2(i,jp1)**2-
     |       op(k,i,jm1)/bmod2(i,jm1)**2)))
          explicit(k,i) =
     |      explicit(k,i)+vdotn_h(k,i,lat)+bdotdh_bvel(k,i,lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0+2,lon1-2
!
! Periodic points for explicit terms.
! This is apparently unnecessary:
!     call periodic_f2d(explicit,lon0,lon1,nlevs)

!     call addfld('EXPLIC1',' ',' ',explicit(lev0:lev1-1,:),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      do i=lon0,lon1
        dvb(i,lat) = dvb(i,lat)/bz(i,lat)
      enddo ! i=lon0,lon1
!      
      do i=lon0,lon1
        do k=lev0,lev1-1
          hdz(k,i) = 1._rp/(hj(k,i,lat)*dz)
          tp1(k,i) = 0.5_rp*(ti(k,i,lat)+te(k,i,lat))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('TP1',' ',' ',tp1(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('HDZ',' ',' ',hdz(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! gmr = G*M(O+)/(2.*R)
      gmr = grav*rmass_op/(2._rp*gask)
      do i=lon0,lon1
        do k=lev0,lev1-2
          tphdz1(k+1,i) = 2.*tp1(k+1,i)*(0.5_rp*
     |      (hdz(k,i)+hdz(k+1,i)))+
     |      gmr ! s13
          tphdz0(k+1,i) = 2.*tp1(k  ,i)*(0.5_rp*
     |      (hdz(k,i)+hdz(k+1,i)))-
     |      gmr ! s12
        enddo ! k=lev0,lev1-2
      enddo ! i=lon0,lon1
!
! Upper and lower boundaries:
      do i=lon0,lon1
        tphdz1(lev0,i) = 2._rp*tp1(lev0,i)*
     |                   (1.5_rp*hdz(lev0,i)-
     |                   0.5_rp*hdz(lev0+1,i))+gmr
        tphdz1(lev1,i) = 2._rp*(2._rp*tp1(lev1-1,i)-
     |                   tp1(lev1-2,i))*
     |                   (1.5_rp*hdz(lev1-1,i)-
     |                   0.5_rp*hdz(lev1-2,i))+gmr
        tphdz0(lev0,i) = 2._rp*(2._rp*tp1(lev0,i)-
     |                   tp1(lev0+1,i))*
     |                   (1.5_rp*hdz(lev0,i)-
     |                   0.5_rp*hdz(lev0+1,i))-gmr
        tphdz0(lev1,i) = 2._rp*tp1(lev1-1,i)*
     |                   (1.5_rp*hdz(lev1-1,i)-
     |                   0.5_rp*hdz(lev1-2,i))-gmr
      enddo ! i=lon0,lon1

!     call addfld('TPHDZ1',' ',' ',tphdz1,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TPHDZ0',' ',' ',tphdz0,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! djint = dj at interfaces:
      do i=lon0,lon1
        do k=lev0,lev1-2
          djint(k+1,i) = 0.5_rp*(dj(k,i,lat)+dj(k+1,i,lat))
        enddo ! k=lev0,lev1-2
        djint(lev0,i) = (1.5_rp*dj(lev0  ,i,lat)-
     |                  0.5_rp*dj(lev0+1,i,lat))
        djint(lev1,i) = (1.5_rp*dj(lev1-1,i,lat)-
     |                  0.5_rp*dj(lev1-2,i,lat))
      enddo ! i=lon0,lon1
!     call addfld('DJINT' ,' ',' ',djint,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      call bdotdh(
     |  djbz(:,lon0:lon1,jm1),djbz(:,:,lat),djbz(:,lon0:lon1,jp1),
     |  bdotdh_djbz,lon0,lon1,lev0,lev1,lat)
!
! divbz = (DIV(B)+(DH*D*BZ)/(D*BZ) (was s7)
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          divbz(k,i) = 
     |      dvb(i,lat)+bdotdh_djbz(k,i)/(dj(k,i,lat)*bz(i,lat)**2)
        enddo ! k=lev0,lev1-1
      enddo ! i=lonbeg,lonend

! Periodic points for divbz apparently not necessary:
!     call periodic_f2d(divbz,lon0,lon1,nlevs)
!
! Set periodic points to zero to avoid NaNS trap:
      if (lon0==1) divbz(:,lon0:lon0+1) = 0._rp
      if (lon1==nlonp4) divbz(:,lon1-1:lon1) = 0._rp
!     call addfld('DIVBZ' ,' ',' ',divbz(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! hdzmbz = (1./(H*DZ)-(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2 (was s10)
! hdzpbz = (1./(H*DZ)+(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2 (was s9 )
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          hdzmbz(k,i) = (hdz(k,i)-0.5_rp*divbz(k,i))*bz(i,lat)**2
          hdzpbz(k,i) = (hdz(k,i)+0.5_rp*divbz(k,i))*bz(i,lat)**2
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('HDZMBZ' ,' ',' ',hdzmbz(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('HDZPBZ' ,' ',' ',hdzpbz(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Sum O+ at time n-1 to explicit terms: N(O+)/(2*DT) (N-1) (was s4)
! Boundary longitudes for optm1_smooth were obtained after first
! latitude scan above.
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          optm1_smooth2(k,i,lat) =
     |      optm1_smooth(k,i,lat)-shapiro/nstep_sub*
     |      (optm1_smooth(k,i+2,lat)+optm1_smooth(k,i-2,lat)-4._rp*
     |      (optm1_smooth(k,i+1,lat)+optm1_smooth(k,i-1,lat))+6._rp*
     |       optm1_smooth(k,i,lat))
        enddo
      enddo
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          explicit(k,i) = explicit(k,i)-
     |      optm1_smooth2(k,i,lat)*dtx2inv*nstep_sub
        enddo ! k=lev0,lev1
      enddo ! i=lonbeg,lonend

!     call addfld('OPTM1_SM1' ,' ',' ',
!    |  optm1_smooth(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('EXPLIC2' ,' ',' ',explicit(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
! 
! Begin coefficients p_coeff, q_coeff, r_coeff  (s1,s2,s3)
      do i=lon0,lon1
        do k=lev0,lev1-1
          diffp(k,i,lat) =   hdzmbz(k,i)*djint(k  ,i)*tphdz0(k  ,i)
          diffq(k,i,lat) = -(hdzpbz(k,i)*djint(k+1,i)*tphdz0(k+1,i)+
     |                       hdzmbz(k,i)*djint(k  ,i)*tphdz1(k  ,i))
          diffr(k,i,lat) =   hdzpbz(k,i)*djint(k+1,i)*tphdz1(k+1,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('P_COEFF0',' ',' ',diffp(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEFF0',' ',' ',diffq(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEFF0',' ',' ',diffr(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Continue coefficients with vertical ion velocity:
      do i=lon0,lon1
        do k=lev0,lev1-2
          driftp(k+1,i,lat) =
     |      0.5_rp*(wi(k+1,i,lat)+wi(k+2,i,lat))*0.5*hdz(k+1,i)
          driftq(k,i,lat) =
     |      -0.5_rp*(wi(k,i,lat)+wi(k+1,i,lat))*6._rp/re
          driftr(k,i,lat) =
     |      -0.5_rp*(wi(k,i,lat)+wi(k+1,i,lat))*0.5_rp*hdz(k,i)
        enddo ! k=lev0,lev1-2
      enddo ! i=lon0,lon1
!
! Upper and lower boundaries:
      do i=lon0,lon1
        driftp(lev0,i,lat) = 
     |    0.5_rp*(wi(lev0,i,lat)+wi(lev0+1,i,lat))*
     |    0.5_rp*hdz(lev0,i)
        driftq(lev1-1,i,lat) =
     |    -0.5_rp*(wi(lev1,i,lat)+wi(lev1-1,i,lat))*
     |    6._rp/re
        driftr(lev1-1,i,lat) =
     |    -0.5_rp*(wi(lev1,i,lat)+wi(lev1-1,i,lat))*
     |    0.5_rp*hdz(lev1-1,i)
      enddo ! i=lon0,lon1

!     call addfld('P_COEFF1',' ',' ',driftp(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEFF1',' ',' ',driftq(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEFF1',' ',' ',driftr(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Continue coefficients with vertical neutral wind:
      do i=lon0,lon1
        do k=lev0,lev1-2
          windp(k+1,i,lat) = bz(i,lat)*bdotu(k,i,lat)*
     |         0.5_rp*hdz(k+1,i)
          windq(k,i,lat) = -bdotu(k,i,lat)*dvb(i,lat)*bz(i,lat)
          windr(k,i,lat) = -bz(i,lat)*bdotu(k+1,i,lat)*
     |         0.5_rp*hdz(k,i)
        enddo ! k=lev0,lev1-2
      enddo ! i=lon0,lon1
!
! Upper and lower boundaries:
      do i=lon0,lon1
        windp(lev0,i,lat) = bz(i,lat)*
     |    (2._rp*bdotu(lev0,i,lat)-
     |    bdotu(lev0+1,i,lat))*0.5_rp*hdz(lev0,i)
        windq(lev1-1,i,lat) = -bdotu(lev1-1,i,lat)*dvb(i,lat)*bz(i,lat)
        windr(lev1-1,i,lat) = bz(i,lat)*
     |   (2._rp*bdotu(lev1-1,i,lat)-bdotu(lev1-2,i,lat))*
     |    0.5_rp*hdz(lev1-1,i)
      enddo ! i=lon0,lon1

!     call addfld('P_COEFF2',' ',' ',windp(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEFF2',' ',' ',windq(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEFF2',' ',' ',windr(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Additions to Q coefficients:
      do i=lon0,lon1
        do k=lev0,lev1-1
          p_coeff(k,i) = diffp(k,i,lat)+driftp(k,i,lat)+windp(k,i,lat)
          q_coeff(k,i) = diffq(k,i,lat)+driftq(k,i,lat)+windq(k,i,lat)
     |      -dtx2inv*nstep_sub
          r_coeff(k,i) = diffr(k,i,lat)+driftr(k,i,lat)+windr(k,i,lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! Upper boundary condition for O+:
      do i=lon0,lon1
        ubca(i) = 0._rp
        ubcb(i) = -bz(i,lat)**2*djint(lev1,i)*tphdz0(lev1,i)-ubca(i) ! t3
        ubca(i) = -bz(i,lat)**2*djint(lev1,i)*tphdz1(lev1,i)+ubca(i) ! t2
!
! Q = Q+B/A*R
        q_coeff(lev1-1,i) = q_coeff(lev1-1,i)+ubcb(i)/ubca(i)*
     |    r_coeff(lev1-1,i) 
!
! F = F -R/A*PHI
        explicit(lev1-1,i) = explicit(lev1-1,i)-opflux(i,lat)*
     |    r_coeff(lev1-1,i)/ubca(i)
        r_coeff(lev1-1,i) = 0._rp
      enddo ! i=lon0,lon1

!     call addfld('EXPLIC3',' ',' ',explicit(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('P_COEFF2',' ',' ',p_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEFF2',' ',' ',q_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEFF2',' ',' ',r_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

!     call addfld('QOP2P_OP',' ',' ',qop2p(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP2D_OP',' ',' ',qop2d(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RK20',' ',' ',rk20(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('RK25',' ',' ',rk25(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Sources and sinks (xiop2p and xiop2d are outputs):
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xiop2p(k,i,lat) = 
     |      0.5_rp*(qop2p(k,i,lat)+qop2p(k+1,i,lat))/((xnmbar(k,i,lat)*
     |      ((rk16+rk17)*n2(k,i,lat)*rmassinv_n2+
     |      rk18*o1(k,i,lat)*rmassinv_o1))+
     |      (rk19(k,i,lat)+rk20(k,i,lat))*ne(k,i,lat)+rk21+rk22)
! 
          xiop2d(k,i,lat) = 
     |      (0.5_rp*(qop2d(k,i,lat)+qop2d(k+1,i,lat))+(rk20(k,i,lat)*
     |      ne(k,i,lat)+rk22)*xiop2p(k,i,lat))/((xnmbar(k,i,lat)*
     |      (rk23*n2(k,i,lat)*rmassinv_n2+rk24*
     |      o1(k,i,lat)*rmassinv_o1+rk26*o2(k,i,lat)*rmassinv_o2))+
     |      rk25(k,i,lat)*ne(k,i,lat)+rk27)
! 
          op_loss(k,i,lat) = 
     |      xnmbar(k,i,lat)*(rk1(k,i,lat)*o2(k,i,lat)*rmassinv_o2+
     |      rk2(k,i,lat)*n2(k,i,lat)*rmassinv_n2+rk10*
     |      n2d(k,i,lat)*rmassinv_n2d)
!
          q_coeff(k,i) = q_coeff(k,i)-op_loss(k,i,lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('XIOP2P',' ',' ',xiop2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('XIOP2D',' ',' ',xiop2d(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_LOSS',' ',' ',op_loss(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('OP_QOP',' ',' ',qop(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_NE' ,' ',' ',ne(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_O1' ,' ',' ',o1(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_TN' ,' ',' ',tn(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_RK19' ,' ',' ',rk19(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('OP_RK25' ,' ',' ',rk25(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Add source term to RHS (explicit terms):
      do i=lon0,lon1
        do k=lev0,lev1-1
          op_prod(k,i) =
     |      0.5_rp*(qop(k,i,lat)+qop(k+1,i,lat))+(rk19(k,i,lat)*
     |      ne(k,i,lat)+rk21)*
     |      xiop2p(k,i,lat)+(rk25(k,i,lat)*ne(k,i,lat)+rk27)*
     |      xiop2d(k,i,lat)+(rk18*xiop2p(k,i,lat)+rk24*xiop2d(k,i,lat))*
     |      o1(k,i,lat)*rmassinv_o1*xnmbar(k,i,lat)
          explicit(k,i) = explicit(k,i)-op_prod(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('EXPLIC4',' ',' ',explicit(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Lower boundary condition N(O+) = Q/L:
      do i=lon0,lon1
        q_coeff(lev0,i) = q_coeff(lev0,i)-p_coeff(lev0,i) 
        explicit(lev0,i) = explicit(lev0,i)-2._rp*p_coeff(lev0,i)*
     |    qop(lev0,i,lat)/
     |    (1.5_rp*op_loss(lev0,i,lat)-0.5_rp*op_loss(lev0+1,i,lat))
        p_coeff(lev0,i) = 0._rp
      enddo ! i=lon0,lon1

!     call addfld('P_COEFF',' ',' ',p_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEFF',' ',' ',q_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEFF',' ',' ',r_coeff(lev0:lev1-2,:),
!    |  'lev',lev0,lev1-2,'lon',lon0,lon1,lat)
!     call addfld('EXPLIC5',' ',' ',explicit(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Tridiagonal solver returns updated O+ in opout (all other args are input):
!     subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1,lonmax,lat,
!    |  idebug)
!
      call trsolv(p_coeff,q_coeff,r_coeff,explicit,
     |  opout(:,lon0:lon1,lat),lev0,lev1,lev0,lev1-1,lon0,lon1,nlonp4,
     |  lat,0)
!
!     call addfld('OP_SOLV',' ',' ',opout(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

        if (debug) 
     |    write(6,"('oplus end third lat scan: lat=',i3)") lat
      enddo ! lat=lat0,lat1
!------------------------ End third latitude scan ---------------------
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            dopdt(k,i) = dtx2inv*nstep_sub*
     |        (opout(k,i,lat)-optm1_smooth2(k,i,lat))
            op_loss_out(k,i) = op_loss(k,i,lat)*opout(k,i,lat)
          enddo
          do k=lev0+1,lev1-1
            diffsum (k,i) = diffexp(k,i,lat)+
     |        diffp (k,i,lat)*opout(k-1,i,lat)+
     |        diffq (k,i,lat)*opout(k  ,i,lat)+
     |        diffr (k,i,lat)*opout(k+1,i,lat)
            driftsum(k,i) = -vdotn_h(k,i,lat)+
     |        driftp(k,i,lat)*opout(k-1,i,lat)+
     |        driftq(k,i,lat)*opout(k  ,i,lat)+
     |        driftr(k,i,lat)*opout(k+1,i,lat)
            windsum (k,i) = -bdotdh_bvel(k,i,lat)+
     |        windp (k,i,lat)*opout(k-1,i,lat)+
     |        windq (k,i,lat)*opout(k  ,i,lat)+
     |        windr (k,i,lat)*opout(k+1,i,lat)
          enddo
        enddo
        call addfld('DOPDT','O+ changing rate','cm-3s-1',
     |    dopdt(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('OP_PROD','O+ production','cm-3s-1',
     |    op_prod(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('OP_LOSS','O+ loss','cm-3s-1',
     |    -op_loss_out(lev0:lev1-1,:),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('OP_DIFF',
     |    'O+ transport due to ambipolar diffusion','cm-3s-1',
     |    diffsum(lev0+1:lev1-1,:),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('OP_DRIFT',
     |    'O+ transport due to ion drift','cm-3s-1',
     |     driftsum(lev0+1:lev1-1,:),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('OP_WIND',
     |    'O+ transport due to neutral wind','cm-3s-1',
     |    windsum(lev0+1:lev1-1,:),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
      enddo
!
! Filter updated O+:
! 
      call filter_op(opout(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,'OPLUS')
!
!----------------------- Begin fourth latitude scan ---------------------
      do lat=lat0,lat1
        if (debug) 
     |    write(6,"('oplus begin fourth lat scan: lat=',i3)") lat
!
!     call addfld('OP_FILT',' ',' ',opout(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Time smoothing:
!
! optm1out(k,i,lat): New O+ at current latitude and time n-1.
! op(k,i,lat)      : O+ at current latitude and time.
! optm1(k,i,lat)   : O+ at current latitude and time n-1.
! opout(k,i,lat)   : New O+ at current latitude and time.
!
      do i=lon0,lon1
        do k=lev0,lev1-1      
          optm1out(k,i,lat) = dtsmooth*op(k,i,lat)+dtsmooth_div2*
     |      (optm1(k,i,lat)+opout(k,i,lat))
        enddo ! k=lev0,lev1-1      
      enddo ! i=lon0,lon1
!
! Upper boundary:
      opout(lev1,:,lat) = 0._rp
      optm1out(lev1,:,lat) = 0._rp
!
! Insure global non-negative O+:
      do i=lon0,lon1
        do k=lev0,lev1-1      
          if (opout(k,i,lat)    < 0._rp) opout(k,i,lat)    = 0._rp
          if (optm1out(k,i,lat) < 0._rp) optm1out(k,i,lat) = 0._rp
        enddo ! k=lev0,lev1-1      
      enddo ! i=lon0,lon1
!
! 12/4/14 btf: Enforce O+ minimum.
! Opfloor is Stan's "smooth floor" (product of two Gaussians, 
!   dependent on latitude and pressure level) (opmin=1000.0):
!
! 2024/08 Haonan Wu:
!   Change the altitude distribution of O+ minimum from Gaussian to logistic,
!   also change the latitude distribution based on magnetic latitudes
!   instead of geographic latitudes
!
      if (opfloor/=0. .and. oprate/=0. .and.
     |    oplev/=spval .and. oplatwidth/=0.) then
        do k=lev0,lev1-1
          do i=lon0,lon1
            opmin = opfloor*exp(-(rlatm(i,lat)*rtd/oplatwidth)**2)
     |        /(1+exp(-oprate*(zpmid(k)-oplev)))
            if (opout(k,i,lat) < opmin) then
              opout(k,i,lat) = opmin
            endif ! opout < opmin
          enddo ! i=lon0,lon1
        enddo ! k=lev0,lev1-1
      endif
!
! End fourth and final latitude scan:
      enddo ! lat=lat0,lat1
!
! Periodic points for outputs:
      call mp_periodic_f3d(opout(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
      call mp_periodic_f3d(optm1out(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
!
      call mp_polelats_f3d(opout(:,lon0:lon1,:),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1,(/1./))
      call mp_bndlats_f3d(opout,nlevs,lon0,lon1,lat0,lat1,1)
      call mp_bndlons_f3d(opout,nlevs,lon0,lon1,lat0,lat1,1,0)
!
      call mp_polelats_f3d(optm1out(:,lon0:lon1,:),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1,(/1./))
      call mp_bndlats_f3d(optm1out,nlevs,lon0,lon1,lat0,lat1,1)
      call mp_bndlons_f3d(optm1out,nlevs,lon0,lon1,lat0,lat1,1,0)
!
! Save outputs on secondary history for diagnostics:
!     do lat=lat0,lat1
!       call addfld('OPOUT',' ',' ',opout(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('OPOUTM1',' ',' ',optm1out(lev0:lev1-1,lon0:lon1,
!    |    lat),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
#ifdef VT
!     code = 113 ; state = 'oplus' ; activity='ModelCode'
      call vtend(113,ier)
#endif
      if (debug) 
     |  write(6,"('oplus returning.')")
      end subroutine oplus
!-----------------------------------------------------------------------
      subroutine oplus_flux(opflux,lon0,lon1,lat0,lat1)
!
! Calculate O+ number flux in opflux for sub oplus (was sub opflux).
!
      use cons_module,only: pi,rtd
      use chapman_module,only: chi     ! was t2 in old sub opflux
      use magfield_module,only: rlatm
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
      real,intent(out) :: opflux(lon0:lon1,lat0:lat1)
!
! Local:
      integer :: i,lat
      real,parameter ::
     |  phid =  2.0e8_rp,
     |  phin = -2.0e8_rp,
     |  ppolar = 0._rp
      real :: a(lon0:lon1)   ! was t3 ("a" needs a better name)
      real :: fed(lon0:lon1) ! was t4
      real :: fen(lon0:lon1) ! was t5
!
! Latitude scan:
      do lat=lat0,lat1
!
! Longitude loop:
        do i=lon0,lon1
!
! Change upper boundary flux as electron heat flux in settei
          if (abs(rlatm(i,lat)) >= pi/4.5_rp) then
            a(i) = 1._rp
          elseif (abs(rlatm(i,lat)) <= pi/18._rp) then
            a(i) = 0._rp
          else
            a(i) = .5_rp*(1._rp+cos(abs(rlatm(i,lat))*
     |             6._rp-pi*4._rp/3._rp))
          endif
          fed(i) = phid*a(i)
          fen(i) = phin*a(i)
          if (chi(i,lat)-0.5_rp*pi >= 0._rp) then
            opflux(i,lat) = fen(i)
          else
            opflux(i,lat) = fed(i)
          endif
          if ((chi(i,lat)*rtd-80._rp)*(chi(i,lat)*rtd-100._rp) < 0._rp)
     |      opflux(i,lat) = .5_rp*(fed(i)+fen(i))+
     |      .5_rp*(fed(i)-fen(i))*
     |        cos(pi*(chi(i,lat)*rtd-80._rp)/20._rp)
!
! Add ppolar if magnetic latitude >= 60 degrees:
          if (abs(rlatm(i,lat))-pi/3._rp >= 0._rp)
     |      opflux(i,lat) = opflux(i,lat)+ppolar
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
      end subroutine oplus_flux
!-----------------------------------------------------------------------
      subroutine divb(dvb,lon0,lon1,lat0,lat1)
      use cons_module,only: cs,dlamda,dphi,re  
!
! Evaluate divergence of B, the unit magnetic field vector.
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1 
      real,intent(out) :: dvb(lon0:lon1,lat0:lat1)
!
! Local:
      integer :: lonbeg,lonend,i,lat,jm1,jp1
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2
!
      do lat=lat0,lat1
        jm1 = lat-1
        jp1 = lat+1
        dvb(:,lat) = 0._rp
        do i=lonbeg,lonend
          dvb(i,lat) = (((bx(i+1,lat)-bx(i-1,lat))/(2._rp*dlamda)+
     |      (cs(jp1)*by(i,jp1)-cs(jm1)*by(i,jm1))/(2._rp*dphi))/
     |      cs(lat)+2._rp*bz(i,lat))/re
        enddo ! i=lonbeg,lonend
      enddo ! lat=lat0,lat1
!
      end subroutine divb
!-----------------------------------------------------------------------
      subroutine rrk(rms,ps1,ps2,he,n2,tr,ans,vni,
     |  lon0,lon1,lev0,lev1)
!
! Returns ambipolar diffusion coefficient in ans and vni.
!
      use params_module,only: zpmid,spval
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,
     |  rmassinv_he
      use input_module,only: colfac,opdiffcap,opdiffrate,opdifflev
      use init_module,only: istep
!
! Args:
      integer,intent(in) :: lon0,lon1,lev0,lev1
      real,dimension(lev0:lev1,lon0:lon1),intent(in) ::
     |  rms,ps1,ps2,he,n2,tr
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: ans,vni
!
! Local:
      integer :: k,i
      real :: opdiffmax
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! 8/28/13 btf: Use n2=1-o2-o-he, and include O+/He collision rate from Wenbin:
!
          vni(k,i) = 18.1_rp*ps1(k,i)*rmassinv_o2+
     |      ps2(k,i)*rmassinv_o1*sqrt(tr(k,i))*
     |      (1._rp-0.064_rp*alog10(tr(k,i)))**2*colfac+
     |      3.6_rp*he(k,i)*rmassinv_he+18.6_rp*n2(k,i)*rmassinv_n2
          ans(k,i) = 1.42E17_rp/(rms(k,i)*vni(k,i))
          vni(k,i) = 16._rp*3.53E-11_rp*vni(k,i)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
! 1/30/16 btf:
! Cap ambipolar diffusion coefficient. Namelist parameter OPDIFFCAP.
! This was tested with various values (1.5e8, 3e8, 6e8, 8e8),
! and was found to improve numerical stability in some storm cases, 
! for example the November, 2003 and July, 2000 storms, with Weimer
! potential model and IMF/OMNI data. Both of these cases generally 
! will not complete if timestep is longer than 10 sec. The Nov, 2003
! may succeed with opdiffcap turned off, but the July, 2000 
! "Bastille day storm" will succeed only with step=10 and opdiffcap=6.e8.
!
! 2024/08 Haonan Wu:
!   Change the altitude distribution of maximum O+ diffusion flux
!   from Gaussian to logistic based on input parameters
!   this may help stabilize the model in some storm cases
!   which needs a steeper cap during the main phase
!   but also needs a smoother cap during the recovery phase
!
      if (opdiffcap /= 0. .and. opdiffrate /= 0 .and.
     |    opdifflev /= spval) then ! default is off
!        if (istep==1) write(6,"('oplus rrk: opdiffcap = ',es12.4)")
!     |    opdiffcap
!       where(ans(:,:) > opdiffcap)
!         ans(:,:) = opdiffcap
!       endwhere
        do k=lev0,lev1-1
          opdiffmax = opdiffcap/
     |      (1+exp(opdiffrate*(zpmid(k)-opdifflev)))
          do i=lon0,lon1
            if (ans(k,i) > opdiffmax) ans(k,i) = opdiffmax
          enddo
        enddo
      endif

      end subroutine rrk
!-----------------------------------------------------------------------
      subroutine diffus(tp,en,hj,ans,lon0,lon1,lev0,lev1)
!
! Evaluates ans = (d/(h*dz)*tp+m*g/r)*en
!
      use cons_module,only: rmass_op,grav,gask
!
! Args:
      integer :: lon0,lon1,lev0,lev1
      real,dimension(lev0:lev1,lon0:lon1),intent(in) :: tp,en,hj
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: ans
!
! Local:
      integer :: k,i,k1,k2
      real :: mgr
!
      mgr = rmass_op*grav/gask
      do i=lon0,lon1
        do k=lev0,lev1-2
          k1 = k+1
          k2 = k+2
          ans(k1,i) = 1._rp/(2._rp*
     |      hj(k1,i)*dlev)*(tp(k2,i)*en(k2,i)-
     |      tp(k,i)*en(k,i))+mgr*en(k1,i)
        enddo
      enddo
!
! Upper and lower boundaries:
      k1 = lev1-1
      k2 = lev1-2
      do i=lon0,lon1
        ans(k1,i) = 1._rp/(hj(k1,i)*dlev)*(tp(k1,i)*en(k1,i)-
     |    tp(k2,i)*en(k2,i))+mgr*en(k1,i)
        ans(1,i) = 1._rp/(hj(1,i)*dlev)*(tp(2,i)*en(2,i)-
     |    tp(1,i)*en(1,i))+mgr*en(1,i)
      enddo
      end subroutine diffus
!-----------------------------------------------------------------------
      subroutine bdotdh(phijm1,phij,phijp1,ans,lon0,lon1,lev0,lev1,lat)
      use cons_module,only: re,dphi,dlamda,cs
!
! Evaluates ans = (b(h)*del(h))*phi
!
! Args:
      integer :: lon0,lon1,lev0,lev1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in)    :: phijm1,phijp1
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: phij
      real,dimension(lev0:lev1,lon0:lon1),intent(out)   :: ans
!
! Local:
      integer :: k,i,lonbeg,lonend
!
      lonbeg = lon0
      if (lon0==1) then
        lonbeg = 3
        ans(:,lon0:lon1) = 0._rp ! set periodic points to zero to avoid NaNS trap
      endif
      lonend = lon1
      if (lon1==nlonp4) then
        lonend = lon1-2
        ans(:,lon1-1:lon1) = 0._rp ! set periodic points to zero to avoid NaNS trap
      endif
!
! Note phij longitude dimension is lon0-2:lon1+2 (only i-1 and i+1 are used).
! Boundary longitudes i-1 and i+1 must have been set before this routine is
! called (e.g., call mp_bndlons_f3d).
!
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          ans(k,i) = 1._rp/re*(bx(i,lat)/(cs(lat)*2._rp*dlamda)*
     |      (phij(k,i+1)-phij(k,i-1))+by(i,lat)*
     |      (phijp1(k,i)-phijm1(k,i))/(2._rp*dphi))
        enddo ! k=lev0,lev1
      enddo ! i=lonbeg,lonend
!
      end subroutine bdotdh
!-----------------------------------------------------------------------
      subroutine bdzdvb(phi,dvb,h,ans,lev0,lev1,lon0,lon1,lat)
!
! Evaluates  ans = (bz*d/(h*dz)+divb)*phi
!
! Args:
      integer :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: dvb(lon0:lon1)
      real,dimension(lev0:lev1,lon0:lon1),intent(in)    :: phi,h
      real,dimension(lev0:lev1,lon0:lon1),intent(out)   :: ans
!
! Local:
      integer :: k,i
!
      do i=lon0,lon1
        do k=lev0+1,lev1-2
          ans(k,i) = bz(i,lat)/
     |      (2._rp*h(k,i)*dz)*(phi(k+1,i)-phi(k-1,i))+
     |      dvb(i)*phi(k,i)
        enddo ! k=lev0+1,lev1-1
      enddo ! i=lon0,lon1
!
! Upper and lower boundaries:
      do i=lon0,lon1
        ans(lev1-1,i) = bz(i,lat)/(h(lev1-1,i)*dz)*(phi(lev1-1,i)-
     |    phi(lev1-2,i))+dvb(i)*phi(lev1-1,i)
        ans(lev0,i) = bz(i,lat)/(h(lev0,i)*dz)*
     |    (phi(lev0+1,i)-phi(lev0,i))+dvb(i)*phi(lev0,i)
      enddo ! i=lon0,lon1
      end subroutine bdzdvb
!-----------------------------------------------------------------------
      subroutine filter_op(opout,lev0,lev1,lon0,lon1,lat0,lat1,name)
!
! Filter updated O+. This is called from outside latitude loop,
! i.e., once per timestep.
!
      use mpi_module,only: mytidi,mp_gatherlons_f3d,mp_scatterlons_f3d
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(inout) :: opout(lev0:lev1,lon0:lon1,lat0:lat1)
      character(len=*),intent(in) :: name
!
! Local:
      integer :: i,j,ier
      real :: op_ik(nlonp4,lev0:lev1),op_kij(lev0:lev1,nlonp4,lat0:lat1)
!
#ifdef VT
!     code = 124 ; state = 'filter_op' ; activity='Filtering'
      call vtbegin(124,ier)
#endif
!
! Define lons in op_kij from current task subdomain opout:
      op_kij = 0._rp
      do j=lat0,lat1
        do i=lon0,lon1
          op_kij(:,i,j) = opout(:,i,j)
        enddo
      enddo ! j=lat0,lat1
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0 
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(op_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,
     |  name)
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
        if (mytidi==0) then
!
! Loop through subdomain latitudes:
          latscan: do j=lat0,lat1
!
! Define 2d array at all longitudes for filter:
            do i=1,nlonp4
              op_ik(i,:) = op_kij(:,i,j)
            enddo ! i=1,nlonp4
!
! Do the filtering:
            call ringfilter(op_ik,lev0,lev1,3,name,j)
!
! Return filtered array to op_kij:
            do i=1,nlonp4
              op_kij(:,i,j) = op_ik(i,:)
            enddo ! i=1,nlonp4
          enddo latscan ! j=lat0,lat1
        endif ! mytidi==0
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(op_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,
     |  name)
!
! Return filtered array to opout at task subdomain:
      do j=lat0,lat1
        do i=lon0,lon1
          opout(:,i,j) = op_kij(:,i,j)
        enddo
      enddo
!
#ifdef VT
!     code = 124 ; state = 'filter_op' ; activity='Filtering'
      call vtend(124,ier)
#endif
      end subroutine filter_op
!-----------------------------------------------------------------------
      end module oplus_module
