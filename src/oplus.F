!
      module oplus_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlat,nlonp4,dz,nlon,dlev
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
! 1/8/08 btf: add option to call filter() and/or filter2() for oplus.
! Default for tiegcm is callfilt1=.false. and callfilt2=.true.
      logical :: callfilt1=.false.  ! if true, call filter()
      logical :: callfilt2=.true.   ! if true, call filter2()
      contains
!-----------------------------------------------------------------------
      subroutine oplus(tn,te,ti,o2,o1,n2,he,n2d,ne,u,v,w,barm,ui,vi,wi,
     |  xnmbar,op,optm1,opout,optm1out,xiop2p,xiop2d,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Update O+ ion at 3d task subdomain.
! Outputs are opout, optm1out, xiop2p, and xiop2d, all other args
!   are input.
! There are 4 latitude scans, with 3d mpi calls in between the loops.
!   (see also 3d gather/scatter calls in sub filter_op).
!
      use params_module,only: zpmid,glat,glon  ! for O+ minimum
      use init_module,only: secs,istep         ! for O+ minimum
      use cons_module,only: rmass_op,gask,grav,re,cs,dphi,dlamda,
     |  shapiro,dtx2inv,boltz,expz,rmassinv_o2,rmassinv_o1,
     |  rmassinv_n2,rmassinv_n2d,p0,dtsmooth,dtsmooth_div2,kut
      use qrj_module,only: 
     |  qop2p, ! O+(2p) ionization from qrj, used in xiop2p
     |  qop2d, ! O+(2d) ionization from qrj, used in xiop2d
     |  qop    ! O+ ionization from qrj
      use chemrates_module,only: ! needed chemical reaction rates
     |  rk1 ,rk2 ,rk10,rk16,rk17,rk18,rk19,rk20,
     |  rk21,rk22,rk23,rk24,rk25,rk26,rk27
#ifdef MPI
      use mpi_module,only: mp_bndlons_f3d, mp_periodic_f3d,
     |  mp_polelats_f3d,mp_geo_halos_f3d
#endif
!
! Args:
      integer,intent(in) :: 
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn, te, ti, ! neutral, electron, and ion temperatures (deg K)
     |  o2, o1,     ! o2, o mass mixing ratios
     |  he,         ! he mass mixing ratio
     |  n2d,        ! n2d
     |  ne,         ! electron density
     |  u,v,w,      ! neutral wind velocities (zonal, meridional, omega)
     |  barm,       ! mean molecular mass
     |  optm1,      ! O+ at time n-1
     |  op,         ! O+ ion
     |  ui,vi,wi,   ! zonal, meridional, and vertical ion velocities
     |  xnmbar      ! p0*e(-z)*barm/kT
!
! N2 is intent(inout) for setting halo points (mp_geo_halos_f3d call below)
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(inout) :: n2
!
! Output fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  opout,    ! O+ output for next timestep
     |  optm1out, ! O+ output for time n-1
     |  xiop2p,xiop2d
!
! Local:
      integer :: k,i,lonbeg,lonend,lat,ier,nlevs
      integer :: jm2,jm1,j0,jp1,jp2 ! lat-2, lat-1, lat, lat+1, lat+2 
      real,dimension(lon0:lon1,lat0:lat1) :: 
     |  opflux,    ! upward number flux of O+ (returned by sub oplus_flux) (t7)
     |  dvb        ! output of sub divb
      real,dimension(lon0:lon1) :: 
     |  ubca, ubcb ! O+ upper boundary condition (were t2,t3)
      real :: explic = 1., gmr
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  psi_n2,           ! n2 mass mixing ratio (1-o2-o-he)
     |  bdzdvb_op,        ! was s7
     |  explicit,         ! was s4
     |  hdz,              ! was s15
     |  tphdz1,tphdz0,    ! were s13,s12 (using gmr)
     |  djint,            ! was s11
     |  divbz,            ! was s7 (DIV(B)+(DH*D*BZ)/(D*BZ)
     |  hdzmbz,hdzpbz,    ! were s10,s9
     |  p_coeff,q_coeff,r_coeff, ! coefficients for tridiagonal solver (s1,s2,s3)
     |  bdotu,            ! was s7 (B.U)
     |  op_loss,          ! was s13
     |  tp                ! te+ti
      real,dimension(lev0:lev1,lon0:lon1,lat0-1:lat1+1) :: hj ! (s10-s12)
!
! For O+ minimum:
      real :: opmin=1000., opmin_lat0=-15., opmin_lat1=15., opmin_zp=1.0
      real :: slt,dum ! only for print to stdout
!
! Local fields at 3d subdomain (must be 3d to bridge latitude scans):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) :: 
     |  bvel,
     |  diffj,       ! (D/(H*DZ)*2.*TP+M*G/R)*N(O+) (s7,s8,s9)
     |  tpj,         ! N(O+)*(te+ti) (s3-s7)
     |  bdotdh_op,   ! (b(h)*del(h))*phi
     |  bdotdh_opj,  ! (b(h)*del(h))*phi
     |  bdotdh_diff, ! (b(h)*del(h))*phi
     |  dj,          ! diffusion coefficients (s13,s14,s15)
     |  optm1_smooth ! op at time n-1, with shapiro smoother (was s1)
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,5) :: f5
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

#ifdef MPI
!
! Make latitudes j=-1,0 from j=2,1, and j=nlat+1,nlat+2 from nlat-1,nlat for n2, 
! then define 2d halos points.  Note that mp_geo_halos_f3d and mp_polelats_f3d 
! do not define "outside corner" halo points, but this should be ok since 
! cross-derivatives are not performed.
!
      call mp_polelats_f3d(n2(:,lon0:lon1,:),   ! 3rd dim is lat0-2:lat1+2
     |  lev0,lev1,lon0,lon1,lat0,lat1,1,(/1./)) ! last arg means no change in polesign
      call mp_geo_halos_f3d(n2,lev0,lev1,lon0,lon1,lat0,lat1,1)      
#endif
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
      call addfld('OPFLUX',' ',' ',opflux,
     |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Divergence is returned in dvb(lon0:lon1,lat0:lat1) by sub divb:
      call divb(dvb,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('oplus after divb.')")
!     call addfld('OP_DIVB',' ',' ',dvb,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!
!----------------------- Begin first latitude scan ---------------------
      do lat=lat0,lat1
        if (debug) write(6,"('oplus begin first lat scan: lat=',i3)")lat
        jm2 = lat-2
        jm1 = lat-1
        j0  = lat
        jp1 = lat+1
        jp2 = lat+2

        call mkdiag_BXYZ('BX',    bx (lon0:lon1,lat),lon0,lon1,lat)
        call mkdiag_BXYZ('BY',    by (lon0:lon1,lat),lon0,lon1,lat)
        call mkdiag_BXYZ('BZ',    bz (lon0:lon1,lat),lon0,lon1,lat)
        call mkdiag_BMAG('BMAG',bmod2(lon0:lon1,lat),lon0,lon1,lat)
!
! Set tp's:
! 8/28/13: Wenbin says this is wrong, it should be using tn instead of
!          te, that he presented this in 2011, but it was not changed
!          in the repository. We will get back to this later, after He 
!          is fixed.
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          tpj(k,i,jm1) = 0.5*(te(k,i,jm1)+ti(k,i,jm1))
          tpj(k,i,lat) = 0.5*(te(k,i,lat)+ti(k,i,lat))
          tpj(k,i,jp1) = 0.5*(te(k,i,jp1)+ti(k,i,jp1))
        enddo
      enddo

!     call addfld('TP_JM1',' ',' ',tpj(:,lon0:lon1,jm1),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TP_J',' ',' ',tpj(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TP_JP1',' ',' ',tpj(:,lon0:lon1,jp1),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! rrk returns djm1,dj,djp1:
!
      call rrk(
     |  tn (:,lon0:lon1,jm1),barm(:,lon0:lon1,jm1),
     |  o2 (:,lon0:lon1,jm1),o1(:,lon0:lon1,jm1),
     |  n2 (:,lon0:lon1,jm1),he(:,lon0:lon1,jm1),
     |  tpj(:,lon0:lon1,jm1),dj(:,lon0:lon1,jm1),
     |  lon0,lon1,lev0,lev1,lat)

      call rrk(
     |  tn (:,lon0:lon1,lat),barm(:,lon0:lon1,lat),
     |  o2 (:,lon0:lon1,lat),o1(:,lon0:lon1,lat),
     |  n2 (:,lon0:lon1,lat),he(:,lon0:lon1,lat),
     |  tpj(:,lon0:lon1,lat),dj(:,lon0:lon1,lat),
     |  lon0,lon1,lev0,lev1,lat)

      call rrk(
     |  tn (:,lon0:lon1,jp1),barm(:,lon0:lon1,jp1),
     |  o2 (:,lon0:lon1,jp1),o1(:,lon0:lon1,jp1),
     |  n2 (:,lon0:lon1,jp1),he(:,lon0:lon1,jp1),
     |  tpj(:,lon0:lon1,jp1),dj(:,lon0:lon1,jp1),
     |  lon0,lon1,lev0,lev1,lat)

      if (debug) write(6,"('oplus after rrk: lat=',i3)") lat
!
! Top level is invalid, so exclude it:
!     call addfld('DJM1',' ',' ',dj(lev0:lev1-1,lon0:lon1,jm1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('LOG10DJ  ',' ',' ',
!    |    log10(dj(lev0:lev1-1,lon0:lon1,lat)),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('DJP1',' ',' ',dj(lev0:lev1-1,lon0:lon1,jp1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      do i=lon0,lon1
        do k=lev0,lev1-1
          tpj(k,i,jm1) = 2.*tpj(k,i,jm1) 
          tpj(k,i,lat) = 2.*tpj(k,i,lat)
          tpj(k,i,jp1) = 2.*tpj(k,i,jp1) 
        enddo
      enddo
      if (debug) write(6,"('oplus after tpj: lat=',i3)") lat
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          hj(k,i,jm1) = gask*tn(k,i,jm1)/
     |      (0.5*(barm(k,i,jm1)+barm(k+1,i,jm1))*grav)
          hj(k,i,lat  ) = gask*tn(k,i,j0 )/
     |      (0.5*(barm(k,i,j0 )+barm(k+1,i,j0 ))*grav)
          hj(k,i,jp1) = gask*tn(k,i,jp1)/
     |      (0.5*(barm(k,i,jp1)+barm(k+1,i,jp1))*grav)
        enddo
      enddo
      if (debug) write(6,"('oplus after hj: lat=',i3)") lat

!     call addfld('HJM1',' ',' ',hj(lev0:lev1-1,:,jm1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('HJ  ',' ',' ',hj(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('HJP1',' ',' ',hj(lev0:lev1-1,:,jp1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! bvel @ jm1 = (B.U)*N(O+)    (J-1) (was s6)
! bvel @ j   = (B.U)*N(O+)      (J) (was s7)
! bvel @ jp1 = (B.U)*N(O+)    (J+1) (was s8)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          bvel(k,i,jm1) = 
     |      (bx(i,jm1)*u(k,i,jm1)+by(i,jm1)*v(k,i,jm1)+ 
     |      hj(k,i,jm1)*bz(i,jm1)*0.5*(w(k,i,jm1)+w(k+1,i,jm1)))*
     |      op(k,i,jm1) 
          bvel(k,i,lat) = 
     |      (bx(i,lat)*u(k,i,j0)+by(i,lat)*v(k,i,j0)+
     |      hj(k,i,lat)*bz(i,lat)*0.5*(w(k,i,j0)+w(k+1,i,j0)))*
     |      op(k,i,j0)
          bvel(k,i,jp1) = 
     |      (bx(i,jp1)*u(k,i,jp1)+by(i,jp1)*v(k,i,jp1)+ 
     |      hj(k,i,jp1)*bz(i,jp1)*0.5*(w(k,i,jp1)+w(k+1,i,jp1)))*
     |      op(k,i,jp1) 
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
      if (debug) write(6,"('oplus after bvel: lat=',i3)") lat

!     call addfld('BVEL_JM1',' ',' ',bvel(lev0:lev1-1,lon0:lon1,jm1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('BVEL_J'  ,' ',' ',bvel(lev0:lev1-1,lon0:lon1,j0),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('BVEL_JP1',' ',' ',bvel(lev0:lev1-1,lon0:lon1,jp1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

!     subroutine diffus(tp,en,hj,ans,lon0,lon1,lev0,lev1,lat)
!     real,dimension(lev0:lev1,lon0:lon1),intent(in) :: tp,en,hj
!     real,dimension(lev0:lev1,lon0:lon1),intent(out) :: ans
!
      tpj(lev1,:,jm1:jp1) = 0.
      call diffus(
     |  tpj  (:,lon0:lon1,jm1),op(:,lon0:lon1,jm1),hj(:,:,jm1),
     |  diffj(:,lon0:lon1,jm1),lon0,lon1,lev0,lev1,lat)
      call diffus(
     |  tpj  (:,lon0:lon1,lat),op(:,lon0:lon1,lat),hj(:,:,lat),
     |  diffj(:,lon0:lon1,lat),lon0,lon1,lev0,lev1,lat)
      call diffus(
     |  tpj  (:,lon0:lon1,jp1),op(:,lon0:lon1,jp1),hj(:,:,jp1),
     |  diffj(:,lon0:lon1,jp1),lon0,lon1,lev0,lev1,lat)

      if (debug) write(6,"('oplus after diffus: lat=',i3)") lat
!     call addfld('DIFF_JM1',' ',' ',diffj(:,lon0:lon1,jm1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('DIFF_J  ',' ',' ',diffj(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('DIFF_JP1',' ',' ',diffj(:,lon0:lon1,jp1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! tpjm2 = 2.*TP*N(O+)    (J-2)
! tpjm1 = 2.*TP*N(O+)    (J-1)
! tp    = 2.*TP*N(O+)    (J)
! tpjp1 = 2.*TP*N(O+)    (J+1)
! tpjp2 = 2.*TP*N(O+)    (J+2)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          tpj(k,i,jm2) = op(k,i,jm2)*(te(k,i,jm2)+ti(k,i,jm2))
          tpj(k,i,jm1) = tpj(k,i,jm1)*op(k,i,jm1)
          tpj(k,i,lat  ) = tpj(k,i,lat  )*op(k,i,lat)
          tpj(k,i,jp1) = tpj(k,i,jp1)*op(k,i,jp1)
          tpj(k,i,jp2) = op(k,i,jp2)*(te(k,i,jp2)+ti(k,i,jp2))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
      if (debug) write(6,"('oplus after tpj2: lat=',i3)") lat

!     call addfld('TPJM2',' ',' ',tpj(lev0:lev1-1,lon0:lon1,jm2),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('TPJM1',' ',' ',tpj(lev0:lev1-1,lon0:lon1,jm1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('TP   ',' ',' ',tpj(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('TPJP1',' ',' ',tpj(lev0:lev1-1,lon0:lon1,jp1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('TPJP2',' ',' ',tpj(lev0:lev1-1,lon0:lon1,jp2),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Shapiro smoother: optm1 is O+ at time n-1 (optm1_smooth was s1)
! optm1_smooth will be used in explicit terms below.
      do i=lon0,lon1
        do k=lev0,lev1-1
          optm1_smooth(k,i,lat) = optm1(k,i,j0)-shapiro*
     |      (optm1(k,i,jp2)+optm1(k,i,jm2)-4.*
     |      (optm1(k,i,jp1)+optm1(k,i,jm1))+6.*
     |       optm1(k,i,j0))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('OPTM1_SM0' ,' ',' ',
!    |  optm1_smooth(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
        if (debug) write(6,"('oplus end first lat scan: lat=',i3)") lat
      enddo ! lat=lat0,lat1
!------------------------- End first latitude scan ---------------------
#ifdef MPI
!
!
! Boundary longitudes:
      f5(:,:,:,1) = dj(:,:,:)
      f5(:,:,:,2) = bvel(:,:,:)
      f5(:,:,:,3) = diffj(:,:,:)
      f5(:,:,:,4) = tpj(:,:,:)
      f5(:,:,:,5) = optm1_smooth(:,:,:)

      call mp_bndlons_f3d(f5,nlevs,lon0,lon1,lat0,lat1,5,0)

      dj(:,:,:)           = f5(:,:,:,1)
      bvel(:,:,:)         = f5(:,:,:,2)
      diffj(:,:,:)        = f5(:,:,:,3)
      tpj(:,:,:)          = f5(:,:,:,4)
      optm1_smooth(:,:,:) = f5(:,:,:,5)
!
! Periodic points for dj:
      call mp_periodic_f3d(dj(:,lon0:lon1,lat0-1:lat1+1),lev0,lev1,
     |  lon0,lon1,lat0,lat1,1)
!
#endif
!----------------------- Begin second latitude scan --------------------
      do lat=lat0,lat1
        if (debug) 
     |    write(6,"('oplus begin first lat scan: lat=',i3)") lat
        jm2 = lat-2
        jm1 = lat-1
        j0  = lat
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
          bdotdh_op(k,i,lat) = dj(k,i,lat)*bz(i,lat)*bdotdh_op(k,i,lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('BDOTDH',' ',' ',bdotdh_op(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! bdotdh_opjm1 = (B(H).DEL(H))*2.*TP*N(O+)    (J-1)
! bdotdh_opj   = (B(H).DEL(H))*2.*TP*N(O+)      (J)
! bdotdh_opjp1 = (B(H).DEL(H))*2.*TP*N(O+)    (J+1)
!
      call bdotdh(
     |  tpj(:,lon0:lon1,jm2),
     |  tpj(:,:,jm1),
     |  tpj(:,lon0:lon1,lat  ),
     |  bdotdh_opj(:,lon0:lon1,jm1),lon0,lon1,lev0,lev1,jm1)
      call bdotdh(
     |  tpj(:,lon0:lon1,jm1),
     |  tpj(:,:,lat  ),
     |  tpj(:,lon0:lon1,jp1),
     |  bdotdh_opj(:,lon0:lon1,lat),lon0,lon1,lev0,lev1,lat)
      call bdotdh(
     |  tpj(:,lon0:lon1,lat  ),
     |  tpj(:,:,jp1),
     |  tpj(:,lon0:lon1,jp2),
     |  bdotdh_opj(:,lon0:lon1,jp1),lon0,lon1,lev0,lev1,jp1)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          bdotdh_opj(k,i,jm1) = bdotdh_opj(k,i,jm1)*dj(k,i,jm1)
          bdotdh_opj(k,i,lat) = bdotdh_opj(k,i,lat)*dj(k,i,lat)
          bdotdh_opj(k,i,jp1) = bdotdh_opj(k,i,jp1)*dj(k,i,jp1)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
        if (debug) 
     |    write(6,"('oplus end second lat scan: lat=',i3)") lat
      enddo ! lat=lat0,lat1
!------------------------ End second latitude scan ---------------------
!
#ifdef MPI
!
! Periodic points for bdotdh_opj (output from bdotdh above):
      call mp_periodic_f3d(bdotdh_opj(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
!
! Boundary longitudes for bdotdh_opj (input to below call to bdotdh):
      call mp_bndlons_f3d(bdotdh_opj,nlevs,lon0,lon1,lat0,lat1,1,0)
#endif
!
!----------------------- Begin third latitude scan ---------------------
      do lat=lat0,lat1
        if (debug) 
     |    write(6,"('oplus begin third lat scan: lat=',i3)") lat
        jm2 = lat-2
        jm1 = lat-1
        j0  = lat
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

!     call addfld('BDOT_DIF',' ',' ',bdotdh_diff(lev0:lev1-1,lon0:lon1,
!    |  lat),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('BDOT_JM1',' ',' ',bdotdh_opj(lev0:lev1-1,lon0:lon1,
!    |  jm1),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('BDOT_J'  ,' ',' ',bdotdh_opj(lev0:lev1-1,lon0:lon1,
!    |  lat),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('BDOT_JP1',' ',' ',bdotdh_opj(lev0:lev1-1,lon0:lon1,
!    |  jp1),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! bdzdvb_op = (BZ*D/(H*DZ)+DIV(*B))*S2
! bdzdvb returns bdzdvb_op.
!
      call bdzdvb(bdotdh_opj(:,lon0:lon1,lat),dvb(:,lat),hj(:,:,lat),
     |  bdzdvb_op,lev0,lev1,lon0,lon1,lat)

!     call addfld('BDZDVB',' ',' ',bdzdvb_op(lev0:lev1-1,:),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Collect explicit terms:
      do i=lon0,lon1
        do k=lev0,lev1-1
          explicit(k,i) = -explic*(bdzdvb_op(k,i)+bdotdh_diff(k,i,lat)+
     |      bdotdh_op(k,i,lat))
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
! Note if input flag DYNAMO<=0, then ui,vi,wi velocities will be zero.
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          explicit(k,i) = explicit(k,i)+1./(2.*re)*
     |      (1./(cs(lat)*dlamda)*(bx(i,lat)*
     |      (bvel(k,i+1,lat)-bvel(k,i-1,lat))+
     |      0.5*(ui(k,i,lat)+ui(k+1,i,lat))*bmod2(i,lat)**2*
     |      (op(k,i+1,j0)/bmod2(i+1,lat)**2-
     |       op(k,i-1,j0)/bmod2(i-1,lat)**2))+
!
     |      1./dphi*(by(i,lat)*(bvel(k,i,jp1)-bvel(k,i,jm1))+
     |      0.5*(vi(k,i,lat)+vi(k+1,i,lat))*bmod2(i,lat)**2*
     |      (op(k,i,jp1)/bmod2(i,jp1)**2-
     |       op(k,i,jm1)/bmod2(i,jm1)**2)))
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
          hdz(k,i) = 1./(hj(k,i,lat)*dz)            ! s15
          tp(k,i) = 0.5*(ti(k,i,j0)+te(k,i,j0))     ! s14
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('HDZ',' ',' ',hdz(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! gmr = G*M(O+)/(2.*R)
      gmr = grav*rmass_op/(2.*gask)
      do i=lon0,lon1
        do k=lev0,lev1-2
          tphdz1(k+1,i) = 2.*tp(k+1,i)*(0.5*(hdz(k,i)+hdz(k+1,i)))+
     |      gmr ! s13
          tphdz0(k+1,i) = 2.*tp(k  ,i)*(0.5*(hdz(k,i)+hdz(k+1,i)))-
     |      gmr ! s12
        enddo ! k=lev0,lev1-2
      enddo ! i=lon0,lon1
!
! Upper and lower boundaries:
      do i=lon0,lon1
        tphdz1(lev0,i) = 2.*tp(lev0,i)*
     |                   (1.5*hdz(lev0,i)-0.5*hdz(lev0+1,i))+gmr
        tphdz1(lev1,i) = 2.*(2.*tp(lev1-1,i)-tp(lev1-2,i))*
     |                   (1.5*hdz(lev1-1,i)-0.5*hdz(lev1-2,i))+gmr
        tphdz0(lev0,i) = 2.*(2.*tp(lev0,i)-tp(lev0+1,i))*
     |                   (1.5*hdz(lev0,i)-0.5*hdz(lev0+1,i))-gmr
        tphdz0(lev1,i) = 2.*tp(lev1-1,i)*
     |                   (1.5*hdz(lev1-1,i)-0.5*hdz(lev1-2,i))-gmr
      enddo ! i=lon0,lon1

!     call addfld('TPHDZ1',' ',' ',tphdz1,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TPHDZ0',' ',' ',tphdz0,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! djint = dj at interfaces:
      do i=lon0,lon1
        do k=lev0,lev1-2
          djint(k+1,i) = 0.5*(dj(k,i,lat)+dj(k+1,i,lat)) ! s11=.5*(s8(k)+s8(k+1))
        enddo ! k=lev0,lev1-2
        djint(lev0,i) = (1.5*dj(lev0  ,i,lat)-0.5*dj(lev0+1,i,lat))
        djint(lev1,i) = (1.5*dj(lev1-1,i,lat)-0.5*dj(lev1-2,i,lat))
      enddo ! i=lon0,lon1
!     call addfld('DJINT' ,' ',' ',djint,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! divbz = (DIV(B)+(DH*D*BZ)/(D*BZ) (was s7)
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          divbz(k,i) = 
     |      dvb(i,lat)+1./(re*dj(k,i,lat)*bz(i,lat)**2)*(bx(i,lat)/
     |      cs(lat)*(dj(k,i+1,lat)*bz(i+1,lat)-dj(k,i-1,lat)*
     |      bz(i-1,lat))/(2.*dlamda)+by(i,lat)*(dj(k,i,jp1)*
     |      bz(i,jp1)-dj(k,i,jm1)*bz(i,jm1))/(2.*dphi))
        enddo ! k=lev0,lev1-1
      enddo ! i=lonbeg,lonend

! Periodic points for divbz apparently not necessary:
!     call periodic_f2d(divbz,lon0,lon1,nlevs)
!
! Set periodic points to zero to avoid NaNS trap:
      if (lon0==1) divbz(:,lon0:lon0+1) = 0.
      if (lon1==nlonp4) divbz(:,lon1-1:lon1) = 0.
!     call addfld('DIVBZ' ,' ',' ',divbz(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! hdzmbz = (1./(H*DZ)-(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2 (was s10)
! hdzpbz = (1./(H*DZ)+(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2 (was s9 )
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          hdzmbz(k,i) = (hdz(k,i)-0.5*divbz(k,i))*bz(i,lat)**2
          hdzpbz(k,i) = (hdz(k,i)+0.5*divbz(k,i))*bz(i,lat)**2
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
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          explicit(k,i) = explicit(k,i)-(optm1_smooth(k,i,lat)-
     |      shapiro*
     |      (optm1_smooth(k,i+2,lat)+optm1_smooth(k,i-2,lat)-4.*
     |      (optm1_smooth(k,i+1,lat)+optm1_smooth(k,i-1,lat))+6.*
     |       optm1_smooth(k,i,lat)))*dtx2inv
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
          p_coeff(k,i) =   hdzmbz(k,i)*djint(k  ,i)*tphdz0(k  ,i)
          q_coeff(k,i) = -(hdzpbz(k,i)*djint(k+1,i)*tphdz0(k+1,i)+
     |                     hdzmbz(k,i)*djint(k  ,i)*tphdz1(k  ,i))
          r_coeff(k,i) =   hdzpbz(k,i)*djint(k+1,i)*tphdz1(k+1,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('P_COEFF0',' ',' ',p_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEFF0',' ',' ',q_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEFF0',' ',' ',r_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! bdotu = B.U (s7)
      do i=lon0,lon1
        do k=lev0,lev1-1
          bdotu(k,i) = bx(i,lat)*u(k,i,j0)+by(i,lat)*v(k,i,j0)+
     |      hj(k,i,lat)*bz(i,lat)*0.5*(w(k,i,j0)+w(k+1,i,j0))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('BDOTU' ,' ',' ',bdotu(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Continue coefficients with vertical ion velocity:
      do i=lon0,lon1
        do k=lev0,lev1-2
          p_coeff(k+1,i) = p_coeff(k+1,i)+
     |      (bz(i,lat)*bdotu(k,i)+0.5*(wi(k+1,i,lat)+wi(k+2,i,lat)))*
     |      0.5*hdz(k+1,i)
          q_coeff(k,i) = q_coeff(k,i)-0.5*(wi(k,i,lat)+wi(k+1,i,lat))*
     |      6./re
          r_coeff(k,i) = r_coeff(k,i)-(bz(i,lat)*bdotu(k+1,i)+
     |      0.5*(wi(k,i,lat)+wi(k+1,i,lat)))*0.5*hdz(k,i) 
        enddo ! k=lev0,lev1-2
      enddo ! i=lon0,lon1
!
!     call addfld('P_COEFF1',' ',' ',p_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEFF1',' ',' ',q_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEFF1',' ',' ',r_coeff(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Upper and lower boundaries:
      do i=lon0,lon1
        p_coeff(lev0,i) = p_coeff(lev0,i)+(bz(i,lat)*
     |    (2.*bdotu(lev0,i)-bdotu(lev0+1,i))+0.5*
     |    (wi(lev0,i,lat)+wi(lev0+1,i,lat)))*0.5*hdz(lev0,i)
        q_coeff(lev1-1,i) = q_coeff(lev1-1,i)-
     |    0.5*(wi(lev1,i,lat)+wi(lev1-1,i,lat))*6./re
        r_coeff(lev1-1,i) = r_coeff(lev1-1,i)-
     |    (bz(i,lat)*(2.*bdotu(lev1-1,i)-bdotu(lev1-2,i))+
     |    0.5*(wi(lev1,i,lat)+wi(lev1-1,i,lat)))*0.5*hdz(lev1-1,i)
      enddo ! i=lon0,lon1
!
! Additions to Q coefficients:
      do i=lon0,lon1
        do k=lev0,lev1-1
          q_coeff(k,i) = q_coeff(k,i)-bdotu(k,i)*dvb(i,lat)*bz(i,lat)-
     |      dtx2inv
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! Upper boundary condition for O+:
      do i=lon0,lon1
        ubca(i) = 0.
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
        r_coeff(lev1-1,i) = 0.
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
!
          psi_n2(k,i) = n2(k,i,j0) ! n2 mixing ratio
!
          xiop2p(k,i,lat) = 
     |      0.5*(qop2p(k,i,lat)+qop2p(k+1,i,lat))/((xnmbar(k,i,lat)*
     |      ((rk16+rk17)*psi_n2(k,i)*rmassinv_n2+
     |      rk18*o1(k,i,j0)*rmassinv_o1))+(rk19(k,i,lat)+rk20(k,i,lat))*
     |      ne(k,i,j0)+rk21+rk22)
! 
          xiop2d(k,i,lat) = 
     |      (0.5*(qop2d(k,i,lat)+qop2d(k+1,i,lat))+(rk20(k,i,lat)*
     |      ne(k,i,j0)+rk22)*xiop2p(k,i,lat))/((xnmbar(k,i,lat)*
     |      (rk23*psi_n2(k,i)*rmassinv_n2+rk24*
     |      o1(k,i,j0)*rmassinv_o1+rk26*o2(k,i,j0)*rmassinv_o2))+
     |      rk25(k,i,lat)*ne(k,i,j0)+rk27)
! 
          op_loss(k,i) = 
     |      xnmbar(k,i,lat)*(rk1(k,i,lat)*o2(k,i,j0)*rmassinv_o2+
     |      rk2(k,i,lat)*psi_n2(k,i)*rmassinv_n2+rk10*
     |      n2d(k,i,j0)*rmassinv_n2d)
!
          q_coeff(k,i) = q_coeff(k,i)-op_loss(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('XIOP2P',' ',' ',xiop2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('XIOP2D',' ',' ',xiop2d(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_LOSS',' ',' ',op_loss(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('OP_QOP',' ',' ',qop(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_NE' ,' ',' ',ne(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_O1' ,' ',' ',o1(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_TN' ,' ',' ',tn(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_BARM' ,' ',' ',barm(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_RK19' ,' ',' ',rk19(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('OP_RK25' ,' ',' ',rk25(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Add source term to RHS (explicit terms):
      do i=lon0,lon1
        do k=lev0,lev1-1
          explicit(k,i) = explicit(k,i)-
     |      0.5*(qop(k,i,lat)+qop(k+1,i,lat))-(rk19(k,i,lat)*
     |      ne(k,i,j0)+rk21)*
     |      xiop2p(k,i,lat)-(rk25(k,i,lat)*ne(k,i,j0)+rk27)*
     |      xiop2d(k,i,lat)-(rk18*xiop2p(k,i,lat)+rk24*xiop2d(k,i,lat))*
     |      o1(k,i,j0)*rmassinv_o1*p0*expz(k)*0.5*(barm(k,i,j0)+
     |      barm(k+1,i,j0))/(boltz*tn(k,i,j0))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('EXPLIC4',' ',' ',explicit(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Lower boundary condition N(O+) = Q/L:
      do i=lon0,lon1
        q_coeff(lev0,i) = q_coeff(lev0,i)-p_coeff(lev0,i) 
        explicit(lev0,i) = explicit(lev0,i)-2.*p_coeff(lev0,i)*
     |    qop(lev0,i,lat)/(1.5*op_loss(lev0,i)-0.5*op_loss(lev0+1,i))
        p_coeff(lev0,i) = 0.
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
!
! Filter updated O+:
! 
      call filter_op(opout(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,kut,'OPLUS')
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
      opout(lev1,:,lat) = 0.
      optm1out(lev1,:,lat) = 0.
!
! Insure global non-negative O+:
      do i=lon0,lon1
        do k=lev0,lev1-1      
          if (opout(k,i,lat)    < 1.e-5) opout(k,i,lat)    = 1.e-5
          if (optm1out(k,i,lat) < 1.e-5) optm1out(k,i,lat) = 1.e-5
        enddo ! k=lev0,lev1-1      
      enddo ! i=lon0,lon1
!
! Report magnitude and spatial extent of imposed O+ minimum at
! first timestep (see opmin settings above):
!
      if (istep==1.and.lat==lat0)
     |  write(6,"('OPLUS: opmin=',e12.4,' opmin_lat0,1=',2e12.4,
     |    ' opmin_zp=',f8.2)") opmin,opmin_lat0,opmin_lat1,opmin_zp
!
! May-June, 2014 (B.Foster and Wenbin): Limit O+ to minimum opmin at 
!   low-latitudes and above given zp level (see opmin settings above):
!
      if (glat(lat) <= opmin_lat1.and.glat(lat) >= opmin_lat0) then
        do i=lon0,lon1
          do k=lev0,lev1-1      
            if (zpmid(k) >= opmin_zp .and. opout(k,i,lat) < opmin) then
!             if (opout(k,i,lat) /= 1.e-5) then ! don't report global min
!               slt = fslt(dum,secs/3600.,glon(i),1)
!               write(6,"('OPLUS: apply O+ min: glat=',f8.2,' glon=',
!    |            f8.2,' zpmid=',f8.2,' opout=',e12.3,' uthr=',f8.2,
!    |            ' slt=',f8.2)") glat(lat),glon(i),zpmid(k),
!    |            opout(k,i,lat),secs/3600.,slt
!             endif ! report to stdout
              opout(k,i,lat) = opmin
            endif ! apply minimum
          enddo
        enddo  
      endif
!
! End fourth and final latitude scan:
      enddo ! lat=lat0,lat1
!
! Periodic points for outputs:
#ifdef MPI
      call mp_periodic_f3d(opout(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
      call mp_periodic_f3d(optm1out(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
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
     |  phid =  2.0e8,
     |  phin = -2.0e8,
     |  ppolar = 0.
      real :: a(lon0:lon1)   ! was t3 ("a" needs a better name)
      real :: fed(lon0:lon1) ! was t4
      real :: fen(lon0:lon1) ! was t5
!
! Latitude scan:
      do lat=lat0,lat1
!
! Longitude loop:
        do i=lon0,lon1
          if (abs(rlatm(i,lat))-pi/24.>=0.) then
            a(i) = 1.
          else
            a(i)=.5*(1.+sin(pi*(abs(rlatm(i,lat))-pi/48.)/(pi/24.)))
            if (a(i) < 0.05) a(i) = 0.05
          endif
          fed(i) = phid*a(i)
          fen(i) = phin*a(i)
          if (chi(i,lat)-0.5*pi >= 0.) then
            opflux(i,lat) = fen(i)
          else
            opflux(i,lat) = fed(i)
          endif
          if ((chi(i,lat)*rtd-80.)*(chi(i,lat)*rtd-100.) < 0.) 
     |      opflux(i,lat) = .5*(fed(i)+fen(i))+.5*(fed(i)-fen(i))*
     |        cos(pi*(chi(i,lat)*rtd-80.)/20.)
!
! Add ppolar if magnetic latitude >= 60 degrees:
          if (abs(rlatm(i,lat))-pi/3. >= 0.) 
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
        dvb(:,lat) = 0.
        do i=lonbeg,lonend
          dvb(i,lat) = (((bx(i+1,lat)-bx(i-1,lat))/(2.*dlamda)+
     |      (cs(jp1)*by(i,jp1)-cs(jm1)*by(i,jm1))/(2.*dphi))/
     |      cs(lat)+2.*bz(i,lat))/re
        enddo ! i=lonbeg,lonend
      enddo ! lat=lat0,lat1
!
      end subroutine divb
!-----------------------------------------------------------------------
      subroutine rrk(t,rms,ps1,ps2,n2,he,tp,ans,lon0,lon1,lev0,lev1,lat)
!
! Returns ambipolar diffusion coefficient in ans.
!
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,boltz,
     |  p0,expz,rmassinv_he
      use input_module,only: colfac
!
! Args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in) ::
     |  t,rms,ps1,ps2,n2,he,tp
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: ans
!
! Local:
      integer :: k,i
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! 8/28/13 btf: Use n2=1-o2-o-he, and include O+/He collision rate from Wenbin:
!
          ans(k,i) = 1.42E17*boltz*t(k,i)/(p0*expz(k)*.5*(rms(k,i)+
     |      rms(k+1,i))*(ps2(k,i)*rmassinv_o1*sqrt(tp(k,i))*(1.-0.064*
     |      alog10(tp(k,i)))**2*colfac+18.6*n2(k,i)*rmassinv_n2+18.1*
     |      ps1(k,i)*rmassinv_o2+3.6*he(k,i)*rmassinv_he))

        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
      end subroutine rrk
!-----------------------------------------------------------------------
      subroutine diffus(tp,en,hj,ans,lon0,lon1,lev0,lev1,lat)
!
! Evaluates ans = (d/(h*dz)*tp+m*g/r)*en
!
      use cons_module,only: rmass_op,grav,gask
!
! Args:
      integer :: lon0,lon1,lev0,lev1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in) ::
     |  tp,en,hj
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
          ans(k1,i) = 1./(2.*hj(k1,i)*dlev)*(tp(k2,i)*en(k2,i)-
     |      tp(k,i)*en(k,i))+mgr*en(k1,i)
        enddo
      enddo
!
! Upper and lower boundaries:
      k1 = lev1-1
      k2 = lev1-2
      do i=lon0,lon1
        ans(k1,i) = 1./(hj(k1,i)*dlev)*(tp(k1,i)*en(k1,i)-
     |    tp(k2,i)*en(k2,i))+mgr*en(k1,i)
        ans(1,i) = 1./(hj(1,i)*dlev)*(tp(2,i)*en(2,i)-
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
      integer :: iperiodic,lon0,lon1,lev0,lev1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in)    :: phijm1,phijp1
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(inout) :: phij
      real,dimension(lev0:lev1,lon0:lon1),intent(out)   :: ans
!
! Local:
      integer :: k,i,lonbeg,lonend
!
      lonbeg = lon0
      if (lon0==1) then
        lonbeg = 3
        ans(:,lon0:lon1) = 0. ! set periodic points to zero to avoid NaNS trap
      endif
      lonend = lon1
      if (lon1==nlonp4) then
        lonend = lon1-2
        ans(:,lon1-1:lon1) = 0. ! set periodic points to zero to avoid NaNS trap
      endif
!
! Note phij longitude dimension is lon0-2:lon1+2 (only i-1 and i+1 are used).
! Boundary longitudes i-1 and i+1 must have been set before this routine is
! called (e.g., call mp_bndlons_f3d).
!
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          ans(k,i) = 1./re*(bx(i,lat)/(cs(lat)*2.*dlamda)*
     |      (phij(k,i+1)-phij(k,i-1))+by(i,lat)*
     |      (phijp1(k,i)-phijm1(k,i))/(2.*dphi))
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
          ans(k,i) = bz(i,lat)/(2.*h(k,i)*dz)*(phi(k+1,i)-phi(k-1,i))+
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
      subroutine filter_op(opout,lev0,lev1,lon0,lon1,lat0,lat1,kut,name)
!
! Filter updated O+. This is called from outside latitude loop,
! i.e., once per timestep.
!
      use filter_module,only: filter,filter2
#ifdef MPI
      use mpi_module,only: mytidi,mp_gatherlons_f3d,mp_scatterlons_f3d
#else
      integer :: mytidi=0
#endif
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      integer,intent(in) :: kut(nlat)
      real,intent(inout) :: opout(lev0:lev1,lon0:lon1,lat0:lat1)
      character(len=*),intent(in) :: name
!
! Local:
      integer :: i,j,ier
      real :: op_ik(nlonp4,lev0:lev1),op_kij(lev0:lev1,nlonp4,lat0:lat1)
      real :: fmin,fmax
!
#ifdef VT
!     code = 124 ; state = 'filter_op' ; activity='Filtering'
      call vtbegin(124,ier)
#endif
!
! Define lons in op_kij from current task subdomain opout:
      op_kij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          op_kij(:,i,j) = opout(:,i,j)
        enddo
      enddo ! j=lat0,lat1
!
#ifdef MPI
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0 
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(op_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,
     |  name)
#endif
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
        if (mytidi==0) then
!
! Loop through subdomain latitudes:
          latscan: do j=lat0,lat1
!
! Define 2d array at all longitudes for filter2:
            do i=1,nlonp4
              op_ik(i,:) = op_kij(:,i,j)
            enddo ! i=1,nlonp4
!
! Do the filtering (requires fft):
! call filter and/or filter2 according to module logicals callfilt1,callfilt2:
!
! If requested, call filter2(), as in composition comp_major and minor:
            if (callfilt2) call filter2(op_ik,lev0,lev1,kut(j),j)
!
! If requested, call filter(), as in dt,duv,etc:
!
            if (callfilt1) call filter(op_ik,lev0,lev1,kut(j),j)
!
! Return filtered array to op_kij:
            do i=1,nlonp4
              op_kij(:,i,j) = op_ik(i,:)
            enddo ! i=1,nlonp4
          enddo latscan ! j=lat0,lat1
        endif ! mytidi==0
#ifdef MPI
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(op_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,
     |  name)
#endif
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
