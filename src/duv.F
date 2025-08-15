!
      subroutine duv(tn,tn_upd,tn_nm,un,vn,un_nm,vn_nm,w_upd,mbar,
     |  scht,schti,z,hdu,hdv,ui,vi,lxx,lyy,lxy,lyx,km,Fe,Fn,
     |  un_upd,unm_upd,vn_upd,vnm_upd,lev0,lev1,lon0,lon1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Advance neutral velocities for current time step:
!
      use params_module,only: nlonp4,dz,spval
      use init_module,only: iday
      use lbc,only: u_lbc,v_lbc
      use cons_module,only: dtx2inv,expz,p0,grav,
     |  xmue,cor,tanphi=>tn,re,dtsmooth_div2,dtsmooth,shapiro
      use fields_module,only: ulbc,vlbc,ulbc_nm,vlbc_nm
      use addfld_module,only: addfld
      use mpi_module,only: mp_periodic_f3d,mp_periodic_f2d,
     |  mp_bndlons_f2d,mp_bndlats_f2d
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Inputs at full subdomain:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,      ! neutral temperature (deg K)
     |  tn_upd,  ! updated neutral temperature (from dt.F)
     |  tn_nm,   ! neutral temperature from time n-1
     |  un,      ! neutral zonal velocity
     |  vn,      ! neutral meridional velocity
     |  un_nm,   ! neutral zonal velocity at time n-1
     |  vn_nm,   ! neutral meridional velocity at time n-1
     |  w_upd,   ! updated vertical velocity (swdot.F)
     |  mbar,    ! mean molecular weight
     |  scht,    ! scale height at midpoints
     |  schti,   ! scale height at interfaces
     |  z,       ! geopotential height
     |  hdu,     ! horizontal diffusion of U (hdif3 in hdif.F)
     |  hdv,     ! horizontal diffusion of V (hdif3 in hdif.F)
     |  ui,      ! zonal ion velocity input
     |  vi,      ! meridional ion velocity input
     |  lxx,     ! xx ion drag coefficient
     |  lyy,     ! yy ion drag coefficient
     |  lxy,     ! xy ion drag coefficient
     |  lyx,     ! yx ion drag coefficient
     |  km,      ! molecular viscosity (cpktkm.F)
     |  Fe,Fn    ! field-aligned ion drag
!
! Outputs at full subdomain:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out)::
     |  un_upd,  ! updated zonal velocity
     |  unm_upd, ! updated zonal velocity at time n-1
     |  vn_upd,  ! updated meridional velocity
     |  vnm_upd  ! updated meridional velocity at time n-1
!
! Local:
      integer :: k,i,lat,lonbeg,lonend,m,n
      real :: wi
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  g,g1,        ! g**2*(kM+XMUE)*MBAR*/(P0*R*T*DS**2)
     |  dwdz,        ! exp(-s)*w(k+1/2)/(2.*Ds)
     |  cen
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  advec_un,    ! zonal advection (output of sub advec)
     |  advec_vn,    ! meridional advection (output of sub advec)
     |  zl,zp,       ! output from glp (and dldp)
     |  unm_smooth,  ! un at time n-1, smoothed
     |  vnm_smooth   ! vn at time n-1, smoothed
      real :: flbc(lon0-2:lon1+2,lat0-2:lat1+2,4) ! ulbc,vlbc,ulbc_nm,vlbc_n
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  eddyvisc,eddyvisci     ! eddy viscosity (?)
!
! For diagnostic, debug:
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  dudt,dvdt,du_cor,dv_cor,du_cen,dv_cen,du_drag,dv_drag,
     |  un_drag_coef,vn_drag_coef,du_visc,dv_visc,du_veradv,dv_veradv
!
! pp, qq, rr, beta AND gamma are (2 x 2) matrices defined for 
!   k = 3/2, K-1/2, 1 and i = 3, IMAX+2, 1
      real,dimension(2,2,lev0:lev1,lon0:lon1) :: 
     |  pp,qq,rr
!
! ss, xx and yy are two component vectors similarly defined.
      real,dimension(2,lev0:lev1,lon0:lon1) :: 
     |  ss,xx
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Horizontal advection:
!     subroutine advec(un,vn,f,hadvec,lev0,lev1,lon0,lon1,lat0,lat1)
!     real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in) :: 
!    |  un,vn,f ! input field with ghost cells for finite differencing
!     real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: hadvec
!      
      call advec(un,vn,un,advec_un,lev0,lev1,lon0,lon1,lat0,lat1)
      call advec(un,vn,vn,advec_vn,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Horizontal pressure forcing output of glp is in zl,zp. Pass full 
! subdomains for derivatives in lat and lon (see sub dldp, called by glp):
!
      call glp(tn,tn_upd,tn_nm,mbar,z,zl,zp,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Smooth un,vn at time n-1 with 2-part shapiro smoother.
! Input un_nm,vn_nm at full subdomain.
! Output in unm_smooth,vnm_smooth at subdomain.
!
      call smooth(un_nm,unm_smooth,shapiro,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
      call smooth(vn_nm,vnm_smooth,shapiro,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Start third latitude scan:
      do lat=lat0,lat1

!       call addfld('ADVEC_U0','ADVEC_U0',' ',advec_un(:,:,lat),
!    |    'ilev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('ADVEC_V0','ADVEC_V0',' ',advec_vn(:,:,lat),
!    |    'ilev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('UNSMOOTH','UNSMOOTH',' ',unm_smooth(:,:,lat),
!    |     'ilev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VNSMOOTH','VNSMOOTH',' ',vnm_smooth(:,:,lat),
!    |     'ilev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('UN_NM_DUV','UN_NM_DUV',' ',un_nm(:,lon0:lon1,lat),
!    |     'ilev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VN_NM_DUV','VN_NM_DUV',' ',vn_nm(:,lon0:lon1,lat),
!    |     'ilev',lev0,lev1,'lon',lon0,lon1,lat)

        do i=lon0,lon1
          do k=lev0,lev1-1
            un_drag_coef(k,i) = Fe(k,i,lat)+.25*
     |       ((lxx(k,i,lat)+lxx(k+1,i,lat))*(ui(k,i,lat)+ui(k+1,i,lat))+
     |        (lxy(k,i,lat)+lxy(k+1,i,lat))*(vi(k,i,lat)+vi(k+1,i,lat)))
            vn_drag_coef(k,i) = Fn(k,i,lat)+.25*
     |       ((lyy(k,i,lat)+lyy(k+1,i,lat))*(vi(k,i,lat)+vi(k+1,i,lat))-
     |        (lyx(k,i,lat)+lyx(k+1,i,lat))*(ui(k,i,lat)+ui(k+1,i,lat)))
            ss(1,k,i) = expz(k)*(dtx2inv*unm_smooth(k,i,lat)+
     |        hdu(k,i,lat)+un_drag_coef(k,i)-
     |        advec_un(k,i,lat)-zl(k,i,lat))
            ss(2,k,i) = expz(k)*(dtx2inv*vnm_smooth(k,i,lat)+
     |        hdv(k,i,lat)+vn_drag_coef(k,i)-
     |        advec_vn(k,i,lat)-zp(k,i,lat))
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
!
        do i=lon0,lon1
          do k=lev0,lev1-1          
            eddyvisc(k,i) = p0*expz(k)*scht(k,i,lat)*xmue(k,iday)/grav
          enddo ! k=lev0+1,lev1-1          
!
! Upper boundary:
          eddyvisci(lev1,i) =
     |      sqrt(eddyvisc(lev1-1,i)**3/eddyvisc(lev1-2,i))            ! s3 top
!            
          do k=lev0+1,lev1-1
            eddyvisci(k,i) = sqrt(eddyvisc(k-1,i)*eddyvisc(k,i))      ! s3
          enddo ! k=lev0+1,lev1-1
!
! Bottom boundary in tlbc (from t_lbc of lbc.F):
          eddyvisci(lev0,i) =
     |      sqrt(eddyvisc(lev0,i)**3/eddyvisc(lev0+1,i))              ! s3 bottom
        enddo ! i=lon0,lon1
!
        do i=lon0,lon1

! G = g*KM/(P0*H*Ds**2) = g**2*(kM+XMUE)*MBAR*/(P0*R*T*DS**2)
          do k=lev0,lev1
            g(k,i) = grav*(km(k,i,lat)+eddyvisci(k,i))/
     |        (p0*schti(k,i,lat)*dz**2)
          enddo ! k=lev0,lev1
          do k=lev0,lev1-1
            g1(k,i) = g(k+1,i)
          enddo
          g1(lev1,i) = 2*g(lev1,i)-g(lev1-1,i)
        enddo ! i=lon0,lon1

!       call addfld('DUV_G',' ',' ',g  ,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

!
! dwdz = exp(-s)*w(k+1/2)/(2.*Ds)
        do i=lon0,lon1
          do k=lev0,lev1-1
            wi = 0.5*(w_upd(k,i,lat)+w_upd(k+1,i,lat))
            dwdz(k,i) = expz(k)*wi/(2.*dz)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('DWDZ',' ',' ',dwdz,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

!
! P = -(G(k-1/2,n) + exp(-s)*w(k,n)/   |1., 0.|
!                            (2.*Ds)) *|0., 1.|
!                                      |_    _|
!
! R = -(G(k+1/2,n) - exp(-s)*w(k,n)/   |1., 0.|
!                            (2.*Ds)) *|0., 1.|
!                                      |_    _|
!  for k = 3/2, K-1/2, 1
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            pp(1,1,k,i) = -(g(k  ,i)+dwdz(k,i))
            rr(1,1,k,i) = -(g1(k ,i)-dwdz(k,i))
            pp(2,2,k,i) = pp(1,1,k,i)
            rr(2,2,k,i) = rr(1,1,k,i)
            pp(1,2,k,i) = 0.
            rr(1,2,k,i) = 0.
            pp(2,1,k,i) = 0.
            rr(2,1,k,i) = 0.
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
!
! Now calculate Q from (5) for k = 3/2, K-1/2, 1
!  Q(1,1,k) = G(k-1/2) + G(k+1/2) + exp(-s)*(1./(2.*Dt) + lamxx(k))
!    and
!  Q(2,2,k) = G(k-1/2) + G(k+1/2) + exp(-s)*(1./(2.*Dt) + lamyy(k))
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            qq(1,1,k,i) = g(k,i) + g1(k,i) + 
     |        expz(k)*(dtx2inv+.5*(lxx(k,i,lat)+lxx(k+1,i,lat)))
            qq(2,2,k,i) = g(k,i) + g1(k,i) +
     |        expz(k)*(dtx2inv+.5*(lyy(k,i,lat)+lyy(k+1,i,lat)))
!
!     Q(1,2) = -exp(-s) * (f+u/r*tan(theta)-lamxy)
!      and 
!     Q(2,1) = exp(-s) * (f+u/r*tan(theta)-lamyx)
!
            cen(k,i) = un(k,i,lat)/re*tanphi(lat)
            qq(1,2,k,i) = -expz(k)*(cor(lat)+cen(k,i)-
     |        .5*(lxy(k,i,lat)+lxy(k+1,i,lat)))
            qq(2,1,k,i) =  expz(k)*(cor(lat)+cen(k,i)-
     |        .5*(lyx(k,i,lat)+lyx(k+1,i,lat)))
          enddo ! k=lev0,lev1-1
!
! Boundaries at k=3/2 and k=k-1/2
! Lower boundary:
!     Q(3/2) = Q(3/2) - P(3/2)
!     S(3/2) = S(3/2) -2.*P(3/2)*Vb
!     P(3/2) = 0.
!
          do n=1,2
            do m=1,2
              qq(m,n,lev0,i) = qq(m,n,lev0,i)-pp(m,n,lev0,i)
            enddo
          enddo
          do n=1,2
            ss(n,lev0,i) = ss(n,lev0,i)-2.*(pp(n,1,lev0,i)*u_lbc(i,lat)+
     |                                      pp(n,2,lev0,i)*v_lbc(i,lat))
          enddo
          do n=1,2
            do m=1,2
              pp(m,n,lev0,i) = 0.
            enddo
          enddo
!
! Upper boundary:
!     Q(K-1/2) = Q(K-1/2) + R(K-1/2)
!     R(K-1/2) = 0.
!
          do n=1,2
            do m=1,2
              qq(m,n,lev1-1,i) = qq(m,n,lev1-1,i)+rr(m,n,lev1-1,i)
            enddo
          enddo
          do n=1,2
            do m=1,2
              rr(m,n,lev1-1,i) = 0.
            enddo
          enddo
        enddo ! i=lon0,lon1

!       call addfld('PP11',' ',' ',pp(1,1,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('PP12',' ',' ',pp(1,2,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('PP21',' ',' ',pp(2,1,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('PP22',' ',' ',pp(2,2,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QQ11',' ',' ',qq(1,1,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QQ12',' ',' ',qq(1,2,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QQ21',' ',' ',qq(2,1,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QQ22',' ',' ',qq(2,2,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RR11',' ',' ',rr(1,1,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RR12',' ',' ',rr(1,2,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RR21',' ',' ',rr(2,1,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RR22',' ',' ',rr(2,2,:,:),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Invoke tridiagonal (2x2) block matrix solver for U and V: 
! pp,qq,rr,ss are input, beta,gamma,yy,xx are output.
! (final solution is in xx)
!
        call blktri(pp,qq,rr,ss,xx,lev0,lev1,lon0,lon1)
!
! Save updated U and V:
        do i=lonbeg,lonend 
          do k=lev0,lev1-1
            un_upd(k,i,lat) = xx(1,k,i)
            vn_upd(k,i,lat) = xx(2,k,i)
          enddo ! k=lev0,lev1-1
!
! Put spval in top nlevp1 level:
          un_upd(lev1,i,lat) = spval
          vn_upd(lev1,i,lat) = spval
!
! Lower boundary is in ulbc,vlbc (fields.F):
          ulbc_nm(i,lat) = ulbc(i,lat)  ! LB for un_nm(itc) = un(itp) 
          vlbc_nm(i,lat) = vlbc(i,lat)  ! LB for vn_nm(itc) = vn(itp) 
          ulbc(i,lat) = u_lbc(i,lat)    ! LB for un(itc)
          vlbc(i,lat) = v_lbc(i,lat)    ! LB for vn(itc) 
        enddo ! i=lonbeg,lonend 

!       call addfld('UN_SOLV',' ',' ',un_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VN_SOLV',' ',' ',vn_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            dudt(k,i) = dtx2inv*(un_upd(k,i,lat)-unm_smooth(k,i,lat))
            dvdt(k,i) = dtx2inv*(vn_upd(k,i,lat)-vnm_smooth(k,i,lat))
            du_cor(k,i) = cor(lat)*vn_upd(k,i,lat)
            dv_cor(k,i) = -cor(lat)*un_upd(k,i,lat)
            du_cen(k,i) = cen(k,i)*vn_upd(k,i,lat)
            dv_cen(k,i) = -cen(k,i)*un_upd(k,i,lat)
            du_drag(k,i) = un_drag_coef(k,i)-
     |        .5*(lxx(k,i,lat)+lxx(k+1,i,lat))*un_upd(k,i,lat)-
     |        .5*(lxy(k,i,lat)+lxy(k+1,i,lat))*vn_upd(k,i,lat)
            dv_drag(k,i) = vn_drag_coef(k,i)-
     |        .5*(lyy(k,i,lat)+lyy(k+1,i,lat))*vn_upd(k,i,lat)+
     |        .5*(lyx(k,i,lat)+lyx(k+1,i,lat))*un_upd(k,i,lat)
          enddo
          do k=lev0+1,lev1-2
            wi = 0.5*(w_upd(k,i,lat)+w_upd(k+1,i,lat))
            du_visc(k,i) = 1/expz(k)*
     |        (g (k,i)        *un_upd(k-1,i,lat)-
     |        (g1(k,i)+g(k,i))*un_upd(k  ,i,lat)+
     |         g1(k,i)        *un_upd(k+1,i,lat))
            dv_visc(k,i) = 1/expz(k)*
     |        (g (k,i)        *vn_upd(k-1,i,lat)-
     |        (g1(k,i)+g(k,i))*vn_upd(k  ,i,lat)+
     |         g1(k,i)        *vn_upd(k+1,i,lat))
            du_veradv(k,i) = wi/(2*dz)*
     |        (un_upd(k+1,i,lat)-un_upd(k-1,i,lat))
            dv_veradv(k,i) = wi/(2*dz)*
     |        (vn_upd(k+1,i,lat)-vn_upd(k-1,i,lat))
          enddo
        enddo
        call addfld('DUDT','U changing rate','cm/s2',
     |    dudt(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DVDT','V changing rate','cm/s2',
     |    dvdt(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DU_CORIOLIS','U Coriolis force','cm/s2',
     |    du_cor(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DV_CORIOLIS','V Coriolis force','cm/s2',
     |    dv_cor(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DU_CENTRIFUGAL','U centrifugal force','cm/s2',
     |    du_cen(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DV_CENTRIFUGAL','V centrifugal force','cm/s2',
     |    dv_cen(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DU_IONDRAG','U ion drag','cm/s2',
     |    du_drag(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DV_IONDRAG','V ion drag','cm/s2',
     |    dv_drag(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DU_PRESSURE','U pressure gradient force','cm/s2',
     |    -zl(lev0:lev1-1,:,lat),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DV_PRESSURE','V pressure gradient force','cm/s2',
     |    -zp(lev0:lev1-1,:,lat),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DU_HORDIF','U horizontal diffusion','cm/s2',
     |    hdu(lev0:lev1-1,lon0:lon1,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DV_HORDIF','V horizontal diffusion','cm/s2',
     |    hdv(lev0:lev1-1,lon0:lon1,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DU_VISCOUS','U viscous force','cm/s2',
     |    du_visc(lev0+1:lev1-2,:),
     |    'lev',lev0+1,lev1-2,'lon',lon0,lon1,lat)
        call addfld('DV_VISCOUS','V viscous force','cm/s2',
     |    dv_visc(lev0+1:lev1-2,:),
     |    'lev',lev0+1,lev1-2,'lon',lon0,lon1,lat)
        call addfld('DU_HORADV','U horizontal advection','cm/s2',
     |    -advec_un(lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DV_HORADV','V horizontal advection','cm/s2',
     |    -advec_vn(lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('DU_VERADV','U vertical advection','cm/s2',
     |    -du_veradv(lev0+1:lev1-2,:),
     |    'lev',lev0+1,lev1-2,'lon',lon0,lon1,lat)
        call addfld('DV_VERADV','V vertical advection','cm/s2',
     |    -dv_veradv(lev0+1:lev1-2,:),
     |    'lev',lev0+1,lev1-2,'lon',lon0,lon1,lat)
!
! End third latitude scan:
      enddo ! lat=lat0,lat1
!
! Fourier smoothing of U and V:
      call filter_uv(un_upd(lev0:lev1-1,:,:),lev0,lev1-1,lon0,lon1,
     |  lat0,lat1,'UN')
      call filter_uv(vn_upd(lev0:lev1-1,:,:),lev0,lev1-1,lon0,lon1,
     |  lat0,lat1,'VN')
!
! Fourth latitude scan:
      do lat=lat0,lat1
!       call addfld('UN_FILT',' ',' ',un_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VN_FILT',' ',' ',vn_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Smooth updated un,vn:
        do i=lon0,lon1 
          do k=lev0,lev1-1
            unm_upd(k,i,lat) = dtsmooth_div2*(un_nm(k,i,lat)+
     |        un_upd(k,i,lat)) + dtsmooth*un(k,i,lat)
            vnm_upd(k,i,lat) = dtsmooth_div2*(vn_nm(k,i,lat)+
     |        vn_upd(k,i,lat)) + dtsmooth*vn(k,i,lat)
          enddo ! k=lev0,lev1-1
!
! put spval in top nlevp1 level:
          unm_upd(lev1,i,lat) = spval
          vnm_upd(lev1,i,lat) = spval
        enddo ! i=lon0,lon1 
!       call addfld('UN_NMOUT',' ',' ',unm_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VN_NMOUT',' ',' ',vnm_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End fourth latitude scan:
      enddo ! lat=lat0,lat1
!
! Boundary (halo) latitudes:
      flbc(:,:,1) = ulbc(:,:)
      flbc(:,:,2) = vlbc(:,:)
      flbc(:,:,3) = ulbc_nm(:,:)
      flbc(:,:,4) = vlbc_nm(:,:)
      call mp_bndlats_f2d(flbc,lon0,lon1,lat0,lat1,4)
      ulbc(:,:)    = flbc(:,:,1)
      vlbc(:,:)    = flbc(:,:,2)
      ulbc_nm(:,:) = flbc(:,:,3)
      vlbc_nm(:,:) = flbc(:,:,4)
!
! Boundary (halo) longitudes:
      flbc(:,:,1) = ulbc(:,:)
      flbc(:,:,2) = vlbc(:,:)
      flbc(:,:,3) = ulbc_nm(:,:)
      flbc(:,:,4) = vlbc_nm(:,:)
      call mp_bndlons_f2d(flbc,lon0,lon1,lat0,lat1,4)
      ulbc(:,:)    = flbc(:,:,1)
      vlbc(:,:)    = flbc(:,:,2)
      ulbc_nm(:,:) = flbc(:,:,3)
      vlbc_nm(:,:) = flbc(:,:,4)
!
! Periodic points for lbc:
      flbc(:,:,1) = ulbc(:,:)
      flbc(:,:,2) = vlbc(:,:)
      flbc(:,:,3) = ulbc_nm(:,:)
      flbc(:,:,4) = vlbc_nm(:,:)
      call mp_periodic_f2d(flbc(lon0:lon1,lat0:lat1,1:4),
     |       lon0,lon1,lat0,lat1,4)
      ulbc(:,:)    = flbc(:,:,1)
      vlbc(:,:)    = flbc(:,:,2)
      ulbc_nm(:,:) = flbc(:,:,3)
      vlbc_nm(:,:) = flbc(:,:,4)
!
!     do lat=lat0,lat1
!       call addfld('UN_FINAL',' ',' ',un_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VN_FINAL',' ',' ',vn_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
      end subroutine duv
!-----------------------------------------------------------------------
      subroutine glp(tn,tn_upd,tn_nm,mbar,z,zl,zp,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Horizontal pressure forcing for un,vn: 
!
      use params_module,only: nlonp4
      use cons_module,only: dz,gask,grav,racs,re_inv,dlamda,dphi
      use addfld_module,only: addfld
      use mpi_module,only: mp_bndlons_f3d, mp_periodic_f3d,
     |  mp_polelats_f3d,mp_bndlats_f3d
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Inputs at full subdomain:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,      ! neutral temperature (deg K)
     |  tn_upd,  ! updated neutral temperature (from dt.F)
     |  tn_nm,   ! neutral temperature from time n-1
     |  mbar,    ! mean molecular weight
     |  z        ! geopotential
!
! Outputs at subdomain:
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: 
     |  zl,zp
!
! Local:
      integer :: k,i,nlevs,lat,lonbeg,lonend
      real,parameter :: wt=0.225
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  ztmp
      real,dimension(lev0:lev1,lon0:lon1) :: tbar,scht
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) :: dzdx,dzdy,zli,zpi
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1,2) :: ftmp
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2

      nlevs = lev1-lev0+1
!
! Latitude scan:
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            tbar(k,i) = tn(k,i,lat)+wt*
     |        (tn_nm(k,i,lat)+tn_upd(k,i,lat)-2.*tn(k,i,lat))
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('TBAR',' ',' ',tbar(:,lon0:lon1),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

        do i=lon0,lon1
          do k=lev0,lev1-1
            scht(k,i) = gask*tbar(k,i)/(mbar(k,i,lat)*grav)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('SCHT_TBAR',' ',' ',scht(:,lon0:lon1),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! In earlier model versions, geopotential bottom boundary was recalculated
! here. In this code, it is just copied from z, where it is already available:
        do i=lon0,lon1
          ztmp(lev0,i,lat) = z(lev0,i,lat)
        enddo ! i=lon0,lon1
!
! This will differ from z(..[itp,itc]) because tn_upd is in tbar.
        do i=lon0,lon1
          do k=lev0,lev1-1
            ztmp(k+1,i,lat) = dz*scht(k,i)+ztmp(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo
!
! Exchange boundary longitudes in ztmp for 4th order differencing in dldp. 
! Not necessary if ztmp has been defined at full subdomain, including lon 
! halos, but as of 4/02, this was not working (s.a. comments above)
!
      call mp_periodic_f3d(ztmp(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
      call mp_polelats_f3d(ztmp(:,lon0:lon1,:),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1,(/1./))
      call mp_bndlats_f3d(ztmp,nlevs,lon0,lon1,lat0,lat1,1)
      call mp_bndlons_f3d(ztmp,nlevs,lon0,lon1,lat0,lat1,1,0)
!
!     do lat=lat0,lat1
!       call addfld('ZTMP',' ',' ',ztmp(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! Calculate latitudinal and longitudinal derivatives in zl,zp:
! (dldp input ztmp is defined at full subdomain)
!
      call divergence(ztmp,dzdx,-1,lev0,lev1,lon0,lon1,lat0,lat1)
      call divergence(z   ,dzdy, 1,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Complete zl,zp at subdomain:
!
      do lat=lat0,lat1
        do i=lonbeg,lonend
          do k=lev0,lev1
            zli(k,i,lat) = dzdx(k,i,lat)/dlamda*grav*racs(lat)
            zpi(k,i,lat) = dzdy(k,i,lat)/dphi*grav*re_inv
          enddo
          do k=lev0,lev1-1
            zl(k,i,lat) = (zli(k,i,lat)+zli(k+1,i,lat))*.5
            zp(k,i,lat) = (zpi(k,i,lat)+zpi(k+1,i,lat))*.5
          enddo ! k=lev0,lev1-1
          zl(lev1,i,lat) = zli(lev1,i,lat)
          zp(lev1,i,lat) = zpi(lev1,i,lat)
        enddo ! i=lonbeg,lonend
!
! Set periodic points zero to avoid NaNS fpe in duv.
! (earlier versions actually set periodic points here)
!        if (lon0==1) then
!          zl(:,lon0:lon0+1,lat) = 0.
!          zp(:,lon0:lon0+1,lat) = 0.
!        endif
!        if (lon1==nlonp4) then
!          zl(:,lon1-1:lon1,lat) = 0.
!          zp(:,lon1-1:lon1,lat) = 0.
!        endif
      enddo ! lat=lat0,lat1
!
! This part is modified to restore the periodic boundary conditions
! for zl and zp. zl and zp are horizontal pressure gradients. NaNs occurred before
! because values at periodic points (lon=1,2 and 75,76) of tn_upd in dt.F were not given. 
! They are used in subroutine dldp to preform fouth order finite difference. 
! Now the halo longtudinal periodic boundary for tn_upd is set in dt.F, 
! so the forth order differencing in dldp can be made betwen lon=3 to 74 
! (for 5 degree resolution), and so the values of periodic 
! points of zl and zp are now no longer NaNs.
!                                                   Wenbin  04/24/08
      ftmp(:,:,:,1) = zl(:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,2) = zp(:,lon0:lon1,lat0:lat1)

      call mp_periodic_f3d(ftmp,lev0,lev1,lon0,lon1,lat0,lat1,2)

      zl(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,1)
      zp(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,2)
!
      end subroutine glp
!-----------------------------------------------------------------------
      subroutine blktri(a,b,c,f,x,lev0,lev1,lon0,lon1)
!
! This procedure solves (I2-I1+1) tridiagonal block matrix
! systems in which all blocks are 2 x 2 matrices.
!
      use matutil_module,only: matinv2
      implicit none
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1
      real,dimension(2,2,lev0:lev1,lon0:lon1),intent(in) ::
     |  a,b,c
      real,dimension(2,lev0:lev1,lon0:lon1),intent(in) ::
     |  f
!
! Output args:
      real,dimension(2,lev0:lev1,lon0:lon1),intent(out) ::
     |  x
!
! Local:
      integer :: k,i
      real,dimension(2,2,lev0:lev1,lon0:lon1) ::
     |  beta,gama
      real,dimension(2,lev0:lev1,lon0:lon1) :: y
!
! Lower boundary:
      do i=lon0,lon1
!
! BETA(K1) = B(K1)**(-1)
        beta(:,:,lev0,i) = matinv2(b(:,:,lev0,i))
!
! Y(K1) = BETA(K1)*F(K1)
        y(:,lev0,i) = matmul(beta(:,:,lev0,i),f(:,lev0,i))
!
! Now deal with levels (K1+1),K2,1
        do k=lev0+1,lev1-1
!
! GAMMA(K-1) = BETA(K-1)*C(K-1)
          gama(:,:,k-1,i) = matmul(beta(:,:,k-1,i),c(:,:,k-1,i))
!
! GAMMA(K) = B(K) - A(K)*GAMMA(K-1)
          gama(:,:,k,i) = b(:,:,k,i)-matmul(a(:,:,k,i),gama(:,:,k-1,i))
!
! BETA(K) = GAMMA(K)**(-1)
          beta(:,:,k,i) = matinv2(gama(:,:,k,i))
!
!  X(K) = F(K) - A(K)*Y(K-1)
          x(:,k,i) = f(:,k,i)-matmul(a(:,:,k,i),y(:,k-1,i))
!
! Y(K) = BETA(K)*X(K)
          y(:,k,i) = matmul(beta(:,:,k,i),x(:,k,i))
        enddo ! k=lev0+1,lev1-1
      enddo ! i=lon0,lon1
!
      do i=lon0,lon1
!
!  X(K2) = Y(K2)
        x(:,lev1-1,i) = y(:,lev1-1,i)
!
! Backward sweep to determine final solution, X(K) for k=k2,k1,-1
        do k=lev1-2,lev0,-1
          x(:,k,i) = y(:,k,i)-matmul(gama(:,:,k,i),x(:,k+1,i))
        enddo ! k=lev1-2,lev1-1,-1
      enddo ! i=lon0,lon1
      end subroutine blktri
!-----------------------------------------------------------------------
      subroutine filter_uv(fout,lev0,lev1,lon0,lon1,lat0,lat1,name)
!
! Filter updated W omega:
!
      use params_module,only: nlonp4
      use addfld_module,only: addfld
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(inout) :: fout(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
      character(len=*),intent(in) :: name
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j,ier
      real :: fik(nlonp4,lev0:lev1),fkij(lev0:lev1,nlonp4,lat0:lat1)
!
#ifdef VT
!     code = 131 ; state = 'filter_uv' ; activity='Filtering'
      call vtbegin(131,ier)
#endif
!
! Save subdomain of current task in fkij:
      fkij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          fkij(:,i,j) = fout(:,i,j)
        enddo
      enddo ! j=lat0,lat1
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0 
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1,name)
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
      if (mytidi==0) then
!
! Define 2d array with all longitudes for filter at each latitude:
        latscan: do j=lat0,lat1
          do i=1,nlonp4
            fik(i,:) = fkij(:,i,j)
          enddo ! i=1,nlonp4
!
          call ringfilter(fik,lev0,lev1,3,name,j)
!
! Return filtered array to fkij:
          do i=1,nlonp4
            fkij(:,i,j) = fik(i,:)
          enddo ! i=1,nlonp4
        enddo latscan ! j=lat0,lat1
      endif ! mytidi==0
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1,name)

!     write(6,"('filter uv after mp_scatterlons: lon0,1=',2i4,' lat0,1='
!    |  ,2i4,' fkij(:,lon0:lon1,:) min,max=',2e12.4)") lon0,lon1,
!    |  lat0,lat1,minval(fkij(:,lon0:lon1,:)),
!    |            maxval(fkij(:,lon0:lon1,:))
!
! Return filtered array to fout at current task longitudes and latitudes:
      do j=lat0,lat1
        do i=lon0,lon1
          fout(:,i,j) = fkij(:,i,j)
        enddo
      enddo
!
#ifdef VT
!     code = 131 ; state = 'filter_uv' ; activity='Filtering'
      call vtend(131,ier)
#endif
      end subroutine filter_uv
