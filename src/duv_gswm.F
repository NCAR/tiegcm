! modified from subroutine duv to include the option of setting
! the lower boundary condition by the GSWM output (subroutine lbc_gswm_duv) 
!
      subroutine duv_gswm(tn,tn_upd,tn_nm,un,vn,un_nm,vn_nm,w_upd,
     |  barm,z,hdu,hdv,ui,vi,lxx,lyy,lxy,lyx,km,un_upd,unm_upd,vn_upd,
     |  vnm_upd,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Advance neutral velocities for current time step:
!
      use params_module,only: nlonp4,dz,nlat
!     use input_module,only: iuivi ! no longer in use (see dynamo input flag)
      use init_module,only: iter
      use bndry_module,only: ub,ub2,uba,vb,vb2,vba,bnd,bnd2,bnda,ci
      use cons_module,only: freq_semidi,dt,dtx2inv,expz,p0,gask,grav,
     |  xmue,cor,tanphi=>tn,re,kut_5,set_wave_filter,dtsmooth_div2,
     |  dtsmooth
#ifdef MPI
      use mpi_module,only: mp_periodic_f3d
#endif
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
     |  barm,    ! mean molecular weight
     |  z,       ! geopotential height
     |  hdu,     ! horizontal diffusion of U (hdif3 in hdif.F)
     |  hdv,     ! horizontal diffusion of V (hdif3 in hdif.F)
     |  ui,      ! zonal ion velocity input
     |  vi,      ! meridional ion velocity input
     |  lxx,     ! xx ion drag coefficient
     |  lyy,     ! yy ion drag coefficient
     |  lxy,     ! xy ion drag coefficient
     |  lyx,     ! yx ion drag coefficient
     |  km       ! molecular viscosity (cpktkm.F)
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
      integer :: k,i,lat,lonbeg,lonend,kutt(nlat)
      integer :: nk,nkm1,nlevs
      complex :: expt,expt2,expta
      real :: 
     |  unlbc(lon0:lon1), ! un lower boundary (t1)
     |  vnlbc(lon0:lon1)  ! vn lower boundary (t2)
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  g,           ! g**2*(kM+XMUE)*MBAR*/(P0*R*T*DS**2)
     |  dwdz,        ! exp(-s)*w(k+1/2)/(2.*Ds)
     |  tni          ! tn at interfaces
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  advec_un,    ! zonal advection (output of sub advec)
     |  advec_vn,    ! meridional advection (output of sub advec)
     |  zl,zp,       ! output from glp (and dldp)
     |  unm_smooth,  ! un at time n-1, smoothed
     |  vnm_smooth   ! vn at time n-1, smoothed
      real :: rtxmue(lev0:lev1) ! sqrt(xmue)
! For diagnostic, debug:
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  unlbc_diag,  ! un lbc redundant in vertical, for diagnostic ki
     |  vnlbc_diag,  ! vn lbc redundant in vertical, for diagnostic ki
     |  ss_un,ss_vn, ! also for diagnostics
     |  pp11,pp12,pp21,pp22, qq11,qq12,qq21,qq22, rr11,rr12,rr21,rr22
      real,dimension(lon0:lon1,lat0:lat1) ::
     |  unlbc_diag_ij,! un lbc , for diagnostic ij
     |  vnlbc_diag_ij ! vn lbc , for diagnostic ij
!
! pp, qq, rr, beta AND gamma are (2 x 2) matrices defined for 
!   k = 3/2, K-1/2, 1 and i = 3, IMAX+2, 1
      real,dimension(2,2,lev0:lev1,lon0:lon1) :: 
     |  pp,qq,rr,beta,gamma 
!
! ss, xx and yy are two component vectors similarly defined.
      real,dimension(2,lev0:lev1,lon0:lon1) :: 
     |  ss,xx,yy
!
! For addfsech:
      nk = lev1-lev0+1
      nkm1 = nk-1
      nlevs = nk
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Complex constants:
      expt  = cexp(ci*freq_semidi*dt*iter)
      expt2 = cexp(ci*.5*freq_semidi*dt*iter)
      expta = 1.
!
! First latitude scan:
      do lat=lat0,lat1
!
! Horizontal advection:
!     subroutine advec(f,hadvec,lev0,lev1,lon0,lon1,lat)
!     real,dimension(lev0:lev1,lon0-2:lon1+2,lat-2:lat+2),intent(in) :: 
!    |  f ! input field with ghost cells for finite differencing
!     real,dimension(lev0:lev1,lon0:lon1),intent(out) :: hadvec
!      
        call advec(un(:,:,lat-2:lat+2),advec_un(:,:,lat), 
     |    lev0,lev1,lon0,lon1,lat)
        call advec(vn(:,:,lat-2:lat+2),advec_vn(:,:,lat), 
     |    lev0,lev1,lon0,lon1,lat)

!       call addfsech('ADVEC_U0',' ',' ',advec_un(:,:,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!       call addfsech('ADVEC_V0',' ',' ',advec_vn(:,:,lat),
!    |    lon0,lon1,nk,nkm1,lat)

!
! End first latitude scan:
      enddo ! lat=lat0,lat1
!
! Horizontal pressure forcing output of glp is in zl,zp. Pass full 
! subdomains for derivatives in lat and lon (see sub dldp, called by glp):
!
      call glp(tn,tn_upd,tn_nm,barm,z,zl,zp,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Second latitude scan:
      do lat=lat0,lat1

!       call addfsech('ZL',' ',' ',zl(:,:,lat),lon0,lon1,nk,nkm1,lat)
!       call addfsech('ZP',' ',' ',zp(:,:,lat),lon0,lon1,nk,nkm1,lat)

        do i=lon0,lon1
          do k=lev0,lev1-1
            advec_un(k,i,lat) = advec_un(k,i,lat)+
     |        zl(k,i,lat)-hdu(k,i,lat)
            advec_vn(k,i,lat) = advec_vn(k,i,lat)+
     |        zp(k,i,lat)-hdv(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfsech('ADVEC_U1',' ',' ',advec_un(:,:,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!       call addfsech('ADVEC_V1',' ',' ',advec_vn(:,:,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!
! End second lat scan:
      enddo ! lat=lat0,lat1
!
! Smooth un,vn at time n-1 with 2-part shapiro smoother.
! Input un_nm,vn_nm at full subdomain.
! Output in unm_smooth,vnm_smooth at subdomain.
!
      call smooth(un_nm,unm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
      call smooth(vn_nm,vnm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
!
     
! Start third latitude scan:
      do lat=lat0,lat1
      
!       call addfsech('UNSMOOTH',' ',' ',unm_smooth(:,:,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!       call addfsech('VNSMOOTH',' ',' ',vnm_smooth(:,:,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!
! Setting lower boundary condition (check for GSWM input)
        call lbc_gswm_duv(unlbc,vnlbc,unlbc_diag,vnlbc_diag,
     |    lon0,lon1,lat,lev0,lev1,expt,expt2,expta)

!      call addfsech('UNLBC_IK','UNLBC_IK','[cm/s]',unlbc_diag,lon0,
!    |    lon1,nk,nkm1,lat)
!      call addfsech('VNLBC_IK','VNLBC_IK','[cm/s]',vnlbc_diag,lon0,
!    |    lon1,nk,nkm1,lat)
     
	unlbc_diag_ij(:,lat) = unlbc_diag(1,:)
	vnlbc_diag_ij(:,lat) = vnlbc_diag(1,:)
!
        do i=lonbeg,lonend
          do k=lev0,lev1-1
            advec_un(k,i,lat) = advec_un(k,i,lat)-
     |        dtx2inv*unm_smooth(k,i,lat)
            advec_vn(k,i,lat) = advec_vn(k,i,lat)-
     |        dtx2inv*vnm_smooth(k,i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

!       call addfsech('ADVEC_U2',' ',' ',advec_un(:,:,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!       call addfsech('ADVEC_V2',' ',' ',advec_vn(:,:,lat),
!    |    lon0,lon1,nk,nkm1,lat)

        do i=lon0,lon1
          do k=lev0,lev1-1
            ss(1,k,i) = expz(k)*(.25*((lxx(k,i,lat)+lxx(k+1,i,lat))*
     |        (ui(k,i,lat)+ui(k+1,i,lat))+(lxy(k,i,lat)+lxy(k+1,i,lat))*
     |        (vi(k,i,lat)+vi(k+1,i,lat)))-advec_un(k,i,lat))
            ss(2,k,i) = expz(k)*(.25*((lyy(k,i,lat)+lyy(k+1,i,lat))*
     |        (vi(k,i,lat)+vi(k+1,i,lat))-(lyx(k,i,lat)+lyx(k+1,i,lat))*
     |        (ui(k,i,lat)+ui(k+1,i,lat)))-advec_vn(k,i,lat))
            ss_un(k,i) = ss(1,k,i)
            ss_vn(k,i) = ss(2,k,i)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfsech('SS_UN',' ',' ',ss_un,lon0,lon1,nk,nkm1,lat)
!       call addfsech('SS_VN',' ',' ',ss_vn,lon0,lon1,nk,nkm1,lat)

!
! xmue = eddy viscosity (cons.F). rtxmue = sqrt(xmue)
        do k=lev0+1,lev1-1          
          rtxmue(k) = sqrt(xmue(k-1)*xmue(k))
        enddo ! k=lev0+1,lev1-1          
!
! tni = tn at interfaces:
        do i=lon0,lon1
          do k=lev0+1,lev1-1          
            tni(k,i) = .5*(tn(k-1,i,lat)+tn(k,i,lat))
          enddo ! k=lev0+1,lev1-1          
        enddo ! i=lon0,lon1
!
! Boundaries:
        rtxmue(lev0) = sqrt(xmue(lev0)**3/xmue(lev0+1))
        rtxmue(lev1) = sqrt(xmue(lev1-1)**3/xmue(lev1-2))
!       write(6,"('duv: lat=',i2,' rtxmue=',/,(6e12.4))") lat,rtxmue
        do i=lon0,lon1
          tni(lev0,i) = tn(lev1,i,lat) ! lbc of tn is stored in top slot
          tni(lev1,i) = tn(lev1-1,i,lat)
!
! G = g*KM/(P0*H*Ds**2) = g**2*(kM+XMUE)*MBAR*/(P0*R*T*DS**2)
          do k=lev0,lev1
            g(k,i) = grav**2*(km(k,i,lat)+rtxmue(k))*barm(k,i,lat)/
     |        (p0*gask*tni(k,i)*dz**2)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

!       call addfsech('TNI'  ,' ',' ',tni,lon0,lon1,nk,nk,lat)
!       call addfsech('DUV_G',' ',' ',g  ,lon0,lon1,nk,nk,lat)
!
! dwdz = exp(-s)*w(k+1/2)/(2.*Ds)
        do i=lon0,lon1
          do k=lev0,lev1-1
            dwdz(k,i) = 0.5*expz(k)*(w_upd(k,i,lat)+w_upd(k+1,i,lat))/
     |        (2.*dz)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfsech('DWDZ',' ',' ',dwdz,lon0,lon1,nk,nkm1,lat)
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
            rr(1,1,k,i) = -(g(k+1,i)-dwdz(k,i))
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
            qq(1,1,k,i) = g(k,i) + g(k+1,i)
            qq(2,2,k,i) = qq(1,1,k,i)
            qq(1,1,k,i) = qq(1,1,k,i) + 
     |        expz(k)*(dtx2inv+.5*(lxx(k,i,lat)+lxx(k+1,i,lat)))
            qq(2,2,k,i) = qq(2,2,k,i) +
     |        expz(k)*(dtx2inv+.5*(lyy(k,i,lat)+lyy(k+1,i,lat)))
!
!     Q(1,2) = -exp(-s) * (f+u/r*tan(theta)-lamxy)
!      and 
!     Q(2,1) = exp(-s) * (f+u/r*tan(theta)-lamyx)
!
            qq(1,2,k,i) = cor(lat)+un(k,i,lat)/re*tanphi(lat)
            qq(2,1,k,i) = qq(1,2,k,i)
            qq(1,2,k,i) = -expz(k)*(qq(1,2,k,i)-
     |        .5*(lxy(k,i,lat)+lxy(k+1,i,lat)))
            qq(2,1,k,i) =  expz(k)*(qq(2,1,k,i)-
     |        .5*(lyx(k,i,lat)+lyx(k+1,i,lat)))
          enddo ! k=lev0,lev1-1
!
! Boundaries at k=3/2 and k=k-1/2
! Lower boundary:
!     Q(3/2) = Q(3/2) - P(3/2)
!     S(3/2) = S(3/2) -2.*P(3/2)*Vb
!     P(3/2) = 0.
!
          qq(1,1,lev0,i) = qq(1,1,lev0,i)-pp(1,1,lev0,i)
          qq(1,2,lev0,i) = qq(1,2,lev0,i)-pp(1,2,lev0,i)
          qq(2,1,lev0,i) = qq(2,1,lev0,i)-pp(1,2,lev0,i)
          qq(2,2,lev0,i) = qq(2,2,lev0,i)-pp(2,2,lev0,i)
          ss(1,lev0,i) = ss(1,lev0,i)-2.*(pp(1,1,lev0,i)*unlbc(i)+
     |                                    pp(1,2,lev0,i)*vnlbc(i))
          ss(2,lev0,i) = ss(2,lev0,i)-2.*(pp(2,1,lev0,i)*unlbc(i)+
     |                                    pp(2,2,lev0,i)*vnlbc(i))
          ss_un(lev0,i) = ss(1,lev0,i)
          ss_vn(lev0,i) = ss(2,lev0,i)
          pp(1,1,lev0,i) = 0.
          pp(1,2,lev0,i) = 0.
          pp(2,1,lev0,i) = 0.
          pp(2,2,lev0,i) = 0.
!
! Upper boundary:
!     Q(K-1/2) = Q(K-1/2) + R(K-1/2)
!     R(K-1/2) = 0.
!
          qq(1,1,lev1-1,i) = qq(1,1,lev1-1,i)+rr(1,1,lev1-1,i)
          qq(1,2,lev1-1,i) = qq(1,2,lev1-1,i)+rr(1,2,lev1-1,i)
          qq(2,1,lev1-1,i) = qq(2,1,lev1-1,i)+rr(2,1,lev1-1,i)
          qq(2,2,lev1-1,i) = qq(2,2,lev1-1,i)+rr(2,2,lev1-1,i)
          rr(1,1,lev1-1,i) = 0.
          rr(1,2,lev1-1,i) = 0.
          rr(2,1,lev1-1,i) = 0.
          rr(2,2,lev1-1,i) = 0.
        enddo ! i=lon0,lon1

!       pp11(:,:) = pp(1,1,:,:) ; pp12(:,:) = pp(1,2,:,:)
!       pp21(:,:) = pp(2,1,:,:) ; pp22(:,:) = pp(2,2,:,:)
!       qq11(:,:) = qq(1,1,:,:) ; qq12(:,:) = qq(1,2,:,:)
!       qq21(:,:) = qq(2,1,:,:) ; qq22(:,:) = qq(2,2,:,:)
!       rr11(:,:) = rr(1,1,:,:) ; rr12(:,:) = rr(1,2,:,:)
!       rr21(:,:) = rr(2,1,:,:) ; rr22(:,:) = rr(2,2,:,:)
!       call addfsech('PP11',' ',' ',pp11,lon0,lon1,nk,nkm1,lat)
!       call addfsech('PP12',' ',' ',pp12,lon0,lon1,nk,nkm1,lat)
!       call addfsech('PP21',' ',' ',pp21,lon0,lon1,nk,nkm1,lat)
!       call addfsech('PP22',' ',' ',pp22,lon0,lon1,nk,nkm1,lat)
!       call addfsech('QQ11',' ',' ',qq11,lon0,lon1,nk,nkm1,lat)
!       call addfsech('QQ12',' ',' ',qq12,lon0,lon1,nk,nkm1,lat)
!       call addfsech('QQ21',' ',' ',qq21,lon0,lon1,nk,nkm1,lat)
!       call addfsech('QQ22',' ',' ',qq22,lon0,lon1,nk,nkm1,lat)
!       call addfsech('RR11',' ',' ',rr11,lon0,lon1,nk,nkm1,lat)
!       call addfsech('RR12',' ',' ',rr12,lon0,lon1,nk,nkm1,lat)
!       call addfsech('RR21',' ',' ',rr21,lon0,lon1,nk,nkm1,lat)
!       call addfsech('RR22',' ',' ',rr22,lon0,lon1,nk,nkm1,lat)
!       call addfsech('SS_UN',' ',' ',ss_un,lon0,lon1,nk,nkm1,lat)
!       call addfsech('SS_VN',' ',' ',ss_vn,lon0,lon1,nk,nkm1,lat)
!
! Invoke tridiagonal (2x2) block matrix solver for U and V: 
! pp,qq,rr,ss are input, beta,gamma,yy,xx are output.
! (final solution is in xx)
!
        call blktri(pp,qq,rr,ss, beta,gamma,yy,xx,
     |    lev0,lev1,lon0,lon1,lat)
!
! Save updated U and V:
        do i=lonbeg,lonend 
          do k=lev0,lev1-1
            un_upd(k,i,lat) = xx(1,k,i)
            vn_upd(k,i,lat) = xx(2,k,i)
          enddo ! k=lev0,lev1-1
!
! Store lower boundary conditions for U and V in top slot:
          un_upd(lev1,i,lat) = unlbc(i)
          vn_upd(lev1,i,lat) = vnlbc(i)

!         write(6,"('duv: i=',i3,' lat=',i3,' unlbc=',e12.4,' vnlbc=',
!    |      e12.4)") i,lat,unlbc(i),vnlbc(i)

        enddo ! i=lonbeg,lonend 

!       call addfsech('UNLBC',' ',' ',unlbc_diag,lon0,lon1,nk,nkm1,lat)
!       call addfsech('VNLBC',' ',' ',vnlbc_diag,lon0,lon1,nk,nkm1,lat)

!       call addfsech('UN_SOLV',' ',' ',un_upd(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nk,lat)
!       call addfsech('VN_SOLV',' ',' ',vn_upd(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nk,lat)
!
! End third latitude scan:
      enddo ! lat=lat0,lat1

!      call addfsech_ij('UNLBC_IJ','UNLBC_IJ','[cm/s]',unlbc_diag_ij,
!    |    lon0,lon1,lat0,lat1)
!      call addfsech_ij('VNLBC_IJ','VNLBC_IJ','[cm/s]',vnlbc_diag_ij,
!    |    lon0,lon1,lat0,lat1)
!
! Set kut for wave filtering according to dlat (2.5 or 5.0):
      call set_wave_filter(36,kut_5,nlat,kutt)
!
! Fourier smoothing of U and V:
      call filter_uv(un_upd(lev0:lev1-1,:,:),lev0,lev1-1,lon0,lon1,
     |  lat0,lat1,kutt)
      call filter_uv(vn_upd(lev0:lev1-1,:,:),lev0,lev1-1,lon0,lon1,
     |  lat0,lat1,kutt)
!
! Fourth latitude scan:
      do lat=lat0,lat1

!       call addfsech('UN_FILT',' ',' ',un_upd(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nk,lat)
!       call addfsech('VN_FILT',' ',' ',vn_upd(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nk,lat)

!
! Store lower boundary conditions for U and V in top slot:
!       do i=lon0,lon1 
!         un_upd(lev1,i,lat) = unlbc(i)
!         vn_upd(lev1,i,lat) = vnlbc(i)
!       enddo ! i=lon0,lon1 
!
! Smooth updated un,vn:
        do i=lon0,lon1 
          do k=lev0,lev1-1
            unm_upd(k,i,lat) = dtsmooth_div2*(un_nm(k,i,lat)+
     |        un_upd(k,i,lat)) + dtsmooth*un(k,i,lat)
            vnm_upd(k,i,lat) = dtsmooth_div2*(vn_nm(k,i,lat)+
     |        vn_upd(k,i,lat)) + dtsmooth*vn(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1 

!       call addfsech('UN_NMOUT',' ',' ',unm_upd(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!       call addfsech('VN_NMOUT',' ',' ',vnm_upd(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!
! End fourth latitude scan:
      enddo ! lat=lat0,lat1
!
#ifdef MPI
!
! Periodic points:
!     call mp_periodic_f3d(un_upd(:,lon0:lon1,lat0-1:lat1+1),
!    |  lev0,lev1,lon0,lon1,lat0-1,lat1+1)
!     call mp_periodic_f3d(vn_upd(:,lon0:lon1,lat0-1:lat1+1),
!    |  lev0,lev1,lon0,lon1,lat0-1,lat1+1)
#endif
!     do lat=lat0,lat1
!       call addfsech('UN_FINAL',' ',' ',un_upd(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nk,lat)
!       call addfsech('VN_FINAL',' ',' ',vn_upd(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nk,lat)
!     enddo ! lat=lat0,lat1
      end subroutine duv_gswm
