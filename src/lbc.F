
      module lbc
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate lower boundary conditions for T,U,V,Z
!
      use params_module,only: nlonp4,nlat,nlev,dz
      use cons_module,only: pi,atm_amu,gask,grav,freq_semidi,
     |  dt,re,dlamda,tgrad,cs,cor,tn
      use cons_module,only: tbound,zbound
      use addfld_module,only: addfld
      implicit none
!
! Total lower boundary conditions returned by this module.
! (dimensioned at full global grid, but defined at subdomains only)
!
      real,dimension(nlonp4,nlat) :: t_lbc, u_lbc, v_lbc, z_lbc
!
! Diurnal and semi-diurnal tidal perturbations using Hough functions: 
!
      complex,dimension(nlat) :: 
     |  t_di , u_di , v_di , z_di, t_sdi, u_sdi, v_sdi, z_sdi
      complex,parameter :: ci=(0.,1.), expta=1.
      complex :: bnd_sdi(nlonp4), bnd_di(nlonp4)
!
! For bndcmp:
      real :: b(nlonp4,3,3),fb(nlonp4,3)
!
! This t0 is different than the t0 in cons.
      real :: t0(nlev+1) 
!
! Lower boundary for helium (mmr):
      real,parameter :: pshelb=0.1154E-5
      contains
!-----------------------------------------------------------------------
      subroutine init_lbc
!
! Called once per run from tgcm.
!
      use cons_module,only: t0cons=>t0
      implicit none
!
! t0 in cons is zeroed out:
      t0cons(:) = 0.
!
! t0 local to this module (different than the t0 in cons) is for
! use in Hough functions:
      t0(:) = 0.
      t0(1) = tbound
      t0(2) = tbound+dz*tgrad
!
      call bndry_diurnal      ! t_di ,u_di ,v_di ,z_di
      call bndry_semidiurnal  ! t_sdi,u_sdi,v_sdi,z_sdi
      call bndcmp
      end subroutine init_lbc
!-----------------------------------------------------------------------
      subroutine tuvz_lbc(istep,iyear,iday,secs,itp)
!
! Update lower boundary subdomains of T,U,V,Z
! This is called every timestep from advance.
!
      use mpi_module,only: lon0,lon1,lat0,lat1
      use init_module,only: iter,
     |  ibgrddata           ! 0/1 flag to get background data
      use input_module,only: step,
     |  gswm_mi_di_ncfile,  ! gswm migrating diurnal data file
     |  gswm_mi_sdi_ncfile, ! gswm migrating semi-diurnal data file
     |  gswm_nm_di_ncfile,  ! gswm non-migrating diurnal data file
     |  gswm_nm_sdi_ncfile  ! gswm non-migrating semi-diurnal data file
      use gswm_module,only: getgswm, ! (nlonp4,nlat)
     |  gswm_mi_di_z, gswm_mi_sdi_z, gswm_nm_di_z, gswm_nm_sdi_z,
     |  gswm_mi_di_t, gswm_mi_sdi_t, gswm_nm_di_t, gswm_nm_sdi_t,
     |  gswm_mi_di_u, gswm_mi_sdi_u, gswm_nm_di_u, gswm_nm_sdi_u,
     |  gswm_mi_di_v, gswm_mi_sdi_v, gswm_nm_di_v, gswm_nm_sdi_v
      use bgrd_data_module,only: get_bgrddata,! (nlonp4,nlat)
     |  bgrddata_z, bgrddata_t, bgrddata_u, bgrddata_v
      use input_module,only: saber_ncfile, tidi_ncfile, ctmt_ncfile
      use saber_tidi,only: get_saber_tidi,saber_t,saber_z,
     |  tidi_u,tidi_v
      use ctmt_module,only: get_ctmt,ctmt_t,ctmt_z,ctmt_u,ctmt_v
      use hist_module,only: modeltime

      implicit none
!
! Args:
      integer,intent(in) :: istep,iyear,iday,itp
      real,intent(in) :: secs
!
! Local:
      integer :: i,j
      real :: rstep
      complex :: t_expt_sdi, t_expt_di, uvz_expt_sdi, uvz_expt_di
!
! Calculate exponentials
      rstep = float(step)
      t_expt_sdi = cexp(ci*freq_semidi*rstep*iter)
      t_expt_di  = cexp(ci*.5*freq_semidi*rstep*iter)
      uvz_expt_sdi = cexp(ci*freq_semidi*dt*iter)
      uvz_expt_di  = cexp(ci*.5*freq_semidi*dt*iter)
!
! Set background constants 
      if (ibgrddata > 0 ) then ! background based on data
        call get_bgrddata(istep,iyear,iday,secs)
        t_lbc(lon0:lon1,lat0:lat1) = bgrddata_t(lon0:lon1,lat0:lat1)
        u_lbc(lon0:lon1,lat0:lat1) = bgrddata_u(lon0:lon1,lat0:lat1)
        v_lbc(lon0:lon1,lat0:lat1) = bgrddata_v(lon0:lon1,lat0:lat1)
        z_lbc(lon0:lon1,lat0:lat1) = bgrddata_z(lon0:lon1,lat0:lat1)
      else !(default see cons module):
        t_lbc(lon0:lon1,lat0:lat1) = tbound
        u_lbc(lon0:lon1,lat0:lat1) = 0.
        v_lbc(lon0:lon1,lat0:lat1) = 0.
        z_lbc(lon0:lon1,lat0:lat1) = zbound
      endif

      if (len_trim(saber_ncfile) > 0 .or.
     |    len_trim(tidi_ncfile) > 0 .or.
     |    len_trim(ctmt_ncfile) > 0) goto 100
!
! Add gswm perturbations:
! 12/2/15 btf: moved from advance.F to lbc.F
!
      call getgswm(istep,iyear,iday,secs)
!
! Add gswm or Hough mode perturbations:
!
! GSWM migrating diurnal:
      if (len_trim(gswm_mi_di_ncfile) > 0) then
        do j=lat0,lat1
          do i=lon0,lon1
            t_lbc(i,j) = t_lbc(i,j)+gswm_mi_di_t(i,j)
            u_lbc(i,j) = u_lbc(i,j)+gswm_mi_di_u(i,j)
            v_lbc(i,j) = v_lbc(i,j)+gswm_mi_di_v(i,j)
            z_lbc(i,j) = z_lbc(i,j)+gswm_mi_di_z(i,j)
          enddo
        enddo
!
! Hough mode diurnal:
      else ! use Hough functions for diurnal tide
        do j=lat0,lat1
          do i=lon0,lon1
            t_lbc(i,j) = t_lbc(i,j)+real(t_di(j)*bnd_di(i)*t_expt_di)
            u_lbc(i,j) = u_lbc(i,j)+real(u_di(j)*bnd_di(i)*uvz_expt_di)
            v_lbc(i,j) = v_lbc(i,j)+real(v_di(j)*bnd_di(i)*uvz_expt_di)
            z_lbc(i,j) = z_lbc(i,j)+real(z_di(j)*bnd_di(i)*uvz_expt_di)
          enddo
        enddo
      endif
!
! GSWM migrating semi-diurnal:
      if (len_trim(gswm_mi_sdi_ncfile) > 0) then
        do j=lat0,lat1
          do i=lon0,lon1
            t_lbc(i,j) = t_lbc(i,j)+gswm_mi_sdi_t(i,j)
            u_lbc(i,j) = u_lbc(i,j)+gswm_mi_sdi_u(i,j)
            v_lbc(i,j) = v_lbc(i,j)+gswm_mi_sdi_v(i,j)
            z_lbc(i,j) = z_lbc(i,j)+gswm_mi_sdi_z(i,j)
          enddo
        enddo
!
! Hough mode semi-diurnal:
      else ! use Hough functions for semi-diurnal tide
        do j=lat0,lat1
          do i=lon0,lon1
            t_lbc(i,j)=t_lbc(i,j)+real(t_sdi(j)*bnd_sdi(i)*t_expt_sdi)
            u_lbc(i,j)=u_lbc(i,j)+real(u_sdi(j)*bnd_sdi(i)*uvz_expt_sdi)
            v_lbc(i,j)=v_lbc(i,j)+real(v_sdi(j)*bnd_sdi(i)*uvz_expt_sdi)
            z_lbc(i,j)=z_lbc(i,j)+real(z_sdi(j)*bnd_sdi(i)*uvz_expt_sdi)
          enddo
        enddo
      endif
!
! GSWM non-migrating diurnal:
      if (len_trim(gswm_nm_di_ncfile) > 0) then
        do j=lat0,lat1
          do i=lon0,lon1
            t_lbc(i,j) = t_lbc(i,j)+gswm_nm_di_t(i,j)
            u_lbc(i,j) = u_lbc(i,j)+gswm_nm_di_u(i,j)
            v_lbc(i,j) = v_lbc(i,j)+gswm_nm_di_v(i,j)
            z_lbc(i,j) = z_lbc(i,j)+gswm_nm_di_z(i,j)
          enddo
        enddo
      endif
!
! GSWM non-migrating semi-diurnal:
      if (len_trim(gswm_nm_sdi_ncfile) > 0) then
        do j=lat0,lat1
          do i=lon0,lon1
            t_lbc(i,j) = t_lbc(i,j)+gswm_nm_sdi_t(i,j)
            u_lbc(i,j) = u_lbc(i,j)+gswm_nm_sdi_u(i,j)
            v_lbc(i,j) = v_lbc(i,j)+gswm_nm_sdi_v(i,j)
            z_lbc(i,j) = z_lbc(i,j)+gswm_nm_sdi_z(i,j)
          enddo
        enddo
      endif
!
! Add SABER and/or TIDI perturbations:
 100  continue
      if (len_trim(saber_ncfile) > 0 .or.
     |    len_trim(tidi_ncfile) > 0)
     |  call get_saber_tidi(modeltime)

      if (len_trim(saber_ncfile) > 0) then
        do j=lat0,lat1
          do i=lon0,lon1
            t_lbc(i,j) = t_lbc(i,j)+saber_t(i,j)
            z_lbc(i,j) = z_lbc(i,j)+saber_z(i,j)
          enddo
        enddo

!       write(6,"('tuvz_lbc: added saber_t,z to t_lbc,z_lbc')")
!       write(6,"('tuvz_lbc: saber_z min,max=',2e12.4)") 
!    |    minval(saber_z),maxval(saber_z)
!       write(6,"('tuvz_lbc: z_lbc min,max=',2e12.4)") 
!    |    minval(z_lbc),maxval(z_lbc)

      endif
      if (len_trim(tidi_ncfile) > 0) then
        do j=lat0,lat1
          do i=lon0,lon1
            u_lbc(i,j) = u_lbc(i,j)+tidi_u(i,j)
            v_lbc(i,j) = v_lbc(i,j)+tidi_v(i,j)
          enddo
        enddo
!       write(6,"('tuvz_lbc: added tidi_u,v to u_lbc,v_lbc')")
      endif
!
! Add CTMT perturbations:
      if (len_trim(ctmt_ncfile) > 0 )
     |  call get_ctmt(istep,iyear,iday,secs)
!
      if (len_trim(ctmt_ncfile) > 0) then
        do j=lat0,lat1
          do i=lon0,lon1
            t_lbc(i,j) = t_lbc(i,j)+ctmt_t(i,j)
            z_lbc(i,j) = z_lbc(i,j)+ctmt_z(i,j)
            u_lbc(i,j) = u_lbc(i,j)+ctmt_u(i,j)
            v_lbc(i,j) = v_lbc(i,j)+ctmt_v(i,j)
          enddo
        enddo
      endif
!
! Save 2d boundaries to secondary histories:
!     call addfld('T_LBC','T_LBC',' ',t_lbc(lon0:lon1,lat0:lat1),
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('U_LBC','U_LBC',' ',u_lbc(lon0:lon1,lat0:lat1),
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('V_LBC','V_LBC',' ',v_lbc(lon0:lon1,lat0:lat1),
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('Z_LBC','Z_LBC',' ',z_lbc(lon0:lon1,lat0:lat1),
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
      end subroutine tuvz_lbc
!-----------------------------------------------------------------------
      subroutine bndry_semidiurnal
!
! Lower boundary conditions for semi-diurnal tide, using Hough functions.
! This is called once per run from init, and returns t_sdi, u_sdi, v_sdi, 
!   z_sdi at nlat latitudes. 
!
      use input_module,only: tide
      implicit none
!
! Local:
      integer,parameter :: nalf=19, malf=2
      real :: p(nlat,nalf,malf),hough(nlat,5,malf),cp(nalf/2+1)
      complex :: dzb(nlat)
      real :: b(5,19),rl(5),bhour(5),rlamda,xdot(19),ydot(19),
     |  ptscal,theta,ptjm(2*nlat+1)
      integer :: n,jm,l,lm1,m,mm1,j,ld,i,nm1
!
      complex zee(5),cl(5),expdlm
      data b/
     | 0.969152, 0.0     , 0.216046, 0.0     , 0.093838,
     | 0.0     , 0.909763, 0.0     , 0.342113, 0.0     ,
     |-0.245226, 0.0     , 0.798445, 0.0     , 0.421218,
     | 0.0     ,-0.408934, 0.0     , 0.645517, 0.0     ,
     | 0.024633, 0.0     ,-0.543993, 0.0     , 0.464159,
     | 0.0     , 0.071127, 0.0     ,-0.643189, 0.0     ,
     |-0.001292, 0.0     , 0.139613, 0.0     ,-0.699495,
     | 0.0     ,-0.006673, 0.0     , 0.225090, 0.0     ,
     | 0.000042, 0.0     ,-0.019654, 0.0     , 0.320141,
     | 0.0     , 0.000394, 0.0     ,-0.043345, 0.0     ,
     |-0.000001, 0.0     , 0.001772, 0.0     ,-0.079831,
     | 0.0     ,-0.000016, 0.0     , 0.005401, 0.0     ,
     | 0.0     , 0.0     ,-0.000112, 0.0     , 0.012932,
     | 0.0     , 0.0     , 0.0     ,-0.000476, 0.0     ,
     | 0.0     , 0.0     , 0.000005, 0.0     ,-0.001490,
     | 0.0     , 0.0     , 0.0     , 0.000031, 0.0     ,
     | 0.0     , 0.0     , 0.0     , 0.0     , 0.000129,
     | 0.0     , 0.0     , 0.0     ,-0.000002, 0.0     ,
     | 0.0     , 0.0     , 0.0     , 0.0     ,-0.000009/
      data rl/7.8519E5, 3.6665E5, 2.1098E5, 1.3671E5, 0.9565E5/
      real,external :: sddot ! util.F
!
! t0 is local to this module (different than the t0 in cons)
      t0(:) = 0.
      t0(1) = tbound
      t0(2) = tbound+dz*tgrad
!
! Longitudinal structure:
      rlamda = -2.*dlamda
      bnd_sdi(1)=cexp(ci*2.*rlamda)
      expdlm=cexp(ci*2.*dlamda)
      do i=2,nlonp4
        bnd_sdi(i)=bnd_sdi(i-1)*expdlm
      enddo
!
! Zero out if user did not provide amp/phase:
      if (all(tide==0.)) then
        t_sdi = 0.
        u_sdi = 0.
        v_sdi = 0.
        z_sdi = 0.
        return
      endif
!
      bhour = tide(6:10)
      do n=1,5
        zee(n)=tide(n)*cexp(ci*pi*bhour(n)/6.)
        cl(n)=csqrt(cmplx(gask/(atm_amu*grav*rl(n))*
     |    (t0(1)*2./7.+(t0(2)-t0(1))/dz)-.25))-.5*ci
      enddo
      jm=2*nlat+1
!
! Set up hough functions (see sphpac.F)
      do n=2,nalf+1
        nm1 = n-1
        do m=2,malf+1
          mm1=m-1
          call alfk(n,m,cp)
          do j=1,jm
            theta = float(j-1)*pi/float(jm-1)
            call lfpt(n,m,theta,cp,ptscal)
            ptjm(j) = ptscal
          enddo
          do j=1,nlat
            p(j,nm1,mm1) = ptjm(2*(nlat+1-j))
          enddo
        enddo
        do j=1,nlat
          p(j,nm1,2)=sqrt(float(n*(n+1)-6))*p(j,nm1,2)-2.*tn(j)*
     |               p(j,nm1,1)
        enddo
      enddo
!
! util.F: real function sddot(n,x,y)
      do l=1,5
        do ld=1,2
          do j=1,nlat
            xdot(:) = p(j,:,ld)
            ydot(:) = b(l,:)
            hough(j,l,ld)=sddot(19,xdot,ydot)
          enddo
        enddo
      enddo
!
! Define module data:
      do j=1,nlat
        t_sdi(j)=0.
        z_sdi(j)=0.
        dzb(j)=0.
      enddo
      do l=1,5
        do j=1,nlat
          z_sdi(j)=z_sdi(j)+zee(l)*hough(j,l,1)
          dzb(j)=dzb(j)+zee(l)*hough(j,l,2)
          t_sdi(j)=t_sdi(j)+ci*atm_amu*grav/gask*zee(l)*cl(l)*
     |      hough(j,l,1)
        enddo
      enddo
      do j=1,nlat
        u_sdi(j)=freq_semidi*re*(1.-(cor(j)/freq_semidi)**2)
        v_sdi(j)=ci*grav*(dzb(j)-2.*cor(j)/(freq_semidi*cs(j))*
     |    z_sdi(j))/u_sdi(j)
        u_sdi(j)=grav*(cor(j)/freq_semidi*dzb(j)-2./cs(j)*
     |    z_sdi(j))/u_sdi(j)
      enddo
!     write(6,"('bndry_semidiurnal: t_sdi min,max=',2e12.4)")
!    |  minval(real(t_sdi)),maxval(real(t_sdi))
!     write(6,"('bndry_semidiurnal: u_sdi min,max=',2e12.4)")
!    |  minval(real(u_sdi)),maxval(real(u_sdi))
!     write(6,"('bndry_semidiurnal: v_sdi min,max=',2e12.4)")
!    |  minval(real(v_sdi)),maxval(real(v_sdi))
!     write(6,"('bndry_semidiurnal: z_sdi min,max=',2e12.4)")
!    |  minval(real(z_sdi)),maxval(real(z_sdi))
      end subroutine bndry_semidiurnal
!-----------------------------------------------------------------------
      subroutine bndry_diurnal
!
! Lower boundary conditions for diurnal tide, using Hough functions.
! This is called once per run from init, and returns t_di, u_di, v_di, 
!   z_di at nlat latitudes. 
!
      use input_module,only: tide2
!
! Local:
      integer,parameter :: nalf=19, malf=2
      real :: p(nlat,nalf,malf),hough(nlat,5,malf),cp(nalf/2+1)
      complex :: dzb(nlat)
      real :: b(1,19),rl (1),bhour(1),rlamda,xdot(19),ydot(19),
     |  ptscal,theta,ptjm(2*nlat+1),pik
      integer :: l,m,j,n,jm,ld,i
      complex zee(1),cl(1),expdlm
!
      data b/
     | 0.282710,
     | 0.0     ,
     |-0.638229,
     | 0.0     ,
     | 0.620521,
     | 0.0     ,
     |-0.336408,
     | 0.0     ,
     | 0.117021,
     | 0.0     ,
     |-0.028332,
     | 0.0     ,
     | 0.005042,
     | 0.0     ,
     |-0.000686,
     | 0.0     ,
     | 0.000074,
     | 0.0     ,
     |-0.000006/
      data rl/0.6909E5/
      real,external :: sddot ! in util.F
!
! t0 is local to this module (different than the t0 in cons)
      t0(:) = 0.
      t0(1) = tbound
      t0(2) = tbound+dz*tgrad
!
! Calculate longitudinal structure
      rlamda = -2.*dlamda
      bnd_di(1)=cexp(ci*rlamda)
      expdlm=cexp(ci*dlamda)
      do i=2,nlonp4
        bnd_di(i)=bnd_di(i-1)*expdlm
      enddo
!
! Zero out if user did not provide amp/phase:
      if (all(tide2==0.)) then
        t_di = 0.
        u_di = 0.
        v_di = 0.
        z_di = 0.
        return
      endif
      bhour(1) = tide2(2)
      pik = 3.14159265358979312 
      do n=1,1
        zee(n)=tide2(n)*cexp(ci*pik*bhour(n)/12.)
        cl(n)=csqrt(cmplx(gask/(atm_amu*grav*rl(n))*
     |  (t0(1)*2./7.+(t0(2)-t0(1))/dz)-.25))-.5*ci
      enddo
      jm=2*nlat+1
!
! Set up hough functions:
!
      do n=1,19
        do m=1,2
          call alfk(n,m,cp)
          do j=1,jm
            theta = float(j-1)*pi/float(jm-1)
            call lfpt(n,m,theta,cp,ptscal)
            ptjm(j)=ptscal
          enddo
          do j=1,nlat
            p(j,n,m) = ptjm(2*(nlat+1-j))
          enddo
        enddo
        do j=1,nlat
          p(j,n,2)=sqrt(float(n*(n+1)-2))*p(j,n,2)-tn(j)*p(j,n,1)
        enddo
      enddo
!
! util.F: real function sddot(n,x,y)
      do l=1,1
        do ld=1,2
          do j=1,nlat
            xdot(:) = p(j,:,ld)
            ydot(:) = b(l,:)
            hough(j,l,ld)=sddot(19,xdot,ydot)
          enddo
        enddo
      enddo
!
! Generate t_di, u_di, v_di, z_di:
      do j=1,nlat
        t_di(j)=0.
        z_di(j)=0.
        dzb(j)=0.
      enddo
      do l=1,1
        do j=1,nlat
          z_di(j)=z_di(j)+zee(l)*hough(j,l,1)
          dzb(j)=dzb(j)+zee(l)*hough(j,l,2)
          t_di(j)=t_di(j)+ci*atm_amu*grav/gask*zee(l)*cl(l)*hough(j,l,1)
        enddo
      enddo
      do j=1,nlat
        u_di(j)=.5*freq_semidi*re*(1.-(cor(j)/(.5*freq_semidi))**2)
        v_di(j)=ci*grav*(dzb(j)-cor(j)/(.5*freq_semidi*cs(j))*z_di(j))/
     |    u_di(j)
        u_di(j)=grav*(cor(j)/(.5*freq_semidi)*dzb(j)-1./cs(j)*z_di(j))/
     |    u_di(j)
      enddo
!     write(6,"('bndry_diurnal: t_di min,max=',2e12.4)")
!    |  minval(real(t_di)),maxval(real(t_di))
!     write(6,"('bndry_diurnal: u_di min,max=',2e12.4)")
!    |  minval(real(u_di)),maxval(real(u_di))
!     write(6,"('bndry_diurnal: v_di min,max=',2e12.4)")
!    |  minval(real(v_di)),maxval(real(v_di))
!     write(6,"('bndry_diurnal: z_di min,max=',2e12.4)")
!    |  minval(real(z_di)),maxval(real(z_di))
      end subroutine bndry_diurnal
!-----------------------------------------------------------------------
      subroutine bndcmp
C     ****
C     ****     CALCULATE MATRICES B(nlonp4,2,2) AND VECTORS FB(nlonp4,2)
C     ****       REPRESENTING THE LOWER BOUBNDARY CONDITION IN COMP,
C     ****       WHERE PSI1 AND PSI2 ARE CALCULATED:
C     ****
C     ****         PSI(K=-1/2) = B * PSI(K=1/2) + FB
C     ****
C     ****     BNDCMP CALLS THE SUBROUTINE BNDEF TO DEFINE THE 2 X 2
C     ****       MATRICES E, F AND THE 2 VECTOR G IN THE GENERAL
C     ****       LOWER BOUNDARY CONDITION:
C     ****
C     ****         E * D(PSI)/DS + F * PSI + G = 0.
C     ****
C     ****         WHERE:
C     ****           PSI = |PSI1| AND THE BOUNDARY CONDITION IS APPLIED
C     ****                 |    |
C     ****                 |PSI2|
C     ****
C     ****           AT LEVEL ZERO
C     ****
C     ****     THIS SUBROUTINE THEN EVALUATES B AND FB FROM:
C     ****
C     ****       B = (E/DS - F/2.)**(-1) * (E/DS + F/2.)
C     ****
C     ****       FB = (E/DS - F/2.)**(-1) * G
C     ****
!
! Local:
      real :: ee(nlonp4,3,3),ff(nlonp4,3,3),gg(nlonp4,3),
     |  wm1(nlonp4,3,3),wm2(nlonp4,3,3),wm3(nlonp4,3,3),
     |  ws1(nlonp4)
      integer :: l,m,i
C     ****
C     ****     CALL BNDEF TO DEFINE E, F AND G IN S1, S2 AND S3
C     ****
      call bndef(ee,ff,gg)
C     ****
C     ****     WM1 = (E/DS - F/2.)
C     ****
C     ****     WM2 = (E/DS + F/2.)
C     ****
      do l = 1,3
        do m = 1,3
          do i = 1,nlonp4
            wm1(i,l,m) = ee(i,l,m)/dz-ff(i,l,m)/2.
            wm2(i,l,m) = ee(i,l,m)/dz+ff(i,l,m)/2.
          enddo
        enddo
      enddo
C     ****
C     ****     WM3 = WM1**(-1)
C     ****
C     ****     WS1 = DET(WM1)
C     ****
      do i = 1,nlonp4
        ws1(i) = 
     |    wm1(i,1,1)*(wm1(i,2,2)*wm1(i,3,3)-wm1(i,2,3)*wm1(i,3,2))+
     |    wm1(i,1,2)*(wm1(i,2,3)*wm1(i,3,1)-wm1(i,2,1)*wm1(i,3,3))+
     |    wm1(i,1,3)*(wm1(i,2,1)*wm1(i,3,2)-wm1(i,2,2)*wm1(i,3,1))
      enddo
C     ****
C     ****     NOW INVERSE OF WM1 IN WM3
C     ****
      do i = 1,nlonp4
        wm3(i,1,1) =(wm1(i,2,2)*wm1(i,3,3)-wm1(i,2,3)*wm1(i,3,2))/ws1(i)
        wm3(i,1,2) =(wm1(i,1,3)*wm1(i,3,2)-wm1(i,1,2)*wm1(i,3,3))/ws1(i)
        wm3(i,1,3) =(wm1(i,1,2)*wm1(i,2,3)-wm1(i,1,3)*wm1(i,2,2))/ws1(i)
        wm3(i,2,1) =(wm1(i,2,3)*wm1(i,3,1)-wm1(i,2,1)*wm1(i,3,3))/ws1(i)
        wm3(i,2,2) =(wm1(i,1,1)*wm1(i,3,3)-wm1(i,1,3)*wm1(i,3,1))/ws1(i)
        wm3(i,2,3) =(wm1(i,1,3)*wm1(i,2,1)-wm1(i,1,1)*wm1(i,2,3))/ws1(i)
        wm3(i,3,1) =(wm1(i,2,1)*wm1(i,3,2)-wm1(i,2,2)*wm1(i,3,1))/ws1(i)
        wm3(i,3,2) =(wm1(i,1,2)*wm1(i,3,1)-wm1(i,1,1)*wm1(i,3,2))/ws1(i)
        wm3(i,3,3) =(wm1(i,1,1)*wm1(i,2,2)-wm1(i,1,2)*wm1(i,2,1))/ws1(i)
      enddo
C     ****
C     ****     B = WM3 * WM2
C     ****
! b and fb are bndry_module module data.
      do l = 1,3
        do m = 1,3
          do i = 1,nlonp4
            b(i,l,m) = wm3(i,l,1)*wm2(i,1,m)+wm3(i,l,2)*wm2(i,2,m)+
     |                 wm3(i,l,3)*wm2(i,3,m)
          enddo
        enddo
      enddo
C     ****
C     ****     FB = WM3 * G
C     ****
      do l = 1,3
        do i = 1,nlonp4
          fb(i,l) = wm3(i,l,1)*gg(i,1)+wm3(i,l,2)*gg(i,2)+wm3(i,l,3)*
     |              gg(i,3)
        enddo
      enddo
      end subroutine bndcmp
!-----------------------------------------------------------------------
      subroutine bndef(ee,ff,gg)
C     ****
C     ****     BNDEF DEFINES THE LOWER BOUNDARY CONDITION FOR THIS
C     ****       VERSION OF THE MODEL
C     ****
C     ****     THE LOWER BOUNDARY CONDITION FOR COMP IS:
C     ****
C     ****       E * D(PSI)/DS + F * PSI +G = 0.
C     ****
C     ****         WHERE:
C     ****           PSI = VECTOR(PSI1,PSI2)
C     ****           E AND F ARE 2 X 2 MATRICES
C     ****           G = VECTOR(G1,G2)
C     ****
C     ****           E, F AND G MAY BE FUNCTIONS OF LATITUDE & LONGITUDE
C     ****
C     ****       THIS SUBROUTINE DEFINES E, F AND G for BNDCMP
C     ****
!
! Args:
      real,intent(out) :: ee(nlonp4,3,3),ff(nlonp4,3,3),gg(nlonp4,3)
!
! Local:
      real :: alfa
      integer :: i
C     ****
C     ****     IN TIGCM AND TIEGCM:
C     ****
C     ****       E = |0.  0.  0.|
C     ****           |          |
C     ****           |0.  1.  0.|
C     ****           |          |
C     ****           |0.  0.  0.|
C     ****
C     ****       F = |1.  1.  0.|
C     ****           |          |
C     ****           |0. -1.  0.|
C     ****           |          |
C     ****           |0.  0.  1.|
C     ****
C     ****       G = |-ALFA|
C     ****           |     |
C     ****           |   0.|
C     ****           |     |
C     ****           | -HEB|
C     ****
C     ****       WHERE:
C     ****         ALFA = 0.22 + 0.014 = 0.234
C     ****
      data alfa/0.234/
      do i = 1,nlonp4
        ee(i,1,1) = 0.
        ee(i,1,2) = 0.
        ee(i,1,3) = 0.
        ee(i,2,1) = 0.
        ee(i,2,2) = 1.
        ee(i,2,3) = 0.
        ee(i,3,1) = 0.
        ee(i,3,2) = 0.
        ee(i,3,3) = 0.
        ff(i,1,1) = 1.
        ff(i,1,2) = 1.
        ff(i,1,3) = 0.
        ff(i,2,1) = 0.
        ff(i,2,2) = -1.
        ff(i,2,3) = 0.
        ff(i,3,1) = 0.
        ff(i,3,2) = 0.
        ff(i,3,3) = 1.
        gg(i,1) = -alfa
        gg(i,2) = 0.
        gg(i,3) = -pshelb
      enddo
      end subroutine bndef
!-----------------------------------------------------------------------
      end module lbc
