!
      module chapman_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlonp4,nlat,dz
      use addfld_module,only: addfld
      implicit none
      real,dimension(nlonp4,nlat) :: 
     |  chi,         ! solar zenith angle
     |  sin_chi,     ! sin(chi)
     |  cos_chi,     ! cos(chi)
     |  rt_sinchi,   ! sqrt(sin(chi))
     |  slt          ! local time
!
      contains
!-----------------------------------------------------------------------
      subroutine chapman(z,tn,o2,o1,n2,barm,vo2,vo1,vn2,sco2,sco1,scn2,
     |  lev0,lev1,lon0,lon1,lat)
      use params_module,only: dz
      use init_module,only: secs,sin_sundec,cos_sundec,idn,uthr
      use cons_module,only: pi,dphi,dlamda,re,rmass_o2,rmass_o1,rmass_n2
      use hist_module,only: modeltime ! for diag
      use fields_module,only: tlbc
!
! Calculate line integrals for o2,o n2:
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  z ,   ! geopotential height (cm)
     |  tn,   ! neutral temperature (deg K)
     |  o2,   ! molecular oxygen (mmr)
     |  o1,   ! atomic oxygen (mmr)
     |  n2,   ! molecular nitrogen (mmr)
     |  barm  ! mean molecular weight
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  vo2,  ! o2 vertical integration
     |  vo1,  ! o1 vertical integration
     |  vn2,  ! n2 vertical integration
     |  sco2, ! o2 slant column integration
     |  sco1, ! o1 slant column integration
     |  scn2  ! n2 slant column integration
!
! Local:
      real :: rlat,coslat,sinlat,factor,rtpi,fmin,fmax
      integer :: k,i
      real :: rp(lev0:lev1,lon0:lon1), ! z+re
     |        ti(lev0:lev1,lon0:lon1)  ! tn at interfaces
      real,dimension(lev0:lev1,lon0:lon1) :: ! for diagnostic plotting
     |  chi_plt, slt_plt, idn_plt
!
! Set up sun related geometry according to current model date and time.
! Also set day/night index idn (init_module).
!
      rtpi = sqrt(pi)
      rlat = -.5*pi+(float(lat-1)+.5)*dphi
      coslat=cos(rlat)
      sinlat=sin(rlat)
      idn(lon0:lon1) = 1
      do i=lon0,lon1
        slt(i,lat)=amod(secs/3600.+(float(i-3)*dlamda+pi)*12./pi,24.)
        chi(i,lat)=acos(sin_sundec*sinlat+cos_sundec*coslat*
     |    cos(pi*(slt(i,lat)-12.)/12.))
        sin_chi(i,lat) = sin(chi(i,lat))
        cos_chi(i,lat) = cos(chi(i,lat))
        rt_sinchi(i,lat) = sqrt(sin_chi(i,lat))
        if (chi(i,lat) > 1.8326) idn(i) = 0

        chi_plt(:,i) = chi(i,lat)
        slt_plt(:,i) = slt(i,lat)
        idn_plt(:,i) = idn(i)
      enddo                

!     write(6,"('chapman: lat=',i3,' modeltime=',4i4,
!    |  ' uthr=',f5.2,/,'slt=',(6e12.4))") lat,modeltime,uthr
!     write(6,"('chi=',(12f5.2))") chi(:,lat)
!     write(6,"('slt=',(12f5.2))") slt(:,lat)
!     write(6,"('idn=',(12i2))") idn(:)

!     call addfld('CHI',' ',' ',chi_plt,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SLT',' ',' ',slt_plt,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('IDN',' ',' ',idn_plt,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! rp = z+re
      do k=lev0,lev1
        rp(k,:) = z(k,lon0:lon1)+re
      enddo
!     call addfld('RP',' ',' ',rp,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ti is tn at interfaces:
!     ti(1,:) = tn(lev1,lon0:lon1) ! tn bottom boundary is stored in top slot
      ti(1,:) = tlbc(lon0:lon1,lat) ! Lower boundary is in tlbc
      do k=lev0+1,lev1-1
        ti(k,:) = .5*(tn(k-1,lon0:lon1)+tn(k,lon0:lon1))
      enddo
      ti(lev1,:) = tn(lev1-1,lon0:lon1) ! nlevp1 <- nlev
!     call addfld('TNI',' ',' ',ti,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
!     call addfld('CHAP_O2',' ',' ',o2,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CHAP_O1',' ',' ',o1,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CHAP_N2',' ',' ',n2,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Calculate line integrals (vo[x] and sco[x] are output):
!
      call line_integ(o2,rmass_o2,vo2,sco2,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,1)
      call line_integ(o1,rmass_o1,vo1,sco1,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,2)
      call line_integ(n2,rmass_n2,vn2,scn2,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,3)
!
! Save column number densities:
!     call addfld('VO2',' ',' ',vo2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('VO1',' ',' ',vo1(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('VN2',' ',' ',vn2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Save slant line integrals:
!     call addfld('SCO2',' ',' ',sco2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCO1',' ',' ',sco1(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCN2',' ',' ',scn2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      end subroutine chapman
!-----------------------------------------------------------------------
      subroutine line_integ(f,fmass,v,s,ti,barm,cos_chi,sin_chi,
     |  rt_sinchi,rp,idn,lon0,lon1,lev0,lev1,lat,if)
      use params_module,only: dz
      use cons_module,only: expz,grav,p0,avo,gask,pi,re
      implicit none
!
! Args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat,if
      real,intent(in) :: fmass
      real,dimension(lev0:lev1,lon0:lon1),intent(in) :: 
     |  ti,        ! tn at interfaces
     |  rp         ! z+re         (S9)
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  f,         ! input species density (mmr)
     |  barm       ! mean molecular mass
      real,dimension(lon0:lon1),intent(in) ::
     |  cos_chi,   ! cosine solar zenith angle
     |  sin_chi,   ! sine solar zenith angle
     |  rt_sinchi  ! sqrt(sin_chi)
      integer,intent(in) :: 
     |  idn(lon0:lon1) ! day-night index
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  v,         ! output column density
     |  s          ! output slant-line integral
!
! Local:
      real,parameter :: big=1.e80
      integer :: k,i
      real :: factor,rtpi,exparg
      real :: 
     |  rtrp(lev0:lev1,lon0:lon1),      ! sqrt(rp/2hp) (S6)
     |  yp  (lev0:lev1,lon0:lon1),      ! yp and ip    (S5)
     |  r2ig(lev0:lev1,lon0:lon1)       ! 2.*ig        (S4)
      real,external :: expo
!
! Top:
      factor = avo*p0*expz(lev1-1)*exp(-.5*dz)/(fmass**2*grav)
      do i=lon0,lon1
        v(:,i) = f(:,i)
      enddo
      where(v < 0.) v = 0.
      do i=lon0,lon1
        v(lev1,i) = factor*.5*(v(lev1-1,i)+v(lev1,i))*
     |                    barm(lev1,i)
      enddo
!
! Integrate from top down:
      do i=lon0,lon1
        do k=lev1-1,lev0,-1
          factor = avo*p0*expz(k)/(fmass*grav)*dz
          v(k,i) = v(k+1,i)+factor*v(k,i)
        enddo
      enddo
!
! Set up for slant-line integral:
      factor = fmass*grav/(2.*gask)
      rtpi = sqrt(pi)
      do i=lon0,lon1
        do k=lev0,lev1
          rtrp(k,i) = sqrt(rp(k,i)*factor/ti(k,i)) ! S6
          yp(k,i) = rtrp(k,i)*abs(cos_chi(i))      ! S5
          if (yp(k,i) >= 8.) then
            yp(k,i) = v(k,i)*rtpi*rtrp(k,i)*
     |        (0.56498823/(0.06651874+yp(k,i)))
          else
            yp(k,i) = v(k,i)*rtpi*rtrp(k,i)*
     |        ((1.0606963+0.5564383*yp(k,i))/((yp(k,i)+1.7245609)*
     |          yp(k,i)+1.0619896))
          endif
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
      factor = grav*fmass/gask
      do i=lon0,lon1
        do k=lev0,lev1
!
! Avoid exceeding max arg to exp on non-unicos platform:
          exparg = rp(k,i)*(1.-sin_chi(i))*factor/ti(k,i)

          if (idn(i)==1.and.exparg < 650.) then ! daytime
            r2ig(k,i) = 2.*v(k,i)*exp(exparg)*
     |        rtpi*rt_sinchi(i)*rtrp(k,i)
          else
            r2ig(k,i) = big
          endif

!         r2ig(k,i) = 2.*v(k,i)*expo(exparg,0)*
!    |      rtpi*rt_sinchi(i)*rtrp(k,i)
        enddo
      enddo
!
! Slant line integrals (0 if obscured by earth):
      do i=lon0,lon1
        do k=lev0,lev1
          if (cos_chi(i) >= 0.) then
            s(k,i) = big
          else
            s(k,i) = r2ig(k,i)-yp(k,i)
          endif
          if (rp(k,i)*sin_chi(i)-re < 0.) s(k,i) = big
          if (cos_chi(i) >= 0.) s(k,i) = yp(k,i)
        enddo
      enddo
      end subroutine line_integ
      end module chapman_module
