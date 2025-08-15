!
      module chapman_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlonp4,nlat
      implicit none
      real,dimension(nlonp4,nlat) :: 
     |  chi         ! solar zenith angle
!
      contains
!-----------------------------------------------------------------------
      subroutine chapman(mbar_ubc,xnmbar_ubc,z,tn,o2,o1,n2,xnmbar,scht,
     |  vo2,vo1,vn2,sco2,sco1,scn2,lev0,lev1,lon0,lon1,lat)
      use params_module,only: glat,glon0
      use init_module,only: secs,sin_sundec,cos_sundec
      use cons_module,only: pi,dtr,rmass_o2,rmass_o1,rmass_n2
      use fields_module,only: tlbc
      use addfld_module,only: addfld
!
! Calculate line integrals for o2,o,n2:
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
      real,dimension(lon0-2:lon1+2),intent(in) ::
     |  mbar_ubc, ! mean molecular weight
     |  xnmbar_ubc
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  z ,   ! geopotential height (cm)
     |  tn,   ! neutral temperature (deg K)
     |  o2,   ! molecular oxygen (mmr)
     |  o1,   ! atomic oxygen (mmr)
     |  n2,   ! molecular nitrogen (mmr)
     |  xnmbar,scht
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  vo2,  ! o2 vertical integration
     |  vo1,  ! o1 vertical integration
     |  vn2,  ! n2 vertical integration
     |  sco2, ! o2 slant column integration
     |  sco1, ! o1 slant column integration
     |  scn2  ! n2 slant column integration
!
! Local:
      real :: rlat
      integer :: k,i
      real :: ti(lev0:lev1,lon0:lon1)  ! tn at interfaces
      real,dimension(lon0:lon1) ::
     |  slt          ! local time
!
! Set up sun related geometry according to current model date and time.
!
      rlat = glat(lat)*dtr
      do i=lon0,lon1
        slt(i)=amod(secs/3600.+glon0(i)/15.,24.)
        chi(i,lat)=acos(sin_sundec*sin(rlat)+cos_sundec*cos(rlat)*
     |    cos(pi*(slt(i)-12.)/12.))
      enddo                

!     write(6,"('chapman: lat=',i3) lat
!     write(6,"('chi=',(12f5.2))") chi(:,lat)
!     write(6,"('slt=',(12f5.2))") slt(:)
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
! Calculate line integrals (v[x] and sc[x] are output):
!
      call line_integ(o2,rmass_o2,vo2,sco2,ti,xnmbar,scht,z,
     |  chi(lon0:lon1,lat),mbar_ubc(lon0:lon1),xnmbar_ubc(lon0:lon1),
     |  lon0,lon1,lev0,lev1)
      call line_integ(o1,rmass_o1,vo1,sco1,ti,xnmbar,scht,z,
     |  chi(lon0:lon1,lat),mbar_ubc(lon0:lon1),xnmbar_ubc(lon0:lon1),
     |  lon0,lon1,lev0,lev1)
      call line_integ(n2,rmass_n2,vn2,scn2,ti,xnmbar,scht,z,
     |  chi(lon0:lon1,lat),mbar_ubc(lon0:lon1),xnmbar_ubc(lon0:lon1),
     |  lon0,lon1,lev0,lev1)
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
      subroutine line_integ(f,fmass,v,s,ti,xnmbar,scht,z,
     |  chi,mbar_ubc,xnmbar_ubc,lon0,lon1,lev0,lev1)
      use params_module,only: dz
      use cons_module,only: expzmid,expz,grav,p0,avo,gask,pi,re,dtr
!
! Args:
      integer,intent(in) :: lon0,lon1,lev0,lev1
      real,intent(in) :: fmass
      real,dimension(lev0:lev1,lon0:lon1),intent(in) :: 
     |  ti         ! tn at interfaces
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  f,         ! input species density (mmr)
     |  xnmbar,scht,z
      real,dimension(lon0:lon1),intent(in) ::
     |  chi,       ! solar zenith angle
     |  mbar_ubc,  ! mean molecular weight
     |  xnmbar_ubc
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  v,         ! output column density
     |  s          ! output slant-line integral
!
! Local:
      real,parameter :: big=1.e80, logbig=log(big)
      integer :: k,i
      real :: yparg,yppart,exparg
      real,dimension(lon0:lon1) :: f_ubc_cm3
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  rp,        ! z+re         (S9)
     |  hp,rp_hp,f_cm3,
     |  yp  ,      ! yp and ip    (S5)
     |  r2ig       ! 2.*ig        (S4)
!
      do i=lon0,lon1
        do k=lev0,lev1
          rp(k,i) = z(k,i)+re
          hp(k,i) = gask*ti(k,i)/(fmass*grav)
          rp_hp(k,i) = rp(k,i)/hp(k,i)
          f_cm3(k,i) = f(k,i)*xnmbar(k,i)/fmass
        enddo
        f_ubc_cm3(i) = .5*(f(lev1-1,i)+f(lev1,i))*xnmbar_ubc(i)/fmass
      enddo
!
! Top:
      do i=lon0,lon1
        v(lev1,i) = avo*p0*expz(lev1-1)*expzmid/(fmass**2*grav)*
     |    .5*(f(lev1-1,i)+f(lev1,i))*mbar_ubc(i)
!       v(lev1,i) = f_ubc_cm3(i)*hp(lev1,i)
      enddo
!
! Integrate from top down:
      do i=lon0,lon1
        do k=lev1-1,lev0,-1
          v(k,i) = v(k+1,i)+avo*p0*expz(k)/(fmass*grav)*dz*f(k,i)
!         v(k,i) = v(k+1,i)+f_cm3(k,i)*scht(k,i)*dz
        enddo
      enddo
!
! Set up for slant-line integral:
      do i=lon0,lon1
        do k=lev0,lev1
          yparg = sqrt(rp_hp(k,i)/2.)*abs(cos(chi(i)))       ! S5
          if (yparg >= 8.) then
            yppart = 0.56498823/(0.06651874+yparg)
          else
            yppart = (1.0606963+0.5564383*yparg)/
     |        ((yparg+1.7245609)*yparg+1.0619896)
          endif
          yp(k,i) = v(k,i)*sqrt(pi*rp_hp(k,i)/2.)*yppart
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
      do i=lon0,lon1
        do k=lev0,lev1
!
! Avoid exceeding max arg to exp on non-unicos platform:
          exparg = rp_hp(k,i)*(1.-sin(chi(i)))

          if (chi(i)<=105*dtr.and.exparg < logbig) then ! daytime
            r2ig(k,i) = 2.*v(k,i)*exp(exparg)*
     |        sqrt(pi*sin(chi(i))*rp_hp(k,i)/2.)
          else
            r2ig(k,i) = big
          endif
        enddo
      enddo
!
! Slant line integrals (0 if obscured by earth):
      do i=lon0,lon1
        do k=lev0,lev1
          if (cos(chi(i)) < 0.) s(k,i) = r2ig(k,i)-yp(k,i)
          if (rp(k,i)*sin(chi(i))-re < 0.) s(k,i) = big
          if (cos(chi(i)) >= 0.) s(k,i) = yp(k,i)
        enddo
      enddo
      end subroutine line_integ
      end module chapman_module
