#include <defs.h>
!
      module getapex_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
!
      use params_module,only:
     |  nlon,   ! number of geographic longitudes (at 5 deg, nlon=72)
     |  nlonp1, ! nlon+1 
     |  nlonp2, ! nlon+2
     |  nlat,   ! number of geographic latitudes (at 5 deg, nlat==36)
     |  nlatp1, ! nlat+1
     |  nmlon,  ! number of geomagnetic longitudes
     |  nmlonp1,! nmlon+1
     |  nmlat   ! number of geomagnetic latitudes
      use apex,only: apex_mka, apex_mall, apex_q2g

      implicit none
      integer,parameter :: 
     |  nalt=2
!
! (formerly in common /dvecddd/)
      real,dimension(nlonp1,nlat,3,2) :: dvec
      real,dimension(nlonp1,nlat) :: dddarr,be3arr
      real,dimension(nlatp1) :: gplat
      real,dimension(nlonp2) :: gplon
      real,dimension(nalt)   :: gpalt
!
! Geographic coordinates of each magnetic grid point, saved for use
! in esmf magnetic/geographic regridding.
      real,dimension(nmlonp1,nmlat) :: gdlatdeg,gdlondeg
!
      contains
!-----------------------------------------------------------------------
      subroutine apxparm(date)
      use cons_module,only: re_dyn,h0,hs,rtd,dtr,ylatg,ylong,dlong,
     |  dlonm,ylonm,ylatm
      use magfield_module,only:
     |  alatm,   ! geomagnetic latitude at each geographic grid point (radians)
     |  alonm,   ! geomagnetic longitude at each geographic grid point (radians)
     |  xb,      ! northward component of magnetic field
     |  yb,      ! eastward component of magnetic field
     |  zb,      ! downward component of magnetic field
     |  bmod,    ! magnitude of magnetic field (gauss?)
     |  rjac,    ! (nlonp1,0:nlatp1,2,2) scaled derivatives of mag coords wrt geog
     |  rjacd,   ! determinant of rjac
     |  p,
     |  cslatg,  ! (0:nlatp1) cos(theta)
     |  snlatg,  ! (0:nlatp1) sin(theta)
     |  im,      ! (nlonp1,0:nlatp1) geomag lon containing each geog point
     |  jm,      ! (nlonp1,0:nlatp1) geomag lon containing each geog point
     |  ig,      ! (nmlonp1,nmlat) geog lon grid containing each geomag point
     |  jg,      ! (nmlonp1,nmlat) geog lat grid containing each geomag point
     |  djm,     ! (nlonp1,0:nlatp1) fraction in lat for grid interp
     |  dim,     ! (nlonp1,0:nlatp1) fraction in lon for grid interp
     |  cslatm,  ! (nlonp1,nlat) cos(thetas)
     |  snlatm,  ! (nlonp1,nlat) sin(thetas)
     |  cslonm,  ! (nlonp1,nlat) cos(lamdas)
     |  snlonm,  ! (nlonp1,nlat) sin(lamdas)
     |  cslong,  ! (nlonp1) cos(lamda)
     |  snlong,  ! (nlonp1) sin(lamda)
     |  av,      ! (nlonp1,0:nlatp1,3,2) the two magnetic vectors a1,a2
     |  wt       ! (4,nmlonp1,nmlat) weights for geo2mag interpolation
!
! Original file for this subroutine 
! ~richmond/prog/tgcmtst/modsrc/apxparm.mk was copied on 2/25/00.
!
! 5/02 B.Foster: adapted for tiegcm1.
!
! Args:
      real,intent(in) :: date
!
! Local:
!
! msgun = Fortran unit number for diagnostics
! ngrf  = Number of epochs in the current DGRF/IGRF; see COFRM in
!         file magfld.f
!
      integer,parameter :: 
     |  msgun=6, 
     |  ngrf=8
!
! Local:
      integer :: j,i,ii,ist,jjm,jjg
      integer :: iplot=0
      real :: rekm,h0km,alt,hr,ror03,glat,glon,dellat,bmag,alon,xlatm,
     |  vmp,w,d,be3,sim,xlatqd,f,xlonmi,qdlon,qdlat,gdlat,xlongi,frki,
     |  frkj,dellon,si,gdlon,date1(1),re
      character(len=80) :: title
!
! Non-scalar arguments returned by APXMALL:
      real :: 
     |  b(3),bhat(3),
     |  d1(3),d2(3),d3(3), 
     |  e1(3),e2(3),e3(3), 
     |  f1(2),f2(2)
!
! Specify grid values
! Center min, max altitudes about 130 km
!
      gpalt(1) =  90.
      gpalt(2) = 170.
      dellat = 180./float(nlatp1-1)
      do j=1,nlatp1
     	gplat(j) = (j-1)*dellat - 90.
      enddo
      dellon = 360./float(nlonp2-2)
      do i=1,nlonp2
     	gplon(i) = (float(i)-1.5)*dellon - 180.
      enddo

!     write(6,"('apxparm: mlon=',i3,' dellon=',f8.3,' gplon=',/,
!    |  (6e12.4))") mlon,dellon,gplon
!     call shutdown('apxparm')

!
      call apex_mka(date,gplat,gplon,gpalt,nlatp1,nlonp2,nalt,ist)
      if (ist /= 0) call shutdown('apex_mka')
!
! Compute dvec, ddd
! re_dyn, h0, hs, rtd, dtr, ylatg, ylong are in constants module (cons.F)
! cslatg, snlatg, alatm, alonm are in magfield module (magfield.F)
!
      re = re_dyn
      rekm = re*1.e-5
      h0km = h0*1.e-5
      alt  = hs*1.e-5
      hr = alt
      ror03= ((rekm + alt)/(rekm + h0km))**3
      do j = 0,nlatp1
	glat = ylatg(j)*rtd
	cslatg(j) = cos(ylatg(j))
	snlatg(j) = sin(ylatg(j))
	do i = 1,nlonp1
	  glon = ylong(i)*rtd
          call apex_mall(
     |      glat,glon,alt,hr,                         !Inputs
     |      b,bhat,bmag,si,                           !Mag Fld
     |      alon,                                     !Apx Lon
     |      xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3,  !Mod Apx
     |      xlatqd,f,f1,f2 , ist)                     !Qsi-Dpl
          if (ist .ne. 0) call shutdown('apex_mall')
!
! Define several fields in magfield module (magfield.F) from apxmall
! output:
!
! alatm,   ! geomagnetic latitude at each geographic grid point (radians)
! alonm,   ! geomagnetic longitude at each geographic grid point (radians)
! im,      ! (nlonp1,0:nlatp1) geomag lon containing each geog point
! jm,      ! (nlonp1,0:nlatp1) geomag lon containing each geog point
!
          alatm(i,j)    =  xlatm*dtr
          alonm(i,j)    =  alon *dtr
          xb   (i,j)    =  b(2)*1.e-5 ! nT -> gauss
          yb   (i,j)    =  b(1)*1.e-5 ! nT -> gauss
          zb   (i,j)    = -b(3)*1.e-5 ! nT -> gauss
          bmod (i,j)    =  bmag*1.e-5 ! nT -> gauss
!
          p    (i,j)    = 1.e+5*abs(sim/si)/(bmag*d)
	  av   (i,j,1,1) = d1(1)*P(I,J)*si
	  av   (i,j,2,1) = d1(2)*P(I,J)*si
	  av   (i,j,3,1) = d1(3)*P(I,J)*si
	  av   (i,j,1,2) =-d2(1)*P(I,J)*si/sim
	  av   (i,j,2,2) =-d2(2)*P(I,J)*si/sim
	  av   (i,j,3,2) =-d2(3)*P(I,J)*si/sim
	  rjac (i,j,1,1) =  f2(2)
	  rjac (i,j,1,2) = -f2(1)
	  rjac (i,j,2,1) = -f1(2)
	  rjac (i,j,2,2) =  f1(1)
	  rjacd(i,j)     = f 
!
! Set up parameters for magnetic to geographic interpolation 
	  xlonmi = (alonm(i,j) - ylonm(1))/dlonm
	  if (xlonmi < 0.) xlonmi = xlonmi + float(nmlon)
	  im(i,j) = xlonmi
	  dim(i,j) = xlonmi - float(im(i,j))
	  im(i,j) = im(i,j) + 1
	  if (im(i,j) >= nmlonp1) im(i,j) = im(i,j) - nmlon
	  alatm(i,j) = min(alatm(i,j),ylatm(nmlat))
          do jjm=2,nmlat
            if (alatm(i,j) > ylatm(jjm)) cycle
            jm(i,j) = jjm - 1
	    djm(i,j) = (alatm(i,j) - ylatm(jm(i,j)))/
     |                 (ylatm(jjm) - ylatm(jm(i,j)))
            exit
          enddo
	  if (j /= 0 .and. j /= nlatp1) then
	    dvec(i,j,1,1) = d1(1)
	    dvec(i,j,2,1) = d1(2)
	    dvec(i,j,3,1) = d1(3)
	    dvec(i,j,1,2) = d2(1)
	    dvec(i,j,2,2) = d2(2)
	    dvec(i,j,3,2) = d2(3)
	    dddarr(i,j)   = d
!
! Scale be3 from 130 km to a reference height of 90 km.
	    be3arr(I,J)   = be3*ror03
!
            cslatm(i,j) = cos(alatm(i,j))
            snlatm(i,j) = sin(alatm(i,j))
            cslonm(i,j) = cos(alonm(i,j))
            snlonm(i,j) = sin(alonm(i,j))
	  endif
        enddo ! i=1,nlonp1

!       write(6,"('apexmag: j=',i3,' bnorth (nT)=',/,(6e12.4))") 
!    |    j,yb(:,j)*1.e5
!       write(6,"('apexmag: j=',i3,' beast  (nT)=',/,(6e12.4))") 
!    |    j,xb(:,j)*1.e5
!       write(6,"('apexmag: j=',i3,' bdown  (nT)=',/,(6e12.4))") 
!    |    j,zb(:,j)*1.e5

      enddo   ! j=0,nlatp1
!
      do i = 1,nlonp1
	cslong(i) = cos(ylong(i))
	snlong(i) = sin(ylong(i))
      enddo
!
! Set up parameters for geographic to magnetic interpolation 
      do i=1,nmlonp1
	qdlon = ylonm(i)*rtd
	do J=1,nmlat
	  qdlat = ylatm(j)*rtd
            call apex_q2g(qdlat,qdlon,alt, gdlat,gdlon, ist)
            if (ist /= 0) then
              write(6,"('apxparm: error from apex_q2g: ist=',i3,
     |          ' i=',i3,' j=',i3)") ist,i,j
              call shutdown('apex_q2g error')
            endif
	  gdlat = gdlat*dtr
	  gdlon = gdlon*dtr
	  xlongi = (gdlon - ylong(1))/dlong
	  if (xlongi < 0.) xlongi = xlongi + float(nlon)
	  ig(i,j) = xlongi
	  frki = xlongi - float(ig(i,j))
	  ig(i,j) = ig(i,j) + 1
	  if (ig(i,j) >= nlonp1) ig(i,j) = ig(i,j) - nlon
	  gdlat = min(gdlat,ylatg(nlatp1))
	  do jjg=1,nlatp1
	    if (gdlat > ylatg(jjg)) cycle
	    jg(i,j) = jjg - 1
	    frkj = (gdlat - ylatg(jg(i,j)))/
     1             (ylatg(jjg) - ylatg(jg(i,j)))
!
! 99/2/25b Add one to JG to account for the fact that AG in geo2mag has
!  a second (J) index starting at 1, while the second index of the
!  array in the calling arguments begins at 0.
!
	    jg(i,j) = jg(i,j) + 1
            exit
          enddo
	  wt(1,i,j) = (1. - frki)*(1. - frkj)
	  wt(2,i,j) =        frki*(1. - frkj)
	  wt(3,i,j) =        frki*frkj
	  wt(4,i,j) = (1. - frki)*frkj
          gdlatdeg(i,j) = gdlat*rtd
          gdlondeg(i,j) = gdlon*rtd
        enddo
      enddo
!
! Contour some results:
! subroutine mkcon(f,i0,i1,cint,log,xlab,ylab,title)
!
!     if (iplot > 0) then
!
! av(nlonp1,0:nlatp1,3,2): the two magnetic vectors a1,a2
! AV(:,:,1:3,1:2)
!     do i=1,2
!       do ii=1,3
!         write(title,"('AV(:,:,',i1,',',i1,')')") ii,i
!         call mkcon(av(:,:,ii,i),imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |      'Geog Latitude',trim(title))
!       enddo
!     enddo
!
! RJAC(:,:,1:2,1:2) 
!     do i=1,2 
!       do ii=1,2
!         write(title,"('RJAC(:,:,',i1,',',i1,')')") ii,i
!         call mkcon(rjac(:,:,ii,i),imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |      'Geog Latitude',trim(title))
!       enddo
!     enddo
!     write(title,"('RJACD')")
!     call mkcon(rjacd,imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!
! alatm,alonm:
!     write(title,"('ALATM')")
!     call mkcon(alatm,imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!     write(title,"('ALONM')")
!     call mkcon(alonm,imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!
! XB,YB,ZZB,bmod:
!     write(title,"('XB')")
!     call mkcon(xb,imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!     write(title,"('YB')")
!     call mkcon(yb,imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!     write(title,"('ZZB')")
!     call mkcon(zzb,imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!     write(title,"('BMOD')")
!     call mkcon(bmod,imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!
! im,jm,dim,djm:
!     write(title,"('IM')")
!     call mkcon(float(im),imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!     write(title,"('JM')")
!     call mkcon(float(jm),imaxgp,jmaxgp,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!!    write(title,"('DIM')")
!!    call mkcon(dim,imaxgp,jmaxgp,0.,0,'Geog Longitude',
!!   |  'Geog Latitude',trim(title))
!!    write(title,"('DJM')")
!!    call mkcon(djm,imaxgp,jmaxgp,0.,0,'Geog Longitude',
!!   |  'Geog Latitude',trim(title))
!
! dvec:
!     do i=1,2
!       do ii=1,3
!         write(title,"('DVEC(:,:,',i1,',',i1,')')") ii,i
!         call mkcon(dvec(:,:,ii,i),imaxgp,jmaxg,0.,0,'Geog Longitude',
!    |      'Geog Latitude',trim(title))
!       enddo
!     enddo
!
! dddarr, be3arr:
!     write(title,"('DDDARR')")
!     call mkcon(dddarr,imaxgp,jmaxg,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!     write(title,"('BE3ARR')")
!     call mkcon(be3arr,imaxgp,jmaxg,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!
! cslatm,snlatm,cslonm,snlonm:
!     write(title,"('CSLATM')")
!     call mkcon(cslatm,imaxgp,jmaxg,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!     write(title,"('SNLATM')")
!     call mkcon(snlatm,imaxgp,jmaxg,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!     write(title,"('CSLONM')")
!     call mkcon(cslonm,imaxgp,jmaxg,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!     write(title,"('SNLONM')")
!     call mkcon(snlonm,imaxgp,jmaxg,0.,0,'Geog Longitude',
!    |  'Geog Latitude',trim(title))
!      
! IG,JG:
!     write(title,"('IG')")
!     call mkcon(float(ig),imaxmp,jmaxm,0.,0,'Geomag Longitude',
!    |  'Geomag Latitude',trim(title))
!     write(title,"('JG')")
!     call mkcon(float(jg),imaxmp,jmaxm,0.,0,'Geomag Longitude',
!    |  'Geomag Latitude',trim(title))
!
!     endif ! iplot
!
      end subroutine apxparm
      end module getapex_module
