!
      module magfield_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate magnetic field data used for transformation 
!   between geographic and geomagnetic coordinate systems. 
!   Sub magdyn is called once per run from main tgcm.F.
! Magdyn:    Establish magnetic field and grid (bx,by,bz).
! Sub sunloc is called once per timestep from advance to determine
!   sun's longitudes for current ut model time.
!
      use params_module,only: nmlonp1,nmlat,nlat,nlonp4,nlatp2,nlonp2,
     |  nlon,nlonp1,nlevp1,nmlon,nlatp1
      use nchist_module,only: nc_open,nc_close,handle_ncerr
      use addfld_module,only: addfld
      implicit none
!
! Netcdf header file:
#ifdef SUN
#include <netcdf3.inc>
#else
#include <netcdf.inc>
#endif
!
      real,dimension(nlonp1,0:nlatp1) ::
     |  alatm,   ! geomagnetic latitude at each geographic grid point (radians)
     |  alonm,   ! geomagnetic longitude at each geographic grid point (radians)
     |  xb,      ! northward component of magnetic field
     |  yb,      ! eastward component of magnetic field
     |  zb,      ! downward component of magnetic field
! apex.F: bmod (i,j)    =  bmag*1.e-5
     |  bmod,    ! magnitude of magnetic field (gauss?)
     |  dmlat,   ! dipole latitude corresponding to apex of field line
     |  rmag11,  ! (a1.a1)/p*sin(i)*cos(thetas)
     |  rmagc,   ! (a1.a2)/p*sin(i)
     |  rmag2,   ! 1./bmod
     |  rmag22,  ! (a2.a2)/p*sin(i)/cos(thetas)
     |  rjacd,   ! determinant of rjac
     |  p
!
! rjac: scaled derivatives of geomagnetic coords wrt geographic coordinates.
! rjac(1,1) = cos(thetas)/cos(theta)*d(lamdas)/d(lamda)
! rjac(1,2) = cos(thetas)*d(lamdas)/d(theta)
! rjac(2,1) = 1./cos(theta)*d(thetas)/d(lamda)
! rjac(2,2) = d(thetas)/d(theta)
! where (lamda,theta) are geographic coordinates
!       (lamdas,thetas) are geomagnetic coordinates
!
      real :: rjac(nlonp1,0:nlatp1,2,2) 
!
! av  = the two magnetic vectors a1 and a2
! av1 = a1
! av2 = a2/cos(thetas)
!
      real :: av  (nlonp1,0:nlatp1,3,2)
!
! Quantities needed to transform scalar fields between geographic and
! geomagnetic coordinate systems (these were in old version cterp.h)
!
      integer :: 
     |  ig(nmlonp1,nmlat),   ! geog lon grid containing each geomag point
     |  jg(nmlonp1,nmlat),   ! geog lat grid containing each geomag point
     |  im(nlonp1,0:nlatp1), ! geomag lon grid containing each geog point
     |  jm(nlonp1,0:nlatp1)  ! geomag lat grid containing each geog point
!
! wt(4) are interpolation weights to be applied to function values at 4 
!   corners of geographic grid element (ig,jg)
! dim and djm are fractions in i and j directions, used for bilinear 
!   interpolation in geomagnetic grid element.
!
      real ::
     |  wt(4,nmlonp1,nmlat),   
     |  dim(nlonp1,0:nlatp1),
     |  djm(nlonp1,0:nlatp1)
!
! Trigonometric factors needed in the calculation of the derivatives
! of the geomagnetic coordinates wrt the geographic coordinates
! (these were in old versions header file trig.h) 
!
      real,dimension(nlonp1,nlat) ::
     |  cslatm, ! cos(thetas)
     |  snlatm, ! sin(thetas)
     |  cslonm, ! cos(lamdas)
     |  snlonm  ! sin(lamdas)
!
! 10/17/02 btf: cslatg and snlatg dim changed from nlat to 0:nlatp1
!   for apxparm (apex.F).
      real ::
     |  cslatg(0:nlatp1), ! cos(theta)
     |  snlatg(0:nlatp1), ! sin(theta)
     |  cslong(nlonp1),   ! cos(lamda)
     |  snlong(nlonp1)    ! sin(lamda)
!
! Magnetic field at geographic grid. These are defined by magdyn after 
!   reading magdat file (formerly in common /MAGFLD/ and header file trgm.h)
! (note bx,by,bz,bmod2 are dimensioned nlonp4 rather than -1:nlonp2, as
!  in earlier versions, so these fields can be referenced at i rather than
!  i-2, e.g. in sub oplus)
!
      real,dimension(nlonp4,-1:nlatp2) :: bx,by,bz,bmod2
      real,dimension(nlonp4,nlat) ::
     |  rlatm,rlonm,dipmag,decmag,sndec,csdec,sn2dec,sncsdc
!
! sunlons: sun's longitude in dipole coordinates (see sub sunloc)
! (this was dlons in earlier versions)
!
      real :: sunlons(nlat)
!
      contains
!-----------------------------------------------------------------------
#include <defs.h>
      subroutine magdyn
!
! Define magnetic field and related quantities:
!
! Local:
!
! sin10 should be the same as dipmin (see cons.F)
!
#if (NLAT==36 && NLON==72)
      real,parameter :: sin10=0.17   ! 5.0 deg horizontal resolution
#elif (NLAT==72 && NLON==144)
      real,parameter :: sin10=0.24   ! 2.5 deg horizontal resolution
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
      real :: cos10
!
! Field are calculated in local arrays dimensioned -1:nlonp2, then transferred 
!   to module data bx,by,bz,bmod2 at 1:nlonp4 for use in the model.
!
      real,dimension(-1:nlonp2,-1:nlatp2) :: bxtmp,bytmp,bztmp,bmod2tmp
      integer :: i,j
!
      cos10 = sqrt(1.-sin10**2)
      do j = 1,nlat
        do i = 1,nlon
          rlatm(i+2,j) = alatm(i,j)
          rlonm(i+2,j) = alonm(i,j)
          dipmag(i+2,j) = atan(zb(i,j)/sqrt(xb(i,j)**2+yb(i,j)**2))
          decmag(i+2,j) = atan2(yb(i,j),xb(i,j))
          sndec(i+2,j) = sin(decmag(i+2,j))
          csdec(i+2,j) = cos(decmag(i+2,j))
          sn2dec(i+2,j) = sndec(i+2,j)**2
          sncsdc(i+2,j) = sndec(i+2,j)*csdec(i+2,j)
          bxtmp(i,j) = yb(i,j)/bmod(i,j)
          bytmp(i,j) = xb(i,j)/bmod(i,j)
          bztmp(i,j) = -zb(i,j)/bmod(i,j)
          bmod2tmp(i,j) = bmod(i,j)
!
! Set minimum dip to 10 degrees
          if (abs(bztmp(i,j))-sin10 < 0.) then
            bxtmp(i,j)=bxtmp(i,j)*(cos10/sqrt(1.-bztmp(i,j)**2))
            bytmp(i,j)=bytmp(i,j)*(cos10/sqrt(1.-bztmp(i,j)**2))
            bztmp(i,j)=sign(sin10,bztmp(i,j))
          endif
        enddo
      enddo
!
! Values at j = -1, 0, nlatp1, nlatp2:
      do j = 1,2
        do i = 1,nlon
          bxtmp(i,j-2) = -bxtmp(1+mod(i-1+nlon/2,nlon),3-j)
          bytmp(i,j-2) = -bytmp(1+mod(i-1+nlon/2,nlon),3-j)
          bztmp(i,j-2) =  bztmp(1+mod(i-1+nlon/2,nlon),3-j)
          bmod2tmp(i,j-2) = bmod2tmp(1+mod(i-1+nlon/2,nlon),3-j)
          bxtmp(i,nlat+j) = -bxtmp(1+mod(i-1+nlon/2,nlon),nlat+1-j)
          bytmp(i,nlat+j) = -bytmp(1+mod(i-1+nlon/2,nlon),nlat+1-j)
          bztmp(i,nlat+j) =  bztmp(1+mod(i-1+nlon/2,nlon),nlat+1-j)
          bmod2tmp(i,nlat+j) = bmod2tmp(1+mod(i-1+nlon/2,nlon),nlat+1-j)
        enddo
      enddo
!
! Periodic points:
      do i = 1,2
        do j = 1,nlat
          rlatm (i,j) = rlatm (i+nlon,j)
          rlonm (i,j) = rlonm (i+nlon,j)
          dipmag(i,j) = dipmag(i+nlon,j)
          decmag(i,j) = decmag(i+nlon,j)
          sndec (i,j) = sndec (i+nlon,j)
          csdec (i,j) = csdec (i+nlon,j)
          sn2dec(i,j) = sn2dec(i+nlon,j)
          sncsdc(i,j) = sncsdc(i+nlon,j)
!
          rlatm (i+nlonp2,j) = rlatm (i+2,j)
          rlonm (i+nlonp2,j) = rlonm (i+2,j)
          dipmag(i+nlonp2,j) = dipmag(i+2,j)
          decmag(i+nlonp2,j) = decmag(i+2,j)
          sndec (i+nlonp2,j) = sndec (i+2,j)
          csdec (i+nlonp2,j) = csdec (i+2,j)
          sn2dec(i+nlonp2,j) = sn2dec(i+2,j)
          sncsdc(i+nlonp2,j) = sncsdc(i+2,j)
        enddo
      enddo
      do i = 1,2
        do j = -1,nlatp2
          bxtmp(i-2,j) = bxtmp(i-2+nlon,j)
          bytmp(i-2,j) = bytmp(i-2+nlon,j)
          bztmp(i-2,j) = bztmp(i-2+nlon,j)
          bmod2tmp(i-2,j) = bmod2tmp(i-2+nlon,j)

          bxtmp(i+nlon,j) = bxtmp(i,j)
          bytmp(i+nlon,j) = bytmp(i,j)
          bztmp(i+nlon,j) = bztmp(i,j)
          bmod2tmp(i+nlon,j) = bmod2tmp(i,j)
        enddo
      enddo
!
! Save local fields -1:nlonp2 to module data at 1:nlonp4.
! This way, model routines (e.g., oplus) can reference i rather than i-2.
      do i=1,nlonp4
        bx(i,:) = bxtmp(i-2,:)
        by(i,:) = bytmp(i-2,:)
        bz(i,:) = bztmp(i-2,:)
        bmod2(i,:) = bmod2tmp(i-2,:)
      enddo ! i=1,nlonp4

!     do j=-1,nlatp2
!       write(6,"('magdyn: j=',i2)") j
!       write(6,"('bx(:,j)=',/,(6e12.4))") bx(:,j)
!       write(6,"('by(:,j)=',/,(6e12.4))") by(:,j)
!       write(6,"('bz(:,j)=',/,(6e12.4))") bz(:,j)
!     enddo ! j=-1,nlatp2
!
      end subroutine magdyn
!-----------------------------------------------------------------------
      subroutine sunloc_apex(iyr,iday,secs)

      use cons_module,only: dtr  ! degrees to radians (pi/180)
      use apex,only: subsol,cofrm,dypol,solgmlon
! am 10/04
! tiegcm uses the sun's longitude in dipole coordinates
! we changed the approximation of the sun's location in geographic
! coodinates from
!      glats=asin(.398749*sin(2.*PI*(iday-80)/365.))
!      glons=pi*(1.-2.*secs/86400.)
! to use the apex routines based on formulas in Astronomical Almanac
! difference is around 6/7 min
! This is called every timestep from advance.
!
! Args:
      integer,intent(in) :: iyr,  ! year
     |   iday ! day of year
      real,intent(in) :: secs    ! ut in seconds
!
! Local:
      integer :: ihr,imn,j
      real :: sec,date,vp,xmlon, ! apex magnetic longitude
     |  sbsllat,    ! geographic latitude of subsolar point (degrees)
     |  sbsllon,    ! geographic longitude of subsolar point (degrees)
     |  colat,      ! Geocentric colatitude of geomagnetic dipole north pole (deg)
     |  elon        ! East longitude of geomagnetic dipole north pole (deg)

      ihr = int(secs/3600.)
      imn = int((secs - float(ihr)*3600.)/60.)
      sec = secs - float(ihr)*3600. - float(imn)*60.

!  calculate subsol point: given universal time
!          input: iyr,iday,ihr,imn,sec
!          output: sbsllat,sbsllon
!
      call subsol(iyr,iday,ihr,imn,sec ,sbsllat,sbsllon)

      date = float(iyr) + float(iday)/365. + float(ihr)/24./365. +
     |  float(imn)/60./24./365.+ sec/60./60./24./365.
      call cofrm(date)
      call dypol(colat,elon,vp)

! calculate geomagn. diploe longitude
!        input: aloni,sbsllat,sbsllon,colat,elon
!        output: xmlon
      call solgmlon(sbsllat,sbsllon,colat,elon,xmlon)
!
      sunlons(1) = xmlon*dtr
      do j = 2,nlat
        sunlons(j) = sunlons(1)
      enddo
!      write(6,*) 'sunloc_apex ',xmlon
!    
      end subroutine sunloc_apex
!-----------------------------------------------------------------------
      end module magfield_module
