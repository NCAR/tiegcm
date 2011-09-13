!
      module efield_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlon,nlev,nlat,nlonp4,nmlon,nmlat,nmlonp1,
     |  nmlonp2,nlevp1,nlonp2,nmlatp1,nmlev
      use mpi_module,only: mytid ! for condx on addfld call only
      use addfld_module,only: addfld
      implicit none
!
! Electric field components ex,ey,ez on geographic grid are allocated at 
!   task subdomains (nlevp1,lon0-2:lon1+2,lat0:lat1) by subroutine alloc_e 
!   (called from allocdata). 
! The electric field (on geographic grid) is used by sub ionvel to calculate 
!   ExB ion drifts ui,vi,wi.
!
      real,allocatable,dimension(:,:,:) :: 
     |  ex, ! zonal component of electric field (geographic)
     |  ey, ! meridional component of electric field (geographic)
     |  ez  ! vertical component of electric field (geographic)
!
! Electric field on geomagnetic grid:
      real ::
     |  emx(nmlonp1,nmlat,nlevp1), ! zonal e-field on magnetic grid
     |  emy(nmlonp1,nmlat,nlevp1), ! meridional e-field on magnetic grid
     |  emz(nmlonp1,nmlat,nlevp1)  ! vertical e-field on magnetic grid
      contains
!-----------------------------------------------------------------------
      subroutine efield(lev0,lev1,lon0,lon1,lat0,lat1)
!
! Calculate 3-d electric field by numerical differentiation of 3-d dynamo 
!   potential on geomagnetic grid. Save 3-d electric field on both magnetic 
!   (emx,emy,emz) and geographic grids (ex,ey,ez).
!
! On input:
!   Phim3d(nmlonp1,nmlat,-2:nlevp1) is 3d electric potential on magnetic 
!   grid (was transformed from source history potential on geographic grid
!   to the magnetic grid by sub dynpotmag in magfield.F).
! On output:
!   Ex,ey,ez: 3-d electric field on geographic grid.
!   Emx,emy,emz: 3-d electric field on magnetic grid.
!
! This routine is called once per time step from sub advance.
! This routine is a rewritten version of the old sub vdrift.
!
      use cons_module,only: pi,dlatm,dlonm,dt0dts,rcos0s,table,rtd,
     |  ylatg,ylong
      use magfield_module,only: alatm,alonm
      use fields_module,only: 
     |  phim3d,    ! 3d electric potential on magnetic grid
     |  emphi3d,   ! 3d eastward electric field magnetic
     |  emlam3d,   ! 3d equatorw. electric field magnetic
     |  emz3d      ! 3d upward (?) electric field magnetic
      use init_module,only: istep
      use params_module,only: 
     |  nmlev,  ! number of geomagnetic pressure levels (nmlev==nlevp1+3)
     |  nmlevp1 ! number of geomagnetic midpoint levels
#ifdef MPI
      use mpi_module,only: mp_periodic_f3d
#endif
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,ii,j,k,lonbeg,lonend,ier
      integer,parameter :: nmlonh = nmlon/2
      real :: 
     |  phi3d(0:nmlonp2,0:nmlatp1,nlevp1), ! electric potential magnetic grid
     |  thetam(nlonp4),
     |  qslot(nlonp4), 
     |  pslot(nlonp4)
      real :: csth0
      integer :: 
     |  islot(nlonp4),
     |  jslot(nlonp4)
      real :: ftmp(lev0:lev1,lon0:lon1,lat0:lat1,3)
!
#ifdef VT
!     code = 117 ; state = 'efield' ; activity='ModelCode'
      call vtbegin(117,ier)
#endif
!     write(6,"('enter efield: lon0,1=',2i3,' lat0,1=',2i3)")
!    |  lon0,lon1,lat0,lat1
!
! Electric field emx,y,z is calculated from phi at the first timestep.
! Thereafter, emphi3d,emlam3d,emz3d are calculated in dynamo on the
!   mag grid, and copied into emx,y,z below.
!
      if(istep == 1) then ! do only for the first timestep at the moment
!
! Copy phim3d to local phi3d. Phim3d is dimensioned in fields.F, and is 
! defined from electric potential dynpot (geographic electric potential
! from source history), by sub dynpotmag in magfield.F.
!
        do k=1,nlevp1
          do j=1,nmlat
            do i=1,nmlonp1 
              phi3d(i,j,k) = phim3d(i,j,k)
            enddo ! i=1,nmlonp1 
          enddo ! j=1,nmlat
        enddo ! k=1,nlevp1
!
! Add points corresponding to j=0 and j=nmlat+1 performing extrapolation
!   across the poles.
!
        do k=1,nlevp1
          do i=1,nmlonp1
            phi3d(i,0,k) = phi3d(1+mod(i-1+nmlonh,nmlon),2,k)
            phi3d(i,nmlatp1,k)=phi3d(1+mod(i-1+nmlonh,nmlon),nmlat-1,k)
          enddo ! i=1,nmlonp1
!
! Periodic points for phi3d: 
          do j=0,nmlatp1
            phi3d(0,j,k) = phi3d(nmlon,j,k)
            phi3d(nmlonp2,j,k) = phi3d(2,j,k)
          enddo ! j=0,nmlatp1
        enddo ! k=1,nlevp1
!
! Emy is the meridional component of the electric field, at all 
!   geomagnetic grid points. The factor dt0dts converts from the 
!   latitudinally distorted geomagnetic grid to standard apex latitudes.
!
        do k=1,nlevp1 
          do j=1,nmlat
            do i=1,nmlonp1
              emy(i,j,k) = -(phi3d(i,j+1,k)-phi3d(i,j-1,k))/
     |          (2.*dlatm)*dt0dts(j)
            enddo ! i=1,nmlonp1
          enddo ! j=1,nmlat
!
! Emx is the zonal electric field for all geomagnetic grid points.
!   The factor rcos0s converts from the latitudinally distorted 
!   geomagnetic grid to standard apex latitudes.
!
          do j=2,nmlat-1
            csth0 = cos(-pi/2.+float(j-1)*dlatm)
            do i=1,nmlonp1
              emx(i,j,k) = -(phi3d(i+1,j,k)-phi3d(i-1,j,k))/
     |          (2.*dlonm*csth0)*rcos0s(j)
            enddo ! i=1,nmlonp1
          enddo ! j=2,nmlat-1
!
! Polar values for emx:
          do i=1,nmlonp1
            emx(i,1,k) = emy(1+mod(i-1+(nmlon/4),nmlon),1,k)
            emx(i,nmlat,k) = emy(1+mod(i-1+((3*nmlon)/4),nmlon),nmlat,k)
          enddo ! i=1,nmlonp1
        enddo ! k=1,nlevp1 
!
! Emz = d(phi)/dz
        do k=2,nlev
          do j=1,nmlat
            do i=1,nmlonp1
              emz(i,j,k) = -(phim3d(i,j,k+1)-phi3d(i,j,k-1))
            enddo ! i=1,nmlonp1
          enddo ! j=1,nmlat
        enddo ! k=2,nlev
!
! istep > 1: copy to emx,y,z:
      else
!
! Copy fields emphi3d, emlam3d, emz3d into emx, emy, emz
        do k=1,nlevp1
          do j=1,nmlat
            do i=1,nmlonp1
              emx(i,j,k) = emphi3d(i,j,k)
              emy(i,j,k) = emlam3d(i,j,k)
              emz(i,j,k) = emz3d(i,j,k)
            enddo ! i=1,nmlonp1
          enddo ! j=1,nmlat
        enddo ! k=1,nlevp1
      endif ! istep
! am 3/10 start with lev=4 since not all mlec are copied over      
      do j=1,nmlat
        call addfld('EMX',' ',' ',emx(:,j,:),
     |	 'mlon',1,nmlonp1,'imlev',4,nmlevp1,j)
        call addfld('EMPHI_EF',' ',' ',emphi3d(:,j,:),
     |	 'mlon',1,nmlonp1,'imlev',1,nmlevp1,j)
        call addfld('EMY',' ',' ',emy(:,j,:),
     |	 'mlon',1,nmlonp1,'mlev',4,nmlevp1,j)
        call addfld('EMZ',' ',' ',emz(:,j,:),
     |	 'mlon',1,nmlonp1,'mlev',4,nmlevp1,j)
      enddo ! j=1,nmlat
!
! Transform emx,emy,emz to geographic space. For each geographic grid 
!   point, determine the corresponding theta0 in the distorted magnetic
!   latitude grid.
! alatm(i,j) is the thetas (standard apex latitude) corresponding to 
!   geographic grid point (i,j).  Find the slot for this value in the 
!   table of theta0 vs thetas, which has 1 degree spacing from 0 to 90
!   degrees.
! Define ex,ey,ez at task subdomain:
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2
!
      do j=lat0,lat1      
        do i=lonbeg,lonend
          ii = i-2
          pslot(i) = abs(alatm(ii,j))*180./pi+1.
          islot(i) = int(pslot(i))
          pslot(i) = pslot(i)-float(islot(i))
!
! Interpolate for theta0 in table using linear interpolation:
          thetam(i) = sign((1.-pslot(i))*table(islot(i),2)+
     |      pslot(i)*table(islot(i)+1,2),alatm(ii,j))
!
! Locate magnetic grid element containing geographic grid point (i,j):
          pslot(i) = (alonm(ii,j)+pi)/dlonm+1.
          islot(i) = int(pslot(i))
          pslot(i) = pslot(i)-float(islot(i))
          qslot(i) = (thetam(i)+pi/2.)/dlatm+1.
          jslot(i) = int(qslot(i))
          qslot(i) = qslot(i)-float(jslot(i))
        enddo ! i=lon0,lon1
!
! Transform emx,emy,emz to geographic frame as ex,ey,ez:
        do i=lonbeg,lonend
          do k=1,nlevp1
            ex(k,i,j) = 
     |        (1.-pslot(i))*(1.-qslot(i))*emx(islot(i)  ,jslot(i)  ,k)+
     |            pslot(i) *(1.-qslot(i))*emx(islot(i)+1,jslot(i)  ,k)+
     |            pslot(i) *    qslot(i) *emx(islot(i)+1,jslot(i)+1,k)+
     |        (1.-pslot(i))*    qslot(i) *emx(islot(i)  ,jslot(i)+1,k)

            ey(k,i,j) = 
     |        (1.-pslot(i))*(1.-qslot(i))*emy(islot(i)  ,jslot(i)  ,k)+
     |            pslot(i) *(1.-qslot(i))*emy(islot(i)+1,jslot(i)  ,k)+
     |            pslot(i) *    qslot(i) *emy(islot(i)+1,jslot(i)+1,k)+
     |        (1.-pslot(i))*    qslot(i) *emy(islot(i)  ,jslot(i)+1,k)
     
          enddo ! k=1,nlevp1
        enddo ! i=lon0,lon1

!       if (j <= 2.or.j >= nlat-1) then
!         k=19 ! zp = +2
!         write(6,"('efield: j=',i3,' i=',i3,' k=',i3,' (lat0,1=',2i3,
!    |      ' lonbeg,end=',2i3,')')") j,i,k,lat0,lat1,lonbeg,lonend
!         write(6,"('  pslot=',/,(8f8.2))")  pslot(lonbeg:lonend)
!         write(6,"('  qslot=',/,(8f8.2))")  qslot(lonbeg:lonend)
!         write(6,"('  islot=',/,(8i8))")  islot(lonbeg:lonend)
!         write(6,"('  jslot=',/,(8i8))")  jslot(lonbeg:lonend)
!         do i=lonbeg,lonend
!           write(6,"('  i=',i3,' lonbeg,end=',2i3,
!    |        ' emx(islot(i),jslot(i),k)=',e12.4)") i,lonbeg,lonend,
!    |          emx(islot(i),jslot(i),k)
!         enddo
!         write(6,"('  ex(k,lonbeg:lonend,j)=',/,(6e12.4))")  
!    |      ex(k,lonbeg:lonend,j)
!       endif

        do i=lonbeg,lonend
          do k=1,nlevp1
            ez(k,i,j) = 
     |        (1.-pslot(i))*(1.-qslot(i))*emz(islot(i)  ,jslot(i)  ,k)+
     |            pslot(i) *(1.-qslot(i))*emz(islot(i)+1,jslot(i)  ,k)+
     |            pslot(i) *    qslot(i) *emz(islot(i)+1,jslot(i)+1,k)+
     |        (1.-pslot(i))*    qslot(i) *emz(islot(i)  ,jslot(i)+1,k)
          enddo ! i=1,nlon
        enddo ! k=2,nlev
      enddo ! j=lat0,lat1
!
! Periodic points:

#ifdef MPI
! real :: ftmp(lev0:lev1,lon0:lon1,lat0:lat1,3)
      ftmp(:,:,:,1) = ex(:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,2) = ey(:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,3) = ez(:,lon0:lon1,lat0:lat1)

      call mp_periodic_f3d(ftmp,lev0,lev1,lon0,lon1,lat0,lat1,3)

      ex(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,1)
      ey(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,2)
      ez(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,3)
#else
      do j=lat0,lat1
        do i=1,2
          ex(:,i,j) = ex(:,nlon+i,j)
          ex(:,nlonp2+i,j) = ex(:,i+2,j)
          ey(:,i,j) = ey(:,nlon+i,j)
          ey(:,nlonp2+i,j) = ey(:,i+2,j)
          ez(:,i,j) = ez(:,nlon+i,j)
          ez(:,nlonp2+i,j) = ez(:,i+2,j)
        enddo
      enddo
#endif
!     do j=lat0,lat1
!       call addfld('EX',' ',' ',ex(lev0:lev1-1,lon0:lon1,j),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!       call addfld('EY',' ',' ',ex(lev0:lev1-1,lon0:lon1,j),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!       call addfld('EZ',' ',' ',ex(lev0:lev1-1,lon0:lon1,j),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!     enddo ! j=lat0,lat1
!
#ifdef VT
!     code = 117 ; state = 'efield' ; activity='ModelCode'
      call vtend(117,ier)
#endif
      end subroutine efield
!-----------------------------------------------------------------------
      subroutine alloc_e(lon0,lon1,lat0,lat1)
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
      allocate(ex(nlevp1,lon0-2:lon1+2,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_e: error allocating',
     |  ' ex: stat=',i3)") istat
      allocate(ey(nlevp1,lon0-2:lon1+2,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_e: error allocating',
     |  ' ey: stat=',i3)") istat
      allocate(ez(nlevp1,lon0-2:lon1+2,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_e: error allocating',
     |  ' ez: stat=',i3)") istat
      end subroutine alloc_e
!-----------------------------------------------------------------------
      end module efield_module
