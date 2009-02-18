!
      subroutine ionvel(z,ui,vi,wi,lev0,lev1,lon0,lon1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate ExB ion velocities from electric field.
! (electric field ex,ey,ez was calculated in efield module, efield.F)
! (this was old sub vdrift2)
!
      use params_module,only: nlon,nlonp2,nlonp4
      use cons_module,only: re
      use efield_module,only: ex,ey,ez ! (nlevp1,lon0-2:lon1+2,lat0:lat1)

      use magfield_module,only: rjac,xb,yb,zb,bmod
      use addfld_module,only: addfld
#ifdef MPI
      use mpi_module,only: mp_periodic_f3d
#else
      use params_module,only: nlon,nlonp2
#endif
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) :: z
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) :: ui,vi,wi
!
! Local:
      integer :: k,i,nlevs,lonbeg,lonend,lat,ier
      real,dimension(lev0:lev1,lon0:lon1) :: eex,eey,eez
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
#ifdef VT
!     code = 121 ; state = 'ionvel' ; activity='ModelCode'
      call vtbegin(121,ier)
#endif
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2
!
! Latitude scan:
      do lat=lat0,lat1
!       call addfld('EX_IONV',' ',' ',ex(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('EY_IONV',' ',' ',ey(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('EZ_IONV',' ',' ',ez(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
! 
! For latitude lat, rotate ex and ey to geographic orientation using 
!   Jacobian.  Divide by distance from center of earth.
! eex = rotated ex
! eey = rotated ey
!
        do i=lonbeg,lonend
          do k=lev0,lev1
            eex(k,i) = (rjac(i-2,lat,1,1)*ex(k,i,lat)+
     |                  rjac(i-2,lat,2,1)*ey(k,i,lat))/(re+z(k,i,lat))

!           if (mod(nlat,lat)==6) 
!    |      write(6,"(/,'ionvel: 'lat=',i2,' k=',i2,' i=',i2,' rjac=',
!    |        2e12.4,' ex=',e12.4,' ey=',e12.4,' z=',e12.4,' eex=',
!    |        e12.4)") lat,k,i-2,rjac(i-2,lat,1,1),rjac(i-2,lat,2,1),
!    |        ex(k,i,lat),ey(k,i,lat),z(k,i,lat),eex(k,i)

            eey(k,i) = (rjac(i-2,lat,1,2)*ex(k,i,lat)+
     |                  rjac(i-2,lat,2,2)*ey(k,i,lat))/(re+z(k,i,lat))

!           if (mod(nlat,lat)==6) 
!    |      write(6,"(/,'ionvel: 'lat=',i2,' k=',i2,' i=',i2,' rjac=',
!    |        2e12.4,' ex=',e12.4,' ey=',e12.4,' z=',e12.4,' eey=',
!    |        e12.4)") lat,k,i-2,rjac(i-2,lat,1,1),rjac(i-2,lat,2,1),
!    |        ex(k,i,lat),ey(k,i,lat),z(k,i,lat),eey(k,i)

          enddo ! k=lev0,lev1
        enddo ! i=longeg,lonend

!       call addfld('Z_IONV',' ',' ', z(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

        do i=lonbeg,lonend
          do k=lev0+1,lev1-1
            eez(k,i) = ez(k,i,lat)/(z(k+1,i,lat)-z(k-1,i,lat))

!           write(6,"('ionvel: lat=',i2,' k=',i2,' i=',i2,
!    |        ' ez(i,lat,k)=',e14.7,' z(k+1)=',e14.7,' z(k-1)=',
!    |        e14.7)") lat,k,i-2,ez(i,lat,k),z(k+1,i,lat),z(k-1,i,lat)

          enddo ! k=lev0+1,lev1-1
        enddo ! i=lonbeg,lonend
!
! Extrapolate for lower and upper boundaries:
        do i=lonbeg,lonend
          eez(1,i) = 2.*eez(2,i)-eez(3,i)
          eez(lev1,i) = 2.*eez(lev1-1,i)-eez(lev1-2,i)
        enddo

!        call addfld('EEX','E_x','V/cm',eex,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!        call addfld('EEY','E_y','V/cm',eey,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!        call addfld('EEZ','E_z','V/cm',eez,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

!
! ion velocities = (e x b/b**2) (x 1.e6 for m/sec)
! ui = zonal, vi = meridional, wi = vertical
        do k=lev0,lev1
          do i=lonbeg,lonend
            ui(k,i,lat) = -(eey(k,i)*zb(i-2,lat)+eez(k,i)*xb(i-2,lat))*
     |        1.e6/bmod(i-2,lat)**2

! diffs in eez
!           write(6,"('ionvel: lat=',i2,' k=',i2,' i=',i2,
!    |        ' eey=',e14.7,' eez=',e14.7)") lat,k,i-2,eey(k,i),eez(k,i)
!           write(6,"('  xb=',e14.7,' zb=',e14.7)")
!    |        xb(i-2,lat),zb(i-2,lat)
!           write(6,"('  bmod=',e14.7,' ui=',e14.7)") 
!    |        bmod(i-2,lat),ui(k,i,lat)

            vi(k,i,lat) =  (eez(k,i)*yb(i-2,lat)+eex(k,i)*zb(i-2,lat))*
     |        1.e6/bmod(i-2,lat)**2
            wi(k,i,lat) =  (eex(k,i)*xb(i-2,lat)-eey(k,i)*yb(i-2,lat))*
     |        1.e6/bmod(i-2,lat)**2
          enddo ! i=lon0,lon1
        enddo ! k=lev0,lev1
!
! Convert ion velocities from meters to cm for rest of the model:
! (this was done at beginning of heelis.F in earlier versions).
        do i=lon0,lon1
          ui(:,i,lat) = ui(:,i,lat)*100.
          vi(:,i,lat) = vi(:,i,lat)*100.
          wi(:,i,lat) = wi(:,i,lat)*100.
        enddo 
!
! End latitude scan:
      enddo ! lat=lat0,lat1
!
! Periodic points for ion velocities:
#ifdef MPI
      call mp_periodic_f3d(ui(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
      call mp_periodic_f3d(vi(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
      call mp_periodic_f3d(wi(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
#else
      do lat=lat0,lat1
        do i=1,2
          ui(:,i,lat) = ui(:,nlon+i,lat)
          ui(:,nlonp2+i,lat) = ui(:,i+2,lat)
          vi(:,i,lat) = vi(:,nlon+i,lat)
          vi(:,nlonp2+i,lat) = vi(:,i+2,lat)
          wi(:,i,lat) = wi(:,nlon+i,lat)
          wi(:,nlonp2+i,lat) = wi(:,i+2,lat)
        enddo
      enddo
#endif
!
      do lat=lat0,lat1
       call addfld('UI_VEL','UI','cm/s',ui(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
       call addfld('VI_VEL','VI','cm/s',vi(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
       call addfld('WI_VEL','WI','cm/s',wi(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
#ifdef VT
!     code = 121 ; state = 'ionvel' ; activity='ModelCode'
      call vtend(121,ier)
#endif
      end subroutine ionvel
