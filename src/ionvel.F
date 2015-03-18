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
      use pdynamo_module,only: ex,ey,ez ! (nlevp1,lon0-2:lon1+2,lat0:lat1)

      use magfield_module,only: rjac,xb,yb,zb,bmod
      use addfld_module,only: addfld
      use diags_module,only: mkdiag_UI,mkdiag_VI,mkdiag_WI
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
      real :: ftmp(lev0:lev1,lon0:lon1,lat0:lat1,3)
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
!       call addfld('Z_IONV',' ',' ',z(lev0:lev1,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
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
            eey(k,i) = (rjac(i-2,lat,1,2)*ex(k,i,lat)+
     |                  rjac(i-2,lat,2,2)*ey(k,i,lat))/(re+z(k,i,lat))
          enddo ! k=lev0,lev1
        enddo ! i=longeg,lonend

        do i=lonbeg,lonend
          do k=lev0+1,lev1-1
            eez(k,i) = ez(k,i,lat)/(z(k+1,i,lat)-z(k-1,i,lat))
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
            vi(k,i,lat) =  (eez(k,i)*yb(i-2,lat)+eex(k,i)*zb(i-2,lat))*
     |        1.e6/bmod(i-2,lat)**2
            wi(k,i,lat) =  (eex(k,i)*xb(i-2,lat)-eey(k,i)*yb(i-2,lat))*
     |        1.e6/bmod(i-2,lat)**2

#if defined(INTERCOMM) || defined(CISMAH)
            ! CMIT modification
            if (vi(k,i,lat) > 2000.0) vi(k,i,lat) = 2000.0
            if (ui(k,i,lat) > 2000.0) ui(k,i,lat) = 2000.0
            if (vi(k,i,lat) < -2000.0) vi(k,i,lat) = -2000.0
            if (ui(k,i,lat) < -2000.0) ui(k,i,lat) = -2000.0
#endif

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
      ftmp(:,:,:,1) = ui(:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,2) = vi(:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,3) = wi(:,lon0:lon1,lat0:lat1)

      call mp_periodic_f3d(ftmp,lev0,lev1,lon0,lon1,lat0,lat1,3)

      ui(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,1)
      vi(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,2)
      wi(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,3)
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

! real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
!     write(6,"('end ionvel: ui=',2e12.4,' vi=',2e12.4,' wi=',2e12.4)")
!    |  minval(ui(:,lon0:lon1,lat0:lat1)),
!    |  maxval(ui(:,lon0:lon1,lat0:lat1)),
!    |  minval(vi(:,lon0:lon1,lat0:lat1)),
!    |  maxval(vi(:,lon0:lon1,lat0:lat1)),
!    |  minval(wi(:,lon0:lon1,lat0:lat1)),
!    |  maxval(wi(:,lon0:lon1,lat0:lat1))
!
! Save ion drifts if requested:
      call mkdiag_UI('UI_ExB',ui(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
      call mkdiag_VI('VI_ExB',vi(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
      call mkdiag_WI('WI_ExB',wi(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
#ifdef VT
!     code = 121 ; state = 'ionvel' ; activity='ModelCode'
      call vtend(121,ier)
#endif
      end subroutine ionvel
