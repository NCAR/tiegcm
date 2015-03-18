!
      subroutine swdot(un,vc,z,w,lev0,lev1,lon0,lon1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate omega for vertical velocity W (s.a., divrg.F):
!
      use cons_module,only: expzmid,dz
      use addfld_module,only: addfld
      use diags_module,only: mkdiag_WN
      use mpi_module,only: mytid
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1

      real,intent(in),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)::
     |  un,vc,z
      real,intent(out) :: 
     |  w(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
!
! Local:
      integer :: k,i,lat
      real :: w_divrg(lev0:lev1,lon0:lon1)
!
! Latitude scan:
      do lat=lat0,lat1
!       call addfld('UN_SWDOT',' ',' ',un(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VC_SWDOT',' ',' ',vc(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! W=DIVRG(V)
!
        call divrg(un(:,:,lat),vc,w_divrg,lon0,lon1,lev0,lev1,lat0,lat1,
     |    lat)
!
! nlevp1 <- 1:
        do i=lon0,lon1
          w(lev1,i,lat) = w_divrg(lev1-1,i)
        enddo
!
! W(K)=expzmid*(expzmid*W(K+1)+dz*S1(K))
        do i=lon0,lon1
          do k=lev1-1,lev0,-1
            w(k,i,lat) = expzmid*(expzmid*w(k+1,i,lat)+dz*w_divrg(k,i))
          enddo
        enddo
!       call addfld('W_SWDOT',' ',' ',w(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Filter W (3-d):
      call filter_w(w,lev0,lev1,lon0,lon1,lat0,lat1,'WN')
!
! Pass omega (filtered w) to diags module for WN:
! (mkdiag_WN is expecting geopotential z in cm)
      do lat=lat0,lat1
        call mkdiag_WN('WN',w(:,lon0:lon1,lat),z(:,lon0:lon1,lat),
     |    lev0,lev1,lon0,lon1,lat)
      enddo

      end subroutine swdot
!-----------------------------------------------------------------------
      subroutine filter_w(wout,lev0,lev1,lon0,lon1,lat0,lat1,name)
!
! Filter updated W omega:
!
      use params_module,only: nlat,nlonp4
      use cons_module,only: kut
      use filter_module,only: filter
#ifdef MPI
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi
#else
      integer :: mytidi=0
#endif
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(inout) :: wout(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
      character(len=*),intent(in) :: name
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j
      real :: w_ik(nlonp4,lev0:lev1),w_kij(lev0:lev1,nlonp4,lat0:lat1)
      real :: fmin,fmax
!
#ifdef VT
!     code = 125 ; state = 'filter_w' ; activity='Filtering'
      call vtbegin(125,ier)
#endif
!
! Define lons in w_kij from current task:
      w_kij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          w_kij(:,i,j) = wout(:,i,j)
        enddo
      enddo ! j=lat0,lat1
!
#ifdef MPI
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0 
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(w_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,name)
#endif
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
      if (mytidi==0) then
!
! Define 2d array with all longitudes for filter at each latitude:
        do j=lat0,lat1
          do i=1,nlonp4
            w_ik(i,:) = w_kij(:,i,j)
          enddo ! i=1,nlonp4
!
! Remove wave numbers > kut(lat):
          call filter(w_ik,lev0,lev1,kut(j),j)
!
! Return filtered array to w_kij:
          do i=1,nlonp4
            w_kij(:,i,j) = w_ik(i,:)
          enddo ! i=1,nlonp4
        enddo ! j=lat0,lat1
      endif ! mytidi==0
#ifdef MPI
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(w_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,
     |  name)

!     write(6,"('filter w  after mp_scatterlons: lon0,1=',2i4,' lat0,1='
!    |  ,2i4,' fkij(:,lon0:lon1,:) min,max=',2e12.4)") lon0,lon1,
!    |  lat0,lat1,minval(w_kij(:,lon0:lon1,:)),
!    |            maxval(w_kij(:,lon0:lon1,:))
#endif
!
! Return filtered array to wout at current task longitudes and latitudes:
      do j=lat0,lat1
        do i=lon0,lon1
          wout(:,i,j) = w_kij(:,i,j)
        enddo
      enddo
!
#ifdef VT
!     code = 125 ; state = 'filter_w' ; activity='Filtering'
      call vtend(125,ier)
#endif
      end subroutine filter_w
!-----------------------------------------------------------------------
