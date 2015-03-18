!
      subroutine divrg(un,vc,w,lon0,lon1,lev0,lev1,lat0,lat1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Called from swdot for vertical velocity W:
!
      use cons_module,only: dlamda_2div3,dlamda_1div12,
     |  dphi_2div3,dphi_1div12,racs
      use params_module,only: nlonp4
      use addfld_module,only: addfld
#ifdef MPI
      use mpi_module,only: mp_periodic_f2d
#endif
      implicit none
!
! Args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat0,lat1,lat
      real,intent(in) ::
     |  un(lev0:lev1,lon0-2:lon1+2),
     |  vc(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
      real,intent(out) :: w(lev0:lev1,lon0:lon1)
!
! Local:
      real,dimension(lev0:lev1,lon0:lon1) :: d2u,d4u,d2vc,d4vc
      integer :: k,i,lonbeg,lonend

!     do i=lon0,lon1
!       write(6,"('divrg: lat=',i3,' i=',i3,' un(:,i)=',/,(6e12.4))")
!    |    lat,i,un(:,i)
!     enddo ! i=lon0,lon1

!     call addfld('DVG_UN',' ',' ',un(lev0:lev1,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('DVG_VC',' ',' ',vc(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Do not use un,vc at lons -1, 0, nlonp4+1, or nlonp4+2:
! (define lons 3->nlonp4-2; periodic points are defined later)
      lonbeg = lon0
      if (lon0==1) then
        lonbeg = 3
        w(:,1:2) = 0.
      endif
      lonend = lon1
      if (lon1==nlonp4) then
        lonend = nlonp4-2
        w(:,lonend+1:nlonp4) = 0.
      endif
!
! Calculate divergence for omega (W):
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          d2u(k,i) = (un(k,i+1)-un(k,i-1))*dlamda_2div3    
          d4u(k,i) = (un(k,i+2)-un(k,i-2))*dlamda_1div12
          d2vc(k,i) = (vc(k,i,lat+1)-vc(k,i,lat-1))*dphi_2div3
          d4vc(k,i) = (vc(k,i,lat+2)-vc(k,i,lat-2))*dphi_1div12
          w(k,i) = ((d2u(k,i)-d4u(k,i))+(d2vc(k,i)-d4vc(k,i)))*racs(lat)
        enddo
        
!       write(6,"(/,'divrg: lat=',i3,' i=',i3)") lat,i
!       write(6,"('  un(:,i+1)=',/,(6e12.4))") un(:,i+1)
!       write(6,"('  un(:,i-1)=',/,(6e12.4))") un(:,i-1)
!       write(6,"('  un(:,i+2)=',/,(6e12.4))") un(:,i+2)
!       write(6,"('  un(:,i-2)=',/,(6e12.4))") un(:,i-2)
!       write(6,"('  d2u(:,i)=',/,(6e12.4))") d2u(:,i)
!       write(6,"('  d4u(:,i)=',/,(6e12.4))") d4u(:,i)

      enddo
!
! Set periodic points in fields that need them 
! subroutine mp_periodic_f2d(f,lon0,lon1,lat0,lat1)
#ifdef MPI
! Oct, 2010: Comment these calls to prevent divergence in T,U,V 
!   between runs w/ different numbers of processors (see also minor.F):
!     call mp_periodic_f2d(d2u,lon0,lon1,lat,lat,1)
!     call mp_periodic_f2d(d4u,lon0,lon1,lat,lat,1)
!     call mp_periodic_f2d(w  ,lon0,lon1,lat,lat,1)
#endif
!     call addfld('DVG_D2U' ,' ',' ',d2u(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('DVG_D4U' ,' ',' ',d4u(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('DVG_W' ,' ',' ',w(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      end subroutine divrg
