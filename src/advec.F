!
      subroutine advec(un,vn,f,hadvec,lev0,lev1,lon0,lon1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Horizontal advection for field f. Return advection in hadvec.
!
      use cons_module,only: dlamda_2div3,dlamda_1div12,racs,
     |  dphi_1div12,dphi_2div3,re_inv
      use params_module,only: nlonp4,rp
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) :: un,vn,f ! input field with ghost cells for finite differencing
      real(rp),dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) ::
     |  hadvec
!
! Local:
      integer :: k,i,lat,lonbeg,lonend
      real(rp),dimension(lev0:lev1,lon0:lon1) ::
     |  ubarl,d2x,d4x,wk1,vbarp,wk2
!
      do lat=lat0,lat1

!     call addfld('F_IN',' ',' ',f(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2
!
! Finite difference in longitude:
      do i=lonbeg,lonend
        do k=lev0,lev1-1
!
! UBARL = (U(I+1)+U(I-1))/2
          ubarl(k,i) = (un(k,i+1,lat)+
     |                  un(k,i-1,lat))*0.5_rp
!
! D2X = (dlamda_2div3*(X(I+1)-X(I-1)))*ubarl
          d2x(k,i) = ((f(k,i+1,lat)-f(k,i-1,lat))*dlamda_2div3)*
     |      ubarl(k,i)
!
! UBARL = (U(I+2)+U(I-2))/2
          ubarl(k,i) = (un(k,i+2,lat)+
     |                  un(k,i-2,lat))*0.5_rp
!
! D4X = (dlamda_1div12*(X(I+2)-X(I-2)))*ubarl
          d4x(k,i) = ((f(k,i+2,lat)-f(k,i-2,lat))*dlamda_1div12)*
     |      ubarl(k,i)
!
! X = (D2X-D4X)*RACS
          wk1(k,i) = (d2x(k,i)-d4x(k,i))*racs(lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lonbeg,lonend

!     call addfld('UBARL',' ',' ',ubarl(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('D2X'  ,' ',' ',d2x(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('D4X'  ,' ',' ',d4x(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('WK   ',' ',' ',wk1(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
!
! Finite difference in latitude:
      do i=lonbeg,lonend
        do k=lev0,lev1-1
!
! 2nd order finite difference in latitude (j+1,j-1):
!
! (10/21/03 btf: 2nd order diffs (next 2 statements) were missing 
!    in tgcm13,14,15, and tiegcm1 until this date.)
!
! VBARP = (V(J+1)+V(J-1))/2
          vbarp(k,i) = (vn(k,i,lat+1)+
     |                  vn(k,i,lat-1))*0.5_rp
!
! D2X = dphi_2div3*(X(J+1)-X(J-1))
          d2x(k,i) = ((f(k,i,lat+1)-f(k,i,lat-1))*dphi_2div3)*
     |      vbarp(k,i)
!
! Fourth order finite difference in latitude (j+2,j-2):
!
! VBARP = (V(J+2)+V(J-2))/2 
          vbarp(k,i) = (vn(k,i,lat+2)+
     |                  vn(k,i,lat-2))*0.5_rp
!
! D4X = dphi_1div12*(X(J+2)-X(J-2))
          d4x(k,i) = ((f(k,i,lat+2)-f(k,i,lat-2))*dphi_1div12)*
     |      vbarp(k,i)
!
! Y = (D2X-D4X)*re_inv
          wk2(k,i) = (d2x(k,i)-d4x(k,i))*re_inv
!
! S = X+Y = ADVEC(X,Y)
          hadvec(k,i,lat) = wk1(k,i)+wk2(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lonbeg,lonend
!
! Set hadvec periodic points to zero to avoid NaNS fpe in advecv.
! Earlier versions apparently assumed this, since advecv references
!   the periodic points. (This output from advec becomes input and
!   output in advecv, see dt.F).
!
      if (lon0==1) hadvec(:,lon0:lon0+1,lat) = 0._rp
      if (lon1==nlonp4) hadvec(:,lon1-1:lon1,lat) = 0._rp

!     call addfld('VBARP' ,' ',' ',vbarp(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('D4XJ'  ,' ',' ',d4x(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('D2XJ'  ,' ',' ',d2x(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('HADVEC',' ',' ',hadvec(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      enddo

      end subroutine advec
