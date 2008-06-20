!
      subroutine advec(f,hadvec,lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Horizontal advection for field f. Return advection in hadvec.
!
      use fields_module,only: un,vn,itp
      use cons_module,only: dlamda_2div3,dlamda_1div12,racs,
     |  dphi_1div12,dphi_2div3,re_inv
      use params_module,only: nlonp4
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat-2:lat+2),intent(in) :: 
     |  f ! input field with ghost cells for finite differencing
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: hadvec
!
! Local:
      integer :: k,i,lonbeg,lonend
      integer :: nk,nkm1
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  ubarl,d2x,d4x,wk,vbarp
!
      nk = lev1-lev0+1
      nkm1 = nk-1

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
          ubarl(k,i) = (un(k,i+1,lat,itp)+
     |                  un(k,i-1,lat,itp))*0.5
!
! D2X = (dlamda_2div3*(X(I+1)-X(I-1)))*ubarl
          d2x(k,i) = ((f(k,i+1,lat)-f(k,i-1,lat))*dlamda_2div3)*
     |      ubarl(k,i)
!
! UBARL = (U(I+2)+U(I-2))/2
          ubarl(k,i) = (un(k,i+2,lat,itp)+
     |                  un(k,i-2,lat,itp))*0.5
!
! D4X = (dlamda_1div12*(X(I+2)-X(I-2)))*ubarl
          d4x(k,i) = ((f(k,i+2,lat)-f(k,i-2,lat))*dlamda_1div12)*
     |      ubarl(k,i)
!
! X = (D2X-D4X)*RACS
          wk(k,i) = (d2x(k,i)-d4x(k,i))*racs(lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lonbeg,lonend

!     call addfld('UBARL',' ',' ',ubarl(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('D2X'  ,' ',' ',d2x(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('D4X'  ,' ',' ',d4x(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('WK   ',' ',' ',wk (lev0:lev1-1,:),
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
          vbarp(k,i) = (vn(k,i,lat+1,itp)+
     |                  vn(k,i,lat-1,itp))*0.5
!
! D2X = dphi_2div3*(X(J+1)-X(J-1))
          d2x(k,i) = ((f(k,i,lat+1)-f(k,i,lat-1))*dphi_2div3)*
     |      vbarp(k,i)
!
! Fourth order finite difference in latitude (j+2,j-2):
!
! VBARP = (V(J+2)+V(J-2))/2 
          vbarp(k,i) = (vn(k,i,lat+2,itp)+
     |                  vn(k,i,lat-2,itp))*0.5
!
! D4X = dphi_1div12*(X(J+2)-X(J-2))
          d4x(k,i) = ((f(k,i,lat+2)-f(k,i,lat-2))*dphi_1div12)*
     |      vbarp(k,i)
!
! D2X = (D2X-D4X)*re_inv
          d2x(k,i) = (d2x(k,i)-d4x(k,i))*re_inv
!
! S = X+D2X = ADVEC(X)
          hadvec(k,i) = wk(k,i)+d2x(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lonbeg,lonend
!
! Set hadvec periodic points to zero to avoid NaNS fpe in advecv.
! Earlier versions apparently assumed this, since advecv references
!   the periodic points. (This output from advec becomes input and
!   output in advecv, see dt.F).
!
      if (lon0==1) hadvec(:,lon0:lon0+1) = 0.
      if (lon1==nlonp4) hadvec(:,lon1-1:lon1) = 0.

!     call addfld('VBARP' ,' ',' ',vbarp(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('D4XJ'  ,' ',' ',d4x(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('D2XJ'  ,' ',' ',d2x(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('HADVEC',' ',' ',hadvec(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      end subroutine advec
!-----------------------------------------------------------------------
      subroutine advecv(f,lbc,advecv_out,lev0,lev1,lon0,lon1,lat)
!
! Vertical advection for field f. Return advection in advecv_out.
!
      use params_module,only: dz,nlonp4
      use fields_module,only: w,itc
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  f ! input field with ghost cells for finite differencing
      real,dimension(lon0:lon1),intent(in) :: lbc ! lower boundary
!
! advecv_out is intent(inout) because it is passed in as output from
! previous call to horizontal advection (see dt.F).
      real,dimension(lev0:lev1,lon0:lon1),intent(inout) :: advecv_out
!
! Local:
      integer :: k,i,nk
      real :: advec(lev0:lev1,lon0:lon1) ! s3
      real :: dsig
!
      nk = lev1-lev0+1
!
! Lower boundary:
      do i=lon0,lon1 
        advec(1,i) = (f(1,i)-lbc(i))*w(1,i,lat,itc)*2.
      enddo ! i=lon0,lon1 
      dsig = .5*(1./dz)
!
! Upper boundary:
      do i=lon0,lon1 
        advec(lev1,i) = 0.
!
! Loop through column:
        do k=lev0+1,lev1-1
          advec(k,i) = (f(k,i)-f(k-1,i))*w(k,i,lat,itc)
        enddo ! k=lev0,lev1
        do k=lev0,lev1-1
          advec(k,i) = (advec(k,i)+advec(k+1,i))*dsig
        enddo ! k=lev0,lev1-1
!
! Return vertical advection. Horizontal advection is already in advecv_out
! on input (see call from sub dt).
        do k=lev0,lev1-1
          advecv_out(k,i) = advecv_out(k,i)+advec(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1 
!
! Set periodic points to zero to avoid NaNS fpe in dt.F.
! (earlier versions actually did periodic points here)
      if (lon0==1) advecv_out(:,lon0:lon0+1) = 0.
      if (lon1==nlonp4) advecv_out(:,lon1-1:lon1) = 0.

!     call addfld('ADVECV',' ',' ',advecv_out(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      end subroutine advecv
