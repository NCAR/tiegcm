!
      subroutine qjoule_ti(un,vn,ui,vi,lxx,lyy,lxy,lyx,qji_ti,
     |  lev0,lev1,lon0,lon1,lat)
!
! Calculate ion joule heating for use in ion temperature.
! This is called from sub dynamics before settei. Qji_ti output from
!   here is passed to settei for ion temperature calculation.
! In earlier model versions, this calculation was at the end of settei.
! See also sub qjoule_tn below for neutral temperature.
!
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  un,vn,          ! zonal, meridional neutral velocity
     |  ui,vi,          ! zonal, meridional ion velocity
     |  lxx,lyy,lxy,lyx ! lambda ion drag coefficients
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  qji_ti          ! ion joule heating for ti (output)
!
! Local:
      integer :: k,i,nlevs
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  uii,vii ! ion velocities at interfaces
!
      nlevs = lev1-lev0+1

!     call addfsech('QJI_UI',' ',' ',ui(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QJI_VI',' ',' ',vi(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QJI_LXX',' ',' ',lxx(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QJI_LYY',' ',' ',lyy(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QJI_LXY',' ',' ',lxy(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QJI_LYX',' ',' ',lyx(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!
! Task subdomain:
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! Ion velocities at interfaces (k+1/2)
          uii(k,i) = .5*(ui(k,i)+ui(k+1,i))  ! s6
          vii(k,i) = .5*(vi(k,i)+vi(k+1,i))  ! s5
!
! Joule heating (s7):
! (note qji_ti at lev1 not defined)
          qji_ti(k,i) = 
     |      .5*((lxx(k,i)+lxx(k+1,i))*uii(k,i)*(uii(k,i)-un(k,i))+
     |          (lxy(k,i)+lxy(k+1,i))*uii(k,i)*(vii(k,i)-vn(k,i))-
     |          (lyx(k,i)+lyx(k+1,i))*vii(k,i)*(uii(k,i)-un(k,i))+
     |          (lyy(k,i)+lyy(k+1,i))*vii(k,i)*(vii(k,i)-vn(k,i)))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('QJI_TI',' ',' ',qji_ti(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)

      end subroutine qjoule_ti
!-----------------------------------------------------------------------
      subroutine qjoule_tn(un,vn,ui,vi,lxx,lyy,lxy,lyx,qji_tn,
     |  lev0,lev1,lon0,lon1,lat)
!
! Calculate ion joule heating for neutral temperature.
! This is called from dynamics before dt.
! (in earlier model versions, this was in dt.F)
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  un,vn,          ! zonal, meridional neutral velocity
     |  ui,vi,          ! zonal, meridional ion velocity
     |  lxx,lyy,lxy,lyx ! lambda ion drag coefficients
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  qji_tn          ! ion joule heating for tn (output)
!
! Local:
      integer :: k,i,nlevs
      real,dimension(lev0:lev1,lon0:lon1) :: vel_zonal,vel_merid
      integer :: iuivi=1
!
! iuivi used to be an input flag, but has been replaced by the dynamo
!   input flag. iuivi is set locally here, so ui,vi are always used in
!   joule heating calculation. But if dynamo==0 (no dynamo), then ui,vi 
!   will be zero, resulting in the same result as if iuivi==0.
!
      if (iuivi==0) then ! never executed, assuming local iuivi==1.
        do i=lon0,lon1
          do k=lev0,lev1-1
            vel_zonal(k,i) = -un(k,i)   ! s2
            vel_merid(k,i) = -vn(k,i)   ! s3
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      else
        do i=lon0,lon1
          do k=lev0,lev1-1
            vel_zonal(k,i) = .5*(ui(k,i)+ui(k+1,i))-un(k,i)  ! s2
            vel_merid(k,i) = .5*(vi(k,i)+vi(k+1,i))-vn(k,i)  ! s3
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      endif
!
!       S4(I,1) = .5*(S2(I,1)**2*(F(I,NLXXK)+F(I,NLXXK+1))+S2(I,1)*
!    1    S3(I,1)*(F(I,NLXYK)+F(I,NLXYK+1)-F(I,NLYXK)-F(I,NLYXK+1))+
!    2    S3(I,1)**2*(F(I,NLYYK)+F(I,NLYYK+1)))

      do i=lon0,lon1
        do k=lev0,lev1-1
          qji_tn(k,i) = 
     |      .5*(vel_zonal(k,i)**2*(lxx(k,i)+lxx(k+1,i))+
     |          vel_zonal(k,i)*vel_merid(k,i)*
     |         (lxy(k,i)+lxy(k+1,i)-lyx(k,i)-lyx(k+1,i))+
     |          vel_merid(k,i)**2*(lyy(k,i)+lyy(k+1,i)))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

      nlevs = lev1-lev0+1
!
! xlf8.1 bug: QJI_LXX has spike at geog pole in 2nd time step
! (tgcm15 does not)
!
!     call addfsech('QJI_UN',' ',' ',un(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('QJI_VN',' ',' ',vn(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('QJI_UI',' ',' ',ui(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('QJI_VI',' ',' ',vi(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('QJI_LXX',' ',' ',lxx(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('QJI_LYY',' ',' ',lyy(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('QJI_LXY',' ',' ',lxy(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('QJI_LYX',' ',' ',lyx(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('QJI_TN',' ',' ',qji_tn(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs-1,lat)

      end subroutine qjoule_tn
