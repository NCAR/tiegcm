      module qjoule
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      contains
!-----------------------------------------------------------------------
      subroutine qjoule_ti(un,vn,w,ui,vi,wi,lam1,op,nplus,
     |    n2p,nop,o2p,barm,tn,qji_ti,lev0,lev1,lon0,lon1,lat)
!
! Calculate ion joule heating for use in ion temperature.
! This is called from sub dynamics before settei. Qji_ti output from
!   here is passed to settei for ion temperature calculation.
! In earlier model versions, this calculation was at the end of settei.
! See also sub qjoule_tn below for neutral temperature.
!
! am 6/05: correct calculation of Joule heating going into the ions
!   - include also vertical velocity by using lambda_1 (mag. direction)
!      and vertical velocity difference
!   - using the ration between the mean molecular weights
!
!     before
!       Q_JI = [ lambda_xx * u_i(u_i-u_n) + lambda_yy*v_i*(v_i-v_n)]
!     now
!       Q_JI = m_n/(m_n+m_i) * lambda_1 * [ (u_i-u_n)^2 + (v_i-v_n)^2+ (w_i-w_n)^2]
!
      use cons_module,only: avo,        ! avogadro number 6.023e23 particles / mole
     |  grav,   ! acceler. due to gravity at lower boundary [cm/s^2]
     |  gask,   ! gas konstant 8.314e7 [ergs/K/mol]
     |  rmass_o2 , rmass_o1, rmass_n2,  ! 32,16,28 [g/mol]
     |  rmass_n4s, rmass_no, rmass_op   ! 14,30,16 [g/mol]
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  un,vn,          ! zonal, meridional  neutral velocity [cm/s]
     |  w,              ! dimensionless vertical  neutral velocity omega [1/s] (interface)
     |  ui,vi,wi,       ! zonal, meridional, vertical ion velocity (ExB) [cm/s] (interface)
     |  lam1,           ! lambda_1 ion drag coefficients [1/s] (interface)
     |  op,        	! O+ 1/cm^3  at half levels
     |  nplus,     	! N+  1/cm^3 at half levels
     |  n2p,       	! N2+ 1/cm^3 at half levels
     |  nop,       	! NO+  1/cm^3 at half levels
     |  o2p,        	! O2+  1/cm^3 at half levels
     |  barm,           ! mean molecular weight (interface)
     |  tn              ! neutral temperature (deg K)
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  qji_ti          ! ion joule heating for ti (output) [ergs/s/g]
!
! Local:
      integer :: k,i
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  vert_vel,        ! scale height  
     |  uii,vii,wii      ! ion velocities at interfaces
      real :: m_i,m_n,mfac
!
! Task subdomain:
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! mean molecular weight of ion at half level [g/mol]
          m_i = op(k,i)*rmass_o1 + o2p(k,i)*rmass_o2+  
     |          nplus(k,i)*rmass_n4s+ n2p(k,i)*rmass_n2+
     |          nop(k,i)*rmass_no
          m_i = m_i/(op(k,i) + o2p(k,i)+nplus(k,i)+ n2p(k,i)+
     |          nop(k,i))
! mean molecular weight of neutrals at half level [g/mol]
          m_n = .5*(barm(k,i)+barm(k+1,i))
! ratio m_n/(m_n+m_i)	
	  mfac = m_i+m_n
	  mfac = m_n/mfac
	  
! scale height H [cm] to convert omega [1/s] to vertical neutral wind [cm/s]
          vert_vel(k,i) = 0.5*(w(k,i)+w(k+1,i))*gask*tn(k,i)/
     |        (.5*(barm(k,i)+barm(k+1,i))*grav)
     
! Ion velocities at half levels
          uii(k,i) = .5*(ui(k,i)+ui(k+1,i))  ! s6
          vii(k,i) = .5*(vi(k,i)+vi(k+1,i))  ! s5
          wii(k,i) = .5*(wi(k,i)+wi(k+1,i))  ! s5
!
! Joule heating (s7): [ergs/s/g]
! (note qji_ti at lev1 not defined)
          qji_ti(k,i) = (uii(k,i)-un(k,i))**2+(vii(k,i)-vn(k,i))**2 +
     |      (wii(k,i)-vert_vel(k,i))**2
          qji_ti(k,i) = mfac*.5*(lam1(k,i)+lam1(k+1,i))*qji_ti(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('QJI_TI','QJI_TI ',' ',qji_ti(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      end subroutine qjoule_ti
!-----------------------------------------------------------------------
      subroutine qjoule_tn(un,vn,w,ui,vi,wi,lam1,barm,tn,qji_tn,
     |  lev0,lev1,lon0,lon1,lat)
!
! Calculate ion joule heating for neutral temperature.
! This is called from dynamics before dt.
! (in earlier model versions, this was in dt.F)
!
! am 6/05: include the vertical component of the velocities
!   - include also vertical velocity by using lambda_1 (mag. direction)
!      and vertical velocity difference
!     change from
!       Q_J = [ lambda_xx * (u_i-u_n)^2 + lambda_yy*(v_i-v_n)^2
!             (lambda_xy-lambda_yx)* (u_i-u_n)*(v_i-v_n)]
!     to
!       Q_J = lambda_1 * [ (u_i-u_n)^2 + (v_i-v_n)^2 + (w_i-w_n)^2]
!
      use cons_module,only:
     |  grav,   ! acceler. due to gravity at lower boundary [cm/s^2]
     |  gask    ! gas konstant 8.314e7 [ergs/K/mol]
      use addfld_module,only: addfld
      use diags_module,only: mkdiag_QJOULE,mkdiag_QJOULE_INTEG
      use input_module,only: joulefac_input => joulefac
      implicit none
!     
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  un,vn,          ! zonal, meridional  neutral velocity [cm/s]
     |  w,              ! dimensionless vertical  neutral velocity omega [1/s] (interface)
     |  ui,vi,wi,       ! zonal, meridional, vertical ion velocity (interface)
     |  lam1,           ! lambda_1 ion drag coefficient
     |  barm,           ! mean molecular weight (interface)
     |  tn              ! neutral temperature (deg K)
       real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  qji_tn          ! ion joule heating for tn (output)
!
! Local:
      integer :: k,i
      real,dimension(lev0:lev1,lon0:lon1) ::vel_zonal,vel_merid,
     |  vel_vert,scheight        ! scale height  
      real :: joulefac
!
#if defined(INTERCOMM) || defined(CISMAH)
!
! CMIT: When copled to the LFM, the currents are strong
!       so we need to reduce joulefac from 1.5 to 1.0:
!
      joulefac = 1.0 ! joule heating multiplication factor
#else
!
! Assign joulefac from value read from namelist (or its default):
!
      joulefac = joulefac_input
#endif
! 
! scheight: scale height H [cm] to convert omega [1/s] to vertical 
! neutral wind [cm/s] on half level	  
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          scheight(k,i) = gask*tn(k,i)/
     |        (.5*(barm(k,i)+barm(k+1,i))*grav)
          vel_zonal(k,i) = .5*(ui(k,i)+ui(k+1,i))-un(k,i)  ! s2
          vel_merid(k,i) = .5*(vi(k,i)+vi(k+1,i))-vn(k,i)  ! s3
          vel_vert(k,i)  = .5*(wi(k,i)+wi(k+1,i)-scheight(k,i)*
     |         ( w(k,i)-w(k+1,i)) )
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
      do i=lon0,lon1
        do k=lev0,lev1-1
          qji_tn(k,i) = .5*(lam1(k,i)+lam1(k+1,i))*
     |      (vel_zonal(k,i)**2 + vel_merid(k,i)**2 + 
     |       vel_vert(k,i)**2)
!
! Apply joule heating factor (formerly in dt.F):
          qji_tn(k,i) = qji_tn(k,i) * joulefac
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!     call addfld('QJI_UN','QJI_UN','cm/s',un(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QJI_VN','QJI_VN','cm/s',vn(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QJI_W','QJI_W','cm/s',vn(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QJI_UI','QJI_UI','cm/s',ui(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QJI_VI','QJI_VI','cm/s',vi(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QJI_WI','QJI_WI','cm/s',wi(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QJI_LAM1','QJI_LAM','1/s',lam1(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QJI_TN','QJI_TN','ergs/s/g',
!    |  qji_tn(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Save diagnostic Total Joule Heating:
      call mkdiag_QJOULE('QJOULE',qji_tn,lev0,lev1,lon0,lon1,lat)
!
! Save diagnostic Height-integrated Joule Heating:
      call mkdiag_QJOULE_INTEG('QJOULE_INTEG',qji_tn(:,lon0:lon1),
     |  lev0,lev1,lon0,lon1,lat)

      end subroutine qjoule_tn
!-----------------------------------------------------------------------
      end module qjoule
