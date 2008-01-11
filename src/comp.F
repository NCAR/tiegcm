!
      subroutine comp(tn,o2,o2_nm,o1,o1_nm,un,vn,w,hdo2,hdo1,
     |                o2_upd,o2nm_upd,o1_upd,o1nm_upd,
     |                lev0,lev1,lon0,lon1,lat0,lat1)
!
! Advance major species O2 and O.
!
      use params_module,only: dz,nlonp4,nlat,dlat
      use input_module,only: difhor
      use init_module,only: glat
      use cons_module,only: pi,hor,dtr,rmassinv_o2,rmassinv_o1,
     |  rmassinv_n2,rmass_o2,rmass_o1,expz,expzmid,expzmid_inv,
     |  difk,dtx2inv,dtsmooth,dtsmooth_div2
      use chemrates_module,only: fs ! from sub comp_o2o
      use bndry_module,only: b,fb
      use fields_module,only: tlbc
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,   ! neutral temperature
     |  o2,   ! O2 (mmr) at current timestep
     |  o1,   ! O  (mmr) at current timestep
     |  o2_nm,! O2 (mmr) at time n-1
     |  o1_nm,! O  (mmr) at time n-1
     |  un,   ! zonal wind velocity at current timestep
     |  vn,   ! meridional wind velocity at current timestep
     |  w,    ! vertical velocity at current timestep
     |  hdo2, ! O2 horizontal diffusion (hdif3)
     |  hdo1  ! O  horizontal diffusion (hdif3)
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out)::
     |  o2_upd   ,o1_upd,  ! output: O2,O updated for next timestep
     |  o2nm_upd ,o1nm_upd ! output: O2,O updated for previous timestep
!
! Local:
      integer :: k,kk,i,lat,isp,km,kp,ktmp,m,lonbeg,lonend
      integer :: nk,nkm1,i0,i1,k0,k1
      integer,parameter :: io2=1,io1=2 ! indices to O2, O, respectively
      real,dimension(lon0:lon1,lev0:lev1,2,2) :: gama
      real,dimension(lon0:lon1,lev0:lev1,2)   :: zz
      real,dimension(lon0:lon1,lev0:lev1)     :: embar
      real,dimension(lon0:lon1,2,2,2)         :: ak
      real,dimension(lon0:lon1,2,2)           :: ep,pk,qk,rk,wkm1,wkm2
      real,dimension(lon0:lon1,2)             :: fk,wkv1,wkv2,ps0
      real,dimension(lon0:lon1)               :: wks1,wks2,wks3,wks4,
     |                                           embar0,dfactor
      real :: rlat
      real :: ak0(2,2),phi(2,3),delta(2,2),tau,t00,small
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  o2nm_smooth, o1nm_smooth,  ! smoothed at time n-1
     |  o2_advec   , o1_advec      ! horizontal advection
      real,dimension(lev0:lev1,lon0:lon1,2) :: upd
! For diagnostics:
      real,dimension(lev0:lev1,lon0:lon1,2) :: zz_ki
      real,dimension(lev0:lev1,lon0:lon1,2,2) :: gama_ki
!
      phi(:,1)=(/0.  ,0.673/)
      phi(:,2)=(/1.35,0.   /)
      phi(:,3)=(/1.11,0.769/)
      tau=1.86e+3
      delta(:,1)=(/1.,0./)
      delta(:,2)=(/0.,1./)
      t00=273.
      small = 1.e-6
      i0=lon0 ; i1=lon1 ; k0=lev0 ; k1=lev1
!
      nk = lev1-lev0+1
      nkm1 = nk-1
!
! Calculate and save horizontal advection in o2_advec, o1_advec:
!
      do lat=lat0,lat1
        call advecl(o2,o1,un,vn,o2_advec,o1_advec,
     |    lev0,lev1,lon0,lon1,lat0,lat1,lat)
      enddo ! lat=lat0,lat1
!
! Save smoothed o2,o at time n-1:
!
      call smooth(o2_nm,o2nm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
      call smooth(o1_nm,o1nm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)

!     do lat=lat0,lat1
!       call addfld('O2SMOOTH','O2SMOOTH',' ',o2nm_smooth(:,i0:i1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!       call addfld('O1SMOOTH','O1SMOOTH',' ',o1nm_smooth(:,i0:i1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!     enddo ! lat=lat0,lat1
!
! Begin latitude scan:
      do lat=lat0,lat1
!
      do i=lon0,lon1
        dfactor(i) = 1.
      enddo
!
! If difhor flag was by user (see input_module):
! dfactor = 
!   .5*(1.+SIN(PI*(ABS(RLATM)-PI/6.)/(PI/3.))) FOR ABS(RLATM).LT.PI/3.
! dfactor = 1. FOR  ABS(RLATM).GE.PI/3
! (dfactor was in sub dfact in earlier versions)
!
      if (difhor > 0) then
        rlat = (glat(1)+(lat-1)*dlat)*dtr
        if (abs(rlat)-pi/4.5 >= 0.) then
          dfactor(:) = hor(lat)+1.
        else
          dfactor(:) = hor(lat)+.5*(1.+sin(pi*(abs(rlat)-pi/9.)/
     |      (pi/4.5)))
        endif
      else
        dfactor(:) = 1.
      endif
!     write(6,"('comp: lat=',i2,' dfactor=',/,(6e12.4))")
!    |  lat,dfactor
!
! Embar:
      do i=lon0,lon1
        do k=lev0,lev1
          embar(i,k) = 1./(o2(k,i,lat)*rmassinv_o2 + 
     |                     o1(k,i,lat)*rmassinv_o1 +
     |     (1.-o2(k,i,lat)-o1(k,i,lat))*rmassinv_n2)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

! real,dimension(lon0:lon1,lev0:lev1)     :: embar
!     call addfld('EMBAR',' ',' ',embar(i0:i1,:),
!    |  'lon',i0,i1,'lev',k0,k1,lat)
!
! ps0 and embar0:
      do i=lon0,lon1
        ps0(i,io2) = 
     |    b(i,1,1)*o2(lev0,i,lat)+b(i,1,2)*o1(lev0,i,lat)+fb(i,1)
        ps0(i,io1) = 
     |    b(i,2,1)*o2(lev0,i,lat)+b(i,2,2)*o1(lev0,i,lat)+fb(i,2)
        embar0(i) = 1./(ps0(i,io2)*rmassinv_o2+ps0(i,io1)*rmassinv_o1+
     |              (1.-ps0(i,io2)-ps0(i,io1))*rmassinv_n2)
! WKS4 = .5*(DMBAR/DZ)/MBAR
        wks4(i) = (embar(i,lev0)-embar0(i))/
     |            (dz*(embar0(i)+embar(i,lev0)))
      enddo ! i=lon0,lon1

!     write(6,"(/,'comp: lat=',i2,' ps0(:,io2)=',/,(6e12.4))") 
!    |  lat,ps0(:,io2)
!     write(6,"('ps0(:,io1)=',/,(6e12.4))") ps0(:,io1)
!     write(6,"('embar0(:)=',/,(6e12.4))") embar0(:)
!     write(6,"('wks4(:)=',/,(6e12.4))") wks4(:)
!
! ep, ak at level 1/2:
      km = 1 
      kp = 2
      do i=lon0,lon1
        ep(i,io2,kp) = 1.-(2./(embar0(i)+embar(i,lev0)))*
     |               (rmass_o2+(embar(i,lev0)-embar0(i))/dz)
        ep(i,io1,kp) = 1.-(2./(embar0(i)+embar(i,lev0)))*
     |               (rmass_o1+(embar(i,lev0)-embar0(i))/dz)
        zz(i,1,:) = 0.
      enddo ! i=lon0,lon1

!     write(6,"('ep(:,1,kp)=',/,(6e12.4))") ep(:,1,kp)
!     write(6,"('ep(:,2,kp)=',/,(6e12.4))") ep(:,2,kp)

      do m=1,2
        do i=lon0,lon1
          ak(i,io2,m,kp) = 
     |      -delta(io2,m)*(phi(io1,3)+(phi(io1,io2)-phi(io1,3))*
     |      .5*(ps0(i,io2)+o2(lev0,i,lat)))-(1.-delta(io2,m))*
     |      (phi(io2,m)-phi(io2,3))*.5*(ps0(i,io2)+o2(lev0,i,lat))
          ak(i,io1,m,kp) = 
     |      -delta(io1,m)*(phi(io2,3)+(phi(io2,io1)-phi(io2,3))*
     |      .5*(ps0(i,io1)+o1(lev0,i,lat)))-(1.-delta(io1,m))*
     |      (phi(io1,m)-phi(io1,3))*.5*(ps0(i,io1)+o1(lev0,i,lat))
        enddo ! i=lon0,lon1

!       write(6,"('lat=',i2,' m=',i2,' ak(io2)=',/,(6e12.4))")
!    |    lat,m,ak(:,io2,m,kp)
!       write(6,"('lat=',i2,' m=',i2,' ak(io1)=',/,(6e12.4))")
!    |    lat,m,ak(:,io1,m,kp)

      enddo ! m=1,2
!
! WKS1=MBAR/M3*(T00/(T0+T))*0.25/(TAU*DET) AT LEVEL 1/2
! tn lower boundary is stored in top slot tn(lev1..). 
      do i=lon0,lon1
!       wks1(i) = 0.5*(embar0(i)+embar(i,lev0))*rmassinv_n2*
!    |    (t00/tn(lev1,i,lat))**0.25/(tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-
!    |    ak(i,1,2,kp)*ak(i,2,1,kp)))
! Lower boundary of tn is now in tlbc:
        wks1(i) = 0.5*(embar0(i)+embar(i,lev0))*rmassinv_n2*
     |    (t00/tlbc(i,lat))**0.25/(tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-
     |    ak(i,1,2,kp)*ak(i,2,1,kp)))
      enddo ! i=lon0,lon1
!     write(6,"('comp: lat=',i3,' wks1=',/,(6e12.4))") lat,wks1
!
! Complete claculation of ak(1/2)
      do m=1,2
        do i=lon0,lon1
          ak(i,io2,m,kp) = ak(i,io2,m,kp)*wks1(i)
          ak(i,io1,m,kp) = ak(i,io1,m,kp)*wks1(i)
          gama(i,lev0,:,m) = 0.
        enddo ! i=lon0,lon1

!       write(6,"('lat=',i2,' m=',i2,' ak(io2)=',/,(6e12.4))")
!    |    lat,m,ak(:,io2,m,kp)
!       write(6,"('lat=',i2,' m=',i2,' ak(io1)=',/,(6e12.4))")
!    |    lat,m,ak(:,io1,m,kp)

      enddo ! m=1,2
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Height (pressure) loop:
! For now (4/02), put k-loop on outside even tho embar and input
!   fields are (k,i), for convenience in verification with tgcm15.
!
      km = 1  ! alternates 2,1,2,1,... during k-loop
      kp = 2  ! alternates 1,2,1,2,... during k-loop
      levloop: do k=lev0,lev1-1                               ! DO 6
        ktmp = km
        km   = kp
        kp   = ktmp
        do i=lon0,lon1
          ep(i,io2,kp) = 1.-(2./(embar(i,k)+embar(i,k+1)))*(rmass_o2+
     |      (embar(i,k+1)-embar(i,k))/dz)
          ep(i,io1,kp) = 1.-(2./(embar(i,k)+embar(i,k+1)))*(rmass_o1+
     |      (embar(i,k+1)-embar(i,k))/dz)
        enddo ! i=lon0,lon1

!       write(6,"('lat=',i3,' k=',i3,' ep(io2)=',/,(6e12.4))") lat,k,
!    |    ep(:,io2,kp)
!       write(6,"('ep(io1)=',/,(6e12.4))") ep(:,io1,kp)

        do m=1,2
          do i=lon0,lon1
!
! AK(K+1/2)
            ak(i,io2,m,kp) = 
     |        -delta(io2,m)*(phi(io1,3)+(phi(io1,io2)-phi(io1,3))*
     |        .5*(o2(k,i,lat)+o2(k+1,i,lat)))-
     |        (1.-delta(io2,m))*(phi(io2,m)-phi(io2,3))*
     |        .5*(o2(k,i,lat)+o2(k+1,i,lat))

            ak(i,io1,m,kp) = 
     |        -delta(io1,m)*(phi(io2,3)+(phi(io2,io1)-phi(io2,3))*
     |        .5*(o1(k,i,lat)+o1(k+1,i,lat)))-
     |        (1.-delta(io1,m))*(phi(io1,m)-phi(io1,3))*
     |        .5*(o1(k,i,lat)+o1(k+1,i,lat))

          enddo ! i=lon0,lon1
        enddo ! m=1,2
!
! WKS1=MBAR/M3*(T00/(T0+T))**0.25/(TAU*DET(ALFA))
        do i=lon0,lon1
          wks1(i) = 0.5*(embar(i,k)+embar(i,k+1))*rmassinv_n2*
     |      (t00/(.5*(tn(k,i,lat)+tn(k+1,i,lat))))**0.25/
     |      (tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-ak(i,1,2,kp)*
     |      ak(i,2,1,kp)))
!
! EDDY DIFFUSION TERMS IN WKS3 AND WKS4 
          wks3(i) = wks4(i)
          wks4(i) = (embar(i,k+1)-embar(i,k))/
     |      (dz*(embar(i,k)+embar(i,k+1)))
        enddo ! i=lon0,lon1
!
! FINISH CALCULATING AK(K+1/2) AND GENERATE PK, QK, RK
        do m=1,2
          do isp=io2,io1
            do i=lon0,lon1
              ak(i,isp,m,kp) = ak(i,isp,m,kp)*wks1(i)

              pk(i,isp,m) = (ak(i,isp,m,km)*(1./dz+ep(i,m,km)/2.)-
     |          expz(k)*(expzmid_inv*difk(k)*dfactor(i)*(1./dz-
     |          wks3(i))+0.25*(w(k,i,lat)+w(k+1,i,lat)))*
     |          delta(isp,m))/dz

              rk(i,isp,m) = (ak(i,isp,m,kp)*(1./dz-ep(i,m,kp)/2.)-
     |          expz(k)*(expzmid*difk(k+1)*dfactor(i)*(1./dz+
     |          wks4(i))-0.25*(w(k,i,lat)+w(k+1,i,lat)))*
     |          delta(isp,m))/dz
   
              qk(i,isp,m) = -(ak(i,isp,m,km)*(1./dz-ep(i,m,km)/2.)+
     |          ak(i,isp,m,kp)*(1./dz+ep(i,m,kp)/2.))/dz+expz(k)*
     |          (((expzmid*difk(k+1)*(1./dz-wks4(i))+expzmid_inv*
     |          difk(k)*(1./dz+wks3(i)))*dfactor(i)/dz+dtx2inv)*
     |          delta(isp,m)-fs(i,k,isp,m,lat))

            enddo ! i=lon0,lon1

!           write(6,"(/,'comp: lat=',i3,' m=',i2,' isp=',i2,' k=',i3)") 
!    |        lat,m,isp,k
!           write(6,"('ak=',/,(6e12.4))") ak(:,isp,m,kp)
!           write(6,"('pk=',/,(6e12.4))") pk(:,isp,m)
!           write(6,"('rk=',/,(6e12.4))") rk(:,isp,m)
!           write(6,"('qk=',/,(6e12.4))") qk(:,isp,m)

          enddo ! isp=io2,io1
        enddo ! m=1,2
!
! Use advection saved from advecl calls at beginning of routine:
!
        do i=lon0,lon1
          fk(i,io2) = o2_advec(k,i,lat)
          fk(i,io1) = o1_advec(k,i,lat)
        enddo ! i=lonbeg,lonend

!       write(6,"(/,'comp: lat=',i3,' k=',i2)") lat,k
!       write(6,"('fk(:,io2)=',/,(6e12.4))") fk(lonbeg:lonend,io2)
!       write(6,"('fk(:,io1)=',/,(6e12.4))") fk(lonbeg:lonend,io1)
!
! Add explicit source terms to fk:
        do i=lon0,lon1
          fk(i,io2) = fk(i,io2)-fs(i,k,io2,0,lat)
          fk(i,io1) = fk(i,io1)-fs(i,k,io1,0,lat)
        enddo ! i=lon0,lon1

!       write(6,"(/,'comp: lat=',i3,' k=',i2)") lat,k
!       write(6,"('fk(:,io2)=',/,(6e12.4))") fk(lonbeg:lonend,io2)
!       write(6,"('fk(:,io1)=',/,(6e12.4))") fk(lonbeg:lonend,io1)
!
! Complete calculation of rhs in fk:
        do i=lonbeg,lonend
          fk(i,io2) = expz(k)*(o2nm_smooth(k,i,lat)*dtx2inv-fk(i,io2)+
     |      hdo2(k,i,lat))
          fk(i,io1) = expz(k)*(o1nm_smooth(k,i,lat)*dtx2inv-fk(i,io1)+
     |      hdo1(k,i,lat))
        enddo ! i=lonbeg,lonend

!       write(6,"(/,'comp: lat=',i3,' k=',i2)") lat,k
!       write(6,"('fk(:,io2)=',/,(6e12.4))") fk(lonbeg:lonend,io2)
!       write(6,"('fk(:,io1)=',/,(6e12.4))") fk(lonbeg:lonend,io1)

!
! fk is ok up to this point.
! In earlier version, periodic points for fk were taken here.
! For now, ignore periodic points.


!
! Lower boundary:
        if (k==lev0) then
          do m=1,2                                                ! DO 16
            do kk=1,2
              do i=lon0,lon1
                qk(i,io2,m) = qk(i,io2,m)+pk(i,io2,kk)*b(i,kk,m) 
                qk(i,io1,m) = qk(i,io1,m)+pk(i,io1,kk)*b(i,kk,m) 
              enddo ! i=lon0,lon1
            enddo ! kk=1,2
          enddo ! m=1,2
          do m=1,2
            do i=lon0,lon1
              fk(i,io2) = fk(i,io2)-pk(i,io2,m)*fb(i,m)
              fk(i,io1) = fk(i,io1)-pk(i,io1,m)*fb(i,m)
              pk(i,:,m) = 0.
            enddo ! i=lon0,lon1
          enddo ! m=1,2

!         do m=1,2
!           write(6,"('comp lbc: m=',i2,' lat=',i2)") m,lat
!           write(6,"('qk(io2)=',/,(6e12.4))") qk(:,io2,m)
!           write(6,"('qk(io1)=',/,(6e12.4))") qk(:,io1,m)
!         enddo ! m=1,2
!         write(6,"('fk(io2)=',/,(6e12.4))") fk(:,io2)
!         write(6,"('fk(io1)=',/,(6e12.4))") fk(:,io1)

!
! Upper boundary:
        elseif (k==lev1-1) then
          do m=1,2
            do i=lon0,lon1
              qk(i,io2,m) = qk(i,io2,m)+(1.+.5*ep(i,m,kp)*dz)/
     |          (1.-.5*ep(i,m,kp)*dz)*rk(i,io2,m)
              qk(i,io1,m) = qk(i,io1,m)+(1.+.5*ep(i,m,kp)*dz)/
     |          (1.-.5*ep(i,m,kp)*dz)*rk(i,io1,m)
              rk(i,:,m) = 0.
            enddo ! i=lon0,lon1
          enddo ! m=1,2

!         do m=1,2
!           write(6,"('comp ubc: m=',i2,' lat=',i2)") m,lat
!           write(6,"('qk(io2)=',/,(6e12.4))") qk(:,io2,m)
!           write(6,"('qk(io1)=',/,(6e12.4))") qk(:,io1,m)
!         enddo ! m=1,2

        endif ! lbc or ubc
!
! QK=ALFAK=QK-PK*GAMA(K-1)
        do m=1,2                                          ! DO 18
          do kk=1,2
            do i=lon0,lon1

!             write(6,"('comp: i=',i2,' kk=',i2,' m=',i2,' k=',i2,
!    |          ' lat=',i2)") i,kk,m,k,lat
!             write(6,"('qk=',e12.4,' pk=',e12.4,' gama=',e12.4))") 
!    |          qk(i,io2,m),pk(i,io2,kk),gama(i,k,kk,m)

              qk(i,io2,m) = qk(i,io2,m)-pk(i,io2,kk)*gama(i,k,kk,m)
              qk(i,io1,m) = qk(i,io1,m)-pk(i,io1,kk)*gama(i,k,kk,m)
            enddo ! i=lon0,lon1
          enddo ! kk=1,2
        enddo ! m=1,2

! Testing will not work here until gama is incremented below..
!       do m=1,2
!         write(6,"('comp: m=',i2,' k=',i2,' lat=',i2)") m,k,lat
!         write(6,"('qk(io2)=',/,(6e12.4))") qk(:,io2,m)
!         write(6,"('qk(io1)=',/,(6e12.4))") qk(:,io1,m)
!       enddo ! m=1,2

!
! WKS1=DET(ALFA)
        do i=lon0,lon1
          wks1(i) = qk(i,1,1)*qk(i,2,2)-qk(i,1,2)*qk(i,2,1)
        enddo ! i=lon0,lon1
!
! WKM1=ALFAI
        do m=1,2
          do i=lon0,lon1
            wkm1(i,io2,m) = (delta(io2,m)*qk(i,io1,io1)-
     |        (1.-delta(io2,m))*qk(i,io2,m))/wks1(i)
            wkm1(i,io1,m) = (delta(io1,m)*qk(i,io2,io2)-
     |        (1.-delta(io1,m))*qk(i,io1,m))/wks1(i)
          enddo ! i=lon0,lon1
        enddo ! m=1,2
!
! WKV1=FK-PK*Z(K)
        do i=lon0,lon1
          wkv1(i,io2) = fk(i,io2)
          wkv1(i,io1) = fk(i,io1)
        enddo ! i=lon0,lon1
!
! GAMA(K+1)=ALFAI*RK
        do m=1,2
          do i=lon0,lon1
            gama(i,k+1,io2,m) = 0.
            gama(i,k+1,io1,m) = 0.
            wkv1(i,io2) = wkv1(i,io2)-pk(i,io2,m)*zz(i,k,m)
            wkv1(i,io1) = wkv1(i,io1)-pk(i,io1,m)*zz(i,k,m)
          enddo ! i=lon0,lon1
          do kk=1,2
            do i=lon0,lon1
              gama(i,k+1,io2,m) = gama(i,k+1,io2,m)+wkm1(i,io2,kk)*
     |          rk(i,kk,m)
              gama(i,k+1,io1,m) = gama(i,k+1,io1,m)+wkm1(i,io1,kk)*
     |          rk(i,kk,m)
            enddo ! i=lon0,lon1
          enddo ! kk=1,2
        enddo ! m=1,2
!
! Z(K+1)=WKM1*WKV1
        do i=lon0,lon1
          zz(i,k+1,:) = 0.
        enddo ! i=lon0,lon1
        do m=1,2
          do i=lon0,lon1
            zz(i,k+1,io2) = zz(i,k+1,io2)+wkm1(i,io2,m)*wkv1(i,m)
            zz(i,k+1,io1) = zz(i,k+1,io1)+wkm1(i,io1,m)*wkv1(i,m)
          enddo ! i=lon0,lon1
        enddo ! m=1,2

!       do m=1,2
!         write(6,"('comp: m=',i2,' k=',i2,' lat=',i2)") m,k,lat
!         write(6,"('gama(k+1,io2)=',/,(6e12.4))") gama(:,k+1,io2,m)
!         write(6,"('gama(k+1,io1)=',/,(6e12.4))") gama(:,k+1,io1,m)
!       enddo
!       write(6,"('zz(k+1,io2)=',/,(6e12.4))") zz(:,k+1,io2)
!       write(6,"('zz(k+1,io1)=',/,(6e12.4))") zz(:,k+1,io1)

!
! End main pressure loop:
      enddo levloop ! k=lev0,lev1-1
!
! Save diagnostics:
!     real,dimension(lev0:lev1,lon0:lon1,2) :: zz_ki
!     real,dimension(lev0:lev1,lon0:lon1,2,2) :: gama_ki
!
      do k=lev0,lev1
        zz_ki(k,:,io2) = zz(:,k,io2)
        zz_ki(k,:,io1) = zz(:,k,io1)
        do m=1,2
          gama_ki(k,:,io2,m) = gama(:,k,io2,m)
          gama_ki(k,:,io1,m) = gama(:,k,io1,m)
        enddo
      enddo
!       call addfld('ZZ_O2',' ',' ',zz_ki(k0:k1-1,:,io2),
!    |    'lev',k0,k1-1,'lon',i0,i1,lat)
!       call addfld('ZZ_O1',' ',' ',zz_ki(k0:k1-1,:,io1),
!    |    'lev',k0,k1-1,'lon',i0,i1,lat)
!       call addfld('GAMAO2M1',' ',' ',gama_ki(k0:k1-1,:,io2,1),
!    |    'lev',k0,k1-1,'lon',i0,i1,lat)
!       call addfld('GAMAO2M2',' ',' ',gama_ki(k0:k1-1,:,io2,2),
!    |    'lev',k0,k1-1,'lon',i0,i1,lat)
!       call addfld('GAMAO1M1',' ',' ',gama_ki(k0:k1-1,:,io1,1),
!    |    'lev',k0,k1-1,'lon',i0,i1,lat)
!       call addfld('GAMAO1M2',' ',' ',gama_ki(k0:k1-1,:,io1,2),
!    |    'lev',k0,k1-1,'lon',i0,i1,lat)
!
! Set upper boundary to zero:
      do i=lon0,lon1
        o2_upd(lev1,i,lat) = 0.
        o1_upd(lev1,i,lat) = 0.
        upd(lev1,i,:) = 0.
      enddo ! i=lon0,lon1
!
! Downward sweep:
      do k=lev1-1,lev0,-1
        do isp=io2,io1
          do i=lon0,lon1
            upd(k,i,isp) = zz(i,k+1,isp)
          enddo ! i=lon0,lon1
          do m=1,2
            do i=lon0,lon1
              upd(k,i,isp) = upd(k,i,isp)-gama(i,k+1,isp,m)*
     |                       upd(k+1,i,m)
            enddo
          enddo ! m=1,2
        enddo ! isp=io2,io1
      enddo ! k=lev1-1,lev0,-1
!
! Transfer to output arrays:
      do k=lev0,lev1
        o2_upd(k,lon0:lon1,lat) = upd(k,:,io2)
        o1_upd(k,lon0:lon1,lat) = upd(k,:,io1)
      enddo ! k=lev0,lev1
!
! Upper boundary:
! kp is carried forward from the last iteration of levloop above.
      do i=lon0,lon1
        o2_upd(lev1,i,lat) = 
     |    (1.+.5*ep(i,io2,kp)*dz)/
     |    (1.-.5*ep(i,io2,kp)*dz)*o2_upd(lev1-1,i,lat)
        o1_upd(lev1,i,lat) = 
     |    (1.+.5*ep(i,io1,kp)*dz)/
     |    (1.-.5*ep(i,io1,kp)*dz)*o1_upd(lev1-1,i,lat)
      enddo ! i=lon0,lon1
!
!       call addfld('O2_SOLV',' ',' ',o2_upd(:,lon0:lon1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!       call addfld('O1_SOLV',' ',' ',o1_upd(:,lon0:lon1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!
! End latitude scan:
      enddo ! lat=lat0,lat1
!
! Filter the new composition species:
!
! Fourier smoothing of O2 and O:
      call filter_o2o(o2_upd,lev0,lev1,lon0,lon1,lat0,lat1)
      call filter_o2o(o1_upd,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Resume latitude scan:
      do lat=lat0,lat1
!       call addfld('O2_FILT',' ',' ',o2_upd(:,lon0:lon1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!       call addfld('O1_FILT',' ',' ',o1_upd(:,lon0:lon1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!
! Time smoothing:
        do i=lon0,lon1 
          do k=lev0,lev1
            o2nm_upd(k,i,lat) = dtsmooth_div2*(o2_nm(k,i,lat)+
     |        o2_upd(k,i,lat)) + dtsmooth*o2(k,i,lat)
            o1nm_upd(k,i,lat) = dtsmooth_div2*(o1_nm(k,i,lat)+
     |        o1_upd(k,i,lat)) + dtsmooth*o1(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1 

!       call addfld('O2NM_OUT',' ',' ',o2nm_upd(:,lon0:lon1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!       call addfld('O1NM_OUT',' ',' ',o1nm_upd(:,lon0:lon1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)

      enddo ! lat=lat0,lat1
#ifdef MPI
!
! Periodic points:
!     call mp_periodic_f3d(o2_upd(:,lon0:lon1,lat0-1:lat1+1),
!    |  lev0,lev1,lon0,lon1,lat0-1,lat1+1)
!     call mp_periodic_f3d(o1_upd(:,lon0:lon1,lat0-1:lat1+1),
!    |  lev0,lev1,lon0,lon1,lat0-1,lat1+1)
#endif
!
! Insure non-negative O2,O:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            if (o2_upd(k,i,lat) < small) o2_upd(k,i,lat) = small
            if (o1_upd(k,i,lat) < small) o1_upd(k,i,lat) = small
            if (o2nm_upd(k,i,lat) < small) o2nm_upd(k,i,lat) = small
            if (o1nm_upd(k,i,lat) < small) o1nm_upd(k,i,lat) = small
            if (1.-small-o2_upd(k,i,lat)-o1_upd(k,i,lat) < 0.) then
              o2_upd(k,i,lat) = o2_upd(k,i,lat)*((1.-small)/
     |                         (o2_upd(k,i,lat)+o1_upd(k,i,lat)))
              o1_upd(k,i,lat) = o1_upd(k,i,lat)*((1.-small)/
     |                         (o2_upd(k,i,lat)+o1_upd(k,i,lat)))
            endif
            if (1.-small-o2nm_upd(k,i,lat)-o1nm_upd(k,i,lat) < 0.) then
              o2nm_upd(k,i,lat) = o2nm_upd(k,i,lat)*((1.-small)/
     |                           (o2nm_upd(k,i,lat)+o1nm_upd(k,i,lat)))
              o1nm_upd(k,i,lat) = o1nm_upd(k,i,lat)*((1.-small)/
     |                           (o2nm_upd(k,i,lat)+o1nm_upd(k,i,lat)))
            endif
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

!       call addfld('O2_OUT',' ',' ',o2_upd(:,i0:i1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!       call addfld('O1_OUT',' ',' ',o1_upd(:,i0:i1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!       call addfld('O2_NMOUT',' ',' ',o2nm_upd(:,i0:i1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!       call addfld('O1_NMOUT',' ',' ',o1nm_upd(:,i0:i1,lat),
!    |    'lev',k0,k1,'lon',i0,i1,lat)
!
      enddo ! lat=lat0,lat1
      end subroutine comp
!-----------------------------------------------------------------------
      subroutine advecl(o2,o1,un,vn,o2_advec,o1_advec,
     |  lev0,lev1,lon0,lon1,lat0,lat1,lat)
!
! Horizontal advection for O2,O. 
! In previous versions, this was sub advecl (inline.F), called from 
!   k-loop in comp.F. Here it is called from beginning of comp.F at 
!   all latitudes and saved for later use (fk) inside comp.F k-loop.
! O2,o1,un,vn already have i-2,i-1,i+1,i+2, and j-1,j-2,j+1,j+2 for 
!   finite differencing.
!
      use cons_module,only: dlamda_2div3 ,dlamda_1div12, dphi_2div3,
     |  dphi_1div12,re_inv,racs
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  o2,o1,un,vn
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) ::
     |  o2_advec, o1_advec
!
! Local:
      integer :: k,i
!
      do k=lev0,lev1-1
        do i=lon0,lon1
          o2_advec(k,i,lat) = .5*racs(lat)*
     |      (dlamda_2div3*(o2(k,i+1,lat)-o2(k,i-1,lat))*
     |      (un(k,i+1,lat)+un(k,i-1,lat))-
     |      dlamda_1div12*(o2(k,i+2,lat)-o2(k,i-2,lat))*
     |      (un(k,i+2,lat)+un(k,i-2,lat)))+
     |      .5*re_inv*
     |      (dphi_2div3*(o2(k,i,lat+1)-o2(k,i,lat-1))*
     |      (vn(k,i,lat+1)+vn(k,i,lat-1))-
     |      dphi_1div12*(o2(k,i,lat+2)-o2(k,i,lat-2))*
     |      (vn(k,i,lat+2)+vn(k,i,lat-2)))

          o1_advec(k,i,lat) = .5*racs(lat)*
     |      (dlamda_2div3*(o1(k,i+1,lat)-o1(k,i-1,lat))*
     |      (un(k,i+1,lat)+un(k,i-1,lat))-
     |      dlamda_1div12*(o1(k,i+2,lat)-o1(k,i-2,lat))*
     |      (un(k,i+2,lat)+un(k,i-2,lat)))+
     |      .5*re_inv*
     |      (dphi_2div3*(o1(k,i,lat+1)-o1(k,i,lat-1))*
     |      (vn(k,i,lat+1)+vn(k,i,lat-1))-
     |      dphi_1div12*(o1(k,i,lat+2)-o1(k,i,lat-2))*
     |      (vn(k,i,lat+2)+vn(k,i,lat-2)))
        enddo ! i=lon0,lon1
      enddo ! k=lev0,lev1-1

      end subroutine advecl
!-----------------------------------------------------------------------
      subroutine filter_o2o(fout,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Filter updated W omega:
!
      use params_module,only: nlat,nlonp4,nlon
      use filter_module,only: filter
      use cons_module,only: kut ! kut(nlat)
#ifdef MPI
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi
      implicit none
#else
      implicit none
      integer :: mytidi=0
#endif
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(inout) :: fout(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j,nlevs,nlons,nlats
      real :: fik(nlonp4,lev0:lev1),fkij(lev0:lev1,nlonp4,lat0:lat1)
      real :: fmin,fmax
!
#ifdef VT
!     code = 131 ; state = 'filter_o2o' ; activity='Filtering'
      call vtbegin(131,ier)
#endif
!
      nlevs = lev1-lev0+1
      nlons = lon1-lon0+1
      nlats = lat1-lat0+1
!
! Define lons in w_ki from current task:
      fkij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          fkij(:,i,j) = fout(:,i,j)
        enddo
      enddo ! j=lat0,lat1
!
#ifdef MPI
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0 
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
      if (mytidi==0) then
!
! Define 2d array with all longitudes for filter at each latitude:
        latscan: do j=lat0,lat1
          do i=1,nlonp4
            fik(i,:) = fkij(:,i,j)
          enddo ! i=1,nlonp4
!
! Remove wave numbers > kut(lat):
          call filter(fik,lev0,lev1,kut(j),j) ! tiegcm
!         call filter2(fik,lev0,lev1,j)       ! timegcm
!
! Return filtered array to fkij:
          do i=1,nlonp4
            fkij(:,i,j) = fik(i,:)
          enddo ! i=1,nlonp4
        enddo latscan ! j=lat0,lat1
      endif ! mytidi==0
#ifdef MPI
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
!
! Return filtered array to fout at current task longitudes and latitudes:
      do j=lat0,lat1
        do i=lon0,lon1
          fout(:,i,j) = fkij(:,i,j)
        enddo
      enddo
!
#ifdef VT
!     code = 131 ; state = 'filter_o2o' ; activity='Filtering'
      call vtend(131,ier)
#endif
      end subroutine filter_o2o
