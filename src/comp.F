!
      subroutine comp(tn,o2,o2_nm,o1,o1_nm,he,he_nm,w,mbar,xnmbar,
     |                hdo2,hdo1,hdhe,o2_upd,o2nm_upd,o1_upd,o1nm_upd,
     |                he_upd,henm_upd,lev0,lev1,lon0,lon1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the
! Open Source Academic Research License Agreement contained in the file
! tiegcmlicense.txt.
!
! Advance major species O2, O and He.
!
      use params_module,only: dz,glat,zpmid
      use init_module,only: iday
      use cons_module,only: pi,hor,dtr,rmassinv_o2,rmassinv_o1,
     |  rmassinv_he,rmassinv_n2,rmass_o2,rmass_o1,rmass_he,expz,expzmid,
     |  expzmid_inv,difk,dtx2inv,dtsmooth,dtsmooth_div2,difhor,
     |  gask,grav,re,p0
      use chemrates_module,only: fs ! from sub comp_o2o
      use lbc,only: b,fb,pshelb
      use fields_module,only: tlbc
      use addfld_module,only: addfld
      use input_module,only: calc_helium,enforce_n2
      use matutil_module,only: matinv3
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,   ! neutral temperature
     |  o2,   ! O2 (mmr) at current timestep
     |  o1,   ! O  (mmr) at current timestep
     |  he,   ! He (mmr) at current timestep
     |  o2_nm,! O2 (mmr) at time n-1
     |  o1_nm,! O  (mmr) at time n-1
     |  he_nm,! He (mmr) at time n-1
     |  w,    ! vertical velocity at current timestep
     |  mbar, ! mean molecular weight
     |  xnmbar,
     |  hdo2, ! O2 horizontal diffusion (hdif3)
     |  hdo1, ! O  horizontal diffusion (hdif3)
     |  hdhe  ! He horizontal diffusion (hdif3)
!
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out)::
     |  o2_upd,  o1_upd,  he_upd,  ! output: O2,O,He updated for next timestep
     |  o2nm_upd,o1nm_upd,henm_upd ! output: O2,O,He updated for previous timestep
!
! Local:
      integer,dimension(3,3),parameter ::
     |  delta = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      real,parameter ::
     |  tau = 1.86e3, t00 = 273., small = 1.e-9,
     |  thdiffalpha = -0.38 ! Thermal diff alpha coefficient for Helium
      real,dimension(3),parameter ::
     |  ss = (/1.710,1.749,1.718/) ! exponent factor for diff_fac (O2,O,N2)
      real,dimension(3,4),parameter ::
     |  phi = reshape(
     |    (/0.  ,0.673,0.270,
     |      1.35,0.   ,0.404,
     |      2.16,1.616,0.   ,
     |      1.11,0.769,0.322/),(/3,4/))
      integer :: k,i,lat,m,n,ktop
      real :: phi34,zmtop,n2top,n2top1,dn2dz,normalize
      real,dimension(3) :: wkv1
      real,dimension(3,3) :: invalpha,wkm1
      real,dimension(lat0:lat1) :: dfactor
      real,dimension(lon0:lon1,lat0:lat1) ::
     |  bo2,bo1,bhe,bn2,bmbar,o1_ub,he_ub,
     |  flx_he,flx_arg,flx00 ! Helium Mass Flux at upper boundary
      real,dimension(3,lon0:lon1,lat0:lat1) :: epep
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  o2_advec,   o1_advec,   he_advec,    ! horizontal advection
     |  o2nm_smooth,o1nm_smooth,henm_smooth, ! smoothed at time n-1
     |  o2i,o1i,hei,barm,dmdz,tni,dtdz,wi,wks1,
     |  eddyppart,eddyrpart,eddyp1part,eddyr1part,
     |  eddyp,eddyr,eddyq,ppart,rpart,qpart,n2,n2nm
      real,dimension(3,lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  ep,fk,zz,upd,dpdt,prod,loss,moldif,eddydif,veradv,
     |  diff_fac ! Correction factor for diffusion coeff betweeen helium and O2,O,N2
      real,dimension(3,3,lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  alpha,ak,molp,molr,molq,molp1,molr1,pk,rk,qk,expzfs,gama
!
      if (difhor > 0) then
        do lat=lat0,lat1
          if (abs(glat(lat)) >= 40.) then
            dfactor(lat) = hor(lat)+1.
          else
            dfactor(lat) = hor(lat)+.5*(1.-cos(abs(glat(lat))*dtr*4.5))
          endif
        enddo
      else
        dfactor = 1.
      endif
!
      do lat=lat0,lat1
        do i=lon0,lon1
          bo2(i,lat) = fb(1)+
     |      b(1,1)*o2(lev0,i,lat)+
     |      b(1,2)*o1(lev0,i,lat)+
     |      b(1,3)*he(lev0,i,lat)
          bo1(i,lat) = fb(2)+
     |      b(2,1)*o2(lev0,i,lat)+
     |      b(2,2)*o1(lev0,i,lat)+
     |      b(2,3)*he(lev0,i,lat)
          bhe(i,lat) = fb(3)+
     |      b(3,1)*o2(lev0,i,lat)+
     |      b(3,2)*o1(lev0,i,lat)+
     |      b(3,3)*he(lev0,i,lat)
          bn2(i,lat) = max(1.-bo2(i,lat)-bo1(i,lat)-bhe(i,lat),0.)
!
          bmbar(i,lat) = 1./
     |      (bo2(i,lat)*rmassinv_o2+
     |       bo1(i,lat)*rmassinv_o1+
     |       bhe(i,lat)*rmassinv_he+
     |       bn2(i,lat)*rmassinv_n2)
!
          do k=lev0+1,lev1
            o2i(k,i,lat) = 0.5*(o2(k-1,i,lat)+o2(k,i,lat))
            o1i(k,i,lat) = 0.5*(o1(k-1,i,lat)+o1(k,i,lat))
            hei(k,i,lat) = 0.5*(he(k-1,i,lat)+he(k,i,lat))
            barm(k,i,lat) = 0.5*(mbar(k-1,i,lat)+mbar(k,i,lat))
            dmdz(k,i,lat) = (mbar(k,i,lat)-mbar(k-1,i,lat))/dz
          enddo
          o2i(lev0,i,lat) = 0.5*(bo2(i,lat)+o2(lev0,i,lat))
          o1i(lev0,i,lat) = 0.5*(bo1(i,lat)+o1(lev0,i,lat))
          hei(lev0,i,lat) = 0.5*(bhe(i,lat)+he(lev0,i,lat))
          barm(lev0,i,lat) = 0.5*(bmbar(i,lat)+mbar(lev0,i,lat))
          dmdz(lev0,i,lat) = (mbar(lev0,i,lat)-bmbar(i,lat))/dz
!
          do k=lev0+1,lev1-1
            tni(k,i,lat) = 0.5*(tn(k-1,i,lat)+tn(k,i,lat))
            dtdz(k,i,lat) = (tn(k,i,lat)-tn(k-1,i,lat))/dz
          enddo
          tni(lev0,i,lat) = tlbc(i,lat)
          tni(lev1,i,lat) = tn(lev1-1,i,lat)
          dtdz(lev0,i,lat) = 2*(tn(lev0,i,lat)-tlbc(i,lat))/dz
          dtdz(lev1,i,lat) = dtdz(lev1-1,i,lat)
!
          do k=lev0,lev1-1
            wi(k,i,lat) = (w(k,i,lat)+w(k+1,i,lat))/2
          enddo
          wi(lev1,i,lat) = 1.5*w(lev1,i,lat)-0.5*w(lev1-1,i,lat)
        enddo
      enddo
!
! Calculate and save Non Escaping Helium mass flux:
!   Hodges [1973], Eq. 9, with P2 set to (1+tn/3330)
!   similar to Hodges and Johnson [1968], Eq. 46.
!   added by EKS
! 4*p*(kT/mg)^3/2*barm*(1+T/3330)*He/(R^2*(2*pi*g)^1/2*m)
      do lat=lat0,lat1
        do i=lon0,lon1
          flx_arg(i,lat) = -4*p0*expz(lev1-1)*expzmid*
     |      sqrt((gask*tni(lev1,i,lat)/(rmass_he*grav))**3)*
     |      barm(lev1,i,lat)*(1+tni(lev1,i,lat)/3330)*hei(lev1,i,lat)/
     |      (re**2*sqrt(2*pi*grav)*rmass_he)
        enddo
      enddo
!
! Zero out helium flux if requested:
      if (calc_helium == 1) then
        call laplacian(flx_arg,flx_he,lon0,lon1,lat0,lat1)
!
        call addfld('FLX_ARG','Argument of Laplacian','g/cm^2/s',
     |    flx_arg,'lon',lon0,lon1,'lat',lat0,lat1,0)
        call addfld('FLX_HE','UB Mass Flux of Helium','g/cm^2/s',
     |    flx_he,'lon',lon0,lon1,'lat',lat0,lat1,0)
!
      else
        flx_he = 0.
      endif
!
! Calculate and save horizontal advection in o2_advec, o1_advec, he_advec:
      do lat=lat0,lat1
        call advec(o2(:,:,lat-2:lat+2),o2_advec(:,:,lat),
     |    lev0,lev1,lon0,lon1,lat)
        call advec(o1(:,:,lat-2:lat+2),o1_advec(:,:,lat),
     |    lev0,lev1,lon0,lon1,lat)
        call advec(he(:,:,lat-2:lat+2),he_advec(:,:,lat),
     |    lev0,lev1,lon0,lon1,lat)
      enddo
!
! Save smoothed O2, O and He at time n-1:
      call smooth(o2_nm,o2nm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
      call smooth(o1_nm,o1nm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
      call smooth(he_nm,henm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
!
! Correction factors for mutual diffusion of helium with O2,O,N2:
! (use tn at interfaces)
! Haonan Wu 2024/03: correct TN lower boundary to t_lbc (not (tn(lev0)+tlbc)/2)
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            do n=1,3 ! O2,O,N2
              diff_fac(n,k,i,lat) = (tni(k,i,lat)/t00)**(1.75-ss(n))
            enddo
          enddo
        enddo
      enddo
!
! WKS1 = MBAR/M4*(T00/T)**0.25/TAU
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            wks1(k,i,lat) = barm(k,i,lat)*rmassinv_n2*
     |        (t00/tni(k,i,lat))**0.25/tau
          enddo
        enddo
      enddo
!
! EP = 1-(M+DMBAR/DZ)/MBAR
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            ep(1,k,i,lat) = 1.-(rmass_o2+dmdz(k,i,lat))/barm(k,i,lat)
            ep(2,k,i,lat) = 1.-(rmass_o1+dmdz(k,i,lat))/barm(k,i,lat)
            ep(3,k,i,lat) = 1.-(rmass_he+dmdz(k,i,lat))/barm(k,i,lat)-
     |        thdiffalpha*dtdz(k,i,lat)/tni(k,i,lat)
          enddo
        enddo
      enddo
!
! ALFA
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            phi34 = diff_fac(3,k,i,lat)*phi(3,4)
!
            alpha(1,1,k,i,lat) = -phi(1,4)-
     |        (phi(1,2)-phi(1,4))*o1i(k,i,lat)-
     |        (diff_fac(1,k,i,lat)*phi(1,3)-phi(1,4))*hei(k,i,lat)
            alpha(2,2,k,i,lat) = -phi(2,4)-
     |        (phi(2,1)-phi(2,4))*o2i(k,i,lat)-
     |        (diff_fac(2,k,i,lat)*phi(2,3)-phi(2,4))*hei(k,i,lat)
            alpha(3,3,k,i,lat) = -phi34-
     |        (diff_fac(1,k,i,lat)*phi(3,1)-phi34)*o2i(k,i,lat)-
     |        (diff_fac(2,k,i,lat)*phi(3,2)-phi34)*o1i(k,i,lat)
            alpha(1,2,k,i,lat) = (phi(1,2)-phi(1,4))*o2i(k,i,lat)
            alpha(1,3,k,i,lat) =
     |        (diff_fac(1,k,i,lat)*phi(1,3)-phi(1,4))*o2i(k,i,lat)
            alpha(2,1,k,i,lat) = (phi(2,1)-phi(2,4))*o1i(k,i,lat)
            alpha(2,3,k,i,lat) =
     |        (diff_fac(2,k,i,lat)*phi(2,3)-phi(2,4))*o1i(k,i,lat)
            alpha(3,1,k,i,lat) =
     |        (diff_fac(1,k,i,lat)*phi(3,1)-phi34)*hei(k,i,lat)
            alpha(3,2,k,i,lat) =
     |        (diff_fac(2,k,i,lat)*phi(3,2)-phi34)*hei(k,i,lat)
          enddo
        enddo
      enddo
!
! AK = MBAR/M4*(T00/T)**0.25/TAU*INV(ALFA)
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            invalpha = matinv3(alpha(:,:,k,i,lat))
            do n=1,3
              do m=1,3
                ak(m,n,k,i,lat) = invalpha(m,n)*wks1(k,i,lat)
              enddo
            enddo
          enddo
        enddo
      enddo
!
! Calculate molecular diffusion coefficients of O2,O1,He
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            do n=1,3
              do m=1,3
                molp (m,n,k,i,lat) =
     |            ak(m,n,k,i,lat)*(1./dz+ep(n,k,i,lat)/2.)
                molr1(m,n,k,i,lat) =
     |            ak(m,n,k,i,lat)*(1./dz-ep(n,k,i,lat)/2.)
              enddo
            enddo
          enddo
        enddo
      enddo
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            do n=1,3
              do m=1,3
                molp1(m,n,k,i,lat) = molp (m,n,k+1,i,lat)
                molr (m,n,k,i,lat) = molr1(m,n,k+1,i,lat)
              enddo
            enddo
          enddo
          do n=1,3
            do m=1,3
              molp1(m,n,lev1,i,lat) =
     |          2*molp (m,n,lev1,i,lat)-molp (m,n,lev1-1,i,lat)
              molr (m,n,lev1,i,lat) =
     |          2*molr1(m,n,lev1,i,lat)-molr1(m,n,lev1-1,i,lat)
            enddo
          enddo
        enddo
      enddo
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            do n=1,3
              do m=1,3
                molq(m,n,k,i,lat) =
     |            molp1(m,n,k,i,lat)+molr1(m,n,k,i,lat)
              enddo
            enddo
          enddo
        enddo
      enddo
!
! Calculate eddy diffusion coefficients
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            eddyppart (k,i,lat) = difk(k,iday)*
     |        (1./dz-dmdz(k,i,lat)/(2*barm(k,i,lat)))
            eddyr1part(k,i,lat) = difk(k,iday)*
     |        (1./dz+dmdz(k,i,lat)/(2*barm(k,i,lat)))
          enddo
        enddo
      enddo
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            eddyp1part(k,i,lat) = eddyppart (k+1,i,lat)
            eddyrpart (k,i,lat) = eddyr1part(k+1,i,lat)
          enddo
          eddyp1part(lev1,i,lat) =
     |      2*eddyppart (lev1,i,lat)-eddyppart (lev1-1,i,lat)
          eddyrpart (lev1,i,lat) =
     |      2*eddyr1part(lev1,i,lat)-eddyr1part(lev1-1,i,lat)
        enddo
      enddo
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            eddyp(k,i,lat) = dfactor(lat)*
     |        expzmid_inv*eddyppart(k,i,lat)
            eddyr(k,i,lat) = dfactor(lat)*
     |        expzmid    *eddyrpart(k,i,lat)
            eddyq(k,i,lat) = dfactor(lat)*
     |        (expzmid    *eddyp1part(k,i,lat)+
     |         expzmid_inv*eddyr1part(k,i,lat))
          enddo
        enddo
      enddo
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            ppart(k,i,lat) = expz(k)*(eddyp(k,i,lat)+0.5*wi(k,i,lat))/dz
            rpart(k,i,lat) = expz(k)*(eddyr(k,i,lat)-0.5*wi(k,i,lat))/dz
            qpart(k,i,lat) = expz(k)*(eddyq(k,i,lat)/dz+dtx2inv)
            do n=1,3
              do m=1,3
                expzfs(m,n,k,i,lat) = expz(k)*fs(i,k,m,n,lat)
              enddo
            enddo
          enddo
        enddo
      enddo
!
! GENERATE PK, QK, RK
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            do n=1,3
              do m=1,3
                pk(m,n,k,i,lat) =
     |            molp(m,n,k,i,lat)/dz-delta(m,n)*ppart(k,i,lat)
                rk(m,n,k,i,lat) =
     |            molr(m,n,k,i,lat)/dz-delta(m,n)*rpart(k,i,lat)
                qk(m,n,k,i,lat) =
     |            -molq(m,n,k,i,lat)/dz+delta(m,n)*qpart(k,i,lat)-
     |            expzfs(m,n,k,i,lat)
              enddo
            enddo
          enddo
        enddo
      enddo
!
! Use advection saved from advecl calls at beginning of routine:
! Add explicit source terms to fk:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            fk(1,k,i,lat) = expz(k)*(o2nm_smooth(k,i,lat)*dtx2inv-
     |         o2_advec(k,i,lat)+fs(i,k,1,0,lat)+hdo2(k,i,lat))
            fk(2,k,i,lat) = expz(k)*(o1nm_smooth(k,i,lat)*dtx2inv-
     |         o1_advec(k,i,lat)+fs(i,k,2,0,lat)+hdo1(k,i,lat))
            fk(3,k,i,lat) = expz(k)*(henm_smooth(k,i,lat)*dtx2inv-
     |         he_advec(k,i,lat)+fs(i,k,3,0,lat)+hdhe(k,i,lat))
          enddo
        enddo
      enddo
!
! Lower boundary:
      do lat=lat0,lat1
        do i=lon0,lon1
          qk(:,:,lev0,i,lat) = qk(:,:,lev0,i,lat)+
     |      matmul(pk(:,:,lev0,i,lat),b)
          fk(:,lev0,i,lat) = fk(:,lev0,i,lat)-
     |      matmul(pk(:,:,lev0,i,lat),fb)
          do n=1,3
            do m=1,3
              pk(m,n,lev0,i,lat) = 0.
            enddo
          enddo
        enddo
      enddo
!
! Upper boundary:
      do lat=lat0,lat1
        do i=lon0,lon1
          do n=1,3
            epep(n,i,lat) =
     |        (1.+.5*ep(n,lev1,i,lat)*dz)/(1.-.5*ep(n,lev1,i,lat)*dz)
          enddo
          do n=1,3
            do m=1,3
              qk(m,n,lev1-1,i,lat) = qk(m,n,lev1-1,i,lat)+
     |          epep(n,i,lat)*rk(m,n,lev1-1,i,lat)
            enddo
          enddo
!
! Separate longitude loop to calculate helium lateral exospheric transport flux UBC (Added by EKS)
! Add lateral exospheric transport mass flux
          flx00(i,lat) = wks1(lev1,i,lat)*p0/grav
          o1_ub(i,lat) = (alpha(2,3,lev1,i,lat)-alpha(2,2,lev1,i,lat))*
     |      flx_he(i,lat)/(flx00(i,lat)*(1./dz-0.5*ep(2,lev1,i,lat)))
          he_ub(i,lat) = (alpha(3,3,lev1,i,lat)-alpha(3,2,lev1,i,lat))*
     |      flx_he(i,lat)/(flx00(i,lat)*(1./dz-0.5*ep(3,lev1,i,lat)))
          do n=1,3
            fk(n,lev1-1,i,lat) = fk(n,lev1-1,i,lat)-
     |        rk(n,2,lev1-1,i,lat)*o1_ub(i,lat)-
     |        rk(n,3,lev1-1,i,lat)*he_ub(i,lat)
          enddo
          do n=1,3
            do m=1,3
              rk(m,n,lev1-1,i,lat) = 0.
            enddo
          enddo
        enddo
      enddo
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do n=1,3
            zz(n,lev0,i,lat) = 0.
          enddo
          do n=1,3
            do m=1,3
              gama(m,n,lev0,i,lat) = 0.
            enddo
          enddo
!
          do k=lev0,lev1-1
!
! ALFA = QK-PK*GAMA(K-1)
! WKM1 = INV(ALFA)
! ALFA refers to the block diagonal matrices,
!   and GAMA to the upper block diagonal matrices
!   in the Thomas algorithm solution
!   to the block tridiagonal system of equations
            wkm1 = matinv3(qk(:,:,k,i,lat)-
     |        matmul(pk(:,:,k,i,lat),gama(:,:,k,i,lat)))
!
! WKV1 = FK-PK*Z(K)
            wkv1 = fk(:,k,i,lat)-matmul(pk(:,:,k,i,lat),zz(:,k,i,lat))
!
! GAMA(K+1) = INV(ALFA)*RK
            gama(:,:,k+1,i,lat) = matmul(wkm1,rk(:,:,k,i,lat))
!
! Z(K+1) = WKM1*WKV1
            zz(:,k+1,i,lat) = matmul(wkm1,wkv1)
          enddo
        enddo
      enddo
!
      do lat=lat0,lat1
        do i=lon0,lon1
!
! Set upper boundary to zero:
          do n=1,3
            upd(n,lev1,i,lat) = 0.
          enddo
!
! Downward sweep:
          do k=lev1-1,lev0,-1
            upd(:,k,i,lat) = zz(:,k+1,i,lat)-
     |        matmul(gama(:,:,k+1,i,lat),upd(:,k+1,i,lat))
          enddo
        enddo
      enddo
!
! Upper boundary:
      do lat=lat0,lat1
        do i=lon0,lon1
          do n=1,3
            upd(n,lev1,i,lat) = epep(n,i,lat)*upd(n,lev1-1,i,lat)
          enddo
          upd(2,lev1,i,lat) = upd(2,lev1,i,lat)+o1_ub(i,lat)
          upd(3,lev1,i,lat) = upd(3,lev1,i,lat)+he_ub(i,lat)
        enddo
      enddo
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            dpdt(1,k,i,lat) =
     |        dtx2inv*(upd(1,k,i,lat)-o2nm_smooth(k,i,lat))
            dpdt(2,k,i,lat) =
     |        dtx2inv*(upd(2,k,i,lat)-o1nm_smooth(k,i,lat))
            dpdt(3,k,i,lat) =
     |        dtx2inv*(upd(3,k,i,lat)-henm_smooth(k,i,lat))
          enddo
!
          do k=lev0,lev1-1
            do n=1,3
              prod(n,k,i,lat) = fs(i,k,n,0,lat)
            enddo
            loss(:,k,i,lat) = matmul(fs(i,k,:,1:3,lat),upd(:,k,i,lat))
          enddo
!
          do k=lev0+1,lev1-1
            moldif(:,k,i,lat) =
     |        (matmul(molp(:,:,k,i,lat),upd(:,k-1,i,lat))-
     |         matmul(molq(:,:,k,i,lat),upd(:,k  ,i,lat))+
     |         matmul(molr(:,:,k,i,lat),upd(:,k+1,i,lat)))/dz/expz(k)
            do n=1,3
              eddydif(n,k,i,lat) =
     |          (eddyp(k,i,lat)*upd(n,k-1,i,lat)-
     |           eddyq(k,i,lat)*upd(n,k  ,i,lat)+
     |           eddyr(k,i,lat)*upd(n,k+1,i,lat))/dz
              veradv(n,k,i,lat) = wi(k,i,lat)/(2*dz)*
     |          (upd(n,k+1,i,lat)-upd(n,k-1,i,lat))
            enddo
          enddo
        enddo
      enddo
!
      do lat=lat0,lat1
        call addfld('DO2DT','O2 changing rate','s-1',
     |    dpdt(1,:,:,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DO1DT','O1 changing rate','s-1',
     |    dpdt(2,:,:,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DHEDT','He changing rate','s-1',
     |    dpdt(3,:,:,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('O2_PROD','O2 production','s-1',
     |    prod(1,lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O1_PROD','O1 production','s-1',
     |    prod(2,lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('HE_PROD','He production','s-1',
     |    prod(3,lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O2_LOSS','O2 loss','s-1',
     |    loss(1,lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O1_LOSS','O1 loss','s-1',
     |    loss(2,lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('HE_LOSS','He loss','s-1',
     |    loss(3,lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O2_HORDIF','O2 horizontal diffusion','s-1',
     |    hdo2(lev0:lev1-1,lon0:lon1,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O1_HORDIF','O1 horizontal diffusion','s-1',
     |    hdo1(lev0:lev1-1,lon0:lon1,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('HE_HORDIF','He horizontal diffusion','s-1',
     |    hdhe(lev0:lev1-1,lon0:lon1,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O2_MOLDIF','O2 molecular diffusion','s-1',
     |    -moldif(1,lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O1_MOLDIF','O1 molecular diffusion','s-1',
     |    -moldif(2,lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('HE_MOLDIF','He molecular diffusion','s-1',
     |    -moldif(3,lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O2_EDDYDIF','O2 eddy diffusion','s-1',
     |    eddydif(1,lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O1_EDDYDIF','O1 eddy diffusion','s-1',
     |    eddydif(2,lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('HE_EDDYDIF','He eddy diffusion','s-1',
     |    eddydif(3,lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O2_VERADV','O2 vertical advection','s-1',
     |    -veradv(1,lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O1_VERADV','O1 vertical advection','s-1',
     |    -veradv(2,lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('HE_VERADV','He vertical advection','s-1',
     |    -veradv(3,lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O2_HORADV','O2 horizontal advection','s-1',
     |    -o2_advec(lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('O1_HORADV','O1 horizontal advection','s-1',
     |    -o1_advec(lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('HE_HORADV','He horizontal advection','s-1',
     |    -he_advec(lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo
!
! Transfer to output arrays:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            o2_upd(k,i,lat) = upd(1,k,i,lat)
            o1_upd(k,i,lat) = upd(2,k,i,lat)
            he_upd(k,i,lat) = upd(3,k,i,lat)
          enddo
        enddo
!       call addfld('O2_SOLV',' ',' ',o2_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('O1_SOLV',' ',' ',o1_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('HE_SOLV',' ',' ',he_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo
!
      if (calc_helium == 0) then
        he_upd = pshelb
        henm_upd = pshelb
      endif
!
! Filter the new composition species:
! Fourier smoothing of O2, O and He:
      call filter_comp(o2_upd,lev0,lev1,lon0,lon1,lat0,lat1)
      call filter_comp(o1_upd,lev0,lev1,lon0,lon1,lat0,lat1)
      call filter_comp(he_upd,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Time smoothing:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            o2nm_upd(k,i,lat) = dtsmooth*o2(k,i,lat)+
     |        dtsmooth_div2*(o2_nm(k,i,lat)+o2_upd(k,i,lat))
            o1nm_upd(k,i,lat) = dtsmooth*o1(k,i,lat)+
     |        dtsmooth_div2*(o1_nm(k,i,lat)+o1_upd(k,i,lat))
            henm_upd(k,i,lat) = dtsmooth*he(k,i,lat)+
     |        dtsmooth_div2*(he_nm(k,i,lat)+he_upd(k,i,lat))
          enddo
        enddo
      enddo
!
! Ensure non-negative O2,O,He,N2:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            if (o2_upd(k,i,lat) < small) o2_upd(k,i,lat) = small
            if (o1_upd(k,i,lat) < small) o1_upd(k,i,lat) = small
            if (he_upd(k,i,lat) < small) he_upd(k,i,lat) = small
!
            if (o2nm_upd(k,i,lat) < small) o2nm_upd(k,i,lat) = small
            if (o1nm_upd(k,i,lat) < small) o1nm_upd(k,i,lat) = small
            if (henm_upd(k,i,lat) < small) henm_upd(k,i,lat) = small
!
            n2(k,i,lat) =
     |        1.-o2_upd(k,i,lat)-o1_upd(k,i,lat)-he_upd(k,i,lat)
            n2nm(k,i,lat) =
     |        1.-o2nm_upd(k,i,lat)-o1nm_upd(k,i,lat)-henm_upd(k,i,lat)
!
            if (n2(k,i,lat) < small) n2(k,i,lat) = small
            if (n2nm(k,i,lat) < small) n2nm(k,i,lat) = small
          enddo
        enddo
      enddo
!
! When O2+O+He~=1 (typically around z=9), N2 will not scale properly with altitudes
! This is an ad-hoc correction to let N2 have the right scaling at very high altitudes
      if (enforce_n2) then
        do lat=lat0,lat1
          do i=lon0,lon1
!
! Find the level where N2 is no longer a major species
            do ktop=lev0,lev1
              if (n2(ktop,i,lat) < 0.01) exit ! Typically around z=7
            enddo
!
! Use extrapolation above this height to replace the original values
! The scaling uses number density instead of mass mixing ratio
            if (ktop>1 .and. ktop<lev1) then
              zmtop = zpmid(ktop)
              n2top = log(n2(ktop,i,lat)*xnmbar(ktop,i,lat)*rmassinv_n2)
              n2top1 =
     |          log(n2(ktop-1,i,lat)*xnmbar(ktop-1,i,lat)*rmassinv_n2)
              dn2dz = (n2top-n2top1)/dz
              do k=ktop+1,lev1
                n2(k,i,lat) = exp(n2top+dn2dz*(zpmid(k)-zmtop))/
     |            xnmbar(k,i,lat)/rmassinv_n2
              enddo
            endif
          enddo
        enddo
      endif
!
! EKS comment: divisor should be the same
! when normalizing O2, O and He to ensure non-negativity
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            normalize = o2_upd(k,i,lat)+o1_upd(k,i,lat)+
     |        he_upd(k,i,lat)+n2(k,i,lat)
            o2_upd(k,i,lat) = o2_upd(k,i,lat)/normalize
            o1_upd(k,i,lat) = o1_upd(k,i,lat)/normalize
            he_upd(k,i,lat) = he_upd(k,i,lat)/normalize
!
            normalize = o2nm_upd(k,i,lat)+o1nm_upd(k,i,lat)+
     |        henm_upd(k,i,lat)+n2nm(k,i,lat)
            o2nm_upd(k,i,lat) = o2nm_upd(k,i,lat)/normalize
            o1nm_upd(k,i,lat) = o1nm_upd(k,i,lat)/normalize
            henm_upd(k,i,lat) = henm_upd(k,i,lat)/normalize
          enddo
        enddo
      enddo
!
      end subroutine comp
!-----------------------------------------------------------------------
      subroutine filter_comp(fout,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Filter updated O2,O:
!
      use params_module,only: nlonp4
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi
      implicit none
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
      integer :: i,j,nlevs,nlons,nlats,ier
      real :: fik(nlonp4,lev0:lev1),fkij(lev0:lev1,nlonp4,lat0:lat1)
!
#ifdef VT
!     code = 131 ; state = 'filter_comp' ; activity='Filtering'
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
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1,'O2O')
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
          call ringfilter(fik,lev0,lev1,1,'O2O',j)
!
! Return filtered array to fkij:
          do i=1,nlonp4
            fkij(:,i,j) = fik(i,:)
          enddo ! i=1,nlonp4
        enddo latscan ! j=lat0,lat1
      endif ! mytidi==0
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1,
     |  'O2O')
!
! Return filtered array to fout at current task longitudes and latitudes:
      do j=lat0,lat1
        do i=lon0,lon1
          fout(:,i,j) = fkij(:,i,j)
        enddo
      enddo
!
#ifdef VT
!     code = 131 ; state = 'filter_comp' ; activity='Filtering'
      call vtend(131,ier)
#endif
      end subroutine filter_comp
