!
      subroutine minor(tn,o2,o1,he,mbar,barm,xnmbar,w,fcomp,fcomp_tm1,
     |  fcomp_out,fcomp_tm1_out,sloss,sprod,flbc,fubc,rmx,phix,alfax,
     |  lev0,lev1,lon0,lon1,lat0,lat1,idebug,name)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Advance minor species fcomp by one time step.
! This is called from comp_no and comp_n4s.
!
      use params_module,only: nlonp4,glat
      use lbc,only: b,fb
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,dzp,
     |  p0,expz,rmass_o2,rmass_o1,rmass_n2,hor,dtr,pi,shapiro,
     |  difk,expzmid_inv,expzmid,dtx2inv,grav,avo,dtsmooth,
     |  dtsmooth_div2,difhor,rmass_he,rmassinv_he
      use init_module,only: iday
      use fields_module,only: tlbc
      use matutil_module,only: matinv3
      use addfld_module,only: addfld
      use mpi_module,only: mp_bndlons_f3d, mp_periodic_f3d,
     |  mp_periodic_f2d
!
! VT means vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,idebug
      character(len=*),intent(in) :: name
!
! Input fields are at full task subdomain:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,       ! neutral temperature (deg K)
     |  o2,       ! molecular oxygen (mmr)
     |  o1,       ! atomic oxygen (mmr)
     |  he,       ! helium (mmr)
     |  mbar,     ! mean molecular weight
     |  barm,     ! mean molecular weight
     |  xnmbar,   ! p0*e(-z)*mbar/kT
     |  w,        ! vertical wind
     |  fcomp,    ! input species (for finite diffs in advec)
     |  fcomp_tm1 ! input species at time n-1 
!
! flbc(:,1)=A, flbc(:,2)=B, flbc(:,3)=C define lower boundary condition, 
!   where: A*DPSX/DZ + B*PSX + C = 0.
! Boundary conditions are allocated at task subdomains, without ghost cells.
!
      real,dimension(lon0:lon1,3,lat0:lat1),intent(in) :: flbc
!
! fubc = diffusive upward number flux at upper boundary.
      real,dimension(lon0:lon1,lat0:lat1),intent(in) :: fubc
!
      real,intent(in) :: 
     |  phix(4), ! diffusion vector
     |  rmx,     ! molecular weight of fcomp minor species
     |  alfax    ! thermal diffusion coefficient
!
! Input production and loss:
! sloss: sx/n(x) where sx is portion of number density source
!   proportional to n(x), the minor species number density.
! sprod: s0, portion of number density source independent of n(x)
!
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) ::
     |  sloss, sprod
!
! Output minor species (allocated at full subdomains):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  fcomp_out,     ! output species 
     |  fcomp_tm1_out  ! output species at time n-1
!
! Local:
      integer :: k,i,lonbeg,lonend,lat,n,ier
!
      integer :: nk,nkm1,nlevs
      real,parameter :: small=1.e-12, tau=1.86e+3, t00=273.
      real :: rlat,wi
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  hadvec,   ! horizontal advection (output of sub advec) (s13)
     |  do2dz,    ! do2/dz(k)    (s6)
     |  do1dz,    ! do1/dz(k)    (s7)
     |  dhedz,    ! dhe/dz(k)
     |  pso2,     ! o2           (s8)
     |  pso1,     ! o1           (s9)
     |  pshe,     ! he
     |  dmdz,     ! dm/dz(k)     (s10)
     |  xmbar_k,  ! barm (k)     (s11)
     |  xmbar_kh, ! mbar (k)     (s12)
     |  tni,      ! tn (k)       (s5)
     |  s0prod,   ! sprod*mx/nmbar (s15)
     |  ex,ax,    ! (s12,s11)
     |  thdiff,   ! thermal diffusion term (s12)
     |  p_coef,q_coef,r_coef,f_rhs, ! coefficients for tridiagonal solver
     |  molp,molq,molr,eddyp,eddyq,eddyr,moldif,eddydif,veradv,
     |  loss_out,dpdt
!
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  ftm1_smooth  ! time-smoothed field at time n-1 (s10)
!
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  ftm1_jdif ! lat diffs (includes lon boundaries) (s8)
!
      integer :: jm2,jm1,j0,jp1,jp2 ! lat-2,lat-1,lat,lat+1,lat+2
!
      real,dimension(lon0:lon1) :: 
     |  xmbari,       ! (t6)
     |  bo2,bo1,bhe,  ! (t4,t5)
     |  dfactor       ! (t7) (formerly output of sub dfact)
!
! Following diffusion related variables are from comp
      real :: phi(3,4) ! mutual thermal diffusion coefficients among major species
      real,dimension(3) ::
     |  ss,       ! exponent factor for diff_fac
     |  diff_fac, ! Correction factor for diffusion coeff betweeen helium and o2,o,n2
     |  salfax,salfaxinvalfa
      real,dimension(3,3) :: alfa,invalfa
!
#ifdef VT
!     code = 114 ; state = 'minor' ; activity='ModelCode'
      call vtbegin(114,ier)
#endif
!
!     write(6,"('Enter minor')")

      phi(:,1)=(/0.  ,0.673,0.270/)
      phi(:,2)=(/1.35,0.   ,0.404/)
      phi(:,3)=(/2.16,1.616,0.   /)
      phi(:,4)=(/1.11,0.769,0.322/)

      ss = (/1.710, 1.749, 1.718/) ! for o2,o,n2

      salfax(1)=phix(1)-phix(4)
      salfax(2)=phix(2)-phix(4)
      salfax(3)=phix(3)-phix(4)

      nk = lev1-lev0+1
      nkm1 = nk-1
      nlevs = nk

      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2 ! nlonp4==nlonp4
!
! Latitude scan over task subdomain:
      do lat=lat0,lat1
!       write(6,"('Minor first latitude scan: lat=',i3)") lat
        jm2 = lat-2
        jm1 = lat-1
        j0  = lat
        jp1 = lat+1
        jp2 = lat+2
!
! Shapiro time smoother:
!
! ftm1_jdif = 4th order diffs in latitude at time n-1:
! f = f(j)-shapiro * ( f(j+2)+f(j-2) - 4.*(f(j+1)+f(j-1)) + 6.*f(j) )
!
        do i=lon0,lon1
          do k=lev0,lev1
            ftm1_jdif(k,i,lat) = fcomp_tm1(k,i,j0) - shapiro * 
     |           (fcomp_tm1(k,i,jp2) + fcomp_tm1(k,i,jm2)  - 
     |        4.*(fcomp_tm1(k,i,jp1) + fcomp_tm1(k,i,jm1)) + 
     |        6.*fcomp_tm1(k,i,j0))
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

      enddo ! lat=lat0,lat1
!------------------------- End first latitude scan ---------------------
!
! Get subdomain boundary longitudes for ftm1_jdif, to complete shapiro
! time smoother. (mp_periodic_f3d call opened up 8/31/11)
!
      call mp_periodic_f3d(ftm1_jdif(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
      call mp_bndlons_f3d(ftm1_jdif,nlevs,lon0,lon1,lat0,lat1,1,0)
!----------------------- Begin second latitude scan --------------------
      do lat=lat0,lat1
!       write(6,"('Minor second latitude scan: lat=',i3)") lat
!
! ftm1_smooth = zonally-smoothed field at time n-1:
! f = f(i)-shapiro * ( f(i+2)+f(i-2) - 4.*(f(i+1)+f(i-1)) + 6.*f(i) )
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          ftm1_smooth(k,i,lat) = ftm1_jdif(k,i,lat) - shapiro *
     |         (ftm1_jdif(k,i+2,lat) + ftm1_jdif(k,i-2,lat)  -
     |      4.*(ftm1_jdif(k,i+1,lat) + ftm1_jdif(k,i-1,lat)) +
     |      6.* ftm1_jdif(k,i,lat))
        enddo ! k=lev0,lev1-1
      enddo ! i=lonbeg,lonend
!
! Set periodic points to zero:
      if (lon0==1) ftm1_smooth(:,lon0:lon0+1,lat) = 0.
      if (lon1==nlonp4) ftm1_smooth(:,lon1-1:lon1,lat) = 0.
!
! Horizontal advection (pass k vs i slices at full task subdomain 
! longitudes, and the 5 latitudes centered over the current latitude).
!
      call advec(fcomp(:,:,lat-2:lat+2),hadvec,lev0,lev1,lon0,lon1,lat)
!
! Periodic points for advection:
!
! Oct, 2010: Comment this call to prevent divergence in N4S and NO
!   between runs w/ different numbers of processors (see also divrg.F):
!     call mp_periodic_f2d(hadvec,lon0,lon1,lat,lat)
!
! Set periodic points to zero:
!     if (lon0==1) hadvec(:,lon0:lon0+1) = 0.
!     if (lon1==nlonp4) hadvec(:,lon1-1:lon1) = 0.

      if (idebug > 0) then
        call addfld('HADVEC',' ',' ',hadvec(lev0:lev1-1,:),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('SLOSS' ,' ',' ',sloss(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('SPROD' ,' ',' ',sprod(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      endif
!
      do i=lon0,lon1
        bo2(i) = b(1,1)*o2(lev0,i,lat)+b(1,2)*o1(lev0,i,lat)+
     |    b(1,3)*he(lev0,i,lat)+fb(1)
        bo1(i) = b(2,1)*o2(lev0,i,lat)+b(2,2)*o1(lev0,i,lat)+
     |    b(2,3)*he(lev0,i,lat)+fb(2)
        bhe(i) = b(3,1)*o2(lev0,i,lat)+b(3,2)*o1(lev0,i,lat)+
     |    b(3,3)*he(lev0,i,lat)+fb(3)
        xmbari(i) = 1./
     |    (bo2(i)*rmassinv_o2+
     |     bo1(i)*rmassinv_o1+
     |     bhe(i)*rmassinv_he+
     |     (1.-bo2(i)-bo1(i)-bhe(i))*rmassinv_n2)
      enddo ! i=lon0,lon1
!
! xmbar_kh = mbar, xmbar_k = barm:
      do i=lon0,lon1
        do k=lev0,lev1
          xmbar_kh(k,i) = mbar(k,i,lat)
          xmbar_k(k,i) = barm(k,i,lat)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
! To be completely correct and consistent, need to
!   add pshe and dhedz below (including calculation of ex), to 
!   include helium in diffusion of minor species.
! Ben Foster 2014/06/16 pointed it out, Haonan Wu 2024/03/20 implemented
!
! Lower boundary:
      do i=lon0,lon1
        dmdz(lev0,i) = (xmbar_kh(lev0,i)-xmbari(i))/dzp
        pso1(lev0,i) = .5*(bo1(i)+o1(lev0,i,lat))  
        pso2(lev0,i) = .5*(bo2(i)+o2(lev0,i,lat))  
        pshe(lev0,i) = .5*(bhe(i)+he(lev0,i,lat))
        do1dz(lev0,i) = (o1(lev0,i,lat)-bo1(i))/dzp
        do2dz(lev0,i) = (o2(lev0,i,lat)-bo2(i))/dzp
        dhedz(lev0,i) = (he(lev0,i,lat)-bhe(i))/dzp
      enddo ! i=lon0,lon1
!
! Levels 2 -> lev1:
      do i=lon0,lon1
        do k=lev0+1,lev1
          dmdz(k,i) = (xmbar_kh(k,i)-xmbar_kh(k-1,i))/dzp
          pso1(k,i) = .5*(o1(k,i,lat)+o1(k-1,i,lat))
          pso2(k,i) = .5*(o2(k,i,lat)+o2(k-1,i,lat))
          pshe(k,i) = .5*(he(k,i,lat)+he(k-1,i,lat))
          do1dz(k,i) = (o1(k,i,lat)-o1(k-1,i,lat))/dzp
          do2dz(k,i) = (o2(k,i,lat)-o2(k-1,i,lat))/dzp
          dhedz(k,i) = (he(k,i,lat)-he(k-1,i,lat))/dzp
        enddo ! k=lev0+1,lev1
      enddo ! i=lon0,lon1
!
! tni = tn at interfaces:
      do i=lon0,lon1
        tni(lev0,i) = tlbc(i,lat) 
        tni(lev1,i) = tn(lev1-1,i,lat)
        do k=lev0+1,lev1-1
          tni(k,i) = .5*(tn(k,i,lat)+tn(k-1,i,lat))
        enddo ! k=lev0+1,lev1-1
      enddo ! i=lon0,lon1

      if (idebug > 0) then
        call addfld('XMBAR_KH',' ',' ',xmbar_kh,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('XMBAR_K' ,' ',' ',xmbar_k ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DMDZ0'   ,' ',' ',dmdz    ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('PSO1'    ,' ',' ',pso1    ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('PSO2'    ,' ',' ',pso2    ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('PSHE'    ,' ',' ',pshe    ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DO1DZ'   ,' ',' ',do1dz   ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DO2DZ'   ,' ',' ',do2dz   ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DHEDZ'   ,' ',' ',dhedz   ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('TNI'     ,' ',' ',tni     ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      endif
 
      do i=lon0,lon1
        do k=lev0,lev1-1
          s0prod(k,i) = sprod(k,i,lat)*rmx/xnmbar(k,i,lat)  ! s15
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

      do i=lon0,lon1
        do k=lev0,lev1
!
! Correction factors for mutual diffusion of helium with o2,o,n2
          do n = 1,3
            diff_fac(n) = (tni(k,i)/t00)**(1.75-ss(n))
          enddo
!
! alpha matrix (also see comp)
          alfa(1,1) = -(phi(1,4)+
     |      (phi(1,2)-phi(1,4))*pso1(k,i)+
     |      (diff_fac(1)*phi(1,3)-phi(1,4))*pshe(k,i))
          alfa(2,2) = -(phi(2,4)+
     |      (phi(2,1)-phi(2,4))*pso2(k,i)+
     |      (diff_fac(2)*phi(2,3)-phi(2,4))*pshe(k,i))
          alfa(3,3) = -(diff_fac(3)*phi(3,4)+
     |      (diff_fac(1)*phi(3,1)-diff_fac(3)*phi(3,4))*pso2(k,i)+
     |      (diff_fac(2)*phi(3,2)-diff_fac(3)*phi(3,4))*pso1(k,i))
          alfa(1,2) = (phi(1,2)-phi(1,4))*pso2(k,i)
          alfa(1,3) = (diff_fac(1)*phi(1,3)-phi(1,4))*pso2(k,i)
          alfa(2,1) = (phi(2,1)-phi(2,4))*pso1(k,i)
          alfa(2,3) = (diff_fac(2)*phi(2,3)-phi(2,4))*pso1(k,i)
          alfa(3,1) =
     |      (diff_fac(1)*phi(3,1)-diff_fac(3)*phi(3,4))*pshe(k,i)
          alfa(3,2) =
     |      (diff_fac(2)*phi(3,2)-diff_fac(3)*phi(3,4))*pshe(k,i)
!
! alpha_x * inv(alpha)
          invalfa = matinv3(alfa)
          do n = 1,3
            salfaxinvalfa(n) = dot_product(salfax,invalfa(:,n))
          enddo
!
! 1-(m+dm/dz/m) + alpha_x*inv(alpha)*(dps/dz-(1-(mi+dm/dz/m))*ps)
          ex(k,i) = 1.-(rmx+dmdz(k,i))/xmbar_k(k,i)+
     |      salfaxinvalfa(1)*(do2dz(k,i)-
     |      (1.-(rmass_o2+dmdz(k,i))/xmbar_k(k,i))*pso2(k,i))+
     |      salfaxinvalfa(2)*(do1dz(k,i)-
     |      (1.-(rmass_o1+dmdz(k,i))/xmbar_k(k,i))*pso1(k,i))+
     |      salfaxinvalfa(3)*(dhedz(k,i)-
     |      (1.-(rmass_he+dmdz(k,i))/xmbar_k(k,i))*pshe(k,i))
!
          dmdz(k,i) = dmdz(k,i)/xmbar_k(k,i)          ! s10
!
          ax(k,i) = 
     |      -xmbar_k(k,i)/(tau*rmass_n2)*(t00/tni(k,i))**0.25/
     |      (phix(4)+salfax(1)*pso2(k,i)+
     |      salfax(2)*pso1(k,i)+salfax(3)*pshe(k,i))

        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

      if (idebug > 0) then
        call addfld('S0PROD',' ',' ',s0prod,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('EX_MNR'    ,' ',' ',ex,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DMDZ1' ,' ',' ',dmdz  ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('AX0'   ,' ',' ',ax    ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      endif
!
! thdiff = EX-ALFAX*D/DS(LN(T(TOT))  (Thermal diffusion term) (s12)
      do i=lon0,lon1
        do k=lev0+1,lev1-1
          thdiff(k,i) = ex(k,i)-alfax*(tni(k+1,i)-tni(k-1,i))/
     |      (2.*dzp*tni(k,i))
        enddo ! k=lev0+1,lev1-1
        thdiff(lev0,i) = ex(lev0,i)-alfax*(tni(lev0+1,i)-tni(lev0,i))/
     |    (dzp*tni(lev0,i))
        thdiff(lev1,i) = ex(lev1,i)-alfax*(tni(lev1,i)-tni(lev1-1,i))/
     |    (dzp*tni(lev1,i))
      enddo ! i=lon0,lon1

      if (idebug > 0) 
     |  call addfld('THDIFF0',' ',' ',thdiff,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! difhor flag is a parameter flag set to 1 in cons.F
! dfactor = 
!   .5*(1.+SIN(PI*(ABS(RLATM)-PI/6.)/(PI/3.))) FOR ABS(RLATM).LT.PI/3.
! dfactor = 1. FOR  ABS(RLATM).GE.PI/3
! (dfactor was in sub dfact in earlier versions)
!
      if (difhor > 0) then
        rlat = glat(lat)*dtr
        if (abs(rlat)-pi/4.5 >= 0.) then
          dfactor(:) = hor(lat)+1.
        else
          dfactor(:) = hor(lat)+.5*(1.+sin(pi*(abs(rlat)-pi/9.)/
     |      (pi/4.5)))
        endif
      else
        dfactor(:) = 1.
      endif
!
! Set up coefficients for tridiagonal solver:
! (W (f4s(w)) is use-associated from fields module)
! p_coef (s1), q_coef (s2), r_coef (s3), and f_rhs (s4)
! (difk = eddy diffusion)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! molecular diffusion coefficients
          molp(k,i) = ax(k  ,i)*(1./dzp+.5*thdiff(k  ,i))
          molr(k,i) = ax(k+1,i)*(1./dzp-.5*thdiff(k+1,i))
          molq(k,i) =
     |      ax(k  ,i)*(1./dzp-.5*thdiff(k  ,i))+
     |      ax(k+1,i)*(1./dzp+.5*thdiff(k+1,i))
!
! eddy diffusion coefficients
          eddyp(k,i) = expzmid_inv*difk(k,iday)*(1./dzp-.5*dmdz(k,i))
          eddyr(k,i) = expzmid*difk(k+1,iday)*(1./dzp+.5*dmdz(k+1,i))
          eddyq(k,i) =
     |      expzmid_inv*difk(k,iday)*(1./dzp+.5*dmdz(k,i))+
     |      expzmid*difk(k+1,iday)*(1./dzp-.5*dmdz(k+1,i))
!
! vertical advection
          wi = 0.5*(w(k,i,lat)+w(k+1,i,lat))
! s1
          p_coef(k,i) = molp(k,i)/dzp-expz(k)*
     |      (eddyp(k,i)*dfactor(i)+0.5*wi)/dzp
! s3
          r_coef(k,i) = molr(k,i)/dzp-expz(k)*
     |      (eddyr(k,i)*dfactor(i)-0.5*wi)/dzp
! s2
          q_coef(k,i) = -molq(k,i)/dzp+expz(k)*
     |      (eddyq(k,i)*dfactor(i)/dzp+dtx2inv-sloss(k,i,lat))
! s4
          f_rhs(k,i) = expz(k)*(ftm1_smooth(k,i,lat)*
     |      dtx2inv-hadvec(k,i)+s0prod(k,i))
        enddo ! k=lev0,lev1-1
!
! Lower boundary (use flbc(:,1-3):
        q_coef(lev0,i) = q_coef(lev0,i)+p_coef(lev0,i)*(flbc(i,1,lat)+
     |    .5*flbc(i,2,lat)*dzp)/(flbc(i,1,lat)-.5*flbc(i,2,lat)*dzp)
        f_rhs(lev0,i) = f_rhs(lev0,i)-p_coef(lev0,i)*flbc(i,3,lat)*dzp/
     |    (flbc(i,1,lat)-.5*flbc(i,2,lat)*dzp)
        p_coef(lev0,i) = 0.
!
! Upper boundary (use fubc(:):
        p_coef(lev1,i) = 1.+.5*dzp*thdiff(lev1,i)
        q_coef(lev1,i) = p_coef(lev1,i)-2.
        r_coef(lev1,i) = 0.
        f_rhs (lev1,i) = -grav*rmx*fubc(i,lat)*dzp/(p0*ax(lev1,i)*avo)
      enddo ! i=lon0,lon1

      if (idebug > 0) then
        call addfld('AX1'    ,' ',' ',ax         ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DMDZ2'  ,' ',' ',dmdz       ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('THDIFF1',' ',' ',thdiff     ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('FSMOOTH',' ',' ',
     |    ftm1_smooth(lev0:lev1-1,lon0:lon1,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('W_OMEGA',' ',' ',w(lev0:lev1,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('P_COEF' ,' ',' ',p_coef(lev0:lev1-1,:),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('R_COEF' ,' ',' ',r_coef(lev0:lev1-1,:),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('Q_COEF' ,' ',' ',q_coef(lev0:lev1-1,:),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('F_RHS'  ,' ',' ',f_rhs      ,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      endif
!
! Solve tridiagonal system (fcomp_out is allocated at full subdomain):
!
!     subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1,lonmax,lat,
!    |  idebug)
!
      call trsolv(p_coef,q_coef,r_coef,f_rhs,fcomp_out(:,lon0:lon1,lat),
     |  lev0,lev1,lev0,lev1,lon0,lon1,nlonp4,lat,0)

      if (idebug > 0)
     |  call addfld('MNR_SOLV'  ,' ',' ',fcomp_out(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1
          dpdt(k,i) = dtx2inv*
     |      (fcomp_out(k,i,lat)-ftm1_smooth(k,i,lat))
        enddo
        do k=lev0,lev1-1
          loss_out(k,i) = sloss(k,i,lat)*fcomp_out(k,i,lat)
        enddo
        do k=lev0+1,lev1-1
          wi = 0.5*(w(k,i,lat)+w(k+1,i,lat))
          moldif(k,i) = -1/expz(k)/dzp*
     |      (molp(k,i)*fcomp_out(k-1,i,lat)-
     |       molq(k,i)*fcomp_out(k  ,i,lat)+
     |       molr(k,i)*fcomp_out(k+1,i,lat))
          eddydif(k,i) = dfactor(i)/dzp*
     |      (eddyp(k,i)*fcomp_out(k-1,i,lat)-
     |       eddyq(k,i)*fcomp_out(k  ,i,lat)+
     |       eddyr(k,i)*fcomp_out(k+1,i,lat))
          veradv(k,i) = wi/(2*dzp)*
     |      (fcomp_out(k+1,i,lat)-fcomp_out(k-1,i,lat))
        enddo
      enddo
!
      call addfld('D'//name//'DT',name//' changing rate','s-1',
     |  dpdt,'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld(name//'_PROD',name//' production','s-1',
     |  s0prod(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      call addfld(name//'_LOSS',name//' loss','s-1',
     |  loss_out(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      call addfld(name//'_MOLDIF',name//' molecular diffusion','s-1',
     |  moldif(lev0+1:lev1-1,:),
     |  'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
      call addfld(name//'_EDDYDIF',name//' eddy diffusion','s-1',
     |  eddydif(lev0+1:lev1-1,:),
     |  'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
      call addfld(name//'_VERADV',name//' vertical advection','s-1',
     |  -veradv(lev0+1:lev1-1,:),
     |  'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
      call addfld(name//'_HORADV',name//' horizontal advection','s-1',
     |  -hadvec(lev0:lev1-1,:),'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      enddo ! lat=lat0,lat1
!----------------------- End second latitude scan ----------------------
!
! Wave filter the minor species output. Filter_minor does 3d gather/scatter
! for fft, so it is isolated from the latitude loop. Pass task subdomain,
! not including ghost cells.
!
      call filter_minor(fcomp_out(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,name)
!
!----------------------- Begin third latitude scan --------------------
      do lat=lat0,lat1
!       write(6,"('Minor third latitude scan: lat=',i3)") lat

        if (idebug > 0)
     |    call addfld('MNR_FILT'  ,' ',' ',fcomp_out(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Time smoothing:
        do i=lonbeg,lonend
          do k=lev0,lev1
            fcomp_tm1_out(k,i,lat) = dtsmooth*fcomp(k,i,lat)+
     |        dtsmooth_div2*(fcomp_tm1(k,i,lat)+fcomp_out(k,i,lat))
          enddo ! k=lev0,lev1
        enddo ! i=lonbeg,lonend

        if (idebug > 0)
     |    call addfld('MNR_SMOO'  ,' ',' ',
     |      fcomp_tm1_out(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      enddo ! lat=lat0,lat1
!----------------------- End third latitude scan ----------------------
!
! Periodic points for fcomp_out and fcomp_tm1_out:
      call mp_periodic_f3d(fcomp_out(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
      call mp_periodic_f3d(fcomp_tm1_out(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
!
!----------------------- Begin fourth latitude scan --------------------
      do lat=lat0,lat1
!       write(6,"('Minor fourth latitude scan: lat=',i3)") lat
!
! Insure density > 0.
! There is a much more elaborate check in earlier versions, which
! is not used here. This results in very small "diamond diffs" in
! N4S, and no diffs at all in NO (wrt tgcm15).
!
      do i=lon0,lon1
        do k=lev0,lev1
          if (fcomp_out(k,i,lat) < small) then
!           write(6,"('minor: fcomp_out < small: lat=',i2,' i=',
!    |        i2,' k=',i2,' fcomp_out=',e12.4)") 
!    |        lat,i,k,fcomp_out(k,i,lat)
            fcomp_out(k,i,lat) = small
          endif
          if (fcomp_tm1_out(k,i,lat) < small) then
!           write(6,"('minor: fcomp_tm1_out < small: lat=',i2,' i=',
!    |        i2,' k=',i2,' fcomp_out=',e12.4)") 
!    |        lat,i,k,fcomp_tm1_out(k,i,lat)
            fcomp_tm1_out(k,i,lat) = small
          endif
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

      if (idebug > 0) then
        call addfld('MNR_OUT'  ,' ',' ',fcomp_out(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('MNR_TM1'  ,' ',' ',fcomp_tm1_out(:,lon0:lon1,lat)
     |    ,'lev',lev0,lev1,'lon',lon0,lon1,lat)
      endif
!
! End fourth and final latitude loop:
      enddo ! lat=lat0,lat1
#ifdef VT
!     code = 114 ; state = 'minor' ; activity='ModelCode'
      call vtend(114,ier)
#endif
!
      end subroutine minor
!-----------------------------------------------------------------------
      subroutine filter_minor(fout,lev0,lev1,lon0,lon1,lat0,lat1,name)
!
! Filter minor species (called from minor after trsolv)
!
      use params_module,only: nlonp4
!
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(inout) :: fout(lev0:lev1,lon0:lon1,lat0:lat1)
      character(len=*),intent(in) :: name
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j,ier
      real :: f_ik(nlonp4,lev0:lev1),f_kij(lev0:lev1,nlonp4,lat0:lat1)
!
#ifdef VT
!     code = 123 ; state = 'filter_minor' ; activity='Filtering'
      call vtbegin(123,ier)
#endif
!     write(6,"('Enter filter_minor')")
!
! Define lons in f_ki from current task:
      f_kij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          f_kij(:,i,j) = fout(:,i,j)
        enddo
      enddo
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(f_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,name)
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
      if (mytidi==0) then
!
! Loop through subdomain latitudes:
        latscan: do j=lat0,lat1
!
! Define 2d array with all longitudes for filter:
          do i=1,nlonp4
            f_ik(i,:) = f_kij(:,i,j)
          enddo ! i=1,nlonp4
!
          call ringfilter(f_ik,lev0,lev1,3,name,j)
!
! Return filtered array to f_kij:
          do i=1,nlonp4
            f_kij(:,i,j) = f_ik(i,:)
          enddo ! i=1,nlonp4
        enddo latscan ! j=lat0,lat1
      endif ! mytidi==0
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(f_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,
     |  name)
!
! Return filtered array to fout at current task longitudes:
      do j=lat0,lat1
        do i=lon0,lon1
          fout(:,i,j) = f_kij(:,i,j)
        enddo
      enddo ! j=lat0,lat1
!
#ifdef VT
!     code = 123 ; state = 'filter_minor' ; activity='Filtering'
      call vtend(123,ier)
#endif
      end subroutine filter_minor
