!
      subroutine minor(tn,un,vn,o2,o1,he,mbar,barm,xnmbar,w,
     |  fcomp,fcomp_tm1,fcomp_out,fcomp_tm1_out,
     |  sloss,sprod,flbc,fubc,rmx,phix,alfax,
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
     |  p0,expz,rmass_o2,rmass_o1,rmass_n2,hor,dtr,shapiro,
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
     |  un,       ! zonal wind velocity
     |  vn,       ! meridional wind velocity
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
      real,parameter :: small=1.e-12, tau=1.86e+3, t00=273.
      real :: wi
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  do2dz,    ! do2/dz(k)    (s6)
     |  do1dz,    ! do1/dz(k)    (s7)
     |  dhedz,    ! dhe/dz(k)
     |  pso2,     ! o2           (s8)
     |  pso1,     ! o1           (s9)
     |  pshe,     ! he
     |  dmdz,     ! dm/dz(k)     (s10)
     |  tni,      ! tn (k)       (s5)
     |  s0prod,   ! sprod*mx/nmbar (s15)
     |  ex,ax,    ! (s12,s11)
     |  thdiff,   ! thermal diffusion term (s12)
     |  p_coef,q_coef,r_coef,f_rhs, ! coefficients for tridiagonal solver
     |  molp,molr,molp1,molr1,molq,
     |  eddyppart,eddyrpart,eddyp1part,eddyr1part,
     |  eddyp,eddyr,eddyq,ppart,rpart,qpart,
     |  dpdt,loss_out,moldif,eddydif,veradv
!
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  ftm1_smooth, ! time-smoothed field at time n-1 (s10)
     |  hadvec       ! horizontal advection (output of sub advec) (s13)
!
      real,dimension(lon0:lon1) ::
     |  bmbar,        ! (t6)
     |  bo2,bo1,bhe   ! (t4,t5)
      real,dimension(lat0:lat1) ::
     |  dfactor       ! (t7) (formerly output of sub dfact)
!
! Following diffusion related variables are from comp
      real,dimension(3),parameter ::
     |  ss = (/1.710,1.749,1.718/) ! exponent factor for diff_fac (o2,o,n2)
      real,dimension(3,4),parameter ::
     |  phi = reshape( ! mutual thermal diffusion coefficients among major species
     |    (/0.  ,0.673,0.270,
     |      1.35,0.   ,0.404,
     |      2.16,1.616,0.   ,
     |      1.11,0.769,0.322/),(/3,4/))
      real,dimension(3) ::
     |  diff_fac, ! Correction factor for diffusion coeff betweeen helium and o2,o,n2
     |  salfax,ep
      real,dimension(3,3) :: alfa,invalfa
!
#ifdef VT
!     code = 114 ; state = 'minor' ; activity='ModelCode'
      call vtbegin(114,ier)
#endif
!
!     write(6,"('Enter minor')")

      salfax(1)=phix(1)-phix(4)
      salfax(2)=phix(2)-phix(4)
      salfax(3)=phix(3)-phix(4)

      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2 ! nlonp4==nlonp4
!
! Shapiro time smoother:
!
      call smooth(fcomp_tm1,ftm1_smooth,shapiro,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Horizontal advection (pass k vs i slices at full task subdomain 
! longitudes, and the 5 latitudes centered over the current latitude).
!
      call advec(un,vn,fcomp,hadvec,lev0,lev1,lon0,lon1,lat0,lat1)
!
!----------------------- Begin second latitude scan --------------------
      do lat=lat0,lat1
!       write(6,"('Minor second latitude scan: lat=',i3)") lat

      if (idebug > 0) then
        call addfld('HADVEC',' ',' ',hadvec(lev0:lev1-1,:,lat),
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
        bmbar(i) = 1./
     |    (bo2(i)*rmassinv_o2+
     |     bo1(i)*rmassinv_o1+
     |     bhe(i)*rmassinv_he+
     |     (1.-bo2(i)-bo1(i)-bhe(i))*rmassinv_n2)
      enddo ! i=lon0,lon1
!
! To be completely correct and consistent, need to
!   add pshe and dhedz below (including calculation of ex), to 
!   include helium in diffusion of minor species.
! Ben Foster 2014/06/16 pointed it out, Haonan Wu 2024/03/20 implemented
!
! Lower boundary:
      do i=lon0,lon1
        dmdz(lev0,i) = (mbar(lev0,i,lat)-bmbar(i))/dzp
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
          dmdz(k,i) = (mbar(k,i,lat)-mbar(k-1,i,lat))/dzp
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
          ep(1) = do2dz(k,i)-pso2(k,i)*
     |      (1-(rmass_o2+dmdz(k,i))/barm(k,i,lat))
          ep(2) = do1dz(k,i)-pso1(k,i)*
     |      (1-(rmass_o1+dmdz(k,i))/barm(k,i,lat))
          ep(3) = dhedz(k,i)-pshe(k,i)*
     |      (1-(rmass_he+dmdz(k,i))/barm(k,i,lat))
!
! 1-(m+dm/dz/m) + alpha_x*inv(alpha)*(dps/dz-(1-(mi+dm/dz/m))*ps)
          ex(k,i) = 1.-(rmx+dmdz(k,i))/barm(k,i,lat)
!
          invalfa = matinv3(alfa)
          do n = 1,3
            ex(k,i) = ex(k,i)+dot_product(salfax,invalfa(:,n))*ep(n)
          enddo
!
          ax(k,i) = 
     |      -barm(k,i,lat)*(t00/tni(k,i))**0.25/(tau*rmass_n2*
     |      (phix(4)+salfax(1)*pso2(k,i)+
     |      salfax(2)*pso1(k,i)+salfax(3)*pshe(k,i)))

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
      if (difhor > 0) then
        if (abs(glat(lat)) >= 40.) then
          dfactor(lat) = hor(lat)+1.
        else
          dfactor(lat) = hor(lat)+.5*(1.-cos(abs(glat(lat))*dtr*4.5))
        endif
      else
        dfactor(lat) = 1.
      endif
!
! Set up coefficients for tridiagonal solver:
! (W (f4s(w)) is use-associated from fields module)
! p_coef (s1), q_coef (s2), r_coef (s3), and f_rhs (s4)
! (difk = eddy diffusion)
!
      do i=lon0,lon1
!
! molecular diffusion coefficients
        do k=lev0,lev1
          molp (k,i) = ax(k,i)*(1./dzp+.5*thdiff(k,i))
          molr1(k,i) = ax(k,i)*(1./dzp-.5*thdiff(k,i))
        enddo
!
        do k=lev0,lev1-1
          molp1(k,i) = molp (k+1,i)
          molr (k,i) = molr1(k+1,i)
        enddo
        molp1(lev1,i) = 2*molp (lev1,i)-molp (lev1-1,i)
        molr (lev1,i) = 2*molr1(lev1,i)-molr1(lev1-1,i)
!
        do k=lev0,lev1
          molq(k,i) = molp1(k,i)+molr1(k,i)
        enddo
!
! eddy diffusion coefficients
        do k=lev0,lev1
          eddyppart (k,i) = difk(k,iday)*
     |      (1./dzp-.5*dmdz(k,i)/barm(k,i,lat))
          eddyr1part(k,i) = difk(k,iday)*
     |      (1./dzp+.5*dmdz(k,i)/barm(k,i,lat))
        enddo
!
        do k=lev0,lev1-1
          eddyp1part(k,i) = eddyppart (k+1,i)
          eddyrpart (k,i) = eddyr1part(k+1,i)
        enddo
        eddyp1part(lev1,i) = 2*eddyppart (lev1,i)-eddyppart (lev1-1,i)
        eddyrpart (lev1,i) = 2*eddyr1part(lev1,i)-eddyr1part(lev1-1,i)
!
        do k=lev0,lev1
          eddyp(k,i) = dfactor(lat)*expzmid_inv*eddyppart(k,i)
          eddyr(k,i) = dfactor(lat)*expzmid    *eddyrpart(k,i)
          eddyq(k,i) = dfactor(lat)*
     |      (expzmid*eddyp1part(k,i)+expzmid_inv*eddyr1part(k,i))
        enddo
!
        do k=lev0,lev1-1
          wi = 0.5*(w(k,i,lat)+w(k+1,i,lat))
          ppart(k,i) = expz(k)*(eddyp(k,i)+0.5*wi)/dzp
          rpart(k,i) = expz(k)*(eddyr(k,i)-0.5*wi)/dzp
          qpart(k,i) = expz(k)*(eddyq(k,i)/dzp+dtx2inv-sloss(k,i,lat))
        enddo
!
! vertical advection
        do k=lev0,lev1-1
! s1
          p_coef(k,i) = molp(k,i)/dzp-ppart(k,i)
! s3
          r_coef(k,i) = molr(k,i)/dzp-rpart(k,i)
! s2
          q_coef(k,i) = -molq(k,i)/dzp+qpart(k,i)
! s4
          f_rhs(k,i) = expz(k)*(ftm1_smooth(k,i,lat)*
     |      dtx2inv-hadvec(k,i,lat)+s0prod(k,i))
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
!     subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1,lonmax)
!
      call trsolv(p_coef,q_coef,r_coef,f_rhs,fcomp_out(:,lon0:lon1,lat),
     |  lev0,lev1,lev0,lev1,lon0,lon1,nlonp4)

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
          moldif(k,i) =
     |      (molp(k,i)*fcomp_out(k-1,i,lat)-
     |       molq(k,i)*fcomp_out(k  ,i,lat)+
     |       molr(k,i)*fcomp_out(k+1,i,lat))/dzp/expz(k)
          eddydif(k,i) =
     |      (eddyp(k,i)*fcomp_out(k-1,i,lat)-
     |       eddyq(k,i)*fcomp_out(k  ,i,lat)+
     |       eddyr(k,i)*fcomp_out(k+1,i,lat))/dzp
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
     |  -moldif(lev0+1:lev1-1,:),
     |  'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
      call addfld(name//'_EDDYDIF',name//' eddy diffusion','s-1',
     |  eddydif(lev0+1:lev1-1,:),
     |  'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
      call addfld(name//'_VERADV',name//' vertical advection','s-1',
     |  -veradv(lev0+1:lev1-1,:),
     |  'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
      call addfld(name//'_HORADV',name//' horizontal advection','s-1',
     |  -hadvec(lev0:lev1-1,:,lat),
     |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
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
