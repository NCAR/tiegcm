! 
      subroutine settei(tn,o2,o1,n2,ne,te,ti,op,o2p,nop,barm,qji_ti,
     |  te_out,ti_out,lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate electron and ion temperatures.
!
      use params_module,only: dz,nlonp4,spval
      use cons_module,only: pi,rtd,evergs,p0,expz,expzmid,expzmid_inv,
     |  boltz,rmassinv_o2,rmassinv_o1,rmassinv_n2,gask,grav,dipmin,avo
      use input_module,only: f107
      use chapman_module,only: chi   ! solar zenith angle (nlonp4,nlat)
      use magfield_module,only: rlatm,dipmag
      use aurora_module,only: qteaur ! (nlonp4,nlat)
      use qrj_module,only: ! Q is modified, all others are input.
     |  qtotal,! total heating 
     |  qop2p, ! o+(2p)
     |  qop2d, ! o+(2d)
     |  qo2p,  ! o2+ ionization
     |  qop,   ! o+  ionization
     |  qn2p,  ! n2+ ionization
     |  qnp,   ! n+  ionization
     |  qnop   ! no+ ionization
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  ne,      ! electron density (cm3)
     |  te,      ! electron temperature (from previous time step)
     |  ti,      ! ion temperature (from previous time step)
     |  op,      ! O+
     |  o2p,     ! O2+
     |  nop,     ! NO+
     |  barm,    ! mean molecular weight
     |  qji_ti   ! joule heating from qjoule_ti (used ui,vi)
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  te_out,  ! output electron temperature (deg K) 
     |  ti_out   ! output ion temperature (deg K)
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: k,i,n,ier
      integer :: nk,nkm1
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  te_int,    ! electron temperature (interfaces)
     |  tn_int,    ! neutral temperature (interfaces)
     |  o2n,       ! O2 number density (midpoints or interfaces)
     |  o1n,       ! O1 number density (midpoints or interfaces)
     |  n2n,       ! N2 number density (midpoints or interfaces)
     |  xnmbar,    ! p0*e(-z)*barm/kT  (minpoints or interfaces)
     |  root_te,   ! sqrt(te)
     |  root_tn,   ! sqrt(tn)
     |  root_ne,   ! sqrt(ne)
     |  tek0,      ! ke/te**2.5 (s15)
     |  h_mid,h_int,
     |  p_coef,    ! coefficient for trisolv     (s1)
     |  q_coef,    ! coefficient for trisolv     (s2)
     |  r_coef,    ! coefficient for trisolv     (s3)
     |  rhs,       ! right-hand-side for trisolv (s4)
     |  qtot,      ! total ionization rate       (s11)
     |  fki,       ! work array
     |  qe,        ! source term                 (s10)
     |  q_eni,     ! heating from electron/neutral and electron/ion collisions
     |  coll_en2v, ! electron/N2vib collision    (s9)
!
! Cooling rates (heat loss):
     |  loss_en2v, ! electron/N2vib loss term    (s10)
     |  loss_en2,  ! electron/N2 loss
     |  loss_eo2,  ! electron/O2 loss
     |  loss_eo1d, ! electron/O(1d) loss
     |  loss_eo1,  ! electron/O loss
     |  loss_xen,  ! L0*(E,N) (s8)
     |  loss_en,   ! electrons/neutrals loss     (s11)
     |  loss_ei,   ! electron/ion loss           (s10)
     |  loss_in    ! ion/neutral loss            (s9)
      real,parameter ::
     |  fpolar = -3.0e+9,  ! polar te flux
     |  del    = 1.e-6  ,
!
! Correction factors for neutral heating due to L(E,O1D)
     |  alam = 0.0069   ,
     |  ad   = 0.0091   ,
     |  sd   = 2.3e-11
      real :: 
     |  root2,  ! sqrt(2)
     |  f107te  ! solar flux
!
! a,fed,fen,fe,sindipmag have a z dimension only for diagnostic plotting:
      real,dimension(lon0:lon1) ::
     |  a,fed,fen, ! day/night
     |  fe,        ! heat flux at upper boundary
     |  sindipmag  ! sin(dipmag)
!
! For diagnostic plotting:
      real,dimension(lev0:lev1-1,lon0:lon1) ::
     |  a_ki,        ! for diagnostic plotting of a
     |  fed_ki,      ! for diagnostic plotting of fed
     |  fen_ki,      ! for diagnostic plotting of fen
     |  fe_ki,       ! for diagnostic plotting of fe
     |  dipmag_ki,   ! for diagnostic plotting of dipmag
     |  chi_ki,      ! for diagnostic plotting of chi
     |  qteaur_ki,   ! for diagnostic plotting of qteaur
     |  sindipmag_ki ! for diagnostic plotting of sindipmag
!
      real,dimension(lev0:lev1) :: expzp ! e(-z) (midpoints or interfaces)
!
#ifdef VT
!     code = 126 ; state = 'settei' ; activity='ModelCode'
      call vtbegin(126,ier)
#endif
!
      root2 = sqrt(2.)
      f107te = f107
      if (f107te > 235.) f107te = 235.
      nk = lev1-lev0+1
      nkm1 = nk-1

!     call addfld('QJI_TI',' ',' ',qji_ti(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        if (abs(rlatm(i,lat))-pi/4.5 >= 0.) then
          a(i) = 1.
        else
          a(i) = .5*(1.+sin(pi*(abs(rlatm(i,lat))-pi/9.)/(pi/4.5)))
        endif 
!
! Increased heat flux for TE fom protonosphere.
        fed(i) = ( -5.0e+7*f107te*a(i)-4.0e+7*f107te)*1.2
        fen(i) = fed(i)/2.
        fed(i) = fed(i)+qteaur(i,lat)     ! t4
        fen(i) = fen(i)+qteaur(i,lat)     ! t5
        if (chi(i,lat)-.5*pi >= 0.) then  ! chi==t2
          fe(i) = fen(i)                ! t1
        else
          fe(i) = fed(i)
        endif
        if ((chi(i,lat)*rtd-80.)*(chi(i,lat)*rtd-100.)>=0.) then
          fe(i) = fe(i)*evergs
        else
          fe(i) = (.5*(fed(i)+fen(i))+.5*(fed(i)-fen(i))*
     |      cos(pi*(chi(i,lat)*rtd-80.)/20.))*evergs
        endif
!
! Add fpolar if magnetic latitude >= 60 degrees:
        if (abs(rlatm(i,lat))-pi/3.>=0.) fe(i) = fe(i)+fpolar*evergs
!
! For plotting (first dimension is lev0:lev1-1):
        a_ki  (:,i) = a(i)
        chi_ki(:,i) = chi(i,lat)
        qteaur_ki(:,i) = qteaur(i,lat)
        fed_ki(:,i) = fed(i)
        fen_ki(:,i) = fen(i)
        fe_ki (:,i) = fe(i)
      enddo ! i=lon0,lon1

!     call addfld('MAGLAT',' ',' ',a_ki     ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('CHI'   ,' ',' ',chi_ki   ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('QTEAUR',' ',' ',qteaur_ki,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('FED'   ,' ',' ',fed_ki   ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('FEN'   ,' ',' ',fen_ki   ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('FE'    ,' ',' ',fe_ki    ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

!     write(6,"('settei: lat=',i2,' fed=',/,(6e12.4))") lat,fed
!     write(6,"('settei: lat=',i2,' fen=',/,(6e12.4))") lat,fen
!     write(6,"('settei: lat=',i2,' fe =',/,(6e12.4))") lat,fe

!
! te,o2,o,n2,tn at interfaces: 
      do i=lon0,lon1
        do k=lev0+1,lev1-1
          te_int(k,i) = .5*(te(k,i)+te(k-1,i))
          o2n(k,i)    = .5*(o2(k,i)+o2(k-1,i))
          o1n(k,i)    = .5*(o1(k,i)+o1(k-1,i))
          n2n(k,i)    = .5*(n2(k,i)+n2(k-1,i))
          tn_int(k,i) = .5*(tn(k,i)+tn(k-1,i))
        enddo ! k=lev0+1,lev1-2
!
! Bottom:
        te_int(lev0,i) = 1.5*te(lev0,i)-.5*te(lev0+1,i)
        o2n(lev0,i)    = 1.5*o2(lev0,i)-.5*o2(lev0+1,i)
        o1n(lev0,i)    = 1.5*o1(lev0,i)-.5*o1(lev0+1,i)
        n2n(lev0,i)    = 1.5*n2(lev0,i)-.5*n2(lev0+1,i)
        tn_int(lev0,i) = 1.5*tn(lev0,i)-.5*tn(lev0+1,i)
!
! Top:
        te_int(lev1,i) = 1.5*te(lev1-1,i)-.5*te(lev1-2,i)
        o2n(lev1,i)    = 1.5*o2(lev1-1,i)-.5*o2(lev1-2,i)
        o1n(lev1,i)    = 1.5*o1(lev1-1,i)-.5*o1(lev1-2,i)
        n2n(lev1,i)    = 1.5*n2(lev1-1,i)-.5*n2(lev1-2,i)
        tn_int(lev1,i) = 1.5*tn(lev1-1,i)-.5*tn(lev1-2,i)
!
! N2:
        do k=lev0,lev1
          if (n2n(k,i) < 0.) n2n(k,i) = 0.
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     call addfld('TE_INT'  ,' ',' ',te_int,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O2_INT'  ,' ',' ',o2n   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O1_INT'  ,' ',' ',o1n   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TN_INT'  ,' ',' ',tn_int,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2_INT'  ,' ',' ',n2n   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Convert o2,o,n2 to number density (interfaces):
      do k=lev0,lev1-1
        expzp(k) = expzmid_inv*expz(k)
      enddo ! k=lev0,lev1
      expzp(lev1) = expzmid*expz(lev1-1)
      do i=lon0,lon1
        do k=lev0,lev1
! xnmbar at interfaces:
          xnmbar(k,i) = p0*expzp(k)*barm(k,i)/(boltz*tn_int(k,i)) ! s8
          o2n(k,i) = xnmbar(k,i)*o2n(k,i)*rmassinv_o2       ! s13
          o1n(k,i) = xnmbar(k,i)*o1n(k,i)*rmassinv_o1       ! s12
          n2n(k,i) = xnmbar(k,i)*n2n(k,i)*rmassinv_n2       ! s11
          root_te(k,i) = sqrt(te_int(k,i))
!
          tek0(k,i) = 7.5e5/
     |      (1.+3.22e4*te_int(k,i)**2/
     |      ne(k,i)*(root_te(k,i)*
     |      (2.82e-17 - 3.41e-21   * te_int (k,i))*n2n(k,i)+
     |      (2.20e-16 + 7.92e-18   * root_te(k,i))*o2n(k,i)+
     |       1.10e-16 * (1.+5.7e-4 * te_int (k,i))*o1n(k,i)))*evergs
!
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     call addfld('XNMBARI',' ',' ',xnmbar,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TEK0'   ,' ',' ',tek0  ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      do i=lon0,lon1
        do k=lev0,lev1-1
          h_mid(k,i) = gask*tn(k,i)/(.5*(barm(k,i)+barm(k+1,i))*grav) ! s7
        enddo ! k=lev0,lev1-1
        do k=lev0,lev1
          h_int(k,i) = gask*tn_int(k,i)/(barm(k,i)*grav)              ! s6
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     call addfld('H_MID'   ,' ',' ',h_mid(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('H_INT'   ,' ',' ',h_int,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      do i=lon0,lon1
        if (abs(dipmag(i,lat)) >= dipmin) then
          sindipmag(i) = (sin(dipmag(i,lat)))**2 ! t2,s2
        else
          sindipmag(i) = (sin(dipmin))**2
        endif
        if (sindipmag(i) < .10) sindipmag(i) = .10
!
! Start coefficients and rhs for trsolv:
        do k=lev0,lev1-1
          p_coef(k,i) = 2./7.*sindipmag(i)/(h_mid(k,i)*dz**2) ! s1
          r_coef(k,i) = p_coef(k,i)*tek0(k+1,i)/h_int(k+1,i)  ! s3
          p_coef(k,i) = p_coef(k,i)*tek0(k  ,i)/h_int(k  ,i)  ! s1
          q_coef(k,i) = -(p_coef(k,i)+r_coef(k,i))            ! s2
          rhs(k,i) = 0.                                       ! s4
        enddo ! k=lev0,lev1-1
!
! Bottom boundary:
        q_coef(lev0,i) = q_coef(lev0,i)-p_coef(lev0,i)
        rhs(lev0,i) = rhs(lev0,i)-2.*p_coef(lev0,i)*tn_int(lev0,i)**3.5
        p_coef(lev0,i) = 0.
!
! Upper boundary:
        q_coef(lev1-1,i) = q_coef(lev1-1,i)+r_coef(lev1-1,i)
        rhs(lev1-1,i) = rhs(lev1-1,i)+r_coef(lev1-1,i)*dz*3.5*
     |    h_int(lev1,i)*fe(i)/tek0(lev1,i)
        r_coef(lev1-1,i) = 0.
      enddo ! i=lon0,lon1

      do i=lon0,lon1
        dipmag_ki(:,i) = dipmag(i,lat)
        sindipmag_ki(:,i) = sindipmag(i)
      enddo

!     call addfld('DIPMAG' ,' ',' ',dipmag_ki,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('SINDIPM',' ',' ',sindipmag_ki,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

!     call addfld('P_COEFa' ,' ',' ',p_coef(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEFa' ,' ',' ',q_coef(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEFa' ,' ',' ',r_coef(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('RHS0'   ,' ',' ',rhs   (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

!
! qtot = total ionization rate = sum(Qxx) = 
! (QO2+) + (QO+) + (QN2+) + (QNO+) + (QN+) + (QO+(2D)) + (QO+(2P))
!
!     call addfld('QO2P' ,' ',' ',qo2p (:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP'  ,' ',' ',qop  (:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QN2P' ,' ',' ',qn2p (:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QNOP' ,' ',' ',qnop (:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QNP'  ,' ',' ',qnp  (:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP2D',' ',' ',qop2d(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP2P',' ',' ',qop2p(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Note re comparison with tgcm15: Because the check for qtot < 1.e-20
! is inside the reduction operation loop, tiegcm1 must sum separately
! for each ion species, as in tgcm15. If tiegcm1 sums all species at
! each grid point (qtot(k,i)=qo2p(k,i)+qop(k,i)+...), there are diffs
! at the bottom boundary (which is where qtot < 1.e-20). By summing
! each species separately, as below, there are no diffs. Ions are also
! summed in qjion.F, but there is no check for < 1.e-20 there, so summing
! all species at each grid point works fine (see qtot in qjion.F).
!
      qtot = 0. ! whole array init
      do n=1,7
        select case (n)
          case (1) 
            fki(:,:) = qo2p(:,:,lat)
          case (2) 
            fki(:,:) = qop(:,:,lat)
          case (3) 
            fki(:,:) = qn2p(:,:,lat)
          case (4) 
            fki(:,:) = qnop(:,:,lat)
          case (5) 
            fki(:,:) = qnp(:,:,lat)
          case (6) 
            fki(:,:) = qop2d(:,:,lat)
          case (7) 
            fki(:,:) = qop2p(:,:,lat)
        end select
        do i=lon0,lon1
          do k=lev0,lev1
            qtot(k,i) = qtot(k,i)+fki(k,i)
            if (qtot(k,i) < 1.e-20) qtot(k,i) = 1.e-20
          enddo
        enddo
      enddo ! n=1,7

!     call addfld('QTOT_SUM',' ',' ',qtot,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      do i=lon0,lon1
        do k=lev0,lev1-1
          qtot(k,i) = sqrt(qtot(k,i)*qtot(k+1,i))
        enddo
        qtot(lev1,i) = 0.
      enddo ! i=lon0,lon1
!     call addfld('QTOT',' ',' ',qtot   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Set Ne (midpoints "(K+1/2)"):
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          root_ne(k,i) = ne(k,i)*ne(k+1,i)
          if (root_ne(k,i) < 1.e4) root_ne(k,i) = 1.e4
          root_ne(k,i) = sqrt(root_ne(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('ROOT_NE',' ',' ',root_ne(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Set up o2,o,n2 number densities at midpoints:
!
      do i=lon0,lon1 ! (DO 20)
        do k=lev0,lev1-1
          xnmbar(k,i) = p0*expz(k)*.5*(barm(k,i)+barm(k+1,i))/
     |      (boltz*tn(k,i))
          o2n(k,i) = xnmbar(k,i)*o2(k,i)*rmassinv_o2  ! s14
          o1n(k,i) = xnmbar(k,i)*o1(k,i)*rmassinv_o1  ! s13
          n2n(k,i) = n2(k,i)
          if (n2n(k,i) < 0.) n2n(k,i) = 0.
          n2n(k,i) = xnmbar(k,i)*n2n(k,i)*rmassinv_n2 ! s12
!
! Calculate source term qe (s10)
! Comment from earlier version (maybe the *1.0 below was once *2.0):
!   "Correction facor of 2 increase in TE heating rate"
!
          qe(k,i) = alog(root_ne(k,i)/(o2n(k,i)+n2n(k,i)+0.1*o1n(k,i)))
          qe(k,i) = exp(-((((0.001996*qe(k,i)+0.08034)*qe(k,i)+1.166)*
     |              qe(k,i)+6.941)*qe(k,i)+12.75))*1.0
!
! Subtract qe from right-hand-side:
          rhs(k,i) = rhs(k,i)-qe(k,i)*qtot(k,i)*evergs
          root_te(k,i) = sqrt(te(k,i))
!
! Electron/N2 collision A(E,N2,VIB) (s9):
!
        if (te(k,i) >= 1000.) then
          coll_en2v(k,i) = 2.e-7*exp(-4605.2/te(k,i))
        else
          coll_en2v(k,i) = 5.71e-8*exp(-3352.6/te(k,i))
        endif
        if (te(k,i) > 2000.)
     |    coll_en2v(k,i) = 2.53e-6*root_te(k,i)*exp(-17620./te(k,i))
!
! Loss due to electron/n2 collision L0(E,N2,VIB)/(NE*N(N2)) (s10)
!
          loss_en2v(k,i) = 3200.*(1./te(k,i)-1./tn(k,i))
          loss_en2v(k,i) = sign(abs(loss_en2v(k,i))+del,loss_en2v(k,i))
          loss_en2v(k,i) = -3200./(te(k,i)*tn(k,i))*
     |      (1.-exp(loss_en2v(k,i)))/loss_en2v(k,i)
          loss_en2v(k,i) = 1.3e-4*loss_en2v(k,i)*coll_en2v(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1 (end DO 20)

!     call addfld('XNMBARM',' ',' ',xnmbar(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('O2N'    ,' ',' ',o2n(lev0:lev1-1,:)   ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('O1N'    ,' ',' ',o1n(lev0:lev1-1,:)   ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('N2N'    ,' ',' ',n2n(lev0:lev1-1,:)   ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('C_EN2V'  ,' ',' ',coll_en2v(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EN2V'  ,' ',' ',loss_en2v(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Calculate and sum cooling rates (heat loss) due to interactions between
! electrons/neutrals, electrons/ions, ions/neutrals
!
      do i=lon0,lon1 ! (DO 500)
        do k=lev0,lev1-1
!
! Electron/N2 loss rate:
! loss_en2 = (L0(E,N2)+L0(E,N2,ROT)+L0(E,N2,VIB))/NE (s11)
!
        loss_en2(k,i) = n2n(k,i)*(1.77E-19*(1.-1.21E-4*te(k,i))*
     |    te(k,i) + 2.9e-14/root_te(k,i) + loss_en2v(k,i))
!
! Start total of electron/neutral loss rate (s11):
!
        loss_en(k,i) = loss_en2(k,i)
!
! Electron/O2 loss rates: (L0(E,O2)+L0(E,O2,ROT)+L0(E,O2,VIB)/NE
!
          loss_eo2(k,i) = o2n(k,i)*(1.21e-18*(1.+3.6e-2*root_te(k,i))*
     |      root_te(k,i)+6.9e-14/root_te(k,i)+3.125e-21*te(k,i)**2)
          loss_en(k,i) = loss_en(k,i)+loss_eo2(k,i)
!
! Electron/O(1d) loss rates: L0(E,O,1D)/(NE*N(O))
!
          loss_eo1d(k,i) = 22713.*(1./te(k,i)-1./tn(k,i))
          loss_eo1d(k,i) = sign(abs(loss_eo1d(k,i))+del,loss_eo1d(k,i))
          loss_eo1d(k,i) = 22713./(te(k,i)*tn(k,i))*
     |      (1.-exp(loss_eo1d(k,i)))/loss_eo1d(k,i)
!
! loss_eo1d function often fails here with bad argument to exp()
! due to high te and/or high loss_eo1d from above.
!         loss_eo1d(k,i) = 1.57e-12*exp((2.4e4+0.3*(te(k,i)-1500.)-
!    |      1.947e-5*(te(k,i)-1500.)*(te(k,i)-4000.))*(te(k,i)-3000.)/
!    |      (3000.*te(k,i)))*loss_eo1d(k,i)
          loss_eo1d(k,i) = 0.
!
! Electron/O1 loss rates: (L0(E,O)+L0(E,O,F))/NE
!
          loss_eo1(k,i) = o1n(k,i)*(7.9e-19*(1.+5.7e-4*te(k,i))*
     |      root_te(k,i)+3.4e-12*(1.-7.e-5*te(k,i))/tn(k,i)*
     |      (150./te(k,i)+0.4))

          loss_en(k,i) = loss_en(k,i)+loss_eo1(k,i)
!
! loss_xen = L0*(E,N) (s8)
!
          loss_xen(k,i) = (loss_en(k,i)+o1n(k,i)*(1.-alam/(ad+sd*
     |      n2n(k,i)))*loss_eo1d(k,i))*root_ne(k,i)*evergs
!
! Complete total electron/neutral loss rate L0(E,N) (s11):
!
          loss_en(k,i) = (loss_en(k,i)+o1n(k,i)*loss_eo1d(k,i))*
     |      root_ne(k,i)*evergs
!
! Calculate L0(E,I) = L(E,I)/(TE-TI), where L(E,I) is loss due to
! interactions between electrons and ions.
!
          loss_ei(k,i) = 3.2e-8*root_ne(k,i)/(root_te(k,i)*te(k,i))*
     |      15.*(op(k,i)+0.5*o2p(k,i)+0.53*nop(k,i))*evergs

          root_tn(k,i) = sqrt(tn(k,i))
! 
! loss_in = ion/neutral cooling = L0(I,N) =L(I,N)/(TI-TN)
!
          loss_in(k,i) = ((6.6e-14*n2n(k,i)+5.8e-14*o2n(k,i)+0.21e-14*
     |      o1n(k,i)*root2*root_tn(k,i))*op(k,i)+(5.45e-14*o2n(k,i)+
     |      5.9e-14*n2n(k,i)+4.5e-14*o1n(k,i))*nop(k,i)+(5.8e-14*
     |      n2n(k,i)+4.4e-14*o1n(k,i)+0.14e-14*o2n(k,i)*root_tn(k,i))*
     |      o2p(k,i))*evergs
!
! Complete tridiagonal matrix coefficients and rhs:
!
! q_coef = q_coef-(L0(E,N)+L0(E,I))/TE**2.5 = Q
!
          q_coef(k,i) = q_coef(k,i)-(loss_en(k,i)+loss_ei(k,i))/
     |      te(k,i)**2.5
!          
! rhs = rhs-L0(E,N)*TN-L0(E,I)*TI
!
          rhs(k,i) = rhs(k,i)-loss_en(k,i)*tn(k,i)-loss_ei(k,i)*ti(k,i)

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1 (end DO 500)

!     call addfld('L_EN2'  ,' ',' ',loss_en2 (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EO2'  ,' ',' ',loss_eo2 (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EO1D' ,' ',' ',loss_eo1d(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EO1'  ,' ',' ',loss_eo1 (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_XEN'  ,' ',' ',loss_xen (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EN'   ,' ',' ',loss_en  (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EI'   ,' ',' ',loss_ei  (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_IN'   ,' ',' ',loss_in  (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEFb' ,' ',' ',q_coef   (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('RHS1'   ,' ',' ',rhs      (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Calculate heating due to electron/neutral and electron/ion collisions
! (ergs/sec/gm):
!
      do i=lon0,lon1 ! (DO 24)
        do k=lev0,lev1-1
          if (te(k,i)-ti(k,i) >= 0.) then
            q_eni(k,i)=loss_ei(k,i)*(te(k,i)-ti(k,i))
          else
            q_eni(k,i) = 0.
          endif
          q_eni(k,i) = (loss_xen(k,i)*(te(k,i)-tn(k,i))+q_eni(k,i))
     |      *avo/xnmbar(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1 (end DO 24)

!     call addfld('Q_ENI',' ',' ',q_eni(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Add collisional heating to Q for use in thermodynamic equation.
      do i=lon0,lon1 ! (DO 27)
        do k=lev0,lev1-2
          qtotal(k+1,i,lat) = qtotal(k+1,i,lat)+
     |      .5*(q_eni(k,i)+q_eni(k+1,i))
        enddo ! k=lev0,lev1-2
!
! Upper and lower boundaries:
        qtotal(lev0,i,lat) = qtotal(lev0,i,lat)+1.5*q_eni(lev0,i)-
     |    0.5*q_eni(lev0+1,i)
        qtotal(lev1,i,lat) = qtotal(lev1,i,lat)+1.5*q_eni(lev1-1,i)-
     |    0.5*q_eni(lev1-2,i)
      enddo ! i=lon0,lon1 ! (DO 27)

!     call addfld('Q_TOT',' ',' ',qtotal(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Solve tridiagonal system:
!
!     call addfld('P_COEF' ,' ',' ',p_coef(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEF' ,' ',' ',q_coef(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEF' ,' ',' ',r_coef(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('RHS2'   ,' ',' ',rhs   (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
!     subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1,lonmax,lat,
!    |  idebug)

      call trsolv(p_coef,q_coef,r_coef,rhs,te_out(:,lon0:lon1),
     |  lev0,lev1,lev0,lev1-1,lon0,lon1,nlonp4,lat,0)
!
! Periodic points:
!     call periodic_f2d(te_out(:,lon0:lon1),lon0,lon1,lev1-lev0+1)

!     call addfld('TE_SOLV',' ',' ',te_out(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Te = Te**(2./7.):
      do i=lon0,lon1
        do k=lev0,lev1-1
          te_out(k,i) = te_out(k,i)**(2./7.)
        enddo
      enddo
!
! 10/21/03 btf: make this check after te*(2/7), rather than before.
!
! Te must be >= Tn:
      do i=lon0,lon1
        do k=lev0,lev1-1
          if (te_out(k,i) < tn(k,i)) te_out(k,i) = tn(k,i)
        enddo
      enddo
!
! 1/9/08 btf: put spval in top level of te:
      te_out(lev1,:) = spval
!
! Te is not defined at lev1 (only up to lev1-1)
!     call addfld('TE_OUT',' ',' ',te_out(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Set ion temperature output. Use joule heating qji_ti from sub 
! qjoule_ti (see qjoule.F). lev1 not calculated.
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          ti_out(k,i) = (qji_ti(k,i)*(xnmbar(k,i)/avo)+
     |      loss_ei(k,i)*te_out(k,i)+loss_in(k,i)*tn(k,i))/
     |      (loss_ei(k,i)+loss_in(k,i))
!
! ti must be at least as large as tn:
          if (ti_out(k,i) < tn(k,i)) ti_out(k,i) = tn(k,i)
        enddo
      enddo
!
! 1/9/08 btf: put spval in top level of ti:
      ti_out(lev1,:) = spval
!
!     call addfld('TI_OUT',' ',' ',ti_out(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
#ifdef VT
!     code = 126 ; state = 'settei' ; activity='ModelCode'
      call vtend(126,ier)
#endif
      end subroutine settei
