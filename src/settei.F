! 
      subroutine settei(tn,o2,o1,he,n2,ne,te,ti,op,o2p,nplus,n2p,nop,
     |  xnmbar,xnmbari,scht,schti,qji_ti,Q1,Q2,te_out,ti_out,
     |  lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate electron and ion temperatures.
!
      use params_module,only: dz,nlonp4,spval,rp
      use cons_module,only: pi,rtd,evergs,
     |  rmassinv_o2,rmassinv_o1,rmassinv_he,rmassinv_n2,dipmin,avo,
     |  rmass_o1,rmassinv_n4s,rmassinv_no
      use input_module,only: f107,do_et,electron_heating, ti_cap, te_cap
      use chapman_module,only: chi   ! solar zenith angle (nlonp4,nlat)
      use magfield_module,only: rlatm,dipmag
      use fields_module,only: tlbc
      use lbc,only: fb,b
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
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  he,      ! helium (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  ne,      ! electron density (cm3)
     |  te,      ! electron temperature (from previous time step)
     |  ti,      ! ion temperature (from previous time step)
     |  op,      ! O+
     |  o2p,     ! O2+
     |  nplus,   ! N+
     |  n2p,     ! N2+
     |  nop,     ! NO+
     |  xnmbar,  ! p0*e(-z)*mbar/kT  (midpoints)
     |  xnmbari, ! p0*e(-z)*barm/kT  (interfaces)
     |  scht,    ! scale height (midpoints)
     |  schti,   ! scale height (interfaces)
     |  qji_ti,  ! joule heating from qjoule_ti (used ui,vi)
     |  Q1,      ! electron ohmic heating
     |  Q2       ! electrojet turbulent heating
!
! Output args:
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
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
      integer :: k,i,ier
      integer :: nk,nkm1
      real(rp),dimension(lev0:lev1,lon0:lon1) ::
     |  te_int,    ! electron temperature (interfaces)
     |  tn_int,    ! neutral temperature (interfaces)
     |  o2n,       ! O2 number density (midpoints or interfaces)
     |  o1n,       ! O1 number density (midpoints or interfaces)
     |  hen,       ! He number density (midpoints or interfaces)
     |  n2n,       ! N2 number density (midpoints or interfaces)
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
     |  qe,        ! source term                 (s10)
     |  q_eni,     ! heating from electron/neutral and electron/ion collisions
     |  coll_en2v, ! electron/N2vib collision    (s9)
!
! Cooling rates (heat loss):
     |  loss_en2v, ! electron/N2vib loss term    (s10)
     |  loss_eo2,  ! electron/O2 loss
     |  loss_eo1d, ! electron/O(1d) loss
     |  loss_eo1,  ! electron/O loss
     |  loss_ehe,  ! electron/He loss
     |  loss_en2,  ! electron/N2 loss
     |  loss_xen,  ! L0*(E,N) (s8)
     |  loss_en,   ! electrons/neutrals loss     (s11)
     |  loss_ei,   ! electron/ion loss           (s10)
     |  loss_in    ! ion/neutral loss            (s9)
      real(rp),parameter ::
     |  fpolar = -3.0e+9,  ! polar te flux
     |  del    = 1.e-6  ,
!
! Correction factors for neutral heating due to L(E,O1D)
     |  alam = 0.0069   ,
     |  ad   = 0.0091   ,
     |  sd   = 2.3e-11
      real(rp) ::
     |  f107te  ! solar flux
!
! a,fed,fen,fe,sindipmag have a z dimension only for diagnostic plotting:
      real(rp),dimension(lon0:lon1) ::
     |  a,fed,fen, ! day/night
     |  fe,        ! heat flux at upper boundary
     |  sindipmag  ! sin(dipmag)
!
! For diagnostic plotting:
      real(rp),dimension(lev0:lev1-1,lon0:lon1) ::
     |  a_ki,        ! for diagnostic plotting of a
     |  fed_ki,      ! for diagnostic plotting of fed
     |  fen_ki,      ! for diagnostic plotting of fen
     |  fe_ki,       ! for diagnostic plotting of fe
     |  dipmag_ki,   ! for diagnostic plotting of dipmag
     |  chi_ki,      ! for diagnostic plotting of chi
     |  qteaur_ki,   ! for diagnostic plotting of qteaur
     |  sindipmag_ki ! for diagnostic plotting of sindipmag
!
#ifdef VT
!     code = 126 ; state = 'settei' ; activity='ModelCode'
      call vtbegin(126,ier)
#endif
!
      f107te = f107
      if (f107te > 235._rp) f107te = 235._rp
      nk = lev1-lev0+1
      nkm1 = nk-1

!     call addfld('QJI_TI',' ',' ',qji_ti(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        if (abs(rlatm(i,lat)) >= pi/4.5_rp) then
          a(i) = 1._rp
!       else
!         a(i) = .5*(1.+sin(pi*(abs(rlatm(i,lat))-pi/9.)/(pi/4.5)))
! Dang, 2019, set a(i)=0 between +10 and -10 magnetic latitude          
        elseif (abs(rlatm(i,lat)) <= pi/18._rp) then
          a(i) = 0._rp
        else 
! Dang, 2019, change parameter: pi/9.->pi/12., pi/4.5->6, per Hanli and Wenbin
!         a(i) = .5*(1.+sin(pi*(abs(rlatm(i,lat))-pi/12.)/(pi/6.)))
! Pham, 2021, try to fix huge gradients in low lat 
          a(i) = .5_rp*(1._rp+cos(abs(rlatm(i,lat))*
     |           6._rp-pi*4._rp/3._rp))
        endif 
!
! Increased heat flux for TE fom protonosphere.
!       fed(i) = ( -5.0e+7*f107te*a(i)-4.0e+7*f107te)*1.2
! Dang, 2019, FeDCoef2 yields a non-zero downward flux at the equator.
! Remove per Wenbin's suggestion, and FeDCoef1 is changed from -5e7 to -9e7.
        fed(i) = -9.0e+7_rp*f107te*a(i)
        fen(i) = fed(i)/2._rp
        fed(i) = fed(i)+qteaur(i,lat)     ! t4
        fen(i) = fen(i)+qteaur(i,lat)     ! t5
        if (chi(i,lat)-.5_rp*pi >= 0._rp) then  ! chi==t2
          fe(i) = fen(i)                ! t1
        else
          fe(i) = fed(i)
        endif
        if ((chi(i,lat)*rtd-80._rp)*
     |    (chi(i,lat)*rtd-100._rp)>=0._rp) then
          fe(i) = fe(i)*evergs
        else
          fe(i) = (.5*(fed(i)+fen(i))+.5_rp*(fed(i)-fen(i))*
     |      cos(pi*(chi(i,lat)*rtd-80._rp)/20._rp))*evergs
        endif
!
! Add fpolar if magnetic latitude >= 60 degrees:
        if (abs(rlatm(i,lat))-pi/3._rp>=0._rp)
     |      fe(i) = fe(i)+fpolar*evergs
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
          te_int(k,i) = .5_rp*(te(k,i)+te(k-1,i))
          o2n(k,i)    = .5_rp*(o2(k,i)+o2(k-1,i))
          o1n(k,i)    = .5_rp*(o1(k,i)+o1(k-1,i))
          hen(k,i)    = .5_rp*(he(k,i)+he(k-1,i))
          n2n(k,i)    = .5_rp*(n2(k,i)+n2(k-1,i))
          tn_int(k,i) = .5_rp*(tn(k,i)+tn(k-1,i))
        enddo ! k=lev0+1,lev1-2
!
! Bottom:
! 2023/02 Dong Lin: added minimum cap for te (Te>=Tn)
! 2024/04 Haonan Wu: changed the lower boundary calculation of major species
        te_int(lev0,i) =
     |    max(1.5_rp*te(lev0,i)-.5_rp*te(lev0+1,i),tlbc(i,lat))
        o2n(lev0,i)    = .5*(fb(1)+
     |    (b(1,1)+1._rp)*o2(lev0,i)+
     |     b(1,2)    *o1(lev0,i)+
     |     b(1,3)    *he(lev0,i))
        o1n(lev0,i)    = .5_rp*(fb(2)+
     |     b(2,1)    *o2(lev0,i)+
     |    (b(2,2)+1._rp)*o1(lev0,i)+
     |     b(2,3)    *he(lev0,i))
        hen(lev0,i)    = .5_rp*(fb(3)+
     |     b(3,1)    *o2(lev0,i)+
     |     b(3,2)    *o1(lev0,i)+
     |    (b(3,3)+1._rp)*he(lev0,i))
        n2n(lev0,i)    =
     |    max(1._rp-o2n(lev0,i)-o1n(lev0,i)-hen(lev0,i),0._rp)
        tn_int(lev0,i) = tlbc(i,lat)
!
! Top:
        te_int(lev1,i) = 1.5_rp*te(lev1-1,i)-.5*te(lev1-2,i)
        o2n(lev1,i)    = .5_rp*(o2(lev1,i)+o2(lev1-1,i))
        o1n(lev1,i)    = .5_rp*(o1(lev1,i)+o1(lev1-1,i))
        hen(lev1,i)    = .5_rp*(he(lev1,i)+he(lev1-1,i))
        n2n(lev1,i)    = .5_rp*(n2(lev1,i)+n2(lev1-1,i))
        tn_int(lev1,i) = 1.5_rp*tn(lev1-1,i)-.5_rp*tn(lev1-2,i)
      enddo ! i=lon0,lon1

!     call addfld('TE_INT'  ,' ',' ',te_int,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TN_INT'  ,' ',' ',tn_int,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O2_INT'  ,' ',' ',o2n   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O1_INT'  ,' ',' ',o1n   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('HE_INT'  ,' ',' ',hen   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2_INT'  ,' ',' ',n2n   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Convert o2,o,he,n2 to number density (interfaces):
      do i=lon0,lon1
        do k=lev0,lev1
          o2n(k,i) = xnmbari(k,i)*o2n(k,i)*rmassinv_o2       ! s13
          o1n(k,i) = xnmbari(k,i)*o1n(k,i)*rmassinv_o1       ! s12
          hen(k,i) = xnmbari(k,i)*hen(k,i)*rmassinv_he
          n2n(k,i) = xnmbari(k,i)*n2n(k,i)*rmassinv_n2       ! s11
          root_te(k,i) = sqrt(te_int(k,i))
!
          tek0(k,i) = 7.5e5_rp/
     |      (1._rp+3.22e4_rp*te_int(k,i)**2/ne(k,i)*
     |      ((2.20e-16_rp + 7.92e-18_rp   * root_te(k,i))*o2n(k,i)+
     |        1.10e-16_rp * (1._rp+5.7e-4_rp * te_int (k,i))*o1n(k,i)+
     |        5.60e-16_rp * hen(k,i)+
     |       (2.82e-17_rp - 3.41e-21_rp   * te_int (k,i))*root_te(k,i)*
     |        n2n(k,i)))*evergs
!
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     call addfld('XNMBARI',' ',' ',xnmbari(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TEK0'   ,' ',' ',tek0  ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      do i=lon0,lon1
        do k=lev0,lev1-1
          h_mid(k,i) = scht(k,i)  ! s7
        enddo ! k=lev0,lev1-1
        do k=lev0,lev1
          h_int(k,i) = schti(k,i) ! s6
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
        if (sindipmag(i) < .10_rp) sindipmag(i) = .10_rp
!
! Start coefficients and rhs for trsolv:
        do k=lev0,lev1-1
          p_coef(k,i) = 2._rp/7._rp*sindipmag(i)/(h_mid(k,i)*dz**2) ! s1
          r_coef(k,i) = p_coef(k,i)*tek0(k+1,i)/h_int(k+1,i)  ! s3
          p_coef(k,i) = p_coef(k,i)*tek0(k  ,i)/h_int(k  ,i)  ! s1
          q_coef(k,i) = -(p_coef(k,i)+r_coef(k,i))            ! s2
          rhs(k,i) = 0._rp                                       ! s4
        enddo ! k=lev0,lev1-1
!
! Bottom boundary:
        q_coef(lev0,i) = q_coef(lev0,i)-p_coef(lev0,i)
        rhs(lev0,i) = rhs(lev0,i)-2._rp*p_coef(lev0,i)*
     |                tn_int(lev0,i)**3.5_rp
        p_coef(lev0,i) = 0._rp
!
! Upper boundary:
        q_coef(lev1-1,i) = q_coef(lev1-1,i)+r_coef(lev1-1,i)
        rhs(lev1-1,i) = rhs(lev1-1,i)+r_coef(lev1-1,i)*dz*3.5_rp*
     |    h_int(lev1,i)*fe(i)/tek0(lev1,i)
        r_coef(lev1-1,i) = 0._rp
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
      do i=lon0,lon1
        do k=lev0,lev1
          qtot(k,i) = qo2p(k,i,lat)+qop(k,i,lat)+qn2p(k,i,lat)+
     |      qnop(k,i,lat)+qnp(k,i,lat)+qop2d(k,i,lat)+qop2p(k,i,lat)
            if (qtot(k,i) < 1.e-20_rp) qtot(k,i) = 1.e-20_rp
        enddo
      enddo

!     call addfld('QTOT_SUM',' ',' ',qtot,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      do i=lon0,lon1
        do k=lev0,lev1-1
          qtot(k,i) = sqrt(qtot(k,i)*qtot(k+1,i))
        enddo
        qtot(lev1,i) = 0._rp
      enddo ! i=lon0,lon1
!     call addfld('QTOT',' ',' ',qtot   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Set Ne (midpoints "(K+1/2)"):
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          root_ne(k,i) = ne(k,i)*ne(k+1,i)
          if (root_ne(k,i) < 1.e4_rp) root_ne(k,i) = 1.e4_rp
          root_ne(k,i) = sqrt(root_ne(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('ROOT_NE',' ',' ',root_ne(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Set up o2,o,he,n2 number densities at midpoints:
!
      do i=lon0,lon1 ! (DO 20)
        do k=lev0,lev1-1
          o2n(k,i) = xnmbar(k,i)*o2(k,i)*rmassinv_o2  ! s14
          o1n(k,i) = xnmbar(k,i)*o1(k,i)*rmassinv_o1  ! s13
          hen(k,i) = xnmbar(k,i)*he(k,i)*rmassinv_he
          n2n(k,i) = xnmbar(k,i)*n2(k,i)*rmassinv_n2  ! s12
!
! Calculate source term qe (s10)
          if (electron_heating == 6) then
! Modified according to Smithtro & Solomon (2008), by Yihui Cai, 2021/10
            qe(k,i) = log(root_ne(k,i)/(o2n(k,i)+n2n(k,i)+o1n(k,i)))
            qe(k,i) = exp((((((-1.249e-5_rp*qe(k,i)-5.755e-4_rp)*qe(k,i)
     |                -9.346e-3_rp)*qe(k,i)-5.900e-2_rp)*qe(k,i)-
     |                4.392e-2_rp)*qe(k,i)+1.056_rp)*qe(k,i)+5.342_rp)
          else
! Comment from earlier version (maybe the *1.0 below was once *2.0):
!   "Correction facor of 2 increase in TE heating rate"
            qe(k,i) = log(root_ne(k,i)/(o2n(k,i)+
     |                n2n(k,i)+0.1_rp*o1n(k,i)))
            qe(k,i) = exp(-((((0.001996_rp*qe(k,i)+0.08034_rp)*
     |                qe(k,i)+1.166_rp)*
     |                qe(k,i)+6.941_rp)*qe(k,i)+12.75_rp))*1.0_rp
          endif
!
! Subtract qe from right-hand-side:
          rhs(k,i) = rhs(k,i)-qe(k,i)*qtot(k,i)*evergs
!
! Subtract electrojet turbulent heating from rhs
          ! Q1/Q2 are electron ohmic and turbulent heating rates,
          ! updated in lamdas, called before settei in dynamics.
          if (do_et) rhs(k,i) = rhs(k,i)-(Q1(k,i)+Q2(k,i))*10.0_rp
!
          root_te(k,i) = sqrt(te(k,i))
!
! Electron/N2 collision A(E,N2,VIB) (s9):
!
        if (te(k,i) >= 1000._rp) then
          coll_en2v(k,i) = 2.e-7_rp*exp(-4605.2_rp/te(k,i))
        else
          coll_en2v(k,i) = 5.71e-8_rp*exp(-3352.6_rp/te(k,i))
        endif
        if (te(k,i) > 2000._rp)
     |    coll_en2v(k,i) = 2.53e-6_rp*root_te(k,i)*
     |    exp(-17620._rp/te(k,i))
!
! Loss due to electron/n2 collision L0(E,N2,VIB)/(NE*N(N2)) (s10)
!
! Separate the case when te is approaching tn to avoid calculating 0/0
          if (abs(te(k,i)-tn(k,i)) < del) then
            loss_en2v(k,i) = 3200._rp/tn(k,i)**2
          else
            loss_en2v(k,i) = 1._rp/(te(k,i)-tn(k,i))*
     |        (1._rp-exp(-3200._rp*(te(k,i)-tn(k,i))/(te(k,i)*tn(k,i))))
          endif
          loss_en2v(k,i) = 1.3e-4_rp*loss_en2v(k,i)*coll_en2v(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1 (end DO 20)

!     call addfld('XNMBAR',' ',' ',xnmbar(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('O2N'    ,' ',' ',o2n(lev0:lev1-1,:)   ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('O1N'    ,' ',' ',o1n(lev0:lev1-1,:)   ,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('HEN'    ,' ',' ',hen(lev0:lev1-1,:)   ,
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
! Electron/O2 loss rates: (L0(E,O2)+L0(E,O2,ROT)+L0(E,O2,VIB)/NE
!
          loss_eo2(k,i) = o2n(k,i)*(1.21e-18_rp*(1._rp+
     |      3.6e-2_rp*root_te(k,i))*
     |      root_te(k,i)+6.9e-14_rp/root_te(k,i)+
     |      3.125e-21_rp*te(k,i)**2)
!
! Electron/O(1d) loss rates: L0(E,O,1D)/(NE*N(O))
!
!         if (abs(te(k,i)-tn(k,i)) < del) then
!           loss_eo1d(k,i) = 22713./tn(k,i)**2
!         else
!           loss_eo1d(k,i) = 1./(te(k,i)-tn(k,i))*
!    |        (1.-exp(-22713.*(te(k,i)-tn(k,i))/(te(k,i)*tn(k,i))))
!         endif
!
! loss_eo1d function often fails here with bad argument to exp()
! due to high te and/or high loss_eo1d from above.
!         loss_eo1d(k,i) = 1.57e-12*exp((2.4e4+0.3*(te(k,i)-1500.)-
!    |      1.947e-5*(te(k,i)-1500.)*(te(k,i)-4000.))*(te(k,i)-3000.)/
!    |      (3000.*te(k,i)))*loss_eo1d(k,i)
          loss_eo1d(k,i) = 0._rp
!
! Electron/O1 loss rates: (L0(E,O)+L0(E,O,F))/NE
!
          loss_eo1(k,i) = o1n(k,i)*(7.9e-19_rp*
     |      (1._rp+5.7e-4_rp*te(k,i))*
     |      root_te(k,i)+3.4e-12_rp*(1._rp-7.e-5_rp*te(k,i))/tn(k,i)*
     |      (150._rp/te(k,i)+0.4_rp))
!
! Electron/He loss rate: L0(E,HE)/NE
!
          loss_ehe(k,i) = hen(k,i)*2.46e-17_rp*root_te(k,i)
!
! Electron/N2 loss rate: (L0(E,N2)+L0(E,N2,ROT)+L0(E,N2,VIB))/NE
!
          loss_en2(k,i) = n2n(k,i)*(1.77E-19_rp*(1.-1.21E-4_rp*te(k,i))*
     |      te(k,i) + 2.9e-14_rp/root_te(k,i) + loss_en2v(k,i))
!
! Total electron/neutral loss rate (s11):
!
          loss_en(k,i) = loss_eo2(k,i)+loss_eo1(k,i)+
     |      loss_ehe(k,i)+loss_en2(k,i)
!
! Cooling correction from electrojet turbulent heating
          if (do_et .and. te(k,i)>500.0_rp .and. Q2(k,i)>0.0_rp)
     |      loss_en(k,i)=loss_en(k,i)*
     |      exp(-7.54E-4_rp*(te(k,i)-500.0_rp))
!
! loss_xen = L0*(E,N) (s8)
!
          loss_xen(k,i) = (loss_en(k,i)+o1n(k,i)*(1._rp-alam/(ad+sd*
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
          loss_ei(k,i) = 3.2e-8_rp*root_ne(k,i)/(root_te(k,i)*te(k,i))*
     |      15._rp*evergs*rmass_o1*
     |      (op(k,i)*rmassinv_o1+o2p(k,i)*rmassinv_o2+
     |      nplus(k,i)*rmassinv_n4s+n2p(k,i)*rmassinv_n2+
     |      nop(k,i)*rmassinv_no)

          root_tn(k,i) = sqrt(2*tn(k,i))
! 
! loss_in = ion/neutral cooling = L0(I,N) =L(I,N)/(TI-TN)
! (rates not on Rees & Roble 1975 were calculated based on Banks 1966)
          loss_in(k,i) = 1e-14_rp*evergs*
     |      (op(k,i)*(5.8_rp*o2n(k,i)+0.21_rp*o1n(k,i)*root_tn(k,i)+
     |      2.8_rp*hen(k,i)+6.6_rp*n2n(k,i))+
     |      o2p(k,i)*(0.14_rp*o2n(k,i)*root_tn(k,i)+
     |      4.36_rp*o1n(k,i)+1.63_rp*hen(k,i)+5.81_rp*n2n(k,i))+
     |      nplus(k,i)*(5.84_rp*o2n(k,i)+5.84_rp*o1n(k,i)+
     |      3.05_rp*hen(k,i)+6.56_rp*n2n(k,i))+
     |      n2p(k,i)*(5.54*o2n(k,i)+4.65*o1n(k,i)+
     |      1.82_rp*hen(k,i)+0.27_rp*n2n(k,i)*root_tn(k,i))+
     |      nop(k,i)*(5.45_rp*o2n(k,i)+4.5_rp*o1n(k,i)+
     |      1.72_rp*hen(k,i)+5.92_rp*n2n(k,i)))
!
! Complete tridiagonal matrix coefficients and rhs:
!
! q_coef = q_coef-(L0(E,N)+L0(E,I))/TE**2.5 = Q
!
          q_coef(k,i) = q_coef(k,i)-(loss_en(k,i)+loss_ei(k,i))/
     |      te(k,i)**2.5_rp
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
          if (te(k,i)-ti(k,i) >= 0._rp) then
            q_eni(k,i)=loss_ei(k,i)*(te(k,i)-ti(k,i))
          else
            q_eni(k,i) = 0._rp
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
     |      .5_rp*(q_eni(k,i)+q_eni(k+1,i))
        enddo ! k=lev0,lev1-2
!
! Upper and lower boundaries:
        qtotal(lev0,i,lat) = qtotal(lev0,i,lat)+1.5_rp*q_eni(lev0,i)-
     |    0.5_rp*q_eni(lev0+1,i)
        qtotal(lev1,i,lat) = qtotal(lev1,i,lat)+1.5_rp*q_eni(lev1-1,i)-
     |    0.5_rp*q_eni(lev1-2,i)
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
          te_out(k,i) = te_out(k,i)**(2._rp/7._rp)
        enddo
      enddo
!
! 10/21/03 btf: make this check after te*(2/7), rather than before.
!
! Te must be >= Tn:
      do i=lon0,lon1
        do k=lev0,lev1-1
          if (te_out(k,i) < tn(k,i)) te_out(k,i) = tn(k,i)
          ! Apply cap to prevent superstorm from exploding
          if (te_out(k,i) > te_cap) te_out(k,i) = te_cap
!
! 2023/02 Dong Lin: N2 cross section becomes negative
! when te>2.82e-17/3.41e-21=8269.7K. With ET enabled,
! Te can go very high in the E region causing problems.
! The following Te cap can help prevent from such problem
!         if (te_out(k,i) > 6000.0) te_out(k,i) = 6000.0
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
          ! Apply cap to prevent superstorm from exploding
          if (ti_out(k,i) > ti_cap) ti_out(k,i) = ti_cap
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
