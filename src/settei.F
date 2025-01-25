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
      use cons_module,only: pi,dtr,evergs,
     |  rmassinv_o2,rmassinv_o1,rmassinv_he,rmassinv_n2,dipmin,avo,
     |  rmass_o1,rmassinv_n4s,rmassinv_no
      use input_module,only: f107,et,electron_heating, ti_cap, te_cap
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
      real(rp),parameter ::
     |  fpolar = -3e9,   ! polar te flux
     |  del    = 1e-6
      integer :: k,i,ier
      real(rp) :: r,
     |  a,fed,fen,       ! day/night
     |  sindipmag        ! sin(dipmag)
      real(rp),dimension(lon0:lon1) :: fe ! heat flux at upper boundary
      real(rp),dimension(lev0:lev1,lon0:lon1) ::
     |  o2i,o1i,hei,n2i, ! mass mixing ratios at interfaces
     |  o2_cm3,o2i_cm3,  ! O2 number density (midpoints or interfaces)
     |  o1_cm3,o1i_cm3,  ! O1 number density (midpoints or interfaces)
     |  he_cm3,hei_cm3,  ! He number density (midpoints or interfaces)
     |  n2_cm3,n2i_cm3,  ! N2 number density (midpoints or interfaces)
     |  tei,             ! electron temperature at interfaces
     |  nei,             ! eletron density at midpoints
     |  qe,              ! source term
     |  root_te,         ! sqrt(te)
!
! Cooling rates (heat loss):
     |  loss_eo2,        ! electron/O2 loss
     |  loss_eo1,        ! electron/O loss
     |  loss_ehe,        ! electron/He loss
     |  coll_en2v,       ! electron/N2 vibrational collision
     |  loss_en2v,       ! electron/N2 vibrational loss
     |  loss_en2,        ! electron/N2 loss
     |  loss_en,         ! electron/neutral loss
     |  loss_ei,         ! electron/ion loss
!
     |  tek0,            ! ke/te**2.5
     |  tek0_h,tek0_h_1, ! tek0/scht
     |  qtot,qtoti,      ! total ionization rate
!
     |  coef_part,
     |  p_coef,          ! coefficient for trisolv
     |  q_coef,          ! coefficient for trisolv
     |  r_coef,          ! coefficient for trisolv
     |  rhs,             ! right-hand-side for trisolv
!
     |  root_tn,         ! sqrt(2.*tn)
     |  loss_in,         ! ion/neutral loss
     |  q_eni,q_eni_i    ! heating from electron/neutral and electron/ion collisions
!
#ifdef VT
!     code = 126 ; state = 'settei' ; activity='ModelCode'
      call vtbegin(126,ier)
#endif
!
      do i=lon0,lon1
        if (abs(rlatm(i,lat)) >= pi/4.5_rp) then
          a = 1._rp
!
! Dang, 2019, set a=0 between +10 and -10 magnetic latitude
        elseif (abs(rlatm(i,lat)) <= pi/18._rp) then
          a = 0._rp
        else 
!
! Pham, 2021, fix huge gradients in low latitude
          a = .5_rp*(1._rp-cos(abs(rlatm(i,lat))*6._rp-pi/3._rp))
        endif
!
! Dang, 2019, FeDCoef2 yields a non-zero downward flux at the equator.
! Remove per Wenbin's suggestion, and FeDCoef1 is changed from -5e7 to -9e7.
        fed = -9.0e+7_rp*min(f107,235._rp)*a
        fen = fed/2._rp
        fed = fed+qteaur(i,lat)
        fen = fen+qteaur(i,lat)
!
        if (chi(i,lat) >= 100._rp*dtr) then
          fe(i) = fen
        elseif (chi(i,lat) <= 80._rp*dtr) then
          fe(i) = fed
        else
          fe(i) = .5_rp*(fed+fen+(fed-fen)*cos(chi(i,lat)*9._rp))
        endif
!
! Add fpolar if magnetic latitude >= 60 degrees:
        if (abs(rlatm(i,lat)) >= pi/3._rp) fe(i) = fe(i)+fpolar
!
        fe(i) = fe(i)*evergs
      enddo ! i=lon0,lon1
!
! te,o2,o,he,n2 at interfaces: 
      do i=lon0,lon1
        do k=lev0+1,lev1-1
          tei(k,i) = .5_rp*(te(k,i)+te(k-1,i))
        enddo ! k=lev0+1,lev1-1
!
! Bottom:
! 2023/02 Dong Lin: add floor for te (Te>=Tn)
        tei(lev0,i) = 1.5_rp*te(lev0,i)-.5_rp*te(lev0+1,i)
        if (tei(lev0,i) < tlbc(i,lat)) tei(lev0,i) = tlbc(i,lat)
!
! Top:
        tei(lev1,i) = 1.5_rp*te(lev1-1,i)-.5_rp*te(lev1-2,i)
!
        do k=lev0+1,lev1
          o2i(k,i) = .5_rp*(o2(k,i)+o2(k-1,i))
          o1i(k,i) = .5_rp*(o1(k,i)+o1(k-1,i))
          hei(k,i) = .5_rp*(he(k,i)+he(k-1,i))
          n2i(k,i) = .5_rp*(n2(k,i)+n2(k-1,i))
        enddo ! k=lev0+1,lev1
!
! 2024/04 Haonan Wu: change lower boundary calculation of major species
        o2i(lev0,i)    = .5_rp*(fb(1)+
     |    (b(1,1)+1._rp)*o2(lev0,i)+
     |     b(1,2)    *o1(lev0,i)+
     |     b(1,3)    *he(lev0,i))
        o1i(lev0,i)    = .5_rp*(fb(2)+
     |     b(2,1)    *o2(lev0,i)+
     |    (b(2,2)+1._rp)*o1(lev0,i)+
     |     b(2,3)    *he(lev0,i))
        hei(lev0,i)    = .5_rp*(fb(3)+
     |     b(3,1)    *o2(lev0,i)+
     |     b(3,2)    *o1(lev0,i)+
     |    (b(3,3)+1._rp)*he(lev0,i))
        n2i(lev0,i)    =
     |    max(1._rp-o2i(lev0,i)-o1i(lev0,i)-hei(lev0,i),0._rp)
      enddo ! i=lon0,lon1
!
!     call addfld('TE_INT'  ,' ',' ',tei,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O2_INT'  ,' ',' ',o2i   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O1_INT'  ,' ',' ',o1i   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('HE_INT'  ,' ',' ',hei   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2_INT'  ,' ',' ',n2i   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Convert o2,o,he,n2 to number density (midpoints and interfaces):
      do i=lon0,lon1
        do k=lev0,lev1
          o2_cm3(k,i) = xnmbar(k,i)*o2(k,i)*rmassinv_o2
          o1_cm3(k,i) = xnmbar(k,i)*o1(k,i)*rmassinv_o1
          he_cm3(k,i) = xnmbar(k,i)*he(k,i)*rmassinv_he
          n2_cm3(k,i) = xnmbar(k,i)*n2(k,i)*rmassinv_n2
          o2i_cm3(k,i) = xnmbari(k,i)*o2i(k,i)*rmassinv_o2
          o1i_cm3(k,i) = xnmbari(k,i)*o1i(k,i)*rmassinv_o1
          hei_cm3(k,i) = xnmbari(k,i)*hei(k,i)*rmassinv_he
          n2i_cm3(k,i) = xnmbari(k,i)*n2i(k,i)*rmassinv_n2
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
!     call addfld('XNMBAR',' ',' ',xnmbar(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('XNMBARI',' ',' ',xnmbari(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O2N',' ',' ',o2_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O1N',' ',' ',o1_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('HEN',' ',' ',he_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2N',' ',' ',n2_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1
          root_te(k,i) = sqrt(tei(k,i))
!
          tek0(k,i) = 7.5e5_rp*evergs/
     |      (1._rp + 3.22e4_rp*tei(k,i)**2/ne(k,i)*
     |      ((2.20e-16_rp + 7.92e-18_rp*root_te(k,i))*o2i_cm3(k,i)+
     |        1.10e-16_rp * (1._rp + 5.7e-4_rp*tei(k,i))*o1i_cm3(k,i)+
     |        5.60e-16_rp * hei_cm3(k,i)+
     |       (2.82e-17_rp - 3.41e-21_rp*tei(k,i))*root_te(k,i)*
     |       n2i_cm3(k,i)))
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
!     call addfld('TEK0'   ,' ',' ',tek0  ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1
          tek0_h(k,i) = tek0(k,i)/schti(k,i)
        enddo ! k=lev0,lev1
        do k=lev0,lev1-1
          tek0_h_1(k,i) = tek0_h(k+1,i)
        enddo ! k=lev0,lev1-1
        tek0_h_1(lev1,i) = 2*tek0_h(lev1,i)-tek0_h(lev1-1,i)
      enddo ! i=lon0,lon1
!
!     call addfld('TEK0_H',' ',' ',tek0_h,
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('TEK0_H_1',' ',' ',tek0_h_1,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        sindipmag = max(sin(dipmag(i,lat))**2,sin(dipmin)**2,.10_rp)
!
! Start coefficients and rhs for trsolv:
        do k=lev0,lev1-1
          coef_part(k,i) = 2._rp*sindipmag/(7._rp*scht(k,i)*dz**2)
          p_coef(k,i) = coef_part(k,i)*tek0_h(k,i)
          r_coef(k,i) = coef_part(k,i)*tek0_h_1(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! qtot = total ionization rate = sum(Qxx) = 
! (QO2+) + (QO+) + (QN2+) + (QNO+) + (QN+) + (QO+(2D)) + (QO+(2P))
!
!     call addfld('QO2P' ,' ',' ',qo2p (:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP'  ,' ',' ',qop  (:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QN2P' ,' ',' ',qn2p (:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QNOP' ,' ',' ',qnop (:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QNP'  ,' ',' ',qnp  (:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP2D',' ',' ',qop2d(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP2P',' ',' ',qop2p(:,lon0:lon1,lat),
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
      do i=lon0,lon1
        do k=lev0,lev1
          qtot(k,i) = qo2p(k,i,lat)+qop(k,i,lat)+qn2p(k,i,lat)+
     |      qnop(k,i,lat)+qnp(k,i,lat)+qop2d(k,i,lat)+qop2p(k,i,lat)
          if (qtot(k,i) < 1.e-20_rp) qtot(k,i) = 1.e-20_rp
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
!     call addfld('QTOT_SUM',' ',' ',qtot,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          qtoti(k,i) = sqrt(qtot(k,i)*qtot(k+1,i))
        enddo ! k=lev0,lev1-1
        qtoti(lev1,i) = sqrt(qtot(lev1,i)**3/qtot(lev1-1,i))
      enddo ! i=lon0,lon1
!
!     call addfld('QTOT',' ',' ',qtoti   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Set Ne (midpoints "(K+1/2)"):
      do i=lon0,lon1
        do k=lev0,lev1-1
          nei(k,i) = sqrt(ne(k,i)*ne(k+1,i))
        enddo ! k=lev0,lev1-1
        nei(lev1,i) = sqrt(ne(lev1,i)**3/ne(lev1-1,i))
        do k=lev0,lev1
          if (nei(k,i) < 1.e2_rp) nei(k,i) = 1.e2_rp
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
!     call addfld('NE_MID',' ',' ',nei,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! Calculate source term qe
          if (electron_heating == 6) then
!
! Modified according to Smithtro & Solomon (2008), by Yihui Cai, 2021/10
            r = log(nei(k,i)/(o2_cm3(k,i)+n2_cm3(k,i)+o1_cm3(k,i)))
            qe(k,i) = exp((((((-1.249e-5_rp*r-5.755e-4_rp)*r-
     |        9.346e-3_rp)*r-5.900e-2_rp)*r-
     |        4.392e-2_rp)*r+1.056_rp)*r+5.342_rp)
          else
!
! Comment from earlier version (maybe the *1.0 below was once *2.0):
!   "Correction facor of 2 increase in TE heating rate"
            r = log(nei(k,i)/(o2_cm3(k,i)+
     |        n2_cm3(k,i)+0.1_rp*o1_cm3(k,i)))
            qe(k,i) =
     |        exp(-((((0.001996_rp*r+0.08034_rp)*r+
     |        1.166_rp)*r+6.941_rp)*r+12.75_rp))*1.0_rp
          endif
!
! Calculate and sum cooling rates (heat loss) due to interactions between
! electrons/neutral, electrons/ions
          root_te(k,i) = sqrt(te(k,i))
!
! Electron/O2 loss rates: (L0(E,O2)+L0(E,O2,ROT)+L0(E,O2,VIB)/NE
          loss_eo2(k,i) = 1.21e-18_rp*(root_te(k,i)+3.6e-2_rp*te(k,i))+
     |      6.9e-14_rp/root_te(k,i)+3.125e-21_rp*te(k,i)**2
!
! Electron/O(1d) loss rates: L0(E,O,1D)/(NE*N(O))
! skipped (not used)
!
! Electron/O1 loss rates: (L0(E,O)+L0(E,O,F))/NE
          loss_eo1(k,i) =
     |      7.9e-19_rp*(1._rp+5.7e-4_rp*te(k,i))*root_te(k,i)+
     |      3.4e-12_rp*(1._rp-7.e-5_rp*te(k,i))/tn(k,i)*
     |      (150._rp/te(k,i)+0.4_rp)
!
! Electron/He loss rate: L0(E,HE)/NE
          loss_ehe(k,i) = 2.46e-17_rp*root_te(k,i)
!
! Electron/N2 collision A(E,N2,VIB):
          if (te(k,i) > 2000._rp) then
            coll_en2v(k,i) =
     |        2.53e-6_rp*root_te(k,i)*exp(-17620._rp/te(k,i))
          elseif (te(k,i) >= 1000._rp) then
            coll_en2v(k,i) = 2.e-7_rp*exp(-4605.2_rp/te(k,i))
          else
            coll_en2v(k,i) = 5.71e-8_rp*exp(-3352.6_rp/te(k,i))
          endif
!
! Loss due to electron/N2 collision L0(E,N2,VIB)/(NE*N(N2))
!
! Separate the case when te is approaching tn to avoid calculating 0/0
          if (abs(te(k,i)-tn(k,i)) < del) then
            loss_en2v(k,i) = 3200._rp/tn(k,i)**2
          else
            loss_en2v(k,i) = 1._rp/(te(k,i)-tn(k,i))*
     |        (1._rp-exp(-3200._rp*(te(k,i)-tn(k,i))/(te(k,i)*tn(k,i))))
          endif
!
! Electron/N2 loss rate: (L0(E,N2)+L0(E,N2,ROT)+L0(E,N2,VIB))/NE
          loss_en2(k,i) =
     |      1.77E-19_rp*(1._rp-1.21E-4_rp*te(k,i))*te(k,i)+
     |      2.9e-14_rp/root_te(k,i)+
     |      1.3e-4_rp*loss_en2v(k,i)*coll_en2v(k,i)
!
! 2024/03 Haonan Wu: there were calculations in previous versions on
! electron/O(1D) vibrational cooling, but it was zeroed out,
! no documentation was found in the model description either,
! therefore all calculations regarding electron/O(1D) are discarded.
!
! Total electron/neutral loss rate:
          loss_en(k,i) = nei(k,i)*evergs*
     |      (o2_cm3(k,i)*loss_eo2(k,i)+o1_cm3(k,i)*loss_eo1(k,i)+
     |       he_cm3(k,i)*loss_ehe(k,i)+n2_cm3(k,i)*loss_en2(k,i))
!
! Cooling correction from electrojet turbulent heating
          if (et .and. te(k,i)>500.0_rp .and. Q2(k,i)>0.0_rp)
     |      loss_en(k,i) = loss_en(k,i)*
     |      exp(-7.54E-4_rp*(te(k,i)-500.0_rp))
!
! Calculate L0(E,I) = L(E,I)/(TE-TI), where L(E,I) is loss due to
! interactions between electrons and ions.
          loss_ei(k,i) = 3.2e-8_rp*nei(k,i)/(root_te(k,i)*te(k,i))*
     |      15._rp*evergs*rmass_o1*
     |      (op   (k,i)*rmassinv_o1+
     |       o2p  (k,i)*rmassinv_o2+
     |       nplus(k,i)*rmassinv_n4s+
     |       n2p  (k,i)*rmassinv_n2+
     |       nop  (k,i)*rmassinv_no)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!     call addfld('L_EO2',' ',' ',loss_eo2(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EO1',' ',' ',loss_eo1(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EHE',' ',' ',loss_ehe(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('C_EN2V',' ',' ',coll_en2v(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EN2V',' ',' ',loss_en2v(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EN2',' ',' ',loss_en2(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EN',' ',' ',loss_en(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_EI',' ',' ',loss_ei(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Complete tridiagonal matrix coefficients and rhs:
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! q_coef = -p_coef-r_coef-(L0(E,N)+L0(E,I))/TE**2.5 = Q
          q_coef(k,i) = -p_coef(k,i)-r_coef(k,i)-
     |      (loss_en(k,i)+loss_ei(k,i))/te(k,i)**2.5_rp
!
! rhs = -qe-L0(E,N)*TN-L0(E,I)*TI
          rhs(k,i) = -qe(k,i)*qtoti(k,i)*evergs-
     |      loss_en(k,i)*tn(k,i)-loss_ei(k,i)*ti(k,i)
!
! Subtract electrojet turbulent heating from rhs
          if (et) rhs(k,i) = rhs(k,i)-(Q1(k,i)+Q2(k,i))*10.0_rp
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
      do i=lon0,lon1
!
! Lower boundary:
        q_coef(lev0,i) = q_coef(lev0,i)-p_coef(lev0,i)
        rhs(lev0,i) = rhs(lev0,i)-
     |    2._rp*p_coef(lev0,i)*tlbc(i,lat)**3.5_rp
        p_coef(lev0,i) = 0._rp
!
! Upper boundary:
        q_coef(lev1-1,i) = q_coef(lev1-1,i)+r_coef(lev1-1,i)
        rhs(lev1-1,i) = rhs(lev1-1,i)+
     |    r_coef(lev1-1,i)*dz*3.5_rp*fe(i)/tek0_h(lev1,i)
        r_coef(lev1-1,i) = 0._rp
      enddo ! i=lon0,lon1
!
!     call addfld('P_COEF',' ',' ',p_coef(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEF',' ',' ',q_coef(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('R_COEF',' ',' ',r_coef(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('RHS'   ,' ',' ',rhs   (lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Solve tridiagonal system:
      call trsolv(p_coef,q_coef,r_coef,rhs,te_out(:,lon0:lon1),
     |  lev0,lev1,lev0,lev1-1,lon0,lon1,nlonp4,lat,0)
!
!     call addfld('TE_SOLV',' ',' ',te_out(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Te = Te**(2./7.):
      do i=lon0,lon1
        do k=lev0,lev1-1
          te_out(k,i) = te_out(k,i)**(2._rp/7._rp)
!
! 10/21/03 btf: make this check after te*(2/7), rather than before.
!
! Te must be >= Tn:
          if (te_out(k,i) < tn(k,i)) te_out(k,i) = tn(k,i)
!
! 2023/02 Dong Lin: N2 cross section becomes negative
! when te>2.82e-17/3.41e-21=8269.7K. With ET enabled,
! Te can go very high in the E region causing problems.
! The following Te cap can help prevent from such problem
          if (te_out(k,i) > te_cap) te_out(k,i) = te_cap
        enddo ! k=lev0,lev1-1
!
! 1/9/08 btf: put spval in top level of te:
        te_out(lev1,i) = spval
      enddo ! i=lon0,lon1
!
! Te is not defined at lev1 (only up to lev1-1)
!     call addfld('TE_OUT',' ',' ',te_out(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          root_tn(k,i) = sqrt(2*tn(k,i))
!
! loss_in = ion/neutral cooling = L0(I,N) = L(I,N)/(TI-TN)
! (rates not on Rees & Roble 1975 were calculated based on Banks 1966)
          loss_in(k,i) = 1e-14_rp*evergs*
     |      (op(k,i)*(5.8_rp*o2_cm3(k,i)+
     |       0.21_rp*o1_cm3(k,i)*root_tn(k,i)+
     |       2.8_rp*he_cm3(k,i)+
     |       6.6_rp*n2_cm3(k,i))+
     |      o2p(k,i)*(0.14_rp*o2_cm3(k,i)*root_tn(k,i)+
     |       4.36_rp*o1_cm3(k,i)+
     |       1.63_rp*he_cm3(k,i)+
     |       5.81_rp*n2_cm3(k,i))+
     |      nplus(k,i)*(5.84_rp*o2_cm3(k,i)+
     |       5.84_rp*o1_cm3(k,i)+
     |       3.05_rp*he_cm3(k,i)+
     |       6.56_rp*n2_cm3(k,i))+
     |      n2p(k,i)*(5.54_rp*o2_cm3(k,i)+
     |       4.65_rp*o1_cm3(k,i)+
     |       1.82_rp*he_cm3(k,i)+
     |       0.27_rp*n2_cm3(k,i)*root_tn(k,i))+
     |      nop(k,i)*(5.45_rp*o2_cm3(k,i)+
     |       4.5_rp*o1_cm3(k,i)+
     |       1.72_rp*he_cm3(k,i)+
     |       5.92_rp*n2_cm3(k,i)))
!
! Set ion temperature output. Use joule heating qji_ti from sub 
! qjoule_ti (see qjoule.F). lev1 not calculated.
          ti_out(k,i) = (qji_ti(k,i)*xnmbar(k,i)/avo+
     |      loss_ei(k,i)*te_out(k,i)+loss_in(k,i)*tn(k,i))/
     |      (loss_ei(k,i)+loss_in(k,i))
!
! ti must be at least as large as tn:
          if (ti_out(k,i) < tn(k,i)) ti_out(k,i) = tn(k,i)
!
! 2024/09 Haonan Wu: add Ti cap in a similar way as Te cap
! This helps the model to stay stable in some extreme cases
            if (ti_out(k,i) > ti_cap) ti_out(k,i) = ti_cap
        enddo ! k=lev0,lev1-1
!
! 1/9/08 btf: put spval in top level of ti:
        ti_out(lev1,i) = spval
      enddo ! i=lon0,lon1
!
!     call addfld('QJI_TI',' ',' ',qji_ti(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('L_IN',' ',' ',loss_in(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('TI_OUT',' ',' ',ti_out(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Calculate heating due to electron/neutral and electron/ion collisions
! (ergs/sec/gm):
      do i=lon0,lon1
        do k=lev0,lev1-1
          q_eni(k,i) = (loss_en(k,i)*(te(k,i)-tn(k,i))+
     |      loss_ei(k,i)*max(te(k,i)-ti(k,i),0._rp))*avo/xnmbar(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!     call addfld('Q_ENI',' ',' ',q_eni(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1-2
          q_eni_i(k+1,i) = .5_rp*(q_eni(k,i)+q_eni(k+1,i))
        enddo ! k=lev0,lev1-2
!
! Upper and lower boundaries:
        q_eni_i(lev0,i) = 1.5_rp*q_eni(lev0  ,i)-0.5_rp*q_eni(lev0+1,i)
        q_eni_i(lev1,i) = 1.5_rp*q_eni(lev1-1,i)-0.5_rp*q_eni(lev1-2,i)
!
! Add collisional heating to Q for use in thermodynamic equation.
        do k=lev0,lev1
          qtotal(k,i,lat) = qtotal(k,i,lat)+q_eni_i(k,i)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
!     call addfld('Q_TOT',' ',' ',qtotal(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
#ifdef VT
!     code = 126 ; state = 'settei' ; activity='ModelCode'
      call vtend(126,ier)
#endif
!
      end subroutine settei
