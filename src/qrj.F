!
      module qrj_module
      use params_module,only: nlevp1,nlonp4
      implicit none
!
! Coefficients for wavelength scans in qrj.
! These are referenced in routines qrj, init_sflux, init_sigmas, and 
! set_qrj_coeff (all in qrj_module.F).
!
      integer,parameter :: lmax=59, l1=16
      real :: 
     |  euveff(nlevp1),  !
     |  sigeuv(3,lmax),  !
     |  rlmeuv(lmax),    !
     |  feuv(lmax),      !
     |  fsrc(l1-1),      !
     |  sigsrc(l1-1),    !
     |  rlmsrc(l1-1),    !
     |  sigmas(6,lmax),  !
     |  quench(4),       !
     |  sflux(lmax),     !
     |  brn2(lmax),      !
     |  bro2(lmax)       !
!
! For euvac, init_sigmas
      integer,parameter :: neuv=37
      real,dimension(neuv) ::
     |  wleuv1,wleuv2,
     |  sigao, sigao2, sigan2,
     |  sigio, sigio2, sigin2,
     |  sigop2p, sigop2d, sigop4s, sigin,
     |  brop4s, brop2d, brop2p, brn2np, bro2op 
!
! Heating and ionization terms set by qrj, and used by other routines.
! These are allocated for task subdomains by alloc_q (called from allocdata)
!
      real,dimension(:,:,:),allocatable ::
     |  rj,     ! total o2 dissociation
     |  qtef,   ! solar photodissociation of n2
     |  qtotal, ! total heating
     |  qop2p,  ! o+(2p)
     |  qop2d,  ! o+(2d)
     |  qo2p,   ! o2+
     |  qop,    ! o+
     |  qn2p,   ! n2+
     |  qnp,    ! n+
     |  qnop    ! no+
!
      contains
!-----------------------------------------------------------------------
      subroutine qrj(sco2,sco1,scn2,tn,no,o2,o1,n4s,xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
!
! Calculate heating and dissociation rates.
!
      use input_module,only: f107
      use init_module,only: sfeps,idn
      use cons_module,only: t0,expz,avo,rmassinv_n4s,rmassinv_no,
     |  rmassinv_o2,rmassinv_o1,rmassinv_n2,expzmid,gask,
     |  expzmid_inv,p0,check_exp
      use bndry_module,only: fb,b
      use chemrates_module,only: disn2p,beta9
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  sco2,sco1,scn2,   ! chapman integrals
     |  tn,no,o2,o1,n4s,  ! tn and species mass mixing ratios
     |  xnmbari           ! p0*e(-z)*barm/kT at interfaces
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: k,i,l,nlevs,ier
      real,parameter ::
     |  do2=8.203E-12   ,
     |  do22=1.1407E-11 ,   
     |  aband=0.143     , ! shumann-runge
     |  bband=9.64E8    , ! shumann-runge
     |  cband=9.03E-19  , ! shumann-runge
     |  e3=0.33         ,
     |  hc = 1.9845E-16   ! C(60)
      real :: rlmeuvinv(lmax),rlmsrcinv(l1-1),factor,fmin,fmax
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  r,      ! s14
     |  rp,     ! s15
     |  taur,   ! s1
     |  tauq,   ! s2
     |  tau1,   ! s3
     |  tau2,   ! s4
     |  tau3,   ! s5
     |  etaur,  ! s6
     |  etauq   ! s7
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  o2i,    ! o2  at interfaces (s1)
     |  o1i,    ! o   at interfaces (s2)
     |  n2i,    ! n2  at interfaces (s3)
     |  n4si,   ! n4s at interfaces (s4)
     |  tni     ! tn  at interfaces (s6)
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  zo2,       ! o2  ionization (s8)
     |  zo1,       ! o   ionization (s9)
     |  zn2,       ! n2  ionization (s10)
     |  zn4s,      ! n4s ionization (s11)
     |  zo2_bro2,  ! o2  ionization * bro2 (s12)
     |  zn2_brn2,  ! n2  ionization * brn2 (s13)
     |  qop2pd,    ! (s4)
     |  quenchfac, ! (s8)
     |  sigchap,   ! (s9)
     |  p3f        ! (s7)
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  sum1,   ! sum(o2,o,n2)(sigma*chapman)   (s5)
     |  sum2,   ! sum(o2,o,n2)(sigma*psi/rmass) (s6)
     |  sum3    ! sum(o2,o,n2,n4s)(sigmas)      (s7)
      logical,parameter :: 
     |  debug=.false.     ! insert print statements
!
! expo() (util.F) is used only if check_exp is true. This will avoid
! NaNS fpe, but will degrade performance. Check_exp is in cons.F.
!
      real,external :: expo ! used only when check_exp is set (util.F)
!
      if (debug) write(6,"('Enter qrj: lat=',i3,' lon0,1=',2i3)") 
     |  lat,lon0,lon1
!
#ifdef VT
!     code = 118 ; state = 'qrj' ; activity='ModelCode'
      call vtbegin(118,ier)
#endif
!
! Exec:
      nlevs = lev1-lev0+1

!     write(6,"('qrj: lat=',i2,' lon0,1=',2i3,' idn(lon0:lon1)=',/,
!    |  (15i3))") lat,lon0,lon1,idn(lon0:lon1)
!     call addfsech('SCO2',' ',' ',sco2(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('SCO1',' ',' ',sco1(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('SCN2',' ',' ',scn2(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('XNMBARI',' ',' ',xnmbari(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!
! Will multiply by inverses:
      do i=1,lmax
        rlmeuvinv(i) = 1./rlmeuv(i)
      enddo
      do i=1,l1-1 ! 1,15
        rlmsrcinv(i) = 1./rlmsrc(i)
      enddo
!
! COMPUTE S14 = RSQ, S15 = RSP
! S1 = TAU(R),  S2 = TAU(Q)
!
!     write(6,"('qrj: lat=',i2,' lon0,1=',2i3,' idn(lon0:lon1)=',/,
!    |  (12i3))") lat,lon0,lon1,idn(lon0:lon1)
!
      do i=lon0,lon1
        do k=lev0,lev1
          taur(k,i) = sigeuv(1,49)*sco2(k,i)+sigeuv(2,49)*sco1(k,i)+
     |                sigeuv(3,49)*scn2(k,i)
          tauq(k,i) = sigeuv(1,20)*sco2(k,i)+sigeuv(2,20)*sco1(k,i)+
     |                sigeuv(3,20)*scn2(k,i)
          tau1(k,i) = 1.3*taur(k,i)
          tau2(k,i) = 2.0*taur(k,i)
          tau3(k,i) = 2.5*taur(k,i)
          if (taur(k,i) > 9.) taur(k,i) = 9.
          if (tauq(k,i) > 9.) tauq(k,i) = 9.
          if (tau1(k,i) > 9.) tau1(k,i) = 9.
          if (tau2(k,i) > 9.) tau2(k,i) = 9.
          if (tau3(k,i) > 9.) tau3(k,i) = 9.
!
! TAU(N) = EXP(-TAU(N)) FOR N = 1,3,1
!
          tau1(k,i) = exp(-tau1(k,i))
          tau2(k,i) = exp(-tau2(k,i))
          tau3(k,i) = exp(-tau3(k,i))
!
! taur = EXP(-TAU(R)), tauq = EXP(-TAU(Q))
!
          etaur(k,i) = exp(-taur(k,i))
          etauq(k,i) = exp(-tauq(k,i))
!
! rp = r(o2p), r = r(op,n2p,np)
!
          r(k,i) = etaur(k,i)+2.*(tau1(k,i)+tau2(k,i)+tau3(k,i))
          rp(k,i) = 1.5*etaur(k,i)/(r(k,i)+tauq(k,i)/taur(k,i)*
     |      etauq(k,i))
          r(k,i) = 2.4*etaur(k,i)/r(k,i)
        enddo
      enddo
      if (debug) write(6,"(/,'qrj after tau: lat=',i3)") lat
!
!     call addfsech('TAUQ',' ',' ',tauq,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('TAUR',' ',' ',taur,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('R'   ,' ',' ',r   ,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('RP'  ,' ',' ',rp  ,lon0,lon1,nlevs,nlevs,lat)
!
! O2,O,N4S at interface levels:
      do i=lon0,lon1
        do k=lev0,lev1-1
          o2i (k+1,i) = 0.5*(o2(k,i)+o2(k+1,i))
          o1i (k+1,i) = 0.5*(o1(k,i)+o1(k+1,i))
          n4si(k+1,i) = 0.
        enddo
      enddo
!
! Bottom boundary:
      do i=lon0,lon1
        o2i(1,i) = .5*((b(i,1,1)+1.)*o2(1,i)+b(i,1,2)*o1(1,i)+fb(i,1))
        o1i(1,i) = .5*(b(i,2,1)*o2(1,i)+(b(i,2,2)+1.)*o1(1,i)+fb(i,2))
        n4si(1,i) = 0.
      enddo
!
! N2 at interfaces:
      do k=lev0,lev1
        n2i(k,:) = 1.-o2i(k,:)-o1i(k,:)
      enddo 
!     call addfsech('O2I'  ,' ',' ',o2i,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('O1I'  ,' ',' ',o1i,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('N2I'  ,' ',' ',n2i,lon0,lon1,nlevs,nlevs,lat)
!
! Initialize heating terms on current processor domain:
! (global module data for use by other routines)
      rj    (lev0:lev1,lon0:lon1,lat) = 0.
      qtef  (lev0:lev1,lon0:lon1,lat) = 0.  
      qtotal(lev0:lev1,lon0:lon1,lat) = 0.
      qop2p (lev0:lev1,lon0:lon1,lat) = 0.
      qop2d (lev0:lev1,lon0:lon1,lat) = 0.
      qo2p  (lev0:lev1,lon0:lon1,lat) = 0.
      qop   (lev0:lev1,lon0:lon1,lat) = 0.
      qn2p  (lev0:lev1,lon0:lon1,lat) = 0.
      qnp   (lev0:lev1,lon0:lon1,lat) = 0.
      qnop  (lev0:lev1,lon0:lon1,lat) = 0.
!
! Initialize ionization rates:
! (array ops on local arrays dimensioned (lon0:lon1,lev0:lev1)):
      zo2      = 0.
      zo1      = 0.
      zn2      = 0.
      zn4s     = 0.
      zo2_bro2 = 0.
      zn2_brn2 = 0.
!
! Summation over wavelength:
      do l=l1,lmax     ! 16,59 (44 wavelengths)
        sum1(:,:) = 0. ! sum(o2,o,n2)(sigma*chapman)
        sum2(:,:) = 0. ! sum(o2,o,n2)(sigma*psi/rmass)
        sum3(:,:) = 0. ! sum(o2,o,n2,n4s)(sigmas)
        do i=lon0,lon1
          do k=lev0,lev1
            sum1(k,i) = sum1(k,i)+sigeuv(1,l)*sco2(k,i)+
     |                            sigeuv(2,l)*sco1(k,i)+
     |                            sigeuv(3,l)*scn2(k,i)
            if (.not.check_exp) then
              sum1(k,i) = feuv(l)*exp(-sum1(k,i))
            else
              sum1(k,i) = feuv(l)*expo(-sum1(k,i),0)
            endif
            sum2(k,i) = sum2(k,i)+sigeuv(1,l)*o2i(k,i)*rmassinv_o2+
     |                            sigeuv(2,l)*o1i(k,i)*rmassinv_o1+
     |                            sigeuv(3,l)*n2i(k,i)*rmassinv_n2
            sum3(k,i) = sum3(k,i)+sigmas(1,l)+sigmas(2,l)+
     |                            sigmas(3,l)+sigmas(4,l)
          enddo
        enddo
!
! Longitude and column domain of current process:
        do i=lon0,lon1
          do k=lev0,lev1
            qtotal(k,i,lat) = qtotal(k,i,lat)+hc*rlmeuvinv(l)*
     |        sum1(k,i)*sum2(k,i)
            rj(k,i,lat) = rj(k,i,lat)+sum1(k,i)*sum3(k,i)
            qop2p(k,i,lat) = qop2p(k,i,lat)+sum1(k,i)*sigmas(6,l)*
     |        o1i(k,i)*rmassinv_o1
            qop2d(k,i,lat) = qop2d(k,i,lat)+sum1(k,i)*sigmas(5,l)*
     |        o1i(k,i)*rmassinv_o1
!
! Add ionization for o2,o,n2:
            zo2(k,i) = zo2(k,i)+sigmas(1,l)*sum1(k,i)*o2i(k,i)*
     |        rmassinv_o2
            zo1(k,i) = zo1(k,i)+sigmas(2,l)*sum1(k,i)*o1i(k,i)*
     |        rmassinv_o1
            zn2(k,i) = zn2(k,i)+sigmas(3,l)*sum1(k,i)*n2i(k,i)*
     |        rmassinv_n2
!
! Photo dissoc of n2:
            qtef(k,i,lat) = qtef(k,i,lat)+(sigeuv(3,l)-sigmas(3,l))*
     |        sum1(k,i)*n2i(k,i)*rmassinv_n2*2. ! *2 is a snoe mod
!
! Ionization for n4s, o2*bro2, and n2*brn2:
            zn4s(k,i) = zn4s(k,i)+sigmas(4,l)*sum1(k,i)*n4si(k,i)*
     |        rmassinv_n4s
            zo2_bro2(k,i) = zo2_bro2(k,i)+sigmas(1,l)*sum1(k,i)*o2i(k,i)
     |        *rmassinv_o2*bro2(l)
            zn2_brn2(k,i) = zn2_brn2(k,i)+sigmas(3,l)*sum1(k,i)*n2i(k,i)
     |        *rmassinv_n2*brn2(l)
          enddo
        enddo
      enddo ! l=l1,lmax
!
!     call addfsech('QOP2P' ,' ',' ', qop2p(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QOP2D' ,' ',' ', qop2d(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('SUM1' ,' ',' ', sum1,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('Q0'  ,' ',' ',qtotal(:,:,lat),lon0,lon1,
!    |  nlevs,nlevs,lat)
!     call addfsech('ZO2' ,' ',' ',zo2 ,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('ZO1' ,' ',' ',zo1 ,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('ZN2' ,' ',' ',zn2 ,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('ZN4S',' ',' ',zn4s,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('ZO2_BRO2',' ',' ',zo2_bro2,lon0,lon1,nlevs,nlevs,
!    |  lat)
!     call addfsech('ZN2_BRN2',' ',' ',zn2_brn2,lon0,lon1,nlevs,nlevs,
!    |  lat)
!     call addfsech('R_QRJ',' ',' ',r,lon0,lon1,nlevs,nlevs,lat)
!
! Multiply Q by efficiency factor:
      do i=lon0,lon1
        do k=lev0,lev1
          qtotal(k,i,lat) = qtotal(k,i,lat)*euveff(k)*avo
        enddo
      enddo
!     call addfsech('Q1'  ,' ',' ',qtotal(:,:,lat),lon0,lon1,
!    |  nlevs,nlevs,lat)
!
! Contributions to O2+, O+, N2+, N+, qtef:
      do i=lon0,lon1
        do k=lev0,lev1
          disn2p(k,i,lat) = zn2(k,i)*r(k,i)*xnmbari(k,i)
          qo2p(k,i,lat) = qo2p(k,i,lat)+(zo2(k,i)*(1.+rp(k,i))-
     |      zo2_bro2(k,i))*xnmbari(k,i)
          qn2p(k,i,lat) = qn2p(k,i,lat)+(zn2(k,i)*(1.+r(k,i))-
     |      zn2_brn2(k,i))*xnmbari(k,i) 
          qnp(k,i,lat)  = qnp(k,i,lat)+(zn4s(k,i)+zn2_brn2(k,i))*
     |      xnmbari(k,i)
          qtef(k,i,lat) = qtef(k,i,lat)*xnmbari(k,i)
          qop2pd(k,i) = qop2p(k,i,lat)+qop2d(k,i,lat)+zo1(k,i)+
     |      zo2_bro2(k,i)
          qop2p(k,i,lat) = (qop2p(k,i,lat)+qop2pd(k,i)*r(k,i)*0.22)*
     |      xnmbari(k,i)
          qop2d(k,i,lat) = (qop2d(k,i,lat)+qop2pd(k,i)*r(k,i)*0.24)*
     |      xnmbari(k,i)
          qop(k,i,lat) = qop(k,i,lat)+(zo1(k,i)+zo2_bro2(k,i)+
     |      qop2pd(k,i)*r(k,i)*0.56)*xnmbari(k,i)
        enddo
      enddo

!     call addfsech('QOP2PD',' ',' ',qop2pd,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QOP',' ',' ',qop(:,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
! 
! Add no ionization to qnop:
      do i=lon0,lon1
        qnop(lev0,i,lat) = qnop(lev0,i,lat)+
     |    beta9(lev0,i,lat)*no(lev0,i)*xnmbari(lev0,i)*rmassinv_no
      enddo
      do i=lon0,lon1
        do k=lev0+1,lev1
          qnop(k,i,lat) = qnop(k,i,lat)+beta9(k,i,lat)*.5*(no(k,i)+
     |      no(k-1,i))*xnmbari(k,i)*rmassinv_no
        enddo
      enddo
!
! tn at interfaces:
      do i=lon0,lon1
        tni(1,i) = tn(lev1,i) ! tn bottom boundary is stored in top slot
        do k=lev0+1,lev1-1
          tni(k,i) = .5*(tn(k-1,i)+tn(k,i))
        enddo
        tni(lev1,i) = tn(lev1-1,i) ! nlevp1 <- nlev
      enddo
!
! Quench:
      do i=lon0,lon1
        do k=lev0,lev1
          factor = avo*p0/gask*expz(1)*expzmid**(2*k-3)
          quenchfac(k,i) = factor/((o2i(k,i)*rmassinv_o2+
     |                              o1i(k,i)*rmassinv_o1+
     |                              n2i(k,i)*rmassinv_n2)*tni(k,i))
          quenchfac(k,i) = quenchfac(k,i)*
     |      (quench(1)*n2i(k,i)*rmassinv_n2+
     |       quench(2)*o2i(k,i)*rmassinv_o2)
          quenchfac(k,i) = quench(3)*quenchfac(k,i)/
     |                    (quench(4)+quenchfac(k,i))
        enddo
      enddo
!
! Summation over wave length:
      sum1(:,:) = 0.
      do l=1,l1-1 ! 1->15
        do i=lon0,lon1
          do k=lev0,lev1
!
! Note check_exp should be true for debug runs only:
            if (.not.check_exp) then
              sigchap(k,i) = sigsrc(l)*fsrc(l)*exp(-sigsrc(l)*sco2(k,i))
            else
              sigchap(k,i) = sigsrc(l)*fsrc(l)*
     |          expo(-sigsrc(l)*sco2(k,i),0)
            endif
!
            sum1(k,i) = sum1(k,i)+sigchap(k,i)*
     |        (hc*rlmsrcinv(l)-do22+quenchfac(k,i))
            rj(k,i,lat) = rj(k,i,lat)+sigchap(k,i)
          enddo
        enddo
      enddo
!
! Update q:
      do i=lon0,lon1
        do k=lev0,lev1
          qtotal(k,i,lat) = qtotal(k,i,lat)+sum1(k,i)*avo*
     |      o2i(k,i)*rmassinv_o2
        enddo
      enddo
!     call addfsech('Q2'  ,' ',' ',qtotal(:,:,lat),lon0,lon1,
!    |  nlevs,nlevs,lat)
!
! Contributions from schumann runge bands:
      do i=lon0,lon1
        do k=lev0,lev1
          if (sco2(k,i) >= 1.e18) then
            p3f(k,i) = (1./(aband*sco2(k,i)+bband*sqrt(sco2(k,i))))*
     |                 (1.+0.11*(f107-65.)/165.)*sfeps
          else
            p3f(k,i) = cband*(1.+0.11*(f107-65.)/165.)*sfeps
          endif
          qtotal(k,i,lat) = qtotal(k,i,lat)+p3f(k,i)*avo*
     |      o2i(k,i)*rmassinv_o2*e3
          rj(k,i,lat) = rj(k,i,lat)+p3f(k,i)/do2
        enddo
      enddo

!     call addfsech('Q'     ,' ',' ',qtotal(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('RJ'    ,' ',' ',rj(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QOP2P' ,' ',' ', qop2p(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QOP2D' ,' ',' ', qop2d(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QTEF'  ,' ',' ',  qtef(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QO2P'  ,' ',' ',  qo2p(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QN2P'  ,' ',' ',  qn2p(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QNP'   ,' ',' ',   qnp(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QOP'   ,' ',' ',   qop(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QNOP'  ,' ',' ',  qnop(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('DISN2P',' ',' ',disn2p(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!
#ifdef VT
!     code = 118 ; state = 'qrj' ; activity='ModelCode'
      call vtend(118,ier)
#endif
      end subroutine qrj
!---------------------------------------------------------------
      subroutine init_sflux
      use input_module,only: f107,f107a
      use init_module,only: sfeps
!
! Flux initialization once per time step, called from advance.
!
! Local:
      real :: euvflx(neuv),wave1(lmax),wave2(lmax)
      real :: flya,hlybr,fexvir,hlya,heiew,xuvfac
      integer :: iscale,n,nn
!
! Use f107,f107a from input_mod (these are either provided by
! the user, or obtained from GPI database, see input_module.F and
! gpi_module.F).
!
      hlybr = 0.
      fexvir = 0.
      hlya = 3.e+11+0.4E+10*(f107-70.)
      heiew = 0.
      xuvfac = 4.0 - (f107-68.0) / (243.0-68.0)
      if (xuvfac .lt. 1.0) xuvfac = 1.0
!
! ISCALE   =0 for Hinteregger contrast ratio method
!          =1 for Hinteregger linear interpolation
!          =2 for Tobiska EUV91 model
!          =3 for Woods & Rottman 10 Nov. 1988 measurement
!          =4 for Woods & Rottman 20 Jun. 1989 measurement
! F107     daily 10.7 cm flux (1.E-22 W m-2 Hz-1)
! F107A    81-day centered average 10.7 cm flux
! HLYBR    ratio of H Ly-b 1026A flux to solar minimum value (optional)
! FEXVIR   ratio of Fe XVI 335A flux to solar minimum value (optional)
! HLYA     H Lyman-alpha flux (photons cm-2 s-1) (optional)
! HEIEW    He I 10830A equivalent width, (milliAngstroms) (optional)
! XUVFAC   factor for scaling flux 16-250A (optional)
! WAVE1    longwave bound of spectral intervals (Angstroms)
! WAVE2    shortwave bound of intervals (= WAVE1 for indiv. lines)
! SFLUX    scaled solar flux returned by subroutine (photons cm-2 s-1)
!
      iscale = 0
      call ssflux(iscale,f107,f107a,hlybr,fexvir,hlya,
     |            heiew, xuvfac, wave1, wave2, sflux)
      call euvac(f107,f107a,euvflx)
!
! Transfer values of sflux to appropriate slots
!
      do n = l1,lmax
        feuv(n) = sflux(n)*sfeps
      enddo
      do n = 1,l1-1
        fsrc(n)   = sflux(n)*sfeps
      enddo
      do n=1,neuv
        nn = n+15
        feuv(nn) = euvflx(n)*sfeps
      enddo
!
!     write(6,"('  n   wave1   wave2       sflux      ',
!    |  'euvflx        fsrc        feuv')")
!     do n=1,lmax                   ! 1->59
!       if (n < l1) then            ! 1->15
!         write(6,"(i3,2f8.2,e12.4,12x,e12.4,e12.4)") 
!    |      n,wave1(n),wave2(n),sflux(n),fsrc(n),feuv(n)
!       elseif (n < l1+neuv-1) then ! 16->52
!         write(6,"(i3,2f8.2,e12.4,e12.4,12x,e12.4)") 
!    |      n,wave1(n),wave2(n),sflux(n),euvflx(n-l1+1),feuv(n)
!       else                        ! 53->59
!         write(6,"(i3,2f8.2,e12.4,12x,12x,e12.4)") 
!    |      n,wave1(n),wave2(n),sflux(n),feuv(n)
!       endif 
!     enddo
      end subroutine init_sflux
!-----------------------------------------------------------------------
      subroutine init_qrj
!
! Called once per run, from init_module.
!
! TOTAL ABSORPTION COEFFICIENTS
!   sigeuv(1,N) = SIGA(O2)
!   sigeuv(2,N) = SIGA(O)
!   sigeuv(3,N) = SIGA(N2)
!
! TOTAL IONIZATION COEFFICIENTS
!   SIGMAS(1,N) = SIGI(O2)
!   SIGMAS(2,N) = SIGI(O+(4S))
!   SIGMAS(3,N) = SIGI(N2)
!   SIGMAS(4,N) = SIGI(N)
!   SIGMAS(5,N) = SIGI(O+(2D))
!   SIGMAS(6,N) = SIGI(O+(2P))
!
      RLMEUV = (/0.17250E-04, 0.16750E-04, 0.16250E-04, 0.15750E-04,
     |           0.15250E-04, 0.14750E-04, 0.14250E-04, 0.13750E-04,
     |           0.13250E-04, 0.12750E-04, 0.12250E-04, 0.12157E-04,
     |           0.11750E-04, 0.11250E-04, 0.10750E-04, 0.10250E-04,
     |           0.10319E-04, 0.10257E-04, 0.97500E-05, 0.97702E-05,
     |           0.92500E-05, 0.87500E-05, 0.82500E-05, 0.77500E-05,
     |           0.78936E-05, 0.77041E-05, 0.76515E-05, 0.72500E-05,
     |           0.70331E-05, 0.67500E-05, 0.62500E-05, 0.62973E-05,
     |           0.60976E-05, 0.57500E-05, 0.58433E-05, 0.55437E-05,
     |           0.52500E-05, 0.47500E-05, 0.46522E-05, 0.42500E-05,
     |           0.37500E-05, 0.36807E-05, 0.32500E-05, 0.30378E-05,
     |           0.30331E-05, 0.27500E-05, 0.28415E-05, 0.25630E-05,
     |           0.22500E-05, 0.17500E-05, 0.12500E-05, 0.75000E-06,
     |           0.41000E-06, 0.27500E-06, 0.19500E-06, 0.12000E-06,
     |           0.60000E-07, 0.30000E-07, 0.15000E-07/)
      sigeuv(1,:) = (/
     |             0.50, 1.50, 3.40, 6.00,10.00,13.00,
     |            15.00,12.00, 2.20, 0.30, 3.00, 0.01,
     |             0.30, 0.10, 1.00, 1.10, 1.00, 1.60,
     |            16.53, 4.00,15.54, 9.85,20.87,27.09,
     |            26.66,25.18,21.96,29.05,25.00,26.27,
     |            26.02,25.80,26.10,25.04,22.00,25.59,
     |            24.06,21.59,20.40,19.39,18.17,18.40,
     |            17.19,16.80,16.80,15.10,15.70,13.20,
     |            10.60, 7.10, 4.00, 1.18, 0.32, 0.10,
     |             1.02, 0.14, .024, .004, .0004/)
      sigeuv(2,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,
     |             0.00, 0.00, 2.12, 4.18, 4.38, 4.23,
     |             4.28, 4.18, 4.18, 8.00,11.35,10.04,
     |            12.21,12.22,12.23,11.90,12.17,12.13,
     |            11.91,11.64,11.25,11.21, 9.64, 9.95,
     |             8.67, 7.70, 7.68, 6.61, 7.13, 6.05,
     |             5.30, 2.90, 1.60, 0.59, 0.16, 0.05,
     |             0.51, 0.07, .012, .002, .0002/)
      sigeuv(3,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,
     |            36.16, 0.70,16.99,46.63,15.05,30.71,
     |            19.26,26.88,35.46,30.94,26.30,29.75,
     |            23.22,23.20,23.10,22.38,23.20,24.69,
     |            24.53,21.85,21.80,21.07,17.51,18.00,
     |            13.00,11.60,11.60,10.30,10.60, 9.70,
     |             8.00, 4.40, 1.90, 0.60, 0.24, 1.16,
     |             0.48, 0.09, .015, .003, .0003/)
      SIGMAS(1,:) = (/
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.27, 0.00, 1.00,
     |            12.22, 2.50, 9.34, 4.69, 6.12, 9.39,
     |            11.05, 9.69, 8.59,23.81,23.00,22.05,
     |            25.94,25.80,26.10,25.04,22.00,25.59,
     |            24.06,21.59,20.40,19.39,18.17,18.40,
     |            17.19,16.80,16.80,15.10,15.70,13.20,
     |            10.60, 7.10, 4.00, 1.18, 0.32, 0.10,
     |             1.02, 0.14, .024, .004, .0004/)
      SIGMAS(2,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             2.12, 4.18, 4.38, 4.23,
     |             4.28, 4.18, 4.18, 4.20, 4.91, 4.01,
     |             3.78, 3.79, 3.67, 3.45, 3.53, 3.52,
     |             3.45, 3.26, 3.15, 3.03, 2.51, 2.59,
     |             2.25, 1.93, 1.92, 1.65, 1.78, 1.51,
     |             1.38, 0.78, 0.46, 0.18, 0.05, 0.015,
     |             0.015, 0.02, 0.004, 0.0006, 0.00006/)
      SIGMAS(3,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,
     |             0.00, 0.00, 0.00, 0.00, 0.00,16.75,
     |            10.18,18.39,23.77,23.20,23.00,25.06,
     |            23.22,23.20,23.10,22.38,23.20,24.69,
     |            24.53,21.85,21.80,21.07,17.51,18.00,
     |            13.00,11.60,11.60,10.30,10.60, 9.70,
     |             8.00, 4.40, 1.90, 0.60, 0.24, 1.16,
     |             0.48, 0.09, .015, .003, .0003/)
      SIGMAS(4,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             1.00, 1.00, 1.00, 1.00, 1.10, 1.10,
     |             1.10, 1.20, 1.20, 1.20, 1.20, 1.20,
     |             1.10, 1.20, 1.15, 1.10, 1.00, 1.00,
     |             1.00, 0.70, 0.80, 0.65, 0.60, 0.60,
     |             0.50, 0.50, 0.40, 0.35, 0.25, 0.20,
     |             0.10, 0.10, 0.10, 0.05, 0.01, 0.,0./)
      SIGMAS(5,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,
     |             3.80, 6.44, 5.52, 5.49,
     |             5.50, 5.50, 5.36, 5.48, 5.46, 5.36,
     |             5.24, 5.06, 4.71, 3.86, 3.98, 3.47,
     |             2.85, 2.84, 2.38, 2.64, 2.12, 1.86,
     |             0.99, 0.51, 0.19, 0.05, 0.015, 0.015,
     |             0.02, 0.004, 0.0006, 0.00006/)
      SIGMAS(6,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,
     |             0.50, 2.93, 2.93, 3.06,
     |             3.09, 3.16, 3.15, 3.10, 3.14, 3.04,
     |             3.48, 3.28, 3.38, 2.95, 2.93, 2.92,
     |             2.58, 2.71, 2.42, 2.07, 1.13, 0.62,
     |             0.22, 0.06, 0.02, 0.02, 0.03, 0.004,
     |             0.0007, 0.00007/)
      BRN2 = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.01, 0.04, 0.04, 0.03, 0.05, 0.05,
     |             0.15, 0.20, 0.20, 0.25, 0.32, 0.34,
     |             0.36, 0.36, 0.36, 0.36, 0.36, 0.36,
     |             0.36, 0.36, 0.36, 0.36, 0.36/)
      BRO2 = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,
     |             .025, .036, .045, .120, .155, .189,
     |             .230, 0.20, 0.20, 0.20, 0.23, 0.25,
     |             0.29, 0.33, 0.33, 0.33, 0.33, 0.33,
     |             0.33, 0.33, 0.33, 0.33, 0.33, 0.33,
     |             0.33, 0.33, 0.33, 0.33, 0.33/)
      euveff(:) = 0.05
      quench = (/7.E-11,5.E-11,3.1401E-12,9.1E-3/)
      end subroutine init_qrj
!---------------------------------------------------------------
      subroutine init_sigmas
!
! Called once per run from init_module.
!
      integer :: m,n,nn
      do m=1,3
        do n=1,lmax
          sigeuv(m,n)=sigeuv(m,n)*1.E-18
          sigmas(m,n)=sigmas(m,n)*1.E-18
        enddo
      enddo
      do n = 1,lmax
        sigmas(4,n) = sigmas(4,n)*1.E-17
        sigmas(5,n) = sigmas(5,n)*1.E-18
        sigmas(6,n) = sigmas(6,n)*1.E-18
      enddo
      do n = 1,l1-1
        rlmsrc(n) = rlmeuv(n)
        sigsrc(n) = sigeuv(1,n)
      enddo
      do n=1,neuv
        brop2p(n) = 0.
        if (n > 14) brop2p(n) = 1.-brop2d(n)-brop4s(n)
        sigop2p(n)=sigio(n)*brop2p(n)
        sigop2d(n)=sigio(n)*brop2d(n)
        sigop4s(n)=sigio(n)*brop4s(n)
      enddo
      do n = 1,neuv
        nn = n+15    ! 16:52
        sigeuv(1,nn) = sigao2(n)
        sigeuv(2,nn) = sigao(n)
        sigeuv(3,nn) = sigan2(n)
        sigmas(1,nn) = sigio2(n)
        sigmas(2,nn) = sigop4s(n)
        sigmas(3,nn) = sigin2(n)
        sigmas(4,nn) = sigin(n)
        sigmas(5,nn) = sigop2d(n)
        sigmas(6,nn) = sigop2p(n)
        brn2(nn) = brn2np(n)
        bro2(nn) = bro2op(n)
      enddo
      end subroutine init_sigmas
!-----------------------------------------------------------------------
      subroutine init_euvac
!
! Called once per run from init_module.
!
! lambdas:
      wleuv1 = 
     +    (/1000.00, 1031.91, 1025.72,  950.00,  977.02,  900.00,
     +       850.00,  800.00,  750.00,  789.36,  770.41,  765.15,
     +       700.00,  703.36,  650.00,  600.00,  629.73,  609.76,
     +       550.00,  584.33,  554.31,  500.00,  450.00,  465.22,
     +       400.00,  350.00,  368.07,  300.00,  303.78,  303.31,
     +       250.00,  284.15,  256.30,  200.00,  150.00,  100.00,
     +        50.00/)
      wleuv2 =
     +    (/1050.00, 1031.91, 1025.72, 1000.00,  977.02,  950.00,
     +       900.00,  850.00,  800.00,  789.36,  770.41,  765.15,
     +       750.00,  703.36,  700.00,  650.00,  629.73,  609.76,
     +       600.00,  584.33,  554.31,  550.00,  500.00,  465.22,
     +       450.00,  400.00,  368.07,  350.00,  303.78,  303.31,
     +       300.00,  284.15,  256.30,  250.00,  200.00,  150.00,
     +       100.00/)
c
c sigao
      sigao =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 1.32e-18,
     +      4.55e-18, 3.50e-18, 5.09e-18, 3.75e-18, 3.89e-18, 4.00e-18,
     +      1.07e-17, 1.15e-17, 1.72e-17, 1.34e-17, 1.34e-17, 1.34e-17,
     +      1.30e-17, 1.31e-17, 1.26e-17, 1.21e-17, 1.21e-17, 1.19e-17,
     +      1.15e-17, 9.69e-18, 9.84e-18, 8.69e-18, 7.70e-18, 7.68e-18,
     +      6.46e-18, 7.08e-18, 6.05e-18, 5.20e-18, 3.73e-18, 1.84e-18,
     +      7.30e-19/)
c
c sigao2
      sigao2 =
     +    (/1.35e-18, 1.00e-18, 1.63e-18, 2.11e-17, 1.87e-17, 1.28e-17,
     +      8.56e-18, 1.66e-17, 2.21e-17, 2.67e-17, 1.89e-17, 2.08e-17,
     +      2.85e-17, 2.74e-17, 2.19e-17, 2.60e-17, 3.21e-17, 2.81e-17,
     +      2.66e-17, 2.28e-17, 2.60e-17, 2.46e-17, 2.31e-17, 2.19e-17,
     +      2.03e-17, 1.81e-17, 1.83e-17, 1.74e-17, 1.68e-17, 1.68e-17,
     +      1.44e-17, 1.58e-17, 1.34e-17, 1.09e-17, 7.51e-18, 3.81e-18,
     +      1.32e-18/)
c
c sigan2
      sigan2 =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 5.10e-17, 2.24e-18, 9.68e-18,
     +      2.02e-17, 1.70e-17, 3.36e-17, 1.65e-17, 1.42e-17, 1.20e-16,
     +      2.47e-17, 2.65e-17, 3.18e-17, 2.33e-17, 2.34e-17, 2.28e-17,
     +      2.28e-17, 2.24e-17, 2.41e-17, 2.45e-17, 2.35e-17, 2.32e-17,
     +      2.17e-17, 1.64e-17, 1.69e-17, 1.39e-17, 1.17e-17, 1.17e-17,
     +      1.05e-17, 1.09e-17, 1.02e-17, 8.39e-18, 4.96e-18, 2.26e-18,
     +      7.20e-19/)
c
c sigio
      sigio =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 1.32e-18,
     +      4.55e-18, 3.50e-18, 5.09e-18, 3.75e-18, 3.89e-18, 4.00e-18,
     +      1.07e-17, 1.15e-17, 1.72e-17, 1.34e-17, 1.34e-17, 1.34e-17,
     +      1.30e-17, 1.31e-17, 1.26e-17, 1.21e-17, 1.21e-17, 1.19e-17,
     +      1.15e-17, 9.69e-18, 9.84e-18, 8.69e-18, 7.70e-18, 7.68e-18,
     +      6.46e-18, 7.08e-18, 6.05e-18, 5.20e-18, 3.73e-18, 1.84e-18,
     +      7.30e-19/)
c
c sigio2
      sigio2 =
     +    (/2.59e-19, 0.00e+00, 1.05e-18, 1.39e-17, 1.55e-17, 9.37e-18,
     +      5.49e-18, 6.41e-18, 1.06e-17, 1.02e-17, 8.47e-18, 1.17e-17,
     +      2.38e-17, 2.38e-17, 2.13e-17, 2.49e-17, 3.11e-17, 2.64e-17,
     +      2.66e-17, 2.28e-17, 2.60e-17, 2.46e-17, 2.31e-17, 2.19e-17,
     +      2.03e-17, 1.81e-17, 1.83e-17, 1.74e-17, 1.68e-17, 1.68e-17,
     +      1.44e-17, 1.58e-17, 1.34e-17, 1.09e-17, 7.51e-18, 3.81e-18,
     +      1.32e-18/)
c
c sigin2
      sigin2 =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 1.43e-17, 8.86e-18, 8.50e-18, 6.58e-17,
     +      1.51e-17, 2.55e-17, 2.92e-17, 2.33e-17, 2.34e-17, 2.28e-17,
     +      2.28e-17, 2.24e-17, 2.41e-17, 2.45e-17, 2.35e-17, 2.32e-17,
     +      2.17e-17, 1.64e-17, 1.69e-17, 1.39e-17, 1.17e-17, 1.17e-17,
     +      1.05e-17, 1.09e-17, 1.02e-17, 8.39e-18, 4.96e-18, 2.26e-18,
     +      7.20e-19/)
c
c sigin
      sigin =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.21e-18,10.29e-18,11.71e-18,10.96e-18,11.24e-18,11.32e-18,
     +     12.10e-18,13.26e-18,12.42e-18,11.95e-18,11.21e-18,11.80e-18,
     +     11.76e-18,11.78e-18,11.77e-18,11.50e-18,11.02e-18,10.58e-18,
     +      9.56e-18, 8.15e-18, 8.30e-18, 7.30e-18, 6.41e-18, 6.40e-18,
     +      5.24e-18, 5.73e-18, 4.87e-18, 3.95e-18, 2.49e-18, 0.99e-18,
     +      0.33e-18/)
c
c O 4S
      brop4s =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 1.00e+00,
     +      1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     +      6.30e-01, 4.30e-01, 3.00e-01, 3.20e-01, 2.90e-01, 2.70e-01,
     +      2.80e-01, 3.00e-01, 2.90e-01, 2.80e-01, 2.80e-01, 2.80e-01,
     +      2.70e-01, 2.60e-01, 2.60e-01, 2.60e-01, 2.50e-01, 2.50e-01,
     +      2.50e-01, 2.50e-01, 2.50e-01, 2.60e-01, 2.70e-01, 2.90e-01,
     +      3.00e-01/)
c O 2D
      brop2d =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      3.70e-01, 5.70e-01, 6.60e-01, 4.60e-01, 4.70e-01, 4.60e-01,
     +      4.50e-01, 4.60e-01, 4.50e-01, 4.50e-01, 4.50e-01, 4.50e-01,
     +      4.30e-01, 4.00e-01, 4.00e-01, 4.00e-01, 3.70e-01, 3.70e-01,
     +      3.60e-01, 3.70e-01, 3.50e-01, 3.50e-01, 3.30e-01, 3.20e-01,
     +      3.20e-01/)
c O 2P
c     brop2p =
c    +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
c    +      0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
c    +      0.00e+00, 0.00e+00, 4.00e-02, 2.20e-01, 2.40e-01, 2.70e-01,
c    +      2.70e-01, 2.40e-01, 2.60e-01, 2.70e-01, 2.70e-01, 2.70e-01,
c    +      2.60e-01, 2.50e-01, 2.60e-01, 2.50e-01, 2.50e-01, 2.50e-01,
c    +      2.30e-01, 2.30e-01, 2.30e-01, 2.20e-01, 2.20e-01, 2.10e-01,
c    +      2.10e-01/)
c N2 -> N+
      brn2np =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 0.00e+00, 1.30e-03, 3.00e-02, 4.60e-02,
     +      4.50e-02, 1.05e-01, 9.00e-02, 1.63e-01, 2.13e-01, 2.13e-01,
     +      3.00e-01, 2.57e-01, 3.35e-01, 3.77e-01, 3.64e-01, 3.46e-01,
     +      3.85e-01/)
c O2 -> O+
      bro2op =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 8.91e-03, 3.86e-02, 3.31e-02, 7.01e-02,
     +      1.44e-01, 1.71e-01, 1.71e-01, 2.08e-01, 2.31e-01, 2.25e-01,
     +      2.37e-01, 2.59e-01, 2.49e-01, 2.99e-01, 3.53e-01, 3.53e-01,
     +      3.71e-01, 3.73e-01, 3.66e-01, 3.93e-01, 4.48e-01, 3.84e-01,
     +      0.00e+00/)
      end subroutine init_euvac
!-----------------------------------------------------------------------
      subroutine euvac(F107,F107A,EUVFLX)
!
! This EUV flux model uses the F74113 solar reference spectrum and
! ratios determined from Hinteregger's SERF1 model. It uses the daily
! F10.7 flux (F107) and the 81 day mean (F107A) as a proxy for scaling
! The fluxes are returned in EUVFLX and correspond to the 37 wavelength
! bins of Torr et al. [1979] Geophys. Res. Lett. p771.
! See Richards et al. [1994] J. Geophys. Res. p8981 for details.
!
! Called from init_sflux, once per timestep.
!
! F107   = input daily 10.7 cm flux index. (e.g. 74)
! F107A  = input 81 day average of daily F10.7 centered on current day
! EUVFLX = output array for EUV flux in units of photons/cm2/sec.
!
      real :: wleuv1,wleuv2,sigao,sigao2,sigan2,sigio,sigio2,
     |  sigin2,brop4s,brop2d,brop2p,sigop2p,sigop2d,sigop4s,
     |  sigin,brn2np,bro2op
!
! Args:
      real,intent(in) :: f107,f107a
      real,intent(out) :: euvflx(neuv)
!
! Local:
      integer :: i
      real :: AFAC(neuv),F74113(neuv),flxfac
      real :: afac_tmp(neuv),f74113_tmp(neuv)
!
! F74113 reference spectrum (doubled below 150-250 A, tripled <150)
! Will be multiplied by 1.0E9 later
!
      F74113=
     | (/1.200,0.450,4.800,3.100,0.460,0.210,1.679,0.800,6.900,
     |   0.965,0.650,0.314,0.383,0.290,0.285,0.452,0.720,1.270,
     |   0.357,0.530,1.590,0.342,0.230,0.360,0.141,0.170,0.260,
     |   0.702,0.758,1.625,3.537,3.000,4.400,1.475,3.500,2.100,
     |   2.467/)
!
! Scaling factors(Ai) for the EUV flux
      AFAC=
     | (/1.0017E-02,7.1250E-03,1.3375E-02,1.9450E-02,2.7750E-03,
     |   1.3768E-01,2.6467E-02,2.5000E-02,3.3333E-03,2.2450E-02,
     |   6.5917E-03,3.6542E-02,7.4083E-03,7.4917E-03,2.0225E-02,
     |   8.7583E-03,3.2667E-03,5.1583E-03,3.6583E-03,1.6175E-02,
     |   3.3250E-03,1.1800E-02,4.2667E-03,3.0417E-03,4.7500E-03,
     |   3.8500E-03,1.2808E-02,3.2750E-03,4.7667E-03,4.8167E-03,
     |   5.6750E-03,4.9833E-03,3.9417E-03,4.4167E-03,5.1833E-03,
     |   5.2833E-03,4.3750E-03/)
!
! 6/20/03 bf: Reverse order of F74113 and AFAC to correspond to descending
! wavelength bins (as in wleuv1,2 and WAVEL,WAVES):
      do i=1,neuv
        afac_tmp(i) = afac(neuv-i+1)  
        f74113_tmp(i) = f74113(neuv-i+1)  
      enddo
      afac = afac_tmp      ! whole array
      f74113 = f74113_tmp
!
! Loop through the wavelengths calculating the scaling factors and
!   the resulting solar flux.
! The scaling factors are restricted to be greater than 0.8
!
      DO 50 I=1,neuv
        FLXFAC=(1.0 + AFAC(I) * (0.5*(F107+F107A) - 80.0))
        IF(FLXFAC.LT.0.8) FLXFAC=0.8
        EUVFLX(I)=F74113(I) * FLXFAC * 1.0E9
 50   CONTINUE
      end subroutine euvac
!-----------------------------------------------------------------------
      subroutine ssflux (iscale, f107, f107a, hlybr, fexvir, hlya,
     |                   heiew, xuvfac, wave1, wave2, sflux)
      integer,parameter :: lmax=59
!
! Args:
      integer,intent(in) :: iscale
      real,intent(in) :: f107,f107a,hlybr,fexvir,hlya,heiew,xuvfac
      real,intent(out) :: wave1(lmax),wave2(lmax),sflux(lmax) 
!
! Local:
      real ::
     |          wavel(lmax), waves(lmax), rflux(lmax), xflux(lmax),
     |          scale1(lmax), scale2(lmax),
     |          tchr0(lmax), tchr1(lmax), tchr2(lmax),
     |          tcor0(lmax), tcor1(lmax), tcor2(lmax),
     |          war1(lmax), war2(lmax),
     |          b1(3), b2(3)
      real :: frat,r1,r2,hlymod,heimod,xuvf
      integer :: l
!
! regression coefficients which reduce to solar min. spectrum:
      b1 = (/1.0, 0.0138, 0.005/) 
      b2 = (/1.0, 0.59425, 0.3811/)
!
! 'best fit' regression coefficients, commented out, for reference:
!     b1 = (/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/)
!
      WAVEL =   (/1750.00, 1700.00, 1650.00, 1600.00, 1550.00, 1500.00,
     |            1450.00, 1400.00, 1350.00, 1300.00, 1250.00, 1215.67,
     |            1200.00, 1150.00, 1100.00, 1050.00, 1031.91, 1025.72,
     |            1000.00,  977.02,  950.00,  900.00,  850.00,  800.00,
     |             789.36,  770.41,  765.15,  750.00,  703.31,  700.00,
     |             650.00,  629.73,  609.76,  600.00,  584.33,  554.37,
     |             550.00,  500.00,  465.22,  450.00,  400.00,  368.07,
     |             350.00,  303.78,  303.31,  300.00,  284.15,  256.30,
     |             250.00,  200.00,  150.00,  100.00,   50.00,   32.00,
     |              23.00,   16.00,    8.00,    4.00,    2.00/)
      WAVES =   (/1700.00, 1650.00, 1600.00, 1550.00, 1500.00, 1450.00,
     |            1400.00, 1350.00, 1300.00, 1250.00, 1200.00, 1215.67,
     |            1150.00, 1100.00, 1050.00, 1000.00, 1031.91, 1025.72,
     |             950.00,  977.02,  900.00,  850.00,  800.00,  750.00,
     |             789.36,  770.41,  765.15,  700.00,  703.31,  650.00,
     |             600.00,  629.73,  609.76,  550.00,  584.33,  554.37,
     |             500.00,  450.00,  465.22,  400.00,  350.00,  368.07,
     |             300.00,  303.78,  303.31,  250.00,  284.15,  256.30,
     |             200.00,  150.00,  100.00,   50.00,   32.00,   23.00,
     |              16.00,    8.00,    4.00,    2.00,    1.00/)
      RFLUX =    (/322.00,  168.00,   95.00,   62.00,   44.00,   25.00,
     |              16.90,   11.80,   19.50,    4.10,   11.10,  249.00,
     |               2.78,    0.70,    3.07,    3.64,    3.18,    4.38,
     |               1.78,    5.96,    4.22,    4.43,    1.93,    0.87,
     |               0.79,    0.24,    0.20,    0.17,    0.39,    0.22,
     |               0.17,    1.50,    0.45,    0.48,    1.58,    0.80,
     |               0.51,    0.31,    0.18,    0.39,    0.21,    0.74,
     |               0.87,    6.00,    0.24,    0.84,    0.10,    0.27,
     |               0.92,    1.84,    0.13,    0.38,  0.0215,  0.0067,
     |              1.E-3,   2.E-3,   1.E-5,   5.E-8,   1.E-10/)
      XFLUX =    (/354.00,  191.00,  110.00,   76.00,   55.00,   28.00,
     |              19.60,   14.30,   25.30,    5.00,   17.20,  401.00,
     |               6.26,    1.51,    6.11,    8.66,    9.04,   13.12,
     |               4.42,   13.18,   12.03,   13.29,    5.01,    2.18,
     |               1.59,    0.67,    0.43,    0.43,    0.72,    0.46,
     |               0.48,    3.02,    1.46,    1.02,    4.86,    1.59,
     |               1.57,    1.67,    0.36,    0.99,    2.20,    1.39,
     |               5.63,   11.28,    2.50,    4.14,    3.16,    0.59,
     |               3.70,    4.85,    0.34,    1.15,    0.18,    0.08,
     |             2.5E-2,   5.E-2,   8.E-4,   3.E-5,   5.E-7/)
      SCALE1 = (/     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     |                0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     |            1692.09,  405.95, 1516.20, 2731.70, 3314.57, 4375.00,
     |            1316.91, 3621.91, 3908.56, 4432.54, 1541.21,  531.73,
     |             364.83,    0.00,  116.00,  129.41,  162.48,   94.07,
     |              41.29,  709.50,    0.00,  268.47, 1561.05,  367.64,
     |             290.06,  184.36,    0.00,   86.15,    7.50,    0.00,
     |               0.00, 2220.00,    0.00,   61.00,    0.00,   86.95,
     |             206.00,  135.89,   60.35,  157.12,    7.06,    0.75,
     |               0.00,    0.00,    0.00,    0.00,    0.00/)
      SCALE2 =  (/   0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    5.34,    0.00,    0.00,    0.00,    0.54,
     |               3.30,    0.00,   12.60,    0.00,    0.00,    0.00,
     |               5.34,   11.63,    2.28,    5.56,   24.93,    8.16,
     |              60.69,    0.00,   28.20,   45.90,   40.80,    1.27,
     |              35.47,   42.80,    1.12,    6.19,    1.26,    0.69,
     |               0.23,    0.46,  7.6E-3,  2.9E-4,  4.8E-6/)
      TCHR0 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,-4.290E+00,-5.709E+00,-8.493E+00,
     |-1.161E+00,-3.429E+00,-5.464E+00,-6.502E+00,-1.912E+00,-4.034E-01,
     |-1.448E-01, 0.000E+00,-9.702E-02,-6.591E-02,-2.338E-02,-1.273E-01,
     |-2.406E-01,-3.351E-01, 0.000E+00,-1.465E+00,-2.405E+00,-7.975E-02,
     |-4.197E-01,-1.971E-01, 0.000E+00,-5.895E-02,-5.815E-03, 0.000E+00,
     | 0.000E+00, 2.138E-01, 0.000E+00,-7.713E-02, 0.000E+00,-3.035E-02,
     |-2.039E-01,-1.749E-01,-1.041E-01,-2.638E-01,-1.094E-02, 0.000E+00,
     |       0.0,       0.0,       0.0,       0.0,       0.0/)
      TCHR1 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,-3.023E-13,-3.745E-13,-5.385E-13,
     |-1.211E-13,-3.868E-13,-3.646E-13,-4.125E-13,-1.527E-13,-4.753E-14,
     |-3.411E-14, 0.000E+00,-1.190E-14,-1.034E-14,-1.343E-14,-1.539E-14,
     |-5.174E-14,-6.934E-14, 0.000E+00,-1.215E-13,-1.537E-13,-2.024E-14,
     |-4.596E-14,-1.562E-14, 0.000E+00,-1.221E-14,-1.123E-15, 0.000E+00,
     | 0.000E+00,-2.263E-13, 0.000E+00,-1.508E-14, 0.000E+00,-1.744E-14,
     |-2.100E-14,-1.805E-14,-8.224E-15,-1.919E-14,-7.944E-16, 0.000E+00,
     |       0.0,       0.0,       0.0,       0.0,       0.0/)
      TCHR2 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0, 3.275E-11, 4.057E-11, 6.160E-11,
     | 1.312E-11, 4.189E-11, 4.167E-11, 4.716E-11, 1.654E-11, 5.150E-12,
     | 3.901E-12, 0.000E+00, 1.289E-12, 1.120E-12, 1.455E-12, 1.667E-12,
     | 5.604E-12, 7.931E-12, 0.000E+00, 1.317E-11, 1.757E-11, 2.194E-12,
     | 4.978E-12, 1.693E-12, 0.000E+00, 1.324E-12, 1.285E-13, 0.000E+00,
     | 0.000E+00, 2.586E-11, 0.000E+00, 1.724E-12, 0.000E+00, 1.889E-12,
     | 2.400E-12, 2.063E-12, 8.911E-13, 2.193E-12, 9.090E-14, 0.000E+00,
     |       0.0,       0.0,       0.0,       0.0,       0.0/)
      TCOR0 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0, 0.000E+00, 0.000E+00, 0.000E+00,
     | 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,-6.060E-02,
     | 0.000E+00,-3.399E-02, 0.000E+00, 0.000E+00, 0.000E+00, 4.866E-02,
     |-1.762E-01, 0.000E+00,-2.412E-01, 0.000E+00, 0.000E+00, 0.000E+00,
     |-4.743E-01,-9.713E-01, 5.891E-02,-1.263E-01,-1.246E+00, 2.870E-01,
     |-4.659E+00, 0.000E+00,-1.058E+00,-3.821E+00,-1.874E+00, 0.000E+00,
     |-1.896E+00,-8.505E-01,-2.101E-04,-2.012E-01,-6.097E-02,-2.925E-02,
     |-4.875E-03,       0.0,       0.0,       0.0,       0.0/)
      TCOR1 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0, 0.000E+00, 0.000E+00, 0.000E+00,
     | 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.877E-03,
     | 0.000E+00, 1.760E-03, 0.000E+00, 0.000E+00, 0.000E+00, 3.313E-04,
     | 3.643E-03, 0.000E+00, 5.225E-03, 0.000E+00, 0.000E+00, 0.000E+00,
     | 4.085E-03, 1.088E-02, 8.447E-04, 3.237E-03, 1.907E-02, 2.796E-03,
     | 4.460E-02, 0.000E+00, 1.007E-02, 3.481E-02, 1.604E-02, 0.000E+00,
     | 2.029E-02, 2.160E-02, 6.342E-04, 3.594E-03, 5.503E-04, 2.687E-04,
     | 4.479E-05,       0.0,       0.0,       0.0,       0.0/)
      TCOR2 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0, 0.000E+00, 0.000E+00, 0.000E+00,
     | 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.846E-03,
     | 0.000E+00, 1.127E-03, 0.000E+00, 0.000E+00, 0.000E+00, 1.891E-04,
     | 2.326E-03, 0.000E+00, 2.801E-03, 0.000E+00, 0.000E+00, 0.000E+00,
     | 2.446E-03, 7.121E-03, 5.204E-04, 1.983E-03, 1.204E-02, 1.721E-03,
     | 2.911E-02, 0.000E+00, 7.177E-03, 2.272E-02, 9.436E-03, 0.000E+00,
     | 1.316E-02, 1.398E-02, 4.098E-04, 2.328E-03, 3.574E-04, 1.745E-04,
     | 2.909E-05,       0.0,       0.0,       0.0,       0.0/)
      WAR1 =  (/     0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    3.80,    6.25,    4.93,    6.06,
     |               2.70,    7.07,    8.62,    9.60,    4.54,    2.37,
     |               0.82,    0.33,    0.24,    0.67,    0.28,    0.55,
     |               1.56,    1.11,    0.77,    1.32,    1.71,    0.44,
     |               1.11,    0.95,    0.39,    0.81,    2.00,    1.49,
     |               6.81,    5.07,    1.63,    5.62,    2.08,    0.59,
     |               3.89,    5.19,    0.35,    1.18,    0.099,   0.04,
     |               0.007,   0.00,    0.00,    0.00,    0.00/)
      WAR2 =  (/     0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,   20.80,   17.90,    9.30,   14.30,
     |               6.90,   12.00,   15.60,   18.60,   10.10,    4.30,
     |              12.40,    8.00,    3.60,    1.80,    0.50,    1.40,
     |               3.90,    2.60,    1.60,    3.40,    4.10,    0.70,
     |               4.30,    4.30,    3.80,    2.60,    6.08,    1.35,
     |              12.60,    9.78,    2.96,   10.20,    4.11,    6.68,
     |               6.62,    8.07,    0.47,    1.73,    0.17,    0.075,
     |               0.012,   0.00,    0.00,    0.00,    0.00/)
C
C Linear Interpolation between SC#21REFW and F79050:
C
      FRAT = (F107-68.) / (243.-68.)
      DO 200 L=1,LMAX
        SFLUX(L) = RFLUX(L) + (XFLUX(L)-RFLUX(L)) * FRAT
  200 CONTINUE
C
C Hinteregger contrast ratio method:
C
      IF (ISCALE .EQ. 0) THEN
        IF (HLYBR .GT. 0.001) THEN
          R1 = HLYBR
        ELSE
          R1 =  B1(1) + B1(2)*(F107A-71.5) + B1(3)*(F107-F107A+3.9)
        ENDIF
        IF (FEXVIR .GT. 0.001) THEN
          R2 = FEXVIR
        ELSE
          R2 =  B2(1) + B2(2)*(F107A-71.5) + B2(3)*(F107-F107A+3.9)
        ENDIF
        DO 100 L=13,LMAX
          SFLUX(L) = (RFLUX(L) + ((R1-1.)*SCALE1(L)
     |                         +  (R2-1.)*SCALE2(L)) / 1000.)
  100   CONTINUE
      ENDIF
C
C Tobiska EUV91 Method:
C
      IF (ISCALE .EQ. 2) THEN
        IF (HLYA .GT. 0.001) THEN
          HLYMOD = HLYA
        ELSE
          IF (HEIEW .GT. 0.001) THEN
             HLYMOD = HEIEW * 3.77847E9 + 8.40317E10
          ELSE
             HLYMOD = 8.70E8 * F107 + 1.90E11
             HLYMOD = 8.70E8 * F107 + 1.90E11
          ENDIF
        ENDIF
        IF (HEIEW .GT. 0.001) THEN
          HEIMOD = HEIEW * 3.77847E9 + 8.40317E10
        ELSE
          HEIMOD = HLYMOD
        ENDIF
        DO 500 L=16,55
          SFLUX(L) = TCHR0(L) + TCHR1(L)*HLYMOD + TCHR2(L)*HEIMOD
     |             + TCOR0(L) + TCOR1(L)*F107 + TCOR2(L)*F107A
 500    CONTINUE
      ENDIF
C
C Woods and Rottman (10 Nov. 1988) spectrum:
C
      IF (ISCALE .EQ. 3) THEN
        DO 550 L=15,55
        SFLUX(L) = WAR1(L)
  550   CONTINUE
      ENDIF
C
C Woods and Rottman (20 June 1989) spectrum:
C
      IF (ISCALE .EQ. 4) THEN
        DO 560 L=15,55
        SFLUX(L) = WAR2(L)
  560   CONTINUE
      ENDIF
C
C Substitute in H Lyman-alpha and XUVFAC if provided:
C
      IF (HLYA .GT. 0.001) SFLUX(12) = HLYA / 1.E9
      IF (XUVFAC .GT. 0.001) THEN
        XUVF = XUVFAC
      ELSE
        XUVF = 1.0
      ENDIF
C
C Convert from gigaphotons to photons, etc.:
C
      DO 600 L=1,LMAX
        WAVE1(L) = WAVEL(L)
        WAVE2(L) = WAVES(L)
        IF (SFLUX(L) .LT. 0.0) SFLUX(L) = 0.0
        IF (WAVEL(L).LT.251.0 .AND. WAVES(L).GT.15.0)
     |     SFLUX(L)=SFLUX(L)*XUVF
        SFLUX(L) = SFLUX(L) * 1.E9
  600 CONTINUE
C
      end subroutine ssflux
!-----------------------------------------------------------------------
      subroutine alloc_q(lon0,lon1,lat0,lat1)
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
      allocate(rj(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' rj: stat=',i3)") istat
      allocate(qtef(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qtef: stat=',i3)") istat
      allocate(qtotal(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qtotal: stat=',i3)") istat
      allocate(qop2p(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qop2p: stat=',i3)") istat
      allocate(qop2d(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qop2d: stat=',i3)") istat
      allocate(qo2p(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qo2p: stat=',i3)") istat
      allocate(qop(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qop: stat=',i3)") istat
      allocate(qn2p(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qn2p: stat=',i3)") istat
      allocate(qnp(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qnp: stat=',i3)") istat
      allocate(qnop(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qnop: stat=',i3)") istat
      end subroutine alloc_q
!-----------------------------------------------------------------------
      end module qrj_module
