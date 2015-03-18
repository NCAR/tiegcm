!
      module qrj_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlevp1,nlonp4
      use addfld_module,only: addfld
      implicit none
!
      integer,parameter :: lmax=37  ! number of bins from 0.5A-1750A
      integer,parameter :: l1=15    ! number of bins in (1050-1750A)
      real :: 
     |  euveff(nlevp1),
     |  sigeuv(3,lmax),
     |  rlmeuv(lmax),
     |  feuv(lmax),
     |  fsrc(l1),
     |  sigsrc(l1),
     |  rlmsrc(l1),
     |  sigin4s(lmax),
     |  quench(4),
     |  wave1(lmax),     !       short bound of wave bins
     |  wave2(lmax),     !       long bound of wave bins 
     |  sfmin(lmax),     !       reference solar minimum flux of EUVAC model
     |  afac(lmax),      !       The A factor of EUVAC model
     |  sflux(lmax)      !       Solar flux for each time step
!
! Branching ratios for photon are branching ratios from photon absorption rate.
! Branching ratios for photolectron are branching ratios from photoionization rate.

      real::
     |  BPhotonI(3,lmax),     ! Photoionization branching ratio for major species
     |  BElectronI(3,lmax),   ! Photoelectron ionization branching ratio
                              ! for major species
     |  brop2pPh(lmax),   ! photoionization branching ratio for O+(2p)
     |  brop2dPh(lmax),   ! photoionization branching ratio for O+(2d)
     |  brop4sPh(lmax),   ! photoionization branching ratio for O+(4s)
     |  bro2DPh(lmax),    ! photodissociation branching ratio for O2
     |  brn2DPh(lmax),    ! photodissociation branching ratio for N2
     |  bro2DIPh(lmax),   ! Photon dissociative ionization branching ratio for O2
     |  brn2DIPh(lmax),   ! Photon dissociative ionization branching ratio for N2
     |  brop2pEl(lmax),   ! electron impact ionization branching ratio for O+(2p) 
     |  brop2dEl(lmax),   ! electron impact ionization branching ratio for O+(2d)
     |  brop4sEl(lmax),   ! electron impact ionization branching ratio for O+(4s)
     |  bro2DIEl(lmax),   ! Photoelectron dissociative ionization branching ratio for O2
     |  brn2DIEl(lmax),   ! Photoelectron dissociative ionization branching ratio for N2
     |  brn2DEl(lmax),    ! Photoelectron dissociation branching ratio for N2
     |  bro2DEl(lmax)     ! Photoelectron dissociation branching ratio for O2
!
! Heating and ionization terms set by qrj, and used by other routines.
! These are allocated for task subdomains by alloc_q (called from allocdata)
!
      real,dimension(:,:,:),allocatable ::
     |  rj,     ! total o2 dissociation frequency (s^-1)
     |  qtef,   ! total N dissociative production rate (cm^-3 s^-1)
     |  qtotal, ! total heating rate
     |  qop2p,  ! o+(2p) production rate
     |  qop2d,  ! o+(2d) production rate
     |  qo2p,   ! o2+ production rate
     |  qop,    ! o+ production rate
     |  qn2p,   ! n2+ production rate
     |  qnp,    ! n+ production rate
     |  qnop    ! no+ production rate
!
      contains
!-----------------------------------------------------------------------
      subroutine qrj(sco2,sco1,scn2,tn,no,o2,o1,n2,n4s,xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
!
! Calculate heating and dissociation rates.
!
      use input_module,only: f107,f107a
      use init_module,only: sfeps
      use cons_module,only: t0,expz,avo,rmassinv_n4s,rmassinv_no,
     |  rmassinv_o2,rmassinv_o1,rmassinv_n2,expzmid,gask,
     |  expzmid_inv,p0,check_exp
      use lbc,only: fb,b
      use chemrates_module,only: beta9
      use fields_module,only: tlbc
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  sco2,sco1,scn2,     ! chapman integrals
     |  tn,no,o2,o1,n2,n4s, ! tn and species mass mixing ratios
     |  xnmbari             ! p0*e(-z)*barm/kT at interfaces
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: k,i,l,ier
      real,parameter ::
     |  do2=8.203E-12   ,
     |  do22=1.1407E-11 ,   
     |  aband=0.143     , ! shumann-runge
     |  bband=9.64E8    , ! shumann-runge
     |  cband=9.03E-19  , ! shumann-runge
     |  e3=0.33         ,
     |  hc = 1.9845E-16   ! C(60)
      real :: rlmeuvinv(lmax),rlmsrcinv(l1),factor,fmin,fmax
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  o2i,    ! o2  at interfaces (s1)
     |  o1i,    ! o   at interfaces (s2)
     |  n2i,    ! n2  at interfaces (s3)
     |  n4si,   ! n4s at interfaces (s4)
     |  tni     ! tn  at interfaces (s6)
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  quenchfac, ! (s8)
     |  sigchap,   ! (s9)
     |  p3f        ! (s7)

! temporary loop variables 
      real::
     |  absorp_o,         ! photoabsorption frequency of O
     |  absorp_o2,        ! photoabsoption frequency of o2
     |  absorp_n2,        ! photoabsoption frequency of n2
     |  ioniz_o,          ! photoionization frequency of o
     |  ioniz_o2,         ! photoionization frequency of o2
     |  ioniz_n2,         ! photoionization frequency of n2
     |  htfac             ! hc/wavelength, for calculation of heating
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  di_o2,          ! total dissociative ionization frequency of o2
     |  di_n2,          ! total dissociative ionization frequency of n2
     |  mn_o2,          ! transfer mass density to number density (O2)
     |  mn_o1,          ! transfer mass density to number density (O)
     |  mn_n2,          ! transfer mass density to number density (N2)
     |  mn_n            ! transfer mass density to number density (N)

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
      call vtbegin(118,ier)
#endif
!
! Exec:
!     call addfld('XNMBARI',' ',' ',xnmbari(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfsech('XNMBARI',' ',' ',xnmbari(:,lon0:lon1),
!    |  lon0,lon1,nlevs,nlevs,lat)
!
!  calculate inverse of wave length
!
      do i=1,lmax
        rlmeuvinv(i) = 1./rlmeuv(i)
      enddo
      do i=1,l1 
        rlmsrcinv(i) = 1./rlmsrc(i)
      enddo
!
! O2,O,N4S at interface levels:
      do i=lon0,lon1
        do k=lev0,lev1-1
          o2i (k+1,i) = 0.5*(o2(k,i)+o2(k+1,i))
          o1i (k+1,i) = 0.5*(o1(k,i)+o1(k+1,i))
          n2i (k+1,i) = 0.5*(n2(k,i)+n2(k+1,i))
          n4si(k+1,i) = 0.
        enddo
      enddo
!
! Bottom boundary:
      do i=lon0,lon1
        o2i(1,i) = .5*((b(i,1,1)+1.)*o2(1,i)+b(i,1,2)*o1(1,i)+
     |    fb(i,1))
        o1i(1,i) = .5*(b(i,2,1)*o2(1,i)+(b(i,2,2)+1.)*o1(1,i)+
     |    fb(i,2))
        n4si(1,i) = 0.
!
! btf tiegcm_he 8/21/13: Not sure about n2i bottom boundary. 
! For now, use 1-o2i-o1i since helium is negligible at this height.
!
        n2i(1,i) = 1.-o2i(1,i)-o1i(1,i)
      enddo
!
! N2 at interfaces:
      do k=lev0,lev1
        do i=lon0,lon1                       !!! temporal fix CISM
          if(n2i(k,i)<=0.)n2i(k,i)=1.e-5     !!! temporal fix CISM
        enddo                                !!! temporal fix CISM
      enddo
!
! calculate variables for transferring mass density to number density
      do i=lon0,lon1
        do k=lev0,lev1
           mn_o2(k,i)=o2i(k,i)*rmassinv_o2
           mn_o1(k,i)=o1i(k,i)*rmassinv_o1
           mn_n2(k,i)=n2i(k,i)*rmassinv_n2
           mn_n(k,i)=n4si(k,i)*rmassinv_n4s
        enddo
      enddo
!	 
!     call addfld('O2I'  ,' ',' ',o2i,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O1I'  ,' ',' ',o1i,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2I'  ,' ',' ',n2i,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Initialize ionization/dissociation arrays on current processor domain:
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
! Initialize local arrays
      di_o2=0.
      di_n2=0.
!
! Summation over wavelength:
      do l=l1+1,lmax   !  from 0.5A to 1050A
        sum1(:,:) = 0. ! sum(o2,o,n2)(sigma*chapman)
        sum2(:,:) = 0. ! sum(o2,o,n2)(sigma*psi/rmass)
        sum3(:,:) = 0. ! sum(o2,o,n2,n4s)(sigmas)
        htfac=hc*rlmeuvinv(l)
        do i=lon0,lon1
          do k=lev0,lev1
            sum1(k,i) = sum1(k,i)+sigeuv(1,l)*sco2(k,i)+
     |                            sigeuv(2,l)*sco1(k,i)+
     |                            sigeuv(3,l)*scn2(k,i)
!
            if (.not.check_exp) then
              sum1(k,i) = feuv(l)*exp(-sum1(k,i))
            else
              sum1(k,i) = feuv(l)*expo(-sum1(k,i),0)
            endif

            sum2(k,i) = sum2(k,i)+sigeuv(1,l)*mn_o2(k,i)+
     |                            sigeuv(2,l)*mn_o1(k,i)+
     |                            sigeuv(3,l)*mn_n2(k,i)
          enddo
        enddo
!
! Longitude and column domain of current process:
        do i=lon0,lon1
          do k=lev0,lev1
            ! absorption/ionization frequency for the three major species (O2, O, and N2)
            absorp_o2= sum1(k,i)*sigeuv(1,l)
            absorp_o= sum1(k,i)*sigeuv(2,l)
            absorp_n2= sum1(k,i)*sigeuv(3,l)
            ioniz_o2=absorp_o2*BPhotonI(1,l)
            ioniz_o=absorp_o*BPhotonI(2,l)
            ioniz_n2=absorp_n2*BPhotonI(3,l)
!
!  ionization/dissociative ionization frequency (s^-1)
            di_o2(k,i)=di_o2(k,i)+absorp_o2*bro2DIPh(l)+
     |                            ioniz_o2*bro2DIEl(l)
            di_n2(k,i)=di_n2(k,i)+absorp_n2*brn2DIPh(l)+
     |                            ioniz_n2*brn2DIEl(l)
            qnp(k,i,lat) = qnp(k,i,lat)+sigin4s(l)*sum1(k,i)
            qn2p(k,i,lat) = qn2p(k,i,lat)+ioniz_n2+
     |                                ioniz_n2*BElectronI(3,l)
            qo2p(k,i,lat) = qo2p(k,i,lat)+ioniz_o2+
     |                                ioniz_o2*BElectronI(1,l)
            qop2p(k,i,lat) = qop2p(k,i,lat)+absorp_o*brop2pPh(l)+
     |                                      ioniz_o*brop2pEl(l) 
            qop2d(k,i,lat) = qop2d(k,i,lat)+absorp_o*brop2dPh(l)+
     |                                      ioniz_o*brop2dEl(l)
            qop(k,i,lat) = qop(k,i,lat)+absorp_o*brop4sPh(l)+
     |                                  ioniz_o*brop4sEl(l)
!
! total dissociation and EUV heating
            rj(k,i,lat) = rj(k,i,lat)+(absorp_o2*bro2DPh(l)
     |                      +ioniz_o2*bro2DEl(l))
            qtef(k,i,lat) = qtef(k,i,lat)+absorp_n2*brn2DPh(l)
     |                      +ioniz_n2*brn2DEl(l)
            qtotal(k,i,lat) = qtotal(k,i,lat)+htfac*
     |                      sum1(k,i)*sum2(k,i)
          enddo
        enddo
      enddo ! l=l1+1,lmax
!
!     call addfld('QOP2P' ,' ',' ', qop2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP2D' ,' ',' ', qop2d(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SUM1' ,' ',' ', sum1,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('Q0'  ,' ',' ',qtotal(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RJ_QRJ',' ',' ',rj(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Multiply Q by efficiency factor:
      do i=lon0,lon1
        do k=lev0,lev1
          qtotal(k,i,lat) = qtotal(k,i,lat)*euveff(k)*avo
        enddo
      enddo
!     call addfld('Q1'  ,' ',' ',qtotal(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! transfer frequency (s^-1) to production rates (cm^-3 s^-1)
      do i=lon0,lon1
        do k=lev0,lev1
           mn_o2(k,i)=mn_o2(k,i)*xnmbari(k,i)
           mn_n2(k,i)=mn_n2(k,i)*xnmbari(k,i)
           mn_o1(k,i)=mn_o1(k,i)*xnmbari(k,i)
           mn_n(k,i)=mn_n(k,i)*xnmbari(k,i)
           di_o2(k,i)=di_o2(k,i)*mn_o2(k,i)
           di_n2(k,i)=di_n2(k,i)*mn_n2(k,i)
           qo2p(k,i,lat) = qo2p(k,i,lat)*mn_o2(k,i)-di_o2(k,i)
           qn2p(k,i,lat) = qn2p(k,i,lat)*mn_n2(k,i)-di_n2(k,i)
           qnp(k,i,lat)  = qnp(k,i,lat)*mn_n(k,i)+di_n2(k,i)
           qop(k,i,lat) = qop(k,i,lat)*mn_o1(k,i)+0.54*di_o2(k,i)
           qop2p(k,i,lat) = qop2p(k,i,lat)*mn_o1(k,i)+0.22*di_o2(k,i)
           qop2d(k,i,lat) = qop2d(k,i,lat)*mn_o1(k,i)+0.24*di_o2(k,i)
           qtef(k,i,lat) = 2.*qtef(k,i,lat)*mn_n2(k,i)+di_n2(k,i)
        enddo
      enddo
!
!     call addfld('QOP',' ',' ',qop(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
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
        tni(lev0,i) = tlbc(i,lat)
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
      do l=1,l1       ! from 1050A to 1750A
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
!
!     call addfld('Q2'  ,' ',' ',qtotal(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
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
!     call addfld('Q'     ,' ',' ',qtotal(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RJ'    ,' ',' ',rj(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP2Pa' ,' ',' ',qop2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP2Da' ,' ',' ', qop2d(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QTEF'  ,' ',' ',  qtef(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QO2P'  ,' ',' ',  qo2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QN2P'  ,' ',' ',  qn2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QNP'   ,' ',' ',   qnp(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOPa'   ,' ',' ',  qop(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QNOP'  ,' ',' ',  qnop(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
#ifdef VT
!     code = 118 ; state = 'qrj' ; activity='ModelCode'
      call vtend(118,ier)
#endif
!
      end subroutine qrj
!---------------------------------------------------------------
      subroutine init_sflux
      use input_module,only: f107,f107a,see_ncfile
      use init_module,only: sfeps,iday,iyear,secs,istep
      use soldata_module,only: get_soldata,soldata,nwave
!
! Flux initialization once per time step, called from advance.
!
! Local:
      integer :: n
! sflu: scaled solar flux returned by subroutine ssflux() (photons cm-2 s-1)
!
      if (len_trim(see_ncfile) == 0) then
         call ssflux(f107,f107a)
         do n = l1+1,lmax
              feuv(n) = sflux(n)*sfeps
         enddo
      
         do n = 1,l1
           fsrc(n)   = sflux(n)*sfeps
         enddo
      else
         call get_soldata(iyear,iday,int(secs),istep)
         if (nwave /= lmax) then
            write(6,"('init_sflux(): wave bins mismatch: nwave=',i4,
     |         ' lmax=',i4)") nwave,lmax
            call shutdown('init_sflux')
         endif

         do n=1,l1
            fsrc(n)=soldata(n)
         enddo
         do n=l1+1,lmax
            feuv(n)=soldata(n)
         enddo
      endif
      end subroutine init_sflux
!-----------------------------------------------------------------------
      subroutine init_qrj

      integer :: m,n
!
! Called once per run, from tgcm.F
!
! Initialize bins (37 bins)
!
      wave1 = (/1700.00, 1650.00, 1600.00, 1550.00, 1500.00,
     |          1450.00, 1400.00, 1350.00, 1300.00, 1250.00,
     |          1200.00, 1215.67, 1150.00, 1100.00, 1050.00,
     |          1027.00,  987.00,  975.00,  913.00,  913.00,
     |           913.00,  798.00,  798.00,  798.00,  650.00,
     |           650.00,  540.00,  320.00,  290.00,  224.00,
     |           155.00,   70.00,   32.00,   18.00,    8.00,
     |             4.00,    0.50/)
      wave2 = (/1750.00, 1700.00, 1650.00, 1600.00, 1550.00,
     |          1500.00, 1450.00, 1400.00, 1350.00, 1300.00,
     |          1250.00, 1215.67, 1200.00, 1150.00, 1100.00,
     |          1050.00, 1027.00,  987.00,  975.00,  975.00,
     |           975.00,  913.00,  913.00,  913.00,  798.00,
     |           798.00,  650.00,  540.00,  320.00,  290.00,
     |           224.00,  155.00,   70.00,   32.00,   18.00,
     |             8.00,    4.00/)
!
! Solar spectrum based on EUVAC and glow for wave length less than 1050 A
! and Woods for wavelength greater than 1050 A
!
! solar minimum flux (when P_index=80, unit:photon cm^-2 S^-1)
!
      sfmin=(/3.397e+11, 1.998e+11, 1.055e+11, 7.260e+10,
     |        5.080e+10, 2.802e+10, 1.824e+10, 1.387e+10,
     |        2.659e+10, 7.790e+09, 1.509e+10, 3.940e+11,
     |        8.399e+09, 3.200e+09, 3.298e+09, 4.235e+09,
     |        4.419e+09, 4.482e+09, 7.156e+08, 1.028e+09,
     |        3.818e+08, 8.448e+08, 3.655e+09, 2.364e+09,
     |        1.142e+09, 1.459e+09, 4.830e+09, 2.861e+09,
     |        8.380e+09, 4.342e+09, 5.612e+09, 1.270e+09,
     |        5.326e+08, 2.850e+07, 2.000e+06, 1.000e+04,
     |        5.010e+01/)
!
! scaling factor A as defined in EUVAC model
!
      afac=(/5.937e-04, 6.089e-04, 1.043e-03, 1.125e-03,
     |       1.531e-03, 1.202e-03, 1.873e-03, 2.632e-03,
     |       2.877e-03, 2.610e-03, 3.739e-03, 4.230e-03,
     |       2.541e-03, 2.099e-03, 3.007e-03, 4.825e-03,
     |       5.021e-03, 3.950e-03, 4.422e-03, 4.955e-03,
     |       4.915e-03, 5.437e-03, 5.261e-03, 5.310e-03,
     |       3.680e-03, 5.719e-03, 5.857e-03, 1.458e-02,
     |       7.059e-03, 2.575e-02, 1.433e-02, 9.182e-03,
     |       1.343e-02, 6.247e-02, 2.000e-01, 3.710e-01,
     |       6.240e-01/)
!
!   0.5*(wave1+wave2) in centimeter
!
      rlmeuv=(/1.725e-05, 1.675e-05, 1.625e-05, 1.575e-05,
     |         1.525e-05, 1.475e-05, 1.425e-05, 1.375e-05,
     |         1.325e-05, 1.275e-05, 1.225e-05, 1.216e-05,
     |         1.175e-05, 1.125e-05, 1.075e-05, 1.038e-05,
     |         1.007e-05, 9.810e-06, 9.440e-06, 9.440e-06,
     |         9.440e-06, 8.555e-06, 8.555e-06, 8.555e-06,
     |         7.240e-06, 7.240e-06, 5.950e-06, 4.300e-06,
     |         3.050e-06, 2.570e-06, 1.895e-06, 1.125e-06,
     |         5.100e-07, 2.500e-07, 1.300e-07, 6.000e-08,
     |         2.250e-08/)
!
! O2 absorption coefficient:
!
      sigeuv(1,:) = (/
     |          5.00e-01, 1.50e+00, 3.40e+00, 6.00e+00, 1.00e+01,
     |          1.30e+01, 1.50e+01, 1.20e+01, 2.20e+00, 4.00e-01,
     |          1.30e+01, 1.00e-02, 1.40e+00, 4.00e-01, 1.00e+00,
     |          1.15e+00, 1.63e+00, 1.87e+01, 3.25e+01, 1.44e+01,
     |          1.34e+01, 1.33e+01, 1.09e+01, 1.05e+01, 2.49e+01,
     |          2.36e+01, 2.70e+01, 2.03e+01, 1.68e+01, 1.32e+01,
     |          7.63e+00, 2.63e+00, 6.46e-01, 2.10e-01, 2.25e-01,
     |          3.40e-02, 4.54e-03/)
!
! O absorption coefficient:
!
      sigeuv(2,:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 3.79e+00, 4.10e+00, 3.00e+00, 4.79e+00,
     |          8.52e+00, 1.31e+01, 1.07e+01, 7.72e+00, 6.02e+00,
     |          3.78e+00, 1.32e+00, 3.25e-01, 1.05e-01, 1.13e-01,
     |          1.70e-02, 2.27e-03/)
!
! N2 absorption coefficient:
!
      sigeuv(3,:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 2.55e+00, 1.15e+02, 1.44e+01,
     |          2.18e+00, 7.17e+01, 1.31e+01, 2.14e+00, 5.45e+01,
     |          2.30e+01, 2.31e+01, 1.97e+01, 1.17e+01, 9.94e+00,
     |          5.09e+00, 1.53e+00, 3.46e-01, 1.14e+00, 1.41e-01,
     |          2.01e-02, 2.53e-03/)
!
! The three major species' ionization branching ratio (off absorption):
!
! O2 
!
      BPhotonI(1,:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 6.13e-01, 8.30e-01, 6.20e-01, 7.86e-01,
     |          7.56e-01, 5.34e-01, 5.74e-01, 5.49e-01, 4.76e-01,
     |          6.73e-01, 9.83e-01, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00/)
!
! O
!
      BPhotonI(2,:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00/)
!
! N2
!
      BPhotonI(3,:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 4.29e-01,
     |          6.80e-01, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00/)
!
! photon ionization branching ratio for O+(2p),o+(2d),O+(4s) 
! (off O photon ionization)
!
! O+(2p)
!
      brop2pPh(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          8.56e-03, 2.52e-01, 2.60e-01, 2.46e-01, 2.41e-01,
     |          2.33e-01, 2.27e-01, 2.26e-01, 2.24e-01, 2.24e-01,
     |          2.24e-01, 2.24e-01/)
!
! O+(2d)
!
      brop2dPh(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 6.98e-02,
     |          3.37e-01, 4.51e-01, 4.24e-01, 4.03e-01, 4.02e-01,
     |          3.92e-01, 3.77e-01, 3.74e-01, 3.78e-01, 3.78e-01,
     |          3.78e-01, 3.78e-01/)
!
! O+(4s)
!
      brop4sPh(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 9.30e-01,
     |          6.55e-01, 2.98e-01, 3.17e-01, 3.46e-01, 3.50e-01,
     |          3.67e-01, 3.89e-01, 3.93e-01, 3.90e-01, 3.90e-01,
     |          3.90e-01, 3.90e-01/)
!
! O2 photon dissociation braching ratio
!
      bro2DPh(:) = (/
     |          1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 3.87e-01, 1.70e-01, 3.80e-01, 2.14e-01,
     |          2.44e-01, 4.66e-01, 4.26e-01, 4.51e-01, 5.24e-01,
     |          3.27e-01, 1.74e-02, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00/)
!
! n2 photon dissociation braching ratio for n(2d)
!
      brn2DPh(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 5.71e-01,
     |          3.20e-01, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00/)
!
! O2 photon dissociative ionization braching ratio
!
      bro2DIPh(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          5.35e-04, 1.08e-01, 2.40e-01, 3.51e-01, 3.76e-01,
     |          4.47e-01, 6.53e-01, 8.92e-01, 1.00e+00, 1.00e+00,
     |          1.00e+00, 1.00e+00/)
!
! n2 photon dissociative ionization braching ratio
!
      brn2DIPh(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 4.72e-03, 9.27e-02, 2.46e-01,
     |          2.53e-01, 2.49e-01, 2.82e-01, 9.60e-01, 9.60e-01,
     |          9.60e-01, 9.60e-01/)
!
! n(4s) photoionization cross section
!
      sigin4s(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 3.24e+00, 2.48e+00, 2.11e+00, 1.12e+01,
     |          1.18e+01, 1.16e+01, 9.35e+00, 6.43e+00, 4.80e+00,
     |          2.55e+00, 6.81e-01, 1.66e-01, 5.68e-01, 7.05e-02,
     |          1.00e-02, 1.27e-03/)
!
! The three major species' electron impact ionization branching ratio
! off its photon ionization rate
!
! O2
!
      BElectronI(1,:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 2.38e-02, 1.05e-01, 2.42e-01,
     |          5.79e-01, 1.61e+00, 4.27e+00, 6.00e+01, 2.03e+01,
     |          5.02e+01, 2.11e+02/)
!
! O
!
      BElectronI(2,:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 1.27e-01, 4.18e-01, 6.94e-01,
     |          1.09e+00, 2.19e+00, 4.99e+00, 7.14e+01, 2.36e+01,
     |          5.06e+01, 2.17e+02/)
!
! N2
!
      BElectronI(3,:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 3.07e-02, 1.78e-01, 3.61e-01,
     |          9.33e-01, 2.86e+00, 7.79e+00, 1.08e+01, 3.22e+01,
     |          8.09e+01, 3.43e+02/)
!
! electron impact ionization branching ratio for O+(2p),o+(2d),O+(4s) 
! (off O photon ionization)
!
! O+(2p)
!
      brop2pEl(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 9.15e-03, 6.12e-02, 1.16e-01,
     |          2.03e-01, 4.36e-01, 1.01e+00, 1.46e+01, 4.77e+00,
     |          1.10e+01, 4.74e+01/)
!
! O+(2d)      
!
      brop2dEl(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 3.39e-02, 1.48e-01, 2.53e-01,
     |          4.18e-01, 8.53e-01, 1.96e+00, 2.82e+01, 9.36e+00,
     |          2.07e+01, 8.85e+01/)
!
! O+(4s)      
!
      brop4sEl(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 8.36e-02, 2.09e-01, 3.25e-01,
     |          4.70e-01, 9.02e-01, 2.02e+00, 2.86e+01, 9.42e+00,
     |          1.89e+01, 8.12e+01/)
!
! photoelectron dissociative ionization branching ratio of O2
!
      bro2DIEl(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 9.60e-04, 1.43e-02, 5.23e-02,
     |          1.63e-01, 5.21e-01, 1.44e+00, 2.03e+01, 6.98e+00,
     |          1.79e+01, 7.61e+01/)
!
! photoelectron dissociative ionization branching ratio of N2
!
      brn2DIEl(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 1.84e-04, 8.49e-03, 3.66e-02,
     |          1.46e-01, 5.71e-01, 1.65e+00, 2.29e+00, 6.95e+00,
     |          1.83e+01, 7.87e+01/)
!
! photoelectron dissociation branching ratio of N2
!
      brn2DEl(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 1.57e-01, 5.15e-01, 7.64e-01,
     |          1.37e+00, 2.91e+00, 6.53e+00, 9.05e+00, 2.53e+01,
     |          5.21e+01, 2.45e+02/)
!
! photoelectron dissociation branching ratio of o2
!
      bro2DEl(:) = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 1.10e-02, 6.53e-01, 7.62e-01, 9.96e-01,
     |          1.27e+00, 2.04e+00, 4.11e+00, 5.70e+01, 1.78e+01,
     |          2.03e+01, 8.79e+01/)
 
 
! 
! transfer units of cross section to cm^2
      do m=1,3
        do n=1,lmax
          sigeuv(m,n)=sigeuv(m,n)*1.E-18
        enddo
      enddo
      do n = 1,lmax
        sigin4s(n) = sigin4s(n)*1.E-18
      enddo
! 
      do n = 1,l1
        rlmsrc(n) = rlmeuv(n)
        sigsrc(n) = sigeuv(1,n)
      enddo
!
      euveff(:) = 0.05
      quench = (/7.E-11,5.E-11,3.1401E-12,9.1E-3/)

      end subroutine init_qrj
!---------------------------------------------------------------
      subroutine ssflux (f107, f107a)
!
! Args:
      real,intent(in) :: f107,f107a
!
! Local:
      real ::  pind
      integer :: l
!
! solar model is the same as EUVAC, i.e., solar flux at 
! (f107,f107a) is scale as: sflux=sfmin(1+afac(P-80.)) 
! where P=0.5(f107+f107a)

      pind=0.5*(f107+f107a)
      do l=1,lmax
         sflux(l)=sfmin(l)*(1+afac(l)*(pind-80.))
         ! set solar flux  to be 80% of the value when pind=80
         ! if it becomes negative
         if (sflux(l) .le. 0.1*sfmin(l)) sflux(l) = 0.1*sfmin(l)
      enddo
!
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
