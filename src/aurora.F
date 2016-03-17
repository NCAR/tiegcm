!
      module aurora_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Auroral oval parameterization. See reference:
! R.G. Roble, E.C. Ridley
! An auroral model for the NCAR thermospheric general circulation model (TGCM)
! Annales Geophysicae,5A, (6), 369-382, 1987. 
!
! The aurora oval is a circle in auroral circle coordinates.  Auroral circle
!  coordinates are offset from magnetic coordinates by offa degrees (radians)
!  towards 0 MLT and by dskofa degrees (radians) towards dusk (18 MLT).
! The aurora assumes a Maxwellian in energy, so that the characteristic
!  energy is half of the mean energy (or mean energy = 2*alfa, where alfa
!  is the characteristic energy).  The Maxwellian is approximated in the
!  aion and bion subroutines.
! The aurora oval is assumed to be a Gaussian in auroral latitude, with
!  peak values on the day (=1) and night (=2) sides that change from one to
!  the other using cosines of the auroral longitude coordinate.
! There is provision for a low energy (~75 eV) aurora at the location of the
!  regular (~1-6 keV) aurora in order to simulate the energy flux found
!  at higher altitudes that is non-Maxwellian, but the flux is usually
!  set to zero (1.e-80).
! There is provision for a proton (MeV) aurora, but the flux is usually
!  set to zero (1.e-20).
! The drizzle is a constant low energy electron flux over the polar cap,
!  which goes to 1/e over twice the 1/e-width of the aurora at the
!  radius of the aurora.
! The cusp is a low energy electron flux centered over the dayside convection
!  entrance at phid at the convection reversal boundary theta0.  The cusp
!  falls off over 5 degrees in latitude and over 20 degrees in longitude
!  to 1/e values of the peak at the center.
! 1.e-20 and 1.e-80 are used to give a near zero answer.
!
! The polar drizzle and cusp electron energies are low, and soft particles
!  have great influence on the over-all thermospheric and ionospheric
!  structure, especially on the electron density profiles at mid-latitudes
!  and in winter since low energy electrons produce ionization at high
!  altitudes where loss rates are very low.  (Comment by Wenbin Wang.)
! The original energies for drizzle and cusp were alfad=0.75, alfac=0.5 keV.
! The original guess at energy fluxes were: ed=0.1+2.0*power/100.,ec=0.1+0.9*power/100.
! The next guess at energy fluxes were: ed=0.01+0.2*power/100., ec=0.01+0.09*power/100.
! The values below reflect higher estimates for the electron energy (lower alt)
!
! Calling sequence (all subs in aurora_module, aurora.F):
!   1) sub aurora_cons called once per time step from advance.
!   2) sub aurora called from dynamics, inside parallel latitude scan.
!   3) subs aurora_cusp and aurora_heat called from sub aurora.
!   4) sub aurora_ions called from sub aurora. 
!
! 1/4/08 btf: Updating with Liying and Wenbin's mods, see e1,e2, h1,h2,
!             and mods to gpi.F. 
! 1/24/08 btf:Changed alfad,ed, alfac,ec as per Wenbin.
!
      use cons_module,only: ! see cons.F
     |  dphi, ! delta latitude geographic (radians)
     |  dtr,  ! degrees to radians (pi/180)
     |  pi,    ! 4.*atan(1.)
     |  crit  ! may be changed from default if amie run
      use params_module,only: nlon,nlonp4,nlevp1,nlat,glat
      use addfld_module,only: addfld
      use amie_module,only: crad,phida,ekvg,efxg,hpi_sh_amie,
     |  hpi_nh_amie,pcp_sh_amie,pcp_nh_amie
      use input_module,only: iamie
      implicit none
      save
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! 11/08 EMERY added eflux=flux*fac_p2e*alfa   (fac_p2e=2*1.602e-9, so 2*alfa in this factor)
      real, parameter:: fac_p2e = 3.204e-9 ! convert from particle to energy flux

      real,parameter :: h2deg = 15.     ! convert from hours to degrees

! Grid dimensions:
      logical :: do_aurora ! flag to calculate aurora at current lat
      integer :: ihem      ! hemisphere at current latitude (1=SH,2=NH)
      integer,parameter :: isouth=1, inorth=2
!
! Polar drizzle parameters:
!   alfad: Characteristic Maxwellian energy of drizzle electrons (keV)
!   ed   : Column energy input of drizzle electrons (ergs/cm**2/s)
!   fd   : Electron particle flux of drizzle electrons (particles/cm**2/s)
! 1/24/08 btf: reduce alfad to 0.5, as per Wenbin:
!
      real,parameter :: alfad = 0.5
      real :: ed ! set in sub aurora_cons as function of hem power
      real :: fd ! set in sub aurora_ions
!
! Polar cusp parameters:
!   alfac: Characteristic Maxwellian energy of polar cusp electons (keV)
!   ec   : Column energy input of polar cusp electrons (ergs/cm**2/s)
!   fc   : Electron particle flux of polar cusp electrons (particles/cm**2/s)
!
! 1/24/08 btf: As per Wenbin: alfac = 0.1 (old tiegcm1.8 value was 1.0), 
!              and ec is set as function of hemispheric power in aurora_cons 
!              (old tiegcm1.8 value for ec was 0.5)
!
      real,parameter :: alfac = 0.1
      real :: ec ! set in sub aurora_cons (function of hem power)
!
      real :: fc ! set in sub aurora_ions
!
! e1: Peak energy flux in noon sector of the aurora (ergs/cm**2/s)
! e2: Peak energy flux in midnight sector of the aurora (ergs/cm**2/s)
! h1: Gaussian 1/e-width of the noon auroral oval in degrees
! h2: Gaussian 1/e-width of the midnight auroral oval in degrees
!
      real ::
     |  e1,e2, ! set in sub aurora_cons (function of hem power)
     |  h1,h2  ! set in sub aurora_cons (function of hem power)
!
! Solar proton parameters for the polar cap (time-gcm only)
!   alfa_sp: Characteristic Maxwellian energy of solar protons (MeV) (was alfad2)
!   e_sp   : Column energy input of solar protons (ergs/cm**2/s) (was ed2)
!   flx_sp : e_sp/1.602e-6, for input to sub bion                (was fd2)
! Add solar protons to ionization if add_sproton is true (time-gcm only)
!
      logical :: add_sproton = .false.
      real,parameter ::
     |  alfa_sp = 10.,
     |  e_sp = 1.e-20
      real :: flx_sp
!
! High energy electron parameters in the auroral oval (time-gcm only):
      logical :: add_helectron = .false.
      real,parameter :: 
     |  alfa30 = 40., ! Characteristic energy of auroral electrons (Mev)
     |  e30 = .05     ! Column energy of auroral electrons (ergs/cm**2/s)
!
! 1/02 bf: low energy protons not included at this time.
!     real,parameter :: ! low energy protons:
!    |  alfap1 = 5.,
!    |  alfap2 = 15.,
!    |  pe1 = 1.e-20,
!    |  pe2 = 1.e-20
!     real :: alfap0,e0p
!
! Additional auroral parameters (see sub aurora_cons):
! (dimension 2 is for south, north hemispheres)
      real ::
     |  theta0(2),    ! convection reversal boundary in radians
     |  offa(2),      ! offset of oval towards 0 MLT relative to magnetic pole (rad)
     |  dskofa(2),    ! offset of oval in radians towards 18 MLT (f(By))
     |  phid(2),      ! dayside convection entrance in MLT (subtract 12h since 0=noon) converted to radians (f(By))
     |  phin(2),      ! night convection entrance in MLT (subtract 12h since 0=noon) converted to radians (f(By))
     |  rrad(2),      ! radius of auroral circle in radians
     |  rradp(2),     ! radius of polar-cap circle in radians for drizls
     |  alfa0,        ! average of noon and midnight characteristic Maxw energies
     |  ralfa,ralfa2, ! difference ratios of characteristic energies
     |  rrote,        ! MLT of maximum auroral energy flux converted to radians
     |  rroth,        ! MLT of maximum auroral Gaussian 1/e-width converted to radians
     |  h0,           ! average of noon and midnight Gaussian 1/e-widths
     |  rh,           ! difference ratio of 1/e-widths (rh=(h2-h1)/(h2+h1))
     |  e0,e20,       ! e0 = average of noon and midnight electrons
     |  ree,re2,      ! difference ratios of peak energy fluxes (ree=(e2-e1)/(e2+e1))
     |  alfa20        ! average of noon and midnight char energies for high alt aurora
!
! The following parameters are used only by heelis module for dynamo:
      real ::
     |  offc(2),      ! offset of convection towards 0 MLT relative to magnetic pole (rad)
     |  dskofc(2),    ! offset of convection in radians towards 18 MLT (f(By))
     |  psim(2),      ! maximum potential in the morning cell (V)
     |  psie(2),      ! minimum potential in the evening cell (V)
     |  pcen(2),      ! potential at the center (offc,dskofc) of the convection pattern (V, f(By))
     |  phidp0(2),    ! angle curvature of convection on plus side of dayside entrance (rad)
     |  phidm0(2),    ! angle curvature of convection on minus side of dayside entrance (rad)
     |  phinp0(2),    ! angle curvature of convection on plus side of nightside entrance (rad)
     |  phinm0(2),    ! angle curvature of convection on minus side of nightside entrance (rad)
     |  rr1(2)        ! exponential fall-off of convection from convection radius
!
! Auroral circle coordinates (see sub aurora):
      real ::
     |  dlat_aur(nlon), ! = rlatm = magnetic latitude (radians)
     |  dlon_aur(nlon), ! = rlonm-sunlons = magnetic longitude (or MLT-12 where 0=noon) converted to radians
     |  colat(nlon),    ! auroral magnetic co-latitude from offset center (radians)
     |  sinlat(nlon), coslat(nlon),     ! sin,cos of dlat_aur (rlatm) magnetic latitude
     |  coslon(nlon), sinlon(nlon),     ! sin,cos of dlon_aur+offset
     |  alon(nlon)                      ! auroral magnetic longitude where 0=noon shifted by -offset (radians)
!
! Local-time dependent energy and particle fluxes:
! (see subs aurora_heat and aurora_cusp)
      real ::
     |  alfa(nlonp4), alfa2(nlonp4), alfa3(nlonp4), ! characteristic energy
     |  flux(nlonp4), flux2(nlonp4), flux3(nlonp4), ! particle flux
     |  drizl(nlonp4),cusp(nlonp4),                 ! polar drizzle and cusp
     |  qteaur(nlonp4,nlat),                        ! for electron temperature
     |  eflux(nlonp4),eflux2(nlonp4),eflux3(nlonp4) ! energy flux

! 11/08 EMERY added energy flux eflux and eflux2d
      real :: cusp2d(nlonp4,nlat),drzl2d(nlonp4,nlat),
     |        alfa2d(nlonp4,nlat),nflx2d(nlonp4,nlat),
     |        eflux2d(nlonp4,nlat)
!
! aureff = vertical profile for auroral heating efficiency (used by qjion)
! 
      real :: aureff(nlevp1)
      real :: bdriz(nlevp1)   ! vertical profile for background drizzle
!
! Subroutines internal to aurora module.
      contains
!-----------------------------------------------------------------------
#include <defs.h>
      subroutine aurora_cons(iprint,iyear,iday,secs)
!
! Set auroral constants (called from sub advance once per time-step).
!
      use params_module,only: dz,nlev
      use hist_module,only:   modeltime
      use input_module,only:  ! from user input
     |  byimf,   ! BY component of IMF (nT)      (e.g., 0.)
     |  ctpoten, ! cross-cap potential (kV)      (e.g., 45.)
     |  power,   ! hemispheric power   (GW)      (e.g., 16.)
     |  iamie,
     |  potential_model  ! electric potential model used
!
! Args:
      integer,intent(in) :: iprint,iday,iyear
      real,intent(in) :: secs
!
! Local:
      real,parameter :: convert=3.1211e+8
      real :: alfa_1,alfa_2,plevel,roth,rote,rcp,rhp,
     |  alfa21,alfa22,e21,e22
      real :: arad(2)
      real :: hp_amie,c25,c35
      real :: byloc   ! local by
      real,parameter :: e=1.e-10
!
! Auroral heating efficiency (nlevp1) according to vertical resolution 
! (NLEV==28, dz==.5, or NLEV==56, dz==.25)
!
#if (NLEV==28)
      aureff =
     |(/0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,
     |  0.55,0.50,0.45,0.37,0.30,0.20,0.16,0.13,0.10,0.08,0.07,0.06,
     |  0.05,0.05,0.05,0.05,0.05/)
      bdriz =
     |(/2.17E-4,4.48E-3,3.78E-2,1.55E-1,3.26E-1,3.85E-1,3.50E-1,2.81E-1,
     |  2.04E-1,1.38E-1,8.74E-2,5.23E-2,3.01E-2,1.68E-2,9.27E-3,5.08E-3,
     |  2.79E-3,1.55E-3,8.62E-4,4.85E-4,2.76E-4,1.58E-4,9.11E-5,5.28E-5,
     |  3.08E-5,1.79E-5,1.04E-5,6.01E-6,3.49E-6/)
#elif (NLEV==56)
        aureff =
     |  (/0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00,
     |    0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00,
     |    0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00,
     |    0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00,
     |    0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00, 0.5500E+00,
     |    0.5250E+00, 0.5000E+00, 0.4750E+00, 0.4500E+00, 0.4100E+00,
     |    0.3700E+00, 0.3350E+00, 0.3000E+00, 0.2500E+00, 0.2000E+00,
     |    0.1800E+00, 0.1600E+00, 0.1450E+00, 0.1300E+00, 0.1150E+00,
     |    0.1000E+00, 0.9000E-01, 0.8000E-01, 0.7500E-01, 0.7000E-01,
     |    0.6500E-01, 0.6000E-01, 0.5500E-01, 0.5000E-01, 0.5000E-01,
     |    0.5000E-01, 0.5000E-01, 0.5000E-01, 0.5000E-01, 0.5000E-01,
     |    0.5000E-01, 0.5000E-01/)
        bdriz =
     |  (/0.2170E-03, 0.2349E-02, 0.4480E-02, 0.2114E-01, 0.3780E-01,
     |    0.9640E-01, 0.1550E+00, 0.2405E+00, 0.3260E+00, 0.3555E+00,
     |    0.3850E+00, 0.3675E+00, 0.3500E+00, 0.3155E+00, 0.2810E+00,
     |    0.2425E+00, 0.2040E+00, 0.1710E+00, 0.1380E+00, 0.1127E+00,
     |    0.8740E-01, 0.6985E-01, 0.5230E-01, 0.4120E-01, 0.3010E-01,
     |    0.2345E-01, 0.1680E-01, 0.1303E-01, 0.9270E-02, 0.7175E-02,
     |    0.5080E-02, 0.3935E-02, 0.2790E-02, 0.2170E-02, 0.1550E-02,
     |    0.1206E-02, 0.8620E-03, 0.6735E-03, 0.4850E-03, 0.3805E-03,
     |    0.2760E-03, 0.2170E-03, 0.1580E-03, 0.1246E-03, 0.9110E-04,
     |    0.7195E-04, 0.5280E-04, 0.4180E-04, 0.3080E-04, 0.2435E-04,
     |    0.1790E-04, 0.1415E-04, 0.1040E-04, 0.8205E-05, 0.6010E-05,
     |    0.4750E-05, 0.3490E-05/)
#else
        write(6,"('>>> aurora: unsupported NLEV=',i4,' dz=',f10.4)")
     |    nlev,dz
        call shutdown('aurora dz')
#endif
!
! Add limits to byimf if use the Heelis convection pattern,this is to have 
! asymmetric dawn and dusk convection cells and By effect. --Wenbin Wang 12/02/2008
!
      byloc = byimf ! init local from original namelist input
      If (potential_model == 'HEELIS') then
        if (byloc > 7.0) byloc = 7.0
        if (byloc < -11.0) byloc = -11.0         
      endif

      plevel = 0.
! 8/1/06 btf: change conditional as per Richmond suggestion:
!     if (power >= 0.01) plevel = 2.09*alog(power)
      if (power >= 1.00) plevel = 2.09*alog(power)
!
! e1 = energy flux in the noon sector of the aurora (ergs/cm**2/s)
! e2 = energy flux in the midnight sector of the aurora (ergs/cm**2/s)
! modified by LQIAN, 2007
! produce realistic oval compared to NOAA empirical auroral oval and TIMED/GUVI
! e1 formula given by Wenbin base on POLARVIS image;
! e2 formula based on Emery et al original auroral parameterization report
      e1 = max(0.50, -2.15 + 0.62*plevel)
      e2=1.+0.11*power
      e0 = 0.5 * (e1 + e2)
      ree = (e2 - e1) / (e1 + e2)
!
! 1/24/08 btf: Set ec and ed as function of power, as per Wenbin:
      ec=(0.24+0.0067*power)/5. 
      ed=0.0012+0.0006*power
!
! h1 = Gaussian 1/e-width of the noon auroral oval in degrees
! h2 = Gaussian 1/e-width of the midnight auroral oval in degrees
! modified by LQIAN, 2007
! produce realistic oval compared to NOAA empirical auroral oval and TIMED/GUVI
! h1 formula given by Wenbin base on POLARVIS image;
! h2 formula based on Emery et al original auroral parameterization report
      h1 = min(2.35, 0.83 + 0.33*plevel)
      h2=2.5+0.025*max(power,55.)+0.01*min(0.,power-55.)
      rh = (h2 - h1) / (h1 + h2)
!
      theta0(isouth) = (-3.80+8.48*(ctpoten**0.1875))*dtr
#if defined(INTERCOMM) || defined(CISMAH)
! Set theta0 = 10 deg so crit(1,2)=15,30 (old values in cons.F) in colath.F for CISM
! This is commented as per coupled tiegcm1.94.2:
!     theta0(isouth) = 10.*dtr
#endif
      theta0(inorth) = theta0(isouth)
      offa(isouth) = 1.0*dtr
      offa(inorth) = 1.0*dtr
      dskofa(isouth) = 0.
      dskofa(inorth) = 0.
!
      rhp = 14.20 + 0.96*plevel
      rcp = -0.43 + 9.69 * (ctpoten**0.1875)
      arad(isouth) = max(rcp,rhp)
      arad(inorth) = max(rcp,rhp)
      rrad(isouth) = arad(isouth)*dtr 
      rrad(inorth) = arad(inorth)*dtr
!
! 1/4/08 btf: New values for alfa_1 and alfa_2, as per Wenbin
      alfa_1 = 1.5
      alfa_2 = 2.
      alfa0 = 0.5*(alfa_1+alfa_2)
      ralfa = (alfa_2-alfa_1) / (alfa_1+alfa_2+e)

! 12/8/08 scs: In-lined rotation units and cleaned up comments
!              Old rotation angles were:  roth = 12.18 - 0.89 * plevel
!                                         rote =  2.62 - 0.55 * plevel

! roth = MLT of max width of aurora in hours
! rote = MLT of max energy flux of aurora in hours

      roth = 0.81 - 0.06 * plevel
      rote = 0.17 - 0.04 * plevel

!  Convert MLT from hours to degrees to radians

      rroth = roth * h2deg * dtr
      rrote = rote * h2deg * dtr

      h0 = 0.5 * (h1 + h2) * dtr
      alfa21 = 0.075
      alfa22 = 0.075
      alfa20 = 0.5 * (alfa21 + alfa22)
      ralfa2 = (alfa22 - alfa21) / (alfa21 + alfa22 + e)
      e21 = 1.e-80
      e22 = 1.e-80
      e20 = 0.5 * (e21 + e22)
      re2 = (e22 - e21) / (e21 + e22)
!
! Set cusp and drizzle parameters:
! (conversion between particle number density and characteristic
!  energy and column energy input)
!
      fc = convert * ec / alfac
      fd = convert * ed / alfad
!
! Solar proton flux:
      flx_sp = e_sp/1.602e-6
!
! The following parameters (offc through rr1) are used only in heelis 
!   potential calculation for the dynamo (see heelis.F)
!
!  tiegcm1.9
!     offc(isouth) = 1.*dtr
!     offc(inorth) = 1.*dtr
!     dskofc(isouth) = 0.
!     dskofc(inorth) = 0.
!     phid(isouth) = 0.
!     phid(inorth) = 0.
!     phin(isouth) = 180.*dtr
!     phin(inorth) = 180.*dtr
!     psim(:) =  0.50 * ctpoten * 1000.
!     psie(:) = -0.50 * ctpoten * 1000.
!     pcen(isouth) = 0.
!     pcen(inorth) = 0.
!     phidp0(:) = 90.*dtr
!     phidm0(:) = 90.*dtr
!     phinp0(:) = 90.*dtr
!     phinm0(:) = 90.*dtr
!     rr1(:) = -2.6
!
! tiegcm original with assymmetry and By effect
!     write(6,"(/,'byloc (local by)=',f8.3)") byloc

      offc(isouth) = 1.1*dtr
      offc(inorth) = 1.1*dtr
      dskofc(isouth) = (-0.08 + 0.15*byloc)*dtr
      dskofc(inorth) = (-0.08 - 0.15*byloc)*dtr
      phid(isouth) = (9.39 + 0.21*byloc - 12.) * h2deg * dtr
      phid(inorth) = (9.39 - 0.21*byloc - 12.) * h2deg * dtr
      phin(isouth) = (23.50 + 0.15*byloc - 12.) * h2deg * dtr
      phin(inorth) = (23.50 - 0.15*byloc - 12.) * h2deg * dtr
      psim(:) =  0.44 * ctpoten * 1000.
      psie(:) = -0.56 * ctpoten * 1000.
      pcen(isouth) = (-0.168 + 0.027*byloc) * ctpoten * 1000.
      pcen(inorth) = (-0.168 - 0.027*byloc) * ctpoten * 1000.
      phidp0(:) = 90.*dtr
      phidm0(:) = 90.*dtr
      phinp0(:) = 90.*dtr
      phinm0(:) = 90.*dtr
      rr1(:) = -2.6
!
! Recalculate some parameters if an amie run:
      if (iamie > 0) then
        theta0(:) = crad(:) ! radians from amie module
        phid(:) = phida(:)  ! radians from amie module
        hp_amie = max(hpi_sh_amie,hpi_nh_amie)
        hp_amie = 0.5*(hpi_sh_amie+hpi_nh_amie)
        plevel = 2.09*alog(hp_amie)
        h1 = min(2.35,0.83+0.33*plevel)
        h2 = 2.87 + 0.15 * plevel
        h0 = 0.5 * (h1 + h2) * dtr
        roth = (12.18- 0.89*plevel)
        rote = (2.62 - 0.55*plevel)
        rroth = roth * dtr
        rrote = rote * dtr
! For Sep 10 (Day 253) 2005 storm interval
        h2 = 4.5 + 0.11*hp_amie
! For October 2003 storm
!       h2 = 4.5 + 0.15*hp_amie
!       rradp(isouth) = crad(isouth)-1.5*dtr
!       rradp(inorth) = crad(inorth)-2.5*dtr
        rradp(isouth) = crad(isouth) + h2*dtr
        rradp(inorth) = crad(inorth) + h2*dtr
        if (rradp(isouth) > 31.*dtr) rradp(isouth) = 31.*dtr
        if (rradp(inorth) > 31.*dtr) rradp(inorth) = 31.*dtr
!       write(6,"('crad =',2f6.2,'  rrad=',2f6.2,'  rradp=',2f6.2,
!    |  ' power =',f5.1,' plev =',f5.2)")
!    |  crad/dtr,rrad/dtr,rradp/dtr,hp_amie,plevel
        rrad(isouth) = rradp(isouth)+0.5*dtr
        rrad(inorth) = rradp(inorth)+0.5*dtr
!
! Set critical colatitude crit(2) 40 deg -- G. Lu 5/11/98
! Therefore, crit2 = 40, crit1 = 25-30, depending on theta0
!
        c35 = min(30.,(theta0(1)+theta0(2))*0.5/dtr + 5.0)
!       c25 = max(25.,c35)
!       c25 = max(30.,c35)
        c25 = max(35.,c35)
        crit(1) = c25 * dtr
        crit(2) = 40. * dtr
        if (crit(1) > 35.*dtr) crit(1) = 35.*dtr
!       write(6,"('aurora: crit(1:2) = ',2f7.1)")crit(:)/dtr

      endif
!
! Report to stdout:
      if (iprint > 0) then
        write(6,"(/,'aurora_cons:')")
        write(6,"('  cusp:    alfac=',f8.3,' ec=',f8.3,' fc=',e12.4)")
     |    alfac,ec,fc
        write(6,"('  drizzle: alfad=',f8.3,' ed=',f8.3,' fd=',e12.4)")
     |    alfad,ed,fd
	write(6,"('  auroral radius = max of rhp,rcp=',2f10.3)") rhp,rcp
	write(6,"('  roth, rote (MLT) =',2f10.3)") roth,rote
	write(6,"('  1/e-widths = h1,h2=',2f10.3)") h1,h2
        write(6,"('  energy flux = e1,e2=',2f10.3)") e1,e2
        write(6,"('  add_sproton = ',l1)") add_sproton
        if (add_sproton) write(6,"('  solar protons: alfa_sp=',f10.3,
     |    ' e_sp = ',e10.3,' flx_sp = ',e10.3)") alfa_sp,e_sp,flx_sp
        if (add_helectron)
     |    write(6,"('  high-energy electrons: alfa30 = ',f10.3,
     |    ' e30 = ',f10.3)") alfa30,e30
        write(6,"(' ')")
      endif
      end subroutine aurora_cons
!-----------------------------------------------------------------------
      subroutine aurora(tn,o2,o1,n2,barm,lev0,lev1,lon0,lon1,lat)
!
! Auroral parameterization driver.
! (called by sub dynamics, from inside parallel latitude scan)
!
      use magfield_module,only: 
     |  rlatm,rlonm, ! magnetic latitude,longitude grid (radians)
     |  sunlons      ! sun's mag longitudes (see sub sunloc, magfield.F) 
      use params_module,only: nlat,zibot,dz
      use mpi_module,only: lat0,lat1
!
      use qrj_module,only:       ! for diag only
     |  qo2p,  ! o2+ ionization
     |  qop,   ! o+  ionization
     |  qn2p,  ! n2+ ionization
     |  qnp    ! n+  ionization
!
! Args:
      integer,intent(in) :: 
     |  lev0,lev1,  ! bottom,top of pressure column (1,nlevp1) 
     |  lon0,lon1,  ! first,last longitude indices of current domain
     |  lat         ! current latitude index 
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,    ! neutral gas temperature (deg K)
     |  o2,    ! O2 mass mixing ratio
     |  o1,    ! O mass mixing ratio
     |  n2,    ! N2 mass mixing ratio (1-o2-o-he)
     |  barm   ! mean molecular weight (see sub addiag)
!
! Local:
      integer :: i,ier
      integer,parameter :: nlond2 = max(1,nlon/2)
      real :: aurlat  ! boundary lat below which auroral is not calculated
                      ! (32.5 degrees at 5 degree resolution)
      real :: ofda,cosofa,sinofa,aslona
!
#ifdef VT
!     code = 130 ; state = 'aurora' ; activity='ModelCode'
      call vtbegin(130,ier)
#endif

!     write(6,"('Enter aurora: lat=',i4,' glat=',f8.2,' lon0,1=',2i4)")
!    |  lat,glat(lat),lon0,lon1
!
! Save ionization rates on entry to secondary history:
!     call addfld('QO2P_PRE',' ',' ',qo2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP_PRE' ,' ',' ',qop(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QN2P_PRE',' ',' ',qn2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Initialize diagnostics to insure they have values before the
! below return when lat < 32.5.
!
      cusp2d(:,lat) = 0.
      drzl2d(:,lat) = 0.
      alfa2d(:,lat) = 0.
      nflx2d(:,lat) = 0.
      eflux2d(:,lat) = 0.
!
! Check latitude, and return if below aurlat (32.5 deg at 5 deg res):
!
      aurlat = (float(lat-nlat)-0.5)*dphi+pi/2.
      do_aurora = .true.
      if (abs(aurlat) < pi/6.) do_aurora = .false.
      if (.not.do_aurora) then
#ifdef VT
        call vtend(130,ier)
#endif
        return
      endif
!
      do i=1,nlon
        dlat_aur(i) = rlatm(i+2,lat)
        dlon_aur(i) = rlonm(i+2,lat)-sunlons(lat)
      enddo
!
! Establish which hemisphere current latitude is in:
! (ihem will be either isouth (1) or inorth (2))
      ihem = int(dlat_aur(nlond2)*2./pi+2.)
!
! Find auroral circle coordinates (formerly in flowv3):
!
      ofda = sqrt(offa(ihem)**2+dskofa(ihem)**2)
      cosofa = cos(ofda)
      sinofa = sin(ofda)
      aslona = asin(dskofa(ihem)/ofda)
      do i=1,nlon
        sinlat(i) = sin(abs(dlat_aur(i)))
        coslat(i) = cos(dlat_aur(i))
        sinlon(i) = sin(dlon_aur(i)+aslona)
        coslon(i) = cos(dlon_aur(i)+aslona)
        colat(i) = acos(cosofa*sinlat(i)-sinofa*coslat(i)*coslon(i))
        alon(i) = amod(atan2(sinlon(i)*coslat(i),sinlat(i)*sinofa+
     |    cosofa*coslat(i)*coslon(i))-aslona+3.*pi,(pi*2.))-pi
      enddo
!
!     write(6,"('aurora: lat=',i2,' sinlat=',/,(6e12.4))") lat,sinlat
!     write(6,"('aurora: lat=',i2,' coslat=',/,(6e12.4))") lat,coslat
!     write(6,"('aurora: lat=',i2,' coslon=',/,(6e12.4))") lat,coslon
!     write(6,"('aurora: lat=',i2,' colat=',/,(6e12.4))") lat,colat
!     write(6,"('aurora: lat=',i2,' alon =',/,(6e12.4))") lat,alon
!
! Make cusp:
      call aurora_cusp(lat)
!
! Make alfa, flux, and drizzle:
      call aurora_heat(lat)

! Save aurora oval on exit to full arrays so can save later to secondary history
      do i=lon0,lon1
        cusp2d(i,lat)=cusp(i)
        drzl2d(i,lat)=drizl(i)
        alfa2d(i,lat)=alfa(i)
        nflx2d(i,lat)=flux(i)
        eflux2d(i,lat)=eflux(i)
      enddo
!
! Auroral additions to ionization rates:
! Sub aurora_ions is a standalone subroutine, i.e., not in aurora_module.
! Pass module data to aurora_ions through argument list.
!
!     subroutine aurora_ions(drizl,cusp,alfa1,alfa2,flux1,flux2,
!    |  tn,o2,o1,barm,zpbot,dzp,lev0,lev1,lon0,lon1,lat)
!
! 11/08 EMERY added eflux=flux*fac_p2e*alfa - (better for falfa1 if alfa quite variable)
      call aurora_ions(fac_p2e,
     |  drizl(lon0:lon1), cusp (lon0:lon1),
     |  alfa (lon0:lon1), alfa2(lon0:lon1), alfa3(lon0:lon1),
     |  eflux(lon0:lon1), eflux2(lon0:lon1), eflux3(lon0:lon1),
     |  flux (lon0:lon1), flux2(lon0:lon1), flux3(lon0:lon1),
     |  tn,o2,o1,n2,barm,zibot,dz,lev0,lev1,lon0,lon1,lat)
!
#ifdef VT
!     code = 130 ; state = 'aurora' ; activity='ModelCode'
      call vtend(130,ier)
#endif

      end subroutine aurora
!-----------------------------------------------------------------------
      subroutine aurora_cusp(lat)
      use params_module,only: nlonp2
      use cons_module,only: check_exp
!
! Calculate horizontal variation of polar cusp heating (was polht.F)
!
! Args:
      integer,intent(in) :: lat
!
! Local:
      integer :: i
      real,parameter :: 
     |  s5=.08726646, 
     |  s20=.34906585,
     |  pi_cusp = 3.1415927 ! hard-wired by tgcm15
!
! expo() (util.F) is used only if check_exp is true. This will avoid
! NaNS fpe, but will degrade performance. Check_exp is in cons.F.
!
      real,external :: expo
!
      do i=1,nlon
        if (.not.check_exp) then
          cusp(i+2) = (exp(-((theta0(ihem)-colat(i))/s5)**2)+
     |                 exp(-((pi_cusp-theta0(ihem)-colat(i))/s5)**2))*
     |                 exp(-(atan2(sin(alon(i)-phid(ihem)),
     |                             cos(alon(i)-phid(ihem)))/s20)**2)

!         write(6,"('cusp: lat=',i3,' i=',i3,' ihem=',i2,' theta0=',
!    |      e12.4,' colat=',e12.4,' alon=',e12.4,' phid=',e12.4,
!    |      ' cusp=',e12.4)") lat,i,ihem,theta0(ihem),colat(i),
!    |      alon(i),phid(ihem),cusp(i+2)

        else
          cusp(i+2) = (expo(-((theta0(ihem)-colat(i))/s5)**2,0)+
     |                expo(-((pi_cusp-theta0(ihem)-colat(i))/s5)**2,0))*
     |                 expo(-(atan2(sin(alon(i)-phid(ihem)),
     |                             cos(alon(i)-phid(ihem)))/s20)**2,0)
        endif
      enddo ! i=3,nlonp2
!
! Periodic points:
      do i=1,2
        cusp(i) = cusp(nlon+i)     ! 1:2 <- 73:74 at 5 deg res
        cusp(nlonp2+i) = cusp(i+2) ! 75:76 <- 3:4
      enddo
!     write(6,"('aurora_cusp: lat=',i2,' cusp=',/,(6e12.4))") 
!    |  lat,cusp
      end subroutine aurora_cusp 
!-----------------------------------------------------------------------
      subroutine aurora_heat(lat)
!
! Calculate alfa, flux, and drizzle.
!
      use params_module,only: nlon,nlonp2
#if defined(INTERCOMM) || defined(CISMAH)
      use cism_coupling_module,only: geng,gflx     !CISM
#endif
!
! Args:
      integer,intent(in) :: lat
!
! Local:
      integer :: i
      real,dimension(nlon) :: 
     |  coslamda,  ! cos(angle from throat)
     |  halfwidth, ! oval 1/e-width
     |  dtheta     ! latitudinal variation (Gaussian)
      real,parameter :: s10=0.174532925
      real :: clat

!
! Low-energy protons:
!     alfap0 = 0.5*(alfap1+alfap2)
!     e0p = 0.5*(pe1+pe2)
!
! coslamda  = cos(lamda)
! halfwidth = auroral half width
! dtheta    = colat-theta0(ihem)
! alfa      = electron energy
!
      if (alfa0 > 0.01) then
        do i=1,nlon
          coslamda(i) = cos(atan2(sin(alon(i)-rrote),
     |                            cos(alon(i)-rrote)))
!
! Auroral oval 1/e-width (equation (1) in Roble,1987):
          halfwidth(i)=h0*(1.-rh*cos(atan2(sin(alon(i)-rroth),
     |      cos(alon(i)-rroth))))
          dtheta(i) = colat(i)-rrad(ihem)
!
! Characteristic energy (equation (2) in Roble,1987):
          alfa(i+2) =  alfa0*(1.-ralfa*coslamda(i))
        enddo ! i=1,nlon
      else
        do i=1,nlon
          coslamda(i)=cos(atan2(sin(alon(i)-rrote),cos(alon(i)-rrote))) ! wk1
          halfwidth(i)=h0*(1.-rh*cos(atan2(sin(alon(i)-rroth),          ! wk2
     |      cos(alon(i)-rroth))))
          dtheta(i) = colat(i)-rrad(ihem)                               ! wk3
!
! 1/02 bf: this alfa not currently in use:
!         alfa(i+2) = alfk+
!    |      alf6*exp(-((dtheta(i)-rd6-rd6v*cos(alon(i)) )/rh6)**2)*
!    |      exp(-(atan2(sin(alon(i)-rrot6),cos(alon(i)-rrot6))/rt6)**2)+
!    |      alf21*exp(-((dtheta(i)-rd21)/rh21)**2)*
!    |      exp(-(atan2(sin(alon(i)-rrot21),cos(alon(i)-rrot21))/
!    |      rt21)**2)
        enddo ! i=1,nlon
      endif
!     write(6,"('aurora_heat: lat=',i3,' alfa=',/,(6e12.4))")
!    |  lat,alfa
!
! Flux, drizzle, alfa2, flux2:
      do i=1,nlon
!
! Energy flux (equation (3) in Roble,1987):
        flux(i+2) = e0*(1.-ree*coslamda(i))*exp(-(dtheta(i)/
     |    halfwidth(i))**2) / (2.*alfa(i+2)*1.602e-9)
! 11/08 EMERY added eflux=flux*fac_p2e*alfa
	eflux(i+2) = flux(i+2)*fac_p2e*alfa(i+2)
        drizl(i+2) = exp(-((dtheta(i)+abs(dtheta(i)))/(2.*h0))**2)
        alfa2(i+2) = alfa20*(1.-ralfa2*coslamda(i))
        flux2(i+2) = e20*(1.-re2*coslamda(i))*exp(-(dtheta(i)/
     |    halfwidth(i))**2) / (2.*alfa2(i+2)*1.602E-9)
	eflux2(i+2) = flux2(i+2)*fac_p2e*alfa2(i+2)
!
! Alfa3, flux3 for high energy electrons:
        alfa3(i+2) = alfa30
        flux3(i+2) = e30*exp(-(dtheta(i)/halfwidth(i))**2)/1.602e-6
	eflux3(i+2) = flux3(i+2)*fac_p2e*alfa3(i+2)
!
! For electron temperature (used in settei):  
        qteaur(i+2,lat) = -7.e+8*exp(-(dtheta(i)/halfwidth(i))**2)
!
! Recalculate alfa, flux, drizl if an amie run (as in old sub amiepa):
        if (iamie > 0) then
! Insure ekvg mean energy >= 1.:
!         if (lat==5) write(6,"('Using AMIE aurora flux >>>',4g10.3)")
!    |      minval(ekvg),maxval(ekvg),minval(efxg),maxval(efxg)
          ekvg(i,lat) = max(ekvg(i,lat),1.)
          alfa(i+2) = ekvg(i,lat)/2.
          alfa2(i+2) = alfa20
          clat = acos(sin(abs(dlat_aur(i))))
!         drizl(i+2) = exp(-((clat-crad(ihem)+
!    |      abs(clat-crad(ihem)))/(2.*s10))**2)
!         drizl(i+2) = exp(-((colat(i)-crad(ihem)+
!    |      abs(colat(i)-crad(ihem)))/s10)**2)
          flux(i+2) = max(efxg(i,lat)/(2.*alfa(i+2)*1.602e-9),1.e-20)
!          flux2(i+2) = e20*(1.-re2*coslamda(i))*
!     |      exp(-((clat-crad(ihem))/
!!    |      exp(-((colat(i)-crad(ihem))/
!     |      halfwidth(i))**2) / (2.*alfa2(i+2)*1.602E-9)
! New drizal formulation - G. Lu /08/17/2006
          drizl(i+2) = exp(-((clat-rradp(ihem)+
     |      abs(clat-rradp(ihem)))/(1*s10))**2)
          flux2(i+2) = e20*(1.-re2*coslamda(i))*
     |      exp(-((clat-rrad(ihem))/
     |      halfwidth(i))**2) / (2.*alfa2(i+2)*1.602E-9)
        endif

      enddo ! i=1,nlon
!
! Periodic points:
      do i=1,2
        alfa  (i) = alfa (nlon+i)
        flux  (i) = flux (nlon+i)
	eflux (i) = eflux(nlon+i)
        drizl (i) = drizl(nlon+i)
        alfa2 (i) = alfa2(nlon+i)
        alfa3 (i) = alfa3(nlon+i)
        flux2 (i) = flux2(nlon+i)
        flux3 (i) = flux3(nlon+i)
	eflux2 (i) = eflux2(nlon+i)
	eflux3 (i) = eflux3(nlon+i)
        qteaur(i,lat)= qteaur(nlon+i,lat)
        alfa  (nlonp2+i) = alfa (i+2)
        flux  (nlonp2+i) = flux (i+2)
	eflux (nlonp2+i) = eflux(i+2)
        drizl (nlonp2+i) = drizl(i+2)
        alfa2 (nlonp2+i) = alfa2(i+2)
        alfa3 (nlonp2+i) = alfa3(i+2)
	flux2 (nlonp2+i) = flux2(i+2)
        flux3 (nlonp2+i) = flux3(i+2)
	eflux2 (nlonp2+i) = eflux2(i+2)
	eflux3 (nlonp2+i) = eflux3(i+2)
        qteaur(nlonp2+i,lat)= qteaur(i+2,lat)
      enddo
!
! Replacing with M-I coupler values. Warning: alfa and flux can not be
! zero, otherwise the model will blow up. Set to small values here. This
! will not change simulations.
!                                                        Wenbin  wang
!
#if defined(INTERCOMM) || defined(CISMAH)
     do i=1,nlonp4
      alfa(i)=geng(i,lat)+0.01
      flux(i)=gflx(i,lat)+0.01
      if(alfa(i)<=0.)alfa(i)=0.01
      if(flux(i)<=0.)flux(i)=0.01
      eflux(i) = flux(i)*fac_p2e*alfa(i)
!     write(6,*)alfa(i),flux(i),i,lat
     enddo
#endif

      end subroutine aurora_heat
      end module aurora_module
!-----------------------------------------------------------------------
! 11/08 EMERY added eflux=flux*fac_p2e*alfa (better when alfa quite variable)
      subroutine aurora_ions (fac_p2e,drizl,cusp,alfa1,alfa2,alfa3,
     |  eflux1,eflux2,eflux3,flux1,flux2,flux3,
     |  tn,o2,o1,n2,barm,zpbot,dzp,lev0,lev1,lon0,lon1,lat)
!
! Calculate aroral additions to ionization rates.
! Note this routine is standalone, not in aurora_module.
!
      use cons_module,only: ! see cons.F
     |  p0,    ! standard pressure (5.0e-4)
     |  expz,  ! exp(-zp) where zp is log pressure (ln(p0/p))
     |  grav,  ! gravity
     |  gask,  ! gas constant (8.314e7)
     |  boltz, ! boltzman's constant (1.38e-16)
     |  rmass_o2, rmass_o1, rmass_n2,        ! molelcular weights of o2,o,n2
     |  rmassinv_o2,rmassinv_o1,rmassinv_n2, ! 1./rmass
     |  avo    ! avogadro number (6.023e23)
!
! Total ionizations. Auroral ionization is added to these:
      use qrj_module,only:  ! (nlevp1,nlonp4) see qrj.F
     |  qo2p,  ! o2+ ionization
     |  qop,   ! o+  ionization
     |  qn2p,  ! n2+ ionization
     |  qnp,   ! n+  ionization
     |  qtef   ! n  production rate (cm^-3s-1)
!
! Auroral parameters:
      use aurora_module,only: alfac, alfad, alfa_sp, ec,ed, fc,fd,
     |  flx_sp, e30, add_sproton, add_helectron
      use addfld_module,only: addfld
      use init_module,only: istep
      implicit none
!
! Args:
      real,intent(in) :: fac_p2e                ! convert from particle to energy flux
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lon0:lon1) :: drizl,cusp,
     |  alfa1,alfa2,alfa3,eflux1,eflux2,eflux3,flux1,flux2,flux3
! 11/08 EMERY added eflux=flux*fac_p2e*alfa as eflux1
      real,intent(in) :: zpbot,dzp
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,o2,o1,n2,barm
!
! Local:
      integer :: i,k,nlevs
      real :: zlev,klev,xn2
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  p0ez, xalfa1, xalfa2, xcusp, xdrizl,       ! input to sub aion
     |  xalfa_sp, xalfa3,
     |  flux1_ion, flux2_ion, cusp_ion, drizl_ion, ! output from sub aion
     |  alfa1_ion, alfa2_ion, alfa3_ion,           ! output from sub aion
     |  alfasp_bion,                               ! output from sub bion
     |  barm_t, qsum, denom, p0ez_mbar, tk_mbar,
     |  qo2p_aur, qop_aur, qn2p_aur ! auroral ionization for O2+, O+, N2+
      real,dimension(lon0:lon1) ::
     |  falfa1, falfa2, fcusp, fdrizl, falfa_sp, falfa3
      real :: qia(5) ! low energy proton source (not in use, 1/02)
      real,dimension(lev0:lev1,lon0:lon1) :: qaurora
! External:
      logical,external :: time2print
      logical :: use_cion = .false.      ! switch between AION and CION
!
      qia(:) = 0.
!
!     write(6,"('aurora_ions: lat=',i2,' flux1=',/,(10e12.4))")
!    |  lat,flux1
!     write(6,"('aurora_ions: lat=',i2,' alfa1=',/,(10e12.4))")
!    |  lat,alfa1
!     write(6,"('aurora_ions: lat=',i2,' alfa2=',/,(6e12.4))")
!    |  lat,alfa2
!     write(6,"('aurora_ions: lat=',i2,' cusp =',/,(6e12.4))")
!    |  lat,cusp 
!     write(6,"('aurora_ions: lat=',i2,' drizl=',/,(6e12.4))")
!    |  lat,drizl 
!
      do i=lon0,lon1
        do k=lev0,lev1
          p0ez(k,i) = (p0*expz(k)/(grav*4.e-6))**0.606
          xalfa1(k,i) = p0ez(k,i)/alfa1(i)
          xalfa2(k,i) = p0ez(k,i)/alfa2(i)
          xalfa3(k,i) = p0ez(k,i)/alfa3(i) ! alfa3(:)==alfa30
          xcusp (k,i) = p0ez(k,i)/alfac
          xdrizl(k,i) = p0ez(k,i)/alfad
!
          p0ez_mbar(k,i) = p0*expz(k)/(boltz*tn(k,i)*
     |      (o2(k,i)*rmassinv_o2+
     |       o1(k,i)*rmassinv_o1+
     |       n2(k,i)*rmassinv_n2))/avo
          tk_mbar(k,i) = gask*tn(k,i)/(grav*barm(k,i))*p0ez_mbar(k,i)
          xalfa_sp(k,i) = ((tk_mbar(k,i)/0.00271)**0.58140)/alfa_sp
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
! Init (whole array operations):
      flux1_ion = 1.e-20
      flux2_ion = 1.e-20
      alfa1_ion = 1.e-20
      alfa2_ion = 1.e-20
      cusp_ion  = 1.e-20
      drizl_ion = 1.e-20
!
! Sub aion calculates auroral electrons: 
!
!     zlev = -7.5 ! time-gcm only
!     klev = int((zlev-zpbot)/dzp)+1
!     write(6,"('zlev=',f8.2,' klev=',i3)") zlev,klev
!
      if (.not.use_cion) then
        call aion(xalfa1,alfa1_ion,lev0,lev1,lon0,lon1,lat)
        call aion(xalfa2,alfa2_ion,lev0,lev1,lon0,lon1,lat)
        call aion(xcusp,cusp_ion,lev0,lev1,lon0,lon1,lat)
        call aion(xdrizl,drizl_ion,lev0,lev1,lon0,lon1,lat)
      else
        call cion(xalfa1,alfa1_ion,alfa1,lev0,lev1,lon0,lon1,lat)
        call cion(xalfa2,alfa2_ion,alfa1,lev0,lev1,lon0,lon1,lat)
        call cion(xcusp,cusp_ion,alfa1,lev0,lev1,lon0,lon1,lat)
        call cion(xdrizl,drizl_ion,alfa1,lev0,lev1,lon0,lon1,lat)
      endif

!     call addfld('ALFA1',' ',' ',alfa1_ion,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('ALFA2',' ',' ',alfa2_ion,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('CUSP',' ',' ',cusp_ion,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('DRIZL',' ',' ',drizl_ion,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!
! Solar protons:
      if (flx_sp > 1.e-19)
     |  call bion(xalfa_sp,alfasp_bion,lev0,lev1,lon0,lon1,lat)
!
! High energy electrons:
      if (e30 > 1.e-19)
     |  call aion(xalfa3, alfa3_ion,lev0,lev1,lon0,lon1,lat)
!
      qaurora = 0.
      do i=lon0,lon1
! 11/08 EMERY added eflux=flux*fac_p2e*alfa;  eflux1 better than flux1 if alfa quite variable
!       falfa1(i) = eflux1(i)/fac_p2e  ! s7  This should be equivalent to alfa1*flux1
        falfa1(i) = alfa1(i)*flux1(i)  ! s7
!       falfa2(i) = eflux2(i)/fac_p2e  ! s8  This should be equivalent to alfa2*flux2
        falfa2(i) = alfa2(i)*flux2(i)  ! s8

        fcusp (i) = cusp(i)*alfac*fc   ! s9
        fdrizl(i) = drizl(i)*alfad*fd  ! s10
!	falfa3(i) = eflux3(i)/fac_p2e  ! s13 This should be equivalent to alfa3*flux3
        falfa3(i) = alfa3(i)*flux3(i)  ! s13 (high energy electrons)
!
! Pressure column, bottom to top:
        do k=lev0,lev1-1
          falfa_sp(i) = drizl(i)*p0ez_mbar(k,i)*flx_sp*1.e6/
     |      (tk_mbar(k,i)*35.)
          barm_t(k,i) = grav*0.5*(barm(k,i)+barm(k+1,i))/(35.e-3*gask* ! s11
     |      tn(k,i))
          qsum(k,i) = 
     |      falfa1(i)*alfa1_ion(k,i)+  ! s7*s3
     |      falfa2(i)*alfa2_ion(k,i)+  ! s8*s4
     |      fcusp (i)*cusp_ion (k,i)+  ! s9*s5
     |      fdrizl(i)*drizl_ion(k,i)   ! s10*s6
          qaurora(k,i) = barm_t(k,i)*(falfa1(i)*alfa1_ion(k,i)+
     |                   falfa2(i)*alfa2_ion(k,i))
!
! Include solar protons if add_sproton is set, and high energy electrons
! if add_helectron is set:
          if (add_helectron) 
     |      qsum(k,i) = qsum(k,i) + falfa3(i)*alfa3_ion(k,i)
          if (add_sproton)
     |      qsum(k,i) = qsum(k,i) + falfa_sp(i)*alfasp_bion(k,i)
!
          qsum(k,i) = qsum(k,i)*barm_t(k,i) ! s1 = s1*s11
!
! Denominator of equations (13-16) in Roble,1987.
          denom(k,i) = 
     |      0.92 * n2(k,i)*rmassinv_n2 +
     |      1.5  * o2(k,i)*rmassinv_o2 + 
     |      0.56 * o1(k,i)*rmassinv_o1
!
! Production of O2+ (equation (15) in Roble,1987):
          qo2p_aur(k,i) = qsum(k,i)*o2(k,i)/(rmass_o2*denom(k,i))+qia(2)
!
! Production of O+ (equation (16) in Roble,1987):
          qop_aur(k,i) = qsum(k,i)*
     |      (0.5 *o2(k,i)*rmassinv_o2+
     |       0.56*o1(k,i)*rmassinv_o1)/denom(k,i) + qia(3)
!
! Production of N2+ (equation (13) in Roble,1987)
          xn2 = n2(k,i)
          if (xn2 <= 1.e-8) xn2 = 1.e-8
          qn2p_aur(k,i) = qsum(k,i)*0.7*xn2/(rmass_n2*denom(k,i))+qia(1)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     write(6,"('aurora: lat=',i3,' lon0,1=',2i3,' falfa1=',/,
!    |  (10e12.4))") lat,lon0,lon1,falfa1
!     write(6,"('aurora: lat=',i3,' lon0,1=',2i3,' qsum=',/,
!    |  (8e12.4))") lat,lon0,lon1,qsum
!     write(6,"('aurora: lat=',i3,' lon0,1=',2i3,' qo2p_aur=',/,
!    |  (9e12.4))") lat,lon0,lon1,qo2p_aur
!     write(6,"('aurora: lat=',i3,' lon0,1=',2i3,' qop_aur=',/,
!    |  (9e12.4))") lat,lon0,lon1,qop_aur
!     write(6,"('aurora: lat=',i3,' lon0,1=',2i3,' qn2p_aur=',/,
!    |  (9e12.4))") lat,lon0,lon1,qn2p_aur
!     write(6,"('aurora: lat=',i3,' lon0,1=',2i3,' falfa2=',/,
!    |  (6e12.4))") lat,lon0,lon1,falfa2
!     write(6,"('aurora: lat=',i3,' lon0,1=',2i3,' fcusp=',/,
!    |  (6e12.4))") lat,lon0,lon1,fcusp
!     write(6,"('aurora: lat=',i3,' lon0,1=',2i3,' fdrizl=',/,
!    |  (6e12.4))") lat,lon0,lon1,fdrizl

!     call addfld('BARMT'   ,' ',' ',barm_t(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('QSUM'    ,' ',' ',qsum    ,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('DENOM'   ,' ',' ',denom(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('QO2P_AUR',' ',' ',qo2p_aur(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('QOP_AUR' ,' ',' ',qop_aur(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('QN2P_AUR',' ',' ',qn2p_aur(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      call addfld('QAURORA',' ','cm-3 s-1',qaurora,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Add auroral to total ionization rates for o2+, o+, and n2+
! (see qrj_module)
      do i=lon0,lon1
        do k=lev0+1,lev1-1
          qo2p(k,i,lat) = qo2p(k,i,lat)+sqrt(qo2p_aur(k,i)*
     |      qo2p_aur(k-1,i))
          qop (k,i,lat) = qop (k,i,lat)+sqrt(qop_aur (k,i)*
     |      qop_aur(k-1,i))
          qn2p(k,i,lat) = qn2p(k,i,lat)+sqrt(qn2p_aur(k,i)*
     |      qn2p_aur(k-1,i))
          qnp (k,i,lat) = qnp (k,i,lat)+
     |      .22/.7 * sqrt(qn2p_aur(k,i)*qn2p_aur(k-1,i))
          qtef(k,i,lat) = qtef(k,i,lat)+1.57*sqrt(qn2p_aur(k,i)*
     |      qn2p_aur(k-1,i))
        enddo ! k=lev0+2,lev1-1
!
! Extrapolate lower boundary:
! TIEGCM-1.9x bug fix from Wenbin's old LTR/CMIT-2.5 code
#if defined(INTERCOMM) || defined(CISMAH)
        if(1.5*qo2p_aur(lev0,i) > 0.5*qo2p_aur(lev0+1,i))then
         qo2p(lev0,i,lat) = qo2p(lev0,i,lat)+1.5*qo2p_aur(lev0,i)-
     |     0.5*qo2p_aur(lev0+1,i)
        endif                         !!! CISM fix
        if(1.5*qop_aur (lev0,i) > 0.5*qop_aur (lev0+1,i))then
         qop (lev0,i,lat) = qop (lev0,i,lat)+1.5*qop_aur (lev0,i)-
     |     0.5*qop_aur (lev0+1,i)
        endif
        if(1.5*qn2p_aur(lev0,i) > 0.5*qn2p_aur(lev0+1,i))then
         qn2p(lev0,i,lat) = qn2p(lev0,i,lat)+1.5*qn2p_aur(lev0,i)-
     |     0.5*qn2p_aur(lev0+1,i)
        endif
        if(1.5*qn2p_aur(lev0,i) > 0.5*qn2p_aur(lev0+1,i))then
         qnp (lev0,i,lat) = qnp (lev0,i,lat)+ .22/.7 *
     |     (1.5*qn2p_aur(lev0,i)-0.5*qn2p_aur(lev0+1,i))
        endif
        if(1.5*qn2p_aur(lev0,i) > 0.5*qn2p_aur(lev0+1,i))then
         qtef(lev0,i,lat) = qtef(lev0,i,lat)+1.57*(1.5*qn2p_aur(lev0,i)-
     |     0.5*qn2p_aur(lev0+1,i))
        endif
#else
        qo2p(lev0,i,lat) = qo2p(lev0,i,lat)+1.5*qo2p_aur(lev0,i)-
     |    0.5*qo2p_aur(lev0+1,i)        
        qop (lev0,i,lat) = qop (lev0,i,lat)+1.5*qop_aur (lev0,i)-
     |    0.5*qop_aur (lev0+1,i)        
        qn2p(lev0,i,lat) = qn2p(lev0,i,lat)+1.5*qn2p_aur(lev0,i)-
     |    0.5*qn2p_aur(lev0+1,i)        
        qnp (lev0,i,lat) = qnp (lev0,i,lat)+ .22/.7 * 
     |    (1.5*qn2p_aur(lev0,i)-0.5*qn2p_aur(lev0+1,i))
        qtef(lev0,i,lat) = qtef(lev0,i,lat)+1.57*(1.5*qn2p_aur(lev0,i)-
     |    0.5*qn2p_aur(lev0+1,i))        
#endif
!
! Extrapolate upper boundary:
        qo2p(lev1,i,lat)= qo2p(lev1,i,lat)+1.5*qo2p_aur(lev1-1,i)-
     |    0.5*qo2p_aur(lev1-2,i)
        qop (lev1,i,lat)= qop (lev1,i,lat)+1.5*qop_aur (lev1-1,i)-
     |    0.5*qop_aur (lev1-2,i)
        qn2p(lev1,i,lat)= qn2p(lev1,i,lat)+1.5*qn2p_aur(lev1-1,i)-
     |    0.5*qn2p_aur(lev1-2,i)
        qnp (lev1,i,lat)= qnp (lev1,i,lat)+ .22/.7 * 
     |    (1.5*qn2p_aur(lev1-1,i)-0.5*qn2p_aur(lev1-2,i))
        qtef(lev1,i,lat) = qtef(lev1,i,lat)+1.57*
     |    (1.5*qn2p_aur(lev1-1,i)-0.5*qn2p_aur(lev1-2,i))
      enddo ! i=lon0,lon1

!     call addfld('QO2P',' ',' ',qo2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QOP' ,' ',' ',qop(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QN2P' ,' ',' ',qn2p(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QNP' ,' ',' ',qnp(lev0:lev1,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      end subroutine aurora_ions
!-----------------------------------------------------------------------
      subroutine aion(si,so,lev0,lev1,lon0,lon1,lat)
      implicit none
!
! Calculates integrated f(x) needed for total auroral ionization.
! See equations (10-12) in Roble,1987.
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in)  :: si(lev0:lev1,lon0:lon1)
      real,intent(out) :: so(lev0:lev1,lon0:lon1)
!
! Local:
      integer :: k,i
      real :: xlog
!
! Coefficients for equation (12) of Roble,1987, giving total 
! ionization rate (revised since 1987):
      real,parameter :: cc(8) =
     | (/3.2333134511131 ,  2.5658873458085 ,  2.2540957232641 ,
     |   0.72971983372673,  1.1069072431948 ,  1.7134937681128 ,
     |   1.8835442312993 ,  0.86472135072090/)
!
! Use the identity x**y = exp(y*ln(x)) for performance 
! (fewer (1/2) trancendental functions are required).
!
      do k=lev0,lev1-1
        do i=lon0,lon1
          xlog = log(si(k,i))
          so(k,i) = cc(1)*exp(cc(2)*xlog-cc(3)*exp(cc(4)*xlog))+
     |              cc(5)*exp(cc(6)*xlog-cc(7)*exp(cc(8)*xlog))
          enddo ! i=lon0,lon1
      enddo ! k=lev0,lev1
      end subroutine aion
!-----------------------------------------------------------------------
      subroutine bion(si,so,lev0,lev1,lon0,lon1,lat)
      implicit none
!
! Calculates integrated f(x) needed for total auroral ionization.
! See equations (10-12) in Roble,1987.
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in)  :: si(lev0:lev1,lon0:lon1)
      real,intent(out) :: so(lev0:lev1,lon0:lon1)
!
! Local:
      integer :: k,i
      real :: xlog
      real,parameter :: cc(8) =
     |  (/0.12718, 4.9119, 1.8429, 0.99336, 0.52472, 1.5565,
     |    0.85732, 1.4116/)
!
! Use the identity x**y = exp(y*ln(x)) for performance 
! (fewer (1/2) trancendental functions are required).
!
      do k=lev0,lev1-1
        do i=lon0,lon1
          xlog = log(si(k,i))
          so(k,i) = cc(1)*exp(cc(2)*xlog-cc(3)*exp(cc(4)*xlog))+
     |              cc(5)*exp(cc(6)*xlog-cc(7)*exp(cc(8)*xlog))
          enddo ! i=lon0,lon1
      enddo ! k=lev0,lev1
      end subroutine bion
!-----------------------------------------------------------------------
      subroutine cion(si,so,alpha,lev0,lev1,lon0,lon1,lat)
!
! Calculate a normalized energy deposition function, which is produced by
! electron precipitation with a characteristic energy from 100 eV to 1 MeV.
! This parameterization is designed to improve the subroutine aion, extending
! the working energy range from around 1-100 keV to 0.1-1000 keV.
! -- Xiaohua Fang, LASP, University of Colorado, May 2008.
! Reference:
!  Fang, X., C. E. Randall, D. Lummerzheim, S. C. Solomon, M. J. Mills,
!  D. Marsh, C. H. Jackman, W. Wang, and G. Lu (2008), Electron impact
!  ionization: A new parameterization for 100 eV to 1 MeV electrons,
!  J. Geophys. Res., doi:10.1029/2008JA013384.
!
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in)  :: si(lev0:lev1,lon0:lon1),alpha(lon0:lon1)
      real,intent(out) :: so(lev0:lev1,lon0:lon1)
!
! Locals:
      real,parameter :: P(4, 8) = RESHAPE(
     | (/ 3.49979e-1, -6.18200e-2, -4.08124e-2,  1.65414e-2,
     |    5.85425e-1, -5.00793e-2,  5.69309e-2, -4.02491e-3,
     |    1.69692e-1, -2.58981e-2,  1.96822e-2,  1.20505e-3,
     |   -1.22271e-1, -1.15532e-2,  5.37951e-6,  1.20189e-3,
     |    1.57018,     2.87896e-1, -4.14857e-1,  5.18158e-2,
     |    8.83195e-1,  4.31402e-2, -8.33599e-2,  1.02515e-2,
     |    1.90953,    -4.74704e-2, -1.80200e-1,  2.46652e-2,
     |   -1.29566,    -2.10952e-1,  2.73106e-1, -2.92752e-2 /),
     | (/4, 8/) )
      real :: cc(8), logE, logY
      integer  :: i,ii,k
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          if (alpha(i) >= 0.1) then
!
! calculate the energy-dependent coefficients
            logE = log(alpha(i))      ! alpha=E0 (in keV)
            do ii=1, 8
              cc(ii)=exp( P(1,ii) + P(2,ii)*logE + P(3,ii)*(logE**2.) +
     |          P(4,ii)*(logE**3.) )
            end do
!
! calculate the energy deposition function
            logY = log( si(k,i) )
            so(k,i) =  cc(1)*exp( cc(2)*logY - cc(3)*exp( cc(4)*logY ) )
     |               + cc(5)*exp( cc(6)*logY - cc(7)*exp( cc(8)*logY ) )
          else
            so(k,i)=0.
          endif
        enddo
      enddo
      end subroutine cion

