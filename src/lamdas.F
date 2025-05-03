!
      subroutine lamdas(tn,xnmbar,o2,o1,he,n2,ti,te,o2p,op,nplus,
     |  n2p,nop,Etot,lxx,lyy,lxy,lyx,lamda1,ped_out,hall_out,Q1,Q2,KncP,
     |  KncH,lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Compute ion drag coefficients 
!      in geographic direction lxx, lyy, lxy, and lyx [1/s] (full levels)
!      in magnetic direction lamda1 [1/s] (full levels)
! Pedersen/ Hall conductivities [S/m] (half levels)
!
      use magfield_module,only: bmod2,sn2dec,csdec,sndec
      use cons_module,only: dipmin,avo,rtd,boltz,
     |  rmass_o2,rmass_o1,rmass_he,rmass_n2,rmass_no,rmass_n4s,
     |  rmassinv_o2,rmassinv_o1,rmassinv_he,rmassinv_n2,
     |  rmassinv_no,rmassinv_n4s
      use magfield_module,only: dipmag,rlatm
      use input_module,only: colfac
      use addfld_module,only: addfld
      use diags_module,only: mkdiag_SIGMAPED,mkdiag_SIGMAHAL,
     |                       mkdiag_LAMDAPED,mkdiag_LAMDAHAL
      implicit none
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,   ! neutral temperature (deg K)
     |  xnmbar,! for conversion to volume density
     |  o2,   ! molecular oxygen (mmr)
     |  o1,   ! atomic oxygen    (mmr)
     |  he,   ! helium           (mmr)
     |  n2,   ! molecular nitrogen (mmr)
     |  ti,   ! ion temperature (deg K)
     |  te,   ! electron temperature (deg K)
     |  o2p,  ! O2+ number density (1/cm3)
     |  op,   ! O+  number density (1/cm3)
     |  nplus,! N+ number density (1/cm3)
     |  n2p,  ! N2+ number density (1/cm3)
     |  nop,  ! NO+ number density (1/cm3)
     |  Etot  ! Total electric field
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  lxx,     ! lamda XX term [1/s] 
     |  lyy,     ! lamda YY term [1/s] 
     |  lxy,     ! lamda XY term [1/s] 
     |  lyx,     ! lamda YX term [1/s] 
     |  lamda1,  ! lamda_1 (sigma_p * B^2) / rho [1/s] magnetic eastward direction
     |  ped_out, ! pedersen conductivity (siemens/m) (will be input to dynamo)
     |  hall_out,! hall conductivity (siemens/m) (will be input to dynamo)
     |  KncP,    ! Nonlinear current coefficient for Pedersen conductacne
     |  KncH,    ! Nonlinear current coefficient for Hall conductacne
     |  Q1,      ! electron ohmic heating
     |  Q2       ! electrojet turbulent heating
!
! Local:
      real,parameter ::
     |  qe = 1.602e-19,         ! electron charge (coulomb)
     |  me = qe/1.7588028E7/10, ! electron mass [g]
     |  mp = qe/9.6489E3/10,    ! proton mass [g]
!
! Ion-neutral momentum transfer collision frequencies (non-resonant):
! 2024/03/21 Haonan Wu: There used to be only O2+, O+, and NO+ in the calculation.
! To be more self-consistent, include N+ and N2+ in the calculation as well.
!
     |  nu_op_o2  = 6.64E-10, ! O+  ~ O2
     |  nu_np_o2  = 7.25E-10, ! N+  ~ O2
     |  nu_n2p_o2 = 4.49E-10, ! N2+ ~ O2
     |  nu_nop_o2 = 4.27E-10, ! NO+ ~ O2
     |  nu_o2p_o  = 2.31E-10, ! O2+ ~ O
     |  nu_np_o   = 4.42E-10, ! N+  ~ O
     |  nu_n2p_o  = 4.42E-10, ! N2+ ~ O
     |  nu_nop_o  = 2.44E-10, ! NO+ ~ O
     |  nu_o2p_he = 0.70E-10, ! O2+ ~ He
     |  nu_op_he  = 1.32E-10, ! O+  ~ He
     |  nu_np_he  = 1.49E-10, ! N+  ~ He
     |  nu_n2p_he = 0.79E-10, ! N2+ ~ He
     |  nu_nop_he = 0.74E-10, ! NO+ ~ He
     |  nu_o2p_n2 = 4.13E-10, ! O2+ ~ N2
     |  nu_op_n2  = 6.82E-10, ! O+  ~ N2
     |  nu_np_n2  = 7.47E-10, ! N+  ~ N2
     |  nu_nop_n2 = 4.34E-10  ! NO+ ~ N2
!
      integer :: k,i
      real :: sqrt_te,        ! sqrt(te)
     |  E1,E2                 ! turbulent eletric field
!
! Local (lon):
      real,dimension(lon0:lon1) :: 
     |  dip,rot_xx,rot_yy,rot_xy,
     |  omega_o2p,omega_op ,omega_np,
     |  omega_n2p,omega_nop,omega_e

! reference pedersen conductivity (siemens/m) without ET effects.
! The name of _ref may be in other studies such as eclipse, flare etc.
! size should be consistent with ped_out.
      real,dimension(lev0:lev1,lon0-2:lon1+2) :: ped_ref_out
      real,dimension(lev0:lev1,lon0-2:lon1+2) :: hall_ref_out
!
! Local (lev,lon):
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  tr,         ! average of tn and ti
     |  o2_cm3, o1_cm3, he_cm3, n2_cm3, ! major species number densities (cm3)
     |  sigma_ped_ref,  ! pedersen conductivity without ET effects (siemens/m)
     |  sigma_ped,  ! pedersen conductivity (siemens/m)
     |  sigma_hall, ! hall conductivity (siemens/m)
     |  ne,         ! electron density (assume o2p+op+nplus+n2p+nop)
     |  lamda2,     ! sighal*b**2/rho
     |  lamda1tmp,  ! temporary lamda1
     |  lamda2tmp,  ! temporary lamda2
!
! Ion-neutral momentum transfer collision frequencies:
! (resonant, temperature dependent)
!
     |  nu_o2p_o2, ! O2+ ~ O2
     |  nu_op_o,   ! O+  ~ O
     |  nu_n2p_n2, ! N2+ ~ N2
!
     |  nu_o2p,    ! [o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~he]n(he)+[o2p~n2]n(n2)
     |  nu_op,     ! [op ~o2]n(o2)+[op ~o]n(o)+[op ~he]n(he)+[op ~n2]n(n2)
     |  nu_np,     ! [np ~o2]n(o2)+[np ~o]n(o)+[np ~he]n(he)+[np ~n2]n(n2)
     |  nu_n2p,    ! [n2p~o2]n(o2)+[n2p~o]n(o)+[n2p~he]n(he)+[n2p~n2]n(n2)
     |  nu_nop,    ! [nop~o2]n(o2)+[nop~o]n(o)+[nop~he]n(he)+[nop~n2]n(n2)
! 
     |  nu_e,      ! electron~neutral
     |  nu_i,      ! ion~neutral
     |  mi,        ! mean ion molecular mass
     |  ke,ki      ! ratio between gyro frequency and neutral collision frequency
!
! Set local needs:
      do i=lon0,lon1
!
!  e/B  [C/T 10^6 cm^3/m^3]
! 1.e10 for SI units = 1.e6 (cm3->m3) * 1.e4 (gauss->tesla)
!
! gyrofrequencies: omega_i = eB/m_i  [1/s]
!                  omega_e = eB/m_e  [1/s]
! 1/10 for conversion from Gauss/g to T/kg
!
        omega_op (i) = qe*bmod2(i,lat)*rmassinv_o1 /(mp*10)
        omega_o2p(i) = qe*bmod2(i,lat)*rmassinv_o2 /(mp*10)
        omega_np (i) = qe*bmod2(i,lat)*rmassinv_n4s/(mp*10)
        omega_n2p(i) = qe*bmod2(i,lat)*rmassinv_n2 /(mp*10)
        omega_nop(i) = qe*bmod2(i,lat)*rmassinv_no /(mp*10)
        omega_e(i)   = qe*bmod2(i,lat)/(me*10)
!
! Sin and cos for rotation of lamdas:
!   dipmag (magnetic dip angle) is in magfield module (magfield.F), 
!   dipmin (minimum dip angle) is in cons module (cons.F) approx. 10 deg.
        if (abs(dipmag(i,lat)) >= dipmin) then
          dip(i) = dipmag(i,lat)
        else
          dip(i) = sign(dipmin,dipmag(i,lat))
        endif
        rot_xx(i) = csdec(i,lat)**2+sin(dip(i))**2*sn2dec(i,lat)
        rot_yy(i) = sn2dec(i,lat)  +sin(dip(i))**2*csdec(i,lat)**2
        rot_xy(i) = cos(dip(i))**2*sndec(i,lat)*csdec(i,lat)
      enddo ! i=lon0,lon1
!
! Ion-neutral momentum transfer collision frequencies [cm^3/s]:
! (Defines "NU" matrix, formerly sub new.F)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          tr(k,i) = 0.5*(ti(k,i)+tn(k,i))            ! ave of tn & ti
          nu_o2p_o2(k,i) = 2.59E-11*sqrt(tr(k,i))*   ! O2+ ~ O2 (resonant)
     |      (1.-0.073*log10(tr(k,i)))**2
          nu_op_o  (k,i) = 3.67e-11*sqrt(tr(k,i))*   ! O+  ~ O  (resonant)
     |      (1.-0.064*log10(tr(k,i)))**2*colfac
          nu_n2p_n2(k,i) = 5.14E-11*sqrt(tr(k,i))*   ! N2+ ~ N2 (resonant)
     |      (1.-0.069*log10(tr(k,i)))**2
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! Major species number densities [1/cm^3]:
! mass density [g/cm^3/mol]: xnmbar half levels
! pressure at half level: p0*expz(k) [dyn/cm^2]
! Boltzmann constant [erg/K]
! neutral temperature [K]
! mass mixing rations[-]: o2, o1, he, n2
! inverse of mean molecular weight [mol/g]:rmassinv_o2, rmassinv_o1, rmassinv_he, rmassinv_n2
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          o2_cm3(k,i) = o2(k,i)*xnmbar(k,i)*rmassinv_o2
          o1_cm3(k,i) = o1(k,i)*xnmbar(k,i)*rmassinv_o1
          he_cm3(k,i) = he(k,i)*xnmbar(k,i)*rmassinv_he
          n2_cm3(k,i) = n2(k,i)*xnmbar(k,i)*rmassinv_n2 ! this part varies from timegcm
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('O2_CM3',' ',' ',o2_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O_CM3' ,' ',' ',o1_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('HE_CM3',' ',' ',he_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2_CM3',' ',' ',n2_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! collision frequency nu_in for each ion [1/s]
!    by multiplying with neutral number density [1/cm^3] and sum over neutrals
! nu_o2p = [o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~he]n(he)+[o2p~n2]n(n2)
! nu_op  = [op ~o2]n(o2)+[op ~o]n(o)+[op ~he]n(he)+[op ~n2]n(n2)
! nu_nop = [nop~o2]n(o2)+[nop~o]n(o)+[nop~he]n(he)+[nop~n2]n(n2)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          nu_o2p(k,i) = nu_o2p_o2(k,i)*o2_cm3(k,i) +
     |                  nu_o2p_o      *o1_cm3(k,i) +
     |                  nu_o2p_he     *he_cm3(k,i) +
     |                  nu_o2p_n2     *n2_cm3(k,i)
          nu_op (k,i) = nu_op_o2      *o2_cm3(k,i) +
     |                  nu_op_o  (k,i)*o1_cm3(k,i) +
     |                  nu_op_he      *he_cm3(k,i) +
     |                  nu_op_n2      *n2_cm3(k,i)
          nu_np (k,i) = nu_np_o2      *o2_cm3(k,i) +
     |                  nu_np_o       *o1_cm3(k,i) +
     |                  nu_np_he      *he_cm3(k,i) +
     |                  nu_np_n2      *n2_cm3(k,i)
          nu_n2p(k,i) = nu_n2p_o2     *o2_cm3(k,i) +
     |                  nu_n2p_o      *o1_cm3(k,i) +
     |                  nu_n2p_he     *he_cm3(k,i) +
     |                  nu_n2p_n2(k,i)*n2_cm3(k,i)
          nu_nop(k,i) = nu_nop_o2     *o2_cm3(k,i) +
     |                  nu_nop_o      *o1_cm3(k,i) +
     |                  nu_nop_he     *he_cm3(k,i) +
     |                  nu_nop_n2     *n2_cm3(k,i)
!
! neutral~electron collision frequency (from Banks & Kockards) nu_en
!
          sqrt_te = sqrt(te(k,i))
          nu_e(k,i) = 
     |      1.82e-10*o2_cm3(k,i)*(sqrt_te+3.60e-2*te(k,i))+
     |      8.90e-11*o1_cm3(k,i)*sqrt_te*(1.+5.70e-4*te(k,i))+
     |      4.60e-10*he_cm3(k,i)*sqrt_te+
     |      2.33e-11*n2_cm3(k,i)*te(k,i)*(1.-1.21e-4*te(k,i))
!
! 6/2/06 btf: Multiply nu_e by 4, as per Richmond:
! The effective electron-neutral collision frequency is increased in 
! an an hoc manner by a factor of 4 in order for the model to produce
! electric fields and currents below 105 km that agree better with
! observations, as recommended by Gagnepain et al. (J. Atmos. Terr. 
! Phys., 39, 1119-1124, 1977).
!
          nu_e(k,i) = nu_e(k,i)*4.

!         write(6,"('lamdas: lat=',i3,' k=',i3,' i=',i3,' te=',
!    |      e12.4,' o2=',e12.4,' o1=',e12.4,' he=',e12.4,' n2=',e12.4,
!    |      ' omega_e=',e12.4,' nu_e=',e12.4)") lat,k,i,te(k,i),
!    |      o2_cm3(k,i),o1_cm3(k,i),he_cm3(k,i),n2_cm3(k,i),
!    |      omega_e(i),nu_e(k,i)

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('NU_O2P',' ',' ',nu_o2p,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('NU_OP' ,' ',' ',nu_op ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('NU_NOP',' ',' ',nu_nop,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('NU_NE' ,' ',' ',nu_e  ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Pedersen and Hall conductivities (siemens/m):
! Qe_fac includes conversion from CGS to SI units 
!  -> e/B [C/T 10^6 m^3/cm^3], see above.
! number densities [1/cm^3]
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! ne = electron density assuming charge equilibrium [1/cm3]:
! (Using ne from fields module will fail, not sure about the reason)
          ne(k,i) = op(k,i)+o2p(k,i)+nplus(k,i)+n2p(k,i)+nop(k,i)
!
! Pedersen conductivity [S/m] (half level):
          sigma_ped(k,i) = qe*1.e10/bmod2(i,lat)*
     |      (ne   (k,i)*omega_e  (i)*nu_e  (k,i)/
     |       (omega_e  (i)**2+nu_e  (k,i)**2)+
     |       op   (k,i)*omega_op (i)*nu_op (k,i)/
     |       (omega_op (i)**2+nu_op (k,i)**2)+
     |       o2p  (k,i)*omega_o2p(i)*nu_o2p(k,i)/
     |       (omega_o2p(i)**2+nu_o2p(k,i)**2)+
     |       nplus(k,i)*omega_np (i)*nu_np(k,i)/
     |       (omega_np (i)**2+nu_np (k,i)**2)+
     |       n2p  (k,i)*omega_n2p(i)*nu_n2p(k,i)/
     |       (omega_n2p(i)**2+nu_n2p(k,i)**2)+
     |       nop  (k,i)*omega_nop(i)*nu_nop(k,i)/
     |       (omega_nop(i)**2+nu_nop(k,i)**2))
!
! Hall conductivity [S/m] (half level):
          sigma_hall(k,i) = qe*1.e10/bmod2(i,lat)*
     |      (ne   (k,i)*omega_e  (i)**2/
     |       (omega_e  (i)**2+nu_e  (k,i)**2)-
     |       op   (k,i)*omega_op (i)**2/
     |       (omega_op (i)**2+nu_op (k,i)**2)-
     |       o2p  (k,i)*omega_o2p(i)**2/
     |       (omega_o2p(i)**2+nu_o2p(k,i)**2)-
     |       nplus(k,i)*omega_np (i)**2/
     |       (omega_np (i)**2+nu_np (k,i)**2)-
     |       n2p  (k,i)*omega_n2p(i)**2/
     |       (omega_n2p(i)**2+nu_n2p(k,i)**2)-
     |       nop  (k,i)*omega_nop(i)**2/
     |       (omega_nop(i)**2+nu_nop(k,i)**2))

          if (sigma_hall(k,i) < 1e-20) sigma_hall(k,i) = 1e-20
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('ELECDEN' ,' ',' ',ne,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ion drag coefficients lamda1,2 [1/s] (full level)
! Pedersen/ Hall conductivity [S/m]: sigma_ped, sigma_hall (half level)
! mass density [g/cm^3/mol]: xnmbar (half level)
! magnetic field [Gauss]: bgauss
! Na Avoagdro number [1/mol]
! conversion: tesla = 1.e-4*bgauss
!             kg/m3 = 1.e-3*g/(1.e-6*cm^3) = 1.e3*g/cm3
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          lamda1tmp(k,i) = sigma_ped (k,i)*(1.e-4*bmod2(i,lat))**2*avo/
     |                     (1.e3*xnmbar(k,i))
          lamda2tmp(k,i) = sigma_hall(k,i)*(1.e-4*bmod2(i,lat))**2*avo/
     |                     (1.e3*xnmbar(k,i))
        enddo ! k=lev0,lev1-1
!
! Lamda1,2 to full (interface) levels:
        do k=lev0,lev1-2
!
! 5/17/04 btf: FP invalid here on SGI dataproc due to negative lamda1tmp
!              at k==lev1 (reference to lamda1tmp at k+1):
!         write(6,"('lamdas: lat=',i3,' i=',i3,' k=',i3,' lamda1tmp=',
!    |      e12.4,' lamda1tmp(k+1,i)=',e12.4)") lat,i,k,
!    |      lamda1tmp(k,i),lamda1tmp(k+1,i)         

          lamda1(k+1,i) = sqrt(lamda1tmp(k,i)*lamda1tmp(k+1,i))
          lamda2(k+1,i) = sqrt(lamda2tmp(k,i)*lamda2tmp(k+1,i))
        enddo ! k=lev0,lev1-1
! Bottom boundary:
        lamda1(lev0,i) = sqrt(lamda1tmp(lev0,i)**3/lamda1tmp(lev0+1,i))
        lamda2(lev0,i) = sqrt(lamda2tmp(lev0,i)**3/lamda2tmp(lev0+1,i))
! Top boundary:
        lamda1(lev1,i)= sqrt(lamda1tmp(lev1-1,i)**3/lamda1tmp(lev1-2,i))
        lamda2(lev1,i)= sqrt(lamda2tmp(lev1-1,i)**3/lamda2tmp(lev1-2,i))
      enddo ! i=lon0,lon1
!
! Save diagnostic ion drag coefficients:
      call mkdiag_LAMDAPED('LAMDA_PED',lamda1(:,lon0:lon1),
     |  lev0,lev1,lon0,lon1,lat)
      call mkdiag_LAMDAHAL('LAMDA_HAL',lamda2(:,lon0:lon1),
     |  lev0,lev1,lon0,lon1,lat)
!
! Rotate lamdas for displacement of geomagnetic and geographic poles:
! (full levels)
!
      do i=lon0,lon1
        do k=lev0,lev1
          lxx(k,i)= lamda1(k,i)*rot_xx(i)
          lyy(k,i)= lamda1(k,i)*rot_yy(i)
          lyx(k,i)= lamda2(k,i)*sin(dip(i))+lamda1(k,i)*rot_xy(i)
          lxy(k,i)= lamda2(k,i)*sin(dip(i))-lamda1(k,i)*rot_xy(i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

      Q1 = 0.0
      Q2 = 0.0
      KncP = 0.0
      KncH = 0.0

      do i=lon0,lon1
        do k=lev0,lev1-1

! mean ion-neutral collision frequency
          nu_i(k,i) =
     |      op   (k,i)*nu_op (k,i)+
     |      o2p  (k,i)*nu_o2p(k,i)+
     |      nplus(k,i)*nu_np (k,i)+
     |      n2p  (k,i)*nu_n2p(k,i)+
     |      nop  (k,i)*nu_nop(k,i)

! ratio between gyro frequency and neutral collision frequency
          ke(k,i) = omega_e(i)/nu_e(k,i)
          ki(k,i) = 1/ne(k,i)*
     |      (op   (k,i)*omega_op (i)/nu_op (k,i)+
     |       o2p  (k,i)*omega_o2p(i)/nu_o2p(k,i)+
     |       nplus(k,i)*omega_np (i)/nu_np (k,i)+
     |       n2p  (k,i)*omega_n2p(i)/nu_n2p(k,i)+
     |       nop  (k,i)*omega_nop(i)/nu_nop(k,i))

! mean ion molecular mass [g/mol]
          mi(k,i) = 1/ne(k,i)*
     |      (op   (k,i)*rmass_o1+
     |       o2p  (k,i)*rmass_o2+
     |       nplus(k,i)*rmass_n4s+
     |       n2p  (k,i)*rmass_n2+
     |       nop  (k,i)*rmass_no)

! Calculate electrojet turbulent heating, added by Jing Liu
! ne: 1/cm^3*1E6 -> 1/m^3
! bmod2: G*1E-4 -> T
! Q1,Q2: J*1E-7 -> erg
! E1,E2: 1E-2 -> V/m 
!    (boltz=1.38E-16 erg/K from cons needs x1e-7 to be in SI unit.
!     mp=1.67e-24 g, needs x1e-3 to be in SI unit.
!     together in sqrt(1e-7/1e-3) gives 1E-2 for E2 to be in SI unit [V/m]
!     and consistent with Etot [V/m]).
          if (abs(rlatm(i,lat))*rtd>50.0 .and. ki(k,i)<1.0) then
            E2 = 1E-2*sqrt(boltz*(te(k,i)+ti(k,i))*(1.0+ki(k,i)**2)/
     |        (mp*mi(k,i)*(1.0-ki(k,i)**2)))*bmod2(i,lat)*1E-4
            E1 = E2*(1.0+1.0/(ke(k,i)*ki(k,i)))

            if (E1 < Etot(k,i)) then
              Q1(k,i) = me*ne(k,i)*1E3*nu_e(k,i)*
     |          (Etot(k,i)/(bmod2(i,lat)*1E-4))**2
              Q2(k,i) = mp*mi(k,i)*1E3*nu_i(k,i)*
     |          ((Etot(k,i)-E1)/(bmod2(i,lat)*1E-4))**2*
     |          (Etot(k,i)/E2-1.0)*ki(k,i)**2/(1.0+ki(k,i)**2)
              KncP(k,i) = (1.0-ki(k,i)**2)/(1.0+ki(k,i)**2)
     |          *(Etot(k,i)-E1)**2/(Etot(k,i)*E1)
              KncH(k,i) = (2.0*ki(k,i))/(1.0+ki(k,i)**2)
     |          *(1.0+1.0/(ke(k,i)*ki(k,i)))
     |          *(Etot(k,i)-E1)**2/(Etot(k,i)*E1)
            endif
          endif
        enddo
      enddo
!
! Output conductivities [S/m] (half levels):
      do i=lon0,lon1
        do k=lev0,lev1
          ! default ped_out/hall_out has ET effects when do_et is true.
          ! use ped_ref_out/hall_ref_out to store the ped without ET.
          ! Note ped_ref_out = ped_out/(1+KncP)
          ! and hall_ref_out = hall_out/(1+KncH)
          ped_out(k,i)  = sigma_ped(k,i)
          ped_ref_out(k,i)  = sigma_ped(k,i)*(1.0+KncP(k,i))
          hall_out(k,i) = sigma_hall(k,i)
          hall_ref_out(k,i) = sigma_hall(k,i)*(1.0+KncH(k,i))
        enddo
      enddo
!
! Save diagnostic conductivities (midpoints):
      call mkdiag_SIGMAPED('SIGMA_PED',ped_out,lev0,lev1,lon0,lon1,lat)
      call mkdiag_SIGMAHAL('SIGMA_HAL',hall_out,lev0,lev1,lon0,lon1,lat)

      call addfld('Q1','Electron Ohmic Heating',' ',
     |  Q1(lev0:lev1-1,lon0:lon1),
     |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      call addfld('Q2','Electrojet Turbulent Heating',' ',
     |  Q2(lev0:lev1-1,lon0:lon1),
     |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      call addfld('KncP','Pedersen nonlinear current coefficient',' ',
     |  KncP(lev0:lev1-1,lon0:lon1),
     |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      call addfld('SIGMA_PED_REF','Pedersen conductivity with ET',
     |  'S/m',ped_ref_out(lev0:lev1-1,lon0:lon1),
     |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      call addfld('KncH','Hall nonlinear current coefficient',' ',
     |  KncH(lev0:lev1-1,lon0:lon1),
     |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      call addfld('SIGMA_HALL_REF','Hall conductivity with ET',
     |  'S/m',hall_ref_out(lev0:lev1-1,lon0:lon1),
     |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      end subroutine lamdas
