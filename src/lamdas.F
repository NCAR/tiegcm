!
      subroutine lamdas(tn,barm,o2,o1,n2,ti,te,o2p,op,nop, lxx,lyy,lxy,
     |  lyx,lamda1,ped_out,hall_out,lev0,lev1,lon0,lon1,lat)
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
      use params_module,only: nlonp4
      use magfield_module,only: bmod2,sn2dec,sncsdc,csdec,sndec
      use cons_module,only: dipmin,rmassinv_o2,rmassinv_o1,rmassinv_n2,
     |  p0,expz,rmass_o2,rmass_o1,boltz,avo
      use magfield_module,only: dipmag
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
     |  barm, ! mean molecular mass
     |  o2,   ! molecular oxygen (mmr)
     |  o1,   ! atomic oxygen    (mmr)
     |  n2,   ! molecular nitrogen (mmr)
     |  ti,   ! ion temperature (deg K)
     |  te,   ! electron temperature (deg K)
     |  o2p,  ! O2+ number density (1/cm3)
     |  op,   ! O+  number density (1/cm3)
     |  nop   ! NO+ number density (1/cm3)
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  lxx,     ! lamda XX term [1/s] 
     |  lyy,     ! lamda YY term [1/s] 
     |  lxy,     ! lamda XY term [1/s] 
     |  lyx,     ! lamda YX term [1/s] 
     |  lamda1,  ! lamda_1 (sigma_p * B^2) / rho [1/s] magnetic eastward direction
     |  ped_out, ! pedersen conductivity (siemens/m) (will be input to dynamo)
     |  hall_out ! hall conductivity (siemens/m) (will be input to dynamo)
!
! Local:
      integer :: k,i,l,lonbeg,lonend,i0,i1
      real :: sqrt_te         ! sqrt(te)
      real,parameter :: qe=1.602e-19,  ! electron charge (coulomb)
     |    qeomeo10 =  1.7588028E7,     ! qe/m_e/10 [C/g]
     |    qeoNao10 =  9.6489E3         ! qe/N_a/10 [C/mol]
      real,parameter :: rmass_nop = 30.! [g/mol]
      real,parameter :: rmassinv_nop=1./rmass_nop ! inverted rmass
      real,dimension(lev0:lev1,lon0:lon1) :: ped_plt,hall_plt
!
! Local (lon):
      real,dimension(lon0:lon1) :: 
     |  bgauss, qe_fac, sindip, cosdip, cos2dip, sin2dip, cos2dec,
     |  omega_o2p    ,omega_op    ,omega_nop,
     |  omega_o2p_inv,omega_op_inv,omega_nop_inv,
     |  omega_e      ,omega_e_inv
!
! Local (lev,lon):
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  xnmbar,     ! for conversion to volume density
     |  xn2,        ! N2 (1-o2-o-he mmr)
     |  tnti,       ! average of tn and ti
     |  o2_cm3, o1_cm3, n2_cm3, ! major species number densities (cm3)
     |  sigma_ped,  ! pedersen conductivity (siemens/m)
     |  sigma_hall, ! hall conductivity (siemens/m)
     |  ne,         ! electron density (assume o2p+op+nop)
     |  lamda2,     ! sighal*b**2/rho
     |  lamda1tmp,  ! temporary lamda1
     |  lamda2tmp,  ! temporary lamda2
     |  lxxnorot,   ! XX before rotation
     |  lyynorot,   ! YY before rotation
     |  lxynorot,   ! XY before rotation
     |  lyxnorot,   ! YX before rotation
!
! Ion-neutral momentum transfer collision frequencies:
!
     |  rnu_o2p_o2, ! O2+ ~ O2 collision freq (resonant, temperature dependent)
     |  rnu_op_o2 , ! O+  ~ O2 collision freq (non-resonant)
     |  rnu_nop_o2, ! NO+ ~ O2 collision freq (non-resonant)
!
     |  rnu_o2p_o,  ! O2+ ~ O  collision freq (non-resonant)
     |  rnu_op_o ,  ! O+  ~ O  collision freq (resonant, temperature dependent)
     |  rnu_nop_o,  ! NO+ ~ O  collision freq (non-resonant)
!
     |  rnu_o2p_n2, ! O2+ ~ N2 collision freq (non-resonant)
     |  rnu_op_n2 , ! O+  ~ N2 collision freq (non-resonant)
     |  rnu_nop_n2, ! NO+ ~ N2 collision freq (non-resonant)
!
     |  rnu_o2p,    ! [[o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~n2]n(n2)]/w(o2p)
     |  rnu_op,     ! [[op ~o2]n(o2)+[op ~o]n(o)+[op ~n2]n(n2)]/w(op )
     |  rnu_nop,    ! [[nop~o2]n(o2)+[nop~o]n(o)+[nop~n2]n(n2)]/w(nop)
! 
     |  rnu_ne      ! electron~neutral
!
! Save input args to secondary history:
!     real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
!    |  tn,   ! neutral temperature (deg K)
!    |  barm, ! mean molecular mass
!    |  o2,   ! molecular oxygen (mmr)
!    |  o1,   ! atomic oxygen    (mmr)
!    |  n2,   ! molecular nitrogen (mmr)
!    |  ti,   ! ion temperature (deg K)
!    |  te,   ! electron temperature (deg K)
!    |  o2p,  ! O2+ number density (1/cm3)
!    |  op,   ! O+  number density (1/cm3)
!    |  nop   ! NO+ number density (1/cm3)

      i0 = lon0 ; i1 = lon1

!     call addfld('TN_lam',' ',' ',tn(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('BARM_lam',' ',' ',barm(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O2_lam',' ',' ',o2(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O1_lam',' ',' ',o1(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2_lam',' ',' ',n2(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TI_lam',' ',' ',ti(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TE_lam',' ',' ',te(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Save input ions to secondary history:
!
!     call addfld('O2P_LAM',' ',' ',o2p(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP_LAM' ,' ',' ',op (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('NOP_LAM',' ',' ',nop(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Set local needs:
      do i=lon0,lon1
        bgauss(i) = bmod2(i,lat) ! magnetic field strength [Gauss]
!
!  e/B  [C/T 10^6 cm^3/m^3]
! 1.e10 for SI units = 1.e6 (cm3->m3) * 1.e4 (gauss->tesla)
!
        qe_fac(i) = qe*1.e10/bgauss(i) 
! 
! gyrofrequencies: omega_i = eB/m_i  [1/s]
!                  omega_e = eB/m_e  [1/s]
! with qeoNao10 = e/Na [C/mol T/kg g/Gauss]
!      qeomeo10 = e/m_e [C/g T/kg g/Gauss]
! 1/10 in qeoNao10 and qeomeo10 for conversion from Gauss/g to T/kg
!
        omega_op (i)    = qeoNao10*bgauss(i)*rmassinv_o1
        omega_o2p(i)    = qeoNao10*bgauss(i)*rmassinv_o2
        omega_nop(i)    = qeoNao10*bgauss(i)*rmassinv_nop
        omega_op_inv (i)= 1./omega_op(i)
        omega_o2p_inv(i)= 1./omega_o2p(i)
        omega_nop_inv(i)= 1./omega_nop(i)
        omega_e(i)      = qeomeo10*bgauss(i)
        omega_e_inv(i)  = 1./omega_e(i)
!
! Sin and cos for rotation of lamdas:
!   dipmag (magnetic dip angle) is in magfield module (magfield.F), 
!   dipmin (minimum dip angle) is in cons module (cons.F) approx. 10 deg.
        if (abs(dipmag(i,lat)) >= dipmin) then
          sindip(i) = sin(dipmag(i,lat))
          cosdip(i) = cos(dipmag(i,lat))
        else
          if (dipmag(i,lat) >= 0.) then
            sindip(i) = sin(dipmin)
            cosdip(i) = cos(dipmin)
          else
            sindip(i) = sin(-dipmin)
            cosdip(i) = cos(-dipmin)
          endif
        endif
        cos2dip(i) = cosdip(i)**2
        sin2dip(i) = sindip(i)**2
        cos2dec(i) = csdec(i,lat)**2
      enddo ! i=lon0,lon1
!
! Ion-neutral momentum transfer collision frequencies [cm^3/s]:
! (Defines "NU" matrix, formerly sub new.F)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          tnti(k,i) = 0.5*(ti(k,i)+tn(k,i))           ! ave of tn & ti
!
! O2 collision frequencies:
          rnu_o2p_o2(k,i) = 2.59E-11*sqrt(tnti(k,i))* ! O2+ ~ O2 (resonant)
     |      (1.-0.073*alog10(tnti(k,i)))**2
          rnu_op_o2 (k,i) = 6.64E-10                  ! O+  ~ O2
          rnu_nop_o2(k,i) = 4.27E-10                  ! NO+ ~ O2
!
! O collision frequencies:
          rnu_o2p_o(k,i) = 2.31E-10                   ! O2+ ~ O
          rnu_op_o (k,i) = 3.67e-11*sqrt(tnti(k,i))*  ! O+  ~ O (resonant)
     |      (1.-0.064*alog10(tnti(k,i)))**2*colfac    
          rnu_nop_o(k,i) = 2.44E-10                   ! NO+ ~ O
!
! N2 collision frequencies:
          rnu_o2p_n2(k,i) = 4.13E-10                  ! O2+ ~ N2
          rnu_op_n2 (k,i) = 6.82E-10                  ! O+  ~ N2
          rnu_nop_n2(k,i) = 4.34E-10                  ! NO+ ~ N2
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! Major species number densities [1/cm^3]:
! mass density [g/cm^3/mol]: xnmbar half levels
! pressure at half level: p0*expz(k) [dyn/cm^2]
! mean molecular weight at full level: barm [g/mol]
! Boltzmann constant [erg/K]
! neutral temperature [K]
! mass mixing rations[-]: o2, o1, xn2
! inverse of mean molecular weight [mol/g]:rmassinv_o2, rmassinv_o1, rmassinv_n2
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xnmbar(k,i) = p0*expz(k)*.5*(barm(k,i)+barm(k+1,i))/
     |      (boltz*tn(k,i))
          o2_cm3(k,i) = o2(k,i)*xnmbar(k,i)*rmassinv_o2
          o1_cm3(k,i) = o1(k,i)*xnmbar(k,i)*rmassinv_o1
          xn2(k,i) = n2(k,i)
          xn2(k,i) = xn2(k,i)*rmassinv_n2            ! this part varies from timegcm
          if (xn2(k,i) < 1.e-20) xn2(k,i) = 1.e-20   ! done as in tiegcm
          n2_cm3(k,i) = xn2(k,i)*xnmbar(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('O2_CM3',' ',' ',o2_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O_CM3' ,' ',' ',o1_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2_CM3',' ',' ',n2_cm3,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! collision frequency nu_in for each ion [1/s]
!    by multiplying with neutral number density [1/cm^3] and sum over neutrals
! nu_in is divided by gyrofrequency omega_i
! nu_in/omega_i [-]:
! rnu_o2p = [[o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~n2]n(n2)]/w(o2p)
! rnu_op  = [[op ~o2]n(o2)+[op ~o]n(o)+[op ~n2]n(n2)]/w(op )
! rnu_nop = [[nop~o2]n(o2)+[nop~o]n(o)+[nop~n2]n(n2)]/w(nop)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          rnu_o2p(k,i) = (rnu_o2p_o2(k,i)*o2_cm3(k,i) +
     |                    rnu_o2p_o (k,i)*o1_cm3(k,i) +
     |                    rnu_o2p_n2(k,i)*n2_cm3(k,i))*omega_o2p_inv(i)
          rnu_op (k,i) = (rnu_op_o2 (k,i)*o2_cm3(k,i) +
     |                    rnu_op_o  (k,i)*o1_cm3(k,i) +
     |                    rnu_op_n2 (k,i)*n2_cm3(k,i))*omega_op_inv(i)
          rnu_nop(k,i) = (rnu_nop_o2(k,i)*o2_cm3(k,i) +
     |                    rnu_nop_o (k,i)*o1_cm3(k,i) +
     |                    rnu_nop_n2(k,i)*n2_cm3(k,i))*omega_nop_inv(i)
!
! neutral~electron collision frequency (from Banks & Kockards) nu_en
! divided by gyrofrequency omega_2:
! nu_en/omega_e [-]
!
          sqrt_te = sqrt(te(k,i))
          rnu_ne(k,i) = 
     |      (2.33e-11*n2_cm3(k,i)*te(k,i)*(1.-1.21e-4*te(k,i))+
     |       1.82e-10*o2_cm3(k,i)*sqrt_te*(1.+3.60e-2*sqrt_te)+
     |       8.90e-11*o1_cm3(k,i)*sqrt_te*(1.+5.70e-4*te(k,i)))*
     |       omega_e_inv(i)
!
! 6/2/06 btf: Multiply rnu_ne by 4, as per Richmond:
! The effective electron-neutral collision frequency is increased in 
! an an hoc manner by a factor of 4 in order for the model to produce
! electric fields and currents below 105 km that agree better with
! observations, as recommended by Gagnepain et al. (J. Atmos. Terr. 
! Phys., 39, 1119-1124, 1977).
!
          rnu_ne(k,i) = rnu_ne(k,i)*4.

!         write(6,"('lamdas: lat=',i3,' k=',i3,' i=',i3,' te=',
!    |      e12.4,' n2=',e12.4,' o2=',e12.4,' o1=',e12.4,' omega_e_inv='
!    |      ,e12.4,' rnu_ne=',e12.4)") lat,k,i,te(k,i),n2_cm3(k,i),
!    |      o2_cm3(k,i),o1_cm3(k,i),omega_e_inv(i),rnu_ne(k,i)

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('RNU_O2P',' ',' ',rnu_o2p,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RNU_OP' ,' ',' ',rnu_op ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RNU_NOP',' ',' ',rnu_nop,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RNU_NE' ,' ',' ',rnu_ne ,
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
          ne(k,i) = op(k,i)+o2p(k,i)+nop(k,i) 
!
! Pedersen conductivity [S/m] (half level):
          sigma_ped(k,i) = qe_fac(i)*
     |      ((op (k,i)*rnu_op (k,i)/(1.+rnu_op (k,i)**2))+
     |       (o2p(k,i)*rnu_o2p(k,i)/(1.+rnu_o2p(k,i)**2))+
     |       (nop(k,i)*rnu_nop(k,i)/(1.+rnu_nop(k,i)**2))+
     |       (ne (k,i)*rnu_ne (k,i)/(1.+rnu_ne (k,i)**2)))
!
! Hall conductivity [S/m] (half level):
          sigma_hall(k,i) = qe_fac(i)*
     |      (ne (k,i)/(1.+rnu_ne (k,i)**2)-
     |       op (k,i)/(1.+rnu_op (k,i)**2)-
     |       o2p(k,i)/(1.+rnu_o2p(k,i)**2)-
     |       nop(k,i)/(1.+rnu_nop(k,i)**2))

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('ELECDEN' ,' ',' ',ne,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SIGPEDin',' ',' ',sigma_ped(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('SIGHALin',' ',' ',sigma_hall(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
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
          lamda1tmp(k,i) = (sigma_ped(k,i)*(1.e-4*bgauss(i))**2)*avo/
     |                     (1.e3*xnmbar(k,i))
          lamda2tmp(k,i) = (sigma_hall(k,i)*(1.e-4*bgauss(i))**2)*avo/
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
!
! Non-rotated lamdas:
        do k=lev0,lev1
          lxxnorot(k,i) = lamda1(k,i)
          lyynorot(k,i) = lamda1(k,i)*sin2dip(i)
          lxynorot(k,i) = lamda2(k,i)*sindip(i)
          lyxnorot(k,i) = lxynorot(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! Save diagnostic ion drag coefficients:
      call mkdiag_LAMDAPED('LAMDA_PED',lamda1(:,lon0:lon1),
     |  lev0,lev1,lon0,lon1,lat)
      call mkdiag_LAMDAHAL('LAMDA_HAL',lamda2(:,lon0:lon1),
     |  lev0,lev1,lon0,lon1,lat)

!     call addfld('LAMDA1',' ',' ',lamda1,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('LAMDA2',' ',' ',lamda2,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('LXXNOROT',' ',' ',lxxnorot,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('LYYNOROT',' ',' ',lyynorot,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('LXYNOROT',' ',' ',lxynorot,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('LYXNOROT',' ',' ',lyxnorot,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Rotate lamdas for displacement of geomagnetic and geographic poles:
! (full levels)
!
      do i=lon0,lon1
        do k=lev0,lev1
          lxx(k,i)= lxxnorot(k,i)*cos2dec(i)+lyynorot(k,i)*sn2dec(i,lat)
          lyy(k,i)= lyynorot(k,i)*cos2dec(i)+lxxnorot(k,i)*sn2dec(i,lat)
          lyx(k,i)= lxynorot(k,i)-(lyynorot(k,i)-lxxnorot(k,i))*
     |      sndec(i,lat)*csdec(i,lat)
          lxy(k,i)= lxynorot(k,i)+(lyynorot(k,i)-lxxnorot(k,i))*
     |      sndec(i,lat)*csdec(i,lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
 
!     call addfld('LXX','LXX','Hz',lxx(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('LYY','LYY','Hz',lyy(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('LXY','LXY','Hz',lxy(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('LYX','LYX','Hz',lyx(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Output conductivities [S/m] (half levels):
      do i=lonbeg,lonend
        do k=lev0,lev1
          ped_out(k,i)  = sigma_ped(k,i)
          hall_out(k,i) = sigma_hall(k,i)
        enddo
      enddo
!
! Save diagnostic conductivities (midpoints):
      call mkdiag_SIGMAPED('SIGMA_PED',ped_out,lev0,lev1,lon0,lon1,lat)
      call mkdiag_SIGMAHAL('SIGMA_HAL',hall_out,lev0,lev1,lon0,lon1,lat)

      end subroutine lamdas
