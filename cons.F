!
      module cons_module
      use params_module,only: dlat,dz,nlon,nlonp1,nlonp4,nlat,nlatp1,
     |  dlev,nlev,nmlat,nmlon,nmlonp1,nmagphrlat,nmagphrlon,
     |  magphrlat1,magphrlat2,magphrlon1,plev1,zst
      implicit none
!
! Define model constants. 
! Parameter constants are cons_module module data and are accessed 
!   in subroutines via use-association. 
! Derived constants are cons_module module data, and are calculated
!   in sub init_cons (contained in cons_module). 
!   Sub init_cons is called by sub init (init_mod.F).
! Parameters referenced here are in params.h (s.a., dims.h)
!
! Real parameter constants:
!
      real,parameter ::
     |  dzp  = dz,           ! alias for dz (also dlev)
     |  re   = 6.37122e8,    ! earth radius (cm)                  C(51)
     |  re_inv = 1./re,      ! inverse of earth radius            C(52)
     |  avo  = 6.023e23,     ! avogadro number                    C(85)
     |  boltz = 1.38E-16,    ! boltzman's constant                C(84)
     |  p0   = 5.0e-4,       ! standard pressure                  C(81)
     |  gask = 8.314e7,      ! gas constant                       C(57)
!    |  dipmin = .005,       ! minimum mag dip angle (tiegcm) (tgcm15=.005)
     |  dipmin = 0.17,       ! minimum mag dip angle (tiegcm)
     |  brn2d = 0.6,         !
     |  tbound = 181.0,      ! background tn at lower boundary (tiegcm)
     |  atm_amu = 28.9,      ! mean mass of surface atmosphere    C(24)
     |  shapiro = 3.0E-2,    ! shapiro smoother constant          C(26)
     |  dtsmooth = 0.95,     ! time smoothing constant            C(30)
     |  dtsmooth_div2 = 0.5*(1.-dtsmooth), !                      C(31)
     |  tgrad = 6.,          ! TN gradient from old /RUNMDS/ (tiegcm)
     |  nob(nlat) = 4.e6,    ! N(NO) LBC from old /RUNMDS/ 
     |  avto = 4.0e-12,      ! from old /RUNMDS/ (not used)
     |  hor(nlat) = .25,     ! horizontal variation of eddy diffusion and
                             !  eddy thrmal conductivity. If unity, value
                             !  of KE at bottom is 5.e-6 (old /RUNMDS/)
     |  prndtl = 1.,         ! prandtl number
     |  evergs = 1.602e-12,  ! 1 eV = 1.602e-12 ergs
     |  tsurplus=5.11*evergs ! surplus heat per event (ergs)    C(45)
!
      integer :: nlonper=nlonp4 ! nlon + periodic points (alias for nlonp4)
!
! Many expressions require x/rmass, but its more efficient on some
! platforms to multiply rather than divide, so set rmassinv = 1./rmass 
! here, and use x*rmassinv in the code.
!
      real,parameter :: 
     |  rmass(3) = (/32.,16.,28./), ! o2,o,n2
     |  rmass_o2 = 32., rmass_o1  = 16., rmass_n2  = 28., 
     |  rmass_o3 = 48., rmass_n4s = 14., rmass_n2d = 14., 
     |  rmass_no = 30., rmass_op  = 16., rmass_co2 = 44.
      real,parameter ::
     |  rmassinv_o2  = 1./rmass_o2,
     |  rmassinv_o1  = 1./rmass_o1,
     |  rmassinv_n2  = 1./rmass_n2,
     |  rmassinv_o3  = 1./rmass_o3,
     |  rmassinv_n4s = 1./rmass_n4s,
     |  rmassinv_n2d = 1./rmass_n2d,
     |  rmassinv_no  = 1./rmass_no,
     |  rmassinv_op  = 1./rmass_op
!
! 2/00: these were in modsrc.snoe (tgcm13mt), but were unused.
! Low-energy protons:
!     real,parameter ::
!       alfalp = 10.,
!       efluxlp = 1.e-20
!
! Model derived constants (see sub init_cons in this module):
!
      real ::
     |  pi,             ! set with 4*atan(1)    C(110)
     |  rtd,            ! radians-to-degrees (180./pi)
     |  dtr,            ! degrees-to-radians (pi/180.)
     |  dphi,           ! delta lat (pi/nlat)   C(2)
     |  dphi_2div3,     ! 2./(3.*dphi)          C(12)
     |  dphi_1div12,    ! 1./(12.*dphi)         C(13)
     |  dlamda,         ! delta lon (2pi/nlon)  C(1)
     |  dlamda_2div3,   ! 2./(3.*dlamda)        C(10)
     |  dlamda_1div12,  ! 1./(12.*dlamda)       C(11)
     |  dt,             ! time step (secs)      C(4)
     |  dtx2,           ! 2*dt                  C(6)
     |  dtx2inv,        ! 1./(2*dt)             C(7)
     |  freq_3m3,       ! frequency of 2-day wave (rad/sec)       C(21)
     |  freq_semidi,    ! frequency of semidiurnal tide (rad/sec) C(23)
     |  freq_ann,       ! frequency of annual tide                C(25)
     |  expz(nlev+1),   ! exp(-z) at midpoints
     |  expzmid,        ! exp(-.5*dz)                             C(86)
     |  expzmid_inv,    ! 1./expzmid                              C(87)
     |  rmassinv(3),    ! inverse of rmass
     |  t0(nlev+1),     ! set by sub lowbound (bndry_mod.F)
     |  racs(nlat),     ! 1./(re*cs(lat))
     |  cs(-1:nlat+2),  ! cos(phi)
     |  sn(nlat),       ! sin(phi)
     |  tn(nlat),       ! tan(phi)
     |  cor(nlat),
     |  grav,           ! accel due to gravity (dependent on lower boundary)
     |  dzgrav,         ! grav/gask C(65)
     |  difk(nlev+1),   ! background eddy diffusion
     |  dift(nlev+1),   ! background thermal conductivity
     |  xmue(nlev+1)    ! eddy viscosity (?)
!
! Constants for dynamo and electric field calculations:
      real,parameter :: re_dyn = 6.378165e8  ! earth radius (cm)
      real,parameter :: h00=9.7e6, r00=re+h00	! use mean earth radius
      real,parameter :: h0 =9.0e6, r0 =re+h0	! use mean earth radius
      real,parameter :: hs=1.3e7
      real :: 
     |  dlatg, dlong, dlatm, dlonm,dmagphrlon,
     |  ylatm(nmlat),    ! magnetic latitudes (radians)
     |  ylonm(nmlonp1),  ! magnetic longitudes (radians)
     |  ylatmagphr(nmagphrlat),  ! magnetosphere latitudes (radians)
     |  ylonmagphr(nmagphrlon),  ! magnetosphere longitudes (radians)
     |  rcos0s(nmlat),   ! cos(theta0)/cos(thetas)
     |  dt0dts(nmlat),   ! d(theta0)/d(thetas)
     |  dt1dts(nmlat),   ! dt0dts/abs(sinim) (non-zero at equator)
     |  table(91,2)
!
! Geographic grid in radians:
      real ::
     |  ylatg(0:nlatp1), ! geographic latitudes (radians)
     |  ylong(nlonp1)    ! geographic longitudes (radians)
!
! Critical colatitude limits for use of Heelis potential in dynamo:
      real,parameter :: 
     |  crit(2) = (/0.261799387, 0.523598775/)
!
! For filtering longitudinal waves (see filter.F):  
!     integer,parameter :: kut(nlat) =
!    |  (/1,2,3,5,6,7,9,10,11,13,14,15,17,17,17,17,17,17,17,17,17,17,17,
!    |    17,15,14,13,11,10,9,7,6,5,3,2,1/)
!
! kut for wave numbers to filter:
      integer :: kut(nlat)
!
! kut for tiegcm at dlat 5.0 degrees:
      integer,parameter :: kut_5(36) =
     |  (/1,2,3,5,6,7,9,10,11,13,14,15,17,17,17,17,17,17,17,17,17,17,17,
     |    17,15,14,13,11,10,9,7,6,5,3,2,1/)
!
! If check_exp is set true, certain routines will use expo() (util.F)
! instead of exp(). expo checks for out of range arguments to the
! exponential, substituting large or small results if the argument
! is out of range. This avoids NaNS fpe's, but degrades performance.
! It will also produce slightly different results.
!
!     logical,parameter :: check_exp = .true.
      logical,parameter :: check_exp = .false.
!
! Special pi for mag field calculations. If pi=4.*atan(1.) and code is
! linked with -lmass lib, then the last 2 digits (16th and 17th) of pi
! are different (56 instead of 12), resulting in theta0(j=49)==0., which 
! is wrong (should be .1110e-15).
!
      real,parameter ::
     |  pi_dyn=3.14159265358979312
!
      contains
!-----------------------------------------------------------------------
      subroutine init_cons
      use input_module,only: step
!
! Set derived constants (this is called from sub init in init_module)
!
! Local:
      real :: z,expdz,phi
      real :: omega = 7.292E-5
      integer :: k,i,j,nfsech,js
!
      pi = 4.*atan(1.)                ! C(110)

!     print *,' init_cons: pi=',pi

      rtd = 180./pi                   ! radians to degrees
      dtr = pi/180.                   ! degrees to radians
      dphi = pi/float(nlat)           ! C(2)
      dphi_2div3 = 2./(3.*dphi)       ! C(12)
      dphi_1div12 = 1./(12.*dphi)     ! C(13)
      dlamda = 2.*pi/float(nlon)      ! C(1)
      dlamda_2div3 = 2./(3.*dlamda)   ! C(10)
      dlamda_1div12 = 1./(12.*dlamda) ! C(11)
!
! expz(kmax) is exp(-zp) at midpoints:
! expz (will replace EXPS) (expz(nlev+1) not used).
!
      expz(:) = 0. ! init
!
! bottom midpoint z = plev1 + 1/2 deltaz (deltaz==dz==0.5 or 0.25)
! (plev1 and dz are in params.h)
      z = plev1+.5*dlev     
      expz(1) = exp(-z)
      expdz = exp(-dlev)
      difk(1) = 5.0e-6
      dift(1) = 5.0e-6/prndtl
      xmue(1) = 5.0e-6
      do k=2,nlev
        expz(k) = expz(k-1)*expdz
        difk(k) = difk(k-1)*expdz
        xmue(k) = difk(k)
        dift(k) = dift(k-1)*expdz
      enddo
      difk(nlev+1) = difk(nlev)*expdz
      dift(nlev+1) = dift(nlev)*expdz
      xmue(nlev+1) = difk(nlev)
      expzmid = exp(-.5*dlev)
      expzmid_inv = 1./expzmid
      do i=1,3
        rmassinv(i) = 1./rmass(i)
      enddo
      js=-(nlat/2)
      do j=1,nlat
        phi=(j+js-.5)*dphi
        cs(j)=cos(phi)
        sn(j)=sin(phi)
        tn(j)=tan(phi)
        cor(j)=2.*omega*sn(j)
        racs(j) = 1./(re*cs(j))
      enddo
!
! cs at 0, -1, nlat+1, and nlat+2 replace the old cssp and csnp:
      cs(-1) = -cs(2)
      cs(0)  = -cs(1)
      cs(nlat+1) = -cs(nlat)
      cs(nlat+2) = -cs(nlat-1)
      dt = float(step)                   ! was C(4)
      dtx2 = 2.*dt                       ! was C(6)
      dtx2inv = 1./dtx2                  ! was C(7)
      freq_3m3 = 2.*pi/(49.7789*60.*60.) ! was C(21)
      freq_semidi = 4.*pi/(24.*60.*60.)  ! was C(23)
      freq_ann = freq_semidi/(2.*365.25) ! was C(25)
!
! Set gravity according to lower boundary:
      if (plev1==-17.) then     ! time-gcm
        grav = 945.
      elseif (plev1==-7.) then  ! tiegcm
        grav = 870.
      else
        write(6,"(/,'>>> WARNING: do not know how to assign gravity',
     |    ' constant with plev1=',f10.2)") plev1
        grav = 945.
        write(6,"('  Default to grav=',f10.2,/)") grav
      endif
      dzgrav = grav/gask   ! C(65)
!
! Set kut for wave filtering according to dlat (2.5 or 5.0):
      call set_wave_filter(36,kut_5,nlat,kut)
      write(6,"('init_cons: dlat=',f6.2,' nlat=',i3,' kut=',/,(12i4))")
     |  dlat,nlat,kut
!
! Set dynamo constants:
      call consdyn
!
! Report to stdout:
!     write(6,"(/,'Model name    = ',a)") tgcm_name
!     write(6,"(  'Model version = ',a)") tgcm_version
      write(6,"(/,'Set constants:')")
      write(6,"('  nlat=',i3,' nlon=',i3,' nlev=',i3)") nlat,nlon,nlev
      write(6,"('  zst=',f8.2,' plev1=',f8.2,' dlev=',f5.2)") 
     |  zst,plev1,dlev
      write(6,"('  dt = ',f8.2,' secs')") dt
      write(6,"('  grav = ',f10.2)") grav
      write(6,"('  freq_3m3 = ',e10.4,'  freq_semidi=',e10.4,
     |  '  freq_ann=',e12.4)") freq_3m3,freq_semidi,freq_ann
      write(6,"('  dipmin = ',f8.3)") dipmin

!     write(6,"('kut     = ',36i3)") kut
!     write(6,"('2*kut+3 = ',36i3)") 2*kut(:)+3
!
      end subroutine init_cons
!-----------------------------------------------------------------------
      subroutine set_wave_filter(nlat5,kut5,nlat,kut)
!
! Args:
      integer,intent(in) :: nlat5,nlat
      integer,intent(in) :: kut5(nlat5)
      integer,intent(out) :: kut(nlat)
!
! Local:
      integer :: j
!
      if (nlat==nlat5) then ! nlat==nlat5==36 (5x5 degree res)
        do j=1,nlat
          kut(j) = kut5(j)
        enddo
      elseif (nlat==nlat5*2) then ! nlat==72 (2.5x2.5 degree res)
        do j=1,nlat5-1
          kut(j*2-1) = kut5(j) ! 1,3,5,...,65,67,69
          kut(j*2)   = kut5(j) ! 2,4,6,...,66,68,70
        enddo
        kut(nlat) = kut5(nlat5)
        kut(nlat-1) = kut5(nlat5)
      else
        write(6,"('set_wave_filter: nlat=',i3,' dlat=',f8.3,
     |    ' not supported.')") nlat,dlat
        stop 'dlat'
      endif
!     write(6,"('set_wave_filter: nlat=',i3,' kut=',/,(12i4))") nlat,kut
      end subroutine set_wave_filter
!-----------------------------------------------------------------------
      subroutine consdyn
      use input_module,only: dynamo
!
! Set derived constants used in dynamo.
!
! Local:
      integer :: j,i,n
      real,parameter :: e=1.e-6, r1=1.06e7, alfa=1.668
      real :: 
     |  tanth0(nmlat),
     |  tanths(nmlat),
     |  theta0(nmlat),
     |  hamh0(nmlat)
      real :: dtheta,table2(91,3:5),tanths2
      real :: fac,rmin,rmax,rmag
!
! Set grid deltas:
      dlatg = pi/float(nlat)
      dlong = 2.*pi/float(nlon)
      dlatm = pi_dyn/float(nmlat-1) ! note use of pi_dyn
      dlonm = 2.*pi_dyn/float(nmlon) 
      dmagphrlon = 360./float(nmagphrlon)
!
! Set geographic latitude array ylatg:
      do j=1,nlat
        ylatg(j) = -0.5*(pi-dlatg)+float(j-1)*dlatg
      enddo ! j=1,nlat
      ylatg(0) = -pi/2.+e
      ylatg(nlatp1) = pi/2.-e
!
! Set geographic longitude array ylong:
      do i=1,nlonp1
        ylong(i) = -pi+float(i-1)*dlong
      enddo ! i=1,nmlonp1
!
! Set magnetic latitudes ylatm and magnetic longitudes ylonm:
!
! ylatm is equally spaced in theta0, but holds corresponding value
!   of thetas.
      do j=1,nmlat
        theta0(j) = -pi_dyn/2.+float(j-1)*dlatm ! note use of pi_dyn
      enddo ! j=1,nmlat
      do j=2,nmlat-1
        tanth0(j) = abs(tan(theta0(j)))
        hamh0(j) = r1*tanth0(j)+r0*tanth0(j)**(2.+2.*alfa)/
     |    (1.+tanth0(j)**2)**alfa
        tanths(j) = sqrt(hamh0(j)/r0)
        ylatm(j) = sign(atan(tanths(j)),theta0(j))
        rcos0s(j) = sqrt((1.+tanths(j)**2)/(1.+tanth0(j)**2))
!
! If dynamo==0 -> dynamo is not calculated
! If dynamo==1 -> old dynamo is called
! If dynamo==2 -> new dynamo is called
!
        if (dynamo /= 2) then ! old dynamo (or no dynamo)
          dt0dts(j) = (2.*r0*tanths(j)*(1.+tanths(j)**2))/
     |      (r1*(1.+tanth0(j)**2)+2.*r0*tanth0(j)**(2.*alfa+1.)*
     |      (1.+alfa+tanth0(j)**2)/(1.+tanth0(j)**2)**alfa)
        else ! dynamo==2 -> new dynamo
          tanths2  = tanths(j)**2
          dt1dts(j) = 
     |      (r0*sqrt(1.+4.*tanths2)*(1.+tanths2))/
     |      (r1*(1.+tanth0(j)**2)+2.*r0*tanth0(j)**(2.*alfa+1.)*
     |      (1.+alfa+tanth0(j)**2)/(1.+tanth0(j)**2)**alfa)
          dt0dts(j) = dt1dts(j)*2.*tanths(j)/sqrt(1.+4.*tanths2)
        endif
      enddo ! j=2,nmlat-1
!
! Magnetic poles:
      ylatm(1) = theta0(1)
      ylatm(nmlat) = theta0(nmlat)
      rcos0s(1) = 1.
      rcos0s(nmlat) = 1.
      dt0dts(1) = 1.
      dt0dts(nmlat) = 1.
!
! Magnetic longitudes:
      do i=1,nmlonp1
        ylonm(i) = -pi+float(i-1)*dlonm
      enddo ! i=1,nmlonp1
      dtheta = pi/(2.*90.)
!
      table(1,1) = 0. 
      table(1,2) = 0. 
      do i=2,91
        table(i,1) = table(i-1,1)+dtheta
      enddo
      do i=2,90
        table2(i,4) = tan(table(i,1))
        table(i,2) = table(i,1)
      enddo ! i=2,90
      table(91,2) = table(91,1)
!     table(91,2) = pi/2.
      do n=1,7
        do i=2,90
          table2(i,3) = table(i,2)
          table(i,2) = tan(table2(i,3))
          table2(i,5) = sqrt(r1/r0*table(i,2)+table(i,2)**(2.*(1.+alfa))
     |      /(1.+table(i,2)**2)**alfa)
          table(i,2) = table2(i,3)-(table2(i,5)-table2(i,4))*2.*
     |      table2(i,5)/(r1/r0*(1.+table(i,2)**2)+2.*table(i,2)**
     |      (2.*alfa+1.)*(1.+alfa+table(i,2)**2)/(1.+table(i,2)**2)**
     |      alfa)
        enddo ! i=2,90
      enddo ! n=1,7

! Define magnetospheric grid vars
      do i=1,nmagphrlon
        ylonmagphr(i) = magphrlon1+(i-1)*dmagphrlon
      enddo
      fac = pi/180.
      rmax = 1./(cos(magphrlat1*fac))**2
      rmin = 1./(cos(magphrlat2*fac))**2
      do i=1,nmagphrlat
        rmag = (rmax-(rmax-rmin)/(nmagphrlat-1)*real(i-1))
        ylatmagphr(i) = acos(sqrt(1./rmag))/fac
      enddo
!
!     write(6,"(/,'consdyn: ylatg =',/,(6e12.4))") ylatg
!     write(6,"(  'consdyn: ylong =',/,(6e12.4))") ylong
!     write(6,"(  'consdyn: ylatm =',/,(6e12.4))") ylatm
!     write(6,"(  'consdyn: rcos0s=',/,(6e12.4))") rcos0s
!     write(6,"(  'consdyn: dt0dts=',/,(6e12.4))") dt0dts
!
      end subroutine consdyn
!-----------------------------------------------------------------------
      end module cons_module
