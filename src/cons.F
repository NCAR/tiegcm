!
#include <defs.h>
      module cons_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: dlat,dz,nlon,nlonp1,nlonp4,nlat,nlatp1,
     |  dlev,nlev,nmlat,nmlon,nmlonp1,zmbot,zmtop,zibot,zitop,dlon
      implicit none
!
! Define model constants. 
! Parameter constants are cons_module module data and are accessed 
!   in subroutines via use-association. 
! Derived constants are cons_module module data, and are calculated
!   in sub init_cons (contained in cons_module). 
!   Sub init_cons is called by sub init (init_mod.F).
! Parameter constants:
!
      integer,parameter :: ndays =366  ! maximum number of days in a year
      real,parameter ::
     |  dzp  = dz,           ! alias for dz (also dlev)
     |  re   = 6.37122e8,    ! earth radius (cm)                  C(51)
     |  re_inv = 1./re,      ! inverse of earth radius            C(52)
     |  avo  = 6.023e23,     ! avogadro number                    C(85)
     |  boltz = 1.38E-16,    ! boltzman's constant                C(84)
     |  p0   = 5.0e-4,       ! standard pressure                  C(81)
     |  gask = 8.314e7,      ! gas constant                       C(57)
!
! dipmin should be same as sin10 (see magfield.F):
!
#if (NLAT==36 && NLON==72)
     |  dipmin = 0.17,       ! minimum mag dip angle (5.0 deg horizontal res)
#elif (NLAT==72 && NLON==144)
     |  dipmin = 0.24,       ! minimum mag dip angle (2.5 deg horizontal res)
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
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
     |  difk(nlev+1,ndays),   ! background eddy diffusion
     |  dift(nlev+1,ndays),   ! background thermal conductivity
     |  xmue(nlev+1,ndays),   ! eddy viscosity (?)
     |  zbound          ! background low bound of Z (formerly ZBA in annual tide)
!
! Constants for dynamo and electric field calculations:
      real,parameter :: h00=9.7e6, r00=re+h00	! use mean earth radius
      real,parameter :: h0 =9.0e6, r0 =re+h0	! use mean earth radius
      real :: 
     |  dlatg, dlong, dlatm, dlonm,
     |  ylatm(nmlat),    ! magnetic latitudes (radians)
     |  ylonm(nmlonp1),  ! magnetic longitudes (radians)
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
! Critical colatitude limits (15,30 deg) for use of Heelis potential in dynamo:
! These may change with "dynamic crit" in Weimer run, or if reading amie data 
      real :: crit(2) = (/0.261799387, 0.523598775/)
!
! Kut is used in filtering longitudinal waves (see filter.F):
!
#if (NLAT==36 && NLON==72)                 /* 5.0 deg horizontal resolution */
      integer,parameter :: kut(nlat) =
     |  (/1,2,3,5,6,7,9,10,11,13,14,15,17,17,17,17,17,17,17,17,17,17,17,
     |    17,15,14,13,11,10,9,7,6,5,3,2,1/)
#elif (NLAT==72 && NLON==144)              /* 2.5 deg horizontal resolution */
       integer,parameter :: kut(nlat) =
     |   (/1  ,1  ,2  ,2  ,4  ,4  ,8  ,8  ,10 ,10 ,12 ,12,
     |     15 ,15 ,18 ,18 ,22 ,22 ,26 ,26 ,30 ,30 ,32 ,32,
     |     34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34,
     |     34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34,
     |     32 ,32 ,30 ,30 ,26 ,26 ,22 ,22 ,18 ,18 ,15 ,15,
     |     12 ,12 ,10 ,10 ,8  ,8  ,4  ,4  ,2  ,2  ,1  ,1/)
!      integer,parameter :: kut(nlat) =
!    |   (/0  ,0  ,1  ,2  ,3  ,4  ,5  ,6  , 7 , 8 , 9 ,10,
!    |     15 ,15 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17,
!    |     17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17,
!    |     17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17,
!    |     17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,15 ,15,
!    |     10 , 9 , 8 , 7 ,6  ,5  ,4  ,3  ,2  ,1  ,0  ,0/)
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
!
! If check_exp is set true, certain routines will use expo() (util.F)
! instead of exp(). expo checks for out of range arguments to the
! exponential, substituting large or small results if the argument
! is out of range. This avoids NaNS fpe's, but degrades performance.
! It will also produce slightly different results.
!
#ifdef DEBUG
      logical,parameter :: check_exp = .true.
#else
      logical,parameter :: check_exp = .false.
#endif
!
! Special pi for mag field calculations. If pi=4.*atan(1.) and code is
! linked with -lmass lib, then the last 2 digits (16th and 17th) of pi
! are different (56 instead of 12), resulting in theta0(j=49)==0., which 
! is wrong (should be .1110e-15).
!
      real,parameter ::
     |  pi_dyn=3.14159265358979312

      integer,parameter :: difhor=1
!
      contains
!-----------------------------------------------------------------------
      subroutine init_cons
      use input_module,only: step,eddy_dif
      use params_module,only: tgcm_version
!
! Set derived constants (this is called from sub init in init_module)
!
! Local:
      real :: z,expdz,phi
      real :: omega = 7.292E-5
      integer :: m,k,i,j,nfsech,js
      integer :: calday(ndays)         ! calenday day in a year
      real :: coeff(9)
      real :: theta(ndays)
!
! Begin calculation of eddy diffusivity:
      if (eddy_dif > 0) then ! use DOY-dependent eddy diffusion
!
! The 9 coefficient of the semiannual eddy diffusivity parameterization:
        coeff(:)=(/4.06e-06, -8.77e-07, -2.28e-06, 1.77e-06, 2.15e-06,
     |            -3.05e-07, -2.66e-07,  4.08e-07, 1.59e-07/)
!
! set up calendar days
        do m=1,ndays
           calday(m)=m
        enddo
!
! set up the theta array
        do m=1,ndays
           theta(m)=2.*3.14*(calday(m)-1.)/(ndays-1.)
        enddo
!      
! calculate difk, dift, and xmue at first model level(bottom) in s^-1 for all
! the days in a year
!
        do m=1,ndays
           difk(1,m) = coeff(1) 
     |             + coeff(2)*sin(theta(m))+coeff(3)*cos(theta(m))
     |             + coeff(4)*sin(2.*theta(m))+coeff(5)*cos(2.*theta(m)) 
     |             + coeff(6)*sin(3.*theta(m))+coeff(7)*cos(3.*theta(m)) 
     |             + coeff(8)*sin(4.*theta(m))+coeff(9)*cos(4.*theta(m)) 
           dift(1,m) = difk(1,m)/prndtl
           xmue(1,m) = difk(1,m)
        enddo
!
!  calculate the rest of the levels
!
        expdz = exp(-dlev)
        do m=1,ndays
           do k=2,nlev+1
              difk(k,m) = difk(k-1,m)*expdz
              dift(k,m) = dift(k-1,m)*expdz
              if (k == nlev+1) then
                 xmue(k,m) = difk(k-1,m)
              else
                 xmue(k,m) = difk(k,m)
              endif
           enddo
        enddo
!
! Use constants if eddy_dif <= 0:
!
      else ! eddy_dif is off
        expz(:) = 0. ! init
        z = zibot+.5*dlev
        expz(1) = exp(-z)
        expdz = exp(-dlev)

!       difk(1,:) = 5.0e-6
!       dift(1,:) = 5.0e-6/prndtl
!       xmue(1,:) = 5.0e-6
!
! 4/23/09 btf: amie mods as per Gang:
! Note this amie code assumes eddy_dif==0 (i.e., DOY eddy diffusion is off)
!
! reduce Eddy diffusion by a factor of 2 - G. Lu, Jul 24, 2006
!     difk(1,:) = 2.5e-6
! Increase Eddy diffusion by a factor of 2 - G. Lu, Dec 20, 2006 - Don't seem to affect Nmf2 very much
!     difk(1,:) = 10.0e-6
!     dift(1,:) = 5.0e-6/prndtl
!     xmue(1,:) = 5.0e-6
! Change made by Roble on May 13, 2007 for Charley Barth
!     difk(1,:) = 2.0e-6
!     dift(1,:) = 2.0e-6/prndtl
!     xmue(1,:) = 2.0e-6
! Change made for April 2002 (when the modeled Rho was smaller than CHAMP)
!     difk(1,:) = 1.0e-6
!     dift(1,:) = 1.0e-6/prndtl
!     xmue(1,:) = 1.0e-6
! Change made for Sep 2005 (when the modeled Rho was larger than CHAMP)
      difk(1,:) = 2.0e-5
      dift(1,:) = 2.0e-5/prndtl
      xmue(1,:) = 2.0e-5

        do k=2,nlev
          expz(k) = expz(k-1)*expdz
          difk(k,:) = difk(k-1,:)*expdz
          dift(k,:) = dift(k-1,:)*expdz
          xmue(k,:) = difk(k,:)
        enddo
        difk(nlev+1,:) = difk(nlev,:)*expdz
        dift(nlev+1,:) = dift(nlev,:)*expdz
        xmue(nlev+1,:) = difk(nlev,:)
      endif  ! eddy_dif on or off
!
      if (eddy_dif > 0) then
        write(6,"('init_cons: eddy_dif=',i3)") eddy_dif
        write(6,"('min,max: difk=',2e12.4,' dift=',2e12.4,' xmue=',
     |    2e12.4)") minval(difk),maxval(difk),minval(dift),maxval(dift),
     |              minval(xmue),maxval(xmue)
        write(6,"('difk(:,1)=',/,(6e12.4))") difk(:,1)
        write(6,"('dift(:,1)=',/,(6e12.4))") dift(:,1)
        write(6,"('xmue(:,1)=',/,(6e12.4))") xmue(:,1)
      endif
!
! end calculation of eddy diffusivity
!
      pi = 4.*atan(1.)                ! C(110)
      rtd = 180./pi                   ! radians to degrees
      dtr = pi/180.                   ! degrees to radians
      dphi = pi/float(nlat)           ! C(2)
      dphi_2div3 = 2./(3.*dphi)       ! C(12)
      dphi_1div12 = 1./(12.*dphi)     ! C(13)
      dlamda = 2.*pi/float(nlon)      ! C(1)
      dlamda_2div3 = 2./(3.*dlamda)   ! C(10)
      dlamda_1div12 = 1./(12.*dlamda) ! C(11)
      zbound = 136.291/sqrt(2.)*1.e5  ! background lower boundary of Z (cm)
!
! expz(kmax) is exp(-zp) at midpoints:
! expz (will replace EXPS) (expz(nlev+1) not used).
!
      expz(:) = 0. ! init
!
! bottom midpoint z = zibot + 1/2 deltaz (deltaz==dz==0.5 or 0.25)
! (zibot and dz are in params.h)
      z = zibot+.5*dlev
      expz(1) = exp(-z)
      do k=2,nlev
        expz(k) = expz(k-1)*expdz
      enddo
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
!
! Set gravity according to lower boundary:
      grav = 870. ! (is 945. in time-gcm)
      dzgrav = grav/gask   ! C(65)
!
! Set dynamo constants:
      call consdyn
!
! Report to stdout:
!     write(6,"(/,'Model name    = ',a)") tgcm_name
      write(6,"(  'Model version = ',a)") tgcm_version
      write(6,"(/,'Set constants:')")
      write(6,"('  nlat=',i3,' nlon=',i3,' nlev=',i3)") nlat,nlon,nlev
      write(6,"('  dz= ',f5.2)") dz
      write(6,"('  dlat=',f6.2,' dlon=',f6.2)") dlat,dlon
      write(6,"('  zbound (cm) = ',e12.4)") zbound
      write(6,"('  zmbot, zmtop = ',2f8.3,
     |  ' (bottom,top midpoint levels)')") zmbot,zmtop
      write(6,"('  zibot, zitop = ',2f8.3,
     |  ' (bottom,top interface levels)')") zibot,zitop
      write(6,"('  dt = ',f8.2,' secs')") dt
      write(6,"('  grav = ',f10.2)") grav
      write(6,"('  freq_3m3 = ',e12.4,'  freq_semidi=',e12.4)")
     |  freq_3m3,freq_semidi
      write(6,"('  dipmin = ',f8.3)") dipmin
      write(6,"('  check_exp = ',l1)") check_exp
      write(6,"('  kut (for filtering) = ',36i3)") kut
!
      end subroutine init_cons
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
      real :: rmin,rmax,rmag
!
! Set grid deltas:
      dlatg = pi/float(nlat)
      dlong = 2.*pi/float(nlon)
      dlatm = pi_dyn/float(nmlat-1) ! note use of pi_dyn
      dlonm = 2.*pi_dyn/float(nmlon) 
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
! If dynamo <= 0 -> no dynamo
! If dynamo >= 1 -> dynamo is called
!
        if (dynamo > 0) then
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
