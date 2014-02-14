!
!  Data assimilation/nudging module for the TIEGCM
!
! Author: Douglas P. Drob
!         Upper Atmospheric Modeling Section
!         Space Science Division
!         Naval Research Laboratory
!         4555 Overlook Ave
!         Washington, DC 20375
!
!  Date: July 2011
!
!  Description:
!
!        Uses cubic spline interpolation over a split timestep
!        to avoid evaluation of MSIS/MSIS emperical funcions at
!        every TIEGCM timestep.
!
!        This module faciliates emperical corrections
!        to the physics of the model.
!
!        The onset of the corrections of the model
!        are gradual so as not to shock the system.
!
! ====================================================

      module etimod

      use params_module,only: nlonp4,nlat,nlev,dz
      use init_module,only: istep,uthr,iter,secs,iday,iyear,glat,glon

      implicit none

      ! for now hardwire

      real,parameter  :: gtheta = 0.00015  ! 1/s  ~2 hours
      
      real,parameter  :: alphaZ = 1.0      ! vertical transition rate
      real,parameter  :: betaZ = -2.0      ! midpoint level in zp (model level)
      real,parameter  :: alphaT = 1.0      ! temporal transition rate
      real,parameter  :: betaT = 5.0	   ! transition midpoint time (days after inittime0)

      logical         :: usezlbc = .true.

      logical         :: nudgeT = .true.
      logical         :: nudgeU = .true.
      logical         :: nudgeV = .true.

      logical         :: nudgezmonly = .true.

      logical         :: useoxf = .false.

      logical         :: t0init = .true.
      integer         :: nudgeyear0 = 2008
      
      real            :: time0

      real,parameter    :: p0 = 5.e-7 !0.00054831657921422918 ! p0*exp(-Z) Hpa,Z=-7,p0=5.e-7
      integer,parameter :: nlbc = 8
      real,parameter    :: Z(0:nlbc) = (/-7.0,-6.75,-6.25,-5.75,-5.25,
     |                                   -4.75,-4.25,-3.75,-3.25/)
      real,parameter    :: zwght(0:nlbc) = (/1.0, 1.0, 1.0, 1.0, 1.0,
     |                                    1.0,1.0,0.5,0.25/)
      real              :: timeramp

      real,dimension(nlonp4,nlat,0:nlbc) :: msistlbc, msiszlbc
      real,dimension(nlonp4,nlat,0:nlbc) :: hwmulbc, hwmvlbc
      real,dimension(nlonp4,nlat,0:nlbc) :: msiso1f, msiso2f

      real,dimension(nlonp4,nlat,0:nlbc,-1:2) :: msistfx, msiszfx
      real,dimension(nlonp4,nlat,0:nlbc,-1:2) :: hwmufx, hwmvfx
      real,dimension(nlonp4,nlat,0:nlbc,-1:2) :: msiso1fx,msiso2fx

      logical           :: switchinit = .true.
      real(4)           :: sw(25) = 1.0   ! switches for MSIS harmonics
      integer           :: lastihour = 0  ! assume we start at 0 UT

      contains

      function tfilter(time)   ! real function for the temporal part of the filter
        implicit none
        real           :: tfilter
        real           :: time
        real,parameter :: sec2days = 1./3600.0/24.0
        tfilter = 1./(1. + exp(-alphaT*((time-time0)*sec2days-betaT)))  ! 0 -> 1
        tfilter = max(tfilter,0.0)
        tfilter = min(tfilter,1.0)
        return
      end function tfilter

      function zfilter(zp)      ! real function for the spatial part of the filter
        implicit none
        real           :: zfilter
        real           :: zp
        zfilter = 1. - 1./(1. + exp(-alphaZ*(zp-betaZ)))      ! 1 -> 0
        zfilter = max(zfilter,0.0)
        zfilter = min(zfilter,1.0)
        return
      end function zfilter

      ! ---------------------------------------------
      ! A subroutine to calculate the time dependent
      ! MSIS lower boundary condition
      ! ---------------------------------------------

      subroutine calclbc()

      use mpi_module,only: lon0,lon1,lat0,lat1,mytid
      implicit none

      integer     :: ccyyddd
      real(4)     :: utsec         ! universal time seconds

      integer     :: i,j,k,n
      integer     :: ihour

      real(4)     :: x1,x2,x3,x4

      real(4)     :: f107a,f107
      real(4)     :: ap(7)

      ! --------------------------------------------------
      ! On the first call set the model feature switches
      ! and build the moving window
      ! --------------------------------------------------

      if (switchinit) then

         if (nudgezmonly) then
            ! MSIS/HWM
            !sw(7:8) = 0.0          ! turn off diurnal and semidiurnal tides
            !sw(14) = 0.0           ! turn off terdiurnal tide
            !sw(10) = 0.0           ! all ut/longitude effects off
            !sw(11) = 0.0
            sw(7:14) = 0.0
            ! NRLMSISE00 
            call tselec7(sw)               
            ! HWM93
            call tselec(sw)
         endif

         ! Compute at hours -1,+0,+1,+2  (assumes start is at ut = 0 hr)
         ! Note HWM and MSIS are peridioc functions in time so uthour = -1 is ok.

         ccyyddd = iyear*1000 + iday
         do n = -1,2
             ihour = int(uthr) + n
             utsec = float(ihour)*3600.0
             call evaletmod(n,ccyyddd,utsec)
         enddo

         switchinit = .false.

      endif

      ! ------------------------------------------------
      ! Check for new interval, if so shift and append
      ! ------------------------------------------------

      ihour = int(uthr)
      if (ihour .ne. lastihour) then
          do n = -1,1
              do k = 0,nlbc
                  do j=lat0,lat1
                      do i=lon0,lon1
                          msistfx(i,j,k,n) = msistfx(i,j,k,n+1)
                          msiszfx(i,j,k,n) = msiszfx(i,j,k,n+1)
                           hwmufx(i,j,k,n) =  hwmufx(i,j,k,n+1)
                           hwmvfx(i,j,k,n) =  hwmvfx(i,j,k,n+1)
                         msiso1fx(i,j,k,n) = msiso1fx(i,j,k,n+1)
                         msiso2fx(i,j,k,n) = msiso2fx(i,j,k,n+1)
                      enddo
                  enddo
              enddo
          enddo

          utsec = float(int(uthr) + 2)*3600.0
          ccyyddd = iyear*1000 + iday
          call evaletmod(2,ccyyddd,utsec)
          if (mytid .eq. 0)
     .         print '(i7,5f10.3)',
     .            ccyyddd,utsec/3600.0,f107,f107a,ap(1),ap(2)
          lastihour = ihour
      endif

      ! ------------------------------
      ! now just interpolate in time
      ! ------------------------------

      x1 = abs(mod(uthr,1.0))
      x2 = x1 + 1.0
      x3 = 1.0 - x1
      x4 = 2.0 - x1
      x1 =  1.5*x1*x1*x1 - 2.5*x1*x1 + 1.0
      x2 = -0.5*x2*x2*x2 + 2.5*x2*x2 - 4.0*x2 + 2.0
      x3 =  1.5*x3*x3*x3 - 2.5*x3*x3 + 1.0
      x4 = -0.5*x4*x4*x4 + 2.5*x4*x4 - 4.0*x4 + 2.0

      do k = 0,nlbc
         do j=lat0,lat1
             do i=lon0,lon1
                  msistlbc(i,j,k) = x1*msistfx(i,j,k,0)
     .                            + x2*msistfx(i,j,k,-1)
     .                            + x3*msistfx(i,j,k,1)
     .                            + x4*msistfx(i,j,k,2)
                  msiszlbc(i,j,k) = x1*msiszfx(i,j,k,0)
     .                            + x2*msiszfx(i,j,k,-1)
     .                            + x3*msiszfx(i,j,k,1)
     .                            + x4*msiszfx(i,j,k,2)
                  hwmulbc(i,j,k) = x1*hwmufx(i,j,k,0)
     .                           + x2*hwmufx(i,j,k,-1)
     .                           + x3*hwmufx(i,j,k,1)
     .                           + x4*hwmufx(i,j,k,2)
                  hwmvlbc(i,j,k) = x1*hwmvfx(i,j,k,0)
     .                           + x2*hwmvfx(i,j,k,-1)
     .                           + x3*hwmvfx(i,j,k,1)
     .                           + x4*hwmvfx(i,j,k,2)
                  msiso1f(i,j,k) = x1*msiso1fx(i,j,k,0)
     .                           + x2*msiso1fx(i,j,k,-1)
     .                           + x3*msiso1fx(i,j,k,1)
     .                           + x4*msiso1fx(i,j,k,2)
                  msiso2f(i,j,k) = x1*msiso2fx(i,j,k,0)
     .                           + x2*msiso2fx(i,j,k,-1)
     .                           + x3*msiso2fx(i,j,k,1)
     .                           + x4*msiso2fx(i,j,k,2)
              enddo
          enddo
      enddo

      return

      contains

      ! ---------------------------------------------
      ! Compute the temperature and altitude of the
      ! lower TIEGCM pressure surface
      ! ---------------------------------------------

      subroutine evaletmod(n,ccyyddd,utsec)

          use oxfmod
          implicit none

          integer,intent(in)  :: n
          integer,intent(in)  :: ccyyddd
          real(4),intent(in)  :: utsec

          integer  :: i,j,k
          integer  :: ii

          real(4)     :: alt           ! geodetic altitude
          real(4)     :: slt           ! solar local time
          real(4)     :: d(9)          ! composition at altitude
          real(4)     :: t(2)          ! temperature at altitude and exosphere

          real(4)     :: apw(2) = -1.0 ! wind model Ap, -1 no DWM
          real(4)     :: w(2)          ! winds at altitude


          real(4)     :: plbc
          real(4)     :: glatin
          real(4)     :: glonin
          !real(4)    :: zm,c2,gv,reff
          real(4)     :: tmass
          real(4)     :: psio1,psio2

          real(8)           :: day,lat,zalt
          real(8)           :: mox,sbrln10o,o1sbr,o2corr
          logical,parameter :: etmodgpi = .false.

!   Commented getindex() out which is just gpi availble elsewhere in the model.
!   Assume constant F10.7(a), Ap lower boundary forcing for now.
!
!     if (etmodgpi) then
!         call getindex(ccyyddd,utsec,1,f107a,f107,ap)  
!     else
          f107a = 70.
          f107 = 70.
          ap(:) = 4.0
!     endif
          do j=lat0,lat1

              glatin = glat(j)

              do i=lon0,lon1

                  if (i .lt. 3) then
                      ii = 70 + i
                  else if (i .gt. 74) then
                      ii = i - 74
                  else
                      ii = i - 2
                  endif

                  glonin = glon(ii)
                  slt = utsec/3600. + glonin/15. ! solar local time

                  do k = 0,nlbc

                      plbc = p0*exp(-Z(k))

                      call ghp7(ccyyddd,utsec,alt,glatin,glonin,slt,
     |                    f107a,f107,ap,d,t,plbc)

                      msistfx(i,j,k,n) = t(2)
                      msiszfx(i,j,k,n) = 100.0 *1000.0 * alt

                      mox = d(2) + 2*d(4)
                      tmass = 16.0*mox + 28.0*d(3)
                      day = mod(float(ccyyddd),1000.)
                      lat = glatin
                      zalt = alt
                      
                      if (useoxf) then
                      
                        call oxf(day,lat,zalt,sbrln10o)

                        o1sbr = 10**sbrln10o
                        o2corr = (mox - o1sbr)/2.0

                        psio1 = 16.0*o1sbr/tmass
                        psio2 = 32.0*o2corr/tmass

                        msiso1fx(i,j,k,n) = psio1
                        msiso2fx(i,j,k,n) = psio2
                        
                      endif
                      
                      ! zm = alt*1000.0*100.0
                      ! c2 = cosd(2.*glatin)
                      ! gv = 980.616*(1.-.0026373*c2)
                      ! reff = 2.*gv/(3.085462e-6 + 2.27e-9*c2) !*1.0e-5
                      ! msiszlbc(i,j,k) = (reff*zm)/(reff + zm)*(gv/980.665)

                      call gws5(ccyyddd,utsec,alt,glatin,glonin,slt,
     .                    f107a,f107,apw,w)

                      hwmufx(i,j,k,n) = w(2)*100.
                      hwmvfx(i,j,k,n) = w(1)*100. ! cm/s

                  enddo
              enddo
          enddo
          return
      end subroutine evaletmod

      end subroutine calclbc

      end module etimod


