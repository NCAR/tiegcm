#include "defs.h"F
      module amie_module
!      
! Module used to read data from the AMIE outputs (POT,mean energy, 
!   and energy flux).
!
      use params_module,only: nlon,nlat,nlev,nlonp1,nlonp4,nlatp1,
     |  nlevp1,nmlon,nmlonp1,nmlat
      use addfld_module,only: addfld
      use input_module,only: iamie
      use aurora_module,only: gekev,gmwpm2	! (nlonp1,0:nlatp1) in geog coords
      use cons_module,only: pi,dtr,rtd,ylonm,ylatm
      implicit none
!

! Define parameters for AMIE input data file:
!  iamie in namelist: iamie=1 for Gang Lu HAO AMIE, iamie=2 for Geoff Crowley ASTRA AMIE,
!   iamie=3 for Aaron Ridley U MI AMIE, iamie=4 for SWMF in U MI AMIE form 
!   iamie=5 for Dartmouth SuperDARN (SD) for 20 50-88mlat and 36 40min MLT (0-2320MLT) ~10 min
!   iamie=6 for VT SuperDARN (VTSD) for 40 50-89mlat and 180 2deg mlon (0-358mlon) ~10 min

! colath.F restricts crit(1) between 15-30 colat so crit(2) is between 30 and 45)

! FIXME!  Want to change these so do not have to revise them for each run
!  mxtimes only used for ASTRA and for SuperDARN, and inmin used only for SuperDARN
!  Do not know why UMI worked with ASTRA (larger) parameters for first day 5deg grid
      integer,parameter :: 
!  SWMF in 3 1-day files of 289 times from 0-0(24)UT every 15 min
!     |  mxtimes = 289,  ! maximum number of times of AMIE data per file (Dec 13,14,15 5 min) + 1 for 0UT 
!     |  lonmx = 180,    ! maximum number of longitudes of AMIE data (1 less than grid), but SWMF is lonp1!
!     |  ithtrns = 89,   ! corresponding to trans lat 1 mlat number of AMIE lats minus 1
!     |  ithmx = 91,     ! max number of AMIE mlats (90/(5/3)=54+1) to equator
!     |  ifpole = 1,	 ! ifpole=1,0 if do,not have pole as endpoints
!     |  inmin = 15,     ! number of minutes between inputs, used only for SD
! Geoff Crowley's amie files:
!  ASTRA AMIE is 40min MLTS 00-24(dup) or 37 times, 31 lats 90-40 1.67 deg, every 5 min
!     |  mxtimes = 1153, ! maximum number of times of AMIE data per file (Dec 13-16 5 min) + 1 for 0UT Dec 17
!     |  lonmx = 36,     ! maximum number of longitudes of AMIE data (1 less than grid), but ASTRA is lonp1!
!     |  ithtrns = 30,   ! corresponding to trans lat 40 mlat number of AMIE lats minus 1
!     |  ithmx = 55,     ! max number of AMIE mlats (90/(5/3)=54+1) to equator
!     |  ifpole = 1,	 ! ifpole=1,0 if do,not have pole as endpoints
!     |  inmin = 5,      ! number of minutes between inputs, used only for SD
! Aaron Ridley's amie files:
!  Ridley AMIE is 2deg 90-44 or 24 lats + 15 lats extended to 14 mlat where potential is zero
!     |  mxtimes = 5881, ! maximum number of times of AMIE data per file (Dec 13-16 1min)
!     |  lonmx = 24,     ! maximum number of longitudes of AMIE data (1 less than grid)
!     |  ithtrns = 24,   ! Aaron: corresponding to trans lat 42 mlat (ie, 90 to 44 dmlat=2) 
!     |  ithmx = 46,     ! max number of AMIE mlats (90/2=45+1) to equator
!     |  ifpole = 1,	 ! ifpole=1,0 if do,not have pole as endpoints
!     |  inmin = 1,      ! number of minutes between inputs, used only for SD
! Gang Lu's AMIE files (usually 1 per day, so need to run each day separately)
!     |  mxtimes = 1440 ! maximum number of times of AMIE data per daily file
!     |  lonmx = 36,    ! maximum number of longitudes of AMIE data
!     |  ithtrns = 30,  ! corresponding to trans lat 40-deg
!     |  ithmx = 55,    ! maximum number of latitudes of AMIE data
!     |  ifpole = 1,	 ! ifpole=1,0 if do,not have pole as endpoints
!     |  inmin = 1,      ! number of minutes between inputs, used only for SD
!  Dartmouth SuperDARN 200612130000-200612170000 north and south 4 days 10min UT 6*24*4=576, 
!    36 40min MLT 00.00-23.33, 20 mlat 2 deg 50-88deg
!     |  mxtimes = 578,  ! maximum number of times of AMIE data per file (Dec 13-16 10 min) + 2 for endpts
!     |  lonmx = 36,     ! maximum number of longitudes of AMIE data (1 less than grid)
!     |  ithtrns = 20,   ! corresponding to trans lat 50 mlat number of AMIE lats minus 1 + pole
!     |  ithmx = 46,     ! max number of AMIE mlats (90/2deg=45+1) to equator
!     |  ifpole = 0,	 ! ifpole=1,0 if do,not have pole as endpoints
!     |  inmin = 10,     ! number of minutes between inputs, used only for SD
!  Virginia Tech SuperDARN north and south 4 days 10min UT 142*4=568 missing 20 min 2350-0010!
!   180 2 deg 0-358mlon, 40 mlat 1 deg 50-89deg
! After 23 Apr 2014, VT SD has begin and end times of 0002(-0004) and 2356(-2358) since
!   inmin cadence is really the default 2 min chunks every inmin minutes
!   Sometimes, especially in the SH, the end 2min chunk of time is 2354(-2356) instead of 2356.
!  If SH missing 2356UT for inmin=2, then should duplicate 2354ut, revise UT to 2356ut, 
!   and cat it on end of daily file to have the same number of UT times as the NH
!     |  mxtimes = 570,  ! maximum number of times of AMIE data per file (Dec 13-16 10 min - 20min ~0UT) + 2 for endpts
!     |  mxtimes = 722,  ! maximum number of times of AMIE data per file 72_20m*10d=720+2
!     |  mxtimes = 1438,  ! maximum number of times of AMIE data per file 718_2m*2d=1436+2
!     |  mxtimes = 2874,  ! maximum number of times of AMIE data per file 718_2m*4d=2872+2
     |  mxtimes = 218,  ! maximum number of times of AMIE data per file 72_20m*3d=216+2
     |  lonmx = 180,    ! maximum number of longitudes of AMIE data (1 less than grid)
     |  ithtrns = 40,   ! corresponding to trans lat 50 mlat number of AMIE lats minus 1 + pole
     |  ithmx = 91,     ! max number of AMIE mlats (90/1deg=90+1) to equator
     |  ifpole = 0,	! ifpole=1,0 if do,not have pole as endpoints
!     |  inmin = 10,     ! number of minutes between inputs, used only for SD
     |  inmin = 20,     ! number of minutes between inputs, used only for SD
!     |  inmin = 2,     ! number of minutes between inputs, used only for SD
! jmxm is odd for ifpole=1, even for ifpole=0 (SD) which includes 2 polar values also
     |  jmxm = 2*ithmx-ifpole ! maximum number of global latitudes 
!                             !  (only one 0 mlat ifpole=1, or inbetween 2 pts for ifpole=0)

      real :: alatm(2*ithmx-ifpole),alonm(lonmx+1)
!  lonp1 and latp1 are AMIE dim of lon (37 or lonmx+1) and lat (31 or ithmx+1 where 30 is 90 to 45.67 in dmlat=1.67)
      integer :: lonp1,latp1
! Define AMIE output fields
      real :: 
     |  tieefx(nmlonp1,nmlat),tieekv(nmlonp1,nmlat)
!
! Define fields for AMIE input data file:
! electric potential in Volt
! mean energy in KeV
! energy flux in W/m^2 (wrong for U MI AMIE - is mW/m2 or ergs/cm2-s)

! NOT USED from Gang Lu's HAO AMIE for aurora.F
      real :: cusplat_nh_amie, cuspmlt_nh_amie, cusplat_sh_amie,
     |  cuspmlt_sh_amie
!
! 6/12: added by Emery for binary Ridley AMIE
      real    :: del,xmlt,dmlat,dlatm,dlonm,dmltm,rot,mltsec
      real, dimension(ithtrns*2+2) :: degmlat
      real, dimension(lonmx+1,ithtrns*2+2) :: mlts, lats, 
     |  TempPotential, AveEOut, EfluxOut

! 4/13: added by Emery for ascii ASTRA AMIE
      real,dimension(mxtimes,2*ithtrns+2,lonmx+1),private :: 
     |    epkv,wpm2,ekev
! 4/14 bae: save daynum amie_day, and minutes in day amie_min 
      integer,dimension(mxtimes),private :: amie_day
      real,dimension(mxtimes),private :: amie_min
 


      contains 
!-----------------------------------------------------------------------
      subroutine init_amie
!
! Was called from tgcm.F, but now called in getamie for istep<=1
! (this is not in init.F to avoid circular dependencies)
!
      use input_module,only: amienh,amiesh
      use aurora_module,only: calc_hp

      integer,external :: nextlu
!
! 6/12 added by Emery for Ridley binary AMIE files
      character (len=100), dimension(100) :: Lines
      character (len=90) :: char90
      integer :: iError, i,j, lu_amie,ls_amie,nt,ihem,jt
      integer :: iy,imo,ida,ih,imin,isec
      real :: mlt,mlat,epot,eflx,qj,energy

! 3/14 added by Emery for SuperDARN files
      integer :: iy2,imo2,ida2,ih2,imin2,isec2,ntim,nt1,nt2,nda,jutsec
      integer :: nl,ifstop,j1,j2
      real :: lblat,glat,glon,mlon,depot,polekvs,polekvn,delmin
      character*1 :: char1
! Local dimension for temporary print-outs
      real :: epot_nsh,epot_xsh,epot_nnh,epot_xnh,hpsh,hpnh
      real :: mlat_nsh,mlat_xsh,mlat_nnh,mlat_xnh
      real :: mlt_nsh,mlt_xsh,mlt_nnh,mlt_xnh
      real :: read_mlt(lonmx+1),read_mlat(2*ithtrns+2)

! 6/12: Set AMIE NH or SH grid dimensions lonp1 and latp1 (to go 90 to 40 mlat in dmlat=2 for Aaron)
!  latp1 includes poles for SD which have to be defined from highest latitude
      latp1 = ithtrns + 1
      lonp1 = lonmx + 1

C     ****     SET GRID SPACING DLATM, DLONG, DLONM
C  DMLAT=lat spacing in degrees of AMIE apex grid 
C 6/12: from amie.nc.out in ~ganglu/amie, have 31 lats w dmlat=1.67 from 90 to 40 mlat
C       and 37 mlts w dmltm=0.67h from 0 to 24MLT, where lonmx=36 (NOT 37!)
C   18*3=54 so 90/54=1.67=dmlat (ithtrns=30 to 40 mlat, ithmx=55 to equator, jmxm=2*ithmx-1)
      DMLAT = 180. / FLOAT(jmxm-2+ifpole)	! AMIE grid spacing in deg (should be 2 deg Aaron, 1.67 deg ASTRA)
      DLATM = DMLAT*dtr
      DLONM = 2.*pi/FLOAT(lonmx)
      DMLTM = 24./FLOAT(lonmx)		! AMIE grid spacing in MLT (should be 1h Aaron, 0.67h ASTRA)
      mltsec = 86400/lonmx
  
      write(6,"(1x,'dmlat dmlon dmlt =',3f8.2)") dmlat,dlonm/dtr,dmltm
       j2 = 1
       j1 = 1
       if (ifpole==0) then
        j2 = 2
        j1 = 2        
        degmlat(1) = -90.
        degmlat(latp1+1) = 90.
       do i=1,lonp1
        lats(i,1) = -90.0
        lats(i,1+latp1) = 90.0
        mlts(i,1) = (i-1)*dmltm
        mlts(i,1+latp1) = (i-1)*dmltm
       enddo
       endif
      do j=j2,latp1
! In 'mirror', the min/max potentials are flipped from NH to SH 
!  Want lats from pole to 90-ithtrns*dmlat in both hem ifpole=1
       degmlat(j+latp1) = 90.0 - (1-ifpole)*dmlat/2. - (j-j1)*dmlat
       degmlat(j) = -90.0 + (1-ifpole)*dmlat/2. + (j-j1)*dmlat
! lats,mlts required for U MI AMIE, and used as a check for others
       do i=1,lonp1
       lats(i,j+latp1) = 90.0 - (1-ifpole)*dmlat/2. - (j-j1)*dmlat		
       lats(i,j) = -90.0  + (1-ifpole)*dmlat/2. + (j-j1)*dmlat
! mlts NOT used for iamie=6 (SD_VT) since only have mlon!
        mlts(i,j) = (i-1)*dmltm
        mlts(i,j+latp1) = (i-1)*dmltm
       enddo
      enddo   
!  Printout messed up if nproc>1
      write (6,"(1x,'degmlat=',10f7.2)") degmlat
      write (6,"(1x,'mlts=',10f7.2)") (mlts(i,1),i=1,lonp1)   

! Set whole globe grid for rdgrd2.F, including 2 polar values for ifpole=0
      alatm(1) = -PI/2. 
      alatm(jmxm) = PI/2.
      j1 = 2
      if (ifpole==0) then
       j1 = 3
       alatm(2) = -PI/2. + dlatm/2. 
       alatm(jmxm-1) = PI/2. - dlatm/2.
      endif
      do j=j1,ithmx
        alatm(j) = alatm(j-1)+dlatm
        alatm(jmxm+1-j) = alatm(jmxm+2-j)-dlatm
      enddo
! FIX so have 1 more SD alonm at -180 then -179,-177 etc
      alonm(1) = -pi + (1-ifpole)*dlonm/2.
      alonm(lonp1) = pi
      do i=2,lonmx
        alonm(i) = alonm(i-1) + dlonm
      enddo
! Ensure alatm(1,jmxm),alonm(1,lonmx=1) same as distorted yaltm,ylonm TIEGCM mag grid
      if (alatm(1) > ylatm(1)) alatm(1) = ylatm(1)
      if (alatm(jmxm) < ylatm(nmlat)) alatm(jmxm) = ylatm(nmlat)
      if (alonm(1) > ylonm(1)) alonm(1) = ylonm(1)
      if (alonm(lonp1) < ylonm(nmlonp1)) alonm(lonp1) = ylonm(nmlonp1)
      write (6,"(1x,'alatm=',10f7.2)") alatm*rtd
      write (6,"(1x,'alonm=',10f8.2)") alonm*rtd
      write (6,"(1x,'TIEGCM ylatm=',10f7.2)") ylatm*rtd
      write (6,"(1x,'TIEGCM ylonm=',10f8.2)") ylonm*rtd

      if (len_trim(amienh) > 0) then
        write(6,"('Reading AMIENH file ',a)") trim(amienh)
        Lines(2) = trim(amienh)
      endif

! 04/13: Do different things for U MI and ASTRA AMIE

! ASTRA AMIE:
      if (iamie==2) then
          lu_amie = nextlu()
        write (6,"(1x,'lu_amie=',i3)") lu_amie
          open(lu_amie,file=trim(amienh),status='old')
! Skip header line
        read (lu_amie,"(a)") char90
! Read every 5 min 4 days of data over 0000to2400 mlts from 90to40, -40to-90 mlat
!  Store in order SH -90 to ~-40 mlat, NH +90 to ~+40mlat
        do nt=1,mxtimes-1
! ASTRA has both MLTS=0 and 24 or i=1 and lonp1 which are dups
         do i=1,lonp1
          do ihem=1,2
           do j=1,latp1
             read(lu_amie,"(6i4,2f8.1,1x,4f12.4)") iy,imo,ida,ih,imin,
     |        isec,mlt,mlat,epot,eflx,qj,energy
             read_mlt(i) = mlt
	     if (ihem==1 .and. j==1 .and. i==1) then
              call cvt2dn (6,iy,nda,imo,ida)
              amie_day(nt) = nda
              amie_min(nt) = ih*60 + imin + isec/60.
!              write (6,"(1x,'nt dayn hr min sec minpd = ',5i5,f9.2)") 
!     |         nt,nda,ih,imin,isec,amie_min(2)
             endif
             if (ihem==1) then
              epkv(nt,j+latp1,i) = epot
              wpm2(nt,j+latp1,i) = eflx
              ekev(nt,j+latp1,i) = energy
              read_mlat(j+latp1) = mlat
             else
              epkv(nt,latp1+1-j,i) = epot
              wpm2(nt,latp1+1-j,i) = eflx
              ekev(nt,latp1+1-j,i) = energy
              read_mlat(latp1+1-j) = mlat
             endif
           enddo
          enddo
         enddo		! do i=1,lonp1
! TEMP find HP and min-max and print out CP and time to fort.52
          epot_nsh = 0.
          epot_xsh = 0.
          epot_nnh = 0.
          epot_xnh = 0.
          do j=1,latp1
           do i=1,lonmx
             epot_nsh = min(epkv(nt,j,i),epot_nsh)
             epot_xsh = max(epkv(nt,j,i),epot_xsh)
             epot_nnh = min(epkv(nt,j+latp1,i),epot_nnh)
             epot_xnh = max(epkv(nt,j+latp1,i),epot_xnh)
             EFluxOut(i,j) = wpm2(nt,j,i) 
             EFluxOut(i,j+latp1) = wpm2(nt,j+latp1,i) 
           enddo
          enddo
!  Calculate hpi_s,nh_amie in GW
       call calc_hp (latp1*2,lonmx,lonmx+1,degmlat,EFluxOut,hpsh,hpnh) 
         if (nt==1) then
        write (52,"('ASTRA AMIE Dec06 CP and HP'/'  nt year',
     | ' nda mnpd cp_sh_kV cp_nh_kV ctpoten hp_sh_GW hp_nh_GW power')")
         endif
          write (52,"(i4,2i5,7f9.2)") nt,iy,amie_day(nt),
     |     amie_min(nt),epot_xsh-epot_nsh,
     |     epot_xnh-epot_nnh,0.5*(epot_xsh-epot_nsh+epot_xnh-epot_nnh),
     |     hpsh,hpnh,0.5*(hpsh+hpnh)
! wrap-around loc 1=lonp1 (mlt=0 and 24) not necessary since MLT=24 is also read in. 
!  Check that read_mlt and read_mlat same as mlts and degmlat
         if (nt==1) then
          ifstop = 0
          do j=1,latp1*2
           if (abs(lats(1,j)-read_mlat(j)) .gt. 0.1) ifstop=1
          enddo
          do i=1,lonmx
           if (abs(mlts(i,1)-read_mlt(i)) .gt. 0.1) ifstop=1
          enddo
          if (ifstop==1) then
           write (6,"(1x,'STOP: lats.ne.read_mlat or mlts.ne.read_mlt'
     |      /(1x,10f8.2))") (lats(1,j),read_mlat(j),j=1,latp1*2),
     |      (mlts(i,1),read_mlt(i),i=1,lonmx)
           stop
          endif
         endif   	! if (nt==1) then    
        enddo		! do nt=1,mxtimes-1
        delmin = amie_min(2) - amie_min(1)
!  TEMP!! Extrapolate to 24 UT at end of time (ie to 17 Dec 2006 from 16 Dec) for ASTRA
        amie_day(mxtimes) = amie_day(mxtimes-1)+1
        amie_min(mxtimes) = 0
        do i=1,lonmx+1
         do j=1,latp1*2
          epkv(mxtimes,j,i) = epkv(mxtimes-1,j,i) 
     |      + (epkv(mxtimes-1,j,i)-epkv(mxtimes-2,j,i))
          wpm2(mxtimes,j,i) = wpm2(mxtimes-1,j,i) 
     |      + (wpm2(mxtimes-1,j,i)-wpm2(mxtimes-2,j,i))
          ekev(mxtimes,j,i) = ekev(mxtimes-1,j,i) 
     |      + (ekev(mxtimes-1,j,i)-ekev(mxtimes-2,j,i))
         enddo
        enddo
       write (6,"(1x,'Reading AMIE times from 1 to mxtimes (dn,mpd)', 
     |   ' every delmin minutes =',i4,3(1x,i4,f9.2),f7.2)")
     |    mxtimes,amie_day(1),amie_min(1),amie_day(mxtimes-1),
     |    amie_min(mxtimes-1),amie_day(mxtimes), amie_min(mxtimes),
     |    delmin
      endif		! if (iamie==2)

! U MI AMIE:
      if (iamie==3 .or. iamie==4) then
 
      if (len_trim(amiesh) > 0) then
        write(6,"('Reading AMIESH file ',a)") trim(amiesh)
        Lines(3) = trim(amiesh)
      else
	Lines(3) = "mirror"
      endif
       
      Lines(1) = "#AMIEFILES"
      Lines(4) = ""
      Lines(5) = ""
      Lines(6) = ""

      Lines(7) = "#DEBUG"
      Lines(8) = "2"
      Lines(9) = "0"
      Lines(10) = ""
      Lines(11) = "#END"

      write(*,*) "Calling IE_set_inputs"

      call IE_set_inputs(Lines)

      write(*,*) "=> Setting up Ionospheric Electrodynamics"

      call IE_Initialize(iError)

      write(*,*) "Setting nmlts and nlats"
! Aaron Ridley says the original AMIE comes as 1min NH arrays of 25 MLTs (0-24 - NOT 1-25!!!)
!  every 2 deg from 90N to 44N (24!)- However, the code interpolates to whatever grid is 
!  set below as lats and mlts for input with dimensions set in IO_SetnMLT(LAT)s
! Can find the dimensions of AMIE in readAMIEoutput.f90

      call IO_SetnMLTs(lonmx+1)

      call IO_SetnLats(ithtrns*2+2)

      write(*,*) "Setting grid"

      call IO_SetGrid(mlts,lats, iError)
      endif	!       if (iamie==3 .or. iamie==4) then

! SuperDARN Dartmouth and Virginia Tech:
      if (iamie==5 .or. iamie==6) then
          lu_amie = nextlu()
          open(lu_amie,file=trim(amienh),status='old')
          ls_amie = nextlu()
          open(ls_amie,file=trim(amiesh),status='old')
        write (6,"(1x,'iamie and unit numbers for SuperDARN N,S=',
     |    3i3)") iamie,lu_amie,ls_amie
! Read 2 to mxtimes-1 data 
        do nt=2,mxtimes-1 
         if (iamie==5) then
! Read time interval year month day hr min sec begin and end
          read(lu_amie,*) iy,imo,ida,ih,imin,isec,
     |      iy2,imo2,ida2,ih2,imin2,isec2,lblat
          read(lu_amie,*) ntim
          read(ls_amie,*) iy,imo,ida,ih,imin,isec,
     |      iy2,imo2,ida2,ih2,imin2,isec2,lblat
          read(ls_amie,*) ntim
          if (ntim .ne. lonmx*latp1) stop  
         else
! iamie=6 for VT SD:  Skip over 12 lines, where 5th line gives first time, 9th gives min,max,CP
          do nl=1,4
           read (lu_amie,"(a1)") char1
           read (ls_amie,"(a1)") char1
          enddo
!  Read begin time for VT SD
          read (ls_amie,"(2x,i4,5(1x,i2))") iy,imo,ida,ih,imin,isec
! Read NH last since VT first time 0002 or last time 2356 can be different for SH w less data
          read (lu_amie,"(2x,i4,5(1x,i2))") iy,imo,ida,ih,imin,isec
          do nl=6,12
           read (lu_amie,"(a1)") char1
           read (ls_amie,"(a1)") char1
          enddo
         endif	! if (iamie==5) then,else
!  use midpoint time and convert to daynumber
              call cvt2dn (6,iy,nda,imo,ida)
              amie_day(nt) = nda
              amie_min(nt) = ih*60 + imin + isec/60. + inmin/2.	!  Corrected for inmin/2 midpt
             if (iamie==6) then
!  Correct for endpoints of day at 0002-0004 and 2356-2358 for any inmin but 2 min (the default chunk size)
!  If SH missing 2356UT for inmin=2, then should duplicate 2354ut, revise UT to 2356ut, 
!   and cat it on end of daily file to have the same number of UT times as the NH
              if (ih==0 .and. imin==2 .and. inmin.ne.2)  
     |         amie_min(nt) = ih*60 + imin + (inmin-2)/2.
              if (ih==23 .and. imin==60-2*inmin .and. inmin.ne.2)
     |         amie_min(nt) = ih*60 + imin + (56-imin)/2.
              if (ih==23 .and. imin==56 .and. inmin.ne.2) 
     |         amie_min(nt) = ih*60 + imin + 1.
             endif
             if (amie_min(nt).ge.1440.) then
              write (6,"(1x,'nt dayn hr min sec minpd(>=1440) = ',5i5,
     |         f9.2)") nt,nda,ih,imin,isec,amie_min(nt)
              stop
             endif
        if (iamie==5) then
! Read Dartmouth SD  
         do i=1,lonmx
           do j=1,latp1-1
! NH read 50 to 88 mlat, 0 to 23.33 mlt
             read(lu_amie,"(f10.1,6f11.1)") epot,depot,
     |        glat,glon,mlat,mlon,mlt	
!  Save 88 to 50 mlat and allow for polar value at j=latp1 (1=latp1+1-j, latp1+1=2*latp1+1-latp1)
             epkv(nt,2*latp1+1-j,i) = epot*1.e-3
!  NOTE:  Need to revise mlat,mlt from begin point to mid-point!
	     read_mlat(2*latp1+1-j) = mlat + dmlat/2.
             read_mlt(i) = mlt + dmltm/2. 
! SH read -50 to -88 mlat, 0 to 23.33 mlt
             read(ls_amie,"(f10.1,6f11.1)") epot,depot,
     |        glat,glon,mlat,mlon,mlt	     
              epkv(nt,latp1+1-j,i) = epot*1.e-3
	     read_mlat(latp1+1-j) = mlat - dmlat/2.
             read_mlt(i) = mlt + dmltm/2. 
           enddo
         enddo		! do i=1,lonmx
        else
! Read VT SD from 50-89mlat and save 89.5 to 50.5 with j=latp1 (1=latp1+1-latp1,latp1+1=2*latp1=1-latp1) for pole
         do j=1,latp1-1
          do i=1,lonmx
!  Read NH 0-358 mlon
           read (lu_amie,"(17x,f13.4,f18.3,72x,f21.4,8x,i4,5(1x,i2))")
     |      mlat,mlon,epot,iy2,imo2,ida2,ih2,imin2,isec2
           jutsec = ih2*3600 + imin2*60 + isec2 + mltsec/2
! rot=15*ihr+min/4 - 70 is mlon at local midnight 
           rot = jutsec*15./3600. - 70.
!  Save 89 to 60 mlat and -180 to +180 mlon (not 1:2:359 mlon)
	     read_mlat(2*latp1+1-j) = mlat + dmlat/2.
             if (i>lonmx/2) then
              read_mlt(i-lonmx/2) = mlon+dlonm*rtd/2. - 360.
              epkv(nt,2*latp1+1-j,i-lonmx/2) = epot*1.e-3
             else
              read_mlt(i+lonmx/2) = mlon+dlonm*rtd/2.
              epkv(nt,2*latp1+1-j,i+lonmx/2) = epot*1.e-3
             endif
!  Read SH 0-358 mlon
           read (ls_amie,"(17x,f13.4,f18.3,72x,f21.4,8x,i4,5(1x,i2))")
     |      mlat,mlon,epot,iy2,imo2,ida2,ih2,imin2,isec2
!  Save -89 to -60 mlat and -180 to +180 mlon (not 1:2:359 mlon)
	     read_mlat(latp1+1-j) = mlat - dmlat/2.
             if (i>lonmx/2) then
              epkv(nt,latp1+1-j,i-lonmx/2) = epot*1.e-3
             else
              epkv(nt,latp1+1-j,i+lonmx/2) = epot*1.e-3
             endif
           enddo	! do i=1,lonmx
         enddo		! do j=1,latp1-1
        endif		! if (iamie==5) then,else
! Calculate polar values for SD
         polekvs = 0. 
         polekvn = 0.
         do i=1,lonmx
           polekvs = polekvs + epkv(nt,2,i)
           polekvn = polekvn + epkv(nt,latp1+2,i)
         enddo
         polekvs = polekvs/lonmx
         polekvn = polekvn/lonmx
	  read_mlat(1) = -90.
	  read_mlat(latp1+1) = 90.
         do i=1,lonmx
          epkv(nt,1,i) = polekvs
          epkv(nt,latp1+1,i) = polekvn
         enddo
! TEMP find min-max and print out CP and time to fort.52
          epot_nsh = 0.
          epot_xsh = 0.
          epot_nnh = 0.
          epot_xnh = 0.
          do j=2,latp1
           do i=1,lonmx
             if (epot_nsh>epkv(nt,j,i)) then
                  epot_nsh = epkv(nt,j,i)
                  mlat_nsh = lats(i,j)
                  if (iamie==5) then
		    mlt_nsh = mlts(i,j)
                  else
                    mlt_nsh = amod((read_mlt(i)+rot)/15.+24.,24.)
                  endif
             endif
             if (epot_xsh<epkv(nt,j,i)) then
                  epot_xsh = epkv(nt,j,i)
                  mlat_xsh = lats(i,j) 
                  if (iamie==5) then
		    mlt_xsh = mlts(i,j)
                  else
                    mlt_xsh = amod((read_mlt(i)+rot)/15.+24.,24.)
                  endif
             endif
             if (epot_nnh>epkv(nt,j+latp1,i)) then
                  epot_nnh = epkv(nt,j+latp1,i)
                  mlat_nnh = lats(i,j+latp1) 
                  if (iamie==5) then
		    mlt_nnh = mlts(i,j+latp1)
                  else
                    mlt_nnh = amod((read_mlt(i)+rot)/15.+24.,24.)
                  endif
             endif
             if (epot_xnh<epkv(nt,j+latp1,i)) then
                  epot_xnh = epkv(nt,j+latp1,i)
                  mlat_xnh = lats(i,j+latp1) 
                  if (iamie==5) then
		    mlt_xnh = mlts(i,j+latp1)
                  else
                    mlt_xnh = amod((read_mlt(i)+rot)/15.+24.,24.)
                  endif
             endif
!             epot_nsh = min(epkv(nt,j,i),epot_nsh)
!             epot_xsh = max(epkv(nt,j,i),epot_xsh)
!             epot_nnh = min(epkv(nt,j+latp1,i),epot_nnh)
!             epot_xnh = max(epkv(nt,j+latp1,i),epot_xnh)
           enddo
          enddo
          if (nt==2) write (52,"(i1,' iamie Superdarn CP'/
     |     '  nt year nda mnpd cp_sh_kV cp_nh_kV ctpoten',
     |     ' nsh:epot,mlt,mlat xsh:epot,mlt,mlat nnh:epot,mlt,mlat',
     |     ' xnh:epot,mlt,mlat')") iamie
           write (52,"(i4,2i5,16f8.2)") nt,iy,amie_day(nt),
     |     amie_min(nt),epot_xsh-epot_nsh,
     |     epot_xnh-epot_nnh,0.5*(epot_xsh-epot_nsh+epot_xnh-epot_nnh),
     |     epot_nsh,mlt_nsh,mlat_nsh,epot_xsh,mlt_xsh,mlat_xsh,
     |     epot_nnh,mlt_nnh,mlat_nnh,epot_xnh,mlt_xnh,mlat_xnh

!  Check that read_mlt and read_mlat same as mlts (or alon) and degmlat
         if (nt==2) then
          ifstop = 0
          do j=1,latp1*2
           if (abs(lats(1,j)-read_mlat(j)) .gt. 0.1) ifstop=1
          enddo
          if (ifstop==1) then
            write (6,"(1x,'STOP: lats.ne.read_mlat'/(1x,10f8.2))") 
     |      (lats(1,j),read_mlat(j),j=1,latp1*2)
          stop
          endif
          if (iamie==5) then
          do i=2,lonmx-1
           if (abs(mlts(i,1)-read_mlt(i)) .gt. 0.1) ifstop=1
          enddo
          if (ifstop==1) then
            write (6,"(1x,'STOP: mlts.ne.read_mlt',(1x,10f8.2))")
     |       (mlts(i,1),read_mlt(i),i=1,lonmx)
            stop
          endif
          else
          do i=2,lonmx-1
           if (abs(alonm(i)*rtd-read_mlt(i)) .gt. 0.1) ifstop=1
          enddo
          if (ifstop==1) then
            write (6,"(1x,'STOP: alon.ne.read_mlt',(1x,10f8.2))")
     |      (alonm(i)*rtd,read_mlt(i),i=1,lonmx)
            stop
          endif
          endif         ! if (iamie==5) then,else
         endif   	! if (nt==1) then    

! wrap-around loc 1=lonp1 (mlt=0 and 24) is necessary since MLT=24 is not read in.
          do j=1,latp1*2
           epkv(nt,j,lonp1) = epkv(nt,j,1)
          enddo
        enddo		! do nt=2,mxtimes-1
! Add times 1,mxtimes assuming 2min SD chunks so add 1 min
        nt1 = 1
        nt2 = mxtimes
        amie_day(nt1) = amie_day(nt1+1)
        amie_min(nt1) = 0
        amie_day(nt2) = amie_day(nt2-1) + 1
        amie_min(nt2) = 0
        do i=1,lonp1
         do j=1,latp1*2
          epkv(nt1,j,i) = epkv(nt1+1,j,i)
          epkv(nt2,j,i) = epkv(nt2-1,j,i)
         enddo
        enddo
       write (6,"(1x,'Reading SuperDARN iamie from 2 to mxtimes-1 ',
     |   '(daynum,minpd) every delmin minutes =',2i4,2(1x,i4,f9.2), 
     |   i10)") iamie,mxtimes,amie_day(2),amie_min(2),
     |   amie_day(mxtimes-1),amie_min(mxtimes-1),inmin
      endif		! if (iamie==5 .or. iamie==6)

      end subroutine init_amie
!-----------------------------------------------------------------------
      subroutine getamie(iyear,iday,iutsec,iprint,istep,model_hpi)
!
!  04/13 bae:  ASTRA AMIE is every 5 min (interp), U MI AMIE is every 1 min (no interp in time)
!
! Read AMIE outputs from amie_ncfile file, returning electric potential,
! auroral mean energy and energy flux at current date and time,
! and the data is linearly interpolated to the model time
! gl - 12/07/2002
!
      use ModKind			! 6/12:  for AMIE U MI 
      use aurora_module,only: calc_hp
      use hist_module,only: modeltime
      use input_module,only: amiesh,amienh
      use magfield_module,only: im,jm,dim,djm
      use dynamo_module,only: mag2geo,nmlat0,phihm
      use input_module,only:  ! from user input
     |  irdpot,  !  flag (2,1,0) if have outside 2 (epot and aurora gekev+gmwpm2), 
!                   1 (epot only), or 0 (internal W05 or Heelis epot) for the 
!                   electric potential (epot) from the potential_model used
     |  ctpoten,  ! cross-cap potential (kV)      (e.g., 45.)
     |  power     ! hemispheric power   (GW)      (e.g., 16.)
!
! Args:
      integer,intent(in) :: iyear,iday,iutsec,iprint,istep
      real,intent(out) :: model_hpi(2)
!
! Local:
      real :: 
     |   potm(lonmx+1,jmxm),efxm(lonmx+1,jmxm),ekvm(lonmx+1,jmxm)

      integer :: ier,lw,liw,isign,intpol(2)
      integer,allocatable :: iw(:)
      real,allocatable :: w(:)
      integer :: i,j,ndmon_nl(13),ndmon_lp(13),nonleap
      integer :: nn, iset, iset1, m, mp1, n, amie_debug,amieday
      real    :: maxefx
! 6/12: added by Emery for Ridley binary AMIE files
      integer, dimension(7) :: itime
      real (Real8_) :: rtime
!      real :: rtime
      integer :: iError, imlt,ilat,jl
      real :: hpi
		
      integer :: nt,ihem
      real :: frac1,frac2,rmodelminpd,amieminpd,amieminpdm1,delmin	

      integer :: ip,ipp
      real :: polekv,polekev,polemwpm2
      real :: alonmi,addhalf

!   
! Data: GSWM data given at the middle of each month 
!       -> day of year for this data assuming non-leap year
!         J  F  M  A   M   J   J   A   S   O   N   D   J
      data ndmon_nl                                ! non_leap
     |  /31,59,90,120,151,181,212,243,273,304,334,365,396/ 
      data ndmon_lp                                ! leap year
     |  /31,60,91,121,152,182,213,244,274,305,335,366,397/ 
!
!

      if (iprint > 0) then
        write(6,"(/,72('-'))")
        write(6,"('GETAMIE:')")
        write(6,"('Modeltime(1:4,dayn,UTh,m,s)= ',4i5)") modeltime(1:4)
        write(6,"('Initial requested iyear=',i4,' iday=',i3,' iutsec=',
     |    i10)") iyear,iday,iutsec
      endif

      if (istep<=1) then
       call init_amie
      endif	! if (istep<=1) then

      maxefx=50.
! 
! Check for leap year
      nonleap = int(mod(real(iyear),4.))  ! nonleap year /= 0
!     if (iyear.eq.2000) nonleap = 1
! 

! Check times
      if (modeltime(3).lt.0.or.modeltime(3).gt.59) then
        write(6,"('>>> getamie: bad imin=',3i4)") modeltime(3)
        stop 'getamie'
      endif
      if (modeltime(4).lt.0.or.modeltime(4).gt.59) then
        write(6,"('>>> getamie: bad isec=',3i4)") modeltime(4)
        stop 'getamie'
      endif

      if (iamie==2 .or. iamie==5 .or. iamie==6) then
! 04/13 bae: ASTRA AMIE (modeltime=daynum,hr,min,sec)
! 04/14 bae: SuperDARN also stored in arrays like ASTRA AMIE - now use minutes per day 
       do nt=2,mxtimes
        rmodelminpd = modeltime(2)*60 + modeltime(3) + 
     |     modeltime(4)/60.
        amieminpdm1 = amie_min(nt-1)
        amieminpd = amie_min(nt)
! Check if go over a day boundary (will not work for Dec 31 to Jan 01!)
        if (amie_day(nt-1)==amie_day(nt)-1) then
          amieminpd = amieminpd + 1440.
          amieminpdm1 = amieminpdm1 - 1440.
        endif
        delmin = amieminpd - amie_min(nt-1)
        if (amie_day(nt-1)==modeltime(1) .and. amieminpd .ge.
     |    rmodelminpd .and. amie_min(nt-1) .le. rmodelminpd) go to 100
        if (amie_day(nt)==modeltime(1) .and. amieminpdm1 .le.   
     |   rmodelminpd .and. amie_min(nt) .ge. rmodelminpd) go to 101
       enddo
        write (6,"(1x,'Cannot find modeltime in amie_day,min; ',
     |   'modeltime rmodelminpd amie_day,min(1,mxtimes) =',4i5,1x,f9.2,
     |   1x,i4,f9.2,i6,1x,i4,f9.2)") modeltime,rmodelminpd,
     |   amie_day(1),amie_min(1),mxtimes,amie_day(mxtimes),
     |   amie_min(mxtimes)
       stop
!  Use nt and nt-1 for interpolation over delmin or amie_min(nt-1,nt)
100    frac1 = (amieminpd-rmodelminpd) / delmin
       go to 102
101    frac1 = (amie_min(nt)-rmodelminpd) / delmin 
102    frac2 = 1. - frac1
       if (frac1<0. .or. frac1>1.) then
        write (6,"(1x,'bad: frac1,2 nt-1,nt modeltime rmodelminpd ',
     |   'amie_day,min(nt-1,nt) amieminpdm1 amieminpd delmin =',2f7.2,
     |   2i5,1x,4i5,1x,f9.2,2(i5,f9.2),3f9.2)") frac1,frac2,nt-1,nt,
     |   modeltime,rmodelminpd,amie_day(nt-1),amie_min(nt-1),
     |   amie_day(nt),amie_min(nt),amieminpdm1,amieminpd,delmin
        stop
       endif
      if (iprint > 0) then
        write(6,"('Initial frac1,2,nt-1,nt,amie_day,min =',2f7.2,2i5,
     |   2(i5,f9.2))") frac1,frac2,nt-1,nt,amie_day(nt-1),
     |   amie_min(nt-1),amie_day(nt),amie_min(nt)
      endif
!  ASTRA AMIE arrays (kV,mw/m2,keV) are read over 0000to2320 mlts from 90to40, -40to-90 mlat
!    but are stored in init_amie in order SH -90 to ~-40 mlat, NH +90 to ~+40mlat
!  Store in U MI AMIE arrays (V,mw/m2,keV) in order SH -90 to ~-40 mlat, NH +90 to ~+40mlat
        do i=1,lonmx+1
         do j=1,latp1*2
          TempPotential(i,j) = frac1*epkv(nt-1,j,i)*1.e+3 
     |     + frac2*epkv(nt,j,i)*1.e+3
         enddo
        enddo
! Only ASTRA has aurora also.  SuperDARN only has ion convection
       if (iamie==2) then
        do i=1,lonmx+1
         do j=1,latp1*2
          AveEOut(i,j) = frac1*ekev(nt-1,j,i) 
     |     + frac2*ekev(nt,j,i)
          EFluxOut(i,j) = frac1*wpm2(nt-1,j,i) 
     |     + frac2*wpm2(nt,j,i)
         enddo
        enddo
       endif	! if (iamie==2) then
      endif	!       if (iamie==2 .or. iamie==5 .or. iamie==6) then

      if (iamie==3 .or. iamie==4) then
! 6/12: added by Emery for Ridley binary AMIE
!   modeltime(1) = daynumber of year 
      itime(1) = iyear
      do i=1,12
       if (nonleap == 0) then
	if (ndmon_lp(i) <= modeltime(1)) then
          itime(2) = i+1
        endif
       else
	if (ndmon_nl(i) <= modeltime(1)) then
          itime(2) = i+1
        endif
       endif
      enddo
       if (nonleap == 0) then
        itime(3) = modeltime(1)-ndmon_lp(itime(2)-1)
       else
        itime(3) = modeltime(1)-ndmon_nl(itime(2)-1)
       endif
      itime(4) = modeltime(2)
      itime(5) = modeltime(3)
      itime(6) = modeltime(4) 
      itime(7) = 0
! Calling time conversion for Ridley AMIE
!      write(*,*) nonleap,itime
      call time_int_to_real(itime, rtime)
!      write(*,*) "rtime =",rtime
      call IO_SetTime(rtime)
!      write(*,*) "after SetTime rtime =",rtime
! 8/13 bae:  Changed in timing_IO_UMI_AMIE.F Method from IE_Closest_ to IE_Interpolate_ for Eflux and AveE
      call get_AMIE_values(rtime)
!      write(*,*) "after get_AMIE rtime =",rtime
! 
! AveEOut in keV (from ~1-9 keV)
      call IO_GetAveE(AveEOut, iError)
      if (iError /= 0) then
       write(*,*) "Error : ",iError
      else
!      write(*,*) "AveEOut itime=",itime
!      write(*,*) AveEOut
      endif

! EfluxOut in mW/m2? (~24 mW/m2 max 06349 0030UT - has to be mW/m2!)
!  6/12: Emery - U MI AMIE EFluxOut is in erg/cm2-s or mW/m2 (NOT in W/m2!)
      call IO_GetEFlux(EfluxOut, iError)
      if (iError /= 0) then
       write(*,*) "Error : ",iError
      else
!      write(*,*) "EfluxOut itime=",itime
!      write(*,*) EfluxOut
      endif

      call IO_GetPotential(TempPotential, iError) 

      if (iError /= 0) then
       write(*,*) "Error : ",iError
       stop
      else
      endif		! for if read potential w error

!     if (itime(4).eq.0 .and. itime(5).eq.2 .and. itime(3).eq.13) then
!       write(*,*) "i,rtime=",itime,rtime
! mlts = 10x 1-25 (same as 25 IEi_HavenMLTS)
!    write(*,*) "mlts=",mlts
! lats = 25x 86,82,78,74,70,66,62,58,54,50 (10 mlats - diff from 39 IEi_HavenLats - expected 88,87-50)
!    write(*,*) "lats=",lats
! TempPotential in Volts 
!  For SWMF, have zeros at poles and at i=1 for each j
!       write(*,*) "SH pot in V, -90 to -42 mlat"
!       write(*,*) ((TempPotential(i,j),i=1,lonmx),j=1,latp1)
!       write(*,*) "NH pot in V, 90 to 42 mlat"
!       write(*,*) ((TempPotential(i,j),i=1,lonmx),j=latp1+1,latp1*2)
!      endif

      if (iamie==4) then
! 8/13 bae found SWMF 1=181 constant at all mlats from 2-180, so had to be revised as ave of #2+180!
       do j=1,latp1+latp1
        TempPotential(1,j) = 0.5*(TempPotential(2,j) + 
     |    TempPotential(lonmx,j))
        AveEOut(1,j) = 0.5*(AveEOut(2,j)+AveEOut(lonmx,j))
        EfluxOut(1,j) = 0.5*(EfluxOut(2,j)+EfluxOut(lonmx,j))
        TempPotential(lonmx+1,j) = TempPotential(1,j)
        AveEOut(lonmx+1,j) = AveEOut(1,j)
        EfluxOut(lonmx+1,j) = EfluxOut(1,j)
       enddo
! TEMP printout - now looks good (except poles still 0!)
!      do ipp=1,latp1
!       write (6,"(1x,'i=1,2 pot,kev,eflx at j=',i3,6e12.4)") ipp,
!     | TempPotential(1,ipp),AveEOut(1,ipp),EfluxOut(1,ipp),
!     | TempPotential(2,ipp),AveEOut(2,ipp),EfluxOut(2,ipp)
!      enddo

!  Pole values are trash for SWMF, so compute values at j=1 and j=latp1+1
! Find pole value as average of closest latitude for SWMF (and for AMIE also)
         ip = 1
         ipp = 2
         polekv = 0. 
         polekev = 0. 
         polemwpm2 = 0. 
         do i=1,lonmx
           polekv = polekv + TempPotential(i,ipp)
           polekev = polekev + AveEOut(i,ipp)
           polemwpm2 = polemwpm2 + EfluxOut(i,ipp)
         enddo
         polekv = polekv/lonmx
         polekev = polekev/lonmx
         polemwpm2 = polemwpm2/lonmx
         do i=1,lonmx+1
           TempPotential(i,ip) = polekv
           AveEOut(i,ip) = polekev           
           EfluxOut(i,ip) = polemwpm2
         enddo
         ip = latp1+1
         ipp = latp1+2
         polekv = 0. 
         polekev = 0. 
         polemwpm2 = 0. 
         do i=1,lonmx
           polekv = polekv + TempPotential(i,ipp)
           polekev = polekev + AveEOut(i,ipp)
           polemwpm2 = polemwpm2 + EfluxOut(i,ipp)
         enddo
         polekv = polekv/lonmx
         polekev = polekev/lonmx
         polemwpm2 = polemwpm2/lonmx
         do i=1,lonmx+1
           TempPotential(i,ip) = polekv
           AveEOut(i,ip) = polekev           
           EfluxOut(i,ip) = polemwpm2
         enddo
	if (istep<3) write (6,"(1x,'polekv,kev,mwpm2 =',3e12.4)") 
     |   polekv,polekev,polemwpm2
      endif	!       if (iamie==4) then
      endif	!       if (iamie==3 .or. iamie==4) then

! get rid of outragous flux values > maxefx

      if (irdpot==2) then
!  Calculate hpi_s,nh_amie in GW
       call calc_hp (latp1*2,lonmx,lonmx+1,degmlat,EfluxOut,
     |   model_hpi(1),model_hpi(2)) 
      else
       model_hpi(1) = power
       model_hpi(2) = power
      endif
!
! gl - 14/07/2002 changed from Emery's original amieterp.f
!
!  Change mag coords so latitude extends to the equator in into the opposite
!   pole, and so MLT is converted to apex longitues between -180 and +180.
!  Put EPOTDUP, EMKEV and WPM2 into S. Hem and N. Hem alike in EPOTM, EKEVM,
!   and EFLXM, dimensioned (IMXMP,JMNH) where MLAT goes from -90 to ~-2
!   for ihsoln=1, and from +90 to ~+2 for ihsoln=2 since GRDSET is for
!   -90 to the equator.  (Will work for +90 to the equator also.)
!   AMIE grid goes from -90 to ~-40 for SH and +90 to ~+40 for NH (6/12: w dmlat=1.67deg))
!   TGCM grid goes from -87.5 to -2.5 for SH and +2.5 to +87.5 for NH (6/12: w dglat=5deg)
!   Hence, at very end, need to reverse order of NH in whole globe arrays.

      if (iamie .ne. 6) then
!  All models but VT SD in MLT start with 0 MLT (NOTE Dartmouth SD starts at (i-1+0.5)*dmltm
! mlongitude starts from -180 degree
! rot=15*ihr+min/4 - 70 is mlon at local midnight for mlt=(alon+rot)/15
! rot=250-iutsec*15/3600 is mlon at local noon for mlt(0 start)=alon(-180 start)-rot
      rot = 250./15. - iutsec/3600.
      addhalf = 0.
      if (iamie==5) addhalf = 0.5
      do i=1,lonp1
        xmlt = (i-1+addhalf)*dmltm - rot + 24.
        xmlt = AMOD(xmlt,24.)
        m = IFIX(xmlt/dmltm + 1.01)
        mp1 = m + 1
        IF (mp1 > lonp1) mp1 = 2
        del = xmlt - (m-1)*dmltm

C  Initialize arrays around equator
        do j=latp1+1,ithmx
          potm(i,j) = 0.
          potm(i,jmxm+1-j) = 0.
         if (irdpot==2) then
          ekvm(i,j) = (1.-del)*AveEOut(m,latp1) + 
     |                        del*AveEOut(mp1,latp1)
          ekvm(i,jmxm+1-j) = (1.-del)*AveEOut(m,latp1+latp1) +
     |                        del*AveEOut(mp1,latp1+latp1)
          efxm(i,j) = 0.
          efxm(i,jmxm+1-j) = 0.
         endif
        enddo
C  Put in AMIE arrays from pole to latp1
        do j=1,latp1
          potm(i,j) = (1.-del)*TempPotential(m,j) + 
     |                       del*TempPotential(mp1,j)
          potm(i,jmxm+1-j) = (1.-del)*TempPotential(m,j+latp1) + 
     |                       del*TempPotential(mp1,j+latp1)
         if (irdpot==2) then
          ekvm(i,j) = (1.-del)*AveEOut(m,j) + del*AveEOut(mp1,j)
          ekvm(i,jmxm+1-j) = (1.-del)*AveEOut(m,j+latp1) + 
     |                       del*AveEOut(mp1,j+latp1)
          efxm(i,j) = (1.-del)*EfluxOut(m,j) + del*EfluxOut(mp1,j)
          efxm(i,jmxm+1-j) = (1.-del)*EfluxOut(m,j+latp1) + 
     |                       del*EfluxOut(mp1,j+latp1)
         endif
        enddo
     
      enddo		!       do i=1,lonp1

!  VT_SD stores in alonm -180 to +180 already
       else
        do i=1,lonp1
         do j=1,latp1
          potm(i,j) = TempPotential(i,j)
          potm(i,jmxm+1-j) = TempPotential(i,j+latp1) 
         enddo
C  Initialize arrays around equator
         do j=latp1+1,ithmx
          potm(i,j) = 0.
          potm(i,jmxm+1-j) = 0.
         enddo
        enddo
       endif		! if (iamie .ne. 6) then,else

C  Set up coeffs to go between POTM(IMXMP,JMNH) and phihm(IMAXM,JMAXMH)in dynamo

! ylatm and ylonm are arrays of latitudes and longitudes of the
! distored magnetic grids in radian - from consdyn in cons.F
!  Convert from amie magnetic grid to distorted magnetic grid
!
! Allocate workspace for regrid routine rgrd2.f:
      lw = nmlonp1+nmlat+2*nmlonp1
      if (.not. allocated(w)) allocate(w(lw),stat=ier)
      if (ier /= 0) write(6,"('>>> horizontal_interp: error allocating',
     |  ' w(lw): lw=',i6,' ier=',i4)") lw,ier
      liw = nmlonp1 + nmlat
      if (.not. allocated(iw)) allocate(iw(liw),stat=ier)
      if (ier /= 0) write(6,"('>>> horzontal_interp: error allocating',
     |  ' iw(liw): liw=',i6,' ier=',i4)") liw,ier
      intpol(:) = 1 ! linear (not cubic) interp in both dimensions

      call rgrd2(lonp1,jmxm,alonm,alatm,potm,nmlonp1,nmlat,
     |           ylonm,ylatm,phihm,intpol,w,lw,iw,liw,ier)

      if (irdpot==2) then
      call rgrd2(lonp1,jmxm,alonm,alatm,ekvm,nmlonp1,nmlat,
     |           ylonm,ylatm,tieekv,intpol,w,lw,iw,liw,ier)
      call rgrd2(lonp1,jmxm,alonm,alatm,efxm,nmlonp1,nmlat,
     |           ylonm,ylatm,tieefx,intpol,w,lw,iw,liw,ier)

C  Convert from TGCM distorted magnetic grid to geographic one
       call mag2geo(tieekv(1,1),gekev(1,0),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
       call mag2geo(tieefx(1,1),gmwpm2(1,0),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
C     ****
C     ****     INSERT PERIODIC POINTS
C     ****
C 4/13 bae: use i=1,nlon of gekev and gmwpm2 in aurora.F and below, BUT nlonp1 is needed for alfa,eflux!
      DO j = 1,nlat
        gekev(nlonp1,j) = gekev(1,j)
        gmwpm2(nlonp1,j) = gmwpm2(1,j)
      ENDDO
      endif	! if (irdpot==2)

! temp - save the AMIE output in magnetic grid
      amie_debug = 1
      if (amie_debug > 0) then
        call addfld('PHIHM2D','2D AMIE-TYPE ELECTRIC POTENTIAL','V',
     |    phihm,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
       if (irdpot==2) then
        call addfld('TIEEKV',' ',' ',tieekv,'mlon',1,nmlonp1,
     |    'mlat',1,nmlat,0)
        call addfld('TIEEFX',' ',' ',tieefx,'mlon',1,nmlonp1,
     |    'mlat',1,nmlat,0)
       endif
      endif
!      
      if (iprint > 0) then
        write(6,"('getamie: AMIE data interpolated to date and time')")
        write(6,"(72('-'),/)")
      endif
      end subroutine getamie    
!-------------------------------------------------------------------

      end module amie_module

