      subroutine kameleon_tiegcm (istep,model_hpi)

c  Ja Soon Shim:  This is an example program mainly to test the time interpolator
c  04/13 bae:  Revised to read Dec 13-15, 2006 SWMF 5 min files of epot,eflux,eave

c Variables to be used for interpolation and data extraction
c  ep = electric potential in kV
c  eave = electron average energy in keV
c  eflux = electron energy flux in W/m2
c  ex,y,z = electric field in x(phi??), y(theta??), z(up??) directions in mV/m

      use input_module,only:  ! from user input
     |  ctpoten  ! cross-cap potential (kV)      (e.g., 45.)
!!   |  ,power   ! hemispheric power   (GW)      (e.g., 16.)
      use input_module,only: kameleon_path,facmw,lbmlat
      use params_module,only: nlonp1,nmlonp1,nlat,nmlon,nmlat,
     |    nlatp1,nlon
      use cons_module,only: pi,dtr,rtd,ylonm,ylatm
      use dynamo_module,only: nmlat0,phihm
      use magfield_module,only: sunlons 
      use magfield_module,only: im,jm,dim,djm
      use init_module,only: iyear,iday,uthr
      use hist_module,only: modeltime
      use aurora_module,only: theta0,offc,dskofc,calc_hp
      use aurora_module,only: gekev,gmwpm2	! (nlonp1,0:nlatp1) in geog coords
      use dynamo_module,only: mag2geo
      use addfld_module,only: addfld
      implicit none		

!
! Args:
      integer,intent(in) :: istep
      real,intent(out) :: model_hpi(2)
!
! Local:
      integer :: iyr,nda,imon,ida,ihh,imm
      integer :: iyear1,iyear2,imon1,imon2,iday1,iday2,ihh1,ihh2,
     |           imm1,imm2
      integer :: j,i,iv,ih,ip,ipp,j1,j2
      real :: rmlt(nmlonp1),mlat,vnx(2,2)
      real :: xshmwpm2,xnhmwpm2,polekv,polekev,polemwpm2
! Define auroral output fields in TIEGCM distorted magnetic grid (same as phihm)
      real ::  tieefx(nmlonp1,nmlat),tieekv(nmlonp1,nmlat)
!      integer,save :: iprint3 
!      real,save :: da(nmlat),Re,areatot,dmlat,dlonm
      real :: rhon,rhox,thetandeg,thetaxdeg,xn,yn,xx,yx,xcen,ycen
      real :: dskofcen(2),offcen(2),radcen(2)

      character*300 cdf_file_path1
      character*300 cdf_file_path2
      character*24 time
       
      real :: time_epoch
      integer :: tid
      integer :: iskp2

! External functions 
      external f_timeinterp_create
      external f_timeinterp_addtimestep
      external f_timeinterp_managemem,f_timeinterp_interpolate
      external f_timeinterp_timetoepoch
      external f_timeinterp_timestrtoepoch
      external f_timeinterp_closeandcleanupmemory

! Find da for HPI estimation
      if (istep<=1) then
! Set up printout for istep<3 for polar values (messed up if nproc>1)
!       iprint3 = 1
! Re=radius of the earth in m
!      Re = 6371.e+3
! ylatm (97 values) and ylonm (80 values 4.5 deg apart) are arrays of latitudes 
!   and longitudes of the distorted tiegcm magnetic grid in radians - from consdyn.h
!  nmlat0 = 49 for mag equator, nmlat = 97, ylatm from -/-90 to +/-0.94, with dellat>3deg 35-48mlat
!   #77-97 or 21 values for 41.9mlat to 90mlat;  Use kameleon only for mlat>=lbmlat
!      dlonm = 2.*pi/float(nmlon)
!      areatot = 0.
!      write (6,"(1x,'nmlon,dlonm =',i4,f7.2)") nmlon,dlonm/dtr
! Set da for NH first from eq to pole (so set j-1 to abs for SH at j=nmlat0)
!      do j=nmlat0,nmlat
!       if (j==nmlat) then
!        dmlat = 90.*dtr - ylatm(j) + 0.5*(ylatm(j)-ylatm(j-1))
!       else
!        dmlat = 0.5*(ylatm(j+1)-abs(ylatm(j-1)))
!       endif
!       da(j) = Re*cos(ylatm(j))*dlonm*(Re*dmlat)
!       write (6,"(1x,'j dmlat mlat da =',i3,3e12.4)") j,dmlat/dtr,
!     |  ylatm(j)/dtr,da(j)
!       areatot = areatot + da(j)*nmlon
!      enddo
! Set SH from NH da
!      do j=1,nmlat0-1
!       da(j) = da(nmlat+1-j)
!      enddo
!       write (6,"(1x,'istep<=1 area of hem = ',e12.4)") areatot
      endif	! if (istep<=1) then

! Redefine RMLT with new sunlons(1-nmlonp1) all the same but change each time step
        do i=1,nmlonp1  
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc) in rad
          rmlt(i) = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
          if (rmlt(i)>24.) rmlt(i) = rmlt(i)-24.
          if (rmlt(i)<0.) rmlt(i) = rmlt(i)+24.
        enddo

! Convert from day-of-year to month,day to get names for SWMF 5-min files
!  sub cvt2md in wei05sc.F
!
      iyr   = iyear
      nda = iday
      call cvt2md(6,iyr,nda,imon,ida)
      ihh = int(uthr)
      imm = int((uthr - ihh)*60.)
      write(time,'(i4.4,2("-",i2.2),"T",2(i2.2,":"),"00.000Z")')   
     &   iyear,imon,ida,ihh,imm 
      if (istep<=1) write (6,"(1x,'time = ',a)") time
C Find times to nearest 5 minutes
       imm1=imm/5*5
       imm2=imm1+5
       if(imm2.eq.60) then
          imm2=0         
          ihh1=ihh
          ihh2=ihh1+1
        else
          ihh1=ihh
          ihh2=ihh1 
       endif 
       if(ihh2.eq.24) then
          ihh2=0         
          iday1=ida
          iday2=iday1+1
        else
          iday1=ida
          iday2=iday1
        endif 
         
         imon1=imon
         imon2=imon1  
         iyear1=iyear
         iyear2=iyear1

C Have no 0005UT 16 Dec, so have special case
       iskp2 = 0
       if (iday2.eq.16 .and.ihh.eq.0 .and.imm.lt.5) iskp2 = 1

! Write 2 cdf file paths
      write(cdf_file_path1,'(a,i4.4,2i2.2,"-",2i2.2,"00-000.cdf")')   
     &   trim(kameleon_path),iyear1,imon1,iday1,ihh1,imm1 
      write(cdf_file_path2,'(a,i4.4,2i2.2,"-",2i2.2,"00-000.cdf")')   
     &   trim(kameleon_path),iyear2,imon2,iday2,ihh2,imm2 
      if (istep<=1) write (6,"(1x,'path1,2 = ',a,1x,a)") 
     |   trim(cdf_file_path1),trim(cdf_file_path2)

c     tid will be set when calling f_timeinterp_create. Just initialize here.
      tid=0
      time_epoch = 0.

      call f_timeinterp_create(tid)
      call f_timeinterp_timestrtoepoch(time, time_epoch)
      if (istep<3) write (6,"(1x,'tid time_epoch =',i5,f18.1)")
     |  tid,time_epoch

      call f_timeinterp_addtimestep(tid, cdf_file_path1)
      if (iskp2==0) call f_timeinterp_addtimestep(tid, cdf_file_path2)

C Find 'ep'
C  NOTE in TimeInterpolator.cpp, could put 3 variables in as:
C   manageMemory(tid, time_epoch, 'ep, 'eave', 'eflux', 3)
C  But fortran wrapper (managemem) has only 1 variable passed in at a time
      call f_timeinterp_managemem(tid, time_epoch, 'ep')
! Get SH (j=1,nmlat0) and NH (j=nmlat0+1,nmlat)
      j1 = 1
      j2 = nmlat
!  Pole values are trash for SWMF, so skip j=1 and nmlat and compute later
!      do j=1,nmlat
!      do j=2,nmlat-1
      do j=j1,j2
       mlat = ylatm(j)*rtd
       ih = 1
       if (j>nmlat0) ih = 2
        do i=1,nmlonp1  
         if (abs(mlat)<lbmlat) then
          phihm(i,j) = 0.
         else
! Obtain electric potential in kV and convert to V
          call f_timeinterp_interpolate(tid, time_epoch, 'ep', 0.0,
     |     mlat, rmlt(i), phihm(i,j))
          phihm(i,j) = phihm(i,j)*1000.
         endif	!  if (abs(mlat)<lbmlat) then,else
        enddo ! i=1,nmlonp1
      enddo ! j=2,nmlat-1

C Find 'eave'
      call f_timeinterp_managemem(tid, time_epoch, 'eave')
!      do j=2,nmlat-1
      do j=j1,j2
       mlat = ylatm(j)*rtd
        do i=1,nmlonp1  
         if (abs(mlat)<lbmlat) then
          tieekv(i,j) = 0.5
         else
! Obtain electron mean energy in keV
         call f_timeinterp_interpolate(tid, time_epoch, 'eave', 0.0,
     |    mlat, rmlt(i), tieekv(i,j)) 
         endif	!  if (abs(mlat)<lbmlat) then,else
        enddo ! i=1,nmlonp1
      enddo ! j=2,nmlat-1

C Find 'eflux' and convert to mW/m2 if necessary
      call f_timeinterp_managemem(tid, time_epoch, 'eflux')
!      do j=2,nmlat-1
      do j=j1,j2
       mlat = ylatm(j)*rtd
        do i=1,nmlonp1  
         if (abs(mlat)<lbmlat) then
          tieefx(i,j) = 0.
         else
! Obtain electron energy flux in w/m2 (SWMF) or in mw/m2 (ASTRA_AMIE)    
          call f_timeinterp_interpolate(tid, time_epoch, 'eflux', 0.0,
     |     mlat, rmlt(i), tieefx(i,j))
          tieefx(i,j) = tieefx(i,j)*facmw 
         endif	!  if (abs(mlat)<lbmlat) then,else
        enddo ! i=1,nmlonp1
      enddo ! j=2,nmlat-1

! TEMP
! printout is messed up if nproc>1
	if (istep<3) then
         write (6,"(1x,'istep nmlon rmlt=',2i3/(10f7.1))") 
     |    istep,nmlon,rmlt
         write (6,"(1x,'nmlat mlat=',i3/(10f7.1))") nmlat,
     |    ((ylatm(i)*rtd),i=1,nmlat)
         ip = 1
         write (6,"(1x,'ip,mlat phihm =',i3,f6.0/(10e12.4))") ip,
     |     ylatm(ip)*rtd,(phihm(i,ip),i=1,nmlon)  
         write (6,"(1x,'ip tiekv =',i3/(10e12.4))") ip, 
     |     (tieekv(i,ip),i=1,nmlon)
         write (6,"(1x,'ip tieefx =',i3/(10e12.4))") ip,
     |     (tieefx(i,ip),i=1,nmlon)
         ip = nmlat
         write (6,"(1x,'ip,mlat phihm =',i3,f6.0/(10e12.4))") ip,
     |     ylatm(ip)*rtd,(phihm(i,ip),i=1,nmlon)  
         write (6,"(1x,'ip tiekv =',i3/(10e12.4))") ip, 
     |     (tieekv(i,ip),i=1,nmlon)
         write (6,"(1x,'ip tieefx =',i3/(10e12.4))") ip,
     |     (tieefx(i,ip),i=1,nmlon)
       endif
!	if (j1.gt. 1 .and. j2.lt.nmlat) then 
! Find pole value as average of closest latitude for SWMF (and for AMIE also)
!  March 2014 B. Emery: kam2 package has bad keV and mW/m2 at poles at 24MLT+/-1to3h; 
!  SH epot varies, NH epot const, but still died after 1.2h only keeping NH epot at pole.  
!  THEREFORE, MUST REPLACE ALL POLAR VALUES! - still bad for kam3 _24mar14 libs
         ip = 1
         ipp = 2
         polekv = 0. 
         polekev = 0. 
         polemwpm2 = 0. 
         do i=1,nmlon
           polekv = polekv + phihm(i,ipp)
           polekev = polekev + tieekv(i,ipp)
           polemwpm2 = polemwpm2 + tieefx(i,ipp)
         enddo
         polekv = polekv/nmlon
         polekev = polekev/nmlon
         polemwpm2 = polemwpm2/nmlon
         do i=1,nmlonp1
           phihm(i,ip) = polekv
           tieekv(i,ip) = polekev           
           tieefx(i,ip) = polemwpm2
         enddo
         ip = nmlat
         ipp = nmlat-1
         polekv = 0. 
         polekev = 0. 
         polemwpm2 = 0. 
         do i=1,nmlon
           polekv = polekv + phihm(i,ipp)
           polekev = polekev + tieekv(i,ipp)
           polemwpm2 = polemwpm2 + tieefx(i,ipp)
         enddo
         polekv = polekv/nmlon
         polekev = polekev/nmlon
         polemwpm2 = polemwpm2/nmlon
         do i=1,nmlonp1
           phihm(i,ip) = polekv
           tieekv(i,ip) = polekev           
           tieefx(i,ip) = polemwpm2
         enddo
! TEMP
	if (istep<3) write (6,"(1x,
     |   'polekv,kev,mwpm2 =',3e12.4)") polekv,polekev,polemwpm2
!       endif

!  Calculate hem HP in GW in distorted magnetic grid
       call calc_hp (nmlat,nmlon,nmlonp1,ylatm/dtr,tieefx,model_hpi(1),
     |   model_hpi(2))
 
C  Convert aurora from TGCM distorted magnetic grid to geographic one for aurora.F
       call mag2geo(tieekv(1,1),gekev(1,0),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
       call mag2geo(tieefx(1,1),gmwpm2(1,0),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
C     ****
C     ****     INSERT PERIODIC POINTS
C     ****
C 4/13 bae: use i=1,nlon of gekev and gmwpm2  in aurora.F and below, BUT nlonp1 is needed for alfa,eflux!
      DO j = 1,nlat
        gekev(nlonp1,j) = gekev(1,j)
        gmwpm2(nlonp1,j) = gmwpm2(1,j)
      ENDDO
! Could calculate theta0 (cusplat and colath) and set MLT for cusp also from calccloc
!  or could use defaults in aurora.F or could revise as in amie.F
! 3/14 bae:   Make offc,dskofc const in revloc, and set default theta0 for both ifbad=1 in input.F from iamie

      call f_timeinterp_closeandcleanupmemory(tid)

        call addfld('PHIHM2D','2D KAMELEON-TYPE ELECTRIC POTENTIAL',
     |    'V',phihm,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
        call addfld('TIEEKV',' ',' ',tieekv,'mlon',1,nmlonp1,
     |    'mlat',1,nmlat,0)
        call addfld('TIEEFX',' ',' ',tieefx,'mlon',1,nmlonp1,
     |    'mlat',1,nmlat,0)
        call addfld('GEKEV',' ',' ',gekev,'lon',1,nlonp1,
     |    'lat',0,nlatp1,0)
        call addfld('GMWPM2',' ',' ',gmwpm2,'lon',1,nlonp1,
     |    'lat',0,nlatp1,0)

      end subroutine
