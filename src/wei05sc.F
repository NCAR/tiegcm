!-----------------------------------------------------------------------
      module wei05sc
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! The Weimer model of high-latitude potential is the intellectual property of 
! Daniel Weimer and may not be extracted, distributed, or used for any purpose 
! other than as implemented in the TIE-GCM.  For further information concerning 
! this model, please contact Dan Weimer.
!
! 2005 Version of the electric and magnetic potential (FAC) models
! by Dan Weimer.  Uses Spherical Cap Harmonic Analysis (SCHA) functions.
! Model description is in:
! Weimer, D. R., Predicting Surface Geomagnetic Variations Using Ionospheric
!    Electrodynamic Models, J. Geophys. Res., 110, A12307, doi:10.1029/
!    2005JA011270, 2005.
! Some information about the model (such as outer boundary calculation)
! is also in the earlier paper:
! Weimer, D. R. (2005), Improved ionospheric electrodynamic models and
!    application to calculating Joule heating rates,  J. Geophys. Res., 110,
!    A05306, doi:10.1029/2004JA010884.
!
! For information about the SCHA, see the paper:
! Haines, G. V., Spherical cap harmonic analysis, J. Geophys. Res., 90, B3,
!  2583, 1985.  (Note that this is in JGR-B, "Solid Earth", rather than JGR-A)
!
! April, 2008:
! This f90 module of the Electric Potential model was translated
!   from the original IDL by Ben Foster (NCAR, foster@ucar.edu)
! Netcdf data file wei05sc.nc was written from original IDL save files
!   W05scBndy.xdr, W05scEpot.xdr, W05scBpot.xdr, and SCHAtable.dat
!
! May 10-12, 2008 bae and btf: 
!   Correct angle by/bz for different hemispheres (change sign)
!   (see bt, angle, and setmodel and setboundary calls)
!
      use params_module,only: nmlat,nmlon,nmlonp1
      use mpi_module,only: mytid ! for stdout print 
      use hist_module,only: modeltime
      implicit none
!
! Coefficients read from netcdf data file wei05sc.nc:
!
      integer,parameter :: 
     |  na=6, nb=7, nex=2, n1_scha=19, n2_scha=7, n3_scha=68,
     |  csize=28, n_schfits=15, n_alschfits=18
      integer :: maxk_scha, maxm_scha, maxl_pot, maxm_pot
      real :: bndya(na), bndyb(nb), ex_bndy(nex), ex_epot(nex),
     |  ex_bpot(nex)
      real :: th0s(n3_scha), allnkm(n1_scha,n2_scha,n3_scha)
      integer :: ab(csize), ls(csize), ms(csize)
      real :: epot_alschfits(n_alschfits,csize),
     |        bpot_alschfits(n_alschfits,csize)
      real :: bpot_schfits(n_schfits,csize), 
     |        epot_schfits(n_schfits,csize)
!
! Intermediate calculations:
!
      real :: rad2deg,deg2rad           ! set by setmodel
      real :: bndyfitr                  ! calculated by setboundary
      real :: esphc(csize),bsphc(csize) ! calculated by setmodel
      real :: tmat(3,3),ttmat(3,3)      ! from setboundary
      integer,parameter :: mxtablesize=500
      real :: plmtable(mxtablesize,csize),colattable(mxtablesize)
      real :: nlms(csize)
      real :: wei05sc_fac(nmlonp1,nmlat)  ! field-aligned current output

! 05/08 bae:  Have ctpoten from both hemispheres from Weimer
      real :: weictpoten(2)

      contains
!-----------------------------------------------------------------------
      subroutine weimer05(by,bz,swvel,swden,wei05sc_ncfile,istep)
!
! 4/2/08 btf:
! Driver to call Weimer 2005 model for tiegcm.
! 7/2/08 bae: Weimer setboundary 2005 constant location and bndyfitr/2 radius to 0 poten
! 07/08:  Added ctpoten to store average weictpoten from wei05loc
!
      use input_module,only: ctpoten
      use cons_module,only: rtd,ylonm,ylatm
      use pdynamo_module,only: nmlat0,phihm
      use magfield_module,only: sunlons 
      use init_module,only: iyear,iday,uthr
      implicit none
!
! Args:
      real,intent(in) :: by,swvel,swden
      real,intent(inout) :: bz
      character(len=*),intent(in) :: wei05sc_ncfile
      integer,intent(in) :: istep
!
! Local:
      real :: angl,angle,bt
      integer :: i,j
      real :: rmlt,mlat,tilt,htilt,hem,ut
      real,parameter :: fill=0.
      integer :: iyr,nda,imo,ida
!
! At least one of by,bz must be non-zero: 
      if (by==0..and.bz==0.) then
        write(6,"(/,'>>> WARNING: by and bz cannot both be zero',
     |    ' when calling the Weimer model: am setting bz=0.01')")
        bz = 0.01 ! note intent(inout) above
      endif
!
      bt = sqrt(by**2+bz**2)
      angl = atan2(by,bz)*rtd
!
! Read netcdf data file if this is first step:
      if (istep==1) call read_wei05sc_ncfile(wei05sc_ncfile)
!
! Convert from day-of-year to month,day and get tilt from date and ut:
! sub cvt2md and function get_tilt below were copied from the 
!   Weimer 2001 code.
!
      iyr   = iyear
      nda = iday
      ut = uthr
      call cvt2md(6,iyr,nda,imo,ida)
      tilt = get_tilt(iyr,imo,ida,ut)
!
! Call Weimer model for southern hemisphere epot:
      hem = -1.
      htilt = hem * tilt
      angle = hem * angl
      call setmodel(angle,bt,htilt,swvel,swden,'epot')
      do j=1,nmlat0
        do i=1,nmlon
!
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc) in rad
          rmlt = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
          mlat = abs(ylatm(j))*rtd
!
! Obtain electric potential and convert from kV to V
          call epotval(mlat,rmlt,fill,phihm(i,j))
          phihm(i,j) = phihm(i,j)*1000.
        enddo ! i=1,nmlon
      enddo ! j=1,nmlat0
!
!  Re-calculate SH values of offa, dskofa, arad, and Heelis phid and phin from
!   Weimer 2005 setboundary values of offc, dskofc, and theta0
	call wei05loc (1)
!
! Call Weimer model for southern hemisphere fac:
      call setmodel(angle,bt,htilt,swvel,swden,'bpot')
      do j=1,nmlat0
        do i=1,nmlon
          rmlt = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
          mlat = abs(ylatm(j))*rtd
          call mpfac(mlat,rmlt,fill,wei05sc_fac(i,j))
        enddo ! i=1,nmlon
      enddo ! j=1,nmlat0
!
! Call Weimer model for northern hemisphere epot:
      hem = 1.
      htilt = hem * tilt
      angle = hem * angl
      call setmodel(angle,bt,htilt,swvel,swden,'epot')
      do j=nmlat0+1,nmlat
        do i=1,nmlon
!
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc) in rad
          rmlt = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
          mlat = abs(ylatm(j))*rtd
!
! Obtain electric potential and convert from kV to V
          call epotval(mlat,rmlt,fill,phihm(i,j))
          phihm(i,j) = phihm(i,j)*1000.
        enddo ! i=1,nmlon
      enddo ! j=1,nmlat0+1,nmlat
!
!  Re-calculate NH values of offa, dskofa, arad, and Heelis phid and phin from
!   Weimer 2005 setboundary values of offc, dskofc, and theta0
	 call wei05loc (2)
!
! Call Weimer model for northern hemisphere fac:
      call setmodel(angle,bt,htilt,swvel,swden,'bpot')
      do j=nmlat0+1,nmlat
        do i=1,nmlon
          rmlt = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
          mlat = abs(ylatm(j))*rtd
          call mpfac(mlat,rmlt,fill,wei05sc_fac(i,j))
        enddo ! i=1,nmlon
      enddo ! j=1,nmlat0
!
! Periodic points:
      do j=1,nmlat
        phihm(nmlonp1,j) = phihm(1,j)
        wei05sc_fac(nmlonp1,j) = wei05sc_fac(1,j)
      enddo ! j=1,nmlat
!
! 07/08 bae:  Return the average of the weictpoten from the SH and NH in ctpoten
      ctpoten = 0.5*(weictpoten(1)+weictpoten(2))
!      write (6,"(1x,'weimer05: CPS,N ctpoten=',3f8.1)") weictpoten,
!     |    ctpoten
!      write(6,"(i4,':weimer05: istep=',i8,' mtime=',3i4,' ctpoten=',
!    |  f8.3)") mytid,istep,modeltime(1:3),ctpoten
!
! Note calculation of pfrac for dynamo is done by sub colath (colath.F),
! called from advance.
!
      end subroutine weimer05
!-----------------------------------------------------------------------
      subroutine read_wei05sc_ncfile(file)
!
! Read coefficients and other data from netcdf data file.
!
      use nchist_module,only: handle_ncerr
      implicit none
#include <netcdf.inc>
!
! Arg:
      character(len=*),intent(in) :: file
!
! Local:
      integer :: istat,ncid
      character(len=120) :: char120
      integer :: rd_na,rd_nb,rd_nex,rd_n1_scha,rd_n2_scha,rd_n3_scha,
     |  rd_csize,rd_n_schfits,rd_n_alschfits
      integer :: id
!
! Open netcdf file for reading:
      istat = nf_open(file,NF_NOWRITE,ncid)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_open of netcdf file ',a)") 
     |    trim(file)
        call handle_ncerr(istat,char120)
      else
!       write(6,"('wei05sc: opened netcdf data file',a)") trim(file)
      endif
!
! Read and check dimensions:
!
! na=6
      istat = nf_inq_dimid(ncid,"na",id) 
      istat = nf_inq_dimlen(ncid,id,rd_na)
      if (rd_na /= na) then
        write(6,"(/,'>>> wei05sc: rd_na /= na: rd_na=',i4,' na=',i4)")
     |    rd_na,na
        call shutdown('wei05sc')
      endif
!
! nb=7
!
      istat = nf_inq_dimid(ncid,"nb",id) 
      istat = nf_inq_dimlen(ncid,id,rd_nb)
      if (rd_na /= na) then
        write(6,"(/,'>>> wei05sc: rd_nb /= nb: rd_nb=',i4,' nb=',i4)")
     |    rd_nb,nb
        call shutdown('wei05sc')
      endif
!
! nex=2
!
      istat = nf_inq_dimid(ncid,"nex",id) 
      istat = nf_inq_dimlen(ncid,id,rd_nex)
      if (rd_nex /= nex) then
        write(6,"(/,'>>> wei05sc: rd_nex /= nex: rd_nex=',i4,' nex=',
     |    i4)") rd_nex,nex
        call shutdown('wei05sc')
      endif
!
! n1_scha=19
!
      istat = nf_inq_dimid(ncid,"n1_scha",id) 
      istat = nf_inq_dimlen(ncid,id,rd_n1_scha)
      if (rd_n1_scha /= n1_scha) then
        write(6,"(/,'>>> wei05sc: rd_n1_scha /= n1_scha: rd_n1_scha=',
     |    i4,' n1_scha=',i4)") rd_n1_scha,n1_scha
        call shutdown('wei05sc')
      endif
!
! n2_scha=7
!
      istat = nf_inq_dimid(ncid,"n2_scha",id) 
      istat = nf_inq_dimlen(ncid,id,rd_n2_scha)
      if (rd_n2_scha /= n2_scha) then
        write(6,"(/,'>>> wei05sc: rd_n2_scha /= n2_scha: rd_n2_scha=',
     |    i4,' n2_scha=',i4)") rd_n2_scha,n2_scha
        call shutdown('wei05sc')
      endif
!
! n3_scha=68
!
      istat = nf_inq_dimid(ncid,"n3_scha",id) 
      istat = nf_inq_dimlen(ncid,id,rd_n3_scha)
      if (rd_n3_scha /= n3_scha) then
        write(6,"(/,'>>> wei05sc: rd_n3_scha /= n3_scha: rd_n3_scha=',
     |    i4,' n3_scha=',i4)") rd_n3_scha,n3_scha
        call shutdown('wei05sc')
      endif
!
! csize=28
!
      istat = nf_inq_dimid(ncid,"csize",id) 
      istat = nf_inq_dimlen(ncid,id,rd_csize)
      if (rd_csize /= csize) then
        write(6,"(/,'>>> wei05sc: rd_csize /= csize: rd_csize=',i4,
     |    ' csize=',i4)") rd_csize,csize
        call shutdown('wei05sc')
      endif
!
! n_schfits=15
!
      istat = nf_inq_dimid(ncid,"n_schfits",id) 
      istat = nf_inq_dimlen(ncid,id,rd_n_schfits)
      if (rd_n_schfits /= n_schfits) then
        write(6,"(/,'>>> wei05sc: rd_n_schfits /= n_schfits: ',
     |    'rd_n_schfits=',i4,' n_schfits=',i4)") rd_n_schfits,n_schfits
        call shutdown('wei05sc')
      endif
!
! n_alschfits=18
!
      istat = nf_inq_dimid(ncid,"n_alschfits",id) 
      istat = nf_inq_dimlen(ncid,id,rd_n_alschfits)
      if (rd_n_alschfits /= n_alschfits) then
        write(6,"(/,'>>> wei05sc: rd_n_alschfits /= n_alschfits: ',
     |    'rd_n_alschfits=',i4,' n_alschfits=',i4)") 
     |    rd_n_alschfits,n_alschfits
        call shutdown('wei05sc')
      endif
!
! integer :: maxk_scha, maxm_scha, maxl_pot, maxm_pot
! maxk_scha = 18 ;
! maxm_scha = 6 ;
! maxl_pot = 12 ;
! maxm_pot = 2 ;
!
      istat = nf_inq_dimid(ncid,"maxk_scha",id) 
      istat = nf_inq_dimlen(ncid,id,maxk_scha)
      istat = nf_inq_dimid(ncid,"maxm_scha",id) 
      istat = nf_inq_dimlen(ncid,id,maxm_scha)
      istat = nf_inq_dimid(ncid,"maxl_pot",id) 
      istat = nf_inq_dimlen(ncid,id,maxl_pot)
      istat = nf_inq_dimid(ncid,"maxm_pot",id) 
      istat = nf_inq_dimlen(ncid,id,maxm_pot)
!     write(6,"('wei05sc: maxk_scha=',i3,' maxm_scha=',i3)") 
!    |  maxk_scha,maxm_scha
!     write(6,"('wei05sc: maxl_pot=',i3,' maxm_pot=',i3)") 
!    |  maxl_pot,maxm_pot
!
! Read variables:
!
! double bndya(na):
      istat = nf_inq_varid(ncid,'bndya',id)
      istat = nf_get_var_double(ncid,id,bndya)
!     write(6,"('wei05sc: bndya=',/,(8f8.3))") bndya
!
! double bndyb(nb):
      istat = nf_inq_varid(ncid,'bndyb',id)
      istat = nf_get_var_double(ncid,id,bndyb)
!     write(6,"('wei05sc: bndyb=',/,(8f8.3))") bndyb
!
! double ex_bndy(nex):
      istat = nf_inq_varid(ncid,'ex_bndy',id)
      istat = nf_get_var_double(ncid,id,ex_bndy)
!     write(6,"('wei05sc: ex_bndy=',/,(8f8.3))") ex_bndy
!
! double th0s(n3_scha):
      istat = nf_inq_varid(ncid,'th0s',id)
      istat = nf_get_var_double(ncid,id,th0s)
!     write(6,"('wei05sc: th0s=',/,(8f8.3))") th0s
!
! double allnkm(n1_scha,n2_scha,n3_scha):
      istat = nf_inq_varid(ncid,'allnkm',id)
      istat = nf_get_var_double(ncid,id,allnkm)
!     write(6,"('wei05sc: allnkm min,max=',2e12.4)")
!    |  minval(allnkm),maxval(allnkm)
!
! int ab(csize):
      istat = nf_inq_varid(ncid,'ab',id)
      istat = nf_get_var_int(ncid,id,ab)
!     write(6,"('wei05sc: ab=',/,(10i4))") ab
!
! int ls(csize):
      istat = nf_inq_varid(ncid,'ls',id)
      istat = nf_get_var_int(ncid,id,ls)
!     write(6,"('wei05sc: ls=',/,(10i4))") ls
!
! int ms(csize):
      istat = nf_inq_varid(ncid,'ms',id)
      istat = nf_get_var_int(ncid,id,ms)
!     write(6,"('wei05sc: ms=',/,(10i4))") ms
!
! double ex_epot(nex):
      istat = nf_inq_varid(ncid,'ex_epot',id)
      istat = nf_get_var_double(ncid,id,ex_epot)
!     write(6,"('wei05sc: ex_epot=',/,(8f8.3))") ex_epot
!
! double ex_bpot(nex):
      istat = nf_inq_varid(ncid,'ex_bpot',id)
      istat = nf_get_var_double(ncid,id,ex_bpot)
!     write(6,"('wei05sc: ex_bpot=',/,(8f8.3))") ex_bpot
!
! double epot_schfits(csize,n_schfits):
      istat = nf_inq_varid(ncid,'epot_schfits',id)
      istat = nf_get_var_double(ncid,id,epot_schfits)
!     write(6,"('wei05sc: epot_schfits min,max=',2e12.4)")
!    |  minval(epot_schfits),maxval(epot_schfits)
!
! double bpot_schfits(csize,n_schfits):
      istat = nf_inq_varid(ncid,'bpot_schfits',id)
      istat = nf_get_var_double(ncid,id,bpot_schfits)
!     write(6,"('wei05sc: bpot_schfits min,max=',2e12.4)")
!    |  minval(bpot_schfits),maxval(bpot_schfits)
!
! double epot_alschfits(csize,n_alschfits):
      istat = nf_inq_varid(ncid,'epot_alschfits',id)
      istat = nf_get_var_double(ncid,id,epot_alschfits)
!     write(6,"('wei05sc: epot_alschfits min,max=',2e12.4)")
!    |  minval(epot_alschfits),maxval(epot_alschfits)
!
! double bpot_alschfits(csize,n_alschfits):
      istat = nf_inq_varid(ncid,'bpot_alschfits',id)
      istat = nf_get_var_double(ncid,id,bpot_alschfits)
!     write(6,"('wei05sc: bpot_alschfits min,max=',2e12.4)")
!    |  minval(bpot_alschfits),maxval(bpot_alschfits)
!
! Close file:
      istat = nf_close(ncid)
      write(6,"('wei05sc: completed read of file ',a)") trim(file)

      end subroutine read_wei05sc_ncfile
!-----------------------------------------------------------------------
      subroutine setmodel(angle,bt,tilt,swvel,swden,model)
!
! Calculate the complete set of the models' SCHA coeficients,
!   given an aribitrary IMF angle (degrees from northward toward +Y),
!   given byimf, bzimf, solar wind velocity (km/sec), and density.
!
      implicit none
!
! Args:
      real,intent(in) :: angle,bt,tilt,swvel,swden
      character(len=*),intent(in) :: model
!
! Local:
      integer :: i,j
      real :: pi,stilt,stilt2,sw,swp,swe,c0,
     |  rang,cosa,sina,cos2a,sin2a
      real :: a(n_schfits)
!
      if (trim(model) /= 'epot'.and.trim(model) /= 'bpot') then
        write(6,"('>>> model=',a)") trim(model)
        write(6,"('>>> setmodel: model must be either',
     |    '''epot'' or ''bpot''')")
        call shutdown('setmodel')
      endif
!
      pi = 4.*atan(1.)
      rad2deg = 180./pi
      deg2rad = pi/180.
!
!     write(6,"('setmodel call setboundary: model=',a,' swvel=',e12.4)")
!    |  model, swvel

      call setboundary(angle,bt,tilt,swvel,swden)
!
      stilt = sin(tilt*deg2rad)
      stilt2 = stilt**2
      sw = bt*swvel/1000.
      if (trim(model) == 'epot') then
        swe = (1.-exp(-sw*ex_epot(2)))*sw**ex_epot(1)
      else
        swe = (1.-exp(-sw*ex_bpot(2)))*sw**ex_bpot(1)
      endif
      c0 = 1.
      swp = swvel**2 * swden*1.6726e-6
      rang = angle*deg2rad
      cosa = cos(rang)
      sina = sin(rang)
      cos2a = cos(2.*rang)
      sin2a = sin(2.*rang)
      if (bt < 1.) then ! remove angle dependency for IMF under 1 nT
        cosa = -1.+bt*(cosa+1.)
        cos2a = 1.+bt*(cos2a-1.)
        sina = bt*sina
        sin2a = bt*sin2a
      endif
      a = (/c0      , swe       , stilt      , stilt2     , swp,  
     |      swe*cosa, stilt*cosa, stilt2*cosa, swp*cosa,          
     |      swe*sina, stilt*sina, stilt2*sina, swp*sina,          
     |      swe*cos2a,swe*sin2a/)
      if (trim(model) == 'epot') then
        esphc(:) = 0.
        do j=1,csize
          do i=1,n_schfits
            esphc(j) = esphc(j)+epot_schfits(i,j)*a(i)
          enddo
        enddo
!       write(6,"('setmodel: esphc=',/,(6e12.4))") esphc
      else
        bsphc(:) = 0.
        do j=1,csize
          do i=1,n_schfits
            bsphc(j) = bsphc(j)+bpot_schfits(i,j)*a(i)
          enddo
        enddo
!       write(6,"('setmodel: bsphc=',/,(6e12.4))") bsphc
      endif
      end subroutine setmodel
!-----------------------------------------------------------------------
      subroutine setboundary(angle,bt,tilt,swvel,swden)
!
! Sets the coefficients that define the low-latitude boundary model,
!   given the IMF and solar wind values.
!
      implicit none
!
! Args:
      real,intent(in) :: angle,bt,tilt,swvel,swden
!
! Local:
      integer :: i
      real :: swp,xc,theta,ct,st,tilt2,cosa,btx,x(na),c(na)
!
! Calculate the transformation matrix to the coordinate system
! of the offset pole.
!
      xc = 4.2
      theta = xc*(deg2rad)
      ct = cos(theta)
      st = sin(theta)
!
      tmat(1,:) = (/ ct, 0., st/) 
      tmat(2,:) = (/ 0., 1., 0./) 
      tmat(3,:) = (/-st, 0., ct/)
!
      ttmat(1,:) = (/ct, 0.,-st/)
      ttmat(2,:) = (/ 0.,1., 0./)
      ttmat(3,:) = (/st, 0., ct/)
!
      swp = swden*swvel**2*1.6726e-6 ! pressure
      tilt2 = tilt**2
      cosa = cos(angle*deg2rad)
      btx = 1.-exp(-bt*ex_bndy(1))
      if (bt > 1.) then
        btx = btx*bt**ex_bndy(2)
      else
        cosa = 1.+bt*(cosa-1.) ! remove angle dependency for IMF under 1 nT
      endif
      x = (/1., cosa, btx, btx*cosa, swvel, swp/)
      c = bndya
      bndyfitr = 0.
      do i=1,na
        bndyfitr = bndyfitr+x(i)*c(i)

!       write(6,"('setboundry: i=',i3,' bndyfitr=',e12.4)") i,bndyfitr

      enddo
      end subroutine setboundary
!-----------------------------------------------------------------------
      subroutine epotval(lat,mlt,fill,epot)
!
! Return the Potential (in kV) at given combination of def. latitude 
!   (lat) and MLT, in geomagnetic apex coordinates (practically identical 
!   to AACGM).  
! If the location is outside of the model's low-latitude boundary, then 
!   the value "fill" is returned.
!
      implicit none
!
! Args:
      real,intent(in) :: lat,mlt,fill
      real,intent(out) :: epot
!
! Local:
      integer :: inside,j,m,mm,skip
      real :: z,phir,plm,colat,nlm
      real :: phim(2),cospm(2),sinpm(2)
!
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
!
      call checkinputs(lat,mlt,inside,phir,colat)
      if (inside == 0) then
        epot = fill
        return
      endif
!
! IDL code: 
! phim=phir # replicate(1,maxm) * ((indgen(maxm)+1) ## replicate(1,n_elements(phir)))
!   where the '#' operator multiplies columns of first array by rows of second array,
!   and the '##' operator multiplies rows of first array by columns of second array.
! Here, maxm == maxm_pot == 2, and phir is a scalar. The above IDL statement then 
!   becomes: phim = ([phir] # [1,1]) * ([1,2] ## [phir]) where phim will be 
!   dimensioned [1,2]
!
      phim(1) = phir
      phim(2) = phir*2.
      cospm(:) = cos(phim(:))
      sinpm(:) = sin(phim(:))
!
      z = 0.
      skip = 0 ! Added by B.Foster, 4/23/14
      do j=1,csize
        if (skip == 1) then
          skip = 0
          cycle
        endif
        m = ms(j)
        if (ab(j)==1) then
          plm = scplm(j,colat,nlm) ! scplm function is in this module
          skip = 0
          if (m == 0) then
            z = z+plm*esphc(j)
          else
            z = z+plm*(esphc(j)*cospm(m)+esphc(j+1)*sinpm(m))
            skip = 1
          endif
        endif ! ab(j)
      enddo
      epot = z 
      end subroutine epotval
!-----------------------------------------------------------------------
      subroutine mpfac(lat,mlt,fill,fac)
      implicit none
!
! Args:
      real,intent(in) :: lat,mlt,fill
      real,intent(out) :: fac
!
! Local:
      integer :: j,m,inside,skip
      real :: phim(2),cospm(2),sinpm(2),cfactor
      real :: re,z,phir,plm,colat,nlm,pi
!
      re = 6371.2 + 110. ! km radius (allow default ht=110)
!
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
!
      call checkinputs(lat,mlt,inside,phir,colat)
      if (inside == 0) then
        fac = fill
        return
      endif
!
      phim(1) = phir
      phim(2) = phir*2.
      cospm(:) = cos(phim(:))
      sinpm(:) = sin(phim(:))
!
      z = 0.
      skip = 0 ! Added by B.Foster, 4/23/14
      jloop: do j=1,csize
        if (skip == 1) then
          skip = 0
          cycle
        endif
        if (ls(j) >= 11) exit jloop
        m = ms(j)
        if (ab(j) == 1) then
          plm = scplm(j,colat,nlm) ! colat and nlm are returned (both reals)
          plm = plm*(nlm*(nlm+1.))
!
! bsphc was calculated in setmodel (when setmodel called with 'bpot')
          if (m==0) then
            z = z-plm*bsphc(j)
          else
            z = z-(plm*(bsphc(j)*cospm(m)+bsphc(j+1)*sinpm(m)))
            skip = 1
          endif
        endif
      enddo jloop ! j=1,csize
      pi = 4.*atan(1.)
      cfactor = -1.e5/(4.*pi*re**2) ! convert to uA/m2
      z = z*cfactor
      fac = z
!     write(6,"('mpfac: lat=',f8.3,' mlt=',f8.3,' fac=',1pe12.4)") 
!    |  lat,mlt,fac
      end subroutine mpfac
!-----------------------------------------------------------------------
      real function scplm(index,colat,nlm)
!
! Return Spherical Cap Harmonic Associated Legendre values, given colat 
!   values and index i into array of L and M values.
!
      implicit none
!
! Args:
      integer,intent(in) :: index
      real,intent(in) :: colat
      real,intent(out) :: nlm
!
! Local:
      integer :: istat,i,j,l,m,skip
      real :: th0,out(1),colata(1),plm1
      real :: cth(mxtablesize)
      real,save :: prevth0=1.e36
      integer,save :: tablesize
!
      scplm = 0.
      skip = 0 ! Added by B.Foster, 4/23/14
      th0 = bndyfitr
      if (prevth0 /= th0) then
        tablesize = 3*nint(th0)
        if (tablesize > mxtablesize) then 
          write(6,"('>>> tablesize > mxtablesize: tablesize=',i8,
     |       ' mxtablesize=',i8,' th0=',e12.4)") tablesize,mxtablesize,
     |       th0
          call shutdown('tablesize')
        endif
        do i=1,tablesize
          colattable(i) = real(i-1)*(th0/real(tablesize-1))
          cth(i) = cos(colattable(i)*deg2rad)
        enddo
        prevth0 = th0
        nlms = 0. ! whole array init 
        do j=1,csize
          if (skip == 1) then
            skip = 0
            cycle
          endif
          l = ls(j)
          m = ms(j)
          nlms(j) = nkmlookup(l,m,th0) ! nkmlookup in this module

! real :: plmtable(mxtablesize,csize)
          call pm_n(m,nlms(j),cth,plmtable(1:tablesize,j),tablesize)
          skip = 0
          if (m /= 0 .and. ab(j) > 0) then
            plmtable(1,j+1) = plmtable(1,j)
            nlms(j+1) = nlms(j)
            skip = 1
          endif
        enddo ! j=1,csize
      endif ! prevth0
      nlm = nlms(index)
      colata(1) = colat
      call interpol_quad(plmtable(1:tablesize,index),
     |  colattable(1:tablesize),colata,out)
      scplm = out(1)
      end function scplm
!-----------------------------------------------------------------------
      subroutine pm_n(m,r,cth,plmtable,tablesize)
!
! Another SCHA function, returns the SCHA version of the associated 
! Legendre Polynomial, Pmn
!
      implicit none
!
! Args:
      integer,intent(in) :: m,tablesize
      real,intent(in) :: r
      real,intent(in) :: cth(tablesize)
      real,intent(out) :: plmtable(tablesize)
!
! Local:
      integer :: i,k,ii
      real :: rm,rk,div,ans,xn
      real,dimension(tablesize) :: a,x,tmp,table
!
      if (m == 0) then 
        a = 1. ! whole array op
      else
        do i=1,tablesize
          a(i) = sqrt(1.-cth(i)**2)**m
        enddo
      endif
      xn = r*(r+1.)
      x(:) = (1.-cth(:))/2.
      table = a ! whole array init
      k = 1
      pmn_loop: do         ! repeat-until loop in idl code
        do i=1,tablesize
          rm = real(m)
          rk = real(k)
          a(i) = a(i)*(x(i)*((rk+rm-1.)*(rk+rm)-xn)/(rk*(rk+rm)))
          table(i) = table(i)+a(i) ! "result" in idl code
        enddo
        k = k+1
        do i=1,tablesize
          div = abs(table(i))
          if (div <= 1.e-6) div = 1.e-6
          tmp(i) = abs(a(i)) / div
        enddo
        if (maxval(tmp) < 1.e-6) exit pmn_loop
      enddo pmn_loop
      ans = km_n(m,r)

      plmtable(:) = table(:)*ans
      end subroutine pm_n
!-----------------------------------------------------------------------
      real function km_n(m,rn)
!
! A normalization function used by the SCHA routines.  See Haines.
!
      implicit none
!
! Args:
      integer,intent(in) :: m
      real,intent(in) :: rn
!
! Local:
      integer :: i,n
      real :: rm
!
      if (m == 0) then 
        km_n = 1.
        return
      endif
      rm = real(m)
      km_n = sqrt(2.*exp(lngamma(rn+rm+1.)-lngamma(rn-rm+1.))) / 
     |  (2.**m*factorial(m))
      end function km_n
!-----------------------------------------------------------------------
      real function nkmlookup(k,m,th0)
!
! Given the size of a spherical cap, defined by the polar cap angle, th0, 
!   and also the values of integers k and m, returns the value of n, a 
!   real number (see Haines).
! It uses interpolation from a lookup table that had been precomputed, 
!   in order to reduce the computation time.
!
      implicit none
!
! Args:
      integer,intent(in) :: k,m
      real,intent(in) :: th0
!
! Local:
      integer :: kk,mm
      real :: th0a(1),out(1)

      if (th0 == 90.) then
        nkmlookup = real(k)
        return
      endif
      th0a(1) = th0
      kk = k+1
      mm = m+1
      if (kk > maxk_scha) then
        call interpol_quad(allnkm(maxk_scha,mm,:),th0s,th0a,out)
      endif
      if (mm > maxm_scha) then
        call interpol_quad(allnkm(kk,maxm_scha,:),th0s,th0a,out)
      endif
      if (th0 < th0s(1)) then
        write(6,"('>>> nkmlookup: th0 < th0s(1): th0=',e12.4,
     |    ' th0s(1)=',e12.4)") th0,th0s(1)
      endif
      call interpol_quad(allnkm(kk,mm,:),th0s,th0a,out)
      nkmlookup = out(1)
      end function nkmlookup
!-----------------------------------------------------------------------
      subroutine checkinputs(lat,mlt,inside,phir,colat)
      implicit none
!
! Args:
      real,intent(in) :: lat,mlt
      integer,intent(out) :: inside
      real,intent(out) :: phir,colat
!
! Local:
      real :: lon,tlat,tlon,radii
!
      lon = mlt*15.
      call dorotation(lat,lon,tlat,tlon)
      radii = 90.-tlat
      inside = 0
      if (radii <= bndyfitr) inside = 1 ! bndyfitr from setboundary
      phir = tlon*deg2rad
      colat = radii
      end subroutine checkinputs
!-----------------------------------------------------------------------
      subroutine dorotation(latin,lonin,latout,lonout)
!
! Uses transformation matrices tmat and ttmat, to convert between
!   the given geomagnetic latatud/longitude, and the coordinate 
!   system that is used within the model,that is offset from the pole.
!
! Rotate Lat/Lon spherical coordinates with the transformation given 
!   by saved matrix. The coordinates are assumed to be on a sphere of 
!   Radius=1. Uses cartesian coordinates as an intermediate step.
!
      implicit none
!
! Args:
      real,intent(in) :: latin,lonin
      real,intent(out) :: latout,lonout
!
! Local:
      real :: latr,lonr,stc,ctc,sf,cf,a,b,pos(3)
      integer :: i
!
      latr = latin*deg2rad
      lonr = lonin*deg2rad
      stc = sin(latr)
      ctc = cos(latr)
      sf = sin(lonr)
      cf = cos(lonr)
      a = ctc*cf
      b = ctc*sf
!
! IDL code: Pos= TM ## [[A],[B],[STC]]
! The ## operator multiplies rows of first array by columns of second array.
! Currently, TM(3,3) = Tmat (or TTmat if "reversed" was set)
! If called w/ single lat,lon, then a,b,stc are dimensioned (1), and
!   Pos is then (1,3)
!
      do i=1,3
        pos(i) = tmat(1,i)*a + tmat(2,i)*b + tmat(3,i)*stc
      enddo
      latout = asin(pos(3))*rad2deg
      lonout = atan2(pos(2),pos(1))*rad2deg
      end subroutine dorotation
!-----------------------------------------------------------------------
      subroutine interpol_quad(v,x,u,p)
!
! f90 translation of IDL function interpol(v,x,u,/quadratic)
!
      implicit none
!
! Args:
      real,intent(in) :: v(:),x(:),u(:)
      real,intent(out) :: p(:)
!
! Local:
      integer :: nv,nx,nu,i,ix
      real :: x0,x1,x2
!
      nv = size(v)
      nx = size(x)
      nu = size(u)
      if (nx /= nv) then
        p(:) = 0.
        return
      endif
      do i=1,nu
        ix = value_locate(x,u(i))
! 01/14 bae: interpol_quad in wei05sc.F is called when inside=1 or radii<bndryfit 
!  for Weimer 2005. The fix to ix<=1 and ix>=nx assures epot is non-zero near 
!  the pole (85.8mlat,0MLT) and the boundary (bndryfit).
        if (ix <=1) ix = 2	      ! bug fix by bae 01/28/14
        if (ix >=nx) ix = nx-1	      ! bug fix by bae 01/29/14
!       if (ix <= 1.or.ix >= nx) then ! bug fix by btf 12/23/09
!          p(i) = 0.
!          cycle                       ! bug fix by btf 12/23/09
!        endif
        x1 = x(ix)
        x0 = x(ix-1)
        x2 = x(ix+1)
        p(i) = v(ix-1) * (u(i)-x1) * (u(i)-x2) / ((x0-x1) * (x0-x2)) +  
     |         v(ix)   * (u(i)-x0) * (u(i)-x2) / ((x1-x0) * (x1-x2)) +  
     |         v(ix+1) * (u(i)-x0) * (u(i)-x1) / ((x2-x0) * (x2-x1))
      enddo
      end subroutine interpol_quad
!-----------------------------------------------------------------------
      integer function value_locate(vec,val)
!
! f90 translation of IDL function value_locate
! Return index i into vec for which vec(i) <= val >= vec(i+1)
! Input vec must be monotonically increasing
!
      implicit none
!
! Args:
      real,intent(in) :: vec(:),val
!
! Local:
      integer :: n,i
!
      value_locate = 0
      n = size(vec)
      if (val < vec(1)) return
      if (val > vec(n)) then
        value_locate = n
        return
      endif
      do i=1,n-1
        if (val >= vec(i) .and. val <= vec(i+1)) then
          value_locate = i
          return
        endif
      enddo
      end function value_locate
!-----------------------------------------------------------------------
      real function lngamma(xx)
!
! This is an f90 translation from C code copied from 
! www.fizyka.umk.pl/nrbook/c6-1.pdf (numerical recipes gammln)
!
      implicit none
      real,intent(in) :: xx
      real :: x,y,tmp,ser
      real :: cof(6) = (/76.18009172947146, -86.50532032941677, 
     |  24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, 
     |  -0.5395239384953e-5/)
      integer :: j
!
      y = xx
      x = xx
      tmp = x+5.5
      tmp = tmp-(x+0.5)*log(tmp)
      ser = 1.000000000190015
      do j=1,5
        y = y+1
        ser = ser+cof(j)/y
      enddo
      lngamma = -tmp+log(2.5066282746310005*ser/x)
      end function lngamma
!-----------------------------------------------------------------------
      real function factorial(n)
      implicit none
      integer,intent(in) :: n
      integer :: m
      if (n <= 0) then
        factorial = 0.
        return
      endif
      if (n == 1) then
        factorial = 1.
        return
      endif
      factorial = real(n)
      do m = n-1,1,-1
        factorial = factorial * real(m)
      enddo
      end function factorial
!-----------------------------------------------------------------------
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
* COORDINATE TRANSFORMATION UTILITIES

CNCAR      Feb 01:  Changed TRANS to GET_TILT s.t. the dipole tilt angle is
C          returned.

	FUNCTION GET_TILT (YEAR,MONTH,DAY,HOUR)
C       SUBROUTINE TRANS(YEAR,MONTH,DAY,HOUR,IDBUG)
	implicit none
	real :: GET_TILT,B3,B32,B33
	integer :: IYR,JD,MJD,I,J,K
CNCAR

	INTEGER YEAR,MONTH,DAY,IDBUG
	REAL HOUR
C         
C      THIS SUBROUTINE DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
C      TRANSFORMATIONS, IDENTIFIED BY K.
C          K=1 TRANSFORMS GSE to GEO
C          K=2     "      GEO to MAG
C          K=3     "      GSE to MAG
C          K=4     "      GSE to GSM
C          K=5     "      GEO to GSM
C          K=6     "      GSM to MAG
C          K=7     "      GSE to GEI
C          K=8     "      GEI to GEO
C          K=9     "      GSM to SM 
C	   K=10    "      GEO to SM 
C	   K=11    "      MAG to SM 
C
C      IF IDBUG IS NOT 0, THEN OUTPUTS DIAGNOSTIC INFORMATION TO
C      FILE UNIT=IDBUG
C	
	INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
	INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
	INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
	INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
	INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG

	PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
	PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
	PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
	PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
	PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
	PARAMETER (MAGSM =11,SMMAG =-11)
C
C      The formal names of the coordinate systems are:
C	GSE - Geocentric Solar Ecliptic
C	GEO - Geographic
C	MAG - Geomagnetic
C	GSM - Geocentric Solar Magnetospheric
C	SM  - Solar Magnetic
C	
C      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
C      ST(I) AND CT(I) ARE SINES & COSINES.       
C
C      Program author:  D. R. Weimer
C
C      Some of this code has been copied from subroutines which had been
C      obtained from D. Stern, NASA/GSFC.  Other formulas are from "Space 
C      Physics Coordinate Transformations: A User Guide" by M. Hapgood (1991).
C
C      The formulas for the calculation of Greenwich mean sidereal time (GMST)
C      and the sun's location are from "Almanac for Computers 1990",
C      U.S. Naval Observatory.
C
	REAL UT,T0,GMSTD,GMSTH,ECLIP,MA,LAMD,SUNLON

CNCAR      Feb 01:  Eliminate unused routines from translib.for: ROTATE,
C          ROTATEV, FROMCART, TOCART, MLT, MAGLONG, SUNLOC.  Remaining
C          are ADJUST and JULDAY
CNCAR      Nov 02: Commons MFIELD and TRANSDAT now only in TRANS (GET_TILT)
C                  So eliminate them as commons.  For Fortran 90, eliminate
C                  the DATA statement for assignments (not block_data)
C       COMMON/MFIELD/EPOCH,TH0,PH0,DIPOLE
C       COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
C         
C       DATA EPOCH,TH0,PH0,DIPOLE/1980.,11.19,-70.76,.30574/
	REAL EPOCH,TH0,PH0,DIPOLE
	REAL CX(9),ST(6),CT(6),AM(3,3,11)
C         
C  TH0 = geog co-lat of NH magnetic pole
C  PH0 = geog longitude of NH magnetic pole
C  DIPOLE = magnitude of the B field in gauss at the equator
	EPOCH = 1980.
	TH0 = 11.19
	PH0 = -70.76
	DIPOLE = 0.30574
CNCAR

CNCAR      Feb 01:  Prevent debug printing to a file
	IDBUG = 0
CNCAR

	IF(YEAR.LT.1900)THEN
	  IYR=1900+YEAR
	ELSE
	  IYR=YEAR
	ENDIF
	UT=HOUR
	JD=JULDAY(MONTH,DAY,IYR)
	MJD=JD-2400001
	T0=(real(MJD)-51544.5)/36525.0
	GMSTD=100.4606184 + 36000.770*T0 + 3.87933E-4*T0*T0 +
     $        15.0410686*UT
	CALL ADJUST(GMSTD)
	GMSTH=GMSTD*24./360.
	ECLIP=23.439 - 0.013*T0
	MA=357.528 + 35999.050*T0 + 0.041066678*UT
	CALL ADJUST(MA)
	LAMD=280.460 + 36000.772*T0 + 0.041068642*UT
	CALL ADJUST(LAMD)
	SUNLON=LAMD + (1.915-0.0048*T0)*SIND(MA) + 0.020*SIND(2.*MA)
	CALL ADJUST(SUNLON)
	IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) YEAR,MONTH,DAY,HOUR
	  WRITE(IDBUG,*) 'MJD=',MJD
	  WRITE(IDBUG,*) 'T0=',T0
	  WRITE(IDBUG,*) 'GMSTH=',GMSTH
	  WRITE(IDBUG,*) 'ECLIPTIC OBLIQUITY=',ECLIP
	  WRITE(IDBUG,*) 'MEAN ANOMALY=',MA
	  WRITE(IDBUG,*) 'MEAN LONGITUDE=',LAMD
	  WRITE(IDBUG,*) 'TRUE LONGITUDE=',SUNLON
	ENDIF

	CX(1)= GMSTD
	CX(2) = ECLIP
	CX(3) = SUNLON
	CX(4) = TH0
	CX(5) = PH0
c Derived later:
c       CX(6) = Dipole tilt angle  
c       CX(7) = Angle between sun and magnetic pole
c       CX(8) = Subsolar point latitude
c       CX(9) = Subsolar point longitude

	DO I=1,5
	  ST(I) = SIND(CX(I))
	  CT(I) = COSD(CX(I))
	ENDDO
C         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0.         
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)      
C         
      AM(1,1,GEIGEO) = CT(1)      
      AM(1,2,GEIGEO) = ST(1)      
      AM(1,3,GEIGEO) = 0.         
      AM(2,1,GEIGEO) = -ST(1)     
      AM(2,2,GEIGEO) = CT(1)      
      AM(2,3,GEIGEO) = 0.         
      AM(3,1,GEIGEO) = 0.         
      AM(3,2,GEIGEO) = 0.         
      AM(3,3,GEIGEO) = 1.         
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) +
     $    AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +  AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
C         
      AM(1,1,GEOMAG) = CT(4)*CT(5) 
      AM(1,2,GEOMAG) = CT(4)*ST(5) 
      AM(1,3,GEOMAG) =-ST(4)       
      AM(2,1,GEOMAG) =-ST(5)       
      AM(2,2,GEOMAG) = CT(5)       
      AM(2,3,GEOMAG) = 0.
      AM(3,1,GEOMAG) = ST(4)*CT(5) 
      AM(3,2,GEOMAG) = ST(4)*ST(5) 
      AM(3,3,GEOMAG) = CT(4)       
C         
      DO I=1,3   
      DO J=1,3   
       AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) +
     $   AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +  AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
C         
      B32 = AM(3,2,GSEMAG)         
      B33 = AM(3,3,GSEMAG)         
      B3  = SQRT(B32*B32+B33*B33)       
      IF (B33.LE.0.) B3 = -B3    
C         
      AM(2,2,GSEGSM) = B33/B3      
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)   
      AM(3,2,GSEGSM) = B32/B3      
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)   
      AM(1,1,GSEGSM) = 1.
      AM(1,2,GSEGSM) = 0.
      AM(1,3,GSEGSM) = 0.
      AM(2,1,GSEGSM) = 0.
      AM(3,1,GSEGSM) = 0.
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) +
     $    AM(I,2,GSEGSM)*AM(J,2,GSEGEO) + AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) +
     $   AM(I,2,GEOMAG)*AM(J,2,GEOGSM) + AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO
C
	ST(6) = AM(3,1,GSEMAG)       
	CT(6) = SQRT(1.-ST(6)*ST(6))      
	CX(6) = ASIND(ST(6))     

        AM(1,1,GSMSM) = CT(6)
        AM(1,2,GSMSM) = 0.
        AM(1,3,GSMSM) = -ST(6)
        AM(2,1,GSMSM) = 0.
        AM(2,2,GSMSM) = 1.
        AM(2,3,GSMSM) = 0.
        AM(3,1,GSMSM) = ST(6)
        AM(3,2,GSMSM) = 0.
        AM(3,3,GSMSM) = CT(6)
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOSM) = AM(I,1,GSMSM)*AM(1,J,GEOGSM) +
     $    AM(I,2,GSMSM)*AM(2,J,GEOGSM) +  AM(I,3,GSMSM)*AM(3,J,GEOGSM)
      ENDDO
      ENDDO
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,MAGSM) = AM(I,1,GSMSM)*AM(J,1,GSMMAG) +
     $   AM(I,2,GSMSM)*AM(J,2,GSMMAG) + AM(I,3,GSMSM)*AM(J,3,GSMMAG)
      ENDDO
      ENDDO
C
      CX(7)=ATAN2D( AM(2,1,11) , AM(1,1,11) )
      CX(8)=ASIND( AM(3,1,1) )
      CX(9)=ATAN2D( AM(2,1,1) , AM(1,1,1) )

      IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) 'Dipole tilt angle=',CX(6)
	  WRITE(IDBUG,*) 'Angle between sun and magnetic pole=',CX(7)
	  WRITE(IDBUG,*) 'Subsolar point latitude=',CX(8)
	  WRITE(IDBUG,*) 'Subsolar point longitude=',CX(9)

        DO K=1,11
         WRITE(IDBUG,1001) K
         DO I=1,3
           WRITE(IDBUG,1002) (AM(I,J,K),J=1,3)
         ENDDO
        ENDDO
 1001   FORMAT(' ROTATION MATRIX ',I2)
 1002   FORMAT(3F9.5)
      ENDIF

CNCAR      Mar 96: return the dipole tilt from this function call.
      GET_TILT = CX(6)
CNCAR

      RETURN
!     END
	end function get_tilt
!-----------------------------------------------------------------------
******************************************************************************
CNCAR      Feb 01:  Eliminate unused routines from translib.for: ROTATE,
C          ROTATEV, FROMCART, TOCART, MLT, MAGLONG, SUNLOC.  Remaining
C          are ADJUST and JULDAY
CNCAR
	SUBROUTINE ADJUST(ANGLE)
	implicit none
	real :: angle
C	ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
 10	CONTINUE
	IF(ANGLE.LT.0.)THEN
	  ANGLE=ANGLE+360.
	  GOTO 10
	ENDIF
 20	CONTINUE
 	IF(ANGLE.GE.360.)THEN
	  ANGLE=ANGLE-360.
	  GOTO 20
	ENDIF
	RETURN
!       END
	end subroutine adjust
******************************************************************************
      FUNCTION JULDAY(MM,ID,IYYY)
      implicit none
      integer :: igreg, iyyy, mm, id, jy, jm, ja, julday
      PARAMETER (IGREG=15+31*(10+12*1582))
      IF (IYYY.EQ.0) call shutdown('There is no Year Zero.')
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF
      RETURN
!     END
	end function julday
!-----------------------------------------------------------------------
      SUBROUTINE CVT2MD(MSGUN,IYEAR,NDA,MON,DAY)
C     SUBROUTINE CVT2MD(MSGUN,IYEAR,NDA,MMDD)
C          This sub converts NDA, the day number of the year, IYEAR,
C          into the appropriate month and day of month (integers)
      implicit none
      integer :: msgun,iyear,nda,mon,miss,numd,i,mmdd
      INTEGER DAY
      INTEGER LMON(12)
      PARAMETER (MISS=-32767)
      SAVE LMON
      DATA LMON/31,28,31,30,31,30,31,31,30,31,30,31/
 
      LMON(2)=28
      IF(MOD(IYEAR,4) .EQ. 0)LMON(2)=29
 
      NUMD=0
      DO 100 I=1,12
      IF(NDA.GT.NUMD .AND. NDA.LE.NUMD+LMON(I))GO TO 200
      NUMD=NUMD+LMON(I)
  100 CONTINUE
      WRITE(MSGUN,'('' CVT2MD:  Unable to convert year & day of year'',
     +              I5,'','',I5,''to month & day of month'')')IYEAR,NDA
      MMDD = MISS
      MON = MISS
      DAY = MISS
      RETURN
  200 MON=I
      DAY=NDA-NUMD
      MMDD = MON*100 + DAY
 
      RETURN
!     END
	end subroutine cvt2md

CNCAR      Routines added to work around non-ANSI trig functions which
C          input degrees instead of radians:  SIND, COSD, ASIND, ATAN2D

      FUNCTION SIND (DEG)
      implicit none
      real :: sind,d2r,r2d,deg
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      SIND = SIN (DEG * D2R)
      RETURN
!     END
	end function sind

      FUNCTION COSD (DEG)
      implicit none
      real :: cosd,d2r,r2d,deg
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)

      COSD = COS (DEG * D2R)
      RETURN
!     END
	end function cosd

      FUNCTION ASIND (RNUM)
      implicit none
      real :: asind,d2r,r2d,rnum
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      ASIND = R2D * ASIN (RNUM)
      RETURN
!     END
	end function asind

      FUNCTION ATAN2D (RNUM1,RNUM2)
      implicit none
      real :: atan2d,d2r,r2d,rnum1,rnum2
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      ATAN2D = R2D * ATAN2 (RNUM1,RNUM2)
      RETURN
!     END
	end function atan2d

!-----------------------------------------------------------------------
      subroutine wei05loc (ih)
! ih=1,2 for SH,NH called from weimer05

!  Jan 2012 bae:  added byloc to limit size of byimf (probably incorrect assymmetry in hems)
!  NCAR July 08 bae:  Write wei05loc to use values from Weimer 2005
!       setboundary for convection offc=4.2 deg, dskofc=0, theta0=bndyfitr
!       where deg should be in radians.  wei05loc replaces calccloc05.
!
! (dimension 2 is for south, north hemispheres)
!  Calculate offa, dskofa, rrad, phid, and phin from Weimer 2005 offc, dskofc, theta0
!   Use Fig 8 of Heelis et al. [JGR, 85, 3315-3324, 1980]
!     This shows:  arad = 18.7 deg, crad = 16.7 deg (so arad = crad + 2 deg)
!           offa = offc = 3 deg (so offa = offc)
!           dskofc = 2 deg, dskofa = -0.5 deg  (so dskofa = dskofc - 2.5 deg)
!   Parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                             and phin (phin(MLT)=23.50 +/- 0.15By - 12)
!   (In aurora_cons, phid=0., phin=180.*rtd)
!     (For zero By, should be phid=21.39MLT*15*rtd, phin=11.5*15*rtd)
!  05/08:  But formulae for ra-rc using IMF CP between 5-7 deg (not 2) so use
!          difference of ra(max IMF CP or HP)-rc(IMF CP) as in aurora.F
! These are the dimensions and descriptions (corrected phid,n) from aurora.F:
!    |  theta0(2), ! convection reversal boundary in radians
!    |  offa(2),   ! offset of oval towards 0 MLT relative to magnetic pole (rad)
!    |  dskofa(2), ! offset of oval in radians towards 18 MLT (f(By))
!    |  phid(2),   ! dayside convection entrance in MLT-12 converted to radians (f(By))
!                      phid is the MLT-12 location of the cusp on the dayside
!    |  phin(2),   ! night convection entrance in MLT-12 converted to radians (f(By))
!    |  rrad(2),   ! radius of auroral circle in radians
!    |  offc(2),   ! offset of convection towards 0 MLT relative to mag pole (rad)
!    |  dskofc(2)  ! offset of convection in radians towards 18 MLT (f(By))
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc)
!
      use pdynamo_module,only: nmlat0,phihm
! Additional auroral parameters (see sub aurora_cons):
      use aurora_module,only: theta0,offa,dskofa,phid,phin,rrad,offc,
     |                        dskofc
      use magfield_module,only: sunlons 
      use input_module,only:  ! from user input
     |  byimf    ! By component of IMF (nT)      (e.g., 0.)
! 05/08:  Added bzimf,ctpoten,power
     | ,bzimf    ! Bz component of IMF (nT)      (e.g., 0.)
     | ,ctpoten  ! cross-cap potential (kV)      (e.g., 45.)
     | ,power    ! hemispheric power   (GW)      (e.g., 16.)
      use cons_module,only: rtd,
     |  ylonm,ylatm  ! magnetic grid lons, lats in radians
      implicit none
!
! Args:
      integer,intent(in) :: ih
!
! Local:
! 05/08:  Added rccp etc
      real :: rccp,racp,rahp,ramx,diffrac,plevel,tmltmin,tmltmax
      real :: diffrac2
      real :: offcdegp(2),dskofcp(2),radegp(2)
      integer :: i,i1,i2,j,j1,j2
      real :: vnx(2,2),hem,mltd,mltn
      integer :: jnx(2,2),inx(2,2)
      real :: offcdeg,dskof,arad,crad
      real :: cosl06,sinl06,colat06,cosm06,cosl18,sinl18,colat18,cosm18
      real :: byloc

! Limit size of byimf in phin and phid calculations (as in aurora.F) 
!  NOTE:  This byloc is assymetric in hemisphere, which is probably not correct
      byloc = byimf
      if (byloc .gt. 7.) byloc = 7.
      if (byloc .lt. -11.) byloc = -11.
!
!  ih=1 is SH, ih=2 is NH
	if (ih .eq. 1) then
	  j1 = 1
	  j2 = nmlat0
	  hem = -1.
	else
	  j1 = nmlat0 + 1
	  j2 = nmlat
	  hem = 1.
	endif
! Print out un-revised values:
!       write (6,"(1x,'Original convection/oval params (hem,By,off,dsk',
!    |    ',rad,phid,n=',10f9.4)") hem,byimf,offc(ih)*rtd,offa(ih)*rtd,
!    |    dskofc(ih)*rtd,dskofa(ih)*rtd,theta0(ih)*rtd,rrad(ih)*rtd,
!    |    phid(ih)*rtd/15.+12.,phin(ih)*rtd/15.+12.
!  Find min/max
	vnx(ih,1) = 0.
	vnx(ih,2) = 0.
	do j=j1,j2
	  do i=1,nmlonp1-1
	    if (phihm(i,j) .gt. vnx(ih,2)) then
	      vnx(ih,2) = phihm(i,j)
	      jnx(ih,2) = j
	      inx(ih,2) = i
	    endif
	    if (phihm(i,j) .lt. vnx(ih,1)) then
	      vnx(ih,1) = phihm(i,j)
	      jnx(ih,1) = j
	      inx(ih,1) = i
	    endif
	  enddo  !  i=1,nmlonp1-1
	enddo  !  j=j1,j2
! 05/08: Calculate weictpoten in kV from Weimer model min/max in V
	weictpoten(ih) = 0.001 * (vnx(ih,2) - vnx(ih,1))
	tmltmin = (ylonm(inx(ih,1))-sunlons(1)) * rtd/15. + 12.
	if (tmltmin .gt. 24.) tmltmin = tmltmin - 24.
	tmltmax = (ylonm(inx(ih,2))-sunlons(1)) * rtd/15. + 12.
	if (tmltmax .gt. 24.) tmltmax = tmltmax - 24.
!       write (6,"('ih Bz By Hp ctpoten,wei min/max potV,lat,mlt=',i2,
!    |    5f8.2,2x,e12.4,2f8.2,2x,e12.4,2f8.2))") ih,bzimf,byimf,power,
!    |    ctpoten,weictpoten(ih),
!    |    vnx(ih,1),ylatm(jnx(ih,1))*rtd,tmltmin,
!    |    vnx(ih,2),ylatm(jnx(ih,2))*rtd,tmltmax
! 05/08: From aurora_cons, calculate convection and aurora radii using IMF convection
!   and power (plevel);  racp (DMSP/NOAA) - rccp (AMIE) = 5.32 (Bz>0) to 6.62 (Bz<0) deg
!  Heelis et al [1980, JGR, 85, pp 3315-3324] Fig 8: ra=rc+2deg, and is 2.5 deg to dusk
	rccp = -3.80+8.48*(weictpoten(ih)**0.1875)
	racp = -0.43+9.69*(weictpoten(ih)**0.1875)
	plevel = 0.
	if (power >=1.00) plevel = 2.09*alog(power)
	rahp = 14.20 + 0.96*plevel
	ramx = max(racp,rahp)
	diffrac = ramx - rccp

!  Set default values
!  Use parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                             and phin (phin(MLT)=23.50 +/- 0.15By - 12)
	mltd = 9.39 - hem*0.21*byloc
	mltn = 23.50 - hem*0.15*byloc
	phid(ih) = (mltd-12.) * 15. / rtd
	phin(ih) = (mltn-12.) * 15. / rtd
! 05/18/08:  Note that phid,phin are only for Heelis and are irrelevant for Weimer
!       write (6,"(1x,'mltd mltn phid,n =',4f8.2)")
!    |   mltd,mltn,phid(ih)*rtd/15.,phin(ih)*rtd/15.
!  Use default constant value of offcdegp from setboundary in Weimer 2005
	offcdeg = 4.2
	offcdegp(ih) = offcdeg
	offc(ih) = offcdegp(ih) / rtd
	offa(ih) = offcdegp(ih) / rtd
!       write (6,"(1x,'offcdeg,rad =',2e12.4)") offcdeg,offc(ih)
	dskof = 0.
	dskofcp(ih) = dskof
	dskofc(ih) = dskof / rtd
!  oval offset is 2.5 deg towards dawn (more neg dskof)
	dskofa(ih) = (dskof-2.5) / rtd
!       write (6,"(1x,'dskof,c,a=',3f8.2)")
!    |    dskof,dskofc(ih)*rtd,dskofa(ih)*rtd
! Set crad from bndyfitr/2 of setboundary of Weimer 2005
	crad = bndyfitr/2.
!      write (6,"(1x,'wei05loc: ih,bz,y,crad =',i2,3f8.2)") 
!    |    ih,bzimf,byimf,crad
!  Fig 8 Heelis et al [1980]: ra=rc+2deg, and shifted 2.5 deg to dusk
	arad = crad + 2.
! 05/08:  Make ra=rc+diffrac(=ramx-rccp) - same difference as in aurora.F
! Choose to have arad=crad(Weimer) + diffrac(same diff as in aurora.F)
	arad = crad + diffrac
! 08/08: OR make ra=ramx=max(racp,rahp) so diffrac=arad-crad
	diffrac2 = ramx - crad
! Choose to have arad=ramx (same as in aurora.F as determined by P/CP)
!       arad = ramx
	theta0(ih) = crad / rtd
	radegp(ih) = arad
	rrad(ih) = arad / rtd
!       write (6,"(1x,'radius: crad,rccp,racp,rahp diffa-c',
!    |   '(aurF,ramx-Weic) ramx,Weic+d,arad deg=',9f8.2)") crad,rccp,
!    |   racp,rahp,diffrac,diffrac2,ramx,crad+diffrac,arad

! Print out revised values (revised 05/08):
!       write (6,"(1x,'Revised convection/oval params (off,dsk,',
!    |    'rad,phid,n=',8f9.4)")offc(ih)*rtd,offa(ih)*rtd,
!    |    dskofc(ih)*rtd,dskofa(ih)*rtd,theta0(ih)*rtd,rrad(ih)*rtd,
!    |    phid(ih)*rtd/15.+12.,phin(ih)*rtd/15.+12.

      return
      end subroutine wei05loc
!-----------------------------------------------------------------------
      end module wei05sc
