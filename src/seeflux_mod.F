!
      module seeflux_module
!
! Use TIMED-SEE EUV flux data
      implicit none
! Or using fitted GOES netCDF SPE flux Data
      integer,allocatable,dimension(:) :: year,jday
      real,allocatable,dimension(:) :: hour,sp_wave1,sp_wave2
      real,allocatable,dimension(:,:) ::
     |  sp_flux   ! photon flux in photons/cm2/s1
      contains
!-----------------------------------------------------------------------
      subroutine init_seeflux
!
! Called from tgcm.F
! (this is not in init.F to avoid circular dependencies)
!
      use input_module,only: seeflux
!
! Read SEE flux data
      if (len_trim(seeflux) > 0) then
        write(6,"('Reading  SEEFLUX = ',a)") trim(seeflux)
        call rdseeflux
      endif
      end subroutine init_seeflux
!-----------------------------------------------------------------------
      subroutine rdseeflux
!
! Read TIMED-SEE EUV flux data
!
      use input_module,only: tempdir,seeflux
      use amie_module,only: rpt_ncerr
!     use netcdf_module,only:nc_open,nc_close,rpt_ncerr
      implicit none
#include "netcdf.inc"
!
! Local:
      character(len=80) :: dskfile
      integer :: istat,n_times,ncid,ndims,nvars,ngatts,idunlim,ier
      integer :: idv_wave1,idv_wave2,id_wave,id_time,nwaves,
     |    idv_year,idv_jday,idv_iter,idv_ut,
     |    idv_spflux
!
! Acquire mss file:
      dskfile = ' '
      call getfile(seeflux,dskfile)
!
      write(6,"(/,72('-'))")
      write(6,"('RDSEEFLUX: read SEE EUV flux data:')")
!
! Open netcdf file:
      istat = nf_open(dskfile,NF_NOWRITE,ncid)
      if (ncid<=0) then
        write(6,"(/,'>>> rdseeflux: error opening SEE flux data ',
     |    'diurnal file ',a)") trim(dskfile)
        stop 'rdseeflux'
      else
        write(6,"('rdseeflux: opened file ',a)") trim(dskfile)
      endif
!
! Get time dimension:
!     istat = nf_inq_unlimdim(ncid,id_time)
      istat = nf_inq_dimid(ncid,'dim1_YEAR',id_time)
      istat = nf_inq_dimlen(ncid,id_time,n_times)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting time dimension')
      istat = nf_inq_dimid(ncid,'dim1_WAVE1',id_wave)
      istat = nf_inq_dimlen(ncid,id_wave,nwaves)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting wave1 dimension')
      write(6,"('rdseeflux: n_times, nwaves =',2i7)") n_times,nwaves
!
! Search for requested output fields
      istat = nf_inq(ncid,ndims,nvars,ngatts,idunlim)
!
! Allocate 1-d fields:
      if (.not. allocated(year)) allocate(year(n_times),stat=ier)
      if (ier /= 0) write(6,"('>>> rdseeflux: error allocating',
     |  ' year: n_times=',i3)")n_times
      if (.not. allocated(jday)) allocate(jday(n_times),stat=ier)
      if (ier /= 0) write(6,"('>>> rdseeflux: error allocating',
     |  ' jday: n_times=',i3)")n_times
      if (.not. allocated(hour)) allocate(hour(n_times),stat=ier)
      if (ier /= 0) write(6,"('>>> rdseeflux: error allocating',
     |  ' hour: n_times=',i3)")n_times
      if (.not. allocated(sp_wave1)) allocate(sp_wave1(n_times),
     |    stat=ier)
      if (ier /= 0) write(6,"('>>> rdseeflux: error allocating',
     |  ' sp_wave1: n_times=',i3)")n_times
      if (.not. allocated(sp_wave2)) allocate(sp_wave2(n_times),
     |    stat=ier)
      if (ier /= 0) write(6,"('>>> rdseeflux: error allocating',
     |  ' sp_wave2: n_times=',i3)")n_times
      if (.not. allocated(sp_flux)) allocate(sp_flux(nwaves,n_times),
     |    stat=ier)
      if (ier /= 0) write(6,"('>>> rdseeflux: error allocating',
     |  ' sp_flux: n_times=',i3)")n_times
!
! Get 1-D fields (n_times)
      istat = nf_inq_varid(ncid,'YEAR',idv_year)
      istat = nf_get_var_int(ncid,idv_year,year)
      istat = nf_inq_varid(ncid,'JDAY',idv_jday)
      istat = nf_get_var_int(ncid,idv_jday,jday)
! Get ut
      istat = nf_inq_varid(ncid,'UT',idv_ut)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting UT id')
      istat = nf_get_var_double(ncid,idv_ut,hour)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting SPE variable ut')
! Get wave1, wave2, sp_flux
      istat = nf_inq_varid(ncid,'WAVE1',idv_wave1)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting sp_wave1 id')
      istat = nf_get_var_double(ncid,idv_wave1,sp_wave1)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting SPE variable sp_wave1')
      istat = nf_inq_varid(ncid,'WAVE2',idv_wave2)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting sp_wave2 id')
      istat = nf_get_var_double(ncid,idv_wave2,sp_wave2)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting SPE variable sp_wave2')
      istat = nf_inq_varid(ncid,'SP_FLUX',idv_spflux)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting sp_flux id')
      istat = nf_get_var_double(ncid,idv_spflux,sp_flux)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdseeflux: Error getting SPE variable sp_flux')
!
! Close the file:
      istat = nf_close(ncid)
      write(6,"('Completed read from SPE flux data file ',a)")
     |   trim(dskfile)
      write(6,"(72('-'),/)")

      end subroutine rdseeflux
!-----------------------------------------------------------------------
      subroutine get_seeflux(iyear,mday,secs,wave1,wave2,sflux)
!
      implicit none
! Args:
      integer, parameter :: lmax=37
      integer,intent(in) :: iyear,mday
      real,intent(in) :: secs
      real,intent(out) :: wave1(lmax),wave2(lmax),sflux(lmax)
!
! Locals:
      integer :: i,i0,i1,ii,msecs,dsecs0,dsecs1,nn
!
! msecs = current model time in seconds (includes day)
      msecs = mday*86400+int(secs)
!Using DAY 299 of 2003 for the background run
!     msecs = 299*86400+int(0)
      nn = size(hour)
! dsecs0,dsecs1 = beginning,ending data times (includes day)
      dsecs0 = jday(1)*86400+int(hour(1)*3600.)
      dsecs1 = jday(nn)  *86400+int(hour(nn)*3600.)
      if (msecs <= dsecs0 .or. msecs >= dsecs1) then
        write(6,"('>>> get_seeflux: cannot find SEE flux data for ',
     |    'model mday=',i4,' secs=',f10.2)") mday,secs
        write(6,"('>>> See seeflux_mod.F for available data.')")
        stop 'get_seeflux'
      endif
!
! Bracket model time:
      i0 = 0
      i1 = 0
       time_loop: do i=1,nn-1
        dsecs0 = jday(i)*86400+int(hour(i)*3600.)
        dsecs1 = jday(i+1)*86400+int(hour(i+1)*3600.)
        if (iyear==year(i) .and. msecs>=dsecs0 .and. msecs<=dsecs1) then
!         write(6,"('SEEFLUX>Found time: msecs=',i8,' dsecs0,1=',2i8)")
!    |      msecs,dsecs0,dsecs1
          i0 = i
          i1 = i+1
! Do simple linear interpolation:
          do ii=1,59
            sflux(ii) = 
     |        (msecs-dsecs0)*(sp_flux(ii,i1)-sp_flux(ii,i0))/
     |        (dsecs1-dsecs0)+sp_flux(ii,i0)
            wave1(ii) = sp_wave1(ii)
            wave2(ii) = sp_wave2(ii)
          enddo
          exit time_loop
        endif
      enddo  time_loop
!
! Stop if unable to bracket model time:
      if (i0==0.or.i1==0) then
        write(6,"(/,'>>> WARNING get_seeflux: could not',
     |    ' bracket model time year',i4,' jday=',i4,' secs=',i5)")
     |    iyear, mday,int(secs)
        stop 'get_seeflux'
      endif
 
!     write(6,"('seeflux: mday=',i3,' ut=',f5.2)")mday,secs/3600.
!     write(6,"('seeflux: see_flux= ',/(8e12.4))")(sflux(i),i=1,59)
!
      end subroutine get_seeflux
!-----------------------------------------------------------------------
      end module seeflux_module
