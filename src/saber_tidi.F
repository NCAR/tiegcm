      module saber_tidi
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Read and import netcdf files containing SABER (T,Z) and TIDI (U,V) data.
!
      use input_module,only: mxlen_filename
      use nchist_module,only:nc_open,nc_close,handle_ncerr
      use params_module,only: model_nlon=>nlon, model_nlat=>nlat
      use init_module,only: iyear,istep
      use hist_module,only: nstep
      implicit none
#include <netcdf.inc>
      private
      public saber_t,saber_z,tidi_u,tidi_v,get_saber_tidi
!
! Subdomains interpolated to current model time:
! (assigned from globals at every step)
!
      real,allocatable,save,dimension(:,:) :: ! (lon0:lon1,lat0:lat1)
     |  saber_t, saber_z, tidi_u, tidi_v
!
! Global data from model run start time to run end time (nlon,nlat,ntime):
! (one-time read at model init, see sub read_data_times).
!
      integer :: nlon,nlat
      integer :: ntime_saber,ntime_tidi
      real,allocatable,save,dimension(:,:,:) :: ! (nlon,nlat,ntime)
     |  saber_data_t, saber_data_z, tidi_data_u, tidi_data_v
      real,allocatable,save,dimension(:) ::     ! fractional year (ntime)
     |  saber_yfrac, tidi_yfrac
      real :: spval_t,spval_z,spval_u,spval_v   ! fill-value (missing data)
!
      contains
!-----------------------------------------------------------------------
      subroutine get_saber_tidi(modeltime)
!
! Driver for obtaining SABER and/or TIDI lbc data, called once per 
! timestep from advance.
!
      use input_module,only: saber_ncfile, tidi_ncfile, 
     |  start_day,calendar_advance
      use init_module,only: start_mtime,stop_mtime
!
! Args:
      integer,intent(in) :: modeltime(4)
!
! Local:
      integer :: mtime0(3),mtime1(3),mtime(3)
!
! Return silently if data files were not provided:
      if (len_trim(saber_ncfile) == 0 .and.
     |    len_trim(tidi_ncfile)  == 0) return
!
! Get model run start and stop times, and current modeltime:
      mtime0 = start_mtime    ! start model time
      mtime1 = stop_mtime     ! stop model time
      mtime = modeltime(1:3)  ! current model time
!
! For perpetual run, use data for the first day (24 hours), 
! then repeat that day of data throughout the perpetual run. 
! (use only hour and minute from the current model time)
!
      if (calendar_advance <= 0) then
        mtime0 = (/start_day,  0,0/)
        mtime1 = (/start_day+1,0,0/)
        mtime(1) = start_day
      endif

      if (istep==1)
     |  write(6,"('get_saber_tidi: calendar_advance=',i2,
     |    ' mtime0=',3i4,' mtime1=',3i4)") calendar_advance,
     |    mtime0,mtime1
!
! Saber data is requested:
      if (len_trim(saber_ncfile) > 0) then
!
! Read saber data file (once per run):
        if (istep==1) call read_saber(saber_ncfile,mtime0,mtime1)
!
! Interpolate saber data to current model time, and define subdomains:
        call interp_saber(mtime)
      endif
!
! Tidi data is requested:
      if (len_trim(tidi_ncfile) > 0) then
!
! Read tidi data file (once per run):
        if (istep==1) call read_tidi(tidi_ncfile,mtime0,mtime1)
!
! Interpolate tidi data to current model time, and define subdomains:
        call interp_tidi(mtime)
      endif
      end subroutine get_saber_tidi
!-----------------------------------------------------------------------
      subroutine read_data_times(ncid,dtype,mtime0,mtime1,istart,istop)
!
! Get indices to start and stop times in data corresponding to
! start and stop times of the model run. Called once per run from 
! read_saber and read_tidi. Also validate nlon,nlat dimensions.
!
! Args:
      integer,intent(in)  :: ncid,mtime0(3),mtime1(3)
      integer,intent(out) :: istart,istop
      character(len=*) :: dtype  ! 'SABER' or 'TIDI'
!
! Local:
      integer :: istat,i
      integer :: id_lon, id_lat, id_time     ! dimension ids
      integer :: idv_year, idv_day, idv_hour ! variable ids
      integer :: ndata,ntime
      real :: start_yfrac, stop_yfrac
      integer,allocatable,dimension(:) :: year,day,hour
      real,allocatable,dimension(:) :: yfrac
!
! Read and validate lat,lon dimensions:
      istat = nf_inq_dimid(ncid,'lon' ,id_lon)
      istat = nf_inq_dimlen(ncid,id_lon,nlon)  ! defines module data nlon
      if (nlon /= model_nlon) then
        write(6,"('>>> Fatal read_data_times: ',a,' nlon=',i5,
     |    ' model_nlon=',i5)") dtype,nlon,model_nlon
        call shutdown('read_data_times')
      endif

      istat = nf_inq_dimid(ncid,'lat' ,id_lat)
      istat = nf_inq_dimlen(ncid,id_lat,nlat)  ! defines module data nlat
      if (nlat /= model_nlat) then
        write(6,"('>>> Fatal read_data_times: nlat=',i5,
     |    ' model_nlon=',i5)") dtype,nlat,model_nlat
        call shutdown('read_data_times')
      endif
!
! Read number of times on data file:
      istat = nf_inq_dimid(ncid,'time',id_time)
      istat = nf_inq_dimlen(ncid,id_time,ndata)
!
! Allocate data date/time variables:
      allocate(year(ndata),stat=istat)
      allocate(day (ndata),stat=istat)
      allocate(hour(ndata),stat=istat)
      year = 0. ; hour = 0. ; day = 0.
      allocate(yfrac(ndata),stat=istat)
      yfrac = 0.
!
! Read data date/time variables, and calculate year-fraction:
      istat = nf_inq_varid(ncid,'year',idv_year)
      istat = nf_inq_varid(ncid,'day' ,idv_day)
      istat = nf_inq_varid(ncid,'hour',idv_hour)

      istat = nf_get_var_int(ncid,idv_year,year)
      istat = nf_get_var_int(ncid,idv_day ,day)
      istat = nf_get_var_int(ncid,idv_hour,hour)

      write(6,"('Reading ',a,' data file.')") dtype
      write(6,"('nlon=',i4,' nlat=',i4,' ndata=',i6)") nlon,nlat,ndata
      write(6,"('Hourly data starts at yyyy/ddd/hh =',i4,'/',i3,'/',
     |  i2)") year(1),day(1),hour(1)
      write(6,"('Hourly data ends   at yyyy/ddd/hh =',i4,'/',i3,'/',
     |  i2)") year(ndata),day(ndata),hour(ndata)
      write(6,"('Total number of data times = ',i6)") ndata

      do i=1,ndata
        if (mod(year(i),4) == 0) then
          yfrac(i) = year(i)+(day(i)-1+(hour(i)/24.))/366.
        else
          yfrac(i) = year(i)+(day(i)-1+(hour(i)/24.))/365.
        endif
      enddo
!
! Calculate fractional year for start and stop times:
      if (mod(iyear,4) == 0) then
        start_yfrac = iyear+(mtime0(1)-1+(mtime0(2)/24.)+
     |    (mtime0(3)/(24.*60.)))/366.
        stop_yfrac = iyear+(mtime1(1)-1+(mtime1(2)/24.)+
     |    (mtime1(3)/(24.*60.)))/366.
      else
        start_yfrac = iyear+(mtime0(1)-1+(mtime0(2)/24.)+
     |    (mtime0(3)/(24.*60.)))/365.
        stop_yfrac = iyear+(mtime1(1)-1+(mtime1(2)/24.)+
     |    (mtime1(3)/(24.*60.)))/365.
      endif
!
! Search for start time in data:
      istart = 0
      do i=1,ndata-1
        if (yfrac(i) <= start_yfrac .and. yfrac(i+1) >= start_yfrac)then
          istart = i
        endif
      enddo
      if (istart==0) then
        write(6,"('>>> ',a,' Could not bracket model start',
     |    ' time ',3i4,' start_yfrac=',f15.6,' dtype=',a)") 
     |  dtype,mtime0,start_yfrac
        call shutdown('read_data_times') 
      endif
!
! Search for stop time in data:
      istop = 0
      do i=1,ndata-1
        if (yfrac(i) <= stop_yfrac .and. yfrac(i+1) >= stop_yfrac)then
          istop = i+1
        endif
      enddo
      if (istop==0) then
        write(6,"('>>> ',a,' Could not bracket model stop time ',3i4,
     |    ' stop_yfrac=',f15.6)") dtype,mtime1,stop_yfrac 
        call shutdown('read_data_times') 
      endif
      if (istart > istop) then
        write(6,"('>>> read_data_times: dtype=',a,' istart must be ',
     |    '<= istop: istart=',i6,' istop=',i6)") 
     |    trim(dtype),istart,istop
        call shutdown('read_data_times')
      endif
      ntime = istop-istart+1 ! defines module data ntime

      write(6,"('Read ',a,' data from start mtime ',3i4,
     |  ' to stop mtime ',3i4,' ntimes=',i6)") 
     |  dtype,mtime0,mtime1,ntime
!
! Define module data year-fractions and ntimes:
      if (dtype=='SABER') then
        allocate(saber_yfrac(ntime),stat=istat)
        saber_yfrac(:) = yfrac(istart:istop)
        ntime_saber = ntime
      else
        allocate(tidi_yfrac(ntime),stat=istat)
        tidi_yfrac(:) = yfrac(istart:istop)
        ntime_tidi = ntime
      endif

      deallocate(year)
      deallocate(day)
      deallocate(hour)
      deallocate(yfrac)
 
      end subroutine read_data_times
!-----------------------------------------------------------------------
      subroutine read_saber(ncfile,mtime0,mtime1)
!
! Read saber data from model run start to end times,
! setting module data saber_data_x. This is called
! once per run from get_saber_tidi.
!
!
! Args:
      integer,intent(in) :: mtime0(3),mtime1(3)
      character(len=*),intent(in) :: ncfile
!
! Local:
      integer :: ncid,istat,i,ntime,istart_saber,istop_saber
      integer :: idv_t,idv_z
      integer :: start(3),count(3)
      character(len=mxlen_filename) :: dskfile
      character(len=80) :: units_t,units_z
!
      dskfile = ' '
      call getfile(ncfile,dskfile)
      write(6,"(/,72('-'))")
!
! Open netcdf file:
      call nc_open(ncid,dskfile,'OLD','READ')
      if (ncid==0) then
        write(6,"(/,'>>> read_saber: error opening netcdf saber ',
     |    'file ',a,' dskfile ',a)") trim(ncfile),trim(dskfile)
        call shutdown('read_saber')
      endif
!
! Get indices to saber data for model start,stop times:
! (istart_saber,istop_saber are returned)
!
      call read_data_times(ncid,'SABER',mtime0,mtime1,
     |  istart_saber,istop_saber)
      ntime = istop_saber-istart_saber+1
!
! Allocate global data for ntime times:
!
      allocate(saber_data_t(nlon,nlat,ntime),stat=istat)
      if (istat /= 0) call shutdown('allocate saber_data_t')
      saber_data_t = 0.
      allocate(saber_data_z(nlon,nlat,ntime),stat=istat)
      if (istat /= 0) call shutdown('allocate saber_data_z')
      saber_data_z = 0.

      start = (/1,1,istart_saber/)
      count = (/nlon,nlat,ntime/)
!
! Read T and check for missing data:
      istat = nf_inq_varid(ncid,'T',idv_t)
      istat = nf_get_vara_double(ncid,idv_t,start,count,saber_data_t) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_saber: Error getting variable T')
!
      istat = nf_get_att_double(ncid,idv_t,'missing_data',spval_t)
      if (any(saber_data_t==spval_t)) then
        write(6,"('>>> Missing data found in SABER T data.')")
        write(6,"('    Between file time indicies ',i6,' and ',i6)") 
     |    istart_saber,istop_saber
        call shutdown('Missing SABER T data')
      endif
!
! Check T units:
      units_t = ' '
      istat = nf_get_att_text(ncid,idv_t,"units",units_t)
      if (units_t(1:5) /= 'deg K') write(6,"('>>> read_saber: ',
     |  'Unknown units of T = ',a)") trim(units_t)
!
! Read Z and check for missing data:
      istat = nf_inq_varid(ncid,'Z',idv_z)
      istat = nf_get_vara_double(ncid,idv_z,start,count,saber_data_z)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_saber: Error getting variable Z')
!
      istat = nf_get_att_double(ncid,idv_z,'missing_data',spval_z)
      if (any(saber_data_z==spval_z)) then
        write(6,"('>>> Missing data found in SABER Z data.')")
        write(6,"('    Between file time indicies ',i6,' and ',i6)") 
     |    istart_saber,istop_saber
        call shutdown('Missing SABER Z data')
      endif
!
! Check Z units:
      units_z = ' '
      istat = nf_get_att_text(ncid,idv_z,"units",units_z)
      if (units_z(1:2) == 'km') then
        saber_data_z = saber_data_z * 1.e5  ! km to cm for model
        write(6,"('Converted saber_data_z units from km to cm')")
        units_z = 'cm'
      elseif (units_z(1:2) == 'cm') then
        write(6,"('>>> read_saber: Unexpected units of Z = ',
     |    a)") trim(units_z)
      endif

      write(6,"('T perturbations data min,max=',2e12.4,' (',a,')')") 
     |  minval(saber_data_t),maxval(saber_data_t),trim(units_t)
      write(6,"('Z perturbations data min,max=',2e12.4,' (',a,')')") 
     |  minval(saber_data_z),maxval(saber_data_z),trim(units_z)

      write(6,"(/,72('-'))")
      end subroutine read_saber
!-----------------------------------------------------------------------
      subroutine read_tidi(ncfile,mtime0,mtime1)
!
! Read tidi data from model run start to end times,
! setting module data tidi_data_x. This is called
! once per run from get_saber_tidi.
!
! Args:
      integer,intent(in) :: mtime0(3),mtime1(3)
      character(len=*),intent(in) :: ncfile
!
! Local:
      integer :: ncid,istat,i,ntime,istart_tidi,istop_tidi
      integer :: idv_u,idv_v,itime,j
      integer :: start(3),count(3)
      character(len=mxlen_filename) :: dskfile
      character(len=80) :: units_u,units_v
!
      dskfile = ' '
      call getfile(ncfile,dskfile)
      write(6,"(/,72('-'))")
!
! Open netcdf file:
      call nc_open(ncid,dskfile,'OLD','READ')
      if (ncid==0) then
        write(6,"(/,'>>> read_tidi: error opening netcdf tidi ',
     |    'file ',a,' dskfile ',a)") trim(ncfile),trim(dskfile)
        call shutdown('read_tidi')
      endif
!
! Get tidi data for model start,stop times:
! (istart_tidi,istop_tidi are returned)
!
      call read_data_times(ncid,'TIDI',mtime0,mtime1,
     |  istart_tidi,istop_tidi)
      ntime = istop_tidi-istart_tidi+1
!
! Allocate global data for ntime times:
!     real,allocatable,save,dimension(:,:,:) :: ! nlon,nlat,ndata
!    |  saber_data_t, saber_data_z, tidi_data_u, tidi_data_v
!
      allocate(tidi_data_u(nlon,nlat,ntime),stat=istat)
      if (istat /= 0) call shutdown('allocate tidi_data_u')
      tidi_data_u = 0.
      allocate(tidi_data_v(nlon,nlat,ntime),stat=istat)
      if (istat /= 0) call shutdown('allocate tidi_data_v')
      tidi_data_v = 0.

      start = (/1,1,istart_tidi/)
      count = (/nlon,nlat,ntime/)
!
! Read U:
      istat = nf_inq_varid(ncid,'U',idv_u)
      istat = nf_get_vara_double(ncid,idv_u,start,count,tidi_data_u) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_tidi: Error getting variable U')
!
      istat = nf_get_att_double(ncid,idv_u,'missing_data',spval_u)
      if (istat == NF_NOERR) then
        if (any(tidi_data_u==spval_u)) then
          write(6,"('>>> Missing data found in TIDI U data')")
          write(6,"('    Between file time indicies ',i6,' and ',i6)") 
     |      istart_tidi,istop_tidi
!         call shutdown('Missing TIDI U data')
        endif
      endif
!
! Check units of U:
      units_u = ' '
      istat = nf_get_att_text(ncid,idv_u,"units",units_u)
      if (units_u(1:3) == 'm/s') then
        tidi_data_u = tidi_data_u * 100.
        write(6,"('Converted units of U from m/s to cm/s')")
        units_u = 'cm/s'
      elseif (trim(units_u) /= 'cm/s') then
        write(6,"('>>> read_tidi: Unexpected units of U = ',
     |    a)") trim(units_u)
      endif
!
! Read V and check for missing data:
      istat = nf_inq_varid(ncid,'V',idv_v)
      istat = nf_get_vara_double(ncid,idv_v,start,count,tidi_data_v) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_tidi: Error getting variable V')
!
      istat = nf_get_att_double(ncid,idv_v,'missing_data',spval_v)
      if (istat == NF_NOERR) then
        if (any(tidi_data_v==spval_v)) then
          write(6,"('>>> Missing data found in TIDI V data')")
          write(6,"('    Between file time indicies ',i6,' and ',i6)") 
     |      istart_tidi,istop_tidi
          call shutdown('Missing TIDI V data')
        endif
      endif
!
! Check units of V:
      units_v = ' '
      istat = nf_get_att_text(ncid,idv_v,"units",units_v)
      if (units_v(1:3) == 'm/s') then
        tidi_data_v = tidi_data_v * 100.
        units_v = 'cm/s'
        write(6,"('Converted units of V from m/s to cm/s')")
      elseif (units_v(1:4) /= 'cm/s') then
        write(6,"('>>> read_tidi: Unexpected units of V = ',
     |    a)") trim(units_v)
      endif

      write(6,"('U perturbations data min,max=',2e12.4,' (',a,')')") 
     |  minval(tidi_data_u),maxval(tidi_data_u),trim(units_u)
      write(6,"('V perturbations data min,max=',2e12.4,' (',a,')')") 
     |  minval(tidi_data_v),maxval(tidi_data_v),trim(units_v)

      write(6,"(72('-'))")
      end subroutine read_tidi
!-----------------------------------------------------------------------
      subroutine interp_saber(mtime)
!
! Interpolate saber data to current model time, defining module data
! saber_t, saber_z (nlon,nlat) at subdomains.
!
      use mpi_module,only: lon0,lon1,lat0,lat1
!
! Args:
      integer,intent(in) :: mtime(3)
!
! Local:
      integer :: ndays,i,j,i0,i1,istat,ibeg,iend
      real :: yfrac ! current model year-fraction
      logical,external :: time2print
!
! Calculate model year-fraction:
      ndays = 365
      if (mod(iyear,4) == 0) ndays = 366
      yfrac = iyear+(mtime(1)-1+(mtime(2)/24.)+
     |    (mtime(3)/(24.*60.)))/float(ndays)
!
! Bracket model time between data i0,i1 times:
      i0 = 0 ; i1 = 0
      do i=1,ntime_saber-1
        if (yfrac >= saber_yfrac(i) .and. yfrac <= saber_yfrac(i+1))
     |    then
          i0 = i
          i1 = i+1
        endif
      enddo
      if (i0==0.and.i1==0) then
        write(6,"('>>> interp_saber: could not bracket model yfrac',
     |    ' =',f15.6,' mtime=',3i5)") yfrac,mtime
        call shutdown('interp_saber')
      endif
!
! Allocate subdomain data (1st step only):
      if (.not.allocated(saber_t)) then
        allocate(saber_t(lon0:lon1,lat0:lat1),stat=istat)
        if (istat /= 0) call shutdown('Error allocating saber_t')
      endif
      if (.not.allocated(saber_z)) then
        allocate(saber_z(lon0:lon1,lat0:lat1),stat=istat)
        if (istat /= 0) call shutdown('Error allocating saber_z')
      endif
!
! Do linear interpolation to model time.
! Note saber_data_t is global, saber_t is subdomain only.
! Subdomains saber_t(3-74,j) <= global saber_data_t(1-72,j)
!
      ibeg = lon0 ; if (lon0==1) ibeg = 3
      iend = lon1 ; if (lon1==nlon+4) iend = lon1-2 ! iend=74
      do j=lat0,lat1
        do i=ibeg,iend
          saber_t(i,j) = finterp1(saber_data_t(i-2,j,i0),
     |      saber_data_t(i-2,j,i1),saber_yfrac(i0),saber_yfrac(i1),
     |      yfrac)
          saber_z(i,j) = finterp1(saber_data_z(i-2,j,i0),
     |      saber_data_z(i-2,j,i1),saber_yfrac(i0),saber_yfrac(i1),
     |      yfrac)
        enddo ! i=lon0,lon1

!       write(6,"('before set_per: j=',i4,' ibeg,iend=',2i4,
!    |    ' saber_t(ibeg:iend,j)=',/,(6e12.4))") j,ibeg,iend,
!    |    saber_t(ibeg:iend,j)
      enddo ! j=lat0,lat1
!
! Set periodic points:
      call set_periodic(saber_t,lon0,lon1,lat0,lat1)
      call set_periodic(saber_z,lon0,lon1,lat0,lat1)

      if (time2print(nstep,istep))
     |  write(6,"('interp_saber: mtime=',3i4,' saber_t min,max=',
     |    2e12.4,' saber_z min,max=',2e12.4)") mtime,
     |    minval(saber_t),maxval(saber_t),
     |    minval(saber_z),maxval(saber_z)

      end subroutine  interp_saber
!-----------------------------------------------------------------------
      subroutine interp_tidi(mtime)
!
! Interpolate tidi data to current model time, defining module data
! tidi_u, tidi_v (nlon,nlat) at subdomains.
!
      use mpi_module,only: lon0,lon1,lat0,lat1
!
! Args:
      integer,intent(in) :: mtime(3)
!
! Local:
      integer :: ndays,i,j,i0,i1,istat,ibeg,iend
      real :: yfrac ! current model year-fraction
      logical,external :: time2print
!
! Calculate model year-fraction:
      ndays = 365
      if (mod(iyear,4) == 0) ndays = 366
      yfrac = iyear+(mtime(1)-1+(mtime(2)/24.)+
     |    (mtime(3)/(24.*60.)))/float(ndays)
!
! Bracket model time between data i0,i1 times:
      i0 = 0 ; i1 = 0
      do i=1,ntime_tidi-1
        if (yfrac >= tidi_yfrac(i) .and. yfrac <= tidi_yfrac(i+1))
     |    then
          i0 = i
          i1 = i+1
        endif
      enddo
      if (i0==0.and.i1==0) then
        write(6,"('>>> interp_tidi: could not bracket model yfrac',
     |    ' =',f15.6,' mtime=',3i5)") yfrac,mtime
        call shutdown('interp_tidi')
      endif
!
! Allocate subdomain data (1st step only):
      if (.not.allocated(tidi_u)) then
        allocate(tidi_u(lon0:lon1,lat0:lat1),stat=istat)
        if (istat /= 0) call shutdown('Error allocating tidi_u')
      endif
      if (.not.allocated(tidi_v)) then
        allocate(tidi_v(lon0:lon1,lat0:lat1),stat=istat)
        if (istat /= 0) call shutdown('Error allocating tidi_v')
      endif
!
! Do linear interpolation to model time:
! Note tidi_data_u is global, tidi_u is subdomain only.
! Subdomains tidi_u(3-74,j) <= global tidi_data_u(1-72,j)
!
      ibeg = lon0 ; if (lon0==1) ibeg = 3
      iend = lon1 ; if (lon1==nlon+4) iend = lon1-2
      do j=lat0,lat1
        do i=ibeg,iend
          tidi_u(i,j) = finterp1(tidi_data_u(i-2,j,i0),
     |      tidi_data_u(i-2,j,i1),tidi_yfrac(i0),tidi_yfrac(i1),
     |      yfrac)
          tidi_v(i,j) = finterp1(tidi_data_v(i-2,j,i0),
     |      tidi_data_v(i-2,j,i1),tidi_yfrac(i0),tidi_yfrac(i1),
     |      yfrac)
        enddo
      enddo
!
! Set periodic points:
      call set_periodic(tidi_u,lon0,lon1,lat0,lat1)
      call set_periodic(tidi_v,lon0,lon1,lat0,lat1)

      if (time2print(nstep,istep))
     |  write(6,"('interp_tidi:  mtime=',3i4,' tidi_u min,max =',
     |    2e12.4,' tidi_v min,max =',2e12.4)") mtime,
     |    minval(tidi_u),maxval(tidi_u),
     |    minval(tidi_v),maxval(tidi_v)

      end subroutine  interp_tidi
!-----------------------------------------------------------------------
      real function finterp1(f0,f1,frac0,frac1,frac)
!
! Args:
      real,intent(in) :: f0,f1,frac0,frac1,frac
!
      finterp1 = f0+(f1-f0)*(frac1-frac0)/(frac1-frac0)
      end function finterp1
!-----------------------------------------------------------------------
      subroutine set_periodic(f,lon0,lon1,lat0,lat1)
      use params_module,only: nlonp4
#ifdef MPI
      use mpi_module,only: mp_periodic_f2d
#endif
!
! lons 1,2 <- nlonp4-3,nlonp4-2 and nlonp4-1,nlonp4 <- 3,4
! lons 1,2 <- 72,73 and lons 75,76 <- 3,4
!
      integer,intent(in) :: lon0,lon1,lat0,lat1
      real,intent(inout) :: f(lon0:lon1,lat0:lat1)
#ifdef MPI
      call mp_periodic_f2d(f,lon0,lon1,lat0,lat1,1)
#else
      f(1:2,:)              = f(nlonp4-3:nlonp4-2,:)
      f(nlonp4-1:nlonp4,:)  = f(3:4,:)
#endif
      end subroutine set_periodic
!-----------------------------------------------------------------------
      end module saber_tidi

