!
      module nchist_module
      use params_module,only: nlon,nlev,nmlev,nlat,nlevp1,spval,nlonp4,
     |  nlonp1,nmlonp1,nmlat,nmlon,nmagphrlon,nmagphrlat
      use fields_module,only: nf4d,nf4d_hist,f4d,fakeflds,shortname_len,
     |  fsech,fsechmag,fsech2d,fsechmag2d,fsechmagphr2d,
     |  itc,itp,fields_4d,fields_3d,fields_2d,poten,foutput
      implicit none
!
! Module to create, read and write netcdf history files:
!
#ifdef SUN
#include "netcdf3.inc"
#else
#include "netcdf.inc"
#endif
!
      integer,private ::
     |  id_time,                   ! dimension id for time (unlimited)
     |  id_mtime,                  ! dimension id for model time (3)
     |  id_lon,id_lat,id_lev,      ! dimension ids for 3d model grid
     |  id_mlon,id_mlat,id_mlev,   ! dimension ids for 3d magnetic grid
     |  id_magphrlon,id_magphrlat, ! dimension ids for 2d magnetosphere grid
     |  id_latlon,                 ! dimension id for coords (2)
     |  id_dtide,                  ! dimension id for dtide (2)
     |  id_sdtide,                 ! dimension id for sdtide (10)
     |  id_fake,                   ! fake dimension id (1)
     |  ids1(1),ids2(2),ids3(3),ids4(4),	! vectors of dim id's
     |  start_2d(2),count_2d(2),	! start loc and count for 2d vars
     |  start_3d(3),count_3d(3),        ! start loc and count for 3d vars
     |  start_4d(4),count_4d(4),	! start loc and count for 4d vars
!
! Variable id's (idv_xxxx). Note lon, lat, lev, and time are
!   coordinate variables (lon(lon), lat(lat), lev(lev), time(time)).
! Variables starting at mtime are "primary" variables.
! idv_f4d(nf4d) are id's for prognostic 4d fields.
!
     |  idv_time,idv_mtime,idv_ut,      ! variable id's
     |  idv_lon,idv_lat,idv_lev,
     |  idv_mlon,idv_mlat,idv_mlev,
     |  idv_magphrlon,idv_magphrlat,
     |  idv_step,idv_iter,idv_year,idv_day,idv_p0,
     |  idv_hpower,idv_ctpoten,idv_byimf,idv_colfac,idv_ncep,
     |  idv_gpi,idv_gswmdi,idv_gswmsdi,idv_gswmnmidi,idv_gswmnmisdi,
     |  idv_f107d,idv_f107a,
     |  idv_mag,idv_dtide,idv_sdtide,idv_nflds,idv_f4d(nf4d),
     |  idv_alfa30,idv_e30,idv_ed2,idv_alfad2
      integer,allocatable :: idv_sech(:)
!
      contains
!-----------------------------------------------------------------------
      subroutine nc_open(ncid,flnm,status,rdwr)
!
! Open or create a netcdf dataset flnm, returning netcdf file
!   id handle in ncid.
! If status=='NEW' or 'REPLACE', a new netcdf file is created,
!   otherwise (status=='OLD'), open existing dataset).
! If an existing file is opened and rdwr=='READ', the dataset is
!   opened read-only (NF_NOWRITE), otherwise (rdwr=='WRITE' or 
!   'APPEND') the dataset is opened read/write (NF_WRITE).
! If the open fails, print error msg and return ncid == 0.
!
! Args:
      integer,intent(out) :: ncid    ! netcdf file id
      character(len=*),intent(in) :: flnm,status,rdwr
!
! Local:
      integer :: istat
      character(len=120) :: char120
!
      ncid = 0
      if (status=='NEW'.or.status=='REPLACE') then
!
! Create new dataset (clobber old one if it pre-exists):
!
        istat = nf_create(flnm,NF_CLOBBER,ncid)
        if (istat /= NF_NOERR) then
          write(char120,"('Error return from nf_create for netcdf ',
     |      'file ',a,' status=',a)") trim(flnm),status
          call handle_ncerr(istat,char120)
          ncid = 0
        else
          write(6,"('Created netcdf file ',a,' ncid=',i8)") 
     |      trim(flnm),ncid
        endif
!
! Open pre-existing dataset:
! NF_WRITE opens with read-write, where write means append, change, etc
! NF_NOWRITE is read-only.
!
      elseif (status=='OLD') then
        if (rdwr=='WRITE'.or.rdwr=='APPEND') then
          istat = nf_open(flnm,NF_WRITE,ncid)
        elseif (rdwr=='READ') then
          istat = nf_open(flnm,NF_NOWRITE,ncid)
        else
          write(6,"('>>> nc_open: unrecognized rdwr=',a)") trim(rdwr)
          write(6,"('  Will open NF_WRITE (read/write).')")
          istat = nf_open(flnm,NF_WRITE,ncid)
        endif
        if (istat /= NF_NOERR) then
          write(char120,"('Error return from nf_open for netcdf ',
     |      'file ',a,' status=',a)") trim(flnm),status
          call handle_ncerr(istat,char120)
          ncid = 0
        else
!         write(6,"('Opened existing netcdf file ',a,' ncid=',i8)") 
!    |      trim(flnm),ncid
        endif
!
! Unknown status request:
      else
        write(6,"(/,'>>> nc_open: unrecognized status = ',a)")
     |    status
        ncid = 0
      endif
      end subroutine nc_open
!-------------------------------------------------------------------
      subroutine nc_close(ncid)
!
! Close a netcdf dataset file:
!
! Arg:
      integer,intent(in) :: ncid
!
! Local:
      integer :: istat
!
      istat = nf_close(ncid)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error return from nf_close')
      end subroutine nc_close
!-----------------------------------------------------------------------
      subroutine nc_rdhist(ncid,diskfile,mtime,itime,ier)
!
! Read history at model time mtime from netcdf file diskfile,
!   currently open under ncid. Define output history in global
!   history structure h (module hist_module in history.f).
!   Prognostic fields are transferred to global fg-array (wrf4d).
! If there is a problem (e.g., mtime not found), return ier > 0
!
      use hist_module,only: h,hist_initype,hist_print,nsource,
     |  nsecsource
      use init_module,only: istep,start_mtime,glon,glat,plev,
     |  gmlon,gmlat,pmlev
      use input_module,only: start_day,start_year
!
! Args:
      integer,intent(in) :: ncid,mtime(3)
      character(len=*) :: diskfile
      integer,intent(out) :: itime,ier
!
! Local:
      integer,parameter :: mxdims=20 ! assume <= mxdims dimensions
      integer :: istat,id_time,id_mtime,mtime_read(3),ntimes,
     |  i,ix,istart2(2),icount2(2),j,irdf4d(nf4d_hist)
      integer :: ndims,nvars,ngatts,idunlim,itype,natts,len,
     |  iddims(mxdims),iscalar
      character(len=80) :: varname
      real :: scalar,var1(1),var22(2,2),var2(2),var10(10)
!
      write(6,"('Reading source history from diskfile ',a,':')")
     |  trim(diskfile)
      ier = 0
!
! Get id for mtime variable:
      istat = nf_inq_varid(ncid,"mtime",idv_mtime) 
      if (istat /= NF_NOERR) then
        call handle_ncerr(istat,'Error getting id of variable mtime')
        ier = 1
        return
      endif
!
! Get number of histories (times) on the file:
      istat = nf_inq_unlimdim(ncid,id_time) ! id of unlimited record var
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error getting unlimited dimension id')
      istat = nf_inq_dimlen(ncid,id_time,ntimes) ! length of time var
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error getting length of time record variable')
!
! Search for requested model time. When found, retain time index itime:
      istart2(1) = 1
      icount2(1) = 3
      icount2(2) = 1
      itime = 0
      do i=1,ntimes
        istart2(2) = i
        istat = nf_get_vara_int(ncid,idv_mtime,istart2,icount2,
     |    mtime_read)
        if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error getting mtime values')
        write(6,"('nc_rdhist: seeking ',3i4,' found ',3i4,'  n=',i3)")
     |    mtime,mtime_read,i
        if (all(mtime_read==mtime)) then
          itime = i
          exit
        endif
      enddo
      if (itime==0) then
        write(6,"(/,'>>> nc_rdhist: mtime ',3i4,' NOT found on ',
     |    'diskfile ',a)") mtime,trim(diskfile)
        ier = 1
        return
      endif
!
! Define history structure from netcdf vars:
!
! Initialize history structure:
      call hist_initype(h,istep)  ! initialize history structure
!
! Define vars already found above:
      h%ihist = itime
      h%modeltime(1:3) = mtime
      h%modeltime(4) = 0
!
! Read global attributes:
      istat = nf_get_att_text(ncid,NF_GLOBAL,"rundate",h%rundate)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"logname",h%logname)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"host",h%host)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"system",h%system)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"model_name",h%model_name)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"model_version",
     |  h%model_version)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"mss_path",h%mss_path)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"mss_source",h%mss_source)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"mss_secsource",
     |  h%mss_secsource)
!
! If this is a continuation run, read initial time from the history.
! Otherwise (this is an initial run) the initial time was set to the 
! model start time by hist_initype.
!
      if (nsource <= 0) then
        istat=nf_get_att_int(ncid,NF_GLOBAL,"initial_year",
     |    h%initial_year)
        if (istat /= NF_NOERR) istat=nf_get_att_int(ncid,NF_GLOBAL,
     |    "start_year" ,h%initial_year)

        istat=nf_get_att_int(ncid,NF_GLOBAL,"initial_day",h%initial_day)
        if (istat /= NF_NOERR) istat=nf_get_att_int(ncid,NF_GLOBAL,
     |    "start_day" ,h%initial_day)

        istat=nf_get_att_int(ncid,NF_GLOBAL,"initial_mtime",
     |    h%initial_mtime)
        if (istat /= NF_NOERR) istat=nf_get_att_int(ncid,NF_GLOBAL,
     |    "start_mtime",h%initial_mtime)

        istat=nf_get_att_int(ncid,NF_GLOBAL,"source_mtime",
     |    h%source_mtime)
      endif
!
! History_type is either "primary" or "secondary":
      istat = nf_get_att_text(ncid,NF_GLOBAL,"history_type",h%type)
!
! Get number of dims, vars, atts, and id of unlimited dimension:
      istat = nf_inq(ncid,ndims,nvars,ngatts,idunlim)
!
! Get dimensions and check against current model parameters:
! (i.e., current model run at resolution A should not read a 
!  source history from a model that was run with resolution B) 
!
      do i=1,ndims
        istat = nf_inq_dim(ncid,i,varname,len)
        if (trim(varname)=='lon') then
          if (len /= nlon) then
            write(6,"(/,'>>> nc_rdhist: length of lon read from',
     |        ' history does not match parameter nlon.')")
            write(6,"('  len=',i4,' nlon=',i4)") len,nlon
            stop 'nlon'
          endif
          h%nlon = len
        endif
        if (trim(varname)=='lat') then
          if (len /= nlat) then
            write(6,"(/,'>>> nc_rdhist: length of lat read from',
     |        ' history does not match parameter nlat.')")
            write(6,"('  len=',i4,' nlat=',i4)") len,nlat
            stop 'nlat'
          endif
          h%nlat = len
        endif
        if (trim(varname)=='lev') then
!
! If number of levels read from source history (len) is less than the
! number of levels allocated in the model (nlev), this is fatal:
!
          if (len < nlevp1) then
            write(6,"(/,'>>> nc_rdhist: length of lev dimension',
     |        ' of history read is less than model parameter nlevp1.')")
            write(6,"('  len=',i4,' nlevp1=',i4)") len,nlevp1
            stop 'nlevp1'
!
! If number of levels read from source history (len) is greater than
! the number of levels allocated in the model (nlevp1), then we will
! read the first nlevp1 levels (e.g., if nlevp1==88 and len==96, then
! read levels 1-88 from the history i.e., the bottom 88 levels).
!
          elseif (len > nlevp1) then
            write(6,"(/,'NOTE nc_rdhist: number of levels in source',
     |        ' history (len=',i3,') exceeds the ',/,2x,
     |        'number of levels allocated in the model (nlevp1=',
     |        i3,')')") len,nlevp1
            write(6,"('  Will read the first ',i3,' levels from the',
     |        ' source into the history buffer.')") nlevp1
          endif
          h%nlev = len
        endif
      enddo
!
! Get variables (should not need to get coord vars lat,lon,lev)
      do i=1,nvars
        istat = nf_inq_var(ncid,i,varname,itype,ndims,iddims,natts)
        select case(trim(varname))
          case('lon')
            istat = nc_get_var_real(ncid,i,glon)
          case('lat') 
            istat = nc_get_var_real(ncid,i,glat)
          case('lev')
            istat = nc_get_var_real(ncid,i,plev)
          case('mlon')
            istat = nc_get_var_real(ncid,i,gmlon)
          case('mlat') 
            istat = nc_get_var_real(ncid,i,gmlat)
          case('mlev')
            istat = nc_get_var_real(ncid,i,pmlev)
          case('time')
            istat = nc_get_var1_real(ncid,i,itime,scalar)
            h%time = scalar
          case('mtime') ! already got it from hist search above
          case('year')
            istat = nf_get_var1_int(ncid,i,itime,iscalar)
            h%year = iscalar
          case('day')
            istat = nf_get_var1_int(ncid,i,itime,iscalar)
            h%day = iscalar
          case('iter')
            istat = nf_get_var1_int(ncid,i,itime,iscalar)
            h%iter = iscalar
          case('ut')
            istat = nc_get_var1_real(ncid,i,itime,scalar)
            h%ut = scalar
          case('mag')
            istat = nc_get_var_real(ncid,i,var22)
            h%mag(:,:) = var22(:,:)
          case('dtide')
            istart2(1) = 1
            istart2(2) = h%ihist
            icount2(1) = 2 
            icount2(2) = 1
            istat = nf_get_vara_double(ncid,i,istart2,icount2,var2)
            h%dtide(:) = var2(:)
          case('sdtide')
            istart2(1) = 1
            istart2(2) = h%ihist
            icount2(1) = 10
            icount2(2) = 1
            istat = nf_get_vara_double(ncid,i,istart2,icount2,var10)
            h%sdtide(:) = var10(:)
          case('f107d')
            istat = nc_get_var_real(ncid,i,var1)
            h%f107d = var1(1)
          case('f107a')
            istat = nc_get_var_real(ncid,i,var1)
            h%f107a = var1(1)
          case('hpower')
            istat = nc_get_var_real(ncid,i,var1)
            h%hpower = var1(1)
          case('ctpoten')
            istat = nc_get_var_real(ncid,i,var1)
            h%ctpoten = var1(1)
          case('byimf')
            istat = nc_get_var_real(ncid,i,var1)
            h%byimf = var1(1)
          case('colfac')
            istat = nc_get_var_real(ncid,i,var1)
            h%colfac = var1(1)
          case('alfa30')
            istat = nc_get_var_real(ncid,i,var1)
            h%alfa30 = var1(1)
          case('e30')
            istat = nc_get_var_real(ncid,i,var1)
            h%e30 = var1(1)
          case('ed2')
            istat = nc_get_var_real(ncid,i,var1)
            h%ed2 = var1(1)
          case('alfad2')
            istat = nc_get_var_real(ncid,i,var1)
            h%alfad2 = var1(1)
          case('p0')
            istat = nc_get_var_real(ncid,i,var1)
            h%p0 = var1(1)
          case('timestep')
            istat = nf_get_var_int(ncid,i,iscalar)
            h%step = iscalar
          case('nflds')
            istat = nf_get_var_int(ncid,i,iscalar)
            h%nflds = iscalar
          case('ncep')
            istat = nf_get_var_int(ncid,i,iscalar)
            h%ncep = iscalar
          case('gpi')
            istat = nf_get_var_int(ncid,i,iscalar)
            h%gpi = iscalar
          case('gswmdi')
            istat = nf_get_var_int(ncid,i,iscalar)
            h%gswmdi = iscalar
          case('gswmsdi')
            istat = nf_get_var_int(ncid,i,iscalar)
            h%gswmsdi = iscalar
          case('gswmnmidi')
            istat = nf_get_var_int(ncid,i,iscalar)
            h%gswmnmidi = iscalar
          case('gswmnmisdi')
            istat = nf_get_var_int(ncid,i,iscalar)
            h%gswmnmisdi = iscalar
          case('maglon')

          case('maglat')

          case default
            if (any(f4d(1:nf4d_hist)%short_name==varname)) then
              do j=1,nf4d_hist
                if (trim(varname)==f4d(j)%short_name) then
                  irdf4d(j) = 1
                  call rdf4d(ncid,varname,itime,j)
                endif
              enddo
            else
              write(6,"('Note nc_rdhist: unknown variable: ',a)") 
     |          trim(varname)
            endif
        end select
      enddo ! i=1,nvars
      h%zptop = plev(nlevp1)
      h%zpbot = plev(1)
!
! Report to stdout:
      call hist_print(h,'READ')
!     
      end subroutine nc_rdhist
!-----------------------------------------------------------------------
      subroutine rdf4d(ncid,name,itime,ix)
!
! Read prognostic field "name" from netcdf file on ncid into into local
! f3drd at full domain. Then subdomains are transferred from f3rd into
! 4-d arrays for use in the model. It would be interesting to try reading
! only the subdomain directly from the file rather than the full domain.
!
      use input_module,only: dynamo
      use hist_module,only: h
      use init_module,only: glat
      use mpi_module,only: lon0,lon1,lat0,lat1
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: ncid,itime,ix
!
! Local:
      integer :: i,ii,j,k,istat,idvar,itype,ndims,iddims(4),natts,
     |  iscalar,idimsizes(4),nx,nxk,lonbeg,lonend
      character(len=16) :: rdname,tmpname
      character(len=120) :: char120
      real :: fakevar(1,1,1),fmin,fmax,fmin1,fmax1
      real :: f3drd(nlon,nlat,nlevp1) ! (i,j,k) as on history 
      real :: zmj2,zmj1,fac
!
! Fields whose lower boundary condition is stored in the top
! pressure slot (k==nlev==nlevp1) (see dt.F and duv.F in model):
      character(len=shortname_len) :: names_lbc_intop(6) = 
     |  (/'TN              ','UN              ','VN              ',
     |    'TN_NM           ','UN_NM           ','VN_NM           '/)
!
! External:
      integer,external :: strloc
!
! Get info about the field:
      istat = nf_inq_varid(ncid,name,idvar)
      istat = nf_inq_var(ncid,idvar,rdname,itype,ndims,iddims,natts)
      if (ndims /= 4) then
        write(6,"(/,'>>> WARNING rdf4d: bad ndims=',i3,
     |    ' (every prognostic should have 4 dimensions)')") ndims
      endif
!
! Get info about the dimensions:
      do i=1,ndims
        istat = nf_inq_dim(ncid,iddims(i),tmpname,idimsizes(i))
      enddo 
!
! Check dimension sizes against current model grid resolution:
      if (.not.fakeflds) then
        if (idimsizes(1) /= nlon .or. idimsizes(2) /= nlat) then
          write(6,"(/,'>>> WARNING rdf4d: bad dimension sizes',
     |      ' for prognostic ',a)") trim(rdname)
          write(6,"('  idimsizes(1:2)=',2i4,' but should be ',
     |      'nlon,nlat=',2i4)") idimsizes(1:2),nlon,nlat
        endif
!
! Error if number of levels in source history is less than number
! of levels in model run (not an error if number of levels in source
! history is greater than nlevp1).
!
        if (idimsizes(3) < nlevp1) then
          write(6,"(/,'>>> WARNING rdf4d: bad idimsizes(3)=',
     |      i4,' nlevp1=',i4,' field ',a)") idimsizes(3),nlevp1,
     |      trim(rdname)
        endif
      endif
      if (fakeflds) then
        start_4d(1:3) = 1
        start_4d(4) = itime
        count_4d(1) = 1
        count_4d(2) = 1
        count_4d(3) = 1
        count_4d(4) = 1
        istat = nf_get_vara_double(ncid,idvar,start_4d,count_4d,fakevar)
        if (istat /= NF_NOERR) then
          write(char120,"('Error return from nf_get_vara_double',
     +      ' for fake field var ',a,' itime=',i2)") name,itime
          call handle_ncerr(istat,char120)
        endif
        return
      endif
!
! Zero out electric potential and return if DYNAMO flag is not set:
      if (trim(name)=='POTEN'.and.dynamo <= 0) then
        f4d(ix)%data = 0.
        return
      endif
!
! Read field into data component of field structure:
      start_4d(1:3) = 1
      start_4d(4) = itime
      count_4d(1) = nlon
      count_4d(2) = nlat
      count_4d(3) = nlevp1
      count_4d(4) = 1
      istat = nf_get_vara_double(ncid,idvar,start_4d,count_4d,f3drd)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_get_vara_double',
     +    ' for field var ',a,' itime=',i2)") name,itime
        call handle_ncerr(istat,char120)
      endif
!
! Transform from f3drd(i,j,k) to f4d(k,i,j):
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2

!     call fminmax(f3drd(lonbeg-2:lonend-2,lat0:lat1,:),
!    |  nlevp1*((lonend-2)-(lonbeg-2)+1)*(lat1-lat0+1),fmin,fmax)
!     write(6,"('rdf4d: lonbeg,end=',2i3,' lat0,1=',2i3,
!    |  ' f3drd min,max=',2e12.4)") lonbeg,lonend,lat0,lat1,fmin,fmax

      ii = 1
      do j=lat0,lat1
        do i=lonbeg,lonend
          do k=1,nlevp1
!
! Local : real :: f3drd(nlon,nlat,nlevp1) ! (i,j,k) as on history 
! Global: real,pointer :: f4d(ix)%data(levd0:levd1,lond0:lond1,latd0:latd1) 
!
            f4d(ix)%data(k,i,j,itp) = f3drd(i-2,j,k)
          enddo ! k=1,nlevp1
          ii = ii+1
        enddo ! i=lonbeg,lonend  
      enddo ! j=lat0,lat1

      call fminmax(f4d(ix)%data(:,lonbeg:lonend,lat0:lat1,itp),
     |  nlevp1*(lonend-lonbeg+1)*(lat1-lat0+1),fmin,fmax)

!     write(6,"('lonbeg,lonend=',2i3)") lonbeg,lonend
      write(6,"('Read field ',a,' 3d min,max=',2e12.4)") 
     |  f4d(ix)%short_name,fmin,fmax
!
! Do periodic points later, when all fields can be loaded into
! mpi messages, e.g., after call nc_rdhist, see rdsource.F
!
! Periodic points for data(k,i,j):
!      f4d(ix)%data(:,1     ,1:nlat,itp) =   ! 1 <- nlon+1 (73)
!    | f4d(ix)%data(:,nlon+1,1:nlat,itp)  
!      f4d(ix)%data(:,2     ,1:nlat,itp) =   ! 2 <- nlon+2 (74)
!    | f4d(ix)%data(:,nlon+2,1:nlat,itp)
!      f4d(ix)%data(:,nlon+3,1:nlat,itp) =   ! nlon+3 (75) <- 3
!    | f4d(ix)%data(:,3     ,1:nlat,itp)
!      f4d(ix)%data(:,nlon+4,1:nlat,itp) =   ! nlon+4 (76) <- 4
!    | f4d(ix)%data(:,4     ,1:nlat,itp)
!

!
! If number of levels in source history (h%nlev) is > number of
! levels in model (nlevp1), then store t,u,v bottom boundaries in
! top (nlevp1) slot. This is important esp for u,v, so that the
! model tuvbnd.F will calculate u,v bottom boundaries correctly.
!
      if (h%nlev > nlevp1 .and. any(names_lbc_intop==name)) then
        write(6,"('rdf4d: nlevp1=',i3,' h%nlev=',i3,' field=',
     |    a,' -- storing bottom boundary in top level.')")
     |    nlevp1,h%nlev,trim(name)
! data(i,k,j):
!       f4d(ix)%data(nlevp1,:,1:nlat,itp) = 
!    |  f4d(ix)%data(1     ,:,1:nlat,itp)
! data(k,i,j):
        f4d(ix)%data(nlevp1,:,1:nlat,itp) = 
     |  f4d(ix)%data(1     ,:,1:nlat,itp)
      endif
!
! New models (tiegcmxx or timegcmxx) do not store bottom boundary
! of t,u,v in top slot. If reading an "old" model history (tgcmxx),
! put top slot into bottom boundary: 
!     if (h%model_version(1:4) == 'tgcm'.and.any(names_lbc_intop==name)) 
!    |  then
!!       f4d(ix)%data(:,1     ,1:nlat,itp) = 
!!    |  f4d(ix)%data(:,nlevp1,1:nlat,itp)
!       f4d(ix)%data(1     ,:,1:nlat,itp) = 
!    |  f4d(ix)%data(nlevp1,:,1:nlat,itp)
!     endif
!
! In old cray-blocked histories, NOPNM was not defined at the
! top level. Define it here with spval for compatability with
! old sources (s.a., read_oldsrc). This is not strictly necessary,
! but is convenient for debug and comparisons between old and new 
! histories with fminmax, etc.
!
!     if (trim(name)=='OP_NM') then
!       f4d(ix)%data(nlevp1,:,1:nlat,itp) = spval
!     endif

!     call fminmax(f4d(ix)%data(:,lon0:lon1,lat0:lat1,itp),
!    |  nlevp1*(lon1-lon0+1)*(lat1-lat0+1),fmin,fmax)
!     write(6,"('Read field ',a,' 3d min,max=',2e12.4)") 
!    |  f4d(ix)%short_name,fmin,fmax
!
! In earlier versions, full-domain dynpot was defined here from poten.
! In this version, subdomains of poten are read from the source history
! by this routine, then mp_dynpot is called (from readsource) to define 
! dynpot at the full domain for each task.
!
      end subroutine rdf4d
!-----------------------------------------------------------------------
      subroutine nc_define(ncid)
!
! Define dimensions, variables and attributes in a new history
!   netcdf dataset (file has been opened, but no histories written yet).
! On input:
!   ncid = file id for open netcdf dataset (in define mode).
!   h = defined history structure (in hist_module)
! On output:
!   Variables and dimensions for the netcdf file are defined.
!   Values are given to the dimension variables, and the file
!     is taken out of define mode (ready to receive histories)
!
      use hist_module,only: h,nsource,nsecsource
      use init_module,only: glon,glat,plev,gmlon,gmlat,pmlev
      use cons_module,only: pi,ylatmagphr,ylonmagphr
!
! Args:
      integer,intent(in) :: ncid
!
! Local:
      integer :: i,istat,idum,ivar1(1),imo,ida,startmtime(4)
      character(len=120) :: char120
      character(len=80) :: char80
      real :: var1(1),rmins
      real,external :: mtime_to_datestr
!
! Define dimensions (field vars will be (time,lev,lat,lon) in
! CDL notation. Mag fields will be (time,mlev,mlat,mlon): 
!
      istat = nf_def_dim(ncid,"time",NF_UNLIMITED,id_time)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining time dimension')
!
! Geographic grid dimensions:
      istat = nf_def_dim(ncid,"lon",nlon,id_lon)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining longitude dimension')
      istat = nf_def_dim(ncid,"lat",nlat,id_lat)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error defining latitude dimension')
!     istat = nf_def_dim(ncid,"lev",nlev,id_lev)
      istat = nf_def_dim(ncid,"lev",nlevp1,id_lev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining levels dimension')
!
! Magnetic grid dimensions:
      istat = nf_def_dim(ncid,"mlon",nmlonp1,id_mlon)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining magnetic longitude dimension')
      istat = nf_def_dim(ncid,"mlat",nmlat,id_mlat)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error defining magnetic latitude dimension')
      istat = nf_def_dim(ncid,"mlev",nlevp1+3,id_mlev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining magnetic levels dimension')
!
! am_09/02 magnetosphere grid dimensions:
      istat = nf_def_dim(ncid,"magphrlon",nmagphrlon,id_magphrlon)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining magnetosphere longitude dimension')
      istat = nf_def_dim(ncid,"magphrlat",nmagphrlat,id_magphrlat)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error defining magnetosphere latitude dimension')
!
      istat = nf_def_dim(ncid,"mtime",3,id_mtime)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining mtime dimension')
      istat = nf_def_dim(ncid,"latlon",2,id_latlon)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining latlon dimension')
      istat = nf_def_dim(ncid,"dtidedim",2,id_dtide)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining dtide dimension')
      istat = nf_def_dim(ncid,"sdtidedim",10,id_sdtide)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining sdtide dimension')
      if (fakeflds) then ! fake dimension for testing
        istat = nf_def_dim(ncid,"fakedim",1,id_fake)
        if (istat /= NF_NOERR) call handle_ncerr(istat,
     +    'Error defining fakedim dimension')
      endif
!
! Define dimension (coordinate) variables and attributes:
!
! Geographic longitude (deg) (coordinate variable lon(lon)):
      ids1(1) = id_lon
      istat = nf_def_var(ncid,"lon",NF_DOUBLE,1,ids1,idv_lon)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining longitude dimension variable')
      write(char80,"('geographic longitude (-west, +east)')")
      istat = nf_put_att_text(ncid,idv_lon,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of longitude dimension variable')
      istat = nf_put_att_text(ncid,idv_lon,"units",12,'degrees_east')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of longitude dimension variable')
!
! Geographic latitude (deg) (coordinate variable lat(lat)):
      ids1(1) = id_lat
      istat = nf_def_var(ncid,"lat",NF_DOUBLE,1,ids1,idv_lat)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining latitude dimension variable')
      write(char80,"('geographic latitude (-south, +north)')")
      istat = nf_put_att_text(ncid,idv_lat,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of latitude dimension variable')
      istat = nf_put_att_text(ncid,idv_lat,"units",13,'degrees_north')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of latitude dimension variable')
!
! Vertical levels (log pressure) (coordinate variable lev(lev)):
      ids1(1) = id_lev
      istat = nf_def_var(ncid,"lev",NF_DOUBLE,1,ids1,idv_lev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining levels dimension variable')
      write(char80,"('log pressure levels')")
      istat = nf_put_att_text(ncid,idv_lev,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of levels dimension variable')
      istat = nf_put_att_text(ncid,idv_lev,"units",8,'ln(p0/p)')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of levels dimension variable')
!
! Magnetic longitude (deg) (coordinate variable mlon(mlon)):
      ids1(1) = id_mlon
      istat = nf_def_var(ncid,"mlon",NF_DOUBLE,1,ids1,idv_mlon)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining magnetic longitude dimension variable')
      write(char80,"('magnetic longitude (-west, +east)')")
      istat = nf_put_att_text(ncid,idv_mlon,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of mag longitude dimension variable')
      istat = nf_put_att_text(ncid,idv_lon,"units",12,'degrees_east')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of magnetic longitude dimension variable')
!
! Magnetic latitude (deg) (coordinate variable mlat(mlat)):
      ids1(1) = id_mlat
      istat = nf_def_var(ncid,"mlat",NF_DOUBLE,1,ids1,idv_mlat)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining magnetic latitude dimension variable')
      write(char80,"('magnetic latitude (-south, +north)')")
      istat = nf_put_att_text(ncid,idv_mlat,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of mag latitude dimension variable')
      istat = nf_put_att_text(ncid,idv_mlat,"units",13,'degrees_north')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of magnetic latitude dimension variable')
!
! Magnetic levels (log pressure) (coordinate variable mlev(mlev)):
      ids1(1) = id_mlev
      istat = nf_def_var(ncid,"mlev",NF_DOUBLE,1,ids1,idv_mlev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining magnetic levels dimension variable')
      write(char80,"('log pressure levels (magnetic grid)')")
      istat = nf_put_att_text(ncid,idv_mlev,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of mag levels dimension variable')
      istat = nf_put_att_text(ncid,idv_mlev,"units",8,'ln(p0/p)')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of magnetic levels dimension variable')
!
! Magnetospheric longitude (deg) (coordinate variable magphrlon(magphrlon)):
      ids1(1) = id_magphrlon
      istat = nf_def_var(ncid,"maglon",NF_DOUBLE,1,ids1,idv_magphrlon)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining magnetospheric longitude dimension variable')
      write(char80,"('magnetospheric longitude (-west, +east)')")
      istat = nf_put_att_text(ncid,idv_magphrlon,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     + 'Error defining long_name of magnetospheric long. dim. var.')
      istat = nf_put_att_text(ncid,idv_magphrlon,"units",7,'degrees')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     + 'Error defining units of magnetospheric longitude dim. var.')
!
! Magnetospheric latitude (deg) (coordinate variable magphrlat(magphrlat)):
      ids1(1) = id_magphrlat
      istat = nf_def_var(ncid,"magphrlat",NF_DOUBLE,1,ids1,
     |   idv_magphrlat)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining magnetospheric latitude dimension variable')
      write(char80,"('magnetospheric latitude (-south, +north)')")
      istat = nf_put_att_text(ncid,idv_magphrlat,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     + 'Error defining long_name of magnetospheric latitude dim. var.')
      istat = nf_put_att_text(ncid,idv_magphrlat,"units",7,'degrees')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of magnetospheric latitude dim. var.')
!
! ids1(1) is id of time dimension for the following defines:
      ids1(1) = id_time	! for 1d time-unlimited vars
!
! Time (coordinate variable time(time)). This is minutes since
! the run's source history start time. The units string is: yyyy-m-d, 
! where yyyy is the year, m is month, and d is day of the source
! start time.
!
      istat = nf_def_var(ncid,"time",NF_DOUBLE,1,ids1,idv_time)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining time dimension variable')
      startmtime(1) = h%initial_day
      startmtime(2:3) = h%initial_mtime(2:3) ; startmtime(4) = 0
      rmins = mtime_to_datestr(h%initial_year,startmtime,imo,ida,char80)
      istat = nf_put_att_text(ncid,idv_time,"long_name",4,"time")
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of time dimension variable')
      istat = nf_put_att_text(ncid,idv_time,"units",len_trim(char80),
     |  trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of time dimension variable')
!
! Source date and time are attributes of time var 
! (they are also global file attributes):
      istat = nf_put_att_int(ncid,idv_time,"initial_year",NF_INT,1,
     |  h%initial_year)
      istat = nf_put_att_int(ncid,idv_time,"initial_day",NF_INT,1,
     |  h%initial_day)
      istat = nf_put_att_int(ncid,idv_time,"initial_mtime",NF_INT,3,
     |  h%initial_mtime)
!
! Define "primary" (non-coordinate) time dependent variables:
! (these variables get a value for each history) 
!
! Model time variable mtime(3,time), where 3 is the mtime dimension 
! (id_mtime) (day,hour,minute), and time is the unlimited dimension 
! (id_time):
!
      ids2(1) = id_mtime
      ids2(2) = id_time
      istat = nf_def_var(ncid,"mtime",NF_INT,2,ids2,idv_mtime)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining mtime variable')
      write(char80,"('model times (day, hour, minute)')")
      istat = nf_put_att_text(ncid,idv_mtime,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name attribute of mtime variable')
      istat = nf_put_att_text(ncid,idv_mtime,"units",17,
     |  'day, hour, minute')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of mtime variable')  
!
! Calendar year and day (this will be advanced if model is advancing
! in calendar time, otherwise will be constant from input)
      istat = nf_def_var(ncid,"year",NF_INT,1,ids1,idv_year)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining calendar year variable')
      istat = nf_put_att_text(ncid,idv_year,"long_name",
     |  13,'calendar year')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of calendar year variable')
!
      istat = nf_def_var(ncid,"day",NF_INT,1,ids1,idv_day)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining calendar day variable')
      istat = nf_put_att_text(ncid,idv_day,"long_name",
     |  12,'calendar day')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of calendar day variable')
!
! Number of time steps from model time 0,0,0 to current model time
! (iter(time)):
      istat = nf_def_var(ncid,"iter",NF_INT,1,ids1,idv_iter)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining iter variable')
      istat = nf_put_att_text(ncid,idv_iter,"long_name",
     |  42,'number of time steps from model time 0,0,0')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of iter variable')
!
! Decimal ut(time) (ids1==id_time is unlimited time dimension):
      call defvar_time_dbl(ncid,"ut",ids1,idv_ut,
     |  "universal time (from model time hour and minute)",
     |  "hours")
!
! Time step (seconds) (integer scalar) (time-independent):
      istat = nf_def_var(ncid,"timestep",NF_INT,0,idum,idv_step)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining timestep variable')
      istat = nf_put_att_text(ncid,idv_step,"long_name",8,'timestep')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of timestep variable')
      istat = nf_put_att_text(ncid,idv_step,"units",7,'seconds')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of timestep variable')
!
! f107 daily and f107 average solar fluxes:
      call defvar_time_dbl(ncid,"f107d",ids1,idv_f107d,
     |  "Daily 10.7 cm solar flux","W/m^2Hz*1.e-22")
      call defvar_time_dbl(ncid,"f107a",ids1,idv_f107a,
     |  "81-day average 10.7 cm solar flux","W/m^2Hz*1.e-22")
!
! Hemispheric power: 
      call defvar_time_dbl(ncid,"hpower",ids1,idv_hpower,
     |  "Hemispheric power","Gw")
!
! Cross-tail potential:
      call defvar_time_dbl(ncid,"ctpoten",ids1,idv_ctpoten,
     |  "Cross-tail potential","Volts")
!
! Byimf:
      call defvar_time_dbl(ncid,"byimf",ids1,idv_byimf,
     |  "BY component of IMF"," ")
!
! Collision factor:
      call defvar_time_dbl(ncid,"colfac",ids1,idv_colfac,
     |  "ion/neutral collision factor"," ")
!
! 1/31/01: Making tides time-dependent (dtide(time,2) and sdtide(time,10):
! Diurnal tide amp/phase:
      ids2(1) = id_dtide
      ids2(2) = id_time
      istat = nf_def_var(ncid,"dtide",NF_DOUBLE,2,ids2,idv_dtide)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining variable dtide')
      write(char80,"('amplitude and phase of diurnal tide')")
      istat = nf_put_att_text(ncid,idv_dtide,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name attribute of variable dtide')
!
! Semi-diurnal tide amp/phase:
      ids2(1) = id_sdtide
      ids2(2) = id_time
      istat = nf_def_var(ncid,"sdtide",NF_DOUBLE,2,ids2,idv_sdtide)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining variable sdtide')
      write(char80,"('amplitude and phase of semi-diurnal tide')")
      istat = nf_put_att_text(ncid,idv_sdtide,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name attribute of variable sdtide')
!
! 1/31/01: Adding ncep flag:
      istat = nf_def_var(ncid,"ncep",NF_INT,1,ids1,idv_ncep)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining ncep variable')
      write(char80,
     |  "('NCEP/NMC flag (1 if using ncep boundaries for t and z)')")
      istat = nf_put_att_text(ncid,idv_ncep,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of ncep variable')
!
! 2/08/01: Adding gpi flag:
      istat = nf_def_var(ncid,"gpi",NF_INT,1,ids1,idv_gpi)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining gpi variable')
      write(char120,
     |  "('GPI flag (1 if using geophysical indices database for ',
     |  'f107d, f107a, hpower, and ctpoten)')") 
      istat = nf_put_att_text(ncid,idv_gpi,"long_name",
     |  len_trim(char120),trim(char120))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of gpi variable')
!
! 7/04/02: Adding gswmdi flag:
      istat = nf_def_var(ncid,"gswmdi",NF_INT,1,ids1,idv_gswmdi)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining gswm variable')
      write(char120,
     |  "('GSWMDI flag (1 if using GSWM diurnal tides as boundary ',
     |  'for Z,TN,UN,VN)')") 
      istat = nf_put_att_text(ncid,idv_gswmdi,"long_name",
     |  len_trim(char120),trim(char120))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining long_name of gswmdi variable')
!
! 7/04/02: Adding gswmsdi flag:
      istat = nf_def_var(ncid,"gswmsdi",NF_INT,1,ids1,idv_gswmsdi)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining gswmsdi variable')
      write(char120,
     |  "('GSWMSDI flag (1 if using GSWM semidiurnal tides as ',
     |  'boundary for Z,TN,UN,VN)')") 
      istat = nf_put_att_text(ncid,idv_gswmsdi,"long_name",
     |  len_trim(char120),trim(char120))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining long_name of gswmsdi variable')
!
! 7/18/02: Adding gswmnmidi flag:
      istat = nf_def_var(ncid,"gswmnmidi",NF_INT,1,ids1,idv_gswmnmidi)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining gswmnmidi variable')
      write(char120,
     |  "('GSWMNMIDI flag (1 if using GSWM nonmigrating diurnal ',
     |  'tides as boundary for Z,TN,UN,VN)')") 
      istat = nf_put_att_text(ncid,idv_gswmnmidi,"long_name",
     |  len_trim(char120),trim(char120))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining long_name of gswmnmidi variable')
!
! 7/19/02: Adding gswmnmisdi flag:
      istat=nf_def_var(ncid,"gswmnmisdi",NF_INT,1,ids1,idv_gswmnmisdi)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining gswmnmisdi variable')
      write(char120,
     |  "('GSWMNMISDI flag (1 if using GSWM nonmigrating semidiurnal ',
     |  'tides as boundary for Z,TN,UN,VN)')") 
      istat = nf_put_att_text(ncid,idv_gswmnmisdi,"long_name",
     |  len_trim(char120),trim(char120))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining long_name of gswmnmisdi variable')
!
! 2/14/01: Add auroral parameters (from input):
!
! alfa30:
      call defvar_time_dbl(ncid,"alfa30",ids1,idv_alfa30,
     |  "Characteristic energy of high-energy auroral electrons",
     |  "KeV")
!
! e30:
      call defvar_time_dbl(ncid,"e30",ids1,idv_e30,
     |  "Energy flux of high-energy auroral electrons",
     |  "ergs/cm2/s")
!
! alfad2:
      call defvar_time_dbl(ncid,"alfad2",ids1,idv_alfad2,
     |  "Characteristic energy of solar protons in the polar cap",
     |  "KeV")
!
! ed2:
      call defvar_time_dbl(ncid,"ed2",ids1,idv_ed2,
     |  "Energy flux of solar protons in the polar cap",
     |  "ergs/cm2/s")
!
! Time-independent "primary" variables:
!
! Magnetic pole coordinates (real(2,2)):
      ids2(1) = id_latlon
      ids2(2) = id_latlon
      istat = nf_def_var(ncid,"mag",NF_DOUBLE,2,ids2,idv_mag)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining variable mag')
      write(char80,"('lat,lon coordinates of S,N magnetic poles')")
      istat = nf_put_att_text(ncid,idv_mag,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name attribute of variable mag')
      istat = nf_put_att_text(ncid,idv_mag,"units",7,
     |  'degrees')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of variable mag')  
!
! Standard pressure:
      istat = nf_def_var(ncid,"p0",NF_DOUBLE,0,idum,idv_p0)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining variable p0')
      istat = nf_put_att_text(ncid,idv_p0,"long_name",
     |  17,"Standard pressure")
!
! Number of fields (will be nf4d for primary histories):
      istat = nf_def_var(ncid,"nflds",NF_INT,0,idum,idv_nflds)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining nflds variable')
      write(char80,"('number of 3-d model fields')")
      istat = nf_put_att_text(ncid,idv_nflds,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name attribute of variable nflds')
!
! All prognostic fields are on primary history:
      if (h%type(1:3)=='pri') then
        do i=1,nf4d
          call defvar_f4d(ncid,f4d(i),idv_f4d(i))
        enddo
!
! Define secondary history fields (geographic, magnetic and magnetospheric):
      elseif (h%type(1:3)=='sec') then
        call defvar_sech(ncid,h%nfgeo,h%nfmag,h%nfgeo2d,h%nfmag2d,
     |    h%nfmagphr)
      else
        write(6,"('>>> nc_define: unknown h%type=',a)") h%type
      endif
!
! Define global file attributes:
!
! Run date and time:
      istat = nf_put_att_text(ncid,NF_GLOBAL,"rundate",
     +  len_trim(h%rundate),trim(h%rundate))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for rundate global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! User login name:
      istat = nf_put_att_text(ncid,NF_GLOBAL,"logname",
     +  len_trim(h%logname),trim(h%logname))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for logname global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Machine host name:
      istat = nf_put_att_text(ncid,NF_GLOBAL,"host",
     +  len_trim(h%host),trim(h%host)) 
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for host global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Machine operating system (from pre-proc macro)
      istat = nf_put_att_text(ncid,NF_GLOBAL,"system",
     +  len_trim(h%system),trim(h%system)) 
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for system global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Model name string:
      istat = nf_put_att_text(ncid,NF_GLOBAL,"model_name",
     +  len_trim(h%model_name),trim(h%model_name)) 
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for model name global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Model version string:
      istat = nf_put_att_text(ncid,NF_GLOBAL,"model_version",
     +  len_trim(h%model_version),trim(h%model_version)) 
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for model version global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Mss path to history file:
      istat = nf_put_att_text(ncid,NF_GLOBAL,"mss_path",
     |  len_trim(h%mss_path),trim(h%mss_path))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for mss_path global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Primary or Secondary histories:
      if (h%type(1:3)=='pri') then
        istat = nf_put_att_text(ncid,NF_GLOBAL,"history_type",
     |    7,'primary')
      elseif (h%type(1:3)=='sec') then
        istat = nf_put_att_text(ncid,NF_GLOBAL,"history_type",
     |    9,'secondary')
      else
        write(6,"('>>> nc_define: unknown h%type=',a)") h%type
        istat = nf_put_att_text(ncid,NF_GLOBAL,"history_type",
     |    7,'unknown')
      endif
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for history_type global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Mss path to source history file:
      if (nsource > 0) then ! initial source file
        write(6,"('h%mss_source=',a)") h%mss_source
        istat = nf_put_att_text(ncid,NF_GLOBAL,"mss_source",
     |    len_trim(h%mss_source)+10,trim(h%mss_source)//' (initial)')
      else ! continuation source file
        istat = nf_put_att_text(ncid,NF_GLOBAL,"mss_source",
     |    len_trim(h%mss_source)+15,
     |    trim(h%mss_source)//' (continuation)')
      endif
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for mss_source global attribute')")
        call handle_ncerr(istat,char120)
      endif
      istat = nf_put_att_int(ncid,idv_time,"source_mtime",NF_INT,3,
     |  h%source_mtime)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_int ',
     |    'for source_mtime=',3i4,' global attribute')") h%source_mtime
        call handle_ncerr(istat,char120)
      endif
!
! Mss path to secsource history file (magnetospheric model only):
      if (nsecsource > 0) then ! initial secsource file
        istat = nf_put_att_text(ncid,NF_GLOBAL,"mss_secsource",
     |    len_trim(h%mss_secsource)+10,trim(h%mss_secsource)//
     |    ' (initial)')
      endif
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for mss_secsource global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Start date and time of initial run:
      istat = nf_put_att_int(ncid,NF_GLOBAL,"initial_year",NF_INT,1,
     |  h%initial_year)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_att_int ',
     |    'for start_year global attribute: h%initial_year=',i4)")
     |    h%initial_year
        call handle_ncerr(istat,char120)
      endif
      istat = nf_put_att_int(ncid,NF_GLOBAL,"initial_day",NF_INT,1,
     |  h%initial_day)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_att_int ',
     |    'for start_day global attribute: h%initial_day=',i4)")
     |    h%initial_day
        call handle_ncerr(istat,char120)
      endif
      istat = nf_put_att_int(ncid,NF_GLOBAL,"initial_mtime",NF_INT,3,
     |  h%initial_mtime)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_att_int ',
     |    'for start_mtime global attribute: h%initial_mtime=',i4)")
     |    h%initial_mtime
        call handle_ncerr(istat,char120)
      endif
!
! Missing value (nc_put_att_real is in util.F):
      var1(1) = spval
      istat = nc_put_att_real(ncid,NF_GLOBAL,"missing_value",1,var1)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_att_real ',
     |    'for missing_value global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Conversion from lev (ln(p0/p)) to pressure mb:
      write(char80,"('p0*e(-lev(k))*1.e-3')")
      istat = nf_put_att_text(ncid,NF_GLOBAL,"lev_to_mb",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for lev_to_mb global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Number of histories on the file:
      ivar1(1) = h%ihist ! this will be updated in nc_wrhist
      istat = nf_put_att_int(ncid,NF_GLOBAL,"nhist",NF_INT,1,ivar1)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_att_int ',
     |    'for nhist global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Delta minutes between histories on the file:
      ivar1(1) = h%delhmins
      istat = nf_put_att_int(ncid,NF_GLOBAL,"delhist_mins",NF_INT,
     |  1,ivar1)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_att_int ',
     |    'for delhist_mins global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Exit define mode:
      istat = nf_enddef(ncid)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error return from nf_enddef')
!
! Mag grid coordinates (consdyn.h) were set in init_consdyn (cons_mod.F)
! (convert to degrees):
!     do i=1,nmlonp1
!       gmlon(i) = ylonm(i)*rtdeg
!     enddo
!     do i=1,nmlat
!       gmlat(i) = ylatm(i)*rtdeg
!     enddo
!     do i=1,nlev
!       pmlev(i+3) = plev(i)
!     enddo
!     do i=3,1,-1
!       pmlev(i) = pmlev(i+1)-dlev
!     enddo
!
! Give values to dimension variables:
!
! Geographic:
      istat = nc_put_var_real(ncid,idv_lon,glon)
      istat = nc_put_var_real(ncid,idv_lat,glat)
      istat = nc_put_var_real(ncid,idv_lev,plev)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_var_real ',
     +    ' to assign values to geographic dimension vars')")
        call handle_ncerr(istat,char120)
      endif
!
! Magnetic:
      istat = nc_put_var_real(ncid,idv_mlon,gmlon)
      istat = nc_put_var_real(ncid,idv_mlat,gmlat)
      istat = nc_put_var_real(ncid,idv_mlev,pmlev)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_var_real ',
     +    ' to assign values to magnetic dimension vars')")
        call handle_ncerr(istat,char120)
      endif
!
! Magnetospheric:
      istat = nc_put_var_real(ncid,idv_magphrlon,ylonmagphr)
      istat = nc_put_var_real(ncid,idv_magphrlat,ylatmagphr)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_var_real ',
     +    ' to assign values to magnetospheric dimension vars')")
        call handle_ncerr(istat,char120)
      endif
      end subroutine nc_define
!-------------------------------------------------------------------
      subroutine defvar_time_dbl(ncid,name,idtime,idvar,
     |  longname,units)

!     call defvar_time_dbl(ncid,"ut",ids1,idv_ut,
!    |  "universal time (from model time hour and minute)",
!    |  "hours")
!
! Define time-dependent variable "name" on open netcdf file ncid.
! Variable has a single dimension in time (e.g., ut(time)).
! Also define long_name and units attributes.
!
! Args:
      integer,intent(in)  :: ncid,idtime(1)
      integer,intent(out) :: idvar
      character(len=*),intent(in) :: name,longname,units
!
! Local:
      integer :: istat
      character(len=80) :: char80
!
! Define the variable:
      istat = nf_def_var(ncid,name,NF_DOUBLE,1,idtime,idvar)
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining double variable ',a)") name
        call handle_ncerr(istat,trim(char80))
      endif
!
! Give long name attribute:
      if (len_trim(longname) > 0) then
        istat = nf_put_att_text(ncid,idvar,"long_name",
     |    len_trim(longname),longname)
        if (istat /= NF_NOERR) then
          write(char80,"('Error defining long_name of double variable ',
     |      a)") name
          call handle_ncerr(istat,trim(char80))
        endif
      endif
!
! Give units attribute:
      if (len_trim(units) > 0) then
        istat = nf_put_att_text(ncid,idvar,"units",len_trim(units),
     |    units)
        if (istat /= NF_NOERR) then
          write(char80,"('Error defining units of variable ',a)") name
          call handle_ncerr(istat,trim(char80))
        endif
      endif
      end subroutine defvar_time_dbl
!-----------------------------------------------------------------------
      subroutine defvar_sech(ncid,nfgeo,nfmag,nfgeo2d,nfmag2d,nfmagphr)
!
! Define secondary history fields on currently open ncid
! (not data -- see wrf3d for data write)
!
! Args:
      integer,intent(in) :: 
     |  ncid,    ! netcdf file id
     |  nfgeo,   ! number of fields on geographic grid
     |  nfmag,   ! number of fields on magnetic grid
     |  nfgeo2d, ! number of fields on geographic grid 2d
     |  nfmag2d, ! number of fields on magnetic grid 2d
     |  nfmagphr ! number of fields on magnetospheric grid
!
! Local:
      integer :: i,ii,ier,istat,iprog
      character(len=80) :: char80
!
! External:
       integer,external :: strloc
!
      if (allocated(idv_sech)) deallocate(idv_sech)
      allocate(idv_sech(nfgeo+nfmag+nfgeo2d+nfmag2d+nfmagphr),stat=ier)
      if (ier /= 0) then
        write(6,"(/,'>>> WARNING defvar_sech: ',
     |  'error return from allocate for idv_sech: nfgeo=',i3,
     |  ' nfmag=',i3,' nfgeo2dr=',i3,
     |  ' nfmag2d=',i3,' nfmagphr=',i3,' ier=',i4)")
     |  nfgeo,nfmag,nfgeo2d,nfmag2d,nfmagphr,ier
      endif
      if (nfgeo > 0) then
        do i=1,nfgeo
          iprog = strloc(f4d%short_name,size(f4d),fsech(i)%short_name)
!
! Define prognostic secondary history field:
          if (iprog > 0) then
            call defvar_f4d(ncid,f4d(iprog),idv_sech(i))
!
! Define secondary history field:
          else
            call defvar_f3d(ncid,fsech(i),idv_sech(i))
          endif
        enddo
      endif
      if (nfmag > 0) then
        do i=1,nfmag
          call defvar_f3d(ncid,fsechmag(i),idv_sech(nfgeo+i))
        enddo
      endif
!     subroutine defvar_f2d(ncid,f,idvar)
      if (nfgeo2d > 0) then
        do i=1,nfgeo2d
          call defvar_f2d(ncid,fsech2d(i),idv_sech(nfgeo+nfmag+i))
        enddo
      endif
      if (nfmag2d > 0) then
        do i=1,nfmag2d
          call defvar_f2d(ncid,fsechmag2d(i),idv_sech(nfgeo+nfmag+
     |       nfgeo2d+i))
        enddo
      endif
      if (nfmagphr > 0) then
        do i=1,nfmagphr
          call defvar_f2d(ncid,fsechmagphr2d(i),idv_sech(nfgeo+nfmag+
     |       nfgeo2d+nfmag2d+i))
        enddo
      endif
      end subroutine defvar_sech
!-------------------------------------------------------------------
      subroutine defvar_f4d(ncid,f,idvar)
!
! Define a prognostic variable on current netcdf history file:
! (not data -- see wrf4d for data write)
!
! Args:
      integer,intent(in)  :: ncid
      integer,intent(out) :: idvar
      type(fields_4d),intent(in) :: f
!
! Local:
      integer :: istat,idimids(4)
      character(len=80) :: char80
      character(len=160) :: char160
!
! Prognostics are dimensioned (lon,lat,lev,time):
      if (.not.fakeflds) then
        idimids(1) = id_lon
        idimids(2) = id_lat
        idimids(3) = id_lev
        idimids(4) = id_time
      else
        idimids(1:3) = id_fake
        idimids(4) = id_time
      endif
! This call returns the variable id, idvar.
      istat = nf_def_var(ncid,f%short_name,NF_DOUBLE,4,idimids,idvar)
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining 4-d field ',a)") 
     |    f%short_name
        call handle_ncerr(istat,trim(char80))
      endif
!
! Add long_name attribute:
      istat = nf_put_att_text(ncid,idvar,"long_name",
     |  len_trim(f%long_name),trim(f%long_name))
!
! Cray apparently crashed here after 1380 steps 
! (char80 was used)
      if (istat /= NF_NOERR) then
        write(char160,"('Error defining long_name of 4-d',
     |    ' variable ',a,': long_name = ',/,a)") trim(f%short_name),
     |    trim(f%long_name)
        call handle_ncerr(istat,trim(char160))
      endif
!
! Add units attribute:
      istat = nf_put_att_text(ncid,idvar,"units",len_trim(f%units),
     |  f%units)
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining units of 4-d',
     |    ' variable ',a,': units = ',a)") trim(f%short_name),
     |    trim(f%units)
        call handle_ncerr(istat,trim(char80))
      endif
      end subroutine defvar_f4d
!-------------------------------------------------------------------
      subroutine defvar_f3d(ncid,f,idvar)
!
! Define a diagnostic variable (geographic and magnetic)
! on the current netcdf history file:
!
! Args:
      integer,intent(in)  :: ncid
      integer,intent(out) :: idvar
      type(fields_3d),intent(inout) :: f
!
! Local:
      integer :: istat,idimids(4),ndims
      character(len=80) :: char80
!
      if (.not.fakeflds) then
        if (.not.f%magnetos.and..not.f%magnetic) then ! field is on geographic grid
          ndims = 4
          idimids(1) = id_lon
          idimids(2) = id_lat
          idimids(3) = id_lev
          idimids(4) = id_time
        elseif (.not.f%magnetos) then    ! field is on magnetic grid
          ndims = 4
          idimids(1) = id_mlon
          idimids(2) = id_mlat
          idimids(3) = id_mlev
          idimids(4) = id_time
        endif
      else
        idimids(1:3) = id_fake
        idimids(4) = id_time
      endif
!
      istat = nf_def_var(ncid,f%short_name,NF_DOUBLE,ndims,
     |   idimids,idvar)
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining diagnostic field ',a)") 
     |    f%short_name
        call handle_ncerr(istat,trim(char80))
!     else
!       write(6,"('defvar_f3d: defined diagnostic field ',a)")
!    |    f%short_name
      endif
      end subroutine defvar_f3d
!-------------------------------------------------------------------
      subroutine defvar_f2d(ncid,f,idvar)
!
! Define a diagnostic variable (geographic, magnetic or magnetosphere)
! on the current netcdf history file:
!
! Args:  
      integer,intent(in)  :: ncid
      integer,intent(out) :: idvar
      type(fields_2d),intent(inout) :: f
!
! Local: 
      integer :: istat,idimids(3),ndims
      character(len=80) :: char80
!
      if (.not.fakeflds) then
        if (.not.f%magnetos.and..not.f%magnetic) then ! 2d geographic
          idimids(1) = id_lon
          idimids(2) = id_lat
        elseif (.not.f%magnetos) then    ! field is on magnetic 2d grid
          idimids(1) = id_mlon
          idimids(2) = id_mlat
        else                             ! field is on magnetospheric 2d grid
          idimids(1) = id_magphrlon
          idimids(2) = id_magphrlat
         endif
         ndims = 3
         idimids(3) = id_time
      else
        idimids(1:2) = id_fake
        idimids(3) = id_time
      endif
!
      istat = nf_def_var(ncid,f%short_name,NF_DOUBLE,ndims,
     |   idimids,idvar)
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining diagnostic field ',a)")
     |    f%short_name
        call handle_ncerr(istat,trim(char80))
!     else
!       write(6,"('defvar_f3d: defined diagnostic field ',a)")
!    |    f%short_name
      endif
      end subroutine defvar_f2d
!-------------------------------------------------------------------
      subroutine nc_wrhist(ncid,iprint)
!
! Write current history structure to open netcdf output file:
!
      use hist_module,only: h,hist_print
      use input_module,only: start,secstart
!
! Args:
      integer,intent(in) :: ncid,iprint
!
! Local:
      integer :: mins,istat,ivar1(1),i,j,iprog,imo,ida,mtimetmp(4)
      character(len=120) :: char120
      character(len=80) :: char80
      integer :: ndims,nvars,ngatts,unlimdimid
      real :: var1(1),var22(2,2),var2(2),var10(10),rmins,rmins1(1)
!
! External:
      integer,external :: strloc
      real,external :: mtime_delta
!
! Total model time in minutes:
! Normally time is current model time (mins) since start_day,mtime
!   (i.e., since start time of the initial run), however if model time 
!   is earlier than initial start time, use time since start time of 
!   current run.
! 
      mtimetmp(1) = h%initial_day
      mtimetmp(2:3) = h%initial_mtime(2:3) ; mtimetmp(4) = 0
      rmins = mtime_delta(mtimetmp,h%modeltime)
      if (rmins < 0.) then ! use time since start time of the run
        if (h%type(1:3)=='pri') then
          mtimetmp(1:3) = start(:,1) ; mtimetmp(4) = 0
          rmins = mtime_delta(mtimetmp,h%modeltime)
        elseif (h%type(1:3)=='sec') then
          mtimetmp(1:3) = secstart(:,1) ; mtimetmp(4) = 0
          rmins = mtime_delta(mtimetmp,h%modeltime)
        else
          write(6,"(/,'>>> WARNING nc_wrhist: unknown h%type=',a)") 
     |      h%type
          rmins = 0.
        endif
      endif
      if (idv_time <= 0) istat = nf_inq_varid(ncid,"time",idv_time)
      istat = nc_put_var1_real(ncid,idv_time,h%ihist,rmins,'time')
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nc_put_var_real defining ',
     |    'time at h%modeltime=',4i4,' rmins=',e12.4,
     |    ' h%initial_mtime=',3i4,' start(:,1)=',3i4,' secstart(:,1)=',
     |    3i4)") h%modeltime,rmins,h%initial_mtime,start(:,1),
     |    secstart(:,1)
        call handle_ncerr(istat,char120)
      endif
!
! Model time (day,hour,minute):
      start_2d(1) = 1
      start_2d(2) = h%ihist
      count_2d(1) = 3
      count_2d(2) = 1
      if (idv_mtime <= 0) istat = nf_inq_varid(ncid,"mtime",idv_mtime)
      istat = nf_put_vara_int(ncid,idv_mtime,start_2d,count_2d,
     |  h%modeltime(1:3))
      if (istat /= NF_NOERR) then
        write(char120,"('Error defining mtime at h%ihist=',i2,
     +    ' h%modeltime=',4i4,' istat=',i5)") h%ihist,h%modeltime,istat
        call handle_ncerr(istat,char120)
      endif
!
! Calendar year (initially from input, then updated if model is
! advancing calendar time):
      if (idv_year <= 0) istat = nf_inq_varid(ncid,"year",idv_year)
      istat = nf_put_var1_int(ncid,idv_year,h%ihist,h%year)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     |    'calendar year at h%ihist=',i3,': h%year=',i5)") 
     |    h%ihist,h%year
        call handle_ncerr(istat,char120)
      endif
!
! Calendar day (initially from input, then updated if model is
! advancing calendar time):
      if (idv_day <= 0) istat = nf_inq_varid(ncid,"day",idv_day)
      istat = nf_put_var1_int(ncid,idv_day,h%ihist,h%day)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     |    'calendar day at h%ihist=',i3,': h%day=',i5)") 
     |    h%ihist,h%day
        call handle_ncerr(istat,char120)
      endif
!
! Iteration number (number of steps from 0,0,0 to current mtime):
      if (idv_iter <= 0) istat = nf_inq_varid(ncid,"iter",idv_iter)
      istat = nf_put_var1_int(ncid,idv_iter,h%ihist,h%iter)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     +    'iter at h%ihist=',i3,': h%iter=',i8)") h%ihist,h%iter
        call handle_ncerr(istat,char120)
      endif
!
! Universal time (decimal hours):
      if (idv_ut <= 0) istat = nf_inq_varid(ncid,"ut",idv_ut)
      istat = nc_put_var1_real(ncid,idv_ut,h%ihist,h%ut,'ut')
!
! Magnetic pole coords:
      if (idv_mag <= 0) istat = nf_inq_varid(ncid,"mag",idv_mag)
      var22(:,:) = h%mag(:,:)
      istat = nc_put_var_real(ncid,idv_mag,var22)
!
! 1/31/01: tides are now time dependent 2d vars:
! Diurnal tide (dtide(time,2):
      start_2d(1) = 1
      start_2d(2) = h%ihist
      count_2d(1) = 2
      count_2d(2) = 1
      if (idv_dtide <= 0) istat = nf_inq_varid(ncid,"dtide",idv_dtide)
      var2(:) = h%dtide(:)
      istat = nf_put_vara_double(ncid,idv_dtide,start_2d,count_2d,var2)
!
! Semi-diurnal tide (sdtide(time,10):
      start_2d(1) = 1
      start_2d(2) = h%ihist
      count_2d(1) = 10
      count_2d(2) = 1
      if (idv_sdtide <= 0) istat = nf_inq_varid(ncid,"sdtide",
     |  idv_sdtide)
      var10(:) = h%sdtide(:)
      istat = nf_put_vara_double(ncid,idv_sdtide,start_2d,count_2d,
     |  var10)
!
! ncep/nmc flag:
      if (idv_ncep <= 0) istat = nf_inq_varid(ncid,"ncep",idv_ncep)
      istat = nf_put_var1_int(ncid,idv_ncep,h%ihist,h%ncep)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     +    'ncep at h%ihist=',i3,': h%ncep=',i3)") h%ihist,h%ncep
        call handle_ncerr(istat,char120)
      endif
!
! gpi flag:
      if (idv_gpi <= 0) istat = nf_inq_varid(ncid,"gpi",idv_gpi)
      istat = nf_put_var1_int(ncid,idv_gpi,h%ihist,h%gpi)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     +    'gpi at h%ihist=',i3,': h%gpi=',i3)") h%ihist,h%gpi
        call handle_ncerr(istat,char120)
      endif
!
! gswmdi flag for diurnal tides:
      if (idv_gswmdi <= 0) istat = nf_inq_varid(ncid,"gswmdi",
     |  idv_gswmdi)
      istat = nf_put_var1_int(ncid,idv_gswmdi,h%ihist,h%gswmdi)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     |   'gswmdi at h%ihist=',i3,': h%gswmdi=',i3)") h%ihist,h%gswmdi
        call handle_ncerr(istat,char120)
      endif
! gswmsdi flag for semidiurnal tides:
      if (idv_gswmsdi <= 0) istat = nf_inq_varid(ncid,"gswmsdi",
     |  idv_gswmsdi)
      istat = nf_put_var1_int(ncid,idv_gswmsdi,h%ihist,h%gswmsdi)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     |   'gswmsdi at h%ihist=',i3,': h%gswmsdi=',i3)") h%ihist,
     |   h%gswmsdi
        call handle_ncerr(istat,char120)
      endif
! gswmnmidi flag for nonmigrating diurnal tides:
      if (idv_gswmnmidi <= 0) istat = nf_inq_varid(ncid,"gswmnmidi",
     |  idv_gswmnmidi)
      istat = nf_put_var1_int(ncid,idv_gswmnmidi,h%ihist,h%gswmnmidi)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     |   'gswmnmidi at h%ihist=',i3,': h%gswmnmidi=',i3)") h%ihist,
     |   h%gswmnmidi
        call handle_ncerr(istat,char120)
      endif
! gswmnmisdi flag for nonmigrating semidiurnal tides:
      if (idv_gswmnmisdi <= 0) istat = nf_inq_varid(ncid,"gswmnmisdi",
     |  idv_gswmnmisdi)
      istat = nf_put_var1_int(ncid,idv_gswmnmisdi,h%ihist,h%gswmnmisdi)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     |   'gswmnmisdi at h%ihist=',i3,': h%gswmnmisdi=',i3)") h%ihist,
     |   h%gswmnmisdi
        call handle_ncerr(istat,char120)
      endif
!
! f107d and f107a (daily and average solar fluxes):
      if (idv_f107d <= 0) istat = nf_inq_varid(ncid,"f107d",
     |  idv_f107d)
      istat = nc_put_var1_real(ncid,idv_f107d,h%ihist,h%f107d,'f107d')
      if (idv_f107a <= 0) istat = nf_inq_varid(ncid,"f107a",
     |  idv_f107a)
      istat = nc_put_var1_real(ncid,idv_f107a,h%ihist,h%f107a,'f107a')
!
! Hemispheric power:
      if (idv_hpower <= 0) istat = nf_inq_varid(ncid,"hpower",
     |  idv_hpower)
      istat = nc_put_var1_real(ncid,idv_hpower,h%ihist,h%hpower,
     |  'hpower')
!
! Cross-tail potential:
      if (idv_ctpoten <= 0) istat = nf_inq_varid(ncid,"ctpoten",
     |  idv_ctpoten)
      istat = nc_put_var1_real(ncid,idv_ctpoten,h%ihist,h%ctpoten,
     |  'ctpoten')
!
! BY imf:
      if (idv_byimf <= 0) istat = nf_inq_varid(ncid,"byimf",idv_byimf)
      istat = nc_put_var1_real(ncid,idv_byimf,h%ihist,h%byimf,
     |  'byimf')
!
! Collision factor:
      if (idv_colfac <= 0) istat = nf_inq_varid(ncid,"colfac",
     |  idv_colfac)
      istat = nc_put_var1_real(ncid,idv_colfac,h%ihist,h%colfac,
     |  'colfac')
!
! alfa30:
      if (idv_alfa30 <= 0) istat = nf_inq_varid(ncid,"alfa30",
     |  idv_alfa30)
      istat = nc_put_var1_real(ncid,idv_alfa30,h%ihist,h%alfa30,
     |  'alfa30')
!
! e30:
      if (idv_e30 <= 0) istat = nf_inq_varid(ncid,"e30",idv_e30)
      istat = nc_put_var1_real(ncid,idv_e30,h%ihist,h%e30,'e30')
!
! alfad2:
      if (idv_alfad2 <= 0) istat = nf_inq_varid(ncid,"alfad2",
     |  idv_alfad2)
      istat = nc_put_var1_real(ncid,idv_alfad2,h%ihist,h%alfad2,
     |  'alfad2')
!
! ed2:
      if (idv_ed2 <= 0) istat = nf_inq_varid(ncid,"ed2",idv_ed2)
      istat = nc_put_var1_real(ncid,idv_ed2,h%ihist,h%ed2,'ed2')
!
! Standard pressure:
      if (idv_p0 <= 0) istat = nf_inq_varid(ncid,"p0",idv_p0)
      istat = nc_put_var1_real(ncid,idv_p0,h%ihist,h%p0,'p0')
!
! Timestep (seconds):
      if (idv_step <= 0) istat = nf_inq_varid(ncid,"timestep",idv_step)
      istat = nf_put_var1_int(ncid,idv_step,h%ihist,h%step)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     +    'ut at h%ihist=',i3,': h%step=',i4)") h%ihist,h%step
        call handle_ncerr(istat,char120)
      endif
!
! Update number of histories on the file (global attribute):
      ivar1(1) = h%ihist
      istat = nf_put_att_int(ncid,NF_GLOBAL,"nhist",NF_INT,1,ivar1)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nc_put_att_int ',
     +    'for nhist global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Number of fields (scalar):
      if (idv_nflds <= 0) istat = nf_inq_varid(ncid,"nflds",idv_nflds)
!     istat = nf_put_var1_int(ncid,idv_nflds,h%ihist,h%nflds)
      istat = nf_put_var_int(ncid,idv_nflds,h%nflds)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var_int defining ',
     +    'nflds: h%nflds=',i4)") h%nflds
        call handle_ncerr(istat,char120)
      endif
!
! Add fields data to primary history:
! Currently (2/00), only prognostic fields are written to the
!   primary histories:
!
      if (h%type(1:3)=='pri') then
        do i=1,nf4d_hist
          call wrf4d(ncid,f4d(i)%short_name,h%ihist,fakeflds,i)
        enddo
!
! Add fields data to secondary history:
! If a prognostic field, write field from fg to the history 
!   (sub wrf4d).
! If a diagnostic field, write field from fsech%data
!   (fsech%data was init to spval by sub set_fsech. If istep > 0,
!    then fsech%data should have been defined by addfsech.
!
      elseif (h%type(1:3)=='sec') then
        if (h%nfgeo > 0) then
          do i=1,h%nfgeo ! geographic fields
            iprog = strloc(f4d%short_name,nf4d,fsech(i)%short_name)
            if (iprog > 0) then
              call wrf4d(ncid,fsech(i)%short_name,h%ihist,fakeflds,
     |          iprog)
            else
              call wrf3d(ncid,fsech(i),h%ihist,fakeflds)
            endif
          enddo
        endif
        if (h%nfmag > 0) then
          do i=1,h%nfmag ! magnetic fields
            call wrf3d(ncid,fsechmag(i),h%ihist,fakeflds)
          enddo
        endif
        if (h%nfgeo2d > 0) then
          do i=1,h%nfgeo2d ! geographic 2d fields
            call wrf2d(ncid,fsech2d(i),h%ihist,fakeflds)
          enddo
        endif
        if (h%nfmag2d > 0) then
          do i=1,h%nfmag2d ! magnetic 2d fields
            call wrf2d(ncid,fsechmag2d(i),h%ihist,fakeflds)
          enddo
        endif
        if (h%nfmagphr > 0) then
          do i=1,h%nfmagphr ! magnetospheric 2d fields
            call wrf2d(ncid,fsechmagphr2d(i),h%ihist,fakeflds)
          enddo
        endif
      else
        write(6,"(/,'>>> WARNING nc_wrhist: unknown h%type=',a)") 
     |    h%type
      endif
!
! Report to stdout:
      if (iprint > 0) call hist_print(h,'WRITE')
      end subroutine nc_wrhist
!-------------------------------------------------------------------
      subroutine wrf2d(ncid,f,itime,fake)
!
! Write diagnostic field to current open history file:
!
! f%data should be allocated and defined.
! If the model has taken at least one timestep (istep > 0), then 
!   f%data should have been defined by user-called sub addfsech. 
! If istep==0, then f%data was init to spval by set_fsech at beginning 
!   of run.
!
      use init_module,only: istep
!
! Args:
      integer,intent(in) :: ncid,itime
      type(fields_2d),intent(in) :: f 
      logical,intent(in) :: fake
!
! Local:
      integer :: k,i,j,istat,idvar,itype,ndims,iddims(2),natts,
     |  idimsizes(2)
      character(len=120) :: char120
      character(len=80) :: char80
      character(len=16) :: rdname
      real :: f2dgeo(nlon,nlat),
     |   f2dmag(nmlonp1,nmlat),
     |   f2dmagphr(nmagphrlon,nmagphrlat) ! note i,j
      real :: fmin,fmax
      real :: fakevar(1,1)
!
! Fake means 2d fields are dimensioned (1,1) for testing:
! (see fakeflds in fields.F)
      if (fake) then
        istat = nf_inq_varid(ncid,f%short_name,idvar)
        if (istat /= NF_NOERR) then
          write(char120,"('wrf3d: Error getting id of field var ',a)")
     |      trim(f%short_name)
          call handle_ncerr(istat,char120)
        endif
        start_3d(1:2) = 1
        start_3d(3) = itime
        count_3d(1) = 1
        count_3d(2) = 1
        count_3d(3) = 1
        fakevar(1,1) = 0.
        istat = nf_put_vara_double(ncid,idvar,start_4d,count_4d,fakevar)
        if (istat /= NF_NOERR) then
          write(char120,"('Error return from nf_put_vara_double',
     +      ' for fake field var ',a,' itime=',i2)") f%short_name,itime
          call handle_ncerr(istat,char120)
        endif
        return
      endif
!
      if (istep > 0) then
        if (any(f%data(:,:) /= spval)) then ! was defined

        else
          write(6,"(/,'>>> WARNING wrf2d: field ',a,' apparently not ',
     |      'defined by addfsech.')") trim(f%short_name) 
        endif
      else  ! istep == 0

      endif
!
! Get field id:
      istat = nf_inq_varid(ncid,f%short_name,idvar)
      if (istat /= NF_NOERR) then
        write(char80,"('wrf2d: Error getting id of field var ',a)")
     |    trim(f%short_name)
        call handle_ncerr(istat,char80)
      endif
!
! Get info about the field: 
      istat = nf_inq_var(ncid,idvar,rdname,itype,ndims,iddims,natts)
      if (ndims /= 3) then
        write(6,"(/,'>>> WARNING wrf2d: bad ndims=',i3,
     |    ' (every diagnostic should have 3 dimensions)')") ndims
      endif
!
! Get info about dimensions:
      do i=1,ndims
        istat = nf_inq_dim(ncid,iddims(i),rdname,idimsizes(i))
      enddo
!
! Check dimension sizes:
      if (.not.fakeflds) then
        if (.not.f%magnetic.and..not.f%magnetos) then ! geographic grid
          if (idimsizes(1) /= nlon .or. idimsizes(2) /= nlat) then
            write(6,"(/,'>>> WARNING wrf2d: bad dimension sizes',
     |        ' for geographic 2d diagnostic field ',a)") trim(rdname)
            write(6,"('  dim sizes=',2i4,' but should be ',
     |        'nlon,nlat=',2i4)") idimsizes(1:2),nlon,nlat
          endif
        elseif (.not.f%magnetos) then ! magnetic grid
          if (idimsizes(1) /= nmlonp1 .or. idimsizes(2) /= nmlat) then
            write(6,"(/,'>>> WARNING wrf2d: bad dimension sizes',
     |        ' for magnetic 2d diagnostic field ',a)") trim(rdname)
            write(6,"('  dim sizes=',2i4,' but should be ',
     |        'nmlonp1,nmlat=',2i4)") idimsizes(1:2),nmlonp1,nmlat
          endif
        else                      ! magnetospheric grid
          if (idimsizes(1) /= nmagphrlon .or.
     |        idimsizes(2) /= nmagphrlat) then
            write(6,"(/,'>>> WARNING wrf2d: bad dimension sizes',
     |        ' for magnetospheric diagnostic field ',a)") trim(rdname)
            write(6,"('  dim sizes=',2i4,' but should be ',
     |        'nmagphrlon,nmagphrlat=',3i4)") idimsizes(1:2),
     |        nmagphrlon,nmagphrlat
          endif

        endif
      endif
!
! Define f2d from f%data (netcdf does not like the pointer f%data):
!     local: real :: f2dgeo(nlonp4,nlat)
!     local: real :: f2dmag(nmlonp1,nmlat)
!     local: real :: f2dmagphr(nmagphrlon,nmagphrlat)
!     field.F:  allocate(fsech2d(i)%data(nlonp4,nlat))
!               allocate(fsechmag2d(i)%data(nmlonp1,nmlat))
!               allocate(fsechmagphr2d(i)%data(nmagphrlon,nmagphrlat))
!
      if (.not.f%magnetic.and..not.f%magnetos) then ! geographic grid
          do i=1,nlon
            f2dgeo(i,:) = f%data(i+2,:)  ! 1->nlon <= 3->nlon+2
          enddo

      elseif (.not.f%magnetos) then ! magnetic grid
          do i=1,nmlonp1
            f2dmag(i,:) = f%data(i,:) 
          enddo
      else                      ! magnetospheric grid
          do i=1,nmagphrlon
            f2dmag(i,:) = f%data(i,:)
          enddo
      endif
!
! Write data to the file:
      start_3d(1:2) = 1
      start_3d(3) = itime
      count_3d(3) = 1
      if (.not.f%magnetic.and..not.f%magnetos) then ! geographic grid
        count_3d(1) = nlon
        count_3d(2) = nlat
        istat = nf_put_vara_double(ncid,idvar,start_3d,count_3d,f2dgeo)
        if (istat /= NF_NOERR) then
          write(char80,"(/,'>>> wrf2d: error return from ',
     |      'nf_put_vara_double for geo2d field var ',a,' itime=',i2)") 
     |      trim(f%short_name),itime
          call handle_ncerr(istat,char80)
        endif
      elseif (.not.f%magnetos) then ! magnetic grid
        count_3d(1) = nmlonp1
        count_3d(2) = nmlat
        istat = nf_put_vara_double(ncid,idvar,start_3d,count_3d,f2dmag)
        if (istat /= NF_NOERR) then
          write(char80,"(/,'>>> wrf2d: error return from ',
     |      'nf_put_vara_double for mag2d field var ',a,' itime=',i2)") 
     |      trim(f%short_name),itime
          call handle_ncerr(istat,char80)
        endif
      else                      ! magnetospheric grid
        count_3d(1) = nmagphrlon
        count_3d(2) = nmagphrlat
        istat = nf_put_vara_double(ncid,idvar,start_3d,count_3d,f2dmag)
        if (istat /= NF_NOERR) then
          write(char80,"(/,'>>> wrf2d: error return from ',
     |      'nf_put_vara_double for magphr2d field var ',a,
     |      ' itime=',i2)") 
     |      trim(f%short_name),itime
          call handle_ncerr(istat,char80)
        endif
      endif
!
! If model has executed one or more time steps, the long_name and
! units of the diagnostic sech field may have been defined by the
! user-called addfsech routine:
!
      if (istep > 0) then
        if (len_trim(f%long_name) > 0 .or. len_trim(f%units) > 0) then
          istat = nf_redef(ncid) ! put dataset in define mode
          if (istat /= NF_NOERR) call handle_ncerr(istat,
     +      'Error return from nf_redef')
          if (len_trim(f%long_name) > 0) then
            istat = nf_put_att_text(ncid,idvar,"long_name",
     |        len_trim(f%long_name),f%long_name)
            if (istat /= NF_NOERR) then
              write(char80,"('Error defining long_name of diagnostic',
     |          ' variable ',a,': long_name = ',a)") 
     |          trim(f%short_name),trim(f%long_name)
              call handle_ncerr(istat,trim(char80))
            endif
          endif
          if (len_trim(f%units) > 0) then
            istat = nf_put_att_text(ncid,idvar,"units",
     |        len_trim(f%units),f%units)
            if (istat /= NF_NOERR) then
              write(char80,"('Error defining units of diagnostic',
     |          ' variable ',a,' units = ',a)") 
     |          trim(f%short_name),trim(f%units)
              call handle_ncerr(istat,trim(char80))
            endif
          endif
          istat = nf_enddef(ncid) ! end define mode
          if (istat /= NF_NOERR) call handle_ncerr(istat,
     |      'Error return from nf_enddef')
        endif ! long_name or units defined
      endif   ! istep > 0
!
! Write min,max of each field:
!     if (istep /= 0) then
!       call fminmaxspv(f3diag,nlon*nlat*nlev,fmin,fmax,spval)
!       write(6,"('wrf3d: Wrote field ',a,' istep=',i3,
!    |    ' 3d min,max=',2e12.4)") f%short_name(1:8),istep,fmin,fmax
!     endif
      end subroutine wrf2d
!-------------------------------------------------------------------
      subroutine wrf3d(ncid,f,itime,fake)
!
! Write diagnostic field to current open history file:
!
! f%data should be allocated and defined.
! If the model has taken at least one timestep (istep > 0), then 
!   f%data should have been defined by user-called sub addfsech. 
! If istep==0, then f%data was init to spval by set_fsech at beginning 
!   of run.
!
      use init_module,only: istep
!
! Args:
      integer,intent(in) :: ncid,itime
      type(fields_3d),intent(in) :: f 
      logical,intent(in) :: fake
!
! Local:
      integer :: k,i,j,istat,idvar,itype,ndims,iddims(4),natts,
     |  idimsizes(4)
      character(len=120) :: char120
      character(len=80) :: char80
      character(len=16) :: rdname
      real :: f3diag(nlon,nlat,nlevp1),f3dmag(nmlonp1,nmlat,nmlev) ! note i,j,k
      real :: fmin,fmax
      real :: fakevar(1,1,1)
!
! Fake means 3d fields are dimensioned (1,1,1) for testing:
! (see fakeflds in flds_mod.f)
      if (fake) then
        istat = nf_inq_varid(ncid,f%short_name,idvar)
        if (istat /= NF_NOERR) then
          write(char120,"('wrf3d: Error getting id of field var ',a)")
     |      trim(f%short_name)
          call handle_ncerr(istat,char120)
        endif
        start_4d(1:3) = 1
        start_4d(4) = itime
        count_4d(1) = 1
        count_4d(2) = 1
        count_4d(3) = 1
        count_4d(4) = 1
        fakevar(1,1,1) = 0.
        istat = nf_put_vara_double(ncid,idvar,start_4d,count_4d,fakevar)
        if (istat /= NF_NOERR) then
          write(char120,"('Error return from nf_put_vara_double',
     +      ' for fake field var ',a,' itime=',i2)") f%short_name,itime
          call handle_ncerr(istat,char120)
        endif
        return
      endif
!
! If the following flags are set, the any() statement can fail with NaNS.
! FFLAGS= -qinitauto=7FF7FFFF -qfloat=nans -qflttrap=invalid:enable 
!   -qsigtrap=xl__trcedump
!
      if (istep > 0) then
        if (any(f%data(:,:,:) /= spval)) then ! was defined
 
        else
          write(6,"(/,'>>> WARNING wrf3d: field ',a,' apparently not ',
     |      'defined by addfsech.')") trim(f%short_name) 
        endif
      endif
!
! Get field id:
      istat = nf_inq_varid(ncid,f%short_name,idvar)
      if (istat /= NF_NOERR) then
        write(char80,"('wrf3d: Error getting id of field var ',a)")
     |    trim(f%short_name)
        call handle_ncerr(istat,char80)
      endif
!
! Get info about the field: 
      istat = nf_inq_var(ncid,idvar,rdname,itype,ndims,iddims,natts)
      if (ndims /= 4) then
        write(6,"(/,'>>> WARNING wrf3d: bad ndims=',i3,
     |    ' (every diagnostic should have 4 dimensions)')") ndims
      endif
!
! Get info about dimensions:
      do i=1,ndims
        istat = nf_inq_dim(ncid,iddims(i),rdname,idimsizes(i))
      enddo
!
! Check dimension sizes:
      if (.not.fakeflds) then
        if (.not.f%magnetic) then ! geographic grid
          if (idimsizes(1) /= nlon .or. idimsizes(2) /= nlat) then
            write(6,"(/,'>>> WARNING wrf3d: bad dimension sizes',
     |        ' for geographic diagnostic field ',a)") trim(rdname)
            write(6,"('  dim sizes=',3i4,' but should be ',
     |        'nlon,nlat,nlev=',3i4)") idimsizes(1:3),nlon,nlat,nlev
          endif
        else                      ! magnetic grid
          if (idimsizes(1) /= nmlonp1 .or. idimsizes(2) /= nmlat) then
            write(6,"(/,'>>> WARNING wrf3d: bad dimension sizes',
     |        ' for magnetic diagnostic field ',a)") trim(rdname)
            write(6,"('  dim sizes=',3i4,' but should be ',
     |        'nmlonp1,nmlat,nmlev=',3i4)") idimsizes(1:3),nmlonp1,
     |        nmlat,nmlev
          endif

        endif
        if (idimsizes(3) < nlev) then
          write(6,"(/,'>>> WARNING wrf3d: dimsize(3) < nlev:',
     |      ' idimsizes(3)=',i4,' nlev=',i4)") idimsizes(3),nlev
        endif
      endif
!
! Define f3diag from f%data (netcdf does not like the pointer f%data):
!
      if (.not.f%magnetic) then ! geographic grid
        do k=1,nlevp1
          do i=1,nlon
            f3diag(i,:,k) = f%data(k,i+2,:) ! 1->nlon <= 3->nlon+2
          enddo
          
!         if (itime==2) then
!         do j=1,nlat
!           write(6,"('wrf3d: j=',i3,' k=',i3,' itime=',i3,' f3diag=',
!    |        /,(6e12.4))") j,k,itime,f3diag(:,j,k)
!         enddo
!         endif

        enddo
 
!       call fminmax(f3diag(:,:,:),nlevp1*nlon*nlat,fmin,fmax)
!       write(6,"('wrf3d: mag field ',a,' fmin,max=',2e12.4)") 
!    |    f%short_name,fmin,fmax

      else                      ! magnetic grid
!
! Secondary history array: fsechmag(n)%data(nlevp1+3,nmlonp1,nmlat)
!   (see allocation in set_fsech in fields.F)
! Local f3dmag(nmlonp1,nmlat,nmlev) ! note i,j,k
! Note nmlev==nlevp1+3 for extra levels at bottom, see transf in dynamo.F.
!   (mag fields in dynamo are (nmlonp1,-2:nlev))
! Transform from (k,i,j) to (i,j,k) for secondary history.
!
        do k=1,nmlev  ! nmlev = nlevp1+3
          do i=1,nmlonp1
            f3dmag(i,:,k) = f%data(k,i,:)
          enddo
        enddo
      endif
!
! Write data to the file:
      start_4d(1:3) = 1
      start_4d(4) = itime
      count_4d(4) = 1
      if (.not.f%magnetic) then ! geographic grid
        count_4d(1) = nlon
        count_4d(2) = nlat
        count_4d(3) = nlevp1
        istat = nf_put_vara_double(ncid,idvar,start_4d,count_4d,f3diag)
        if (istat /= NF_NOERR) then
          write(char80,"(/,'>>> wrf3d: error return from ',
     |      'nf_put_vara_double for geo field var ',a,' itime=',i2)") 
     |      trim(f%short_name),itime
          call handle_ncerr(istat,char80)
        endif
      else                      ! magnetic grid
        count_4d(1) = nmlonp1
        count_4d(2) = nmlat
        count_4d(3) = nmlev
        istat = nf_put_vara_double(ncid,idvar,start_4d,count_4d,f3dmag)
        if (istat /= NF_NOERR) then
          write(char80,"(/,'>>> wrf3d: error return from ',
     |      'nf_put_vara_double for mag field var ',a,' itime=',i2)") 
     |      trim(f%short_name),itime
          call handle_ncerr(istat,char80)
        endif
      endif
!
! If model has executed one or more time steps, the long_name and
! units of the diagnostic sech field may have been defined by the
! user-called addfsech routine:
!
      if (istep > 0) then
        if (len_trim(f%long_name) > 0 .or. len_trim(f%units) > 0) then
          istat = nf_redef(ncid) ! put dataset in define mode
          if (istat /= NF_NOERR) call handle_ncerr(istat,
     |      'Error return from nf_redef')
          if (len_trim(f%long_name) > 0) then
            istat = nf_put_att_text(ncid,idvar,"long_name",
     |        len_trim(f%long_name),f%long_name)
            if (istat /= NF_NOERR) then
              write(char80,"('Error defining long_name of diagnostic',
     |          ' variable ',a,': long_name = ',a)") 
     |          trim(f%short_name),trim(f%long_name)
              call handle_ncerr(istat,trim(char80))
            endif
          endif
          if (len_trim(f%units) > 0) then
            istat = nf_put_att_text(ncid,idvar,"units",
     |        len_trim(f%units),f%units)
            if (istat /= NF_NOERR) then
              write(char80,"('Error defining units of diagnostic',
     |          ' variable ',a,' units = ',a)") 
     |          trim(f%short_name),trim(f%units)
              call handle_ncerr(istat,trim(char80))
            endif
          endif
          istat = nf_enddef(ncid) ! end define mode
          if (istat /= NF_NOERR) call handle_ncerr(istat,
     |      'Error return from nf_enddef')
        endif ! long_name or units defined
      endif   ! istep > 0
!
! Write min,max of each field:
!     if (istep /= 0) then
!       if (.not.f%magnetic) then
!         call fminmaxspv(f3diag,nlon*nlat*nlev,fmin,fmax,spval)
!         write(6,"('wrf3d: Wrote geo field ',a,' istep=',i3,
!    |      ' 3d min,max=',2e12.4)") 
!    |      f%short_name(1:8),istep,fmin,fmax
!       else
!         call fminmaxspv(f3diag,nlon*nlat*nlev,fmin,fmax,spval)
!         write(6,"('wrf3d: Wrote mag field ',a,' istep=',i3,
!    |      ' 3d min,max=',2e12.4)") 
!    |      f%short_name(1:8),istep,fmin,fmax
!       endif
!     endif
      end subroutine wrf3d
!-------------------------------------------------------------------
      subroutine wrf4d(ncid,name,itime,fake,ix)
!
! Write prognostic field type f%data to current netcdf output 
! history file attached to ncid:
!
      use input_module,only: dynamo
      use init_module,only: istep
      use mpi_module,only: lon0,lon1,lat0,lat1
      use fields_module,only: dynpot
!
! Args:
      integer,intent(in) :: ncid,itime,ix
      character(len=*),intent(in) :: name
      logical,intent(in) :: fake
!
! Local:
      integer :: nxk,k,i,j,idvar,istat,lonbeg,lonend
      real :: f3d(nlon,nlat,nlevp1) ! note i,j,k for history
      character(len=120) :: char120
      real :: fakevar(1,1,1),fmin,fmax
!
! Fake means 3d fields are dimensioned (1,1,1):
! (see fakeflds in flds_mod.f)
      if (fake) then
        istat = nf_inq_varid(ncid,name,idvar)
        if (istat /= NF_NOERR) then
          write(char120,"('wrf4d: Error getting id of field var ',a)")
     |      trim(name)
          call handle_ncerr(istat,char120)
        endif
        start_4d(1:3) = 1
        start_4d(4) = itime
        count_4d(1) = 1
        count_4d(2) = 1
        count_4d(3) = 1
        count_4d(4) = 1
        fakevar(1,1,1) = 0.
        istat = nf_put_vara_double(ncid,idvar,start_4d,count_4d,fakevar)
        if (istat /= NF_NOERR) then
          write(char120,"('Error return from nf_put_vara_double',
     +      ' for fake field var ',a,' itime=',i2)") name,itime
          call handle_ncerr(istat,char120)
        endif
        return
      endif
!
! Transform from (k,i,j) to (i,j,k) for the history.
! This routine is executed by root task only. The root task has
! gathered all subdomain data from other tasks into foutput
! (see call mp_gather2root in advance). Foutput is a pointer to 
! all fields at 3d grid, declared in fields.F and allocated in 
! allocdata.F). 
!
! f3d(nlon,nlat,nlevp1)            ! note i,j,k for history
! dynpot(nlonp1,0:nlatp1,nlevp1)   ! 3d electric potential geographic
!
      if (trim(name) /= "POTEN") then
        do j=1,nlat
          do k=1,nlevp1
            f3d(:,j,k) = foutput(k,3:nlon+2,j,ix)
          enddo
        enddo
      else               
!
! Save geographic potential dynpot from dynamo.
!
        do j=1,nlat
          do k=1,nlevp1
            f3d(:,j,k) = dynpot(1:nlon,j,k)
          enddo
        enddo
      endif
!
! Get field id:
      istat = nf_inq_varid(ncid,name,idvar)
      if (istat /= NF_NOERR) then
        write(char120,"('wrf4d: Error getting id of field var ',a)")
     |    trim(name)
        call handle_ncerr(istat,char120)
      endif
!
! Put data onto netcdf file (this is where the majority of output
!   i/o happens). 
! (Should probably make a nc_put_vara_real for SUN vs non-SUN)
!
      start_4d(1:3) = 1
      start_4d(4) = itime
      count_4d(1) = nlon
      count_4d(2) = nlat
      count_4d(3) = nlevp1
      count_4d(4) = 1
      istat = nf_put_vara_double(ncid,idvar,start_4d,count_4d,f3d)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_vara_double',
     +    ' for field var ',a,' itime=',i2)") name,itime
        call handle_ncerr(istat,char120)
      endif
!
! Write min,max of each field if not echoing source history:
!     if (istep /= 0) then
!       call fminmaxspv(f3d,nlon*nlev*nlat,fmin,fmax,spval)
!       write(6,"('wrf4d: Wrote field ',a,' istep=',i3,
!    |    ' 3d min,max=',2e12.4)") name(1:8),istep,fmin,fmax
!     endif
      end subroutine wrf4d
!-------------------------------------------------------------------
      subroutine nc_fileinfo(ncid,fileinfo)
!
! A netcdf file on ncid has been written and is ready to be
! closed and disposed. This routine builds a fileinfo string
! giving beginning and ending history dates, model version string,
! and "primary" or "secondary". This is added as a global
! file attribute, and is included in the mss comment field
! if the file is disposed to the mss, and is therefore accessible 
! via msls -x. Example fileinfo string:
!
! 96080  80, 0, 0 to 96084  84, 0, 0 by  720 tgcm22 primary
!
! Args:
      integer,intent(in) :: ncid
      character(len=80),intent(out) :: fileinfo
!
! Local:
      integer :: istat,idvmtime,idtime,ntimes,idvyear,idvday,
     |  mtime_beg(3),mtime_end(3),istart2(2),icount2(2),
     |  iyear_beg,iyear_end,iday_beg,iday_end,idelhist,
     |  ibegiyd,iendiyd
      character(len=16) :: mversion,htype
      character(len=80) :: contents_desc
!
! Init:
      fileinfo = ' '
!
! Get number of histories (times) on the file:
      istat = nf_inq_unlimdim(ncid,idtime) ! id of unlimited record var
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'nc_fileinfo: Error getting unlimited dimension id')
      istat = nf_inq_dimlen(ncid,idtime,ntimes) ! length of time var
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'nc_fileinfo: Error getting length of time record variable')
!
! Get id for mtime variable:
      istat = nf_inq_varid(ncid,"mtime",idvmtime) 
      if (istat /= NF_NOERR)
     |  call handle_ncerr(istat,
     |    'nc_fileinfo: Error getting id of variable mtime')
!
! Get beginning and ending mtimes:
      istart2(1) = 1
      istart2(2) = 1
      icount2(1) = 3
      icount2(2) = 1
      istat = nf_get_vara_int(ncid,idvmtime,istart2,icount2,
     |  mtime_beg)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'nc_fileinfo: Error getting beginning mtime')
      istart2(2) = ntimes
      istat = nf_get_vara_int(ncid,idvmtime,istart2,icount2,
     |  mtime_end)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'nc_fileinfo: Error getting ending mtime')
!
! Get beginning and ending year:
      istat = nf_inq_varid(ncid,"year",idvyear) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'nc_fileinfo: Error getting id of variable year')
      istat = nf_get_var1_int(ncid,idvyear,1,iyear_beg)
      istat = nf_get_var1_int(ncid,idvyear,ntimes,iyear_end)
!
! Get beginning and ending day:
      istat = nf_inq_varid(ncid,"day",idvday) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'nc_fileinfo: Error getting id of variable day')
      istat = nf_get_var1_int(ncid,idvday,1,iday_beg)
      istat = nf_get_var1_int(ncid,idvday,ntimes,iday_end)
!
! Get delta time between histories (minutes), model version,
! and history type:
      istat = nf_get_att_int (ncid,NF_GLOBAL,"delhist_mins",idelhist)
      mversion = ' '
      istat = nf_get_att_text(ncid,NF_GLOBAL,"model_version",mversion)
      htype = ' '
      istat = nf_get_att_text(ncid,NF_GLOBAL,"history_type",htype)
!
! Construct file info string:
      ibegiyd = (iyear_beg-iyear_beg/100*100)*1000+iday_beg
      iendiyd = (iyear_end-iyear_end/100*100)*1000+iday_end
      write(fileinfo,"(i5.5,' ',i3,' ',i2,' ',i2,' to ',
     |                 i5.5,' ',i3,' ',i2,' ',i2,' by ',i4,
     |  ' ',a,' ',a)") ibegiyd, mtime_beg,iendiyd, mtime_end, 
     |  idelhist,trim(mversion),trim(htype)
!
! Add file info string and description of file info global attributes:
      istat = nf_redef(ncid)
      istat = nf_put_att_text(ncid,NF_GLOBAL,"contents",
     |  len_trim(fileinfo),trim(fileinfo))
      write(contents_desc,"('yyddd day hour min',
     |  ' to yyddd day hour min by delta_mins')")
      istat = nf_put_att_text(ncid,NF_GLOBAL,"contents_desc",
     |  len_trim(contents_desc),trim(contents_desc))
      istat = nf_enddef(ncid)
      end subroutine nc_fileinfo
!-------------------------------------------------------------------
      subroutine handle_ncerr(istat,msg)
!
! Handle a netcdf lib error:
!
      integer,intent(in) :: istat
      character(len=*),intent(in) :: msg
!
      write(6,"(/72('-'))")
      write(6,"('>>> Error from netcdf library:')")
      write(6,"(a)") trim(msg)
      write(6,"('istat=',i5)") istat
      write(6,"(a)") nf_strerror(istat)
      write(6,"(72('-')/)")
      return
      end subroutine handle_ncerr
!-------------------------------------------------------------------
!
! Remaining functions are netcdf wrappers for 32-bit vs 64-bit platforms.
!
!-------------------------------------------------------------------
      integer function nc_put_att_real(ncid,varid,name,len,vals)
!
! Args:
      integer,intent(in) :: ncid,varid,len
      character(len=*),intent(in) :: name
      real,intent(in) :: vals(len)
!
! Local:
      integer :: istat
!
#if defined(SUN)
!     istat = nf_put_att_real(ncid,varid,name,NF_FLOAT,
!    |  len,vals)
      istat = nf_put_att_double(ncid,varid,name,NF_FLOAT,
     |  len,vals)
#else
      istat = nf_put_att_double(ncid,varid,name,NF_DOUBLE,
     |  len,vals)
#endif
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'nc_put_att_real: error from nf_put_att_real/double')
      nc_put_att_real = istat
      end function nc_put_att_real
!-------------------------------------------------------------------
      integer function nc_put_var_real(ncid,idvar,vals)
!
! Args:
      integer,intent(in) :: ncid,idvar
      real,intent(in) :: vals(*)
!
! Local:
      integer :: istat
!
#if defined(SUN)
!     istat = nf_put_var_real(ncid,idvar,vals)
      istat = nf_put_var_double(ncid,idvar,vals)
#else
      istat = nf_put_var_double(ncid,idvar,vals)
#endif
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'nc_put_var_real: error from nf_put_var_real/double')
      nc_put_var_real = istat
      end function nc_put_var_real
!-------------------------------------------------------------------
      integer function nc_get_var_real(ncid,idvar,vals)
!
! Args:
      integer,intent(in) :: ncid,idvar
      real,intent(out) :: vals(*)
!
! Local:
      integer :: istat
!
#if defined(SUN)
!     istat = nf_get_var_real(ncid,idvar,vals)
      istat = nf_get_var_double(ncid,idvar,vals)
#else
      istat = nf_get_var_double(ncid,idvar,vals)
#endif
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'nc_get_var_real: error from nf_get_var_real/double')
      nc_get_var_real = istat
      end function nc_get_var_real
!-------------------------------------------------------------------
      integer function nc_get_var1_real(ncid,idvar,index,val)
!
! Args:
      integer,intent(in) :: ncid,idvar,index
      real,intent(out) :: val
!
! Local:
      integer :: istat
!
#if defined(SUN)
!     istat = nf_get_var1_real(ncid,idvar,index,val)
      istat = nf_get_var1_double(ncid,idvar,index,val)
#else
      istat = nf_get_var1_double(ncid,idvar,index,val)
#endif
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'nc_get_var1_real: error from nf_get_var1_real/double')
      nc_get_var1_real = istat
      end function nc_get_var1_real
!-------------------------------------------------------------------
      integer function nc_put_var1_real(ncid,idvar,index,val,name)
!
! Args:
      integer,intent(in) :: ncid,idvar,index
      real,intent(in) :: val
      character(len=*),intent(in) :: name
!
! Local:
      integer :: istat
      character(len=80) :: msg
!
#if defined(SUN)
!     istat = nf_put_var1_real(ncid,idvar,index,val)
      istat = nf_put_var1_double(ncid,idvar,index,val)
#else
      istat = nf_put_var1_double(ncid,idvar,index,val)
#endif
      if (istat /= NF_NOERR) then
        write(msg,"('nc_put_var1_real: error from nf_put_var1_real',
     |    '/double for variable ',a)") trim(name)
        call handle_ncerr(istat,msg)
      endif
      nc_put_var1_real = istat
      end function nc_put_var1_real
      end module nchist_module
