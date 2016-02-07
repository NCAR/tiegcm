!
      module nchist_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Module to create, read and write netcdf history files:
!
      use params_module,only: nlon,nlev,nmlevp1,nimlevp1,
     |  nlat,nlevp1,spval,nlonp4,nlonp2,nmlonp1,nmlat,nmlon,zmbot,zmtop,
     |  zibot,zitop,mxfsech,zpmag,zpimag,gmlon,gmlat,glon,glat,zpmid,
     |  zpint
      use fields_module,only: nf4d,nf4d_hist,f4d,fakeflds,shortname_len,
     |  fsechist,itc,itp,fields_4d,fields_3d,fields_2d,poten,foutput,
     |  tlbc,ulbc,vlbc,tlbc_glb,ulbc_glb,vlbc_glb,
     |  tlbc_nm,ulbc_nm,vlbc_nm,tlbc_nm_glb,ulbc_nm_glb,vlbc_nm_glb,
     |  nf3d,fzg
      use input_module,only: mxlen_filename,
     |  gswm_mi_di_ncfile, gswm_mi_sdi_ncfile, 
     |  gswm_nm_di_ncfile, gswm_nm_sdi_ncfile,ctmt_ncfile
      use init_module,only: istep,start_mtime
      implicit none
!
! Included file svn_revision.inc contains the svn revision number used in the build.
! (character(len=8) svn_revision)
!
#include "svn_revision.inc"
!
#ifdef SUN
#include <netcdf3.inc>
#else
#include <netcdf.inc>
#endif
!
      integer,private ::
     |  id_time,                   ! dimension id for time (unlimited)
     |  id_mtimedim,               ! dimension id for model time (3)
     |  id_lon,id_lat,             ! horizontal dimension ids
     |  id_lev,id_ilev,            ! vertical dimension ids (midpoints,interfaces)
     |  id_mlon,id_mlat,           ! horizontal dimension ids for magnetic grid
     |  id_mlev,id_imlev,          ! vertical dimension ids for magnetic grid
     |  id_latlon,                 ! dimension id for coords (2)
     |  id_dtide,                  ! dimension id for dtide (2)
     |  id_sdtide,                 ! dimension id for sdtide (10)
     |  id_fake,                   ! fake dimension id (1)
     |  id_datelen,                ! length of date and time
     |  id_filelen,                ! max length of a file name
     |  ids1(1),ids2(2),ids3(3),ids4(4),! vectors of dim id's
     |  start_2d(2),count_2d(2),	! start loc and count for 2d vars
     |  start_3d(3),count_3d(3),        ! start loc and count for 3d vars
     |  start_4d(4),count_4d(4),	! start loc and count for 4d vars
!
! Variable id's (idv_xxxx). Note lon, lat, lev, and time are
!   coordinate variables (lon(lon), lat(lat), lev(lev), time(time)).
! Variables starting at mtime are "primary" variables.
! idv_f4d(nf4d) are id's for prognostic 4d fields.
!
     |  idv_time,idv_mtime,idv_ut,idv_calendar_adv,
     |  idv_lon,idv_lat,idv_lev,idv_ilev,
     |  idv_lbc,idv_tlbc,idv_ulbc,idv_vlbc,
     |  idv_tlbc_nm,idv_ulbc_nm,idv_vlbc_nm,
     |  idv_mlon,idv_mlat,idv_mlev,idv_imlev,idv_p0_model,
     |  idv_step,idv_iter,idv_ntask,idv_year,idv_day,idv_p0,idv_grav,
     |  idv_zg,idv_hpower,idv_ctpoten,idv_bximf,idv_byimf,idv_bzimf,
     |  idv_colfac,idv_swden,idv_swvel,idv_al,idv_kp,
     |  idv_gswm_mi_di_ncfile,idv_gswm_mi_sdi_ncfile,
     |  idv_gswm_nm_di_ncfile,idv_gswm_nm_sdi_ncfile,
     |  idv_see_ncfile,idv_gpi_ncfile,idv_ncep_ncfile,
     |  idv_f107d,idv_f107a,idv_imf_ncfile,idv_mag,idv_dtide,
     |  idv_sdtide,idv_f4d(nf4d),idv_e1,idv_e2,idv_h1,idv_h2,idv_alfac,
     |  idv_ec,idv_alfad,idv_ed,idv_writedate,idv_saber_ncfile,
     |  idv_tidi_ncfile,idv_crit1,idv_crit2,idv_coupled_cmit,
     |  idv_joulefac,idv_bgrddata_ncfile,idv_ctmt_ncfile
      integer,allocatable :: idv_sech(:)
      integer :: idv_fsech(mxfsech)
!
! If check_nan is set, check for presence of NaN's in data when writing 
!   history. 
!
#ifdef AIX
      logical :: check_nan = .true.
#else
      logical :: check_nan = .false.
#endif
      character(len=16) :: datasource_kp=' '
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
!-----------------------------------------------------------------------
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
      if (istat /= NF_NOERR) then
        call handle_ncerr(istat,
     |    'Error return from nf_close')
      else
!       write(6,"('nc_close: closed ncid ',i3)") ncid
      endif
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
      use hist_module,only: sh,hist_initype,hist_print,nsource,
     |  nsecsource
      use input_module,only: start_day
!
! Args:
      integer,intent(in) :: ncid,mtime(3)
      character(len=*),intent(in) :: diskfile
      integer,intent(out) :: itime,ier
!
! Local:
      integer,parameter :: mxdims=20 ! assume <= mxdims dimensions
      integer,parameter :: mxvars=100 ! max # 4d vars
      integer :: istat,id_time,id_mtimedim,mtime_read(3),ntimes,
     |  i,ix,istart2(2),icount2(2),j,i4dims(4)
      integer :: ndims,nvars,ngatts,idunlim,itype,natts,len,
     |  iddims(mxdims),iscalar,id
      character(len=80) :: varname
      character(len=80) :: varnames(mxvars)=' '
      character(len=mxlen_filename) :: filename
      real :: scalar,var1(1),var22(2,2),var2(2),var10(10),varlbc
      real,dimension(nlevp1) :: zpmid_rd,zpint_rd
!
      integer,parameter :: mxf4d=100
      character(len=16) :: f4dnames(mxf4d) ! list of 4d fields read
      logical :: found 
      integer :: filelen ! file name length dimension on the file
!
! External:
      logical,external :: arrayeq ! util.F
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
      call hist_initype(sh,istep)  ! initialize history structure
!
! Define vars already found above:
      sh%ihist = itime
      sh%modeltime(1:3) = mtime
      sh%modeltime(4) = 0
!
! Read global attributes:
      istat = nf_get_att_text(ncid,NF_GLOBAL,"label",sh%label)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"create_date",
     |  sh%create_date)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"logname",sh%logname)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"host",sh%host)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"system",sh%system)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"model_name",sh%model_name)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"model_version",
     |  sh%model_version)
      istat = nf_get_att_int(ncid,NF_GLOBAL,"source_mtime",
     |  sh%source_mtime)
      istat = nf_get_att_text(ncid,NF_GLOBAL,"initial_file",
     |  sh%initial_file)
      if (istat /= NF_NOERR) ! was not on "old" histories
     |  sh%initial_file = ' '
!
! History_type is either "primary" or "secondary":
      istat=nf_get_att_text(ncid,NF_GLOBAL,"history_type",sh%hist_type)
!
! Check for tuv_lbc_intop (lbc of t,u,v stored in top slot)
!   tuv_lbc_intop == 0 --> new history (tuv lbc *not* stored in top slot)
!   tuv_lbc_intop == 1 --> old history (tuv lbc *is*  stored in top slot)
! If tuv_lbc_intop global attribute is not on the file, assume
! "old" history, i.e., tuv_lbc_intop==1:
!
      istat=nf_get_att_int(ncid,NF_GLOBAL,"tuv_lbc_intop",
     |  sh%tuv_lbc_intop)
      if (istat /= NF_NOERR) then ! old hist
        write(6,"('nc_rdhist: did not find tuv_lbc_intop in source',
     |    ' history ',a,' --> assuming tuv_lbc_intop==1')") 
     |    trim(diskfile)
        sh%tuv_lbc_intop = 1
      endif
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
            call shutdown('nlon')
          endif
          sh%nlon = len
        endif
        if (trim(varname)=='lat') then
          if (len /= nlat) then
            write(6,"(/,'>>> nc_rdhist: length of lat read from',
     |        ' history does not match parameter nlat.')")
            write(6,"('  len=',i4,' nlat=',i4)") len,nlat
            call shutdown('nlat')
          endif
          sh%nlat = len
        endif
        if (trim(varname)=='lev') then
!
! If number of levels read from source history (len) is less than the
! number of levels allocated in the model (nlev), this is fatal:
!
          if (len < nlevp1) then
            write(6,"(/,'>>> FATAL: source history vertical dimension',
     |        ' does not match model levels dimension.')")
            write(6,"(/,'nc_rdhist: length of lev dimension',
     |        ' on source history is less than model parameter ',
     |        'nlevp1.')")
            write(6,"('  len=',i4,' nlevp1=',i4)") len,nlevp1
            call shutdown('nlevp1')
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
          sh%nlev = len
        endif
      enddo
!
! Get max length of file name for calls to rdfilename:
      istat=nf_inq_dimid(ncid,'filelen',id) 
      if (istat /= NF_NOERR) 
     |  call handle_ncerr(istat,'get id of filelen dimension')
      istat = nf_inq_dimlen(ncid,id,filelen)
      if (istat /= NF_NOERR) 
     |  call handle_ncerr(istat,'get filelen dimension')
!
! Get variables (should not need to get coord vars lat,lon,lev)
      do i=1,nvars
        istat = nf_inq_var(ncid,i,varname,itype,ndims,iddims,natts)
        if (istat /= NF_NOERR) then
          write(6,"('>>> nc_rdhist: error ',
     |      'inquiring about var: i=',i4,' nvars=',i4)") i,nvars
          call handle_ncerr(istat,'inquiring about a var')
        endif
        select case(trim(varname))
          case('lon')
            istat = nf_get_var_double(ncid,i,glon)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading glon')
          case('lat') 
            istat = nf_get_var_double(ncid,i,glat)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading glat')
!
! New hist: lev(lev) is midpoints
          case('lev')
            if (sh%tuv_lbc_intop <= 0) then ! new hist (lev is midpoints)
              istat = nf_get_var_double(ncid,i,zpmid_rd)
              if (istat /= NF_NOERR) 
     |          call handle_ncerr(istat,'reading zpmid (new history)')

!             write(6,"('nc_rdhist sh%nlev=',i3,': read zpmid from new',
!    |          ' hist: zpmid_rd=',/,(10f7.3))") sh%nlev,zpmid_rd

              if (.not.arrayeq(zpmid,zpmid_rd,nlevp1)) then
                write(6,"(/,'>>> nc_rdhist: zpmid_rd (as read from new',
     |            ' history) is not equal to model zpmid')")
                write(6,"('nlevp1=',i3,' zpmid_rd=',/,(10f7.3))") 
     |            nlevp1,zpmid_rd
                write(6,"('nlevp1=',i3,' zpmid   =',/,(10f7.3))")
     |            nlevp1,zpmid
                call shutdown('zpmid')
              endif
!
! Old hist lev(lev) is interfaces:
            else                         ! old hist (lev is interfaces)
              istat = nf_get_var_double(ncid,i,zpint_rd)
              if (istat /= NF_NOERR) 
     |          call handle_ncerr(istat,'reading zpint (old history)')

!             write(6,"('nc_rdhist sh%nlev=',i3,': read zpint from old',
!    |          ' hist: zpint_rd=',/,(10f7.3))") sh%nlev,zpint_rd

              if (.not.arrayeq(zpint,zpint_rd,nlevp1)) then
                write(6,"(/,'>>> nc_rdhist: zpint_rd (as read from old',
     |            ' history) is not equal to model zpint')")
                write(6,"('nlevp1=',i3,' zpint_rd=',/,(10f7.3))") 
     |            nlevp1,zpint_rd
                write(6,"('nlevp1=',i3,' zpint   =',/,(10f7.3))")
     |            nlevp1,zpint
                call shutdown('zpint')
              endif
            endif
!
! ilev(ilev) interface levels are on new histories, or histories written
! by mksrc.ncl (i.e., new lbc from old histories):
          case('ilev')
            if (sh%tuv_lbc_intop <= 0) then ! new hist (lev is midpoints)
              istat = nf_get_var_double(ncid,i,zpint_rd)
              if (istat /= NF_NOERR) 
     |          call handle_ncerr(istat,'reading zpint (new history)')
              if (.not.arrayeq(zpint,zpint_rd,nlevp1)) then
                write(6,"(/,'>>> nc_rdhist: zpint_rd (as read from new',
     |            ' history) is not equal to model zpint')")
                write(6,"('nlevp1=',i3,' zpint_rd=',/,(10f7.3))") 
     |            nlevp1,zpint_rd
                write(6,"('nlevp1=',i3,' zpint   =',/,(10f7.3))")
     |            nlevp1,zpint
                call shutdown('zpint')
              endif
            else  ! old hist ilev interfaces (probably written by mksrc.ncl)
              istat = nf_get_var_double(ncid,i,zpint_rd)
              if (istat /= NF_NOERR) 
     |          call handle_ncerr(istat,'reading zpint (old history)')
              if (.not.arrayeq(zpint,zpint_rd,nlevp1)) then
                write(6,"(/,'>>> nc_rdhist: zpint_rd (as read from new',
     |            ' history) is not equal to model zpint')")
                write(6,"('nlevp1=',i3,' zpint_rd=',/,(10f7.3))") 
     |            nlevp1,zpint_rd
                write(6,"('nlevp1=',i3,' zpint   =',/,(10f7.3))")
     |            nlevp1,zpint
                call shutdown('zpint')
              endif
            endif
!
! Magnetic coordinates:
          case('mlon')
            istat = nf_get_var_double(ncid,i,gmlon)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading gmlon')
          case('mlat') 
            istat = nf_get_var_double(ncid,i,gmlat)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading gmlat')
          case('mlev')  ! mag midpoints

          case('imlev') ! mag interfaces
!
! Time:
          case('time')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading time')
            sh%time = scalar
!
! If this is a continuation run, read initial times
! (If an initial run, initial times were set by sub hist_initype (hist.F))
!
            if (nsource <= 0) then
              istat=nf_get_att_int(ncid,i,"initial_year",
     |          sh%initial_year)
              istat=nf_get_att_int(ncid,i,"initial_day",
     |          sh%initial_day)
              istat=nf_get_att_int(ncid,i,"initial_mtime",
     |          sh%initial_mtime)
            endif
!
          case('mtime') ! already got it from hist search above
          case('year')
            istat = nf_get_var1_int(ncid,i,itime,iscalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading year')
            sh%year = iscalar
          case('day')
            istat = nf_get_var1_int(ncid,i,itime,iscalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading day')
            sh%day = iscalar
          case('iter')
            istat = nf_get_var1_int(ncid,i,itime,iscalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading iter')
            sh%iter = iscalar
          case('ntask_mpi')
            istat = nf_get_var1_int(ncid,i,itime,iscalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading ntask_mpi')
            sh%ntask_mpi = iscalar
          case('coupled_cmit')
            istat = nf_get_var1_int(ncid,i,itime,iscalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading coupled_cmit')
            sh%coupled_cmit = iscalar
          case('ut')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading ut')
            sh%ut = scalar
          case('dtide')
            istart2(1) = 1
            istart2(2) = sh%ihist
            icount2(1) = 2 
            icount2(2) = 1
            istat = nf_get_vara_double(ncid,i,istart2,icount2,var2)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading dtide')
            sh%dtide(:) = var2(:)
          case('sdtide')
            istart2(1) = 1
            istart2(2) = sh%ihist
            icount2(1) = 10
            icount2(2) = 1
            istat = nf_get_vara_double(ncid,i,istart2,icount2,var10)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading sdtide')
            sh%sdtide(:) = var10(:)
          case('f107d')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading f107d')
            sh%f107d = scalar
          case('f107a')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading f107a')
            sh%f107a = scalar
          case('hpower')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading hpower')
            sh%hpower = scalar
          case('ctpoten')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading ctpoten')
            sh%ctpoten = scalar
          case('Kp')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading Kp')
            sh%kp = scalar
          case('bximf')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading bximf')
            sh%bximf = scalar
          case('byimf')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading byimf')
            sh%byimf = scalar
          case('bzimf')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading bzimf')
            sh%bzimf = scalar
          case('swvel')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading swvel')
            sh%swvel = scalar
          case('swden')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading swden')
            sh%swden = scalar
          case('al')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading al')
            sh%al = scalar
          case('colfac')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading colfac')
            sh%colfac = scalar
          case('e1')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) call handle_ncerr(istat,'reading e1')
            sh%e1 = scalar
          case('e2')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) call handle_ncerr(istat,'reading e2')
            sh%e2 = scalar
          case('h1')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) call handle_ncerr(istat,'reading h1')
            sh%h1 = scalar
          case('h2')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) call handle_ncerr(istat,'reading h2')
            sh%h2 = scalar
          case('ec')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) call handle_ncerr(istat,'reading ec')
            sh%ec = scalar
          case('ed')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) call handle_ncerr(istat,'reading ed')
            sh%ed = scalar
          case('alfac')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading alfac')
            sh%alfac = scalar
          case('alfad')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading alfad')
            sh%alfad = scalar
          case('joulefac')
            istat = nf_get_var1_double(ncid,i,itime,scalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading joulefac')
            sh%joulefac = scalar
          case('p0')
            istat = nf_get_var_double(ncid,i,var1)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading p0')
            sh%p0 = var1(1)
          case('p0_model')
            istat = nf_get_var_double(ncid,i,var1)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading p0_model')
            sh%p0_model = var1(1)
          case('grav')
            istat = nf_get_var_double(ncid,i,var1)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading grav')
            sh%grav = var1(1)
          case('timestep')
            istat = nf_get_var_int(ncid,i,iscalar)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading timestep')
            sh%step = iscalar
!
! Read file names on current itime history: 
! File name vars are type char, dimensioned (filelen,ntime)
! See subroutine rdfilename(ncid,idv,vname,filelen,filename)
!
! GPI data file:
          case('gpi_ncfile')
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%gpi_ncfile)
!
! SABER data file:
          case('saber_ncfile')
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%saber_ncfile)
!
! TIDI data file:
          case('tidi_ncfile')
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%tidi_ncfile)
!
! Background data file:
          case('bgrddata_ncfile')
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%bgrddata_ncfile)
!
! CTMT data file:
          case('ctmt_ncfile')
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%ctmt_ncfile)
!
! NCEP data file: ncep_ncfile(time,filelen), where filelen=mxlen_filename
! (ncep files used only by timegcm -- will always be '[none]' for tiegcm)
          case('ncep_ncfile')                               ! new histories
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%ncep_ncfile)
!
! SEE data file: see_ncfile(time,filelen), where filelen=mxlen_filename
          case('see_ncfile')                               ! new histories
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%see_ncfile)
!
! IMF data file: imf_ncfile(time,filelen), where filelen=mxlen_filename
          case('imf_ncfile')                               ! new histories
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%imf_ncfile)
!
! GSWM data files (from namelist input):
!   character(len=mxlen_filename) ::
!   gswm_mi_di_ncfile, ! GSWM migrating diurnal data file
!   gswm_mi_sdi_ncfile,! GSWM migrating semi-diurnal data file
!   gswm_nm_di_ncfile, ! GSWM non-migrating diurnal data file
!   gswm_nm_sdi_ncfile,! GSWM non-migrating semi-diurnal data file
!
! gswm_mi_di_ncfile:
          case('gswm_mi_di_ncfile')
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%gswm_mi_di_ncfile)
!
! gswm_mi_sdi_ncfile:
          case('gswm_mi_sdi_ncfile')
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%gswm_mi_sdi_ncfile)
!
! gswm_nm_di_ncfile:
          case('gswm_nm_di_ncfile')
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%gswm_nm_di_ncfile)
!
! gswm_nm_sdi_ncfile:
          case('gswm_nm_sdi_ncfile')
            call rdfilename(ncid,itime,i,trim(varname),filelen,
     |        sh%gswm_nm_sdi_ncfile)

          case('maglon')
          case('maglat')
!
! Interface level of lower boundary (scalar coord for TLBC,ULBC,VLBC):
          case('LBC')
            istat = nf_get_var1_double(ncid,i,itime,varlbc)
            if (istat /= NF_NOERR) 
     |        call handle_ncerr(istat,'reading LBC')
            if (varlbc /= zibot) then
              write(6,"(/,'>>> WARNING nc_rdhist read LBC: varlbc=',
     |          f8.2,' but should equal zibot=',f8.2)") varlbc,zibot
            else
              sh%lbc = varlbc
              write(6,"('Read LBC from source history: i=',i3,
     |          ' sh%lbc=',f8.2)") i,sh%lbc
            endif
!
! Read lbc's for t,u,v. See fields.F. These will exist on source history 
!   only if sh%tuv_lbc_intop==0, i.e., "new histories"). Read into global 
!   arrays, then store subdomains.
!
          case('TLBC')
            call rdlbc(ncid,i,tlbc,varname,itime)
          case('ULBC')
            call rdlbc(ncid,i,ulbc,varname,itime)
          case('VLBC')
            call rdlbc(ncid,i,vlbc,varname,itime)
          case('TLBC_NM')
            call rdlbc(ncid,i,tlbc_nm,varname,itime)
          case('ULBC_NM')
            call rdlbc(ncid,i,ulbc_nm,varname,itime)
          case('VLBC_NM')
            call rdlbc(ncid,i,vlbc_nm,varname,itime)
!
! Check for variable name in 4-d fields:
          case default
            if (any(f4d(1:nf4d_hist)%short_name==varname)) then
              do j=1,nf4d_hist
                if (trim(varname)==f4d(j)%short_name) then
                  call rdf4d(ncid,varname,itime,j)
                  sh%nflds = sh%nflds+1
                  varnames(sh%nflds) = varname
                endif
              enddo
!
! 10/27/05 btf: This conditional covers "W" from old histories.
! (var name "W" has been changed to "OMEGA" on the new histories)
!
            elseif (varname=='W') then
              do j=1,nf4d_hist
                if (f4d(j)%short_name=='OMEGA') then
                  call rdf4d(ncid,varname,itime,j)
                  sh%nflds = sh%nflds+1
                  varnames(sh%nflds) = varname
                endif
              enddo
!
! Variable is unknown to the model:
            else
              write(6,"('Note nc_rdhist: unused variable: ',a)") 
     |          trim(varname)
            endif
        end select
      enddo ! i=1,nvars
!
      sh%zmtop = zmtop
      sh%zmbot = zmbot
      sh%zitop = zitop
      sh%zibot = zibot
      if (associated(sh%fnames)) deallocate(sh%fnames)
      if (sh%nflds > 0) then
        allocate(sh%fnames(sh%nflds),stat=ier)
        if (ier /= 0) then
          write(6,"('>>> nc_rdhist: error allocating sh%fnames with ',
     |      ' sh%nflds=',i4)") sh%nflds
          call shutdown('source hist nflds')
        endif
      else
        write(6,"('>>> sh%nflds=',i4)") sh%nflds
        call shutdown('sh%nflds')
      endif
      do i=1,sh%nflds
        sh%fnames(i) = varnames(i)
      enddo
!
! Report to stdout:
      call hist_print(sh,'READ',diskfile)
      end subroutine nc_rdhist
!-----------------------------------------------------------------------
      subroutine rdlbc(ncid,idv,flbc,vname,itime)
!
! Read global TLBC, ULBC, or VLBC from current history, and transfer to
!  subdomains tlbc, ulbc, or vlbc (fields.F).
!
      use mpi_module,only: lon0,lon1,lat0,lat1
      use fields_module,only: lond0,lond1,latd0,latd1
! Routine to add fields to secondary histories:
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: ncid,idv,itime
      character(len=*),intent(in) :: vname
      real,intent(out) :: flbc(lond0:lond1,latd0:latd1)
!
! Local:
      integer :: istat,i,j,lonbeg,lonend
      real,dimension(nlon,nlat) :: glbc ! lbc at global domain
      character(len=80) :: char80
!
      start_3d(1:2) = 1
      start_3d(3) = itime
      count_3d(1) = nlon
      count_3d(2) = nlat
      count_3d(3) = 1
!
! Read global domain into local glbc:
      istat = nf_get_vara_double(ncid,idv,start_3d,count_3d,glbc)
      if (istat /= NF_NOERR) then
        write(char80,"('Error reading ',a,' ncid=',i3)") 
     |    trim(vname),ncid
        call handle_ncerr(istat,trim(char80))
      endif
!     write(6,"('rdlbc read ',a,': glbc min,max=',2e12.4)")
!    |  trim(vname),minval(glbc),maxval(glbc)
!
! Store subdomain in flbc (actual arg is tlbc, ulbc, or vlbc):
! (Zero-out periodic points)
      lonbeg = lon0
      if (lon0==1) then
        lonbeg = 3
        flbc(1:2,:) = 0.
      endif
      lonend = lon1
      if (lon1==nlonp4) then
        lonend = lon1-2
        flbc(lon1-1:lon1,:) = 0.
      endif
      do j=lat0,lat1
        do i=lonbeg,lonend
          flbc(i,j) = glbc(i-2,j)
        enddo
      enddo
      
#ifndef MPI
! Periodic points for non-MPI runs:
      flbc(1:2,:)              = flbc(nlonp4-3:nlonp4-2,:)
      flbc(nlonp4-1:nlonp4,:)  = flbc(3:4,:)
#endif
 
!     write(6,"('Read ',a,' subdomain from source history: min,max=',
!    |  2e12.4)")  trim(vname), minval(flbc(lon0:lon1,lat0:lat1)),
!    |                          maxval(flbc(lon0:lon1,lat0:lat1))

!     if (trim(vname)=='TLBC') then
!       do i=lon0,lon1
!         write(6,"('rdlbc: i=',i3,' lat0,1=',2i3,' tlbc(i,lat0:lat1)=',
!    |      /,(6e12.4))") i,lat0,lat1,tlbc(i,lat0:lat1)
!       enddo ! j=lat0,lat1
!     endif
!     if (trim(vname)=='ULBC') then
!       do i=lon0,lon1
!         write(6,"('rdlbc: i=',i3,' lat0,1=',2i3,' ulbc(i,lat0:lat1)=',
!    |      /,(6e12.4))") i,lat0,lat1,ulbc(i,lat0:lat1)
!       enddo ! j=lat0,lat1
!     endif
!     if (trim(vname)=='VLBC') then
!       do i=lon0,lon1
!         write(6,"('rdlbc: i=',i3,' lat0,1=',2i3,' vlbc(i,lat0:lat1)=',
!    |      /,(6e12.4))") i,lat0,lat1,vlbc(i,lat0:lat1)
!       enddo ! j=lat0,lat1
!     endif
 
      end subroutine rdlbc
!-----------------------------------------------------------------------
      subroutine rdf4d(ncid,name,itime,ix)
!
! Read prognostic field "name" from netcdf file on ncid into into local
! f3drd at full domain. Then subdomains are transferred from f3rd into
! 4-d arrays for use in the model. It would be interesting to try reading
! only the subdomain directly from the file rather than the full domain.
!
      use input_module,only: dynamo,calc_helium
      use hist_module,only: sh
      use mpi_module,only: lon0,lon1,lat0,lat1
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: ncid,itime,ix
!
! Local:
      integer :: i,ii,j,k,istat,idvar,itype,ndims,iddims(4),natts,
     |  iscalar,idimsizes(4),nx,nxk,lonbeg,lonend,nnans,nanfatal
      character(len=16) :: rdname,tmpname
      character(len=120) :: char120
      real :: fakevar(1,1,1),fmin,fmax,fmin1,fmax1
      real :: f3drd(nlon,nlat,nlevp1) ! (i,j,k) as on history 
      real :: zmj2,zmj1,fac
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
        if (idimsizes(1) /= nlon .or. idimsizes(2) /= nlat .or.
     |      idimsizes(3) /= nlevp1) then
          write(6,"(/,'>>> WARNING rdf4d: bad dimension sizes',
     |      ' for field ',a)") trim(rdname)
          write(6,"('  idimsizes(1:3)=',3i4,' but should be ',
     |      'nlon,nlat,nlevp1=',3i4)") idimsizes(1:3),nlon,nlat,
     |      nlevp1
          call shutdown('source history coordinate dimensions')
        endif
      endif
!
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
! Zero out helium if requested by namelist calc_helium=0:
      if (trim(name)=='HE'.and.calc_helium <= 0) then
        f4d(ix)%data = 0.
        write(6,"('rdf4d: namelist calc_helium=',i3,' -> am setting',
     |    ' prognostic helium (HE) to zero.')") calc_helium
        return
      endif
      if (trim(name)=='HE_NM'.and.calc_helium <= 0) then
        f4d(ix)%data = 0.
        write(6,"('rdf4d: namelist calc_helium=',i3,' -> am setting',
     |    ' helium at time t-1 (HE_NM) to zero.')") calc_helium
        return
      endif
!
! Read field into local f3drd:
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
! Check for NaNs in the source history:
!     subroutine check_nans(f,id1,id2,id3,name,n_total,ispval,spval,
!    |  iprint,ifatal)
! real :: f3drd(nlon,nlat,nlevp1) ! (i,j,k) as on history 
!
      nanfatal = 0 ! reading nans from source history is non-fatal
      if (check_nan)
     |  call check_nans(f3drd,nlon,nlat,nlevp1,name(1:8),nnans,0,0.,1,
     |    nanfatal)
!
! Transform from f3drd(i,j,k) to f4d(k,i,j):
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2
      do j=lat0,lat1
        do i=lonbeg,lonend
          do k=1,nlevp1
!
! Local : real :: f3drd(nlon,nlat,nlevp1) ! (i,j,k) as on history 
! Global: real,pointer :: f4d(ix)%data(levd0:levd1,lond0:lond1,latd0:latd1) 
!
            f4d(ix)%data(k,i,j,itp) = f3drd(i-2,j,k)
          enddo ! k=1,nlevp1
        enddo ! i=lonbeg,lonend  
      enddo ! j=lat0,lat1
!
! Define t,u,v lbc subdomains from top slot if tuv_lbc_intop==1: 
! If tuv_lbc_intop==0, these arrays were read by nc_rdhist.
! These arrays are in fields.F, allocated by sub init_lbc as
!   (lond0:lond1,latd0:latd1)
! 
!
      if (sh%tuv_lbc_intop > 0) then ! (is "old" history)
        select case(trim(name))
          case('TN')
            do j=lat0,lat1
              tlbc(:,j) = f4d(ix)%data(nlevp1,:,j,itp)
            enddo ! j=lat0,lat1
            write(6,"('rdf4d: tlbc from top slot min,max=',2f10.2)")
     |        minval(tlbc(lon0:lon1,lat0:lat1)),
     |        maxval(tlbc(lon0:lon1,lat0:lat1))
            f4d(ix)%data(nlevp1,:,:,itp) = spval ! top slot is now missing data
          case('UN')
            do j=lat0,lat1
              ulbc(:,j) = f4d(ix)%data(nlevp1,:,j,itp)
            enddo
            write(6,"('rdf4d: ulbc from top slot min,max=',2f10.2)")
     |        minval(ulbc(lon0:lon1,lat0:lat1)),
     |        maxval(ulbc(lon0:lon1,lat0:lat1))
            f4d(ix)%data(nlevp1,:,:,itp) = spval ! top slot is now missing data
          case('VN')
            do j=lat0,lat1
              vlbc(:,j) = f4d(ix)%data(nlevp1,:,j,itp)
            enddo
            write(6,"('rdf4d: vlbc from top slot min,max=',2f10.2)")
     |        minval(vlbc(lon0:lon1,lat0:lat1)),
     |        maxval(vlbc(lon0:lon1,lat0:lat1))
            f4d(ix)%data(nlevp1,:,:,itp) = spval ! top slot is now missing data
	    
          case('TN_NM')
            do j=lat0,lat1
              tlbc_nm(:,j) = f4d(ix)%data(nlevp1,:,j,itp)
            enddo ! j=lat0,lat1
            write(6,"('rdf4d: tlbc_nm from top slot min,max=',2f10.2)")
     |        minval(tlbc_nm(lon0:lon1,lat0:lat1)),
     |        maxval(tlbc_nm(lon0:lon1,lat0:lat1))
            f4d(ix)%data(nlevp1,:,:,itp) = spval ! top slot is now missing data
          case('UN_NM')
            do j=lat0,lat1
              ulbc_nm(:,j) = f4d(ix)%data(nlevp1,:,j,itp)
            enddo
            write(6,"('rdf4d: ulbc_nm from top slot min,max=',2f10.2)")
     |        minval(ulbc_nm(lon0:lon1,lat0:lat1)),
     |        maxval(ulbc_nm(lon0:lon1,lat0:lat1))
            f4d(ix)%data(nlevp1,:,:,itp) = spval ! top slot is now missing data
          case('VN_NM')
            do j=lat0,lat1
              vlbc_nm(:,j) = f4d(ix)%data(nlevp1,:,j,itp)
            enddo
            write(6,"('rdf4d: vlbc_nm from top slot min,max=',2f10.2)")
     |        minval(vlbc_nm(lon0:lon1,lat0:lat1)),
     |        maxval(vlbc_nm(lon0:lon1,lat0:lat1))
            f4d(ix)%data(nlevp1,:,:,itp) = spval ! top slot is now missing data
          case default ! do nothing if not one of these fields
        end select ! case(trim(name))
      endif ! if (sh%tuv_lbc_intop > 0)
!
! Report to stdout:
!
      call fminmaxspv(f4d(ix)%data(:,lonbeg:lonend,lat0:lat1,itp),
     |  nlevp1*(lonend-lonbeg+1)*(lat1-lat0+1),fmin,fmax,spval)
      write(6,"('Read field ',a,' 3d subdomain min,max=',2e12.4)") 
     |  f4d(ix)%short_name,fmin,fmax
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
      use hist_module,only: h,sh,nsource,nsecsource
      use cons_module,only: pi
      use input_module,only: potential_model,ncep_ncfile
      use params_module,only: tgcm_name,tgcm_version
!
! Args:
      integer,intent(in) :: ncid
!
! Local:
      integer :: i,ii,istat,idum,ivar1(1),imo,ida,startmtime(4),nbytes
      character(len=120) :: char120
      character(len=80) :: char80
      character(len=24) :: create_date,create_time
      real :: var1(1),rdays,rmins
      real,external :: mtime_to_datestr
!
! Define dimensions (field vars will be (time,lev,lat,lon) in
! CDL notation. Mag fields will be (time,[mlev or imlev],mlat,mlon): 
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
!
! Midpoint levels dimension "lev"=nlevp1:
      istat = nf_def_dim(ncid,"lev",nlevp1,id_lev)     ! midpoint levels dimension
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining midpoint levels dimension')
!
! Interface levels dimension "ilev"=nlevp1:
      istat = nf_def_dim(ncid,"ilev",nlevp1,id_ilev)   ! interface levels dimension
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining interface levels dimension')
!
! Magnetic grid dimensions:
      istat = nf_def_dim(ncid,"mlon",nmlonp1,id_mlon)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining magnetic longitude dimension')
      istat = nf_def_dim(ncid,"mlat",nmlat,id_mlat)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error defining magnetic latitude dimension')
!
! Mag midpoint levels dimension:
      istat = nf_def_dim(ncid,"mlev",nmlevp1,id_mlev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining magnetic midpoints levels dimension')
!
! Mag interface levels dimension:
      istat = nf_def_dim(ncid,"imlev",nimlevp1,id_imlev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining magnetic interface levels dimension')

      istat = nf_def_dim(ncid,"mtimedim",3,id_mtimedim)
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
! String length of date and time:
      istat = nf_def_dim(ncid,"datelen",len(h%write_date),id_datelen)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining datelen dimension')
!
! String length of file names (mxlen_filename is parameter in input.F):
      istat = nf_def_dim(ncid,"filelen",mxlen_filename,id_filelen)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining filelen dimension')
!
! Define dimension (coordinate) variables and attributes:
!
! Time (coordinate variable time(time)). This is days since
! the initial run's start time. The units string is: yyyy-m-d, 
! where yyyy is the year, m is month, and d is day of the source
! start time.
!
      ids1(1) = id_time	! for 1d time-unlimited vars
      istat = nf_def_var(ncid,"time",NF_DOUBLE,1,ids1,idv_time)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining time dimension variable')
      startmtime(1:3) = sh%initial_mtime(1:3) ; startmtime(4) = 0
!     rdays = mtime_to_datestr(sh%initial_year,startmtime,imo,ida,
!    |  char80)
      rmins = mtime_to_datestr(sh%initial_year,startmtime,imo,ida,
     |  char80)
      istat = nf_put_att_text(ncid,idv_time,"long_name",4,"time")
      istat = nf_put_att_text(ncid,idv_time,"units",len_trim(char80),
     |  trim(char80))
      istat = nf_put_att_int(ncid,idv_time,"initial_year",NF_INT,1,
     |  h%initial_year)
      istat = nf_put_att_int(ncid,idv_time,"initial_day",NF_INT,1,
     |  h%initial_day)
      istat = nf_put_att_int(ncid,idv_time,"initial_mtime",NF_INT,3,
     |  h%initial_mtime)
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
! Midpoint levels coordinate variable lev(lev):
      ids1(1) = id_lev
      istat = nf_def_var(ncid,"lev",NF_DOUBLE,1,ids1,idv_lev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining midpoint levels coordinate variable lev(lev)')
! long name of lev coord var:
      write(char80,"('midpoint levels')")
      istat = nf_put_att_text(ncid,idv_lev,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |'Error defining long_name of midpoint levels coordinate variable')
! short name = ln(p0/p)
      istat = nf_put_att_text(ncid,idv_lev,"short_name",
     |  8,"ln(p0/p)")
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining short_name of midpoint levels coord var')
! lev coord var is unitless:
      istat = nf_put_att_text(ncid,idv_lev,"units",0," ")
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining units of midpoint levels coord var')
! positive='up'
      istat = nf_put_att_text(ncid,idv_lev,"positive",2,'up')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining positive attribute of midpoint levels coord var')
! standard name of midpoints lev coord var:
      write(char80,"('atmosphere_ln_pressure_coordinate')")
      istat = nf_put_att_text(ncid,idv_lev,"standard_name",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining standard_name of midpoint levels coord var')
! formula terms for midpoints lev coord var:
      write(char80,"('p0: p0 lev: lev')")
      istat = nf_put_att_text(ncid,idv_lev,"formula_terms",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining formula terms of midpoint levels coord var')
! formula to obtain pressure from midpoint lev:
      write(char80,"('p(k) = p0 * exp(-lev(k))')")
      istat = nf_put_att_text(ncid,idv_lev,"formula",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining formula for midpoint levels coord var')
!
! Interface levels coordinate array:
      ids1(1) = id_ilev
      istat = nf_def_var(ncid,"ilev",NF_DOUBLE,1,ids1,idv_ilev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining interface levels coordinate variable ilev(ilev)')
! long name of ilev coord var:
      write(char80,"('interface levels')")
      istat = nf_put_att_text(ncid,idv_ilev,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |'Error defining long_name of interface levels coordinate var')
! short name = ln(p0/p)
      istat = nf_put_att_text(ncid,idv_ilev,"short_name",
     |  8,"ln(p0/p)")
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining short_name of interface levels coord var')
! ilev coord var is unitless:
      istat = nf_put_att_text(ncid,idv_ilev,"units",0," ")
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining units of interface levels coord var')
! positive='up'
      istat = nf_put_att_text(ncid,idv_ilev,"positive",2,'up')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |'Error defining positive attribute of interface levels coord var')
! standard name of interface ilev coord var:
      write(char80,"('atmosphere_ln_pressure_coordinate')")
      istat = nf_put_att_text(ncid,idv_ilev,"standard_name",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining standard_name of interface levels coord var')
! formula terms for interface ilev coord var:
      write(char80,"('p0: p0 lev: ilev')")
      istat = nf_put_att_text(ncid,idv_ilev,"formula_terms",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining formula terms of interface levels coord var')
! formula to obtain pressure from interface ilev:
      write(char80,"('p(k) = p0 * exp(-ilev(k))')")
      istat = nf_put_att_text(ncid,idv_ilev,"formula",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining formula for interface levels coord var')
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
      istat = nf_put_att_text(ncid,idv_mlon,"units",12,'degrees_east')
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
! Magnetic midpoint coordinate variable mlev(mlev)):
      ids1(1) = id_mlev
      istat = nf_def_var(ncid,"mlev",NF_DOUBLE,1,ids1,idv_mlev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining magnetic levels coord variable mlev')
! long name of mlev:
      write(char80,"('magnetic midpoint levels')")
      istat = nf_put_att_text(ncid,idv_mlev,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of mlev coordinate variable')
! short name = ln(p0/p)
      istat = nf_put_att_text(ncid,idv_mlev,"short_name",
     |  8,"ln(p0/p)")
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining short_name of mlev levels coord var')
! units of mlev (unitless):
      istat = nf_put_att_text(ncid,idv_mlev,"units",0," ")
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of mlev coordinate variable')
! positive attribute of mlev:
      istat = nf_put_att_text(ncid,idv_mlev,"positive",2,'up')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining positive attribute for mlev coord var')
! standard_name of mlev:
      write(char80,"('atmosphere_ln_pressure_coordinate')")
      istat = nf_put_att_text(ncid,idv_mlev,"standard_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining standard_name attribute for mlev coord var')
! formula terms for mlev:
      write(char80,"('p0: p0 lev: mlev')")
      istat = nf_put_att_text(ncid,idv_mlev,"formula_terms",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining formula terms of mlev coord var')
! formula to obtain pressure from mlev:
      write(char80,"('p(k) = p0 * exp(-mlev(k))')")
      istat = nf_put_att_text(ncid,idv_mlev,"formula",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining formula for mlev coord var')
!
! Magnetic interface coordinate variable imlev(imlev)):
      ids1(1) = id_imlev
      istat = nf_def_var(ncid,"imlev",NF_DOUBLE,1,ids1,idv_imlev)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining magnetic levels coord variable imlev')
! long name of imlev:
      write(char80,"('magnetic interface levels')")
      istat = nf_put_att_text(ncid,idv_imlev,"long_name",
     +  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of imlev coordinate variable')
! short name = ln(p0/p)
      istat = nf_put_att_text(ncid,idv_imlev,"short_name",
     |  8,"ln(p0/p)")
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining short_name of imlev levels coord var')
! units of imlev (unitless):
      istat = nf_put_att_text(ncid,idv_imlev,"units",0," ")
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining units of imlev coordinate variable')
! positive attribute of imlev:
      istat = nf_put_att_text(ncid,idv_imlev,"positive",2,'up')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining positive attribute for imlev coord var')
! standard_name of imlev:
      write(char80,"('atmosphere_ln_pressure_coordinate')")
      istat = nf_put_att_text(ncid,idv_imlev,"standard_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining standard_name attribute for imlev coord var')
! formula terms for imlev:
      write(char80,"('p0: p0 lev: imlev')")
      istat = nf_put_att_text(ncid,idv_imlev,"formula_terms",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining formula terms of imlev coord var')
! formula to obtain pressure from imlev:
      write(char80,"('p(k) = p0 * exp(-imlev(k))')")
      istat = nf_put_att_text(ncid,idv_imlev,"formula",
     |  len_trim(char80),trim(char80)) 
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     | 'Error defining formula for imlev coord var')
!
! ids1(1) is id of time dimension for the following defines:
      ids1(1) = id_time	! for 1d time-unlimited vars
!
! Source date and time are attributes of time var 
! (they are also global file attributes):
      istat = nf_put_att_int(ncid,idv_time,"initial_year",NF_INT,1,
     |  sh%initial_year)
      istat = nf_put_att_int(ncid,idv_time,"initial_day",NF_INT,1,
     |  sh%initial_day)
      istat = nf_put_att_int(ncid,idv_time,"initial_mtime",NF_INT,3,
     |  sh%initial_mtime)
!
! Define "primary" (non-coordinate) time dependent variables:
! (these variables get a value for each history) 
!
! Model time variable mtime(3,time), where 3 is the mtime dimension 
! (id_mtimedim) (day,hour,minute), and time is the unlimited dimension 
! (id_time):
!
      ids2(1) = id_mtimedim
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
! Calendar advance:
      istat = nf_def_var(ncid,"calendar_advance",NF_INT,1,ids1,
     |  idv_calendar_adv)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining integer calendar advance variable')
      write(char80,"('calendar advance flag',
     |  ' (1 if advancing calendar time)')")
      istat = nf_put_att_text(ncid,idv_calendar_adv,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining long_name of calendar advance variable')
!
! Date each history was written:
      ids2(1) = id_datelen
      ids2(2) = id_time
      istat= nf_def_var(ncid,"write_date",NF_CHAR,2,ids2,idv_writedate)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining history write_date variable')
      write(char80,"('Date and time each history was written')")
      istat = nf_put_att_text(ncid,idv_writedate,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of write_date')
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
! Number of MPI tasks:
      istat = nf_def_var(ncid,"ntask_mpi",NF_INT,1,ids1,idv_ntask)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining ntask_mpi variable')
      istat = nf_put_att_text(ncid,idv_ntask,"long_name",
     |  19,'Number of MPI tasks')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name of ntask_mpi variable')
!
! If this run was coupled w/ CISM/CMIT:
      istat = nf_def_var(ncid,"coupled_cmit",NF_INT,1,ids1,
     |  idv_coupled_cmit)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining coupled_cmit variable')
      istat = nf_put_att_text(ncid,idv_coupled_cmit,"long_name",
     |  40,'1 if coupled with CISM/CMIT, 0 otherwise')
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining long_name of coupled_cmit variable')
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
     |  "Daily 10.7 cm solar flux","1.e-22 W/m2/Hz")
      call defvar_time_dbl(ncid,"f107a",ids1,idv_f107a,
     |  "81-day average 10.7 cm solar flux","1.e-22 W/m2/Hz")
!
! Hemispheric power: 
      call defvar_time_dbl(ncid,"hpower",ids1,idv_hpower,
     |  "Hemispheric power","GW")
!
! Cross-tail potential:
      call defvar_time_dbl(ncid,"ctpoten",ids1,idv_ctpoten,
     |  "Cross-tail potential","kV")
!
! Kp index (data_source will be either namelist or gpi_ncfile):
      call defvar_time_dbl(ncid,"Kp",ids1,idv_kp,"Kp index"," ")
      istat = nf_put_att_text(ncid,idv_kp,"data_source",
     |  len_trim(datasource_kp),datasource_kp)
!
! BX,BY,BZ:
      call defvar_time_dbl(ncid,"bximf",ids1,idv_bximf,
     |  "BX component of IMF","nT")
      call defvar_time_dbl(ncid,"byimf",ids1,idv_byimf,
     |  "BY component of IMF","nT")
      call defvar_time_dbl(ncid,"bzimf",ids1,idv_bzimf,
     |  "BZ component of IMF","nT")
!
! swvel,swden:
      call defvar_time_dbl(ncid,"swvel",ids1,idv_swvel,
     |  "Solar wind velocity","km/s")
      call defvar_time_dbl(ncid,"swden",ids1,idv_swden,
     |  "Solar wind density","cm3")
!
! al:
      call defvar_time_dbl(ncid,"al",ids1,idv_al,
     |  "Lower magnetic auroral activity index","nT")
!
! Collision factor:
      call defvar_time_dbl(ncid,"colfac",ids1,idv_colfac,
     |  "ion/neutral collision factor"," ")
!
! Joule heating factor:
      call defvar_time_dbl(ncid,"joulefac",ids1,idv_joulefac,
     |  "joule heating factor"," ")
!
! 1/31/01: Making tides time-dependent (dtide(time,2) and sdtide(time,10):
! Diurnal tide amp/phase:
      ids2(1) = id_dtide
      ids2(2) = id_time
      istat = nf_def_var(ncid,"dtide",NF_DOUBLE,2,ids2,idv_dtide)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining variable dtide')
      write(char80,"('amplitude and phase of diurnal tide mode (1,1)')")
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
      write(char80,"('amplitudes and phases of semi-diurnal tide')")
      istat = nf_put_att_text(ncid,idv_sdtide,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining long_name attribute of variable sdtide')
!
! Define ncep data file (used by timegcm only):
      call define_ncfile(ncid,"ncep_ncfile",idv_ncep_ncfile,
     |  "long_name","ncep data file",
     |  "data","10 mb lower boundary for TN and Z")
!
! Define gpi data file:
      call define_ncfile(ncid,"gpi_ncfile",idv_gpi_ncfile,
     |  "long_name","GeoPhysical Indices data file",
     |  "data","Kp, f107, f107a, ctpoten, power") 
!
! Define saber data file:
      call define_ncfile(ncid,"saber_ncfile",idv_saber_ncfile,
     |  "long_name","SABER lbc data file",
     |  "data","T,Z") 
!
! Define tidi data file:
      call define_ncfile(ncid,"tidi_ncfile",idv_tidi_ncfile,
     |  "long_name","TIDI lbc data file",
     |  "data","U,V") 
!
! Define background data file:
      call define_ncfile(ncid,"bgrddata_ncfile",idv_bgrddata_ncfile,
     |  "long_name","background lbc data file",
     |  "data","T,U,Z")
!
! Define ctmt data file:
      call define_ncfile(ncid,"ctmt_ncfile",idv_ctmt_ncfile,
     |  "long_name","CTMT lbc data file",
     |  "data","U,V, T, Z")
!
! Define TIMED SEE data file:
      call define_ncfile(ncid,"see_ncfile",idv_see_ncfile,
     |  "long_name","TIMED SEE data file",
     |  "data","Solar flux data")
!
! Define imf data file:
      call define_ncfile(ncid,"imf_ncfile",idv_imf_ncfile,
     |  "long_name","IMF data file",
     |  "data","Bx,By,Bz,Swvel,Swden,Kp,f107,f107a")
!
! Define gswm data files:
      call define_ncfile(ncid,"gswm_mi_di_ncfile",
     |  idv_gswm_mi_di_ncfile,
     | "long_name","GSWM migrating diurnal tides data file",
     | "data","lower boundary perturbations for TN, UN, VN, and Z")
      call define_ncfile(ncid,"gswm_mi_sdi_ncfile",
     |  idv_gswm_mi_sdi_ncfile,
     | "long_name","GSWM migrating semi-diurnal tides data file",
     | "data","lower boundary perturbations for TN, UN, VN, and Z")
      call define_ncfile(ncid,"gswm_nm_di_ncfile",
     |  idv_gswm_nm_di_ncfile,
     | "long_name","GSWM non-migrating diurnal tides data file",
     | "data","lower boundary perturbations for TN, UN, VN, and Z")
      call define_ncfile(ncid,"gswm_nm_sdi_ncfile",
     |  idv_gswm_nm_sdi_ncfile,
     | "long_name","GSWM non-migrating semi-diurnal tides data file",
     | "data","lower boundary perturbations for TN, UN, VN, and Z")
!
! Auroral parameters (alfa30, e30, alfad2, ed2 are timegcm only)
!
! alfa30:
!     call defvar_time_dbl(ncid,"alfa30",ids1,idv_alfa30,
!    |  "Characteristic energy of high-energy auroral electrons",
!    |  "KeV")
!
! e30:
!     call defvar_time_dbl(ncid,"e30",ids1,idv_e30,
!    |  "Energy flux of high-energy auroral electrons",
!    |  "ergs/cm2/s")
!
! alfad2:
!     call defvar_time_dbl(ncid,"alfad2",ids1,idv_alfad2,
!    |  "Characteristic energy of solar protons in the polar cap",
!    |  "KeV")
!
! ed2:
!     call defvar_time_dbl(ncid,"ed2",ids1,idv_ed2,
!    |  "Energy flux of solar protons in the polar cap",
!    |  "ergs/cm2/s")
!
! 1/24/08 btf: Auroral parameters:
!
! e1:
      call defvar_time_dbl(ncid,"e1",ids1,idv_e1,
     |  "Peak energy flux in noon sector of the aurora",
     |  "ergs/cm2/s")
! e2:
      call defvar_time_dbl(ncid,"e2",ids1,idv_e2,
     |  "Peak energy flux in midnight sector of the aurora",
     |  "ergs/cm2/s")
! h1:
      call defvar_time_dbl(ncid,"h1",ids1,idv_h1,
     |  "Gaussian half-width of the noon auroral oval",
     |  "degrees")
! h2:
      call defvar_time_dbl(ncid,"h2",ids1,idv_h2,
     |  "Gaussian half-width of the midnight auroral oval",
     |  "degrees")
! alfac:
      call defvar_time_dbl(ncid,"alfac",ids1,idv_alfac,
     |  "Characteristic Maxwellian energy of polar cusp electrons",
     |  "keV")
! ec:
      call defvar_time_dbl(ncid,"ec",ids1,idv_ec,
     |  "Column energy input of polar cusp electrons",
     |  "ergs/cm**2/s")
! alfad:
      call defvar_time_dbl(ncid,"alfad",ids1,idv_alfad,
     |  "Characteristic Maxwellian energy of drizzle electrons",
     |  "keV")
! ed:
      call defvar_time_dbl(ncid,"ed",ids1,idv_ed,
     |  "Column energy input of drizzle electrons",
     |  "ergs/cm**2/s")
! crit(1):
      call defvar_time_dbl(ncid,"crit1",ids1,idv_crit1,
     |  "critical cross-over latitude 1","deg")
! crit(2):
      call defvar_time_dbl(ncid,"crit2",ids1,idv_crit2,
     |  "critical cross-over latitude 2","deg")
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
! Reference pressure used to convert ln(p0/p) to pressure (millibars (hPa)):
      istat = nf_def_var(ncid,"p0",NF_DOUBLE,0,idum,idv_p0)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining variable p0')
      istat = nf_put_att_text(ncid,idv_p0,"long_name",
     |  18,"Reference pressure") ! for conversion of vertical coordinates
!
! 4/26/06 btf: p0 units can be either hPa or millibars. Using the latter
!              now to be consistent w/ microbar units of p0_model
      istat = nf_put_att_text(ncid,idv_p0,"units",9,"millibars")
!     istat = nf_put_att_text(ncid,idv_p0,"units",3,"hPA")
!
! Reference pressure used by the model (dynes/cm2):
! p0_model units can be either dynes/cm2 or microbars. Using the latter 
!   to be consistent w/ 1981 paper by Dickinson, Ridley, and Roble
!
      istat = nf_def_var(ncid,"p0_model",NF_DOUBLE,0,idum,idv_p0_model)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining variable p0')
      char80 = ' ' 
      write(char80,"('Reference pressure (as used by the model)')")
      istat = nf_put_att_text(ncid,idv_p0_model,"long_name",
     |  len_trim(char80),trim(char80))
      istat = nf_put_att_text(ncid,idv_p0_model,"units",9,"microbars")
!     istat = nf_put_att_text(ncid,idv_p0_model,"units",9,"dynes/cm2")
!
! Gravity (constant used in the model, independent of height):
      istat = nf_def_var(ncid,"grav",NF_DOUBLE,0,idum,idv_grav)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error defining variable grav')
      istat = nf_put_att_text(ncid,idv_grav,"long_name",
     |  26,"gravitational acceleration")
      istat = nf_put_att_text(ncid,idv_grav,"units",
     |  4,"cm/s")
      write(char80,"('constant used in the model, independent of ',
     |  'height')")
      istat = nf_put_att_text(ncid,idv_grav,"info",
     |  len_trim(char80),trim(char80))
!
! All prognostic fields are on primary history:
      if (h%hist_type(1:3)=='pri') then
        nbytes = 8
        do i=1,nf4d_hist
          call defvar_f4d(ncid,f4d(i),idv_f4d(i),nbytes)
        enddo
!
! Define secondary history fields (geographic or magnetic):
      elseif (h%hist_type(1:3)=='sec') then
        call def_fsech(ncid)
      else
        write(6,"('>>> nc_define: unknown h%hist_type=',a)") h%hist_type
      endif
!
! Define lower boundaries of t,u,v (both secondary and primary):
      call deflbc(ncid,"TLBC",idv_tlbc,h%hist_type)
      call deflbc(ncid,"ULBC",idv_ulbc,h%hist_type)
      call deflbc(ncid,"VLBC",idv_vlbc,h%hist_type)
      call deflbc(ncid,"TLBC_NM",idv_tlbc_nm,h%hist_type)
      call deflbc(ncid,"ULBC_NM",idv_ulbc_nm,h%hist_type)
      call deflbc(ncid,"VLBC_NM",idv_vlbc_nm,h%hist_type)
!
! Define LBC: bottom interface level (scalar coord of TLBC,ULBC,VLBC):
      istat = nf_def_var(ncid,"LBC",NF_DOUBLE,0,idum,idv_lbc)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining var LBC')
      write(char80,
     |  "('Interface level of t,u,v lower boundary condition')")
      istat = nf_put_att_text(ncid,idv_lbc,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'Error defining long_name attribute of variable LBC')
!     write(6,"('defined LBC: idv_lbc=',i3)") idv_lbc
!
! Define global file attributes:
!
! CF Conventions:
      istat = nf_put_att_text(ncid,NF_GLOBAL,"Conventions",6,"CF-1.0")
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for Conventions global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! User-defined job comment (namelist parameter "label"):
      if (len_trim(h%label) > 0) then
        istat = nf_put_att_text(ncid,NF_GLOBAL,"label",
     |    len_trim(h%label),trim(h%label))
        if (istat /= NF_NOERR) then
          write(char120,"('Error return from nf_put_att_text ',
     |      'for label global attribute')")
          call handle_ncerr(istat,char120)
        endif
      endif

! Create date and time:
      call datetime(create_date,create_time)
      h%create_date = trim(create_date)//' '//trim(create_time)
      istat = nf_put_att_text(ncid,NF_GLOBAL,"create_date",
     +  len_trim(h%create_date),trim(h%create_date))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for create_date global attribute')")
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
     +  len_trim(tgcm_name),trim(tgcm_name)) 
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for model name global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Model version string:
      istat = nf_put_att_text(ncid,NF_GLOBAL,"model_version",
     +  len_trim(tgcm_version),trim(tgcm_version)) 
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for model version global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! SVN version (see included file revision.inc, created by the build):
      write(6,"('nc_define: svn_revision = ',a)") trim(svn_revision)
      istat = nf_put_att_text(ncid,NF_GLOBAL,"svn_revision",
     +  len_trim(svn_revision),trim(svn_revision)) 
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     +    'for SVN version global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Output path (from namelist):
      istat = nf_put_att_text(ncid,NF_GLOBAL,"output_file",
     |  len_trim(h%output_file),trim(h%output_file))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for output_file global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Primary or Secondary histories (h%hist_type):
      istat = nf_put_att_text(ncid,NF_GLOBAL,"history_type",
     |  len_trim(h%hist_type),trim(h%hist_type))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for history_type global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Run type ('initial' or 'continuation')
      istat = nf_put_att_text(ncid,NF_GLOBAL,"run_type",
     |  len_trim(h%run_type),trim(h%run_type))
      if (istat /= NF_NOERR) then
        write(char120,"('nc_define: Error return from nf_put_att_text ',
     |    'for run_type global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Source history file (for this run, either initial or continuation):
      if (nsource > 0) then ! initial source file
        istat = nf_put_att_text(ncid,NF_GLOBAL,"source_file",
     |    len_trim(h%source_file)+10,trim(h%source_file)//' (initial)')
      else ! continuation source file
        istat = nf_put_att_text(ncid,NF_GLOBAL,"source_file",
     |    len_trim(h%source_file)+15,
     |    trim(h%source_file)//' (continuation)')
      endif
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for source_file global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Source model time (will be same as start time for this run,
! whether initial or continuation)
      istat = nf_put_att_int(ncid,NF_GLOBAL,"source_mtime",NF_INT,3,
     |  h%source_mtime)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for source_mtime global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Initial file (source file from the initial run)
      istat = nf_put_att_text(ncid,NF_GLOBAL,"initial_file",
     |  len_trim(h%initial_file),trim(h%initial_file))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for initial_file global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Initial mtime (this is also an attribute of the time coord var)
      istat = nf_put_att_int(ncid,NF_GLOBAL,"initial_mtime",NF_INT,3,
     |  h%initial_mtime)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for initial_mtime global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Method 1 for conversion from lev (ln(p0/p)) to pressure in hPa 
!   (millibars), using p0:
!
      write(char80,"('p0*exp(-lev(k))')")
      istat = nf_put_att_text(ncid,NF_GLOBAL,"lev_to_hPa_method1",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for lev_to_hPa_method1 global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Method 2 for conversion from lev (ln(p0/p)) to pressure in hPa 
!   (millibars), using p0_model:
!
      write(char80,"('p0_model*1.e-3*exp(-lev(k))')")
      istat = nf_put_att_text(ncid,NF_GLOBAL,"lev_to_hPa_method2",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for lev_to_hPa_method2 global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Number of histories on the file:
      ivar1(1) = h%ihist ! this will be updated in nc_wrhist
      istat = nf_put_att_int(ncid,NF_GLOBAL,"nhist",NF_INT,1,ivar1)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_int ',
     |    'for nhist global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Delta minutes between histories on the file:
      ivar1(1) = h%delhmins
      istat = nf_put_att_int(ncid,NF_GLOBAL,"delhist_mins",NF_INT,
     |  1,ivar1)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_int ',
     |    'for delhist_mins global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Missing_value (this can be used by procs if the field does not
! have its own missing_value attribute)
      var1(1) = spval
      istat = nf_put_att_double(ncid,NF_GLOBAL,"missing_value",
     |  NF_DOUBLE,1,var1)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from ',
     |    'nf_put_att_double for missing_value global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Potential model used (weimer, heelis, or none):
      istat = nf_put_att_text(ncid,NF_GLOBAL,"potential_model",
     |  len_trim(potential_model),trim(potential_model))
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for potential_model global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! tuv_lbc_intop==0 (lbc of t,u,v are not stored in top k-slot)
! (only "old" histories stored lbc of t,u,v in top slot)
      ivar1(1) = h%tuv_lbc_intop ! see init in hist.F
      istat = nf_put_att_int(ncid,NF_GLOBAL,"tuv_lbc_intop",
     |  NF_INT,1,ivar1)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_int ',
     |    'for tuv_lbc_intop global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Exit define mode:
      istat = nf_enddef(ncid)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     +  'Error return from nf_enddef')
!
! Give values to coordinate variables:
!
! Geographic horizontal:
      istat = nf_put_var_double(ncid,idv_lon,glon)
      if (istat /= NF_NOERR)
     |  call handle_ncerr(istat,'Error giving values to glon coord var')
      istat = nf_put_var_double(ncid,idv_lat,glat)
      if (istat /= NF_NOERR)
     |  call handle_ncerr(istat,'Error giving values to glat coord var')
!
! Vertical:
      istat = nf_put_var_double(ncid,idv_lev,zpmid)  ! midpoint levels
      if (istat /= NF_NOERR)
     | call handle_ncerr(istat,'Error giving values to zpmid coord var')
      istat = nf_put_var_double(ncid,idv_ilev,zpint) ! interface levels
      if (istat /= NF_NOERR)
     | call handle_ncerr(istat,'Error giving values to zpint coord var')
!
! Bottom interface level LBC for t,u,v lbc:
      istat = nf_put_var_double(ncid,idv_lbc,zibot)
      if (istat /= NF_NOERR)
     | call handle_ncerr(istat,'Error giving value to LBC scalar var')
!
! Magnetic:
      istat = nf_put_var_double(ncid,idv_mlon,gmlon)
      istat = nf_put_var_double(ncid,idv_mlat,gmlat)
!
      istat = nf_put_var_double(ncid,idv_mlev,zpmag) 
      istat = nf_put_var_double(ncid,idv_imlev,zpimag)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_var_double ',
     +    ' to assign values to magnetic dimension vars')")
        call handle_ncerr(istat,char120)
      endif

!     write(6,"('nc_define: mag midpoints zpmag=',/,(8f8.2))")
!    |  zpmag
!     write(6,"('nc_define: mag interfaces zpimag=',/,(8f8.2))") 
!    |  zpimag

      end subroutine nc_define
!-----------------------------------------------------------------------
      subroutine define_ncfile(ncid,name,idv,att1name,att1val,
     |  att2name,att2val)
!
! Args:
      integer,intent(in) :: ncid
      integer,intent(out) :: idv
      character(len=*),intent(in) :: name,att1name,att1val,
     |  att2name,att2val
!
! Local:
      integer :: istat,ids2(2)
      character(len=120) :: char120
!
      if (len_trim(name) <= 0) then
        write(6,"('>>> define_ncfile: empty name! returning..')")
        return
      endif
!
      ids2(1) = id_filelen
      ids2(2) = id_time
!
! Define time-dependent char variable:
      istat= nf_def_var(ncid,name,NF_CHAR,2,ids2,idv)
      if (istat /= NF_NOERR) then
        write(char120,"('Error defining file var ',a)") trim(name)
        call handle_ncerr(istat,char120)
      endif
!
! First attribute:
      if (len_trim(att1name) > 0.and.len_trim(att1val) > 0) then
        istat = nf_put_att_text(ncid,idv,att1name,len_trim(att1val),
     |    trim(att1val))
        if (istat /= NF_NOERR) then
          write(char120,"('Error defining ',a,' attribute to var ',
     |      a,': attribute value=',a)") trim(name),trim(att1name),
     |      trim(att1val)
          call handle_ncerr(istat,char120)
        endif
      endif
!
! Second attribute:
      if (len_trim(att2name) > 0.and.len_trim(att2val) > 0) then
        istat = nf_put_att_text(ncid,idv,att2name,len_trim(att2val),
     |    trim(att2val))
        if (istat /= NF_NOERR) then
          write(char120,"('Error defining ',a,' attribute to var ',
     |      a,': attribute value=',a)") trim(name),trim(att2name),
     |      trim(att2val)
          call handle_ncerr(istat,char120)
        endif
      endif
      end subroutine define_ncfile
!-----------------------------------------------------------------------
      subroutine defvar_time_dbl(ncid,name,idtime,idvar,
     |  longname,units)
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
      subroutine deflbc(ncid,vname,idv,type)
      use input_module,only: sech_nbyte
!
! Define t, u, or v lbc on current history file. Return new varialbe
!   id in idv. (vname should be TLBC, ULBC, or VLBC)
! If defining a secondary history, check sech_nbyte.
!
! Args:
      integer,intent(in) :: ncid
      integer,intent(out) :: idv
      character(len=*),intent(in) :: vname,type
!
! Local:
      integer :: istat,len
      character(len=80) :: char80,units
!
      ids3(1) = id_lon
      ids3(2) = id_lat
      ids3(3) = id_time
!
! Define variable:
! Always write doubles to primary files. For secondary files, write
!   4-byte floats if sech_nbyte=4, or 8-byte doubles if sech_nbyte=8
!
      if (trim(type)=='primary') then
        istat = nf_def_var(ncid,trim(vname),NF_DOUBLE,3,ids3,idv)
      else
        if (sech_nbyte==4) then
          istat = nf_def_var(ncid,trim(vname),NF_REAL,3,ids3,idv)
        else
          istat = nf_def_var(ncid,trim(vname),NF_DOUBLE,3,ids3,idv)
        endif
      endif
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining 2-d var ',a)") trim(vname)
        call handle_ncerr(istat,trim(char80))
      endif
!
! Long name attribute:
! (Variable names ending in '_NM' are at previous time n-1)
!
      write(char80,"('Lower boundary condition of ',a,'N')") vname(1:1)
      len = len_trim(vname)
      if (vname(len-2:len)=='_NM') 
     |  char80 = trim(char80)//' (TIME N-1)'
!
      istat = nf_put_att_text(ncid,idv,"long_name",
     |  len_trim(char80),trim(char80))
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining long_name attribute of ',
     |    'variable ',a)") trim(vname)
        call handle_ncerr(istat,trim(char80))
      endif
!
! Units attribute:
      units = 'cm/s'
      if (vname(1:1)=='T') units = 'K'
      istat = nf_put_att_text(ncid,idv,"units",len_trim(units),
     |  trim(units))
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining units attribute of variable ',
     |    a)") trim(vname)
        call handle_ncerr(istat,trim(char80))
      endif
!
! Coordinates attribute:
      istat = nf_put_att_text(ncid,idv,"coordinates",3,'LBC')
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining coordinates attribute of ',
     |    'variable',a)") trim(vname)
        call handle_ncerr(istat,trim(char80))
      endif
!     write(6,"('deflbc: defined var ',a,' idv=',i3)") trim(vname),idv
      end subroutine deflbc
!-----------------------------------------------------------------------
      subroutine def_fsech(ncid)
!
! Define secondary history fields (geo or mag, 2d+time or 3d+time):
!
      use fields_module,only: fsechist
      use hist_module,only: nsech,isecout
      use input_module,only: secout,sech_nbyte
      implicit none
!
! Args:
      integer,intent(in) :: ncid
!
! Local:
      integer :: ix,istat,iprog,idimids(4),idv,isdef
      character(len=80) :: char80
      real :: var1(1)
!
! External:
      integer,external :: strloc
!
      if (nsech==2) then
        istat = nf_redef(ncid) ! put dataset in define mode
        if (istat /= NF_NOERR)
     |    call handle_ncerr(istat,'error from nf_redef')
      endif
!
! Fields loop:
      floop: do ix=1,mxfsech
        if (len_trim(fsechist(ix)%short_name)==0) cycle floop
!
! Define prognostic only if this is 1st history on this file:
        iprog=strloc(f4d%short_name,size(f4d),fsechist(ix)%short_name)
        if (iprog > 0.and.nsech==1) then
!
! Define prognostic secondary history field:
          call defvar_f4d(ncid,f4d(iprog),idv_fsech(ix),sech_nbyte)
!
! 10/27/05 btf: Assume user-requested 'W' refers to "old hist" W, i.e. OMEGA:
!
        elseif (nsech==1.and.fsechist(ix)%short_name=='W') then
          write(6,"('Note: defining requested secondary history ',
     |      'variable ''W'' as prognostic ''OMEGA''')") 
          iprog = strloc(f4d%short_name,size(f4d),'OMEGA')
          call defvar_f4d(ncid,f4d(iprog),idv_fsech(ix),sech_nbyte)
          cycle floop
        endif ! iprog and nsech
!
! Define diagostics only if not first timestep (because in that case, 
! diags have not yet been allocated by addfld):
        if (istep==0.or.iprog > 0) cycle floop
!
! If data pointer is not associated, addfld was probably never called
! for this field. In this case, do not define the field. (subs def_fsech
! and wrfsech will also check this pointer, see nchist.F, and sub
! mp_gather2root_fsech in mpi.F)
!
        if (.not.associated(fsechist(ix)%data)) then
          write(6,"('>>> WARNING: Field ',a,' will not be defined',
     |      ' on secondary histories because',/,4x,'data was not ',
     |      'allocated. (addfld was probably not called for this ',
     |      'field)')") trim(fsechist(ix)%short_name)
          cycle floop
        endif
!
! Get dim id's for 3 dimensions (geo or mag):
        if (fsechist(ix)%dimsizes(3) /= 0) then ! 3d var (lon,lat,lev,time)
          istat=nf_inq_dimid(ncid,fsechist(ix)%dimnames(1),idimids(1)) 
          if (istat /= NF_NOERR) then
            write(char80,"('def_fsech: Error getting id for dim ',
     |        a,' (1st of 3)')") trim(fsechist(ix)%dimnames(1))
            call handle_ncerr(istat,trim(char80))
          endif 
          istat=nf_inq_dimid(ncid,fsechist(ix)%dimnames(2),idimids(2)) 
          if (istat /= NF_NOERR) then
            write(char80,"('def_fsech: Error getting id for dim ',
     |        a,' (2nd of 3)')") trim(fsechist(ix)%dimnames(2))
            call handle_ncerr(istat,trim(char80))
          endif 
          istat=nf_inq_dimid(ncid,fsechist(ix)%dimnames(3),idimids(3)) 
          if (istat /= NF_NOERR) then
            write(char80,"('def_fsech: Error getting id for dim ',
     |        a,' (3rd of 3)')") trim(fsechist(ix)%dimnames(3))
            call handle_ncerr(istat,trim(char80))
          endif 
          idimids(4) = id_time
!
! Define 3d var (lon,lat,lev,time)
          istat = nf_inq_varid(ncid,fsechist(ix)%short_name,idv)
          isdef = 0                  ! assume var is not yet defined
          if (istat == 0) isdef = 1  ! var is already defined -- cannot redef
          if (isdef == 0) then ! safe to define
            if (sech_nbyte == 4) then ! write 4-byte reals
              istat = nf_def_var(ncid,fsechist(ix)%short_name,NF_REAL,
     |          4,idimids,idv_fsech(ix)) ! returns idv_fsech(ix)
            else                      ! write 8-byte doubles
              istat = nf_def_var(ncid,fsechist(ix)%short_name,NF_DOUBLE,
     |          4,idimids,idv_fsech(ix)) ! returns idv_fsech(ix)
            endif
            if (istat /= NF_NOERR) then
              write(char80,"('def_fsech: Error defining 3d field ',a)")
     |          trim(fsechist(ix)%short_name)
              call handle_ncerr(istat,trim(char80))
            endif
          endif  ! if var not yet defined
!
! Get dim id's for 2 dimensions (geo or mag):
        else                              ! 2d var (lon,lat,time)
          istat=nf_inq_dimid(ncid,fsechist(ix)%dimnames(1),idimids(1)) 
          if (istat /= NF_NOERR) then
            write(char80,"('def_fsech: Error getting id for dim ',
     |        a,' (1st of 2)')") trim(fsechist(ix)%dimnames(1))
            call handle_ncerr(istat,trim(char80))
          endif
          istat=nf_inq_dimid(ncid,fsechist(ix)%dimnames(2),idimids(2)) 
          if (istat /= NF_NOERR) then
            write(char80,"('def_fsech: Error getting id for dim ',
     |        a,' (2nd of 2)')") trim(fsechist(ix)%dimnames(2))
            call handle_ncerr(istat,trim(char80))
          endif 
          idimids(3) = id_time
!
! Define 2d var (lon,lat,time)
          istat = nf_inq_varid(ncid,fsechist(ix)%short_name,idv)
          isdef = 0                  ! assume var is not yet defined
          if (istat == 0) isdef = 1  ! var is already defined -- cannot redef
          if (isdef == 0) then ! safe to define
            if (sech_nbyte == 4) then ! write 4-byte reals
              istat = nf_def_var(ncid,fsechist(ix)%short_name,NF_REAL,
     |          3,idimids(1:3),idv_fsech(ix)) ! returns idv_fsech(ix)
            else
              istat = nf_def_var(ncid,fsechist(ix)%short_name,NF_DOUBLE,
     |          3,idimids(1:3),idv_fsech(ix)) ! returns idv_fsech(ix)
            endif
            if (istat /= NF_NOERR) then
              write(char80,"('def_fsech: Error defining 2d field ',a)")
     |          fsechist(ix)%short_name
              call handle_ncerr(istat,trim(char80))
            endif
          endif ! isdef==0
        endif  ! 2d or 3d

!       write(6,"('def_fsech: defined sech field ',a,': dimnames=',
!    |    a,',',a,',',a,' dimsizes=',3i4)") fsechist(ix)%short_name,
!    |    trim(fsechist(ix)%dimnames(1)),trim(fsechist(ix)%dimnames(2)),
!    |    trim(fsechist(ix)%dimnames(3)),fsechist(ix)%dimsizes
!
! Attributes:
        if (len_trim(fsechist(ix)%long_name) > 0) then
          istat = nf_put_att_text(ncid,idv_fsech(ix),"long_name",
     |    len_trim(fsechist(ix)%long_name),trim(fsechist(ix)%long_name))
        endif
        if (len_trim(fsechist(ix)%units) > 0) then
          istat = nf_put_att_text(ncid,idv_fsech(ix),"units",
     |    len_trim(fsechist(ix)%units),trim(fsechist(ix)%units))
        endif
!
! Missing value:
        var1(1) = spval
        istat = nf_put_att_double(ncid,idv_fsech(ix),"missing_value",
     |    NF_DOUBLE,1,var1)
        if (istat /= NF_NOERR) then
          write(char80,"('def_fsech: Error return from ',
     |      'nf_put_att_double for missing_value var attribute')")
          call handle_ncerr(istat,char80)
        endif
!
! End field loop:
      enddo floop ! i=1,mxfsech
      if (nsech==2) istat = nf_enddef(ncid) ! take out of define mode
      end subroutine def_fsech
!-----------------------------------------------------------------------
      subroutine defvar_f4d(ncid,f,idvar,nbytes)
!
! Define a prognostic variable on current netcdf history file:
! (not data -- see wrf4d for data write)
!
! Args:
      integer,intent(in)  :: ncid,nbytes
      integer,intent(out) :: idvar
      type(fields_4d),intent(in) :: f
!
! Local:
      integer :: istat,idimids(4)
      character(len=80) :: char80
      character(len=160) :: char160
      real :: var1(1)
!
! Prognostics are dimensioned (lon,lat,lev,time):
      if (.not.fakeflds) then
        idimids(1) = id_lon
        idimids(2) = id_lat
        idimids(3) = id_lev                             ! midpoints
        if (f%vcoord(1:3)=='int') idimids(3) = id_ilev  ! interfaces
        idimids(4) = id_time
      else
        idimids(1:3) = id_fake
        idimids(4) = id_time
      endif
!
! This call returns the variable id, idvar.
!   nbytes == 4 means write 4-byte reals (to secondary files), otherwise
!   nbytes == 8 means write 8-byte doubles (either secondary or primary files).
!
      if (nbytes == 4) then     ! write 4-byte reals
        istat = nf_def_var(ncid,f%short_name,NF_REAL,4,idimids,idvar)
      else                      ! write 8-byte doubles
        istat = nf_def_var(ncid,f%short_name,NF_DOUBLE,4,idimids,idvar)
      endif
      if (istat /= NF_NOERR) then
        write(char80,"('defvar_f4d: Error defining 4-d field ',a,
     |    ' nbytes=',i2)") f%short_name,nbytes
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
! 1/18/06 btf: Use f%units for densities) instead of "fraction" or "1",
!              and do not add alternate units.
!
      istat = nf_put_att_text(ncid,idvar,"units",len_trim(f%units),
     |  f%units)
      if (istat /= NF_NOERR) then
        write(char80,"('Error defining units of 4-d',
     |    ' variable ',a,': units = ',a)") trim(f%short_name),
     |    trim(f%units)
        call handle_ncerr(istat,trim(char80))
      endif
!
! 10/27/05 btf: Add former_name attribute to OMEGA to indicate it
!               was previously known as W:
      if (f%short_name=='OMEGA') then
        istat = nf_put_att_text(ncid,idvar,"former_name",1,"W")
      endif
!
! Add missing_value attribute:
      var1(1) = spval
      istat = nf_put_att_double(ncid,idvar,"missing_value",
     |  NF_DOUBLE,1,var1)
      if (istat /= NF_NOERR) then
        write(char160,"('defvar_f4d: Error return from ',
     |    'nf_put_att_double for missing_value var attribute')")
        call handle_ncerr(istat,char160)
      endif
      end subroutine defvar_f4d
!-----------------------------------------------------------------------
      subroutine nc_wrhist(ncid,diskfile,iprint)
!
! Write current history structure to open netcdf output file:
!
      use hist_module,only: sh,h,hist_print
      use input_module,only: start,secstart,start_year
      use cons_module,only: crit,rtd
!
! Args:
      integer,intent(in) :: ncid,iprint
      character(len=*),intent(in) :: diskfile
!
! Local:
      integer :: mins,istat,ivar1(1),i,j,iprog,imo,ida,mtimeinit(4),
     |  ivar2(2),ii,mfinal
      character(len=120) :: char120
      character(len=240) :: char240
      character(len=80) :: char80
      integer :: ndims,nvars,ngatts,unlimdimid
      real :: var1(1),var22(2,2),var2(2),var10(10),rmins,rdays,
     |  crit1,crit2
!
! External:
      integer,external :: strloc
      real,external :: mtime_delta
!
! Time: decimal days since start time of initial run:
! h%initial_mtime was set by sub hist_initype (if initial run),
! or was read from source history (if continuation run). 

!>>> nc_wrhist: time from initial start must be >= 0: mtimeinit= 365   0   0  
!    h%modeltime=   1   1   0 delta (mins)=-524100.00
! 
      mtimeinit(1:3) = sh%initial_mtime(1:3) ; mtimeinit(4) = 0
      rmins = mtime_delta(mtimeinit,h%modeltime)

!     write(6,"('nc_wrhist calc rmins: h%year=',i5,' h%modeltime=',
!    |  4i4,' sh%year=',i5,' mtimeinit=',4i4,' rmins=',f12.2)") 
!    |  h%year,h%modeltime,sh%year,mtimeinit,rmins 

      if (mtimeinit(1) > h%modeltime(1)) then
        mtimeinit(:) = h%modeltime(:)
        rmins = mtime_delta(mtimeinit,h%modeltime)
        write(6,"('nc_wrhist: reset mtimeinit = h%modeltime = ',4i5)")
     |    mtimeinit
      endif
      rdays = rmins/1440.
      
      if (rmins < 0.) then
        write(6,"(/,'>>> nc_wrhist: time from initial start must be',
     |    ' >= 0: mtimeinit=',3i4,'  h%modeltime=',3i4,' delta (mins)=',
     |    f10.2)") mtimeinit(1:3),h%modeltime(1:3),rmins
        call shutdown('delta time from init run')
      else
!       write(6,"('nc_wrhist: mtimeinit=',3i4,' h%modeltime=',3i4,
!    |    ' delta (mins)=',f10.2,' delta (days)=',f8.2)")
!    |    mtimeinit(1:3),h%modeltime(1:3),rmins,rdays
      endif
      if (idv_time <= 0) istat = nf_inq_varid(ncid,"time",idv_time)
! 1/29/08 btf: Use minutes rather than days 
!              (see time units "minutes since.." above)
!     istat = nf_put_var1_double(ncid,idv_time,h%ihist,rdays)
      istat = nf_put_var1_double(ncid,idv_time,h%ihist,rmins)
      if (istat /= NF_NOERR) then
        write(char240,"('Error from nf_put_var1_double defining ',
     |    'time at h%modeltime=',4i4,' mtimeinit=',3i4,
     |    ' delta (days)=',f8.2)") h%modeltime,mtimeinit(1:3),rdays
        call handle_ncerr(istat,char240)
      endif
!
! Model time (day,hour,minute):
      start_2d(1) = 1
      start_2d(2) = h%ihist
      count_2d(1) = 3
      count_2d(2) = 1
      if (idv_mtime <= 0) istat = nf_inq_varid(ncid,"mtime",idv_mtime)
!
! If model day is > 365, write actual day of year, i.e. day-365
! (e.g., if model day is 366, write model day 1 to the history)
! Check for leap-year.
!
      if((mod(start_year,4).eq.0.and.mod(start_year,100).ne.0).or.
     |   (mod(start_year,400).eq.0)) then
         mfinal=366
      else
         mfinal=365
      endif
      if (h%modeltime(1) > mfinal) then
        write(6,"('Note: Changed model day from ',i4,' to ',i4,
     |    ' before writing to history file.')")
     |    h%modeltime(1), h%modeltime(1)-mfinal
        h%modeltime(1) = h%modeltime(1)-mfinal
      endif
!
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
! Calendar advance:
      if (idv_calendar_adv <= 0) 
     |  istat = nf_inq_varid(ncid,"calendar_advance",idv_calendar_adv)
      istat = nf_put_var1_int(ncid,idv_calendar_adv,h%ihist,
     |  h%calendar_advance)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     |    'calendar advance flag at h%ihist=',i3,
     |    ': h%calendar_advance=',i5)") h%ihist,h%calendar_advance
        call handle_ncerr(istat,char120)
      endif
!
! Date and time history was written:
      if (idv_writedate <= 0) istat = nf_inq_varid(ncid,"write_date",
     |  idv_writedate)
      start_2d(1) = 1
      start_2d(2) = h%ihist
      count_2d(1) = len(h%write_date) 
      count_2d(2) = 1
      istat = nf_put_vara_text(ncid,idv_writedate,start_2d,count_2d,
     |  h%write_date)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_vara_text writing ',
     |    'h%write_date=',a)") trim(h%write_date)
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
! Number of mpi tasks:
      if (idv_ntask <= 0) istat = nf_inq_varid(ncid,"ntask_mpi",
     |  idv_ntask)
      istat = nf_put_var1_int(ncid,idv_ntask,h%ihist,h%ntask_mpi)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     |    'ntask_mpi at h%ihist=',i3,': h%ntask_mpi=',i8)") h%ihist,
     |    h%ntask_mpi
        call handle_ncerr(istat,char120)
      endif
!
! If coupled with CISM/CMIT:
      if (idv_coupled_cmit <= 0) istat = nf_inq_varid(ncid,
     |  "coupled_cmit",idv_coupled_cmit)
      istat = nf_put_var1_int(ncid,idv_coupled_cmit,h%ihist,
     |  h%coupled_cmit)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_var1_int defining ',
     |    'coupled_cmit at h%ihist=',i3,': h%coupled_cmit=',i8)") 
     |    h%ihist,h%coupled_cmit
        call handle_ncerr(istat,char120)
      endif
!
! Universal time (decimal hours):
      if (idv_ut <= 0) istat = nf_inq_varid(ncid,"ut",idv_ut)
      istat = nf_put_var1_double(ncid,idv_ut,h%ihist,h%ut)
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
! Update gpi_ncfile:
      call update_ncfile(ncid,"gpi_ncfile",h%gpi_ncfile,idv_gpi_ncfile,
     |  h%ihist)
!
! Update ncep_ncfile:
      call update_ncfile(ncid,"ncep_ncfile",h%ncep_ncfile,
     |  idv_ncep_ncfile,h%ihist)
!
! Update see_ncfile:
      call update_ncfile(ncid,"see_ncfile",h%see_ncfile,idv_see_ncfile,
     |  h%ihist)
!
! Update imf_ncfile:
      call update_ncfile(ncid,"imf_ncfile",h%imf_ncfile,idv_imf_ncfile,
     |  h%ihist)
!
! Update gswm files:
      call update_ncfile(ncid,"gswm_mi_di_ncfile",h%gswm_mi_di_ncfile,
     |  idv_gswm_mi_di_ncfile,h%ihist)
      call update_ncfile(ncid,"gswm_mi_sdi_ncfile",h%gswm_mi_sdi_ncfile,
     |  idv_gswm_mi_sdi_ncfile,h%ihist)
      call update_ncfile(ncid,"gswm_nm_di_ncfile",h%gswm_nm_di_ncfile,
     |  idv_gswm_nm_di_ncfile,h%ihist)
      call update_ncfile(ncid,"gswm_nm_sdi_ncfile",h%gswm_nm_sdi_ncfile,
     |  idv_gswm_nm_sdi_ncfile,h%ihist)
!
! Update saber_ncfile:
      call update_ncfile(ncid,"saber_ncfile",h%saber_ncfile,
     |  idv_saber_ncfile,h%ihist)
!
! Update tidi_ncfile:
      call update_ncfile(ncid,"tidi_ncfile",h%tidi_ncfile,
     |  idv_tidi_ncfile,h%ihist)
!
! Update bgrddata_ncfile:
      call update_ncfile(ncid,"bgrddata_ncfile",h%bgrddata_ncfile,
     |  idv_bgrddata_ncfile,h%ihist)
!
! Update ctmt_ncfile:
      call update_ncfile(ncid,"ctmt_ncfile",h%ctmt_ncfile,
     |  idv_ctmt_ncfile,h%ihist)
!
! f107d and f107a (daily and average solar fluxes):
      if (idv_f107d <= 0) istat = nf_inq_varid(ncid,"f107d",
     |  idv_f107d)
      istat = nf_put_var1_double(ncid,idv_f107d,h%ihist,h%f107d)
      if (idv_f107a <= 0) istat = nf_inq_varid(ncid,"f107a",
     |  idv_f107a)
      istat = nf_put_var1_double(ncid,idv_f107a,h%ihist,h%f107a)
!
! Hemispheric power:
      if (idv_hpower <= 0) istat = nf_inq_varid(ncid,"hpower",
     |  idv_hpower)
      istat = nf_put_var1_double(ncid,idv_hpower,h%ihist,h%hpower)
!
! Cross-tail potential:
      if (idv_ctpoten <= 0) istat = nf_inq_varid(ncid,"ctpoten",
     |  idv_ctpoten)
      istat = nf_put_var1_double(ncid,idv_ctpoten,h%ihist,h%ctpoten)
!
! Kp index:
      if (idv_kp <= 0) istat = nf_inq_varid(ncid,"Kp",idv_kp)
      istat = nf_put_var1_double(ncid,idv_kp,h%ihist,h%kp)
!
! BX,BY,BZ imf:
      if (idv_bximf <= 0) istat = nf_inq_varid(ncid,"bximf",idv_bximf)
      istat = nf_put_var1_double(ncid,idv_bximf,h%ihist,h%bximf)
      if (idv_byimf <= 0) istat = nf_inq_varid(ncid,"byimf",idv_byimf)
      istat = nf_put_var1_double(ncid,idv_byimf,h%ihist,h%byimf)
      if (idv_bzimf <= 0) istat = nf_inq_varid(ncid,"bzimf",idv_bzimf)
      istat = nf_put_var1_double(ncid,idv_bzimf,h%ihist,h%bzimf)
!
! swvel,swden:
      if (idv_swvel <= 0) istat = nf_inq_varid(ncid,"swvel",idv_swvel)
      istat = nf_put_var1_double(ncid,idv_swvel,h%ihist,h%swvel)
      if (idv_swden <= 0) istat = nf_inq_varid(ncid,"swden",idv_swden)
      istat = nf_put_var1_double(ncid,idv_swden,h%ihist,h%swden)
!
! al:
      if (idv_al <= 0) istat = nf_inq_varid(ncid,"swden",idv_al)
      istat = nf_put_var1_double(ncid,idv_al,h%ihist,h%al)
!
! Collision factor:
      if (idv_colfac <= 0) istat = nf_inq_varid(ncid,"colfac",
     |  idv_colfac)
      istat = nf_put_var1_double(ncid,idv_colfac,h%ihist,h%colfac)
!
! Joule heating factor:
      if (idv_joulefac <= 0) istat = nf_inq_varid(ncid,"joulefac",
     |  idv_joulefac)
      istat = nf_put_var1_double(ncid,idv_joulefac,h%ihist,h%joulefac)
 
! alfa30:
!     if (idv_alfa30 <= 0) istat = nf_inq_varid(ncid,"alfa30",
!    |  idv_alfa30)
!     istat = nf_put_var1_double(ncid,idv_alfa30,h%ihist,h%alfa30)
!
! e30:
!     if (idv_e30 <= 0) istat = nf_inq_varid(ncid,"e30",idv_e30)
!     istat = nf_put_var1_double(ncid,idv_e30,h%ihist,h%e30)
!
! alfad2:
!     if (idv_alfad2 <= 0) istat = nf_inq_varid(ncid,"alfad2",
!    |  idv_alfad2)
!     istat = nf_put_var1_double(ncid,idv_alfad2,h%ihist,h%alfad2)
!
! ed2:
!     if (idv_ed2 <= 0) istat = nf_inq_varid(ncid,"ed2",idv_ed2)
!     istat = nf_put_var1_double(ncid,idv_ed2,h%ihist,h%ed2)

! e1,e2:
      if (idv_e1 <= 0) istat = nf_inq_varid(ncid,"e1",idv_e1)
      istat = nf_put_var1_double(ncid,idv_e1,h%ihist,h%e1)
      if (idv_e2 <= 0) istat = nf_inq_varid(ncid,"e2",idv_e2)
      istat = nf_put_var1_double(ncid,idv_e2,h%ihist,h%e2)
! h1,h2:
      if (idv_h1 <= 0) istat = nf_inq_varid(ncid,"h1",idv_h1)
      istat = nf_put_var1_double(ncid,idv_h1,h%ihist,h%h1)
      if (idv_h2 <= 0) istat = nf_inq_varid(ncid,"h2",idv_h2)
      istat = nf_put_var1_double(ncid,idv_h2,h%ihist,h%h2)
! alfac,ec:
      if (idv_alfac <= 0) istat = nf_inq_varid(ncid,"alfac",idv_alfac)
      istat = nf_put_var1_double(ncid,idv_alfac,h%ihist,h%alfac)
      if (idv_ec <= 0) istat = nf_inq_varid(ncid,"ec",idv_ec)
      istat = nf_put_var1_double(ncid,idv_ec,h%ihist,h%ec)
! alfad,ed:
      if (idv_alfad <= 0) istat = nf_inq_varid(ncid,"alfad",idv_alfad)
      istat = nf_put_var1_double(ncid,idv_alfad,h%ihist,h%alfad)
      if (idv_ed <= 0) istat = nf_inq_varid(ncid,"ed",idv_ed)
      istat = nf_put_var1_double(ncid,idv_ed,h%ihist,h%ed)
! crit1,2:
      crit1 = crit(1)*rtd
      if (idv_crit1 <= 0) istat = nf_inq_varid(ncid,"crit1",idv_crit1)
      istat = nf_put_var1_double(ncid,idv_crit1,h%ihist,crit1)
      crit2 = crit(2)*rtd
      if (idv_crit2 <= 0) istat = nf_inq_varid(ncid,"crit2",idv_crit2)
      istat = nf_put_var1_double(ncid,idv_crit2,h%ihist,crit2)
!
! Reference pressure used to convert ln(p0/p) to pressure (hPa):
! (h%p0 was set in sub define_hist, output.F)
      if (idv_p0 <= 0) istat = nf_inq_varid(ncid,"p0",idv_p0)
      istat = nf_put_var1_double(ncid,idv_p0,h%ihist,h%p0)
!
! Reference pressure used by the model (microbars):
! (h%p0_model was set in sub define_hist, output.F)
      if (idv_p0_model <= 0) istat = nf_inq_varid(ncid,"p0_model",
     |  idv_p0_model)
      istat = nf_put_var1_double(ncid,idv_p0_model,h%ihist,h%p0_model)
!
! Gravitational constant used in the model:
      if (idv_grav <= 0) istat = nf_inq_varid(ncid,"grav",idv_grav)
      istat = nf_put_var1_double(ncid,idv_grav,h%ihist,h%grav)
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
! Update run_type (initial or continuation):
      istat = nf_redef(ncid) ! put dataset in define mode

      istat = nf_put_att_text(ncid,NF_GLOBAL,"run_type",
     |  len_trim(h%run_type),trim(h%run_type))
      if (istat /= NF_NOERR) then
        write(char120,"('nc_wrhist: Error return from nf_put_att_text ',
     |    'for run_type global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Update source_file and source mtime:
      if (h%run_type(1:4)=='init') then ! initial source file
        istat = nf_put_att_text(ncid,NF_GLOBAL,"source_file",
     |    len_trim(h%source_file)+10,trim(h%source_file)//' (initial)')
      elseif (h%run_type(1:4)=='cont') then ! continuation source file
        istat = nf_put_att_text(ncid,NF_GLOBAL,"source_file",
     |    len_trim(h%source_file)+15,
     |    trim(h%source_file)//' (continuation)')
      endif
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for source_file global attribute')")
        call handle_ncerr(istat,char120)
      endif
      istat = nf_put_att_int(ncid,NF_GLOBAL,"source_mtime",NF_INT,3,
     |  h%source_mtime)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_text ',
     |    'for source_mtime global attribute')")
        call handle_ncerr(istat,char120)
      endif

      istat = nf_enddef(ncid) ! take out of define mode
!
! Update number of histories on the file (global attribute):
      ivar1(1) = h%ihist
      istat = nf_put_att_int(ncid,NF_GLOBAL,"nhist",NF_INT,1,ivar1)
      if (istat /= NF_NOERR) then
        write(char120,"('Error return from nf_put_att_int ',
     +    'for nhist global attribute')")
        call handle_ncerr(istat,char120)
      endif
!
! Add fields data to primary history:
! Currently (2/00), only prognostic fields are written to the
!   primary histories:
!
      idv_tlbc = 0 ; idv_ulbc = 0 ; idv_vlbc = 0          ! forces wrlbc to inquire
      idv_tlbc_nm = 0 ; idv_ulbc_nm = 0 ; idv_vlbc_nm = 0 ! forces wrlbc to inquire
      if (h%hist_type(1:3)=='pri') then
        do i=1,nf4d_hist
          call wrf4d(ncid,f4d(i)%short_name,h%ihist,fakeflds,i)
!
! Write ZG after Z (zg was calculated by sub calczg in addiag.F,
! then gathered to fzg by mp_gather2root_f3d):
! 3/1/06 btf: do not write ZG to primary histories.
!         if (trim(f4d(i)%short_name)=='Z')
!    |      call wrf3d(ncid,fzg,h%ihist,'ZG',idv_zg)
        enddo ! i=1,nf4d_hist
!
! Update t,u,v lbc to primary history: 
! Global arrays t,u,vlbc_glb were gathered from subdomains for history 
!   output by sub mp_gather2root_lbc (mpi.F).
!
        call wrlbc(ncid,"TLBC",tlbc_glb,h%ihist,idv_tlbc,'pri')
        call wrlbc(ncid,"ULBC",ulbc_glb,h%ihist,idv_ulbc,'pri')
        call wrlbc(ncid,"VLBC",vlbc_glb,h%ihist,idv_vlbc,'pri')
!
        call wrlbc(ncid,"TLBC_NM",tlbc_nm_glb,h%ihist,idv_tlbc_nm,'pri')
        call wrlbc(ncid,"ULBC_NM",ulbc_nm_glb,h%ihist,idv_ulbc_nm,'pri')
        call wrlbc(ncid,"VLBC_NM",vlbc_nm_glb,h%ihist,idv_vlbc_nm,'pri')
!
! Add fields data to secondary history:
      elseif (h%hist_type(1:3)=='sec') then
        call wrfsech(ncid,h%ihist) ! write secondary prog and diag fields
!
! Update t,u,v lbc to secondary history:
        call wrlbc(ncid,"TLBC",tlbc_glb,h%ihist,idv_tlbc,'sec')
        call wrlbc(ncid,"ULBC",ulbc_glb,h%ihist,idv_ulbc,'sec')
        call wrlbc(ncid,"VLBC",vlbc_glb,h%ihist,idv_vlbc,'sec')
        call wrlbc(ncid,"TLBC_NM",tlbc_nm_glb,h%ihist,idv_tlbc_nm,'sec')
        call wrlbc(ncid,"ULBC_NM",ulbc_nm_glb,h%ihist,idv_ulbc_nm,'sec')
        call wrlbc(ncid,"VLBC_NM",vlbc_nm_glb,h%ihist,idv_vlbc_nm,'sec')
      else
        write(6,"(/,'>>> WARNING nc_wrhist: unknown h%hist_type=',a)") 
     |    h%hist_type
      endif
!
! Report to stdout:
!     write(6,"('nc_wrhist: iprint=',i3)") iprint
      if (iprint > 0) call hist_print(h,'WRITE',diskfile)
      end subroutine nc_wrhist
!-----------------------------------------------------------------------
      subroutine wrfsech(ncid,itime)
!
! Write data values of secondary history fields to ncfile.
! (fields were defined on the file by def_fsech)
! Data was gathered to root task by mp_gather2root_sechist (nchist.F)
!
      use fields_module,only: fsechist,shortname_len
      use hist_module,only: nsech,h
      implicit none
!
! Args:
      integer,intent(in) :: ncid,itime
!
! Local:
      integer :: ix,iprog,istat,idvar,iddims(4),itype,ndims,natts,i
      character(len=80) :: char80
      real,allocatable :: f2d(:,:),f3d(:,:,:)
!
! External:
      integer,external :: strloc
!
!     write(6,"('Enter wrfsech: nsech=',i3,' itime=',i3)") nsech,itime
      floop: do ix=1,mxfsech
        if (len_trim(fsechist(ix)%short_name)==0) cycle floop
!
! Write prognostic:
        iprog=strloc(f4d%short_name,size(f4d),fsechist(ix)%short_name)
        if (iprog > 0) then
          call wrf4d(ncid,fsechist(ix)%short_name,h%ihist,fakeflds,
     |      iprog)
          cycle floop
        elseif (fsechist(ix)%short_name=='W') then
          write(6,"('Note: defining requested secondary history ',
     |      'variable ''W'' as prognostic ''OMEGA''')") 
          iprog = strloc(f4d%short_name,size(f4d),'OMEGA')
          call wrf4d(ncid,'OMEGA',h%ihist,fakeflds,iprog)
          cycle floop
        endif
!
! Wrote diagnostic only if not first timestep:
        if (istep==0) cycle floop
!
! If addfld was never called for this field, data will not be
! allocated -- in this case, leave the field defined on the
! history, but do not write data (netcdf will default to _Fillvalue)
!
        if (.not.associated(fsechist(ix)%data)) then
!         write(6,"('>>> WARNING: sechist field ',a,': data not',
!    |      ' allocated.')") trim(fsechist(ix)%short_name)
          cycle floop
        endif
!
! Get field id:
        istat = nf_inq_varid(ncid,fsechist(ix)%short_name,idvar)
        if (istat /= NF_NOERR) then
          write(char80,"('wrfsech: Error getting var id of field ',a)")
     |      trim(fsechist(ix)%short_name)
          call handle_ncerr(istat,char80)
        endif
!
! 3d+time field:
        if (fsechist(ix)%dimsizes(3) > 0) then
!
! Allocate and transfer to local 3d output array:
          allocate(f3d(fsechist(ix)%dimsizes(1),
     |                 fsechist(ix)%dimsizes(2),
     |                 fsechist(ix)%dimsizes(3)),stat=istat)
          if (istat /= 0) then
            write(6,"(/,'>>> wrfsech: Error allocating f3d: ',
     |        'dimsizes=',3i4)") fsechist(ix)%dimsizes
            call shutdown('wrfsech')
          endif
          f3d = fsechist(ix)%data
!
! Write values for 3d+time var (lon,lat,lev,time)
          start_4d(1:3) = 1
          start_4d(4) = itime
          count_4d(2:3) = fsechist(ix)%dimsizes(2:3)
          count_4d(4) = 1
!
! Write 3d+time geo field:
! (write only 3:nlonp2 to exclude periodic points):
          if (fsechist(ix)%geo) then    ! 3d geo
            count_4d(1) = nlon

!         write(6,"('wrfsech: writing 3d sech field ',a,'(',a,'=',i3,
!    |      ',',a,'=',i3,',',a,'=',i3,'): min,max=',2e12.4)") 
!    |      fsechist(ix)%short_name(1:10),
!    |      trim(fsechist(ix)%dimnames(1)),fsechist(ix)%dimsizes(1),
!    |      trim(fsechist(ix)%dimnames(2)),fsechist(ix)%dimsizes(2),
!    |      trim(fsechist(ix)%dimnames(3)),fsechist(ix)%dimsizes(3),
!    |      minval(f3d(3:nlonp2,:,:)),maxval(f3d(3:nlonp2,:,:))

            istat = nf_put_vara_double(ncid,idvar,start_4d,count_4d,
     |        f3d(3:nlonp2,:,:))
!
! Write 3d+time mag field:
          else                          ! 3d mag
            count_4d(1) = nmlonp1

!         write(6,"('wrfsech: writing 3d sech mag field ',a,'(',a,'=',
!    |      i3,' ',a,'=',i3,',',a,'=',i3,'): min,max=',2e12.4)") 
!    |      trim(fsechist(ix)%short_name),
!    |      trim(fsechist(ix)%dimnames(1)),fsechist(ix)%dimsizes(1),
!    |      trim(fsechist(ix)%dimnames(2)),fsechist(ix)%dimsizes(2),
!    |      trim(fsechist(ix)%dimnames(3)),fsechist(ix)%dimsizes(3),
!    |      minval(f3d),maxval(f3d)

            istat = nf_put_vara_double(ncid,idvar,start_4d,count_4d,
     |        f3d)
          endif                         ! geo or mag
          if (istat /= NF_NOERR) then
            write(char80,"('wrfsech: Error writing 3d+time field ',a)")
     |        trim(fsechist(ix)%short_name)
            call handle_ncerr(istat,char80)
          endif
!
! Write 2d+time var (lon,lat,time):
        else
!
! Allocate and transfer to local 2d array:
          allocate(f2d(fsechist(ix)%dimsizes(1),
     |                 fsechist(ix)%dimsizes(2)),stat=istat)
          if (istat /= 0) then
            write(6,"(/,'>>> wrfsech: Error allocating f2d: ',
     |        'dimsizes=',3i4)") fsechist(ix)%dimsizes(1:2)
            call shutdown('wrfsech')
          endif
          f2d = fsechist(ix)%data(:,:,1)
!
! Write values for 2d+time var (lon,lat,time)
          start_3d(1:2) = 1
          start_3d(3) = itime
          count_3d(2) = fsechist(ix)%dimsizes(2)
          count_3d(3) = 1
!
! Write 2d+time geo field:
! (write only 3:nlonp2 to exclude periodic points):
          if (fsechist(ix)%geo) then      ! 2d geo
            count_3d(1) = nlon

!         write(6,"('wrfsech: writing 2d geo sech field ',a,'(',a,'=',
!    |      i3,',',a,'=',i3,'): min,max=',2e12.4)") 
!    |      fsechist(ix)%short_name(1:10),
!    |      trim(fsechist(ix)%dimnames(1)),fsechist(ix)%dimsizes(1),
!    |      trim(fsechist(ix)%dimnames(2)),fsechist(ix)%dimsizes(2),
!    |      minval(f2d(3:nlonp2,:)),maxval(f2d(3:nlonp2,:))

            istat = nf_put_vara_double(ncid,idvar,start_3d,count_3d,
     |        f2d(3:nlonp2,:))
!
! Write 2d+time mag field:
          else                           ! 2d mag
            count_3d(1) = nmlonp1

!         write(6,"('wrfsech: writing 2d mag sech field ',a,'(',a,'=',
!    |      i3,',',a,'=',i3,'): min,max=',2e12.4)") 
!    |      trim(fsechist(ix)%short_name),
!    |      trim(fsechist(ix)%dimnames(1)),fsechist(ix)%dimsizes(1),
!    |      trim(fsechist(ix)%dimnames(2)),fsechist(ix)%dimsizes(2),
!    |      minval(f2d),maxval(f2d)

            istat = nf_put_vara_double(ncid,idvar,start_3d,count_3d,
     |        f2d)
          endif
          if (istat /= NF_NOERR) then
            write(char80,"('wrfsech: Error writing 2d+time field ',a)")
     |        trim(fsechist(ix)%short_name)
            call handle_ncerr(istat,char80)
          endif
        endif
        if (allocated(f2d)) deallocate(f2d)
        if (allocated(f3d)) deallocate(f3d)
      enddo floop ! ix=1,mxfsech
      end subroutine wrfsech
!-----------------------------------------------------------------------
      subroutine update_ncfile(ncid,vname,ncfile,idv,ihist)
!
! Update file text var vname for history ihist. 
! Value of the variable to write is ncfile.
!
      implicit none
!
! Args:
      integer,intent(in) :: ncid,ihist
      integer,intent(inout) :: idv
      character(len=*),intent(in) :: vname,ncfile
!
! Local:
      integer :: istat,start(2),count(2)
      character(len=len(ncfile)) :: file
      character(len=120) :: char120
!
      if (idv <= 0) 
     |  istat = nf_inq_varid(ncid,vname,idv)
      start(1) = 1
      start(2) = ihist
      count(2) = 1
      file = ncfile
      if (len_trim(ncfile) <= 0) file = '[none]'
      count(1) = len(ncfile)
      istat = nf_put_vara_text(ncid,idv,start,count,file)
      if (istat /= NF_NOERR) then
        write(char120,"('Error from nf_put_vara_text to update ncfile ',
     |    a)") trim(file)
        call handle_ncerr(istat,char120)
      else
!       write(6,"('Updated var ',a,' for ihist=',i3,' file=',a)")
!    |    trim(vname),ihist,trim(file) 
      endif

      end subroutine update_ncfile
!-----------------------------------------------------------------------
      subroutine wrlbc(ncid,vname,flbc,itime,idv,type)
      implicit none
!
! Update global lbc of t,u,v (TLBC,ULBC,VLBC (lon,lat,time)) to 
!   current history (may be primary or secondary history).
! This is executed by master task only (subdomains were gathered
!   into flbc by mp_gather2root_lbc)
!
! Args:
      integer,intent(in) :: ncid,itime
      integer,intent(inout) :: idv
      character(len=*),intent(in) :: vname,type
      real,intent(in) :: flbc(nlonp4,nlat) ! data to save on file
!
! Local:
      integer :: istat,j
      character(len=80) :: char80
      character(len=120) :: char120
      real :: fglb(nlon,nlat) ! local for write to history
!
      start_3d(1:2) = 1
      start_3d(3) = itime
      count_3d(1) = nlon
      count_3d(2) = nlat
      count_3d(3) = 1
!
! idv may be 0 for primary continuation file: 
      if (idv <= 0) then
        istat = nf_inq_varid(ncid,trim(vname),idv)
        if (istat /= NF_NOERR) then
          write(char80,"('Error getting var id for ',a,
     |      ' (ncid=',i3,')')") trim(vname),ncid
          call handle_ncerr(istat,char80)
        else
!         write(6,"('wrlbc: vname=',a,' idv=',i3,' (idv was 0 on ',
!    |      'input)')") trim(vname),idv
        endif
      else
!       write(6,"('wrlbc: vname=',a,' idv on input=',i3)") 
!    |    trim(vname),idv
      endif
!
      fglb(:,:) = flbc(3:nlon+2,:)
      istat = nf_put_vara_double(ncid,idv,start_3d,count_3d,fglb)
      if (istat /= NF_NOERR) then
        write(6,"('>>> wrlbc: error writing ',a,': itime=',i4,
     |    ' type=',a,' global (nlon,nlat) min,max=',2f8.2)") 
     |    trim(vname),itime,trim(type),minval(fglb),maxval(fglb)
        write(char80,"('Error return from nf_put_vara_double for ',
     |    'var ',a,': itime=',i2,' type=',a)") trim(vname),itime,
     |    trim(type)
        call handle_ncerr(istat,char80)
      else
!       write(6,"('wrlbc: wrote ',a,': itime=',i4,' type=',a,
!    |    ' global (nlon,nlat) min,max=',2f8.2)") trim(vname),
!    |    itime,trim(type),minval(fglb),maxval(fglb)
      endif

!     if (trim(vname)=='TLBC') then
!       do j=1,nlat
!         write(6,"('wrlbc: j=',i3,' tlbc(:,j)=',/,(6e12.4))") 
!    |      j,fglb(:,j)
!       enddo
!     endif

      end subroutine wrlbc
!-----------------------------------------------------------------------
      subroutine wrf4d(ncid,name,itime,fake,ix)
!
! Write prognostic field type f%data to current netcdf output 
! history file attached to ncid:
!
      use input_module,only: dynamo
!
! Args:
      integer,intent(in) :: ncid,itime,ix
      character(len=*),intent(in) :: name
      logical,intent(in) :: fake
!
! Local:
      integer :: nxk,k,i,j,idvar,istat,lonbeg,lonend,nnans,nanfatal
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
          write(char120,"('wrf4d: Error return from nf_put_vara_double',
     +      ' for fake field var ',a,' itime=',i2)") name,itime
          call handle_ncerr(istat,char120)
        endif
        return
      endif
!
!  Electric potential in f4d(i_poten) has been updated by sub
!  prepare_phig3d in pdynamo.F.
!
      do j=1,nlat
        do k=1,nlevp1
          f3d(:,j,k) = foutput(k,3:nlon+2,j,ix)
        enddo
      enddo
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
!
! Check for NaNs:
      nanfatal = 1 ! abort if nans are found (this is normal default)
!     nanfatal = 0 ! temporary: allow NaNs for debug
      if (check_nan) 
     |  call check_nans(f3d,nlon,nlat,nlevp1,name(1:8),nnans,
     |    0,0.,1,nanfatal)
      end subroutine wrf4d
!-----------------------------------------------------------------------
      subroutine wrf3d(ncid,f,itime,name,idv)
!
! Write global 3d field type f(k,i,j) to current netcdf output history 
! file attached to ncid:
!
! Args:
      integer,intent(in) :: ncid,itime,idv
      real :: f(nlevp1,nlonp4,nlat)
      character(len=*),intent(in) :: name
!
! Local:
      integer :: nxk,k,i,j,istat,lonbeg,lonend,nnans,nanfatal
      real :: f3d(nlon,nlat,nlevp1) ! note i,j,k for history
      character(len=120) :: char120
      real :: fmin,fmax
!
! Transform from (k,i,j) to (i,j,k) for the history.
! f3d(nlon,nlat,nlevp1)            ! note i,j,k for history
!
      do j=1,nlat
        do k=1,nlevp1
          f3d(:,j,k) = f(k,3:nlon+2,j)
        enddo
      enddo
!
! Put data onto netcdf file (this is where the majority of output
!   i/o happens). 
!
      start_4d(1:3) = 1
      start_4d(4) = itime
      count_4d(1) = nlon
      count_4d(2) = nlat
      count_4d(3) = nlevp1
      count_4d(4) = 1
      istat = nf_put_vara_double(ncid,idv,start_4d,count_4d,f3d)
      if (istat /= NF_NOERR) then
        write(char120,"('wrf3d: Error return from nf_put_vara_double',
     +    ' for field var ',a,' itime=',i2)") name,itime
        call handle_ncerr(istat,char120)
      endif
!
! Write min,max of each field if not echoing source history:
!     if (istep /= 0) then
!       call fminmaxspv(f3d,nlon*nlev*nlat,fmin,fmax,spval)
!       write(6,"('wrf3d: Wrote field ',a,' istep=',i3,
!    |    ' 3d min,max=',2e12.4)") name,istep,fmin,fmax
!     endif
!
! Check for NaNs:
      nanfatal = 1 ! abort if nans are found
      if (check_nan) 
     |  call check_nans(f3d,nlon,nlat,nlevp1,name,nnans,
     |    0,0.,1,nanfatal)
      end subroutine wrf3d
!-----------------------------------------------------------------------
      subroutine nc_fileinfo(ncid,fileinfo)
!
! A netcdf file on ncid has been written and is ready to be
! closed and disposed. This routine builds a fileinfo string
! giving beginning and ending history dates, model version string,
! and "primary" or "secondary". This is added as a global
! file attribute, and is included in the mss comment field
! if the file is disposed to the hpss, and is therefore accessible 
! via "hsi ls -lA". Example fileinfo string:
!
! 96080  80, 0, 0 to 96084  84, 0, 0 by  720 tgcm22 primary
!
! Args:
      integer,intent(in) :: ncid
      character(len=*),intent(out) :: fileinfo
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
!     write(6,"('nc_fileinfo: fileinfo=',a)") trim(fileinfo)
      end subroutine nc_fileinfo
!-----------------------------------------------------------------------
      subroutine rdfilename(ncid,itime,idv,vname,filelen,filename)
      use input_module,only: mxlen_filename
!
! Args:
      integer,intent(in) :: ncid,itime,idv,filelen
      character(len=*),intent(in) :: vname
      character(len=*),intent(out) :: filename ! (mxlen_filename)
!
! Local:
      integer :: i,ii,istat,len,id
      integer :: index(2)
      character(len=filelen) :: nameonfile
!
      filename = ' ' ! init
!
! Get name on file at current itime (char (filelen,ntime)):
      nameonfile = ' '
      index(2) = itime
      do ii=1,filelen
        index(1) = ii
        istat = nf_get_var1_text(ncid,idv,index,nameonfile(ii:ii))
      enddo
      if (istat /= NF_NOERR) 
     |  call handle_ncerr(istat,'rdfilename nf_get_var1_text')
!
! Transfer name to output filename:
      len = len_trim(nameonfile)
      if (len <= mxlen_filename) then
        filename = trim(nameonfile)
      else ! len on file is larger than current mxlen_filename
        filename = nameonfile(1:mxlen_filename)
        write(6,"('>>> Warning rdfilename: filelen on file=',i4,
     |    ' len_trim(',a,')=',i4,' but mxlen_filename=',i4,
     |    ' -> file name may be truncated')") 
     |    filelen,vname,len,mxlen_filename
      endif
      if (len_trim(filename)==0) filename = '[none]'
!     write(6,"('rdfilename returning ',a,' = ',a)") 
!    |  vname,trim(filename)
      end subroutine rdfilename
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      end module nchist_module
