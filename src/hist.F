!
      module hist_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use input_module,only: mxlen_filename
      use params_module,only: spval,ispval
      implicit none
!
! History module contains non-user input variables to control
! primary and secondary history i/o, and defines the history
! structure type "history". A global variable "h" is declared
! type(history).
! 
! Contained subroutines:
!   sub hist_initype: initialize a history structure
!   sub hist_print  : print a history structure to stdout
!   sub hist_init   : initialize history control variables
!
! Procedure to add a new variable to the histories:
!   1. Add declarations to the history type below.
!   2. Add initialization of the new history structure component
!      in sub hist_initype in this module.
!   3. Add write of new var to stdout in hist_print in this module.
!   4. Give value to the new structure component in sub define_hist
!      (in output.f). The value may come from another module
!      (e.g., "use input_mod"), or an included header file.
!      (Define_hist is called prior to writing each history to disk).
!   5. Add definition of the new variable to the netcdf history file 
!      in sub nc_define (nchist_mod.F) (long name, units, etc)
!      (declare variable id at top of nchist_mod.F, e.g., idv_ut,
!       and make nc call according to dimensions and data type).
!      (nc_define is called when a new netcdf history file is created)
!   6. Add value to the netcdf variable in nc_wrhist (nchist_mod.F)
!      using the history structure component (e.g., h%ut) defined
!      in step 4 above.
!   7. Add var to case statement in sub nc_rdhist (nchist_mod.F)
!      for reading the variable from the source history.
! Make a short run to write the new variable to a history, which
!   can then be used as a source history in subsequent runs.
!
! History parameters:
!
      integer,parameter :: 
     |  iprint = 0      ! print level flag
      logical :: netcdf ! true if writing netcdf histories
!
! Primary history related variables (not user input):
      integer ::
     |  nseries_prim,   ! number of primary history time series (1-mxseries)
     |  iseries_prim,   ! index to primary history time series
     |  nsource,	! number of source files given (0 or 1)
     |  nfiles_prim,    ! number of history output files to be written
     |  nstep,          ! total number of steps this run
     |  ncid,           ! netcdf file id (ignored if .not.netcdf)
     |  nhist,          ! number of prim histories written to current hist file
     |  nhist_total,    ! number of prim histories written during run
     |  ioutfile,       ! output(ioutfile) is current open hist file name
     |  modeltime(4)    ! current model time (day,hrs,mins,secs)
!
! Secondary history related variables (not user input):
! Primary histories can be written in up to mxseries separate
!   time sequences. Secondary histories can be written in up to 
!   mxseries_sech time sequences. (see params.h)
!
      integer ::
     |  isechist,       ! secondary history flag (global 0 or 1)
     |  nseries_sech,   ! number of secondary hist time series (<=mxseries_sech)
     |  iseries_sech,   ! index to current sech time series
     |  istep_sech,     ! number of steps since beginning of sech series
     |  nfiles_sech,    ! number of history files to be written
     |  ncidsech,       ! netcdf sech file id (ignored if .not.netcdf)
     |  nsech,          ! number of sech histories written to current file
     |  nsecsource,     ! number of secsource files given (for mhd model)
     |  nsech_total,    ! number of sech histories written during run
     |  isecout,        ! secout(isecout) is current open hist file name
     |  nfsech          ! number of secondary history fields (SECFLDS)
!
! History structure type:
!
      type history
        character(len=80) ::
     |    label              ! user-defined job comment (from namelist)
        character(len=16) ::
     |    model_version,     ! version string (e.g., tgcm24)
     |    model_name         ! model name (e.g., time-gcm)
        character(len=16) ::
     |    hist_type,         ! 'primary' or 'secondary'
     |    run_type,          ! 'initial' or 'continuation'
     |    host,              ! host machine (getenv of 'HOST')
     |    system,            ! operating system (pre-proc macro)
     |    logname            ! user login (getenv of 'LOGNAME')
        character(len=24) ::
     |    create_date,       ! date history file was created
     |    write_date         ! date current history was written
        character(len=mxlen_filename) ::
     |    source_file,       ! source file (this run)
     |    initial_file,      ! initial file (initial run)
     |    output_file,       ! output file (from namelist read)
     |    see_ncfile,        ! SEE data file
     |    gpi_ncfile,        ! gpi data file
     |    ncep_ncfile,       ! ncep data file (timegcm only)
     |    imf_ncfile,        ! IMF data file
     |    gswm_mi_di_ncfile, ! GSWM migrating diurnal data file
     |    gswm_mi_sdi_ncfile,! GSWM migrating semi-diurnal data file
     |    gswm_nm_di_ncfile, ! GSWM non-migrating diurnal data file
     |    gswm_nm_sdi_ncfile,! GSWM non-migrating semi-diurnal data file
     |    saber_ncfile,      ! SABER U,Z lbc data file
     |    tidi_ncfile,       ! TIDI U,V lbc data file
     |    ctmt_ncfile,       ! CTMT lbc data file
     |    bgrddata_ncfile    ! Background data file
        integer ::
     |    nlat,              ! number of latitude grid points
     |    nlon,              ! number of longitude grid points
     |    nlev,              ! number of vertical levels
     |    ihist,             ! ihist = nth history on current file
     |    delhmins,          ! delta minutes between histories
     |    modeltime(4),      ! model time (day,hour,minute,seconds)
     |    iter,              ! iteration (# steps from 0,0,0)
     |    year,              ! 4-digit year
     |    day,               ! calendar day (1-367)
     |    calendar_advance,  ! 0/1
     |    step,              ! time step (seconds)
     |    initial_day,       ! initial starting day
     |    initial_year,      ! initial starting year
     |    initial_mtime(3),  ! initial starting model time
     |    source_mtime(3),   ! start time of source history
     |    tuv_lbc_intop,     ! 1 (old histories), or 0 (new histories)
     |    ntask_mpi,         ! number of mpi tasks this run
     |    coupled_cmit       ! 1 if coupled with CISM/CMIT, 0 otherwise
        real :: 
     |    ut,                ! ut (decimal hours, from modeltime)
     |    time,              ! model time ("minutes since yyyy-m-d")
     |    missing_value,     ! missing value (usually spval)
     |    zmtop,             ! top midpoint level
     |    zmbot,             ! bottom midpoint level
     |    zitop,             ! top interface level
     |    zibot,             ! bottom interface level
     |    p0,                ! reference pressure (5.e-7 mb (or hPa))
     |    p0_model,          ! reference pressure used in the model (5.e-4 ubars)
     |    grav,              ! gravity constant used in model
     |    hpower,            ! hemispheric power
     |    ctpoten,           ! cross-cap potential
     |    bximf,             ! bx imf
     |    byimf,             ! by imf
     |    bzimf,             ! bz imf
     |    swvel,             ! solar wind velocity
     |    swden,             ! solar wind density
     |    al,                ! solar activity index
     |    kp,                ! solar activity index
     |    f107d,             ! daily f10.7 cm flux
     |    f107a,             ! 81-day average f10.7 cm flux
     |    dtide(2),          ! diurnal tide (tide2(2) from input)
     |    sdtide(10),        ! semi-diurnal tide (tide(10) from input)
     |    colfac,            ! collision factor
     |    e1,e2,             ! peak auroral energy flux
     |    h1,h2,             ! gaussian half-width for aurora
     |    alfac,ec,          ! aurora cusp energy
     |    alfad,ed,          ! aurora drizzle energy
     |    lbc,               ! lower boundary interface level
     |    joulefac           ! joule heating factor
        logical dynamo       ! if true, is dynamo history
        logical coupled_ccm  ! if true, run was coupled with ccm
        integer :: nflds,    ! total number of fields on history
                             ! nflds = nfgeo+nfmag+nfgeo2d+nfmag2d
     |             nfgeo,    ! number of fields on geographic grid
     |             nfmag,    ! number of fields on magnetic grid
     |             nfgeo2d,  ! number of fields on geographic grid 2d
     |             nfmag2d   ! number of fields on magnetic grid 2d
        character(len=16),
     |    pointer :: fnames(:) ! pointer to field names on this history
      end type history
!
! Global history structure variable:
      type(history) :: 
     |  h,	! output histories (primary or secondary)
     |  sh	! source history (input)
!
! Initialization values:
      real :: real_init = spval ! initialization value for real vars
      contains
!-------------------------------------------------------------------
      subroutine hist_initype(h,istep)
      use params_module,only: zibot
      use input_module,only: start_year,start_day,start,
     |  calendar_advance,source_start,source
!
! Initialize a history structure:
!
      type(history),intent(out) :: h
      integer,intent(in) :: istep
!
      h%label = ' '
      h%model_name = ' '
      h%model_version = ' '
      h%logname = ' '
      h%host = ' '
      h%system = ' '
      h%hist_type = ' '
      h%run_type = ' '
      h%create_date = ' '
      h%write_date = ' '
      h%source_file = ' '    ! set by sub output
      h%output_file = ' '    ! set by sub output
!
! If this is an initial run, set initial calendar time to model 
! start time. Otherwise (continuation run), the initial time has 
! been read by rdsource (see nchist.F)
!
      if (nsource > 0) then
        h%initial_year  = start_year
        h%initial_day   = start_day 
        h%initial_mtime = start(:,1)
        h%initial_file = source
        h%source_mtime  = source_start
      endif

      h%ihist = -1
      h%delhmins = 0
      h%modeltime = (/-1,-1,-1,-1/)
      h%nlat = 0
      h%nlon = 0
      h%nlev = 0
      h%zmtop = 0.
      h%zmbot = 0.
      h%zitop = 0.
      h%zibot = 0.
      h%time = -1
      h%iter = -1
      h%year = -1
      h%day = -1
      h%calendar_advance = calendar_advance
      h%ut   = real_init
      h%step = -1

      h%see_ncfile = '[none]'
      h%gpi_ncfile = '[none]'
      h%ncep_ncfile = '[none]'

      h%gswm_mi_di_ncfile  = '[none]'
      h%gswm_mi_sdi_ncfile = '[none]'
      h%gswm_nm_di_ncfile  = '[none]'
      h%gswm_nm_sdi_ncfile = '[none]'

      h%saber_ncfile = '[none]'
      h%tidi_ncfile = '[none]'
      h%bgrddata_ncfile = '[none]'
      h%ctmt_ncfile = '[none]'

      h%nflds  = 0
      h%dynamo = .false.
      h%coupled_ccm = .false.
      h%tuv_lbc_intop = 0 ! this version writes "new" histories
      h%lbc = zibot
      h%ntask_mpi = 0
#if defined(INTERCOMM) || defined(CISMAH)
      h%coupled_cmit = 1
#else
      h%coupled_cmit = 0
#endif
!
! If istep==0, h%f107, etc were defined from source history: 
      if (istep > 0) then
        h%hpower  = real_init
        h%ctpoten = real_init
        h%kp      = real_init
        h%f107d   = real_init
        h%f107a   = real_init
      endif
      h%byimf = real_init
      h%bzimf = real_init
      h%swvel = real_init
      h%swden = real_init
      h%al = real_init
      h%kp = real_init
      h%dtide = (/real_init,real_init/)
      h%sdtide(:) = real_init
      h%colfac = real_init
      h%joulefac = real_init
      h%p0     = real_init
      h%p0_model = real_init
      h%grav   = real_init

!     h%alfa30 = real_init
!     h%e30    = real_init
!     h%alfad2 = real_init
!     h%ed2    = real_init

      h%e1 = real_init
      h%e2 = real_init
      h%h1 = real_init
      h%h2 = real_init
      h%alfac = real_init
      h%ec = real_init
      h%alfad = real_init
      h%ed = real_init
!
      if (associated(h%fnames)) deallocate(h%fnames)
      end subroutine hist_initype
!-------------------------------------------------------------------
      subroutine hist_print(h,rdwr,diskfile)
!
! Print info about a history structure:
!
! Args:
      type(history),intent(in) :: h
      character(len=*),intent(in) :: rdwr,diskfile
!
! Locals:
      integer :: i
!
      write(6,"(/,72('-'))")
      if (h%hist_type(1:3)=='pri') then
        if (rdwr=='WRITE') then
          write(6,"('Write TGCM PRIMARY HISTORY',/,'Diskfile: ',a)") 
     |      trim(diskfile)
        else
          write(6,"('Read TGCM PRIMARY HISTORY (source history)',
     |       /,'Diskfile: ',a)") trim(diskfile)
        endif
      elseif (h%hist_type(1:3)=='sec') then
        if (rdwr=='WRITE') then
          write(6,"('Write TGCM SECONDARY HISTORY',/,'Diskfile: ',a)")
     |      trim(diskfile)
        else ! don't generally read sech, but you never know..
          write(6,"('Read TGCM SECONDARY HISTORY',/,'Diskfile=',a)")
     |      trim(diskfile)
          write(6,"('WARNING: Why am I reading a secondary history?')")
          write(6,"(9x,'The model will not start from secondary ',
     |      'history.')")
        endif
      else
        write(6,"('>>> hist_print: unknown type (apparently ',
     |    'not primary or secondary.')")
      endif
!
      write(6,"(2x,'label         = ',a)") trim(h%label)
      write(6,"(2x,'model_name    = ',a)") trim(h%model_name)
      write(6,"(2x,'model_version = ',a)") trim(h%model_version)
      write(6,"(2x,'create_date= ',a)") trim(h%create_date)
      write(6,"(2x,'write_date = ',a)") trim(h%write_date)
      write(6,"(2x,'logname    = ',a)") trim(h%logname)
      write(6,"(2x,'host       = ',a)") trim(h%host)
      write(6,"(2x,'system     = ',a)") trim(h%system)
      write(6,"(2x,'run_type   = ',a)") trim(h%run_type)
      write(6,"(2x,'source_file  = ',a)") trim(h%source_file)
      write(6,"(2x,'output_file  = ',a)") trim(h%output_file)
      write(6,"(2x,'source_mtime = ',3i4)") h%source_mtime
      write(6,"(2x,'initial_year = ',i4)") h%initial_year
      write(6,"(2x,'initial_day  = ',i4)") h%initial_day
      write(6,"(2x,'initial_mtime= ',3i4)") h%initial_mtime
      write(6,"(2x,'type       = ',a)") trim(h%hist_type)
      write(6,"(2x,'ihist      = ',i3,
     |  ' (nth history on history file)')") h%ihist
      write(6,"(2x,'delhmins= ',i4,
     |  ' (delta minutes between histories)')") h%delhmins
      write(6,"(2x,'calendar year,day = ',i4,',',i3)") h%year,h%day
      if (h%calendar_advance<=0) then
        write(6,"(4x,'(model is NOT being advanced in calendar time)')")
      else
        write(6,"(4x,'(model IS being advanced in calendar time)')")
      endif
      write(6,"(2x,'modeltime = ',i3,',',i2,',',i2,',',i2,
     |  ' (model time day,hour,minute,seconds)')") h%modeltime
      if (h%time <= 1440.) then
        write(6,"(2x,'time    = ',f10.2,' (minutes in current day)')") 
     |    h%time
      else
        write(6,"(2x,'time    = ')")
      endif
      write(6,"(2x,'ut      = ',f5.2,' (ut hours)')") h%ut
      write(6,"(2x,'step    = ',i4,' (time step in seconds)')")
     |  h%step
      write(6,"(2x,'iter    = ',i8,' (number of steps from 0,0,0)')")
     |  h%iter
      write(6,"(2x,'nlat    = ',i4,' (number of latitudes)')") h%nlat
      write(6,"(2x,'nlon    = ',i4,' (number of longitudes)')") h%nlon
      write(6,"(2x,'nlev    = ',i4,' (number of levels)')") h%nlev
      write(6,"(2x,'zmtop   = ',f7.3,' (top midpoint level)')") h%zmtop
      write(6,"(2x,'zmbot   = ',f7.3,' (bottom midpoint level)')") 
     |  h%zmbot
      write(6,"(2x,'zitop   = ',f7.3,' (top interface level)')") h%zitop
      write(6,"(2x,'zibot   = ',f7.3,' (bottom interface level)')") 
     |  h%zibot
      write(6,"(2x,'dtide   = ',e8.1,' ',f5.1,
     |  ' (amp/phase of diurnal tide)')") h%dtide
      write(6,"(2x,'sdtide  = ',5e8.1,' ',5f5.1,/,
     |  4x,'(amp/phase of semi-diurnal tide)')") h%sdtide
      write(6,"(2x,'f107d   = ',e12.4,' (daily solar flux)')") h%f107d
      write(6,"(2x,'f107a   = ',e12.4,' (average solar flux)')") h%f107a
      write(6,"(2x,'hpower  = ',e12.4,' (Gw)')") h%hpower
      write(6,"(2x,'ctpoten = ',e12.4,' (Volts)')") h%ctpoten
      write(6,"(2x,'kp      = ',e12.4,' ()')") h%kp
      write(6,"(2x,'byimf   = ',e12.4)") h%byimf
      write(6,"(2x,'bzimf   = ',e12.4)") h%bzimf
      write(6,"(2x,'swvel   = ',e12.4)") h%swvel
      write(6,"(2x,'swden   = ',e12.4)") h%swden
      write(6,"(2x,'al      = ',e12.4)") h%al
      write(6,"(2x,'e1,e2   = ',2e12.4,' (ergs/cm2/s)')") h%e1,h%e2
      write(6,"(2x,'alfac,alfad  = ',2e12.4,' (keV)')") h%alfac,h%alfad
      write(6,"(2x,'ec,ed = ',2e12.4,' (ergs/cm2/s)')") h%ec,h%ed

!     if (h%alfa30 /= real_init)
!    |  write(6,"(2x,'alfa30  = ',e8.2,' (KeV) (same as alfa_sp)')") 
!    |    h%alfa30
!     if (h%e30 /= real_init)
!    |  write(6,"(2x,'e30     = ',e8.2,' (ergs/cm2/s) (same as e_sp)')") 
!    |    h%e30
!     if (h%alfad2 /= real_init) 
!    |  write(6,"(2x,'alfad2  = ',e8.2,' (KeV)')") h%alfad2
!     if (h%ed2 /= real_init)
!    |  write(6,"(2x,'ed2     = ',e8.2,' (ergs/cm2/s)')") h%ed2

      write(6,"(2x,'colfac    = ',e10.2)") h%colfac
      write(6,"(2x,'joulefac  = ',e10.2)") h%joulefac
      write(6,"(2x,'p0        = ',e10.2)") h%p0
      if (h%p0_model /= real_init) write(6,"(2x,'p0_model= ',e10.2)")
     |  h%p0_model
      if (h%grav /= real_init)
     |  write(6,"(2x,'grav    = ',e10.2)") h%grav

!     write(6,"(2x,'gswmdi  = ',i2,' (GSWM flag diurnal tides)')") 
!    |  h%gswmdi
!     write(6,"(2x,'gswmsdi = ',i2,' (GSWM flag semidiurnal tides)')") 
!    |  h%gswmsdi
!     write(6,"(2x,'gswmnmdi = ',i2,' (GSWM flag nonmigrating',
!    |  ' diurnal tides)')") 
!    |  h%gswmnmdi
!     write(6,"(2x,'gswmnmsdi = ',i2,' (GSWM flag nonmigrating',
!    |  ' semidiurnal tides)')") 
!    |  h%gswmnmsdi

      write(6,"(2x,'nflds   = ',i3,' (number of model fields)')")
     |  h%nflds

      write(6,"(2x,'gswm_mi_di_ncfile  = ',a)")
     |  trim(h%gswm_mi_di_ncfile)
      write(6,"(2x,'gswm_mi_sdi_ncfile = ',a)")
     |  trim(h%gswm_mi_sdi_ncfile)
      write(6,"(2x,'gswm_nm_di_ncfile  = ',a)")
     |  trim(h%gswm_nm_di_ncfile)
      write(6,"(2x,'gswm_nm_sdi_ncfile = ',a)")
     |  trim(h%gswm_nm_sdi_ncfile)

      write(6,"(2x,'see_ncfile   = ',a)") trim(h%see_ncfile)
      write(6,"(2x,'gpi_ncfile   = ',a)") trim(h%gpi_ncfile)
      write(6,"(2x,'ncep_ncfile  = ',a)") trim(h%ncep_ncfile)
      write(6,"(2x,'imf_ncfile   = ',a)") trim(h%imf_ncfile)
      write(6,"(2x,'saber_ncfile = ',a)") trim(h%saber_ncfile)
      write(6,"(2x,'tidi_ncfile  = ',a)") trim(h%tidi_ncfile)
      write(6,"(2x,'ctmt_ncfile  = ',a)") trim(h%ctmt_ncfile)
      write(6,"(2x,'bgrddata_ncfile = ',a)") trim(h%bgrddata_ncfile)

      write(6,"(2x,'tuv_lbc_intop=',i3,' (if 1, then lbc of t,u,v',
     |  ' are stored in top k slot (old histories))')") h%tuv_lbc_intop
      write(6,"(2x,'LBC = ',f8.2,' (lower boundary interface ',
     |  'level)')") h%lbc
      write(6,"(2x,'ntask_mpi    =',i4,' (number of MPI tasks)')")
     |  h%ntask_mpi
      write(6,"(2x,'coupled_cmit =',i4,' (1 if coupled with CISM/CMIT,', 
     |  ' 0 otherwise)')") h%coupled_cmit
!
      if (associated(h%fnames)) then
        write(6,"('  There are ',i3,' fields on this history,',
     +    ' as follows:')") h%nflds
        do i=1,h%nflds
          write(6,"(2x,a)",advance="NO") h%fnames(i)(1:8)
          if (mod(i,8)==0.or.i==h%nflds) write(6,"(' ')")
        enddo
      endif
      write(6,"(72('-'),/)")
      return
      end subroutine hist_print
!-------------------------------------------------------------------
      subroutine hist_init
!
! Initialize non-input history variables, using validated namelist 
! user inputs from the input module:
! (this routine is called from init module after input validation)
!
      use params_module,only: mxseries,mxhvols
      use input_module,only: start,stop,step,source,output,mkhvols,
     |  hist,secsource,secstart,secstop,sechist,secout,
     |  secflds
!
! Local:
      integer :: nh,nsteps,i,secs_hist,nsteps_hist,n
      integer(kind=8) :: secs_start(mxseries),secs_stop(mxseries)
      character(len=80) :: hvols(mxhvols)
!
! External:
      integer(kind=8),external :: mtime_to_nsec
      integer,external :: numfiles
!
! Number of primary history time series:
      nseries_prim = (size(start)-count(start==ispval))/3
!
      secs_start(:) = 0
      do i=1,nseries_prim
        secs_start(i) = mtime_to_nsec(start(:,i))
        secs_stop(i)  = mtime_to_nsec(stop(:,i))
      enddo
!
! Total number of timesteps this run:
      nstep = (secs_stop(nseries_prim)-secs_start(1))/step
      write(6,"('hist_init: nstep=',i5)") nstep
!
! Number of source files given (0 or 1):
      if (len_trim(source)==0) then
        nsource = 0
      else
        nsource = 1
      endif
!
! Number of primary files to be written (nfiles_prim), and 
! total primary histories to be written (nhist_total):
      nfiles_prim = numfiles('prim',nseries_prim,nsource,0,nhist_total)
!
! Initialize for beginning of run:
      iseries_prim = 1
      nsech = 0
      iseries_sech = 1
!
! Number of secsource files given (0 or 1):
      if (len_trim(secsource)==0) then
        nsecsource = 0
      else
        nsecsource = 1
      endif
!
! Set secondary history non-input variables:
!
! n = total number of secondary history inputs read:
      n = size(secstart)-count(secstart==ispval) +
     |    size(secstop) -count(secstop==ispval)  +
     |    size(sechist) -count(sechist==ispval)  +
     |    size(secflds) -count(len_trim(secflds)==0)
!
! Number of secondary history time series:
      nseries_sech = (size(secstart)-count(secstart==ispval))/3
!
! isechist: secondary histories flag:
! If no secondary history inputs were read, return.
      isechist = 0 
      if (n > 0) isechist = 1
      if (isechist==0) return
!
! Secondary history init:
      secs_start(:) = 0
      secs_start(:) = 0
      do i=1,nseries_sech
        secs_start(i) = mtime_to_nsec(secstart(:,i))
        secs_stop(i)  = mtime_to_nsec(secstop(:,i))
      enddo
!
      nfsech = size(secflds)-count(len_trim(secflds)==0)
      nfiles_sech = numfiles('sech',nseries_sech,1,0,nsech_total)
      end subroutine hist_init
!-------------------------------------------------------------------
      subroutine copyhist(hsrc,hdest)
!
! Copy history structure hsrc to hdest
!
      implicit none
!
! Args:
      type(history),intent(in) :: hsrc    ! input source history
      type(history),intent(out) :: hdest  ! output destination history
!
! Local:
      integer :: istat
!
! Copy structure wholesale:
      hdest = hsrc
!
! Allocate hdest%fnames and assign from hsrc:
      if (hsrc%nflds > 0) then
        allocate(hdest%fnames(hsrc%nflds),stat=istat)
      else
        write(6,"('>>> copyhist: hsrc%nflds=',i3)") hsrc%nflds
      endif
      hdest%fnames = hsrc%fnames
      write(6,"('Copied history structure at modeltime ',
     |  3i4)") hdest%modeltime(1:3)
      end subroutine copyhist
!-------------------------------------------------------------------
      end module hist_module
