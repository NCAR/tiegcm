!
      module hist_module
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
     |  nfsech_mag,     ! number of sech fields on magnetic grid
     |  nfsech_geo,     ! number of sech fields on geographic grid
     |  nfsech_geo2d,   ! number of sech fields on geographic 2d grid
     |  nfsech_mag2d,   ! number of sech fields on magnetic 2d grid
     |  nfsech_magphr   ! number of sech fields on magnetospheric 2d grid
!
! History structure type:
!
      type history
        character(len=16) ::
     |    model_version,     ! version string (e.g., tgcm24)
     |    model_name         ! model name (e.g., time-gcm)
        character(len=16) ::
     |    type,              ! 'primary' or 'secondary'
     |    host,              ! host machine (getenv of 'HOST')
     |    system,            ! operating system (pre-proc macro)
     |    logname            ! user login (getenv of 'LOGNAME')
        character(len=24) ::
     |    rundate,           ! local date of run
     |    createdate,        ! date history file was created
     |    writedate          ! most recent date file was written to
        character(len=80) ::
     |    mss_path,          ! mss path to file containing this history
     |    mss_source,        ! mss path to source history file
     |    mss_secsource      ! mss path to secsource history file
        integer ::
     |    nlat,              ! number of latitude grid points
     |    nlon,              ! number of longitude grid points
     |    nlev,              ! number of vertical levels
     |    ihist,             ! ihist = nth history on current file
     |    delhmins,          ! delta minutes between histories
     |    modeltime(4),      ! model time (day,hour,minute,seconds)
     |    iter,              ! iteration (# steps from 0,0,0)
     |    year,              ! 4-digit year
!    |    calday,            ! calendar day (1-367)
     |    day,               ! calendar day (1-367)
     |    calendar_advance,  ! 0/1
     |    step,              ! time step (seconds)
     |    ncep,              ! nmc/ncep flag (0 or 1)
     |    gpi,               ! gpi flag (0 or 1)
     |    gswmdi,            ! gswm flag diurnal tides(0 or 1)
     |    gswmsdi,           ! gswm flag semidiurnal tides (0 or 1)
     |    gswmnmdi,          ! gswm flag nonmigrationg diurnal tides (0 or 1)
     |    gswmnmsdi,         ! gswm flag nonmigrationg semidiurnal tides (0 or 1)
     |    initial_year,      ! initial starting year
     |    initial_day,       ! initial starting day
     |    initial_mtime(3),  ! initial starting model time
     |    source_mtime(3)    ! start time of source history
        real :: 
     |    ut,                ! ut (decimal hours, from modeltime)
     |    time,              ! model time ("minutes since yyyy-m-d")
     |    missing_value,     ! missing value (usually spval)
     |    zptop,             ! top Zp level
     |    zpbot,             ! bottom Zp level
     |    p0,                ! reference pressure
     |    hpower,            ! hemispheric power
     |    ctpoten,           ! cross-cap potential
     |    byimf,             ! by imf
     |    f107d,             ! daily f10.7 cm flux
     |    f107a,             ! 81-day average f10.7 cm flux
     |    mag(2,2),          ! magnetic pole coordinates
     |    dtide(2),          ! diurnal tide (tide2(2) from input)
     |    sdtide(10),        ! semi-diurnal tide (tide(10) from input)
     |    colfac,            ! collision factor
     |    alfa30,e30,alfad2,ed2 ! auroral parameters from input
        logical dynamo       ! if true, is dynamo history
        logical coupled_ccm  ! if true, run was coupled with ccm
        integer :: nflds,    ! total number of fields on history
                             ! nflds = nfgeo+nfmag+nfgeo2d+nfmag2d+nfmagphr
     |             nfgeo,    ! number of fields on geographic grid
     |             nfmag,    ! number of fields on magnetic grid
     |             nfgeo2d,  ! number of fields on geographic grid 2d
     |             nfmag2d,  ! number of fields on magnetic grid 2d
     |             nfmagphr  ! number of fields on magnetospheric grid
        character(len=16),
     |    pointer :: fnames(:) ! pointer to field names on this history
      end type history
!
! Global history structure variable:
      type(history) :: h
      contains
!-------------------------------------------------------------------
      subroutine hist_initype(h,istep)
      use input_module,only: start_year,start_day,start,
     |  calendar_advance,source_start
!
! Initialize a history structure:
!
      type(history),intent(out) :: h
      integer,intent(in) :: istep
!
      h%model_name = ' '
      h%model_version = ' '
      h%logname = ' '
      h%host = ' '
      h%system = ' '
      h%type = "unknown"
      h%rundate = ' '
      h%createdate = ' '
      h%writedate = ' '
      h%mss_source = ' '
      h%mss_secsource = ' '
      h%mss_path = ' '
!
! If this is an initial run, set initial calendar time to model 
! start time. Otherwise (continuation run), the initial time has 
! been read by rdsource (see nchist.F)
!
      if (nsource > 0) then
        h%initial_year  = start_year
        h%initial_day   = start_day 
        h%initial_mtime = start(:,1)
        h%source_mtime  = source_start
      endif

      h%ihist = -1
      h%delhmins = 0
      h%modeltime = (/-1,-1,-1,-1/)
      h%nlat = 0
      h%nlon = 0
      h%nlev = 0
      h%zptop = 0.
      h%zpbot = 0.
      h%time = -1
      h%iter = -1
      h%year = -1
      h%day = -1
      h%calendar_advance = calendar_advance
      h%ut   = -1.
      h%step = -1
      h%ncep =  0
      h%gpi  =  0
      h%gswmdi    = 0
      h%gswmsdi   = 0
      h%gswmnmdi = 0
      h%gswmnmsdi= 0
      h%nflds  =-1
      h%dynamo = .false.
      h%coupled_ccm = .false.
!
! If istep==0, h%f107, etc were defined from source history: 
      if (istep > 0) then
        h%hpower = -1.
        h%ctpoten = -1.
        h%f107d = -1.
        h%f107a = -1.
      endif
      h%byimf = -1.
      h%mag(:,1) = (/-1.,-1./)
      h%mag(:,2) = (/-1.,-1./)
      h%dtide = (/-1.,-1./)
      h%sdtide = (/-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1./)
      h%colfac = -1.
      h%p0     = -1.
      h%alfa30 = -1.
      h%e30    = -1.
      h%alfad2 = -1.
      h%ed2    = -1.
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
      if (h%type(1:3)=='pri') then
        if (rdwr=='WRITE') then
          write(6,"('Write TGCM PRIMARY HISTORY (diskfile=',a,'):')")
     |      trim(diskfile)
        else
          write(6,"('Read TGCM PRIMARY HISTORY (source history)',
     |       ' (diskfile=',a,')')") trim(diskfile)
        endif
      elseif (h%type(1:3)=='sec') then
        if (rdwr=='WRITE') then
          write(6,"('Write TGCM SECONDARY HISTORY (diskfile=',a,'):')")
     |      trim(diskfile)
        else ! don't generally read sech, but you never know..
          write(6,"('Read TGCM SECONDARY HISTORY (diskfile=',a,'):')")
     |      trim(diskfile)
          write(6,"('WARNING: Why am I reading a secondary history?')")
          write(6,"(9x,'The model will not start from secondary ',
     |      'history.')")
        endif
      else
        write(6,"('>>> print_hist: unknown type (apparently ',
     |    'not primary or secondary.')")
      endif
!
      write(6,"(2x,'model_name    = ',a)") trim(h%model_name)

      write(6,"(2x,'model_version = ',a)") trim(h%model_version)

      write(6,"(2x,'rundate    = ',a)") trim(h%rundate)
!     write(6,"(2x,'createdate = ',a)") trim(h%createdate)
!     write(6,"(2x,'writedate  = ',a)") trim(h%writedate)
      write(6,"(2x,'logname    = ',a)") trim(h%logname)
      write(6,"(2x,'host       = ',a)") trim(h%host)
      write(6,"(2x,'system     = ',a)") trim(h%system)
      write(6,"(2x,'mss_path   = ',a)") trim(h%mss_path)
      write(6,"(2x,'mss_source = ',a)") trim(h%mss_source)
      write(6,"(2x,'mss_secsource = ',a)") trim(h%mss_secsource)
      write(6,"(2x,'source_mtime = ',3i4)") h%source_mtime
      write(6,"(2x,'initial_year = ',i4)") h%initial_year
      write(6,"(2x,'initial_day  = ',i4)") h%initial_day
      write(6,"(2x,'initial_mtime= ',3i4)") h%initial_mtime
      write(6,"(2x,'type       = ',a)") trim(h%type)
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
      write(6,"(2x,'zptop   = ',f6.2,' (top Zp level)')") h%zptop
      write(6,"(2x,'zpbot   = ',f6.2,' (bottom Zp level)')") h%zpbot
      write(6,"(2x,'mag     = ',4f7.2,' (magnetic pole coords)')") 
     |  h%mag
      write(6,"(2x,'dtide   = ',e8.1,' ',f5.1,
     |  ' (amp/phase of diurnal tide)')") h%dtide
      write(6,"(2x,'sdtide  = ',5e8.1,' ',5f5.1,/,
     |  4x,'(amp/phase of semi-diurnal tide)')") h%sdtide
      write(6,"(2x,'f107d   = ',f6.2,' (daily solar flux)')") h%f107d
      write(6,"(2x,'f107a   = ',f6.2,' (average solar flux)')") h%f107a
      write(6,"(2x,'hpower  = ',f5.2,' (Gw)')") h%hpower
      write(6,"(2x,'ctpoten = ',f5.2,' (Volts)')") h%ctpoten
      write(6,"(2x,'byimf   = ',f5.2)") h%byimf
      write(6,"(2x,'alfa30  = ',e8.2,' (KeV) (same as alfa_sp)')") 
     |  h%alfa30
      write(6,"(2x,'e30     = ',e8.2,' (ergs/cm2/s) (same as e_sp)')") 
     |  h%e30
      write(6,"(2x,'alfad2  = ',e8.2,' (KeV)')") h%alfad2
      write(6,"(2x,'ed2     = ',e8.2,' (ergs/cm2/s)')") h%ed2
      write(6,"(2x,'colfac  = ',f5.2)") h%colfac
      write(6,"(2x,'p0      = ',e8.2)") h%p0
      write(6,"(2x,'ncep    = ',i2,' (NMC/NCEP flag)')") h%ncep
      write(6,"(2x,'gpi     = ',i2,' (GPI flag)')") h%gpi
      write(6,"(2x,'gswmdi  = ',i2,' (GSWM flag diurnal tides)')") 
     |  h%gswmdi
      write(6,"(2x,'gswmsdi = ',i2,' (GSWM flag semidiurnal tides)')") 
     |  h%gswmsdi
      write(6,"(2x,'gswmnmdi = ',i2,' (GSWM flag nonmigrating',
     |  ' diurnal tides)')") 
     |  h%gswmnmdi
      write(6,"(2x,'gswmnmsdi = ',i2,' (GSWM flag nonmigrating',
     |  ' semidiurnal tides)')") 
     |  h%gswmnmsdi
      write(6,"(2x,'nflds   = ',i3,' (number of model fields)')")
     |  h%nflds
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
      use params_module,only: mxseries,mxhvols,ispval
      use input_module,only: start,stop,step,source,output,mkhvols,
     |  hist,save,secsource,secstart,secstop,sechist,secsave,secout,
     |  secflds,secfmag,secfgeo2d,secfmag2d,secfmagphr,mxhist_prim,
     |  mxhist_sech
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
      write(6,"('hist_init: nstep=',i3)") nstep
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
     |    size(secsave) -count(secsave==ispval)  +
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
! Number of fields on magnetic grid to be written to secondary histories:
      nfsech_mag = size(secfmag)-count(len_trim(secfmag)==0)
!
! Number of fields on geographic grid to be written to secondary histories:
      nfsech_geo = size(secflds)-count(len_trim(secflds)==0)
!
! Number of fields on geographic 2d grid to be written to secondary histories:
      nfsech_geo2d = size(secfgeo2d)-count(len_trim(secfgeo2d)==0)
!
! Number of fields on magnetic 2d grid to be written to secondary histories:
      nfsech_mag2d = size(secfmag2d)-count(len_trim(secfmag2d)==0)
!
! Number of fields on magnetospheric grid to be written to secondary histories:
      nfsech_magphr = size(secfmagphr)-count(len_trim(secfmagphr)==0)
!
! Number of secondary files to be written (nfiles_sech), and
! total secondary histories to be written (nsech_total):
!
      nfiles_sech = numfiles('sech',nseries_sech,1,0,nsech_total)
!
      end subroutine hist_init
!-------------------------------------------------------------------
      end module hist_module
