!
      module input_module
      use params_module,only: mxhvols,mxseries,mxseries_sech,mxfsech,
     |  nlat,nlon,nlev,glat1,dlat,glon1,dlon,tgcm_version,tgcm_name,
     |  spval,ispval,mxday,nlonp4,mxind_time
      use mk_hvols,only: mkhvols
      implicit none
!
! Read and validate user inputs via namelist.
! Principle interface to the model is via a single call to sub input 
!   (which makes calls to the input module subroutines), and via 
!   "use input_module" statements in model subprograms.
! This module is dependent on the model header file "params.h"
!
! The model may reference either the module variables or the input_type
!   structure inp (e.g., inp%label, inp%f107, etc). Both methods
!   use "use input_module" statements, e.g.: 
!   use input_module                  ! include entire module
!   use input_module,only: f107,f107a ! include only solar fluxes 
!   use input_module,only: inp        ! include input_type structure
!
! Procedure to add new input parameters:
!   0. Add any necessary parameter declarations in params.h
!      (these are generally for dimensioning input variables)
!   1. Declare new module variables in proper category below
!      (e.g., model-wide, primary histories, secondary histories, etc)
!   2. Declare corresponding components of input_type (use same names)
!   3. Add variables to namelist/tgcm_input/
!   4. Initialize variables in inp_init
!   5. Define tgcm_type components in inp_deftyp (from module variables)
!   6. Add to stdout print in inp_print
!   7. Validate values read in appropriate routine (e.g., inp_hist, etc)
!      (validation may include setting non-input variables related
!       to the input values. These may be referenced from included
!       header files or other modules).
!
! Namelist user input variables:
!
      character(len=80) :: 
     |  label,           ! optional generic text label for this run
     |  tempdir,         ! temporary directory 
     |  magvol,          ! file name or mss path to magnetic data file
     |  amievol          ! file or mss path of amie data file (optional)
!
! date and calday are no longer supported, and are replaced by start_day, 
! start_year, and calendar_advance. Date and calday are retained here so 
! error usage statements can be issued if user sets one of them.
!
      integer ::
     |  start_day,       ! starting day of year (integer 0->365)
     |  start_year,      ! starting year (4-digit integer yyyy)
     |  calendar_advance,! if > 0, advance calendar day from start_day
     |  date(3),         ! old: model starting year, day ( 2 ints yyyy,dd)
     |  calday,          ! old: starting calendar day (0-mxday)
     |  step,            ! model time step (integer seconds)
     |  dispose,         ! dispose output files to mss if dispose==1 or 2
     |  difhor,          ! horizontal eddy diffusion flag 0/1
     |  iuivi,           ! ion drifts in momentum flag (replaced by dynamo flag)
     |  dynamo,          ! 0/1 flag for dynamo
     |  nmc,             ! flag to use NMC boundary conditions 0/1
     |  tideann,         ! 0/1 flag for annual tide
     |  aurora,          ! 0/1 flag for aurora
     |  magphr,          ! 0/1 flag for magnetosphere
     |  ntask_lat,       ! number of tasks in latitude  dimension
     |  ntask_lon        ! number of tasks in longitude dimension
      real ::
     |  mag(2,2),        ! lat,lon of south,north magnetic poles
                         ! (this is vestigal -- can be removed)
     |  tide(10),        ! semidiurnal tide amplitudes and phases
     |  tide2(2),        ! diurnal tide amplitude and phase
     |  tide3m3(2),      ! 2-day wave amplitude and phase
     |  f107,            ! 10.7 cm daily solar flux
     |  f107a,           ! 10.7 cm average (81-day) solar flux
     |  power,           ! hemispheric power (gw) (hpower on histories)
     |  ctpoten,         ! cross-cap potential (volts)
     |  power_time(4,mxind_time),   ! time-dependent hemispheric power input
     |  ctpoten_time(4,mxind_time), ! time-dependent cross tail potential input
     |  byimf,           ! BY component of IMF
     |  bzimf,           ! Bz component of IMF in nT
     |  swvel,           ! Solar wind velocity in km/s
     |  swden,           ! Solar wind density in #/cm3
     |  AL,              ! AL lower magnetic auroral activity index in nT
                         ! if present, ALUSE=true; if absent, AL=-20, ALUSE=false
     |  colfac           ! collision factor
        logical :: aluse ! logical to use AL in Weimer model or not
      character(len=80) ::
     |  potential_model,  ! electric potential model used
                          !  Values can be 'HEELIS', 'WEIMER', or 'NONE'
                          !  If absent, the default value is set to 'HEELIS'
     |  weimer_ncfile,    ! mss path or file path to netcdf weimer coef file
                          ! see comments in weimer_mod.f
     |  weimer_ncfiledef, ! default path to weimer coef data file
! hrindices_ncfile NOT MADE YET! 11/15/02:
     |  hrindices_ncfile, ! mss path or file path to netcdf hourly geophysical
                          ! data indices file.  see comments in hrindices_mod.f
     |  hrindices_ncfiledef, ! default path to hrly indices data file
!
     |  gpi_ncfile,       ! mss path or file path to netcdf gpi data file
                          ! see comments in gpi_mod.f
     |  gswm_di_ncfile,   ! mss path or file path to netcdf gswm data file
                          ! for diurnal tides see comments in gswm_mod.f
     |  gswm_sdi_ncfile,  ! mss path or file path to netcdf gswm data file
                          ! for semidirurnal tides see comments in gswm_mod.f
     |  gswm_nmidi_ncfile,! mss path or file path to netcdf gswm data file
                          ! for nonmigrating diurnal tides see comments in gswm_mod.f
     |  gswm_nmisdi_ncfile,! mss path or file path to netcdf gswm data file
                           ! for nonmigrating semidiurnal tides see comments in gswm_mod.f
!
! Non-input gpi variables:
     |  gpi_ncfiledef     ! default path to gpi data file
      integer,parameter :: ngpivars = 4
      real :: gpi_vars(ngpivars) ! f107,f107a,power,ctpoten
      character(len=16) ::
     |  gpi_names(ngpivars)      ! names of gpi_vars
      integer :: 
     |  ntimes_ctpoten,    ! number of times provided in ctpoten_time
     |  ntimes_power       ! number of times provided in power_time
!
! Primary history user input (dimension parameters are in params.h):
      character(len=80) :: 
     |  source,            ! file containing source history (optional)
     |  output(mxhvols)    ! output file(s) (required)
      integer ::
     |  source_start(3),   ! source history model time
     |  start(3,mxseries), ! primary history model start time(s)
     |  stop(3,mxseries),  ! primary history model stop time(s)
     |  hist(3,mxseries),  ! primary history disk write frequency
     |  save(3,mxseries),  ! primary history file save frequency
     |  mxhist_prim,       ! max number of histories per primary file
     |  msreten,           ! retention period for history files
     |  noutput            ! number of output files given
!
! Secondary history user input (dimension parameters are in params.h):
      character(len=80) :: 
     |  secsource,            ! file containing source sec_history (for mhd)
     |  secout(mxhvols)       ! secondary history output file(s)
      character(len=16) ::
     |  secflds(mxfsech),     ! secondary history output fields
     |  secfmag(mxfsech),     ! secondary output fields on magnetic grid
     |  secfgeo2d(mxfsech),   ! secondary output fields on geographic 2d grid
     |  secfmag2d(mxfsech),   ! secondary output fields on magnetic 2d grid
     |  secfmagphr(mxfsech)   ! secondary output fields on magnetospheric 2d grid
      integer ::
     |  secstart(3,mxseries), ! secondary history model start time(s)
     |  secstop(3,mxseries),  ! secondary history model stop time(s)
     |  sechist(3,mxseries),  ! secondary history disk write frequency
     |  secsave(3,mxseries),  ! secondary history file save frequency
     |  mxhist_sech           ! max number of histories per secondary file
!
! Define input structure type: 
! (Use same names as variables above)
!
! Model-wide input:
      type input_type
        character(len=80) :: 
     |    label, tempdir, magvol, amievol, gpi_ncfile,
     |     gswm_di_ncfile, gswm_sdi_ncfile, gswm_nmidi_ncfile,
     |     gswm_nmisdi_ncfile, potential_model, weimer_ncfile,
     |     hrindices_ncfile
        integer ::
     |    step,start_day,start_year,calendar_advance,
     |    dispose, difhor, iuivi, dynamo, nmc, tideann, aurora, magphr,
     |    ntask_lat,ntask_lon
        real ::
     |    mag(2,2), tide(10),tide2(2),tide3m3(2),
     |    f107,f107a, power,ctpoten,byimf,colfac,bzimf,swvel,swden,al,
     |    power_time(4,mxind_time),ctpoten_time(4,mxind_time)
! 
! Primary histories:
        character(len=80) :: 
     |    source, 
     |    output(mxhvols)
        integer ::
     |    source_start(3),
     |    start(3,mxseries),
     |    stop(3,mxseries),
     |    hist(3,mxseries),
     |    save(3,mxseries),
     |    mxhist_prim,
     |    msreten
!
! Secondary histories:
        character(len=80) :: secsource ! (for mhd model)
        character(len=80) :: secout(mxhvols)
        character(len=16) :: secflds(mxfsech)
        character(len=16) :: secfmag(mxfsech)
        character(len=16) :: secfgeo2d(mxfsech)
        character(len=16) :: secfmag2d(mxfsech)
        character(len=16) :: secfmagphr(mxfsech)
        integer ::
     |    secstart(3,mxseries),
     |    secstop(3,mxseries),
     |    sechist(3,mxseries),
     |    secsave(3,mxseries),
     |    mxhist_sech
      end type input_type
!
! inp is module variable of type input_type:
! (is defined from namelist input vars)
      type(input_type) :: inp
!
! Namelist for read:
      namelist/tgcm_input/ 
     |  label,tempdir,magvol,amievol,date,calday,step,dispose,
     |  source,source_start,output,start,stop,hist,save,
     |  secout,secstart,secstop,sechist,secsave,secflds,
     |  secfmag,secfgeo2d,secfmag2d,secfmagphr,potential_model,
     |  mag,difhor,iuivi,dynamo,nmc,tide,tide2,tide3m3,f107,f107a,
     |  power,ctpoten,byimf,bzimf,swvel,swden,al,colfac,tideann,aurora,
     |  magphr,gpi_ncfile,gswm_di_ncfile,gswm_sdi_ncfile,
     |  gswm_nmidi_ncfile,gswm_nmisdi_ncfile,mxhist_prim,mxhist_sech,
     |  msreten,ntask_lat,ntask_lon,start_day,start_year,
     |  calendar_advance,ctpoten_time,power_time
!
! List of fields that are always written to secondary histories:
      character(len=16) :: secflds_mandatory(4) =
     | (/'TN              ',
     |   'O2              ',
     |   'O1              ',
     |   'Z               '/)
      character(len=16) :: secfmag_mandatory(1) =
     | (/'ZMAG            '/)
!
!
! Current working (execution) directory (set by tgcm.F):
      character(len=120) :: cwd  ! current working directory
!
      contains
!-------------------------------------------------------------------
      subroutine input(mytid,ntask)
!
! Read and validate user namelist input (called from main tgcm.F).
!      
! Args:
      integer,intent(in) :: mytid,ntask
!
! Initialize:
      call inp_init
!
! Do namelist read:
      call inp_read(mytid)
!
! Validate model-wide input:
! (this can be split out into separate routines later)
      call inp_model(ntask)
!
! Validate history input:
      call inp_primhist
      call inp_sechist
!
! Define type(input_typ) inp from validated inputs:
      call inp_deftyp
! 
! Print type(input_typ) inp:
      call inp_print(ntask)
!
      end subroutine input
!-------------------------------------------------------------------
      subroutine inp_init
!
! Initialize input variables:
!
      label   = ' '
      tempdir = ' '
      magvol  = ' '
      amievol = ' '
!
      gpi_ncfile = ' '      ! user input gpi file path
      gpi_ncfiledef = ' '   ! default gpi file path
!     write(gpi_ncfiledef,"('/TGCM/data/gpi_1979-2000.nc')")
!     write(gpi_ncfiledef,"('/TGCM/data/gpi_1979001-2001031.nc')")
!     write(gpi_ncfiledef,"('/TGCM/data/gpi_1979001-2002181.nc')")
      write(gpi_ncfiledef,"('/TGCM/data/gpi_1979001-2002365.nc')")
      gpi_vars(:) = spval
!
      weimer_ncfile = ' '     ! user input weimer coef file path
      weimer_ncfiledef = ' '  ! default weimer coef file path
      write(weimer_ncfiledef,"('/TGCM/data/weimer2001_coeffs.nc')")
!
      hrindices_ncfile = ' '    ! user input hrly indices data file path
      hrindices_ncfiledef = ' ' ! default hrly indices data file path
!     write(hrindices_ncfiledef,
!    |  "('/TGCM/data/hrindices_1988001-2002181.nc')")
!
      gswm_di_ncfile = ' '  ! diurnal tide: user input gswm file path
      gswm_sdi_ncfile= ' '  ! semidiurnal tide: user input gswm file path
      gswm_nmidi_ncfile= ' '  ! nonmigrating diurnal tide: user input gswm file path
      gswm_nmisdi_ncfile= ' ' ! nonmigrating semidiurnal tide: user input gswm file path
!
      date(:) = ispval ! old
      calday  = ispval ! old
!
      start_day  = ispval
      start_year = ispval
      calendar_advance = ispval
!
      step    = ispval
      dispose = ispval
      mag(:,:)= spval
      difhor  = ispval
      iuivi   = ispval
      dynamo  = ispval
      nmc     = 0      ! default is no nmc boundaries
      tide(:) = spval
      tide2(:)= spval
      tide3m3(:) = spval
      tideann = ispval
      aurora  = ispval
      magphr  = ispval
      ntask_lat = ispval
      ntask_lon = ispval
      f107    = spval
      f107a   = spval
      power   = spval
      ctpoten = spval
      power_time(:,:) = spval
      ctpoten_time(:,:) = spval
      byimf   = spval
      bzimf   = spval
      swvel   = spval
      swden   = spval
      al      = spval
      colfac  = spval
!
      potential_model = ' '
      source = ' '
      output = ' '
      source_start(:) = ispval
      start(:,:) = ispval
      stop(:,:)  = ispval
      hist(:,:)  = ispval
      save(:,:)  = ispval
      mxhist_prim = 10 ! default max number of histories per primary file
      if (dlat == 2.5) mxhist_prim = 4 ! 2.5 degree horizontal resolution 
      msreten    = 365 ! default retention period for mss history files
!
      secout(:)  = ' '
      secsource  = ' '
      secflds(:) = ' '
      secfmag(:) = ' '
      secfgeo2d(:)  = ' '
      secfmag2d(:)  = ' '
      secfmagphr(:) = ' '
      secstart(:,:) = ispval
      secstop(:,:)  = ispval
      sechist(:,:)  = ispval
      secsave(:,:)  = ispval
      mxhist_sech = 24 ! default max number of histories per sech file
      if (dlat == 2.5) mxhist_sech = 10 ! 2.5 degree horizontal resolution 
      end subroutine inp_init
!-------------------------------------------------------------------
      subroutine inp_model(ntask)
!
! Args:
      integer,intent(in) :: ntask
!
! Local:
      integer :: i,n,ier
      character(len=16) :: logname
      real :: rday,rhour,rmin,rval
!
! Set temporary directory (settmpdir is in util.F):
!
! Get login name:
      logname = ' '
      call getenv('LOGNAME',logname)
      if (len_trim(logname)==0) then
        write(6,"(/,'>>> INPUT inp_model: Cannot get LOGNAME',
     |    ' environment variable.',/)")
        call shutdown('LOGNAME')
      endif
!
! If tempdir is same as execution directory (cwd), turn off tempdir:
      if (trim(tempdir)=='.'.or.trim(tempdir)==trim(cwd)) then
        write(6,"(/,'NOTE: Temporary scratch directory tempdir (',
     |    a,')')") trim(cwd)
        write(6,"('      is same as execution directory -- will',
     |    ' use only cwd.')")
        tempdir = ' '
      endif
      write(6,"('input: tempdir = ',a)") trim(tempdir)
!
! Model time step (secs):
! If step < 60, it is assumed to be in minutes (enforce step <= 10 
!   minutes), otherwise is in seconds (must be multiple of 60, i.e.,
!   even number of minutes).
! 
      if (step==ispval) then
        write(6,"(/,'>>> INPUT: need model time step STEP ',
     |    '(integer seconds)',/)")
        call shutdown('STEP')
      endif
      if (step <= 0) then
        write(6,"(/,'>>> INPUT: bad time step STEP: ',i5,
     |    ' (must be > 0)',/)") step
        call shutdown('STEP')
      endif
!
! Old starting date and calday (no longer supported):
      if (any(date/=ispval)) then
        write(6,"(/,'>>> INPUT: DATE is no longer supported as an ',
     |    'input parameter.')")
        call usage_calendar
        call shutdown('date')
      endif
      if (calday /= ispval) then
        write(6,"(/,'>>> INPUT: CALDAY is no longer supported ',
     |    'as an input parameter.')")
        call usage_calendar
        call shutdown('calday')
      endif
!
! Verify start_year, start_day, calendar_advance:
! start_day is starting day of the year (1-365)
! start_year is starting year (4-digit yyyy)
! calendar_advance is 0/1 for whether or not to advance calendar time.
!
      if (start_year==ispval) then
        write(6,"(/,'>>> INPUT: need START_YEAR (4-digit integer ',
     |    'starting year.')")
        call usage_calendar
        call shutdown('start_year')
      endif
      if (start_year <= 0) then
        write(6,"(/,'>>> INPUT: bad START_YEAR=',i4,' (must be > 0)')")
     |    start_year
        call shutdown('start_year')
      endif
      if (start_day==ispval) then
        write(6,"(/,'>>> INPUT: need START_DAY (integer calendar ',
     |    'starting day of year (0->365))')")
        call usage_calendar
        call shutdown('start_day')
      endif
      if (start_day <= 0.or.start_day > 365) then
        write(6,"(/,'>>> INPUT: bad START_DAY=',i4,' (must be between',
     |    ' 1 and 365)')") start_day
        call shutdown('start_day')
      endif
      if (calendar_advance==ispval) then
        write(6,"(/,'INPUT NOTE: CALENDAR_ADVANCE was not provided',
     |    ' by input.')")
        write(6,"('Will default to 1 (WILL advance calendar day)',/)")
        calendar_advance = 1
      endif
!
! Magnetic field data file for dynamo:
!   Original cray-blocked file: /ECRIDLEY/ECR90/ECRMG6
!   3-record cray-blocked file: /FOSTER/tgcm/mag.dat (1/99)
!   netcdf file: /TGCM/data/magdat.nc (2/00)
!     (see ~foster/tgcm/mkmag for code that wrote the netcdf file
!      from the original cray-blocked file)
! Should be able to eliminate the need for this external file
!   with Richmond's new dynamo code (summer 00?).
!
! 12/01: Dimension names in /TGCM/data/magdat.nc were changed to match
!        parameter names in tiegcm1, writing new file /TGCM/data/magfield.nc
!        (/TGCM/data/magdat.nc was unchanged). This is still for 5.0h resolution.
!
      n = len_trim(magvol)
      if (n <= 0) then
        write(6,"('Input: did not read value for MAGVOL',
     |    ' (magnetic field data file)')")
        write(magvol,"('/TGCM/data/magfield.nc')")
        write(6,"('  Will use default MAGVOL = ',a)") magvol
      else
        write(6,"('Input: will use magnetic data file ',a)") 
     |    trim(magvol)
      endif
! 
! NMC boundaries for T and Z are not supported in tiegcm:
! If nmc boundaries are not used, zero out znmc and tnmc.
! 5/01: eliminating nmc.h.
!
!     if (nmc > 0) then
!       write(6,"(/,'>>> INPUT: NMC boundaries not supported ',
!    |    'in model ',a,/)") trim(tgcm_version)
!       call shutdown('NMC')
!     else 
!       znmc = 0.
!       tnmc = 0.
!     endif
!
! Tide:
      n = size(tide)-count(tide==spval)
      if (n /= 10) then
        write(6,"(/,'>>> INPUT: must have 10 entries for TIDE,',
     |    ' got ',i3)") n
        call shutdown('TIDE')
      endif
      do i=1,5
        if (tide(i) < 0.) then
          write(6,"(/,'>>> INPUT: amplitudes for TIDE(1:5) must ',
     |      'be > 0.: tide=',/,(5e12.4))") tide
          call shutdown('TIDE')
        endif
      enddo
!
! Tide2:
      n = size(tide2)-count(tide2==spval)
      if (n /= 2) then
        write(6,"(/,'>>> INPUT: must have 2 entries for TIDE2,',
     |    ' got ',i3)") n
        call shutdown('TIDE2')
      endif
      if (tide2(1) < 0.) then
        write(6,"(/,'>>> INPUT: amplitude for TIDE2(1) must ',
     |    'be > 0.: tide2=',e12.4)") tide2
        call shutdown('TIDE2')
      endif
!
! Tide3m3 (in time-gcm only, not in tiegcm):
      n = size(tide3m3)-count(tide3m3==spval)
      if (n > 0) then
        if (tgcm_name /= "time-gcm") then
          write(6,"(/,'>>> INPUT: there is no 2-day wave in ',a)")
     |      trim(tgcm_name)
          write(6,"( '           (tide3m3 is ignored)')")
        else
          if (n /= 2) then
            write(6,"(/,'>>> INPUT: need 2 real values for TIDE3M3')")
            write(6,"(' (amplitude and phase of 2-day wave)')")
            call shutdown('TIDE3M3')
          endif
        endif
      endif
      if (n <= 0 .and. tgcm_name=="time-gcm") then
        write(6,"(/,'>>> WARNING INPUT: did not read input parameter',
     |    'TIDE3M3 -- am setting TIDE3M3 = 0. (no 2-day wave)')")
        tide3m3(:) = 0.
      endif
      if (any(tide3m3==spval)) tide3m3(:) = 0.
!
! Annual tide flag:
      if (tideann==ispval) then
        tideann = 1 ! default is on
      elseif (tideann /= 0 .and. tideann /= 1) then
        write(6,"(/,'>>> INPUT: TIDEANN must be either 0 or 1:',
     |    ' tideann=',i5)") tideann
        call shutdown('TIDEANN')
      endif
!
! Aurora flag:
      if (aurora==ispval) then
        write(6,"(/,'Input: setting default AURORA = 0')")
        aurora = 0
      endif
!
! Magnetospheric flag:
      if (magphr==ispval) then
        write(6,"(/,'Input: setting default MAGPHR = 0')")
        magphr = 0
      endif
!
! Number of tasks in lat,lon dimensions (ntask_lat*ntask_lon must == ntask):
#ifdef MPI
      if (ntask_lat==ispval.or.ntask_lon==ispval) then
        call mkntask(ntask,ntask_lat,ntask_lon,ier)
        if (ier /= 0) then
          write(6,"(/,'>>> INPUT: error from mkntask. ntask_lat=',
     |      i3,' ntask_lon=',i3,' ntask=',i3)") ntask_lat,ntask_lon,
     |      ntask
          call shutdown('MKNTASK')
        else
          write(6,"('Input: mkntask chose ntask_lon=',i3,
     |      ' ntask_lat=',i3,' (ntask=',i3,')')") 
     |      ntask_lon,ntask_lat,ntask
        endif
!       write(6,"(/,'>>> INPUT: for MPI runs, you must specify ',
!    |    'NTASK_LAT as the number of',/,4x,'MPI tasks in the ',
!    |    'latitude dimension, and NTASK_LON as the number of',
!    |    /,4x,'MPI tasks in the longitude dimension.',/)")
!       call shutdown("NTASK_LAT|LON")
      endif
      if (ntask_lat*ntask_lon /= ntask) then
        write(6,"(/,'>>> INPUT: NTASK_LAT * NTASK_LON must ',
     |    'equal the total number of tasks for the run.')")
        write(6,"(4x,'ntask_lat=',i3,' ntask_lon=',i3,
     |    ' ntask=',i4,/)") ntask_lat,ntask_lon,ntask
        call shutdown("NTASK_LAT|LON")
      endif
#endif
!
! Geophysical indices (s.a., /ingpi/ in ingpi.h)
! Note: as of 2/99, requiring only 1 value for each, i.e., not
!       yet allowing time-dependent multiple valued input and
!       not yet using sub getgpi to get indices from mss database.
!
! Daily f10.7 cm flux:
      if (f107 /= spval .and. f107 <= 0.) then
        write(6,"(/,'>>> INPUT: f107 must be positive: f107=',
     |    e12.4)") f107
        call shutdown('F107')
      endif
! c(61) is replaced by f107 from input_mod.
!     c(61) = f107 ! this may get reset later by getgpi or tail
!
! 80-day mean f10.7 flux:
      if (f107a /= spval .and. f107a <= 0.) then
        write(6,"(/,'>>> INPUT: f107a must be positive: f107a=',
     |    e12.4)") f107a
        call shutdown('F107A')
      endif
!
! There are 3 input options for ctpoten and power:
!   1) Constants are provided by the user (ctpoten, power)
!   2) Time series are provided by the user (ctpoten_time, power_time)
!   3) GPI database (neither constans nor time series were provided)
!
! Cannot provide both constant and time series:
      if (ctpoten /= spval .and. any(ctpoten_time /= spval)) then
        write(6,"(/,'>>> INPUT: you cannot provide both a constant',
     |    '(CTPOTEN) and a time-series (CTPOTEN_TIME) for cross-',
     |    'tail potential.')")
        call shutdown('CTPOTEN')
      endif
      if (power /= spval .and. any(power_time /= spval)) then
        write(6,"(/,'>>> INPUT: you cannot provide both a constant',
     |    ' (POWER) and a time-series (POWER_TIME) for hemispheric',
     |    ' power.')")
        call shutdown('POWER')
      endif
!
! Validate cross-tail potential constant, if provided:
      if (ctpoten /= spval .and. ctpoten <= 0.) then
        write(6,"(/,'>>> INPUT: ctpoten must be positive: ctpoten=',
     |    e12.4)") ctpoten
        call shutdown('CTPOTEN')
      endif
!
! Validate hemispheric power constant, if provided:
      if (power /= spval .and. power <= 0.) then
        write(6,"(/,'>>> INPUT: power must be positive: power=',
     |    e12.4)") power
        call shutdown('POWER')
      endif
!
! Validate times and values in user-provided time series for
! crosstail potential and hemispheric power:
      call validate_indices_time(ctpoten_time,ntimes_ctpoten,
     |  'CTPOTEN_TIME')
      call validate_indices_time(power_time,ntimes_power,
     |  'POWER_TIME')
!
! gpi_names and gpi_vars must be in the same order:
      gpi_names(1) = "f107            "
      gpi_names(2) = "f107a           "
      gpi_names(3) = "power           "
      gpi_names(4) = "ctpoten         "
!
      gpi_vars(1) = f107
      gpi_vars(2) = f107a
      gpi_vars(3) = power
      gpi_vars(4) = ctpoten
      if (any(power_time /= spval)) gpi_vars(3) = 0.   ! dummy non-spval
      if (any(ctpoten_time /= spval)) gpi_vars(4) = 0. ! dummy non-spval
!
! User did NOT provide one or more gpi_vars, meaning user wants
!   to use GPI database for those missing indices. In this case, 
!   if the user provided gpi_ncfile, use it, otherwise use the 
!   default file gpi_ncfiledef.
! (If the user DID provide some of these indices, then those values
!  will be used for those indices)
!
      if (any(gpi_vars==spval)) then
        write(6,"(/,'Input: The following GPI inputs were NOT ',
     |    'provided by the user:')")
        do i=1,ngpivars
          if (gpi_vars(i)==spval) write(6,"(4x,a)",advance='no')
     |      trim(gpi_names(i))          
        enddo
        write(6,"(/,4x,'The GPI database will be used to provide ',
     |    'values for these indices.')")
        if (len_trim(gpi_ncfile)==0) then
          write(6,"(4x,'Since gpi data file GPI_NCFILE was also NOT',
     |      ' provided,',/,4x,'I will use the default GPI_NCFILE: ',
     |      a)") trim(gpi_ncfiledef)
          gpi_ncfile = gpi_ncfiledef
        else
          write(6,"(4x,'Will use user-provided path to GPI data',
     |      ' GPI_NCFILE = ',a)") trim(gpi_ncfile)
        endif
!
! Model must advance in calendar time for GPI database to be used:
! (see calday verification above)
!
!       if (calday==0) then
        if (calendar_advance <= 0) then
          write(6,"(/,'>>> INPUT: Model must advance in calendar',
     |      ' time if GPI database is to be used.')")
          write(6,"('  To make a GPI run, you must set ',
     |      'CALENDAR_ADVANCE = 1')")
          call shutdown('GPI')
        endif
!
! User provided *all* gpi_vars:
! In this case, use the provided values. If user the also provided a
!   gpi_ncfile, then stop with error message.
!
      else
        write(6,"(/,'Input: User provided all gpi input ',
     |    'variables:')")
        do i=1,ngpivars
          if ((trim(gpi_names(i))=="ctpoten".and.ntimes_ctpoten > 0).or.
     |      (trim(gpi_names(i))=="ctpoten".and.ntimes_ctpoten > 0)) then
            write(6,"('  ',a,' = (time series)')") gpi_names(i)
          else
            write(6,"('  ',a,' = ',f10.2)") gpi_names(i),gpi_vars(i)
          endif
        enddo
        if (len_trim(gpi_ncfile) > 0) then
          write(6,"(/,'>>> INPUT: User provided all gpi input but',
     |      ' ALSO provided gpi data file path:',/,
     |      '  GPI_NCFILE = ',a)") trim(gpi_ncfile)
          write(6,"('  If you want to use the provided gpi values ',
     |      'please do NOT provide GPI_NCFILE.')")
          write(6,"('  If you want to make a GPI run (i.e., use GPI',
     |      ' data file), please comment out',/,4x,'(i.e., do NOT ',
     |      'provide) those gpi vars you want to get from ',/,4x,
     |      'the GPI database.',/)")
          call shutdown('GPI_NCFILE')
        endif
      endif
!
! Check to see if have hrindices_ncfile (when avail! 11/15/02)
!
! Check electric potential model:
      if (potential_model == 'WEIMER') then

! 5/7/04 btf: Weimer is not working. Setmodel01 in weimer01.F produces
!             NaNs for Coef(). (coefficients read from netcdf file appear
!             to be ok)
!
        write(6,"('>>> Potential_model WEIMER is not available -- ',
     |    'sub setmodel01 in weimer01.F produces NaNs for coeffs.')")
        call shutdown('weimer')
!
!       write (6,"(4x,'Will use the Weimer 2001 electric potential',
!    |    ' model')")
!       if (len_trim(weimer_ncfile)==0) then
!         write(6,"(4x,'Since Weimer coef data file WEIMER_NCFILE was ',
!    |      'NOT provided,',/,4x,'I will use the default ',
!    |      'WEIMER_NCFILE: ',a)") trim(weimer_ncfiledef)
!         weimer_ncfile = weimer_ncfiledef
!       else
!         write(6,"(4x,'Will use user-provided path to Weimer coefs',
!    |      ' Weimer_NCFILE = ',a)") trim(weimer_ncfile)
!       endif

      elseif (potential_model == 'HEELIS') then
        write (6,"(4x,'Will use the Heelis param electric potential',
     |    ' model')")
      elseif (potential_model == 'NONE') then
        write (6,"(4x,'Will use NONE (ie zero) electric potential',
     |    ' model')")
      endif
      if (len_trim(potential_model)==0) then
        write (6,"(4x,'Will use default Heelis elecric potential',
     |    ' model')")
        potential_model='HEELIS'
      endif
!
! BY component of solar IMF magnetic field:
      if (byimf==spval) then
        write(6,"(/,'>>> INPUT: need 1 value for byimf',
     |    ' e.g.: BYIMF=0.')")
        write(6,"('NOTE: multiple (time-dependent) values for CTPOTEN',
     |    ' not yet supported')")
        call shutdown('BYIMF')
      endif
!
! Limits on By NOT necessary for Weimer (or for NONE)
      if (potential_model == 'HEELIS') then
        if (byimf > 7.) then
          write(6,"(/,'>>> INPUT: byimf too big (must be <= 7.):',
     |      ' byimf=',e12.4)") byimf
          call shutdown('BYIMF')
        endif
        if (byimf < -11.) then
          write(6,"(/,'>>> INPUT: byimf too small (must be >= -11.):',
     |      ' byimf=',e12.4)") byimf
          call shutdown('BYIMF')
        endif
      endif
!
! Bz component of solar IMF magnetic field:
      if (bzimf==spval) then
        write(6,"('Input: Using default BZIMF=0.')")
        bzimf = 0.
      endif
!
! Solar wind velocity:
      if (swvel==spval) then
        write(6,"('Input: Using default SWVEL=400 km/s')")
        swvel = 400.
      endif
!
! Solar wind density:
      if (swden==spval) then
        write(6,"('Input: Using default SWDEN=4 #/cm3')")
        swden = 4.
      endif
!
! AL, lower auroral magnetic activity index:
      if (AL==spval) then
        write(6,"('Input: Using default AL=-20., ALUSE=FALSE')")
        al = -20.
        aluse = .false.
      else
        write(6,"('Input: AL provided, setting ALUSE=TRUE')")
        write(6,"('       AL = ',e12.4)") al
        aluse = .true.
      endif
!
! Collision factor:
      if (colfac==spval) then
        write(6,"('Input: Using default colfac = 1.5')")
        colfac = 1.5
      endif
!
! iuivi and dynamo flags:
! iuivi is no longer used -- is replaced by dynamo flag:
      if (iuivi /= ispval) then
        write(6,"(/,'Input: IUIVI flag is no longer used.')") 
        write(6,"('Please use DYNAMO flag to turn dynamo on or off.')")
        write(6,"('If DYNAMO=0, then dynamo is NOT called, and ',
     |    'electric potential and ion drifts will be zero.')")
        write(6,"('If DYNAMO=1, then dynamo IS called, and ',
     |    'electric potential and ion drifts are calculated.')")
        write(6,"('Default is DYNAMO=1')")
        call shutdown('IUIVI')
      endif
!
! 1/12/05 btf: removed old dynamo (was in dynamo_old.F)
! dynamo <= 0 -> no dynamo (dynamo routines are not called, poten==0)
! dynamo >  0 -> "new" dynamo (dynamo module in dynamo.F)
!
      if (dynamo==ispval) then
        dynamo = 1
      else
        if (dynamo <= 0) then
          write(6,"('Input: dynamo=',i3,' --> dynamo will NOT be ',
     |      'called.')") dynamo
        else
          write(6,"('Input: dynamo=',i3,' --> new dynamo will be ',
     |      'called.')") dynamo
        endif
      endif
      end subroutine inp_model
!-----------------------------------------------------------------------
      subroutine mkntask(ntask,ntask_lat,ntask_lon,ier)
      integer,intent(in) :: ntask
      integer,intent(out) :: ntask_lat,ntask_lon,ier
      integer :: i,j,ntlat(nlonp4),ntlon(nlonp4),nchoice,ngap,
     |  ngap_prev,ii
      ntask_lat = 0
      ntask_lon = 0
      nchoice=0
      do i=1,nlonp4
        do j=1,nlat
          if (i*j==ntask) then
            nchoice = nchoice+1
            write(6,"('mkntask: i=',i2,' j=',i2,' i*j=',i5,' ntask=',
     |        i5)") i,j,i*j,ntask
            ntlat(nchoice)=j 
            ntlon(nchoice)=i 
          endif
        enddo
      enddo
!
! Choose combinations of ntlat*ntlon==ntask in which ntlon==ntlat,
!   or if that does not exist, use combination with smallest delta
!   ntlon-ntlat.
!
      do i=1,nchoice
        if (ntlon(i) == ntlat(i)) then
          ntask_lat = ntlat(i)
          ntask_lon = ntlon(i)
        endif
      enddo
      if (ntask_lat==0.or.ntask_lon==0) then
        ngap_prev = nlonp4*nlat
        do i=1,nchoice
          ngap = ntlon(i)-ntlat(i)
          if (abs(ngap) < ngap_prev) then
            ngap = ntlon(i)-ntlat(i)
            ngap_prev = abs(ngap)
            ii = i
          endif
        enddo
        ntask_lat = ntlat(ii)
        ntask_lon = ntlon(ii)
!
! If they are not equal, ntlon should be > ntlat because nlon > nlat.
        if (ntask_lon < ntask_lat) then
          i = ntask_lat
          ntask_lat = ntask_lon
          ntask_lon = i
        endif
      endif
      ier =0
      if (ntask_lat==0.or.ntask_lon==0.or.ntask_lat*ntask_lon /= ntask)
     |   ier = 1
      end subroutine mkntask
!-----------------------------------------------------------------------
      subroutine inp_primhist
!
! Validate primary history inputs:
!
! Local:
      integer :: n,i,ii,nhists,nsaves,nsteps,nsteps_hist(mxseries),
     |  nsteps_save(mxseries),nstarts,nstops,nstep_total,nsrc,
     |  nout,nhists_total,modeltime(4),nfiles_prim
      integer(kind=8) :: step8,
     |  nsec_start(mxseries),nsec_stop(mxseries),
     |  nsec_hist (mxseries),nsec_save(mxseries),
     |  nsec_index,nsec0
      character(len=80) :: ch80
      character(len=80) :: hvols(mxhvols)
!
! External:
      integer,external :: numfiles
      integer(kind=8),external :: mtime_to_nsec
!
! 8-byte integer step:
      step8 = step
!
! Dispose flag:
      if (dispose==ispval) then
        write(6,"(/,'Input: dispose flag not read --',
     |    ' will default to dispose=1.',/)")
        dispose = 1
      elseif (dispose /= 0 .and. dispose /= 1 .and. dispose /= 2) then
        write(6,"(/,'>>> INPUT: DISPOSE flag must be 0, 1, or 2:',
     |    ' dispose = ',i4)") dispose
        write(6,"('DISPOSE = 0 -> do not dispose to mss.')")
        write(6,"('DISPOSE = 1 -> dispose to mss during model',
     |    ' execution.')")
        write(6,"('DISPOSE = 2 -> dispose to mss after model',
     |    ' execution.')")
        call shutdown('DISPOSE')
      endif
!
! Model start time(s):
      n = size(start)-count(start==ispval)
      if (mod(n,3) /= 0) then
        write(6,"(/,'>>> INPUT: START must be given as series of ',
     |    '3-integer triplet times',/,11x,'(day,hr,min).',/)")
        call shutdown('START')
      endif 
      nstarts = n/3 ! number of start times given
      if (nstarts < 1 .or. nstarts > mxseries) then
        write(6,"(/,'>>> INPUT: At least one and a maximum of ',i3,
     |    ' 3-integer START times are allowed.',/)") mxseries
        call shutdown('START')
      endif
      do i=1,nstarts
        call validate_mtime(start(:,i),mxday,'START')
      enddo
!
! Start time(s) must be multiple of step, and must increase: 
      nsec_start(:) = 0
      do i=1,nstarts
        nsec_start(i) = mtime_to_nsec(start(:,i))
        if (mod(nsec_start(i),step8) /= 0) then
          write(6,"(/,'>>> INPUT: START time ',i1,' must be a ',
     |      'multiple of step:',/,11x,'START=',3i4,' STEP=',i4,/)")
     |      i,start(:,i),step
          call shutdown('START')
        endif
        if (i > 1) then
          if (nsec_start(i-1) > nsec_start(i)) then
            write(6,"(/,'>>> INPUT: START times must increase.',/,
     |        11x,'START ',i2,' = ',3i4,' START ',i2,' = ',3i4,
     |        /)") i-1,start(:,i-1),i,start(:,i)
            call shutdown('START')
          endif
        endif
      enddo 
!
! If advancing in calendar time (calday > 0), then the starting model 
!   day must be the same as the starting calendar day (calday):
!
!     if (calday /= 0) then
      if (calendar_advance > 0) then
        if (start(1,1) /= start_day) then
          write(6,"(/,'>>> INPUT: Starting calendar day START_DAY =',
     |      i4)") start_day
          write(6,"(11x,'Starting model day START(1)=',i4)")
     |      start(1,1)
          write(6,"(11x,'CALENDAR_ADVANCE = ',i2)") calendar_advance
          write(6,"('If the model is to be advanced in calendar time,',
     |      ' the starting model day',/,'  must be equal to the ',
     |      'starting calendar day.',/)")
          call shutdown('START')
        endif
      endif
!
! Stop time(s):
      n = size(stop)-count(stop==ispval)
      if (mod(n,3) /= 0) then
        write(6,"(/,'>>> INPUT: STOP must be given as series of ',
     |    '3-integer triplet times (day,hr,min).',/)")
        call shutdown('STOP')
      endif 
      nstops = n/3 ! number of stop times given
      if (nstops < 1 .or. nstops > mxseries) then
        write(6,"(/,'>>> INPUT: At least one and a maximum of ',i2,
     |    ' 3-integer STOP times are allowed.',/)") mxseries
        call shutdown('STOP')
      endif
      do i=1,nstops
        call validate_mtime(stop(:,i),mxday,'STOP')
      enddo
!
! Stop time(s) must be multiple of step, and must increase: 
      nsec_stop(:) = 0
      do i=1,nstops
        nsec_stop(i) = mtime_to_nsec(stop(:,i))
        if (mod(nsec_stop(i),step8) /= 0) then
          write(6,"(/,'>>> INPUT: STOP time ',i1,' must be a ',
     |      'multiple of step:',/,11x,'STOP=',3i4,' STEP=',i4,/)")
     |      i,stop(:,i),step
          call shutdown('STOP')
        endif
        if (i > 1) then
          if (nsec_stop(i-1) > nsec_stop(i)) then
            write(6,"(/,'>>> INPUT: STOP times must increase.',/,
     |        11x,'STOP 1 =',3i4,' STOP 2 =',3i4,/)") stop
            call shutdown('STOP')
          endif
        endif
      enddo 
!
! Stop time(s) must be > start times:
! Note: this module does not cross year boundaries in a single run.
!       To cross a year boundary, make a run up to day 365 or 366,
!       then use that file as the source (source_start=365 or 366), 
!       for a new run with START=1,0,0. (note mxday==366)
!
      if (nsec_start(1)==nsec_stop(1)) then
        write(6,"(/,'NOTE input: start(1)==stop(1): no time steps',
     |    ' will be taken this run.')")
      else
        do i=1,nstops
          if (nsec_start(i) >= nsec_stop(i)) then
            write(6,"(/,'>>> INPUT: STOP time ',i1,' must be > ',
     |        'START time',/,11x,'STOP time ',i1,' = ',3i4,
     |        ' START =',3i4,/)") i,i,stop(:,i),start(:,i)
            call shutdown('START/STOP')
          endif
        enddo
      endif
!
! History write frequencies:
      n = size(hist)-count(hist==ispval)
      if (mod(n,3) /= 0) then
        write(6,"(/,'>>> INPUT: HIST must be given as series of ',
     |    '3-integer triplet times (day,hr,min).',/)")
        call shutdown('HIST')
      endif 
      nhists = n/3 ! number of hist times given
      if (nhists < 1 .or. nhists > mxseries) then
        write(6,"(/,'>>> INPUT: At least one and a maximum of ',i2,
     |    ' 3-integer HIST times are allowed.',/)") mxseries
        call shutdown('HIST')
      endif
      do i=1,nhists
        call validate_mtime(hist(:,i),mxday,'HIST')
      enddo
!
! History write frequencies must be multiple of step: 
      nsec_hist(:) = 0
      do i=1,nhists
        nsec_hist(i) = mtime_to_nsec(hist(:,i))
        if (nsec_hist(i)==0) then
          write(6,"(/,'>>> INPUT: HIST write frequency ',i1,
     |      ' must be > 0',/)") i
          call shutdown('HIST')
        endif
        if (mod(nsec_hist(i),step8) /= 0) then
          write(6,"(/,'>>> INPUT: HIST time ',i1,' must be a ',
     |      'multiple of step:',/,11x,'HIST=',3i4,' STEP=',i4,/)")
     |      i,hist(:,i),step
          call shutdown('HIST')
        endif
      enddo 
!
! History save frequencies:
      n = size(save)-count(save==ispval)
      if (mod(n,3) /= 0) then
        write(6,"(/,'>>> INPUT: SAVE must be given as series of ',
     |    '3-integer triplet times (day,hr,min).',/)")
        call shutdown('SAVE')
      endif 
      nsaves = n/3 ! number of save times given
      if (nsaves < 1 .or. nsaves > mxseries) then
        write(6,"(/,'>>> INPUT: At least one and a maximum of ',i2,
     |    ' 3-integer SAVE times are allowed.',/)") mxseries
        call shutdown('SAVE')
      endif
      do i=1,nsaves
        call validate_mtime(save(:,i),mxday,'SAVE')
      enddo
!
! History save frequencies must be multiple of step: 
      nsec_save(:) = 0
      do i=1,nsaves
        nsec_save(i) = mtime_to_nsec(save(:,i))
        if (nsec_save(i)==0) then
          write(6,"(/,'>>> INPUT: SAVE history save frequency ',i1,
     |      ' must be > 0',/)") i
          call shutdown('SAVE')
        endif
        if (mod(nsec_save(i),step8) /= 0) then
          write(6,"(/,'>>> INPUT: SAVE time ',i1,' must be a ',
     |      'multiple of step:',/,11x,'SAVE=',3i4,' STEP=',i4,/)")
     |      i,save(:,i),step
          call shutdown('SAVE')
        endif
      enddo 
!
! Must have same number of time sequences:
      if (nstarts /= nstops .or. nstarts /= nhists .or. 
     |    nstarts /= nsaves .or. nstops  /= nhists .or.
     |    nstops  /= nsaves .or. nhists  /= nsaves) then
        write(6,"(/,'>>> INPUT: must provide same number of times ',
     |    'for',/,11x,'START, STOP, HIST, and SAVE.')")
        write(6,"(11x,'nstarts=',i3,' nstops=',i3,' nhists=',i3,
     |    ' nsaves=',i3,/)") nstarts,nstops,nhists,nsaves
        call shutdown('ntimes')
      endif
!
! SAVEs must be multiples of HISTs:
      do i=1,nstarts
        if (mod(nsec_save(i),nsec_hist(i)) /= 0) then
          write(6,"(/,'>>> INPUT: SAVE frequencies must be multiple',
     |      ' of HIST frequencies.',/,11x,'SAVE ',i1,' =',3i4,
     |      ' HIST ',i1,' =',3i4,/)") i,save(:,i),i,hist(:,i)
          call shutdown('mod(SAVE,HIST)')
        endif
      enddo
!
! Number of steps in each time series must be a multiple of HIST
! and SAVE frequencies for that series:
      do i=1,nstarts
        if (i==1) then
          nsteps = (nsec_stop(i)-nsec_start(i))/step
        else
          nsteps = (nsec_stop(i)-nsec_stop(1))/step
        endif 
        nsteps_hist(i) = nsec_hist(i)/step
        if (mod(nsteps,nsteps_hist(i)) /= 0) then
          write(6,"(/,'>>> INPUT: number of steps in time series ',
     |      i1,' must be multiple of the ',/,11x,'number of steps',
     |      ' in HIST ',i1,'.')") i,i
          write(6,"(11x,'nsteps ',i1,' = ',i6,' nsteps_hist ',i1,
     |      ' = ',i3)") i,nsteps,i,nsteps_hist(i)
          write(6,"(11x,'START',i1,' = ',3i4,' STOP',i1,' = ',3i4,
     |      ' HIST',i1,' = ',3i4,/)") i,start(:,i),i,stop(:,i),
     |      i,hist(:,i)
          call shutdown('HIST')
        endif
        nsteps_save(i) = nsec_save(i)/step
        if (mod(nsteps,nsteps_save(i)) /= 0) then
          write(6,"(/,'>>> INPUT: number of steps in time series ',
     |      i1,' must be multiple of the ',/,11x,'number of steps',
     |      ' in SAVE ',i1,'.')") i,i
          write(6,"(11x,'nsteps ',i1,' = ',i6,' nsteps_save ',i1,
     |      ' = ',i3)") i,nsteps,i,nsteps_save(i)
          write(6,"(11x,'START',i1,' = ',3i4,' STOP',i1,' = ',3i4,
     |      ' SAVE',i1,' = ',3i4,/)") i,start(:,i),i,stop(:,i),
     |      i,save(:,i)
          call shutdown('SAVE')
        endif
      enddo
!
! Time series cannot overlap (However, they can touch, front-to-back):
      if (nstarts > 1) then
        do i=2,nstarts
          if (nsec_start(i) < nsec_stop(i-1)) then
            write(6,"(/,'>>> INPUT: primary history time series',
     |        ' cannot overlap.')")
            write(6,"(11x,'For series ',i2,': START=',3i4,
     |        ' STOP=',3i4)") i-1,start(:,i-1),start(:,i-1)
            write(6,"(11x,'For series ',i2,': START=',3i4,
     |        ' STOP=',3i4,/)") i,start(:,i),stop(:,i)
            call shutdown('START/STOP')
          endif
        enddo
      endif
!
! Total steps this run (nstep is in hist_module.F)
      nstep_total = (nsec_stop(nstarts)-nsec_start(1))/step8
!
! Source history file (optional):
      nsrc = len_trim(source)
!
! Source start time (must be given if SOURCE file was provided):
      n = size(source_start)-count(source_start==ispval)
      if (nsrc > 0 .and. n <= 0) then
        write(6,"(/,'>>> INPUT: If SOURCE is provided, must also',
     |    ' provide SOURCE_START time.',/,11x,'SOURCE=',a,/)")
     |    trim(source)
        call shutdown('SOURCE_START')
      endif
      if (nsrc > 0) then
        if (n /= 3) then
          write(6,"(/,'>>> INPUT: need 3 values for SOURCE_START ',
     |      'time (day,hour,minute),',/,11x,
     |      'e.g.: SOURCE_START=1,0,0',/)")
          call shutdown('START')
        endif
        call validate_mtime(source_start,mxday,'SOURCE_START')
!
! Model start time hr,min must be same as source_start hr,min:
! (days can be different)
        if (start(2,1) /= source_start(2) .or.
     |      start(3,1) /= source_start(3)) then
          write(6,"(/,'>>> INPUT: START time (hr,min) ',
     |      'must be the same as SOURCE_START time.')")
          write(6,"('  START        = ',3i4)") start(:,1)
          write(6,"('  SOURCE_START = ',3i4)") source_start(:)
          write(6,"('  (START and SOURCE_START days can be ',
     |      'different)',/)")
          call shutdown('START')
        endif
      endif
!
! Validate that times given by ctpoten_time and/or power_time are
! within start and stop times (the times themselves were validated 
! in inp_model). Note user must provide a value for each start time.
!
      if (ntimes_ctpoten > 0) then
!
! Must provide first ctpoten_time(:,1) same as START time.
        nsec0 = mtime_to_nsec(int(ctpoten_time(1:3,1)))
        if (nsec0 /= nsec_start(1)) then
          write(6,"(/,'>>> INPUT: Please provide first CTPOTEN_TIME',
     |      ' time equal to first START time.')")
          call shutdown('starting ctpoten_time')
        endif
!
! All ctpoten times must be within start and stop times:
        do i=1,ntimes_ctpoten
          nsec_index = mtime_to_nsec(int(ctpoten_time(1:3,i)))
          if (nsec_index < nsec_start(1) .or. 
     |        nsec_index > nsec_stop(1)) then
            write(6,"(/,'>>> INPUT: CTPOTEN_TIME ',3i4,' is outside ',
     |        'model START/STOP times: START=',3i4,' STOP=',3i4)")
     |        int(ctpoten_time(1:3,i)),start(:,1),stop(:,1)
            call shutdown('ctpoten_time')
          endif 
        enddo
!
! Init:
        ctpoten = ctpoten_time(4,1)
      endif ! ntimes_ctpoten
      if (ntimes_power > 0) then
!
! Must provide first power_time(:,1) same as START time.
        nsec0 = mtime_to_nsec(int(power_time(1:3,1)))
        if (nsec0 /= nsec_start(1)) then
          write(6,"(/,'>>> INPUT: Please provide first POWER_TIME',
     |      ' time equal to first START time.')")
          call shutdown('starting power_time')
        endif
!
! All power times must be within start and stop times:
        do i=1,ntimes_power
          nsec_index = mtime_to_nsec(int(power_time(1:3,i)))
          if (nsec_index < nsec_start(1) .or. 
     |        nsec_index > nsec_stop(1)) then
            write(6,"(/,'>>> INPUT: POWER_TIME ',3i4,' is outside ',
     |        'model START/STOP times: START=',3i4,' STOP=',3i4)")
     |        int(power_time(1:3,i)),start(:,1),stop(:,1)
            call shutdown('power_time')
          endif 
        enddo
!
! Init:
        power = power_time(4,1)
      endif ! ntimes_power
!
! Primary output volumes:
! (integer function mkhvols either echoes histvols to hvols, or if 
!  histvols(2)=='to',then it expands histvols from 'volfirst','to',
!  'vollast','by','n' to hvols)
!
      hvols = ' '
      nout = mkhvols(output,hvols,mxhvols)
      if (nout==0) then
        write(6,"(/,'>>> INPUT: need at least one output volume',/)") 
        call shutdown('OUTPUT')
      endif
      output = hvols
!
! Max number of histories per primary file:
! (mxhist_prim is an input parameter with default = 13).
      if (mxhist_prim < 0) then
        write(6,"('>>> INPUT: maximum number of histories per ',
     |    'primary file must be > 0: mxhist_prim=',i4)") mxhist_prim
        call shutdown('MXHIST_PRIM')
      endif
!
! mss retention period:
      if (msreten < 1 .or. msreten > 32767) then
        write(6,"('>>> INPUT: bad MSRETEN = ',i6)") msreten
        write(6,"('Mss retention period MSRETEN must be >= 1',
     |    ' and <= 32767')")
        call shutdown('MSRETEN')
      endif
!
! No dups of file names allowed:
      ch80 = ' '
      do i=1,nout
        ch80 = output(i)  
        output(i) = 'dummy'
        if (any(output==ch80)) then
          write(6,"(/,'>>> INPUT: Duplicate OUTPUT file names = ',
     |      a,/)") trim(ch80)
          call shutdown('OUTPUT')
        endif
        output(i) = ch80
      enddo
!
! Check that sufficient primary output files have been provided:
! Func numfiles returns number of files that will be needed, and
!   also returns total histories to be written:
!
      nfiles_prim = numfiles('prim',nstarts,nsrc,0,nhists_total)
      if (nout < nfiles_prim) then
        write(6,"(/,'>>> INPUT: Will need ',i3,' OUTPUT files, but',
     |    ' read only ',i3)") nfiles_prim,nout
        write(6,"(11x,'Total number of steps this run = ',i6)") 
     |    nstep_total
        write(6,"(11x,'Total number of primary histories this run = ',
     |    i5)") nhists_total
        write(6,"(11x,'Maximum number of primary histories per file = ',
     |    i3,/)") mxhist_prim
        call shutdown('OUTPUT')
      endif
      end subroutine inp_primhist
!-------------------------------------------------------------------
      subroutine inp_sechist
!
! Validate secondary history inputs:
!
! Local:
      integer :: n,i,ii,nstarts,nstops,nhists,nsaves,nout,nsteps,
     |  nsteps_hist(mxseries_sech),nsteps_save(mxseries_sech),
     |  nsechs_total,nseriesp,nflds_sech,nhists_total,
     |  nfiles_sech,nonblank
      integer(kind=8) :: step8,
     |  nsec_start(mxseries_sech),nsec_stop(mxseries_sech),
     |  nsec_hist (mxseries_sech),nsec_save(mxseries_sech)
      character(len=80) :: ch80
      character(len=80) :: hvols(mxhvols)
      character(len=16) ::
     |  secflds_tmp(mxfsech)
      logical :: found
!
! External:
      integer,external :: numfiles
      integer(kind=8),external :: mtime_to_nsec
!
! 8-byte integer step:
      step8 = step
!
! n = total number of secondary history inputs read:
      n = size(secstart)-count(secstart==ispval) +
     |    size(secstop) -count(secstop==ispval)  +
     |    size(sechist) -count(sechist==ispval)  +
     |    size(secsave) -count(secsave==ispval)  +
     |    size(secflds) -count(len_trim(secflds)==0) +
     |    size(secfmag) -count(len_trim(secfmag)==0) +
     |    size(secfgeo2d) -count(len_trim(secfgeo2d)==0)+
     |    size(secfmag2d) -count(len_trim(secfmag2d)==0)+
     |    size(secfmagphr)-count(len_trim(secfmagphr)==0)
!
! Secondary output volumes:
! (integer function mkhvols either echoes histvols to hvols, or if 
!  histvols(2)=='to',then it expands histvols from 'volfirst','to',
!  'vollast','by','n' to hvols)
!
      hvols = ' '
      nout = mkhvols(secout,hvols,mxhvols)
      if (nout==0.and.n > 0) then
        write(6,"(/,'>>> INPUT: need at least one secondary ',
     |    'history output volume',/)") 
        call shutdown('SECOUT')
      endif
      secout = hvols
      n = n+nout
      if (n <= 0) return
!
! Secondary history start time(s):
      n = size(secstart)-count(secstart==ispval)
      if (mod(n,3) /= 0) then
        write(6,"(/,'>>> INPUT: SECSTART must be given as series of ',
     |    '3-integer triplet times',/,11x,'(day,hr,min).',/,11x,
     |    'A maximum of ',i3,' secstart times are allowed.',/)")
     |    mxseries_sech
        call shutdown('SECSTART')
      endif 
      nstarts = n/3 ! number of start times given
      if (nstarts < 1 .or. nstarts > mxseries_sech) then
        write(6,"(/,'>>> INPUT: At least one and a maximum of ',i3,
     |    ' 3-integer SECSTART times are allowed.',/)") mxseries_sech
        call shutdown('SECSTART')
      endif
      do i=1,nstarts
        call validate_mtime(secstart(:,i),mxday,'SECSTART')
      enddo
!
! Secondary start time(s) must be multiple of step, and must increase: 
      nsec_start(:) = 0
      do i=1,nstarts
        nsec_start(i) = mtime_to_nsec(secstart(:,i))
        if (mod(nsec_start(i),step8) /= 0) then
          write(6,"(/,'>>> INPUT: SECSTART time ',i1,' must be a ',
     |      'multiple of step:',/,11x,'SECSTART=',3i4,' STEP=',i4,/)")
     |      i,secstart(:,i),step
          call shutdown('SECSTART')
        endif
        if (i > 1) then
          if (nsec_start(i-1) > nsec_start(i)) then
            write(6,"(/,'>>> INPUT: SECSTART times must increase.',/,
     |        11x,'SECSTART ',i2,' = ',3i4,' SECSTART ',i2,' = ',3i4,
     |        /)") i-1,secstart(:,i-1),i,secstart(:,i)
            call shutdown('SECSTART')
          endif
        endif
      enddo 
!
! Secondary start times must be >= first primary start times and
! <= last primary stop time:
      nseriesp = (size(start)-count(start==ispval))/3
      do i=1,nstarts
        if (nsec_start(i) < mtime_to_nsec(start(:,1))) then
          write(6,"(/,'>>> INPUT: all secondary start times SECSTART',
     |      ' must be >= first model START time.')")
          write(6,"(11x,'First model START = ',3i4,' SECSTART ',i2,
     |      ' = ',3i4,/)") start(:,1),i,secstart(:,i)
          call shutdown('SECSTART')
        endif
        if (nsec_start(i) > mtime_to_nsec(stop(:,nseriesp))) then
          write(6,"(/,'>>> INPUT: all secondary start times SECSTART',
     |      ' must be <= last model STOP time.')")
          write(6,"(11x,' SECSTART ',i2,' = ',3i4,' last STOP=',
     |      3i4)") secstart(:,1),i,stop(:,nseriesp)
          call shutdown('SECSTART')
        endif
      enddo
!
! Secondary history stop time(s):
      n = size(secstop)-count(secstop==ispval)
      if (mod(n,3) /= 0) then
        write(6,"(/,'>>> INPUT: SECSTOP must be given as series of ',
     |    '3-integer triplet times',/,11x,'(day,hr,min).',/,11x,
     |    'A maximum of ',i3,' secstop times are allowed.',/)")
     |    mxseries_sech
        call shutdown('SECSTOP')
      endif 
      nstops = n/3 ! number of stop times given
      if (nstops < 1 .or. nstops > mxseries_sech) then
        write(6,"(/,'>>> INPUT: At least one and a maximum of ',i3,
     |    ' 3-integer SECSTOP times are allowed.',/)") mxseries_sech
        call shutdown('SECSTOP')
      endif
      do i=1,nstops
        call validate_mtime(secstop(:,i),mxday,'SECSTOP')
      enddo
!
! Stop time(s) must be multiple of step, and must increase: 
      nsec_stop(:) = 0
      do i=1,nstops
        nsec_stop(i) = mtime_to_nsec(secstop(:,i))
        if (mod(nsec_stop(i),step8) /= 0) then
          write(6,"(/,'>>> INPUT: SECSTOP time ',i1,' must be a ',
     |      'multiple of step:',/,11x,'SECSTOP=',3i4,' STEP=',i4,/)")
     |      i,secstop(:,i),step
          call shutdown('SECSTOP')
        endif
        if (i > 1) then
          if (nsec_stop(i-1) > nsec_stop(i)) then
            write(6,"(/,'>>> INPUT: SECSTOP times must increase.',/,
     |        11x,'SECSTOP ',i2,' = ',3i4,' SECSTOP ',i2,' = ',3i4,
     |        /)") i-1,secstop(:,i-1),i,secstop(:,i)
            call shutdown('SECSTOP')
          endif
        endif
      enddo 
!
! Secondary stop times must be > secondary start times:
      do i=1,nstops
        if (nsec_stop(i) <= nsec_start(i)) then
          write(6,"(/,'>>> INPUT: SECSTART must be < SECSTOP for ',
     |      'all time series.')")
          write(6,"('For time series ',i2,': SECSTART=',3i4,' SECSTOP=',
     |      3i4,/)") i,secstart(:,i),secstop(:,i)
          call shutdown('SECSTART/SECSTOP')
        endif
      enddo
!
! Secondary stop times must be <= last primary stop time:
      nstarts = n/3 ! number of start times given
      do i=1,nstarts
        if (nsec_stop(i) > mtime_to_nsec(stop(:,nseriesp))) then
          write(6,"(/,'>>> INPUT: all secondary stop times must be',
     |      ' <= final model stop time.')")
          write(6,"('For sech time series ',i2,': SECSTOP=',3i4,
     |      ' Model STOP = ',3i4,/)") i,secstop(:,i),stop(:,nseriesp)
          call shutdown('SECSTOP')
        endif
      enddo
!
! Secondary history write frequencies:
      n = size(sechist)-count(sechist==ispval)
      if (mod(n,3) /= 0) then
        write(6,"(/,'>>> INPUT: SECHIST must be given as series of ',
     |    '3-integer triplet times (day,hr,min).',/,11x,'A maximum',
     |    ' of ',i3,' SECHIST times are allowed.',/)") mxseries_sech
        call shutdown('SECHIST')
      endif 
      nhists = n/3 ! number of hist times given
      if (nhists < 1 .or. nhists > mxseries_sech) then
        write(6,"(/,'>>> INPUT: At least one and a maximum of ',i3,
     |    ' 3-integer SECHIST times are allowed.',/)")
        call shutdown('SECHIST')
      endif
      do i=1,nhists
        call validate_mtime(sechist(:,i),mxday,'SECHIST')
      enddo
!
! Secondary history write frequencies must be multiples of step: 
      nsec_hist(:) = 0
      do i=1,nhists
        nsec_hist(i) = mtime_to_nsec(sechist(:,i))
        if (nsec_hist(i)==0) then
          write(6,"(/,'>>> INPUT: SECHIST write frequency ',i1,
     |      ' must be > 0',/)") i
          call shutdown('SECHIST')
        endif
        if (mod(nsec_hist(i),step8) /= 0) then
          write(6,"(/,'>>> INPUT: SECHIST time ',i1,' must be a ',
     |      'multiple of step:',/,11x,'SECHIST=',3i4,' STEP=',i4,/)")
     |      i,sechist(:,i),step
          call shutdown('SECHIST')
        endif
      enddo 
!
! Secondary history save frequencies:
      n = size(secsave)-count(secsave==ispval)
      if (mod(n,3) /= 0) then
        write(6,"(/,'>>> INPUT: SECSAVE must be given as series of ',
     |    '3-integer triplet times (day,hr,min).',/,11x,'A maximum',
     |    ' of ',i3,' SECSAVE times are allowed.',/)") mxseries_sech
        call shutdown('SECSAVE')
      endif 
      nsaves = n/3 ! number of save times given
      if (nsaves < 1 .or. nsaves > mxseries_sech) then
        write(6,"(/,'>>> INPUT: At least one and a maximum of ',i3,
     |    ' 3-integer SECSAVE times are allowed.',/)")
        call shutdown('SECSAVE')
      endif
      do i=1,nsaves
        call validate_mtime(secsave(:,i),mxday,'SECSAVE')
      enddo
!
! Secondary history save frequencies must be multiples of step: 
      nsec_save(:) = 0
      do i=1,nsaves
        nsec_save(i) = mtime_to_nsec(secsave(:,i))
        if (nsec_save(i)==0) then
          write(6,"(/,'>>> INPUT: SECSAVE save frequency ',i1,
     |      ' must be > 0',/)") i
          call shutdown('SECSAVE')
        endif
        if (mod(nsec_save(i),step8) /= 0) then
          write(6,"(/,'>>> INPUT: SECSAVE frequency ',i1,' must be a ',
     |      'multiple of step:',/,11x,'SECSAVE=',3i4,' STEP=',i4,/)")
     |      i,sechist(:,i),step
          call shutdown('SECSAVE')
        endif
      enddo 
!
! Must have same number of time sequences:
      if (nstarts /= nstops .or. nstarts /= nhists .or. 
     |    nstarts /= nsaves .or. nstops  /= nhists .or.
     |    nstops  /= nsaves .or. nhists  /= nsaves) then
        write(6,"(/,'>>> INPUT: must provide same number of times ',
     |    'for',/,11x,'SECSTART, SECSTOP, SECHIST, and SECSAVE.')")
        write(6,"(11x,'nstarts=',i3,' nstops=',i3,' nhists=',i3,
     |    ' nsaves=',i3,/)") nstarts,nstops,nhists,nsaves
        call shutdown('nsechtimes')
      endif
!     nseries_sech = nstarts
!
! SAVEs must be multiples of HISTs:
      do i=1,nstarts
        if (mod(nsec_save(i),nsec_hist(i)) /= 0) then
          write(6,"(/,'>>> INPUT: SECSAVE frequencies must be ',
     |      'multiples of SECHIST frequencies.',/,11x,'SECSAVE ',
     |      i1,' =',3i4,' SECHIST ',i1,' =',3i4,/)") i,secsave(:,i),
     |      i,sechist(:,i)
          call shutdown('mod(SECSAVE,SECHIST)')
        endif
      enddo
!
! Number of steps in each time series must be a multiple of SECHIST
! and SECSAVE frequencies for that series:
      do i=1,nstarts
        nsteps = (nsec_stop(i)-nsec_start(i))/step
        nsteps_hist(i) = nsec_hist(i)/step
        if (mod(nsteps,nsteps_hist(i)) /= 0) then
          write(6,"(/,'>>> INPUT: number of steps in time series ',
     |      i1,' must be multiple of the ',/,11x,'number of steps',
     |      ' in SECHIST ',i1,'.')") i,i
          write(6,"(11x,'nsteps ',i1,' = ',i6,' nsteps_hist ',i1,
     |      ' = ',i3)") i,nsteps,i,nsteps_hist(i)
          write(6,"(11x,'Time series ',i2,': SECSTART = ',3i4,
     |      ' SECSTOP = ',3i4,' SECHIST = ',3i4,/)") i,secstart(:,i),
     |      secstop(:,i),sechist(:,i)
          call shutdown('SECHIST')
        endif
        nsteps_save(i) = nsec_save(i)/step
        if (mod(nsteps,nsteps_save(i)) /= 0) then
          write(6,"(/,'>>> INPUT: number of steps in time series ',
     |      i1,' must be multiple of the ',/,11x,'number of steps',
     |      ' in SECSAVE ',i1,'.')") i,i
          write(6,"(11x,'nsteps ',i1,' = ',i6,' nsteps_save ',i1,
     |      ' = ',i3)") i,nsteps,i,nsteps_save(i)
          write(6,"(11x,'Time series ',i2,': SECSTART = ',3i4,
     |      ' SECSTOP = ',3i4,' SECSAVE = ',3i4,/)") i,secstart(:,i),
     |      secstop(:,i),secsave(:,i)
          call shutdown('SECSAVE')
        endif
      enddo
!
! Time series cannot overlap (However, they can touch, front-to-back):
      if (nstarts > 1) then
        do i=2,nstarts
          if (nsec_start(i) < nsec_stop(i-1)) then
            write(6,"(/,'>>> INPUT: secondary history time series',
     |        ' cannot overlap.')")
            write(6,"(11x,'For series ',i2,': SECSTART=',3i4,
     |        ' SECSTOP=',3i4)") i-1,secstart(:,i-1),secstop(:,i-1)
            write(6,"(11x,'For series ',i2,': SECSTART=',3i4,
     |        ' SECSTOP=',3i4,/)") i,secstart(:,i),secstop(:,i)
            call shutdown('SECSTART/SECSTOP')
          endif
        enddo
      endif
!
! Max number of histories per secondary file (optional input):
      if (mxhist_sech < 0) then
        write(6,"('>>> INPUT: maximum number of histories per ',
     |    'secondary file must be > 0: mxhist_sech=',i4)") mxhist_sech
        call shutdown('MXHIST_SECH')
      endif
!
! No dups of secout file names allowed:
      ch80 = ' '
      do i=1,nout
        ch80 = secout(i)  
        secout(i) = 'dummy'
        if (any(secout==ch80)) then
          write(6,"(/,'>>> INPUT: Duplicate SECOUT file names = ',
     |      a,/)") trim(ch80)
          call shutdown('SECOUT')
        endif
        secout(i) = ch80
      enddo
!
! Check that sufficient secondary output files have been provided:
! Func numfiles returns number of files that will be needed, and
!   also returns total histories to be written:
!
      nfiles_sech = numfiles('sech',nstarts,1,0,nhists_total)
      if (nout < nfiles_sech) then
        write(6,"(/,'>>> INPUT: Will need ',i3,' SECOUT files, but',
     |    ' read only ',i3)") nfiles_sech,nout
        write(6,"(11x,'Total number of secondary histories this ',
     |    'run = ',i5)") nhists_total
        write(6,"(11x,'Maximum number of secondary histories per ',
     |    'file = ',i3,/)") mxhist_sech
        call shutdown('SECOUT')
      endif
!
! Pack names so no blank names occur from 1->nflds_sech
      nflds_sech = count(len_trim(secflds) > 0)
      call packstr(secflds,mxfsech,nonblank)
      if (nonblank /= nflds_sech) then
        write(6,"('>>> WARNING: Input after packstr(secflds): ',
     |    'nonblank /= nflds_sech: nonblank=',i3,' nflds_sech=',i3)")
     |    nonblank,nflds_sech
        write(6,"('secflds(mxfsech=',i3,')=')") mxfsech
        do i=1,mxfsech
          write(6,"('secflds(',i2,')=',a)") i,secflds(i)
        enddo
      endif
!
! Secondary history fields:
! Fields that are forced on the secondary histories are listed
!   in string array secflds_mandatory().
! If user does not provide secflds fields, these mandatory fields
!   are written by default. Any mandatory fields NOT listed by
!   the user are added to the end of the user's list.
!
      if (nflds_sech==0) then
        write(6,"(/,'INPUT NOTE: no secondary history fields were ',
     |    'requested (SECFLDS).')")
        write(6,"('I will write the default minimum set of ',
     |    'fields to secondary histories:')")
        write(6,"(4a12)") secflds_mandatory
        do i=1,size(secflds_mandatory)
          secflds(i) = secflds_mandatory(i)
        enddo
      else ! nflds_sech > 0: enforce mandatory fields
        do i=1,size(secflds_mandatory)
          if (.not.any(secflds==secflds_mandatory(i))) then
            nflds_sech = nflds_sech+1
            secflds(nflds_sech) = secflds_mandatory(i)
            write(6,"('INPUT NOTE: adding mandatory field ',
     |        a,' to secondary history fields (field ',i3,')')")
     |        secflds(nflds_sech)(1:8),nflds_sech
          endif
        enddo
      endif

! Check for dups of secflds field names:
      ch80 = ' '
      do i=1,nflds_sech
        ch80 = secflds(i)  
        secflds(i) = 'dummy'
        if (any(secflds==ch80)) then
          write(6,"(/,'>>> INPUT: Duplicate SECFLDS field names = ',
     |      a,/)") trim(ch80)
          call shutdown('SECFLDS')
        endif
        secflds(i) = ch80
      enddo
!
! Secondary history fields on the magnetic grid:
      nflds_sech = count(len_trim(secfmag) > 0)
!
! Pack names so no blank names occur from 1->nflds_sech
      call packstr(secfmag,mxfsech,nonblank)
      if (nonblank /= nflds_sech) then
        write(6,"('>>> WARNING: Input after packstr(secfmag): ',
     |    'nonblank /= nflds_sech: nonblank=',i3,' nflds_sech=',i3)")
     |    nonblank,nflds_sech
        write(6,"('secfmag(mxfsech=',i3,')=')") mxfsech
        do i=1,mxfsech
          write(6,"('secfmag(',i2,')=',a)") i,secfmag(i)
        enddo
      endif
!
! Check for dups of secfmag field names:
      if (nflds_sech > 0) then
        ch80 = ' '
        do i=1,nflds_sech
          ch80 = secfmag(i)
          secfmag(i) = 'dummy'
          if (any(secfmag==ch80)) then
            write(6,"(/,'>>> INPUT: Duplicate SECFMAG field names = ',
     |        a,/)") trim(ch80)
            call shutdown('SECFMAG')
          endif
          secfmag(i) = ch80
        enddo
      endif 
!
! Add mandatory 3d magnetic fields to secondary history:
      do i=1,size(secfmag_mandatory)
        write(6,"('i=',i3,' Mandatory secfmag field ',a)") 
     |    i,trim(secfmag_mandatory(i))
        found = .false.
        if (nflds_sech > 0) then
          do ii=1,nflds_sech
            if (trim(secfmag(ii))==trim(secfmag_mandatory(i))) then
              found = .true.
            endif
          enddo
        endif
        if (.not.found) then
          write(6,"('Input: adding mandatory mag field ',a,
     |      ' to secondary history.')")
          nflds_sech = nflds_sech+1
          secfmag(nflds_sech) = secfmag_mandatory(i)
        endif
      enddo
!
! Secondary history fields on the geographic 2d grid:
      nflds_sech = count(len_trim(secfgeo2d) > 0)
!
! Pack names so no blank names occur from 1->nflds_sech
      call packstr(secfgeo2d,mxfsech,nonblank)
      if (nonblank /= nflds_sech) then
        write(6,"('>>> WARNING: Input after packstr(secfgeo2d): ',
     |    'nonblank /= nflds_sech: nonblank=',i3,' nflds_sech=',i3)")
     |    nonblank,nflds_sech
        write(6,"('secgeo2d(mxfsech=',i3,')=')") mxfsech
        do i=1,mxfsech
          write(6,"('secfgeo2d(',i2,')=',a)") i,secfgeo2d(i)
        enddo
      endif
!
! Check for dups of secfgeo2d field names:
      if (nflds_sech > 0) then
        ch80 = ' '
        do i=1,nflds_sech
          ch80 = secfgeo2d(i)
          secfgeo2d(i) = 'dummy'
          if (any(secfgeo2d==ch80)) then
            write(6,"(/,'>>> INPUT: Duplicate SECFGEO2D field names = ',
     |        a,/)") trim(ch80)
            call shutdown('SECFGEO2D')
          endif
          secfgeo2d(i) = ch80
        enddo
      endif
!
! Secondary history fields on the magnetic 2d grid:
      nflds_sech = count(len_trim(secfmag2d) > 0)
!
! Pack names so no blank names occur from 1->nflds_sech
      call packstr(secfmag2d,mxfsech,nonblank)
      if (nonblank /= nflds_sech) then
        write(6,"('>>> WARNING: Input after packstr(secfmag2d): ',
     |    'nonblank /= nflds_sech: nonblank=',i3,' nflds_sech=',i3)")
     |    nonblank,nflds_sech
        write(6,"('secfmag2d(mxfsech=',i3,')=')") mxfsech
        do i=1,mxfsech
          write(6,"('secfmag2d(',i2,')=',a)") i,secfmag2d(i)
        enddo
      endif
!
! Check for dups of secfmag 2d field names:
      if (nflds_sech > 0) then
        ch80 = ' '
        do i=1,nflds_sech
          ch80 = secfmag2d(i)
          secfmag2d(i) = 'dummy'
          if (any(secfmag2d==ch80)) then
            write(6,"(/,'>>> INPUT: Duplicate SECFMAG2D field names = ',
     |        a,/)") trim(ch80)
            call shutdown('SECFMAG2D')
          endif
          secfmag2d(i) = ch80
        enddo
      endif
!
! Secondary history fields on the magnetospheric grid:
      nflds_sech = count(len_trim(secfmagphr) > 0)
!
! Pack names so no blank names occur from 1->nflds_sech
      call packstr(secfmagphr,mxfsech,nonblank)
      if (nonblank /= nflds_sech) then
        write(6,"('>>> WARNING: Input after packstr(secfmagphr): ',
     |    'nonblank /= nflds_sech: nonblank=',i3,' nflds_sech=',i3)")
     |    nonblank,nflds_sech
        write(6,"('secfmagphr(mxfsech=',i3,')=')") mxfsech
        do i=1,mxfsech
          write(6,"('secfmagphr(',i2,')=',a)") i,secfmagphr(i)
        enddo
      endif
!
! Check for dups of secfmagphr field names:
      if (nflds_sech > 0) then
        ch80 = ' '
        do i=1,nflds_sech
          ch80 = secfmagphr(i)
          secfmagphr(i) = 'dummy'
          if (any(secfmagphr==ch80)) then
            write(6,"(/,'>>> INPUT: Duplicate SECFMAGPHR field
     |        names',' = ', a,/)") trim(ch80)
            call shutdown('SECFMAGPHR')
          endif
          secfmagphr(i) = ch80
        enddo
      endif
      end subroutine inp_sechist
!-----------------------------------------------------------------------
      subroutine inp_deftyp
!
! Define inp (type(input_type)) from validated namelist input 
!   module variables:
!
      inp%label   = label
      inp%tempdir = tempdir
      inp%magvol  = magvol
      inp%amievol = amievol
      inp%potential_model = potential_model
      inp%hrindices_ncfile = hrindices_ncfile
      inp%weimer_ncfile = weimer_ncfile
      inp%gpi_ncfile = gpi_ncfile
      inp%gswm_di_ncfile    = gswm_di_ncfile
      inp%gswm_sdi_ncfile   = gswm_sdi_ncfile
      inp%gswm_nmidi_ncfile = gswm_nmidi_ncfile
      inp%gswm_nmisdi_ncfile= gswm_nmisdi_ncfile
      inp%start_day  = start_day
      inp%start_year = start_year
      inp%calendar_advance = calendar_advance
      inp%step    = step
      inp%dispose = dispose
      inp%ntask_lon = ntask_lon
      inp%ntask_lat = ntask_lat
      inp%mag     = mag
      inp%difhor  = difhor
      inp%iuivi   = iuivi
      inp%dynamo  = dynamo
      inp%nmc     = nmc
      inp%tide(:) = tide(:)
      inp%tide2(:)= tide2(:)
      inp%tide3m3(:) = tide3m3(:)
      inp%tideann = tideann
      inp%aurora  = aurora
      inp%f107    = f107
      inp%f107a   = f107a
      inp%power   = power
      inp%ctpoten = ctpoten
      inp%power_time(:,:) = power_time(:,:)
      inp%ctpoten_time(:,:) = ctpoten_time(:,:)
      inp%byimf   = byimf
      inp%bzimf   = bzimf
      inp%swvel   = swvel
      inp%swden   = swden
      inp%AL      = AL
      inp%colfac  = colfac
!
      inp%source  = source
      inp%output  = output
      inp%source_start  = source_start
      inp%start   = start
      inp%stop    = stop
      inp%hist    = hist
      inp%save    = save
      inp%mxhist_prim = mxhist_prim
      inp%msreten = msreten
!
      inp%secsource= secsource
      inp%secout   = secout
      inp%secstart = secstart
      inp%secstop  = secstop
      inp%sechist  = sechist
      inp%secsave  = secsave
      inp%secflds  = secflds
      inp%secfmag  = secfmag
      inp%secfgeo2d= secfgeo2d
      inp%secfmag2d= secfmag2d
      inp%secfmagphr  = secfmagphr
      inp%mxhist_sech = mxhist_sech
      end subroutine inp_deftyp
!-------------------------------------------------------------------
      subroutine inp_print(ntask)
!
! Print values of inp (input_type):
!
! Args:
      integer,intent(in) :: ntask
! Local:
      integer :: i,n
!
      write(6,"(/,72('-'))")
      write(6,"('USER INPUT PARAMETERS:')")
!
! Model-wide:
      if (len_trim(label) > 0)
     |  write(6,"('  label  = ',a,/,4x,'(optional text label for',
     |    ' current run)')") trim(inp%label)      
      if (len_trim(tempdir) > 0)
     |  write(6,"('  tempdir = ',a,' (runtime temporary directory)')") 
     |    trim(inp%tempdir)      
      if (len_trim(magvol) > 0)
     |  write(6,"('  magvol = ',a,/,4x,
     |    '(file or mss path containing magnetic data)')")
     |    trim(inp%magvol)      
      if (len_trim(potential_model) > 0)
     |  write(6,"('  high-lat electric potential model:  ',
     |        'potential_model = ',a)") trim(potential_model)
      if (len_trim(hrindices_ncfile) > 0)
     |  write(6,"('  hrly indices:  hrindices_ncfile = ',a)")
     |        trim(hrindices_ncfile)
      if (len_trim(weimer_ncfile) > 0)
     |  write(6,"('  weimer coefs: weimer_ncfile = ',a)")
     |        trim(weimer_ncfile)
      if (len_trim(gpi_ncfile) > 0)
     |  write(6,"('  gpi run: gpi_ncfile = ',a)") trim(gpi_ncfile)  
      if (len_trim(gswm_di_ncfile) > 0)
     |  write(6,"('  gswm run diurnal tides: gswm_di_ncfile = ',a)") 
     |        trim(gswm_di_ncfile)
      if (len_trim(gswm_sdi_ncfile) > 0)
     |  write(6,"('  gswm run semidiurnal tides: gswm_sdi_ncfile = ' 
     |        ,a)") trim(gswm_sdi_ncfile)
      if (len_trim(gswm_nmidi_ncfile) > 0)
     |  write(6,"('  gswm run nonmigrating diurnal tides: ',
     |    'gswm_nmidi_ncfile = ',a)") trim(gswm_nmidi_ncfile)
      if (len_trim(gswm_nmisdi_ncfile) > 0)
     |  write(6,"('  gswm run nonmigrating semidiurnal tides: ',
     |    'gswm_nmisdi_ncfile = ',a)") trim(gswm_nmisdi_ncfile)
      if (len_trim(amievol) > 0)
     |  write(6,"('  amievol = ',a,/,4x,
     |    '(file or mss path containing amie data)')") trim(inp%amievol)
      write(6,"('  start_year = ',i4,' (starting calendar day)')")
     |  inp%start_year
      write(6,"('  start_day  = ',i4,' (starting calendar year)')")
     |  inp%start_day
      if (inp%calendar_advance /= 0) then
        write(6,"('  calendar_advance =',i2,' (model will be advanced ',
     |    'in calendar time starting on this day)')") 
     |    inp%calendar_advance
      else
        write(6,"('  calendar_advance  =',i1,' (model will NOT be ',
     |    'advanced in calendar time)')") inp%calendar_advance
      endif

      write(6,"('  step    =',i4,' (model timestep (seconds))',
     |  i4)") inp%step
      write(6,"('  dispose =',i4,' (mss dispose flag 0/1/2)')") 
     |  inp%dispose
#ifdef MPI      
      write(6,"('  ntask_lon = ',i2,' (number of mpi tasks in ',
     |  'longitude dimension.')") ntask_lon
      write(6,"('  ntask_lat = ',i2,' (number of mpi tasks in ',
     |  'latitude  dimension.')") ntask_lat
      write(6,"('  total tasks = ntask_lon*ntask_lat = ',i4)")
     |  ntask
#endif
!
! Primary histories:
      if (len_trim(inp%source) > 0) then
        write(6,"('  source  = ',a,/,4x,'(file or mss path',
     |    ' containing source history)')") trim(inp%source)      
        write(6,"('  source_start = ',(i3,',',i2,',',i2),
     |    ' (model time of source history)')")
     |    inp%source_start
      endif
      n = size(output)-count(output==' ')
      write(6,"('  output (primary history output files) = ',
     |  /,(4x,a,', ',a))") 
     |  (trim(inp%output(i)),i=1,n)
      n = (size(inp%start)-count(inp%start==ispval))/3
      write(6,"('  start (model start times) =',
     |  /,4(4x,i3,',',i2,',',i2))") (inp%start(:,i),i=1,n)
      n = (size(inp%stop)-count(inp%stop==ispval))/3
      write(6,"('  stop (model stop times) =',
     |  /,4(4x,i3,',',i2,',',i2))") (inp%stop(:,i),i=1,n)
      n = (size(inp%hist)-count(inp%hist==ispval))/3
      write(6,"('  hist (primary history disk write frequencies) =',
     |  /,4(4x,i3,',',i2,',',i2))") (inp%hist(:,i),i=1,n)
      n = (size(inp%save)-count(inp%save==ispval))/3
      write(6,"('  save (primary history file save frequencies) =',
     |  /,4(4x,i3,',',i2,',',i2))") (inp%save(:,i),i=1,n)
      write(6,"('  Maxmimum number of histories per primary file = ',
     |  i3)") inp%mxhist_prim
      write(6,"('  Mass store retention period for history files = ',
     |  i5)") inp%msreten
!
! Secondary histories:
      if (len_trim(inp%secsource) > 0) then
        write(6,"('  secsource  = ',a,/,4x,'(file or mss path',
     |    ' containing secsource history)')") trim(inp%secsource)
        write(6,"('  secsource_start = ',(i3,',',i2,',',i2),
     |    ' (model time of secsource history)')")
     |    inp%source_start
      endif
      n = size(secout)-count(secout==' ')
      if (n > 0)
     |  write(6,"('  secout (secondary history output files)=',
     |    /,(4x,a,', ',a))") 
     |    (trim(inp%secout(i)),i=1,n)
      n = (size(inp%secstart)-count(inp%secstart==ispval))/3
      if (n > 0)
     |  write(6,"('  secstart (secondary history start times) =',
     |    /,4(4x,i3,',',i2,',',i2))") (inp%secstart(:,i),i=1,n)
      n = (size(inp%secstop)-count(inp%secstop==ispval))/3
      if (n > 0)
     |  write(6,"('  secstop (secondary history stop times) =',
     |    /,4(4x,i3,',',i2,',',i2))") (inp%secstop(:,i),i=1,n)
      n = (size(inp%sechist)-count(inp%sechist==ispval))/3
      if (n > 0)
     |  write(6,"('  sechist (secondary history disk write',
     |    ' frequencies) =',/,4(4x,i3,',',i2,',',i2))") 
     |    (inp%sechist(:,i),i=1,n)
      n = (size(inp%secsave)-count(inp%secsave==ispval))/3
      if (n > 0)
     |  write(6,"('  secsave (secondary history file save',
     |    ' frequencies) =',/,4(4x,i3,',',i2,',',i2))") 
     |    (inp%secsave(:,i),i=1,n)
      n = (size(inp%secflds)-count(len_trim(inp%secflds)==0))
      if (n > 0)
     |  write(6,"('  secflds (secondary history fields)',
     |    ' =',/,(4x,5a12))") (inp%secflds(i),i=1,n)
      n = (size(inp%secfmag)-count(len_trim(inp%secfmag)==0))
      if (n > 0)
     |  write(6,"('  secfmag (secondary history fields on magnetic ',
     |    'grid) =',/,(4x,5a12))") (inp%secfmag(i),i=1,n)
      n = (size(inp%secfgeo2d)-count(len_trim(inp%secfgeo2d)==0))
      if (n > 0)
     |  write(6,"('  secfgeo2d (secondary history fields on geo. ',
     |    '2d grid) =',/,(4x,5a12))") (inp%secfgeo2d(i),i=1,n)
      n = (size(inp%secfmag2d)-count(len_trim(inp%secfmag2d)==0))
      if (n > 0)
     |  write(6,"('  secfmag2d (secondary history fields on magnetic ',
     |    '2d grid) =',/,(4x,5a12))") (inp%secfmag2d(i),i=1,n)
      n = (size(inp%secfmagphr)-count(len_trim(inp%secfmagphr)==0))
      if (n > 0)
     |  write(6,"('secfmagphr (secondary history fields on magnetphr',
     |    'grid) =',/,(4x,5a12))") (inp%secfmagphr(i),i=1,n)
      write(6,"('  Maximum number of histories per secondary file = ',
     |  i3)") inp%mxhist_sech
!
! More model-wide inputs:
      write(6,"('  mag (lat,lon of south,north magnetic poles) =',
     |  /,4x,4f8.2)") inp%mag
      write(6,"('  difhor = ',i2,' (horizontal eddy diffusion flag)')")
     |  inp%difhor
!
! iuivi flag is replaced by dynamo flag:
!     write(6,"('  iuivi  = ',i2,' (ion drifts momentum flag)')")
!    |  inp%iuivi
      if (inp%dynamo <= 0) then
        write(6,"('  dynamo = ',i2,' (dynamo will NOT be calculated)')")
     |    inp%dynamo
        write(6,"(14x,
     |    '(electric potential and ion drifts will be zero')")
      else
        write(6,"('  dynamo = ',i2,' (dynamo will be calculated)')")
     |    inp%dynamo
      endif
      write(6,"('  nmc    = ',i2,' (flag for NMC boundary',
     |  ' conditions)')") inp%nmc
      write(6,"('  tide (amplitudes and phases of semidiurnal tide) =',
     |  /,4x,5e8.1,5f6.2)") inp%tide
      write(6,"('  tide2 (amplitude and phase of diurnal tide) =',
     |  /,4x,e8.1,f6.2)") inp%tide2
      write(6,"('  tide3m3 (amplitude and phase of 2-day wave)=',
     |  /,4x,e8.1,f6.2)") inp%tide3m3
      write(6,"('  tideann = ',i2,' (0/1 flag for annual tides)')")
     |  inp%tideann
      write(6,"('  aurora  = ',i2,' (0/1 flag for aurora)')") inp%aurora
      write(6,"('  magphr  = ',i2,' (0/1 flag for magnetosphere)')")
     |  inp%magphr
!
! If any of f107,f107a,power,ctpoten are spval, this means GPI database
! will be used for those indices, otherwise will use user-provided
! values:
!
      if (f107 /= spval) then
        write(6,"('  f107    = ',f9.3,' (daily 10.7 cm solar flux)')")
     |    inp%f107
      elseif (len_trim(gpi_ncfile) > 0) then
        write(6,"('  f107 not provided:',
     |    ' will use gpi data file ',a)") trim(gpi_ncfile)
      endif
!
      if (f107a /= spval) then
        write(6,"('  f107a   = ',f9.3,' (81-day ave 10.7 cm solar',
     |    ' flux)')") inp%f107a
      elseif (len_trim(gpi_ncfile) > 0) then
        write(6,"('  f107a not provided:',
     |    ' will use gpi data file ',a)") trim(gpi_ncfile)
      endif
!
      if (power /= spval) then
        if (ntimes_power==0) then
          write(6,"('  power   = ',f9.3,' (hemispheric power (gw)')")
     |      inp%power
        else
          write(6,"('  power   = ',f9.3,' (hemispheric power (gw)',
     |      ' (user provided time series)')") inp%power
        endif
      elseif (len_trim(gpi_ncfile) > 0) then
        write(6,"('  power not provided:',
     |    ' will use gpi data file ',a)") trim(gpi_ncfile)
      endif
!
      if (ctpoten /= spval) then
        if (ntimes_ctpoten==0) then
          write(6,"('  ctpoten = ',f9.3,' (cross-cap potential ',
     |      '(volts)')") inp%ctpoten
        else
          write(6,"('  ctpoten = ',f9.3,' (cross-cap potential ',
     |      '(volts) (user provided time series)')") inp%ctpoten
        endif
      elseif (len_trim(gpi_ncfile) > 0) then
        write(6,"('  ctpoten not provided:',
     |    ' will use gpi data file ',a)") trim(gpi_ncfile)
      endif
!
      write(6,"('  byimf   = ',f9.3,' (BY component of IMF)')")
     |  inp%byimf
      write(6,"('  bzimf   = ',f9.3,' (Bz component of IMF)')")
     |  inp%bzimf
      write(6,"('  swvel   = ',f9.3,' (solar wind velocity)')")
     |  inp%swvel
      write(6,"('  swden   = ',f9.3,' (solar wind density)')")
     |  inp%swden
      write(6,"('  AL      = ',f9.3,' (AL, lower auroral mag index)')")
     |  inp%AL
      write(6,"('  colfac  = ',f9.3,' (collision factor)')")
     |  inp%colfac
!
      write(6,"('END USER INPUT PARAMETERS')")
      write(6,"(72('-'),/)")
      end subroutine inp_print
!-------------------------------------------------------------------
      subroutine validate_mtime(mtime,mxday,label)
!
! Validate a model time (day,hr,min) from input.
! (may be start or stop time, or frequency,
!  e.g., history write frequency)
! If there is a bad value, stop with error message.
! Label is used to print error msg (usually the keyword name
!   from namelist)
!
      integer,intent(in) :: mtime(3),mxday
      character(len=*),intent(in) :: label
      integer :: ier
!
      ier = 0
!
! Day:
      if (mtime(1) < 0 .or. mtime(1) > mxday) then
        write(6,"(/,'>>> input ',a,': bad model day: ',i5,
     |    ' (must be >= 0 and <= mxday)')") label,mtime(1)
        ier = 1
      endif 
!
! Hour:
      if (mtime(2) < 0 .or. mtime(2) > 23) then
        write(6,"(/,'>>> input ',a,': bad model hour: ',i5,
     |    ' (must be >= 0 and <= 23)')") label,mtime(2)
        ier = 1
      endif 
!
! Minute:
      if (mtime(3) < 0 .or. mtime(3) > 59) then
        write(6,"(/,'>>> input ',a,': bad model minute: ',i5,
     |    ' (must be >= 0 and <= 59)')") label,mtime(1)
        ier = 1
      endif 
!
      if (ier > 0) call shutdown('mtime')
      end subroutine validate_mtime
!-----------------------------------------------------------------------
      subroutine usage_calendar
!
! Print usage statement for calendar inputs.
!
      character(len=60) :: char72(10)
      integer :: i
!      123456789-123456789-123456789-123456789-123456789-123456789-
      char72=(/
     |';                                                           ',
     |'; To set calendar start time and control calendar advance,  ',
     |'; please use the following namelist input parameters, e.g., ',
     |';                                                           ',
     |' START_YEAR = 1983  ; starting year (4-digit integer yyyy)  ',
     |' START_DAY  = 80    ; starting day of year (integer 1->365) ',
     |' CALENDAR_ADVANCE=1 ; if 1, advance calendar time           ',
     |'                    ; if 0, do not advance calendar time    ',
     |';                                                           ',
     |'                                                            '/)
      do i=1,10
        if (len_trim(char72(i)) > 0) write(6,"(a)") char72(i)
      enddo 
      end subroutine usage_calendar
!-----------------------------------------------------------------------
      subroutine validate_indices_time(rindex,ntimes,name)
!
! Validate user provided time series for solar index "name". Rindex will
! be either ctpoten_time or power_time, from namelist input.
!
! Args:
      real,intent(in) :: rindex(4,mxind_time)
      integer,intent(out) :: ntimes
      character(len=*),intent(in) :: name
!
! Local:
      integer :: i,ii,n
      integer(kind=8) :: nsec0,nsec1
!
! External:
      integer(kind=8),external :: mtime_to_nsec
!
! Validate times and values in user provided time series:
      ntimes = 0
      if (any(rindex /= spval)) then
!
! Time series must be provided in groups of 4 (day,hr,min,value):
        n = 0
        do i=1,mxind_time
          do ii=1,4
            if (rindex(ii,i) /= spval) n = n+1
          enddo
        enddo
        if (mod(n,4) /= 0) then
          write(6,"('>>> INPUT: must provide ',a,' in groups',
     |      ' of 4: n=',i5,' mod(n,4)=',i5)") name,n,mod(n,4)
          call shutdown('validate_indices_time')
        endif
        do i=1,mxind_time
          if (any(rindex(:,i) /= spval)) then
            call validate_mtime(int(rindex(1:3,i)),mxday,name)
            if (rindex(4,i) <= 0.) then
              write(6,"(/,'>>> INPUT: ',a,' values must be ',
     |          'positive: i=',i3,' rindex(4,i)=',e12.4)")
     |          name,i,rindex(4,i)
              call shutdown('validate_indices_time')
            endif 
            ntimes = ntimes+1
          endif
        enddo ! i=1,mxind_time
!
! Times must be increasing (assume non-spvals are from 1 to ntimes, 
! i.e., not interspersed with spvals):
        if (ntimes > 1) then
          nsec0 = mtime_to_nsec(int(rindex(1:3,1)))
          do i=2,ntimes
            nsec1 = mtime_to_nsec(int(rindex(1:3,i)))
            if (nsec0 >= nsec1) then
              write(6,"(/,'>>> INPUT: ',a,' times must increase',
     |          '. Check time at i=',i3,' day,hr,min=',3i4,' and the',
     |          ' time previous to i.')") name,i,int(rindex(1:3,i)) 
              call shutdown('validate_indices_time')
            endif
            nsec0 = nsec1
          enddo ! i=1,ntimes
        endif ! ntimes > 1
      endif ! ctpoten_time
      end subroutine validate_indices_time
!-----------------------------------------------------------------------
      end module input_module
