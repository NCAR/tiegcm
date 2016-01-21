      module init_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlon,nlat,nlonp4
      implicit none
!
! istep: the current time step index. istep is initialized
! to 0, then incremented before each time step (see advance.f).
!
      integer :: istep     ! time step index
!
! iter (iteration number): the number of timesteps (at the current 
! step length) from model time 0,0,0 to current model time (day,hour,min).
!
      integer :: iter
!
! iyear and iday represent the current model calendar date
! (iyear is 4-digit). Uthr is current decimal hours. These
! are updated per timestep in advnce.
!
      integer :: iyear     ! Current model calendar 4-digit year
      integer :: iday      ! Current model calendar day
      real :: uthr         ! Current ut (decimal hours)
      integer :: 
     |  start_mtime(3),    ! starting model time (day,hr,min)
     |  stop_mtime(3)      ! stop model time (day,hr,min)
!
      integer :: igswm_mi_di  ! 0/1 flag for GSWM data diurnal tide
      integer :: igswm_mi_sdi ! 0/1 flag for GSWM data semidiurnal tide
      integer :: igswm_nm_di  ! 0/1 flag for GSWM data nonmigrating diurnal tide
      integer :: igswm_nm_sdi ! 0/1 flag for GSWM data nonmigrating semidiurnal tide
      integer :: igetgswm     ! will be 1 only if at least one of the above flags are set
      integer :: ictmt        ! 0/1 flag to get CTMT data
      integer :: ibgrddata    ! 0/1 flag to get background data
!
      character(len=16) ::
     |  host,              ! host machine
     |  system,            ! operating system of host (from pre-proc macros)
     |  logname            ! user login name
      character(len=8) ::
     |  rundate,           ! current local date of run
     |  runtime            ! current local time of run
!
! Secs is updated in advnce, and is referenced in opflux, settei,
!   sun, and chapmn.
      real :: 
     |  secs,           ! current ut in seconds
     |  sfeps,          ! solar output change due to orbital eccentricity
!    |  alfalp,efluxlp  ! low-energy protons in aurora
     |  sundec,         ! sun's declination (updated in advnce)
     |  sin_sundec,     ! sin(sundec) (updated in advnceday)
     |  cos_sundec      ! cos(sundec) (updated in advnceday)
!
! Day/night index is set by chapman.F:
      integer :: idn(nlonp4)        ! day/night index
!
      contains
!-----------------------------------------------------------------------
      subroutine init
!
! Initialize (this is called by tgcm.F after input):
! Some init also takes place in inp_model (input_mod.f)
!
      use input_module,only: start,stop,step,secflds,
     |  start_year,start_day,calendar_advance,gpi_ncfile,
     |  gswm_mi_di_ncfile,ctmt_ncfile,
     |  gswm_mi_sdi_ncfile,gswm_nm_di_ncfile,gswm_nm_sdi_ncfile,
     |  mxhist_prim,mxhist_sech,output,secout,mkhvols,source_start,
     |  see_ncfile,bgrddata_ncfile
      use hist_module,only: hist_init,isechist,nstep,nhist_total,
     |  nsech_total,nsource,nseries_prim,nseries_sech,nfiles_prim,
     |  nfiles_sech,nfsech
      use fields_module,only: init_4d,init_3d,init_lbc,init_fsech,
     |  fsechist
      use cons_module,only: pi,init_cons
      use filter_module,only: trigs,ntrigs,ifax
      use mpi_module,only: lon0,lon1,lat0,lat1,mytid
      use soldata_module,only: rd_soldata
!
! External:
      integer,external :: mtime_to_nstep
!
! Local:
      real :: theta0
      integer :: i,iprintf
!
#if (MPI == 0)
      lon0 = 1
      lon1 = nlonp4
      lat0 = 1
      lat1 = nlat
#endif
!
! Initialize derived constants (init_cons is in cons_module):
      call init_cons
!
! Get login name:
      logname = ' '
      call getenv('LOGNAME',logname)
      if (len_trim(logname)==0) then
        write(6,"(/,'>>> init: Cannot get LOGNAME environment ',
     |    'variable.',/)")
        call shutdown('LOGNAME')
      endif
!
! Get host name:
      call gethostsname(host)
!
! Operating system (based on pre-proc macro):
      call setosys(system)
!
! Get run date (current date and time):
      call datetime(rundate,runtime)
!
! Iter is the number of time steps from 0,0,0 to the current model
! time, using the current step length. Iter is incremented once per 
! timestep in advnce.
!
      iter = mtime_to_nstep(start(:,1),step)
!
! iyear and iday are current calendar year and day.
! If calendar_advance > 0, the model is advanced in calendar time, 
!   starting at start_day. If calendar_advance==0, model is NOT 
!   advanced in calendar time (start_day is held constant).
! iyear and iday are incremented in advance if the model is
!   advancing in calendar time.
!
      iyear = start_year ! from input
      iday  = start_day  ! from input
      write(6,"('init: iyear=',i5,' iday=',i5)") iyear,iday
      sfeps = 1.
!
! If model is being advanced in calendar time, initialize orbital 
! eccentricity.
!
      if (calendar_advance > 0) then
        theta0 = 2.*pi*float(iday)/365.
        sfeps = 1.000110+0.034221*cos(theta0)+0.001280*sin(theta0)+
     |          0.000719*cos(2.*theta0)+0.000077*sin(2.*theta0)
      endif
!
! 2/00: these were in modsrc.snoe (tgcm13mt), but were unused.
! Low-energy protons:
!     alfalp = 10.
!     efluxlp = 1.e-20
!
! GSWM integer flags. 11/2/05 btf: these have been taken off
!   the histories (replaced by the 4 file names, see input.F
!   and nchist.F)
! 
! Set GSWM flag for diurnal tides:
      igswm_mi_di = 0
      if (len_trim(gswm_mi_di_ncfile) > 0) igswm_mi_di = 1
      if (igswm_mi_di > 0) 
     |  write(6,"('  gswm_mi_di_ncfile = ',a)") trim(gswm_mi_di_ncfile)
! Set GSWM flag for semidiurnal tides:
      igswm_mi_sdi = 0
      if (len_trim(gswm_mi_sdi_ncfile) > 0) igswm_mi_sdi = 1
      if (igswm_mi_sdi > 0) 
     |  write(6,"('  gswm_mi_sdi_ncfile = ',a)")trim(gswm_mi_sdi_ncfile)
! Set GSWM flag for nonmigrating diurnal tides:
      igswm_nm_di = 0
      if (len_trim(gswm_nm_di_ncfile) > 0) igswm_nm_di = 1
      if (igswm_nm_di > 0) 
     |  write(6,"('  gswm_nm_di_ncfile = ',a)") trim(gswm_nm_di_ncfile)
! Set GSWM flag for nonmigrating semidiurnal tides:
      igswm_nm_sdi = 0
      if (len_trim(gswm_nm_sdi_ncfile) > 0) igswm_nm_sdi = 1
      if (igswm_nm_sdi > 0) 
     |  write(6,"('  gswm_nm_sdi_ncfile = ',a)")  
     |  trim(gswm_nm_sdi_ncfile)
!
! If any gswm tides are used, set igetgswm==1:
      igetgswm = 0
      if (igswm_mi_di>0 .or. igswm_mi_sdi>0 .or.
     |    igswm_nm_di>0 .or. igswm_nm_sdi>0) igetgswm = 1
!
! Background saber/tidi, msis/hwm:
      ibgrddata = 0
      if (len_trim(bgrddata_ncfile) > 0) ibgrddata = 1
      if (ibgrddata > 0)
     |  write(6,"('  bgrddata_ncfile = ',a)") trim(bgrddata_ncfile)
!
! CTMT lbc data file:
      ictmt = 0
      if (len_trim(ctmt_ncfile) > 0) ictmt = 1
      if (ictmt > 0)
     |  write(6,"('  ctmt_ncfile = ',a)") trim(ctmt_ncfile)
!
! ixtimep is 4th dimension index to fg-array for previous time step
! ixtimec is 4th dimension index to fg-array for current time step
! (see fogcm.f)
!
!     ixtimep = 1
!     ixtimec = 1
!
! Initialize amie, and get amie file if necessary:
!     call init_amie
!
! Read SEE data file:
      if (len_trim(see_ncfile) > 0) then
         call rd_soldata
      endif
!
! Initialize non-input history variables for beginning of run:
      call hist_init
!
! setfft calls set99 for fft init. This call returns trigs and ifax,
! in filter_module.F.
!
      call setfft(trigs,ifax,ntrigs,nlon)
!
! Init starting model time:
      start_mtime = start(:,1)
      stop_mtime = stop(:,1)
!
! Initialize field structures:
      iprintf = 0
      call init_4d(lon0,lon1,lat0,lat1,mytid,iprintf)
      call init_3d(lon0,lon1,lat0,lat1,mytid,iprintf)
!
! Initialize secondary history fields:
      if (isechist > 0) call init_fsech ! fields.F
!
! Allocate lbc for t,u,v:
      call init_lbc(lon0,lon1,lat0,lat1)
!
! Initialize sun's declination:
      sundec=atan(tan(23.5*pi/180.)*sin(2.*pi*float(iday-80)/365.))
      sin_sundec = SIN(sundec) ! C(95)
      cos_sundec = COS(sundec) ! C(96)
!
! Report to stdout:
      write(6,"(/,'Model run initialization:')")
      write(6,"('  nstep      = ',i12,4x,
     |  '(Number of time steps this run)')") nstep
      write(6,"('  iter       = ',i12,4x,
     |  '(Initial iteration number)')") iter
      write(6,"('  iyear      = ',i6,4x,
     |  '(Beginning calendar year)')") iyear
      write(6,"('  iday       = ',i6,4x,
     |  '(Beginning calendar day)')") iday
      write(6,"('  igswm_mi_di = ',i6,4x,
     |  '(If > 0, GSWM diurnal tidal database will be used.)')")
     |  igswm_mi_di
      write(6,"('  igswm_mi_sdi= ',i6,4x,
     |  '(If > 0, GSWM semidiurnal tidal database will be used.)')")
     |  igswm_mi_sdi
      write(6,"('  igswm_nm_di= ',i6,4x,
     |  '(If > 0, GSWM nonmigrating diurnal tidal database will',
     |  ' be used.)')") igswm_nm_di
      write(6,"('  igswm_nm_sdi= ',i6,4x,
     |  '(If > 0, GSWM nonmigrating semidiurnal tidal database will',
     |  ' be used.)')") igswm_nm_sdi
      write(6,"('  ibgrddata= ',i6,4x,
     |  '(If > 0, background datafile will',
     |  ' be used.)')") ibgrddata
      write(6,"('  ictmt= ',i6,4x,
     |  '(If > 0, CTMT tidal database will be used.)')") ictmt
!
      if (nsource > 0) then
        write(6,"(/,'This is an initial run:')")
        write(6,"('  start_year = ',i6,5x,
     |    '(Starting year of initial run)')") start_year
        write(6,"('  start_day  = ',i6,5x,
     |    '(Starting day of initial run)')") start_day
        write(6,"('  start_mtime= ',i4,2i3,1x,
     |    '(Starting mtime of initial run)')") start_mtime
      endif
!
! Report re primary histories to stdout:
      write(6,"(/,'Primary Histories:')")
      write(6,"('  nsource       = ',i5,2x,
     |  '(If > 0, a primary source history was provided)')") nsource
      write(6,"('  nseries_prim  = ',i5,2x,
     |  '(Number of primary time series)')") nseries_prim
      write(6,"('  nhist_total   = ',i5,2x,
     |  '(Number of primary histories to be written)')") nhist_total
      write(6,"('  nfiles_prim   = ',i5,2x,
     |  '(Number of primary output files to be written)')") nfiles_prim
      write(6,"('  mxhist_prim   = ',i5,2x,
     |  '(Maximum number of primary histories per file)')") mxhist_prim
!
! Report re secondary histories to stdout:
      if (isechist > 0) then
        write(6,"(/,'Secondary Histories:')")
        write(6,"('  nseries_sech  = ',i5,2x,
     |    '(Number of secondary time series)')") nseries_sech
        write(6,"('  nsech_total   = ',i5,2x,
     |    '(Number of secondary histories to be written)')") nsech_total
        write(6,"('  nfiles_sech   = ',i5,2x,
     |    '(Number of secondary output files to be written)')")
     |    nfiles_sech
        write(6,"('  mxhist_sech   = ',i5,2x,
     |    '(Maximum number of secondary histories per file)')")
     |    mxhist_sech
        write(6,"('  nfsech        = ',i5,2x,
     |    '(Number of requested secondary history fields)')") nfsech
        do i=1,nfsech
          write(6,"('  secondary history field ',i3,': ',a)") 
     |      i,fsechist(i)%short_name
        enddo
      endif
      end subroutine init
!-----------------------------------------------------------------------
      end module init_module
