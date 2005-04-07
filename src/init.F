!
      module init_module
      use params_module,only: nlon,nlat,nlevp1,nlonp4,nmlonp1,nmlat,
     |  nmlev,glon1,dlon,glat1,dlat,nlev,plev1,dlev
      implicit none
!
! Initialize and store model non-input variables and constants:
!
      real :: glon(nlon),glat(nlat),plev(nlevp1)
      real :: gmlon(nmlonp1),gmlat(nmlat),pmlev(nmlev)
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
     |  start_mtime(3)     ! starting model time (day,hr,min)
!
! getgpi flag must go here rather than in gpi_mod to avoid circular 
! module dependency with nchist_mod and init_mod.
!
      integer :: igetgpi       ! 0/1 flag for GPI data (see gpi_ncfile
      integer :: igetgswmdi    ! 0/1 flag for GSWM data diurnal tide
      integer :: igetgswmsdi   ! 0/1 flag for GSWM data semidiurnal tide
      integer :: igetgswmnmdi  ! 0/1 flag for GSWM data nonmigrating diurnal tide
      integer :: igetgswmnmsdi ! 0/1 flag for GSWM data nonmigrating semidiurnal tide
      integer :: igetgswm ! will be 1 only if at least one of the above flags are set
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
      use input_module,only: start,step,secflds,secfmag,secfmagphr,
     |  start_year,start_day,calendar_advance,gpi_ncfile,gswm_di_ncfile,
     |  gswm_sdi_ncfile,gswm_nmdi_ncfile,gswm_nmsdi_ncfile,
     |  mxhist_prim,mxhist_sech,output,secout,mkhvols,nmc,source_start,
     |  sd_ncfile
      use hist_module,only: hist_init,isechist,nfsech_geo,nfsech_mag,
     |  nfsech_geo2d,nfsech_mag2d,nfsech_magphr,nstep,nhist_total,
     |  nsech_total,nsource,nseries_prim,nseries_sech,nfiles_prim,
     |  nfiles_sech
      use fields_module,only: set_fsech,init_4d,init_3d,
     |  fsech,fsechmag,fsech2d,fsechmag2d,fsechmagphr2d
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
! Set GPI flag:
      igetgpi = 0
      if (len_trim(gpi_ncfile) > 0) igetgpi = 1
      if (igetgpi > 0) 
     |  write(6,"('  gpi_ncfile = ',a)") trim(gpi_ncfile)
!
! Set GSWM flag for diurnal tides:
      igetgswmdi = 0
      if (len_trim(gswm_di_ncfile) > 0) igetgswmdi = 1
      if (igetgswmdi > 0) 
     |  write(6,"('  gswm_di_ncfile = ',a)") trim(gswm_di_ncfile)
! Set GSWM flag for semidiurnal tides:
      igetgswmsdi = 0
      if (len_trim(gswm_sdi_ncfile) > 0) igetgswmsdi = 1
      if (igetgswmsdi > 0) 
     |  write(6,"('  gswm_sdi_ncfile = ',a)") trim(gswm_sdi_ncfile)
! Set GSWM flag for nonmigrating diurnal tides:
      igetgswmnmdi = 0
      if (len_trim(gswm_nmdi_ncfile) > 0) igetgswmnmdi = 1
      if (igetgswmnmdi > 0) 
     |  write(6,"('  gswm_nmdi_ncfile = ',a)") trim(gswm_nmdi_ncfile)
! Set GSWM flag for nonmigrating semidiurnal tides:
      igetgswmnmsdi = 0
      if (len_trim(gswm_nmsdi_ncfile) > 0) igetgswmnmsdi = 1
      if (igetgswmnmsdi > 0) 
     |  write(6,"('  gswm_nmsdi_ncfile = ',a)")  
     |  trim(gswm_nmsdi_ncfile)
!
! If any gswm tides are used, set igetgswm==1:
      igetgswm = 0
      if (igetgswmdi>0   .or. igetgswmsdi>0 .or.
     |    igetgswmnmdi>0 .or. igetgswmnmsdi>0) igetgswm = 1

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
      if (len_trim(sd_ncfile) > 0) then
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
!
! Initialize field structures:
      iprintf = 0
      call init_4d(lon0,lon1,lat0,lat1,mytid,iprintf)
      call init_3d(lon0,lon1,lat0,lat1,mytid,iprintf)
!
! Initialize secondary history fields:
      if (isechist > 0) then
        call set_fsech
      endif
!
! Initialize sun's declination:
      sundec=atan(tan(23.5*pi/180.)*sin(2.*pi*float(iday-80)/365.))
      sin_sundec = SIN(sundec) ! C(95)
      cos_sundec = COS(sundec) ! C(96)
!
! Define geographic grid:
      do i=1,nlon
        glon(i) = glon1+(i-1)*dlon
      enddo
      do i=1,nlat
        glat(i) = glat1+(i-1)*dlat
      enddo
      do i=1,nlev
        plev(i) = plev1+(i-1)*dlev
      enddo
!
! Report to stdout:
      write(6,"(/,'Model run initialization:')")
      write(6,"('  nstep      = ',i6,4x,
     |  '(Number of time steps this run)')") nstep
      write(6,"('  iter       = ',i6,4x,
     |  '(Initial iteration number)')") iter
      write(6,"('  iyear      = ',i6,4x,
     |  '(Beginning calendar year)')") iyear
      write(6,"('  iday       = ',i6,4x,
     |  '(Beginning calendar day)')") iday
      write(6,"('  igetgpi = ',i6,4x,
     |  '(If > 0, geophysical indices database will be used.)')")
     |  igetgpi
      write(6,"('  igetgswmdi = ',i6,4x,
     |  '(If > 0, GSWM diurnal tidal database will be used.)')")
     |  igetgswmdi
      write(6,"('  igetgswmsdi= ',i6,4x,
     |  '(If > 0, GSWM semidiurnal tidal database will be used.)')")
     |  igetgswmsdi
      write(6,"('  igetgswmnmdi= ',i6,4x,
     |  '(If > 0, GSWM nonmigrating diurnal tidal database will',
     |  ' be used.)')") igetgswmnmdi
      write(6,"('  igetgswmnmsdi= ',i6,4x,
     |  '(If > 0, GSWM nonmigrating semidiurnal tidal database will',
     |  ' be used.)')") igetgswmnmsdi
!
! ncep/nmc are in time-gcm only:
!     write(6,"('  ncep    = ',i6,4x,
!    |  '(If > 0, use NCEP Z and TN 10 mb lower boundaries.')")
!    |  ncep
!     write(6,"('  nmc     = ',i6,4x,
!    |  '(If > 0, use NMC  Z and TN 10 mb lower boundaries.')")
!    |  nmc
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
        write(6,"('  nfsech_geo    = ',i5,2x,
     |    '(Number of secondary history fields on geographic grid)')") 
     |    nfsech_geo
        write(6,"('  nfsech_mag    = ',i5,2x,
     |    '(Number of secondary history fields on magnetic grid)')") 
     |    nfsech_mag
        write(6,"('  nfsech_geo2d    = ',i5,2x,
     |    '(Number of secondary history fields on geographic ',
     |    '2d grid)')")nfsech_geo2d
        write(6,"('  nfsech_mag2d    = ',i5,2x,
     |    '(Number of secondary history fields on magnetic 2d grid)')")
     |    nfsech_mag2d
        write(6,"('  nfsech_magphr   = ',i5,2x,
     |   '(Number of secondary history fields on magnetoshere grid)')")
     |    nfsech_magphr
!
! Report secondary history fields: 
        write(6,"(/,'Secondary history fields on geographic grid ',
     |    ' (number of fields =',i3,'):')") nfsech_geo
        do i=1,nfsech_geo
          if (len_trim(fsech(i)%long_name)>0) then
            write(6,"('  Field ',a,' (',a,')')") 
     |        fsech(i)%short_name(1:8),trim(fsech(i)%long_name)
          else
            write(6,"('  Field ',a)") fsech(i)%short_name(1:8)
          endif
        enddo
        write(6,"(/,'Secondary history fields on magnetic ',
     |    'grid (number of fields =',i3,'):')") nfsech_mag
        do i=1,nfsech_mag
          write(6,"('  Mag field ',a,' (diagnostic)')") 
     |      fsechmag(i)%short_name
        enddo
        write(6,"(/,'Secondary history fields on geographic 2d grid ',
     |    ' (number of fields =',i3,'):')") nfsech_geo2d
        do i=1,nfsech_geo2d
          write(6,"('  Field ',a,' (diagnostic)')")
     |      fsech2d(i)%short_name
        enddo
        write(6,"(/,'Secondary history fields on magnetic ',
     |    '2d grid (number of fields =',i3,'):')") nfsech_mag2d
        do i=1,nfsech_mag2d
          write(6,"('  Mag 2d field ',a,' (diagnostic)')")
     |      fsechmag2d(i)%short_name
        enddo
        write(6,"(/,'Secondary history fields on magnetospheric ',
     |    '2d grid  (number of fields =',i3,'):')") nfsech_magphr
        do i=1,nfsech_magphr
          write(6,"('  Magphr field ',a,' (diagnostic)')")
     |      fsechmagphr2d(i)%short_name
        enddo
      endif
      end subroutine init
!-----------------------------------------------------------------------
      end module init_module
