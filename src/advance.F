!
      subroutine advance
!
! Advance the model nstep time steps.
!
      use fields_module
      use hist_module,only: nstep,modeltime,nsource
      use heelis_module,only: heelis,colath
      use weimer_module,only: weimer01
      use input_module,only: start,step,iaurora=>aurora,idynamo=>dynamo,
     |  f107,f107a,ctpoten,power,calendar_advance,ntimes_ctpoten,
     |  ctpoten_time,ntimes_power,power_time,potential_model
      use init_module,only: istep,uthr,iter,secs,iday,iyear,igetgpi,
     |  igetgswmdi,igetgswmsdi,igetgswmnmidi,igetgswmnmisdi
      use cons_module,only: dt,re,cs,racs,kut,ylatm,rtd
      use filter_module,only: filter
      use qrj_module,only: init_sflux
      use magfield_module,only: sunloc
      use efield_module,only: efield
      use aurora_module,only: aurora_cons,cusp2d,drzl2d,alfa2d,nflx2d
      use hdif_module,only: hdif1,hdif2
      use timing_module,only: elapsed_steps,elapsed_mpi,
     |  start_timing, end_timing
      use dynamo_module,only: prep_dynamo,prep_dynamo_dyn0,dynamo,
     |  nodynamo,phihm,nmlat0,zpotenm3d
      use gpi_module,only: rdgpi,getgpi
      use gswm_module,only: rdgswmdi,rdgswmsdi,rdgswmnmidi,
     |  rdgswmnmisdi,getgswmdi,getgswmsdi,getgswmnmidi,
     |  getgswmnmisdi
      use params_module,only: nlonp4,nlat
      use output,only: outhist
!
! Timing stats from mpi module:
      use mpi_module,only: 
     |  do_rtc_mpi            ,
     |  rtcmp_gather2root     ,
     |  rtcmp_bndlats         ,
     |  rtcmp_bndlons         ,
     |  rtcmp_bndlons_f3d     ,
     |  rtcmp_polelat         ,
     |  rtcmp_gatherlons      ,
     |  rtcmp_gatherlons_f3d  ,
     |  rtcmp_scatterlons     ,
     |  rtcmp_scatterlons_f3d ,
     |  rtcmp_periodic_f4d    ,
     |  rtcmp_periodic_f3d    ,
     |  rtcmp_periodic_f2d    ,
     |  rtcmp_dynpot
#ifdef MPI
      use mpi_module,only: 
     |  mytid,mp_gather2root,mp_bndlats,ntask,lat0,lat1,lon0,lon1,
     |  mp_bndlons,mp_gatherlons,mp_scatterlons,mp_polelats,
     |  mp_periodic_f4d,mp_updatephi
!
      implicit none
#include <mpif.h>
#else
      implicit none
      integer :: lat0=1,lat1=nlat, lon0=1,lon1=nlonp4
      integer :: mytid=0,ntask=1
#endif
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j,k,lat,itmp,icount_step,nlatlocal,icount_mpi,ier,
     |  iprint
      real :: fmin,fmax,elapsed_step,elapsed_gather2root
      integer(kind=8) :: nsecs ! current model time in seconds,
                               ! including model day
      logical :: time2write,wrprim,wrsech,save_prim,save_sech,
     |  newseries_sech,newseries_prim,iseries_sech,iseries_prim
      real :: w_ik(nlonp4,nlevp1),avtime
!
! If do_rtc is true, do real-time-clock timing. rtc is an AIX function,
!   called by sub timer (util.F). Sub timer brackets rtc calls with
!   mpi_barrier. Do_rtc is passed to sub dynamics.
! Note: non-rtc timing (e.g., elapsed_step, elapsed_steps, etc) use
!   the system_clock, called from the timing module (timing.F).
!
      logical :: do_rtc=.false.
      real :: tsec, rtc_dynamics, rtc_step, rtc_elapstep,
     |  rtc_elapsteps, rtc_polelat, rtc_bndlats, rtc_bndlons, 
     |  rtc_addiag, rtc_hdif1, rtc_hdif2, rtc_gather2root,
     |  rtc_dynamo
      logical :: debug=.false.
!
! External:
      logical,external :: wrhist
!
#ifdef VT
!     code = 115 ; state = 'advance' ; activity='ModelCode'
      call vtbegin(115,ier)
#endif
!
      write(6,"(/,'Enter advance: iter=',i6,' nstep=',i6)") iter,nstep
!
! Init total elapsed timing (sum of elapsed for each timestep):
!     elapsed_steps = 0. ! uses system_clock in timing module (timing.F)
      rtc_elapsteps = 0. ! uses real-time-clock, see timer in util.F (local)
      avtime = 0.
!
! Init timestep index, modeltime, and nsecs:
      istep = 0
      modeltime(1:3) = start(:,1)
      modeltime(4) = 0
      elapsed_step = 0.
      call modeltime_to_nsecs(modeltime,nsecs) ! sets nsecs, includes day
!
! Echo source history to primary history output if an initial run,
!   and write initial secondary history if necessary (istep==0).
! Source history will not be written if a continuation run (i.e., if
!   SOURCE was *not* provided by user input.
! (if mpi, only master writes histories)
!
! Source was read into data(:,:,:,itp), and output always writes from
! data(:,:,:,itc), so swap indices to echo itp source data:
        itmp = itp
        itp = itc
        itc = itmp
#ifdef MPI
!
! If an initial run, echo source history to output before first time step. 
!   Because only subdomains are read by rdhist, and only the root task 
!   outputs histories, the subdomain data must be gathered to the root.
!
        call mp_gather2root(itc,'prim')
        if (mytid==0) call outhist(istep,modeltime)
#else
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,itc)
        enddo
        call outhist(istep,modeltime)
#endif
!
! nstep may be zero if user only wants to copy source history:
      if (nstep==0) then
        write(6,"(/,'ADVANCE: model is not advanced in time because ',
     |    ' start==stop.')")
        return
      endif
!
! Reswap indices, so model reads itp data and updates itc data:
      itmp = itp
      itp = itc
      itc = itmp
!
! Main time loop:
  100 continue
#ifdef VT
!     code = 127 ; state = 'timestep' ; activity='ModelCode'
      call vtbegin(127,ier)
#endif
      call start_timing(icount_step,'time-steps',1)
      if (do_rtc) call timer(rtc_step,tsec,'begin')
      iter=iter+1
      istep = istep+1
!
! Increment model time by one time step (day,hr,min,secs):
!
      nsecs = nsecs+step
      call nsecs_to_modeltime(nsecs,modeltime) ! increments modeltime(4)
      secs = mod(iter*int(dt),86400)
      uthr = secs/3600.

!     write(6,"('step=',i4,' iter=',i5,' int(dt)=',i3,' modeltime=',
!    |  3i4,' secs=',f10.2,' uthr=',f5.2)") istep,iter,int(dt),
!    |  modeltime(1:3),secs,uthr

!
! Update calendar day and sun's location and declination:
! (sub advance_day also updates sfeps)
!
      if (calendar_advance > 0) call advance_day
      call sunloc(iday,secs) ! locate sun's longitude
      if (debug) write(6,"('advance after sunloc')")
!
! Interpolate power and/or ctpoten to current model time, if time-dependent
! values were read from input:
      if (ntimes_ctpoten > 0) 
     |  call set_index(ctpoten_time,ntimes_ctpoten,nsecs,ctpoten,
     |    'ctpoten')
      if (ntimes_power > 0) 
     |  call set_index(power_time,ntimes_power,nsecs,power,'power')
!
! Get gpi data if necessary:
      if (igetgpi > 0) then
        if (istep==1) call rdgpi ! read gpi data
        iprint = 0
        if (istep==1) iprint = 1
        call getgpi(iyear,iday,int(secs),iprint)
      endif
!
! Get gswm diurnal tidal data if necessary:
      if (igetgswmdi > 0) then
        if (istep==1) call rdgswmdi ! read gswm diurnal tidal data
        iprint = 0
        if (istep==1) iprint = 1
        call getgswmdi(iyear,iday,int(secs),iprint)
      endif
!
! Get gswm semidiurnal tidal data if necessary:
      if (igetgswmsdi > 0) then
        if (istep==1) call rdgswmsdi ! read gswm semidiurnal tidal data
        iprint = 0
        if (istep==1) iprint = 1
        call getgswmsdi(iyear,iday,int(secs),iprint)
      endif
!
! Get gswm nonmigrating diurnal tidal data if necessary:
      if (igetgswmnmidi > 0) then
        if (istep==1) call rdgswmnmidi ! read gswm nonmigrating diurnal tidal data
        iprint = 0
        if (istep==1) iprint = 1
        call getgswmnmidi(iyear,iday,int(secs),iprint)
      endif
!
! Get gswm nonmigrating semidiurnal tidal data if necessary:
      if (igetgswmnmisdi > 0) then
        if (istep==1) call rdgswmnmisdi ! read gswm nonmigrating semidiurnal tidal data
        iprint = 0
        if (istep==1) iprint = 1
        call getgswmnmisdi(iyear,iday,int(secs),iprint)
      endif
!
! Report to stdout:
      if (istep == 1) then
        write(6,"(/,'Step ',i6,' of ',i6,' mtime=',2i4,2i3)")
     |    istep,nstep,modeltime
      else
        if (nstep <= 100 .or. (nstep > 100 .and. mod(istep,10)==0)) then
          if (do_rtc) then
            write(6,"('Step ',i6,' of ',i6,' mtime=',2i4,2i3,
     |        ' rtc secs/step=',f6.2,' elapsed mins=',f8.2)") 
     |        istep,nstep,modeltime,rtc_elapstep,rtc_elapsteps/60.
          else
            write(6,"('Step ',i6,' of ',i6,' mtime=',2i4,2i3,
     |        ' secs/step=',f6.2,' elapsed mins=',f8.2)")
     |        istep,nstep,modeltime,elapsed_step,elapsed_steps/60.
          endif
          if (igetgpi > 0) write(6,"('gpi: f107=',f7.2,' f107a=',f7.2,
     |      ' ctpoten=',f7.2,' power=',f7.2)") 
     |      f107,f107a,ctpoten,power
          if (ntimes_ctpoten > 0.or.ntimes_power > 0) 
     |      write(6,"('ctpoten=',f8.2,' power=',f8.2)") ctpoten,power
        endif
      endif
!
! Update aurora parameters:
      iprint = 0
      if (istep==1) iprint = 1
      call aurora_cons(iprint)
      if (debug) write(6,"('advance after aurora_cons')")
!
! init_sflux calls ssflux and euvac (see qrj_module.F):
      call init_sflux
      if (debug) write(6,"('advance after init_sflux')")
!
! Calculate electric field (ex,ey,ez in efield module, efield.F):
      call efield(1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('advance after efield')")
!
! Calculate vc, barm, xnmbar, and z (tn,o2,o1,vn are input).
!
      if (do_rtc) call timer(rtc_addiag,tsec,'begin')
!    
      if(igetgswmdi==0.and.igetgswmsdi==0.and.
     |   igetgswmnmidi==0.and.igetgswmnmisdi==0) then ! no GSWM input used 
        call addiag(
     |    tn	 (levd0,lond0,latd0,itp),     ! in
     |    o2	 (levd0,lond0,latd0,itp),     ! in
     |    o1	 (levd0,lond0,latd0,itp),     ! in
     |    vn	 (levd0,lond0,latd0,itp),     ! in
     |    vc	 (levd0,lond0,latd0,itp),     ! out
     |    barm   (levd0,lond0,latd0,itp),     ! out
     |    xnmbar (levd0,lond0,latd0),	      ! out
     |    xnmbari(levd0,lond0,latd0),	      ! out
     |    xnmbarm(levd0,lond0,latd0),	      ! out
     |    z	 (levd0,lond0,latd0,itp),     ! out (note itp is output)
     |    lon0,lon1,1,nlevp1,lat0,lat1)
!   
      else 				             ! GSWM input used 
        call addiag_gswm(
     |    tn	 (levd0,lond0,latd0,itp),     ! in
     |    o2	 (levd0,lond0,latd0,itp),     ! in
     |    o1	 (levd0,lond0,latd0,itp),     ! in
     |    vn	 (levd0,lond0,latd0,itp),     ! in
     |    vc	 (levd0,lond0,latd0,itp),     ! out
     |    barm   (levd0,lond0,latd0,itp),     ! out
     |    xnmbar (levd0,lond0,latd0),	      ! out
     |    xnmbari(levd0,lond0,latd0),	      ! out
     |    xnmbarm(levd0,lond0,latd0),	      ! out
     |    z	 (levd0,lond0,latd0,itp),     ! out (note itp is output)
     |    lon0,lon1,1,nlevp1,lat0,lat1)
      endif
!      
      if (debug) write(6,"('advance after addiag')")
      if (do_rtc) then
        call timer(rtc_addiag,tsec,'end')
        write(6,"('Advance step ',i4,': rtc time for addiag=',f10.2)")
     |    istep,tsec 
      endif
!
! Define boundary latitudes -1,0 across the south pole, and
!   nlat+1,nlat+2 across the north pole:
!
#ifdef MPI
      if (do_rtc) call timer(rtc_polelat,tsec,'begin')
      call mp_polelats(itp)
      if (do_rtc) then
        call timer(rtc_polelat,tsec,'end')
        write(6,"('Advance step ',i4,': rtc time for pole_lat=',f10.2)")
     |    istep,tsec 
      endif
!
! Update ghost cell boundary lats and lons.
! Bndlats is called first to exchange lat0-1,2 and lat1+1,2 at lon0->lon1, 
!   then bndlons is called  to exchange lon0-1,2 and lon1+1,2 at lat0-2->lat1+2).
!
      f4d(:)%mpi = .true. 
      if (do_rtc) call timer(rtc_bndlats,tsec,'begin')
      call mp_bndlats(f4d,nf4d,itp)
      if (do_rtc) then
        call timer(rtc_bndlats,tsec,'end')
        write(6,"('Advance step ',i4,': rtc time for mp_bndlats=',
     |    f10.2,' (beginning of step)')") istep,tsec
      endif
!
      if (do_rtc) call timer(rtc_bndlons,tsec,'begin')
      call mp_bndlons(f4d,nf4d,itp)
      if (do_rtc) then
        call timer(rtc_bndlons,tsec,'end')
        write(6,"('Advance step ',i4,': rtc time for mp_bndlons=',
     |    f10.2,' (beginning of step)')") istep,tsec
      endif
#else
      call mk_polelat(0     ,1     ,itp)
      call mk_polelat(-1    ,2     ,itp)
      call mk_polelat(lat1+1,lat1  ,itp)
      call mk_polelat(lat1+2,lat1-1,itp)
#endif
!
! For Z, itc==itp (itp was set by addiag):
      z(:,:,:,itc) = z(:,:,:,itp)
!
! Horizontal diffusion, first step:
! hdif1 saves nrh and kmh at lats -2->nlat in nrh and kmh using
!   fields un_nm and vn_nm at j+1 and j+2, and tn_nm and barm at j+1.
! hdif1 calls mp_bnd[lat,lons]_kmh
!
      if (do_rtc) call timer(rtc_hdif1,tsec,'begin')
!
! Call method 1: pass starting address of first 3 dims of arrays,
!   latitude loop is inside the subroutine.
! hdif1 as an external subroutine: average time over 20 steps = .009 sec.
! This works with f4d data pointer only if hdif1 is not contained in a 
!   module.
!
      call hdif1(
     |  tn_nm(levd0,lond0,latd0,itp),
     |  un_nm(levd0,lond0,latd0,itp),
     |  vn_nm(levd0,lond0,latd0,itp),
     |  barm (levd0,lond0,latd0,itp),
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('advance after hdif1')")

      if (do_rtc) then
        call timer(rtc_hdif1,tsec,'end')
        write(6,"('Advance step ',i4,': rtc time for hdif1=',
     |    f10.2)") istep,tsec
      endif
!
! Horizontal diffusion, second step:
! hdif2 saves 3d diagnostics f3d (fkldt,u,v,o2,o) at j+1 (0->37)
! hdif2: WRITE (fkldt,u,v,o2,o1) <- READ (tn_nm,un_nm,vn_nm,o2_nm,o1_nm)
! The kld terms will be used later in hdif3.
!
      if (do_rtc) call timer(rtc_hdif2,tsec,'begin')
!
! Call method 1 (contiguous actual args, lat loop inside):
! Average over 20 steps = .014 sec. 
      call hdif2(
     |  tn_nm(levd0,lond0,latd0,itp),  ! 4d input
     |  un_nm(levd0,lond0,latd0,itp),  ! 4d input
     |  vn_nm(levd0,lond0,latd0,itp),  ! 4d input
     |  o2_nm(levd0,lond0,latd0,itp),  ! 4d input
     |  o1_nm(levd0,lond0,latd0,itp),  ! 4d input
     |  kldt ,      ! 3d output
     |  kldu ,      ! 3d output
     |  kldv ,      ! 3d output
     |  kldo2,      ! 3d output
     |  kldo1,      ! 3d output
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('advance after hdif2')")

      if (do_rtc) then
        call timer(rtc_hdif2,tsec,'end')
        write(6,"('Advance step ',i4,': rtc time for hdif2=',
     |    f10.2)") istep,tsec
      endif
!
! Dynamo calls Heelis (heelis.F), Weimer (wei01gcm.F), or neither
!   for high latitude electric potential, according to user-provided
!   "model_potential".
! Get high latitude (Heelis or other) colatitudes, NH pfrac, and poten phihm.
!  If Weimer is used, then theta0,phid etc is changed before use in aurora
!   in dynamics.
!
      if (potential_model == 'WEIMER') then
        call weimer01
        if (debug) write(6,"('advance after weimer01: istep=',i4)")istep
      elseif (potential_model == 'HEELIS') then
        call heelis
        if (debug) write(6,"('advance after heelis: istep=',i4)") istep
      else  !  potential_model='NONE'
        do j=1,nmlat0
          do i=1,nmlonp1
	    phihm(i,j) = 0.
          enddo ! i=1,nmlonp1
        enddo ! j=1,nmlat0
        call colath
      endif
!
! Model dynamics. Sub dynamics calls physics, chemistry, and dynamics
! routines for current time step:
!
      if (do_rtc) call timer(rtc_dynamics,tsec,'begin')
      call dynamics(istep,do_rtc)
      if (debug) write(6,"('advance after dynamics: istep=',i4)") istep
      if (do_rtc) then
        call timer(rtc_dynamics,tsec,'end')
        write(6,"('Advance step ',i4,': rtc time for dynamics=',
     |    f10.2)") istep,tsec
      endif
!
! Save Cusp, drizzle, alfa and nflux to secondary history to plot later
!     call addfsech_ij ('CUSP','2D CUSP','0-1',cusp2d,1,nlonp4,
!    |  1,nlat)
!     call addfsech_ij ('DRIZZLE','2D DRIZZLE','0-1',drzl2d,1,
!    |  nlonp4,1,nlat)
!     call addfsech_ij ('ALFA','2D ALFA','keV',alfa2d,1,
!    |  nlonp4,1,nlat)
!     call addfsech_ij ('NFLUX','2D NFLUX','#/cm3?',nflx2d,1,
!    |  nlonp4,1,nlat)
!
! Call dynamo if namelist input flag is set:
!
! zpotenm3d is 3-d geopotential height in mag coords. This field is forced
! on the secondary histories. Here init to zero:
!
      zpotenm3d = 0. ! whole array init (this field forced onto secondary hist)
      if (do_rtc.and.idynamo > 0) call timer(rtc_dynamo,tsec,'begin')
!
! idynamo==0 -> Do not call dynamo. However, if potential_model is
!               requested, set electric potential accordingly.
      if (idynamo <= 0) then
        if (trim(potential_model) /= 'NONE') then
          call prep_dynamo_dyn0(
     |      z      (levd0,lond0,latd0,itp),
     |      1,nlevp1,lon0,lon1,lat0,lat1)
            if (debug) write(6,"('advance after prep_dynamo_dyn0')")
          if (mytid==0) then 
            call nodynamo
            if (debug) write(6,"('advance after nodynamo')")
          endif
        endif   ! potential_model /= NONE
!
! idynamo==2 -> All tasks call prep_dynamo, and master task calls dynamo:
!
      else
        call prep_dynamo(
     |    tn     (levd0,lond0,latd0,itp),
     |    un     (levd0,lond0,latd0,itp),
     |    vn     (levd0,lond0,latd0,itp),
     |    w      (levd0,lond0,latd0,itp),
     |    z      (levd0,lond0,latd0,itp),
     |    barm   (levd0,lond0,latd0,itp),
     |    ped    (levd0,lond0,latd0),     ! from lamdas.F
     |    hall   (levd0,lond0,latd0),     ! from lamdas.F
     |    1,nlevp1,lon0,lon1,lat0,lat1)
          if (debug) write(6,"('advance after prep_dynamo')")
        if (mytid==0) call dynamo
        if (debug) write(6,"('advance after dynamo')")
      endif
      if (do_rtc.and.idynamo > 0) call timer(rtc_dynamo,tsec,'end')
!
! Update time-step timing before any history i/o:
      call end_timing(icount_step,elapsed_step,1)
      elapsed_steps = elapsed_steps+elapsed_step
      if (do_rtc) then
        call timer(rtc_step,tsec,'end')
        rtc_elapstep = tsec
        rtc_elapsteps = rtc_elapsteps+tsec
      endif
      rtc_elapstep = tsec
      rtc_elapsteps = rtc_elapsteps+tsec
#ifdef VT
!       code = 127 ; state = 'timestep' ; activity='ModelCode'
        call vtend(127,ier)
#endif
!
! Update slave tasks with new phim3d from the dynamo. Make this call
! after timing, so that timing will show the difference between master
! and slave tasks as a consequence of the serial dynamo. (slaves will
! be blocked here until the master is done with the dynamo)
!
#ifdef MPI
      if (idynamo > 0) call mp_updatephi
#endif
!
! Write output history if necessary (output writes data%(:,:,:,itc)).
!
      time2write = wrhist(istep,modeltime,
     |  wrprim, save_prim, newseries_prim, iseries_prim,
     |  wrsech, save_sech, newseries_sech, iseries_sech)
!
! If its time to write a history, root task must gather subdomain
! data from slave tasks before writing the history. This is done
! by sub mp_gather2root. Only root task writes to the history.
!
#ifdef MPI
      if (time2write.and.ntask > 1) then
        call start_timing(icount_mpi,'mp_gather2root',1)
        if (do_rtc) call timer(rtc_gather2root,tsec,'begin')
        call mp_gather2root(itc,'prim')
        if (debug) write(6,"('advance after mp_gather2root')")
        if (do_rtc) then
          call timer(rtc_gather2root,tsec,'end')
          write(6,"('Advance step ',i4,': rtc time for gather2root=',
     |      f10.2)") istep,tsec
        endif
!
! Comment this mp_gather2root call if sech fields were written by 
!   the master task only (e.g. from serial dynamo).
! This call can also be commented out if only progostics are saved
!   on secondary histories (i.e., addfsech was not called)
!
        if (wrsech) call mp_gather2root(itc,'sech')
        call end_timing(icount_mpi,elapsed_gather2root,1)
        elapsed_mpi = elapsed_mpi+elapsed_gather2root
      endif
#else
!
! Update foutput if serial non-MPI run:
!
      if (time2write) then
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,itc)
        enddo
      endif
#endif
!
! Write the history to output file:
      if (mytid==0) then
        call outhist(istep,modeltime)
        if (debug) write(6,"('advance after outhist')")
      endif
!
! Each mpi task must receive the 2 lats before its first
! updated lat (lat0-1,lat0-2), and the 2 lats after its last
! updated lat (lat1+1,lat2+2). Each task must also send its
! first 2 (lat0,lat0+1) and last 2 (lat1,lat1-1) updated lats
! to the appropriate "adjacent" tasks.
!
#ifdef MPI
      f4d(:)%mpi = .true. 
      if (do_rtc) call timer(rtc_bndlats,tsec,'begin')

      call mp_bndlats(f4d,nf4d,itc)

      if (debug) write(6,"('advance after mp_bndlats')")
      if (do_rtc) then
        call timer(rtc_bndlats,tsec,'end')
        write(6,"('Advance step ',i4,': rtc time for mp_bndlats=',
     |    f10.2,' (end of step)')") istep,tsec
      endif
!
      if (do_rtc) call timer(rtc_bndlons,tsec,'begin')

      call mp_bndlons(f4d,nf4d,itc)

      if (debug) write(6,"('advance after mp_bndlons')")
      if (do_rtc) then
        call timer(rtc_bndlons,tsec,'end')
        write(6,"('Advance step ',i4,': rtc time for mp_bndlons=',
     |    f10.2,' (end of step)')") istep,tsec
      endif
!
! Periodic points for all updated fields:
      call mp_periodic_f4d(itc)
#else
! Non-mpi serial run:
      call mk_polelat(0     ,1     ,itc)
      call mk_polelat(-1    ,2     ,itc)
      call mk_polelat(lat1+1,lat1  ,itc)
      call mk_polelat(lat1+2,lat1-1,itc)
      call set_periodic_f4d(itc)
#endif
!
! Swap field data time indices, so current updated data becomes previous
! step data for next step:
      itmp = itp
      itp = itc
      itc = itmp
!
! Update time-step timing:
!     call end_timing(icount_step,elapsed_step,1)
!     elapsed_steps = elapsed_steps+elapsed_step
!     call timer(rtc_step,tsec,'end')
!     rtc_elapstep = tsec
!     rtc_elapsteps = rtc_elapsteps+tsec
!
! Return for next time step:
      if (istep < nstep) then
#ifdef MPI
        if (ntask > 1) then
          call mpi_barrier(MPI_COMM_WORLD,ier)
          if (debug) write(6,"('advance after end of step barrier')")
        endif
#endif
        if (debug) write(6,"('advance end iteration for step ',i4)")
     |    istep
        goto 100
      endif
      avtime = avtime/nstep
!
! Report timing stats from mpi_module:
!
      if (do_rtc_mpi) then
        write(6,"(/,'Total rtc secs in mp_gather2root     = ',f8.2)") 
     |    rtcmp_gather2root
        write(6,"('Total rtc secs in mp_bndlats         = ',f8.2)") 
     |    rtcmp_bndlats
        write(6,"('Total rtc secs in mp_bndlons         = ',f8.2)") 
     |    rtcmp_bndlons
        write(6,"('Total rtc secs in mp_bndlons_f3d     = ',f8.2)") 
     |    rtcmp_bndlons_f3d
        write(6,"('Total rtc secs in mp_polelat         = ',f8.2)") 
     |    rtcmp_polelat
        write(6,"('Total rtc secs in mp_gatherlons      = ',f8.2)") 
     |    rtcmp_gatherlons
        write(6,"('Total rtc secs in mp_gatherlons_f3d  = ',f8.2)") 
     |    rtcmp_gatherlons_f3d
        write(6,"('Total rtc secs in mp_scatterlons     = ',f8.2)") 
     |    rtcmp_scatterlons
        write(6,"('Total rtc secs in mp_scatterlons_f3d = ',f8.2)") 
     |    rtcmp_scatterlons_f3d
        write(6,"('Total rtc secs in mp_periodic_f2d    = ',f8.2)")
     |    rtcmp_periodic_f2d
        write(6,"('Total rtc secs in mp_periodic_f3d    = ',f8.2)") 
     |    rtcmp_periodic_f3d
        write(6,"('Total rtc secs in mp_periodic_f4d    = ',f8.2)") 
     |    rtcmp_periodic_f4d
        write(6,"('Total rtc secs in mp_dynpot          = ',f8.2)") 
     |    rtcmp_dynpot
      endif
!
      if (do_rtc) then
        write(6,"(/,'Total rtc elapsed time for ',i4,' time ',
     |    ' steps = ',f12.3,' (seconds)')") nstep,rtc_elapsteps
        write(6,"('Average rtc secs per step = ',f8.3)")
     |    rtc_elapsteps/float(nstep)
      endif
#ifdef VT
!     code = 115 ; state = 'advance' ; activity='ModelCode'
      call vtend(115,ier)
#endif
!
      end subroutine advance
!-----------------------------------------------------------------------
      subroutine advance_day
!
! Advance calendar day if needed. Also update sfeps.
!
      use init_module,only: iter,iyear,iday,sfeps,sundec,sin_sundec,
     |  cos_sundec
      use cons_module,only: dt,pi
      implicit none
!
! Local:
      integer :: idayit,idayprev,iyearprev,iyr4,iyr100,lpyr,ienda
      real :: delta,theta0
!
      idayit = iter*int(dt)/86400
      if (idayit*86400 == iter*int(dt)) then
        idayprev = iday
        iyearprev = iyear
        iday = iday + 1
!
! lpyr = 1(0) if is (not) a leap year
        iyr4 = iyear/4
        iyr100 = iyear/100
        lpyr = 0
        if (iyr4*4 == iyear .and. iyr100*100 /= iyear) lpyr=1
        ienda = 365 + lpyr
        if (iday > ienda) then
          iyear = iyear + 1
          iday = iday - ienda
        endif                          !   for past year's end
!
! Recalculate sun's declination
        delta = atan(tan(23.5*pi/180.)*sin(2.*pi*float(iday-80)/365.))
        sin_sundec = sin(delta) ! C(95)
        cos_sundec = cos(delta) ! C(96)
!
! Update sfeps:
! sfeps is 6% variation in solar output over a year
!   caused by the orbital eccentricity.
!
        theta0 = 2.*pi*float(iday)/365.
        sfeps = 1.000110+0.034221*cos(theta0)+0.001280*sin(theta0)
     1        +0.000719*cos(2.*theta0)+0.000077*sin(2.*theta0)
!
        write(6,"('Advancing day (previous,present)=',4i5,' sfeps=',
     |   e12.4)") idayprev,iyearprev,iday,iyear,sfeps
      endif
      end subroutine advance_day
!-----------------------------------------------------------------------
      subroutine set_index(rindex,ntimes,msecs,outindex,name)
!
! User has provided time-dependent ctpoten (ctpoten_time) and/or
!   power (power_time) via namelist input. This routine interpolates
!   these inputs to the current model time msecs, returning outindex.
! This routine is called separately (from advance) for ctpoten and power,
!   (i.e., rindex will be either ctpoten_time or power_time from input).
! Note times in seconds are 8-byte integers.
!
      use params_module,only: 
     |  mxind_time ! max number of time-dependent solar index points
      implicit none
!
! Args:
      real,intent(in) :: 
     |  rindex(4,mxind_time) ! user input times and values (day,hour,min,value)
      integer,intent(in) :: 
     |  ntimes  ! number of valid time/values in rindex(:,1:ntimes)
      integer(kind=8),intent(in) ::
     |  msecs   ! current model time in seconds
      real,intent(out) :: outindex ! output interpolated value
      character(len=*),intent(in) :: name
!
! Local:
      integer :: i
      integer(kind=8) :: nsec0,nsec1
!
! External:
      integer(kind=8),external :: mtime_to_nsec
      real,external :: finterp_bigints
!
! If model time is beyond last rindex time, use last rindex time:
      nsec1 = mtime_to_nsec(int(rindex(1:3,ntimes)))
      if (msecs > nsec1) then
        outindex = rindex(4,ntimes)
        goto 100
      endif 
!
! Bracket model time:
      do i=1,ntimes-1
        nsec0 = mtime_to_nsec(int(rindex(1:3,i)))
        nsec1 = mtime_to_nsec(int(rindex(1:3,i+1)))
!
! If model time is at a provided time, interpolation is not necessary:
        if (nsec0 == msecs) then
          outindex = rindex(4,i)
          goto 100
        endif
        if (nsec1 == msecs) then
          outindex = rindex(4,i+1)
          goto 100
        endif
!
! Interpolate to model time msecs:
        if (msecs >= nsec0 .and. msecs <= nsec1) then
          outindex = finterp_bigints(rindex(4,i),rindex(4,i+1),nsec0,
     |      nsec1,msecs)
          goto 100
        endif
      enddo ! i=1,ntimes-1
!
! Error if model time could not be bracketed. This should not happen,
! but you never know...
      write(6,"(/,'>>> set_index: could not bracket model time ',
     |  i10)") msecs
      write(6,"('ntimes=',i3)") ntimes
      do i=1,ntimes
        write(6,"('  i=',i3,' ntimes=',i3,' day,hr,min=',3f7.2,
     |    ' value=',e12.4)") i,ntimes,rindex(1:3,i),rindex(4,i)
      enddo
      call shutdown('set_index')
!
! Report to stdout:
 100  continue
!     write(6,"('set_index: ',a,' = ',e12.4)") name,outindex
      end subroutine set_index
