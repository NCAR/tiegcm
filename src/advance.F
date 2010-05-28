      module advance_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      implicit none
!
! Advance the model in time.
!
      contains
!-----------------------------------------------------------------------
      subroutine advance
!
! Advance the model nstep time steps.
!
      use fields_module
      use hist_module,only: nstep,modeltime,nsource
      use heelis_module,only: heelis,colath
      use weimer_module,only: weimer01  ! Weimer 2001 model
      use wei05sc,only: weimer05,       ! Weimer 2005 model
     |  wei05sc_fac  ! field-aligned current diagnostic (nmlonp1,nmlat)
      use input_module,only: start,step,idynamo=>dynamo,f107,f107a,
     |  calendar_advance,gpi_ncfile,potential_model,wei05sc_ncfile,
     |  imf_ncfile
      use input_module,only: 
     |  power  ,power_time  ,ntimes_power, rd_power,
     |  ctpoten,ctpoten_time,ntimes_ctpoten, rd_ctpoten,
     |  bximf  ,bximf_time  ,ntimes_bximf,
     |  byimf  ,byimf_time  ,ntimes_byimf,
     |  bzimf  ,bzimf_time  ,ntimes_bzimf,
     |  swden  ,swden_time  ,ntimes_swden,
     |  swvel  ,swvel_time  ,ntimes_swvel,
     |  al     ,al_time     ,ntimes_al,
     |  kp     ,kp_time     ,ntimes_kp
      use init_module,only: istep,uthr,iter,secs,iday,iyear
      use cons_module,only: dt,re,cs,racs,ylatm,rtd
      use filter_module,only: filter
      use qrj_module,only: init_sflux
      use magfield_module,only: sunloc_apex
      use efield_module,only: efield
      use aurora_module,only: aurora_cons,cusp2d,drzl2d,alfa2d,nflx2d,
     |  eflux2d
      use hdif_module,only: hdif1,hdif2
      use timing_module,only: timer,timing
      use dynamo_module,only: prep_dynamo,dynamo,
     |  phihm,nmlat0,zpotenm3d
      use gpi_module,only: rdgpi,getgpi
      use imf_module,only: getimf
      use gswm_module,only: getgswm
      use params_module,only: nlonp4,nlat,nlev
      use output,only: outhist
      use addfld_module,only: addfld
      use lbc,only: tuvz_lbc
#if defined(INTERCOMM) || defined(CISMAH)
      use cism_coupling_module,only: import, export, geng,gflx
#endif

!      
! current due to plasma pressure & gravity -> add to RHS of dynamo    
      use magpres_g_module, only: magpres_grav,j_pg
!
#ifdef MPI
      use mpi_module,only: 
     |  mytid,mp_gather2root,mp_bndlats,ntask,lat0,lat1,lon0,lon1,
     |  mp_bndlons,mp_gatherlons,mp_scatterlons,mp_polelats,
     |  mp_periodic_f4d,mp_updatephi,mp_updateemphi,mp_updateemlam,
     |  mp_updateemz
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
#if defined(INTERCOMM) || defined(CISMAH)
      integer :: cmodeltime(4) !CISM coupling time
#endif
      real :: fmin,fmax
      integer(kind=8) :: nsecs ! current model time in seconds,
                               ! including model day
      logical :: time2write,wrprim,wrsech,save_prim,save_sech,
     |  newseries_sech,newseries_prim,iseries_sech,iseries_prim
      real :: secs_per_step,dday,
     |  time0,          time1,
     |  time0_step,     time1_step,
     |  time0_dynamics, time1_dynamics,
     |  time0_dynamo,   time1_dynamo,
     |  time0_phist,    time1_phist,
     |  time0_shist,    time1_shist,
     |  time0_prep,     time1_prep,
     |  time0_init,     time1_init,
     |  time0_barrier,  time1_barrier
      character(len=3) :: timing_type
!
      logical :: debug=.false.
!
! External:
      logical,external :: wrhist,time2print
      real,external :: hp_from_bz_swvel,hp_from_kp,ctpoten_from_kp
!
#ifdef VT
!     code = 115 ; state = 'advance' ; activity='ModelCode'
      call vtbegin(115,ier)
#endif
!
      write(6,"(/,'Enter advance: iter=',i6,' nstep=',i6)") iter,nstep
      call timer(time0_init,time1_init,'INIT',0,0)
!
! Init timestep index, modeltime, and nsecs:
      istep = 0
      modeltime(1:3) = start(:,1)
      modeltime(4) = 0
      call modeltime_to_nsecs(modeltime,nsecs) ! sets nsecs, includes day
      timing_type = 'sys'
      if (timing%rtc) timing_type = 'rtc'
      if (debug) write(6,"('advance after update modeltime..')")
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
        if (debug) 
     |    write(6,"('advance calling mp_gather2root for prim..')")
        call mp_gather2root(itc,'prim')
        if (debug) write(6,"('advance calling outhist..')")
        if (mytid==0) call outhist(istep,modeltime)
        if (debug) write(6,"('advance after outhist..')")
#else
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,itc)
        enddo

        tlbc_glb(lon0:lon1,lat0:lat1) = tlbc(lon0:lon1,lat0:lat1)
        ulbc_glb(lon0:lon1,lat0:lat1) = ulbc(lon0:lon1,lat0:lat1)
        vlbc_glb(lon0:lon1,lat0:lat1) = vlbc(lon0:lon1,lat0:lat1)

        if (debug) write(6,"('advance calling outhist..')")
        call outhist(istep,modeltime)
        if (debug) write(6,"('advance after outhist..')")
#endif
!
! nstep may be zero if user only wants to copy source history:
      if (nstep==0) then
        write(6,"(/,'ADVANCE: model is not advanced in time because ',
     |    ' start==stop.')")
        call timer(time0_init,time1_init,'INIT',1,0) ! end init timing
        return
      endif
!
! Reswap indices, so model reads itp data and updates itc data:
      itmp = itp
      itp = itc
      itc = itmp
!
      call timer(time0_init,time1_init,'INIT',1,0) ! end init timing
      if (timing%level >= 2)
     |  write(6,"('Time in INIT =          ',
     |    f12.2,' Advance: step ',i5)") time1_init,istep
!
! Main time loop:
  100 continue
#ifdef VT
!     code = 127 ; state = 'timestep' ; activity='ModelCode'
      call vtbegin(127,ier)
#endif
      call timer(time0_step,time1_step,'STEP',0,0) ! start step timing
      call timer(time0_prep,time1_prep,'PREP',0,0) ! start prep timing
      iter=iter+1
      istep = istep+1
      nsecs = nsecs+step
#if defined(INTERCOMM) || defined(CISMAH)
      cmodeltime=modeltime             !!! CISM ad hoc coupling
#endif
      call nsecs_to_modeltime(nsecs,modeltime) ! increments modeltime(4)
      secs = mod(iter*int(dt),86400)
      uthr = secs/3600.
      if (debug) write(6,"('advance after update time params..')")
!
! Update calendar day and sun's location and declination:
! (sub advance_day also updates sfeps)
!
      if (calendar_advance > 0) call advance_day
      call sunloc_apex(iyear,iday,secs)! locate sun's longitude
      if (debug) write(6,"('advance after sunloc')")
      if (calendar_advance > 0) dday = float(iday)+secs/86400.
!
! Interpolate input parameters to current model time, if time-dependent 
! values were read from input. If namelist read parameter indices_interp 
! is zero, interpolation is not done.
!
      if (ntimes_ctpoten > 0) 
     |  call set_index(ctpoten_time,ntimes_ctpoten,nsecs,ctpoten,
     |    'ctpoten')
      if (ntimes_power > 0) 
     |  call set_index(power_time,ntimes_power,nsecs,power,'power')
      if (ntimes_bximf > 0) 
     |  call set_index(bximf_time,ntimes_bximf,nsecs,bximf,'bximf')
      if (ntimes_byimf > 0) 
     |  call set_index(byimf_time,ntimes_byimf,nsecs,byimf,'byimf')
      if (ntimes_bzimf > 0) 
     |  call set_index(bzimf_time,ntimes_bzimf,nsecs,bzimf,'bzimf')
      if (ntimes_swden > 0) 
     |  call set_index(swden_time,ntimes_swden,nsecs,swden,'swden')
      if (ntimes_swvel > 0) 
     |  call set_index(swvel_time,ntimes_swvel,nsecs,swvel,'swvel')
      if (ntimes_al > 0) 
     |  call set_index(al_time,ntimes_al,nsecs,al,'al')
      if (ntimes_kp > 0) 
     |  call set_index(kp_time,ntimes_kp,nsecs,kp,'kp')
!
! Get gpi data if necessary:
      if (len_trim(gpi_ncfile) > 0) then
        if (istep==1) call rdgpi ! read gpi data
        iprint = 0
        if (istep==1) iprint = 1
        call getgpi(iyear,iday,int(secs),iprint)
      endif
!
! Get imf data if necessary:
      if (len_trim(imf_ncfile) > 0) then
        iprint = 0
        if (istep==1) iprint = 1
        call getimf(iyear,iday,int(secs),istep,iprint)
      endif
!
! If the potential model is Weimer and power was not provided
! by namelist, then calculate power from bz and swvel.
      if (potential_model(1:6)=='WEIMER'.and.
     |  rd_power==spval) then
        power = hp_from_bz_swvel(bzimf,swvel)
      endif
!
! If Kp was provided, calculate power and/or ctpoten
! (see also sub check_solar in input.F)
      if (kp /= spval) then
        if (potential_model(1:6)/='WEIMER'.and.rd_power==spval) then
          power = hp_from_kp(kp)
!         write(6,"('Advance calculated power from Kp: ',
!    |      ' KP=',f5.2,' POWER=',f8.2)") kp,power
        endif
        if (rd_ctpoten==spval) then
          ctpoten = ctpoten_from_kp(kp)
!         write(6,"('Advance calculated ctpoten from Kp: ',
!    |      ' KP=',f5.2,' CTPOTEN=',f8.2)") kp,ctpoten
        endif
      endif
!
! Get GSWM tidal boundary data:
!
      call getgswm(istep,iyear,iday,secs)
!
! Update lower boundaries of t,u,v,z:
!
      call tuvz_lbc
!
#if defined(INTERCOMM) || defined(CISMAH)
!
! Receive aurora parameters from the M-I coupler for CISM
!
      if (mytid==0) call import
!
! Scatter data to other tasks if using MPI, periodic points already added
! in ci_receive. Only need to scatter auroral characteristic energy and
! flux, potential is used in dynamo which is a serial code
!
#ifdef MPI
      call cism_scatter(geng,lon0,lon1,lat0,lat1)
      call cism_scatter(gflx,lon0,lon1,lat0,lat1)
#endif
#endif

!
! Report to stdout:
      if (time2print(nstep,istep)) then
        if (istep > 1) then
          write(6,"('Step ',i6,' of ',i6,' mtime=',4i3,
     |      ' secs/step (',a,') =',f6.2)") istep,nstep,modeltime,
     |      timing_type,secs_per_step
        else
          write(6,"('Step ',i6,' of ',i6,' mtime=',3i3)")
     |      istep,nstep,modeltime(1:3)
        endif
        if (len_trim(gpi_ncfile) > 0) then
          write(6,"('GPI run: power=',f7.2,' ctpoten=',f7.2,' f107=',
     |      f7.2,' f107a=',f7.2)") power,ctpoten,f107,f107a
        endif
        if (len_trim(imf_ncfile) > 0) then
          write(6,"('IMF run: bx=',f7.2,' by=',f7.2,' bz=',f7.2,
     |      ' swvel=',f7.2,' swden=',f7.2,' ctpoten=',f7.2,
     |      ' power=',f7.2)") bximf,byimf,bzimf,swvel,swden,
     |      ctpoten,power
        endif
      endif ! time2print
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
!     subroutine addiag(tn,o2,o1,vn,vc,barm,xnmbar,xnmbari,xnmbarm,z,
!    |  zg,lon0,lon1,lev0,lev1,lat0,lat1)
!
      call addiag(
     |  tn	 (levd0,lond0,latd0,itp),   ! in
     |  o2	 (levd0,lond0,latd0,itp),   ! in
     |  o1	 (levd0,lond0,latd0,itp),   ! in
     |  vn	 (levd0,lond0,latd0,itp),   ! in
     |  vc	 (levd0,lond0,latd0,itp),   ! out
     |  barm   (levd0,lond0,latd0,itp),     ! out
     |  xnmbar (levd0,lond0,latd0),	    ! out
     |  xnmbari(levd0,lond0,latd0),	    ! out
     |  xnmbarm(levd0,lond0,latd0),	    ! out
     |  z      (levd0,lond0,latd0,itp),     ! out (note itp is output)
     |  zg     (levd0,lond0,latd0),         ! out (z with varying grav)
     |  lon0,lon1,1,nlevp1,lat0,lat1)
!   
      if (debug) write(6,"('advance after addiag')")
!
! Define boundary latitudes -1,0 across the south pole, and
!   nlat+1,nlat+2 across the north pole:
!
#ifdef MPI
      call mp_polelats(itp)
!
! Update ghost cell boundary lats and lons.
! Bndlats is called first to exchange lat0-1,2 and lat1+1,2 at lon0->lon1, 
!   then bndlons is called  to exchange lon0-1,2 and lon1+1,2 at lat0-2->lat1+2).
!
      f4d(:)%mpi = .true. 
      call mp_bndlats(f4d,nf4d,itp)
      call mp_bndlons(f4d,nf4d,itp)
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
      call hdif1(
     |  tn_nm(levd0,lond0,latd0,itp),
     |  un_nm(levd0,lond0,latd0,itp),
     |  vn_nm(levd0,lond0,latd0,itp),
     |  barm (levd0,lond0,latd0,itp),
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('advance after hdif1')")
!
! Horizontal diffusion, second step:
! hdif2 saves 3d diagnostics f3d (fkldt,u,v,o2,o) at j+1 (0->37)
! hdif2: WRITE (fkldt,u,v,o2,o1) <- READ (tn_nm,un_nm,vn_nm,o2_nm,o1_nm)
! The kld terms will be used later in hdif3.
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
!
! Dynamo calls Heelis (heelis.F), Weimer (wei01gcm.F), or neither
!   for high latitude electric potential, according to user-provided
!   "model_potential".
! Get high latitude (Heelis or Weimer) colatitudes, NH pfrac, and poten phihm.
!  If Weimer is used, then theta0,phid etc is changed before use in aurora
!   in dynamics.
!
      if (potential_model == 'WEIMER01') then
        if (byimf==0..and.bzimf==0.) then
          write(6,"(/,'>>> WARNING: byimf and bzimf are both zero',
     |      ' before call to weimer01.',/,'This will cause the ',
     |      'weimer model to produce NaNs, so I am setting bzimf ',
     |      'to .01')") 
          bzimf = .01
        endif
        call weimer01
        if (debug) write(6,"('advance after weimer01: istep=',i4)")istep

!       write(6,"('wei01: mtime=',i3,',',i2,',',i2,' by,bz=',2f6.2,
!    |    ' swvel=',f8.2,' swden=',f6.2,' al=',f8.2,' phihm min,max=',
!    |    2e12.4)") modeltime(1:3),byimf,bzimf,swvel,swden,al,
!    |    minval(phihm),maxval(phihm)

      elseif (potential_model == 'WEIMER05'.or.
     |        potential_model == 'WEIMER') then
!
! sub weimer05 (wei05sc.F), calculates mag electric potential in phihm.
! (if byimf==bzimf==0, then weimer05 will set bzimf = .01)
!
        call weimer05(byimf,bzimf,swvel,swden,wei05sc_ncfile,istep)
        if (debug) write(6,"('advance after weimer05: istep=',i4)")istep

!       write(6,"('wei05: mtime=',i3,',',i2,',',i2,' by,bz=',2f8.2,
!    |    ' swvel=',f8.2,' swden=',f6.2,' phihm min,max=',
!    |    2e12.4)") modeltime(1:3),byimf,bzimf,swvel,swden,
!    |    minval(phihm),maxval(phihm)

        call addfld ('W05_EPOT','Weimer05 Electric Potential',
     |    'V',phihm,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
        call addfld ('W05_FAC','Weimer05 Field-Aligned Current',
     |    ' ',wei05sc_fac,'mlon',1,nmlonp1,'mlat',1,nmlat,0)

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
      call timer(time0_prep,time1,'PREP',1,0) ! end step-prep timing
      time1_prep = time1
!
! Model dynamics. Sub dynamics calls physics, chemistry, and dynamics
! routines for current time step:
!
      call timer(time0_dynamics,time1_dynamics,'DYNAMICS',0,0)
      call dynamics(nstep,istep)
      call timer(time0_dynamics,time1_dynamics,'DYNAMICS',1,0)
      if (debug) write(6,"('advance after dynamics: istep=',i4)") istep
      if (timing%level >= 2.and.time2print(nstep,istep))
     |  write(6,"('Time in DYNAMICS =     ',f12.2,
     |    ' Advance: step ',i5)") time1_dynamics,istep
!
! Save Cusp, drizzle, alfa and nflux to secondary history to plot later
!     call addfld ('CUSP','2D CUSP','0-1',
!    |  cusp2d(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld ('DRIZZLE','2D DRIZZLE','0-1',
!    |  drzl2d(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld ('ALFA','2D ALFA','keV',
!    |  alfa2d(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld ('NFLUX','2D NFLUX','#/cm2-s',
!    |  nflx2d(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld ('EFLUX','2D EFLUX','mW/m^2',
!    |  eflux2d(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Call dynamo if namelist input flag is set:
      call timer(time0_dynamo,time1_dynamo,'DYNAMO',0,0)
!
! zpotenm3d is 3-d geopotential height in mag coords. This field is forced
! on the secondary histories. Here init to zero:
!
      zpotenm3d = 0. ! whole array init (this field forced onto secondary hist)
!
! Calculate addition to electrodynamo RHS (current due to plasma 
! pressure and gravity)
!
      if(j_pg) then	
        call magpres_grav (
     |    z	(levd0,lond0,latd0,itp),   ! geopotential input
     |    te	(levd0,lond0,latd0,itp),   ! electron temperature K
     |    ti	(levd0,lond0,latd0,itp),   ! ion temperature K
     |    ne	(levd0,lond0,latd0,itp),   ! electron density 1/cm^3
     |    op	(levd0,lond0,latd0,itp),   ! O+ 1/cm^3
     |    nplus (levd0,lond0,latd0),	   ! N+  1/cm^3
     |    n2p	(levd0,lond0,latd0),	   ! N2+ 1/cm^3
     |    nop	(levd0,lond0,latd0),	   ! NO+  1/cm^3
     |    o2p	(levd0,lond0,latd0,itp),   ! O2+  1/cm^3
     |    1,nlevp1,lon0,lon1,lat0,lat1)
          if (debug) write(6,"('advance after magpres_grav: istep=',
     |    i4)") istep
      endif 
!
! Prepare inputs to dynamo and gather to master task:
!  
      call prep_dynamo(
     |  tn     (levd0,lond0,latd0,itp),
     |  un     (levd0,lond0,latd0,itp),
     |  vn     (levd0,lond0,latd0,itp),
     |  w      (levd0,lond0,latd0,itp),
     |  z      (levd0,lond0,latd0,itp),
     |  barm   (levd0,lond0,latd0,itp),
     |  ped    (levd0,lond0,latd0),     ! from lamdas.F
     |  hall   (levd0,lond0,latd0),     ! from lamdas.F
     |  1,nlevp1,lon0,lon1,lat0,lat1)
!
! Only master task calls dynamo: 
!
      if (debug) write(6,"('advance after prep_dynamo')")
      if (mytid==0) call dynamo
      if (debug) write(6,"('advance after dynamo')")
#if defined(INTERCOMM) || defined(CISMAH)
!
! Send conductance and neutral wind current to M-I coupler for CISM
!
      if (mytid==0) call export(cmodeltime)
#endif
#ifdef VT
!       code = 127 ; state = 'timestep' ; activity='ModelCode'
        call vtend(127,ier)
#endif
!
! End timing for dynamo and step (before i/o):
      call timer(time0_dynamo,time1_dynamo,'DYNAMO',1,0)
      if (timing%level >= 2.and.time2print(nstep,istep))
     |  write(6,"('Time in DYNAMO =       ',
     |    f12.2,' Advance: step ',i5)") time1_dynamo,istep
!
! Update slave tasks with new phim3d from the dynamo. Make this call
! after timing, so that timing will show the difference between master
! and slave tasks as a consequence of the serial dynamo. (slaves will
! be blocked here until the master is done with the dynamo)
!
#ifdef MPI
      call mp_updatephi
      call mp_updateemphi  ! am 12/12/05 should be gathered if we stick to this
      call mp_updateemlam
      call mp_updateemz
#endif
!
! Write output history if necessary (output writes data%(:,:,:,itc)).
!
      time2write = wrhist(istep,modeltime,
     |  wrprim, save_prim, newseries_prim, iseries_prim,
     |  wrsech, save_sech, newseries_sech, iseries_sech)
!
! History i/o timing is not included in STEP segment:
      call timer(time0_step,time1,'STEP',1,0) ! suspend step timing for i/o
      time1_step = time1
!
      time1_phist = 0.
      time1_shist = 0.
!
! If its time to write a history, root task must gather subdomain
! data from slave tasks before writing the history. This is done
! by sub mp_gather2root. Only root task writes to the history.
!
#ifdef MPI
      if (time2write.and.ntask > 1) then
        if (wrprim) call timer(time0_phist,time1_phist,'PHIST',0,0)
        call mp_gather2root(itc,'prim')
        if (debug) write(6,"('advance after mp_gather2root')")
        if (wrsech) then
          call timer(time0_shist,time1_shist,'SHIST',0,0)
          call mp_gather2root(itc,'sech')
        endif
      elseif (time2write.and.ntask==1) then
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,itc)
        enddo
        tlbc_glb(lon0:lon1,lat0:lat1) = tlbc(lon0:lon1,lat0:lat1)
        ulbc_glb(lon0:lon1,lat0:lat1) = ulbc(lon0:lon1,lat0:lat1)
        vlbc_glb(lon0:lon1,lat0:lat1) = vlbc(lon0:lon1,lat0:lat1)
        tlbc_nm_glb(lon0:lon1,lat0:lat1) = tlbc_nm(lon0:lon1,lat0:lat1)
        ulbc_nm_glb(lon0:lon1,lat0:lat1) = ulbc_nm(lon0:lon1,lat0:lat1)
        vlbc_nm_glb(lon0:lon1,lat0:lat1) = vlbc_nm(lon0:lon1,lat0:lat1)
      endif
#else
      call timer(time0_phist,time1_phist,'PHIST',0,0) ! non-MPI
      call timer(time0_shist,time1_shist,'SHIST',0,0) ! non-MPI
!
! Update foutput if serial non-MPI run:
!
      if (time2write) then
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,itc)
        enddo
        tlbc_glb(lon0:lon1,lat0:lat1) = tlbc(lon0:lon1,lat0:lat1)
        ulbc_glb(lon0:lon1,lat0:lat1) = ulbc(lon0:lon1,lat0:lat1)
        vlbc_glb(lon0:lon1,lat0:lat1) = vlbc(lon0:lon1,lat0:lat1)
        tlbc_nm_glb(lon0:lon1,lat0:lat1) = tlbc_nm(lon0:lon1,lat0:lat1)
        ulbc_nm_glb(lon0:lon1,lat0:lat1) = ulbc_nm(lon0:lon1,lat0:lat1)
        vlbc_nm_glb(lon0:lon1,lat0:lat1) = vlbc_nm(lon0:lon1,lat0:lat1)
      endif
#endif
!
! Write the history to output file:
      if (mytid==0) then
        call outhist(istep,modeltime)
        if (debug) write(6,"('advance after outhist')")
      endif
#ifdef MPI
      call timer(time0_step,time1,'STEP',0,0) ! resume step timing
!
      if (wrprim) call timer(time0_phist,time1_phist,'PHIST',1,0) 
      if (wrsech) call timer(time0_shist,time1_shist,'SHIST',1,0)
      if (timing%level >= 2.and.time2print(nstep,istep)) then
        write(6,"('Time in PHIST =        ',f12.2,
     |    ' Advance: step ',i5)") time1_phist,istep
        write(6,"('Time in SHIST =        ',f12.2,
     |    ' Advance: step ',i5)") time1_shist,istep
      endif
!
! Each mpi task must receive the 2 lats before its first
! updated lat (lat0-1,lat0-2), and the 2 lats after its last
! updated lat (lat1+1,lat2+2). Each task must also send its
! first 2 (lat0,lat0+1) and last 2 (lat1,lat1-1) updated lats
! to the appropriate "adjacent" tasks.
!
      call timer(time0_prep,time1,'PREP',0,0) ! resume prep timing
      f4d(:)%mpi = .true. 
      call mp_bndlats(f4d,nf4d,itc)
      if (debug) write(6,"('advance after mp_bndlats')")
      call mp_bndlons(f4d,nf4d,itc)
      if (debug) write(6,"('advance after mp_bndlons')")
!
! Periodic points for all updated fields:
      call mp_periodic_f4d(itc)
#else
! Non-mpi serial run:
      call timer(time0_step,time1,'STEP',0,0) ! resume step timing (non-MPI)
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
      call timer(time0_prep,time1,'PREP',1,0) ! end prep timing
      time1_prep = time1_prep+time1
      call timer(time0_step,time1,'STEP',1,0) ! end step timing
      time1_step = time1_step+time1
      if (timing%level >= 2.and.time2print(nstep,istep)) then
        write(6,"('Time in PREP =         ',
     |    f12.2,' Advance: step ',i5)") time1_prep,istep
        write(6,"('Time in STEP =         ',
     |    f12.2,' Advance: step ',i5)") time1_step,istep
      endif
      secs_per_step = time1_step
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
!
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
      integer :: idayit,idayprev,iyearprev,iyr4,iyr100,iyr400,lpyr,ienda
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
        iyr400 = iyear/400
        lpyr = 0
        if((iyr4*4 == iyear .and. iyr100*100 /= iyear).or.
     |     (iyr400*400 == iyear)) lpyr=1
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
      use input_module,only: indices_interp
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
          if (indices_interp > 0) then
            outindex = finterp_bigints(rindex(4,i),rindex(4,i+1),nsec0,
     |        nsec1,msecs)
          else
            outindex = rindex(4,i)
          endif
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
!-----------------------------------------------------------------------
      end module advance_module
