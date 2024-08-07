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
      logical :: cleanup
!
      contains
!-----------------------------------------------------------------------
      subroutine advance
!
! Advance the model nstep time steps.
!
      use fields_module
      use hist_module,only: nstep,modeltime
      use heelis_module,only: heelis
      use wei05sc,only: weimer05        ! Weimer 2005 model
      use input_module,only: pristart,step,idynamo=>dynamo,
     |  calendar_advance,gpi_ncfile,potential_model,wei05sc_ncfile,
     |  imf_ncfile,current_pg,current_kq,iamie,amie_ibkg,
     |  start_day,nudge_lbc,nudge_f4d,nudge_alpha,
     |  doEclipse,eclipse_list
      use input_module,only: 
     |  power  ,power_time  ,ntimes_power, rd_power,
     |  ctpoten,ctpoten_time,ntimes_ctpoten, rd_ctpoten,
     |  bximf  ,bximf_time  ,ntimes_bximf,
     |  byimf  ,byimf_time  ,ntimes_byimf,
     |  bzimf  ,bzimf_time  ,ntimes_bzimf,
     |  swden  ,swden_time  ,ntimes_swden,
     |  swvel  ,swvel_time  ,ntimes_swvel,
     |  al     ,al_time     ,ntimes_al,
     |  kp     ,kp_time     ,ntimes_kp,
     |  f107   ,f107_time   ,ntimes_f107,
     |  f107a  ,f107a_time  ,ntimes_f107a
      use init_module,only: istep,uthr,iter,secs,iday,iyear
      use cons_module,only: dt
      use qrj_module,only: init_sflux
      use magfield_module,only: sunloc_apex
      use pdynamo_module,only: pefield
      use aurora_module,only: aurora_cons
      use hdif_module,only: hdif1,hdif2
      use current,only: noso_crrt,noso_crdens,nscrrt_glb=>nscrrt
      use timing_module,only: timer,timing
      use mpitime_module,only: mpi_timer
      use gpi_module,only: rdgpi,getgpi
      use imf_module,only: getimf
      use gswm_module,only: getgswm
      use params_module,only: nlonp4,nlat,
     |  nlon,nlonp2,dlon,glon,glon0
      use amie_module,only: getamie,tiepot
      use output,only: outhist,output_hist
      use addfld_module,only: addfld
      use lbc,only: tuvz_lbc,t_lbc,u_lbc,v_lbc,z_lbc
      use eclipse_module,only: read_eclipse_list,find_eclipse
      use nudge_module,only:
     |  nudge_check=>check,nudge_init=>init,
     |  nudge_update=>update,nudge_finalize=>finalize,
     |  nudge_nlb=>nlb,nudge_nf4d=>nf4d,
     |  nudge_time=>time,nudge_itime=>itime,
     |  ext_lbc=>lbc,ext_f4d=>f4d,lbc_flag,
     |  f4d_idx,wrap,maxlev,latbeg,latend,lonbeg,lonend,
     |  offbeg,offend,lonbeg1,lonend1,hori_weight,vert_weight
!
! Routines and timing for parallel dynamo (pdynamo.F):
      use pdynamo_module,only: dynamo_inputs,pdynamo,pefield,
     |  prepare_phig3d,nmlat0,phihm 
      use diags_module,only: mkdiag_JQR

      use mage_coupling_module,only: import_mage, export_mage,
     |  process_import

!      
! current due to plasma pressure & gravity -> add to RHS of dynamo    
      use magpres_g, only: magpres_grav
!
      use mpi_module,only: 
     |  mytid,mp_bndlats,ntask,lat0,lat1,lon0,lon1,
     |  mlat0,mlat1,mlon0,mlon1,mp_scatter_phim,
     |  mp_bndlons,mp_polelats,mp_periodic_f4d,TIEGCM_WORLD
      use getapex_module,only: apxparm
!
      implicit none
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j,itmp,ier,iprint,month,day,modelsec,ifld,k
      integer :: cmodeltime(4) !MAGE coupling time
      real :: fmin,fmax,delta,fac1,fac2
      integer(kind=8) :: nsecs ! current model time in seconds,
                               ! including model day
      logical :: time2write,wrprim,wrsech,newseries_sech,newseries_prim
      integer :: iseries_sech,iseries_prim
      real :: secs_per_step,dday,
     |  time0,          time1,
     |  time0_step,     time1_step,
     |  time0_dynamics, time1_dynamics,
     |  time0_pdynamo,  time1_pdynamo,
     |  time0_phist,    time1_phist,
     |  time0_shist,    time1_shist,
     |  time0_prep,     time1_prep,
     |  time0_init,     time1_init
      character(len=3) :: timing_type
      logical :: adv
      real,dimension(mlon0:mlon1,mlat0:mlat1) :: nscrrt
      real,dimension(:,:),allocatable :: wt,ext_fld,tmp_lbc
      real,dimension(:,:,:),allocatable :: tmp_f4d,nclbc
      real,dimension(:,:,:,:),allocatable :: ncf4d
!
      logical :: debug=.false.
      integer :: lev0=1,lev1=nlevp1
!
! External:
      logical,external :: wrhist,time2print
      real,external :: hp_from_bz_swvel,hp_from_kp,ctpoten_from_kp
      integer,external :: mtime_to_secs
!
#ifdef VT
!     code = 115 ; state = 'advance' ; activity='ModelCode'
      call vtbegin(115,ier)
#endif
!
! A hook for handling OS signals (SIGINT/SIGTERM)
      call register_signal_handler
!
      write(6,"(/,'Enter advance: iter=',i10,' nstep=',i10)") iter,nstep
      call timer(time0_init,time1_init,'INIT',0,0)
!
! Init timestep index, modeltime, and nsecs:
      istep = 0
      modeltime = pristart(:,1)
      call modeltime_to_nsecs(modeltime,nsecs) ! sets nsecs, includes day
      timing_type = 'sys'
      if (timing%rtc) timing_type = 'rtc'
      if (debug) write(6,"('advance after update modeltime..')")
!
! Source was read into data(:,:,:,itp), and output always writes from
! data(:,:,:,itc), so swap indices to echo itp source data:
      itmp = itp
      itp = itc
      itc = itmp
!
! Start phist timing (source hist may or may not be written):
      call timer(time0_phist,time1_phist,'PHIST',0,0) ! start phist
!
! Echo source history to primary history output if an initial run,
      time2write = wrhist(istep,modeltime,
     |  wrprim, newseries_prim, iseries_prim,
     |  wrsech, newseries_sech, iseries_sech)

      if (time2write) then
!
! If an initial run, echo source history to output before first time step. 
!
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,itc)
        enddo
        call outhist(istep,modeltime)
!
! nstep may be zero if user only wants to copy source history:
        if (nstep==0) then
          write(6,"(/,'ADVANCE: model is not advanced in time because ',
     |      ' start==stop.')")
          call timer(time0_init,time1_init,'INIT',1,0) ! end init timing
          return
        endif
      endif ! time2write
      call timer(time0_phist,time1,'PHIST',1,0) ! suspend phist timing
      time1_phist = time1
!
! Reswap indices, so model reads itp data and updates itc data:
      itmp = itp
      itp = itc
      itc = itmp
!
! Initialize eclipse list with external data
      if (doEclipse) call read_eclipse_list(eclipse_list)
!
! glon0 corresponds to lon0:lon1 (1,nlonp4)
      glon0(3:nlonp2) = glon
      glon0(1:2) = glon(1)-(/2,1/)*dlon
      glon0(nlonp2+1:nlonp4) = glon(nlon)+(/1,2/)*dlon

! this nudge_lbc is from nudge module, not from input module
      if (any(nudge_lbc) .or. any(nudge_f4d)) then
        if (mytid == 0) call nudge_check
        call nudge_init

        if (latbeg < latend) then
          allocate(wt(lon0:lon1,latbeg:latend))
          allocate(ext_fld(lon0:lon1,latbeg:latend))
          allocate(tmp_lbc(lon0:lon1,latbeg:latend))
          allocate(tmp_f4d(maxlev,lon0:lon1,latbeg:latend))
          allocate(nclbc(lon0:lon1,latbeg:latend,nudge_nlb))
          allocate(ncf4d(maxlev,lon0:lon1,latbeg:latend,nudge_nf4d))
        endif
      endif
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
      cmodeltime=modeltime             !!! MAGE coupling
      call nsecs_to_modeltime(nsecs,modeltime) ! increments modeltime(4)
      secs = mod(iter*int(dt),86400)
      uthr = secs/3600.
      if (debug) write(6,"('advance after update time params..')")
!
! Update calendar day and sun's location and declination:
! (sub advance_day also updates sfeps)
!
      if (calendar_advance > 0) call advance_day(adv)
      call sunloc_apex(iyear,iday,secs)! locate sun's longitude
      if (debug) write(6,"('advance after sunloc')")
      if (calendar_advance > 0) dday = real(iday)+secs/86400.
!
! Interpolate input parameters to current model time, if time-dependent 
! values were read from input. If namelist read parameter indices_interp 
! is zero, interpolation is not done.
!
      if (ntimes_ctpoten > 0) 
     |  call set_index(ctpoten_time,ntimes_ctpoten,nsecs,ctpoten)
      if (ntimes_power > 0) 
     |  call set_index(power_time,ntimes_power,nsecs,power)
      if (ntimes_bximf > 0) 
     |  call set_index(bximf_time,ntimes_bximf,nsecs,bximf)
      if (ntimes_byimf > 0) 
     |  call set_index(byimf_time,ntimes_byimf,nsecs,byimf)
      if (ntimes_bzimf > 0) 
     |  call set_index(bzimf_time,ntimes_bzimf,nsecs,bzimf)
      if (ntimes_swden > 0) 
     |  call set_index(swden_time,ntimes_swden,nsecs,swden)
      if (ntimes_swvel > 0) 
     |  call set_index(swvel_time,ntimes_swvel,nsecs,swvel)
      if (ntimes_al > 0) 
     |  call set_index(al_time,ntimes_al,nsecs,al)
      if (ntimes_kp > 0) 
     |  call set_index(kp_time,ntimes_kp,nsecs,kp)
      if (ntimes_f107 > 0) 
     |  call set_index(f107_time,ntimes_f107,nsecs,f107)
      if (ntimes_f107a > 0) 
     |  call set_index(f107a_time,ntimes_f107a,nsecs,f107a)
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
        call getimf(iyear,iday,int(secs),istep)
      endif
!
! If the potential model is Weimer and power was not provided
! by namelist, then calculate power from bz and swvel.
      if (potential_model(1:6)=='WEIMER'.and.
     |  rd_power==spval) then
        power = hp_from_bz_swvel(bzimf,swvel)
!       write(6,"('Advance calculated power from bz,swvel: ',
!    |    ' bz=',f8.2,' swvel=',f8.2,' power=',f8.2)")
!    |    bzimf,swvel,power
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
! Get AMIE data
      if (iamie==1) then
        iprint = 0
        if (istep==1) iprint = 1
        if (iprint>0) write(6,"('advance calling getamie...')")
        call getamie(iyear,iday,int(secs),amie_ibkg,iprint)
      endif
!
! Get GSWM tidal boundary data:
!
      call getgswm(istep,iday,secs)
!
! Update lower boundaries of t,u,v,z:
!
      call tuvz_lbc(istep,iday,secs)
!
! Update eclipse elements array if necessary:
      if (doEclipse) then
        call to_month_day(iyear,iday,month,day)
        call find_eclipse(iyear,month,day)
      endif
!
      modelsec = mtime_to_secs(modeltime)-start_day*86400

! read external fields and interpolate spatially
! this sets up both ext_lbc and ext_f4d
      if (any(nudge_lbc) .or. any(nudge_f4d)) then
        if (latbeg<latend .and. nudge_time(nudge_itime+1)<modelsec)
     |    call nudge_update(modelsec)
        delta = nudge_time(nudge_itime+1)-nudge_time(nudge_itime)
        fac1 = (nudge_time(nudge_itime+1)-modelsec)/delta
        fac2 = (modelsec-nudge_time(nudge_itime))/delta
      endif

! lower boundary based on nudge data
      if (any(nudge_lbc)) then
        if (latbeg<latend .and.
     |    nudge_time(nudge_itime)<=modelsec .and.
     |    modelsec<=nudge_time(nudge_itime+1)) then

! temporal interpolation
          nclbc = fac1*ext_lbc(:,:,1,:)+fac2*ext_lbc(:,:,2,:)
          wt = nudge_alpha*hori_weight

          do ifld = 1,nudge_nlb
            if (lbc_flag(ifld)) then
              select case (ifld)
                case (1)
                  tmp_lbc = t_lbc(lon0:lon1,latbeg:latend)
                case (2)
                  tmp_lbc = u_lbc(lon0:lon1,latbeg:latend)
                case (3)
                  tmp_lbc = v_lbc(lon0:lon1,latbeg:latend)
                case (4)
                  tmp_lbc = z_lbc(lon0:lon1,latbeg:latend)
              endselect

              ext_fld = wt*nclbc(:,:,ifld)+(1-wt)*tmp_lbc

! find the corresponding grid and replace the model fields
              if (wrap) then
                tmp_lbc = ext_fld
              else
                if (lonbeg<lonend .and. offbeg==offend)
     |            tmp_lbc(lonbeg:lonend,:) = ext_fld(lonbeg:lonend,:)

                if (lonbeg<lonend .and. offbeg+1==offend) then
                  tmp_lbc(lonbeg:lonend1,:) = ext_fld(lonbeg:lonend1,:)
                  tmp_lbc(lonbeg1:lonend,:) = ext_fld(lonbeg1:lonend,:)
                endif

                if (lonbeg>lonend .and. offbeg==offend+1) then
                  tmp_lbc(lonbeg:lon1,:) = ext_fld(lonbeg:lon1,:)
                  tmp_lbc(lon0:lonend,:) = ext_fld(lon0:lonend,:)
                endif
              endif

              select case (ifld)
                case (1)
                  t_lbc(lon0:lon1,latbeg:latend) = tmp_lbc
                case (2)
                  u_lbc(lon0:lon1,latbeg:latend) = tmp_lbc
                case (3)
                  v_lbc(lon0:lon1,latbeg:latend) = tmp_lbc
                case (4)
                  z_lbc(lon0:lon1,latbeg:latend) = tmp_lbc
              endselect
            endif
          enddo
        endif
      endif
!
! Receive aurora parameters from the M-I coupler for MAGE
!
      call import_mage
!
! Send conductance and neutral wind current to M-I coupler for MAGE
!
      call export_mage(cmodeltime)
!
! Process the imported data
!
      call process_import
!
! Report to stdout:
      if (time2print(nstep,istep)) then
        if (istep > 1) then
          write(6,"('Step ',i8,' of ',i8,' mtime=',4i3,
     |      ' secs/step (',a,') =',f6.2)") istep,nstep,modeltime,
     |      timing_type,secs_per_step
        else
          write(6,"('Step ',i8,' of ',i8,' mtime=',4i3)")
     |      istep,nstep,modeltime
        endif
!
! GPI data run:
        if (len_trim(gpi_ncfile) > 0) then
          write(6,"('GPI run: mtime=',4i3,' power=',f8.3,' ctpoten=',
     |      f8.3,' f107=',f8.3,' f107a=',f8.3)")modeltime,
     |      power,ctpoten,f107,f107a
        endif
!
! Weimer/IMF/OMNI data run:
        if (len_trim(imf_ncfile) > 0) then
          write(6,"('IMF run: mtime=',4i3,' bx=',f8.3,' by=',f8.3,
     |      ' bz=',f8.3,' swvel=',f8.3,' swden=',f8.3)")
     |      modeltime,bximf,byimf,bzimf,swvel,swden
        endif
      endif ! time2print
!
! Write output history if necessary (master task writes data%(:,:,:,itc)).
! This is called early in the timestep so time2write is available to other
! modules who may need it.
!
      time2write = wrhist(istep,modeltime,
     |  wrprim, newseries_prim, iseries_prim,
     |  wrsech, newseries_sech, iseries_sech)
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
! Calculate electric field on magnetic subdomains (first 
! timestep only), and regrid ! to geographic subdomains 
! in ex,ey,ez. (after first timestep, pthreed will calculate
! the electric field)
!
      call mpi_timer('pefield',0,0)
      if (debug) write(6,"('advance call pefield')")
      call pefield
      if (debug) write(6,"('advance after pefield')")
      call mpi_timer('pefield',1,0)
!
! Calculate vc, barm, xnmbar, and z (tn,o2,o1,vn are input).
!
      call addiag(
     |  tn     (levd0,lond0,latd0,itp), ! in
     |  o2     (levd0,lond0,latd0,itp), ! in
     |  o1     (levd0,lond0,latd0,itp), ! in
     |  he     (levd0,lond0,latd0,itp), ! in
     |  vn     (levd0,lond0,latd0,itp), ! in
     |  vc     (levd0,lond0,latd0,itp), ! out
     |  mbar   (levd0,lond0,latd0,itp), ! out
     |  barm   (levd0,lond0,latd0,itp), ! out
     |  xnmbar (levd0,lond0,latd0,itp), ! out
     |  xnmbari(levd0,lond0,latd0,itp), ! out
     |  scht   (levd0,lond0,latd0,itp), ! out
     |  schti  (levd0,lond0,latd0,itp), ! out
     |  z      (levd0,lond0,latd0,itp), ! out (note itp is output)
     |  zg     (levd0,lond0,latd0),     ! out (z with varying grav)
     |  n2     (levd0,lond0,latd0),     ! out (1.-o2-o1-he)
     |  lon0,lon1,1,nlevp1,lat0,lat1)
!   
      if (debug) write(6,"('advance after addiag')")
!
! Define boundary latitudes -1,0 across the south pole, and
!   nlat+1,nlat+2 across the north pole:
!
      if (ntask==1) then
        call mk_polelat( 0,1,itp)
        call mk_polelat(-1,2,itp)
        call mk_polelat(nlat+1,nlat  ,itp)
        call mk_polelat(nlat+2,nlat-1,itp)
      else
        call mp_polelats(itp)
        if (debug) write(6,"('advance after mp_polelats')")
!
! Update ghost cell boundary lats and lons.
! Bndlats is called first to exchange lat0-1,2 and lat1+1,2 at lon0->lon1, 
!   then bndlons is called  to exchange lon0-1,2 and lon1+1,2 at lat0-2->lat1+2).
!
        f4d(:)%mpi = .true. 
        call mp_bndlats(f4d,nf4d,itp)
        call mp_bndlons(f4d,nf4d,itp)
      endif ! ntask==1
!
! For Z, itc==itp (itp was set by addiag):
      z(:,:,:,itc) = z(:,:,:,itp)
      if (debug) write(6,"('advance after z')")
!
! Horizontal diffusion, first step:
! hdif1 saves nrh and kmh at lats -2->nlat in nrh and kmh using
!   fields un_nm and vn_nm at j+1 and j+2, and tn_nm and mbar at j+1.
! hdif1 calls mp_bnd[lats,lons]_f3d
!
      call hdif1(
     |  tn_nm(levd0,lond0,latd0,itp),
     |  un_nm(levd0,lond0,latd0,itp),
     |  vn_nm(levd0,lond0,latd0,itp),
     |  mbar (levd0,lond0,latd0,itp),
     |  fnrh (levd0,lond0,latd0),
     |  fkmh (levd0,lond0,latd0),
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('advance after hdif1')")
!
! Horizontal diffusion, second step:
! hdif2 saves 3d diagnostics f3d (fkldt,u,v,o2,o,he) at j+1 (0->37)
! hdif2: WRITE (fkldt,u,v,o2,o1,he) <- READ (tn_nm,un_nm,vn_nm,o2_nm,o1_nm,he_nm)
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
     |  he_nm(levd0,lond0,latd0,itp),  ! 4d input
     |  fnrh (levd0,lond0,latd0),      ! 3d input
     |  fkmh (levd0,lond0,latd0),      ! 3d input
     |  kldt (levd0,lond0,latd0),      ! 3d output
     |  kldu (levd0,lond0,latd0),      ! 3d output
     |  kldv (levd0,lond0,latd0),      ! 3d output
     |  kldo2(levd0,lond0,latd0),      ! 3d output
     |  kldo1(levd0,lond0,latd0),      ! 3d output
     |  kldhe(levd0,lond0,latd0),      ! 3d output
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('advance after hdif2')")
!
! Dynamo calls Heelis (heelis.F), Weimer (wei05sc.F), or neither
!   for high latitude electric potential, according to user-provided
!   "model_potential".
! Get high latitude (Heelis or Weimer) colatitudes, NH pfrac, and poten phihm.
!  If Weimer is used, then theta0,phid etc is changed before use in aurora
!   in dynamics.
!
      if (potential_model == 'WEIMER05'.or.
     |    potential_model == 'WEIMER') then
!
! sub weimer05 (wei05sc.F), calculates mag electric potential in phihm.
! (if byimf==bzimf==0, then weimer05 will set bzimf = .01)
! Weimer will calculate ctpoten (user cannot provide constant
!   namelist ctpoten with the Weimer model)
!
        call weimer05(byimf,bzimf,swvel,swden,wei05sc_ncfile,istep)
        if (debug) write(6,"('advance after weimer05: istep=',i4)")istep

!       call addfld ('W05_EPOT','Weimer05 Electric Potential',
!    |    'V',phihm(mlon0:mlon1,mlat0:mlat1),
!    |    'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)

      elseif (potential_model == 'HEELIS') then
        call heelis
        if (debug) write(6,"('advance after heelis: istep=',i4)") istep
      else  !  potential_model='NONE'
        do j=1,nmlat0
          do i=1,nmlonp1
	    phihm(i,j) = 0.
          enddo ! i=1,nmlonp1
        enddo ! j=1,nmlat0
      endif
!
! If an amie run, use tiepot from amie.F:
! (phihm and tiepot are dimensioned (nmlonp1,nmlat))
      if (iamie > 0) then
        phihm = tiepot ! whole array op
        call fminmax(phihm,nmlonp1*nmlat,fmin,fmax)
        write(6,"('Advance: amie phihm min,max=',2e12.4)") fmin,fmax
      endif
!
! Set pfrac for dynamo, regardless of potential_model (in colath.F): 
      call colath

      call timer(time0_prep,time1,'PREP',1,0) ! suspend step-prep timing
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
! Call dynamo if namelist idynamo > 0:
!
      if (idynamo > 0) then
        if (calendar_advance > 0) then
          if (adv) call apxparm(iyear+iday/366.)
        endif
!
! Calculate addition to electrodynamo RHS (current due to plasma 
! pressure and gravity)
!
        if (current_pg > 0) then	
          call timer(time0,time1,"MAGPRES_GRAV",0,0)
          call magpres_grav (
     |      z	(levd0,lond0,latd0,itp),   ! geopotential input
     |      te	(levd0,lond0,latd0,itp),   ! electron temperature K
     |      ti	(levd0,lond0,latd0,itp),   ! ion temperature K
     |      ne	(levd0,lond0,latd0,itp),   ! electron density 1/cm^3
     |      op	(levd0,lond0,latd0,itp),   ! O+ 1/cm^3
     |      nplus (levd0,lond0,latd0),	   ! N+  1/cm^3
     |      n2p	(levd0,lond0,latd0),	   ! N2+ 1/cm^3
     |      nop	(levd0,lond0,latd0),	   ! NO+  1/cm^3
     |      o2p	(levd0,lond0,latd0,itp),   ! O2+  1/cm^3
     |      1,nlevp1,lon0,lon1,lat0,lat1)
          if (debug) write(6,"('advance after magpres_grav: istep=',
     |      i4)") istep
          call timer(time0,time1,"MAGPRES_GRAV",1,0)
        endif ! current_pg > 0
!
! Prepare neutral inputs for parallel dynamo:
! Calculate vertical velocity wn:
        call calc_wn(
     |    w    (lev0:lev1,lon0:lon1,lat0:lat1,itp),
     |    scht (lev0:lev1,lon0:lon1,lat0:lat1,itp),
     |    wn   (lev0:lev1,lon0:lon1,lat0:lat1),     ! output
     |    1,nlevp1,lon0,lon1,lat0,lat1)

        call mpi_timer('dynamo_inputs',0,0)
        call dynamo_inputs(
     |    un   (lev0:lev1,lon0:lon1,lat0:lat1,itp),
     |    vn   (lev0:lev1,lon0:lon1,lat0:lat1,itp),
     |    wn   (lev0:lev1,lon0:lon1,lat0:lat1),
     |    z    (lev0:lev1,lon0:lon1,lat0:lat1,itp),
     |    ped  (lev0:lev1,lon0:lon1,lat0:lat1),     ! from lamdas.F
     |    hall (lev0:lev1,lon0:lon1,lat0:lat1),     ! from lamdas.F
     |    1,nlevp1,lon0,lon1,lat0,lat1)
        call mpi_timer('dynamo_inputs',1,0)
        if (debug) write(6,"('advance after dynamo_inputs: istep=',
     |    i5)") istep
!
! Call parallel dynamo (pdynamo.F has its own timing):
        call timer(time0_pdynamo,time1_pdynamo,'PDYNAMO',0,0)
        call pdynamo
        call timer(time0_pdynamo,time1_pdynamo,'PDYNAMO',1,0)
        if (debug) write(6,"('advance after pdynamo: istep=',i5)")
     |    istep
        if (timing%level >= 2.and.time2print(nstep,istep))
     |    write(6,"('Time in DYNAMO =       ',
     |      f12.2,' Advance: step ',i5)") time1_pdynamo,istep
!
! Calculate current diagnostics (current.F90):
! (note sub nosocoef is called by the root task from complete_integrals
!  in pdynamo. Sub nosocoef is a stub in current.F90 that calls sub noso_coef)
!
        if (current_kq > 0) then
          if (mytid==0) call noso_crrt
!
! Save JQR: Upward current density (2d)
!
          call mp_scatter_phim(nscrrt_glb,nscrrt)
          call mkdiag_JQR('JQR',nscrrt,mlon0,mlon1,mlat0,mlat1)
          call noso_crdens
        endif
!
! Prepare potential on geographic grid if writing history:
        call prepare_phig3d(wrprim,wrsech)
        if (debug) write(6,"('advance after prepare_phig3d')")
!
! If idynamo <= 0, dynamo is not called:
      else
        if (time2print(nstep,istep))
     |    write(6,"('Note: idynamo=',i3,' -> dynamo was not called')")
     |      idynamo
      endif ! idynamo
!
      if (any(nudge_f4d)) then
        if (latbeg<latend .and.
     |    nudge_time(nudge_itime)<=modelsec .and.
     |    modelsec<=nudge_time(nudge_itime+1)) then

          ncf4d = fac1*ext_f4d(:,:,:,1,:)+fac2*ext_f4d(:,:,:,2,:)

          do ifld = 1,nudge_nf4d
            if (nudge_f4d(ifld)) then
              tmp_f4d = f4d(f4d_idx(ifld))%
     |          data(1:maxlev,lon0:lon1,latbeg:latend,itc)

              do k = 1,maxlev
                wt = nudge_alpha*hori_weight*vert_weight(k,ifld)
                ext_fld = wt*ncf4d(k,:,:,ifld)+(1-wt)*tmp_f4d(k,:,:)

                if (wrap) then
                  tmp_f4d(k,:,:) = ext_fld
                else
                  if (lonbeg<lonend .and. offbeg==offend)
     |              tmp_f4d(k,lonbeg:lonend,:) =
     |              ext_fld(lonbeg:lonend,:)

                  if (lonbeg<lonend .and. offbeg+1==offend) then
                    tmp_f4d(k,lonbeg:lonend1,:) =
     |                ext_fld(lonbeg:lonend1,:)
                    tmp_f4d(k,lonbeg1:lonend,:) =
     |                ext_fld(lonbeg1:lonend,:)
                  endif

                  if (lonbeg>lonend .and. offbeg==offend+1) then
                    tmp_f4d(k,lonbeg:lon1,:) = ext_fld(lonbeg:lon1,:)
                    tmp_f4d(k,lon0:lonend,:) = ext_fld(lon0:lonend,:)
                  endif
                endif
              enddo

              f4d(f4d_idx(ifld))%
     |          data(1:maxlev,lon0:lon1,latbeg:latend,itc) = tmp_f4d
            endif
          enddo
        endif
      endif
!
#ifdef VT
!       code = 127 ; state = 'timestep' ; activity='ModelCode'
        call vtend(127,ier)
#endif
!
! Write output history if necessary (output writes data%(:,:,:,itc)).
!
      time2write = wrhist(istep,modeltime,
     |  wrprim, newseries_prim, iseries_prim,
     |  wrsech, newseries_sech, iseries_sech)
!
! History i/o timing is not included in STEP segment:
      call timer(time0_step,time1,'STEP',1,0) ! suspend step timing for i/o
      time1_step = time1
!
! If its time to write a history, all tasks write to the history.
!
      if (time2write) then
        if (wrprim) call timer(time0_phist,time1,'PHIST',0,0) ! resume phist
        if (wrsech) call timer(time0_shist,time1_shist,'SHIST',0,0)
!
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,itc)
        enddo
!
! Saving amieoutput here - G. Lu, Aug. 27, 2003
        if (iamie > 0) then
          call prep_output
          call amieoutput
        endif
!
! Write the history to output file:
        if (debug) write(6,"('advance call outhist')")
        call outhist(istep,modeltime)
        if (debug) write(6,"('advance after call outhist')")
!
! Update i/o timing:
        if (wrprim) then
          call timer(time0_phist,time1,'PHIST',1,0) ! end phist timing
          time1_phist = time1_phist+time1
        endif
        if (wrsech) call timer(time0_shist,time1_shist,'SHIST',1,0)
      endif ! time2write
      call timer(time0_step,time1,'STEP',0,0) ! resume step timing (non-MPI)
!
! Only Intel and GNU can capture SIGINT/SIGTERM and perform cleanup work
! Other compilers will simply skip this check
      if (cleanup) then
! A new volume will always be created for primary histories (restart)
        call output_hist(istep,.true.)
        call shutdown('Received termination signal, exiting ...')
      endif
!
! Periodic points for all updated fields:
      call mp_periodic_f4d(itc)
      if (debug) write(6,"('advance after mp_periodic_f4d')")
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
        if (ntask > 1) then
          call mpi_barrier(TIEGCM_WORLD,ier) ! is this necessary??
          if (debug) write(6,"('advance after end of step barrier')")
        endif
!
        if (debug) write(6,"('advance end iteration for step ',i4)")
     |    istep
        goto 100
      endif
!
      if (any(nudge_lbc) .or. any(nudge_f4d)) then
        if (latbeg < latend) call nudge_finalize
      endif
!
#ifdef VT
!     code = 115 ; state = 'advance' ; activity='ModelCode'
      call vtend(115,ier)
#endif
!
      end subroutine advance
!-----------------------------------------------------------------------
      subroutine advance_day(adv)
!
! Advance calendar day if needed. Also update sfeps.
!
      use init_module,only: iter,iyear,iday,sfeps,sin_sundec,cos_sundec
      use cons_module,only: dt,pi
      implicit none
!
      logical,intent(out) :: adv
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
        delta = atan(tan(23.5*pi/180.)*sin(2.*pi*real(iday-80)/365.))
        sin_sundec = sin(delta) ! C(95)
        cos_sundec = cos(delta) ! C(96)
!
! Update sfeps:
! sfeps is 6% variation in solar output over a year
!   caused by the orbital eccentricity.
!
        theta0 = 2.*pi*real(iday)/365.
        sfeps = 1.000110+0.034221*cos(theta0)+0.001280*sin(theta0)
     1        +0.000719*cos(2.*theta0)+0.000077*sin(2.*theta0)
!
        write(6,"('Advancing day (previous,present)=',4i5,' sfeps=',
     |   e12.4)") idayprev,iyearprev,iday,iyear,sfeps
        adv = .true.
      else
        adv = .false.
      endif
      end subroutine advance_day
!-----------------------------------------------------------------------
      subroutine set_index(rindex,ntimes,msecs,outindex)
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
     |  rindex(5,mxind_time) ! user input times and values (day,hour,min,sec,value)
      integer,intent(in) :: 
     |  ntimes  ! number of valid time/values in rindex(:,1:ntimes)
      integer(kind=8),intent(in) ::
     |  msecs   ! current model time in seconds
      real,intent(out) :: outindex ! output interpolated value
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
      nsec1 = mtime_to_nsec(int(rindex(1:4,ntimes)))
      if (msecs > nsec1) then
        outindex = rindex(5,ntimes)
        goto 100
      endif 
!
! Bracket model time:
      do i=1,ntimes-1
        nsec0 = mtime_to_nsec(int(rindex(1:4,i)))
        nsec1 = mtime_to_nsec(int(rindex(1:4,i+1)))
!
! If model time is at a provided time, interpolation is not necessary:
        if (nsec0 == msecs) then
          outindex = rindex(5,i)
          goto 100
        endif
        if (nsec1 == msecs) then
          outindex = rindex(5,i+1)
          goto 100
        endif
!
! Interpolate to model time msecs:
        if (msecs >= nsec0 .and. msecs <= nsec1) then
          if (indices_interp > 0) then
            outindex = finterp_bigints(rindex(5,i),rindex(5,i+1),nsec0,
     |        nsec1,msecs)
          else
            outindex = rindex(5,i)
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
        write(6,"('  i=',i3,' ntimes=',i3,' day,hr,min,sec=',4f7.2,
     |    ' value=',e12.4)") i,ntimes,rindex(1:4,i),rindex(5,i)
      enddo
      call shutdown('set_index')
!
! Report to stdout:
 100  continue
!     write(6,"('set_index: ',e12.4)") outindex
      end subroutine set_index
!-----------------------------------------------------------------------
      subroutine calc_wn(omega,scht,wn,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Calculate vertical velocity wn for dynamo inputs:
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) ::
     |  omega,scht
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: wn
!
! Local:
      integer :: j,i,k

!     write(6,"('calc_wn: omega=',2e12.4,' scht=',2e12.4)")
!    |  minval(omega),maxval(omega),minval(scht),maxval(scht)
!
      do j=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            wn(k,i,j)=0.5*(omega(k,i,j)+omega(k+1,i,j))*scht(k,i,j)
          enddo
          wn(lev1,i,j) = 0.
        enddo
      enddo
      end subroutine calc_wn
!-----------------------------------------------------------------------
      subroutine nodynamo
!
! Dynamo was not called: zero out electric field and potential:
!
      use fields_module,only: 
     |  phim3d,  ! (nmlonp1,nmlat,-2:nlevp1) ! 3d electric potential magnetic
     |  emphi3d, ! (nmlonp1,nmlat,-2:nlevp1) ! 3d eastward electric field magnetic
     |  emlam3d, ! (nmlonp1,nmlat,-2:nlevp1) ! 3d equatorw. electric field magnetic
     |  emz3d    ! (nmlonp1,nmlat,-2:nlevp1) ! 3d upward (?) electric field magnetic

      phim3d = 0.
      emphi3d = 0.
      emlam3d = 0.
      emz3d = 0.

!     write(6,"('nodynamo: zeroed out phim3d, emphi3d,',
!    |  ' emlam3d, and emz3d')")

      end subroutine nodynamo
!-----------------------------------------------------------------------
      subroutine register_signal_handler
! Register signal handler to OS to catch SIGINT/SIGTERM
! Intel and GNU have different forms of function calls
! Other compilers are not supported yet

#ifdef INTEL
      use ifport,only: signal,sigint,sigterm
#endif

#ifdef GNU
! Defined by POSIX, hard code here
      integer,parameter :: sigint=2,sigterm=15
#endif

      integer :: ret

#ifdef INTEL
      ret = signal(sigint,signal_handler,-1)
      ret = signal(sigterm,signal_handler,-1)
#endif

#ifdef GNU
      ret = signal(sigint,signal_handler)
      ret = signal(sigterm,signal_handler)
#endif

      cleanup = .false.

      end subroutine register_signal_handler
!-----------------------------------------------------------------------
      function signal_handler(signum) result(ret)
! This handler does only one thing: change the cleanup flag to true
! Note that only with Intel and GNU this function will be called
! Otherwise the cleanup flag stays false all the time

      integer,intent(in) :: signum
      integer :: ret

      cleanup = .true.
      ret = 0

      end function signal_handler
!-----------------------------------------------------------------------
      end module advance_module
