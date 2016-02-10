!
      program tgcm
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Main program for tgcm models.
!
! initial cvs import was tiegcm1
!
      use params_module,only: tgcm_version
      use input_module,only: input,dynamo,step,cwd,pid,hpss_path,
     |  potential_model,start_year
      use dispose_module,only: init_dispose
      use fields_module,only: itp
      use init_module,only: init,iyear
      use hist_module,only: nstep
      use advance_module,only: advance
      use cons_module,only: init_cons
      use lbc,only: bndry_diurnal, bndry_semidiurnal, bndcmp
      use qrj_module,only: init_qrj
      use lbc,only: init_lbc
      use magfield_module,only: magdyn
      use amie_module,only: init_amie
      use getapex_module,only: apxparm
      use timing_module,only: timer,timer_report,init_timer
      use mpi_module,only: mytid,ntask

      use diags_module,only: init_diags
      use cons_module,only: consdyn,set_geogrid
      use mpitime_module,only: report_mpitime
#ifdef MPI
      use mpi_module,only: mp_init,mp_close,mp_distribute_geo,
     |  mp_distribute_mag,mp_exchange_tasks,report_mpi_timing,
     |  mpi_timing,time_totalrun
      use my_esmf,only: esmf_init
#endif
#if defined(INTERCOMM) || defined(CISMAH)
      use cism_coupling_module,only: initialize,finalize
#endif
#ifdef OMP
!     use omp_module,only: init_omp
#endif
      implicit none
!
! Local:
      integer :: ier,icount_tgcm,icount_apex
      real :: cpu1,cpu2
      character(len=8) ::
     |  rundate,           ! current local date
     |  runtime            ! current local time
      real ::
     |  time0,     time1,
     |  time0_apx, time1_apx,
     |  time0_run, time1_run,
     |  time0_init,time1_init
      real :: starttime,endtime
      real :: hours,days
      real,external :: mpi_wtime
!
#ifndef IRIX
      call cpu_time(cpu1)
#endif
!
! Report starting time and other info to stdout:
      call startup_message

#ifdef INTERCOMM
!
! This code block is required for coupling with LFM (CMIT):
!
! InterComm 2.0 requires that IC_Setprogname be called prior to
! MPI_Init.
      call IC_Setprogname("tiegcm")
#endif

!
! Initialize MPI:
#ifdef MPI
      call mp_init
      if (mpi_timing) starttime = mpi_wtime()
#endif

!
! Initialize timing (must be called after mp_init):
      call init_timer
!
! Init timing for the run, get cwd and pid, get user input,
!   set up 2d decomposition.
!
      call timer(time0_run,time1_run,'RUN',0,0)    ! start run timing
      call timer(time0_init,time1_init,'INIT',0,0) ! start init timing
      call getcwdir(cwd)
      write(6,"('Current working directory (cwd) = ',a)") trim(cwd)
!
! Namelist read and input verification:
!
      call input(mytid,ntask) ! get user input (mpi)
      if (mytid==0.and.len_trim(hpss_path) > 0) call init_dispose
!
! set_geogrid and consdyn (both in cons.F) must be called before
! mp_distribute, for use by sub define_mag.
      call set_geogrid
      call consdyn
!
! Call apex code if doing new dynamo
      if (dynamo > 0.or.
     |   (dynamo <= 0 .and.trim(potential_model) /= 'NONE')) then
        call timer(time0_apx,time1_apx,'APXPARM',0,0)
        write(6,"('tgcm: dynamo=',i2,' -- calling apxparm: start_year=',
     |    i5)") dynamo,start_year
        call apxparm(real(start_year))
        call timer(time0_apx,time1_apx,'APXPARM',1,0)
        write(6,"('Time in apxparm = ',f6.3,' (secs)')") time1_apx
      endif
!
#ifdef MPI
!
! Do 2d decomposition of geographic and magnetic grids:
      call mp_distribute_geo
      call mp_distribute_mag
!
! Make all task structures available to all tasks:
      call mp_exchange_tasks(1)
!
! Set up geo and mag ESMF grids (esmf.F):
      call esmf_init
#endif

#if defined(INTERCOMM) || defined(CISMAH)
!
! Initialize code coupling communication framework 
! (i.e. InterComm or AdHoc file exchanges)
!
      call initialize
#endif
!
! Do initialization:
      call init
!
! Initialize diagnostic fields for secondary histories:
      call init_diags(1)
!
! Initialization for qrj:
! (this is not in init_module to avoid circular dependency between
!  init_module and qrj_module)
      call init_qrj
!
! Read amie data files if requested
! (this is not in init_module to avoid circular dependencies)
      call init_amie
!
! Read source history:
      call readsource(ier)
!
! Define magnetic field and related quantities:
      call magdyn
!
! Set several lower boundary parameters (including Hough mode functions):
      call init_lbc ! lbc module

      call timer(time0_init,time1_init,'INIT',1,0) ! end init timing
!
! Advance the model (timing in main time-step loop is done in advance):
      call advance
!
! Report to stdout:
      write(6,"(' ')")
#ifdef MPI
      write(6,"('MPI run with ntask = ',i3)") ntask
#endif
      write(6,"('nstep=',i8,' step=',i8)") nstep,step
      hours = real(nstep*step)/3600. ! seconds to hours
      days  = hours/24.
      write(6,"('Model simulation time: hours=',f10.2,' days=',f10.2)")
     |  hours,days
!
! End timing:
!
#ifndef IRIX
      call cpu_time(cpu2)
      write(6,"('Cpu time for run = ',f10.2)") cpu2-cpu1
#endif
      call timer(time0_run,time1_run,'RUN',1,0) ! end total run timing

#if defined(INTERCOMM) || defined(CISMAH)
!
! Finalize code coupling communication framework 
! (i.e. InterComm or AdHoc file exchanges)
!
      call finalize
#endif
!
! Finalize mpi:
!
#ifdef MPI
!
! Report detailed mpi timing (see mpi.F):
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_totalrun = time_totalrun+(endtime-starttime)
        call report_mpi_timing
      endif
!
! report_mpitime reports on dynamo timing using mpi_wtime function.
      call report_mpitime
      call mp_close ! mpi.F
#endif
      call timer_report
!
      call final_message
      write(6,"('NORMAL EXIT')")
      end program tgcm
!-----------------------------------------------------------------------
      subroutine startup_message
      use params_module,only: tgcm_version
      character(len=8) ::
     |  rundate,           ! current local date
     |  runtime            ! current local time
      character(len=16) ::
     |  host,              ! host machine
     |  system,            ! operating system of host (from pre-proc macros)
     |  logname            ! user login name

      call datetime(rundate,runtime) 
      call gethostsname(host)
      call setosys(system)
      logname = ' '
      call getenv('LOGNAME',logname)
      if (len_trim(logname)==0) logname = "unknown"
      write(6,"(/,72('='))")
      write(6,"('Begin execution of ',a,' at ',a,'  ',a)")
     |  tgcm_version,rundate,runtime
      write(6,"('  Host    = ',a)") trim(host)
      write(6,"('  System  = ',a)") trim(system)
      write(6,"('  Logname = ',a)") trim(logname)
      write(6,"  (72('='),/)")
      end subroutine startup_message
!-----------------------------------------------------------------------
      subroutine final_message
!
! Print end-of-execution message to stdout with date and time:
!
      use params_module,only: tgcm_version
      character(len=8) ::
     |  rundate,           ! current local date
     |  runtime            ! current local time
      call datetime(rundate,runtime)
      write(6,"('End execution of ',a,' at ',a,'  ',a)")
     |  trim(tgcm_version),rundate,runtime
      end subroutine final_message
!-----------------------------------------------------------------------
