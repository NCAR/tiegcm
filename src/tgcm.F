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
      use input_module,only: input,dynamo,step,cwd,pid,hpss_path
      use dispose_module,only: init_dispose
      use fields_module,only: itp
      use init_module,only: init,iyear
      use hist_module,only: nstep
      use advance_module,only: advance
      use cons_module,only: init_cons
      use lbc,only: bndry_diurnal, bndry_semidiurnal, bndcmp
      use qrj_module,only: init_qrj
      use magfield_module,only: magfield
      use apex_module,only: apxparm
      use timing_module,only: timer,timer_report,init_timer
      use mpi_module,only: mytid,ntask
#ifdef MPI
      use mpi_module,only: mp_init,mp_close,mp_distribute
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
      integer :: ier,icount_tgcm,icount_apex,nsecs
      real :: cpu1,cpu2
      character(len=8) ::
     |  rundate,           ! current local date
     |  runtime            ! current local time
      real ::
     |  time0,     time1,
     |  time0_apx, time1_apx,
     |  time0_run, time1_run,
     |  time0_init,time1_init
      real,external :: cputime  ! util.F
!
! Report starting time and other info to stdout:
      call startup_message
!
! Initialize timing for entire run:
#ifndef IRIX
      call cpu_time(cpu1)
#endif
!
! Initialize timing (must be called after mp_init):
      call init_timer
!
! Get user input. If MPI job, init mpi and set up 2-d decomposition
!   across tasks. Start timing contains an mpi_barrier, so must be
!   called after mp_init.
!
#ifdef MPI
      call mp_init
#endif
!
! Init timing for the run, get cwd and pid, get user input,
!   set up 2d decomposition.
!
      call timer(time0_run,time1_run,'RUN',0,0)    ! start run timing
      call timer(time0_init,time1_init,'INIT',0,0) ! start init timing
      call getcwd(cwd)
      write(6,"('Current working directory (cwd) = ',a)") trim(cwd)
      call getpid(pid)
      write(6,"('Process ID (pid) = ',i8)") pid
      call input(mytid,ntask) ! get user input (mpi)
      if (mytid==0.and.len_trim(hpss_path) > 0) call init_dispose
#ifdef MPI
      call mp_distribute
#endif
#if defined(INTERCOMM) || defined(CISMAH)
!
! Initialize code coupling communication framework 
! (i.e. InterComm or AdHoc file exchanges)
!
      if(mytid==0) call initialize
#endif
!
! Do initialization:
      call init
!
! Initialization for qrj:
! (this is not in init_module to avoid circular dependency between
!  init_module and qrj_module)
      call init_qrj
!
! Read source history:
      call readsource(ier)
!
! Call apex code if doing dynamo.
      call apxparm(real(iyear))
!
! Set up magnetic field data (read from magfield file if dynamo<=0) 
      call magfield
!
! Set lower boundary Hough mode functions:
      call bndry_diurnal
      call bndry_semidiurnal
      call bndcmp

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
      write(6,"('nstep=',i5,' step=',i5)") nstep,step
      nsecs = nstep*step
      write(6,"('Model simulation time = ',i8,' secs ',/,
     |  '  (minutes=',f10.2,', hours=',f8.2,', days=',f10.6,')')")
     |  nsecs,float(nsecs)/60.,float(nsecs)/3600.,
     |  float(nsecs)/(24.*3600.)
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
      if (mytid==0)call finalize
#endif

!
! Finalize mpi:
!
#ifdef MPI
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
