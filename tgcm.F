!
      program tgcm
!
! Main program for tgcm models.
!
! initial cvs import was tiegcm1
!
      use params_module,only: tgcm_version
      use input_module,only: input,dynamo,step
      use fields_module,only: itp
      use init_module,only: init,iyear
      use hist_module,only: nstep
      use cons_module,only: init_cons
      use bndry_module,only: lowbound
      use qrj_module,only: init_sigmas,init_qrj,init_euvac
      use magfield_module,only: magfield
      use apex_module,only: apxparm
      use timing_module
#ifdef MPI
      use mpi_module,only: mp_init,mp_close,mytid,ntask,mp_distribute
#endif
#ifdef OMP
!     use omp_module,only: init_omp
#endif
      implicit none
!
! Local:
      integer :: ier,icount_tgcm,icount_apex,nsecs
      real :: elapsed,cpu1,cpu2
      character(len=8) ::
     |  rundate,           ! current local date
     |  runtime            ! current local time
      real,external :: cputime  ! util.F
!
! Report starting time and other info to stdout:
      call startup_message
!
! Initialize timing for entire run:
      call start_timing(icount_tgcm,'total run')
      cpu1 = cputime()
!
! Get user input. If MPI job, init mpi and set up 2-d decomposition
!   across tasks.
!
#ifdef MPI
      call mp_init
      call input(mytid,ntask) ! get user input (mpi)
      call mp_distribute
#else
      call input(-1,0)    ! get user input (non-mpi)
#endif
!
! Init threads if doing shared memory:
#ifdef OMP
      call init_omp
#endif
!
! Do initialization:
      call init
!
! Initialization for qrj:
! (this is not in init_module to avoid circular dependency between
!  init_module and qrj_module)
      call init_qrj
      call init_euvac
      call init_sigmas
!
! Read source history:
      call readsource(ier)
!
! Call apex code if doing new dynamo (should eliminate some of the
! magfield calls if this is set).
      if (dynamo == 2) then
        call start_timing(icount_apex,'apxparm')
        call apxparm(real(iyear))
        call end_timing(icount_apex,elapsed)
        write(6,"('Elapsed secs for apxparm=',f6.3)") elapsed
      endif
!
! Read magnetic data file and set magnetic field parameters:
! (data file not read if dynamo==2, since in this case apxparm
!  was called above)
!
      call magfield(dynamo)
!
! Set lower boundaries:
      call lowbound
!
! Advance the model (timing in main time-step loop is done in advance):
      call advance
!
! If MPI job, finalize mpi:
!
#ifdef MPI
      call mp_close
      write(6,"('MPI run with ntask = ',i3)") ntask
#endif
!
! Report timing stats:
!
      write(6,"(' ')")

      nsecs = nstep*step
      write(6,"('Model simulation time = ',i8,' secs ',/,
     |  '  (minutes=',f8.2,', hours=',f8.2,', days=',f6.2,')')") 
     |  nsecs,float(nsecs)/60.,float(nsecs)/3600.,
     |  float(nsecs)/(24.*3600.)
      write(6,"('Elapsed secs for ',i6,' time-steps = ',f10.2,
     |  ' (ave per step=',f6.2,')')") nstep,elapsed_steps,
     |  elapsed_steps/real(nstep)
#ifdef MPI
      write(6,"('Elapsed secs for gather2root=',f10.2)") elapsed_mpi
#endif
      write(6,"('Elapsed secs for writing primary histories=',f10.2)")
     |  elapsed_prim
      write(6,"('Elapsed secs for writing secondary histories=',f10.2)")
     |  elapsed_sech
!
      call end_timing(icount_tgcm,elapsed)
      cpu2 = cputime()
      write(6,"('Elapsed secs for run = ',f10.2,' (minutes=',f10.2, 
     |  ', hours=',f8.2,')')") elapsed,elapsed/60.,elapsed/3600.
      write(6,"('Cpu time for run = ',f10.2)") cpu2-cpu1
!
      call datetime(rundate,runtime) 
      write(6,"('End execution of ',a,' at ',a,'  ',a)")
     |  tgcm_version,rundate,runtime
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
