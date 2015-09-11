!
      module timing_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use hist_module,only: nstep
      use input_module,only: step
      implicit none
#ifdef AIX
      real,external :: rtc
#endif
#ifdef MPI
#include <mpif.h>
#endif
!
! timing%level == 0 -> no timing
! timing%level == 1 -> timing on broad code segments only
! timing%level == 2 -> detailed selected per-routine timing
!
! timing%rtc == .true. -> use rtc timing (AIX only)
! timing%sys == .true. -> use f90 intrinsic system_clock
!
      type timing_type
        integer :: 
     |    level       ! 0 (no timing), 1 (segment timing), 2 (per routine timing)
        logical :: 
     |    rtc,        ! use real-time-clock timer (AIX only)
     |    sys         ! use system_clock f90 intrinsic
        real ::
     |    run,        ! timing for entire run
     |    init,       ! time spent before time loop
     |    step,       ! per step timing (not including i/o)
     |    prep,       ! preparing and concluding a time step 
     |    dynamics,   ! timing for dynamics (dynamics.F)
     |    pdynamo,    ! timing for pdynamo
     |    phist,      ! timing for primary history i/o
     |    shist,      ! timing for secondary history i/o
     |    io,         ! total i/o (phist+shist)
     |    qrj,        ! timing for qrj
     |    gw,         ! timing for gravity waves (timegcm only)
     |    oplus,      ! time in sub oplus
     |    cmpminor,   ! minor species n2d, n4s, no
     |    cmpmajor,   ! major species o2 and o
     |    dt,         ! time in dt
     |    duv,        ! time in duv
     |    magpres_grav,
     |    apxparm
      end type timing_type
      type (timing_type) :: timing
!
      integer,parameter :: nsegs = 18
      character(len=16) :: timing_segs(nsegs) =
     |  (/'RUN             ','INIT            ','STEP            ',
     |    'PREP            ','DYNAMICS        ','PDYNAMO         ',
     |    'PHIST           ','SHIST           ','IO              ',
     |    'QRJ             ','GW              ','OPLUS           ',
     |    'CMPMINOR        ','CMPMAJOR        ','DT              ',
     |    'DUV             ','MAGPRES_GRAV    ','APXPARM         '/)
      contains
!-----------------------------------------------------------------------
      subroutine init_timer
!
! Real-Time-Clock (rtc) is IBM only, whereas the f90 intrinsic 
!   system_clock is used for cross-platform comparison.
! 
! Warning: do not set both rtc and sys true.
!
!     timing%level    = 2       ! adds timing statements on per-timestep basis
      timing%level    = 1       ! print timing results only at end of run
      timing%rtc      = .false. ! default is no rtc
      timing%sys      = .true.  ! default is to use system_clock
#ifdef AIX
!
! Use rtc timing under AIX (IBM only):
!     timing%rtc      = .true.  
!     timing%sys      = .false.
!
! Use sys_clock timing for comparison with non-IBM platforms: 
      timing%rtc      = .false.
      timing%sys      = .true.  
#endif
!
! sys timing parameters:
      timing%run      = 0.
      timing%init     = 0.
      timing%step     = 0.
      timing%prep     = 0.
      timing%dynamics = 0.
      timing%pdynamo  = 0.
      timing%phist    = 0.
      timing%shist    = 0.
      timing%io       = 0.
      timing%qrj      = 0.
      timing%gw       = 0.
!
      write(6,"('init_timer: level=',i3,' rtc=',l1,' sys=',l1)")
     |  timing%level, timing%rtc, timing%sys

      end subroutine init_timer
!-----------------------------------------------------------------------
      subroutine timer(time0,time1,seg,ibegend,ibarrier)
!
! Args:
      real,intent(inout) :: time0
      real,intent(out) :: time1
      character(len=*) :: seg
      integer,intent(in) :: ibegend,ibarrier
!
! Local:
      real :: esecs 
      integer :: i
!
! Use system_clock:
      if (timing%sys) then                  ! using sys timing
        if (ibegend == 0) then              ! start sys timing
          call start_timing(time0,seg,ibarrier)     ! time0 is output
          time1 = 0.
        else                                ! end sys timing
          call end_timing(time0,seg,esecs,ibarrier) ! time0 in, esecs out
          time1 = esecs
        endif  
      endif
!
! Use rtc:
      if (timing%rtc) then                           ! using rtc timing
        if (ibegend == 0) then
          call rtc_timer(time0,time1,seg,0,ibarrier) ! start rtc timing
          time1 = 0.
        else
          call rtc_timer(time0,time1,seg,1,ibarrier) ! end rtc timing
          esecs = time1
        endif  
      endif
!
! Update module data, if appropriate:
      if (ibegend == 1) then  ! end timing (rtc or sys)
        select case(trim(seg))
          case('RUN')
            timing%run = timing%run+esecs
          case('INIT')
            timing%init = timing%init+esecs
          case('STEP')
            timing%step = timing%step+esecs
          case('PREP')
            timing%prep= timing%prep+esecs
          case('DYNAMICS')
            timing%dynamics= timing%dynamics+esecs
          case('QRJ')
            timing%qrj= timing%qrj+esecs
          case('GW')
            timing%gw= timing%gw+esecs
          case('OPLUS')
            timing%oplus= timing%oplus+esecs
          case('CMPMINOR')
            timing%cmpminor= timing%cmpminor+esecs
          case('CMPMAJOR')
            timing%cmpmajor= timing%cmpmajor+esecs
          case('DT')
            timing%dt= timing%dt+esecs
          case('DUV')
            timing%duv= timing%duv+esecs
          case('PDYNAMO')
            timing%pdynamo= timing%pdynamo+esecs
          case('PHIST')
            timing%phist= timing%phist+esecs
            timing%io= timing%io+esecs
          case('SHIST')
            timing%shist= timing%shist+esecs
            timing%io= timing%io+esecs
          case('MAGPRES_GRAV')
            timing%magpres_grav= timing%magpres_grav+esecs
          case('APXPARM')
            timing%apxparm= timing%apxparm+esecs
          case default
            write(6,"('>>> timer: unknown segment: ',a)") trim(seg)
        end select
      endif ! end rtc timing
      end subroutine timer
!-----------------------------------------------------------------------
      subroutine start_timing(count_out,string,ibarrier)
!
! f90 intrinsic system_clock returns current clock count, number of
! clock counts per second, and max possible value of the count.
!
! Args:
      real,intent(out) :: count_out
      character(len=*),intent(in) :: string
      integer,intent(in) :: ibarrier
!
! Local:
      integer :: icount1,icount_rate,icount_max,ier
!
      call system_clock(icount1,icount_rate,icount_max)
      count_out = float(icount1)
#ifdef MPI
        if (ibarrier > 0) call mpi_barrier(MPI_COMM_WORLD,ier)
#endif

!     write(6,"('start_timing for ',a,': secs=',f10.3)") 
!    |  trim(string),real(icount1)/real(icount_rate)

      end subroutine start_timing
!-----------------------------------------------------------------------
      subroutine end_timing(count_in,string,esecs,ibarrier)
!
! f90 intrinsic system_clock returns current clock count, number of
! clock counts per second, and max possible value of the count.
!
! Args:
      real,intent(in) :: count_in
      real,intent(out) :: esecs
      character(len=*),intent(in) :: string
      integer,intent(in) :: ibarrier
!
! Local:
      integer :: icount1,icount2,icount_rate,icount_max,ier
      real :: secs1,secs2
!
      call system_clock(icount2,icount_rate,icount_max)
      icount1 = int(count_in)
      if (icount2 < icount1) icount2 = icount2+icount_max
      secs1 = real(icount1)/real(icount_rate)
      secs2 = real(icount2)/real(icount_rate)
      esecs = secs2-secs1
#ifdef MPI
        if (ibarrier > 0) call mpi_barrier(MPI_COMM_WORLD,ier)
#endif

!     write(6,"('end_timing for ',a,' : secs1=',f10.3,' secs2=',f10.3,
!    |  ' elapsed=',f10.3)") trim(string),secs1,secs2,esecs

      end subroutine end_timing
!-----------------------------------------------------------------------
      subroutine rtc_timer(time0,tsec,seg,ibegend,ibarrier)
!
! Timer: this routine is called twice for each timing result.
! It uses the AIX real-time-clock function rtc() (not available on
!   non-IBM systems, but more accurate than the f90 intrinsic
!   system_clock, used by the timing module in timing.F)
!
! If an MPI run, mpi_barrier is called to synchronize tasks before 
! timing is started in the first call, and before timing is completed 
! in the second call.
!
      implicit none
!
! Args:
      real,intent(inout) :: time0  ! starting time
      real,intent(out) :: tsec      ! elapsed time in millisecs
      integer,intent(in) :: ibegend,ibarrier
      character(len=*),intent(in) :: seg 
!
! Local:
      integer :: ier
      integer,save :: ncalls=0 ! only for non-AIX warning message.
      real :: time1 ! ending time
!
      ncalls = ncalls+1
#ifndef AIX
      if (ncalls==1)
     |  write(6,"('>>> timer: rtc() not available on non-AIX ',
     |    ' systems')")
      return
#endif
      tsec = 0.
!
! Begin: return time0 as starting time:
      if (ibegend == 0) then
#ifdef MPI
        if (ibarrier > 0) call mpi_barrier(MPI_COMM_WORLD,ier)
#endif
#ifdef AIX
        time0 = rtc()   
!       write(6,"('rtc_timer begin: seg ',a,' time0=',f15.2)") 
!    |    trim(seg),time0
#endif
!
! End: time0 is now input, return time1-time0 in tsec:
      else
#ifdef MPI
        if (ibarrier > 0) call mpi_barrier(MPI_COMM_WORLD,ier)
#endif
#ifdef AIX
        time1 = rtc()   
        tsec = time1-time0
!       write(6,"('rtc_timer end: seg ',a,' elapsed=',f8.2)") 
!    |    trim(seg),tsec
#endif
      endif
      end subroutine rtc_timer
!-----------------------------------------------------------------------
      subroutine timer_report
!
! Local:
      real :: av_secs_per_step
      integer :: nsecs
! 
      write(6,"(/,72('-'))")
      if (timing%sys) then
        write(6,"('TIMER (system_clock):')")
      elseif (timing%rtc) then
        write(6,"('TIMER (AIX rtc):')")
      endif
      write(6,"('Elapsed run time      = ',f10.2,' (secs)  ',
     |  f6.2,' (mins)  ',f5.2,' (hrs)')")
     |  timing%run, timing%run/60., timing%run/3600.
      write(6,"(' ')")

      write(6,"('Elapsed init time     = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%',' (includes source i/o)')") 
     |  timing%init,timing%init/3600.,
     |  (timing%init*100.)/timing%run

      write(6,"('Elapsed step time     = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%step,
     |  timing%step/3600., (timing%step*100.)/timing%run

      write(6,"('  Elapsed prep time   = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%prep,
     |  timing%prep/3600., (timing%prep*100.)/timing%run

      write(6,"('  Elapsed apxparm time= ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%apxparm,
     |  timing%apxparm/3600., (timing%apxparm*100.)/timing%apxparm

      write(6,"('  Elapsed dynamics    = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%dynamics,
     |  timing%dynamics/3600., 
     |  (timing%dynamics*100.)/timing%run
!
! Timing for modules called from dynamics:
!
! mgw is timegcm only:
!     write(6,"('    Elapsed mgw       = ',f8.2,' secs, ',
!    |  f6.3,' hrs, ',f5.1,'%')") timing%gw,
!    |  timing%gw/3600., 
!    |  (timing%gw*100.)/timing%run

      write(6,"('    Elapsed qrj       = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%qrj,
     |  timing%qrj/3600., 
     |  (timing%qrj*100.)/timing%run

      write(6,"('    Elapsed oplus     = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%oplus,
     |  timing%oplus/3600., 
     |  (timing%oplus*100.)/timing%run

      write(6,"('    Elapsed cmpminor  = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%',' (n4s,n2d,no)')") 
     |  timing%cmpminor,timing%cmpminor/3600., 
     |  (timing%cmpminor*100.)/timing%run

      write(6,"('    Elapsed cmpmajor  = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%',' (o2, o, he)')")
     |  timing%cmpmajor,timing%cmpmajor/3600., 
     |  (timing%cmpmajor*100.)/timing%run

      write(6,"('    Elapsed dt        = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%dt,timing%dt/3600., 
     |  (timing%dt*100.)/timing%run

      write(6,"('    Elapsed duv       = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%duv,timing%duv/3600., 
     |  (timing%duv*100.)/timing%run
!
! Add other reports from dynamics here..
!
      write(6,"('  Elapsed pdynamo     = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%pdynamo,
     |  timing%pdynamo/3600., 
     |  (timing%pdynamo*100.)/timing%run
 
      write(6,"('  Elapsed magpres_grav= ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%magpres_grav,
     |  timing%magpres_grav/3600., 
     |  (timing%magpres_grav*100.)/timing%run

      write(6,"('Elapsed i/o           = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%io,
     |  timing%io/3600., 
     |  (timing%io*100.)/timing%run

      write(6,"('  Elapsed Primary io  = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%phist,
     |  timing%phist/3600., 
     |  (timing%phist*100.)/timing%run

      write(6,"('  Elapsed Secondary   = ',f8.2,' secs, ',
     |  f6.3,' hrs, ',f5.1,'%')") timing%shist,
     |  timing%shist/3600., 
     |  (timing%shist*100.)/timing%run

      av_secs_per_step = timing%step/float(nstep)
      write(6,"(/,'Average secs/step = ',f8.2)") av_secs_per_step
      nsecs = nstep*step
      write(6,"('Average mins/simulated day  = ',f8.2)")
     |  (timing%run/(float(nsecs)/(24.*3600.)))/60.

      write(6,"(72('-'),/)")
      end subroutine timer_report
!-----------------------------------------------------------------------
      end module timing_module
