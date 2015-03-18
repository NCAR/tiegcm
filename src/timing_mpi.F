      module mpitime_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! MPI timing module.
!
      implicit none
!
! Intrinsic MPI wall-clock timer:
      real,external :: mpi_wtime
!
! Max number of code segs to be timed:
      integer,parameter :: mxsegs=100
      integer :: nsegs
!
! Timing_mpi type contains timing info for a particular code segment:
      type mpitime_seg_type
        character(len=80) :: name
        real :: time
        real :: total_time
      end type mpitime_seg_type
      type (mpitime_seg_type) :: codesegs(mxsegs)
      contains
!-----------------------------------------------------------------------
      subroutine mpitime_init
      integer :: i
      do i=1,mxsegs
        codesegs(i)%name = ' '
        codesegs(i)%time = 0.
        codesegs(i)%total_time = 0.
      enddo
      nsegs = 0
      write(6,"('mpitime_init: Completed initialization of mpi ',
     |  'timing')")
      end subroutine mpitime_init
!-----------------------------------------------------------------------
      subroutine mpi_timer(name,startstop,iprint)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: startstop,iprint
!
! Local:
      integer :: i,iseg,ncalls=0
      real :: current_time
!
! Get current time:
      current_time = mpi_wtime()
!
! Init if first call:
      ncalls = ncalls+1
      if (ncalls==1) call mpitime_init

      if (len_trim(name)==0) then
        write(6,"('>>> mpi_timer: blank name -- timing will not ',
     |    'be reported')")
        return
      endif
!
! Check for name:
      iseg = 0
      loop1: do i=1,nsegs
        if (trim(name)==trim(codesegs(i)%name)) then
          iseg = i
          exit loop1 
        endif
      enddo loop1
      if (iseg == 0) then ! is a new code segment
        nsegs = nsegs+1
        if (nsegs > mxsegs) then
          write(6,"('>>> mpi_timer: too many code segments: mxsegs=',
     |      i4,' Please increase parameter mxsegs')") mxsegs
          call shutdown('mpi_timer mxsegs')
        endif
        iseg = nsegs
        codesegs(iseg)%name = name

        if (iprint > 0)
     |    write(6,"('mpi_timer: added seg ',a,' nsegs=',i3)")
     |      trim(codesegs(nsegs)%name),nsegs

      endif
!
! Start timing for this segment:
      if (startstop==0) then
        codesegs(iseg)%time = current_time
!
! Stop timing for this segment:
      else
        codesegs(iseg)%time = current_time - codesegs(iseg)%time
        codesegs(iseg)%total_time = codesegs(iseg)%total_time + 
     |    codesegs(iseg)%time
        if (iprint > 0) 
     |    write(6,"('mpi_timer: iseg=',i3,' segment=',a,' time=',e12.4,
     |      ' total_time=',e12.4)") iseg,trim(codesegs(iseg)%name),
     |      codesegs(iseg)%time,codesegs(iseg)%total_time
      endif
      end subroutine mpi_timer
!-----------------------------------------------------------------------
      subroutine report_mpitime
      use init_module,only: istep
      use hist_module,only: nstep
      use timing_module,only: timing
      use mpi_module,only: mytid

      integer :: i
      write(6,"(/,72('-'))")
      write(6,"('Report MPI wall-clock timing (mpi_wtime): ',
     |  ' mytid=',i4,'  istep=',i8,'  nstep=',i8)") mytid,istep,nstep
!
! If total runtime is available (last timestep), then include percent
! of total run time spent in this segment:
!
      if (timing%run > 0.) then
        do i=1,nsegs
          write(6,"('Total Time (mins) in ''',a,''' = ',e12.4,
     |      ' (',f6.2,'% of total run time)')") 
     |      codesegs(i)%name(1:30),codesegs(i)%total_time/60.,
     |      (codesegs(i)%total_time/60.)*100./(timing%run/60.)
        enddo
!
! If total runtime not available (not last timestep), print most recent
! time calculated (for current timestep, or loop iteration), and total
! accumulated time in this segment so far in the run: 
!
      else
        do i=1,nsegs
          write(6,"('Time (mins) in ''',a,''' = ',
     |      e12.4,' Total time this segment = ',e12.4)")
     |      codesegs(i)%name(1:30),codesegs(i)%time/60.,
     |      codesegs(i)%total_time/60.
        enddo
      endif
      write(6,"(72('-'),/)")
      end subroutine report_mpitime
!-----------------------------------------------------------------------
      end module mpitime_module
