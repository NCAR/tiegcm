!
      logical function wrhist(istep,modeltime,
     |  wrprim, newseries_prim, iseries_p,
     |  wrsech, newseries_sech, iseries_s)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Return true if it is time to write a primary or secondary history, 
!   given current istep and modeltime(4).
! Return wrprim true if its time to write a primary history.
! Return wrsech true if its time to write a secondary history.
! Return newseries_prim true if a new primary time series is beginning.
! Return newseries_sech true if a new secondary time series is beginning.
!
! The current time series number is returned in iseries_p for primary
! histories, and iseries_s for secondaries. These are on the stack so
! this routine can be called from single (e.g., from sub outhist) or 
! multiple threads/tasks (e.g. advnce and addfld).
!
      use input_module,only: step,start,stop,hist,secstart,secstop,
     |  sechist
      use hist_module,only: nsource,isechist,nseries_prim,
     |  nseries_sech
      implicit none
!
! Args:
      integer,intent(in) :: 
     |  istep,          ! current step number in model run (0-nstep)
     |  modeltime(4)    ! current model time (day,hr,min,sec)
      integer,intent(out) :: iseries_p,iseries_s
      logical,intent(out) :: 
     |  wrprim,  newseries_prim,
     |  wrsech,  newseries_sech
!
! Local:
      integer :: i, 
     |  nstep_hist,    ! number of steps between history writes
     |  istep_prim,    ! current step in current primary series
     |  istep_sech     ! current step in current secondary series
      integer(kind=8) :: 
     |  nsecs,         ! current modeltime in seconds
     |  nsec_hist,     ! number of seconds between history writes
     |  nsec_start,    ! start time of current series (secs)
     |  nsec_stop,     ! stop time of current series (secs)
     |  nsec_stop_prev ! stop time of previous series (secs)
!
! External:
      integer,external :: mtime_to_nstep 
      integer(kind=8),external :: mtime_to_nsec
!
! Init:
!     write(6,"('wrhist: istep=',i4,' modeltime=',4i4)") istep,modeltime
      wrhist = .false.
      wrprim = .false. 
      wrsech = .false.
      newseries_prim = .false.
      newseries_sech = .false.
      call modeltime_to_nsecs(modeltime,nsecs) ! sets nsecs
!
! Determine current primary history time series:
! (can be zero, i.e., model is in between time series)
!
      iseries_p = 0
      primloop: do i=1,nseries_prim
        nsec_start = mtime_to_nsec(start(:,i))
        nsec_stop  = mtime_to_nsec(stop (:,i))
        if (nsecs >= nsec_start.and.nsecs <= nsec_stop) then
          iseries_p = i
!
! Determine if this is beginning of a new primary time series:
          if (nsecs==nsec_start) then
            newseries_prim = .true.
            if (i > 1) then
              nsec_stop_prev = mtime_to_nsec(stop(:,i-1))
              if (nsecs == nsec_stop_prev) newseries_prim = .false.
            endif
          else
            if (i > 1) then
              nstep_hist = mtime_to_nstep(hist(:,i),step)
              nsec_hist = nstep_hist*step
              nsec_stop_prev = mtime_to_nsec(stop(:,i-1))
              if (nsecs-nsec_hist == nsec_start.and.
     |            nsec_start == nsec_stop_prev) newseries_prim = .true.
            endif
          endif
          exit primloop
        endif
      enddo primloop
!
! Step 0 of run:
      if (istep==0) then
        if (nsource > 0) then
          wrhist = .true.
          wrprim = .true.
          newseries_prim = .true.
        endif
      else
!
! Check for primary history write:
        if (iseries_p > 0) then
          nsec_start = mtime_to_nsec(start(:,iseries_p))
          nsec_stop  = mtime_to_nsec(stop (:,iseries_p))
          nstep_hist = mtime_to_nstep(hist(:,iseries_p),step)

!         write(6,"('wrhist call mtime_to_nstep: iseries_p=',
!    |      i2,' step=',i4)") iseries_p,step

          istep_prim = (nsecs-nsec_start)/step
          if ((nsecs >= nsec_start.and.nsecs <= nsec_stop).and.
     |        mod(istep_prim,nstep_hist)==0) then
            wrhist = .true.
            wrprim = .true.
          endif
        endif
      endif
!
! Determine current secondary history time series index:
! (can be zero, i.e., model is in between time series)
!
      if (isechist > 0) then ! writing secondary histories this run
        iseries_s = 0
        sechloop: do i=1,nseries_sech
          nsec_start = mtime_to_nsec(secstart(:,i))
          nsec_stop  = mtime_to_nsec(secstop (:,i))
          if (nsecs >= nsec_start.and.nsecs <= nsec_stop) then
            iseries_s = i
!
! Determine if this is beginning of a new secondary time series:
            if (nsecs==nsec_start) then
              newseries_sech = .true.
              if (i > 1) then
                nsec_stop_prev = mtime_to_nsec(secstop(:,i-1))
                if (nsecs == nsec_stop_prev) newseries_sech = .false.
              endif
            else
              if (i > 1) then
                nstep_hist = mtime_to_nstep(sechist(:,i),step)
                nsec_hist = nstep_hist*step
                nsec_stop_prev = mtime_to_nsec(secstop(:,i-1))
                if (nsecs-nsec_hist == nsec_start.and.
     |              nsec_start == nsec_stop_prev) newseries_sech=.true.
              endif
            endif
            exit sechloop
          endif
        enddo sechloop
!
! Step 0 of run:
        if (istep==0) then
          if (isechist > 0.and.all(modeltime(1:3)==secstart(:,1))) then
            wrhist = .true.
            wrsech = .true.
            newseries_sech = .true.
          endif
        else
!
! Check secondary histories:
          if (isechist > 0 .and. iseries_s > 0) then
            nsec_start = mtime_to_nsec(secstart(:,iseries_s))
            nsec_stop  = mtime_to_nsec(secstop (:,iseries_s))
            nstep_hist = mtime_to_nstep(sechist(:,iseries_s),step)
            istep_sech = (nsecs-nsec_start)/step
!
! Check for secondary history write:
            if (nsecs >= nsec_start.and.nsecs <= nsec_stop.and.
     |        mod(istep_sech,nstep_hist)==0) then
                wrhist = .true.
                wrsech = .true.
            endif
          endif
        endif
      endif
!     write(6,"(8x,'wrhist=',l1,' wrprim=',l1,
!    |  ' newseries_prim=',l1,' iseries_p=',i3)")
!    |  wrhist,wrprim,newseries_prim,iseries_p
      end function wrhist
