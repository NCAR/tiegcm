!
      integer function numfiles(type,nseries,nsource,nhist_output1,
     |  nhist_total)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Determine number of history files needed (primary if type=='prim',
! or secondary if type=='sech'). A new file will be needed whenever 
! a new time series is started, or a file is full (i.e., has 
! mxhist_[prim,sech] histories). This is called from input.
!
      use input_module,only: step,
     |  mxhist_prim,start   ,stop,   hist,
     |  mxhist_sech,secstart,secstop,sechist
      implicit none
!
! Args:
      character(len=4),intent(in) :: type  ! 'prim' or 'sech'
      integer,intent(in) :: 
     |  nseries,       ! number of time series
     |  nsource,       ! number of source files (0/1)
     |  nhist_output1  ! number of histories on 1st output file
                       ! (known only if nsource==0)
      integer,intent(out) :: nhist_total   ! total histories to be written
!
! Local:
      integer :: i,mtime(3),nhfile,nhseries,modeltime(4)
      integer(kind=8) :: nsecs,idelsecs,nsecstop
      integer(kind=8),external :: mtime_to_nsec
      integer :: mtime_curr(4),mtime_stop(4)
!
! Primary histories:
      if (type=='prim') then
        numfiles = 0      ! return var
        nhist_total = 0   ! total number of histories (output arg)
!
! Iterate through each primary history time series:
        do i=1,nseries
          idelsecs = mtime_to_nsec(hist(:,i)) ! delta secs between histories
          nsecstop = mtime_to_nsec(stop(:,i)) ! stop secs
          numfiles = numfiles+1 ! start new file at each new time series
          nsecs = mtime_to_nsec(start(:,i)) ! time in secs
! debug:
!         call nsecs_to_modeltime(nsecs,mtime_curr)
!         call nsecs_to_modeltime(nsecstop,mtime_stop)
!         write(6,"(/,'numfiles prim: i=',i2,' mtime_curr=',3i4,
!    |      ' mtime_stop=',3i4,' numfiles=',i2)") 
!    |      i,mtime_curr(1:3),mtime_stop(1:3),numfiles
! end debug
!
! If start history of a new time series is the same as last history
! of the previous time series, we do not rewrite the start history 
! on the new file (nsecs was incremented at end of previous series).
! Otherwise, set nsecs to first history in the new series:
!
          if (i > 1) then
            if (all(start(:,i) == stop(:,i-1)))
     |        nsecs = nsecs+idelsecs
          endif
!
! If continuation run, start with number of hist on 1st output file
! (this is known only when rdsource calls numfiles):
!
          nhfile = 0   ! number of histories in current file
          if (nsource==0.and.i==1) nhfile = nhist_output1-1
!
! Iterate over histories in this time series:
!
          nhseries = 0 ! number of histories in current time series
          primloop: do ! loop over histories in this time series
            nhseries = nhseries+1
            nhfile = nhfile+1
            if (nhfile > mxhist_prim) then
              numfiles = numfiles+1
              nhfile = 1
            endif
! debug:
!           call nsecs_to_modeltime(nsecs,mtime_curr)
!           write(6,"(' primhist: i=',i2,' mtime_curr=',3i4,
!    |        ' nhseries=',i2,' nhfile=',i2,' numfiles=',i2)")
!    |        i,mtime_curr(1:3),nhseries,nhfile,numfiles
! end debug
            nsecs = nsecs+idelsecs
            if (nsecs > nsecstop) exit primloop
          enddo primloop
          nhist_total = nhist_total+nhseries
!         write(6,"('prim nhist_total=',i2)") nhist_total
        enddo ! i=1,nseries
!
! When source file not given (i.e., continuation run), the first
!   (source) history is not echoed to the file, so nhist_total is 
!   one less:
!
        if (nsource==0) then
          nhist_total = nhist_total-1
!         write(6,"('prim: nsource==0, so decrement nhist_total',
!    |      ' by 1 to ',i2)") nhist_total
        endif
!
! Secondary histories:
      elseif (type=='sech') then
        numfiles = 0     ! return var
        nhist_total = 0  ! total number of histories (output arg)
        do i=1,nseries
          idelsecs = mtime_to_nsec(sechist(:,i))
          nsecstop = mtime_to_nsec(secstop(:,i))
          numfiles = numfiles+1 ! start new file at each new time series
          nsecs = mtime_to_nsec(secstart(:,i))
! debug:
!         call nsecs_to_modeltime(nsecs,mtime_curr)
!         call nsecs_to_modeltime(nsecstop,mtime_stop)
!         write(6,"(/,'numfiles sech: i=',i2,' mtime_curr=',3i4,
!    |      ' mtime_stop=',3i4,' numfiles=',i2)") 
!    |      i,mtime_curr(1:3),mtime_stop(1:3),numfiles
! end debug
!
! If start history of a new time series is the same as last history
! of the previous time series, we do not rewrite the start history 
! on the new file:
!
          if (i > 1) then
            if (all(secstart(:,i) == secstop(:,i-1)))
     |        nsecs = nsecs+idelsecs
          endif
          nhfile = 0   ! number of histories in current file
          nhseries = 0 ! number of histories in current time series
          sechloop: do ! loop over histories in this time series
            nhseries = nhseries+1
            nhfile = nhfile+1
            if (nhfile > mxhist_sech) then
              numfiles = numfiles+1
              nhfile = 1
            endif
! debug:
!           call nsecs_to_modeltime(nsecs,mtime_curr)
!           write(6,"(' sech: i=',i2,' mtime_curr=',3i4,
!    |        ' nhseries=',i2,' nhfile=',i2,' numfiles=',i2)")
!    |        i,mtime_curr(1:3),nhseries,nhfile,numfiles
! end debug
            nsecs = nsecs+idelsecs
            if (nsecs > nsecstop) exit sechloop
          enddo sechloop
          nhist_total = nhist_total+nhseries
!         write(6,"('sech nhist_total=',i2)") nhist_total
        enddo ! i=1,nseries
!
! Type must be either "prim" or "sech":
      else
        write(6,"('>>> numfiles: unrecognized type = ',a)") trim(type)
        call shutdown('numfiles')
      endif
      end function numfiles
