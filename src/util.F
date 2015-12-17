!
! Utility subprograms for tgcm:
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
!-------------------------------------------------------------------
      integer function nextlu()
      implicit none
!
! Return an unopened fortan logical unit number (not 5 or 6):
! Do not return a previously returned lu.
!
      logical isopen
      integer lu
      integer,save :: lureq(99-7+1)=0   ! lu's given out so far
      do lu=7,99
        inquire(lu,opened=isopen)
        if (.not.isopen) then
!         if (lureq(lu)<=0) then   ! do not use previously returned lu
            nextlu = lu
            lureq(lu) = lu
            return
!         endif
        endif
      enddo
      write(6,"(/'>>> nextlu: all logical units apparently in use')")
      nextlu = 0
      call shutdown('nextlu')
      end function nextlu
!-------------------------------------------------------------------
      subroutine rmnamelistcomments(filename,luout,comcharin,echo)
      implicit none
! Remove comments from namelist file (filename) and write the results
! to file unit luout.  See rmcomments for additional documentation.
!
! This subroutine was created for CMIT:  when TIEGCM is coupled to
! CMIT, namelist cannot be reliably read from stdin.
!
! Args:
      character(len=*),intent(in) :: filename
      integer,intent(in) :: luout
      character(len=1),intent(in) :: comcharin
      integer,intent(in) :: echo
! Local:
      integer :: fileId=88, iErr

      open(fileId, FILE=filename, STATUS='OLD', iostat=iErr)
      if (iErr /= 0) then
         write(6,"('>>> ERROR inp_read: Could not open TIEGCM ',
     +        'input namelist file tiegcm_namelist.inp ',
     +        'with unit luin=',i2,' ios=',i5)") fileId,iErr
         call shutdown('inp_read')
      endif
      call rmcomments(fileId,luout,comcharin,echo)
      close(fileId)

      return
      end subroutine rmnamelistcomments
!-------------------------------------------------------------------
      subroutine rmcomments(luin,luout,comcharin,echo)
      implicit none
!
! Read input lines from unit lu. If current line contains the comment
!   character comcharin, strip the line from position of comchar to end,
!   and write any remaining line to a new unit. If no comment in current 
!   line, write entire line to new unit. 
! Return new unit, rewound (e.g., ready to be read by namelist).
! If echo > 0, echo output lines to stdout.
! If comcharin is ' ', then default comment char is ';'
!
! Args:
      integer,intent(in) :: luin,luout,echo
      character(len=1),intent(in) :: comcharin
! Local:
      character(len=1) :: comchar
      logical isopen
      integer :: i,lens,ios,compos,nline,nlcfields
      character(len=6400) :: line
      character(len=64) :: newcfields(30)
!
      if (luin <= 0) then
        write(6,"('>>> rmcomments: bad input luin=',i5)") luin
        call shutdown('rmcomments')
      endif
      if (luout <= 0) then
        write(6,"('>>> rmcomments: bad input luout=',i5)") luout
        call shutdown('rmcomments')
      endif
      if (len_trim(comcharin) > 0) then
        comchar = comcharin
      else
        comchar = ';'
        write(6,"('rmcomments: using default semicolon as ',
     +    'comment character.')")
      endif
      isopen = .false.
      inquire(unit=luin,opened=isopen)
      if (.not.isopen) then
        open(unit=luin,iostat=ios)
        if (ios /= 0) then
          write(6,"('>>> WARNING rmcomments: error opening input',
     +      ' file with unit luin=',i2,' ios=',i5)") luin,ios
          call shutdown('rmcomments')
        endif
      endif
      nline = 0
      read_loop: do
        line = ' '
        read(luin,"(a)",iostat=ios) line
        if (ios > 0) then
          write(6,"('>>> rmcomments: error reading from input',
     +      ' unit luin=',i3,' at line ',i5)") luin,nline
          return
        endif
        if (ios < 0) exit read_loop	! eof
        nline = nline+1
!
! Remove line if it has only "E" in column 1 (this was an
! old "Echo" directive from f77/cray namelist):
!
        if (line(1:1)=='E'.and.trim(line)=='E') cycle read_loop
!
! Use only non-commented part of line:
!
        compos = index(line,comchar)
        if (compos == 1) cycle read_loop
        if (compos > 0) line = line(1:compos-1)
        if (len_trim(adjustl(line))==0) cycle read_loop
!
! Write to new unit:
        write(luout,"(a)") trim(line)
        if (echo > 0) write(6,"(a)") line(1:len_trim(line))
      enddo read_loop  
      rewind luout
      return      
      end subroutine rmcomments
!-------------------------------------------------------------------
      subroutine fminmax(f,n,fmin,fmax)
!
! Return min and max of fields f(n) (between -1.e36,+1.e36)
!
      implicit none
!
! Args: 
      integer,intent(in) :: n
      real,intent(in) :: f(n)
      real,intent(out) :: fmin,fmax
!
! Local:
      integer :: i
!
      fmin = 1.e36
      fmax = -1.e36
      do i=1,n
        fmin = min(f(i),fmin)
        fmax = max(f(i),fmax)
      enddo
      end subroutine fminmax
!-------------------------------------------------------------------
      subroutine fminmaxspv(f,n,fmin,fmax,spv)
      implicit none
!
! Return min and max of fields f(n) (between -1.e36,+1.e36)
! Ignore any f(i)==spv.
!
! Args:
      integer,intent(in) :: n
      real,intent(in) :: f(n),spv
      real,intent(out) :: fmin,fmax
!
! Local:
      integer :: i
!
      fmin = 1.e36
      fmax = -1.e36
      do i=1,n
        if (f(i) /= spv) then
          fmin = min(f(i),fmin)
          fmax = max(f(i),fmax)
        endif
      enddo
      end subroutine fminmaxspv
!-------------------------------------------------------------------
      integer function ixfind(iarray,idim,itarget,icount)
!
! Search iarray(idim) for itarget, returning first index in iarray 
! where iarray(idim)==target. Also return number of elements of
! iarray that == itarget in icount.
!

!
! Args:
      integer,intent(in) :: idim,itarget
      integer,intent(in) :: iarray(idim)
      integer,intent(out) :: icount
!
! Local:
      integer :: i
!
      ixfind = 0
      icount = 0
      if (.not.any(iarray==itarget)) return
      icount = count(iarray==itarget)
      do i=1,idim
        if (iarray(i)==itarget) then
          ixfind = i
          exit
        endif
      enddo
      end function ixfind
!-------------------------------------------------------------------
      real function finterp(f0,f1,isec0,isec1,isec)
!
! Do linear interpolation between f0 (which is at isec0) and 
! f1 (which is at isec1) to isec.
!
! Args:
      real,intent(in) :: f0,f1
      integer,intent(in) :: isec0,isec1,isec
!
      if (isec==isec0) then
        finterp = f0
      elseif (isec==isec1) then
        finterp = f1
      else
        finterp = f0+(f1-f0)*float(isec-isec0)/float(isec1-isec0)
      endif
      end function finterp
!-------------------------------------------------------------------
      real function finterp_bigints(f0,f1,isec0,isec1,isec)
!
! Do linear interpolation between f0 (which is at isec0) and 
! f1 (which is at isec1) to isec. Same as finterp except integer
! parameters are 8-byte.
!
! Args:
      real,intent(in) :: f0,f1
      integer(kind=8),intent(in) :: isec0,isec1,isec
!
      finterp_bigints = 
     |  f0+(f1-f0)*float(isec-isec0)/float(isec1-isec0)
      end function finterp_bigints
!-------------------------------------------------------------------
      real function sddot(n,x,y)
      implicit none
!
! Call sdot (single precision) if on Cray, or ddot (double precision) 
!   if on SGI. (ddot must be called even if -r8 on sgi compiler command 
!   line). Ddot is from -lblas on the sgi.
! On IBM AIX use dot_product()
!
! 2/10/00: removing incx,incy args (i.e., incx=incy=1 on unicos
!   and irix -- IBM dot_product does not use increment args --
!   this function must be called with stride-1 vectors 
!   (see bndry.f, bndry2.f, bndrya.f, threed.f, transf.f)
!
      integer,intent(in) :: n
      real,intent(in) :: x(n),y(n)
!
#ifdef UNICOS
      real,external :: sdot
      sddot = sdot(n,x,1,y,1)
#elif SUN
      real,external :: sdot
      sddot = dot_product(x,y)
#elif IRIX
      real,external :: ddot
      sddot = ddot(n,x,1,y,1)
#elif AIX
      sddot = dot_product(x,y)
#elif OSF1
      sddot = dot_product(x,y)
#elif LINUX
      sddot = dot_product(x,y)
#else
      write(6,"('>>> WARNING sddot: unresolved OS pre-processor',
     |  ' directive.')")
#endif
      end function sddot
!-------------------------------------------------------------------
      real function vsum(n,v,inc)
!
! Call single precision vector sum "ssum" on Cray/unicos, or
!   double precision "dsum" on SGI/irix (from -lblas on the sgi):
! On IBM AIX use sum
!
      integer,intent(in) :: n,inc
      real,intent(in) :: v(n)
!
#ifdef UNICOS
      vsum = ssum(n,v,inc)
#elif IRIX
      vsum = dsum(n,v,inc)
#elif AIX
      vsum = sum(v,n)
#elif SUN
      vsum = sum(v,n)
#elif OSF1
      vsum = sum(v,n) ! a wild guess
#elif LINUX
      vsum = sum(v,n) ! another guess
#else
      write(6,"('>>> WARNING vsum: OS system cpp directive',
     |  ' not found.')")
#endif
      end function vsum
!-------------------------------------------------------------------
      complex function sdcmplx(x1,x2)
!
! Call single precision cmplx on unicos, or
! double precision dcmplx on irix:
!
      real,intent(in) :: x1,x2
!
#if defined(IRIX) || defined(AIX)
      sdcmplx = dcmplx(x1,x2)
#elif UNICOS
      sdcmplx = cmplx(x1,x2)
#elif SUN
      sdcmplx = cmplx(x1,x2)
#elif OSF1
      sdcmplx = cmplx(x1,x2)
#elif LINUX
      sdcmplx = cmplx(x1,x2)
#else
      write(6,"('>>> WARNING sddot: OS system cpp directive',
     |  ' not found.')")
#endif
      end function sdcmplx
!-------------------------------------------------------------------
      real function expo(x,iprint)
!
! To avoid overflow/underflow on ieee system, argument range to 
!   exp() must be:  -708.3964 < x  < 709.7827
!
      real,intent(in) :: x 
      integer,intent(in) :: iprint
      real,parameter :: xmin=-708., xmax=+709., 
     |  big=.1e305, small=.1e-305
!
#if defined(IRIX) || defined(AIX) || defined(OSF1) || defined(SUN) || defined(LINUX)
      if (x >= xmin .and. x <= xmax) then
        expo = exp(x)
      elseif (x < xmin) then
        if (iprint > 0) write(6,"('expo iprint=',i2,' x=',e12.4,
     |    ' setting expo = 0.')") iprint,x
        expo = 0.
      else
        if (iprint > 0) write(6,"('expo iprint=',i2,' x=',e12.4,
     |    ' setting expo = big')") iprint,x
        expo = big
      endif
#elif UNICOS
      expo = exp(x)
#endif
      end function expo
!-------------------------------------------------------------------
      integer function isystem(command)
      implicit none
!
! Execute command to the shell. UNICOS and IRIX use ishell(),
! AIX uses system().
!
      character(len=*),intent(in) :: command
      call system(trim(command))
      isystem = 0
      end function isystem
!-------------------------------------------------------------------
      real function mtime_delta(mtime0,mtime1)
!
! Return difference in time (minutes) from mtime0 to mtime1:
!
      implicit none
      integer,intent(in) :: mtime0(4),mtime1(4)
      real :: rmins0,rmins1
!
      rmins0 = mtime0(1)*1440.+mtime0(2)*60.+mtime0(3)+mtime0(4)/60.
      rmins1 = mtime1(1)*1440.+mtime1(2)*60.+mtime1(3)+mtime1(4)/60.
      mtime_delta = rmins1-rmins0
!     if (mtime_delta < 0.) then
!       write(6,"('>>> WARNING mtime_delta: negative delta: mtime0=',
!    |    4i4,' mtime1=',4i4,' delta=',f10.2)") 
!    |    mtime0,mtime1,mtime_delta
!     endif
      end function mtime_delta
!-----------------------------------------------------------------------
      real function mtime_to_datestr(iyear,mtime,imo,ida,datestr)
      implicit none
!
! Given model time mtime (day,hr,min,sec) and iyear (4-digit year),
! return imo (2 digit month), ida (2 digit day of the month), and
! a date string "minutes since yyyy-m-d". Function value is real
! minutes into the day from modeltime(2:4) hour, min, and sec.
!
! Args:
      integer,intent(in) :: iyear,mtime(4)
      integer,intent(out) :: imo,ida
      character(len=*),intent(out) :: datestr
!
! Local:
      character(len=2) :: f_mon,f_day,f_hr,f_min
      character(len=120) :: format
      integer :: ndmon(12) =
!          J  F  M  A  M  J  J  A  S  O  N  D
     |  (/31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: m,id,iday
!
      mtime_to_datestr = -1.
      datestr = ' '
      ndmon(2) = 28
      if (mod(iyear,4).eq.0) ndmon(2) = 29
!
      iday = mtime(1)
      id = 0
      do m=1,12
        id = id+ndmon(m)
        if (id.eq.iday) then
          id = ndmon(m)
          goto 100
        endif
        if (id.gt.iday) then
          id = ndmon(m)-id+iday
          goto 100
        endif
      enddo
      write(6,"('>>> iyd2date: could not find date for mtime=',4i4,
     |  ' iyear=',i4)") mtime,iyear
      imo = 0
      ida = 0
      return
 100  continue
      imo = m
      ida = id
!
! Return real minutes:
      mtime_to_datestr = float(mtime(2))*60.+float(mtime(3))+
     |                   float(mtime(4))/60.
!
! Return date string:
      if (imo <= 9) then
        if (ida <= 9) then
          f_day = 'i1'
          f_mon = 'i1'
        else ! ida > 9
          f_day = 'i2'
          f_mon = 'i1'
        endif
      else ! imo > 9
        if (ida <= 9) then
          f_day = 'i1'
          f_mon = 'i2'
        else ! ida > 9
          f_day = 'i2'
          f_mon = 'i2'
        endif
      endif
      if (mtime(2) <= 9) then
        if (mtime(3) <= 9) then
          f_hr = 'i1'
          f_min = 'i1'
        else
          f_hr = 'i1'
          f_min = 'i2'
        endif
      else
        if (mtime(3) <= 9) then
          f_hr = 'i2'
          f_min = 'i1'
        else
          f_hr = 'i2'
          f_min = 'i2'
        endif
      endif
!
! Example format:
!   ('minutes since ',i4,'-',i1,'-',i2,' ',i1,':',i2,':0')
! Example date string:
!   minutes since 1997-3-21 0:0:0
!
      format = '(''minutes since '',i4,''-'','//f_mon//',''-'','//f_day
      format = trim(format)//','' '','//f_hr//','':'','//f_min//','
      format = trim(format)//''':0'')'
      write(datestr,format) iyear,imo,ida,mtime(2),mtime(3)
      end function mtime_to_datestr
!-----------------------------------------------------------------------
      integer function mtime_to_mins(mtime)
      implicit none
!
! Given model time mtime (day,hr,min), return equivalent time 
! in minutes (includes day):
!
! Arg:
      integer,intent(in)  :: mtime(3)
!
      mtime_to_mins = mtime(1)*24*60+mtime(2)*60+mtime(3)
      end function mtime_to_mins
!-----------------------------------------------------------------------
      integer function mtime_to_nstep(mtime,istepsecs)
      implicit none
!
! Return number of steps from time 0,0,0 to given model time 
! mtime(day,hr,min) (istepsecs is time step in seconds)
! (This function used to be itera in earlier versions of tgcm)
!
! Args:
      integer,intent(in) :: mtime(3),istepsecs
!
      mtime_to_nstep=((mtime(1)*24+mtime(2))*60+mtime(3))*60/istepsecs
      end function mtime_to_nstep
!-----------------------------------------------------------------------
      integer(kind=8) function mtime_to_nsec(mtime)
      implicit none
      integer(kind=8),parameter :: nsecperday=24*60*60 ! 86400 seconds/day
      integer,intent(in) :: mtime(3)
      mtime_to_nsec = mtime(1)*nsecperday + mtime(2)*3600 + mtime(3)*60
      end function mtime_to_nsec
!-----------------------------------------------------------------------
      subroutine mins_to_mtime(mins,mtime)
      implicit none
!
! Given minutes mins, return equivalent model time (day,hr,min)
!
! Args:
      integer,intent(in) :: mins
      integer,intent(out) :: mtime(3)
!
! Local:
      integer,parameter :: minperday=24*60
!
      mtime(:) = 0
      if (mins==0) return
      mtime(1) = mins/minperday                        ! days
      mtime(2) = mod(mins,minperday)/60                ! hours
      mtime(3) = mins-(mtime(1)*minperday+mtime(2)*60) ! minutes
      end subroutine mins_to_mtime
!-------------------------------------------------------------------
      subroutine nsecs_to_modeltime(nsecs,modeltime)
      implicit none
!
! Given seconds nsecs, return equivalent model time (day,hr,min,sec)
!
! Args:
      integer(kind=8),intent(in) :: nsecs
      integer,intent(out) :: modeltime(4)
!
! Local:
      integer(kind=8),parameter :: nsecperday=24*60*60 ! 86400 seconds/day
!
      modeltime(:) = 0
      if (nsecs==0) return
      modeltime(1) = nsecs/nsecperday            ! days
      modeltime(2) = mod(nsecs,nsecperday)/3600  ! hours
      modeltime(3) = (nsecs-(modeltime(1)*nsecperday+modeltime(2)*3600))
     |               /60                         ! mins
      modeltime(4) = nsecs-(modeltime(1)*nsecperday+modeltime(2)*3600+
     |  modeltime(3)*60)                         ! seconds
      end subroutine nsecs_to_modeltime
!-------------------------------------------------------------------
      subroutine modeltime_to_nsecs(modeltime,nsecs)
!
! Given modeltime(4) (day,hr,min,sec), return total seconds
! in nsecs.
!
! Args:
      integer,intent(in)  :: modeltime(4)  
      integer(kind=8),intent(out) :: nsecs
! 
! Local:
      integer(kind=8),parameter :: nsecperday=24*60*60 ! 86400 seconds/day
!
      nsecs = modeltime(1)*nsecperday + modeltime(2)*3600 +
     |        modeltime(3)*60 + modeltime(4)
      end subroutine modeltime_to_nsecs
!-------------------------------------------------------------------
      subroutine setosys(system)
      implicit none
      character(len=*),intent(out) :: system
! 
      system = ' ' 
#ifdef UNICOS
      system = 'UNICOS'
#elif IRIX
      system = 'IRIX'
#elif AIX
      system = 'AIX'
#elif OSF1
      system = 'OSF1'
#elif SUN
      system = 'SUN'
#elif LINUX
      system = 'LINUX'
#else
      write(6,"('>>> WARNING setosys: unresolved OS cpp directive.')")
      system = 'unknown'
#endif
      end subroutine setosys
!-------------------------------------------------------------------
      subroutine datetime(curdate,curtime)
!
! Return character*8 values for current date and time.
! (sub date_and_time is an f90 intrinsic)
!
      implicit none
!
! Args:
      character(len=*),intent(out) :: curdate,curtime
!
! Local:
      character(len=8) :: date
      character(len=10) :: time
      character(len=5) :: zone
      integer :: values(8)
!
      curdate = ' '
      curtime = ' '
      call date_and_time(date,time,zone,values)
!
!     write(6,"('datetime: date=',a,' time=',a,' zone=',a)")
!    |  date,time,zone
!     write(6,"('datetime: values=',8i8)") values
!
      curdate(1:2) = date(5:6)
      curdate(3:3) = '/'
      curdate(4:5) = date(7:8)
      curdate(6:6) = '/'
      curdate(7:8) = date(3:4)
!
      curtime(1:2) = time(1:2)
      curtime(3:3) = ':'
      curtime(4:5) = time(3:4)
      curtime(6:6) = ':'
      curtime(7:8) = time(5:6)
!
      end subroutine datetime
!-------------------------------------------------------------------
      subroutine gethostsname(host)
!
! Return a host name, by various methods depending on platform:
!
      implicit none
!
! Args:
      character(len=*),intent(out) :: host
!
! Local:
      integer :: istat
!
! External:
#if defined(UNICOS)
      integer,external :: gethost
#endif
#if defined(AIX)
      integer,external :: gethostname
#endif
!
      host = ' '
#ifdef UNICOS
      istat = gethost(host)
      if (istat <= 0) then
        write(6,"('>>> gethostname UNICOS: error ',i4,' from gethost')") 
     |    istat
        host = 'unknown'
      endif
#elif defined(AIX)
      istat = gethostname(host) ! this yields "bb0001en" on babyblue
      if (istat /= 0) then
        write(6,"('>>> gethostname AIX: error ',i4,
     |    ' from gethostname')") istat
        host = 'unknown'
      endif
#elif defined(OSF1)
!
! There is a gethostname function on prospect (compaq), but apparently 
! there is no fortran binding.
      write(6,"('>>> gethostsname: do not know how to get host name ',
     |  'under OSF1')")
      host = 'unknown'
#elif defined(IRIX)
!
! Under IRIX, the interactive environment (e.g. utefe) defines $HOST, 
! whereas the batch environment (e.g. ute) defines $QSUB_HOST.
! 
      call getenv('HOST',host)                             ! interactive
      if (len_trim(host)==0) call getenv('QSUB_HOST',host) ! batch
      if (len_trim(host)==0) then
        write(6,"(/,'>>> gethost under IRIX: Cannot get HOST or ',
     |    'QSUB_HOST environment variables.',/)")
        host = 'unknown'
      endif
#elif defined(SUN) || defined(LINUX)
      call getenv('HOST',host)
      if (len_trim(host)==0) then
        write(6,"(/,'>>> gethost: Cannot get HOST environment ',
     |    'variable.',/)")
        host = 'unknown'
      endif
#else
      write(6,"('>>> WARNING gethost: unresolved OS cpp directive.')")
      host = 'unknown'
#endif
      end subroutine gethostsname
!-----------------------------------------------------------------------
      integer function strloc(strarray,nstr,substr)
      implicit none
!
! Search for string str in string array strarray(nstr).
! Its important NOT to use comparisons such as "if (trim(str0)==trim(str1))" 
!   to avoid memory leaks by Intel ifort compiler (this is not a problem 
!   w/ PGI and xlf compilers).
!
! Args:
      integer,intent(in) :: nstr
      character(len=*),intent(in) :: strarray(nstr)
      character(len=*) :: substr
!
! Local:
      integer :: i
      character(len=len(substr)) :: subtrim
      character(len=len(strarray(1))) :: strtrim
!
      strloc = 0
      aloop: do i=1,nstr
        if (len_trim(strarray(i)) > 0) then
          subtrim = trim(substr)
          strtrim = trim(strarray(i)) 
          if (subtrim == strtrim) then
            strloc = i
            exit aloop
          endif
        endif
      enddo aloop
      end function strloc
!-----------------------------------------------------------------------
      subroutine packstr(strings,nstrings,nonblank)
      implicit none
!
! Collect non-blank elements of strings(nstrings) at beginning
!   of the array. Return number of non-blank elements in nonblank. 
! On output, elements strings(1->nonblank) are non-blank (unless
!   nonblank==0), and remaining elements are blank.
!
! Args
      integer,intent(in) :: nstrings
      character(len=*),intent(inout) :: strings(nstrings)
      integer,intent(out) :: nonblank
!
! local:
      character(len=len(strings(1))) :: strings_tmp(nstrings)
      integer :: i
!
      strings_tmp(:) = ' '
      nonblank = 0
      do i=1,nstrings
        if (len_trim(strings(i)) > 0) then 
          nonblank = nonblank+1
          strings_tmp(nonblank) = strings(i)    
        endif
      enddo
      strings(:) = strings_tmp(:)
      end subroutine packstr
!-----------------------------------------------------------------------
      subroutine fftrans(a,na,work,nw,trigs,ntrigs,ifax,inc,jump,n,lot,
     |  isign)
!
! Do fft. If unicos, the original library version, FFT991  is used.
! If non-unicos, call FFT999, which is in fft9.f (the cray lib source,
!   with FFT991 changed to FFT999).
! See also setfft below.
! This is called from filter.F, e.g.:
!   call fftrans(fx,wfft,trigs,ifax,1,nlonp4,nlon,nlevs,-1)
!
      implicit none
!
! Args:
      integer,intent(in) :: inc,jump,n,lot,isign,ntrigs,na,nw
      real,intent(inout) :: a(na)
      real,intent(in) :: work(nw),trigs(ntrigs)
      integer,intent(in) :: ifax(13)
!
#ifdef UNICOS
      call fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
#else
      call fft999(a,na,work,nw,trigs,ntrigs,ifax,inc,jump,n,lot,isign)
#endif
      end subroutine fftrans
!-----------------------------------------------------------------------
      subroutine setfft(trigs,ifax,ntrigs,imax)
!
! Set up fft (called once per run from con.f).
! If unicos, the original library version, SET99  is used.
! If non-unicos, call SET999, which is in fft9.f (the cray lib source,
!   with SET99 changed to SET999.
!
      implicit none
!
! Args:
      integer,intent(in) :: ifax(13),ntrigs,imax
      real,intent(in) :: trigs(ntrigs)
!
#ifdef UNICOS
      call set99(trigs,ifax,imax)
#else
      call set999(trigs,ifax,ntrigs,imax)
#endif
      end subroutine setfft
!-----------------------------------------------------------------------
      subroutine timer(time0,tsec,begend)
!
! Timer: this routine is called twice for each timing result.
! It uses the AIX real-time-clock function rtc() (not available on
!   non-IBM systems, but more accurate than the f90 intrinsic
!   system_clock, used by the timing module in timing.F)
!
! When begend=='begin' it is the first call, and time0 is returned
!   as the beginning time.
! When begend=='end' it is the second call, time0 is input (from
!   the first call), and tsec is returned as elapsed time in seconds
!   between the 2 calls (time1-time0).
!
! If an MPI run, mpi_barrier is called to synchronize tasks before 
! timing is started in the first call, and before timing is completed 
! in the second call.
!
#ifdef MPI
      use mpi
#endif
      implicit none
!
! Args:
      real,intent(inout) :: time0  ! starting time
      real,intent(out) :: tsec      ! elapsed time in millisecs
      character(len=*),intent(in) :: 
     |  begend ! 'begin' (first call) or 'end' (second call) 
!
! Local:
      integer :: ier
      integer,save :: ncalls=0 ! only for non-AIX warning message.
      real :: time1 ! ending time
#ifdef AIX
      real,external :: rtc
#endif
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
      if (trim(begend)=='start'.or.trim(begend)=='begin') then
#ifdef MPI
        call mpi_barrier(MPI_COMM_WORLD,ier)
#endif
#ifdef AIX
        time0 = rtc()   
#endif
!
! End: time0 is now input, return time1-time0 in tsec:
      elseif (trim(begend)=='end') then
#ifdef MPI
        call mpi_barrier(MPI_COMM_WORLD,ier)
#endif
#ifdef AIX
        time1 = rtc()   
        tsec = time1-time0
#endif
      endif
      end subroutine timer
!-----------------------------------------------------------------------
      subroutine getcwdir(cwd)
#ifdef MPI
      use mpi_module,only: mytid
      use mpi
#endif
!
! Return current working directory in cwd:
!
      implicit none
      character(len=*),intent(out) :: cwd
      integer,external :: isystem,nextlu
      integer :: istat,lu,ier
!
#ifdef MPI
      if (mytid==0) then
        istat = isystem('pwd > cwd')
        if (istat /= 0) write(6,"('>>> WARNING getcwdir: error return ',
     |    'isystem(''pwd > cwd'')')")
      endif
      call mpi_barrier(MPI_COMM_WORLD,ier)
      lu = nextlu()
      open(lu,file='cwd',status='old')
      read(lu,"(a)") cwd
      close(lu)
#else
      istat = isystem('pwd > cwd')
      if (istat /= 0) write(6,"('>>> WARNING getcwdir: error return ',
     |  'isystem(''pwd > cwd'')')")
      lu = nextlu()
      open(lu,file='cwd',status='old')
      read(lu,"(a)") cwd
      close(lu)
#endif
      end subroutine getcwdir
!-----------------------------------------------------------------------
      character(len=*) function trcase(string)
!
! Translate case of input string string, i.e., return string like string,
! but with lower case character converted to upper case, and vice-versa.
!
      implicit none
!
! Args:
      character(len=*),intent(in) :: string
!
! Local:
      character(len=*),parameter :: ucase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      character(len=*),parameter :: lcase = "abcdefghijklmnopqrstuvwxyz"
      integer :: i,lpos,upos
!
      trcase = ' '
      do i=1,len_trim(string)
        trcase(i:i) = string(i:i)
        lpos = index(lcase,trcase(i:i))
        upos = index(ucase,trcase(i:i))
        if (lpos > 0) then
          trcase(i:i) = ucase(lpos:lpos) ! lower to upper
        elseif (upos > 0) then
          trcase(i:i) = lcase(upos:upos) ! upper to lower
        endif
      enddo
      end function trcase
!-----------------------------------------------------------------------
      subroutine shutdown(msg)
      use my_esmf
#ifdef MPI
      use mpi
#endif
!
! An fatal error has occurred -- shut down the model, including MPI.
!
      implicit none
!
! Args:
      character(len=*) :: msg
!
! Local:
      integer :: ier,errorcode
!     character(len=80) :: errorcode
!
      write(6,"(/,28('>'),' MODEL SHUTDOWN ',28('<'))")
      write(6,"('Shutdown: stop message: ',a)") trim(msg)
#ifdef MPI
!
! ESMF is nested inside MPI, so shutdown ESMF first:
      call ESMF_Finalize(endflag=ESMF_END_KEEPMPI,rc=ier)
      if (ier /= ESMF_SUCCESS) then
        write(6,"('>>> Shutdown: error return from ESMF_Finalize:',
     |    ' ier=',i4)") ier
      else
        write(6,"('Shutdown after ESMF_Finalize')")
      endif
!
! Shutdown MPI:
      write(6,"('Shutdown calling mpi_abort..')")
      write(6,"(72('>'))")
      call mpi_abort(MPI_COMM_WORLD,errorcode,ier)
#endif
      write(6,"(/,'For information, please see TGCM Users Guide at ',
     |  /,'http://www.hao.ucar.edu/modeling/tgcm/doc/userguide')")
      write(6,"(72('>'))")
      stop 'shutdown'
      end subroutine shutdown
!-----------------------------------------------------------------------
      integer function real_bsearch(inarray,first,last,key)
      implicit none
!
! inarray is a sorted array in asending order. The binary search routine search 
! for the beginning index where the value of key fall within in inarray (or the index
! where the value of array[index] equals to key).
!
! Args:
      integer,intent(in)::first,last
      real,dimension(first:last),intent(in) :: inarray
      real,intent(in)::key
!
! Local:
      integer :: nmid,nmin,nmax
!
      if ((key < inarray(first)) .or. (key > inarray(last))) then
          real_bsearch=-1
          return
      endif

      nmin=first
      nmax=last
      do
         if (nmin > nmax) then
            real_bsearch=nmin-1
	    exit
         end if
         nmid=(nmin+nmax)/2
	 if (key > inarray(nmid)) then
	   nmin=nmid+1
         else if (key < inarray(nmid)) then
	   nmax=nmid-1
         else
	   real_bsearch=nmid
	   exit
         end if 
      end do

      end function real_bsearch
!-----------------------------------------------------------------------
      integer function long_bsearch(inarray,first,last,key)
      implicit none
!
! inarray is a sorted array in asending order. The binary search routine search 
! for the beginning index where the value of key fall within in inarray (or the index
! where the value of array[index] equals to key).
!
! Args:
      integer,intent(in)::first,last
      integer(kind=8),dimension(first:last),intent(in) :: inarray
      integer(kind=8),intent(in)::key
!
! Local:
      integer :: nmid,nmin,nmax
!
      if ((key < inarray(first)) .or. (key > inarray(last))) then
          long_bsearch=-1
          return
      endif

      nmin=first
      nmax=last
      do
         if (nmin > nmax) then
            long_bsearch=nmin-1
	    exit
         end if
         nmid=(nmin+nmax)/2
	 if (key > inarray(nmid)) then
	   nmin=nmid+1
         else if (key < inarray(nmid)) then
	   nmax=nmid-1
         else
	   long_bsearch=nmid
	   exit
         end if 
      end do

      end function long_bsearch
!-----------------------------------------------------------------------
      integer function int_bsearch(inarray,first,last,key)
      implicit none
!
! inarray is a sorted array in asending order. The binary search routine search 
! for the beginning index where the value of key fall within in inarray (or the index
! where the value of array[index] equals to key).
!
! Args:
      integer,intent(in)::first,last
      integer,dimension(first:last),intent(in) :: inarray
      integer,intent(in)::key
!
! Local:
      integer :: nmid,nmin,nmax
!
      if ((key < inarray(first)) .or. (key > inarray(last))) then
          int_bsearch=-1
          return
      endif

      nmin=first
      nmax=last
      do
         if (nmin > nmax) then
            int_bsearch=nmin-1
	    exit
         end if
         nmid=(nmin+nmax)/2
	 if (key > inarray(nmid)) then
	   nmin=nmid+1
         else if (key < inarray(nmid)) then
	   nmax=nmid-1
         else
	   int_bsearch=nmid
	   exit
         end if 
      end do

      end function int_bsearch
!-----------------------------------------------------------------------
      subroutine check_nans(f,id1,id2,id3,name,n_total,ispval,spval,
     |  iprint,ifatal)
!
! Check for existence of +/-INF and NaN's in field f(id1,id2,id3).
!   If ispval > 0 -> replace any INF or NaNs with spval
!   If iprint==1  -> print warnings only if INF or NaNs are found
!   If iprint==2  -> always print number of INF and NaNs found
!   If ifatal > 0 -> stop program when first INF or NaNs are found
! Note: Can distinguish between +/-INF (not really NaNs), but cannot 
!       distinguish between types of actual NaNs (+/-NaNQ and NaNS).
! IBM only. See pp318-319 User's Guide Version 8 XL Fortran for AIX
!
      implicit none
!
! Args:
      integer,intent(in) :: id1,id2,id3,iprint,ifatal,ispval
      integer,intent(out) :: n_total ! total number of +/-INF+NaNs
      real,intent(inout) :: f(id1,id2,id3)
      real,intent(in) :: spval
      character(len=*),intent(in) :: name 
!
! Local:
      real :: plus_inf,minus_inf,plus_nanq,minus_nanq,sig_nan
!
! For double precision 8-byte reals (-qrealsize=8):
      data plus_inf   /z'7ff0000000000000'/  ! INF   (overflow)
      data minus_inf  /z'fff0000000000000'/  ! -INF  (underflow)
      data plus_nanq  /z'7ff8000000000000'/  ! NaNQ  (plus quiet NaN)
      data minus_nanq /z'fff8000000000000'/  ! -NaNQ (minus quiet NaN)
      data sig_nan    /z'7ff0000000000001'/  ! NaNS  (signalling NaN)
!
! For single precision (4-byte) reals:
!     data plus_inf   /z'7f800000'/  ! INF   (overflow)
!     data minus_inf  /z'ff800000'/  ! -INF  (underflow)
!     data plus_nanq  /z'7fc00000'/  ! NaNQ  (plus quiet NaN)
!     data minus_nanq /z'ffc00000'/  ! -NaNQ (minus quiet NaN)
!     data sig_nan    /z'7f800001'/  ! NaNS  (signalling NaN)
!
      integer :: i1,i2,i3
      integer :: 
     |  n_plus_inf,   ! number of INF
     |  n_minus_inf,  ! number of -INF
     |  n_nan         ! total number of NaNs (+/-NaNQ and NaNS)
!
! Init:
      n_plus_inf = 0
      n_minus_inf = 0
      n_nan = 0
      n_total = 0
!
! Scan array:
      do i3=1,id3
        do i2=1,id2
!
! +/-INF are detected by simple comparison:
          n_plus_inf   = n_plus_inf   + count(f(:,i2,i3)==plus_inf) 
          n_minus_inf  = n_minus_inf  + count(f(:,i2,i3)==minus_inf) 
!
! NaNs (NaNQ or NaNS) are detected by (a /= a):
          n_nan        = n_nan        + count(f(:,i2,i3)/=f(:,i2,i3))
          n_total = n_plus_inf+n_minus_inf+n_nan
!
!         write(6,"('i3=',i3,' i2=',i3,' n_plus_inf=',i8,' n_minus_inf='
!    |      ,i8,' n_nan=',i8,' n_total=',i8)") i3,i2,n_plus_inf,
!    |      n_minus_inf,n_nan,n_total
!
! Fatal when first INF or NaN is found:
          if (ifatal > 0 .and. n_total > 0) then
            write(6,"(/,'>>> FATAL: Found INF and/or NaNs in field ',
     |        a)") name
            write(6,"('  Dimensions id1,id2,id3=',3i4)") id1,id2,id3
            write(6,"('  First INF or NaN found at id2=',i4,', id3=',
     |        i4)") i2,i3
            write(6,"('  n_plus_inf   = ',i6)") n_plus_inf
            write(6,"('  n_minus_inf  = ',i6)") n_minus_inf
            write(6,"('  n_nan (NaNS or NaNQ) = ',i6)") n_nan
            write(6,"('  data(:,',i3,',',i3,') = ',/,(6e12.4))") 
     |        i2,i3,f(:,i2,i3)
            call shutdown('check_nans')
          endif ! ifatal > 0
!
! Replace any INF or NaNs with spval:
          if (ispval > 0 .and. n_total > 0) then
            do i1=1,id1
              if (f(i1,i2,i3)==plus_inf.or.f(i1,i2,i3)==minus_inf.or.
     |            f(i1,i2,i3)/=f(i1,i2,i3)) f(i1,i2,i3) = spval
            enddo
          endif 
        enddo ! i2=1,id2
      enddo ! i3=1,id3
!
! Print level 1 (print warnings only if INF or NaNs are found):
      if (iprint==1) then
        if (n_plus_inf > 0) write(6,"('>>> WARNING: found ',
     |    i6,' INF values in field ',a,' (id1,2,3=',3i4,')')") 
     |    n_plus_inf,name,id1,id2,id3
        if (n_minus_inf > 0) write(6,"('>>> WARNING: found ',
     |    i6,' -INF values in field ',a,' (id1,2,3=',3i4,')')") 
     |    n_minus_inf,name,id1,id2,id3
        if (n_nan > 0) write(6,"('>>> WARNING: found ',i6,
     |    ' NaNS or NaNQ values in field ',a,' (id1,2,3=',3i4,')')") 
     |    n_nan,name,id1,id2,id3
!       if (ispval > 0 .and. n_total > 0) 
!    |    write(6,"('>>> Replaced ',i8,' values with spval ',e12.4)")
!    |      n_total,spval
!
! Print level 2 (always print number of nans found):
      elseif (iprint==2) then 
        write(6,"('Checking for INF and NaNs in field ',a,' id1,2,3=',
     |    3i4)") name,id1,id2,id3
        print *,'  n_plus_inf   (',plus_inf,  ') = ',n_plus_inf
        print *,'  n_minus_inf  (',minus_inf, ') = ',n_minus_inf
        print *,'  n_nan        (',plus_nanq,'+',sig_nan,') = ',n_nan
        print *,'  n_total      (total INF+NaNs) = ',n_total
!       if (ispval > 0)
!    |  print *,'  Replaced ',n_total,' values with spval ',spval
      endif
      end subroutine check_nans
!-------------------------------------------------------------------
      subroutine finterpb(f0,f1,fout,n,isec0,isec1,isec)
      implicit none
!
! Do linear interpolation between f0 (which is at isec0) and
! f1 (which is at isec1) to isec. f0,f1,and fout are arrays
! dimensions n.
!
! Args:
      integer,intent(in) :: n
      real,intent(in) :: isec0,isec1,isec
      real,intent(in) :: f0(n),f1(n)
      real,intent(out) :: fout(n)
!
! Local:
      integer :: i
!       
      do i=1,n
        fout(i) = f0(i)+(f1(i)-f0(i))*(isec-isec0)/
     |    (isec1-isec0)
      enddo
      end subroutine finterpb
!-------------------------------------------------------------------
      recursive subroutine expand_path(path)
!
! Expand any environment variables imbedded in path, and return 
!   expanded path. 
! Procedure:
!   If '$' is found in input path, then an env var is defined as 
!   that part of path following the '$' up to (not including) the 
!   next delimiter. The value of the env var is substituted in place 
!   of the env var string. If no '$' is found, the routine returns 
!   without changing path. 
! Environment vars can be set (using setenv) in the user's .cshrc file, 
!   in the job script (e.g., setenv from a shell var), or set manually 
!   in the shell before executing the model.
!
! The 7 recognized delimiters (meaning end of env var name) are:
!   '/' (forward slash), 
!   '.' (dot), 
!   '_' (underscore), 
!   '-' (dash), 
!   ':' (colon), 
!   '#' (pound sign), and 
!   '%' (percent sign)
!
! This routine is recursive, so multiple env vars can be used in the
!   same path, and in combination with different delimiters, see 
!   examples below.
!
! Examples:
!   path = '$TGCMDATA/dir1/file.nc'   (the env var is $TGCMDATA)
!   path = '$MYDIR/$MYSUBDIR/file.nc' (env vars are $MYDIR, $MYSUBDIR)
!   path = '$USER.$MODEL_$NUM.nc'     (3 env vars and different delims)
!   path = '$FILEPATH'                (entire path in one env var)
! Last example:
!   In the job script:
!     set model = $tiegcm  ! set a shell var
!     setenv MODEL $model  ! set env var from shell var
!   In the namelist input:
!     histfile = '$TGCMDATA/TGCM.$MODEL.p001-2002-080.nc' or
!     histfile = '$TGCMDATA/TGCM.$MODEL.p001-$YEAR-$DAY.nc'
!
      implicit none
!
! Args:
      character(len=*),intent(inout) :: path
!
! Local:
      character(len=224) :: path_out,envvar_value
      character(len=80) :: envvar_name
      integer,parameter :: ndelim=7
      character(len=1) :: delimiters(ndelim) = 
     |  (/ '/', '.', '-', '_', ':', '#', '%'/)
      integer :: i,idollar,idelim
!
      if (len_trim(path)==0) then
        write(6,"('>>> WARNING expand_path: path is empty.')")
        return
      endif
!
!     write(6,"('Enter expand_path: path=',a)") trim(path)
!
      idollar = index(path,'$')
      if (idollar <= 0) return ! no env var in path
!
! Env var is between idollar and next slash 
! (or end of path if there is no slash after idollar):
!
      idelim = 0
      do i=idollar+1,len_trim(path) ! find next delimiter
        if (any(delimiters==path(i:i))) then
          idelim = i
          exit
        endif
      enddo
      if (idelim <= 0) idelim = len_trim(path)+1
      envvar_name = path(idollar+1:idelim-1)

!     write(6,"('expand_path: path=',a,' idollar=',i3,
!    |  ' idelim=',i3,' envvar_name=',a)") 
!    |  trim(path),idollar,idelim,trim(envvar_name)

!
! Get value of env var (getenv is f90 intrinsic):
      call getenv(trim(envvar_name),envvar_value)
      if (len_trim(envvar_value) <= 0) then
        write(6,"('>>> WARNING expand_path: error retrieving ',
     |    'value for env var ''',a,'''')") trim(envvar_name)
        return
      else
!       write(6,"('expand_path: envvar=',a,' value=',a)")
!    |    trim(envvar_name),trim(envvar_value)
      endif
!
! Put together the expanded output path:
      if (idollar > 1) then
        if (idelim < len_trim(path)) then
          path_out = path(1:idollar-1)//trim(envvar_value)//
     |      path(idelim:len_trim(path))
        else
          path_out = path(1:idollar-1)//trim(envvar_value)
        endif
      else     ! idollar == 1
        if (idelim < len_trim(path)) then
          path_out = trim(envvar_value)//path(idelim:len_trim(path))
        else
          path_out = trim(envvar_value)
        endif
      endif
!
! Return new path, and make recursive call for more env vars:
      path = trim(path_out)
!     write(6,"('expand_path returning path = ''',a,'''')") trim(path)
!
! Recursive call to expand any additional env vars:
      call expand_path(path) ! expand next env var
!
      end subroutine expand_path
!-----------------------------------------------------------------------
      logical function arrayeq(a0,a1,n)
      implicit none
      real,intent(in) :: a0(n),a1(n) 
      integer,intent(in) :: n
      integer :: i
!     
      arrayeq = .true.
      do i=1,n
        if (a0(i) /= a1(i)) then
          arrayeq = .false.
          return 
        endif
      enddo
      end function arrayeq
!-----------------------------------------------------------------------
      real function hp_from_bz_swvel(bz,swvel)
!
! Calculate hemispheric power from bz, swvel:
! Emery, et.al., (2008, in press, JGR)
! 6/3/08: Enforce minimum hp of 4.0 before *fac.
! 6/6/08: Reset minimum hp from 4.0 to 2.5 before *fac,
!         as per Emery email of 6/5/08.
!
      implicit none
      real,intent(in) :: bz,swvel  ! in
      real :: hp                   ! out
      real :: fac = 2.0
!
      if (bz < 0.) then
        hp = 6.0 + 3.3*abs(bz) + (0.05 + 0.003*abs(bz))*
     |    (min(swvel,700.)-300.)
      else
        hp = 5.0 + 0.05 * (min(swvel,700.)-300.)
      endif
!     hp = max(4.0,hp)*fac
      hp = max(2.5,hp)*fac
      hp_from_bz_swvel = hp
      end function hp_from_bz_swvel
!-----------------------------------------------------------------------
      real function hp_from_kp(kp)
!
! Calculate hemispheric power from Kp.
!
      implicit none
      real,intent(in) :: kp
!
      hp_from_kp = 0.
      if (kp < 0..or.kp > 9.) then
        write(6,"('>>> hp_from_kp: Bad kp=',e12.4)")
        call shutdown('hp_from_kp')
      endif
!
! Formula for hemispheric power:
! modified by LQIAN, 2007
! for fkp<=7: formula given by Zhang Yongliang based on TIMED/GUVI
! for fkp>7: power is 153.13 when fkp=7 from Zhang's formula, 
! assume power is 300.(based on NOAA satellites) when fkp=9 
! do linear interporation in between
!
      if (kp <=7.) hp_from_kp = 16.82*exp(0.32*kp)-4.86
      if (kp > 7.) hp_from_kp = 153.13 + (kp-7.)/(9.-7.)*(300.-153.13)
      end function hp_from_kp
!-----------------------------------------------------------------------
      real function ctpoten_from_kp(kp)
!
! Calculate cross-tail potential from Kp.
!
      implicit none
      real,intent(in) :: kp
!
      ctpoten_from_kp = 0.
      if (kp < 0..or.kp > 9.) then
        write(6,"('>>> ctpoten_from_kp: Bad kp=',e12.4)") kp
        call shutdown('ctpoten_from_kp')
      endif
!
! Formula for potential:
! modified by LQIAN, 2007
! formula given by Wenbin based on data fitting
!
      ctpoten_from_kp = 15.+15.*kp + 0.8*kp**2
      end function ctpoten_from_kp
!-----------------------------------------------------------------------
      logical function time2print(nstep,istep)
      implicit none
        integer,intent(in) :: nstep,istep
        time2print = .false.
        if (nstep <= 100 .or. (nstep > 100 .and. mod(istep,10)==0))
     |    time2print = .true.
      end function time2print
!-----------------------------------------------------------------------
      real function fslt(slt,ut,glon,iflag)
      implicit none
!
! If iflag = 1 -> return slt, given ut and glon (slt=ut+glon/15.)
! If iflag = 2 -> return ut, given slt and glon (ut=slt-glon/15.)
! If iflag = 3 -> return glon, given slt and ut (glon=(slt-ut)*15.)
! (slt and ut always in decimal hours, glon in -180 -> 180)
!
! Args:
      real :: slt,ut,glon
      integer,intent(in) :: iflag
!
      if (iflag==1) then
        fslt = ut + glon/15.
        if (fslt < 0.) fslt = fslt+24.
        if (fslt >= 24.) fslt = fslt-24.
      elseif (iflag==2) then
        fslt = slt - glon/15.  
        if (fslt < 0.) fslt = fslt+24.
        if (fslt > 24.) fslt = fslt-24.
      elseif (iflag==3) then
        fslt = (slt-ut) * 15.
        if (fslt >= 180.) fslt = fslt-360.
        if (fslt < -180.) fslt = fslt+360.
      else
        write(6,"('>>> fslt: bad iflag=',i3)") iflag
        stop 'fslt'
      endif
      return
      end
!-------------------------------------------------------------------
