!
! Utility subprograms for tgcm:
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
          if (lureq(lu)<=0) then   ! do not use previously returned lu
            nextlu = lu
            lureq(lu) = lu
            return
          endif
        endif
      enddo
      write(6,"(/'>>> nextlu: all logical units apparently in use')")
      nextlu = 0
      stop 'nextlu'
      end function nextlu
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
      character*80 line
      character(len=64) :: newcfields(30)
!
      if (luin <= 0) then
        write(6,"('>>> rmcomments: bad input luin=',i5)") luin
        stop 'rmcomments'
      endif
      if (luout <= 0) then
        write(6,"('>>> rmcomments: bad input luout=',i5)") luout
        stop 'rmcomments'
      endif
      if (len_trim(comcharin) > 0) then
        comchar = comcharin
      else
        comchar = ';'
        write(6,"('rmcomments: using default semicolon as ',
     +    'comment character.')")
      endif
      inquire(unit=luin,opened=isopen)
      if (.not.isopen) then
        open(unit=luin,iostat=ios)
        if (ios /= 0) then
          write(6,"('>>> WARNING rmcomments: error opening input',
     +      ' file with unit luin=',i2,' ios=',i5)") luin,ios
          stop 'rmcomments'
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
      subroutine mkdiskflnm(msspath,diskname)
      implicit none
!
! Given a mss path (msspath), construct and return a local 
! file name (diskname) which is the same as msspath except
! that slashes ('/') in the msspath are replaced by dots ('.')
! in the diskname. (The initial slash in msspath is ignored)
!
! For example, if on input msspath = '/FOSTER/dir1/file', then
! diskname would be returned as 'FOSTER.dir1.file'
!
! Args:
      character(len=*),intent(in) :: msspath
      character(len=*),intent(out) :: diskname
!
! Local:
      integer :: lmsspath,i,ii
!
      lmsspath = len_trim(msspath)
      if (lmsspath==0) then
        write(6,"('WARNING mkdiskflnm: zero length msspath.',
     +    ' Returning diskname: ''disk.file''')") 
        diskname = 'disk.file'
        return
      endif
      diskname = ' '
      do i=1,len(diskname)
!
! Do not replace leading slash of msspath. 
! If there are no occurrences of '/' in the msspath (i.e., its
!   not really an mss path, but just a file name), return the 
!   msspath unchanged.
!
        ii = i+1 
        if (index(msspath,'/')==0) ii = i
        if (ii > lmsspath) exit
        if (msspath(ii:ii)=='/') then
          diskname(i:i)='.'
        else
          diskname(i:i) = msspath(ii:ii)
        endif
      enddo
      if (len_trim(diskname)==0) then
        write(6,"('WARNING mkdiskflnm: zero length diskname output.')")
        write(6,"('  msspath=',a)") trim(msspath)
        write(6,"('Returning diskname: ''disk.file''')")
        diskname = 'disk.file'
      endif
!     write(6,"('mkdiskflnm returning:',/,'  msspath=',a,
!    +  /,'  diskname=',a)") trim(msspath),trim(diskname)
      end subroutine mkdiskflnm
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
        fmin = amin1(f(i),fmin)
        fmax = amax1(f(i),fmax)
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
          fmin = amin1(f(i),fmin)
          fmax = amax1(f(i),fmax)
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
      finterp = f0+(f1-f0)*float(isec-isec0)/float(isec1-isec0)
      end function finterp
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
#else
      write(6,"('>>> WARNING sddot: OS system cpp directive',
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
#if defined(IRIX) || defined(AIX) || defined(OSF1) || defined(SUN)
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
#if defined(UNICOS) || defined(IRIX)
      integer,external :: ishell
      isystem = ishell(trim(command))
#elif defined (AIX)
      integer,external :: system
      isystem = system(trim(command)//"\0")
#elif defined(OSF1)
      integer,external :: system
      isystem = system(trim(command))
#elif defined(SUN)
      integer,external :: system
      isystem = system(trim(command))
#endif
      end function isystem
!-------------------------------------------------------------------
      real function cputime()
!
! Return user cpu time:
!
      implicit none
!
! Cray unicos and SGI IRIX use the second() function:
!
#if defined(UNICOS) || defined(IRIX)
      real(kind=4),external :: second
      cputime = second()
#elif AIX
!
! IBM AIX uses sub cpu_time
! To use other than default time type use setrteopts, e.g.:
!     call setrteopts('cpu_time_type=usertime')
! (see p.589 of AIX Language Reference)
!
      real :: time
      call cpu_time(time)
      cputime = time
#elif OSF1
      cputime = 0. ! don't know of one on OSF, except in ncaru
#elif SUN
      cputime = 0. ! don't know of one on OSF, except in ncaru
#endif
      end function cputime 
!-------------------------------------------------------------------
      integer function ilink(from,to,iprint)
      implicit none
!
! Args:
      character(len=*),intent(in) :: from,to
      integer,intent(in) :: iprint
!
! External:
      integer,external :: link
!
#if defined(UNICOS) || defined(IRIX) || defined(SUN) || defined(OSF1)
      ilink = link(trim(from),trim(to))
      if (ilink.eq.0) then
        if (iprint > 0) 
     |  write(6,"('Linked file ',a,' to ',a)") trim(from),trim(to)
      else
        write(6,"(/,'>>> WARNING: link command of ',a,' to ',a,
     |    ' failed: ierror=',i3,/)") trim(from),trim(to),ilink
!       stop 'ilink'
      endif
#elif defined(AIX) 
!
! Note: link will fail if executing from different file system
!       than tmpdir (e.g., executing from /home but tmpdir is
!       in /ptmp)
!
      ilink = link(trim(from)//"\0",trim(to)//"\0")
      if (ilink.eq.0) then
        if (iprint > 0) 
     |    write(6,"('Linked file ',a,' to ',a)") trim(from),trim(to)
      else
        write(6,"('>>> WARNING: link command of ',a,' to ',a,
     |    ' failed: ierror=',i3)") trim(from),trim(to),ilink
      endif
#endif
      end function ilink
!-------------------------------------------------------------------
      integer function iunlink(file,iprint)
      implicit none
      character(len=*),intent(in) :: file
      integer,intent(in) :: iprint
      integer,external :: unlink
!
#if defined(UNICOS) || defined(IRIX) || defined(SUN) || defined(OSF1)
      iunlink = unlink(trim(file))
      if (iunlink.eq.0) then
        if (iprint > 0) 
     |    write(6,"('Unlinked file ',a)") trim(file)
      else
        if (iprint > 0) 
     |    write(6,"('Note: unlink of ',a,' failed',
     |      ' (possibly non-existant file).')") trim(file)
      endif
#elif AIX
      iunlink = unlink(trim(file)//"\0")
      if (iunlink.eq.0) then
        if (iprint > 0) 
     |    write(6,"('Unlinked file ',a)") trim(file)
      else
        if (iprint > 0) 
     |    write(6,"('Note: unlink of ',a,' failed',
     |      ' (possibly non-existant file).')") trim(file)
      endif
#endif
      end function iunlink
!-----------------------------------------------------------------------
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
#else
      write(6,"('>>> WARNING setosys: unresolved OS cpp directive.')")
      system = 'unknown'
#endif
      end subroutine setosys
!-------------------------------------------------------------------
      subroutine settmpdir(logname,model_version,tmpdir)
!
! If tempdir is not provided by the user via input, set it here
!   according to platform. These defaults may need to be changed as
!   directory systems on individual machines change.
! If the user does not want to set it, and this routine is not
!   suitable, the user may set tempdir='.' in input. This will
!   work on all systems, but in this case history files will not be 
!   linked to tempdir, and could be lost if not saved to mss.
!   (e.g., if cwd is $TMPDIR, or some such scratch directory)
! Note however, that if the cwd and tmpdir are in different file
!   systems, the link (e.g., from savefile) will fail.
!
      implicit none
      character(len=*),intent(in) :: logname,model_version
      character(len=*),intent(out) :: tmpdir
!
      tmpdir = ' '
#ifdef UNICOS
      tmpdir = '/tmp/'//trim(logname)//'/'//trim(model_version)
#elif IRIX
      tmpdir = '/ptmp/'//trim(logname)//'/'//trim(model_version)
#elif AIX
      tmpdir = '/ptmp/'//trim(logname)//'/'//trim(model_version)
      tmpdir = trim(tmpdir)//'/tmpdir'
!     tmpdir = '/ptmp/'//trim(logname)//'/'//trim(model_version)
!     tmpdir = '~/hist/hist'
!     tmpdir = '/home/babyblue/foster/hist'
!     tmpdir = '.'
#elif OSF1
      tmpdir = '/ptmp/'//trim(logname)//'/'//trim(model_version)
#elif SUN
!     tmpdir = '/home/tgcm/dev/hist/hist'
!     tmpdir = '/e/foster/tgcm15'
      tmpdir = '/e/foster/tiegcm1/tmpdir'
#else
      write(6,"('>>> WARNING settmpdir: unresolved OS cpp directive.')")
      write(6,"('            Will set to cwd (tempdir = ''.'')')")
      tmpdir = '.'
#endif
!     write(6,"('settmpdir: set tempdir = ',a)") trim(tmpdir)
      end subroutine settmpdir
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
#elif defined(SUN)
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
      logical function isadigit(s)
      implicit none
      character(len=1),intent(in) :: s
      character(len=1) :: digits(10)=
     +  (/'0','1','2','3','4','5','6','7','8','9'/)
      integer :: i
      isadigit = .false.
      do i=1,10
        if (s==digits(i)) then
          isadigit = .true.
          exit
        endif
      enddo
      end function isadigit
!-----------------------------------------------------------------------
      integer function decode_int(str,sout)
      implicit none
!
! Return integer decoded from digits in string str.
! Only use digits from end of str. Also return string
! sout, with integers replaced by dollar signs.
!
! Args:
      character(len=*),intent(in) :: str
      character(len=*),intent(out) :: sout
!
! Local:
      character(len=len(str)) :: snum,s
      character(len=8) :: format
      integer :: i,ipos,slen
!
! External:
      logical,external :: isadigit
!
      decode_int = -1
      s = adjustl(str)
      slen = len_trim(s)
      if (slen>9999) then
        write(6,"('>>> decode_int: string too long:',
     +    ' slen=',i10)") slen
        return
      endif
      ipos = 0
      snum = ' '
      sout = trim(s)
      do i=slen,1,-1
        if (.not.isadigit(s(i:i))) exit
        ipos = ipos+1
        snum(slen-ipos+1:slen-ipos+1) = s(i:i)
        sout(i:i) = '$'
      enddo
      if (ipos==0) return
      snum = adjustl(snum)
      write(format,"('(i',i4,')')") len_trim(snum)
!     read(trim(snum),format,err=100) decode_int ! IBM does not like this
      read(snum,format,err=100) decode_int
      return
 100  write(6,"('>>> decode_int: error making integer from ',
     +  'string: ',a)") trim(snum)
      decode_int = -1
      end function decode_int
!-----------------------------------------------------------------------
      subroutine encode_str(str,num)
      implicit none
!
! Given str containing n characters '$', return
! same string, but with the '$' replaced by num.
!
      character(len=*),intent(inout) :: str
      character(len=len(str)) :: newstr
      integer,intent(in) :: num
      integer :: i,ipos,ndollar
      character(len=16) :: format
!
      ndollar = 0
      do i=1,len(str)
        if (str(i:i)=='$') ndollar = ndollar+1
      enddo
      if (ndollar==0) return
      newstr = str
      write(format,"('(i',i4,'.',i4,')')") ndollar,ndollar
      format = trim(adjustl(format))
      ipos = index(str,'$')
      write(newstr(ipos:ipos+ndollar-1),format,err=100) num
      str = newstr
      return
 100  write(6,"('>>> encode_str: error writing integer ',
     +  'num=',i12,' to string.')") num
      return
      end subroutine encode_str
!-----------------------------------------------------------------------
      integer function strloc(strarray,nstr,str)
      implicit none
!
! Search for string str in string array strarray(nstr).
!
! Args:
      integer,intent(in) :: nstr
      character(len=*),intent(in) :: strarray(nstr)
      character(len=*) :: str
!
! Local:
      integer :: i
!
      strloc = 0
      aloop: do i=1,nstr
        if (len_trim(strarray(i)) > 0) then
          if (trim(str) == trim(strarray(i))) then
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
      subroutine fftrans(a,work,trigs,ifax,inc,jump,n,lot,isign)
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
      integer,intent(in) :: inc,jump,n,lot,isign
      real,intent(inout) :: a(n)
      real,intent(in) :: work(n),trigs(n)
      integer,intent(in) :: ifax(13) 
!
#ifdef UNICOS
      call fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
#else
      call fft999(a,work,trigs,ifax,inc,jump,n,lot,isign)
#endif
      end subroutine fftrans
!-----------------------------------------------------------------------
      subroutine setfft(trigs,ifax,imax)
!
! Set up fft (called once per run from con.f). 
! If unicos, the original library version, SET99  is used.
! If non-unicos, call SET999, which is in fft9.f (the cray lib source,
!   with SET99 changed to SET999.
!
      implicit none
!
! Args:
      integer,intent(in) :: ifax(13),imax
      real,intent(in) :: trigs(imax)
!
#ifdef UNICOS
      call set99(trigs,ifax,imax)
#else
      call set999(trigs,ifax,imax)
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
#ifdef MPI
#include "mpif.h"
#endif
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
