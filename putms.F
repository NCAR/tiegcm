!
      subroutine putms(mspath,dskfile,tmpdir,lnkfile,wrpass,mscomment,
     |  msreten,msasync)
      implicit none
!
! Called by sub savefile (output.f):
!       call putms(outpath,diskfile,tempdir,' ','NCARTGCM',' ',msasync)
c
c Dispose existing disk file dskfile to mss as mspath, with optional
c   write password. Also save a link of the file in tmpdir.
c   (if mspath arg is all blanks, do not do the mswrite, i.e., only
c    save link to tmpdir)
c
c On input:
c   character mspath is full mss path and file name to which dskfile
c     is to be disposed
c   character dskfile is existing disk file in cwd
c   character tmpdir is directory in which to save a link to dskfile
c   character wrpass, if non-zero length, is write password for mspath
c   integer msreten is the retention period for the mss file (days)
c   If msasync > 0, asynchronous dispose will be requested.
c On output:
c   dskfile has been disposed to mspath, optionally with write password
c   A link to dskfile has been saved in tmpdir
c
c Notes:
c   Non-existence of dskfile is a fatal error.
c   Ishell is used only if tmpdir has to be created; all other commands
c     are via fortran calls.
c   If wrpass has non-zero length, then it is given as a write password
c     on the mss file. If wrpass length > 8, only the first 8 chars
c     are used (this is an mswrite restriction).
c   By default, the file is disposed to mss in transparent format
c   If mspath already exists on the mss, it will be overwritten,
c     provided the write password given (wrpass(1:8)) is correct.
c   Errors in links or mswrite results in warning message, but is
c     non-fatal.
c   If tmpdir does not exist, it is created with global permissions
c     prior to linking.
c
c B. Foster 2/95
!
! Args:
      character(len=*),intent(in) :: 
     |  mspath,dskfile,tmpdir,lnkfile,wrpass
      character(len=*),intent(in) :: mscomment
      integer,intent(in) :: msreten,msasync
!
! Local:
      integer,parameter :: mxfilelen = 240 ! maximum number chars in file name
      logical exists
      character(len=80) :: errmsg
      character(len=120) :: opts,testfile
      character(len=mxfilelen) :: tmpfile,cmd
      character*10 wpass
      integer :: lmspath,ldskfile,ltmpdir,llnkfile,ltmpfile,
     |  lwrpass,lcomment,ier,istat
      character(len=120),save :: dskfile_prev=' '
      character(len=8),save :: curdate,curtime
!
! External:
      integer,external :: ilink,iunlink,isystem
c
c If any of first 3 arguments have zero length, stop with error msg:
c
      write(6,"(' ',/72('-'),/'PUTMS:')")
      lmspath = len_trim(mspath)
      ldskfile = len_trim(dskfile)
      ltmpdir = len_trim(tmpdir)
      if (ldskfile.le.0) then
        write(6,"('>>> Bad dskfile=',a)") dskfile
        write(6,"(72('-'),/)")
        stop 'putms'
      endif
      if (ltmpdir.le.0) then
        write(6,"('>>> Bad tmpdir=',a)") tmpdir
        write(6,"(72('-'),/)")
        stop 'putms'
      endif
c
c If dskfile does not exist, stop with fatal error message:
c
      inquire(file=dskfile,exist=exists)
      if (.not.exists) then
        write(6,"('>>> Cannot find file ',a)") dskfile(1:ldskfile)
        write(6,"(72('-'),/)")
        stop 'putms'
      endif
c
c Construct full tmp path with file name (tmpfile = tmpdir/dskfile, if
c   lnkfile is zero-length, tmpdir/lnkfile otherwise):
c (tmpfile may be up to mxfilelen characters long)
c
      llnkfile = len_trim(lnkfile)
      if (llnkfile.gt.0) then
        if (ltmpdir+llnkfile+1.gt.mxfilelen) then
          write(6,"('>>> WARNING: lengths of tmpdir+lnkfile too ',
     +      'long (sum must be <= ',i4,'): ltmpdir=',i3,' llnkfile=',
     +      i3)") mxfilelen,ltmpdir,llnkfile
        endif
        write(tmpfile,"(a,'/',a)") tmpdir(1:ltmpdir),lnkfile(1:llnkfile)
      else
        if (ltmpdir+ldskfile+1.gt.mxfilelen) then
          write(6,"('>>> WARNING: lengths of tmpdir+dskfile too ',
     +      'long (sum must be <= ',i4,'): ltmpdir=',i3,' ldskfile=',
     +      i3)") mxfilelen,ltmpdir,ldskfile
        endif
        write(tmpfile,"(a,'/',a)") tmpdir(1:ltmpdir),dskfile(1:ldskfile)
      endif
      ltmpfile = len_trim(tmpfile)
c
c If tmpdir does not exist, create it with global permissions:
c
#if defined(OSF1) 
!
! 12/9/02 bf: Currently, the only OSF system being used is lemieux.psc.edu.
! Lemieux's frontend archive storage is a machine named golem.psc.edu. The
! Below links and rcp attempt do not work from the lemieux compute nodes,
! and may be using unnecessary wallclock time, so am not using them. Storage
! of history output to golem will be via rcp in the job script after model 
! execution (e.g, see job script tiegcm1lx.job)
! 
      write(6,"('PUTMS returning because this is an OSF job.')")
      return
!
! Inquire(file=name,exist=exists) statement under OSF always returns false
! if "name" is a directory (presumably it must be a file name). Therefore,
! to test existence of the directory tmpdir we attempt to create a file in 
! tmpdir with the "touch" command, and then inquire on that:
!
      write(testfile,"(a,'/test.file')") trim(tmpdir)
      cmd = ' '
      write(cmd,"('touch ',a)") trim(testfile)
      istat = isystem(trim(cmd)) 
      inquire(file=trim(testfile),exist=exists)
#else
      inquire(file=tmpdir,exist=exists)
#endif
      if (.not.exists) then
        cmd = ' '
        write(cmd,"('mkdir -m 777 ',a)") tmpdir(1:ltmpdir)
!
! isystem is local (util.F), and uses platform dependent escape
! functions to the shell (for AIX, isystem adds null-termination
! to the command string).
! 
        istat = isystem(cmd(1:len_trim(cmd)))
        if (istat.eq.0) then
          write(6,"('Created ',a,' with global permissions')")
     +      tmpdir(1:ltmpdir)
        else
          write(6,"('>>> WARNING: error creating tmpdir',a,
     +      ' istat=',i3)") tmpdir(1:ltmpdir)
        endif
      endif
c
c Dispose dskfile to mspath:
c
      if (lmspath.gt.0) then
        lwrpass = len_trim(wrpass)
        lcomment = len_trim(mscomment)
        opts = ' '
        if (lcomment.gt.0) 
     +    write(opts,"('CMNT=',a)") mscomment(1:lcomment)
        if (lwrpass.gt.8) then
          write(6,"('>>> WARNING: only first 8-chars of wrpass ',
     +      'will be used',/'    Write password will be: ',a)")
     +      wrpass(1:8)
        endif
#ifdef MSS
!
! If asynchronous dispose is required, wait for previous request
! to complete before initiating current request. Wait only if this is
! not the first request, and if we are disposing the same diskfile as  
! last time. Also set NOMAIL so mail is sent only in case of error.
!
        if (msasync > 0) then
          opts = opts(1:len_trim(opts))//',nowait,nomail'
          if (len_trim(dskfile_prev) > 0.and.          ! not first dispose
     |        trim(dskfile)==trim(dskfile_prev)) then  ! disposing same file
            call datetime(curdate,curtime)
            write(6,"('Begin wait for completion of previous ',
     |        'asynchronous mss write of disk file',/,'  ',a,
     |        ' at ',a,' ',a)") trim(dskfile_prev),curdate,curtime
#ifdef LINUX
            write(6,"('>>> Must set msasync = 0 on Linux systems..')")
            stop 'msasync'
#else
            call msrecall(trim(dskfile_prev))
#endif
            call datetime(curdate,curtime)
            write(6,"('Completed wait at ',a,' ',a)") curdate,curtime
          endif
        endif
!
! Set mss class-of-service reliability to economy. This means only one
! copy of the mss file will be written. The default reliability is
! "standard", meaning 2 copies will be saved. See msclass man page.
!
        opts = opts(1:len_trim(opts))//',reliability=economy'
#if defined(UNICOS) || defined(IRIX) || defined(OSF1) || defined(AIX)
        if (lwrpass.le.0) then
          call mswrite(ier,dskfile(1:ldskfile),mspath(1:lmspath),
     +      ' ',msreten,opts)
        else
          write(wpass,"(' ,',a)") wrpass(1:8)
          call mswrite(ier,dskfile(1:ldskfile),mspath(1:lmspath),
     +      wpass,msreten,opts)
        endif
        if (ier.eq.0) then
          if (msasync > 0) then
            write(6,"('Requested async dispose of file ',a,' to mss',
     |        ' path ',/,'  ',a)") trim(dskfile),trim(mspath)
          else
            write(6,"('Disposed file ',a,' to mss ',a)")
     +        trim(dskfile),trim(mspath)
          endif
          if (lcomment.gt.0)
     +      write(6,"('Mss comment: ',a)") mscomment(1:lcomment)
          write(6,"('Mss retention period is ',i5,' days.')") msreten
        else
          call mserror(errmsg)
          write(6,"('>>> WARNING: error from mswrite of disk file ',
     +      a,' to mspath ',a)") dskfile(1:ldskfile),mspath(1:lmspath)
          if (lwrpass.gt.0) write(6,"('    Write password = ',a)")
     +      wrpass(1:8)
          write(6,"('    Error message: ',a)") errmsg
        endif ! error return from mswrite
#else
!
! 2/00: mswrite appears to be on most machines -- if OS is
!       unresolved, try msrcp:
!
        write(6,"('>>> WARNING putms: unresolved OS cpp directive',
     |    ' -- will try msrcp system call.')")
        cmd = ' '
        cmd = "msrcp "//dskfile(1:ldskfile)//" mss:"//trim(mspath)
        write(6,"('Executing isystem call on AIX: cmd=',/,a)")
     |    trim(cmd)
        istat = isystem(cmd) ! cmd is null-terminated by isystem for aix
        if (istat==0) then
          write(6,"('Disposed file ',a,' to mss ',a)")
     +      trim(dskfile),trim(mspath)
          if (lcomment.gt.0)
     +      write(6,"('Mss comment: ',a)") mscomment(1:lcomment)
          write(6,"('Mss retention period is ',i5,' days.')") msreten
        endif
#endif
!
! MSS not set -- try rcp. For now, assume target is golem.psc.edu
!   (for job running on lemieux.psc.edu)
!
#else
        write(6,"('putms calling rcpfile: dskfile=',a)")
     |    trim(dskfile)
        call rcpfile(dskfile,'golem.psc.edu:'//trim(dskfile))
#endif
      endif ! lmspath > 0
c
c Save a link to tmpdir (first remove any pre-existing tmpfile):
c
      istat = iunlink(tmpfile,0)
      istat = ilink(dskfile,tmpfile,1)
      if (lmspath.le.0) write(6,"('(File NOT disposed to mss)')")
      write(6,"(72('-'),/)")
      if (len_trim(dskfile) > 120) write(6,"('>>> WARNING putms:',
     |  ' dskfile too long for dskfile_prev: len_trim(dskfile)=',
     |  i4)") len_trim(dskfile)
      dskfile_prev(1:len_trim(dskfile)) = trim(dskfile)
      return
      end
