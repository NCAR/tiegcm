!
      subroutine putms(mspath,dskfile,tmpdir,wrpass,mscomment,
     |  msreten)
      use input_module,only: cwd ! current working directory
      implicit none
!
! Called by sub savefile (output.f):
!       call putms(outpath,diskfile,tempdir,' ','NCARTGCM',' ')
! Called by master task 0 only.
!
! Dispose existing disk file dskfile to mss as mspath, with optional
!   write password. Also save a link of the file in tmpdir.
!   (if mspath arg is all blanks, do not do the msrcp, i.e., only
!    save link to tmpdir)
!
! On input:
!   character mspath is full mss path and file name to which dskfile
!     is to be disposed
!   character dskfile is existing disk file in cwd
!   character tmpdir is directory in which to save a link to dskfile
!   character wrpass, if non-zero length, is write password for mspath
!   integer msreten is the retention period for the mss file (days)
! On output:
!   dskfile has been disposed to mspath, optionally with write password
!   A link to dskfile has been saved in tmpdir
!
! Notes:
!   Non-existence of dskfile is a fatal error.
!   System call is made only if tmpdir has to be created; all other commands
!     are via fortran calls.
!   If wrpass has non-zero length, then it is given as a write password
!     on the mss file.
!   By default, the file is disposed to mss in transparent format
!   If mspath already exists on the mss, it will be overwritten,
!     provided the write password given (wrpass(1:8)) is correct.
!   Errors in links or msrcp results in warning message, but is
!     non-fatal.
!   If tmpdir does not exist, it is created with global permissions
!     prior to linking.
!
! B. Foster 2/95
! B. Foster 6/04: Rewritten, using only msrcp (eliminate mswrite).
! B. Foster 7/13/04: 
!   Removed msasync option (async dispose causes problems with mscomment, 
!   and does not work on all platforms. Use DISPOSE=2 instead, to defer
!   mss file disposes until after model execution).
!
! Args:
      character(len=*),intent(in) :: 
     |  mspath,dskfile,tmpdir,wrpass,mscomment
      integer,intent(in) :: msreten
!
! Local:
      integer,parameter :: mxfilelen = 240 ! maximum number chars in file name
      logical exists
      character(len=520) :: opts,testfile,cmd,tempfile
      character(len=8),save :: curdate,curtime
      character(len=len(tmpdir)) :: tempdir
      integer :: ier,istat
!
! External:
      integer,external :: ilink,iunlink,isystem,my_msrcp
!
      write(6,"(' ',/72('-'),/'PUTMS:')")

!     write(6,"('  mspath  = ',a)") trim(mspath) 
!     write(6,"('  dskfile = ',a)") trim(dskfile) 
!     write(6,"('  tmpdir  = ',a)") trim(tmpdir) 
!     write(6,"('  mscomment = ',a)") trim(mscomment) 
!     write(6,"('  msreten = ',i3)") msreten
!     write(6,"('  wrpass  = ',a)") trim(wrpass) 

!
! dskfile must be non-blank:
!
      if (len_trim(dskfile) < 0) then
        write(6,"('>>> putms: dskfile must be non-blank.')")
        write(6,"(72('-'),/)")
        call shutdown('putms')
      endif
!
! If dskfile does not exist, stop with fatal error message:
!
      inquire(file=dskfile,exist=exists)
      if (.not.exists) then
        write(6,"('>>> Cannot find file ',a)") trim(dskfile)
        write(6,"(72('-'),/)")
        call shutdown('putms')
      endif
!
! Transfer tmpdir to local var:
      tempdir = tmpdir
!
! If tempdir is empty or same as cwd, make tempdir cwd:
      if (len_trim(tmpdir) <= 0.or.trim(tempdir)==trim(cwd)) then 
        tempdir = ' '
        tempdir(1:1) = '.'
!       write(6,"('Note putms: using cwd as tmpdir..')")
      endif
!
! If tmpdir does not exist, create it with global permissions:
!
      if (trim(tempdir) /= '.') then ! tempdir is not cwd
#if defined(OSF1) 
!
! Inquire(file=name,exist=exists) statement under OSF always returns false
! if "name" is a directory (presumably it must be a file name). Therefore,
! to test existence of the directory tmpdir we attempt to create a file in 
! tmpdir with the "touch" command, and then inquire on that:
!
        write(testfile,"(a,'/test.file')") trim(tempdir)
        cmd = ' '
        write(cmd,"('touch ',a)") trim(testfile)
        istat = isystem(trim(cmd)) 
        inquire(file=trim(testfile),exist=exists)
#else
        inquire(file=tempdir,exist=exists)
#endif
        if (.not.exists) then
          cmd = ' '
          write(cmd,"('mkdir -m 777 ',a)") trim(tempdir)
!
! isystem is local (util.F), and uses platform dependent escape
! functions to the shell (for AIX, isystem adds null-termination
! to the command string).
! 
          istat = isystem(cmd(1:len_trim(cmd)))
          if (istat.eq.0) then
            write(6,"('Created ',a,' with global permissions')")
     +        trim(tempdir)
          else
            write(6,"('>>> WARNING: error creating tempdir',a,
     +        ' istat=',i3)") trim(tempdir)
          endif
        endif
      endif ! tempdir /= '.'
!
! Dispose dskfile to mspath:
!
      if (len_trim(mspath) > 0) then
#ifdef MSS
!
! Build options to msrcp command:
!
! Set mss class-of-service reliability to economy. 
! This means only one copy of the mss file will be written. 
        opts = '-class reliability=econ'
!
! Write password:
        if (len_trim(wrpass) > 0) 
     |    opts = trim(opts)//' -wpwd "'//trim(wrpass)//'"'
!
! Retention period:
        if (msreten > 0) then
          opts = trim(opts)//' -period '
          write(opts(len_trim(opts)+1:len_trim(opts)+4),"(i4)") 
     |      msreten
        endif
!
!       write(6,"('putms: opts for msrcp=',a)") trim(opts)
!
! Execute msrcp command:
!       write(6,"('putms calling my_msrcp to dispose file ',a)")
!    |    trim(mspath) 
        istat = my_msrcp(opts,dskfile,"mss:"//trim(mspath))
!
! If successful, report to stdout:
        if (istat==0) then
          write(6,"('Disposed file ',a,' to mss ',a)") trim(dskfile),
     |      trim(mspath)
          write(6,"('Mss retention period = ',i5)") msreten
!
! Add comment field (build options to mscomment command):
          if (len_trim(mscomment) > 0) then
            opts = ' '
            if (len_trim(wrpass) > 0) opts = trim(opts)//' -wpwd '//
     |        trim(wrpass)
            opts = trim(opts)//' -c "'//trim(mscomment)//'"' 
            cmd = "mscomment "//trim(opts)//' '//trim(mspath)
            write(6,"(a)") trim(cmd)
            istat = isystem(cmd) ! cmd is null-terminated by isystem for aix
            if (istat /= 0) 
     |        write(6,"('>>> WARNING: Error from mscomment command:',
     |          ' istat=',i8)") istat
          endif
!
! my_msrcp returned error:
        else
          write(6,"('>>> WARNING putms: Error return from my_msrcp: ',
     |      'istat=',i10)") istat
        endif
!
! End #ifdef MSS
#endif
      endif ! len_trim(mspath) > 0
!
! Save a link to tmpdir (remove any pre-existing tempdir/dskfile):
!
      if (trim(tempdir) /= '.') then
        tempfile = trim(tempdir)//'/'//trim(dskfile)
        istat = iunlink(tempfile,1)
        istat = ilink(dskfile,tempfile,1)
      endif
      if (len_trim(mspath) <= 0) 
     |  write(6,"('(File NOT disposed to mss)')")
      write(6,"(72('-'),/)")
      end subroutine putms
