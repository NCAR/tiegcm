!
      subroutine getms(mspath,diskfile,tmpdir,lnkfile)
#ifdef MPI
      use mpi_module,only: mytid
#endif
      implicit none
c
c Get file mspath from mss to diskfile, saving a link in tmpdir.
c If the file already exists in tmpdir, link it to diskfile rather than
c   acquiring from the mss (see procedure below).
c
c On input (all character data types):
c   mspath is full mss path and file name of needed file.
c   diskfile is the disk file name (not a path) needed in the cwd.
c   tmpdir is a temporary directory in which the file is to be
c     searched for and saved (linked) to for possible future use.
c   lnkfile is "in use" if it contains non-blank characters,
c     otherwise it is "not in use". If lnkfile is in use, then it is the
c     file name to use in linking to/from tmpdir.
c On output:
c   diskfile is in cwd, and a link has been saved to tmpdir/lnkfile
c     (if lnkfile is in use), or tmpdir/diskfile (if lnkfile not in use)
c
c Procedure:
c   1. If any of the first 3 arguments are all blanks, stop with err msg
c   2. If directory tmpdir does not exist, create it with global permissions
c   3. If diskfile exists (cwd) AND lnkfile NOT in use, then link diskfile to
c      tmpdir/diskfile and return.
c   4. If diskfile does not exist in cwd OR lnkfile IS in use, then:
c      a) if lnkfile is NOT in use, check for existence of tmpdir/diskfile,
c         and if it exists, link it to diskfile and return
c      b) if lnkfile IS in use, check for existence of tmpdir/lnkfile,
c         and if it exists, link it to diskfile and return
c      c) if neither a nor b above are satisfied, then acquire mspath from
c         the mss to diskfile, and link to tmpdir/diskfile (if lnkfile not in
c         use), or to tmpdir/lnkfile (if lnkfile is in use)
c
c Notes:
c   The use of lnkfile simply allows the file to be looked for and
c     saved in tmpdir with a different name than diskfile. If lnkfile
c     is in use, any pre-existing diskfile is not used, and is removed
c     before msrcp if msrcp is necessary.
c   No special options are used in the msrcp (no passwords, etc).
c   Ishell is used only if tmpdir has to be created; all other commands
c     are via fortran calls.
c   If the file is linked from tmpdir, any pre-existing diskfile is
c     removed before the link is made.
c   If a link to the file is saved in tmpdir, any pre-existing file there
c     is first removed before linking from cwd.
c   A warning is printed if any links fail; all other errors are fatal
c     and program is stopped with error message to stdout.
c   Length of (number of non-blank chars) tmpdir+diskfile+1 or
c     tmpdir+lnkfile+1 (if lnkfile is in use) must be <= mxfilelen
c
c B. Foster 2/95
! 10/99: cpp conditionals for UNICOS, SGI, or AIX 
! 11/01: Using call msrcp (libdcs.a) instead of msread.
!
! Args:
      character(len=*),intent(in) :: mspath,diskfile,tmpdir,lnkfile
!
! Local:
      integer,parameter :: mxfilelen = 240 ! maximum number chars in file name
      logical exists
      character(len=1024) :: dcs_output
!
! tmpfile = tmpdir/diskfile or tmpdir/lnkdir
      character(len=mxfilelen) :: tmpfile,cmd 
      integer :: istat,ier
!
! External:
! (ilink, iunlink, and isystem are in util.F)
      integer,external :: ilink,iunlink,isystem
!
#ifdef MPI
!
! If an MPI job, only master task should execute this routine, and
! slave tasks must wait. Mytid is use-associated from the mpi module. 
!
#include "mpif.h"
      if (mytid /= 0) goto 100
!     write(6,"('getms: MPI master task: mspath = ',a)") trim(mspath) 
#endif 
c
c If any of the first 3 arguments have zero length, stop with error msg:
c
      write(6,"(' ',/72('-'),/'GETMS:')")
!     write(6,"('enter getms: len_trim(diskfile)=',i3,
!    |  ' trim(diskfile)=$',a,'$')") len_trim(diskfile),trim(diskfile)

!     subroutine getms(mspath,diskfile,tmpdir,lnkfile)
!     write(6,"('  mspath  = ',a)") trim(mspath)
!     write(6,"('  diskfile = ',a)") trim(diskfile)
!     write(6,"('  tmpdir  = ',a)") trim(tmpdir)
!     write(6,"('  lnkfile = ',a)") trim(lnkfile)

      if (len_trim(mspath).le.0) then
        write(6,"('>>> Bad mspath=',a)") mspath
        write(6,"(72('-'),/)")
        stop 'getms'
      endif
      if (len_trim(diskfile).le.0) then
        write(6,"('>>> Bad diskfile=',a)") diskfile
        write(6,"(72('-'),/)")
        stop 'getms'
      endif
      if (len_trim(tmpdir).le.0) then
        write(6,"('>>> Bad tmpdir=',a)") tmpdir
        write(6,"(72('-'),/)")
        stop 'getms'
      endif
#if (MSS == 0)
      call get_diskfile(diskfile,tmpdir) 
      goto 100
#endif
c
c Construct full tmp path with file name (tmpfile = tmpdir/diskfile, if
c   lnkfile is zero-length, tmpdir/lnkfile otherwise):
c (tmpfile may be up to mxfilelen characters long)
c
      if (len_trim(lnkfile).gt.0) then
        if (len_trim(tmpdir)+len_trim(lnkfile)+1.gt.mxfilelen) then
          write(6,"('>>> WARNING: lengths of tmpdir+lnkfile too ',
     |      'long (sum must be <= 319): len_trim(tmpdir)=',i3,
     |      ' len_trim(lnkfile)=',i3)") len_trim(tmpdir),
     |      len_trim(lnkfile)
        endif
        write(tmpfile,"(a,'/',a)") trim(tmpdir),trim(lnkfile)
      else
        if (len_trim(tmpdir)+len_trim(diskfile)+1.gt.mxfilelen) then
          write(6,"('>>> WARNING: lengths of tmpdir+diskfile too ',
     |      'long (sum must be <= 319): len_trim(tmpdir)=',i3,
     |      ' len_trim(diskfile)=',i3)") 
     |      len_trim(tmpdir),len_trim(diskfile)
        endif
        write(tmpfile,"(a,'/',a)") trim(tmpdir),trim(diskfile)
      endif
c
c If tmpdir does not exist, create it with global permissions:
c
      inquire(file=tmpdir,exist=exists)
      if (.not.exists) then
        cmd = ' '
        write(cmd,"('mkdir -p -m 777 ',a)") trim(tmpdir)
        write(6,"('getms: mkdir command: ',a)") trim(cmd)
!
! isystem is local (util.F), and uses platform dependent escape
! functions to the shell (for AIX, isystem adds null-termination
! to the command string).
! 
        istat = isystem(trim(cmd))
        if (istat.eq.0) then
          write(6,"('Created ',a,' with global permissions')")
     +      trim(tmpdir)
        else
          write(6,"('>>> WARNING: error creating tmpdir ',a,
     +      ' istat=',i3)") trim(tmpdir),istat
        endif
      endif
c
c If file already exists in cwd and lnkfile not in use, save a link
c   in tmpdir (remove any pre-existing tmpfile):
c   (tmpfile = tmpdir/diskfile in this case)
c
      inquire(file=diskfile,exist=exists)
      if (exists.and.len_trim(lnkfile).le.0) then
        write(6,"('Found file ',a,' in cwd')") trim(diskfile)
!       istat = iunlink(tmpfile)
!       istat = ilink(diskfile,tmpfile)
c
c If file not in cwd (or lnkfile was passed in), check in tmpdir:
c (in this case tmpfile=tmpdir/lnkfile OR tmpdir/diskfile)
c
      else
        inquire(file=tmpfile,exist=exists)
c
c If it exists in tmpdir, link it to cwd:
c
        if (exists) then
          istat = iunlink(diskfile,0)
          istat = ilink(tmpfile,diskfile,1)
c
c If not in cwd or tmpdir, get from mss to cwd:
c (if lnkfile is in use, diskfile may already exist in cwd, 
c  so remove it first)
c
        else
          write(6,"('Doing msrcp of mss file ',a,' to disk file ',
     +      a,'...')") trim(mspath),trim(diskfile)
          istat = iunlink(diskfile,0)
!
! Callable dcs commands are in libdcs:
!
#if defined(SUN) || defined(IRIX)
          call msrcp(ier,'mss:'//trim(mspath)//' '//trim(diskfile),
     |      dcs_output)
#else
      call msrcp_(ier,'mss:'//trim(mspath)//' '//trim(diskfile),
     |  dcs_output)
#endif
!
          if (ier.ne.0) then
            write(6,"('>>> Error from msrcp reading mspath ',a,
     +        ' to diskfile ',a)") trim(mspath),trim(diskfile)
            write(6,"('    ier=',i4,' dcs_output=',a)") 
     |        ier,trim(dcs_output)
            stop 'getms'
c
c If msrcp was successful, save a link in tmpdir:
c (remove any pre-existing copy in tmpdir)
c
          else
            write(6,"('Acquired mss file ',a,' to disk file ',
     +        a)") trim(mspath),trim(diskfile)
            istat = iunlink(tmpfile,0)
            istat = ilink(diskfile,tmpfile,1)
          endif ! msrcp return
        endif   ! is or is not in tmpdir
      endif     ! is or is not in cwd
      write(6,"(72('-'),/)")
#ifdef MPI
!
! If mpi job, slave tasks jump down here and wait for master proc:
 100  continue
      call mpi_barrier(MPI_COMM_WORLD,ier)
!     if (mytid==0) then
!       write(6,"('getms MPI: master passed barrier after ',
!    |    ' acquiring mspath ',a)") trim(mspath)
!     else
!       write(6,"('getms MPI: task ',i2,' passed barrier...')") mytid 
!     endif
#endif
      return
      end
!-----------------------------------------------------------------------
      subroutine get_diskfile(diskfile,tmpdir) 
!
! This routine is called if mss is not available (MSS==0).
! Check for existence of file diskfile in cwd. If not there,
! check for existence of tmpdir/diskfile. If not there either
! stop with error message. If tmpdir/diskfile does exist, copy
! it to the cwd. 
!
      implicit none
      character(len=*) :: diskfile,tmpdir
      character(len=80) :: path,cmd
      logical :: exists
      integer :: istat
      integer,external :: isystem
!
! diskfile must be non-null:
      if (len_trim(diskfile).le.0) then
        write(6,"('>>> get_diskfile: Bad diskfile=',a)") diskfile
        stop 'get_diskfile'
      endif
!
! Check for existence of diskfile in cwd:
      inquire(file=diskfile,exist=exists)
      if (exists) then
        write(6,"('get_diskfile: found ',a,' in cwd.')") trim(diskfile)
        return
      else
        write(6,"('get_diskfile: file ',a,' not found.')") 
     |    trim(diskfile)
      endif
!
! Check for existence of diskfile in tmpdir:
!     if (len_trim(tmpdir).le.0) then
!       write(6,"('>>> get_diskfile: Bad tmpdir=',a)") tmpdir
!       stop 'get_diskfile'
!     endif
      if (len_trim(diskfile)+len_trim(tmpdir) > 800) 
     |  write(6,"('>>> WARNING get_diskfile: length of diskfile+tmpdir',
     |    ' > 800')") 
      if (len_trim(tmpdir) > 0) then
        path = trim(tmpdir)//'/'//trim(diskfile) 
      else
        path = trim(diskfile)
      endif
      inquire(file=path,exist=exists)
      if (.not.exists) then
        write(6,"(/,'>>> get_diskfile: could not find file ',a,/)") 
     |    trim(path)
        stop 'get_diskfile'
      else
        write(cmd,"('cp ',a,' .')") trim(path)
        istat = isystem(trim(cmd))
        if (istat.eq.0) then
          write(6,"('get_diskfile: copied ',a,' to cwd.')") trim(path)
        else
          write(6,"('>>> WARNING: error copying ',a,
     +      ' to cwd: istat=',i3)") trim(path),istat
        endif
      endif
      end subroutine get_diskfile
