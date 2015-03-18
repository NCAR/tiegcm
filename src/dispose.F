!
      module dispose_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Maintain an external csh script containing the hsi command necessary
! to dispose history output files to the NCAR HPSS. The script name is
! dispose_xxxxxx.hsi, where xxxxxx is a PID unique to the model run.
!
      use input_module,only: hpss_path,pid
      use init_module,only: logname
      implicit none
      character(len=32) :: dispose_script
      integer :: lu_dispose
      integer :: lenrec = 1024
      contains
!-----------------------------------------------------------------------
      subroutine init_dispose
!
! Open and initialize cshell script that will dispose history files
! to hpss in separate job step after model execution. This is called
! once per run from input.F.
!
      integer,external :: nextlu,isystem
      integer :: istat
      character(len=32) :: command
!
! Make dispose script name from pid (see getprocessid in util.F):
      write(dispose_script,"('dispose_',i8.8,'.hsi')") pid
!     write(6,"('init_dispose: pid=',i8,' dispose_script=',a)") 
!    |  pid,trim(dispose_script)
!
! Open new file and initialize as csh script:
      lu_dispose = nextlu()
      open(file=dispose_script,unit=lu_dispose,status='REPLACE')
      write(lu_dispose,"('#!/bin/csh')")
      close(lu_dispose)   
      write(command,"('chmod u+x ',a)") trim(dispose_script)
      istat = isystem(command)
      write(6,"('Dispose: initialized hpss dispose script ',a,
     |  ' lu_dispose=',i3)") trim(dispose_script),lu_dispose
      end subroutine init_dispose
!-----------------------------------------------------------------------
      integer function isput_in_script(putfile,script)
!
! Args:
      character(len=*),intent(in) :: putfile,script
!
! Local:
      character(len=lenrec) :: line
      integer :: i,ipos
!
      isput_in_script = 0
      open(file=script,unit=lu_dispose,status='OLD',action='READ')
      i = 1
 100  continue
      read(lu_dispose,"(a)",end=900) line 
      ipos = index(line,'hsi put')
      if (ipos > 0 .and. index(line,trim(putfile)) > 0) then
        isput_in_script = i
        goto 101
      endif
      i = i+1
      goto 100
 900  continue ! end of file
 101  close(lu_dispose)
      end function isput_in_script
!-----------------------------------------------------------------------
      subroutine update_annotation(line,annotate_cmd)
!
! Args:
      integer,intent(in) :: line              ! line number to update
      character(len=*),intent(in) :: annotate_cmd
!
! Local:
      character(len=lenrec) :: rdline,tmpfile
      character(len=120) :: command
      integer :: i,lutmp,istat
      logical :: wrote
!
! External:
      integer,external :: nextlu,isystem
!
      open(file=dispose_script,unit=lu_dispose,status='OLD',
     |  action='READ')
      lutmp = nextlu()
      tmpfile = trim(dispose_script)//'.tmp'
      open(file=tmpfile,unit=lutmp,status='REPLACE')
      i = 1
      wrote = .false.
 100  continue
      if (i/=line) then
        read(lu_dispose,"(a)",end=900) rdline
        write(lutmp,"(a)") trim(rdline)
      else
        read(lu_dispose,"(a)",end=900) rdline
        write(lutmp,"(a)") trim(annotate_cmd)
        wrote = .true.
      endif
      i = i+1
      goto 100
 900  continue
!
! Append to end of file: 
      if (.not.wrote) write(lutmp,"(a)") trim(annotate_cmd)
!
! Close both files, move tmp file to dispose_script file, and make 
!   it executable.
!
      close(lutmp)
      close(lu_dispose)
      write(command,"('mv ',a,' ',a)")trim(tmpfile),trim(dispose_script)
      istat = isystem(command)
      write(command,"('chmod u+x ',a)") trim(dispose_script)
      istat = isystem(command)
      end subroutine update_annotation
!-----------------------------------------------------------------------
      subroutine add_dispose(diskfile,annotation)
!
! Add hsi put and hsi annotation lines for the current dispose request 
!   to the hsi dispose script. This is called from sub savefile (output.F).
!
! Args:
      character(len=*),intent(in) :: diskfile,annotation
!
! Local:
      integer :: i,iline,istat,iprint=1
      character(len=1020) :: hpss_file,string
      character(len=8) :: opt
      character(len=lenrec) :: hsiput_cmd, annotate_cmd
      logical :: exists

!     write(6,"('Enter add_dispose: diskfile=',a,' annotation=',
!    |  a)") trim(diskfile),trim(annotation)
!
! Print warning if diskfile does not exist:
      inquire(file=trim(diskfile),exist=exists)
      if (.not.exists) write(6,"(/,'>>> WARNING add_dispose: ',
     |  ' disk file ',a,' does not exist.')") trim(diskfile)
!
! Check if put command for the current file is already in dispose script:
      iline = isput_in_script(diskfile,dispose_script)
!
! Construct hsi put command (with or without -P depending on hpss_path):
      hpss_file = ' '
      write(hpss_file,"(a,'/',a)") trim(hpss_path),trim(diskfile)
      string = ' ' ; write(string,"('/home/',a)") trim(logname)
      opt = ' '
      if (trim(string)==trim(hpss_path)) then ! hpss_path does not have subdir
        write(hsiput_cmd,"('hsi put ',a,' : ',a)")
     |    trim(diskfile),trim(hpss_file)
      else                                    ! hpss_path has subdir
        write(hsiput_cmd,"('hsi put -P ',a,' : ',a)")
     |    trim(diskfile),trim(hpss_file)
        opt = '-P'
      endif
!
! Write hsi put command to dispose script:
      if (iline==0) then ! hsi put was not found for this file
        open(file=dispose_script,unit=lu_dispose,status='OLD',
     |    position='APPEND')
        write(lu_dispose,"(' ')") ! space for readability
        write(lu_dispose,"(a)") trim(hsiput_cmd)
        close(lu_dispose)
      endif
!
! Check again for the put command, and if found, add annotation command:
      iline = isput_in_script(diskfile,dispose_script)
      if (iline > 0) then
!
! Construct annotation command:
        annotate_cmd = ' '
        annotate_cmd = 'hsi annotate -A'
!
! hsi requires that double quotes around the annotation be escaped w/ a
! single backslash in the hsi script.  xlf90 and PGI compilers require 
! escape of the backslash, ifort does not.
!
#if defined(AIX) || defined(PGI)
        annotate_cmd =
     |    trim(annotate_cmd)//' '//'\\"'//trim(annotation)//'\\"'
#else
        annotate_cmd =
     |    trim(annotate_cmd)//' '//'\"'//trim(annotation)//'\"'
#endif
        annotate_cmd = trim(annotate_cmd)//' '//trim(hpss_file)
!       write(6,"('add_dispose: annotate_cmd = ',a)") 
!    |    trim(annotate_cmd)
!
! Update annotation for this file:
        call update_annotation(iline+1,trim(annotate_cmd))
      endif
      end subroutine add_dispose      
      end module dispose_module
