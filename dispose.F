!
      module dispose_module
!
! DISPOSE=2 in the model namelist input file means the user wants to dispose 
!   history files *after* model execution. In this case, sub input calls sub
!   init_dispose, which opens a csh script "dispose.csh" in the cwd which 
!   will be used to dispose the files to the mss after model execution.
! When it is time to save a file to the mss, sub savefile calls sub
!   add_dispose, which adds two commands to the dispose script: an msrcp 
!   command to dispose the file to mss, and an mscomment to add the fileinfo 
!   comment to the mss file. If the msrcp command is already in the script 
!   for the current file, it is left alone. If the mscomment command is 
!   already there for the current file, it is updated with the current comment. 
!   Thus, in the final dispose script, each file should have a single msrcp 
!   and a single mscomment command.
! Following model execution, the dispose script can be executed in an 
!   interactive loadleveler job step following the model execution step, 
!   or it can be executed manually on the command line from the directory
!   in which the model was run (assuming POSTCLEAN was 0).
!
      implicit none
      character(len=16) :: dispose_script = 'dispose.csh     '
      integer :: lu_dispose
      integer :: lenrec = 240
      contains
!-----------------------------------------------------------------------
      subroutine init_dispose
!
! Open and initialize cshell script that will dispose history files
! to mss in separate job step after model execution. This is called
! once per run from input.F.
!
      integer,external :: nextlu,isystem
      integer :: istat
      character(len=32) :: command
!
      lu_dispose = nextlu()
      open(file=dispose_script,unit=lu_dispose,status='REPLACE')
      write(lu_dispose,"('#!/bin/csh')")
      close(lu_dispose)
      write(command,"('chmod u+x ',a)") trim(dispose_script)
      istat = isystem(command)
      write(6,"('Dispose: initialized mss dispose script ',a,
     |  ' lu_dispose=',i3)") trim(dispose_script),lu_dispose
      end subroutine init_dispose
!-----------------------------------------------------------------------
      integer function line_in_script(mscommand,file)
!
! Args:
      character(len=*),intent(in) :: mscommand,file
!
! Local:
      character(len=lenrec) :: line
      integer :: i,ipos
!
      line_in_script = 0
      open(file=dispose_script,unit=lu_dispose,status='OLD',
     |  action='READ')
      i = 1
 100  continue
      read(lu_dispose,"(a)",end=900) line 
      ipos = index(line,mscommand)
      if (ipos > 0) then
        ipos = index(line,trim(file))
        if (ipos > 0) then
          line_in_script = i
          close(lu_dispose) 
          return
        endif
      endif
      i = i+1
      goto 100
      close(lu_dispose) 
      return
 900  continue ! end of file
      close(lu_dispose) 
      end function line_in_script
!-----------------------------------------------------------------------
      subroutine update_mscomment(line,comment_command)
!
! Args:
      integer,intent(in) :: line              ! line number to update
      character(len=*),intent(in) :: comment_command  ! mscomment command
!
! Local:
      character(len=16) :: tmpfile = 'dispose_tmp.csh '
      character(len=lenrec) :: rdline
      character(len=120) :: command
      integer :: i,lutmp,istat
!
! External:
      integer,external :: nextlu,isystem
!
      open(file=dispose_script,unit=lu_dispose,status='OLD',
     |  action='READ')
      lutmp = nextlu()
      open(file=tmpfile,unit=lutmp,status='REPLACE')
      i = 1
 100  continue
      if (i/=line) then
        read(lu_dispose,"(a)",end=900) rdline
        write(lutmp,"(a)") trim(rdline)
      else
        read(lu_dispose,"(a)",end=900) rdline
        write(lutmp,"(a)") trim(comment_command)
      endif
      i = i+1
      goto 100
 900  continue
!
! Close both files, move tmp file to dispose_script file, and make 
!   it executable.
      close(lutmp)
      close(lu_dispose)
      write(command,"('mv ',a,' ',a)")trim(tmpfile),trim(dispose_script)
      istat = isystem(command)
      write(command,"('chmod u+x ',a)") trim(dispose_script)
      istat = isystem(command)
      end subroutine update_mscomment
!-----------------------------------------------------------------------
      subroutine add_dispose(mspath,diskfile,tmpdir,lnkfile,wrpass,
     |  mscomment,msreten)
!
! Add mswrite and mscomment lines for the current dispose request to the 
!   dispose script. This is called from sub savefile (output.F).
!
! Args:
      character(len=*),intent(in) ::
     |  mspath,diskfile,tmpdir,lnkfile,wrpass
      character(len=*),intent(in) :: mscomment
      integer,intent(in) :: msreten
!
! Local:
      integer :: i,istat,iprint=1
      character(len=120) :: opts,reten
      character(len=lenrec) :: command
      character(len=8) :: wpass
!
! External:
      integer,external :: ilink,iunlink
!
! Link file to tmpdir:
      if (trim(tmpdir) /= '.') then
        istat = iunlink(trim(tmpdir)//'/'//trim(diskfile),iprint)
        istat = ilink(diskfile,trim(tmpdir)//'/'//trim(diskfile),1)
        if (istat /= 0) 
     |    write(6,"('>>> WARNING dispose: error ',i3,' from ilink.')") 
     |      istat
      else
        write(6,"('>>> WARNING dispose: NOT linking file ',a,
     |    ' to tmpdir because tmpdir==cwd')") trim(diskfile)
      endif
!
! Check if file is already in dispose script:
      i = line_in_script('msrcp',diskfile)
      if (i > 0) goto 100 ! msrcp for this file is already in script
!
! Construct msrcp command:
      opts = ' '
!
! Write password:
      wpass = ' '
      if (len_trim(wrpass) > 0) then
        if (len_trim(wrpass) > 8) then
          write(6,"('>>> WARNING: only first 8-chars of wrpass ',
     |      'will be used',/'    Write password will be: ',a)")
     |      wrpass(1:8)
          wpass = wrpass(1:8)
        else
          wpass = trim(wrpass)
        endif
        write(opts,"('-wpwd ',a)") trim(wpass)
      endif
!
! Retention period:
      if (msreten > 0) then
        if (msreten <= 9) then
          write(reten,"(' -period ',i1)") msreten
        elseif (msreten <= 99) then
          write(reten,"(' -period ',i2)") msreten
        elseif (msreten <= 999) then
          write(reten,"(' -period ',i3)") msreten
        elseif (msreten <= 9999) then
          write(reten,"(' -period ',i4)") msreten
        else
          write(6,"('WARNING dispose: retention period ',i8,
     |      ' is too big -- am defaulting to 365 days.')")
          write(reten,"('-period 365')")
        endif
        opts = trim(opts)//trim(reten)
      endif
      command = 'msrcp '//trim(opts)//' '//trim(diskfile)//' mss:'
     |  //trim(mspath)
!
! Write msrcp command to dispose script:
!
      open(file=dispose_script,unit=lu_dispose,status='OLD',
     |  position='APPEND')
      write(lu_dispose,"(' ')") ! space for readability
      write(lu_dispose,"(a)") trim(command)
      close(lu_dispose)
      write(6,"('Dispose: added msrcp of file ',a,
     |  ' to dispose script.')") trim(diskfile)
!
! Jump to here if msrcp command was present in script for current file.
 100  continue
!
! Construct mscomment command:
      if (len_trim(mscomment) > 0) then
        if (len_trim(wpass) > 0) then
          write(command,"('mscomment -wpwd ',a,' -c')") trim(wpass)
          command = trim(command)//' "'//trim(mscomment)//'" '
     |      //trim(mspath) 
        else
          command = 'mscomment -c "'//trim(mscomment)//'" '
     |      //trim(mspath) 
        endif
!
! If mscomment line is already in script for this file, update it:
        i = line_in_script('mscomment',mspath)
        if (i > 0) then
          call update_mscomment(i,command)
          write(6,"('Dispose: updated mscomment on dispose script',
     |      ' for file ',a)") trim(diskfile)
!
! Otherwise (mscomment line not present), add it for current file:
        else
          open(file=dispose_script,unit=lu_dispose,status='OLD',
     |      position='APPEND')
          write(lu_dispose,"(a)") trim(command)
          close(lu_dispose)
          write(6,"('Dispose: added mscomment on dispose script',
     |      ' for file ',a)") trim(diskfile)
        endif
      endif ! len_trim(mscomment) > 0
      end subroutine add_dispose      
      end module dispose_module
