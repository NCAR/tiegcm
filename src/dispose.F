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
!   in which the model was run.
!
      implicit none
      character(len=16) :: dispose_script = 'dispose.csh     '
      integer :: lu_dispose
      integer :: lenrec = 240
      contains
!-----------------------------------------------------------------------
      subroutine init_dispose(idispose)
!
! Open and initialize cshell script that will dispose history files
! to mss in separate job step after model execution. This is called
! once per run from input.F.
!
      integer,intent(in) :: idispose
      integer,external :: nextlu,isystem,iunlink
      integer :: istat
      character(len=32) :: command
      logical :: exists
!
! Remove any pre-existing dispose script (so dispose_script is not
! executed if namelist DISPOSE is changed to 0 for a new run)
!
      inquire(file=trim(dispose_script),exist=exists)
      if (exists) then
        istat = iunlink(dispose_script,0)
        if (istat == 0) then
          write(6,"('Dispose: removed pre-existing script ',a)")
     |      trim(dispose_script)
        else
          write(6,"('>>> Dispose: error removing pre-existing',
     |      ' script ',a)") trim(dispose_script)
        endif
      endif 
!
! Make new dispose script if dispose flag is set to 2:
      if (idispose==2) then 
        lu_dispose = nextlu()
        open(file=dispose_script,unit=lu_dispose,status='REPLACE')
        write(lu_dispose,"('#!/bin/csh')")
        close(lu_dispose)   
        write(command,"('chmod u+x ',a)") trim(dispose_script)
        istat = isystem(command)
        write(6,"('Dispose: initialized mss dispose script ',a,
     |    ' lu_dispose=',i3)") trim(dispose_script),lu_dispose
      endif
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
          goto 101
        endif
      endif
      i = i+1
      goto 100
 101  close(lu_dispose) 
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
      subroutine add_dispose(mspath,diskfile,wrpass,mscomment,msreten)
!
! Add mswrite and mscomment lines for the current dispose request to the 
!   dispose script. This is called from sub savefile (output.F).
!
! Args:
      character(len=*),intent(in) :: mspath,diskfile,wrpass
      character(len=*),intent(in) :: mscomment
      integer,intent(in) :: msreten
!
! Local:
      integer :: i,istat,iprint=1
      character(len=120) :: opts,reten
      character(len=lenrec) :: msrcp_cmd, mscomment_cmd
      character(len=8) :: wpass
      logical :: exists
!
! Print warning if diskfile does not exist:
      inquire(file=trim(diskfile),exist=exists)
      if (.not.exists) write(6,"(/,'>>> WARNING add_dispose: ',
     |  ' disk file ',a,' does not exist.')") trim(diskfile)
!
! Check password:
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
      endif
!     write(6,"('add_dispose: wpass=',a)") wpass
!
! Construct mscomment command:
      if (len_trim(mscomment) > 0) then
        if (len_trim(wpass) > 0) then
          write(mscomment_cmd,"('mscomment -wpwd ',a,' -c')")trim(wpass)
          mscomment_cmd = trim(mscomment_cmd)//' "'//trim(mscomment)//
     |      '" '//trim(mspath) 
        else
          mscomment_cmd = 'mscomment -c "'//trim(mscomment)//'" '
     |      //trim(mspath) 
        endif
      endif
!     write(6,"('add_dispose: mscomment_cmd=',a)") trim(mscomment_cmd)
!
! Check if file is already in dispose script:
      i = line_in_script('msrcp',diskfile)
      if (i > 0) goto 100 ! msrcp for this file is already in script
!
! Construct msrcp command:
      opts = ' '
!
! Write password:
      write(opts,"('-wpwd ',a)") trim(wpass)
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
      msrcp_cmd = 'msrcp '//trim(opts)//' '//trim(diskfile)//' mss:'
     |  //trim(mspath)
!
! Write msrcp command to dispose script:
!
      open(file=dispose_script,unit=lu_dispose,status='OLD',
     |  position='APPEND')
      write(lu_dispose,"(' ')") ! space for readability
      write(lu_dispose,"(a)") trim(msrcp_cmd)
      close(lu_dispose)
!     write(6,"('Dispose: added msrcp of file ',a,
!    |  ' to dispose script.')") trim(diskfile)
!     write(6,"('  msrcp_cmd=',a)") trim(msrcp_cmd)
!
! Jump to here if msrcp command was present in script for current file.
 100  continue
!
! If mscomment line is already in script for this file, update it:
      if (len_trim(mscomment) > 0) then
        i = line_in_script('mscomment',mspath)
        if (i > 0) then
          call update_mscomment(i,mscomment_cmd)
!         write(6,"('Dispose: updated mscomment on dispose script',
!    |      ' for file ',a)") trim(diskfile)
!         write(6,"('  ',a)") trim(mscomment_cmd)
!
! Otherwise (mscomment line not present), add it for current file:
        else
          open(file=dispose_script,unit=lu_dispose,status='OLD',
     |      position='APPEND')
          write(lu_dispose,"(a)") trim(mscomment_cmd)
          close(lu_dispose)
!         write(6,"('Dispose: added mscomment on dispose script',
!    |      ' for file ',a)") trim(diskfile)
!         write(6,"('  ',a)") trim(mscomment_cmd)
        endif
      endif ! len_trim(mscomment) > 0
      end subroutine add_dispose      
      end module dispose_module
