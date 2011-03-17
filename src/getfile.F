!
      subroutine getfile(remote,local)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      implicit none
!
! Args:
      character(len=*),intent(in) :: 
     |  remote    ! file to get (may be relative or absolute path)
      character(len=*),intent(out) :: 
     |  local     ! local disk file name returned
!
! Local:
      integer :: istat
      logical :: exists
      character(len=1024) :: cmd
!
! External (util.F):
      integer,external :: isystem
!
      write(6,"(/,72('-'))")
!
! Remote must be non-blank:
      if (len_trim(remote) <= 0) then
        write(6,"(/,'>>> Getfile: file argument ''remote'' is blank.')")
        call shutdown('Getfile: blank file argument')
      endif
      local = ' '
      write(6,"('Getfile: remote=',a)") trim(remote)
!
! expand any env vars imbedded in remote (util.F)
      call expand_path(remote) 
!
! Check for existence of remote:
      inquire(file=trim(remote),exist=exists)
!
! If remote exists assign local <- remote and return, otherwise
! shutdown with fatal error:
      if (exists) then
        write(6,"('Getfile: Found file ',a)") trim(remote)
        local = remote
      else  ! remote does not exist
        write(6,"(/,'>>> Getfile: could not find file ',a)")
     |    trim(remote)
        call shutdown('Getfile: file not found')
      endif ! existence of remote
      write(6,"(72('-'),/)")
      end subroutine getfile
