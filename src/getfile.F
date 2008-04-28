!
      subroutine getfile(remote,local)
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
      integer,external :: isystem,my_msrcp
!
      write(6,"(/,72('-'))")
!
! Remote must be non-blank:
      if (len_trim(remote) <= 0) then
        write(6,"(/,'>>> Getfile: file argument ''remote'' is blank.')")
        call shutdown('Getfile: blank file argument')
      endif
      local = ' '
!
      write(6,"('Getfile: remote=',a)") trim(remote)

      call expand_path(remote) ! expand any env vars imbedded in remote (util.F)
!     write(6,"('Getfile after expand_path: remote=',a)") trim(remote)
!
! Check for existence of remote as provided:
      inquire(file=trim(remote),exist=exists)
!
! If remote exists and local is blank (or local==remote), just return.
! If local is non-blank, then link or copy remote to local.
!   (try hard link, symbolic link, and copy, in that order)
!
      if (exists) then
        write(6,"('Getfile: Found file ',a)") trim(remote)
        local = remote
        goto 100
!
! Remote as given does not exist. First check for existence of 
! mkdiskflnm(remote,local). If that does not exist, try obtaining 
! from mss (if MSS macro is set), otherwise give up.
!
      else  ! remote does not exist
        call mkdiskflnm(remote,local) ! in util.F
        inquire(file=trim(local),exist=exists)
        if (exists) then
          write(6,"('Getfile: Found file ',a)") trim(local)
          goto 100
        endif
        if (remote(1:1)=='/') then ! is absolute path --> try mss
#ifdef MSS
!
! -noreplace option to msrcp means use local if it already exists:
          istat=my_msrcp('-noreplace','mss:'//trim(remote),trim(local)) 
          if (istat==0) then
            write(6,"('Getfile: copied ',a,' from mss to local ',a)")
     |        trim(remote),trim(local)
          else
            write(6,"(/,'>>> Getfile: msrcp of ',a,' to local ',a,
     |        ' failed: istat=',i6)") trim(remote),trim(local),istat
            call shutdown('Getfile: msrcp failed.')
          endif
#else
          write(6,"(/,'>>> Getfile: could not find file ',a)")
     |      trim(remote)
          write(6,"('    (NCAR mass store not available ',
     |      '(macro MSS is not defined))')")
          call shutdown('Getfile: file not found.')
#endif
        else ! remote does not exist, and is not absolute path
          write(6,"(/,'>>> Getfile: could not find file ',a)")
     |      trim(remote)
          call shutdown('Getfile: file not found')
        endif
      endif ! existence of remote
 100  continue
      write(6,"(72('-'),/)")
      end subroutine getfile
