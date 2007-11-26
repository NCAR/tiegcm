!
      subroutine getfile(remote,local)
      implicit none
!
! Get input (read-only) file "remote" to local disk file "local". 
! Remote must be non-blank. If local is non-blank and not the same as 
! remote, the file will be linked or copied to local. If local is blank 
! or the same as remote, no links or copies are made (but local is returned 
! as either remote itself, or in the case of an mss transfer, as a local 
! file name constructed from remote).
!
! If remote does not exist on the local disk system, and is an absolute 
! path (begins with "/"), this routine will look for remote on the NCAR MSS. 
! If remote exists on the MSS, it is copied to a local disk file name 
! constructed from the mss path (local is returned to the caller).
!
! If this routine is successful, local can be opened and read.
!
! Procedure:
!   0. If remote is blank, shutdown with error message.
!   1. If remote exists and local is blank or is same as remote, then simply 
!      return (no links or copies are made -- file can be read "in situ").
!   2. If remote exists and local is non-blank (and not equal to remote), 
!      then try the following in order, returning when successful:
!      a) Attempt a hard link from remote to local (ln)
!      b) Attempt a symbolic link from remote to local (ln -s)
!      c) Attempt to copy remote to local (cp)
!      If all of a-c above fail, call shutdown.
!   3. If remote does not exist and is not an absolute path, call shutdown.
!   4. If remote does not exist and is an absolute path (begins with "/"),
!      then go to NCAR mss (sub my_msrcp). If remote exists on mss, it is
!      copied to a local file name constructed from the mss path (sub mkdiskflnm).
!      If mss acquire fails, call shutdown.
!
! Examples of namelist read parameter SOURCE (startup history):
!
!   Example 1: SOURCE = '/ptmp/mydir/myfile.nc'
!              If this file exists on the current disk system, it will be 
!              opened and read "in situ" (no links or copy). If it does not 
!              exist, the model will look on the mss (and will shutdown if 
!              it does not exist there).
!
!   Example 2: SOURCE = '/USER/mydir/myfile.nc'
!              Same as example 1 above, but probably will be found on the
!              mss instead of the local disk system. If found on the mss,
!              it will be copied to the execution directory as 
!              USER.mydir.myfile.nc (see sub mkdiskflnm).
!
!   Example 3: SOURCE = '/fis/hao/data/TGCM.timegcm.p001-2002-080.nc'
!              This file exists on SCD systems ($TGCMDATA), so it will be 
!              opened and read "in situ", i.e., no links or copies will be 
!              made.
!
!   Example 4: SOURCE = 'TGCM.timegcm.p001-2002-080.nc'
!              If you use tgcm_config (e.g., from a tgcm job script), 
!              this file is linked from $TGCMDATA to your execution directory 
!              before model execution. The model will find it there.
!
! This routine is called from the following modules (with blank local), to 
! obtain the file specified by the corresponding namelist read input. 
! Example SOURCE files, and all other files in this list are available in 
! the $TGCMDATA directory.
!
!  rdsource.F:  SOURCE
!  gpi.F:       GPI_NCFILE
!  ncep.F:      NCEP_NCFILE
!  gswm.F:      GSWM_MDI_FILE, GSWM_MSDI_FILE, GSWM_NMDI_FILE, GSWM_NMSDI_FILE
!  solgar.F:    SOLGAR_BNDRY_FILE
!  solgar.F:    SOLGAR_IMPORT_FILE
!  wei01gcm.F:  (coefficients file -- no namelist)
!  zatmos.F:    ZATMOS_FILE
!
! Example namelist input that assumes files were linked from $TGCMDATA by
! tgcm_config script (see also tgcm command MkNamelist) (non-migrating GSWM 
! are not available as of 12/04).
!
!   SOURCE = 'TGCM.timegcm.p001-2002-080.nc'
!   SOLGAR_BNDRY_FILE  = 'solgar_bndry.nc'
!   SOLGAR_IMPORT_FILE = 'solgar_import.nc'
!   ZATMOS_FILE        = 'zatmos_bndry.nc'
!   GSWM_MDI_FILE      = 'gswm_diurn_5.0d_32km.nc'
!   GSWM_MSDI_FILE     = 'gswm_semi_5.0d_32km.nc'
!  ;GSWM_NMDI_FILE     = 'gswm_nonmig_diurn_5.0d_32km.nc'
!  ;GSWM_NMSDI_FILE    = 'gswm_nonmig_semi_5.0d_32km.nc'
!   GPI_NCFILE         = 'gpi_1979001-2004060.nc'
!
! This routine may call the following subroutines and functions:
!   subroutine mkdiskflnm (util.F)
!   subroutine shutdown (util.F)
!   functions isystem, ilink, isymlink, my_msrcp (util.F)
!
! WARNING: Because argument "local" is intent(inout) (i.e., may be reset 
!          upon return), DO NOT call this routine with a blank constant as 
!          the "local" argument, e.g. call getfile(path,' ').
!
! Args:
      character(len=*),intent(in) :: 
     |  remote    ! file to get (may be relative or absolute path)
      character(len=*),intent(inout) :: 
     |  local     ! local disk file
!
! Local:
      integer :: istat
      logical :: exists
      character(len=1024) :: cmd
!
! External (util.F):
      integer,external :: isystem,ilink,isymlink,my_msrcp
!
      write(6,"(/,72('-'))")
!
! Remote must be non-blank:
      if (len_trim(remote) <= 0) then
        write(6,"(/,'>>> Getfile: file argument ''remote'' is blank.')")
        call shutdown('Getfile: blank file argument')
      endif
!
      write(6,"('Getfile: remote=',a)") trim(remote)
      write(6,"('         local =',a)") trim(local)

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
        if (len_trim(local)==0) then
          write(6,"('  (this file was not copied or linked)')")
          local = remote
          goto 100
        elseif (trim(local)==trim(remote)) then
          write(6,"('  (this file was not copied or linked)')")
          goto 100
        endif
!
! Local file was given and is different than remote -> 
! copy or link to local (funcs ilink, isymlink, isystem are in util.F):
!
        istat = ilink(remote,local,1)     ! try hard link
        if (istat==0) goto 100            ! hard link succeeded
        istat = isymlink(remote,local,1)  ! try symbolic link
        if (istat==0) goto 100            ! sym link succeeded
        write(cmd,"('cp ',a,' ',a)") trim(remote),trim(local)
        istat = isystem(cmd) 
        if (istat==0) then                ! copy was successful
          write(6,"('Getfile: copied file ',a,' to ',a)") trim(remote),
     |      trim(local)
        else
          write(6,"('>>> Getfile: the following copy command failed:',
     |      a)") trim(cmd)  
          write(6,"('istat = ',i6)") istat
          call shutdown('Getfile: links and copy failed.')
        endif
!
! Remote as given does not exist. First check for existence of 
! mkdiskflnm(remote,local). If that does not exist, try obtaining 
! from mss (if MSS macro is set), otherwise give up.
!
      else  ! remote does not exist
        if (len_trim(local)==0) then
          call mkdiskflnm(remote,local) ! in util.F
        endif
        inquire(file=trim(local),exist=exists)
        if (exists) then
          write(6,"('Getfile: Found file ',a)") trim(local)
          write(6,"('  (this file was not copied or linked)')")
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
          call shutdown('Getfile: msrcp failed.')
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
