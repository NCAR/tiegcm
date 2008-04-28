!
      subroutine putms(mspath,dskfile,wrpass,mscomment,msreten)
      use input_module,only: cwd ! current working directory
      implicit none
!
! Args:
      character(len=*),intent(in) :: mspath,dskfile,wrpass,mscomment
      integer,intent(in) :: msreten
!
! Local:
      integer,parameter :: mxfilelen = 240 ! maximum number chars in file name
      logical exists
      character(len=520) :: opts,testfile,cmd,tempfile
      character(len=8),save :: curdate,curtime
      integer :: ier,istat
!
! External:
      integer,external :: isystem,my_msrcp
!
      write(6,"(' ',/72('-'),/'PUTMS:')")

!     write(6,"('  mspath  = ',a)") trim(mspath) 
!     write(6,"('  dskfile = ',a)") trim(dskfile) 
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
#else
        write(6,"('>>> putms: NCAR Mass Store is not available --',
     |    ' file ',a,' was NOT disposed to mss.')") trim(dskfile)
!
! End #ifdef MSS
#endif
      else  ! mspath is blank
        write(6,"('>>> putms: mspath is blank: file ',a,' NOT disposed',
     |    ' to mss')") trim(mspath)
      endif ! len_trim(mspath)
!
      write(6,"(72('-'),/)")
      end subroutine putms
