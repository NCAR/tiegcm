!
      subroutine inp_read(mytid)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Read input namelist. The ONLY reason this is not in input_mod.f
! is so mpi_mod can be recompiled easily (i.e., avoid MPI dependency
! in input_mod.F.
!
#ifdef SUN
      use input_module
#else
      use input_module,only: tgcm_input
#endif
      implicit none
#ifdef MPI
#include <mpif.h>
#endif
!
! Args:
      integer,intent(in) :: mytid
!
! Local:
      integer :: luin=7,ios,ier
      character(len=80) :: rec80
!
! Namelist read:
! If mpi job, only master can write new comment-stripped input file.
!   (i.e., slaves have to wait for master to call rmcomments)
!
#ifdef MPI
      if (mytid==0) then      
        call rmcomments(5,luin,';',1)
        call mpi_barrier(MPI_COMM_WORLD,ier)
      else
        call mpi_barrier(MPI_COMM_WORLD,ier)
      endif
#else
      call rmcomments(5,luin,';',1)
#endif
      write(6,"(/,'Reading input data...',/)")
      read(luin, nml=tgcm_input,err=900)
!     read(luin, nml=tgcm_input,iostat=ios,err=900)
!     if (ios /= 0) then
!       write(6,"(/,'>>> Error from namelist read: ios=',i5,
!    |    ' luin=',i3)") ios,luin
!       call shutdown('inp_read')
!     endif
      close(luin)
      write(6,"('Completed successful read of namelist inputs.',/)")
      return
 900  continue
!
! Error in namelist read:
! (on the Crays, ios=1324 means unrecognized keyword in namelist input)
!
      write(6,"(/,72('>'),/,'ERROR in namelist read of user inputs: ',
     |  ' lu=',i2,' ios=',i5/)") luin,ios
      backspace(luin)
      read(luin,"(a)") rec80
      write(6,"('This might be the result of an unrecognized ',
     |  'or misspelled keyword in the input file.')")
      write(6,"('Please check your input file in the vicinity ',
     |  'of the following line:')")
      write(6,"(/,a,/)") rec80
      write(6,"(72('<'),/)")
      close(luin)
      call shutdown('inp_read')
      end subroutine inp_read

