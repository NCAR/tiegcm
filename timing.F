!
      module timing_module
      implicit none
!
! f90 intrinsic system_clock returns current clock count, number of
! clock counts per second, and max possible value of the count.
!
      real :: elapsed_prim=0., elapsed_sech=0., elapsed_steps=0.,
     |        elapsed_mpi=0.
      contains
!-----------------------------------------------------------------------
      subroutine start_timing(icount_out,string)
!
! Args:
      integer,intent(out) :: icount_out
      character(len=*),intent(in) :: string
!
! Local:
      integer :: icount1,icount_rate,icount_max
!
      call system_clock(icount1,icount_rate,icount_max)
      icount_out = icount1

!     write(6,"('start_timing for ',a,': secs=',f10.3)") 
!    |  trim(string),real(icount1)/real(icount_rate)

      end subroutine start_timing
!-----------------------------------------------------------------------
      subroutine end_timing(icount_in,esecs)
!
! Args:
      integer,intent(in) :: icount_in
      real,intent(out) :: esecs
!
! Local:
      integer :: icount1,icount2,icount_rate,icount_max
      real :: secs1,secs2
!
      call system_clock(icount2,icount_rate,icount_max)
      icount1 = icount_in
      if (icount2 < icount1) icount2 = icount2+icount_max
      secs1 = real(icount1)/real(icount_rate)
      secs2 = real(icount2)/real(icount_rate)
      esecs = secs2-secs1

!     write(6,"('end_timing: secs1=',f10.3,' secs2=',f10.3,
!    |  ' elapsed=',f10.3)") secs1,secs2,esecs

      end subroutine end_timing
      end module timing_module
