!BOP
!MODULE: ModMpiOrig and ModMpi - the MPI variables and functions
!DESCRIPTION:
! In Fortran 90 it is customary to use a module instead of including files.
! The ModMpiOrig and ModMpi modules provide an interface to the selected 
! MPI header file mpif90.h. The MPI header file is copied from an
! operating system specific version during installation time.
!
! The original MPI header files have been modified such 
! that the MPI\_REAL parameter is set to the value of the 
! MPI\_DOUBLE\_PRECISION parameter if iRealPrec is 1. The iRealPrec
! parameter is set to 0 in ModMpiOrig, while it is set to 0 or 1
! according to the actual precision of the real numbers in ModMpi.
! If ModMpi is used, it should be compiled with the same precision as the F90
! code using it.
!
!REVISION HISTORY:
! 07/02/2003 G. Toth <gtoth@umich.edu> - initial version of ModMpi
! 07/20/2003 G. Toth - change the MPI_REAL definition
! 07/30/2004 G. Toth - updated the description for the modified mpif90.h files.
!INTERFACE:
      module ModMpiOrig
  !EOP
  !BOC
        implicit none
        integer, parameter :: iRealPrec = 0
        include 'mpif90.h'
  !EOC
      end module ModMpiOrig

!BOP
!INTERFACE:
      module ModMpi
  !EOP
  !BOC
        implicit none
  ! iRealPrec = 1 if the code is compiled with 8 byte reals and 0 otherwise
        integer, parameter :: 
     |    iRealPrec = (1.00000000011 - 1.0)*10000000000.0
        include 'mpif90.h'
  !EOC
      end module ModMpi
