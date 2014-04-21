!^CFG COPYRIGHT UM
!BOP
!MODULE: ModKind - define various precisions in a machine independent way

!DESCRIPTION:
! The Fortran 77 style real*4 and real*8 declarations are obsolete,
! and compilers often issue warnings. The real and double precision
! types are machine and compiler flag dependent.
! The Fortran 90 way is to define the {\bf kind} parameter.
! Typical usage:
!\begin{verbatim}
! real(Real8_) :: CpuTime  ! variable declaration
! CpuTime = 0.0_Real8_     ! 8 byte real constant
!\end{verbatim}

!INTERFACE:
      module ModKind

!PUBLIC DATA MEMBERS:
       integer, parameter :: Real4_=selected_real_kind(6,30)
       integer, parameter :: Real8_=selected_real_kind(12,100)

!EOP
      end module ModKind
