      Module ModAMIE_Interface

      character (len=100) :: AMIE_FileName
      integer :: AMIE_nLats, AMIE_nMlts, AMIE_nTimes

  ! For a single file
      real*4, allocatable,dimension(:)     :: AMIE_Lats, AMIE_MLTs
      real*4, allocatable,dimension(:,:,:,:) :: AMIE_Potential,
     |   AMIE_EFlux,AMIE_AveE
      real*4, allocatable,dimension(:,:,:,:) :: AMIE_Value
      real*8, allocatable,dimension(:,:)     :: AMIE_Time

!  ! For a North and South file
!  real, allocatable, dimension(:,:)     :: AMIE_Lats,AMIE_MLTs
!  real, allocatable, dimension(:,:,:,:) :: AMIE_Potential,AMIE_EFlux,AMIE_AveE
!  real, allocatable, dimension(:,:,:,:) :: AMIE_Value
!  real*8, allocatable, dimension(:,:)   :: AMIE_Time

      integer, parameter :: AMIE_Closest_     = 1
      integer, parameter :: AMIE_After_       = 2
      integer, parameter :: AMIE_Interpolate_ = 3

      integer :: AMIE_iDebugLevel = 0

      integer :: AMIE_South_ = 1
      integer :: AMIE_North_ = 2

      end Module ModAMIE_Interface
