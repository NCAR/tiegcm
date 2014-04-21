      module ModIE_Interface
      
        use ModKind

        real, allocatable, dimension(:,:,:) :: IEr3_HaveLats, 
     |   IEr3_HaveMLTs
        real, allocatable, dimension(:,:,:) :: IEr3_HavePotential
        real, allocatable, dimension(:,:,:) :: IEr3_HaveEFlux
        real, allocatable, dimension(:,:,:) :: IEr3_HaveAveE

        real (Real8_) :: IEd_CurrentTime
!        real :: IEd_CurrentTime
        integer             :: IEi_HavenLats
        integer             :: IEi_HavenMLTs
        integer             :: IEi_HavenBLKs
        integer             :: IEi_HavenTimes

        real (Real8_)               :: GMd_NeedTime = -1.0e32
!        real                :: GMd_NeedTime = -1.0e32
        real, allocatable, dimension(:,:) :: GMr2_NeedLats, 
     |   GMr2_NeedMLTs
        real, allocatable, dimension(:,:) :: GMr2_NeedPotential
        real, allocatable, dimension(:,:) :: GMr2_NeedEFlux
        real, allocatable, dimension(:,:) :: GMr2_NeedAveE
        integer                           :: GMi_NeednLats
        integer                           :: GMi_NeednMLTs
        integer                           :: GMi_NeednTimes
        integer, allocatable, dimension(:,:,:) :: 
     |   GMi3_InterpolationIndices
        real, allocatable, dimension(:,:,:)    :: 
     |   GMr3_InterpolationRatios

        real (Real8_)               :: IMd_NeedTime = -1.0e32
!        real               :: IMd_NeedTime = -1.0e32
        real, allocatable, dimension(:,:) :: IMr2_NeedLats, 
     |   IMr2_NeedMLTs
        real, allocatable, dimension(:,:) :: IMr2_NeedPotential
        real, allocatable, dimension(:,:) :: IMr2_NeedEFlux
        real, allocatable, dimension(:,:) :: IMr2_NeedAveE
        integer                           :: IMi_NeednLats
        integer                           :: IMi_NeednMLTs
        integer                           :: IMi_NeednTimes
        integer, allocatable, dimension(:,:,:) :: 
     |   IMi3_InterpolationIndices
        real, allocatable, dimension(:,:,:)    :: 
     |    IMr3_InterpolationRatios

        real (Real8_)               :: IOd_NeedTime = -1.0e32
!        real                :: IOd_NeedTime = -1.0e32
        real, allocatable, dimension(:,:) :: IOr2_NeedLats, 
     |   IOr2_NeedMLTs
        real, allocatable, dimension(:,:) :: IOr2_NeedPotential
        real, allocatable, dimension(:,:) :: IOr2_NeedEFlux
        real, allocatable, dimension(:,:) :: IOr2_NeedAveE
        integer                           :: IOi_NeednLats
        integer                           :: IOi_NeednMLTs
        integer                           :: IOi_NeednTimes
        integer, allocatable, dimension(:,:,:) :: 
     |   IOi3_InterpolationIndices
        real, allocatable, dimension(:,:,:)    :: 
     |   IOr3_InterpolationRatios
        real :: IOr_NeedIMFBz   = -1.0e32
        real :: IOr_NeedIMFBy   = -1.0e32 
        real :: IOr_NeedSWV     = -1.0e32 
        real :: IOr_NeedHPI     = -1.0e32 
        real :: IOr_NeedHPINorm = -1.0e32 
        real :: IOr_NeedKp      = -1.0e32 
        logical :: IOl_IsNorth  = .true.

        integer                           :: iDebugLevel = 0
        integer                           :: iProc = 0

        integer, parameter                :: IE_Closest_     = 1
        integer, parameter                :: IE_After_       = 2
        integer, parameter                :: IE_Interpolate_ = 3

        character (len=100) :: IE_NameOfEFieldModel
        character (len=100) :: IE_NameOfAuroralModel
        character (len=100) :: IE_NameOfSolarModel
        character (len=100) :: IE_NameOfModelDir

        logical :: UseGridBasedIE

      end module ModIE_Interface
