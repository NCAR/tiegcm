!
! defs.h contains c-preprocessor macros.
!
#include <defs.h>
!
      module params_module
!
! Geographic grid parameters:
      integer,parameter :: 
     |  nlat = NLAT,                 ! number of latitudes
     |  nlon = NLON,                 ! number of longitudes
     |  nlev = NLEV,                 ! number of midpoint levels 
     |  nilev = NLEV,                ! number of interface levels
     |  nlonp4 = nlon+4,             ! includes periodic points in lon
     |  nlevp1 = nlev+1,
     |  nilevp1 = nilev+1,
     |  nlonp1=nlon+1, nlatp1=nlat+1,
     |  nlonp2=nlon+2, nlatp2=nlat+2
      real,parameter ::
     |  glat1 = GLAT1,               ! first latitude
     |  dlat  = DLAT,                ! delta latitude
     |  glon1 = GLON1,               ! first longitude
     |  dlon  = DLON,                ! delta longitude
     |  dlev  = (ZITOP-ZIBOT)/NLEV,  ! delta level (same for midpoints and interfaces)
     |  dz    = dlev,                ! delta level (alias for dlev)
     |  zmtop = ZMTOP,               ! top midpoint level
     |  zmbot = ZMBOT,               ! bottom midpoint level
     |  zitop = ZITOP,               ! top interface level
     |  zibot = ZIBOT                ! bottom interface level
      real,parameter ::
     |  spval = 1.e36
      integer,parameter ::
     |  ispval = 999
!
! Magnetic grid:
      integer,parameter ::
     |  nmlat = 97,         ! number of magnetic latitudes
     |  nmlon = 80,         ! number of magnetic longitudes
     |  nmlonp1=nmlon+1,
     |  nmlonp2=nmlon+2,
     |  nmlatp1=nmlat+1,
     |  nmlath=(nmlat+1)/2, ! index to magnetic equator
     |  nmlev = nlev+3,   ! number of vertical levels of mag fields (-2,nlev)
     |  nimlev  = nmlev,    ! number of interface levels of mag fields (-2,nlev)
     |  nmlevp1 = nmlev+1,  ! 
     |  nimlevp1= nimlev+1  ! 
!
! kbotdyn is zp k-index at which upward dynamo integrals begin (dynamo.F):
! For zp -8 at dzp=0.5,  kbotdyn = 5 
! For zp -8 at dzp=0.25, kbotdyn = 9 
!
! -7 to +7 by 0.50 ("normal" vertical resolution)
#if (NLEV==28)
!     integer,parameter :: kbotdyn = 1
!
! -7 to +7 by 0.25 ("double" vertical resolution)
#elif (NLEV==56) 
!     integer,parameter :: kbotdyn = 1
#endif
!
! For modules hist_mod and input_mod:
      character(len=16),parameter :: tgcm_version = 'tiegcm1.81      '
!
! tgcm_name is either "tiegcm" or "time-gcm"
      character(len=8),parameter :: tgcm_name    = 'tiegcm'
!
      integer,parameter ::
     |  mxday=366,         ! max allowed calendar day
     |  mxhvols = 100,     ! max number of output history files
     |  mxseries = 10,     ! max number of time series for primary histories
     |  mxseries_sech = 50,! max number of secondary history time series
     |  mxfsech = 100,     ! max number of fields on secondary histories
     |  mxind_time = 100   ! max number of time-dependent solar index points
!
      end module params_module
