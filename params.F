!
! defs.h contains c-preprocessor macros.
!
#include "defs.h"
!
      module params_module
!
! Geographic grid parameters:
      integer,parameter :: 
     |  nlat = NLAT,                 ! number of latitudes
     |  nlon = NLON,                 ! number of longitudes
     |  nlev = NLEV,                 ! number of midpoint levels
     |  nlevp1 = nlev+1,             ! number of interface levels
     |  nlonp4 = nlon+4,             ! includes periodic points in lon
     |  nlonp1=nlon+1, nlatp1=nlat+1,
     |  nlonp2=nlon+2, nlatp2=nlat+2
      real,parameter ::
     |  glat1 = GLAT1,               ! first latitude
     |  dlat  = DLAT,                ! delta latitude
     |  glon1 = GLON1,               ! first longitude
     |  dlon  = DLON,                ! delta longitude
     |  plev1 = ZBOT,                ! bottom level
     |  dlev  = (ZTOP-ZBOT)/NLEV,    ! delta level
     |  dz    = dlev,                ! delta level (alias for dlev)
     |  zsb   = ZBOT,                ! bottom level (alias for plev1)
     |  zst   = ZTOP                 ! top level
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
     |  nmlev = nlevp1+3    ! number of vertical levels of mag fields (-2,nlev)
!
! Magnetospheric grid:
      integer,parameter ::
     |  nmagphrlat = 31,
     |  nmagphrlon = 40
      real,parameter ::
     |  magphrlat1 = 71.97, magphrlat2 = 10.14,
     |  magphrlon1 = -180.
!
! For modules hist_mod and input_mod:
      character(len=8),parameter :: tgcm_version = 'tiegcm1 '
!
! tgcm_name is either "tiegcm" or "time-gcm"
      character(len=8),parameter :: tgcm_name    = 'tiegcm'
!
      integer,parameter ::
     |  mxday=366,         ! max allowed calendar day
     |  mxhvols = 100,     ! max number of output history files
     |  mxseries = 10,     ! max number of time series for primary histories
     |  mxseries_sech = 50,! max number of secondary history time series
     |  mxfsech = 50       ! max number of fields on secondary histories
!
      end module params_module
