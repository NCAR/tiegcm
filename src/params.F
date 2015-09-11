!
! defs.h contains c-preprocessor macros.
!
#include <defs.h>
!
      module params_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
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
      real :: glon(nlon),glat(nlat),zpmid(nlevp1),zpint(nlevp1)
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
     |  nmlath=(nmlat+1)/2  ! index to magnetic equator

      real,parameter :: zpbot_dyn  = -8.25 ! bottom midpoint boundary of dynamo
      real,parameter :: zpibot_dyn = -8.5  ! bottom interface boundary of dynamo
      real,parameter :: dmlev = dlev
!
! nmlev_diff = number of levels from zibot down to zpibot_dyn (3 at 5-deg, 6 at 2.5-deg)
! nmlevp1 = total number of mag levels (-2 -> 29 = 32 levels at 5-deg)
!
      integer,parameter ::
     |  nmlev_diff = int((zibot-zpibot_dyn)/dlev), 
     |  nmlevp1 = nlevp1+nmlev_diff, 
     |  nimlevp1= nmlevp1
!
! At 5-deg resolution:
!   mlev0 = -2, mlev1 = nlevp1 = 29
! At 2.5-deg resolution:
!   mlev0 = -5, mlev1 = nlevp1 = 57
!
      integer,parameter :: mlev0 = 1 - nmlev_diff ! (-2 at 5-deg, -5 at 2.5-deg)
      integer,parameter :: mlev1 = nlevp1

      real :: gmlon(nmlonp1), ! magnetic longitude (deg)
     |        gmlat(nmlat),   ! magnetic latitude (deg)
     |        zpmag(nmlevp1),   ! magnetic midpoint levels
     |        zpimag(nimlevp1)  ! magnetic interface levels
!
! For modules hist_mod and input_mod:
      character(len=16),parameter :: tgcm_version = 'tiegcm_trunk    '
!
! tgcm_name is either "tiegcm" or "time-gcm"
      character(len=8),parameter :: tgcm_name    = 'tiegcm'
!
      integer,parameter ::
     |  mxhvols = 500,     ! max number of output history files
     |  mxseries = 10,     ! max number of time series for primary histories
     |  mxseries_sech = 50,! max number of secondary history time series
     |  mxfsech = 500,     ! max number of fields on secondary histories
     |  mxind_time = 500   ! max number of time-dependent solar index points
!
      end module params_module
