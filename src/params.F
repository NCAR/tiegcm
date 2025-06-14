!
! defs.h contains c-preprocessor macros.
!
#include <defs.h>
!
      module params_module
      use, intrinsic :: iso_fortran_env
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Geographic grid parameters:
      integer,parameter :: rp = REAL64
      real(rp),parameter ::
     |  dlat  = DLAT,                ! delta latitude
     |  glat1 = -90._rp+dlat/2._rp,          ! first latitude
     |  dlon  = DLON,                ! delta longitude
     |  glon1 = GLON1,               ! first longitude
     |  dlev  = DLEV,                ! delta level (same for midpoints and interfaces)
     |  dz    = dlev,                ! delta level (alias for dlev)
     |  zitop = ZITOP,               ! top interface level
     |  zibot = ZIBOT,               ! bottom interface level
     |  zmtop = zitop+dz/2._rp,          ! top midpoint level
     |  zmbot = zibot+dz/2._rp           ! bottom midpoint level
      integer,parameter :: 
     |  nlat = 180/dlat,             ! number of latitudes
     |  nlon = 360/dlon,             ! number of longitudes
     |  nlev = (zitop-zibot)/dlev,   ! number of midpoint levels 
     |  nilev = nlev,                ! number of interface levels
     |  nlonp4 = nlon+4,             ! includes periodic points in lon
     |  nlevp1 = nlev+1,
     |  nilevp1 = nilev+1,
     |  nlonp1=nlon+1, nlatp1=nlat+1,
     |  nlonp2=nlon+2, nlatp2=nlat+2
      real(rp) :: glon(nlon),glat(nlat),zpmid(nlevp1),zpint(nlevp1)
      real(rp),parameter ::
     |  spval = 1.e36_rp
      integer,parameter ::
     |  ispval = 999
!
! Magnetic grid:
      integer,parameter ::
     |  nres_grid = NRES_GRID,
     |  nmlat = 3*2**nres_grid+1,    ! number of magnetic latitudes
     |  nmlon = 5*2**(nres_grid-1),  ! number of magnetic longitudes
     |  nmlonp1=nmlon+1,
     |  nmlonp2=nmlon+2,
     |  nmlatp1=nmlat+1,
     |  nmlath=(nmlat+1)/2  ! index to magnetic equator

      real(rp),parameter :: zpbot_dyn  = -8.25_rp ! bottom midpoint boundary of dynamo
      real(rp),parameter :: zpibot_dyn = -8.5_rp  ! bottom interface boundary of dynamo
      real(rp),parameter :: dmlev = dlev
!
! nmlev_diff = number of levels from zibot down to zpibot_dyn
! nmlevp1 = total number of mag levels
!
      integer,parameter ::
     |  nmlev_diff = int((zibot-zpibot_dyn)/dlev), 
     |  nmlevp1 = nlevp1+nmlev_diff, 
     |  nimlevp1= nmlevp1
!
      integer,parameter :: mlev0 = 1 - nmlev_diff
      integer,parameter :: mlev1 = nlevp1

      real(rp) :: gmlon(nmlonp1), ! magnetic longitude (deg)
     |        gmlat(nmlat),   ! magnetic latitude (deg)
     |        zpmag(nmlevp1),   ! magnetic midpoint levels
     |        zpimag(nimlevp1)  ! magnetic interface levels
!
      real(rp),dimension(nlonp4) :: glon0 ! include ghost points (1->nlonp4)
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
