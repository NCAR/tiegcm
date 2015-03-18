!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Definitions of grid parameters for pre-processor.
! See parameters.h.
!
!------------------------------------
! 5.0 degree horizontal:
!
! Latitude dimension:
#define NLAT  (36)
#define GLAT1 (-87.5)
#define DLAT  (5.)
!
! Longitude dimension:
#define NLON  (72)
#define GLON1 (-180.)
#define DLON  (5.)
!
!------------------------------------
! 2.5 degree horizontal:
!
! Latitude dimension:
!#define NLAT  (72)
!#define GLAT1 (-88.75)
!#define DLAT  (2.5)
!
! Longitude dimension:
!#define NLON  (144)
!#define GLON1 (-180.)
!#define DLON  (2.5)
!
!------------------------------------
! Vertical column dimension:
! There are 2 supported vertical resolutions:
!
! ZBOT  ZTOP  DZ   NLEV
! -7     5   0.5    28  "normal resolution" 2 grid points per scale height
! -7     5   0.25   56  "double resolution" 4 grid points per scale height
!
! Define interface and midpoint levels
! Vertical column -7 to +7 by 0.50 ("normal")
#define ZIBOT (-7.0)
#define ZITOP (7.0)
#define ZMBOT (-6.75)
#define ZMTOP (7.25)
#define NLEV (28)
!
! Vertical column -7 to +7 by 0.25 ("double")
!#define ZIBOT (-7.0)
!#define ZITOP (7.0)
!#define ZMBOT (-6.875)
!#define ZMTOP (7.125)
!#define NLEV (56)
!
