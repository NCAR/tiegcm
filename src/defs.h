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
! Vertical column -7 to +7 by 0.50 ("normal")
#define ZBOT (-7.)
#define ZTOP (7.)
#define NLEV (28)
!
! Vertical column -7 to +7 by 0.25 ("double")
!#define ZBOT (-7.)
!#define ZTOP (7.)
!#define NLEV (56)
!
