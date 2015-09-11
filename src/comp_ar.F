#include <defs.h>
      module ar_module
!
! Advance argon by one time step.
!
! Boundary conditions, production and loss for argon are defined
! by comp_ar, and referenced by minor_ar. Comp_ar is called
! from a latitude loop in dynamics. After comp_ar, dynamics calls
! minor_ar, which passes this module data to sub minor. Sub
! minor contains 3d mpi calls and its own latitude loops. 
!
      use params_module,only: nlevp1
      use addfld_module,only: addfld
      implicit none
!
! Boundary conditions and production and loss terms are allocated
! on subdomains by sub alloc_ar (called from allocdata.F).
!
      real,allocatable,dimension(:,:)   :: ar_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: ar_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  ar_prod,     ! production of argon (k,i,j)
     |  ar_loss      ! loss of argon       (k,i,j)
!
#if (NLEV==28)
!
! Initial condition for argon at 5.0-deg resolution. This is used 
! to expediate equilibration of argon when it is not on the
! source history (i.e., is zero on startup). These data were 
! obtained from a run of the glbmean model (Aug, 2014).
!
      real,parameter :: ar_glbm(nlevp1) = (/
     | 1.164e-02, 1.074e-02, 9.443e-03, 7.905e-03, 6.342e-03,
     | 4.927e-03, 3.737e-03, 2.776e-03, 2.017e-03, 1.430e-03,
     | 9.855e-04, 6.575e-04, 4.229e-04, 2.613e-04, 1.548e-04,
     | 8.797e-05, 4.806e-05, 2.535e-05, 1.298e-05, 6.488e-06,
     | 3.184e-06, 1.541e-06, 7.381e-07, 3.511e-07, 1.662e-07,
     | 7.840e-08, 3.689e-08, 1.733e-08, 7.436e-09 /)

#elif (NLEV==56)
!
! Initial condition for argon at 2.5-deg resolution. This is used 
! to expediate equilibration of argon when it is not on the
! source history (i.e., is zero on startup). These data were 
! obtained from a run of the glbmean model (Aug, 2014).
!
      real,parameter :: ar_glbm(nlevp1) = (/
     | 1.164e-02, 1.124e-02, 1.074e-02, 1.013e-02, 9.443e-03,
     | 8.691e-03, 7.905e-03, 7.114e-03, 6.342e-03, 5.609e-03,
     | 4.927e-03, 4.302e-03, 3.737e-03, 3.229e-03, 2.776e-03,
     | 2.373e-03, 2.017e-03, 1.704e-03, 1.430e-03, 1.192e-03,
     | 9.855e-04, 8.084e-04, 6.575e-04, 5.298e-04, 4.229e-04,
     | 3.341e-04, 2.613e-04, 2.022e-04, 1.548e-04, 1.173e-04,
     | 8.797e-05, 6.533e-05, 4.806e-05, 3.505e-05, 2.535e-05,
     | 1.820e-05, 1.298e-05, 9.201e-06, 6.488e-06, 4.554e-06,
     | 3.184e-06, 2.218e-06, 1.541e-06, 1.068e-06, 7.381e-07,
     | 5.095e-07, 3.511e-07, 2.417e-07, 1.662e-07, 1.142e-07,
     | 7.840e-08, 5.380e-08, 3.689e-08, 2.529e-08, 1.733e-08,
     | 1.188e-08, 7.436e-09 /)
#endif
!
      contains
!-----------------------------------------------------------------------
      subroutine comp_ar(ar,lev0,lev1,lon0,lon1,lat)
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  ar           ! argon from previous step
!
! Local:
      integer :: i
!
! ar_lbc, ar_ubc, ar_glbm, ar_prod, ar_loss are module data above.
!
      do i=lon0,lon1 
        ar_lbc(i,1,lat) = 0.
        ar_lbc(i,2,lat) = 1.
        ar_lbc(i,3,lat) = -sqrt(ar_glbm(1)*ar_glbm(2))
        ar_ubc(i,lat)   = 0.
      enddo
!     write(6,"('comp_ar: lat=',i3,' ar_lbc(3)=',/,(6e12.4))") 
!    |  lat,ar_lbc(:,3,lat)
!
! Production and loss of argon are zero:
      ar_prod = 0.
      ar_loss = 0.

      end subroutine comp_ar
!-----------------------------------------------------------------------
      subroutine minor_ar(tn,o2,o1,n2,he,ar,ar_nm,ar_out,arnm_out,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
      use cons_module,only: rmass_ar
      implicit none
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn,     ! neutral temperature (deg K)
     |  o2,     ! molecular oxygen (mmr)
     |  o1,     ! oxygen family (mmr)
     |  n2,     ! molecular nitrogen (mmr)
     |  he,     ! helium (mmr)
     |  ar,     ! argon (mmr)
     |  ar_nm   ! argon at time n-1
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  ar_out,    ! ar output
     |  arnm_out   ! ar output at time n-1
!
! Local:
      real,parameter :: phi_ar(3) = (/1.042,1.509,1.176/) ! from tgcm24
      real,parameter :: alfa_ar   = 0.17    ! thermal diffusion coefficient (from tgcm24)
      real,parameter :: xyar      = 1.e-10  ! from tgcm24
!
      call minor(tn,o2,o1,n2,he,ar,ar_nm,ar_out,arnm_out,ar_loss,
     |  ar_prod,ar_lbc,ar_ubc,rmass_ar,phi_ar,alfa_ar,lev0,lev1,
     |  lon0,lon1,lat0,lat1,0,'AR')
!
      end subroutine minor_ar
!-----------------------------------------------------------------------
      subroutine alloc_ar(lon0,lon1,lat0,lat1)
!
! Allocate subdomains (without ghost cells) to module data for boundary
! conditions and production and loss terms. This is called once per run
! from sub allocdata (allocdata.F).
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
! Allocate subdomains to boundary conditions:
      allocate(ar_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_ar: error allocating',
     |  ' ar_ubc: stat=',i3)") istat
      allocate(ar_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_ar: error allocating',
     |  ' ar_lbc: stat=',i3)") istat
      ar_ubc = 0. ; ar_lbc = 0.
!
! Allocate subdomains to production and loss:
      allocate(ar_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_ar: error allocating',
     |  ' ar_prod: stat=',i3)") istat
      allocate(ar_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_ar: error allocating',
     |  ' ar_loss: stat=',i3)") istat
      ar_prod = 0. ; ar_loss = 0.
!
      write(6,"('alloc_ar: allocated module data')")
!
      end subroutine alloc_ar
!-----------------------------------------------------------------------
      end module ar_module
