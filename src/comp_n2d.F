!
      module n2d_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      implicit none
!
! Boundary conditions, production and loss for N2D are defined
! by comp_n2d, and referenced by minor_n2d. Comp_n2d is called
! from a latitude loop in dynamics. After comp_n2d, dynamics calls
! minor_n2d, which passes this module data to sub minor. Sub 
! minor contains 3d mpi calls and its own latitude loops.
!
! Boundary conditions and production and loss terms are allocated
! subdomains by sub alloc_n2d (called from allocdata.F).
!
      real,allocatable,dimension(:,:)   :: n2d_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: n2d_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  n2d_prod,     ! production of n2d (k,i,j)
     |  n2d_loss      ! loss of n2d       (k,i,j)
!
! These numbers are copied from N4S
      real :: phi_n2d(4) = (/0.651, 0.731, 0.627, 0.741/)
      real,parameter :: alfa_n2d = 0. ! thermal diffusion coefficient
!
      contains
!-----------------------------------------------------------------------
      subroutine alloc_n2d(lon0,lon1,lat0,lat1)
!
! Allocate subdomains (without ghost cells) to module data for boundary 
! conditions and production and loss terms. This is called once per run 
! from sub allocdata (allocdata.F).
!
      use params_module,only: nlevp1
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
! Allocate subdomains to boundary conditions:
      allocate(n2d_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n2d: error allocating',
     |  ' n2d_ubc: stat=',i3)") istat
      allocate(n2d_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n2d: error allocating',
     |  ' n2d_lbc: stat=',i3)") istat
!
! Allocate subdomains to production and loss:
      allocate(n2d_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n2d: error allocating',
     |  ' n2d_prod: stat=',i3)") istat
      allocate(n2d_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n2d: error allocating',
     |  ' n2d_prod: stat=',i3)") istat

!     write(6,"('alloc_n2d: allocated module data')")
!
      end subroutine alloc_n2d
!-----------------------------------------------------------------------
      subroutine comp_n2d(o2,o1,no,ne,op,n2p,nop,xnmbar,xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
!
! Advance n2d by one time step. This is called from dynamics at
! each subdomain latitude.
!
      use cons_module,only: brn2d,rmassinv_o2,rmassinv_o1,rmassinv_no,
     |  rmass_n2d
      use chemrates_module,only: rk3,ra1,ra3,beta2,beta4,beta5,beta6,
     |  beta7,rk10
      use qrj_module,only: qtef
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! Input arrays:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  o2,     ! molecular oxygen (mmr)
     |  o1,     ! atomic oxygen (mmr)
     |  no,     ! nitric oxide (mmr)
     |  ne,     ! electron density (cm3)
     |  op,     ! O+ ion
     |  n2p,    ! N2+ ion
     |  nop,    ! NO+ ion
     |  xnmbar, ! p0*e(-z)*mbar/kT
     |  xnmbari ! p0*e(-z)*barm/kT
!
! Local:
      integer :: k,i
      integer ::i0,i1,nk,nkm1
!
! For addfld:
      i0 = lon0
      i1 = lon1
      nk = lev1-lev0+1
      nkm1 = nk-1
!
      do i=lon0,lon1
!
! Value at bottom given by photochemical equilibrium.
        n2d_lbc(i,1,lat) = 0.
        n2d_lbc(i,2,lat) = 1.
        n2d_lbc(i,3,lat) = -rmass_n2d/xnmbari(lev0,i)*
     |    qtef(lev0,i,lat)*brn2d/
     |    (xnmbari(lev0,i)*
     |    (beta2(lev0,i,lat)*o2(lev0,i)*rmassinv_o2+
     |     beta4*o1(lev0,i)*rmassinv_o1+
     |     beta6*no(lev0,i)*rmassinv_no)+
     |     beta7+beta5(lev0,i,lat)*ne(lev0,i))
!
! Zero diffusive flux at top:
        n2d_ubc(i,lat) = 0.
!
! Sources:
        do k=lev0,lev1-1
!
! Total production of N2D:
          n2d_prod(k,i,lat) = .5*(qtef(k,i,lat)+qtef(k+1,i,lat))*brn2d+
     |      rk3(k,i,lat)*n2p(k,i)*xnmbar(k,i)*o1(k,i)*rmassinv_o1+
     |     (ra1(k,i,lat)*nop(k,i)*0.85+ra3(k,i,lat)*n2p(k,i)*0.9)*
     |      sqrt(ne(k,i)*ne(k+1,i))
!
! Total loss of N2D:
          n2d_loss(k,i,lat) = -(xnmbar(k,i)*
     |      (beta2(k,i,lat)*o2(k,i)*rmassinv_o2+
     |       beta4*o1(k,i)*rmassinv_o1+
     |       beta6*no(k,i)*rmassinv_no)+
     |       beta7+beta5(k,i,lat)*sqrt(ne(k,i)*ne(k+1,i))+
     |       rk10*op(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('N2D_PROD' ,' ',' ',n2d_prod(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('N2D_LOSS' ,' ',' ',n2d_loss(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)

      end subroutine comp_n2d
!-----------------------------------------------------------------------
      subroutine minor_n2d(tn,o2,o1,he,mbar,barm,xnmbar,w,n2d,n2d_nm1,
     |  n2d_out,n2d_nm1_out,lev0,lev1,lon0,lon1,lat0,lat1)
      use cons_module,only: rmass_n2d
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  he,      ! helium (mmr)
     |  mbar,    ! mean molecular weight
     |  barm,    ! mean molecular weight
     |  xnmbar,  ! p0*e(-z)*mbar/kT
     |  w,       ! vertical wind
     |  n2d,     ! N2D (mmr)
     |  n2d_nm1  ! N2D at time n-1 
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  n2d_out,    ! N2D output 
     |  n2d_nm1_out ! N2D output at time n-1
!
!     write(6,"('enter minor_n2d')")
!
! Minor returns n2d_out and n2d_nm1_out. Module data n2d_prod,
! n2d_loss, etc, were defined by comp_n2d.  
!
      call minor(tn,o2,o1,he,mbar,barm,xnmbar,w,n2d,n2d_nm1,
     |  n2d_out,n2d_nm1_out,n2d_loss,n2d_prod,n2d_lbc,n2d_ubc,
     |  rmass_n2d,phi_n2d,alfa_n2d,lev0,lev1,lon0,lon1,lat0,lat1,
     |  0,'N2D')

      end subroutine minor_n2d
!-----------------------------------------------------------------------
      end module n2d_module
