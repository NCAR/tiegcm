!
      module n4s_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      implicit none
!
! Boundary conditions, production and loss for N4S are defined
! by comp_n4s, and referenced by minor_n4s. Comp_n4s is called
! from a latitude loop in dynamics. After comp_n4s, dynamics calls
! minor_n4s, which passes this module data to sub minor. Sub 
! minor contains 3d mpi calls and its own latitude loops.
!
! Boundary conditions and production and loss terms are allocated
! subdomains by sub alloc_n4s (called from allocdata.F).
!
      real,allocatable,dimension(:,:)   :: n4s_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: n4s_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  n4s_prod,     ! production of n4s (k,i,j)
     |  n4s_loss      ! loss of n4s       (k,i,j)
!
! N4S~He coefficient (0.627) is guessed based on F~He coefficient from
! http://demonstrations.wolfram.com/BinaryDiffusionCoefficientsForGases
      real :: phi_n4s(4) = (/0.651, 0.731, 0.627, 0.741/)
      real,parameter :: alfa_n4s = 0. ! thermal diffusion coefficient
!
      contains
!-----------------------------------------------------------------------
      subroutine alloc_n4s(lon0,lon1,lat0,lat1)
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
      allocate(n4s_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n4s: error allocating',
     |  ' n4s_ubc: stat=',i3)") istat
      allocate(n4s_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n4s: error allocating',
     |  ' n4s_lbc: stat=',i3)") istat
!
! Allocate subdomains to production and loss:
      allocate(n4s_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n4s: error allocating',
     |  ' n4s_prod: stat=',i3)") istat
      allocate(n4s_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n4s: error allocating',
     |  ' n4s_prod: stat=',i3)") istat

!     write(6,"('alloc_n4s: allocated module data')")
!
      end subroutine alloc_n4s
!-----------------------------------------------------------------------
      subroutine comp_n4s(o2,o1,n2,xnmbar,xnmbari,no,n2d,ne,o2p,op,n2p,
     |  nplus,nop,lev0,lev1,lon0,lon1,lat)
!
! Advance n4s by one time step. This is called from dynamics at
! each subdomain latitude.
!
      use qrj_module,only: qtef
      use cons_module,only: brn2d,rmass_n4s,
     |  rmassinv_n2d,rmassinv_o2,rmassinv_o1,rmassinv_n2,rmassinv_no
      use chemrates_module,only: beta1,beta3,beta4,beta5,beta7,
     |  beta8,beta17,ra1,ra3,rk2,rk4,rk6,rk8
      use addfld_module,only: addfld
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  xnmbar,  ! p0*e(-z)*mbar/kT
     |  xnmbari, ! p0*e(-z)*barm/kT
     |  no,      ! nitric oxide (mmr)
     |  n2d,     ! N2D (mmr) (from sub comp_n2d)
     |  ne,      ! N2D (mmr)
     |  o2p,     ! O2+ ion
     |  op,      ! O+ ion
     |  n2p,     ! N2+ ion
     |  nplus,   ! N+ ion
     |  nop      ! NO+ ion
!
! Local:
      integer :: k,i
      real :: o2_cm3_lbc,o1_cm3_lbc,n2_cm3_lbc,n2d_cm3_lbc,no_cm3_lbc,
     |  o2_cm3,o1_cm3,n2_cm3,n2d_cm3,no_cm3
!
!     write(6,"('enter comp_n4s: lat=',i2)") lat
!
! Lower boundary:
! n4s_lbc(:,1)=A, n4s_lbc(:,2)=B, n4s_lbc(:,3)=C define lower boundary 
!   condition where A*DPSX/DZ + B*PSX + C = 0.
!
      do i=lon0,lon1
        o2_cm3_lbc = xnmbari(lev0,i)*o2(lev0,i)*rmassinv_o2
        o1_cm3_lbc = xnmbari(lev0,i)*o1(lev0,i)*rmassinv_o1
        n2_cm3_lbc = xnmbari(lev0,i)*n2(lev0,i)*rmassinv_n2
        n2d_cm3_lbc = xnmbari(lev0,i)*n2d(lev0,i)*rmassinv_n2d
        no_cm3_lbc = xnmbari(lev0,i)*no(lev0,i)*rmassinv_no
!
! Value at bottom given by photochemical equilibrium.
        n4s_lbc(i,1,lat) = 0.
        n4s_lbc(i,2,lat) = 1.
        n4s_lbc(i,3,lat) = -rmass_n4s/xnmbari(lev0,i)*
     |   (qtef(lev0,i,lat)*(1.-brn2d)+
     |    n2d_cm3_lbc*
     |    (beta4*o1_cm3_lbc+
     |     beta5(lev0,i,lat)*ne(lev0,i)+beta7)+
     |    beta8(lev0,i,lat)*no_cm3_lbc)/
!
     |   (beta1(lev0,i,lat)*o2_cm3_lbc+
     |    beta3(lev0,i,lat)*no_cm3_lbc+
     |    beta17(lev0,i,lat)*o1_cm3_lbc*n2_cm3_lbc)
!
! Zero diffusive flux at top:
        n4s_ubc(i,lat) = 0.
      enddo ! i=lon0,lon1
!
! Sources:
      do i=lon0,lon1
        do k=lev0,lev1-1
          o2_cm3 = xnmbar(k,i)*o2(k,i)*rmassinv_o2
          o1_cm3 = xnmbar(k,i)*o1(k,i)*rmassinv_o1
          n2_cm3 = xnmbar(k,i)*n2(k,i)*rmassinv_n2
          n2d_cm3 = xnmbar(k,i)*n2d(k,i)*rmassinv_n2d
          no_cm3 = xnmbar(k,i)*no(k,i)*rmassinv_no
!
          n4s_prod(k,i,lat) =
     |      .5*(qtef(k,i,lat)+qtef(k+1,i,lat))*(1.-brn2d)+
     |      n2d_cm3*
     |      (beta4*o1_cm3+
     |      beta5(k,i,lat)*.5*(ne(k,i)+ne(k+1,i))+beta7)+
     |      .5*(beta8(k,i,lat)+beta8(k+1,i,lat))*no_cm3+
     |      rk2(k,i,lat)*op(k,i)*n2_cm3+
     |      rk6*nplus(k,i)*o2_cm3+
     |      rk8*nplus(k,i)*o1_cm3+
     |      sqrt(ne(k,i)*ne(k+1,i))*
     |      (ra1(k,i,lat)*nop(k,i)*0.15+ra3(k,i,lat)*n2p(k,i)*1.1)
!
          n4s_loss(k,i,lat) = -(beta1(k,i,lat)*o2_cm3+
     |      beta3(k,i,lat)*no_cm3+
     |      beta17(k,i,lat)*o1_cm3*n2_cm3+
     |      rk4*o2p(k,i))

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('MBAR_N4S'  ,' ',' ',xnmbar(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('N4S_PROD',' ',' ',n4s_prod(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('N4S_LOSS',' ',' ',n4s_loss(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
 
      end subroutine comp_n4s
!-----------------------------------------------------------------------
      subroutine minor_n4s(tn,un,vn,o2,o1,he,mbar,barm,xnmbar,w,
     |  n4s,n4s_nm1,n4s_out,n4s_nm1_out,lev0,lev1,lon0,lon1,lat0,lat1)
      use cons_module,only: rmass_n4s
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  un,      ! zonal wind velocity
     |  vn,      ! meridional wind velocity
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  he,      ! helium (mmr)
     |  mbar,    ! mean molecular weight
     |  barm,    ! mean molecular weight
     |  xnmbar,  ! p0*e(-z)*mbar/kT
     |  w,       ! vertical wind
     |  n4s,     ! N4S (mmr)
     |  n4s_nm1  ! N4S at time n-1 
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  n4s_out,    ! N4S output 
     |  n4s_nm1_out ! N4S output at time n-1
!
!     write(6,"('enter minor_n4s')")
!
! Minor returns n4s_out and n4s_nm1_out. Module data n4s_prod,
! n4s_loss, etc, were defined by comp_n4s.  
!
      call minor(tn,un,vn,o2,o1,he,mbar,barm,xnmbar,w,n4s,n4s_nm1,
     |  n4s_out,n4s_nm1_out,n4s_loss,n4s_prod,n4s_lbc,n4s_ubc,
     |  rmass_n4s,phi_n4s,alfa_n4s,lev0,lev1,lon0,lon1,lat0,lat1,
     |  0,'N4S')

      end subroutine minor_n4s
!
      end module n4s_module
