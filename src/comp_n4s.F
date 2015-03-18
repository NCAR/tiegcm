!
      module n4s_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlevp1,nlonp4,nlat
      use addfld_module,only: addfld
      implicit none
!
! Boundary conditions, production and loss for N4S are defined
! by comp_n4s, and referenced by minor_n4s. Comp_n4s is called
! from a latitude loop in dynamics. After comp_n4s, dynamics calls
! minor_n4s, which passes this module data to sub minor. Sub 
! minor contains 3d mpi calls and its own latitude loops.
!
!     real,dimension(nlonp4,nlat)   :: n4s_ubc ! upper boundary
!     real,dimension(nlonp4,3,nlat) :: n4s_lbc ! lower boundary
!     real,dimension(nlevp1,nlonp4,nlat) ::
!    |  n4s_prod,     ! production of n4s
!    |  n4s_loss      ! loss of n4s
!
! Boundary conditions and production and loss terms are allocated
! subdomains by sub alloc_n4s (called from allocdata.F).
!
      real,allocatable,dimension(:,:)   :: n4s_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: n4s_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  n4s_prod,     ! production of n4s (k,i,j)
     |  n4s_loss      ! loss of n4s       (k,i,j)
      real :: phi_n4s(3) = (/0.651, 0.731, 0.741/)
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
      subroutine comp_n4s(tn,o2,o1,n2,barm,xnmbarm,no,n2d,ne,o2p,op,n2p,
     |  nplus,nop,lev0,lev1,lon0,lon1,lat)
!
! Advance n4s by one time step. This is called from dynamics at
! each subdomain latitude.
!
      use qrj_module,only: qtef
      use cons_module,only: brn2d,p0,expz,expzmid_inv,boltz,
     |  rmass_n4s,rmassinv_n2d,rmassinv_o2,rmassinv_o1,
     |  rmassinv_n2,rmassinv_no
      use chemrates_module,only: beta1,beta3,beta4,beta5,beta7,
     |  beta8,beta17,ra1,ra3,rk2,rk4,rk6,rk8
      use fields_module,only: tlbc ! lower boundary interface level of TN
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  barm,    ! mean molecular weight
     |  xnmbarm, ! p0*e(-z)*barm/kT
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
      integer ::i0,i1,nk,nkm1
      real :: xnmbar_lbc
!
!     write(6,"('enter comp_n4s: lat=',i2)") lat
      i0 = lon0
      i1 = lon1
      nk = lev1-lev0+1
      nkm1 = nk-1
!
! Lower boundary:
! n4s_lbc(:,1)=A, n4s_lbc(:,2)=B, n4s_lbc(:,3)=C define lower boundary 
!   condition where A*DPSX/DZ + B*PSX + C = 0.
! Note tn(lev1,i) assumes tn bottom boundary is stored in top slot.
!
      do i=lon0,lon1
!       xnmbar_lbc = p0*expz(1)*expzmid_inv*barm(lev0,i)/
!    |    (boltz*tn(lev1,i))
        xnmbar_lbc = p0*expz(1)*expzmid_inv*barm(lev0,i)/
     |    (boltz*tlbc(i,lat))
!
! Value at bottom given by photochemical equilibrium.
        n4s_lbc(i,1,lat) = 0.
        n4s_lbc(i,2,lat) = 1.
        n4s_lbc(i,3,lat) = -rmass_n4s/xnmbar_lbc*
     |   (qtef(lev0,i,lat)*(1.-brn2d)/xnmbar_lbc+
     |    n2d(lev0,i)*rmassinv_n2d*
     |    (beta4*xnmbar_lbc*o1(lev0,i)*rmassinv_o1+
     |     beta5(lev0,i,lat)*ne(lev0,i)+beta7)+
     |    beta8(lev0,i,lat)*no(lev0,i)*rmassinv_no)/
!
     |   (beta1(lev0,i,lat)*o2(lev0,i)*rmassinv_o2+
     |    beta3(lev0,i,lat)*no(lev0,i)*rmassinv_no+
     |    xnmbar_lbc*beta17(lev0,i,lat)*
     |    o1(lev0,i)*rmassinv_o1*n2(lev0,i)*rmassinv_n2)
!
! Zero diffusive flux at top:
        n4s_ubc(i,lat) = 0.
      enddo ! i=lon0,lon1
!
! Sources:
      do i=lon0,lon1
        do k=lev0,lev1-1
!
          n4s_prod(k,i,lat) = (.5*(qtef(k,i,lat)+qtef(k+1,i,lat))*
     |      (1.-brn2d))+xnmbarm(k,i)*(n2d(k,i)*rmassinv_n2d*
     |      (xnmbarm(k,i)*beta4*o1(k,i)*rmassinv_o1+
     |      beta5(k,i,lat)*.5*(ne(k,i)+ne(k+1,i))+beta7)+
     |      .5*(beta8(k,i,lat)+beta8(k+1,i,lat))*no(k,i)*rmassinv_no)+
     |      xnmbarm(k,i)*(rk2(k,i,lat)*op(k,i)*n2(k,i)*rmassinv_n2+
     |      rk6*nplus(k,i)*o2(k,i)*rmassinv_o2+
     |      rk8*nplus(k,i)*o1(k,i)*rmassinv_o1)+
     |      sqrt(ne(k,i)*ne(k+1,i))*(ra1(k,i,lat)*nop(k,i)*0.15+
     |      ra3(k,i,lat)*n2p(k,i)*1.1)
!
          n4s_loss(k,i,lat) = -xnmbarm(k,i)*(beta1(k,i,lat)*o2(k,i)*
     |      rmassinv_o2+beta3(k,i,lat)*no(k,i)*rmassinv_no+
     |      xnmbarm(k,i)*beta17(k,i,lat)*o1(k,i)*rmassinv_o1*
     |      n2(k,i)*rmassinv_n2)-rk4*o2p(k,i)

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('MBAR_N4S'  ,' ',' ',xnmbarm(lev0:lev1-1,i0:i1),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('N4S_PROD',' ',' ',n4s_prod(lev0:lev1-1,i0:i1,lat),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('N4S_LOSS',' ',' ',n4s_loss(lev0:lev1-1,i0:i1,lat),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
 
      end subroutine comp_n4s
!-----------------------------------------------------------------------
      subroutine minor_n4s(tn,o2,o1,n2,he,n4s,n4s_nm1,n4s_out,
     |  n4s_nm1_out,lev0,lev1,lon0,lon1,lat0,lat1)
      use cons_module,only: rmass_n4s
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  he,      ! helium (mmr)
     |  n4s,     ! N4S (mmr)
     |  n4s_nm1  ! N4S at time n-1 
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  n4s_out,    ! N4S output 
     |  n4s_nm1_out ! N4S output at time n-1
!
! Local:
      integer :: lat
!
!     write(6,"('enter minor_n4s')")
!
! Minor returns n4s_out and n4s_nm1_out. Module data n4s_prod,
! n4s_loss, etc, were defined by comp_n4s.  
!
      call minor(tn,o2,o1,n2,he,n4s,n4s_nm1,n4s_out,n4s_nm1_out,
     |  n4s_loss,n4s_prod,n4s_lbc,n4s_ubc,rmass_n4s,phi_n4s,
     |  alfa_n4s,lev0,lev1,lon0,lon1,lat0,lat1,0,'N4S')

!     do lat=lat0,lat1
!       call addfld('N4S_OUT' ,' ',' ',
!    |    n4s_out(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('N4S_TM1' ,' ',' ',
!    |    n4s_nm1_out(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
 
      end subroutine minor_n4s
!
      end module n4s_module
