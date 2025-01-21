!
      module no_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlevp1,nlonp4,nlat
      use addfld_module,only: addfld
      implicit none
!
! Boundary conditions and production and loss terms are allocated
! subdomains by sub alloc_no (called from allocdata.F).
!
      real,allocatable,dimension(:,:)   :: no_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: no_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  no_prod,     ! production of no (k,i,j)
     |  no_loss      ! loss of no       (k,i,j)
!
! NO~He coefficient (0.662) is from Langenberg et al. (2020)
! Technical note: Determination of binary gas phase diffusion
! coefficients of unstable and adsorbing atmospheric trace gases at
! low temperature - Arrested Flow and Twin Tube method
      real :: phi_no(4) = (/0.814, 0.866, 0.662, 0.926/)
      real,parameter :: alfa_no = 0. ! thermal diffusion coefficient
!
      contains
!-----------------------------------------------------------------------
      subroutine alloc_no(lon0,lon1,lat0,lat1)
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
! Allocate (i,j) subdomains to boundary conditions:
      allocate(no_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_no: error allocating',
     |  ' no_ubc: stat=',i3)") istat
      allocate(no_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_no: error allocating',
     |  ' no_lbc: stat=',i3)") istat
!
! Allocate (k,i,j) subdomains to production and loss:
      allocate(no_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_no: error allocating',
     |  ' no_prod: stat=',i3)") istat
      allocate(no_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_no: error allocating',
     |  ' no_prod: stat=',i3)") istat

!     write(6,"('alloc_no: allocated module data')")
!
      end subroutine alloc_no
!-----------------------------------------------------------------------
      subroutine comp_no(o2,o1,n2,xnmbar,xnmbari,n4s,n2d,o2p,
     |  lev0,lev1,lon0,lon1,lat)
!
! Advance NO by one time step:
!
      use cons_module,only: nob,rmass_no,
     |  rmassinv_o2,rmassinv_n4s,rmassinv_n2d,rmassinv_o1,rmassinv_n2
      use chemrates_module,only: beta1,beta2,beta3,beta6,beta8,beta9,
     |  beta9n,beta17,rk5
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  o2, o1, n2, ! o2, o, n2 mass mixing ratios
     |  xnmbar,     ! p0*e(-z)*mbar/kT
     |  xnmbari,    ! p0*e(-z)*barm/kT
     |  n4s,        ! N(4S)
     |  n2d,        ! N(2D)
     |  o2p         ! O2+
!
! Local:
      integer :: k,i
      integer ::i0,i1,nk,nkm1
!
      i0 = lon0
      i1 = lon1
      nk = lev1-lev0+1
      nkm1 = nk-1
!
! Diffusive equilibrium at upper boundary:
      do i=lon0,lon1
        no_ubc(i,lat) = 0.
!
! Lower boundary:
!
        no_lbc(i,1,lat) = 0.
        no_lbc(i,2,lat) = 1.
        no_lbc(i,3,lat) = -nob(lat)*rmass_no/xnmbari(lev0,i)
      enddo ! i=lon0,lon1
!
! Sources:
      do i=lon0,lon1
        do k=lev0,lev1-1
!
          no_prod(k,i,lat) = xnmbar(k,i)**2*o2(k,i)*rmassinv_o2*
     |      (beta1(k,i,lat)*n4s(k,i)*rmassinv_n4s+
     |       beta2(k,i,lat)*n2d(k,i)*rmassinv_n2d)
!
          no_prod(k,i,lat) = no_prod(k,i,lat)+xnmbar(k,i)**3*
     |      beta17(k,i,lat)*o1(k,i)*rmassinv_o1*n2(k,i)*
     |      rmassinv_n2*n4s(k,i)*rmassinv_n4s
!
          no_loss(k,i,lat) = 
     |      -xnmbar(k,i)*(beta3(k,i,lat)*n4s(k,i)*rmassinv_n4s+
     |      beta6*n2d(k,i)*rmassinv_n2d)-
     |      .5*(beta8(k,i,lat)+beta8(k+1,i,lat)+
     |      beta9(k,i,lat)+beta9(k+1,i,lat))-rk5*o2p(k,i)-
     |      .5*(beta9n(k,i,lat)+beta9n(k+1,i,lat))

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('MBAR_NO' ,' ',' ',xnmbar(lev0:lev1-1,i0:i1),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('NO_PROD',' ',' ',no_prod(lev0:lev1-1,i0:i1,lat),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('NO_LOSS',' ',' ',no_loss(lev0:lev1-1,i0:i1,lat),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!
      end subroutine comp_no
!-----------------------------------------------------------------------
      subroutine minor_no(tn,o2,o1,he,mbar,barm,xnmbar,w,no,no_nm1,
     |  no_out,no_nm1_out,lev0,lev1,lon0,lon1,lat0,lat1)
      use cons_module,only: rmass_no
!
! Input args at full task subdomains:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,      ! neutral temperature (deg K)
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  he,      ! helium (mmr)
     |  mbar,    ! mean molecular weight
     |  barm,    ! mean molecular weight
     |  xnmbar,  ! p0*e(-z)*mbar/kT
     |  w,       ! vertical wind
     |  no,      ! NO (mmr)
     |  no_nm1   ! NO at time n-1
!
! Output args also at full task subdomains:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  no_out,    ! NO output
     |  no_nm1_out ! NO output at time n-1
!
! Minor returns no_out and no_nm1_out. Module data no_prod,
! no_loss, etc, were defined by comp_no.
!
!     subroutine minor(tn,o2,o1,fcomp,fcomp_tm1,fcomp_out,fcomp_tm1_out,
!    |  sloss,sprod,flbc,fubc,rmx,phix,alfax,lev0,lev1,lon0,lon1,lat0,
!    |  lat1,idebug)

      call minor(tn,o2,o1,he,mbar,barm,xnmbar,w,no,no_nm1,
     |  no_out,no_nm1_out,no_loss,no_prod,no_lbc,no_ubc,
     |  rmass_no,phi_no,alfa_no,lev0,lev1,lon0,lon1,lat0,lat1,0,'NO')

      end subroutine minor_no
!-----------------------------------------------------------------------
      end module no_module
