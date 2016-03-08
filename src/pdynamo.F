      module pdynamo_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Parallel dynamo module.
!
!
      use params_module,only: nmlevp1,nmlon,nmlonp1,nlon,nlonp4,nlat,
     |  nmlat,gmlat,nmlath,nlev,nlevp1,mlev0,mlev1,nmlev_diff,dlat,spval
      use addfld_module,only: addfld
      use mpitime_module,only: mpi_timer
      use mpi_module,only: mlat0,mlat1,mlon0,mlon1,lon0,lon1,lat0,lat1,
     |  mp_mag_periodic_f2d
      use heelis_module,only: phihm
      use init_module,only: istep
      use input_module,only: current_pg,current_kq
      implicit none
!
! 3d pointers to fields regridded to magnetic subdomains (i,j,k):
! (mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)
      real,pointer,dimension(:,:,:),save :: ped_mag,hal_mag,zpot_mag, 
     |  adotv1_mag, adotv2_mag, je1pg_mag, je2pg_mag
!
! 2d pointers to fields on magnetic subdomains (i,j):
! (mlon0:mlon1,mlat0:mlat1)
      real,pointer,dimension(:,:),save :: sini_mag,adota1_mag,
     |  adota2_mag,a1dta2_mag,be3_mag
!
! 2d coefficients and RHS terms for PDE on magnetic subdomains 
! including halo points:
      real,allocatable,dimension(:,:),save :: 
     |  zigm11,     ! sigma11*cos(theta0)
     |  zigmc,      ! sigmac
     |  zigm2,      ! sigma2
     |  zigm22,     ! sigma22/cos(theta0)
     |  rim1,rim2,
     |  rhs,        ! right-hand side
     |  phimsolv,   ! solution direct from solver (nhem only)
     |  phim2d      ! solution with phihm and both nhem and shem
#if defined(INTERCOMM) || defined(CISMAH)
      real,allocatable,dimension(:,:),save :: 
     |  zigm1, nsrhs
#endif
!
! 3d potential and electric field on mag subdomains (see sub pthreed):
! (mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)
      real,allocatable,dimension(:,:,:),save :: 
     |  phim3d,        ! 3d electric potential
     |  ed13d,ed23d,   ! 3d electric field for current calculations
     |  ephi3d,        ! 3d eastward electric field
     |  elam3d,        ! 3d equatorward electric field
     |  emz3d,         ! 3d upward electric field
     |  zpotm3d        ! 3d geopotential (values at all levels)   
!
! 3d electric field on geographic subdomains (see sub pefield):
! (nlevp1,lon0-2,lon1+2,lat0:lat1)
      real,allocatable,dimension(:,:,:),save :: ex,ey,ez
!
! 3d electric potential on geographic subdomains (lon0:lon1,lat0:lat1,nlevp1)
! This will be regridded from phim3d for output to history files.
      real,allocatable,dimension(:,:,:),save :: phig3d
!
! coef_rhs and coef_cofum are allocated without halo points:
      real,allocatable,dimension(:,:),save :: 
     |  coef_rhs    ! c0(:,:,10) at subdomains (10th only, i.e., rhs)
      real,allocatable,dimension(:,:,:),save :: coef_cofum  ! cofum(mlon0:1,mlat0:1,1-9)
!
! Fields at mag equator (mlon0:mlon1,mlev0:mlev1):
      real,allocatable,dimension(:,:),save :: 
     |  ped_meq, hal_meq, adotv1_meq, adotv2_meq,         ! (mlon0:mlon1,mlev0:mlev1)
     |  je1pg_meq, je2pg_meq
      real,allocatable,dimension(:,:,:),save :: fmeq_out  ! (mlon0:mlon1,mlev0:mlev1,6)
      real,allocatable,dimension(:,:,:,:),save :: fmeq_in ! (mlon0:mlon1,mlat0:mlat1,mlev0:mlev1,6)
      integer :: jlateq      ! latitude index of mag equator
      
      real,allocatable,dimension(:,:) :: ressolv ! mag subdomains w/o halos
      real :: l2norm
!
! Global longitude values near mag equator and poles for complete_integrals and rhs.
! The nf2d 6 fields are: zigm11,zigm22,zigmc,zigm2,rim1,rim2
#if defined(INTERCOMM) || defined(CISMAH)
      integer,parameter :: nf2d=7        ! 7 2d fields for CISM
#else
      integer,parameter :: nf2d=6        ! 6 2d fields
#endif
      real :: feq_jpm1(nmlonp1,2,nf2d)   ! 6 fields at 2 lats (eq-1, eq+1)
      real :: fpole_jpm2(nmlonp1,4,nf2d) ! fields at S pole+1,2 and N pole-1,2
!
! Declarations for remaining serial part of dynamo: 
! (see sub stencils, solver, and related subroutines)
!
! Global 2d fields for root task to complete non-parallel part of dynamo.
! Gathered from subdomains by in sub gather_pdyn.
      real,dimension(nmlonp1,nmlat) :: phim_glb, ! phim_glb is used by noso_crrt
     |  zigm11_glb,zigm22_glb,zigmc_glb,zigm2_glb,rhs_glb
      real,dimension(nmlonp1,nmlat,2) :: rim_glb ! pde solver output
      real,dimension(0:nmlonp1,0:nmlat+1) :: phisolv
!
! Dimensions of the 5 grid resolutions for the multi-grid PDE:
      integer,parameter ::
     |  nmlon0=nmlon+1,
     |  nmlat0=(nmlat+1)/2,
     |  nmlon1=(nmlon0+1)/2,
     |  nmlat1=(nmlat0+1)/2,
     |  nmlon2=(nmlon1+1)/2,
     |  nmlat2=(nmlat1+1)/2,
     |  nmlon3=(nmlon2+1)/2,
     |  nmlat3=(nmlat2+1)/2,
     |  nmlon4=(nmlon3+1)/2,
     |  nmlat4=(nmlat3+1)/2
!
! Unmodified coefficients for using modified mudpack:
      real,dimension(nmlon0,nmlat0,9) :: cofum
!
! Space needed for descretized coefficients of of dynamo pde at all
!   5 levels of resolution:
!
      integer,parameter ::
     |  ncee=10*nmlon0*nmlat0+9*(nmlon1*nmlat1+nmlon2*nmlat2+nmlon3*
     |    nmlat3+nmlon4*nmlat4)
!
! Coefficients are stored in 1-d array cee(ncee)
! cee transmits descretized dynamo PDE coefficients to the multi-grid
!   mudpack solver. (cee was formerly in ceee.h)
! The common block /cee_com/ is retained from earlier versions because
!   of the equivalencing below of coefficient arrays c0, c1, etc.
!
      real :: cee(ncee)
      common/cee_com/ cee
!
! The following parameters nc0,nc1,... are pointers to the beginning of
!   the coefficients for each level of resolution.
!
      integer,parameter ::
     | nc0=1,
     | nc1=nc0+10*nmlon0*nmlat0,
     | nc2=nc1+9 *nmlon1*nmlat1,
     | nc3=nc2+9 *nmlon2*nmlat2,
     | nc4=nc3+9 *nmlon3*nmlat3
!
! nc(1:6) are pointers to beginning of coefficient blocks at each of
!   5 levels of resolution:
! nc(1) = nc0, pointer to coefficients for highest resolution.
! nc(2) = nc1, pointer to coefficients at half the resolution of nc0,
!   and so on for nc(3), nc(4), nc(5), etc.
! nc(6) = ncee, the dimension of the entire cee array, containing
!   coefficients for all 5 levels of resolution.
!
      integer :: nc(6)
      integer :: mlevtop ! for saving mag fields to sechist

      real ::
     |  c0(nmlon0,nmlat0,10),
     |  c1(nmlon1,nmlat1,9),
     |  c2(nmlon2,nmlat2,9),
     |  c3(nmlon3,nmlat3,9),
     |  c4(nmlon4,nmlat4,9)
      equivalence
     |  (cee,c0),
     |  (cee(nc1),c1),
     |  (cee(nc2),c2),
     |  (cee(nc3),c3),
     |  (cee(nc4),c4)

      logical :: debug = .false.

      real,dimension(nmlonp1,nmlat0) :: pfrac  ! NH fraction of potential

      contains
!-----------------------------------------------------------------------
      subroutine pdynamo
      use mpi_module,only: mytid,mp_scatter_coeffs,mp_scatter_phim,
     |  mp_mag_halos,mp_gather_pdyn
      use diags_module,only: mkdiag_PHIM2D
!     use solver_module,only: solve_using_pcg
!
! Parallel dynamo. 
! Required fields from neutral atmos have been regridded to magnetic
!   grid by sub dynamo_inputs (called from advance).
!
! Transform needed fields to geomagnetic coordinates
! Perform field-line integrations
! Evaluate PDE coefficients and RHS
! The PDE is divided by 1/ DT0DTS (in dyncal divided by 1/cos(theta_0)
! Sigma_(phi phi) = zigm11/ rcos0s * dt0dts
! Sigma_(lam lam) = zigm22 * rcos0s / dt0dts
! Sigma_(phi lam) = +-(zigm2-zigmc)
! Sigma_(lam phi) = -+(zigm2+zigmc)
! K_(m phi)^D     =   rim(1) * dt0dts
! K_(m lam)^D     = +-rim(2) * rcos0s
!
! (Sub dynamo_inputs has already been called from advance)
!
! Local:
      integer :: i,j,n,ier,lonend
      real :: coefglb_rhs(nmlonp1,nmlat),coefglb_cofum(nmlonp1,nmlat,9)
      character(len=16) :: cofum_name
!
! Latitude index to mag equator (module data)
      jlateq = (nmlat+1)/2   ! index to mag equator (module data)
!
! Fieldline integration: 
      call mpi_timer("pdynamo_fieldline_integrals",0,0) ! name,startstop,iprint
      call fieldline_integrals
      call mpi_timer("pdynamo_fieldline_integrals",1,0)
      if (debug) write(6,"('pdynamo after fieldline_integrals')")
!
! Equatorial and polar values, hemisphere folding: 
      call mpi_timer("pdynamo_complete_integrals",0,0) ! name,startstop,iprint
      call complete_integrals
      call mpi_timer("pdynamo_complete_integrals",1,0) ! name,startstop,iprint
      if (debug) write(6,"('pdynamo after complete_integrals')")
!
! Calculate right-hand side on mag subdomains:
! (mag halos are needed in rim1 for rhs calculation)
      call mp_mag_halos(rim1,mlon0,mlon1,mlat0,mlat1,1,1)
      call mp_mag_halos(rim2,mlon0,mlon1,mlat0,mlat1,1,1) 
      call mpi_timer("pdynamo_rhspde",0,0) ! name,startstop,iprint
      call rhspde
      call mpi_timer("pdynamo_rhspde",1,0) ! name,startstop,iprint
      if (debug) write(6,"('pdynamo after rhspde')")
!
! Gather needed arrays to root task and complete dynamo serially as in 
! original code. Timing is calculated for mp_gather_pdyn, which is 
! called by sub gather_pdyn, and is the largest wallclock consumer
! of the entire parallel dynamo code.
!
      call gather_pdyn
      if (debug) write(6,"('pdynamo after gather_pdyn')")
!
! Define global mag coeffs arrays, with non-zero in shem indices 
! real :: coefglb_rhs(nmlonp1,nmlat),coefglb_cofum(nmlonp1,nmlat,9)
      coefglb_rhs = 0.   ! whole-array init
      coefglb_cofum = 0. ! whole-array init
!
! Complete dynamo serially by root task:
      if (mytid==0) then  ! root task
        call mpi_timer("stencils",0,0)
        call stencils          ! set up stencils
        call mpi_timer("stencils",1,0)
        if (debug) write(6,"('pdynamo: root task after stencils')")
!
        call mpi_timer("solver",0,0)
        call solver(cofum,c0)  ! mudpack PDE solver
        call mpi_timer("solver",1,0)
        if (debug) write(6,"('pdynamo: root task after solver')")

        coefglb_rhs  (:,1:nmlath)   = c0(:,:,10)
        coefglb_cofum(:,1:nmlath,:) = cofum(:,:,:)
      endif
!
! rim1 after solver is needed for highlat_poten. rim_glb is distributed 
! to subdomains as rim1, and mag halos set. This will overwrite rim1 from 
! fieldline_integrals, complete_integrals, etc.
!
      call mp_scatter_phim(rim_glb(:,:,1),rim1(mlon0:mlon1,
     |  mlat0:mlat1))
      if (debug) write(6,"('pdynamo: after mp_scatter_phim')")
      call mp_mag_halos(rim1,mlon0,mlon1,mlat0,mlat1,1,1)
!
! Add high latitude potential from empirical model (heelis or weimer) 
! to solution rim1, defining phim2d on mag subdomains. 
!     allocate(phim2d(mlon00:mlon11,mlat00:mlat11),stat=istat)
!
      call mpi_timer("highlat_poten",0,0)
      call highlat_poten
      call mpi_timer("highlat_poten",1,0)
      if (debug) write(6,"('pdynamo: after highlat_poten')")

!     write(6,"('pdynamo after highlat_poten: phim2d=',2e12.4)")
!    |  minval(phim2d(mlon0:mlon1,mlat0:mlat1)),
!    |  maxval(phim2d(mlon0:mlon1,mlat0:mlat1))
!
! Gather phim2d to global phim_glb on root task, for use by noso_crrt.
!
      if (current_kq > 0) then
        call mp_gather_pdyn(phim2d(mlon0:mlon1,mlat0:mlat1),
     |    mlon0,mlon1,mlat0,mlat1,phim_glb,nmlonp1,nmlat,1) 
        if (debug) write(6,"('pdynamo after mp_gather_pdyn')")
      endif

      call mkdiag_PHIM2D('PHIM2D',phim2d(mlon0:mlon1,mlat0:mlat1),
     |  mlon0,mlon1,mlat0,mlat1)
!
! Expand phim2d to phim3d, first setting mag halos in phim2d from 
! hightlat_poten. phim3d will then be the final potential from pdynamo.
!
      call mp_mag_halos(phim2d,mlon0,mlon1,mlat0,mlat1,1,1)
      if (debug) write(6,
     |  "('pdynamo after mp_mag_halos, before pthreed')")
      call mpi_timer("pthreed",0,0)
      call pthreed
      call mpi_timer("pthreed",1,0)
      if (debug) write(6,"('pdynamo after pthreed')")
!
! Distribute 2d global phisolv (mudpack solution) to subdomains 
! in phimsolv (exclude phim halos in this call).  phimsolv is used 
! for residual testing by productAX below.
!
!     call mp_scatter_phim(phisolv(1:nmlonp1,1:nmlat),
!    |  phimsolv(mlon0:mlon1,mlat0:mlat1))
!
! Scatter coeffs from global to subdomains for productAX:
!     call mp_scatter_coeffs(coefglb_rhs,coefglb_cofum,
!    |  coef_rhs,coef_cofum)
!
! matrix vector product: 
!     ressolv = 0.0         ! mag subdomains w/o halos
!
! Call productAX with phimsolv on subdomains (distributed phisolv), 
! and with halos set:
!     call productAX(coef_cofum,phimsolv,mlon0,mlon1,
!    |     mlat0,mlat1,ressolv)

!     call solve_using_pcg(phisolv,coef_cofum,
!    |  coef_rhs(mlon0:mlon1,mlat0:mlat1),10,1e-6,
!    |  mlon0,mlon1,mlat0,mlat1,nmlon,nmlat)

! calculate residual on subdomain      .
!     l2norm = 0.
!     lonend = mlon1
!     if (mlon1 == nmlon+1) lonend = mlon1-1     .
!     do j = mlat0,mlat1
!       do i = mlon0,lonend
!         ressolv(i,j) = coef_rhs(i,j)-ressolv(i,j)
!         l2norm = l2norm + ressolv(i,j)*ressolv(i,j)	  
!       enddo
!     enddo
!     write(6,*) 'L2norm (task ',mytid,' subdomain) = ',l2norm

      if (debug) write(6,"('pdynamo returning')")
      end subroutine pdynamo
!-----------------------------------------------------------------------
      subroutine dynamo_inputs(un,vn,wn,zpot,ped,hal,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Provide needed field inputs to the dynamo by regridding the fields
! from geo to mag.
!
      use params_module,only: dlev,zpbot_dyn,zibot
      use cons_module,only: h0,dtr
      use magfield_module,only: zb,bmod
      use mpi_module,only: mp_mageq,mp_periodic_f3d,mp_periodic_f2d,
     |  mp_magpole_3d
      use mpi_module,only: tasks,nmagtaski,mytid ! only for check

      use my_esmf,only: 
     |  geo_src_grid, mag_des_grid, esmf_regrid,
     |  esmf_set2d_geo, esmf_get_2dfield,
     |  esmf_set3d_geo, esmf_get_3dfield

      use my_esmf,only: ! esmf 3d (i,j,k) fields on geo grid
     |  geo_ped, geo_hal, geo_zpot, geo_adotv1, geo_adotv2,
     |  geo_je1pg, geo_je2pg
      use my_esmf,only: ! esmf 2d (i,j) fields on geo grid
     |  geo_sini, geo_adota1, geo_adota2, geo_a1dta2, geo_be3

      use my_esmf,only: ! esmf 3d (i,j,k) fields on mag grid
     |  mag_ped, mag_hal, mag_zpot, mag_adotv1, mag_adotv2,
     |  mag_je1pg, mag_je2pg, esmf_3dfields
      use my_esmf,only: ! esmf 2d (i,j) fields on mag grid
     |  mag_sini, mag_adota1, mag_adota2, mag_a1dta2, mag_be3
      use params_module,only: glat,glon

      use magpres_g,only: ! (nlevp1,lon0:lon1,lat0:lat1)
     |  je1oD_geo,  ! J_e1(pg)/D [A/cm^2] 
     |  je2oD_geo   ! J_e2(pg)/D [A/cm^2]
!
! Args:
      integer :: lev0,lev1,lon0,lon1,lat0,lat1,rc,lonbeg,lonend
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in)::
     |  un,   ! neutral zonal velocity (cm/s)
     |  vn,   ! neutral meridional velocity (cm/s)
     |  wn    ! neutral vertical velocity (cm/s)
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),
     |  intent(inout)::
     |  zpot, ! geopotential height
     |  ped,  ! pedersen conductivity (lamdas.F)
     |  hal   ! hall conductivity     (lamdas.F)
!
! Local:
      integer :: i,j,k
!
! 3d locals dimensioned mlev0:mlev1 (mlev0 = -2 at 5-deg, -5 at 2.5-deg)
      real,dimension(mlev0:mlev1,lon0:lon1,lat0:lat1) ::
     |  undyn,vndyn,wndyn,zpotdyn,peddyn,haldyn

      real,dimension(mlev0:mlev1,lon0:lon1,lat0:lat1) :: 
     |  adotv1,adotv2,je1pg,je2pg

      real,dimension(lon0:lon1,lat0:lat1) :: 
     |  sini,adota1,adota2,a1dta2,be3
      real,dimension(lon0:lon1,lat0:lat1,4) :: f2d
      real,parameter :: unitvm(nmlon)=1.
!
! Set constants:
!   rl1,rl2 are rates at which sigma1 and sigma2 decay with height
!     below bottom of model.
!
      real,parameter ::
     |  rl1 = 5.e5,
     |  rl2 = 3.e5
!
! External:
      real,external :: sddot ! util.F
!
! For esmf_set3d_geo call:
      real,dimension(mlev0:mlev1,lon0:lon1,lat0:lat1,8) :: f3d
      character(len=8) :: fnames(8)
!
!     write(6,"('Enter dynamo_inputs: un=',2e12.4,' vn=',2e12.4,
!    |  ' wn=',2e12.4,' zpot=',2e12.4,' ped=',2e12.4,' hal=',2e12.4)")
!    |  minval(un),maxval(un),minval(vn),maxval(vn),
!    |  minval(wn),maxval(wn),minval(zpot),maxval(zpot),
!    |  minval(ped),maxval(ped),minval(hal),maxval(hal)
!
! mlevtop is for saving mag fields with addfld, see below: 
!   (actually mlevtop one level below the top):
!   At 5.0-deg: mlev0,mlev1 = -2,29 mag_field(81,97,mlev0:mlev1)
!   At 2.5-deg: mlev0,mlev1 = -5,57 mag_field(81,97,mlev0:mlev1)
!
      mlevtop = mlev1+2                  ! 5.0-deg (mlev1=57)
      if (dlat == 2.5) mlevtop = mlev1+5 ! 2.5-deg (mlev1=63)
!
! First, transfer lev0:lev1 (and all geo i,j) from inputs to xxxdyn:
      do j=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            undyn(k,i,j)   = un(k,i,j)
            vndyn(k,i,j)   = vn(k,i,j)
            wndyn(k,i,j)   = wn(k,i,j)
            zpotdyn(k,i,j) = zpot(k,i,j)
            peddyn(k,i,j)  = ped(k,i,j)
            haldyn(k,i,j)  = hal(k,i,j)
            je1pg(k,i,j)   = je1oD_geo(k,i,j)
            je2pg(k,i,j)   = je2oD_geo(k,i,j)
          enddo ! k=lev0,lev1
!
! Extend neutral atmosphere inputs down to mlev0 (zpbot_dyn=-8.5)
! At 5-deg res, this will be k=0,-2, and at 2.5-deg res, k=0,-5.
! Set three equally spaced levels for Z, take U, V, and W
! to be constant, and extrapolate sigmas exponentially.
!
! 3/3/16 btf (as per Astrid and Art): Change extrapolation of Pedersen and Hall:
!   from (float(nmlev_diff-1)*rl1)) to (2.*rl1)) for peddyn, and
!   from (float(nmlev_diff-1)*rl2)) to (2.*rl2)) for haldyn.
! According to Astrid This will lead to smaller conductivities between
!   zp -8.5 to zp -7.
!
          do k=lev0-1,mlev0,-1
            zpotdyn(k,i,j) = h0+float(k+nmlev_diff-1)*
     |        (zpotdyn(lev0,i,j)-h0)/float(nmlev_diff)
            peddyn(k,i,j) = peddyn(lev0,i,j)*exp((zpotdyn(k,i,j)+
     |        zpotdyn(k+1,i,j)-zpotdyn(lev0,i,j)-zpotdyn(lev0+1,i,j))/
     |        (2.*rl1))
            haldyn(k,i,j) = haldyn(1,i,j)*exp((zpotdyn(k,i,j)+
     |        zpotdyn(k+1,i,j)-zpotdyn(lev0,i,j)-zpotdyn(lev0+1,i,j))/
     |        (2.*rl2))
            undyn(k,i,j) = undyn(lev0,i,j)
            vndyn(k,i,j) = vndyn(lev0,i,j)
            wndyn(k,i,j) = wndyn(lev0,i,j)
            je1pg(k,i,j) = 0.
            je2pg(k,i,j) = 0.
          enddo
        enddo ! i=lon0,lon1
      enddo ! j=lat0,lat1

!     do k=lev0-1,mlev0,-1
!       write(6,"('dynamo_inputs: k=',i4,' mlev0=',i4,' undyn=',2e12.4,
!    |    ' vndyn=',2e12.4,' wndyn=',2e12.4,' zpotdyn=',2e12.4,
!    |    ' peddyn=',2e12.4,' haldyn=',2e12.4)") k,mlev0,
!    |    minval(undyn(k,:,:)),maxval(undyn(k,:,:)),
!    |    minval(vndyn(k,:,:)),maxval(vndyn(k,:,:)),
!    |    minval(wndyn(k,:,:)),maxval(wndyn(k,:,:)),
!    |    minval(zpotdyn(k,:,:)),maxval(zpotdyn(k,:,:)),
!    |    minval(peddyn(k,:,:)),maxval(peddyn(k,:,:)),
!    |    minval(haldyn(k,:,:)),maxval(haldyn(k,:,:))
!     enddo
!
! real,dimension(mlev0:mlev1,lon0:lon1,lat0:lat1) ::
      do j=lat0,lat1
        call addfld('Z_DYN',' ',' ',zpotdyn(1:nlev,:,j),
     |    'ilev',lev0,lev1-1,'lon',lon0,lon1,j)
        call addfld('U_DYN',' ',' ',undyn(1:nlev,:,j),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
        call addfld('V_DYN',' ',' ',vndyn(1:nlev,:,j),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
        call addfld('W_DYN',' ',' ',wndyn(1:nlev,:,j),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
        call addfld('PED_DYN',' ',' ',peddyn(1:nlev,:,j),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
        call addfld('HAL_DYN',' ',' ',haldyn(1:nlev,:,j),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
      enddo
!
! Calculate some 2d and 3d fields:
      call calc_adotv(zpotdyn,undyn,vndyn,wndyn,
     |  adotv1,adotv2,adota1,adota2,
     |  a1dta2,be3,mlev0,mlev1,lon0,lon1,lat0,lat1)
!
! Calculate sini sin(I_m) (zb and bmod are from magfield)
!     real,dimension(nlonp1,0:nlatp1) :: xb,bmod
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2
      do j=lat0,lat1
        do i=lonbeg,lonend 
          sini(i,j) = zb(i-2,j)/bmod(i-2,j) ! sin(I_m)
        enddo
        if (lon0==1) sini(1:2,j) = 
     |    zb(nlon:nlon+1,j)/bmod(nlon:nlon+1,j) 
        if (lon1==nlonp4) sini(nlon+3:nlonp4,j) = 
     |    zb(1:2,j)/bmod(1:2,j)
      enddo
!
!     call addfld('BMOD','2d from apex_subs',' ',
!    |  bmod(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('ZB','2d from apex_subs',' ',
!    |  zb(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)

! subroutine mp_periodic_f3d(f,lev0,lev1,lon0,lon1,lat0,lat1,nf)
! real,dimension(lev0:lev1,lon0:lon1,lat0:lat1,4) :: f3d
!
      f3d(:,:,:,1) = peddyn(:,lon0:lon1,lat0:lat1)
      f3d(:,:,:,2) = haldyn(:,lon0:lon1,lat0:lat1)
      f3d(:,:,:,3) = adotv1
      f3d(:,:,:,4) = adotv2
      call mp_periodic_f3d(f3d,mlev0,mlev1,lon0,lon1,lat0,lat1,4)
      peddyn(:,lon0:lon1,lat0:lat1) = f3d(:,:,:,1)
      haldyn(:,lon0:lon1,lat0:lat1) = f3d(:,:,:,2)
      adotv1 = f3d(:,:,:,3) 
      adotv2 = f3d(:,:,:,4) 
!
! Set 3d field values on geographic source grid:
!
      fnames = (/'PED     ','HAL     ','ZPOT    ','SCHT    ',
     |           'ADOTV1  ','ADOTV2  ','JE1PG   ','JE2PG   '/)

      f3d(:,:,:,1) = peddyn
      f3d(:,:,:,2) = haldyn
      f3d(:,:,:,3) = zpotdyn
      f3d(:,:,:,4) = adotv1
      f3d(:,:,:,5) = adotv2
      f3d(:,:,:,6) = je1pg
      f3d(:,:,:,7) = je2pg

      esmf_3dfields(1) = geo_ped
      esmf_3dfields(2) = geo_hal
      esmf_3dfields(3) = geo_zpot
      esmf_3dfields(4) = geo_adotv1
      esmf_3dfields(5) = geo_adotv2
      esmf_3dfields(6) = geo_je1pg
      esmf_3dfields(7) = geo_je2pg

      call esmf_set3d_geo(esmf_3dfields,fnames,
     |  f3d,7,mlev0,mlev1,lon0,lon1,lat0,lat1)

      geo_ped    = esmf_3dfields(1)
      geo_hal    = esmf_3dfields(2)
      geo_zpot   = esmf_3dfields(3)
      geo_adotv1 = esmf_3dfields(4)
      geo_adotv2 = esmf_3dfields(5)
      geo_je1pg  = esmf_3dfields(6)
      geo_je2pg  = esmf_3dfields(7)

      peddyn  = f3d(:,:,:,1)
      haldyn  = f3d(:,:,:,2)
      zpotdyn = f3d(:,:,:,3)
      adotv1  = f3d(:,:,:,4)
      adotv2  = f3d(:,:,:,5)
      je1pg   = f3d(:,:,:,6)
      je2pg   = f3d(:,:,:,7)
!
! 2d fields need only be calculated in first timestep:
      if (istep==1) then
        f2d(:,:,1) = adota1
        f2d(:,:,2) = adota2
        f2d(:,:,3) = a1dta2
        f2d(:,:,4) = be3
        call mp_periodic_f2d(f2d,lon0,lon1,lat0,lat1,4)
        adota1 = f2d(:,:,1)
        adota2 = f2d(:,:,2)
        a1dta2 = f2d(:,:,3)
        be3    = f2d(:,:,4)
!
! Set 2d field values on geographic grid:
! (esmf fields on source grid exclude periodic points)
!
        call esmf_set2d_geo(geo_sini,  geo_src_grid,'SINI    ',
     |    sini,  lon0,lon1,lat0,lat1)
        call esmf_set2d_geo(geo_adota1,geo_src_grid,'ADOTA1  ',
     |    adota1,lon0,lon1,lat0,lat1)
        call esmf_set2d_geo(geo_adota2,geo_src_grid,'ADOTA2  ',
     |    adota2,lon0,lon1,lat0,lat1)
        call esmf_set2d_geo(geo_a1dta2,geo_src_grid,'A1DTA2  ',
     |    a1dta2,lon0,lon1,lat0,lat1)
        call esmf_set2d_geo(geo_be3,   geo_src_grid,'BE3     ',
     |    be3,   lon0,lon1,lat0,lat1)
      endif
!
! Save original fields on geographic grid to secondary histories:
!
      do j=lat0,lat1
        call addfld('PED_GEO',' ','S/m',ped(:,:,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
        call addfld('HAL_GEO',' ','S/m',hal(:,:,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
        call addfld('ZPOT_GEO',' ','cm',zpot(:,:,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! Beware, adotv1,2, etc have levels mlev0:mlev1, e.g., -2:nlevp1,
! so must specify lev0:lev1
!
        call addfld('ADOTV1_GEO',' ',' ',adotv1(lev0:lev1,:,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
        call addfld('ADOTV2_GEO',' ',' ',adotv2(lev0:lev1,:,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
        call addfld('JE1PG_GEO',' ',' ',je1pg(lev0:lev1,:,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
        call addfld('JE2PG_GEO',' ',' ',je2pg(lev0:lev1,:,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
      enddo
!     call addfld('SINI_GEO',' ',' ',sini,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('ADOTA1_GEO',' ',' ',adota1,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('ADOTA2_GEO',' ',' ',adota2,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('A1DTA2_GEO',' ',' ',a1dta2,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('BE3_GEO',' ',' ',be3,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Regrid 3d geo fields to mag grid: 
      call esmf_regrid(geo_ped  ,mag_ped,  'geo2mag',3)
      call esmf_regrid(geo_hal  ,mag_hal,  'geo2mag',3)
      call esmf_regrid(geo_zpot ,mag_zpot, 'geo2mag',3)
      call esmf_regrid(geo_adotv1 ,mag_adotv1, 'geo2mag',3)
      call esmf_regrid(geo_adotv2 ,mag_adotv2, 'geo2mag',3)
      call esmf_regrid(geo_je1pg  ,mag_je1pg,  'geo2mag',3)
      call esmf_regrid(geo_je2pg  ,mag_je2pg,  'geo2mag',3)
!
! Regrid 2d geo fields to mag grid:
      if (istep==1) then
        call esmf_regrid(geo_sini   ,mag_sini  , 'geo2mag',2)
        call esmf_regrid(geo_adota1 ,mag_adota1, 'geo2mag',2)
        call esmf_regrid(geo_adota2 ,mag_adota2, 'geo2mag',2)
        call esmf_regrid(geo_a1dta2 ,mag_a1dta2, 'geo2mag',2)
        call esmf_regrid(geo_be3    ,mag_be3   , 'geo2mag',2)
      endif
!
! Define pdynamo module data pointers to the regridded mag fields.
! First arg of esmf_get_field is input esmf field (my_esmf module), 
!   second arg is output data pointer (pdynamo module)
! (These destination grid fields have periodic points allocated and set)
!
! Get regridded 3d mag fields:
! subroutine esmf_get_3dfield(input field, output fptr, name)
! First arg is input ESMF field, second is output pointer:
!
      call esmf_get_3dfield(mag_ped   ,ped_mag,   "PED     ")
      call esmf_get_3dfield(mag_hal   ,hal_mag,   "HAL     ")
      call esmf_get_3dfield(mag_zpot  ,zpot_mag,  "ZPOT    ")
      call esmf_get_3dfield(mag_adotv1,adotv1_mag,"ADOTV1  ")
      call esmf_get_3dfield(mag_adotv2,adotv2_mag,"ADOTV2  ")
      call esmf_get_3dfield(mag_je1pg ,je1pg_mag, "JE1PG   ")
      call esmf_get_3dfield(mag_je2pg ,je2pg_mag, "JE1PG   ")
!
! Get regridded 2d mag fields:
! First arg is input ESMF field, second is output pointer:
!
      if (istep==1) then
        call esmf_get_2dfield(mag_sini  ,sini_mag  , "SINI    ")
        call esmf_get_2dfield(mag_adota1,adota1_mag, "ADOTA1  ")
        call esmf_get_2dfield(mag_adota2,adota2_mag, "ADOTA2  ")
        call esmf_get_2dfield(mag_a1dta2,a1dta2_mag, "A1A2M   ")
        call esmf_get_2dfield(mag_be3   ,be3_mag   , "BE3     ")
      endif
!
! Set mag equator:
! fmeq_in are input fields on 3d mag subdomains. 
!
! allocate(fmeq_in(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1,6)
! allocate(fmeq_out(mlon0:mlon1,mlev0:mlev1,6),stat=istat)
! allocate(ped_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
!
      fmeq_in(:,:,:,1) = ped_mag(:,:,:)
      fmeq_in(:,:,:,2) = hal_mag(:,:,:)
      fmeq_in(:,:,:,3) = adotv1_mag(:,:,:)
      fmeq_in(:,:,:,4) = adotv2_mag(:,:,:)
      fmeq_in(:,:,:,5) = je1pg_mag(:,:,:)
      fmeq_in(:,:,:,6) = je2pg_mag(:,:,:)
!
! Tasks w/ mag equator send eq data(i,k) to other tasks in their tidi:
!
      call mp_mageq(fmeq_in,fmeq_out,6,mlon0,mlon1,mlat0,mlat1,
     |  mlev0,mlev1)
!
! Output arrays now have mag equator data on longitude subdomain
!   and full column (mlon0:mlon1,mlev0:mlev1)
! These will be used in fieldline_integrals.
!
      ped_meq(:,:)    = fmeq_out(:,:,1)
      hal_meq(:,:)    = fmeq_out(:,:,2)
      adotv1_meq(:,:) = fmeq_out(:,:,3)
      adotv2_meq(:,:) = fmeq_out(:,:,4)
      je1pg_meq(:,:)  = fmeq_out(:,:,5)
      je2pg_meq(:,:)  = fmeq_out(:,:,6)
!
! Save geopotential on magnetic grid in zpotm3d, then
! limit max zpot_mag to h0 for use in fieldline integrals
! and pthreed. It is not necessary to set poles of zpotm3d
! (as in old dynamo for zpotenm3d), since pthreed does not
! reference the poles of zpotm3d).
!
! allocate(zpotm3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
!
      do k=mlev0,mlev1
        do j=mlat0,mlat1
          do i=mlon0,mlon1
            zpotm3d(i,j,k) = zpot_mag(i,j,k)
            if (zpot_mag(i,j,k) < h0) zpot_mag(i,j,k)=h0
          enddo
        enddo
      enddo
!
! Save mag fields to secondary histories 
!
! ZMAG is a mandatory mag field on secondary histories, so this
! addfld call should not be commented. 
! (note cannot use mlev0=-2 in addfld call)
!
      do j=mlat0,mlat1
        call addfld('ZMAG','ZMAG from pdynamo','km',
     |    zpot_mag(:,j,mlev0:mlev1-1)*1.e-5,'mlon',mlon0,mlon1,
     |    'imlev',1,mlevtop,j)
      enddo

      do j=mlat0,mlat1
        call addfld('PED_MAG','PED on mag grid','S/m',
     |    ped_mag(:,j,mlev0:mlev1-1),'mlon',mlon0,mlon1,
     |    'mlev',1,mlevtop,j)

        call addfld('HAL_MAG','HALL on mag grid','S/m',
     |    hal_mag(:,j,mlev0:mlev1-1),'mlon',mlon0,mlon1,
     |    'mlev',1,mlevtop,j)

        call addfld('ZPOT_MAG','ZPOT on mag grid','cm',
     |    zpot_mag(:,j,mlev0:mlev1-1),'mlon',mlon0,mlon1,
     |    'imlev',1,mlevtop,j)

        call addfld('ADOTV1_MAG','ADOTV1 on mag grid',' ',
     |    adotv1_mag(:,j,mlev0:mlev1-1),'mlon',mlon0,mlon1,
     |    'mlev',1,mlevtop,j)

        call addfld('ADOTV2_MAG','ADOTV2 on mag grid',' ',
     |    adotv2_mag(:,j,mlev0:mlev1-1),'mlon',mlon0,mlon1,
     |    'mlev',1,mlevtop,j)

        call addfld('JE1PG_MAG','JE1PG on mag grid',' ',
     |    je1pg_mag(:,j,mlev0:mlev1-1),'mlon',mlon0,mlon1,
     |    'mlev',1,mlevtop,j)

        call addfld('JE2PG_MAG','JE2PG on mag grid',' ',
     |    je2pg_mag(:,j,mlev0:mlev1-1),'mlon',mlon0,mlon1,
     |    'mlev',1,mlevtop,j)
      enddo

!     call addfld('SINI_MAG',  'SINI on mag grid',' ',
!    |  sini_mag,'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ADOTA1_MAG','ADOTA1 on mag grid',' ',
!    |  adota1_mag,'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ADOTA2_MAG','ADOTA2 on mag grid',' ',
!    |  adota2_mag,'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('A1DTA2_MAG','A1DTA2 on mag grid',' ',
!    |  a1dta2_mag,'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('BE3_MAG'   ,'BE3 on mag grid',' ',
!    |  be3_mag,'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!
!     if (current_pg > 0) then
!       do j=mlat0,mlat1
!         call addfld('JE1PG',
!    |     'JE1PG: Eastward current due to plasma pressure and gravity',
!    |      'A/cm^2',je1pg_mag(:,j,1:mlev1),'mlon',mlon0,mlon1,'mlev',
!    |      1,mlev1,j)
!         call addfld('JE2PG',
!    |     'JE2PG: Downward current due to plasma pressure and gravity',
!    |      'A/cm^2',je2pg_mag(:,j,1:mlev1),'mlon',mlon0,mlon1,'mlev',
!    |      1,mlev1,j)
!       enddo
!     endif
      end subroutine dynamo_inputs
!-----------------------------------------------------------------------
      subroutine calc_adotv(z,un,vn,wn,adotv1,adotv2,adota1,adota2,
     |  a1dta2,be3,lev0,lev1,lon0,lon1,lat0,lat1)
      use cons_module,only: r0,h0
      use getapex_module,only: dvec       ! (nlonp1,nlat,3,2)
      use getapex_module,only: dddarr     ! (nlonp1,nlat)
      use getapex_module,only: be3arr     ! (nlonp1,nlat)
      use magfield_module,only: alatm  ! (nlonp1,0:nlatp1)
!
! Args:
! (actual level args are mlev0,mlev1, and 3d arrays have these dims)
!
      integer :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in)::
     |  z,    ! geopotential height (cm)
     |  un,   ! neutral zonal velocity (cm/s)
     |  vn    ! neutral meridional velocity (cm/s)
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) ::
     |  wn    ! vertical velocity (cm/s)
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: 
     |  adotv1, adotv2
      real,dimension(lon0:lon1,lat0:lat1),intent(out) :: 
     |  adota1, adota2, a1dta2, be3
!
! Local:
      integer :: k,i,j,lonbeg,lonend,im2
      real :: r0or,rat,sinalat
      real :: clm2(lon0:lon1,lat0:lat1)
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2

      adotv1 = 0. ; adotv2 = 0.
      do j=lat0,lat1
        do i=lonbeg,lonend
          im2 = i-2
          sinalat = sin(alatm(im2,j))     ! sin(lam)
          clm2(i,j) = 1.-sinalat*sinalat  ! cos^2(lam)
          be3(i,j) = 1.e-9*be3arr(im2,j)  ! be3 is in T (be3arr in nT)
          do k=lev0,lev1-1
!
! d_1 = (R_0/R)^1.5
            r0or = r0/(r0 + 0.5*(z(k,i,j)+z(k+1,i,j))-h0)
            rat = 1.e-2*r0or**1.5 ! 1/100 conversion in cm
!
! A_1 dot V = fac( d_1(1) u + d_1(2) v + d_1(3) w
            adotv1(k,i,j) = rat*(
     |        dvec(im2,j,1,1)*un(k,i,j)+
     |        dvec(im2,j,2,1)*vn(k,i,j)+
     |        dvec(im2,j,3,1)*wn(k,i,j))
!
! Note: clm2 is being used here to represent the squared cosine of the
!  quasi-dipole latitude, not of the M(90) latitude, since the wind
!  values are aligned vertically, not along the field line.
!
            rat = rat*sqrt((4.-3.*clm2(i,j))/(4.-3.*r0or*clm2(i,j)))
!
! A_2 dot V = fac( d_2(1) u + d_2(2) v + d_2(3) w
            adotv2(k,i,j) = rat*(
     |        dvec(im2,j,1,2)*un(k,i,j)+
     |        dvec(im2,j,2,2)*vn(k,i,j)+
     |        dvec(im2,j,3,2)*wn(k,i,j))
          enddo ! k=lev0,lev1-1
!
! Calculation of adota(n) = d(n)**2/D
!   a1dta2   = (d(1) dot d(2)) /D
!
          adota1(i,j) = (dvec(im2,j,1,1)**2 + dvec(im2,j,2,1)**2 +
     |                   dvec(im2,j,3,1)**2)/dddarr(im2,j)
          adota2(i,j) = (dvec(im2,j,1,2)**2 + dvec(im2,j,2,2)**2 +
     |                   dvec(im2,j,3,2)**2)/dddarr(im2,j)
          a1dta2(i,j) = (dvec(im2,j,1,1)*dvec(im2,j,1,2) +
     |                   dvec(im2,j,2,1)*dvec(im2,j,2,2) +
     |                   dvec(im2,j,3,1)*dvec(im2,j,3,2))/dddarr(im2,j)
        enddo ! i=lon0,lon1
!
! Must specify levels 1:nlev for addfld because lev0,lev1 are 
! actually mlev0,mlev1.
!
!       call addfld('Z_ADOTV',' ',' ',z(1:nlev,lon0:lon1,j),
!    |    'lev',1,nlev,'lon',lon0,lon1,j)
!
! 3d output fields:
!       call addfld('ADOTV1',' ',' ',adotv1(1:nlev,:,j),
!    |    'lev',1,nlev,'lon',lon0,lon1,j)
!       call addfld('ADOTV2',' ',' ',adotv2(1:nlev,:,j),
!    |    'lev',1,nlev,'lon',lon0,lon1,j)

      enddo ! j=lat0,lat1
!
! 2d output fields:
!     call addfld('ADOTA1',' ',' ',adota1,'lon',lon0,lon1,
!    |  'lat',lat0,lat1,0)
!     call addfld('ADOTA2',' ',' ',adota2,'lon',lon0,lon1,
!    |  'lat',lat0,lat1,0)
!     call addfld('A1DTA2',' ',' ',a1dta2,'lon',lon0,lon1,
!    |  'lat',lat0,lat1,0)
!     call addfld('BE3',' ',' ',be3,'lon',lon0,lon1,
!    |  'lat',lat0,lat1,0)

      end subroutine calc_adotv
!-----------------------------------------------------------------------
      subroutine fieldline_integrals
      use cons_module,only: ylatm,r0,h0,re
!
! Local:
      integer :: ier,i,j,k
      real :: sinlm, clm2, absinim, ra, sqomrra, sqrra, afac, htfac
      real :: rora,del,omdel,sig1,sig2,ue1,ue2,jgr,jmp,fac,je1pg,je2pg
      real,dimension(mlon0:mlon1) :: sindm, cosdm, ram, aam, cosiam,
     |  csthdam, rtadram
      real,dimension(mlon0:mlon1,mlev0:mlev1) :: rrm, sinidm, cosidm,
     |  costhdm, rtramrm, htfunc, htfunc2, je1pg_diag,je2pg_diag
!
! Initialize coefficients:
      zigm11 = 0.
      zigm22 = 0.
      zigm2  = 0.
      zigmc  = 0.
      rim1   = 0.
      rim2   = 0.
#if defined(INTERCOMM) || defined(CISMAH)
      zigm1  = 0.
#endif 

      je1pg_diag  = 0.
      je2pg_diag  = 0.
!
! Subdomain latitude scan for field line integrals:
      do j=mlat0,mlat1
!
! Skip only poles (poles will be left 0):
! (If poles are not skipped, clm2 will be zero and it will 
!  crash at ra = r0/clm2 below)
!       if (j==1.or.j==nmlat) cycle 
!
! Skip poles and equator:
        if (j==1.or.j==nmlat.or.j==nmlat/2+1) cycle 
!
        sinlm   = sin(ylatm(j))                ! sin(lam_m)
        clm2    = 1. - sinlm*sinlm             ! cos^2(lam_m)
        absinim = abs(sinlm)/sqrt(1.-.75*clm2) ! | sin I_m |
        ra      = r0/clm2                      ! (R_E + H_0)/cos^2(lam_m)
        sqomrra = sqrt(1.-r0/ra)               ! sqrt(1/ R_0/R_A) = sin(lam_m)
        sqrra   = sqrt(r0/ra)                  ! sqrt(R_0/R_A)
        afac    = 2.*sqrt(ra-r0)               ! 2*sqrt(R_A-R_0) (afac is NaN at mag poles)
        htfac   = sqrt(ra-.75*r0)              ! sqrt(R_A -3/4*R_0)
!
        do i=mlon0,mlon1
          aam(i) = afac/abs(sini_mag(i,j))
        enddo
!
! 2*sqrt( h_A - h_0 )/ |sin I_m | w.r to reference height A(h_R)
        do k=mlev0,mlev1
          do i=mlon0,mlon1
!
! rr = r0+z-z0 radius of magnetic point
! (Note zpot_mag max was set to h0 in dynamo_inputs)
!           rrm(i,k) = r0+zpot_mag(i,j,k)-h0 ! timegcm
            rrm(i,k) = re+zpot_mag(i,j,k)    ! tiegcm
!
! rtramr = ra-r if +ive, zero otherwise
            rtramrm(i,k) = max(0.,ra-rrm(i,k))
            rtramrm(i,k) = sqrt(rtramrm(i,k))
          enddo ! i=mlon0,mlon1
        enddo ! k=mlev0,mlev1
!
! Interpolate to midpoints:
! htfunc = factor by which to multiply AAM(I)*d(sqrt(ra-r)) = ds
        do k=mlev0,mlev1-1
          do i=mlon0,mlon1
            rrm(i,k)     = .5*(rrm(i,k)+rrm(i,k+1))
            rtramrm(i,k) = rtramrm(i,k)-rtramrm(i,k+1)
            htfunc(i,k)  = sqrt(ra-.75*rrm(i,k))/htfac
            htfunc2(i,k) = htfunc(i,k)**2
          enddo ! i=mlon0,mlon1

!         write(6,"('fieldline_integrals: j=',i3,' k=',i3,' mlat0,1=',
!    |      2i3,' mlon0,1=',2i3,' rrm(:,k)=',/,(6e12.4))") j,k,mlat0,
!    |      mlat1,mlon0,mlon1,rrm(:,k)
!         write(6,"('fieldline_integrals: j=',i3,' k=',i3,' mlat0,1=',
!    |      2i3,' mlon0,1=',2i3,' rtramrm(:,k)=',/,(6e12.4))") j,k,mlat0
!    |      ,mlat1,mlon0,mlon1,rtramrm(:,k)
!         write(6,"('fieldline_integrals: j=',i3,' k=',i3,' mlat0,1=',
!    |      2i3,' mlon0,1=',2i3,' htfunc(:,k)=',/,(6e12.4))") j,k,mlat0,
!    |      mlat1,mlon0,mlon1,htfunc(:,k)
!         write(6,"('fieldline_integrals: j=',i3,' k=',i3,' mlat0,1=',
!    |      2i3,' mlon0,1=',2i3,' htfunc2(:,k)=',/,(6e12.4))") j,k,mlat0
!    |      ,mlat1,mlon0,mlon1,htfunc2(:,k)

        enddo ! k=mlev0,mlev1
!
! Compute integrals: 
        do k=mlev0,mlev1-1
          do i=mlon0,mlon1
!
! (R_E+h)/(R_E+h_A) < 1 -> h_A > h
            rora = min(1.,rrm(i,k)/ra)
!
! (lam_m - lam) / lam_m =
! sqrt(1-r_0/r_A)sqrt(r/r_A) - sqrt(r_0/r_A)sqrt(1-r/r_A)
            del = (sqomrra*sqrt(rora)-sqrra*sqrt(1.-rora))/
     |        abs(ylatm(j))
            omdel = 1. - del
!
! Interpolate conductivities and winds in latitude along field line, assuming
!   linear variation between foot of field line and magnetic equator.
!   (For field lines other than those near the magnetic equator, del is nearly
!   zero, so that the interpolated values are essentially the values for the
!   latitude of the foot of the field line; inaccuracy of the assumption of
!   linear variation is thus unimportant for these field lines.)
!
! Here, mag equator ped_meq, etc. are from mp_mageq, called from dynamo inputs:
            sig1 = omdel*ped_mag(i,j,k) + del*ped_meq(i,k)
            sig2 = omdel*hal_mag(i,j,k) + del*hal_meq(i,k)
            ue1  = omdel*adotv1_mag(i,j,k) + del*adotv1_meq(i,k)
            ue2  = omdel*adotv2_mag(i,j,k) + del*adotv2_meq(i,k)

!           write(6,"('fieldline_integrals: i=',i3,' k=',i3,' j=',i3,
!    |        ' omdel=',e12.4,' del=',e12.4,' sig1,2=',2e12.4,' ue1,2=',
!    |        2e12.4)") i,k,j,omdel,del,sig1,sig2,ue1,ue2
!
! height varying factors: ds = aam*htfunc
!    d_1^2/D   = 1/htfunc * adota1_mag(i,j)
!    d_2^2/D   = htfunc   * adota2_mag(i,j)
!    d_1*d_2/D = 1        * a1dta2_mag(i,j)
!
! zigm11: int (sigma_p*d_1^2/D) ds : d_1^2/D
! zigm22: int (sigma_p*d_2^2/D) ds : d_2^2/D
!
            zigm11(i,j) = zigm11(i,j) + sig1*rtramrm(i,k)
            zigm22(i,j) = zigm22(i,j) + sig1*rtramrm(i,k)*htfunc2(i,k)
!
! zigmc: int (sigma_p*d_1*d_2/D) ds
! zigm2: int (sigma_h) ds
!
            zigmc(i,j) = zigmc(i,j) + sig1*rtramrm(i,k)*htfunc(i,k)
            zigm2(i,j) = zigm2(i,j) + sig2*rtramrm(i,k)*htfunc(i,k)
#if defined(INTERCOMM) || defined(CISMAH)
            zigm1(i,j) = zigm1(i,j) + sig1*rtramrm(i,k)*htfunc(i,k)
#endif

! 
! rim1: int [sigma_p*d_1^2/D u_e2+(sigma_h-(sigma_p*d_1*d_2)/D) u_e1] ds
! rim2: int [(sigma_h+sigma_p*d_1*d_2/D) u_e2-sigma_p*d_2^2/D u_e1 ] ds
! 
            rim1(i,j) = rim1(i,j) + (sig1*adota1_mag(i,j)*ue2 +
     |       (sig2 - sig1*a1dta2_mag(i,j))*htfunc(i,k)*ue1)*rtramrm(i,k)
            rim2(i,j) = rim2(i,j) + (sig1*adota2_mag(i,j)*htfunc2(i,k)*
     |        ue1 - (sig2 + sig1*a1dta2_mag(i,j))*htfunc(i,k)*ue2)*
     |        rtramrm(i,k)
!     
! rim1 add = 1/B_e3 * int [-k(ti+te)*grad Ne + rho*g]xB / B^2*d1/D ds
! rim2 add = 1/B_e3 * int [-k(ti+te)*grad Ne + rho*g]xB / B^2*d2/D ds
! je1pg & je2pg [A/cm^2]
!     
            if (current_pg > 0) then
              fac = be3_mag(i,j) 
              fac = rtramrm(i,k)*htfunc(i,k)/fac !   ds/be3/aam
              je1pg = omdel*je1pg_mag(i,j,k) + del*je1pg_meq(i,k)
              je2pg = omdel*je2pg_mag(i,j,k) + del*je2pg_meq(i,k)
              rim1(i,j) = rim1(i,j)+ 1.e4*je1pg*fac ! 1e4 for A/cm^2 -> A/m^2
              rim2(i,j) = rim2(i,j)+ 1.e4*je2pg*fac ! 1e4 for A/cm^2 -> A/m^2
	    endif
          enddo ! i=mlon0,mlon1
        enddo ! k=mlev0,mlev1
!
!       call addfld('ZIGM11_a','ZIGM11_a integrals before *.01','S',
!    |    zigm11(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('ZIGM22_a','ZIGM22_a integrals before *.01','S',
!    |    zigm22(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('ZIGMC_a','ZIGMC_a integrals before *.01','S',
!    |    zigmc(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('ZIGM2_a','ZIGM2_a integrals before *.01','S',
!    |    zigm2(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('RIM1_a','RIM1_a integrals before *.01','A/m',
!    |    rim1(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('RIM2_a','RIM2_a integrals before *.01','A/m',
!    |    rim2(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!
! Complete calculation and place result in /coefm/ zigm's
!   rim's are in A/m multiply by 1/100 to convert from [cm] to [m]
!
! At this point, 
!   zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
!   zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
!   zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
!   zigm2  is int[sigma_h] ds,         i.e. Sigma_h
!
!   rim1 is int[(sigma_h-sigma_p*d_1*d_2/D)u_e1 + sigma_p*d_1^2/D u_e2] *A(h_r)*
!                B_e3 ds, i.e.  K_(m phi)^D/abs(sin Im)
!   rim2 is int[(sigma_h+sigma_p*d_1*d_2/D)u_e2 - sigma_p*d_2^2/D u_e1] *A(h_r)*
!                B_e3 ds, K_(m lam)^D ( minus in northern hemisphere
!   Change sign of RIM(2) in S. hemisphere to be compatible with transf
! At this point, rim2 is +-K_(m lam)^D
!
        do i = mlon0,mlon1
          zigm11(i,j) = 1.e-2*zigm11(i,j)*aam(i)*adota1_mag(i,j)
          zigm22(i,j) = 1.e-2*zigm22(i,j)*aam(i)*adota2_mag(i,j)
          zigmc(i,j)  = 1.e-2*zigmc (i,j)*aam(i)*a1dta2_mag(i,j)
          zigm2(i,j)  = 1.e-2*zigm2 (i,j)*aam(i)
#if defined(INTERCOMM) || defined(CISMAH)
          zigm1(i,j)  = 1.e-2*zigm1(i,j)*aam(i)
#endif
          rim1(i,j)   = 1.e-2*rim1(i,j)*aam(i)*be3_mag(i,j)
          rim2(i,j)   = 1.e-2*rim2(i,j)*aam(i)*be3_mag(i,j)
        enddo ! i = 1,nmlon
!
! 12/9/11 btf: Similiar diffs w/ odyn as above, with addition of
! some diffs apparently in geo pole positions (presumably from
! diffs in adota1,2,be3) Overall structure is ok, diffs not visible
! in full-field plots.
!
!       call addfld('ZIGM11_b','ZIGM11_b integrals after *.01','S',
!    |    zigm11(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('ZIGM22_b','ZIGM22_b integrals after *.01','S',
!    |    zigm22(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('ZIGMC_b','ZIGMC_b integrals after *.01','S',
!    |    zigmc(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('ZIGM2_b','ZIGM2_b integrals after *.01','S',
!    |    zigm2(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('RIM1_b','RIM1_b integrals after *.01','A/m',
!    |    rim1(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)
!       call addfld('RIM2_b','RIM2_b integrals after *.01','A/m',
!    |    rim2(mlon0:mlon1,j),'mlon',mlon0,mlon1,'mlat',j,j,0)

      enddo ! j=mlat0,mlat1 (without poles)
!
      end subroutine fieldline_integrals
!-----------------------------------------------------------------------
      subroutine alloc_pdyn
!
! Allocate and initialize arrays for parallel dynamo (module data)
! (called once per run from allocdata.F)
!
      integer :: istat,n
      integer :: mlon00,mlon11,mlat00,mlat11
!
      mlon00=mlon0-1 ; mlon11=mlon1+1
      mlat00=mlat0-1 ; mlat11=mlat1+1
!
! 2d fields on mag subdomains (i,j):
! Certain fields are allocated with halos mlon0-1:mlon1+1,mlat0-1:mlat1+1
!
      allocate(zigm11(mlon00:mlon11,mlat00:mlat11),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: zigm11')
      zigm11 = 0.
      allocate(zigmc(mlon00:mlon11,mlat00:mlat11) ,stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: zigmc')
      zigmc = 0.
      allocate(zigm2(mlon00:mlon11,mlat00:mlat11) ,stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: zigm2')
      zigm2 = 0.
      allocate(zigm22(mlon00:mlon11,mlat00:mlat11),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: zigm22')
      zigm22 = 0.
      allocate(rhs(mlon00:mlon11,mlat00:mlat11)   ,stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: rhs')
      rhs = 0.
      allocate(rim1(mlon00:mlon11,mlat00:mlat11)  ,stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: rim1')
      rim1 = 0.
      allocate(rim2(mlon00:mlon11,mlat00:mlat11)  ,stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: rim2')
      rim2 = 0.
      allocate(phimsolv(mlon00:mlon11,mlat00:mlat11),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: phimsolv')
      phimsolv = 0.
      allocate(phim2d(mlon00:mlon11,mlat00:mlat11),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: phim2d')
      phim2d = 0.
#if defined(INTERCOMM) || defined(CISMAH)
      allocate(zigm1(mlon00:mlon11,mlat00:mlat11),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: zigm1')
      zigm1 = 0.
      allocate(nsrhs(mlon00:mlon11,mlat00:mlat11),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: nsrhs')
      nsrhs = 0.
#endif
!
! 3d phim and electric field on mag subdomains:
      allocate(phim3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: phim3d')
      phim3d = 0.
      allocate(ed13d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: ed13d')
      ed13d = 0.
      allocate(ed23d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: ed23d')
      ed23d = 0.
      allocate(ephi3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: ephi3d')
      ephi3d = 0.
      allocate(elam3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: elam3d')
      elam3d = 0.
      allocate(emz3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: emz3d')
      emz3d = 0.
      allocate(zpotm3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: zpotm3d')
      zpotm3d = 0.
!
! 3d electric field on geographic subdomains:
      allocate(ex(nlevp1,lon0-2:lon1+2,lat0:lat1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: ex')
      allocate(ey(nlevp1,lon0-2:lon1+2,lat0:lat1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: ey')
      allocate(ez(nlevp1,lon0-2:lon1+2,lat0:lat1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: ez')
      ex=0. ; ey=0. ; ez=0.
!
! 3d electric potential on geographic subdomains (regridded from phim3d):
      allocate(phig3d(lon0:lon1,lat0:lat1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: phig3d')
      phig3d = 0.

      allocate(coef_rhs(mlon0:mlon1,mlat0:mlat1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: coef_rhs')
      coef_rhs = 0.

      allocate(coef_cofum(mlon0:mlon1,mlat0:mlat1,1:9),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: coef_cofum')
      coef_cofum = 0.
!
      allocate(ped_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: ped_meq')
      ped_meq = 0.
      allocate(hal_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: hal_meq')
      hal_meq = 0.
      allocate(adotv1_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: adotv1_meq')
      adotv1_meq = 0.
      allocate(adotv2_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: adotv2_meq')
      adotv2_meq = 0.
      allocate(je1pg_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: je1pg_meq')
      je1pg_meq = 0.
      allocate(je2pg_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: je2pg_meq')
      je2pg_meq = 0.
!
! Fields input to mp_mageq (6 fields at full mag subdomain i,j,k):

      write(6,"('alloc_pdyn: allocate fmeq_in with mlev0,1=',2i4)")
     |  mlev0,mlev1

      allocate(fmeq_in(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1,6),
     |  stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: fmeq_in')
      fmeq_in = 0.
!
! Fields output by mp_mageq (6 fields at mag subdomain i,k)
      allocate(fmeq_out(mlon0:mlon1,mlev0:mlev1,6),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: fmeq_out')
      fmeq_out = 0.

! am 9/14/2010      
! for L2norm (only for checking not necessary and can be removed later)     
      allocate(ressolv(mlon0:mlon1,mlat0:mlat1),stat=istat)
      if (istat /= 0) call shutdown('alloc_pdyn: ressolv')
      ressolv = 0. 
        
      end subroutine alloc_pdyn
!-----------------------------------------------------------------------
      subroutine complete_integrals
!
! Field line integrals for each hemisphere have been calculated in
!   mag subdomains. Now, complete these arrays with equator and polar
!   values, and sum the integrals from the 2 hemispheres.
!   (most of this was in sub transf in old serial dynamo)
! This is done by obtaining global mag fields via mpi exchange
!   of mag subdomains, completing the global fields, and updating
!   the subdomains.
! The 6 2d fields are: zigm11,zigm22,zigmc,zigm2,rim1,rim2
!
      use cons_module,only: rcos0s,dt1dts
      use mpi_module,only: mp_mageq_jpm1,mp_magpole_2d,mp_mag_foldhem
      use mpi_module,only: mytid
#if defined(INTERCOMM) || defined(CISMAH)
      use mpi_module,only: mp_mag_halos
      use cism_coupling_module,only: gzigm1,gzigm2,gnsrhs,cism_ucurrent
      use my_esmf,only: 
     |  mag_zigm1, mag_zigm2, mag_nsrhs,
     |  geo_zigm1, geo_zigm2, geo_nsrhs
#endif

      implicit none
!
! Local:
      integer :: i,ii,j,ifld,lonend
      real :: fmsub(mlon0:mlon1,mlat0:mlat1,nf2d)
      real :: corfac
      real,parameter :: unitvm(nmlon)=1.
!
! External:
      real,external :: sddot ! util.F
!
! For equatorial values, we need latitudes eq+1 and eq-1:
! Local feq_jpm1(nmlonp1,2,6) is returned by mp_mageq_jpm1,
!   where the 2 dim contains lats nmlath-1, nmlath+1. These
!   are global in lon, even tho each subd uses only its own i's.
! These mag equator values do not show up on plots because
!   of the small factor .06 and .125.
! The 6 2d fields are: zigm11,zigm22,zigmc,zigm2,rim1,rim2
! Must specify mlon0:mlon1,mlat0:mlat1 because these fields
!   were allocated to include single-point halos (these calls 
!   exclude the halo points):
!
      fmsub(:,:,1) = zigm11(mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,2) = zigm22(mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,3) = zigmc (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,4) = zigm2 (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,5) = rim1  (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,6) = rim2  (mlon0:mlon1,mlat0:mlat1)
#if defined(INTERCOMM) || defined(CISMAH)
      fmsub(:,:,7) = zigm1(mlon0:mlon1,mlat0:mlat1)
#endif

      call mp_mageq_jpm1(fmsub,mlon0,mlon1,mlat0,mlat1,nmlonp1,
     |  feq_jpm1,nf2d)
!
! From sub fieldline_integrals:
!   zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
!   zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
!   zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
!   zigm2  is int[sigma_h] ds,         i.e. Sigma_h
!
!   rim1 is int[(sigma_h-sigma_p*d_1*d_2/D)u_e1 + sigma_p*d_1^2/D u_e2] *A(h_r)*
!                B_e3 ds, i.e.  K_(m phi)^D/abs(sin Im)
!   rim2 is int[(sigma_h+sigma_p*d_1*d_2/D)u_e2 - sigma_p*d_2^2/D u_e1] *A(h_r)*
!                B_e3 ds, K_(m lam)^D ( minus in northern hemisphere
!   Change sign of RIM(2) in S. hemisphere to be compatible with transf
! At this point, rim2 is +-K_(m lam)^D
!
! Equatorial values:
! Assume that quantities primarily dependent on Pedersen conductivity
!   have field-line integrals 1/4 as large as the averages for next-higher
!   field lines; quantities primarily dependent on Hall conductivity
!   have field-line integrals 0.12 as large as the averages for next-higher
!   field lines.  Exact values chosen should not be important for potential
!   calculation, as long as they are physically reasonable and not too
!   different from adjacent values.
!
      do j=mlat0,mlat1
        if (j==jlateq) then
          zigm11(mlon0:mlon1,j) = .125*(feq_jpm1(mlon0:mlon1,1,1)+
     |                                  feq_jpm1(mlon0:mlon1,2,1))
          zigm22(mlon0:mlon1,j) = .125*(feq_jpm1(mlon0:mlon1,1,2)+
     |                                  feq_jpm1(mlon0:mlon1,2,2))
          zigmc (mlon0:mlon1,j) = .125*(feq_jpm1(mlon0:mlon1,1,3)+
     |                                  feq_jpm1(mlon0:mlon1,2,3))
          zigm2 (mlon0:mlon1,j) = .060*(feq_jpm1(mlon0:mlon1,1,4)+
     |                                  feq_jpm1(mlon0:mlon1,2,4))
          rim1  (mlon0:mlon1,j) = .060*(feq_jpm1(mlon0:mlon1,1,5)+
     |                                  feq_jpm1(mlon0:mlon1,2,5))
          rim2  (mlon0:mlon1,j) = .060*(feq_jpm1(mlon0:mlon1,1,6)+
     |                                  feq_jpm1(mlon0:mlon1,2,6))
#if defined(INTERCOMM) || defined(CISMAH)
          zigm1(mlon0:mlon1,j) = .060*(feq_jpm1(mlon0:mlon1,1,7)+
     |                                  feq_jpm1(mlon0:mlon1,2,7))
#endif

!
! Include the boundary condition at the equator eq.(5.30) in
! Richmond (1995) Ionospheric Electrodynamics use. Mag. Apex Coord.
!   J.Goemag.Geoelectr. 47,191-212
! Sig_phiphi/abs(sin Im) = 0.5*Sig_cowling/abs(sin Im)
!   = 0.5/abs(sin Im)*(Sig_phiphi - Sig_philam*sig_lamphi/Sig_lamlam)
!   = 0.5/abs(sin Im)*(Sig_phiphi + (Sig_h-sig_c)*(Sig_h+sig_c)/Sig_lamlam)
!  rim(1) / |sin I_m| = I_1 = R/2*(K_mphi - Sig_philam/Sig_lamlam*K_mlam)
!
          do i=mlon0,mlon1
            zigm11(i,j) = zigm11(i,j)+(zigm2(i,j)-zigmc(i,j))*
     |        (zigm2(i,j)+zigmc(i,j))/zigm22(i,j)
            rim1(i,j) = rim1(i,j) - (zigm2(i,j)-zigmc(i,j))/
     |        zigm22(i,j)*rim2(i,j)
!
! 12/9/11 btf: current dynamo.F (tiegcm and timegcm) has the following 
!   2 "do-nothing" lines, but earlier versions (e.g., tag timegcm1-2dev7 
!   has 0.5* for both lines, like the commented lines below.
!
            zigm11(i,j) = zigm11(i,j)      ! current dynamo.F
            rim1(i,j) = rim1(i,j)          ! current dynamo.F

!           zigm11(i,j) = 0.5*zigm11(i,j)  ! old timegcm1-2dev7 tag
!           rim1(i,j) = 0.5*rim1(i,j)      ! old timegcm1-2dev7 tag

          enddo ! i=mlon0,mlon1
        endif ! j at equator
      enddo ! j=mlat0,mlat1
!
! zigm11, etc are allocated to include halo points, so must
! specify subdomains subsections here:
!
!     call addfld('ZIGM11_c','ZIGM11_c after mageq','S',
!    |  zigm11(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('ZIGM22_c','ZIGM22_c after mageq','S',
!    |  zigm22(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('ZIGMC_c','ZIGMC_c after mageq','S',
!    |  zigmc(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('ZIGM2_c','ZIGM2_c after mageq','S',
!    |  zigm2(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('RIM1_c','RIM1_c after mageq','A/m',
!    |  rim1(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('RIM2_c','RIM2_c after mageq','A/m',
!    |  rim2(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!
! Using notation of Richmond (1995) on right-hand side below:
! Sigma_(phi phi) = zigm11*abs(sin I_m)
! Sigma_(lam lam) = zigm22/abs(sin I_m)
! Sigma_(phi lam) = +-(zigm2-zigmc)
! Sigma_(lam phi) = -+(zigm2+zigmc)
! K_(m phi)^D     = rim(1)*abs(sin I_m)
! K_(m lam)^D     = +-rim(2)
!
! Transforming PDE from original apex (theta_a) to new apex grid (theta_0)
!   which is equally spaced in magnetic latitude
! SCALE quantities to modified (0) magnetic latitude system, multiplying or
!   dividing by abs(sin I_m) [inverse contained in DT1DTS] as necessary.
! Sign of K_(m lam)^D in southern hemisphere remains reversed.
! for the mixed terms the transformation from the integration and differentiation
! canceled out (zigmc, zigm2)
! DT1DTS : d theta_0/ d theta_a / abs(sin I_m)
! RCOS0S : cos(theta_0)/ cos(theta_a)
!
! corfac: abs(I_m)*d theta_a/d theta_0 * cos(theta_0)/ cos(theta_a)
! zigm11: abs(I_m)*d theta_a/d theta_0 * cos(theta_0)/ cos(theta_a)
! zigm22: 1/abs(I_m)*d theta_0/d theta_a * cos(theta_a)/ cos(theta_0)
! rim(1): abs(I_m)*d theta_a/d theta_0
! rim(2): cos(theta_a)/ cos(theta_0)
!
      do j=mlat0,mlat1
        if (j==1.or.j==nmlat) cycle   ! skip poles
        corfac = rcos0s(j)/dt1dts(j)
        do i=mlon0,mlon1
          zigm11(i,j) = zigm11(i,j)*corfac
          zigm22(i,j) = zigm22(i,j)/corfac
          rim1(i,j) = rim1(i,j)/dt1dts(j)
          rim2(i,j) = rim2(i,j)/rcos0s(j)
        enddo
      enddo
!
! For polar values, we need south pole plus 1 and 2 (j==2,3),
!   and north pole minus 1 and 2 (j==nmlat-1,nmlat-2). These
!   are returned by sub mp_magpole_jpm2 (mpi.F):
! Must specify (mlon0:mlon1,mlat0:mlat1) because zigmxx and rims
!   are allocated to include halo cells.
!
      fmsub(:,:,1) = zigm11(mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,2) = zigm22(mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,3) = zigmc (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,4) = zigm2 (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,5) = rim1  (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,6) = rim2  (mlon0:mlon1,mlat0:mlat1)
#if defined(INTERCOMM) || defined(CISMAH)
      fmsub(:,:,7) = zigm1 (mlon0:mlon1,mlat0:mlat1)
#endif
!
! mp_magpole_2d returns fpole_jpm2(nmlonp1,1->4,nf) as:
!   1: j = 2       (spole+1)
!   2: j = 3       (spole+2)
!   3: j = nmlat-1 (npole-1)
!   4: j = nmlat-2 (npole-2)
!
! subroutine mp_magpole_2d(f,ilon0,ilon1,ilat0,ilat1,
!   nglblon,jspole,jnpole,fpole_jpm2,nf)

      call mp_magpole_2d(fmsub,mlon0,mlon1,mlat0,mlat1,nmlonp1,
     |  1,nmlat,fpole_jpm2,nf2d)

!     write(6,"('complete_integrals: zigm11 global lons at j==2 = ',/,
!    |  (6e12.4))") fpole_jpm2(:,1,1)
!     write(6,"('complete_integrals: zigm11 global lons at j==nmlat-1 ='
!    |  ,/,(6e12.4))") fpole_jpm2(:,3,1)
!     write(6,"('complete_integrals: zigm22 global lons at j==2 = ',/,
!    |  (6e12.4))") fpole_jpm2(:,1,2)
!     write(6,"('complete_integrals: zigm22 global lons at j==nmlat-1 ='
!    |  ,/,(6e12.4))") fpole_jpm2(:,3,2)
!
! the PDE is divided by 1/ DT0DTS
! Sigma_(phi phi) = zigm11/ rcos0s * dt0dts
! Sigma_(lam lam) = zigm22 * rcos0s / dt0dts
! Sigma_(phi lam) = +-(zigm2-zigmc)
! Sigma_(lam phi) = -+(zigm2+zigmc)
! K_(m phi)^D     =   rim(1) * dt0dts
! K_(m lam)^D     = +-rim(2) * rcos0s
! 
! Compute polar values for the conductances, 4th order interpolation:
!
!
      do j=mlat0,mlat1
!
! South pole:
        if (j==1) then ! south pole (use fpole_jpm2(nmlon,1->2,nf)
          zigm11(mlon0,j)=(4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,1,1))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,2,1)))/
     |      (3.*float(nmlon))
          zigm22(mlon0,j)=(4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,1,2))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,2,2)))/
     |      (3.*float(nmlon))
          zigmc (mlon0,j)=(4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,1,3))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,2,3)))/
     |      (3.*float(nmlon))
          zigm2 (mlon0,j)=(4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,1,4))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,2,4)))/
     |      (3.*float(nmlon))
#if defined(INTERCOMM) || defined(CISMAH)
          zigm1(mlon0,j) = (4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,1,7))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,2,7)))/
     |      (3.*float(nmlon))
#endif
!
! Extend south pole over longitude:
          do i=mlon0+1,mlon1
            zigm11(i,j) = zigm11(mlon0,j)
            zigm22(i,j) = zigm22(mlon0,j)
            zigmc (i,j) = zigmc (mlon0,j)
            zigm2 (i,j) = zigm2 (mlon0,j)
#if defined(INTERCOMM) || defined(CISMAH)
            zigm1 (i,j) = zigm1 (mlon0,j)
#endif
          enddo ! i=mlon0,mlon1
!
! RHS vector (I_1,I_2): average over south pole:
! (use fpole_jpm2(i,1,nf), i.e. j==2, and lons across the pole)
          lonend = mlon1
          if (mlon1==nmlonp1) lonend = mlon1-1
          do i=mlon0,lonend
            ii = 1+mod(i-1+nmlon/2,nmlon)
            rim1(i,j) = 0.5*(fpole_jpm2(i,1,5)-fpole_jpm2(ii,1,5))
            rim2(i,j) = 0.5*(fpole_jpm2(i,1,6)-fpole_jpm2(ii,1,6))
          enddo
!
! North pole:
        elseif (j==nmlat) then ! north pole (use fpole_jpm2(nmlon,3->4,1,nf)
          zigm11(mlon0,j)=(4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,3,1))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,4,1)))/
     |      (3.*float(nmlon))
          zigm22(mlon0,j)=(4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,3,2))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,4,2)))/
     |      (3.*float(nmlon))
          zigmc (mlon0,j)=(4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,3,3))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,4,3)))/
     |      (3.*float(nmlon))
          zigm2 (mlon0,j)=(4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,3,4))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,4,4)))/
     |      (3.*float(nmlon))
#if defined(INTERCOMM) || defined(CISMAH)
          zigm1(mlon0,j) = (4.*
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,3,7))-
     |      sddot(nmlon,unitvm,fpole_jpm2(1:nmlon,4,7)))/
     |      (3.*float(nmlon))
#endif
!
! Extend north pole over longitude:
          do i=mlon0+1,mlon1
            zigm11(i,j) = zigm11(mlon0,j)
            zigm22(i,j) = zigm22(mlon0,j)
            zigmc (i,j) = zigmc (mlon0,j)
            zigm2 (i,j) = zigm2 (mlon0,j)
#if defined(INTERCOMM) || defined(CISMAH)
            zigm1 (i,j) = zigm1 (mlon0,j)
#endif
          enddo ! i=mlon0,mlon1
!
! RHS vector (I_1,I_2): average over north pole:
! (use fpole_jpm2(i,3,nf), i.e. j==nmlat-1, and lons across the pole)
          lonend = mlon1
          if (mlon1==nmlonp1) lonend = mlon1-1
          do i=mlon0,lonend
            ii = 1+mod(i-1+nmlon/2,nmlon)
            rim1(i,j) = 0.5*(fpole_jpm2(i,3,5)-fpole_jpm2(ii,3,5))
            rim2(i,j) = 0.5*(fpole_jpm2(i,3,6)-fpole_jpm2(ii,3,6))
          enddo
        endif ! south or north pole
      enddo ! j=mlat0,mlat1
!
!     call addfld('ZIGM11_d','ZIGM11_d after mag poles','S',
!    |  zigm11(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('ZIGM22_d','ZIGM22_d after mag poles','S',
!    |  zigm22(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('ZIGMC_d','ZIGMC_d after mag poles','S',
!    |  zigmc(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('ZIGM2_d','ZIGM2_d after mag poles','S',
!    |  zigm2(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('RIM1_d','RIM1_d after mag poles','A/m',
!    |  rim1(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('RIM2_d','RIM2_d after mag poles','A/m',
!    |  rim2(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)

#if defined(INTERCOMM) || defined(CISMAH)
!
!
! Calculate global neutral winds generated current to transfer to CISM M-I couplier
!

      call mp_mag_halos(rim1,mlon0,mlon1,mlat0,mlat1,1,1)
      call mp_mag_halos(rim2,mlon0,mlon1,mlat0,mlat1,1,1) 

      call cism_ucurrent(rim1,rim2,mlon0,mlon1,mlat0,mlat1,nsrhs)

!     call cism_ucurrent
!
! Transform height-integrated parameters from APEX coordinates to geographic
! coordinates to be passed to the M-I coupled in the CISM
!

      call mag2geo_2d(zigm1(mlon0:mlon1,mlat0:mlat1),gzigm1,mag_zigm1,
     |     geo_zigm1,'ZIGM1   ')
      call mag2geo_2d(zigm2(mlon0:mlon1,mlat0:mlat1),gzigm2,mag_zigm2,
     |     geo_zigm2,'ZIGM2   ')
      call mag2geo_2d(nsrhs(mlon0:mlon1,mlat0:mlat1),gnsrhs,mag_nsrhs,
     |     geo_nsrhs,'NSRHS   ')

#endif
!
! Set coefficients for current density calculations.
! Note nosocoef is a stub in current.F90 (but not in the
! current module) which calls noso_coef. This is necessary
! to avoid circular dependency between the pdynamo and current
! modules. Subs noso_crrt and noso_crdens (current.F90) are 
! called from advance after pdynamo.
!
      if (current_kq > 0) then
        call gather_pdyn
        if (mytid==0) call nosocoef
      endif
!
! Fold south hemisphere over onto north, and set periodic points:
      fmsub(:,:,1) = zigm11(mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,2) = zigm22(mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,3) = zigmc (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,4) = zigm2 (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,5) = rim1  (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,6) = rim2  (mlon0:mlon1,mlat0:mlat1)
#if defined(INTERCOMM) || defined(CISMAH)
      fmsub(:,:,7) = zigm1 (mlon0:mlon1,mlat0:mlat1)
#endif
      call mp_mag_foldhem(fmsub,mlon0,mlon1,mlat0,mlat1,nf2d)
      call mp_mag_periodic_f2d(fmsub,mlon0,mlon1,mlat0,mlat1,nf2d) 

      zigm11(mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,1)
      zigm22(mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,2)
      zigmc (mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,3)
      zigm2 (mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,4)
      rim1  (mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,5)
      rim2  (mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,6)
#if defined(INTERCOMM) || defined(CISMAH)
      zigm1 (mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,7)
#endif
!
! Reverse sign of zigmc in northern hemisphere.
      do j=mlat0,mlat1
        if (j >= nmlath) then
          zigmc(mlon0:mlon1,j) = -zigmc(mlon0:mlon1,j)
        endif
      enddo

!     do j=mlat0,mlat1
!       if (j==nmlath) write(6,"('j=nmlath=',i3,' mlon0,1=',2i3,
!    |    ' zigmc(mlon0:mlon1,j) after sign reversal =',/,(6e12.4))") 
!    |    j,mlon0,mlon1,zigmc(mlon0:mlon1,j)
!     enddo

!     call addfld('ZIGM11_e','ZIGM11_e after fold hem','S',
!    |  zigm11(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('ZIGM22_e','ZIGM22_e after fold hem','S',
!    |  zigm22(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('ZIGMC_e','ZIGMC_e after fold hem','S',
!    |  zigmc(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('ZIGM2_e','ZIGM2_e after fold hem','S',
!    |  zigm2(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('RIM1_e','RIM1_e after fold hem','A/m',
!    |  rim1(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!     call addfld('RIM2_e','RIM2_e after fold hem','A/m',
!    |  rim2(mlon0:mlon1,mlat0:mlat1),'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
!
      end subroutine complete_integrals
!-----------------------------------------------------------------------
      subroutine rhspde
      use cons_module,only: pi_dyn,dlatm,dlonm,r0
      use cons_module,only: rcos0s,dt1dts
!
! Calculate right-hand side from rim1,2 on mag subdomains.
! Use global longitude arrays for poles and equator obtained
!   by sub complete_integrals.
!
! Local:
      integer :: j,i,ii
      real :: pi
      real,dimension(nmlat) :: tint1
      real,parameter :: unitvm(nmlon)=1.
      real :: rim2_npm1(nmlonp1), ! global rim2 at nmlat-1
     |        rim2_eqp1(nmlonp1), ! global rim2 at meq+1
     |        rim2_meq (nmlonp1), ! global rim2 at mag eq
     |        rim2_tmp (nmlonp1), ! tmp array
     |        rim1_meq (nmlonp1), ! global rim1 at mag eq
     |        zigm2_meq(nmlonp1), ! needed for rim1_meq
     |        zigmc_meq(nmlonp1), ! needed for rim1_meq
     |        zigm22_meq(nmlonp1) ! needed for rim1_meq
!     real :: fsave(mlon0:mlon1,mlat0:mlat1)
!
! External:
      real,external :: sddot ! util.F
!
!     call addfld('RIM1',' ','S',rim1(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('RIM2',' ','S',rim2(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)

      pi = pi_dyn
      do j=1,nmlat
        tint1(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo
!
! Init rhs subdomains:
      rhs(:,:) = 0.
!
! Will need rim2 at npole-1 and mag equator:
! rim2_npm1: global rim2 at nmlat-1: 
      rim2_npm1(:) = fpole_jpm2(:,3,6)+fpole_jpm2(:,1,6)
! rim2_meq: global rim2 at mag equator jlateq:
      rim2_meq(:) = .060*(feq_jpm1(:,1,6)+feq_jpm1(:,2,6))
      rim2_meq(:) = rim2_meq(:)/rcos0s(jlateq)
      rim2_meq(:) = rim2_meq(:)*2. ! fold eq on itself
!
! Perform differentiation of rim(2) w.r.t. lam_0:
!  +/- (d [ K_(m lam)^D * cos(lam_m)]/ d lam_0 ) /cos ( lam_0) =
!  + (d [ K_(m lam)^D * cos(lam_m)]/ d lam_m ) /cos ( lam_m) / (RCOS0S*DT0DTS) =
!  +/- (d [ K_(m lam)^D(0) * cos(lam_0)]/ d lam_0 ) /cos ( lam_0) =
!
! Lat scan to define rhs subdomains:
      do j=mlat0,mlat1
!
! North Pole (redundant in longitude):
        if (j == nmlat) then       ! north pole
          do i=mlon0,mlon1
            rhs(i,j) = -2./float(nmlon)*sddot(nmlon,unitvm,rim2_npm1)/
     |        tint1(nmlat-1)
          enddo
!
! Include the boundary condition at the equator.
! rhs(equator)/R = d (K_mphi^DT(0) - sig_philam/sig_lamlam*K_mlam^DT(0)) / d phi_m
!                + d (cos lam_0 * K_mlam^DT(0))/ d lam_0
! from Cicely's notes:
! I_1 = 0.5*(K_(m phi)^DT(0) - Sig_(phi lam)/Sig_(lam lam)*K_(ml am)^DT(0))
! I_2 = K_(m lam)^DT(0)
! differentiate
! rhs = (I_1(i+1/2,j)-I_1(i-1/2,j))/dlonm +
!       (2*cos(lam_0)_(j+1/2)*I_2(i,j+1/2))/dlat_0
! (first calc global mag equator as in complete_integrals)
!
        elseif (j == jlateq) then  ! mag equator
          rim2_eqp1(:) = feq_jpm1(:,2,6)/rcos0s(j-1)+
     |                   feq_jpm1(:,1,6)/rcos0s(j+1)
          zigm22_meq(:)  = .125*(feq_jpm1(:,1,2)+feq_jpm1(:,2,2))
          zigmc_meq (:)  = .125*(feq_jpm1(:,1,3)+feq_jpm1(:,2,3))
          zigm2_meq (:)  = .060*(feq_jpm1(:,1,4)+feq_jpm1(:,2,4))
          rim1_meq  (:)  = .060*(feq_jpm1(:,1,5)+feq_jpm1(:,2,5))
          rim2_tmp  (:)  = .060*(feq_jpm1(:,1,6)+feq_jpm1(:,2,6))
          do i=1,nmlonp1
            rim1_meq(i) = rim1_meq(i) - (zigm2_meq(i)-zigmc_meq(i))/
     |        zigm22_meq(i)*rim2_tmp(i)
          enddo
          rim1_meq(:) = rim1_meq(:)/dt1dts(j)
          rim1_meq(:) = rim1_meq(:)*2. ! fold eq on itself

          do i=mlon0,mlon1
            if (i==1) then               ! western most lon
              rhs(i,j) = 0.5/dlonm*(rim1(i+1,j)-rim1_meq(nmlon))
              rhs(i,j) = rhs(i,j)+1./dlatm*(tint1(j)*rim2(i,j)+
     |          tint1(j+1)*rim2_eqp1(i))
            elseif (i==nmlonp1) then       ! eastern most lon
              rhs(i,j) = 0.5/dlonm*(rim1_meq(1)-rim1(i-1,j))
              rhs(i,j) = rhs(i,j)+1./dlatm*(tint1(j)*rim2(i,j)+
     |          tint1(j+1)*rim2_eqp1(i))
            else                         ! body of i subdomain
! Note that rim1 halos were set before calling this subroutine.
              rhs(i,j) = 0.5/dlonm*(rim1(i+1,j)-rim1(i-1,j))
              rhs(i,j) = rhs(i,j)+1./dlatm*(tint1(j)*rim2(i,j)+
     |           tint1(j+1)*rim2_eqp1(i))
            endif
          enddo ! i=mlon0,mlon1
!         write(6,"('rhs: j==jlateq=',i4,' rhs(mlon0:mlon1,j)=',/,
!    |      (6e12.4))") j,rhs(mlon0:mlon1,j)
!
! North hemisphere (not npole and not equator):
! (allow south hemisphere to remain 0)
! (use rim1 instead of tint33)
!
        elseif (j > jlateq) then ! north hem only (excluding npole)
          do i=mlon0,mlon1
            rhs(i,j) = 1./(dlonm*tint1(j))*0.5*(rim1(i+1,j)-rim1(i-1,j))
            if (j == jlateq+1) then 
              rhs(i,j) = rhs(i,j)+1./(dlatm*tint1(j))*
     |         0.5*(rim2(i,j+1)*tint1(j+1)-rim2_meq(i)*tint1(j-1))
            else
              rhs(i,j) = rhs(i,j)+1./(dlatm*tint1(j))*
     |         0.5*(rim2(i,j+1)*tint1(j+1)-rim2(i,j-1)*tint1(j-1))
            endif
          enddo
!         write(6,"('rhs: j=nhem=',i4,' rhs(mlon0:mlon1,j)=',/,
!    |      (6e12.4))") j,rhs(mlon0:mlon1,j)
        endif
      enddo ! j=mlat0,mlat1
!
! scale (multiply by earth radius in meter  = R0*1.E-2)
![( d K_(m phi)^D / d phi /(cos(theta_m)?) +
! (d [ K_(m lam)^D * cos(lam_m)]/ d lam_m ) /cos ( lam_m) ] * R / (RCOS0S*DT0DTS)
! ~ J_(Mr)*r^2*cos(theta_m)/cos(theta_0)/DT0DTS
! theta_m = theta_s
!
      do j=mlat0,mlat1
        do i=mlon0,mlon1
          rhs(i,j) = rhs(i,j)*r0*1.e-2
        enddo
      enddo
!
! rhs here cannot be compared w/ odyn because nhem rhs here is defined
! in nhem indices, whereas in odyn nhem rhs is defined in shem indices.
! However, rhs odyn vs pdyn can be compared after sub gather_pdyn, where
! nhem rhs data is transferred to shem indices of rhs_glb.
!
!     call addfld('RHS_N',' ','S',rhs(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)

      end subroutine rhspde
!-----------------------------------------------------------------------
      subroutine gather_pdyn
!
! Gather needed global arrays to root task, so it can finish non-parallel
! part of dynamo (beginning after sub rhspde) as in original code 
!
      use mpi_module,only: mytid,mp_gather_pdyn
!
! Local:
! 7 fields to gather: zigm11,zigm22,zigmc,zigm2,rim1,rim2,rhs
! 
      integer,parameter :: nf = 7
      real :: fmsub(mlon0:mlon1,mlat0:mlat1,nf)
      real :: fmglb(nmlonp1,nmlat,nf)
      real :: rhs_nhem(nmlonp1,nmlat)
      integer :: i,j,jj
!
! These calls exclude halo points in zigm11, etc.
!
      fmsub(:,:,1) = zigm11  (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,2) = zigm22  (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,3) = zigmc   (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,4) = zigm2   (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,5) = rim1    (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,6) = rim2    (mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,7) = rhs     (mlon0:mlon1,mlat0:mlat1)

      call mp_gather_pdyn(fmsub,mlon0,mlon1,mlat0,mlat1,
     |  fmglb,nmlonp1,nmlat,nf)
!
! Gather 3d zpot in separate call (tell mp_gather_pdyn there
! are nmlevp1) fields.  When building under intel/ifort, must have 
! -heap-arrays flag set for this call to work:
!
      if (mytid==0) then
        zigm11_glb(:,:) = fmglb(:,:,1)
        zigm22_glb(:,:) = fmglb(:,:,2)
        zigmc_glb(:,:)  = fmglb(:,:,3)
        zigm2_glb(:,:)  = fmglb(:,:,4)
        rim_glb(:,:,1)  = fmglb(:,:,5)
        rim_glb(:,:,2)  = fmglb(:,:,6)
        rhs_nhem(:,:)   = fmglb(:,:,7)
!
! Transfer local global rhs_nhem (from sub rhspde) to rhs_glb, 
! so the latter has data in south hemisphere.
        rhs_glb= 0. ! init
        do j=1,nmlat
!
! Transfer north pole to equator:
          if (j == nmlat) then
            do i=1,nmlonp1
              rhs_glb(i,jlateq) = rhs_nhem(i,j)
            enddo
! Transfer equator to south pole:
          elseif (j == jlateq) then
            do i=1,nmlonp1
              rhs_glb(i,1) = rhs_nhem(i,j)
            enddo
! Transfer north hem to south hem:
          elseif (j > jlateq) then ! 50 -> 96
            jj = j-jlateq+1        !  2 -> 48 
            do i=1,nmlonp1
              rhs_glb(i,jj) = rhs_nhem(i,j)
            enddo
          endif
        enddo ! j=mlat0,mlat1
!
! This rhs can be compared w/ odyn rhs after rhspde:
!       call addfld('RHS',' ','S',rhs_glb,'mlon',1,nmlonp1,
!    |    'mlat',1,nmlat,0)
 
      endif ! mytid==0
      end subroutine gather_pdyn
!-----------------------------------------------------------------------
      subroutine productAX(cofum,phi,mlon0,mlon1,mlat0,mlat1,ax)
!      
! am 9/14/2010     
! calculates Ax and put it into ax      
!
      implicit none
!
! Args:
      integer, intent(in):: mlon0,mlon1,mlat0,mlat1
      real, intent(in):: cofum(mlon0:mlon1,mlat0:mlat1,9),
     |     phi(mlon0-1:mlon1+1,mlat0-1:mlat1+1) ! when called with phimsolv
!
! Local:
      real, intent(out)::ax(mlon0:mlon1,mlat0:mlat1)
      integer :: i,j,lonend
!      
! 9/21/11 btf: Set lonend to avoid exceeding upper bound of x(i+1,..)      
! (note eastern most tasks have mlon1==nmlonp1)
!
      lonend = mlon1
      if (mlon1 == nmlon+1) lonend = mlon1-1
      do j = mlat0,mlat1
	 do i = mlon0,lonend
	   ax(i,j) = (
     +     cofum(i,j,1)*phi(i+1,j)+
     +     cofum(i,j,2)*phi(i+1,j+1)+
     +     cofum(i,j,3)*phi(i,j+1)+
     +     cofum(i,j,4)*phi(i-1,j+1)+
     +     cofum(i,j,5)*phi(i-1,j)+
     +     cofum(i,j,6)*phi(i-1,j-1)+
     +     cofum(i,j,7)*phi(i,j-1)+
     +     cofum(i,j,8)*phi(i+1,j-1)+
     +     cofum(i,j,9)*phi(i,j))
	enddo
      enddo
!
! Save product to secondary histories:
!     call addfld('AX','AX product',' ',ax,'mlon',mlon0,mlon1,
!    |  'mlat',mlat0,mlat1,0)
      end subroutine productAX
!-----------------------------------------------------------------------
      subroutine highlat_poten
!
! Global PDE solution rim_glb(:,:,1) has been scattered to mag subdomains 
! in rim1, and halos set (this overwrites previous rim1 from fieldline 
! integrations).  Now add high latitude potential from empirical model 
! (heelis or weimer), defining phim2d on mag subdomains. After this pthreed 
! will expand phim2d to phim3d.
! 
! Input:  rim1(mag subdomains)   ! Solution from mudpack solver (nhem only)
!         pfrac(nmlonp1,nmlat0)  ! NH fraction of potential
!         phihm(nmlonp1,nmlat)   ! potential in magnetic
! Output: phim2d(mag subdomains) ! solution with phihm in both nhem and shem
! Both rim1 and phim2d are dimensioned (mlon00:mlon11,mlat00:mlat11)
!
! Both phihm and pfrac have been set by either heelis or weimer.
! phihm is on 2d global mag grid, pfrac is in north hemisphere only
!
      use input_module,only: iamie
      use amie_module,only: pcp_nh_amie,pcp_sh_amie
      use hist_module,only:   modeltime
      use init_module,only:   iyear
      use cons_module,only: dlatm,dlonm,pi_dyn,ylatm,rtd,crit
! Local:
      logical,parameter :: mod_heelis = .false.   ! true == modified
      integer :: i,j,jn,js
      real :: fac
!
! Add empirical model potential at high latitude:
      fac = 1.0
      if (mod_heelis) fac = 0.  ! modified heelis

      do j=mlat0,mlat1
        if (j > nmlath) cycle   ! south only (including equator)
        jn = nmlat-j+1
        js = nmlath-j+1
        do i=mlon0,mlon1
          phim2d(i,j) = rim1(i,j)+fac*(1.-pfrac(i,js))*(phihm(i,j)-
     |    phihm(i,jn))

        enddo
      enddo

      do j=mlat0,mlat1
        if (j <= nmlath) cycle ! north only (excluding equator)
        do i=mlon0,mlon1
          phim2d(i,j) = rim1(i,j)
        enddo
      enddo
!
! There is also a PHIM2D saved at beginning of pthreed below.
!     call addfld('PHIM2D',' ','VOLTS',
!    |  phim2d(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('RIM1',' ','VOLTS',
!    |  rim1(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)

      end subroutine highlat_poten
!-----------------------------------------------------------------------
      subroutine pthreed
!
! phim2d is now 2d electric potential on mag subdomains, 
! with high-latitude potential added from empirical model (see sub 
! highlat_poten), and mag halos set. Now expand phim2d in vertical,
! defining phim3d. Also calculate electric field ed13d, ed23d for
! later current calculations, and ephi3d, elam3d and emz3d for conversion
! to geographic grid (sub pefield), and calculation of ion drifts 
! by sub ionvel.
!
      use cons_module,only: ylatm,h0,re,pi_dyn,table,dlatm,dlonm,rcos0s,
     |  r0,dt1dts,dt0dts
      use mpi_module,only: mp_mageq_jpm3,mp_magpole_2d,mp_mag_jslot,
     |  mp_mag_halos,mp_magpoles
      use diags_module,only: mkdiag_ED12
!
! Local:
      real,parameter :: eps = 1.e-10, unitvm(nmlon)=1.
      integer,parameter :: mxneed=nmlat+2
      integer :: i,j,k,n,mlon00,mlon11,mlat00,mlat11
      real :: csth0,cosltm,sinltm,sym,pi,phims,phimn,siniq
      real,dimension(nmlonp1) :: thetam,pslot,qslot
      integer,dimension(nmlonp1) :: islot,jslot,ip1f,ip2f,ip3f
      real,dimension(mlon0-1:mlon1+1,mlat0-1:mlat1+1) :: ed1,ed2
      real,dimension(mlon0-1:mlon1+1,mlat0-1:mlat1+1) :: ephi,elam
      real :: fpole2d_jpm2(nmlonp1,4,4) ! global lons at S pole+1,2 and N pole-1,2
      real :: fpoles(nmlonp1,2,1) ! global lons at poles (1 field only)
      real :: fmsub(mlon0:mlon1,mlat0:mlat1,4)
      real :: fmsub1(mlon0-1:mlon1+1,mlat0-1:mlat1+1,5)
      real :: feq_jpm3(nmlonp1,-3:3,1) ! global lons at equator +/- 3
      integer :: jneed(mxneed) ! lats needed from other tasks for interp
      integer :: njneed,icount
      real,dimension(mlon0-1:mlon1+1,mxneed) ::
     |  phineed,   ! phim2d at needed latitudes
     |  ed1need,   ! ed1 at needed latitudes
     |  ed2need,   ! ed2 at needed latitudes
     |  ephineed,  ! ephi at needed latitudes
     |  elamneed   ! elam at needed latitudes
      real,dimension(mlon0-1:mlon1+1,mxneed,5) :: fmneed
      real :: phi0j0,phi1j0,phi0j1,phi1j1
      real :: ed1i0j0,ed1i1j0,ed1i0j1,ed1i1j1
      real :: ed2i0j0,ed2i1j0,ed2i0j1,ed2i1j1
      real :: ephi0j0,ephi1j0,ephi0j1,ephi1j1
      real :: elam0j0,elam1j0,elam0j1,elam1j1
!
! External:
      integer,external :: ixfind
      real,external :: sddot

      if (debug)
     |  write(6,"('Enter pthreed: phim2d=',2e12.4)") 
     |    minval(phim2d),maxval(phim2d)

      pi = pi_dyn
      mlon00=mlon0-1 ; mlon11=mlon1+1
      mlat00=mlat0-1 ; mlat11=mlat1+1

!     call addfld('PHIM2D',' ','VOLTS',
!    |  phim2d(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!
! Calculate ed1,ed2 components of electric field:
! phim2d has halos set, so when mlon0==1, i-1 should wrap
!   to value at i==nmlon, and when mlon1==nmlonp1, i+1 should 
!   wrap to value at i==2.
!
      ed1 = 0.
      ephi = 0.
      do j=mlat0,mlat1
        if (j==1.or.j==nmlat) cycle
        csth0 = cos(-pi/2.+(j-1)*dlatm)
        do i=mlon0,mlon1
          ed1(i,j) = -(phim2d(i+1,j)-phim2d(i-1,j))/(2.*dlonm*csth0)*
     |      rcos0s(j)/(r0*1.e-2)
          ephi(i,j) = ed1(i,j)*(r0*1.e-2)
        enddo ! i=mlon0,mlon1
      enddo ! j=mlat0,mlat1
!
! Southern hemisphere (excluding equator):
      ed2 = 0.
      elam = 0.
      do j=mlat0,mlat1
        if (j >= nmlath) cycle 
        do i=mlon0,mlon1
          ed2(i,j)= -(phim2d(i,j+1)-phim2d(i,j-1))/(2.*dlatm)*dt1dts(j)/
     |      (r0*1.e-2)
          elam(i,j)= -(phim2d(i,j+1)-phim2d(i,j-1))/(2.*dlatm)*dt0dts(j)
        enddo ! i=mlon0,mlon1
      enddo ! j=mlat0,mlat1
      if (debug) 
     |  write(6,"('pthreed after shem: ed1=',2e12.4,' ed2=',2e12.4,
     |    ' ephi=',2e12.4,' elam=',2e12.4)") minval(ed1),maxval(ed1),
     |    minval(ed2),maxval(ed2),minval(ephi),maxval(ephi),
     |    minval(elam),maxval(elam)
!
! Northern hemisphere (excluding equator):
      do j=mlat0,mlat1
        if (j <= nmlath) cycle
        do i=mlon0,mlon1
          ed2(i,j) = (phim2d(i,j+1)-phim2d(i,j-1))/(2.*dlatm)*dt1dts(j)/
     |      (r0*1.e-2)
          elam(i,j)= -(phim2d(i,j+1)-phim2d(i,j-1))/(2.*dlatm)*dt0dts(j)
        enddo ! i=mlon0,mlon1
      enddo ! j=mlat0,mlat1

!     call addfld('ED1_a',' ',' ',ed1(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ED2_a',' ',' ',ed2(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('EPHI_a',' ',' ',ephi(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ELAM_a',' ',' ',elam(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!
! Need ed1,2 at global longitudes at j==2 and j==nmlat-1:
! mp_magpole_2d: Return fpole_jpm2(nglblon,1->4,nf) as:
!   1: j = jspole+1 (spole+1)
!   2: j = jspole+2 (spole+2) (unused here)
!   3: j = jnpole-1 (npole-1)
!   4: j = jnpole-2 (npole-2) (unused here)
!
      fmsub(:,:,1) = ed1(mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,2) = ed2(mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,3) = ephi(mlon0:mlon1,mlat0:mlat1)
      fmsub(:,:,4) = elam(mlon0:mlon1,mlat0:mlat1)

! Returns fpole2d_jpm2(nmlonp1,4,4) ! global lons at S pole+1,2 and N pole-1,2
      call mp_magpole_2d(fmsub,mlon0,mlon1,mlat0,mlat1,nmlonp1,
     |  1,nmlat,fpole2d_jpm2,4)
!
! Poles: average over four surrounding points
      do i = 1,nmlonp1
        ip1f(i) = i + nmlon/4
        if (ip1f(i) > nmlonp1) ip1f(i) = ip1f(i) - nmlon
        ip2f(i) = i + nmlon/2
        if (ip2f(i) > nmlonp1) ip2f(i) = ip2f(i) - nmlon
        ip3f(i) = i + 3*nmlon/4
        if (ip3f(i) > nmlonp1) ip3f(i) = ip3f(i) - nmlon
      enddo ! i=1,nmlonp1
!
! S pole:
      do j=mlat0,mlat1
        if (j==1) then
          do i=mlon0,mlon1
            ed1(i,j)=.25*(fpole2d_jpm2(i,1,1)-fpole2d_jpm2(ip2f(i),1,1)+
     |              fpole2d_jpm2(ip1f(i),1,2)-fpole2d_jpm2(ip3f(i),1,2))
            ed2(i,j)=.25*(fpole2d_jpm2(i,1,2)-fpole2d_jpm2(ip2f(i),1,2)-
     |              fpole2d_jpm2(ip1f(i),1,1)+fpole2d_jpm2(ip3f(i),1,1))
           ephi(i,j)=.25*(fpole2d_jpm2(i,1,3)-fpole2d_jpm2(ip2f(i),1,3)+
     |              fpole2d_jpm2(ip1f(i),1,4)-fpole2d_jpm2(ip3f(i),1,4))
           elam(i,j)=.25*(fpole2d_jpm2(i,1,4)-fpole2d_jpm2(ip2f(i),1,4)-
     |              fpole2d_jpm2(ip1f(i),1,3)+fpole2d_jpm2(ip3f(i),1,3))
          enddo ! i=mlon0,mlon1
!
! N pole:
        elseif (j==nmlat) then
          do i=mlon0,mlon1
            ed1(i,j)=.25*(fpole2d_jpm2(i,3,1)-fpole2d_jpm2(ip2f(i),3,1)+
     |              fpole2d_jpm2(ip1f(i),3,2)-fpole2d_jpm2(ip3f(i),3,2))
            ed2(i,j)=.25*(fpole2d_jpm2(i,3,2)-fpole2d_jpm2(ip2f(i),3,2)-
     |              fpole2d_jpm2(ip1f(i),3,1)+fpole2d_jpm2(ip3f(i),3,1))
           ephi(i,j)=.25*(fpole2d_jpm2(i,3,3)-fpole2d_jpm2(ip2f(i),3,3)+
     |              fpole2d_jpm2(ip1f(i),3,4)-fpole2d_jpm2(ip3f(i),3,4))
           elam(i,j)=.25*(fpole2d_jpm2(i,3,4)-fpole2d_jpm2(ip2f(i),3,4)-
     |              fpole2d_jpm2(ip1f(i),3,3)+fpole2d_jpm2(ip3f(i),3,3))
          enddo ! i=mlon0,mlon1
        endif ! S or N pole
      enddo ! j=mlat0,mlat1

!     call addfld('ED1_b',' ',' ',ed1(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ED2_b',' ',' ',ed2(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('EPHI_b',' ',' ',ephi(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ELAM_b',' ',' ',elam(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!
! Equator: derivative of quadratic polynomial (3 given points)
! For equator and equator +/- 1 of ed2, we need equator and 
! equator +/- 3 of phim2d (note feq_jpm3(nmlonp1,-3:3,1)):
!
      fmsub(:,:,1) = phim2d(mlon0:mlon1,mlat0:mlat1)
      call mp_mageq_jpm3(fmsub(:,:,1),mlon0,mlon1,mlat0,mlat1,nmlonp1,
     |  feq_jpm3,1)
!
! Note timegcm does not use dt1dts here, but tiegcm does.
      do j=mlat0,mlat1
        if (j==nmlath-1) then     ! equator-1
          do i=mlon0,mlon1
            ed2(i,j) = (4.*feq_jpm3(i,-2,1)-feq_jpm3(i,-3,1)-
     |       3.*feq_jpm3(i,-1,1))/(2.*dlatm)*dt1dts(nmlath-1)/(r0*1.e-2)
          enddo
        elseif (j==nmlath) then   ! equator
          do i=mlon0,mlon1
            ed2(i,j) = (4.*feq_jpm3(i,1,1)-feq_jpm3(i,2,1)-
     |        3.*feq_jpm3(i,0,1))/(2.*dlatm)*dt1dts(nmlath)/(r0*1.e-2)
           elam(i,j) = (4.*feq_jpm3(i,1,1)-feq_jpm3(i,2,1)-
     |        3.*feq_jpm3(i,0,1))/(2.*dlatm)
          enddo
        elseif (j==nmlath+1) then ! equator+1
          do i=mlon0,mlon1
            ed2(i,j) = (4.*feq_jpm3(i,2,1)-feq_jpm3(i,3,1)-
     |        3.*feq_jpm3(i,1,1))/(2.*dlatm)*dt1dts(nmlath+1)/(r0*1.e-2)
          enddo
        endif ! equator +/- 1
      enddo ! j=mlat0,mlat1

!     call addfld('ED1_c',' ',' ',ed1(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ED2_c',' ',' ',ed2(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('EPHI_c',' ',' ',ephi(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ELAM_c',' ',' ',elam(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!
! Set halos for 3d calculations:
      fmsub1(:,:,1) = ed1
      fmsub1(:,:,2) = ed2
      fmsub1(:,:,3) = ephi
      fmsub1(:,:,4) = elam

      if (debug)
     |  write(6,"('pthreed call mp_mag_halos: fmsub1(ed1)=',2e12.4,
     |  ' fmsub1(ed2)=',2e12.4,' fmsub1(ephi)=',2e12.4,' fmsub1(elam)=',
     |    2e12.4)") minval(fmsub1(:,:,1)),maxval(fmsub1(:,:,1)),
     |             minval(fmsub1(:,:,2)),maxval(fmsub1(:,:,2)),
     |             minval(fmsub1(:,:,3)),maxval(fmsub1(:,:,3)),
     |             minval(fmsub1(:,:,4)),maxval(fmsub1(:,:,4))

      call mp_mag_halos(fmsub1,mlon0,mlon1,mlat0,mlat1,1,4)

      ed1 = fmsub1(:,:,1)
      ed2 = fmsub1(:,:,2)
      ephi = fmsub1(:,:,3)
      elam = fmsub1(:,:,4)

!     call addfld('ED1',' ',' ',ed1(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ED2',' ',' ',ed2(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('EPHI',' ',' ',ephi(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
!     call addfld('ELAM',' ',' ',elam(mlon0:mlon1,mlat0:mlat1),
!    |  'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)

!
! Determine latitudes needed for interpolation that fall
! outside a task's latitudinal subdomain: 
      njneed = 0    ! number of unique latitudes needed
      jneed(:) = -1 ! j-indices of needed latitudes
      do k=mlev0,mlev1
        do j=mlat0,mlat1
          if (j==1.or.j==nmlat) cycle ! exclude poles
          sym = 1.
          if (j < nmlath) sym = -1.
          cosltm = cos(ylatm(j))
          do i=mlon0,mlon1
            if (i==nmlonp1) cycle

            thetam(i)=(re+zpotm3d(i,j,mlev0))/(re+zpotm3d(i,j,k))

!           write(6,"('pthreed: i=',i4,' j=',i4,' k=',i4,
!    |        ' zpotm3d(i,j,mlev0)=',e12.4,' zpotm3d(i,j,k)=',e12.4,
!    |        ' thetam(i)=',e12.4)") i,j,k,zpotm3d(i,j,mlev0),
!    |        zpotm3d(i,j,k),thetam(i)

            thetam(i) = acos(sqrt(thetam(i))*cosltm*(1.-eps))

!           write(6,"('pthreed: i=',i4,' j=',i4,' k=',i4,' cosltm=',
!    |        e12.4,' thetam(i)=',e12.4)") i,j,k,cosltm,thetam(i)

            pslot(i) = thetam(i)*180./pi+1.
            islot(i) = pslot(i)
            pslot(i) = pslot(i)-float(islot(i))

            thetam(i) = ((1.-pslot(i))*table(islot(i),2)+pslot(i)*
     |        table(islot(i)+1,2))*sym ! thetam negative for south hem

            islot(i) = i
            pslot(i) = 0.
            qslot(i) = (thetam(i)+pi/2.)/dlatm+1.
            jslot(i) = qslot(i)
            qslot(i) = qslot(i)-float(jslot(i))
!
! Save j index if outside subdomain w/ halos:
            if ((jslot(i) < mlat00 .or. jslot(i) > mlat11).and.
     |         .not.any(jslot(i)==jneed)) then
              njneed = njneed+1
              if (njneed > mxneed) call shutdown('njneed')
              jneed(njneed) = jslot(i) 
            endif ! jslot is outside subdomain
!
! Save j+1 index if outside subdomain:
            if ((jslot(i)+1 < mlat00 .or. jslot(i)+1 > mlat11).and.
     |           .not.any(jslot(i)+1==jneed)) then
              njneed = njneed+1
              if (njneed > mxneed) call shutdown('njneed')
              jneed(njneed) = jslot(i)+1 
            endif ! jslot(i)+1 is outside subdomain
          enddo ! i=mlon0,mlon1
        enddo ! j=mlat0,mlat1
      enddo ! k=mlev0,mlev1
!
! Get phim2 at needed latitudes (note inclusion of phim2d halos).
!     real,intent(in) :: fin(mlon00:mlon11,mlat00:mlat11,nf) ! data at current subdomain
!     real,intent(out) :: fout(mlon00:mlon11,mxneed,nf) ! returned data at needed lats
!
      fmsub1(:,:,1) = phim2d
      fmsub1(:,:,2) = ed1
      fmsub1(:,:,3) = ed2
      fmsub1(:,:,4) = ephi
      fmsub1(:,:,5) = elam
      if (debug) write(6,"('pthreed call mp_mag_jslot')")
      call mp_mag_jslot(fmsub1,mlon00,mlon11,mlat00,mlat11,
     |  fmneed,jneed,mxneed,5)
      if (debug) write(6,"('pthreed after mp_mag_jslot')")
      phineed = fmneed(:,:,1)
      ed1need = fmneed(:,:,2)
      ed2need = fmneed(:,:,3)
      ephineed= fmneed(:,:,4)
      elamneed= fmneed(:,:,5)
!
      ephi3d = 0.
      elam3d = 0.
      emz3d = 0.

      do k=mlev0,mlev1
        do j=mlat0,mlat1
          if (j==1.or.j==nmlat) cycle ! exclude poles
          sym = 1.
          if (j < nmlath) sym = -1.
          cosltm = cos(ylatm(j))
          sinltm = sin(ylatm(j))
          do i=mlon0,mlon1
            if (i==nmlonp1) cycle
            siniq = 2.*sinltm/sqrt(4.-3.*cosltm**2)

            thetam(i)=(re+zpotm3d(i,j,mlev0))/(re+zpotm3d(i,j,k))
            thetam(i) = acos(sqrt(thetam(i))*cosltm*(1.-eps))

            pslot(i) = thetam(i)*180./pi+1.
            islot(i) = pslot(i)
            pslot(i) = pslot(i)-float(islot(i))

            thetam(i) = ((1.-pslot(i))*table(islot(i),2)+pslot(i)*
     |        table(islot(i)+1,2))*sym ! thetam negative for south hem

            islot(i) = i
            pslot(i) = 0.
            qslot(i) = (thetam(i)+pi/2.)/dlatm+1.
            jslot(i) = qslot(i)
            qslot(i) = qslot(i)-float(jslot(i))
!
! Check for jslot in subdomain:
            if (jslot(i) >= mlat00.and.jslot(i) <= mlat11) then ! within subdomain
              phi0j0 = phim2d(islot(i)  ,jslot(i))
              phi1j0 = phim2d(islot(i)+1,jslot(i))
              ed1i0j0 =   ed1(islot(i)  ,jslot(i))
              ed1i1j0 =   ed1(islot(i)+1,jslot(i))
              ed2i0j0 =   ed2(islot(i)  ,jslot(i))
              ed2i1j0 =   ed2(islot(i)+1,jslot(i))
              ephi0j0 =  ephi(islot(i)  ,jslot(i))
              ephi1j0 =  ephi(islot(i)+1,jslot(i))
              elam0j0 =  elam(islot(i)  ,jslot(i))
              elam1j0 =  elam(islot(i)+1,jslot(i))
            else                                           ! jslot outside subdomain
              n = ixfind(jneed,mxneed,jslot(i),icount)
              if (n==0) then
                write(6,"('>>> pthreed: i=',i4,' j=',i4,' k=',i4,
     |            ' could not find jslot ',i4,' in jneed=',/,
     |            (10i4))") i,j,k,jslot(i),jneed
                call shutdown('jslot(i) not in jneed')
              endif
              phi0j0  = phineed(islot(i)  ,n)
              phi1j0  = phineed(islot(i)+1,n)
              ed1i0j0 = ed1need(islot(i)  ,n)
              ed1i1j0 = ed1need(islot(i)+1,n)
              ed2i0j0 = ed2need(islot(i)  ,n)
              ed2i1j0 = ed2need(islot(i)+1,n)
              ephi0j0 =ephineed(islot(i)  ,n)
              ephi1j0 =ephineed(islot(i)+1,n)
              elam0j0 =elamneed(islot(i)  ,n)
              elam1j0 =elamneed(islot(i)+1,n)
            endif
!
! Check for jslot+1 in subdomain:
            if (jslot(i)+1 >= mlat00.and.jslot(i)+1 <= mlat11) then ! within subdomain
              phi0j1 = phim2d(islot(i)  ,jslot(i)+1)
              phi1j1 = phim2d(islot(i)+1,jslot(i)+1)
              ed1i0j1 =   ed1(islot(i)  ,jslot(i)+1)
              ed1i1j1 =   ed1(islot(i)+1,jslot(i)+1)
              ed2i0j1 =   ed2(islot(i)  ,jslot(i)+1)
              ed2i1j1 =   ed2(islot(i)+1,jslot(i)+1)
              ephi0j1 =  ephi(islot(i)  ,jslot(i)+1)
              ephi1j1 =  ephi(islot(i)+1,jslot(i)+1)
              elam0j1 =  elam(islot(i)  ,jslot(i)+1)
              elam1j1 =  elam(islot(i)+1,jslot(i)+1)
            else                                           ! jslot+1 outside subdomain
              n = ixfind(jneed,mxneed,jslot(i)+1,icount)
              if (n==0) then
                write(6,"('>>> pthreed: i=',i4,' j=',i4,' k=',i4,
     |            ' could not find jslot+1 ',i4,' in jneed=',/,
     |            (10i4))") i,j,k,jslot(i)+1,jneed
                call shutdown('jslot(i)+1 not in jneed')
              endif
              phi0j1  = phineed(islot(i)  ,n)
              phi1j1  = phineed(islot(i)+1,n)
              ed1i0j1 = ed1need(islot(i)  ,n)
              ed1i1j1 = ed1need(islot(i)+1,n)
              ed2i0j1 = ed2need(islot(i)  ,n)
              ed2i1j1 = ed2need(islot(i)+1,n)
              ephi0j1 =ephineed(islot(i)  ,n)
              ephi1j1 =ephineed(islot(i)+1,n)
              elam0j1 =elamneed(islot(i)  ,n)
              elam1j1 =elamneed(islot(i)+1,n)
            endif
!
! Do the interpolation:
            phim3d(i,j,k) = (1.-qslot(i))*((1.-pslot(i))*
     |        phi0j0 + pslot(i) * phi1j0) + qslot(i)*((1.-pslot(i))*
     |        phi0j1 + pslot(i) * phi1j1)

            ed13d(i,j,k) = (1.-qslot(i))*((1.-pslot(i))*
     |        ed1i0j0 + pslot(i) * ed1i1j0) + qslot(i)*((1.-pslot(i))*
     |        ed1i0j1 + pslot(i) * ed1i1j1)
            ed23d(i,j,k) = (1.-qslot(i))*((1.-pslot(i))*
     |        ed2i0j0 + pslot(i) * ed2i1j0) + qslot(i)*((1.-pslot(i))*
     |        ed2i0j1 + pslot(i) * ed2i1j1)
!
            ephi3d(i,j,k) = (1.-qslot(i))*((1.-pslot(i))*
     |        ephi0j0 + pslot(i) * ephi1j0) + qslot(i)*((1.-pslot(i))*
     |        ephi0j1 + pslot(i) * ephi1j1)
            elam3d(i,j,k) = (1.-qslot(i))*((1.-pslot(i))*
     |        elam0j0 + pslot(i) * elam1j0) + qslot(i)*((1.-pslot(i))*
     |        elam0j1 + pslot(i) * elam1j1)

          enddo ! i=mlon0,mlon1
        enddo ! j=mlat0,mlat1
      enddo ! k=mlev0,mlev1
!
! Mag poles for phim:
! mp_magpoles returns global longitudes at S,N poles in fpoles(nglblon,2,nf)
! 
      if (debug) write(6,"('pthreed call mp_magpoles')")

      call mp_magpoles(phim2d(mlon0:mlon1,mlat0:mlat1),
     |  mlon0,mlon1,mlat0,mlat1,nmlonp1,1,nmlat,fpoles,1)

      if (debug) write(6,"('pthreed after mp_magpoles')")

      phims=sddot(nmlon,unitvm,fpoles(1:nmlon,1,1))/float(nmlon)
      phimn=sddot(nmlon,unitvm,fpoles(1:nmlon,2,1))/float(nmlon)

      do k=mlev0,mlev1
        do j=mlat0,mlat1
          if (j==1) then
            do i=mlon0,mlon1
              phim3d(i,j,k) = phims
              ed13d(i,j,k) = ed1(i,j)
              ed23d(i,j,k) = ed2(i,j)
              ephi3d(i,j,k) = ephi(i,j)
              elam3d(i,j,k) = ed2(i,j)*(r0*1.e-2)
            enddo
          elseif (j==nmlat) then
            do i=mlon0,mlon1
              phim3d(i,j,k) = phimn
              ed13d(i,j,k) = ed1(i,j)
              ed23d(i,j,k) = ed2(i,j)
              ephi3d(i,j,k) = ephi(i,j)
              elam3d(i,j,k) = -ed2(i,j)*(r0*1.e-2)
            enddo
          endif ! poles
        enddo ! j=mlat0,mlat1
      enddo ! k=mlev0,mlev1
!
! Upward electric field:
      do k=mlev0+1,mlev1-1
        do j=mlat0,mlat1
          do i=mlon0,mlon1
            emz3d(i,j,k) = -(phim3d(i,j,k+1)-phim3d(i,j,k-1))
          enddo
        enddo
      enddo
      do j=mlat0,mlat1
        do i=mlon0,mlon1
          emz3d(i,j,mlev1) = emz3d(i,j,mlev1-1)
          emz3d(i,j,mlev0) = emz3d(i,j,mlev0+1)
        enddo
      enddo

      if (debug) write(6,"('pthreed before mp_mag_periodic_f2d')")
      do k=mlev0,mlev1
        call mp_mag_periodic_f2d(phim3d(:,:,k),mlon0,mlon1,
     |    mlat0,mlat1,1)
        call mp_mag_periodic_f2d(ed13d(:,:,k),mlon0,mlon1,
     |    mlat0,mlat1,1)
        call mp_mag_periodic_f2d(ed23d(:,:,k),mlon0,mlon1,
     |    mlat0,mlat1,1)
        call mp_mag_periodic_f2d(ephi3d(:,:,k),mlon0,mlon1,
     |    mlat0,mlat1,1)
        call mp_mag_periodic_f2d(ephi3d(:,:,k),mlon0,mlon1,
     |    mlat0,mlat1,1)
        call mp_mag_periodic_f2d(elam3d(:,:,k),mlon0,mlon1,
     |    mlat0,mlat1,1)
        call mp_mag_periodic_f2d(emz3d(:,:,k),mlon0,mlon1,
     |    mlat0,mlat1,1)
      enddo
      if (debug) write(6,"('pthreed after mp_mag_periodic_f2d')")
!
! Cannot use mlev0 (-2) in addfld call, so plot from level 1:
!
      do j=mlat0,mlat1
        call addfld('PHIM3D',' ',' ',phim3d(:,j,mlev0:mlev1-1),
     |    'mlon',mlon0,mlon1,'mlev',1,mlevtop,j)
        call addfld('EMX',' ',' ',ephi3d(:,j,mlev0:mlev1-1),
     |    'mlon',mlon0,mlon1,'mlev',1,mlevtop,j)
        call addfld('EMY',' ',' ',elam3d(:,j,mlev0:mlev1-1),
     |    'mlon',mlon0,mlon1,'mlev',1,mlevtop,j)
        call addfld('EMZ',' ',' ',emz3d(:,j,mlev0:mlev1-1),
     |    'mlon',mlon0,mlon1,'mlev',1,mlevtop,j)
        call mkdiag_ED12('ED1',ed13d(:,j,mlev0:mlev1-1),
     |    mlon0,mlon1,1,mlevtop,j)
        call mkdiag_ED12('ED2',ed23d(:,j,mlev0:mlev1-1),
     |    mlon0,mlon1,1,mlevtop,j)
      enddo
      if (debug) write(6,"('pthreed returning')")
      end subroutine pthreed
!-----------------------------------------------------------------------
      subroutine pefield
!
! Parallel version of efield. Get electric field on geographic subdomains.
! (ex,ey,ez in this module).  This is called from advance, early in the 
! timestep. If first timestep, calculate electric field in mag coords from 
! phim3d, otherwise use electric field on mag subdomains from pthreed of 
! the previous step.
!
      use cons_module,only: dt0dts,dlatm,dlonm,pi,rcos0s,table
      use mpi_module,only: mp_magpole_3d,mp_mag_halos,mp_magpoles,
     |  mytid
      use magfield_module,only: alatm,alonm  ! (nlonp1,0:nlatp1)
      use my_esmf,only: mag_ephi3d,mag_elam3d,mag_emz3d,
     |                  geo_ephi3d,geo_elam3d,geo_emz3d
      use diags_module,only: mkdiag_EXYZ
      implicit none
!
! Local:
      integer :: i,ii,j,k,lonbeg,lonend
      real :: fpole3d_jpm2(nmlonp1,4,mlev0:mlev1,1) ! global lons at S pole+1,2 and N pole-1,2
      real :: fpoles(nmlonp1,2,mlev0:mlev1) ! global lons at poles
      real :: phi3d(mlon0-1:mlon1+1,mlat0-1:mlat1+1,mlev0:mlev1)
      real :: csth0
      integer,dimension(nlonp4) :: islot,jslot
      real,dimension(nlonp4) :: pslot,qslot,thetam
      real,dimension(lon0:lon1,lat0:lat1,mlev0:mlev1) :: 
     |  exgeo,eygeo,ezgeo

      if (debug) write(6,"('Enter pefield: istep=',i3)") istep

!     do j=mlat0,mlat1
!       call addfld('PHIM3D_EFLD',' ','VOLTS',
!    |    phim3d(mlon0:mlon1,j,1:mlev1),
!    |    'mlon',mlon0,mlon1,'imlev',1,mlev1,j)
!     enddo
!
! Note if istep > 1, then 3d electric field ephi3d,elam3d,emz3d was
! calculated by pthreed in previous step.
      if (istep==1) then
        if (debug)
     |  write(6,"('pefield (step=1): phim3d(:,:,:) min,max=',
     |    2e12.4)") minval(phim3d(:,:,:)),maxval(phim3d(:,:,:))
!
! Copy phim3d to local phi3d, and set halo points:
! (phim3d does not have halo points)
!
        phi3d = 0.
        do k=mlev0,mlev1
          do j=mlat0,mlat1
            do i=mlon0,mlon1
              phi3d(i,j,k) = phim3d(i,j,k)
            enddo
          enddo
        enddo
        do k=mlev0,mlev1
          if (debug) write(6,"('pefield (step=1) k=',i4,
     |      ' call mp_mag_halos for phi3d')") k

          call mp_mag_halos(phi3d(:,:,k),mlon0,mlon1,mlat0,mlat1,1,1)

          if (debug) write(6,"('pefield (step=1) k=',i4,
     |      ' after mp_mag_halos for phi3d')") k
        enddo
!
! Return fpole3d_jpm2(nglblon,1->4,mlev0:mlev1,nf) as:
!   1: j = jspole+1 (spole+1)
!   2: j = jspole+2 (spole+2) not used here
!   3: j = jnpole-1 (npole-1)
!   4: j = jnpole-2 (npole-2) not used here
!
        if (debug) write(6,"('pefield (step=1) call mp_magpole_3d',
     |    ' for phi3d')")
        call mp_magpole_3d(phim3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),
     |   mlon0,mlon1,mlat0,mlat1,nmlevp1,nmlonp1,1,nmlat,fpole3d_jpm2,1)
        if (debug) write(6,"('pefield (step=1) after mp_magpole_3d',
     |    ' for phi3d')")
!
! Set j=0 and j=nmlat+1 of local phi3d. This overwrites the far 
! north and south halo points set by mp_mag_halos above. 
!
        do j=mlat0,mlat1
          if (j==1) then
            do i=mlon0,mlon1
              ii = 1+mod(i-1+nmlon/2,nmlon) ! over the south pole
              phi3d(i,j-1,:) = fpole3d_jpm2(ii,1,:,1)
            enddo
          elseif (j==nmlat) then
            do i=mlon0,mlon1
              ii = 1+mod(i-1+nmlon/2,nmlon) ! over the north pole
              phi3d(i,j+1,:) = fpole3d_jpm2(ii,3,:,1)
            enddo
          endif ! poles or not
        enddo ! j=mlat0,mlat1
!
! Meridional component of electric field:
        do k=mlev0,mlev1
          do j=mlat0,mlat1
            do i=mlon0,mlon1
              elam3d(i,j,k) = -(phi3d(i,j+1,k)-phi3d(i,j-1,k))/
     |          (2.*dlatm)*dt0dts(j)
            enddo
          enddo
        enddo
!
! Zonal component of electric field:
        do j=mlat0,mlat1
          if (j==1.or.j==nmlat) cycle
          csth0 = cos(-pi/2.+float(j-1)*dlatm)
          do k=mlev0,mlev1
            do i=mlon0,mlon1
              ephi3d(i,j,k) = -(phi3d(i+1,j,k)-phi3d(i-1,j,k))/
     |          (2.*dlonm*csth0)*rcos0s(j) 
            enddo
          enddo
        enddo
!
! Polar values for ephi3d (need global lons at poles of elam3d):
!
        if (debug) write(6,"('pefield call mp_magpoles for elam3d')")
        call mp_magpoles(elam3d(:,:,mlev0:mlev1),mlon0,mlon1,
     |    mlat0,mlat1,nmlonp1,1,nmlat,fpoles,nmlevp1)
        if (debug) write(6,"('pefield after mp_magpoles for elam3d')")

        do k=mlev0,mlev1
          do j=mlat0,mlat1
            if (j==1) then         ! south pole
              do i=mlon0,mlon1
                ii = 1+mod(i-1+(nmlon/4),nmlon) ! over the south pole
                ephi3d(i,j,k) = fpoles(ii,1,k)
              enddo
            elseif (j==nmlat) then ! north pole
              do i=mlon0,mlon1
                ii = 1+mod(i-1+((3*nmlon)/4),nmlon) ! over the north pole
                ephi3d(i,j,k) = fpoles(ii,2,k)
              enddo
            endif ! poles or not
          enddo ! j=mlat0,mlat1
        enddo ! k=mlev0,mlev1
!
! emz = d(phi)/dz
        do k=mlev0+1,mlev1-1
          do j=mlat0,mlat1
            do i=mlon0,mlon1
              emz3d(i,j,k) = -(phim3d(i,j,k+1)-phi3d(i,j,k-1))
            enddo
          enddo
        enddo ! k=mlev0+1,mlev1-1
      endif ! istep == 1
!
! If istep > 1, then ephi3d, elam3d and emz3d were set by pthreed
! in the previous step.
!
!     do j=mlat0,mlat1
!     call addfld('EMX_EFLD','EMX in pefield',' ',ephi3d(:,j,1:mlev1),
!    |    'mlon',mlon0,mlon1,'mlev',1,mlev1,j)
!     call addfld('EMY_EFLD','EMY in pefield',' ',elam3d(:,j,1:mlev1),
!    |    'mlon',mlon0,mlon1,'mlev',1,mlev1,j)
!     call addfld('EMZ_EFLD','EMZ in pefield',' ',emz3d(:,j,1:mlev1),
!    |    'mlon',mlon0,mlon1,'mlev',1,mlev1,j)
!     enddo ! j=1,nmlat
!
! Use ESMF to regrid the electric field to the geographic grid:      
!     allocate(ephi3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
!     real,dimension(lon0:lon1,lat0:lat1,mlev0:mlev1) :: exgeo,...

      if (debug) write(6,"('pefield call mag2geo_3d')")
      call mag2geo_3d(ephi3d,exgeo,mag_ephi3d,geo_ephi3d,'EPHI3D  ')
      call mag2geo_3d(elam3d,eygeo,mag_elam3d,geo_elam3d,'ELAM3D  ')
      call mag2geo_3d(emz3d, ezgeo,mag_emz3d ,geo_emz3d ,'EMZ3D   ')
      if (debug) write(6,"('pefield after mag2geo_3d')")
!
! Define ex,ey,ez on geographic subdomains for ionvel:
      do k=1,nlevp1
        do j=lat0,lat1
          do i=lon0,lon1
            ex(k,i,j) = exgeo(i,j,k)
            ey(k,i,j) = eygeo(i,j,k)
            ez(k,i,j) = ezgeo(i,j,k)
          enddo
        enddo

        if (debug)
     |  write(6,"('pefield: k=',i4,' ex=',2e12.4,' ey=',2e12.4,
     |    ' ez=',2e12.4)") k,
     |   minval(ex(k,lon0:lon1,lat0:lat1)),
     |   maxval(ex(k,lon0:lon1,lat0:lat1)),
     |   minval(ey(k,lon0:lon1,lat0:lat1)),
     |   maxval(ey(k,lon0:lon1,lat0:lat1)),
     |   minval(ez(k,lon0:lon1,lat0:lat1)),
     |   maxval(ez(k,lon0:lon1,lat0:lat1))

      enddo ! k=1,nlevp1

      do j=lat0,lat1
        call mkdiag_EXYZ('EX',ex(1:nlevp1,lon0:lon1,j),1,nlevp1,
     |    lon0,lon1,j)
        call mkdiag_EXYZ('EY',ey(1:nlevp1,lon0:lon1,j),1,nlevp1,
     |    lon0,lon1,j)
        call mkdiag_EXYZ('EZ',ez(1:nlevp1,lon0:lon1,j),1,nlevp1,
     |    lon0,lon1,j)
      enddo

      if (debug) write(6,"('pefield returning')")
      end subroutine pefield
!-----------------------------------------------------------------------
      subroutine mag2geo_3d(fmag,fgeo,ESMF_mag,ESMF_geo,fname)

! From pefield:
!     call mag2geo_3d(ephi3d,exgeo,mag_ephi3d,geo_ephi3d,'EPHI3D  ')
!     allocate(ephi3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
!     real,dimension(lon0:lon1,lat0:lat1,mlev0:mlev1) :: exgeo,...
!
! Convert field on geomagnetic grid fmag to geographic grid in fgeo.
!
      use my_esmf,only: esmf_set3d_mag,esmf_regrid,esmf_get_3dfield,
     |  ESMF_Field
!
! Args:
      character(len=*) :: fname
      type(ESMF_Field),intent(inout) :: ESMF_mag, ESMF_geo
      real,intent(in) :: fmag(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)
      real,intent(out) :: fgeo(lon0:lon1,lat0:lat1,mlev0:mlev1)
!
! Local:
      integer :: j
      character(len=8) :: fnames(1)
      type(ESMF_Field) :: magfields(1)
      real,pointer,dimension(:,:,:) :: fptr

      fgeo = 0.
      fnames(1) = fname
      magfields(1) = ESMF_mag
!
! Put fmag into ESMF mag field on mag source grid:
      call esmf_set3d_mag(magfields,fnames,fmag,1,mlev0,mlev1,
     |  mlon0,mlon1,mlat0,mlat1)
!
! Regrid to geographic destination grid, defining ESMF_geo:
      call esmf_regrid(ESMF_mag,ESMF_geo,'mag2geo',3)
!
! Put regridded geo field into pointer:
      call esmf_get_3dfield(ESMF_geo,fptr,fname)
!
! Transfer from pointer to output arg:
      do j=lat0,lat1
        fgeo(:,j,:) = fptr(:,j,:)
      enddo

      end subroutine mag2geo_3d
!-----------------------------------------------------------------------
!
#if defined(INTERCOMM) || defined(CISMAH)
      subroutine mag2geo_2d(fmag,fgeo,ESMF_mag,ESMF_geo,fname)
!
! Convert field on geomagnetic grid fmag to geographic grid in fgeo.
!
      use my_esmf,only: esmf_set2d_mag,esmf_regrid,esmf_get_2dfield,
     |  ESMF_Field
!
! Args:
      character(len=*) :: fname
      type(ESMF_Field),intent(inout) :: ESMF_mag, ESMF_geo
      real,intent(in) :: fmag(mlon0:mlon1,mlat0:mlat1)
      real,intent(out) :: fgeo(lon0:lon1,lat0:lat1)
!
! Local:
      integer :: j
      character(len=8) :: fnames(1)
      type(ESMF_Field) :: magfields(1)
      real,pointer,dimension(:,:) :: fptr

      fgeo = 0.
      fnames(1) = fname
      magfields(1) = ESMF_mag
!
! Put fmag into ESMF mag field on mag source grid:
      call esmf_set2d_mag(magfields,fnames,fmag,1,
     |  mlon0,mlon1,mlat0,mlat1)
!
! Regrid to geographic destination grid, defining ESMF_geo:
      call esmf_regrid(ESMF_mag,ESMF_geo,'mag2geo',2)
!
! Put regridded geo field into pointer:
      call esmf_get_2dfield(ESMF_geo,fptr,fname)
!
! Transfer from pointer to output arg:
      do j=lat0,lat1
        fgeo(:,j) = fptr(:,j)
      enddo

      end subroutine mag2geo_2d
#endif
!-----------------------------------------------------------------------
      subroutine geo2mag_3d(fgeo,fmag,ESMF_geo,ESMF_mag,fname)
!
! Convert field on geographic grid fgeo to geomagnetic grid in fmag.
!
      use my_esmf,only: esmf_set3d_geo,esmf_regrid,esmf_get_3dfield,
     |  ESMF_Field
!
! Args:
      character(len=*) :: fname
      type(ESMF_Field),intent(inout) :: ESMF_geo, ESMF_mag
      real,intent(in) :: fgeo(nlevp1,lon0:lon1,lat0:lat1)
      real,intent(out) :: fmag(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)
!
! Local:
      integer :: j,k
      character(len=8) :: fnames(1)
      type(ESMF_Field) :: geofields(1)
      real,pointer,dimension(:,:,:) :: fptr
      real :: fgeotmp(mlev0:mlev1,lon0:lon1,lat0:lat1)
      integer :: lbnd(3),ubnd(3)

!     do j=lat0,lat1
!       write(6,"('geo2mag_3d: j=',i4,' fgeo(:,lon0:lon1,j)=',
!    |    2e12.4)") j,minval(fgeo(:,lon0:lon1,j)),
!    |                maxval(fgeo(:,lon0:lon1,j))
!     enddo

      fmag = 0.
      fnames(1) = fname
      geofields(1) = ESMF_geo
!
! Put fgeo into ESMF geo field on geo source grid:
! subroutine esmf_set3d_geo(fields,fnames,f,nf,ilev0,ilev1,
!   ilon0,ilon1,ilat0,ilat1)
!
! This leaves bottom k-levels (-2 to 0 or -6 to 0 depending on resolution)
! of fgeotmp zeroed out:
!
      fgeotmp = 0.
      fgeotmp(1:nlevp1,:,:) = fgeo(:,:,:)
! 
! Note esmf_set3d_geo will change dimension order from input 
!   fgeotmp(k,i,j) to output ESMF_geo(i,j,k).
!
      call esmf_set3d_geo(geofields,fnames,fgeotmp,1,mlev0,mlev1,
     |  lon0,lon1,lat0,lat1)
!
! Get ESMF_geo and print bounds and min,max of levels 1,nlevp1:
      call esmf_get_3dfield(ESMF_geo,fptr,fname)

!     write(6,"('geo2mag_3d after esmf_get_3dfield: lbnd(fptr)=',3i4,
!    |  ' ubnd(fptr)=',3i4)") lbound(fptr),ubound(fptr)

!     do j=lat0,lat1
!       write(6,"('geo2mag_3d: j=',i4,' fptr to ',
!    |    'ESMF_geo(:,j,1:nlevp1)=',2e12.4)") 
!    |    j,minval(fptr(:,j,1:nlevp1)),maxval(fptr(:,j,1:nlevp1))
!     enddo
!
! Regrid to magnetic destination grid, defining ESMF_mag:
      call esmf_regrid(ESMF_geo,ESMF_mag,'geo2mag',3)
!
! Put regridded mag field into pointer:
      call esmf_get_3dfield(ESMF_mag,fptr,fname)
!
! Transfer from pointer to output arg:
      do j=mlat0,mlat1
        fmag(:,j,:) = fptr(:,j,:)

!       write(6,"('geo2mag_3d: j=',i4,' fmag(:,j,1:mlev1)=',2e12.4)")
!    |    j,minval(fmag(:,j,1:mlev1)),maxval(fmag(:,j,1:mlev1))

      enddo

      end subroutine geo2mag_3d
!-----------------------------------------------------------------------
      subroutine prepare_phig3d(wrprim,wrsech)
!
! If its time to write a history, convert phim3d to geographic grid in
! phig3d, and transfer to appropriate f4d array, where it will be
! gathered to the root task before writing to the history.
!
      use fields_module,only: itc,poten     ! pointer to f4d(i_poten)
      use my_esmf,only: mag_phi3d,geo_phi3d
!
! Args:
      logical,intent(in) :: wrprim,wrsech
!
! Local:
      integer :: i,j
!
      if (.not.wrprim.and..not.wrsech) return
!
! Transform 3d mag phim3d to geo phig3d, for write to histories:
      call mag2geo_3d(phim3d,phig3d,mag_phi3d,geo_phi3d,'PHIM3D  ')
!
!     do j=lat0,lat1
!       call addfld('PHIG3D','3D EPOT on Geographic grid','Volts',
!    |    phig3d(lon0:lon1,j,:),'lon',lon0,lon1,'ilev',1,nlevp1,j)
!     enddo ! j=lat0,lat1
!
! Update poten (pointer to f4d(i_poten)) for gather2root and write to
! history. In fields.F: poten(levd0:levd1,lond0:lond1,latd0:latd1,2)
!
      do j=lat0,lat1
        do i=lon0,lon1
          poten(:,i,j,itc) = phig3d(i,j,:)
        enddo
      enddo

!     do j=lat0,lat1
!       call addfld('POTEN_PHIG3D','3D EPOT on Geographic grid','Volts',
!    |    poten(:,lon0:lon1,j,itc),'ilev',1,nlevp1,'lon',lon0,lon1,j)
!     enddo ! j=lat0,lat1

      end subroutine prepare_phig3d
!-----------------------------------------------------------------------
      subroutine define_phim3d(ixt)
!
! Convert electric potential read from source history on geographic 
! grid to magnetic grid in phim3d. This is called once per run from 
! sub readsource (rdsource.F).
!
      use my_esmf,only: mag_ped,geo_ped
      use fields_module,only: poten
!
! Args:
      integer,intent(in) :: ixt
!
! Local:
      integer :: j
!
! Define phim3d from poten if doing dynamo:
! (poten was read from source history)

!     do j=lat0,lat1
!       write(6,"('define_phim3d: j=',i4,' poten(:,lon0:lon1,j,ixt)=',
!    |    2e12.4)") j,minval(poten(:,lon0:lon1,j,ixt)),
!    |                maxval(poten(:,lon0:lon1,j,ixt))
!     enddo

      call geo2mag_3d(poten(:,lon0:lon1,lat0:lat1,ixt),
     |  phim3d,geo_ped,mag_ped,'PHIG3D  ')

!     do j=mlat0,mlat1
!       write(6,"('define_phim3d: j=',i3,' phim3d(:,j,:) min,max=',
!    |    2e12.4)") j,minval(phim3d(:,j,:)),maxval(phim3d(:,j,:))
!     enddo

      end subroutine define_phim3d
!-----------------------------------------------------------------------
!
! May, 2012 btf:
! From here to end of module are routines from original serial code,
! which operate in global magnetic space, i.e., the remaining serial
! portion of the dynamo code. These routine will either be replaced by
! a new parallel solver, or they will be parallelized in shared memory
! with OpenMP directives.
!
!-----------------------------------------------------------------------
      subroutine stencils
      use cons_module,only: pi_dyn,dlatm,dlonm
!
! Locals:
      integer :: i,j,jj,jjj,j0,n,ncc,nmaglon,nmaglat
      real :: sym
      real :: cs(nmlat0)
      real :: c0_plt(nmlonp1,nmlat,10)
      real :: cofum_plt(nmlonp1,nmlat,9)
!
! Set index array nc and magnetic latitude cosine array:
! nc pointes to the start of the coefficient array for each level
      nc(1) = nc0
      nc(2) = nc1
      nc(3) = nc2
      nc(4) = nc3
      nc(5) = nc4
      nc(6) = ncee

      do j=1,nmlat0
        cs(j) = cos(pi_dyn/2.-(nmlat0-j)*dlatm)
      enddo ! j=1,nmlat0
!
! Set up difference coefficients. Replace zigm11 by A, zigm22 by B,
! zigmc by C, and zigm2 by D.
!
      j0 = nmlat0-nmlath
      do j=1,nmlath       !  1,49 (assuming nmlat=97)
        jj = nmlath+j-1   ! 49,97
        jjj = nmlath-j+1  ! 49,1
!
! factor 4 from 5-point diff. stencil
! Sigma_(phi lam)/( 4*Delta lam* Delta lon )
! Sigma_(phi lam)/( 4*Delta lam* Delta lon )
! Sigma_(lam lam)*cos(lam_m)*DT0DTS/(Delta lam)^2
! -zigmc_north = southern hemis. 49,1 equator-pole
! -zigm2_north = southern hemis. 49,1 equator-pole
!  zigm22 = southern hemis. 49,1 equator-pole
!
        do i=1,nmlonp1
          zigmc_glb(i,jj)   = (zigmc_glb(i,jj) +zigm2_glb(i,jj))/
     |      (4.*dlatm*dlonm)
          zigm2_glb(i,jj)   = zigmc_glb(i,jj)-2.*zigm2_glb(i,jj)/
     |      (4.*dlatm*dlonm)
          zigm22_glb(i,jj)  = zigm22_glb(i,jj)*cs(j0+j)/dlatm**2
          zigmc_glb(i,jjj)  = -zigmc_glb(i,jj)
          zigm2_glb(i,jjj)  = -zigm2_glb(i,jj)
          zigm22_glb(i,jjj) = zigm22_glb(i,jj)
        enddo ! i=1,nmlonp1
        if (j /= nmlath) then
!
! Sigma_(phi phi)/( cos(lam_m)*DT0DTS*(Delta lon)^2 )
! zigm11 = southern hemis. 49,1 equator-pole
!
          do i = 1,nmlonp1
            zigm11_glb(i,jj)  = zigm11_glb(i,jj)/(cs(j0+j)*dlonm**2)
            zigm11_glb(i,jjj) = zigm11_glb(i,jj)
          enddo
        endif
      enddo ! j=1,nmlath
!
! Set zigm11 to zero at megnetic poles to avoid floating exception
! (values at poles are not used):
!
      do i = 1,nmlonp1
        zigm11_glb(i,1) = 0.0
        zigm11_glb(i,nmlat) = 0.0
      enddo
!
! Save global zigm for comparison w/ odyn:
! real,dimension(nmlonp1,nmlat) :: zigm11_glb,zigm22_glb,zigmc_glb,zigm2_glb,rhs_glb
!
!     call addfld('ZIGM11_GLB',' ',' ',zigm11_glb,'mlon',1,nmlonp1,
!    |  'mlat',1,nmlat,0)
!     call addfld('ZIGM22_GLB',' ',' ',zigm22_glb,'mlon',1,nmlonp1,
!    |  'mlat',1,nmlat,0)
!     call addfld('ZIGMC_GLB',' ',' ',zigmc_glb,'mlon',1,nmlonp1,
!    |  'mlat',1,nmlat,0)
!     call addfld('ZIGM2_GLB',' ',' ',zigm2_glb,'mlon',1,nmlonp1,
!    |  'mlat',1,nmlat,0)
!     call addfld('RHS_GLB',' ',' ',rhs_glb,'mlon',1,nmlonp1,
!    |  'mlat',1,nmlat,0)
!
! Clear array for difference stencils at all levels:
      call clearcee(cee,nmlon0,nmlat0)
!
! isolve = 2 -> modified mudpack solver (modified and unmodified coefficients)
! (currently the only option available)
!
! real,dimension(nmlon0,nmlat0,9) :: cofum
      cofum(:,:,:) = 0. ! init
!
! Calculate contribution to stencils from each PDE coefficient
!
! Sigma_(phi phi)/( cos(lam_m)*dt0dts*(Delta lon)^2 ) 
      sym = 1.
      call stencmd(zigm11_glb,cs,nmlon0,nmlat0,sym,cee,1)
!
! Sigma_(lam lam)*cos(lam_m)*dt0dts/(Delta lam)^2
      sym = 1.
      call stencmd(zigm22_glb,cs,nmlon0,nmlat0,sym,cee,4)
!
! Sigma_(phi lam)/( 4*Delta lam* Delta lon )
      sym = -1.
      call stencmd(zigmc_glb,cs,nmlon0,nmlat0,sym,cee,2)
!
! Sigma_(lam phi)/( 4*Delta lam* Delta lon )
      sym = -1.
      call stencmd(zigm2_glb,cs,nmlon0,nmlat0,sym,cee,3)
!
! Insert RHS in finest stencil (formerly sub rths):
      do j = 1,nmlat0
        jj = nmlath-nmlat0+j
        do i = 1,nmlon0
          c0(i,j,10) = rhs_glb(i,jj)
        enddo ! i = 1,nmlon0
      enddo ! j = 1,nmlat0
      c0(nmlonp1,1,10) = c0(1,1,10)
!
! Set boundary condition at the pole:
      call edges(c0,nmlon0,nmlat0)
      call edges(c1,nmlon1,nmlat1)
      call edges(c2,nmlon2,nmlat2)
      call edges(c3,nmlon3,nmlat3)
      call edges(c4,nmlon4,nmlat4)
      call edges(cofum,nmlon0,nmlat0)
!
! Divide stencils by cos(lam_0) (not rhs):
      call divide(c0,nmlon0,nmlat0,nmlon0,nmlat0,cs,1)
      call divide(c1,nmlon1,nmlat1,nmlon0,nmlat0,cs,1)
      call divide(c2,nmlon2,nmlat2,nmlon0,nmlat0,cs,1)
      call divide(c3,nmlon3,nmlat3,nmlon0,nmlat0,cs,1)
      call divide(c4,nmlon4,nmlat4,nmlon0,nmlat0,cs,1)
      call divide(cofum,nmlon0,nmlat0,nmlon0,nmlat0,cs,0)
!
! Set value of solution to 1. at pole:
      do i=1,nmlon0
        c0(i,nmlat0,10) = 1.
      enddo
!
#if defined(INTERCOMM) || defined(CISMAH)
!
! 3/20/14 btf: cism coupling not yet ready in pdynamo (however,
!   this part will not have to be parallelized, since we are in
!   global domain here anyway).
!
! Convert LFM geographic potential to geomagnetic coordinates
!
      call cism_pot2mag
c      do jj=1,nmlat
c       do i=1,nmlonp1
c        write(6,*)phihm(i,jj),i,jj
c       enddo
c      enddo
#endif
!
! Modify stencils and RHS so that the NH high lat potential is inserted at
!  high latitude.  The SH high lat potential will be added back later.
!  pfrac = fraction of dynamo in solution in the NH. = 1 low lat, = 0 hi lat
!    cons_module: crit(1)=15, crit(2)=30 deg colats, or hi-lat > 75 deg,
!      dynamo < 60 deg, and combination between 60-75 mag lat.
! The dynamo is symmetric about the magnetic equator, but the high latitude
!  is anti-symmetric in both hemispheres.  However, since Mudpack uses the
!  NH potential pattern, then the SH potential pattern must be added
!  back into the 2-D phim before the call threed, and before it is
!  transformed to geographic coordinates.
!
! (In original dynamo, this was conditional on .not.mod_heelis)
!
      ncc = 1
      nmaglon = nmlon0
      nmaglat = nmlat0
      do n=1,5
        call stenmd(nmaglon,nmaglat,cee(ncc),phihm(1,nmlat0),pfrac)
        ncc = ncc+9*nmaglon*nmaglat
        if (n==1) ncc = ncc+nmaglon*nmaglat ! rhs is in 10th slot
        nmaglon = (nmaglon+1)/2
        nmaglat = (nmaglat+1)/2
      enddo ! n=1,5

!
! (note nmlat0==nmlath, and nmlon0==nmlon+1)
! real :: c0_plt(nmlonp1,nmlat,10)
! real :: c0(nmlon0,nmlat0,10),
!
      c0_plt = 0.
      do n=1,10
        c0_plt(:,1:nmlath,n) = c0(:,:,n)

!       do j=1,nmlath
!         write(6,"('n=',i2,' j=',i3,' c0(1,j,n)=',e15.7,
!    |      ' c0(nmlon:nmlonp1,j,n)=',2e15.7)") 
!    |      n,j,c0_plt(1,j,n),c0_plt(nmlon:nmlonp1,j,n)
!       enddo

      enddo

!     call addfld('C0_01',' ',' ',c0_plt(:,:,1),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('C0_02',' ',' ',c0_plt(:,:,2),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('C0_03',' ',' ',c0_plt(:,:,3),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('C0_04',' ',' ',c0_plt(:,:,4),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('C0_05',' ',' ',c0_plt(:,:,5),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('C0_06',' ',' ',c0_plt(:,:,6),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('C0_07',' ',' ',c0_plt(:,:,7),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('C0_08',' ',' ',c0_plt(:,:,8),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('C0_09',' ',' ',c0_plt(:,:,9),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('C0_10',' ',' ',c0_plt(:,:,10),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)

      cofum_plt = 0.
      do n=1,9
        cofum_plt(:,1:nmlath,n) = cofum(:,:,n)
      enddo
!     call addfld('COFUM_01',' ',' ',cofum_plt(:,:,1),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('COFUM_02',' ',' ',cofum_plt(:,:,2),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('COFUM_03',' ',' ',cofum_plt(:,:,3),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('COFUM_04',' ',' ',cofum_plt(:,:,4),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('COFUM_05',' ',' ',cofum_plt(:,:,5),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('COFUM_06',' ',' ',cofum_plt(:,:,6),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('COFUM_07',' ',' ',cofum_plt(:,:,7),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('COFUM_08',' ',' ',cofum_plt(:,:,8),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('COFUM_09',' ',' ',cofum_plt(:,:,9),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)

      end subroutine stencils
!-----------------------------------------------------------------------
      subroutine clearcee(cee,nlon0,nlat0)
!
! Zero C arrays for stencil coefficients.
! Cee will contain:
!   c0(nmlon0,nmlat0,10), c1(nmlon1,nmlat1,9), c2(nmlon2,nmlat2,9),
!   c3(nmlon3,nmlat3,9),  c4(nmlon4,nmlat4,9)
!
! Args:
      integer,intent(in) :: nlon0,nlat0
      real,intent(out) :: cee(*)
!
! Local:
      integer :: nlon,nlat,n,m,i
      real :: time0,time1

!     if (dynamo_timing) call timer(time0,time1,"clearcee",0,0)
!
! Compute total size of cee
      nlon = nlon0
      nlat = nlat0
      n = 0
      do m=1,5 ! 5 resolution levels
        n = n+nlon*nlat
        nlon = (nlon+1)/2
        nlat = (nlat+1)/2
      enddo ! m=1,5 ! 5 resolution levels
      n = 9*n+nlon0*nlat0
!
! Clear cee:
      do i=1,n
        cee(i) = 0.
      enddo

      end subroutine clearcee
!-----------------------------------------------------------------------
      subroutine stencmd(zigm,cs,nlon0,nlat0,sym,cee,ncoef)
!
! Calculate contribution fo 3 by 3 stencil from coefficient zigm
! at each grid point and level. 
!
! Args:
      integer,intent(in) :: 
     |  nlon0, ! longitude dimension of finest grid level
     |  nlat0, ! latitude dimension of finest grid level
     |  ncoef  ! integer identifier of coefficient
      real,intent(in) :: 
     |  zigm(nlon0,nlat0), ! coefficients (nlon0+1/2,(nlat0+1)/2) 
     |  sym,               !  1. if zigm symmetric w.r.t. equator 
                           ! -1. for antisymmetry
     |  cs(nlat0)
      real,intent(inout) :: 
     |  cee(*)  ! output stencil array consisting of c0,c1,c2,c3,c4
!
! Local:
      integer :: nc,nlon,nlat,n
      real :: wkarray(-15:nmlon0+16,nmlat0)
!
! Perform half-way interpolation and extend zigm in wkarray:
!
      call htrpex(zigm,nlon0,nlat0,sym,wkarray)
!
! Calculate contribution to stencil for each grid point and level:
!
      nc = 1
      nlon = nlon0
      nlat = nlat0
!
! Calculate modified and unmodified stencil on finest grid
!
      call cnmmod(nlon0,nlat0,nlon,nlat,cee(nc),ncoef,wkarray,cofum)

!     if (dynamo_timing) call timer(time0,time1,"stencmd",0,0)
!
! Stencils on other grid levels remain the same.
      nc = nc+10*nlon*nlat 
      nlon = (nlon+1)/2
      nlat = (nlat+1)/2
!
      do n=2,5
        call cnm(nlon0,nlat0,nlon,nlat,cee(nc),ncoef,wkarray)
        nc = nc+9*nlon*nlat
        if (n==1) nc = nc+nlon*nlat
        nlon = (nlon+1)/2
        nlat = (nlat+1)/2
      enddo ! n=1,5
      end subroutine stencmd
!-----------------------------------------------------------------------
      subroutine htrpex(coeff,nmlon0,nmlat0,sym,wkarray)
!
! Perform half-way interpolation on array coeff and extend over 16 grid
! points. Result returned in wkarray.
!
! Args:
      integer,intent(in) :: nmlon0,nmlat0
      real,intent(in) :: coeff(nmlon0,nmlat0),sym
      real,intent(out) :: wkarray(-15:nmlon0+16,nmlat0)
!
! Local:
      integer :: i,j,jj
!
! Copy coeff into positions in wkarray:
      do j=1,nmlat0
        jj = nmlat0-j+1
        do i=1,nmlon0
          wkarray(i,j) = sym*coeff(i,jj)
        enddo ! i=1,nmlon0
      enddo ! j=1,nmlat0
!
! Extend over 32 grid spaces to allow for a total of 5 grid levels:
      do i=1,16
        do j=1,nmlat0
          wkarray(1-i,j) = wkarray(nmlon0-i,j) 
          wkarray(nmlon0+i,j) = wkarray(1+i,j)
        enddo ! j=1,nmlat0
      enddo ! i=1,16
      end subroutine htrpex
!-----------------------------------------------------------------------
      subroutine cnm(nlon0,nlat0,nlon,nlat,c,ncoef,wkarray)
!
! Compute contribution to stencil from zigm(ncoef) on grid nlon by nlat,
! Finest grid is nlon0 by nlat0.
!
! Args:
      integer,intent(in) :: 
     |  nlon0,nlat0, ! finest grid dimensions
     |  nlon,nlat    ! output grid dimensions
      real,intent(in) :: wkarray(-15:nmlon0+16,nmlat0)
!
! ncoef: integer id of coefficient:
! ncoef = 1 for zigm11
! ncoef = 2 for zigm12 (=zigmc+zigm2)
! ncoef = 3 for zigm21 (=zigmc-zigm2)
! ncoef = 4 for zigm22
!
      integer,intent(in) :: ncoef
      real,intent(inout) :: 
     |  c(nlon,nlat,*)    ! output array for grid point stencils at
                          ! resolution nlon x nlat
!
! Local:
      integer :: i,j,nint,i0,j0
      real,parameter :: pi=3.141592654
      real :: wk(nlon0,3)
      real :: time0,time1

!     if (dynamo_timing) call timer(time0,time1,"cnm",0,0)
!
! Compute separation of grid points of resolution nlon x nlat within
! grid of resolution nlon0,nlat0. Evaluate dlon and dlat, grid spacing 
! of nlon x nlat.
!
      nint = (nlon0-1)/(nlon-1)
!
! Scan wkarray nlon x nlat calculating and adding contributions to stencil
! from zigm(ncoef)
      i0 = 1-nint
      j0 = 1-nint
!
! zigm11:
! am 2001-6-27 include boundary condition at equator
      if (ncoef==1) then
        do j = 1,nlat-1
          do i = 1,nlon
            c(i,j,1) = c(i,j,1)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            c(i,j,5) = c(i,j,5)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            c(i,j,9) = c(i,j,9)-.5*(wkarray(i0+(i+1)*nint,j0+j*nint)+
     |        2.*wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigm12 (=zigmc+zigm2)
      elseif (ncoef==2) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,2) = c(i,j,2)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            c(i,j,4) = c(i,j,4)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            c(i,j,6) = c(i,j,6)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            c(i,j,8) = c(i,j,8)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            wk(i,1) = .5*(wkarray(i0+(i+1)*nint,j0+j*nint)-
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            wk(i,2) = (c(i,j,3)+wk(i,1))*(c(i,j,7)-wk(i,1))
            wk(i,3) = sign(wk(i,1),c(i,j,3)+c(i,j,7))
            if (wk(i,2) >= 0.) wk(i,3) = 0.
            c(i,j,3) = c(i,j,3)+wk(i,1)+wk(i,3)
            c(i,j,7) = c(i,j,7)-wk(i,1)+wk(i,3)
            c(i,j,9) = c(i,j,9)-2.*wk(i,3)
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigm21 (=zigmc-zigm2)
      elseif (ncoef==3) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,2) = c(i,j,2)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j+1)*nint))
            c(i,j,4) = c(i,j,4)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j+1)*nint))
            c(i,j,6) = c(i,j,6)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j-1)*nint))
            c(i,j,8) = c(i,j,8)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j-1)*nint))
            wk(i,1) = .5*(wkarray(i0+i*nint,j0+(j+1)*nint)-
     |        wkarray(i0+i*nint,j0+(j-1)*nint))
            wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
            wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
            if (wk(i,2) >= 0.) wk(i,3) = 0.
            c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
            c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
            c(i,j,9) = c(i,j,9)-2.*wk(i,3)
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
        j = 1
        do i=1,nlon
          c(i,j,2) = c(i,j,2)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |      wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,4) = c(i,j,4)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |      wkarray(i0+i*nint,j0+(j+1)*nint))
          wk(i,1) = .5*(wkarray(i0+i*nint,j0+j*nint)+
     |      wkarray(i0+i*nint,j0+(j+1)*nint))
          wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
          wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
          if (wk(i,2) >= 0.) wk(i,3) = 0.
          c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
          c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
          c(i,j,9) = c(i,j,9)-2.*wk(i,3)
        enddo ! i=1,nlon
!
! zigm22:
      elseif (ncoef==4) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,3) = c(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j+1)*nint))
            c(i,j,7) = c(i,j,7)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j-1)*nint))
            c(i,j,9) = c(i,j,9)-.5*(wkarray(i0+i*nint,j0+(j-1)*nint)
     |        +2.*wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j+1)*nint))
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
        j = 1
        do i=1,nlon
          c(i,j,3) = c(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |      +wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,9) = c(i,j,9)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |                            wkarray(i0+i*nint,j0+(j+1)*nint))
        enddo ! i=1,nlon
      endif ! ncoef
      end subroutine cnm
!-----------------------------------------------------------------------
      subroutine cnmmod(nlon0,nlat0,nlon,nlat,c,ncoef,wkarray,cofum)
!
! Compute contribution to stencil from zigm(ncoef) on grid nlon by nlat,
! Finest grid is nlon0 by nlat0.
!
! Args:
      integer,intent(in) :: 
     |  nlon0,nlat0, ! finest grid dimensions
     |  nlon,nlat    ! output grid dimensions
      real,intent(in) :: wkarray(-15:nmlon0+16,nmlat0)
      real,dimension(nmlon0,nmlat0,9),intent(inout) :: cofum
!
! ncoef: integer id of coefficient:
! ncoef = 1 for zigm11
! ncoef = 2 for zigm12 (=zigmc+zigm2)
! ncoef = 3 for zigm21 (=zigmc-zigm2)
! ncoef = 4 for zigm22
!
      integer,intent(in) :: ncoef
      real,intent(inout) :: 
     |  c(nlon,nlat,*)    ! output array for grid point stencils at
                          ! resolution nlon x nlat
!
! Local:
      integer :: i,j,nint,i0,j0
      real,parameter :: pi=3.141592654
      real :: wk(nlon0,3)
!
! Compute separation of grid points of resolution nlon x nlat within
! grid of resolution nlon0,nlat0. Evaluate dlon and dlat, grid spacing 
! of nlon x nlat.
!
      nint = (nlon0-1)/(nlon-1)
!
! Scan wkarray nlon x nlat calculating and adding contributions to stencil
! from zigm(ncoef)
      i0 = 1-nint
      j0 = 1-nint
!
! zigm11:
! am 2001-6-27 include boundary condition at equator
      if (ncoef==1) then
        do j = 1,nlat-1
          do i = 1,nlon
            c(i,j,1) = c(i,j,1)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            c(i,j,5) = c(i,j,5)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            c(i,j,9) = c(i,j,9)-.5*(wkarray(i0+(i+1)*nint,j0+j*nint)+
     |        2.*wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
!
! Unmodified:
            cofum(i,j,1) = cofum(i,j,1)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+(i+1)*nint,j0+j*nint))
            cofum(i,j,5) = cofum(i,j,5)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+(i-1)*nint,j0+j*nint))
            cofum(i,j,9) = cofum(i,j,9)-
     |        .5*(wkarray(i0+(i+1)*nint,j0+j*nint)+
     |         2.*wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigm12 (=zigmc+zigm2)
      elseif (ncoef==2) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,2) = c(i,j,2)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            c(i,j,4) = c(i,j,4)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            c(i,j,6) = c(i,j,6)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            c(i,j,8) = c(i,j,8)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            wk(i,1) = .5*(wkarray(i0+(i+1)*nint,j0+j*nint)-
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
!
! Unmodified:
            cofum(i,j,2) = c(i,j,2)
            cofum(i,j,4) = c(i,j,4)
            cofum(i,j,6) = c(i,j,6)
            cofum(i,j,8) = c(i,j,8)
            cofum(i,j,3) = cofum(i,j,3)+wk(i,1)
            cofum(i,j,7) = cofum(i,j,7)-wk(i,1)
!
            wk(i,2) = (c(i,j,3)+wk(i,1))*(c(i,j,7)-wk(i,1))
            wk(i,3) = sign(wk(i,1),c(i,j,3)+c(i,j,7))
            if (wk(i,2) >= 0.) wk(i,3) = 0.
            c(i,j,3) = c(i,j,3)+wk(i,1)+wk(i,3)
            c(i,j,7) = c(i,j,7)-wk(i,1)+wk(i,3)
            c(i,j,9) = c(i,j,9)-2.*wk(i,3)
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigm21 (=zigmc-zigm2)
      elseif (ncoef==3) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,2) = c(i,j,2)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j+1)*nint))
            c(i,j,4) = c(i,j,4)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j+1)*nint))
            c(i,j,6) = c(i,j,6)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j-1)*nint))
            c(i,j,8) = c(i,j,8)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j-1)*nint))
            wk(i,1) = .5*(wkarray(i0+i*nint,j0+(j+1)*nint)-
     |        wkarray(i0+i*nint,j0+(j-1)*nint))
!
! Unmodified:
            cofum(i,j,2) = c(i,j,2)
            cofum(i,j,4) = c(i,j,4)
            cofum(i,j,6) = c(i,j,6)
            cofum(i,j,8) = c(i,j,8)
            cofum(i,j,1) = cofum(i,j,1)+wk(i,1)
            cofum(i,j,5) = cofum(i,j,5)-wk(i,1)
!
            wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
            wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
            if (wk(i,2) >= 0.) wk(i,3) = 0.
            c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
            c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
            c(i,j,9) = c(i,j,9)-2.*wk(i,3)
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
        j = 1
        do i=1,nlon
          c(i,j,2) = c(i,j,2)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |      wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,4) = c(i,j,4)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |      wkarray(i0+i*nint,j0+(j+1)*nint))
          wk(i,1) = .5*(wkarray(i0+i*nint,j0+j*nint)+
     |      wkarray(i0+i*nint,j0+(j+1)*nint))

          cofum(i,j,2) = c(i,j,2)
          cofum(i,j,4) = c(i,j,4)
          cofum(i,j,1) = cofum(i,j,1)+wk(i,1)
          cofum(i,j,5) = cofum(i,j,5)-wk(i,1)

          wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
          wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
          if (wk(i,2) >= 0.) wk(i,3) = 0.
          c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
          c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
          c(i,j,9) = c(i,j,9)-2.*wk(i,3)
        enddo ! i=1,nlon
!
! zigm22:
      elseif (ncoef==4) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,3) = c(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j+1)*nint))
            c(i,j,7) = c(i,j,7)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j-1)*nint))
            c(i,j,9) = c(i,j,9)-.5*(wkarray(i0+i*nint,j0+(j-1)*nint)
     |        +2.*wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j+1)*nint))
!
! Unmodified:
            cofum(i,j,3) = cofum(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j+1)*nint))
            cofum(i,j,7) = cofum(i,j,7)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j-1)*nint))
            cofum(i,j,9) = cofum(i,j,9)-.5*(wkarray(i0+i*nint,j0+(j-1)*
     |        nint)+2.*wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j+1)*nint))

          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
        j = 1
        do i=1,nlon
          c(i,j,3) = c(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |      +wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,9) = c(i,j,9)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |                            wkarray(i0+i*nint,j0+(j+1)*nint))
          cofum(i,j,3) = cofum(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |                            wkarray(i0+i*nint,j0+(j+1)*nint))
          cofum(i,j,9) = cofum(i,j,9)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |                            wkarray(i0+i*nint,j0+(j+1)*nint))

        enddo ! i=1,nlon
      endif ! ncoef
      end subroutine cnmmod
!-----------------------------------------------------------------------
      subroutine ceee(cee,nx,ny,cf)
!
! Called from mudpack solvers to transfer coefficients.
!
      save
      integer,intent(in) :: nx,ny
      real,intent(in) :: cee(nx,ny,*)
      real,intent(out) :: cf(nx,ny,*)
      integer :: i,j,n
      real :: time0,time1

      do n = 1,9
        do j = 1,ny
          do i = 1,nx
            cf(i,j,n) = cee(i,j,n)
          enddo
!         write(6,"('ceee: n=',i2,' j=',i3,' nx=',i3,' ny=',i3,
!    |      ' cf(:,j,n)=',/,(6e12.4))") n,j,nx,ny,cf(:,j,n)
        enddo
      enddo
      end subroutine ceee
!-----------------------------------------------------------------------
      subroutine edges(c,nlon,nlat)
!
! Insert equatorial and polar boundary conditions in stencil c(nlon,nlat,9)
!
! Args:
      integer,intent(in) :: nlon,nlat
      real,intent(out) :: c(nlon,nlat,*)
!
! Local:
      integer :: n,i

      do n=1,8
        do i=1,nlon
          c(i,nlat,n) = 0.
        enddo
      enddo
      do i=1,nlon
        c(i,nlat,9) = 1.
      enddo
      end subroutine edges
!-----------------------------------------------------------------------
      subroutine divide(c,nlon,nlat,nlon0,nlat0,cs,igrid)
!
! Divide stencil C by cos(theta(i,j))
!
! Args:
      integer,intent(in) :: nlon,nlat,nlon0,nlat0,igrid
      real,intent(in) :: cs(*)
      real,intent(out) :: c(nlon,nlat,*)
!
! Local:
      integer :: nint,j0,n,j,i
!
      nint = (nlon0-1)/(nlon-1)
      j0 = 1-nint
      do n = 1,9
        do j = 1,nlat-1
          do i = 1,nlon
            c(i,j,n) = c(i,j,n)/(cs(j0+j*nint)*nint**2)
          enddo ! i = 1,nlon
        enddo ! j = 1,nlat-1
      enddo ! n = 1,9
!
! This is unnecessary since cs(1)==1. 
! (this was in serial dynamo of timegcm, but not in tiegcm).
!
!     if (nint==1.and.igrid > 0) then
!       do i = 1,nlon
!         c(i,1,10) = c(i,1,10)/cs(1)
!       enddo ! i = 1,nlon
!     endif
      end subroutine divide
!-----------------------------------------------------------------------
      subroutine stenmd(inlon,inlat,c,phihm,pfrac)
      use cons_module,only: dlatm,dtr
!
! Modify stencil to set potential to heelis value within auroral circle.
!
! Args:
      integer,intent(in) :: inlon,inlat
      real,intent(inout) :: c(inlon,inlat,*)
      real,dimension(nmlon0,nmlat0),intent(in) :: 
     |  phihm,  ! heelis potential (from subs potm, flwv32)
     |  pfrac   ! fractional presence of dynamo (from sub colath)
!
! Local:
      integer :: nint,i0,j0,i,j,n,jj
!
! Compute separation of grid points for this resolution:
      nint = (nmlon0-1)/(inlon-1)
      i0 = 1-nint
      j0 = 1-nint
!
! If nint==1, then we are at the highest resolution.
! Correct RHS, which is in c(10)
!
      if (nint==1) then
        do j=1,inlat
          do i=1,inlon
            c(i,j,10) = pfrac(i,j)*c(i,j,10)+(1.-pfrac(i,j))*c(i,j,9)*
     |        (dlatm/(10.*dtr))**2*phihm(i,j)
          enddo ! i=1,inlon
        enddo ! j=1,inlat
      endif
!
! Modify stencil, c(i,j,n),n=1,9:
!
      if (nint==1) then
        do j=1,inlat
          jj = j0+j*nint
          do n = 1,8
            do i = 1,inlon
              c(i,j,n) = c(i,j,n)*pfrac(i0+i*nint,jj)
              cofum(i,j,n) = cofum(i,j,n)*pfrac(i0+i*nint,jj)
            enddo ! i = 1,inlon
          enddo ! n = 1,8
          do i = 1,inlon
            c(i,j,9) = c(i,j,9)*pfrac(i0+i*nint,jj)+
     |        (1.-pfrac(i0+i*nint,jj))*c(i,j,9)*
     |        (dlatm*float(nint)/(10.*dtr))**2
            cofum(i,j,9) =cofum(i,j,9)*pfrac(i0+i*nint,jj)+
     |        (1.-pfrac(i0+i*nint,jj))*cofum(i,j,9)*
     |        (dlatm*float(nint)/(10.*dtr))**2
          enddo ! i = 1,inlon
        enddo ! j=1,inlat
      else ! nint /= 1
        do j=1,inlat
          jj = j0+j*nint
          do n = 1,8
            do i = 1,inlon
              c(i,j,n) = c(i,j,n)*pfrac(i0+i*nint,jj)
            enddo ! i = 1,inlon
          enddo ! n = 1,8
          do i = 1,inlon
            c(i,j,9) = c(i,j,9)*pfrac(i0+i*nint,jj)+
     |        (1.-pfrac(i0+i*nint,jj))*c(i,j,9)*
     |        (dlatm*float(nint)/(10.*dtr))**2
          enddo ! i = 1,inlon
        enddo ! j=1,inlat
      endif ! nint
      end subroutine stenmd
!-----------------------------------------------------------------------
      subroutine solver(cofum,c0)
!
! Call mudpack to solve PDE. Solution is returned in rim:
!     real,dimension(nmlonp1,nmlat,2) :: rim
! Mudpack solvers:
!   isolve = 0  org. mud v5.      (modified stencils neq direct solution)
!   isolve = 1  muh hybrid solver (no convergence => only as direct solver)
!   isolve = 2  modified mud      (residual calculated with unmodified stencils
!                                  same solution as with direct solver, if same
!                                  coefficient matrix is used)
! Only isolve==2 is supported in pdynamo.
!
! Args:
      real,dimension(nmlon0,nmlat0,9),intent(in) :: cofum
      real,intent(in) :: c0(nmlon0,nmlat0,10)
!
! Locals:
      integer :: i,j,jntl,ier,isolve
      real :: l2norm
      real,dimension(nmlon0,nmlat0)   :: ressolv
      real,dimension(nmlon0,nmlat0,9) :: cofum_solv

! Module data:
! real,dimension(nmlonp1,nmlat,2) :: rim_glb ! pde solver output

!     call addfld('RIM1_PRESOLV',' ','A/m',rim_glb(:,:,1),'mlon',1,
!    |  nmlonp1,'mlat',1,nmlat,0)
!     call addfld('RIM2_PRESOLV',' ','A/m',rim_glb(:,:,2),'mlon',1,
!    |  nmlonp1,'mlat',1,nmlat,0)

      jntl = 0
      ier = 0
      isolve = 2
      call mudmod(rim_glb,phisolv,jntl,isolve,ier)! solver in mudmod.F
      if (ier < 0 ) then       ! not converged
        write(6,*) 'muh: use direct solver'
        call muh(rim_glb,jntl)            ! solver in mud.F
      endif

!     call addfld('RIM1_SOLV',' ','A/m',rim_glb(:,:,1),'mlon',1,nmlonp1,
!    |  'mlat',1,nmlat,0)
!     call addfld('RIM2_SOLV',' ','A/m',rim_glb(:,:,2),'mlon',1,nmlonp1,
!    |  'mlat',1,nmlat,0)

      l2norm=0.
      ressolv = 0.0
      do j = 1,nmlat0
        do i = 1,nmlon0-1
          cofum_solv(i,j,:)=  cofum(i,j,:)
!
! fields: phisolv(0:nmlonp1,0:nmlat+1)       ! 2d solution/ electric potential
!
	  ressolv(i,j) = (
     +    cofum_solv(i,j,1)*phisolv(i+1,j)+
     +    cofum_solv(i,j,2)*phisolv(i+1,j+1)+
     +    cofum_solv(i,j,3)*phisolv(i,j+1)+
     +    cofum_solv(i,j,4)*phisolv(i-1,j+1)+
     +    cofum_solv(i,j,5)*phisolv(i-1,j)+
     +    cofum_solv(i,j,6)*phisolv(i-1,j-1)+
     +    cofum_solv(i,j,7)*phisolv(i,j-1)+
     +    cofum_solv(i,j,8)*phisolv(i+1,j-1)+
     +    cofum_solv(i,j,9)*phisolv(i,j))   

          ressolv(i,j) = c0(i,j,10)-ressolv(i,j)
	  l2norm = l2norm + ressolv(i,j)*ressolv(i,j)
        end do
      end do
!     write(6,*) 'L2norm (global root task) = ',l2norm

      end subroutine solver
!-----------------------------------------------------------------------

#if defined(INTERCOMM) || defined(CISMAH)
      subroutine cism_pot2mag
!
! This subroutine converts high latitude potential in geographic coordinate 
! obtained form the CMIT M-I couplier to geomagnetci coordinate to be used
! in dynamo calculations
!
! Uses
!
      use params_module,only: nlonp1
      use cism_coupling_module, only:gpotm
      use magfield_module,only: ig,jg,wt
!
! Local:
!
      integer :: jj

      do jj=mlat0,mlat1
         call geo2mag(phihm(1,jj),gpotm(3:nlonp1+2,:),
     |        ig,jg,wt,nlonp1,nmlonp1,nmlon,nmlat,jj)
         ! Periodic point
         phihm(nmlonp1,jj) = phihm(1,jj)
      enddo

      end subroutine cism_pot2mag
!-----------------------------------------------------------------------
      subroutine geo2mag(fmag,fgeo,long,latg,wght,nlonp1_geo,nlonp1_mag,
     |  nlon_mag,nlat_mag,lat)
!
! Transform field fgeo on geographic grid to geomagnetic grid using
!   indices long,latg and weights wght. Return field fmag on magnetic
!   grid.
!
! Args:
      integer,intent(in) :: nlonp1_geo,nlonp1_mag,nlon_mag,nlat_mag,lat
      integer,dimension(nlonp1_mag,nlat_mag),intent(in) :: long,latg
      real,intent(in) :: fgeo(nlonp1_geo,*),wght(4,nlonp1_mag,nlat_mag)
      real,intent(out) :: fmag(nlonp1_mag,*)
!
! Local:
      integer :: i
!
      do i=mlon0,mlon1
        fmag(i,1) =
     |    fgeo(long(i,lat)  ,latg(i,lat)  )*wght(1,i,lat)+
     |    fgeo(long(i,lat)+1,latg(i,lat)  )*wght(2,i,lat)+
     |    fgeo(long(i,lat)+1,latg(i,lat)+1)*wght(3,i,lat)+
     |    fgeo(long(i,lat)  ,latg(i,lat)+1)*wght(4,i,lat)
!       if (iprint > 0) write(6,"('geo2mag: i=',i3,' lat=',i3,' long=',
!    |    i3,' latg=',i3,' wght=',4e12.4,' fgeo=',e12.4,' fmag=',
!    |    e12.4)") i,lat,long(i,lat),latg(i,lat),wght(:,i,lat),
!    |    fgeo(long(i,lat),latg(i,lat)),fmag(i,1)
      enddo
      end subroutine geo2mag
#endif
!-----------------------------------------------------------------------
      end module pdynamo_module
