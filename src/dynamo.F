!
      module dynamo_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
!----------------------------------------------------------------------- 
!
! Module for electrodynamic dynamo (the "e" in tiegcm). 
!
! 5/02 btf: This is the "new" dynamo code, adapted by Astrid Maute
! and Art Richmond from the original ("old") dynamo in tgcm15,
! written by E. Cicely Ridley completed in about 1992. 
!
! Sub prep_dynamo is called from advance. Prep_dynamo prepares dynamo 
! input fields and gathers them to the root task (mp_dynamo_gather).
! Next, advance calls sub dynamo (the dynamo driver), which controls
! the dynamo subroutine calls. Prep_dynamo and the dynamo driver are
! called once per timestep, after sub dynamics.
!
! Since the dynamo is serial (as of 4/02), the master task must receive
! subdomain data from other tasks for input fields to the dynamo 
! (on the geographic grid). These fields are exchanged via sub 
! mp_dynamo_gather which is called from prep_dynamo below.
! 
      use params_module,only: 
     |  nlon,   ! number of geographic longitudes (at 5 deg, nlon==72)
     |  nlonp1, ! nlon+1
     |  nlonp4, ! nlon+4
     |  nlonp2, ! nlon+2
     |  nlev,   ! number of pressure levels (e.g., 28)
     |  nlevp1, ! nlev+1
     |  nlat,   ! number of geographic latitudes (at 5 deg, nlat==36)
     |  nlatp1, ! nlat+1
     |  nlatp2, ! nlat+2
     |  nmlon,  ! number of geomagnetic grid longitudes
     |  nmlonp1,! nmlon+1
     |  nmlat,  ! number of geomagnetic grid latitudes
     |  nmlatp1,! nmlat+1
     |  nmlath, ! (nmlat+1)/2 (index to magnetic equator)
     |  nmlev,  ! number of geomagnetic pressure levels (nmlev==nlevp1+3)
     |  nmlevp1,! number of geomagnetic midpoint levels
     |  nimlevp1,! number of geomagnetic interface levels
     |  spval
      use diags_module,only: mkdiag_ED12,mkdiag_PHIM2D
!    
! J(mag.pressure, gravity) to dynamo RHS
      use magpres_g_module,only: j_pg,tpg1,tpg2,je1oD_geo,je2oD_geo
!
! Routine to add fields to secondary histories:
      use addfld_module,only: addfld
      implicit none
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
!
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
!
! Unmodified coefficients for using modified mudpack:
      real,dimension(nmlon0,nmlat0,9) :: cofum
!
! phim: Single level dynamo potential in geomagnetic coordinates, as
!   output by PDE solver mud (formerly in dynphi.h):
! The dynamo is symmetric about the magnetic equator, but the high latitude
!  is anti-symmetric in both hemispheres.  However, since Mudpack uses the
!  NH potential pattern, then the SH potential pattern must be added
!  back into the 2-D phim before the call threed, and before it is
!  transformed to geographic coordinates.
! jn index for NH part of potential (nmlat down to ~nmlat0)
! jp index for NH pfrac (nmlat0 down to 1), the fraction of the dynamo
!
      real,dimension(nmlonp1,nmlat) :: phim
      integer :: jn,jp
!
! Inputs to the dynamo on the geographic grid. (formerly in dynamo.h):
      real,dimension(nlonp4,nlat,nlevp1) ::
     |  sigma_ped,        ! Pedersen conductivity (see lamdas.F)
     |  sigma_hall,       ! Hall conductivity     (see lamdas.F)
     |  zpoten,           ! geopotential height
     |  scheight,         ! scale height RT/(MBAR*G)
     |  unvel,vnvel,wnvel,! u,v,w neutral velocities (cm/s)
     |  je1oD_full,je2oD_full ! J(mag.pres,gravity) on full domain
!
! Fields on geomagnetic grid (at a single latitude):
! (suffix "m" is for "mag")
! (in ealier versions, these were in common /transmag_priv/)
!
      real,dimension(nmlonp1,-2:nlev) ::
     |  sigma_pedm,       ! Pedersen conductivity (was sigma1m)
     |  sigma_hallm,      ! Hall conductivity     (was sigma2m)
     |  wnvelm,           ! vertical velocity     (was wm)
     |  vxbm
      real,dimension(nmlonp1,-2:nlevp1) :: 
     |  zpotenm           ! geopotential          (was zm)
      real,dimension(nmlonp1,-2:nlev,2) :: adotvm
      real,dimension(nmlonp1,2)         :: adotam
      real,dimension(nmlonp1)           :: a1a2m, siniam, bmodm, pm
      
      real,dimension(nmlonp1,-2:nlev) :: ! J(mag.pres,gravity) mag. grid
     |  je1_pg,	 ! J(gravity)
     |  je2_pg	 ! J(mag.pres)
!
! Coefficients and RHS terms for PDE on geomagnetic grid:
! (formerly in coefm.h)
!
      real,dimension(nmlonp1,nmlat) ::
     |  zigm11,  ! sigma11*cos(theta0)
     |  zigmc,   ! sigmac
     |  zigm2,   ! sigma2
     |  zigm22   ! sigma22/cos(theta0)
#if defined(INTERCOMM) || defined(CISMAH)
      real,dimension(nmlonp1,nmlat) :: Zigm1 ! CISM, height-intergrated
!                                            ! Pedersen conductivity
!
! Parameters to be passed to LFM for CISM. Note the dimension of these variables is
! defined to match mag2geo specifications, see detial in cons.F and Apex.F.
!
      real,dimension(nlonp4,0:nlatp1) ::
     |  gzigm2,   ! sigma2    ! CISM height integrated Hall conductivity geographic
     |  gzigm1,   ! sigma1    ! CISM height integrated Pedersen conductivity geographic
     |  gnsrhs    ! nsrsh     ! CISM height integrated neutral wind currents geographic
#endif
!
! rim(1)=id(1), rim(2)=id(2)/cos(theta0)
      real,dimension(nmlonp1,nmlat,2) :: rim 
      real,dimension(nmlonp1,nmlath)  :: rhs ! right-hand side 
      real,dimension(nmlonp1,-2:nlev) :: rhs_plt ! diag
!
! 3-d geopotential heights in geomagnetic space (was zzm in dynphi.h)
      real,dimension(nmlonp1,nmlat,-2:nlevp1) :: zpotenm3d
!
! isolve = 0 -> original mud version 5.
! isolve = 1 -> muh hybrid solver (only as direct solver -- slow)
! isolve = 2 -> modified mudpack solver (mudmod)
! Default is isolve=2 for DYNAMO=2 (new dynamo), however under SGI/IRIX,
!   the hybrid solver is used because SGI cr
!
#ifdef IRIX
      integer,parameter :: isolve = 1
#else
      integer,parameter :: isolve = 2 ! default is 2 for new dynamo
!     integer,parameter :: isolve = 1 
!     integer,parameter :: isolve = 0 ! no longer available
#endif
! 
! Flag for original Heelis contribution to dynamo
!  p*rhs = (1-p)[ sigma(phi,phi)/cos(lam_m)*(delta lam_m/delta phi_m)^2 + 
!  sigma(lam,lam)*cos(lam_m) ]1/w^2* [poten_H - poten]
!
! or simplified
!  add to rhs: - d /d phi_m * [sigma_R /cos lam_m * d (poten-poten_R)/d phi_m]-
!     - d /d lam_m * [sigma_R *cos lam_m * d (poten-poten_R)/d lam_m]
!    
      logical,parameter :: mod_heelis = .false.   ! true == modified
      real ::  zigm_r(nmlonp1,nmlat),! sigma_r
     |  rim_r(nmlonp1,nmlat,2),      ! right-hand side 
     |  nsrim_r1(nmlon0,nmlat,10),   ! difference stencil for rim_r(1)
     |  nsrim_r2(nmlon0,nmlat,10),   ! difference stencil for rim_r(2)
                                     ! 10 dimension added for using nsstencil
     |  rhs_r(nmlonp1,nmlat)         ! right-hand side 
#if defined(INTERCOMM) || defined(CISMAH)
      real :: nsrhs(nmlonp1,nmlat)         ! CISM current
#endif
!
! Flag for J_rR contribution to dynamo
! J_rR = d / d phi_m [ f zigm_phiphi/ cos(lam_m)* d /d |lam_m| poten +
!        f zigm_philam   d /d phi_m poten ] + d / d |lam_m| [ f*zigm_lamphi*
!        d /d phi_m poten + f*zigm_lamlam*cos(lam_m)*d /d |lam_m| poten ]
!    
      logical,parameter :: J_rR = .false.   
!
! Flag for sigma_M (equivalent magnetic cond.) contribution to dynamo
! add to the conductances: Sigma^T_phiphi + Sigma_phiphi^M
!			   Sigma^T_philam + Sigma_H^M
!			   Sigma^T_lamphi - Sigma_H^M
! Sigma_phiphi^M: equivalent magnetospheric zonal Pedersen conductance
! Sigma_H^M: equivalent magnetospheric Hall conductance
! - acts on below the convection reversal boundary
! 
      logical,parameter :: eqMgCnd = .false.
      real,parameter :: eqMg_crit = 0.314159266 ! 18 deg colatitude where C_r,C_i start
      real ::  zigmc_r(nmlonp1,nmlat),  ! C_r/ Sigma_phiphi^M on tiegcm grid
     |   zigmc_i(nmlonp1,nmlat)         ! C_i/ Sigma_H^M on tiegcm grid
     
! values for equivalent magnet. conductances in 1deg increasing colatitude
! [ Peymirat, Richmond 1993: JGR, Vol98, pp.15467, fig. 4 ]    
      real :: c_r(15),c_i(15)
      data c_r/0.,1.4,3.,5.,8.9,14.9,21.9,    ! Sigma_phiphi^M [S]
     |   29.9,29.7,25.5,17.5,8.5,1.,.7,0./
      data c_i/25.8,25.8,23.7,20.9,17.3,      ! Sigma_H^M  [S]
     |   13.6,9.7,6.1,3.5,1.7,.8,0.,0.,0.,0./ 
!
! Flag for calculation of K_(q,lam)
!   icalkqlam = 1 -> will be calculated
!   icalkqlam = 0 -> will NOT be calculated
!
      integer,parameter :: icalkqlam = 0
      
      real ::
     |  je13d(nmlonp1,nmlat,-2:nlevp1),    ! 3d je1 magnetic
     |  je13d_geo(nlonp1,0:nlatp1,nlevp1), ! 3d je1 geographic
     |  je23d(nmlonp1,nmlat,-2:nlevp1),    ! 3d je2 magnetic
     |  je23d_geo(nlonp1,0:nlatp1,nlevp1)  ! 3d je2 geographic
      real ::
     |  je13d_diag(nlonp4,nlat,nlevp1), ! 3d je1 geographic for addfld
     |  je23d_diag(nlonp4,nlat,nlevp1)  ! 3d je1 geographic for addfld
!
! 3d arrays for calculation of height-integrated sheet-current density
! K_(q,lam) (if icalkqlam==1)
!
      real,dimension(nmlonp1,nmlat,-2:nlevp1) :: 
     |  sigma1m3d, sigma2m3d, adotv1m3d, adotv2m3d, ed13d, ed23d,
     |  je1oD_pg3d, ! plasma pressure, gravity (j_pg=true) a/cm^2
     |  je2oD_pg3d  ! plasma pressure, gravity (j_pg=true) a/cm^2
      real,dimension(nmlonp1,nmlat) :: 
     |  bmodm3d,a1a2m3d,adotam3d,sinim3d
!
! Equatorial quantities:
      real,dimension(nmlonp1,-2:nlev)   :: sigma1e,sigma2e,
     |  je1_pge,je2_pge ! j_pg = .true. J(mag.pressure, gravity)
      real,dimension(nmlonp1,-2:nlev,2) :: adotve
!
! For current density code (current.F) (icalkqlam==1):
      real :: nscoef(nmlon0,nmlat,10),nscrrt(nmlon0,nmlat)
!
! For dot products:
      real,parameter :: unitvm(nmlon)=1., unitv(nlon)=1.
!
! Electric potential from heelis or weimer:
      real,dimension(nmlonp1,nmlat0) :: pfrac  ! NH fraction of potential
      real,dimension(nmlonp1,nmlat)  :: phihm  ! potential in magnetic
      real,dimension(nlonp4,nlat,nlevp1) :: dynpot_diag
! 
! Work array (private to this module)
      real :: wkarray(-15:nmlon0+16,nmlat0)
      logical,parameter :: debug = .false. ! add print statements to stdout
      contains
!-----------------------------------------------------------------------
      subroutine prep_dynamo(tn,un,vn,w,z,barm,ped,hall,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Prepare geographic-grid fields for input to the dynamo, and gather them
! to the root task. This is executed by all tasks, and is called from 
! advance before the dynamo itself (which is executed by master task only).
! calculate vertical velocity at half pressure level
! 
      use cons_module,only: gask,grav
      integer :: lev0,lev1,lon0,lon1,lat0,lat1
      real :: fmin,fmax
!
! Input fields at task subdomains. These are from current (not updated)
! fields (itp):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,   ! neutral temperature (deg K)
     |  un,   ! neutral zonal velocity (cm/s)
     |  vn,   ! neutral meridional velocity (cm/s)
     |  w,    ! dimensionless vertical velocity (1/s)
     |  z,    ! geopotential height (cm)
     |  barm, ! mean molecular mass (g/mole)
     |  ped,  ! pedersen conductivity (lamdas.F) (S)
     |  hall  ! hall conductivity     (lamdas.F) (S)	  
!
! Local:
      integer :: k,i,lat,nk,lonbeg,lonend
!
      if (debug) write(6,"('Enter prep_dynamo..')")
      nk = lev1-lev0+1
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-1
!
! Define subdomain part of dynamo input fields (geographic):
! Also transform from (k,i,lat) to (i,lat,k). The latter is used
!   in the dynamo code, whereas the former is used in the rest
!   of the model.
! (In earlier versions, this was in lamdas, where 1-73 <= 3-75.
!  Here, 3-75 <= 3-75)
!
      do lat=lat0,lat1
        do i=lonbeg,lonend
          do k=lev0,lev1-1
            sigma_ped (i,lat,k) = ped (k,i,lat)
            sigma_hall(i,lat,k) = hall(k,i,lat)
            zpoten    (i,lat,k) = z   (k,i,lat)
            scheight  (i,lat,k) = gask*tn(k,i,lat)/
     |        (.5*(barm(k,i,lat)+barm(k+1,i,lat))*grav)
            unvel     (i,lat,k) = un  (k,i,lat)
            vnvel     (i,lat,k) = vn  (k,i,lat)
            wnvel     (i,lat,k) = .5*(w(k,i,lat)+w(k+1,i,lat))*
     |        scheight(i,lat,k)
            je1oD_full(i,lat,k) = je1oD_geo  (k,i,lat) ! [A/cm^2]
            je2oD_full(i,lat,k) = je2oD_geo  (k,i,lat) ! [A/cm^2]
          enddo ! k=lev0,lev1
          zpoten(i,lat,lev1) = z(lev1,i,lat)
        enddo ! i=lon0,lon1
!
!      call addfld('SIGMAPED','sigma-ped','S/m',
!    |   sigma_ped(lon0:lon1,lat,:),'lon',lon0,lon1,'lev',
!    |   lev0,lev1,lat)
!      call addfld('SIGMAHAL','sigma-hall','S/m',
!    |   sigma_hall(lon0:lon1,lat,:),'lon',lon0,lon1,'lev',
!    |   lev0,lev1,lat)
!      call addfld('ZPOTEN','geop.height','cm',
!    |   zpoten(lon0:lon1,lat,:),'lon',lon0,lon1,'ilev',
!    |   lev0,lev1,lat)
!      call addfld('SCHEIGHT','scl.Height','cm',
!    |   scheight(lon0:lon1,lat,:),'lon',lon0,lon1,'lev',
!    |   lev0,lev1,lat)
!      call addfld('UNVEL','zonal neutral wind','cm/s',
!    |   unvel(lon0:lon1,lat,:),'lon',lon0,lon1,'lev',
!    |   lev0,lev1,lat)
!      call addfld('VNVEL','merid. neutral wind','cm/s',
!    |   vnvel(lon0:lon1,lat,:),'lon',lon0,lon1,'lev',
!    |   lev0,lev1,lat)
!      call addfld('WNVEL','upward neutral wind','cm/s',
!    |   wnvel(lon0:lon1,lat,:),'lon',lon0,lon1,'lev',
!    |   lev0,lev1,lat)
!      call addfld('JE1oD_PG','Je1(p,g)','A/cm^2',
!    |   je1oD_full(lon0:lon1,lat,:),'lon',lon0,lon1,'lev',
!    |   lev0,lev1,lat)
!      call addfld('JE2oD_PG','Je1(p,g)','A/cm^2',
!    |   je2oD_full(lon0:lon1,lat,:),'lon',lon0,lon1,'lev',
!    |   lev0,lev1,lat)
!
      enddo ! lat=lat0,lat1
!
#ifdef MPI
!
! Gather dynamo input fields to the root task, defining module data 
! above at the global domain on the root task:
!
      call mp_dynamo_gather
#endif
      if (debug) write(6,"('Exit prep_dynamo..')")
      end subroutine prep_dynamo
!-----------------------------------------------------------------------
      subroutine dynamo
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
      use cons_module,only: dlatm,dlonm,pi_dyn,ylatm,rtd
      use magfield_module,only: im,jm,dim,djm
      use fields_module,only: dynpot,phim3d
      use input_module,only: potential_model
      use init_module,only: istep	! for sigma_R high latitudes
!
! Local:
      integer :: i,j,jj,jjj,j0,jntl,k,n,ncc,nmaglat,nmaglon,ier
      real :: sym,fac
      real :: cs(nmlat0)
      real,dimension(nlonp1,0:nlatp1) :: phih  ! potential in geographic
      real,dimension(nlonp4,nlat) :: phih_diag
!
! First transform fields to geomagnetic grid, and perform field line
! integrations:
      if (debug) write(6,"('Enter dynamo (call transf)..')")
      call transf
      if (debug) write(6,"('dynamo after call transf..')")
!
! values calculated in transf:
! zigm11 = Sigma_(phi phi)(0)= Sigma_(phi phi)(m) * cos lam_0 /cos lam_m * d lam_m /d lam_0
! zigm22 = Sigma_(lam lam)(0)= Sigma_(lam lam)(m) * cos lam_m /cos lam_0 * d lam_0 /d lam_m
! +-(zigm2-zigmc)= Sigma_(phi lam)(0) = Sigma_(phi lam)(m)
! -+(zigm2+zigmc)= Sigma_(lam phi)(0) = Sigma_(lam phi)(m)
! rim(1) = K_(m phi)^D(0)    = K_(m phi)^D(m) * d lam_m /d lam_0   
! rim(2) = +/-K_(m lam)^D(0) = +/-K_(m lam)^D(m) * cos lam_m /cos lam_0 
!
#if defined(INTERCOMM) || defined(CISMAH)
!
! Calculate global neutral winds generated current to transfer to CISM M-I couplier
!
      call cism_ucurrent
!
! Transform height-integrated parameters from APEX coordinates to geographic
! coordinates to be passed to the M-I coupled in the CISM
!
      call mag2geo(zigm1(1,1),gzigm1(3:nlon+3,0:nlatp1),im(1,0),jm(1,0),
     |     dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
      call mag2geo(zigm2(1,1),gzigm2(3:nlon+3,0:nlatp1),im(1,0),jm(1,0),
     |     dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
      call mag2geo(nsrhs(1,1),gnsrhs(3:nlon+3,0:nlatp1),im(1,0),jm(1,0),
     |     dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
#endif
!
! Calculate coefficient of PDE for each hemisphere if K_(q,lam) is 
!   calculated
!
      if (icalkqlam==1) call nosocoef
!
! Fold southern hemisphere onto northern hemisphere 
! -Value at the equator is also folded therefore at the equatorial boundary
!  factor of 1/2 introduced
! -added values in array index 49-97 from equator to poles
! -SH array (1-49) has original SH values
! -reverse sign of zigmc (to be compatible with Cicely's) 
! -sign of K_(m lam)^D in southern hemisphere is reversed, therefore 
!     K_(m lam)^N - (- K_(m lam)^S)
!
! zigm11 = Sigma_(phi phi)(0)^T
! zigm22 = Sigma_(lam lam)(0)^T
! zigmc  =-Sigma_c(0)^T 
! zigm2  = Sigma_h(0)^T 
! rim(1) = K_(m phi)^D(0)^T
! rim(2) = K_(m lam)^D(0)^T
!
      do j=1,nmlath
        do i=1,nmlonp1
          zigm11(i,nmlatp1-j) = (zigm11(i,nmlatp1-j)+zigm11(i,j))
          zigmc(i,nmlatp1-j)  =  zigmc(i,nmlatp1-j)+zigmc(i,j)
          zigmc(i,nmlat+1-j)  = -zigmc(i,nmlat+1-j)
!
          zigm2(i,nmlat+1-j)  =  zigm2(i,nmlat+1-J)+zigm2(i,j)
          zigm22(i,nmlat+1-j) = (zigm22(i,nmlat+1-J)+zigm22(i,j))
          rim(i,nmlatp1-j,1)  = (rim(i,nmlat+1-j,1)+rim(i,j,1))
          rim(i,nmlatp1-j,2)  = (rim(i,nmlat+1-j,2)+rim(i,j,2))  
        enddo ! i=1,nmlonp1
      enddo ! j=1,nmlath
!
! calculate RHS of PDE from rim(1) and rim(2)
!  R_0 * d lam_m/d lam_0 * 1/cos lam_0 [ d K_(m phi)^DT(0) / d phi  +
!  (d [ K_(m lam)^DT(0) * cos(lam_0)]/ d lam_0  ] 
!
      call rhspde
!
! Set index array nc array:
! nc pointes to the start of the coefficient array for each level
!
      nc(1) = nc0
      nc(2) = nc1
      nc(3) = nc2
      nc(4) = nc3
      nc(5) = nc4
      nc(6) = ncee
!
! Use pi_dyn from cons module rather than 4*atan(1) to avoid small
! differences generated by the atan in -lmass lib (-lmass was not
! used in earlier versions)
! set magnetic latitude cosine array: cos lam_0
!
      do j=1,nmlat0
        cs(j) = cos(pi_dyn/2.-(nmlat0-j)*dlatm)
      enddo ! j=1,nmlat0
!
! Set up difference coefficients. 
!
! mixed terms: factor 4 from 5-point diff. stencil
! "NH" array index 49:97 equator to pole
! zigmc  = Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm2  = Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm22 = Sigma_(lam lam)^T(0)*cos(lam_0)/d lam_0^2
! zigm11 = Sigma_(phi phi)^T(0)/ cos(lam_0) / d lon^2 )
!
! "SH" array index 49:1 equator to pole
! zigmc  = -Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm2  = -Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm22 = Sigma_(lam lam)^T(0)*cos(lam_0)/d lam_0^2
! zigm11 = Sigma_(phi phi)^T(0)/ cos(lam_0) / d lon^2 )
! 
      j0 = nmlat0-nmlath
      do j=1,nmlath       !  1,49 
        jj = nmlath+j-1   ! 49,97  added values ()^T 
        jjj = nmlath-j+1  ! 49,1
        do i=1,nmlonp1
          zigmc(i,jj)   = (zigmc(i,jj) +zigm2(i,jj))/(4.*dlatm*dlonm)
          zigm2(i,jj)   = zigmc(i,jj)-2.*zigm2(i,jj)/(4.*dlatm*dlonm)
          zigm22(i,jj)  = zigm22(i,jj)*cs(j0+j)/dlatm**2
          zigmc(i,jjj)  = -zigmc(i,jj)
          zigm2(i,jjj)  = -zigm2(i,jj)
          zigm22(i,jjj) = zigm22(i,jj)
        enddo ! i=1,nmlonp1
        if (j /= nmlath) then
          do i = 1,nmlonp1
            zigm11(i,jj) = zigm11(i,jj)/(cs(j0+j)*dlonm**2)
            zigm11(i,jjj) = zigm11(i,jj)
          enddo
        endif
      enddo ! j=1,nmlath
!
! Set zigm11 to zero at magnetic poles to avoid floating exception 
! (values at poles are not used):
!
      do i = 1,nmlonp1  
        zigm11(i,1)     = 0.0
        zigm11(i,nmlat) = 0.0
      enddo
!      
! for J_rR only contribution to rhs
!      
      if(J_rR) call calrhs_jrr(phihm)
!
! for modified heelis: 
!      
      if(mod_heelis) then
         if(istep == 1) call set_zigmar ! set sigma_R
	 call add_zigmar    		! diff. sigma_R; add to sigmas
	 call diff_rimr(phihm)	    	! differentiate rim_R; calculate rhs_r
	 call add_rimr	                ! add rhs_r of rhs
      endif
!
! for equivalent magnet. conductances: 
      if(eqMgCnd) then
         if(istep == 1) call set_cicr   ! interpolate C_i, C_r to grid
	 call add_cicr    		! diff. C_i,C_r; add to sigmas
      endif
!
! Save 2d (mlon,mlat) fields to secondary histories:
!
!     call addfld('ZIGM11','Sig^T_phiphi/cos lam0/d^2phi','S',
!    |  zigm11,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('ZIGM22','Sig^T_lamlam*cos lam0/d^2lam','S',
!    |  zigm22,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('ZIGMC','Sig^T_philam/(4dlam dphi)','S',
!    |  zigmc,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('ZIGM2','Sig^T_lamphi/(4dlam dphi)','S',
!    |  zigm2,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('RIM1','K_mphi^DT','A/m',
!    |  rim(:,:,1),'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('RIM2','K_mlam^DT','A/m',
!    |  rim(:,:,2),'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!
! Clear array cee for difference stencils at all levels:
      if (debug) write(6,"('dynamo before call clearcee)..')")
      call clearcee(cee,nmlon0,nmlat0)
!
! Calculate contribution to stencils from each PDE coefficient
! isolve = 0 -> original mud version 5.
! isolve = 1 -> muh hybrid solver (only as direct solver -- slow)
! isolve = 2 -> modified mudpack solver
!
! isolve = 2 -> modified mudpack solver (modified and unmodified coefficient
!               stencil)
      if (debug) write(6,"('dynamo call stencmd: isolve=',i3)") isolve
! for Sigma_(phi lam)(0)/( 4*d lam_0* d lon ) & 
!     Sigma_(lam phi)(0)/( 4*d lam_0* d lon )
! at the equator it's the value from SH -> sign change since NH is used
! sign changed back in case value is used somewhere else
!
      if (isolve==2) then
        cofum(:,:,:) = 0. ! init
!
! Sigma_(phi phi)(0)/( cos(lam_0)*(d lon)^2 )
        call stencmd(zigm11(1,nmlat0),nmlon0,nmlat0,cee,1)
!
! Sigma_(lam lam)(0)*cos(lam_0)/(d lam_0)^2
        call stencmd(zigm22(1,nmlat0),nmlon0,nmlat0,cee,4)
!
! Sigma_(phi lam)(0)/( 4*d lam_0* d lon )
	zigmc(:,nmlat0) = -zigmc(:,nmlat0)
        call stencmd(zigmc(1,nmlat0),nmlon0,nmlat0,cee,2)
	zigmc(:,nmlat0) = -zigmc(:,nmlat0)
!
! Sigma_(lam phi)(0)/( 4*d lam_0* d lon )
	zigm2(:,nmlat0) = -zigm2(:,nmlat0)
        call stencmd(zigm2(1,nmlat0),nmlon0,nmlat0,cee,3)
	zigm2(:,nmlat0) = -zigm2(:,nmlat0)
!
! isolve /= 2: original or hybrid solver (only modified stencil).
!
      else
!
! Sigma_(phi phi)(0)/( cos(lam_0)*(d lon)^2 )
        call stencil(zigm11(1,nmlat0),nmlon0,nmlat0,cee,1)
!
! Sigma_(lam lam)(0)*cos(lam_0)*/(d lam_0)^2
        call stencil(zigm22(1,nmlat0),nmlon0,nmlat0,cee,4)
!
! Sigma_(phi lam)(0)/( 4*d lam_0* d lon )
	zigmc(:,nmlat0) = -zigmc(:,nmlat0)
        call stencil(zigmc(1,nmlat0),nmlon0,nmlat0,cee,2)
	zigmc(:,nmlat0) = -zigmc(:,nmlat0)
!
! Sigma_(lam phi)(0)/( 4*d lam_0* d lon )
	zigm2(:,nmlat0) = -zigm2(:,nmlat0)
        call stencil(zigm2(1,nmlat0),nmlon0,nmlat0,cee,3)
	zigm2(:,nmlat0) = -zigm2(:,nmlat0)
!
      endif ! isolve
!
! Insert RHS in finest stencil (formerly sub rths):
      do j = 1,nmlat0
        jj = nmlath-nmlat0+j
        do i = 1,nmlon0
          c0(i,j,10) = rhs(i,jj)
        enddo ! i = 1,nmlon0
      enddo ! j = 1,nmlat0
!
! Set boundary condition at the pole:
      if (debug) write(6,"('dynamo call edges..')")
      call edges(c0,nmlon0,nmlat0)
      call edges(c1,nmlon1,nmlat1)
      call edges(c2,nmlon2,nmlat2)
      call edges(c3,nmlon3,nmlat3)
      call edges(c4,nmlon4,nmlat4)
      if (isolve==2) 
     |  call edges(cofum,nmlon0,nmlat0)
!
! Divide stencils by cos(lam_0) (not rhs already divided):
      if (debug) write(6,"('dynamo call divide..')")
      call divide(c0,nmlon0,nmlat0,nmlon0,nmlat0,cs,1)
      call divide(c1,nmlon1,nmlat1,nmlon0,nmlat0,cs,1)
      call divide(c2,nmlon2,nmlat2,nmlon0,nmlat0,cs,1)
      call divide(c3,nmlon3,nmlat3,nmlon0,nmlat0,cs,1)
      call divide(c4,nmlon4,nmlat4,nmlon0,nmlat0,cs,1)
      if (isolve==2) 
     |  call divide(cofum,nmlon0,nmlat0,nmlon0,nmlat0,cs,0)
!
! Set value of solution to 1. at pole:
      do i=1,nmlon0
        c0(i,nmlat0,10) = 1.
      enddo
#if defined(INTERCOMM) || defined(CISMAH)
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
!  high latitude -> phihm(1,nmlat0). The SH high lat potential will be added back later.
!  pfrac = fraction of dynamo in solution in the NH. = 1 low lat, = 0 hi lat
!    cons_module: crit(1)=15, crit(2)=30 deg colats, or hi-lat > 75 deg,
!      dynamo < 60 deg, and combination between 60-75 mag lat.
! The dynamo is symmetric about the magnetic equator, but the high latitude
!  is anti-symmetric in both hemispheres.  However, since Mudpack uses the
!  NH potential pattern, then the SH potential pattern must be added
!  back into the 2-D phim before the call threed, and before it is
!  transformed to geographic coordinates.
!
      if(.not.mod_heelis.and..not.J_rR) then   ! original heelis
        ncc = 1
        nmaglon = nmlon0
        nmaglat = nmlat0
        do n=1,5
          if (isolve==2) then
            call stenmd(nmaglon,nmaglat,cee(ncc),phihm(1,nmlat0),pfrac)
          else
            call stenmod(nmaglon,nmaglat,cee(ncc),phihm(1,nmlat0),pfrac)
          endif
          ncc = ncc+9*nmaglon*nmaglat
          if (n==1) ncc = ncc+nmaglon*nmaglat ! rhs is in 10th slot
          nmaglon = (nmaglon+1)/2
          nmaglat = (nmaglat+1)/2
        enddo ! n=1,5
      endif
!
      jntl = 0
!
! Call mudpack to solve PDE. Solution is returned in rim:
!     real,dimension(nmlonp1,nmlat,2) :: rim 
! Mudpack solvers:
!   isolve = 0  org. mud v5.      (modified stencils neq direct solution)
!   isolve = 1  muh hybrid solver (no convergence => only as direct solver)
!   isolve = 2  modified mud      (residual calculated with unmodified stencils
!                                  same solution as with direct solver, if same
!                                  coefficient matrix is used)
! ouput rim: electric potential from south to north pole
!       ier: 0 if converged, <0 not converged
!
      if (debug) write(6,"('dynamo call solver: isolve=',i3)") isolve
!
      ier = 0 
      if(isolve==0) then
        call mud(rim,jntl,isolve,ier)	 ! solver in mud.F
        if(ier < 0 ) then 	! not converged
	  write(6,*) 'muh: use direct solver'
	  call muh(rim,jntl) 		! solver in muh2cr.F
	endif
      elseif (isolve==1) then
        call muh(rim,jntl)        	! solver in muh2cr.F
      elseif (isolve==2) then
        call mudmod(rim,jntl,isolve,ier)     	! solver in mudmod.F
        if(ier < 0 ) then 	! not converged
	  write(6,*) 'muh: use direct solver'
	  call muh(rim,jntl) 		! solver in muh2cr.F
	endif
      else
        write(6,*) 'dynamo: solver type ',isolve,' not implemented.'
        call shutdown('isolve')
      endif
      if (debug) write(6,"('dynamo after solver..')")
!
! Copy output potential from rim to phim(nmlonp1,nmlat)
! for weimer - heelis (is hemispherical symmetric-> correction would be not
!        necessary):
!    Correct the SH potential for the anti-symmetric imposed NH high lat potential
!    electrodynamo equation was solved with NH Heelis potential -> call stenmod
!
! modified Heelis: no correction since symmetric
!
!  pfrac = fraction of dynamo in solution in the NH. = 1 low lat, = 0 hi lat
!    cons_module: crit(1)=15, crit(2)=30 deg colats, or hi-lat > 75 deg,
!      dynamo < 60 deg, and combination between 60-75 mag lat.
!  ie, need (1.-pfrac(i,nmlat0)) at phim(i,1), etc
! jn index for NH part of potential (nmlat down to ~nmlat0)
! jp index for NH pfrac (nmlat0 down to 1)
! 
      fac = 1.0
      if(mod_heelis) fac = 0.   ! modified heelis
      do j=1,nmlat0		! SH
        jn = nmlat  - j + 1
        jp = nmlat0 - j + 1
        do i=1,nmlonp1
          phim(i,j)=rim(i,j,1)+fac*(1.-pfrac(i,jp))*(phihm(i,j)-
     |        phihm(i,jn))
        enddo ! i=1,nmlonp1
      enddo ! j=1,nmlat0
      do j=nmlat0+1,nmlat	! NH
        do i=1,nmlonp1
          phim(i,j) = rim(i,j,1)
        enddo ! i=1,nmlonp1
      enddo ! j=1,nmlat
!
! Save single-level high latitude potential:
!      if (potential_model == 'WEIMER01') then
!    	call addfld('PHIHM2D','2D WEIMER01 ELECTRIC POTENTIAL',
!     |    'V', phihm,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!      elseif (potential_model == 'WEIMER05') then
!    	call addfld('PHIHM2D','2D WEIMER05 ELECTRIC POTENTIAL',
!     |    'V', phihm,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!      elseif (potential_model == 'HEELIS') then
!    	call addfld('PHIHM2D','2D HEELIS ELECTRIC POTENTIAL',
!     |    'V',phihm,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!      else
!    	call addfld('PHIHM2D','2D ELECTRIC POTENTIAL',
!     |    'V',phihm,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!      endif
!
      call mkdiag_PHIM2D('PHIM2D',phim,1,nmlonp1,1,nmlat)
!
!  am_02/02: calculate current for both hemisphere:
!                      R**2 * J_mr / dt0dts /rcos0s
       if(icalkqlam.eq.1) call nosocrrt
!
! Call threed to generate 3-d potential array in geomagnetic coordinates
!   from 2-d solver output phim, corrected for the SH potential.
!   phim3d(nmlonp1,nmlat,-2:nlevp1) is in fields.F.
!
      if (debug) write(6,"('dynamo call threed..')")
      call threed
      if (debug) write(6,"('dynamo after threed..')")
!
! Save 3d mag electric potential:
! phim3d(nmlonp1,nmlat,-2:nlevp1)  ! 3d electric potential magnetic
! Note nmlev = nlevp1, so 1:nmlev cannot be used in the addfld call.
!
!     do j=1,nmlat
!       call addfld('PHIM3D','El.Poten','V',
!    |    phim3d(:,j,1:nlevp1),'mlon',1,nmlonp1,'imlev',1,nlevp1,j)
!     enddo
!
! am_02/02: calculate  Je_1,  K_(q,phi), K_(q,lam)
      if(icalkqlam.eq.1) call nosocrdens
!
! Transform phim3d to geographic coordinates in dynpot (fields.F):
!   phim3d(nmlonp1,nmlat,-2:nlevp1)  ! 3d electric potential magnetic
!   dynpot(nlonp1,0:nlatp1,nlevp1),  ! 3d electric potential geographic
!
      do k=1,nlevp1
        call mag2geo(phim3d(1,1,k),dynpot(1,0,k),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
        call mag2geo(je13d(1,1,k),je13d_geo(1,0,k),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
        call mag2geo(je23d(1,1,k),je23d_geo(1,0,k),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
      enddo ! k=1,nlevp1
!
! Periodic point:
      do k=1,nlevp1
        do j=0,nlatp1
          dynpot(nlonp1,j,k)    = dynpot(1,j,k)
          je13d_geo(nlonp1,j,k) = je13d_geo(1,j,k)
          je23d_geo(nlonp1,j,k) = je23d_geo(1,j,k)
        enddo ! j=0,nlatp1
      enddo ! k=1,nlevp1
!
! Save electric potential on geographic coords to secondary history:
! dynpot(nlonp1,0:nlatp1,nlevp1),  ! 3d electric potential geographic
! dynpot_diag(nlonp4,nlat,nlevp1)  ! for addfld
!
!      do j=1,nlat
!        dynpot_diag(:,:,:) = spval
!        dynpot_diag(3:nlonp2,j,:) = dynpot(1:nlon,j,:)         ! 3,74 <= 1,72
!        dynpot_diag(1:2,j,:) = dynpot_diag(nlonp1:nlonp2,j,:)   ! 1,2 <= 73,74
!        dynpot_diag(nlonp4-1:nlonp4,j,:) = dynpot_diag(3:4,j,:) ! 75,76 <= 3,4
!        call addfld('DYNPOT','elec. potential (geo)','V',
!    |     dynpot_diag(:,j,:),'lon',1,nlonp4,'ilev',1,nlevp1,j)
! 
! je13d_diag(nlonp4,nlat,nlevp1), ! 3d je1 geographic for addfld
!        je13d_diag(:,:,:) = spval
!        je13d_diag(3:nlonp2,j,:) = je13d_geo(1:nlon,j,:)      ! 3,74 <= 1,72
!        je13d_diag(1:2,j,:) = je13d_diag(nlonp1:nlonp2,j,:)   ! 1,2 <= 73,74
!        je13d_diag(nlonp4-1:nlonp4,j,:) = je13d_diag(3:4,j,:) ! 75,76 <= 3,4
!        call addfld('JE13D','J_e1(p,g)','A/m^2',
!    |     je13d_diag(:,j,:),'lon',1,nlonp4,'lev',1,nlevp1,j)
!      enddo ! j=1,nlat
!
! Transform single-level heelis magnetic potential phihm to geographic 
! in phih:
      call mag2geo(phihm(1,1),phih(1,0),im(1,0),jm(1,0),dim(1,0),
     |  djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
!
! Periodic point:
      do j=0,nlatp1  
        phih(nlonp1,j) = phih(1,j)
      enddo ! j=0,nlatp1  
!
! real,dimension(nlonp4,nlat) :: phih_diag
!
!     do j=1,nlat
!       phih_diag(3:nlonp4-2,j) = phih(1:nlon,j)
!       phih_diag(1:2,j) = phih_diag(nlon+1:nlon+2,j)
!       phih_diag(nlon+3:nlon+4,j) = phih_diag(3:4,j)
!     enddo
!     call addfld('PHIH2D','elec. potential (geo)','V',
!    |    phih,'lon',1,nlonp4,'lat',1,nlat,0)
!
      if (debug) write(6,"('dynamo returning..')")
      end subroutine dynamo
!-----------------------------------------------------------------------
      subroutine transf
!
! -maps the geographic fields to the geomagnetic coordinate
!  system
! - call the fieldline integration routine to calculate the
!   conductances as well as the height integrated current density for the
!   right hand side
! - include the boundary condition at the equator
! - transformation from lam_m to lam_0 (regular grid)
!
      use cons_module,only: h0,r0,dt1dts,rcos0s
      use magfield_module,only: 
     |  alatm,    ! (nlonp1,0:nlatp1)
     |  zb,       ! (nlonp1,0:nlatp1)
     |  bmod,     ! (nlonp1,0:nlatp1)
     |  ig,       ! (nmlonp1,nmlat) geog lon at each geomag grid point
     |  jg,       ! (nmlonp1,nmlat) geog lat at each geomag grid point
     |  wt        ! (4,nmlonp1,nmlat) interpolation weights
      use getapex_module,only:           ! (see sub apxparm, apex.F)
     |  be3arr,   ! (nlonp1,nlat), 
     |  dvec,     ! (nlonp1,nlat,3,2)
     |  dddarr    ! (nlonp1,nlat)
!
! Local:
      integer :: i,ii,k,kk,j,jj,lat,n
      real :: z0,sinalat,r0or,rat,corfac
      real,dimension(nlonp1,nlat) :: clm2
      real,dimension(nlonp1,0:nlatp1) :: a1dta2,sini,be3
!
! Fields to be transformed to geomagnetic space (formerly in transmag.h):
! (these will be input to the geographic to magnetic transformation). 
      real,dimension(nlonp1,0:nlatp1,-2:nlev) ::
     |  ssigma1,   ! pedersen conductivity
     |  ssigma2,   ! hall conductivity
     |  je1oDtrn,  ! J_ed1/D [A/cm^2]
     |  je2oDtrn   ! J_ed2/D [A/cm^2]
      real,dimension(nlonp1,0:nlatp1,-2:nlevp1) :: zz
      real,dimension(nlonp1,0:nlatp1,-2:nlev,2) :: adotv
      real,dimension(nlonp1,0:nlatp1,2)         :: adota
!
! Diagnostics for plotting:
      real,dimension(nmlonp1,nlevp1) ::
     |  zigm11_plt, zigmc_plt, zigm2_plt, zigm22_plt, rim1_plt, rim2_plt
      real,dimension(nlonp4,nlevp1) ::
     |  ped_plt, hall_plt, z_plt, h_plt, u_plt, v_plt, w_plt,
     |  sigma1_plt,sigma2_plt,zz_plt,adotv1_plt,adotv2_plt
!
! External:
      real,external :: sddot ! in util.F
!
! Set constants:
!   rl1,rl2 are rates at which sigma1 and sigma2 decay with height
!     below bottom of model.
!   z0 is lowest level for start of field line integration set in h0 
!     (h0 in cons module)
!
      real,parameter :: 
     |  rl1 = 5.e5,
     |  rl2 = 3.e5
!
! Pack inputs from 3->nlon+3 to 1->nlon+1, as in tgcm15 (end of sub lamdas):
      do lat=1,nlat
        do i=1,nlon+1
          ii = i+2
          sigma_ped (i,lat,:) = sigma_ped (ii,lat,:)
          sigma_hall(i,lat,:) = sigma_hall(ii,lat,:)
          zpoten    (i,lat,:) = zpoten    (ii,lat,:)
          scheight  (i,lat,:) = scheight  (ii,lat,:)
          unvel     (i,lat,:) = unvel     (ii,lat,:)
          vnvel     (i,lat,:) = vnvel     (ii,lat,:)
          wnvel     (i,lat,:) = wnvel     (ii,lat,:)
          je1oD_full(i,lat,:) = je1oD_full(ii,lat,:)
          je2oD_full(i,lat,:) = je2oD_full(ii,lat,:)
        enddo ! i=1,nlonp4-3
      enddo ! lat=1,nlat
!
      z0 = h0  ! 90 km 
!
      do j=1,nlat
        do i=1,nlonp1
          sinalat = sin(alatm(i,j))         ! sin(lam)
          clm2(i,j) = 1. - sinalat*sinalat  ! cos^2(lam)
          sini(i,j) = zb(i,j)/bmod(i,j)     ! sin(I_m)
          be3(i,j)  = 1.e-9*be3arr(i,j)     ! be3 is in T (be3arr in nT)

        enddo ! i=1,nlonp1
      enddo ! j=1,nlat
!
! Calculate quantities to be transformed to geomagnetic space:
      do k=1,nlev        
        do j=1,nlat
          do i=1,nlonp1
            ssigma1(i,j,k) = sigma_ped (i,j,k)
            ssigma2(i,j,k) = sigma_hall(i,j,k)
            zz(i,j,k) = zpoten(i,j,k)
!           ww(i,j,k) = wnvel (i,j,k)
            je1oDtrn(i,j,k) = je1oD_full(i,j,k)  ! [A/cm^2]
            je2oDtrn(i,j,k) = je2oD_full(i,j,k)  ! [A/cm^2]
          enddo ! i=1,nlonp1
        enddo ! j=1,nlat
      enddo ! k=1,nlev        
      do j=1,nlat
        do i=1,nlonp1
          zz(i,j,nlevp1) = zpoten(i,j,nlevp1)
        enddo ! i=1,nlonp1
      enddo ! j=1,nlat
!
! Extend fields down to 90 km inserting 3 extra levels.
! 90 km is assumed to be the lower boundary of the ionosphere. 
! Set three equally spaced levels for Z, take U, V, and W
!   to be constant, and extrapolate sigmas exponentially.
!
      do k=0,-2,-1
        do j=1,nlat
          do i=1,nlonp1
            zz(i,j,k) = z0+float(k+2)*(zz(i,j,1)-z0)/3.
            ssigma1(i,j,k) = ssigma1(i,j,1)*exp((zz(i,j,k)+zz(i,j,k+1)-
     |        zz(i,j,1)-zz(i,j,2))/(2.*rl1))
            ssigma2(i,j,k) = ssigma2(i,j,1)*exp((zz(i,j,k)+zz(i,j,k+1)-
     |        zz(i,j,1)-zz(i,j,2))/(2.*rl2))
            je1oDtrn(i,j,k) = je1oDtrn(i,j,1)*exp((zz(i,j,k)+zz
     |        (i,j,k+1)-zz(i,j,1)-zz(i,j,2))/(2.*rl2))
            je2oDtrn(i,j,k) = je2oDtrn(i,j,1)*exp((zz(i,j,k)+zz
     |        (i,j,k+1)-zz(i,j,1)-zz(i,j,2))/(2.*rl2))
          enddo ! i=1,nlonp1
        enddo ! j=1,nlat
      enddo ! k=0,-2,-1
!
! Calculation of adotv = ue1,ue2 (m/s):
      do k=-2,nlev
        kk = k
        if (kk < 1) kk = 1
        do j=1,nlat
          do i=1,nlonp1
!
! am 6/04: d_1,D_2,D are referenced to h_0 (90 km) s. apex.F
! in the following they are scaled from h_0 to the actual altitude
! d_1 = (R_0/R)^1.5
            r0or = r0/(r0 + .5*(zz(i,j,k)+zz(i,j,k+1)) - z0)    
            rat = 1.e-2*r0or**1.5 ! 1/100 convertion in cm
!
! A_1 dot V = fac( d_1(1) u + d_1(2) v + d_1(3) w
            adotv(i,j,k,1) = rat*(                              
     |        dvec(i,j,1,1)*unvel(i,j,kk) +                         
     |        dvec(i,j,2,1)*vnvel(i,j,kk) +
     |        dvec(i,j,3,1)*wnvel(i,j,kk))
!
! Note: clm2 is being used here to represent the squared cosine of the
!  quasi-dipole latitude, not of the M(90) latitude, since the wind
!  values are aligned vertically, not along the field line.
!
            rat = rat*sqrt((4.-3.*clm2(i,j))/(4.-3.*r0or*clm2(i,j)))
!
! A_2 dot V = fac( d_2(1) u + d_2(2) v + d_2(3) w
            adotv(i,j,k,2) = rat*(                              
     |        dvec(i,j,1,2)*unvel(i,j,kk) + 
     |        dvec(i,j,2,2)*vnvel(i,j,kk) +
     |        dvec(i,j,3,2)*wnvel(i,j,kk))
          enddo ! i=1,nlonp1
        enddo ! j=1,nlat
      enddo ! k=-2,nlev
!
! Calculation of adota(n) = d(n)**2/D
!                a1dta2   = (d(1) dot d(2)) /D
      do j=1,nlat
        do i=1,nlonp1 
          adota(i,j,1) = (dvec(i,j,1,1)**2 + dvec(i,j,2,1)**2 +
     |                    dvec(i,j,3,1)**2)/dddarr(i,j)
          adota(i,j,2) = (dvec(i,j,1,2)**2 + dvec(i,j,2,2)**2 +
     |                    dvec(i,j,3,2)**2)/dddarr(i,j)
          a1dta2(i,j)  = (dvec(i,j,1,1)*dvec(i,j,1,2) +
     |                    dvec(i,j,2,1)*dvec(i,j,2,2) +
     |                    dvec(i,j,3,1)*dvec(i,j,3,2))/dddarr(i,j)
        enddo ! i=1,nlonp1 
      enddo ! j=1,nlat
!
! Values at poles:
      do k=-2,nlev
        zz(1,0,k) = 
     |    (9.*sddot(nlon,unitv,zz(1,1,k))-
     |        sddot(nlon,unitv,zz(1,2,k)))/(8.*float(nlon))
        zz(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,zz(1,nlat  ,k))-
     |        sddot(nlon,unitv,zz(1,nlat-1,k)))/(8.*float(nlon))

        ssigma1(1,0,k) = 
     |    (9.*sddot(nlon,unitv,ssigma1(1,1,k))-
     |        sddot(nlon,unitv,ssigma1(1,2,k)))/(8.*float(nlon))
        ssigma1(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,ssigma1(1,nlat  ,k))-
     |        sddot(nlon,unitv,ssigma1(1,nlat-1,k)))/(8.*float(nlon))

        ssigma2(1,0,k) = 
     |    (9.*sddot(nlon,unitv,ssigma2(1,1,k))-
     |        sddot(nlon,unitv,ssigma2(1,2,k)))/(8.*float(nlon))
        ssigma2(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,ssigma2(1,nlat  ,k))-
     |        sddot(nlon,unitv,ssigma2(1,nlat-1,k)))/(8.*float(nlon))
!
! Calculate velocity distribution also at the poles.
! Used later for calculating K_{q lam}
!
        adotv(1,0,k,1) =  
     |     (9.*sddot(nlon,unitv,adotv(1,1,k,1)) - 
     |      sddot(nlon,unitv,adotv(1,2,k,1)))/(8.*float(nlon))
        adotv(1,nlatp1,k,1) =   
     |     (9.*sddot(nlon,unitv,adotv(1,nlat  ,k,1)) -
     |      sddot(nlon,unitv,adotv(1,nlat-1,k,1)))/(8.*float(nlon))
        adotv(1,0,k,2) =  
     |     (9.*sddot(nlon,unitv,adotv(1,1,k,2)) - 
     |      sddot(nlon,unitv,adotv(1,2,k,2)))/(8.*float(nlon))
        adotv(1,nlatp1,k,2) =  
     |     (9.*sddot(nlon,unitv,adotv(1,nlat  ,k,2)) -
     |      sddot(nlon,unitv,adotv(1,nlat-1,k,2)))/(8.*float(nlon))
!
! Calculate J(grad p, g) also at the poles.
! Used later for calculating K_{q lam}     
        je1oDtrn(1,0,k) = 
     |    (9.*sddot(nlon,unitv,je1oDtrn(1,1,k))-
     |        sddot(nlon,unitv,je1oDtrn(1,2,k)))/(8.*float(nlon))
        je1oDtrn(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,je1oDtrn(1,nlat  ,k))-
     |        sddot(nlon,unitv,je1oDtrn(1,nlat-1,k)))/(8.*float(nlon))   
        je2oDtrn(1,0,k) = 
     |    (9.*sddot(nlon,unitv,je2oDtrn(1,1,k))-
     |        sddot(nlon,unitv,je2oDtrn(1,2,k)))/(8.*float(nlon))
        je2oDtrn(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,je2oDtrn(1,nlat  ,k))-
     |        sddot(nlon,unitv,je2oDtrn(1,nlat-1,k)))/(8.*float(nlon))

!
! Extend in longitude:
        do i = 2,nlon
          zz(i,0,k)             = zz(1,0,k)
          zz(i,nlatp1,k)        = zz(1,nlatp1,k)
          ssigma1(i,0,k)        = ssigma1(1,0,k)
          ssigma1(i,nlatp1,k)   = ssigma1(1,nlatp1,k)
          ssigma2(i,0,k)        = ssigma2(1,0,k)
          ssigma2(i,nlatp1,k)   = ssigma2(1,nlatp1,k)
	  adotv(i,nlatp1,k,1)   = adotv(1,nlatp1,k,1)
	  adotv(i,nlatp1,k,2)   = adotv(1,nlatp1,k,2)
	  adotv(i,0,k,1)        = adotv(1,0,k,1)
	  adotv(i,0,k,2)        = adotv(1,0,k,2)
          je1oDtrn(i,0,k)	= je1oDtrn(1,0,k)
          je1oDtrn(i,nlatp1,k)  = je1oDtrn(1,nlatp1,k)
          je2oDtrn(i,0,k)	= je2oDtrn(1,0,k)
          je2oDtrn(i,nlatp1,k)  = je2oDtrn(1,nlatp1,k)
        enddo ! i = 2,nlon
      enddo ! k=-2,nlev
!
! Values at the poles:
      zz(1,0,nlevp1)= (9.*sddot(nlon,unitv,zz(1,1,nlevp1))-
     |                  sddot(nlon,unitv,zz(1,2,nlevp1)))/
     |                  (8.*float(nlon))
      zz(1,nlatp1,nlevp1) =  (9.*
     |                  sddot(nlon,unitv,zz(1,nlat,nlevp1))-
     |                  sddot(nlon,unitv,zz(1,nlat-1,nlevp1)))/
     |                  (8.*float(nlon))
      adota(1,0,1) = (9.*sddot(nlon,unitv,adota(1,1,1))-
     |                   sddot(nlon,unitv,adota(1,2,1)))/
     |                   (8.*float(nlon))
      adota(1,nlatp1,1) =  (9.*
     |                  sddot(nlon,unitv,adota(1,nlat,1))-
     |                  sddot(nlon,unitv,adota(1,nlat-1,1)))/
     |                  (8.*float(nlon))
      adota(1,0,2) = (9.*sddot(nlon,unitv,adota(1,1,2))-
     |                   sddot(nlon,unitv,adota(1,2,2)))/
     |                   (8.*float(nlon))
      adota(1,nlatp1,2) =  (9.*
     |                  sddot(nlon,unitv,adota(1,nlat,2))-
     |                  sddot(nlon,unitv,adota(1,nlat-1,2)))/
     |                  (8.*float(nlon))
      a1dta2(1,0)       = (9.*sddot(nlon,unitv,a1dta2(1,1))-
     |                   sddot(nlon,unitv,a1dta2(1,2)))/
     |                   (8.*float(nlon))
      a1dta2(1,nlatp1)  = (9.*
     |                  sddot(nlon,unitv,a1dta2(1,nlat))-
     |                  sddot(nlon,unitv,a1dta2(1,nlat-1)))/
     |                  (8.*float(nlon))
      sini(1,0)       = (9.*sddot(nlon,unitv,sini(1,1))-
     |                   sddot(nlon,unitv,sini(1,2)))/
     |                   (8.*float(nlon))
      sini(1,nlatp1)  = (9.*
     |                   sddot(nlon,unitv,sini(1,nlat))-
     |                   sddot(nlon,unitv,sini(1,nlat-1)))/
     |                   (8.*float(nlon))
      be3(1,0)        = (9.*sddot(nlon,unitv,be3(1,1))-
     |                   sddot(nlon,unitv,be3(1,2)))/
     |                  (8.*float(nlon))
      be3(1,nlatp1)   = (9.*
     |                   sddot(nlon,unitv,be3(1,nlat))-
     |                   sddot(nlon,unitv,be3(1,nlat-1)))/
     |                   (8.*float(nlon))
!
! Extend in longitude:
      do i = 2,nlon
        zz(i,0,nlevp1)     = zz(1,0,nlevp1)
        zz(i,nlatp1,nlevp1)= zz(1,nlatp1,nlevp1)
        adota(i,0,1)       = adota(1,0,1)
        adota(i,nlatp1,1)  = adota(1,nlatp1,1)
        adota(i,0,2)       = adota(1,0,2)
        adota(i,nlatp1,2)  = adota(1,nlatp1,2)
        a1dta2(i,0)        = a1dta2(1,0)
        a1dta2(i,nlatp1)   = a1dta2(1,nlatp1)
        sini(i,0)          = sini(1,0)
        sini(i,nlatp1)     = sini(1,nlatp1)
        be3(i,0)           = be3(1,0)
        be3(i,nlatp1)      = be3(1,nlatp1)
      enddo ! i = 2,nlon
!
! Periodic points:
      do j = 0,nlatp1
        do k = -2,nlev
           zz(nlonp1,j,k)      = zz(1,j,k)
           ssigma1(nlonp1,j,k) = ssigma1(1,j,k)
           ssigma2(nlonp1,j,k) = ssigma2(1,j,k)
           adotv(nlonp1,j,k,1) = adotv(1,j,k,1)
           adotv(nlonp1,j,k,2) = adotv(1,j,k,2)
           je1oDtrn(nlonp1,j,k)= je1oDtrn(1,j,k)
           je2oDtrn(nlonp1,j,k)= je2oDtrn(1,j,k)
        enddo ! k = -2,nlev
        zz(nlonp1,j,nlevp1)= zz(1,j,nlevp1)
        adota(nlonp1,j,1) = adota(1,j,1)
        adota(nlonp1,j,2) = adota(1,j,2)
        a1dta2(nlonp1,j)  = a1dta2(1,j)
        sini(nlonp1,j)    = sini(1,j)
        be3(nlonp1,j)     = be3(1,j)
      enddo ! j = 0,nlatp1
!
! Transform needed fields to geomagnetic coordinate system, one latitude
! at a time. 
!
! Equatorial values are used in sub integrals for interpolating 
!   conductivities and winds along the magnetic field lines:
!
! subroutine geo2mag(fmag,fgeo,long,latg,wght,nlonp1_geo,nlonp1_mag,
!|  nlon_mag,nlat_mag,lat)
!
      jj = nmlat/2+1 
      do k=-2,nlev
        call geo2mag(sigma1e(1,k),ssigma1(1,0,k),ig,jg,wt,nlonp1,
     |    nmlonp1,nmlon,nmlat,jj)
        sigma1e(nmlonp1,k) = sigma1e(1,k) ! periodic point
!
        call geo2mag(sigma2e(1,k),ssigma2(1,0,k),ig,jg,wt,nlonp1,
     |    nmlonp1,nmlon,nmlat,jj)
        sigma2e(nmlonp1,k) = sigma2e(1,k) ! periodic point
!
        call geo2mag(adotve(1,k,1),adotv(1,0,k,1),ig,jg,wt,nlonp1,
     |    nmlonp1,nmlon,nmlat,jj)
        adotve(nmlonp1,k,1) = adotve(1,k,1) ! periodic point
!
        call geo2mag(adotve(1,k,2),adotv(1,0,k,2),ig,jg,wt,nlonp1,
     |    nmlonp1,nmlon,nmlat,jj)
        adotve(nmlonp1,k,2) = adotve(1,k,2) ! periodic point
!	
        if(j_pg) then               ! for magnetic perturbation
          call geo2mag(je1_pge(1,k),je1oDtrn(1,0,k),ig,jg,wt,
     |      nlonp1,nmlonp1,nmlon,nmlat,jj)
          call geo2mag(je2_pge(1,k),je2oDtrn(1,0,k),ig,jg,wt,
     |      nlonp1,nmlonp1,nmlon,nmlat,jj)
          je1_pge(nmlonp1,k) = je1_pge(1,k)  ! periodic point
          je2_pge(nmlonp1,k) = je2_pge(1,k)
        endif
      enddo ! k=-2,nlev
!
! Magnetic latitude loop for field line integration: 
!
      maglat_loop: do j=2,nmlat-1
!      
! mapping from geographic grid to geomagnetic grid      
! height varying: sigma_ped, sigma_hall, u_e1, u_e2, z
! indepent of height: d_1^2/D, d_2^2/D, d_1*d_2/D, sin I_m, B_e3
!
        do k=-2,nlev
          call geo2mag(sigma_pedm(1,k),ssigma1(1,0,k),ig,jg,wt,nlonp1,
     |      nmlonp1,nmlon,nmlat,j)
          call geo2mag(sigma_hallm(1,k),ssigma2(1,0,k),ig,jg,wt,nlonp1,
     |      nmlonp1,nmlon,nmlat,j)
          call geo2mag(zpotenm(1,k),zz(1,0,k),ig,jg,wt,nlonp1,
     |      nmlonp1,nmlon,nmlat,j)
          do n=1,2
            call geo2mag(adotvm(1,k,n),adotv(1,0,k,n),ig,jg,wt,nlonp1,
     |        nmlonp1,nmlon,nmlat,j)
          enddo ! n=1,2
          sigma_pedm(nmlonp1,k)  = sigma_pedm(1,k) ! periodic point
          sigma_hallm(nmlonp1,k) = sigma_hallm(1,k)
!	
!  current due to plasma gradient and gravity je1_pg,je2_pg [A/cm^2] 
          if(j_pg) then
            call geo2mag(je1_pg(1,k),je1oDtrn(1,0,k),ig,jg,wt,
     |        nlonp1,nmlonp1,nmlon,nmlat,j)
            call geo2mag(je2_pg(1,k),je2oDtrn(1,0,k),ig,jg,wt,
     |        nlonp1,nmlonp1,nmlon,nmlat,j)
            je1_pg(nmlonp1,k) = je1_pg(1,k)  ! periodic point
            je2_pg(nmlonp1,k) = je2_pg(1,k)
          endif
!	  
        enddo ! k=-2,nlev
!
        k=nlevp1
        call geo2mag(zpotenm(1,k),zz(1,0,k),ig,jg,wt,nlonp1,
     |    nmlonp1,nmlon,nmlat,j)
!
        do n=1,2
          call geo2mag(adotam(1,n),adota(1,0,n),ig,jg,wt,nlonp1,
     |      nmlonp1,nmlon,nmlat,j)
        enddo ! n=1,2
!
        call geo2mag(a1a2m ,a1dta2,ig,jg,wt,nlonp1,nmlonp1,nmlon,
     |    nmlat,j)
        call geo2mag(siniam,sini  ,ig,jg,wt,nlonp1,nmlonp1,nmlon,
     |    nmlat,j)
        call geo2mag(bmodm ,be3   ,ig,jg,wt,nlonp1,nmlonp1,nmlon,
     |    nmlat,j)
!	
! real,dimension(nmlonp1,-2:nlev) :: ! J(mag.pres,gravity) mag. grid
!  je1_pg,	 ! J(gravity)
!  je2_pg	 ! J(mag.pres)
!
!       call addfld('je1_pg' ,'je1_pg','A/cm^2',je1_pg(:,1:nmlev),
!    |    'mlon',1,nmlonp1,'mlev',1,nmlev,j)
!       call addfld('je2_pg' ,'je2_pg','A/cm^2',je2_pg(:,1:nmlev),
!    |    'mlon',1,nmlonp1,'mlev',1,nmlev,j)
!       call addfld('SIGMA1M','sig_ped','S/m',sigma_pedm(:,1:nmlev),
!    |    'mlon',1,nmlonp1,'mlev',1,nmlev,j)
!       call addfld('SIGMA2M','sig_hall','S/m',sigma_hallm(:,1:nmlev),
!    |    'mlon',1,nmlonp1,'mlev',1,nmlev,j)
!
! Calculation of k_(q,lam)
!
        if (icalkqlam==1) then
          do i=1,nmlon
            bmodm3d(i,j)   = bmodm(i)
            a1a2m3d(i,j)   = a1a2m(i)
            adotam3d(i,j)  = adotam(i,1)
            sinim3d(i,j)   = siniam(i)
          enddo ! i=1,nmlon
!
! Save geopotential in 3d array (both are module data above):
          do k=-2,nlev
            do i=1,nmlon
              if (zpotenm(i,k) < z0) zpotenm(i,k) = z0 
              zpotenm3d(i,j,k) = zpotenm(i,k)
!
! Store 3d arrays in magnetic coordinates for postprocessing of the
! current density k_(q,lam)
!
              sigma1m3d(i,j,k) = sigma_pedm(i,k)
              sigma2m3d(i,j,k) = sigma_hallm(i,k)
              adotv1m3d(i,j,k) = adotvm(i,k,1)
              adotv2m3d(i,j,k) = adotvm(i,k,2)
!	                    ! gravity + plasma pressure [A/cm^2]
              je1oD_pg3d(i,j,k) = je1_pg(i,k)  
              je2oD_pg3d(i,j,k) = je2_pg(i,k)  
            enddo ! i=1,nmlon
!
! Periodic point:
            bmodm3d(nmlonp1,j)  = bmodm3d(1,j)
            a1a2m3d(nmlonp1,j)  = a1a2m3d(1,j)
            adotam3d(nmlonp1,j) = adotam3d(1,j)
            sinim3d(nmlonp1,j)  = sinim3d(1,j)
            sigma1m3d(nmlonp1,j,:) = sigma1m3d(1,j,:)
            sigma2m3d(nmlonp1,j,:) = sigma2m3d(1,j,:)
            adotv1m3d(nmlonp1,j,:) = adotv1m3d(1,j,:)
            adotv2m3d(nmlonp1,j,:) = adotv2m3d(1,j,:)
	                     ! gravity + mag.pressure [A/cm^2]  
            je1oD_pg3d(nmlonp1,j,:)= je1oD_pg3d(1,j,:)  
            je2oD_pg3d(nmlonp1,j,:)= je2oD_pg3d(1,j,:) 
          enddo ! k=-2,nlevp1
!
! Upper boundary:
          do i=1,nmlon
            if (zpotenm(i,nlevp1) < z0) zpotenm(i,nlevp1) = z0
            zpotenm3d(i,j,nlevp1) = zpotenm(i,nlevp1)
            sigma1m3d(i,j,nlevp1) = sigma_pedm(i,nlev)
            sigma2m3d(i,j,nlevp1) = sigma_hallm(i,nlev)
            adotv1m3d(i,j,nlevp1) = adotvm(i,nlev,1)
            adotv2m3d(i,j,nlevp1) = adotvm(i,nlev,2)
	                    ! gravity + mag.pressure [A/cm^2]  
            je1oD_pg3d(i,j,nlevp1)  = je1_pg(i,nlev)
            je2oD_pg3d(i,j,nlevp1)  = je2_pg(i,nlev)
          enddo ! i=1,nmlon
!
! Without calculation of k_(q,lam)
        else  ! icalkqlam /= 1
          do k=-2,nlev
            do i=1,nmlon  
              if (zpotenm(i,k) < z0) zpotenm(i,k) = z0
              zpotenm3d(i,j,k) = zpotenm(i,k)
            enddo ! i=1,nmlon  
          enddo ! k=-2,nlev
          do i=1,nmlon  
            if (zpotenm(i,nlevp1) < z0) zpotenm(i,nlevp1) = z0
            zpotenm3d(i,j,nlevp1) = zpotenm(i,nlevp1)
          enddo ! i=1,nmlon  
        endif ! icalkqlam
!
! Skip the equator:
        if (j==nmlat/2+1) cycle maglat_loop
!
! Call fieldline_integrals to perform field line integrations and 
!   evaluate PDE coefficients and RHS for current latitude j:
!
        call fieldline_integrals(j)
      enddo maglat_loop ! j=2,nmlat-1 Main magnetic latitude loop
!
! Values for the poles j=1 and j=nmlat for the 3d arrays.
!
      if (icalkqlam==1) then
        do jj=1,2
          j = jj
          if (jj==2) j = nmlat
          call geo2mag(a1a2m,a1dta2,ig,jg,wt,nlonp1,nmlonp1,
     |      nmlon,nmlat,j)
          call geo2mag(siniam,sini,ig,jg,wt,nlonp1,nmlonp1,
     |      nmlon,nmlat,j)
          call geo2mag(bmodm,be3,ig,jg,wt,nlonp1,nmlonp1,
     |      nmlon,nmlat,j)
          call geo2mag(adotam(1,1),adota(1,0,1),ig,jg,wt,
     |      nlonp1,nmlonp1,nmlon,nmlat,j)
          do k=-2,nlev
            call geo2mag(sigma_pedm(1,k),ssigma1(1,0,k),ig,jg,wt,
     |        nlonp1,nmlonp1,nmlon,nmlat,j)
            call geo2mag(sigma_hallm(1,k),ssigma2(1,0,k),ig,jg,wt,
     |        nlonp1,nmlonp1,nmlon,nmlat,j)
            do n=1,2
              call geo2mag(adotvm(1,k,n),adotv(1,0,k,n),ig,jg,wt,
     |          nlonp1,nmlonp1,nmlon,nmlat,j)
            enddo ! n=1,2
!
! Periodic points:
            sigma_pedm(nmlonp1,k) = sigma_pedm(1,k)
            sigma_hallm(nmlonp1,k) = sigma_hallm(1,k)
!	  
            if(j_pg) then
              call geo2mag(je1_pg(1,k),je1oDtrn(1,0,k),ig,jg,wt,
     |      	nlonp1,nmlonp1,nmlon,nmlat,j)
              call geo2mag(je2_pg(1,k),je2oDtrn(1,0,k),ig,jg,wt,
     |      	nlonp1,nmlonp1,nmlon,nmlat,j)
              je1_pg(nmlonp1,k) = je1_pg(1,k)  ! periodic point
              je2_pg(nmlonp1,k) = je2_pg(1,k)
            endif
!	  
          enddo ! k=-2,nlev
!
          do i = 1,nmlon
            bmodm3d(i,j)   = bmodm(i)
            a1a2m3d(i,j)   = a1a2m(i)
            adotam3d(i,j)  = adotam(i,1)
            sinim3d(i,j)   = siniam(i)
            do k = -2,nlev
              sigma1m3d(i,j,k) = sigma_pedm(i,k)
              sigma2m3d(i,j,k) = sigma_hallm(i,k)
              adotv1m3d(i,j,k) = adotvm(i,k,1)
              adotv2m3d(i,j,k) = adotvm(i,k,2)
!	      
              je1oD_pg3d(i,j,k)  = je1_pg(i,k)  ! gravity + mag.pressure
              je2oD_pg3d(i,j,k)  = je2_pg(i,k)  ! gravity + mag.pressure
            enddo
            sigma1m3d(i,j,nlevp1) = sigma_pedm(i,nlev)
            sigma2m3d(i,j,nlevp1) = sigma_hallm(i,nlev)
            adotv1m3d(i,j,nlevp1) = adotvm(i,nlev,1)
            adotv2m3d(i,j,nlevp1) = adotvm(i,nlev,2)
          enddo ! i = 1,nmlon
!
! Periodic points:
          bmodm3d(nmlonp1,j)  = bmodm3d(1,j)
          a1a2m3d(nmlonp1,j)  = a1a2m3d(1,j)
          adotam3d(nmlonp1,j) = adotam3d(1,j)
          sinim3d(nmlonp1,j)  = sinim3d(1,j)
!
          sigma1m3d(nmlonp1,j,:) = sigma1m3d(1,j,:)
          sigma2m3d(nmlonp1,j,:) = sigma2m3d(1,j,:)
          adotv1m3d(nmlonp1,j,:) = adotv1m3d(1,j,:)
          adotv2m3d(nmlonp1,j,:) = adotv2m3d(1,j,:)
          je1oD_pg3d(nmlonp1,j,:)  = je1oD_pg3d(1,j,:)  ! gravity + mag.pressure
          je2oD_pg3d(nmlonp1,j,:)  = je2oD_pg3d(1,j,:)  ! gravity + mag.pressure
        enddo ! jj=1,2
      endif ! calkqlam==1
!
! From sub fieldline_integrals:
!   zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
!   zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
!   zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
!   zigm2  is int[sigma_h] ds,         i.e. Sigma_h
!
!   rim1 is int[(sigma_h-sigma_p*d_1*d_2/D)u_e1 + sigma_p*d_1^2/D u_e2]*B_e3 ds,
!           i.e.  K_(m phi)^D/abs(sin Im)
!   rim2 is -/-int[(sigma_h+sigma_p*d_1*d_2/D)u_e2 - sigma_p*d_2^2/D u_e1]*B_e3 ds, 
!          i.e. +/- K_(m lam)^D 
!         (minus sign in northern hemisphere
!          SH same minus sign but should be plus -> take care of this in transf)
!
! Equatorial values:
! Assume that quantities primarily dependent on Pedersen conductivity
!   have field-line integrals 1/4 as large as the averages for next-higher
!   field lines; quantities primarily dependent on Hall conductivity
!   have field-line integrals 0.12 as large as the averages for next-higher
!   field lines.  Exact values chosen should not be important for potential
!   calculation, as long as they are physically reasonable and not too
!   different from adjacent values.
! add factor 1/2 for the equatorial values since both hemispheres are
! added together - two physical points but equator is only one point
!       
      j = nmlat/2+1
      do i = 1,nmlon
        zigm11(i,j)= .125*(zigm11(i,j-1)+ zigm11(i,j+1))
        zigm22(i,j)= .125*(zigm22(i,j-1)+ zigm22(i,j+1))
        zigmc(i,j) = .125*(zigmc(i,j-1) + zigmc(i,j+1))
#if defined(INTERCOMM) || defined(CISMAH)
        zigm1(i,j) = .06 *(zigm1(i,j-1) + zigm1(i,j+1))     ! cism
#endif
        zigm2(i,j) = .06 *(zigm2(i,j-1) + zigm2(i,j+1))
        rim(i,j,1) = .06 *(rim(i,j-1,1) + rim(i,j+1,1))
        rim(i,j,2) = .06 *(rim(i,j-1,2) + rim(i,j+1,2))
      enddo ! i = 1,nmlon
!
! Include the boundary condition at the equator eq.(5.30) in
! Richmond (1995) Ionospheric Electrodynamics use. Mag. Apex Coord. 
!   
! Sigma_(phi phi)/|sin Im| = 0.5*Sigma_cowling/|sin Im|
!        = 0.5/ |sin Im|*[Simag_(phi phi) -
!                        Sigma_(phi lam)*Sigma_(lam phi)/Sigma_(lam lam)]
!        = 0.5/|sin Im|*[Simag_(phi phi) + 
!                       (Sigma_h-Sigma_c)*(Sigma_h + Sigma_c)/Sigma_(lam lam)]
!  K_(m phi)^D / |sin I_m| =  1/2/|sin Im|*(K_(m phi) - 
!                       Sigma_(phi lam)/Sigma_(lam lam)*K_(m lam)^D)
! factor 1/2 is taken care of above when the conductances at the
! equator are set
!
      j = nmlath      ! nmlath = (nmlat+1)/2
      do i = 1,nmlon
        zigm11(i,j) = zigm11(i,j)+ (zigm2(i,j)-zigmc(i,j))*
     |                (zigm2(i,j)+zigmc(i,j))/zigm22(i,j)
        rim(i,j,1)  = rim(i,j,1) - (zigm2(i,j)-zigmc(i,j))/
     |                zigm22(i,j)*rim(i,j,2)
        zigm11(i,j) = zigm11(i,j) 
        rim(i,j,1)  = rim(i,j,1)  
      enddo

!           
! Using notation of Richmond (1995) on right-hand side below:
! Sigma_(phi phi)/ |sin I_m|  = zigm11
! Sigma_(lam lam) * |sin I_m| = zigm22
! Sigma_(phi lam)             = +/-(zigm2-zigmc)
! Sigma_(lam phi)             = -/+(zigm2+zigmc)
! K_(m phi)^D / |sin I_m|     = rim(1)
! K_(m lam)^D                 = +/-rim(2)
!
! Transforming PDE from original apex (lam_m) to new apex grid (lam_0)
!     lam_m is irregular spaced in latitude
!     lam_0 is regular spaced in latitude (used for derivatives)
! the whole PDE is divided by d lam_m /d lam_0
! 
! sign of K_(m lam)^D in southern hemisphere remains reversed.
! for the mixed terms the transformation factors cancel out (zigmc, zigm2)
! DT1DTS : d lam_0/ d lam_m / |sin I_m|
! RCOS0S : cos(lam_0)/ cos(lam_m)
!
! corfac: |sin I_m|*d lam_m/d lam_0 * cos(lam_0)/ cos(lam_m)
! zigm11: |sin I_m|*d lam_m/d lam_0 * cos(lam_0)/ cos(lam_m)
! zigm22: 1/|sin I_m|*d lam_0/d lam_m * cos(lam_m)/ cos(lam_0)
! rim(1): |sin I_m|*d lam_m/d lam_0
! rim(2): cos(lam_m)/ cos(lam_0)
!
      do j=2,nmlat-1
        corfac = rcos0s(j)/dt1dts(j)
        do i=1,nmlon
          zigm11(i,j) = zigm11(i,j)*corfac
          zigm22(i,j) = zigm22(i,j)/corfac
          rim(i,j,1)  = rim(i,j,1)/dt1dts(j) 
          rim(i,j,2)  = rim(i,j,2)/rcos0s(j) 
        enddo ! i,nmlon
      enddo ! j=2,nmlat-1
!
! Periodic points:
      do j=1,nmlat
        zigm11(nmlonp1,j) = zigm11(1,j)
        zigmc (nmlonp1,j) = zigmc (1,j)
#if defined(INTERCOMM) || defined(CISMAH)
        zigm1 (nmlonp1,j) = zigm1 (1,j)     ! cism
#endif
        zigm2 (nmlonp1,j) = zigm2 (1,j)
        zigm22(nmlonp1,j) = zigm22(1,j)
        rim(nmlonp1,j,:)  = rim(1,j,:)
      enddo ! j=1,nmlat
!
! zigm11 = Sigma_(phi phi)(0)= Sigma_(phi phi)(m) * cos lam_0 /cos lam_m * d lam_m /d lam_0
! zigm22 = Sigma_(lam lam)(0)= Sigma_(lam lam)(m) * cos lam_m /cos lam_0 * d lam_0 /d lam_m
! +-(zigm2-zigmc)= Sigma_(phi lam)(0) = Sigma_(phi lam)(m)
! -+(zigm2+zigmc)= Sigma_(lam phi)(0) = Sigma_(lam phi)(m)
! rim(1) = K_(m phi)^D(0)    = K_(m phi)^D(m) * d lam_m /d lam_0   
! rim(2) = +/-K_(m lam)^D(0) = +/-K_(m lam)^D(m) * cos lam_m /cos lam_0 
!
! Polar values for Z:
      do k=-2,nlevp1
        zpotenm3d(1,1,k) = 
     |    (4.*sddot(nmlon,unitvm,zpotenm3d(1,2,k))-
     |        sddot(nmlon,unitvm,zpotenm3d(1,3,k)))/(3.*float(nmlon))
        zpotenm3d(1,nmlat,k) = 
     |    (4.*sddot(nmlon,unitvm,zpotenm3d(1,nmlat-1,k))-
     |        sddot(nmlon,unitvm,zpotenm3d(1,nmlat-2,k)))/
     |    (3.*float(nmlon))
!
! Extend Z over longitude:
        do i=1,nmlon
          zpotenm3d(i,1,k) = zpotenm3d(1,1,k)
          zpotenm3d(i,nmlat,k) = zpotenm3d(1,nmlat,k)
        enddo ! i=1,nmlon
!
! Periodic points:
        do j=1,nmlat
          zpotenm3d(nmlonp1,j,k) = zpotenm3d(1,j,k)
        enddo ! j=1,nmlat
      enddo ! k=-2,nlevp1
!  
! Compute polar values for the conductances, 4th order interpolation:
! 
      zigm11(1,    1) = (4.*sddot(nmlon,unitvm,zigm11(1,      2))-
     1  sddot(nmlon,unitvm,zigm11(1,      3)))/(3.*float(nmlon))
      zigm11(1,nmlat) = (4.*sddot(nmlon,unitvm,zigm11(1,nmlat-1))-
     1  sddot(nmlon,unitvm,zigm11(1,nmlat-2)))/(3.*float(nmlon))
      zigmc(1,    1) = (4.*sddot(nmlon,unitvm,zigmc(1,      2))-
     1  sddot(nmlon,unitvm,zigmc(1,      3)))/(3.*float(nmlon))
      zigmc(1,nmlat) = (4.*sddot(nmlon,unitvm,zigmc(1,nmlat-1))-
     1  sddot(nmlon,unitvm,zigmc(1,nmlat-2)))/(3.*float(nmlon))
#if defined(INTERCOMM) || defined(CISMAH)
      zigm1(1,    1) = (4.*sddot(nmlon,unitvm,zigm1(1,      2))-     ! cism
     1  sddot(nmlon,unitvm,zigm1(1,      3)))/(3.*float(nmlon))
      zigm1(1,nmlat) = (4.*sddot(nmlon,unitvm,zigm1(1,nmlat-1))-      ! cism
     1  sddot(nmlon,unitvm,zigm1(1,nmlat-2)))/(3.*float(nmlon))
#endif
      zigm2(1,    1) = (4.*sddot(nmlon,unitvm,zigm2(1,      2))-
     1  sddot(nmlon,unitvm,zigm2(1,      3)))/(3.*float(nmlon))
      zigm2(1,nmlat) = (4.*sddot(nmlon,unitvm,zigm2(1,nmlat-1))-
     1  sddot(nmlon,unitvm,zigm2(1,nmlat-2)))/(3.*float(nmlon))
      zigm22(1,    1) = (4.*sddot(nmlon,unitvm,zigm22(1,      2))-
     1  sddot(nmlon,unitvm,zigm22(1,      3)))/(3.*float(nmlon))
      zigm22(1,nmlat) = (4.*sddot(nmlon,unitvm,zigm22(1,nmlat-1))-
     1  sddot(nmlon,unitvm,zigm22(1,nmlat-2)))/(3.*float(nmlon))
! 
! Extend over longitude                                        
      do i = 2,nmlon
        zigm11(i,    1)  = zigm11(1,    1)
        zigm11(i,nmlat)  = zigm11(1,nmlat)
        zigmc(i,     1)  = zigmc(1,     1)
        zigmc(i, nmlat)  = zigmc(1, nmlat)
#if defined(INTERCOMM) || defined(CISMAH)
        zigm1(i,     1)  = zigm1(1,     1)     ! cism
        zigm1(i, nmlat)  = zigm1(1, nmlat)     ! cism
#endif
        zigm2(i,     1)  = zigm2(1,     1)
        zigm2(i, nmlat)  = zigm2(1, nmlat)
        zigm22(i,     1) = zigm22(1,    1)
        zigm22(i, nmlat) = zigm22(1,nmlat)
      enddo ! i = 2,nmlon
!
! RHS vector (I_1,I_2): average over poles:
      do i = 1,nmlon
        rim(i,1,1) = .5*(rim(i,2,1)-rim(1+mod(i-1+nmlon/2,nmlon),2,1))
        rim(i,nmlat,1) = .5*(rim(i,nmlat-1,1)-
     |    rim(1+mod(i-1+nmlon/2,nmlon),nmlat-1,1))
        rim(i,1,2) = .5*(rim(i,2,2)-rim(1+mod(i-1+nmlon/2,nmlon),2,2))
        rim(i,nmlat,2) = .5*(rim(i,nmlat-1,2)-
     |    rim(1+mod(i-1+nmlon/2,nmlon),nmlat-1,2))
      enddo ! i = 1,nmlon
! 
! Periodic points:
      do j=1,nmlat
        zigm11(nmlonp1,j) = zigm11(1,j)
        zigmc (nmlonp1,j) = zigmc (1,j)
#if defined(INTERCOMM) || defined(CISMAH)
        zigm1 (nmlonp1,j) = zigm1 (1,j)        ! cism
#endif
        zigm2 (nmlonp1,j) = zigm2 (1,j)
        zigm22(nmlonp1,j) = zigm22(1,j)
        rim(nmlonp1,j,:)  = rim(1,j,:)
      enddo ! j=1,nmlat
!
! Folding south onto northern hemisphere and calculation of RHS has 
!   been moved to subroutine dynamo.
!
! Save 3d potential on magnetic grid to secondary history:
! PLEASE DO NOT COMMENT THIS OUT -- ZMAG is a mandatory mag field for sech
! real,dimension(nmlonp1,nmlat,-2:nlevp1) :: zpotenm3d
!
      do j=1,nmlat
        call addfld('ZMAG','GEOPOTENTIAL HEIGHT (MAGNETIC)','CM',
     |    zpotenm3d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nimlevp1,j)
      enddo ! j=1,nmlat
!
      end subroutine transf
!-----------------------------------------------------------------------
      subroutine fieldline_integrals(latm)
!
! Perform approximated field line integration at geomagnetic grid
! zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
! zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
! zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
! zigm2  is int[sigma_h] ds,         i.e. Sigma_h
#if defined(INTERCOMM) || defined(CISMAH)
!        zigm1  is int[sigma_p] ds,         i.e. Sigma_p !CISM
#endif
!
! rim1 [A/m] is int[(sigma_h-sigma_p*d_1*d_2/D)u_e1 + sigma_p*d_1^2/D u_e2] *A(h_r)*
!         B_e3 ds, i.e.  K_(m phi)^D/abs(sin Im)
! rim2 [A/m] is -/- int[(sigma_h+sigma_p*d_1*d_2/D)u_e2 - sigma_p*d_2^2/D u_e1] *A(h_r)*
!         B_e3 ds, i.e. +/- K_(m lam)^D 
!        
      use cons_module,only: h0,       ! h0 = 90 km
     |      re,    !  earth radius= 6.37122e8 (cm)
     |      r0,    ! re+h0
     |      ylatm  ! magnetic latitudes
      integer,intent(in) :: latm ! latitude index of field line
!
! Local:
      integer :: i,k
      real :: sinlm,clm2,absinim,ra,sqomrra,sqrra,afac,htfac
      real :: rtrmn,rora,del,omdel,sig1,sig2,ue1,ue2,jgr,jmp,fac,
     |  je1pg,je2pg
      real,dimension(nmlonp1) :: sindm, cosdm, ram, aam, cosiam,
     |  csthdam, rtadram
      real,dimension(nmlonp1,-2:nlevp1) :: rrm, sinidm, cosidm,
     |  costhdm, rtramrm
!
! Quantities needed for field-line integration:
      real,dimension(nmlonp1,-2:nlev) :: htfunc, htfunc2
!
! values at reference height h0
      sinlm   = sin(ylatm(latm))             ! sin(lam_m)
      clm2    = 1. - sinlm*sinlm             ! cos^2(lam_m)
      absinim = abs(sinlm)/sqrt(1.-.75*clm2) ! | sin I_m |
      ra      = r0/clm2                      ! R_A = (R_E + h_0)/cos^2(lam_m)
      sqomrra = sqrt(1.-r0/ra)               ! sqrt(1 - R_0/R_A) = sin(lam_m) 
      sqrra   = sqrt(r0/ra)                  ! sqrt(R_0/R_A) = cos(lam_m) 
      afac    = 2.*sqrt(ra-r0)               ! 2*sqrt(R_A-R_0)
      htfac   = sqrt(ra-.75*r0)              ! sqrt(R_A -3/4*R_0) 
!
! A(h)=A(h0)*htfunc height varying factor for ds = -A(h) d sqrt(h_A-h)
! A(h_0) = 2*sqrt(h_A-h_0)/ |sin I_m| 
      do i=1,nmlon
        aam(i) = afac/abs(siniam(i))
      enddo ! i=1,nmlon
!
! Evaluate 3d fields:
      do k=-2,nlevp1
        do i=1,nmlon
          rrm(i,k) = re+zpotenm(i,k)          ! R = R_E + h
          rtramrm(i,k) = max(0.,ra-rrm(i,k))  ! h_A - h
          rtramrm(i,k) = sqrt(rtramrm(i,k))   ! sqrt( h_A - h )
        enddo ! i=1,nmlon
      enddo ! k=-2,nlevp1
!
! rrm: R interpolation to the middle of horizontal slices
! rtramrm: -d sqrt(h_A-h) 
! htfunc = factor by which to multiply ds = -A(h0)*htfunc* d sqrt(h_A-h) 
!
      do k=-2,nlev
        do i=1,nmlon
          rrm(i,k)     = .5*(rrm(i,k)+rrm(i,k+1))
          rtramrm(i,k) = rtramrm(i,k)-rtramrm(i,k+1)
          htfunc(i,k)  = sqrt(ra-.75*RRM(i,k))/htfac
          htfunc2(i,k) = htfunc(i,k)**2
        enddo ! i=1,nmlon
      enddo ! k=-2,nlev
!
! Initialize coefficients:
      do i=1,nmlon
        zigm11(i,latm) = 0.
        zigm22(i,latm) = 0.
        zigm2 (i,latm) = 0.
        zigmc (i,latm) = 0.
#if defined(INTERCOMM) || defined(CISMAH)
        zigm1 (i,latm) = 0.     ! CISM
#endif
        rim(i,latm,1)  = 0.
        rim(i,latm,2)  = 0.
!	
	tpg1(i,latm) = 0.
	tpg2(i,latm) = 0.
      enddo ! i=1,nmlon
!
! Compute integrals:
! rora :  R/R_A always <1 since h_A > h
!
! assumption for linear interpolation to field line value x(fl)
!    x(fl) = wgt_eq*x(eq)+ wgt_m*x(m)
!  (lam_m- lam)/lam_m ~ (sin lam_m - sin lam) / lam_m for (lam_m- lam) small
!  use dipole approximation cos^2 lam = R/R_A
!  leads to 
!  (lam_m- lam)/lam_m ~ [sqrt(1-R_0/R_A)sqrt(R/R_A) - sqrt(R_0/R_A)sqrt(1-R/R_A)]/lam_m
!
! del  : wgt_eq - weighting factor for equator value
! omdel: wgt_m  - weighting factor for fieldline footpoint value
!
      do k=-2,nlev
        do i=1,nmlon
          rora = min(1.,rrm(i,k)/ra) 
          del = (sqomrra*sqrt(rora)-sqrra*sqrt(1.-rora))/
     |           abs(ylatm(latm))
          omdel = 1. - del
!
! Interpolate conductivities and winds in latitude, assuming
!   linear variation between foot of field line and magnetic equator.
!   (For field lines other than those near the magnetic equator, del is nearly
!   zero, so that the interpolated values are essentially the values for the
!   latitude of the foot of the field line; inaccuracy of the assumption of
!   linear variation is thus unimportant for these field lines.)
!
          sig1 = omdel*sigma_pedm(i,k)  + del*sigma1e(i,k) 
          sig2 = omdel*sigma_hallm(i,k) + del*sigma2e(i,k)
          ue1  = omdel*adotvm(i,k,1) + del*adotve(i,k,1)
          ue2  = omdel*adotvm(i,k,2) + del*adotve(i,k,2)
!
! height varying factors: 
! (note that the factors for d and from ds can cancel)
!    ds        = -A(h0)*htfunc*d sqrt(h_A-h) [cm]
!    d_1^2/D   = 1/htfunc * adotam(i,1)
!    d_2^2/D   = htfunc   * adotam(i,2)
!    d_1*d_2/D = 1        * a1a2m(i)
!
! zigm11: int (sigma_p*d_1^2/D) ds : d_1^2/D & A(h0) outside height loop
! zigm22: int (sigma_p*d_2^2/D) ds : d_2^2/D & A(h0) outside height loop
! zigmc:  int (sigma_p*d_1*d_2/D) ds: d_1*d_2/D & A(h0) outside height loop
! zigm2:  int (sigma_h) ds: A(h0) outside height loop
! rim(1): int Be3[sigma_p*d_1^2/D u_e2+(sigma_h-(sigma_p*d_1*d_2)/D) u_e1] ds: 
!             Be3& A(h0) outside height loop
! rim(2): -/- int Be3[(sigma_h+sigma_p*d_1*d_2/D) u_e2-sigma_p*d_2^2/D u_e1 ] ds: 
!              Be3& A(h0) outside height loop
!
          zigm11(i,latm) = zigm11(i,latm) + sig1*rtramrm(i,k)
          zigm22(i,latm) = zigm22(i,latm) + sig1*rtramrm(i,k)*
     |      htfunc2(i,k)
          zigmc(i,latm) = zigmc(i,latm) + sig1*rtramrm(i,k)*htfunc(i,k)
          zigm2(i,latm) = zigm2(i,latm) + sig2*rtramrm(i,k)*htfunc(i,k)
#if defined(INTERCOMM) || defined(CISMAH)
          zigm1(i,latm) = zigm1(i,latm) + sig1*rtramrm(i,k)*htfunc(i,k)     ! CISM
#endif
          rim(i,latm,1) = rim(i,latm,1) + (sig1*adotam(i,1)*ue2 +(sig2 -
     |      sig1*a1a2m(i))*htfunc(i,k)*ue1)*rtramrm(i,k)
          rim(i,latm,2) = rim(i,latm,2)+(sig1*adotam(i,2)*htfunc2(i,k)*
     |      ue1-(sig2 + sig1*a1a2m(i))*htfunc(i,k)*ue2)*rtramrm(i,k)
!     
! rim1 add = 1/B_e3 * int [-k(ti+te)*grad Ne + rho*g]xB / B^2*d1/D ds
! rim2 add = 1/B_e3 * int [-k(ti+te)*grad Ne + rho*g]xB / B^2*d2/D ds
! je1pg & je2pg [A/cm^2]
!     
	  if(j_pg) then
	    fac = bmodm(i) ! 1/be3
	    fac = rtramrm(i,k)*htfunc(i,k)/fac !   ds/be3/aam
!	    
            je1pg  = omdel*je1_pg(i,k) + del*je1_pge(i,k)
            je2pg  = omdel*je2_pg(i,k) + del*je2_pge(i,k)
!	    
            rim(i,latm,1) = rim(i,latm,1)+ 1.e4*je1pg*fac ! 1e4 for A/cm^2 -> A/m^2
            rim(i,latm,2) = rim(i,latm,2)+ 1.e4*je2pg*fac ! 1e4 for A/cm^2 -> A/m^2
! only for testing	    
!	    tpg1(i,latm) = tpg1(i,latm)+1.e2*je1pg*fac*aam(i)*bmodm(i)  ! A/m 1e4*1e-2
!	    tpg2(i,latm) = tpg2(i,latm)+1.e2*je2pg*fac*aam(i)*bmodm(i)  ! A/m 1e4*1e-2	
	  endif
        enddo ! i=1,nmlon
      enddo ! k=-2,nlev
!
! Complete calculation
!
! At this point, 
!   zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
!   zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
!   zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
!   zigm2  is int[sigma_h] ds,         i.e. Sigma_h
!
!   rim1 is int[(sigma_h-sigma_p*d_1*d_2/D)u_e1 + sigma_p*d_1^2/D u_e2] *
!                B_e3 ds, i.e.  K_(m phi)^D/abs(sin Im)
!   rim2 is -/-int[(sigma_h+sigma_p*d_1*d_2/D)u_e2 - sigma_p*d_2^2/D u_e1] *
!                B_e3 ds, i.e. +/- K_(m lam)^D (minus sign in northern hemisphere
!                SH same minus sign but should be plus -> take care of this in transf)
!
!   ds [cm] therefore convert from [cm] to [m] by multiplying with 1/100
!
      do i = 1,nmlon
        zigm11(i,latm) = 1.e-2*zigm11(i,latm)*aam(i)*adotam(i,1)
        zigm22(i,latm) = 1.e-2*zigm22(i,latm)*aam(i)*adotam(i,2)
        zigmc(i,latm)  = 1.e-2*zigmc(i,latm)*aam(i)*a1a2m(i)
        zigm2(i,latm)  = 1.e-2*zigm2(i,latm)*aam(i)
#if defined(INTERCOMM) || defined(CISMAH)
        zigm1(i,latm)  = 1.e-2*zigm1(i,latm)*aam(i)     ! CISM
#endif
        rim(i,latm,1)  = 1.e-2*rim(i,latm,1)*aam(i)*bmodm(i)  ! [A/m]
        rim(i,latm,2)  = 1.e-2*rim(i,latm,2)*aam(i)*bmodm(i)  ! [A/m]
      enddo ! i = 1,nmlon
!      
      end subroutine fieldline_integrals
!-----------------------------------------------------------------------
      subroutine rhspde
!
!  differentiate RHS 
!  R_0 * d lam_m/d lam_0 * 1/cos lam_0 [ d K_(m phi)^DT(0) / d phi  +
!  (d [ K_(m lam)^DT(0) * cos(lam_0)]/ d lam_0  ] 
!
      use cons_module,only: pi_dyn,dlatm,dlonm,r0
!
! Local:
      integer :: i,j,jj
      real :: pi,fac
      real,dimension(nmlat) :: tint1,tint2,tint3
      real,dimension(-1:nmlonp1+1) :: tint33
!
! External:
      real,external :: sddot ! in util.F
!
      pi = pi_dyn
!
! cos lam_0:
      do j=1,nmlat
        tint1(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo ! j=1,nmlat
!
! Perform differentiation (1/ cos lam0 ) * d K_(m phi)^DT(0)/ d phi
! central differencing
! rhs(i,j) = 1/ cos lam0(j) *  [K_(m phi)^DT(0)(i+1,j)-K_(m phi)^DT(0)(i-1,j)]/ (2 d lon)
! rim(1) =  K_(m phi)^DT(0)
! copy rim(1) into tint33 and set up wrap around points
!
      do j=2,nmlath-1   !  2,48                             
        jj = j+nmlath-1 ! 49,96
        do i=1,nmlon
          tint33(i) = rim(i,jj,1)             ! tint33(1:nmlon)
        enddo ! i=1,nmlon
        do i=1,2                              ! wrap around points
          tint33(i-2)     = tint33(i-2+nmlon) ! -1:0 <= nmlon-1:nmlon
          tint33(i+nmlon) = tint33(i)         ! nmlon+1:nmlon+2 <= 1,2
        enddo ! i=1,2
        do i=1,nmlon                                         
          rhs(i,j) = 1./(dlonm*tint1(nmlath+j-1))*
     |      .5*(tint33(i+1)-tint33(i-1))      ! tint33 2:nmlon+1 and 0:nmlon-1
        enddo ! i=1,nmlon
      enddo ! j=2,nmlath-1
!
! Perform differentiation (1/ cos lam0 ) * d [K_(m lam)^DT(0) cos lam_0]/ d lam_0:
! central differencing
! rhs(i,j) = rhs(i,j) + 1/ cos lam0(j) *  
!  [K_(m lam)^DT(0)(i,j+1)cos lam0(j+1)-K_(m lam)^DT(0)(i,j-1)cos lam0(j-1)]/ (2 d lat)
!
      do j=nmlath+1,nmlat-1 ! 50,96                        
        jj = j-nmlath+1     !  2,48
        do i=1,nmlon
          rhs(i,jj) = rhs(i,jj)+1./(dlatm*tint1(j))*
     |      .5*(rim(i,j+1,2)*tint1(j+1)-rim(i,j-1,2)*tint1(j-1))
        enddo ! i=1,nmlon
      enddo ! j=nmlath+1,nmlat-1
!
! polar value:
      rhs(1,nmlath) = -2./float(nmlon)*
     |  sddot(nmlon,unitvm,rim(1,nmlat-1,2))/tint1(nmlat-1)
!
! equatorial value:
! [K_(m phi)^DT(0)(i+1,j)-K_(m phi)^DT(0)(i-1,j)]/ (2 d lon) +
! [K_(m lam)^DT(0)(i,j+1/2)cos lam0(j+1/2)-0]/ ( d lat/2)
!   with   -  K_(m lam)^DT(0)(i,j_eq) = 0
!          -  cos lam0(j_eq) =1
!
      i = 1
      rhs(i,1) = 0.5/dlonm*(rim(i+1,nmlath,1)-rim(nmlon,nmlath,1))
      rhs(i,1) = rhs(i,1)+1./dlatm*(tint1(nmlath)*rim(i,nmlath,2)+
     |                            tint1(nmlath+1)*rim(i,nmlath+1,2))
      do i = 2,nmlon-1
        rhs(i,1) = 0.5/dlonm*(rim(i+1,nmlath,1)-rim(i-1,nmlath,1))
        rhs(i,1) = rhs(i,1)+1./dlatm*(tint1(nmlath)*rim(i,nmlath,2)+
     |                              tint1(nmlath+1)*rim(i,nmlath+1,2))
      enddo   
      i = nmlon
      rhs(i,1) = 0.5/dlonm*(rim(1,nmlath,1)-rim(i-1,nmlath,1))
      rhs(i,1) = rhs(i,1)+1./dlatm*(tint1(nmlath)*rim(i,nmlath,2)+
     |                            tint1(nmlath+1)*rim(i,nmlath+1,2))
!
! Extend over longitude:
      do i=2,nmlon
        rhs(i,nmlath) = rhs(1,nmlath)
      enddo ! i=2,nmlon
!
! Periodic points:
      do j=1,nmlath
        rhs(nmlonp1,j) = rhs(1,j)
      enddo ! j=1,nmlath
!
! multiply by earth radius + h_0 in meter  = R0*1.E-2 
!  R_0 * d lam_m/d lam_0 * 1/cos lam_0 [ d K_(m phi)^DT(0) / d phi  +
!  (d [ K_(m lam)^DT(0) * cos(lam_0)]/ d lam_0  ] 
!
      fac= r0*1.e-2
      do j=1,nmlath
        do i=1,nmlonp1
          rhs(i,j) = rhs(i,j)*fac
        enddo ! i=1,nmlonp1
      enddo ! j=1,nmlath
!
      end subroutine rhspde
!-----------------------------------------------------------------------
      subroutine clearcee(cee,nlon0,nlat0)
!
! Zero C arrays for stencil coefficients.
! Cee will contain:
!   c0(nmlon0,nmlat0,10), c1(nmlon1,nmlat1,9), c2(nmlon2,nmlat2,9),
!   c3(nmlon3,nmlat3,9),  c4(nmlon4,nmlat4,9)
!
      integer,intent(in) :: nlon0,nlat0
      real,intent(out) :: cee(*)
!
! Local:
      integer :: nlon,nlat,n,m,i
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
      subroutine stencil(zigm,nlon0,nlat0,cee,ncoef)
!
! Calculate contribution for 3 by 3 stencil from coefficient zigm
! at each grid point and level.
!
      integer,intent(in) :: 
     |  nlon0, ! longitude dimension of finest grid level
     |  nlat0, ! latitude dimension of finest grid level
     |  ncoef   ! integer identifier of coefficient
      real,intent(in) :: 
     |  zigm(nlon0,nlat0)  ! coefficients (nlon0+1/2,(nlat0+1)/2) 
      real,intent(inout) :: 
     |  cee(*)  ! output stencil array consisting of c0,c1,c2,c3,c4
!
! Local:
      integer :: nc,nlon,nlat,n
!
! copy Sigma into array and extend in longitude
!
      call htrpex(zigm,nlon0,nlat0)
!
! Calculate contribution to stencil for each grid point and level:
!
      nc = 1
      nlon = nlon0
      nlat = nlat0
      do n=1,5 ! 5 levels of resolution
        call cnm(nlon0,nlat0,nlon,nlat,cee(nc),ncoef)
        nc = nc+9*nlon*nlat
        if (n==1) nc = nc+nlon*nlat
        nlon = (nlon+1)/2
        nlat = (nlat+1)/2
      enddo ! n=1,5
      end subroutine stencil
!-----------------------------------------------------------------------
      subroutine stencmd(zigm,nlon0,nlat0,cee,ncoef)
!
! Calculate contribution for 3 by 3 stencil from coefficients Sigma
! at each grid point and grid level. 
!
! Args:
      integer,intent(in) :: 
     |  nlon0, ! longitude dimension of finest grid level
     |  nlat0, ! latitude dimension of finest grid level
     |  ncoef  ! integer identifier of coefficient
      real,intent(in) :: 
     |  zigm(nlon0,nlat0)   ! coefficients (nlon0,(nlat+1)/2) 
      real,intent(inout) :: 
     |  cee(*)  ! output stencil array consisting of c0,c1,c2,c3,c4
!        
! Local:
      integer :: nc,nlon,nlat,n
!
! Perform half-way interpolation and extend zigm in wkarray:
!
      call htrpex(zigm,nlon0,nlat0)
!
! Calculate contribution to stencil for each grid point and level:
!
      nc = 1
      nlon = nlon0
      nlat = nlat0
!
! Calculate modified and unmodified stencil on finest grid
!
      call cnmmod(nlon0,nlat0,nlon,nlat,cee(nc),ncoef)
!
! Stencils on other grid levels remain the same.
      nc = nc+10*nlon*nlat 
      nlon = (nlon+1)/2
      nlat = (nlat+1)/2
!
      do n=2,5
        call cnm(nlon0,nlat0,nlon,nlat,cee(nc),ncoef)
        nc = nc+9*nlon*nlat
        if (n==1) nc = nc+nlon*nlat
        nlon = (nlon+1)/2
        nlat = (nlat+1)/2
      enddo ! n=1,5
      end subroutine stencmd
!-----------------------------------------------------------------------
      subroutine htrpex(coeff,nmlon0,nmlat0)
!
! Copy coefficients into array and extend array over 16 grid points 
! for the wrap around points on the 5 different grid levels. 
! Result returned in array(1:nmlath) - equator to pole.
!
! Args:
      integer,intent(in) :: nmlon0,nmlat0        ! nmlon, nmlath
      real,intent(in) :: coeff(nmlon0,nmlat0)    ! Sigma(SH) pole to equator
!        
! Local:
      integer :: i,j,jj
!
! Copy coeff into positions in array:
! MODIFY: NO SYM. USE NH VALUES
!
      do j=1,nmlat0       ! 1:nmlath
        do i=1,nmlon0
          wkarray(i,j) = coeff(i,j)
        enddo ! i=1,nmlon0
      enddo ! j=1,nmlat0
!
! Extend over 32 grid spaces to allow for wrap around points for
! a total of 5 grid levels:
!
      do i=1,16
        do j=1,nmlat0
          wkarray(1-i,j)      = wkarray(nmlon0-i,j) 
          wkarray(nmlon0+i,j) = wkarray(1+i,j)
        enddo ! j=1,nmlat0
      enddo ! i=1,16
!
      end subroutine htrpex
!-----------------------------------------------------------------------
      subroutine cnm(nlon0,nlat0,nlon,nlat,c,ncoef)
!
! Compute contribution to stencil from Sigma (array) on grid nlon by nlat,
! Finest grid is nlon0 by nlat0.
! uses upwinding method if diagonal dominance is not preserved
!
! ncoef: integer id of Sigma coefficient:
! ncoef = 1 for zigm11 Sigma_(phi phi)
! ncoef = 2 for zigmc  Sigma_(phi lam)
! ncoef = 3 for zigm2  Sigma_(lam phi)
! ncoef = 4 for zigm22 Sigma_(lam lam)
!
! Args:
      integer,intent(in) :: 
     |  nlon0,nlat0,      ! finest grid dimensions
     |  nlon,nlat         ! output grid dimensions
      integer,intent(in) :: ncoef ! flag for type of Sigma
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
! zigm11 Sigma_(phi phi)
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
! zigmc Sigma_(phi lam)
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
! zigm2  Sigma_(lam phi)
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
! zigm22 Sigma_(lam lam)
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
      subroutine cnmmod(nlon0,nlat0,nlon,nlat,c,ncoef)
!
! Compute contribution to stencil from zigm(ncoef) on finest grid nlon0 by nlat0.
! calculates the stencil with (c) and without (cofum)
! the use of upwinding method if diagonal dominance is not preserved
! output: c     ( upwinding)
!         cofum (no upwinding)
!
! ncoef: integer id of Sigma coefficient:
! ncoef = 1 for zigm11 Sigma_(phi phi)
! ncoef = 2 for zigmc  Sigma_(phi lam)
! ncoef = 3 for zigm2  Sigma_(lam phi)
! ncoef = 4 for zigm22 Sigma_(lam lam)
!
! Args:
      integer,intent(in) :: 
     |  nlon0,nlat0, ! finest grid dimensions
     |  nlon,nlat    ! output grid dimensions
      integer,intent(in) :: ncoef  ! flag for type of Sigma coeff.
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
! zigm11: Sigma_(phi phi)
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
            cofum(i,j,1)=cofum(i,j,1)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            cofum(i,j,5)=cofum(i,j,5)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            cofum(i,j,9) = cofum(i,j,9)-
     |        .5*(wkarray(i0+(i+1)*nint,j0+j*nint)+
     |         2.*wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigmc Sigma_(phi lam)
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
! test for upwinding
            wk(i,2) = (c(i,j,3)+wk(i,1))*(c(i,j,7)-wk(i,1))
            wk(i,3) = sign(wk(i,1),c(i,j,3)+c(i,j,7))
            if (wk(i,2) >= 0.) wk(i,3) = 0.
            c(i,j,3) = c(i,j,3)+wk(i,1)+wk(i,3)
            c(i,j,7) = c(i,j,7)-wk(i,1)+wk(i,3)
            c(i,j,9) = c(i,j,9)-2.*wk(i,3)
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigm2 Sigma_(lam phi)
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
! test for upwinding
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
! zigm22: Sigma_(lam lam)
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
     |                                 wkarray(i0+i*nint,j0+(j+1)*nint))
          cofum(i,j,9) = cofum(i,j,9)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |                                 wkarray(i0+i*nint,j0+(j+1)*nint))

        enddo ! i=1,nlon
      endif ! ncoef
!
      end subroutine cnmmod
!-----------------------------------------------------------------------
      subroutine divide(c,nlon,nlat,nlon0,nlat0,cs,igrid)
!
! Divide stencil C by cos(theta(i,j))
!
      integer,intent(in) :: nlon,nlat,   ! dim. of finest grid
     |    nlon0,nlat0,                   ! dim. of actual grid
     |    igrid                          ! flag for divide by cos lam_0
      real,intent(in) :: cs(*)           ! cos lam_0
      real,intent(out) :: c(nlon,nlat,*) ! coefficient array for grid level
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
      end subroutine divide
!-----------------------------------------------------------------------
      subroutine edges(c,nlon,nlat)
!
! Insert polar boundary conditions in stencil c(nlon,nlat,9)
!    c(:,nlat,1:8) = 0 & diagonal term c(:,nlat,9) = 1.
!
      integer,intent(in) :: nlon,nlat    ! dimension of coeff. array
      real,intent(out) :: c(nlon,nlat,*) ! coefficient array
!        
! Local:
      integer :: n,i
!
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
      subroutine stenmod(inlon,inlat,c,phihm,pfrac)
      use cons_module,only: dlatm,dtr
!
! Modify stencil c and cofum to set potential to heelis value within auroral circle.
!    c     is the stencil w   upwinding
! pfrac set up in heelis module: fractional presence of dynamo equation
!        pfrac = 1  for |colam_m| < crit(2)
!        pfrac = 0  for |colam_m| > crit(1)
!        pfrac = (colat_m-crit(1))/(crit(2)-crit(1))  
!                                for crit(2) < |lam_m| < crit(1)
! to ensure diagonal dominance only the diogonal element c(9) is used
! factor of (dlatm/(10.*dtr))**2 to eliminate dependence on grid spacing
! for c(9)
!
      integer,intent(in) :: inlon,inlat
      real,dimension(nmlon0,nmlat0),intent(in) :: 
     |  phihm,  ! heelis potential (from subs potm, flwv32)
     |  pfrac   ! fractional presence of dynamo (from sub colath)
      real,intent(inout) :: c(inlon,inlat,*)
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
      do j=1,inlat
        jj = j0+j*nint
        do n = 1,8
          do i = 1,inlon
            c(i,j,n) = c(i,j,n)*pfrac(i0+i*nint,jj)
          enddo ! i = 1,inlon
        enddo ! n = 1,8
        do i = 1,inlon
          c(i,j,9) = c(i,j,9)*pfrac(i0+i*nint,jj)+
     |      (1.-pfrac(i0+i*nint,jj))*c(i,j,9)*
     |      (dlatm*float(nint)/(10.*dtr))**2
        enddo ! i = 1,inlon
      enddo ! j=1,inlat
      end subroutine stenmod
!-----------------------------------------------------------------------
      subroutine stenmd(inlon,inlat,c,phihm,pfrac)
      use cons_module,only: dlatm,dtr
!
! Modify stencil c and cofum to set potential to heelis value within auroral circle.
!    c     is the stencil w   upwinding
!    cofum is the stencil w/o upwinding
! pfrac set up in heelis module: fractional presence of dynamo equation
!        pfrac = 1  for |colam_m| < crit(2)
!        pfrac = 0  for |colam_m| > crit(1)
!        pfrac = (colat_m-crit(1))/(crit(2)-crit(1))  
!                                for crit(2) < |lam_m| < crit(1)
! to ensure diagonal dominance only the diogonal element c(9) is used
! factor of (dlatm/(10.*dtr))**2 to eliminate dependence on grid spacing
! for c(9)
!
      integer,intent(in) :: inlon,inlat
      real,dimension(nmlon0,nmlat0),intent(in) :: 
     |  phihm,  ! heelis potential (from subs potm, flwv32)
     |  pfrac   ! fractional presence of dynamo (from sub colath)
      real,intent(inout) :: c(inlon,inlat,*)
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
      subroutine ceee(cee,nx,ny,cf)
!
! Called from mudpack solvers to transfer coefficients.
!
      save
      integer,intent(in) :: nx,ny
      real,intent(in) :: cee(nx,ny,9)
      real,intent(out) :: cf(nx,ny,9)
!
      integer :: i,j,n
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
      subroutine threed
      use cons_module,only: ylatm,r0,r00,pi_dyn,table,dlatm,dlonm,
     |  rcos0s,re,dt1dts,dt0dts,rtd,ylonm
      use fields_module,only: phim3d,
     |  emphi3d,   ! 3d eastward electric field magnetic
     |  emlam3d,   ! 3d equatorw. electric field magnetic
     |  emz3d      ! 3d upward (?) electric field magnetic
      use init_module,only: istep	! for debug
!
! calculates electric field 
!      E_d1 = -1/(R cos lam_m) d Phi/ d phi_m
!      E_d2 = 1/ (R |sinI_m|)  d Phi/ d lam_m
! calculate 3D electric potential and electric field assuming
! that the fieldlines are dipolar and with constant electric potential
! along the field line
! 
! calculates the electric field (used in subroutine efield)
!      R E_phi = -1/(cos lam_m) d Phi/ d phi_m
!      R E_lam = d Phi/ d lam_m
!        E_z   = d Phi/ d zp
! with R the Erath radius
!
! Local:
      integer :: k,i,j,ip1f,ip2f,ip3f
      real :: cosltm,pi,phimn,phims,csth0,sym,ctm,sinltm,siniq
      real,dimension(nmlonp1) :: thetam,pslot,qslot
      integer,dimension(nmlonp1) :: islot,jslot
      real,dimension(nmlonp1,nmlat) :: ed1,ed2
      real,dimension(nmlonp1,nmlat) :: ephi,elam
!
! Externals:
      real,external :: sddot ! in util.F
!
      pi = pi_dyn
!
! Calculate E_d1, E_d2 components of electric field:
! E_d1 = -1/(R cos lam_0)* (cos lam_0/cos lam_m) d Phi/ d phi_m
      do j=2,nmlat-1
        csth0 = cos(-pi/2.+(j-1)*dlatm)
        do i=2,nmlon
          ed1(i,j) = -(phim(i+1,j)-phim(i-1,j))/(2.*dlonm*csth0)*
     |      rcos0s(j)/(r0*1.e-2)
          ephi(i,j) = ed1(i,j)*(r0*1.e-2)
        enddo ! i=2,nmlon
        ed1(1,j) = -(phim(2,j)-phim(nmlon,j))/(2.*dlonm*csth0)*
     |      rcos0s(j)/(r0*1.e-2)
        ed1(nmlonp1,j) = ed1(1,j)
        ephi(1,j)       = ed1(1,j)*(r0*1.e-2)
        ephi(nmlonp1,j) = ed1(nmlonp1,j)*(r0*1.e-2)
      enddo ! j=2,nmlat-1
!
! E_d2 = +/- 1/(R |sinI_m|) d Phi/d lam_0* (d lam_0/ (d lam_m*|sin I_m|))
! Southern hemisphere:
      do j=2,nmlath-1
        do i=1,nmlonp1
          ed2(i,j) = -(phim(i,j+1)-phim(i,j-1))/(2.*dlatm)*dt1dts(j)/
     |      (r0*1.e-2)
          elam(i,j) = -(phim(i,j+1)-phim(i,j-1))/(2.*dlatm)*dt0dts(j)
        enddo ! i=1,nmlonp1
      enddo ! j=2,nmlath-1
!
! Northern hemisphere:
      do j=nmlath+1,nmlat-1
        do i=1,nmlonp1
          ed2(i,j) = (phim(i,j+1)-phim(i,j-1))/(2.*dlatm)*dt1dts(j)/
     |      (r0*1.e-2)
          elam(i,j) = -(phim(i,j+1)-phim(i,j-1))/(2.*dlatm)*dt0dts(j)
        enddo ! i=1,nmlonp1
      enddo ! j=nmlath+1,nmlat-1
!
! Poles: average over four surounding points
      do i = 1,nmlonp1
        ip1f = i + nmlon/4
        if (ip1f > nmlonp1) ip1f = ip1f - nmlon
        ip2f = i + nmlon/2
        if (ip2f > nmlonp1) ip2f = ip2f - nmlon
        ip3f = i + 3*nmlon/4
        if (ip3f > nmlonp1) ip3f = ip3f - nmlon
        ed1(i,1) = .25*(ed1(i,2) - ed1(ip2f,2) +
     |                  ed2(ip1f,2) - ed2(ip3f,2))
        ed1(i,nmlat) = .25*(ed1(i,nmlat-1) - ed1(ip2f,nmlat-1) +
     |                      ed2(ip1f,nmlat-1) - ed2(ip3f,nmlat-1))
        ed2(i,1) = .25*(ed2(i,2) - ed2(ip2f,2) -
     |                  ed1(ip1f,2) + ed1(ip3f,2))
        ed2(i,nmlat) = .25*(ed2(i,nmlat-1) - ed2(ip2f,nmlat-1) -
     |                      ed1(ip1f,nmlat-1) + ed1(ip3f,nmlat-1))
        ephi(i,1) = .25*(ephi(i,2) - ephi(ip2f,2) +
     |                  elam(ip1f,2) - elam(ip3f,2))
        ephi(i,nmlat) = .25*(ephi(i,nmlat-1) - ephi(ip2f,nmlat-1) +
     |                      elam(ip1f,nmlat-1) - elam(ip3f,nmlat-1))
        elam(i,1) = .25*(elam(i,2) - elam(ip2f,2) -
     |                  ephi(ip1f,2) + ephi(ip3f,2))
        elam(i,nmlat) = .25*(elam(i,nmlat-1) - elam(ip2f,nmlat-1) -
     |                      ephi(ip1f,nmlat-1) + ephi(ip3f,nmlat-1))
!
! Equator: derivative of quadratic polynomial (3 given points)
        ed2(i,nmlath) = (4.*phim(i,nmlath+1)-phim(i,nmlath+2)
     |    -3.*phim(i,nmlath))/(2.*dlatm)*dt1dts(nmlath)/(R0*1.e-2)
        ed2(i,nmlath+1) = (4.*phim(i,nmlath+2)-phim(i,nmlath+3)
     |    -3.*phim(i,nmlath+1))/(2.*dlatm)*dt1dts(nmlath+1)/(R0*1.e-2)
        ed2(i,nmlath-1) = (4.*phim(i,nmlath-2)-phim(i,nmlath-3)
     |    -3.*phim(i,nmlath-1))/(2.*dlatm)*dt1dts(nmlath-1)/(R0*1.e-2)
        elam(i,nmlath) = (4.*phim(i,nmlath+1)-phim(i,nmlath+2)
     |    -3.*phim(i,nmlath))/(2.*dlatm)
      enddo ! i = 1,nmlonp1
!
! Comment from tgcm15 (threed.F):
! For each geomagnetic grid point over the 3d grid, determine the
!   location of the foot of the magnetic field passing through the
!   grid point. Assume that the field lines are dipole in shape. 
!   Recalling that the model assumes constant electric potential
!   along each field line, interpolation in the 2d array phim gives
!   the potential at the original grid point.
! ar 3/2010
! For calculating electric fields and ExB drifts as a function of
!  height, neglect height variations of R_E+h, except when those height
!  variations cause large height variations in functions like
!  arccos{sqrt[(R_E+h_0)/(R_E+h)]}.  Consistent with this assumption,
!  neglect height variations of the geomagnetic-field strength.
!
      do k=-2,nlevp1
        do j=2,nmlat-1
          cosltm = cos(ylatm(j))     ! dipole latitude          
	  sinltm = sin(ylatm(j))
          sym = 1.                   ! northern hemisphere
          if (j < nmlath) sym = -1.  ! southern hemisphere
          do i=1,nmlon
!
! latitude lam_m of footpoint
! cos lam_m = cos lam_q [(R_E+h_0)/R_E+h]^0.5 
!
            ctm = sqrt(r0/(re+zpotenm3d(i,j,k)))*cosltm 
            thetam(i) = acos(ctm)
	    siniq = 2.*sinltm/sqrt(4.-3.*cosltm**2)
!
! Find position of lam_m in table of lam_m vs lam_0
! table is spaced 1 deg in lam_m
! Determine lam_0 (theta0) by interpolation in table:
            pslot(i) = thetam(i)*180./pi+1.
            islot(i) = int(pslot(i))
            pslot(i) = pslot(i)-float(islot(i))
            thetam(i) = ((1.-pslot(i))*table(islot(i),2)+pslot(i)*
     |        table(islot(i)+1,2))*sym ! thetam negative for south hem
!
! Now identify geomagnetic grid element in which the point lies:
! longitudinal grid point are the same
            islot(i) = i
            qslot(i) = (thetam(i)+pi/2.)/dlatm+1.
            jslot(i) = int(qslot(i))
            qslot(i) = qslot(i)-float(jslot(i))
!
! Interpolate in phim using linear interpolation in latitude:
! no interpolation necessary in longitude
            phim3d(i,j,k) = (1.-qslot(i))*
     |        phim(islot(i)  ,jslot(i)) +qslot(i)*
     |        phim(islot(i)  ,jslot(i)+1)
!
            ed13d(i,j,k) = (1.-qslot(i))*
     |        ed1(islot(i),jslot(i))+qslot(i)*
     |        ed1(islot(i),jslot(i)+1)
!
            ed23d(i,j,k) = (1.-qslot(i))*
     |        ed2(islot(i),jslot(i))+qslot(i)*
     |        ed2(islot(i),jslot(i)+1)
!
            emphi3d(i,j,k) = (1.-qslot(i))*
     |        ephi(islot(i),jslot(i))+qslot(i)*
     |        ephi(islot(i),jslot(i)+ 1)   
            emlam3d(i,j,k) = -siniq*(r0*1.e-2)*( (1.-qslot(i))*
     |        ed2(islot(i),jslot(i))+qslot(i)*
     |        ed2(islot(i),jslot(i)+1) )
          enddo ! i=1,nmlon
        enddo ! j=2,nmlat-1
      enddo ! k=-2,nlevp1
!
! Calculate value of phim3d at geomagnetic poles:
      phims = sddot(nmlon,unitvm,phim(1,1))/float(nmlon)
      phimn = sddot(nmlon,unitvm,phim(1,nmlat))/float(nmlon)
      do k = -2,nlevp1
        do i = 1,nmlon
          phim3d(i,1,k) = phims
          phim3d(i,nmlat,k) = phimn
!
          ed13d(i,1,k)     = ed1(i,1)
          ed13d(i,nmlat,k) = ed1(i,nmlat)
          ed23d(i,1,k)     = ed2(i,1)
          ed23d(i,nmlat,k) = ed2(i,nmlat)
!
          emphi3d(i,1,k)     = ephi(i,1)
          emphi3d(i,nmlat,k) = ephi(i,nmlat)          
	  emlam3d(i,1,k)     =  ed2(i,1)*(r0*1.e-2)
	  emlam3d(i,nmlat,k) = -ed2(i,nmlat)*(r0*1.e-2) 
        enddo ! i = 1,nmlon
      enddo ! k = -2,nlevp1
!          
      do k=-1,nlev
        do j=1,nmlat
          do i=1,nmlon
             emz3d(i,j,k) = -(phim3d(i,j,k+1)-phim3d(i,j,k-1))
          enddo ! i=1,nmlon
        enddo ! j=2,nmlat-1
      enddo ! k=2,nlevp1-1
!
      do j=1,nmlat
        do i=1,nmlon
          emz3d(i,j,nlevp1) = emz3d(i,j,nlev)
          emz3d(i,j,-2) = emz3d(i,j,-1)
        enddo
      enddo
!
! Periodic point:
      do k=-2,nlevp1
        do j=1,nmlat
          phim3d(nmlonp1,j,k) = phim3d(1,j,k)
          ed13d(nmlonp1,j,k)  = ed13d(1,j,k)
          ed23d(nmlonp1,j,k)  = ed23d(1,j,k)
          emphi3d(nmlonp1,j,k)= emphi3d(1,j,k)
          emlam3d(nmlonp1,j,k)= emlam3d(1,j,k)
          emz3d(nmlonp1,j,k)  = emz3d(1,j,k)
        enddo ! j=1,nmlat
      enddo ! k=1,nlevp1
!
! Save 3d electric potential on magnetic grid to secondary history:
! From fields.F:
! Electric potential on geographic and magnetic grids:
!    |  nmlev = nlev+3,     ! number of vertical levels of mag fields (-2,nlev)
!    |  nimlev  = nmlev,    ! number of interface levels of mag fields (-2,nlev)
!    |  nmlevp1 = nmlev+1,  ! 
!    |  nimlevp1= nimlev+1  ! 
!     real ::
!    |  dynpot(nlonp1,0:nlatp1,nlevp1),   ! 3d electric potential geographic
!    |  phim3d(nmlonp1,nmlat,-2:nlevp1),  ! 3d electric potential magnetic
!    |  emphi3d(nmlonp1,nmlat,-2:nlevp1), ! 3d eastward electric field magnetic
!    |  emlam3d(nmlonp1,nmlat,-2:nlevp1), ! 3d equatorw. electric field magnetic
!    |  emz3d(nmlonp1,nmlat,-2:nlevp1)    ! 3d upward (?) electric field magnetic
! Local:
! real,dimension(nmlonp1,nmlat,-2:nlevp1) :: zpotenm3d
!
!     do j=1,nmlat
!       call addfld('EMPHI3D','R E_mphi','[V/m*m]',
!    |	emphi3d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nmlevp1,j)
!      call addfld('EMLAM3D','R E_mlam','[V/m*m]',
!    |	emlam3d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nmlevp1,j)
!      call addfld('EMZ3D','E_mZ','[V/m*m]',
!    |	emz3d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nmlevp1,j)
!      call addfld('ED1M3D','ed_1 (magn.)','[V/m]',
!    |	ed13d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nmlevp1,j)
!      call addfld('ED2M3D'  ,'ed_2 (magn.)','[V/m]',
!    |	ed23d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nmlevp1,j)
!      call addfld('PHIM3D','El.Poten','V',
!    |	 phim3d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nmlevp1,j)
!      call addfld('ZPOTEN3D','ZPOTEN3D','[cm]',
!    |	 zpotenm3d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nmlevp1,j)
!     enddo ! j=1,nmlat

      do j=1,nmlat
        call mkdiag_ED12('ED1',ed13d(:,j,:),1,nmlonp1,1,nmlevp1,j)
        call mkdiag_ED12('ED2',ed23d(:,j,:),1,nmlonp1,1,nmlevp1,j)
      enddo
!      
      end subroutine threed
!-----------------------------------------------------------------------
      subroutine mag2geo(am,ag,im,jm,dim,djm,lg,lm,nlong,nlatg,nlonm,
     |  nlatm)
!
! Transform field am on geomagnetic grid to ag on the geographic grid using
!   indices im,jm and weights am. Return field ag on geographic grid.
!
! Args:
      integer,intent(in) :: lg,lm,nlong,nlatg,nlonm,nlatm
      integer,intent(in) :: im(lg,*),jm(lg,*)
      real,intent(in) :: am(lm,*),dim(lg,*),djm(lg,*)
      real,intent(out) :: ag(lg,*)
!
! Local:
      integer :: ig,jg
!
      do jg=1,nlatg 
        do ig=1,nlong
          ag(ig,jg) =
     |    am(im(ig,jg)  ,jm(ig,jg))  *(1.-dim(ig,jg))*(1.-djm(ig,jg))+
     |    am(im(ig,jg)+1,jm(ig,jg))  *    dim(ig,jg) *(1.-djm(ig,jg))+
     |    am(im(ig,jg)  ,jm(ig,jg)+1)*(1.-dim(ig,jg))*djm(ig,jg)+
     |    am(im(ig,jg)+1,jm(ig,jg)+1)*    dim(ig,jg) *djm(ig,jg)
        enddo ! ig=1,nlong
      enddo ! jg=1,nlatg 
      end subroutine mag2geo
!-----------------------------------------------------------------------
      subroutine geo2mag(fmag,fgeo,long,latg,wght,nlonp1_geo,nlonp1_mag,
     |  nlon_mag,nlat_mag,lat)
!
! Transform field fgeo on geographic grid to geomagnetic grid using
!   indices long,latg and weights wght. Return field fmag on magnetic grid.
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
      do i=1,nlon_mag
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
!-----------------------------------------------------------------------
#ifdef MPI
      subroutine mp_dynamo_gather
      use mpi_module,only: mytid,lon0,lon1,lat0,lat1,mxlat,mxlon,
     |  ntask,mpi_comm_world,mpi_real8,irstat,tasks,handle_mpi_err
!
! Local:
      integer,parameter :: nfld=9
      integer :: k,i,j,n,nn,ier,len,nlons,nlats,ifld,msgtag,ireqrecv,
     |  ireqsend,ilon0,ilon1,ilat0,ilat1
      integer :: idest = 0
      real :: fmin,fmax
      real,allocatable :: 
     |  rcvbuf(:,:,:,:,:), ! Receive buffer for root task (i,j,k,ifld,ntask)
     |  sndbuf(:,:,:,:)    ! Send buffer for slave tasks (i,j,k,ifld)
      character(len=8) :: fnames(nfld) =
     |  (/'SIGMAPED','SIGMAHAL','ZPOTEN  ','SCHEIGHT','UNVEL   ',
     |    'VNVEL   ','WNVEL   ','JE1oD_F ','JE2oD_F '/)
!
! Allocate send buffer for each task, and a receive buffer for the 
! root task. (mxlon,mxlat are max number of lons,lats held by all tasks)
!
      allocate(rcvbuf(mxlon,mxlat,nlevp1,nfld,ntask),stat=ier)
      if (ier /= 0)
     |  write(6,"('>>> mp_dynamo_gather: error allocating rcvbuf:',
     |    ' mxlon=',i3,' mxlat=',i3,' nlevp1=',i3,' nfld=',i3,
     |    ' ntask=',i3,' ier=',i4)") mxlon,mxlat,nlevp1,nfld,ntask,ier
!
      allocate(sndbuf(mxlon,mxlat,nlevp1,nfld),stat=ier)
      if (ier /= 0)
     |  write(6,"('>>> mp_dynamo_gather: error allocating sndbuf:',
     |    ' mxlon=',i3,' mxlat=',i3,' nlevp1=',i3,' nfld=',i3,
     |    ' ier=',i4)") mxlon,mxlat,nlevp1,nfld,ier
      len = mxlon*mxlat*nlevp1*nfld ! buffer length
!
! Non-root tasks send to master:
      if (mytid /= 0) then  ! mytid /= 0
        sndbuf(:,:,:,:) = 0.
        nlons = lon1-lon0+1
        nlats = lat1-lat0+1
        do ifld=1,nfld
          select case (ifld)
          case(1)
            sndbuf(1:nlons,1:nlats,:,ifld) = 
     |        sigma_ped(lon0:lon1,lat0:lat1,:)
          case(2)
            sndbuf(1:nlons,1:nlats,:,ifld) = 
     |        sigma_hall(lon0:lon1,lat0:lat1,:)
          case(3)
            sndbuf(1:nlons,1:nlats,:,ifld) = 
     |        zpoten(lon0:lon1,lat0:lat1,:)
          case(4)
            sndbuf(1:nlons,1:nlats,:,ifld) = 
     |        scheight(lon0:lon1,lat0:lat1,:)
          case(5)
            sndbuf(1:nlons,1:nlats,:,ifld) = 
     |        unvel(lon0:lon1,lat0:lat1,:)
          case(6)
            sndbuf(1:nlons,1:nlats,:,ifld) = 
     |        vnvel(lon0:lon1,lat0:lat1,:)
          case(7)
            sndbuf(1:nlons,1:nlats,:,ifld) = 
     |        wnvel(lon0:lon1,lat0:lat1,:)
          case(8)
            sndbuf(1:nlons,1:nlats,:,ifld) = 
     |        je1oD_full(lon0:lon1,lat0:lat1,:)
          case(9)
            sndbuf(1:nlons,1:nlats,:,ifld) = 
     |        je2oD_full(lon0:lon1,lat0:lat1,:)
          end select
        enddo
      endif
!
! Gather to root:
      call mpi_gather(sndbuf,len,MPI_REAL8,rcvbuf,len,MPI_REAL8,0,
     |  MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_dynamo_gather: mpi_gather to root')
!
! Root sorts data from each task's subdomain into its global arrays:
! (root task already defined its subdomain in prep_dynamo, so skip
!  that part of the receive buffer (it sent zeroes to task 0 anyway,
!  because task0 did not load the send buffer above).
!
      if (mytid==0) then
        do n=1,ntask-1
          nn = n+1
          ilon0 = tasks(n)%lon0
          ilon1 = tasks(n)%lon1
          ilat0 = tasks(n)%lat0
          ilat1 = tasks(n)%lat1
          nlons = ilon1-ilon0+1
          nlats = ilat1-ilat0+1
          do ifld=1,nfld
            select case (ifld)
            case(1)
              sigma_ped(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld,nn)
            case(2)
              sigma_hall(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld,nn)
            case(3)
              zpoten(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld,nn)
            case(4)
              scheight(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld,nn)
            case(5)
              unvel(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld,nn)
            case(6)
              vnvel(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld,nn)
            case(7)
              wnvel(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld,nn)
            case(8)
              je1oD_full(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld,nn)
            case(9)
              je2oD_full(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld,nn)
            end select
          enddo ! ifld=1,nfld
        enddo ! n=1,ntask-1 (receive from slaves only)
      endif ! (mytid==0) then
!
! Release local buffer space:
      deallocate(rcvbuf)
      deallocate(sndbuf)
!
      end subroutine mp_dynamo_gather
#endif
!-----------------------------------------------------------------------
#if defined(INTERCOMM) || defined(CISMAH)
!FIXME:  cism_ucurrent should go in the "cism_coupling" module...
! BOP
! IROUTINE: cism_ucurrent
! INTERFACE:
!
      subroutine cism_ucurrent
!
! DESCRIPTION:
!
! This subroutine calculate height-integrated neutral wind generated
! field-aligned current (dynamo) to be passed to the M-I coupler to
! solve electric potential. This subroutine is based on the subroutine
! 'nosocoef' in 'current.F' written by Astrid Maute.
! The height-integrated neutral wind field-alined current is calculated in a
! Quisi-Dipole Coordinate that is defined in detail in Richmond (1995). This
! coordinate system removes the 1/|sinI_m| factor in the partial differential
! equation for the electric potential in the Modified Apex Coordinate system.
! 1/|sinI_m| is not defined at magnetic equator in the Modified Apex
! Coordinate system, but is well defined in the Quisi-Dipole Coordinate system
! (I still need to see how this works).
! The neutral dynamo currents are already calculted in subroutine
! 'fieldline-integrals' in the 'dynamo.F' as a global variable
! 'rim(nmlonp1,nmlat,2)'. Subroutine 'rshpde' has the formula to calculate
! height-integrated neutral wind current, but the current there is the sum
! of two hemispheres. We want a global distribution of this current for the M-I
! coupler. Thus the code here is an expanded version of that in "rhspde", but a
! stripped version of "nosocoef". 'nosocoef' also calculates other
! coefficients (lhs) for the potential equation to obtain total field-aligned currents
! including both magnetosphere and thermosphere originated currents. We only need
! thermospheric originated currents for the CISM M-I coupler.
!
! This subroutine is called by subroutine 'dynamo' after 'call transfer'
! in 'dynamo.F'
!                                                  ---------- Wenbin Wang 09/20/05
! USES
!
      use cons_module,only: dlonm,dlatm,pi,r0

!
! PARAMETERS:
! RETURN VALUE: nsrhs(nmlonp1,nmlat)   ! defined as a global variable above
!
! !REVISION HISTORY:
!
! EOP
!
! Calculate height-integrated field-aligned neutral wind currents for both hemisphere
!
! Local:
!
      real :: cs(nmlat)
      real :: dfac
      integer :: j,je,jj,i,n
!
! Externals:
!
      real,external :: sddot ! in util.F
!
! Calculate coefficients for dynamo pde for both hemisphere
!
!
! Clear arrays
!
      nsrhs(:,:)    = 0.0
!
! Calculate magnetic latitude cosin array
!
      do j = 1,nmlat     ! -pi/2 to pi/2
        cs(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo
!
! Calculate right hand side of pde from rim(1) and rim(2)
!
      do j = 2,nmlath-1  ! 2,48  south pole-1 to equator-1
        jj = j+nmlath-1     ! 50,96 equator+1 to north pole-1
!
! Differentiate rim(1) w.r.t lamda
!
        do i = 2,nmlon-1
          nsrhs(i,j)  = 1.0/(dlonm*cs(j))*
     |                      0.5*(rim(i+1,j,1)-rim(i-1,j,1))
          nsrhs(i,jj) = 1.0/(dlonm*cs(jj))*
     |                      0.5*(rim(i+1,jj,1)-rim(i-1,jj,1))
       enddo
!
! Values at longitudinal boundaries
!
        nsrhs(1,j)  = 1.0/(dlonm*cs(j))*
     |                      0.5*(rim(2,j,1)-rim(nmlon,j,1))
        nsrhs(1,jj) = 1.0/(dlonm*cs(jj))*
     |                      0.5*(rim(2,jj,1)-rim(nmlon,jj,1))
        nsrhs(nmlon,j)  = 1.0/(dlonm*cs(j))*
     |                      0.5*(rim(1,j,1)-rim(nmlon-1,j,1))
        nsrhs(nmlon,jj) = 1.0/(dlonm*cs(jj))*
     |                      0.5*(rim(1,jj,1)-rim(nmlon-1,jj,1))
      enddo
!
! Differentiate rim(2) w.r.t theta0
!
      do j = 2,nmlath-1  ! 2,48  south pole -1 to equator-1 "negative sigh"
       jj = j+nmlath-1  ! 50,96 equator+1 to north pole-1
       do i = 1,nmlon
        nsrhs(i,j)  = nsrhs(i,j) - 1.0/(dlatm*cs(j))*0.5*
     |               (rim(i,j+1,2)*cs(j+1)-rim(i,j-1,2)*cs(j-1))
        nsrhs(i,jj) = nsrhs(i,jj) + 1.0/(dlatm*cs(jj))*0.5*
     |               (rim(i,jj+1,2)*cs(jj+1)-rim(i,jj-1,2)*cs(jj-1))
       enddo
      enddo
!
! Calculate value at the poles by averaging over i:nmlon
!
      nsrhs(1,nmlat) = -2./float(nmlon)*
     |                  sddot(nmlon,unitvm,rim(1,nmlat-1,2))/cs(nmlat-1)
      nsrhs(1,1)     = -2./float(nmlon)*
     |                  sddot(nmlon,unitvm,rim(1,2,2))/cs(2)
!
! Extend over longitude
!
      nsrhs(:,nmlat) = nsrhs(1,nmlat)
      nsrhs(:,1)     = nsrhs(1,1)
!
! Calculate equator values
!
      je = nmlath
      i = 1
      nsrhs(i,je) = 0.5/dlonm*(rim(i+1,je,1)-rim(nmlon,je,1))
      nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)*
     |                 rim(i,je,2)+ cs(je+1)*rim(i,je+1,2))
!
      do i = 2,nmlon-1
        nsrhs(i,je) = 0.5/dlonm*(rim(i+1,je,1)-rim(i-1,je,1))
       nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)*
     |                    rim(i,je,2)+ cs(je+1)*rim(i,je+1,2))
      enddo
!
      i = nmlon
      nsrhs(i,je) = 0.5/dlonm*(rim(1,je,1)-rim(i-1,je,1))
      nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)*
     |                 rim(i,je,2)+ cs(je+1)*rim(i,je+1,2))
!
! Periodic points
!
      nsrhs(nmlonp1,:) = nsrhs(1,:)
!
! Scale rhs by refernce radius (R_E + H0) in meters dfac = r0*1e-2
!
      dfac = r0*1.0e-2
      nsrhs(:,:) = -1.*nsrhs(:,:)/dfac
!
      end subroutine cism_ucurrent
#endif
!-----------------------------------------------------------------------
      end module dynamo_module
