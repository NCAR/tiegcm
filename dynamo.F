!
      module dynamo_module
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
! Dymamo calls heelis to add heelis potential at high latitudes.
! See heelis routines in heelis.F.
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
     |  nmlev   ! number of geomagnetic pressure levels (nmlev==nlevp1+3)
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
!
      real,dimension(nmlonp1,nmlat) :: phim
!
! Inputs to the dynamo on the geographic grid. (formerly in dynamo.h):
      real,dimension(nlonp4,nlat,nlevp1) ::
     |  sigma_ped,        ! Pedersen conductivity (see lamdas.F)
     |  sigma_hall,       ! Hall conductivity     (see lamdas.F)
     |  zpoten,           ! geopotential height
     |  scheight,         ! scale height RT/(MBAR*G)
     |  unvel,vnvel,wnvel ! u,v,w neutral velocities (cm/s)
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
!
! Coefficients and RHS terms for PDE on geomagnetic grid:
! (formerly in coefm.h)
!
      real,dimension(nmlonp1,nmlat) ::
     |  zigm11,  ! sigma11*cos(theta0)
     |  zigmc,   ! sigmac
     |  zigm2,   ! sigma2
     |  zigm22   ! sigma22/cos(theta0)
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
!   the hybrid solver is used because SGI crashes in mudmod if isolve==2. 
!
#ifdef IRIX
      integer,parameter :: isolve = 1
#else
      integer,parameter :: isolve = 2 ! default is 2 for new dynamo
!     integer,parameter :: isolve = 1 
!     integer,parameter :: isolve = 0
#endif
!
! Flag for calculation of K_(q,lam)
!   icalkqlam = 1 -> will be calculated
!   icalkqlam = 0 -> will NOT be calculated
!
      integer,parameter :: icalkqlam = 0
!
! 3d arrays for calculation of height-integrated sheet-current density
! K_(q,lam) (if icalkqlam==1)
!
      real,dimension(nmlonp1,nmlat,-2:nlevp1) :: 
     |  sigma1m3d, sigma2m3d, adotv1m3d, adotv2m3d, ed13d, ed23d
      real,dimension(nmlonp1,nmlat) :: 
     |  bmodm3d,a1a2m3d,adotam3d,sinim3d
!
! Equatorial quantities:
      real,dimension(nmlonp1,-2:nlev)   :: sigma1e,sigma2e
      real,dimension(nmlonp1,-2:nlev,2) :: adotve
!
! For current density code (current.F) (icalkqlam==1):
      real :: nscoef(nmlon0,nmlat,10),nscrrt(nmlon0,nmlat)
!
      contains
!-----------------------------------------------------------------------
      subroutine prep_dynamo(tn,un,vn,w,z,barm,ped,hall,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Prepare geographic-grid fields for input to the dynamo, and gather them
! to the root task. This is executed by all tasks, and is called from 
! advance before the dynamo itself (which is executed by master task only).
!
      use cons_module,only: gask,grav
!
! Args:
      integer :: lev0,lev1,lon0,lon1,lat0,lat1
      real :: fmin,fmax
!
! Input fields at task subdomains. These are from current (not updated)
! fields (itp):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,   ! neutral temperature (deg K)
     |  un,   ! neutral zonal velocity (cm/s)
     |  vn,   ! neutral meridional velocity (cm/s)
     |  w,    ! vertical velocity
     |  z,    ! geopotential height
     |  barm, ! mean molecular mass
     |  ped,  ! pedersen conductivity (lamdas.F)
     |  hall  ! hall conductivity     (lamdas.F)
!
! Local:
      integer :: k,i,lat,nk,lonbeg,lonend
!
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
          enddo ! k=lev0,lev1
          zpoten(i,lat,lev1) = z(lev1,i,lat)
        enddo ! i=lon0,lon1
!
! (In tgcm15, these calls are in lamdas.F)
! This is executed by all mpi tasks, so the second mp_gather2root call
!   in advance should be executed: call mp_gather2root(itc,'prim')
!
!       call addfsech_ik('SIGMAPED',' ',' ',sigma_ped(lon0:lon1,lat,:),
!    |    lon0,lon1,nk,nk-1,lat)
!       call addfsech_ik('SIGMAHAL',' ',' ',sigma_hall(lon0:lon1,lat,:),
!    |    lon0,lon1,nk,nk-1,lat)
!       call addfsech_ik('ZPOTEN'  ,' ',' ',zpoten(lon0:lon1,lat,:),
!    |    lon0,lon1,nk,nk-1,lat)
!       call addfsech_ik('SCHEIGHT',' ',' ',scheight(lon0:lon1,lat,:),
!    |    lon0,lon1,nk,nk-1,lat)
!       call addfsech_ik('UNVEL'   ,' ',' ',unvel(lon0:lon1,lat,:),
!    |    lon0,lon1,nk,nk-1,lat)
!       call addfsech_ik('VNVEL'   ,' ',' ',vnvel(lon0:lon1,lat,:),
!    |    lon0,lon1,nk,nk-1,lat)
!       call addfsech_ik('WNVEL'   ,' ',' ',wnvel(lon0:lon1,lat,:),
!    |    lon0,lon1,nk,nk-1,lat)

!       do k=lev0,lev1-1
!         write(6,"('prep_dynamo: k=',i3,' lat=',i3,
!    |      ' zpoten(lonbeg:lonend,lat,k)=',
!    |      /,(6e12.4))") k,lat,zpoten(lonbeg:lonend,lat,k)
!       enddo ! k=lev0,lev1-1

      enddo ! lat=lat0,lat1

!     call fminmax(z(1:nlevp1,lon0:lon1,lat0:lat1),
!    |  nlevp1*(lon1-lon0+1)*(lat1-lat0+1),fmin,fmax)
!     write(6,"('prep_dynamo: z subd min,max=',2e12.4)") fmin,fmax

! zpoten(nlonp4,nlat,nlevp1) ::
!     call fminmax(zpoten(lonbeg:lonend,lat0:lat1,:),
!    |  (lonend-lonbeg+1)*(lat1-lat0+1)*nlevp1,fmin,fmax)
!     write(6,"('prep_dynamo: zpoten subdomain min,max=',2e12.4)")
!    |  fmin,fmax

#ifdef MPI
!
! Gather dynamo input fields to the root task, defining module data 
! above at the global domain on the root task:
!
      call mp_dynamo_gather
#endif
      end subroutine prep_dynamo
!-----------------------------------------------------------------------
      subroutine dynamo
      use cons_module,only: dlatm,dlonm,pi_dyn
      use magfield_module,only: im,jm,dim,djm
      use fields_module,only: dynpot,phim3d
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
! Local:
      integer :: i,j,jj,jjj,j0,jntl,k,n,ncc,nmaglat,nmaglon
      real :: sym
      real :: array(-15:nmlon0+16,nmlat0),cs(nmlat0)
!
! Heelis potential, and related:
      real,dimension(nmlonp1,nmlat0)  :: colatc ! colatitudes
      real,dimension(nmlon0 ,nmlat0)  :: pfrac  ! fraction of potential
      real,dimension(nmlonp1,nmlat)   :: phihm  ! potential in magnetic
      real,dimension(nlonp1,0:nlatp1) :: phih   ! potential in geographic
!
! First transform fields to geomagnetic grid, and perform field line
! integrations:
      call transf
!
! Calculate coefficient of PDE for each hemisphere if K_(q,lam) is 
!   calculated
!
      if (icalkqlam==1) call nosocoef
!
! Fold southern hemisphere over on to northern  (was in transf.F version tgcm15)
! Value at the equator is also folded therefore at the equatorial boundary
! Condition factor 1/2 introduced
! Size of arrays 49,97 from equator to poles (added)
! Southern hemisphere 1,49 original values
!
      do j=1,nmlath
        do i=1,nmlonp1
          zigm11(i,nmlatp1-j) = (zigm11(i,nmlatp1-j)+zigm11(i,j))
          zigmc(i,nmlatp1-j)  =  zigmc(i,nmlatp1-j)+zigmc(i,j)
!
! Reverse sign of zigmc to be compatible with Cicely's
          zigmc(i,nmlat+1-j)  = -zigmc(i,nmlat+1-j)
!
          zigm2(i,nmlat+1-j)  =  zigm2(i,nmlat+1-J)+zigm2(i,j)
          zigm22(i,nmlat+1-j) = (zigm22(i,nmlat+1-J)+zigm22(i,j))
!
! Sign of K_(m lam)^D in southern hemisphere remains reversed; 
!   therefore + sign
          rim(i,nmlatp1-j,1)  = (rim(i,nmlat+1-j,1)+rim(i,j,1))
          rim(i,nmlatp1-j,2)  = (rim(i,nmlat+1-j,2)+rim(i,j,2))  
        enddo ! i=1,nmlonp1
      enddo ! j=1,nmlath
!
! Calculate RHS of PDE from rim(1) and rim(2) (was in transf.F version tgcm15)
!  [( d K_(m phi)^D / d phi /(cos(theta_m)?) +
!  (d [ K_(m lam)^D * cos(lam_m)]/ d lam_m ) /cos ( lam_m) ] * R / (RCOS0S*DT0DTS)
! ~ J_(Mr)*r^2*cos(theta_m)/cos(theta_0)/DT0DTS
!
      call rhspde
!
! Set index array nc and magnetic latitude cosine array:
! nc pointes to the start of the coefficient array for each level
      nc(1) = nc0
      nc(2) = nc1
      nc(3) = nc2
      nc(4) = nc3
      nc(5) = nc4
      nc(6) = ncee
!
! Use pi_dyn from cons module rather than 4*atan(1) to avoid small
! differences generated by the atan in -lmass lib (-lmass was not
! used in earlier versions).
!
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
          zigmc(i,jj)   = (zigmc(i,jj) +zigm2(i,jj))/(4.*dlatm*dlonm)
          zigm2(i,jj)   = zigmc(i,jj)-2.*zigm2(i,jj)/(4.*dlatm*dlonm)
          zigm22(i,jj)  = zigm22(i,jj)*cs(j0+j)/dlatm**2
          zigmc(i,jjj)  = -zigmc(i,jj)
          zigm2(i,jjj)  = -zigm2(i,jj)
          zigm22(i,jjj) = zigm22(i,jj)
        enddo ! i=1,nmlonp1
        if (j /= nmlath) then
!
! Sigma_(phi phi)/( cos(lam_m)*DT0DTS*(Delta lon)^2 )
! zigm11 = southern hemis. 49,1 equator-pole
!
          do i = 1,nmlonp1
            zigm11(i,jj) = zigm11(i,jj)/(cs(j0+j)*dlonm**2)
            zigm11(i,jjj) = zigm11(i,jj)
          enddo
        endif
      enddo ! j=1,nmlath
!
! Set zigm11 to zero at megnetic poles to avoid floating exception 
! (values at poles are not used):
!
      do i = 1,nmlonp1
        zigm11(i,1) = 0.0
        zigm11(i,nmlat) = 0.0
      enddo
!
! Save 2d (mlon,mlat) fields to secondary histories:
      call addfsech_ij('ZIGM11',' ',' ',zigm11,1,nmlonp1,1,nmlat)
      call addfsech_ij('ZIGMC' ,' ',' ',zigmc ,1,nmlonp1,1,nmlat)
      call addfsech_ij('ZIGM2' ,' ',' ',zigm2 ,1,nmlonp1,1,nmlat)
      call addfsech_ij('ZIGM22',' ',' ',zigm22,1,nmlonp1,1,nmlat)
      call addfsech_ij('RIM1'  ,' ',' ',rim(:,:,1),1,nmlonp1,1,nmlat)
      call addfsech_ij('RIM2'  ,' ',' ',rim(:,:,2),1,nmlonp1,1,nmlat)
!
! Clear array for difference stencils at all levels:
      call clearcee(cee,nmlon0,nmlat0)
!
! isolve = 0 -> original mud version 5.
! isolve = 1 -> muh hybrid solver (only as direct solver -- slow)
! isolve = 2 -> modified mudpack solver
!
! isolve = 2 -> modified mudpack solver (modified and unmodified coefficients)
      if (isolve==2) then
        cofum(:,:,:) = 0. ! init
!
! Calculate contribution to stencils from each PDE coefficient
!     subroutine stencmd(zigm,cs,nlon0,nlat0,sym,cee,array,ncoef)
!
! Sigma_(phi phi)/( cos(lam_m)*dt0dts*(Delta lon)^2 )
        array(1,1) = 1.
        sym = 1.
        call stencmd(zigm11,cs,nmlon0,nmlat0,sym,cee,array(1,1),1)
!
! Sigma_(lam lam)*cos(lam_m)*dt0dts/(Delta lam)^2
        array(1,1) = 4.                                 
        sym = 1.
        call stencmd(zigm22,cs,nmlon0,nmlat0,sym,cee,array(1,1),4)
!
! Sigma_(phi lam)/( 4*Delta lam* Delta lon )
        array(1,1) = 2.                                 
        sym = -1.
        call stencmd(zigmc,cs,nmlon0,nmlat0,sym,cee,array(1,1),2)
!
! Sigma_(lam phi)/( 4*Delta lam* Delta lon )
        array(1,1) = 3.                                 
        sym = -1.
        call stencmd(zigm2,cs,nmlon0,nmlat0,sym,cee,array(1,1),3)
!
! isolve /= 2: original or hybrid solver.
      else
!
! Sigma_(phi phi)/( cos(lam_m)*DT0DTS*(Delta lon)^2 )
        array(1,1) = 1.
        sym = 1.
        call stencil(zigm11,nmlon0,nmlat0,sym,cee,array(1,1),1)
!
! Sigma_(lam lam)*cos(lam_m)*DT0DTS/(Delta lam)^2
        array(1,1) = 4.
        sym = 1.
        call stencil(zigm22,nmlon0,nmlat0,sym,cee,array(1,1),4)
!
! Sigma_(phi lam)/( 4*Delta lam* Delta lon )
        array(1,1) = 2.
        sym = -1.
        call stencil(zigmc,nmlon0,nmlat0,sym,cee,array(1,1),2)
!
! Sigma_(lam phi)/( 4*Delta lam* Delta lon )
        array(1,1) = 3.
        sym = -1.
        call stencil(zigm2,nmlon0,nmlat0,sym,cee,array(1,1),3)
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
      call edges(c0,nmlon0,nmlat0)
      call edges(c1,nmlon1,nmlat1)
      call edges(c2,nmlon2,nmlat2)
      call edges(c3,nmlon3,nmlat3)
      call edges(c4,nmlon4,nmlat4)
      if (isolve==2) 
     |  call edges(cofum,nmlon0,nmlat0)
!
! Divide stencils by cos(lam_0) (not rhs):
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
!
! Get heelis colatitudes, pfrac, and potential phihm:
!
      call heelis(colatc,pfrac,phihm,nmlat0,nmlon0,nmlat,nmlon,nmlonp1)
!
! Save single-level heelis potential:
      call addfsech_ij('PHIHM2D','2d HEELIS ELECTRIC POTENTIAL','VOLTS',
     |  phihm,1,nmlonp1,1,nmlat)
!
! Modify stencils and RHS so that heelis potential is inserted at
! high latitude.
!
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
      if(isolve==0) then
        call mud(rim,jntl,isolve) ! solver in mud.F
      elseif (isolve==1) then
        call muh(rim,jntl)        ! solver in muh2cr.F
      elseif (isolve==2) then
        call mudmod(rim,jntl)     ! solver in mudmod.F
      else
        write(6,*) 'dynamo: solver type ',isolve,' not implemented.'
        stop 'isolve'
      endif
!
! Copy output potential from rim to phim(nmlonp1,nmlat):
      do j=1,nmlat
        do i=1,nmlonp1
          phim(i,j) = rim(i,j,1)
        enddo ! i=1,nmlonp1
!       write(6,"('dynamo: j=',i3,' phim(:,j)=',/,(6e12.4))")
!    |    j,phim(:,j)
      enddo ! j=1,nmlat
!
! Save single-level potential to secondary histories (mag grid):
      call addfsech_ij('PHIM2D','2d ELECTRIC POTENTIAL','VOLTS',
     |  phim,1,nmlonp1,1,nmlat)
!
!  am_02/02: calculate current for both hemisphere:
!                      R**2 * J_mr / dt0dts /rcos0s
       if(icalkqlam.eq.1) call nosocrrt
!
! Call threed to generate 3-d potential array in geomagnetic coordinates
!   from 2-d solver output phim. phim3d(nmlonp1,nmlat,-2:nlevp1) is in
!   fields.F.
!
      call threed
!
! Save 3d mag electric potential:
! phim3d(nmlonp1,nmlat,-2:nlevp1)  ! 3d electric potential magnetic
!
      do j=1,nmlat
        call addfsech_ik('PHIM3D','ELECTRIC POTENTIAL (MAG)','VOLTS',
     |    phim3d(:,j,:),1,nmlonp1,nmlev,nmlev-1,j)
      enddo
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
      enddo ! k=1,nlevp1
!
! Periodic point:
      do k=1,nlevp1
        do j=0,nlatp1
          dynpot(nlonp1,j,k) = dynpot(1,j,k)
        enddo ! j=0,nlatp1
      enddo ! k=1,nlevp1
!
! Save electric potential on geographic coords to secondary history:
      do j=1,nlat
        call addfsech_ik('DYNPOT',' ',' ',dynpot(:,j,:),
     |    1,nlonp1,nlevp1,nlevp1-1,j)
      enddo ! j=1,nlat
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
! Save 2d heelis potential on geographic grid:
      call addfsech_ij('PHIH2D','2D HEELIS POTENTIAL (GEOG)','VOLTS',
     |  phih,1,nlonp1,1,nlat)
      end subroutine dynamo
!-----------------------------------------------------------------------
      subroutine transf
      use cons_module,only: h0,r0,dt1dts,rcos0s
      use magfield_module,only: 
     |  alatm,    ! (nlonp1,0:nlatp1)
     |  zb,       ! (nlonp1,0:nlatp1)
     |  bmod,     ! (nlonp1,0:nlatp1)
     |  ig,       ! (nmlonp1,nmlat) geog lon at each geomag grid point
     |  jg,       ! (nmlonp1,nmlat) geog lat at each geomag grid point
     |  wt        ! (4,nmlonp1,nmlat) interpolation weights
      use apex_module,only:           ! (see sub apxparm, apex.F)
     |  be3arr,   ! (nlonp1,nlat), 
     |  dvec,     ! (nlonp1,nlat,3,2)
     |  dddarr    ! (nlonp1,nlat)
!
! Local:
      real,parameter :: unitv(nmlon)=1.
      integer :: i,ii,k,kk,j,jj,lat,n
      real :: z0,sinalat,r0or,rat,corfac
      real,dimension(nlonp1,nlat) :: clm2
      real,dimension(nlonp1,0:nlatp1)        :: a1dta2,sini,be3
!
! Fields to be transformed to geomagnetic space (formerly in transmag.h):
! (these will be input to the geographic to magnetic transformation). 
      real,dimension(nlonp1,0:nlatp1,-2:nlev) ::
     |  ssigma1,  ! pedersen conductivity
     |  ssigma2   ! hall conductivity
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
        enddo ! i=1,nlonp4-3
      enddo ! lat=1,nlat
!
! Small diffs with tgcm15 in sigma_ped and sigma_hall
!     do lat=1,nlat
!       do k=1,nlev
!         do i=1,nlonp1
!           ped_plt(i+2,k) = sigma_ped(i,lat,k)
!           hall_plt(i+2,k) = sigma_hall(i,lat,k)
!           z_plt(i+2,k) = zpoten(i,lat,k)
!           h_plt(i+2,k) = scheight(i,lat,k)
!           u_plt(i+2,k) = unvel(i,lat,k)
!           v_plt(i+2,k) = vnvel(i,lat,k)
!           w_plt(i+2,k) = wnvel(i,lat,k)
!         enddo ! i=1,nlonp1
!       enddo
!       call addfsech_ik('SIGMAPED',' ',' ',ped_plt,
!    |    1,nlonp4,nlevp1,nlev,lat)
!       call addfsech_ik('SIGMAHAL',' ',' ',hall_plt,
!    |    1,nlonp4,nlevp1,nlev,lat)
!       call addfsech_ik('ZPOTEN'  ,' ',' ',z_plt,
!    |    1,nlonp4,nlevp1,nlev,lat)
!       call addfsech_ik('SCHEIGHT',' ',' ',h_plt,
!    |    1,nlonp4,nlevp1,nlev,lat)
!       call addfsech_ik('UNVEL'   ,' ',' ',u_plt,
!    |    1,nlonp4,nlevp1,nlev,lat)
!       call addfsech_ik('VNVEL'   ,' ',' ',v_plt,
!    |    1,nlonp4,nlevp1,nlev,lat)
!       call addfsech_ik('WNVEL'   ,' ',' ',w_plt,
!    |    1,nlonp4,nlevp1,nlev,lat)
!     enddo ! lat=1,nlat
!
      z0 = h0
!
      do j=1,nlat
        do i=1,nlonp1
          sinalat = sin(alatm(i,j))         ! sin(lam)
          clm2(i,j) = 1. - sinalat*sinalat  ! cos^2(lam)
          sini(i,j) = zb(i,j)/bmod(i,j)     ! sin(I_m)
          be3(i,j) = 1.e-9*be3arr(i,j)      ! be3 is in T (be3arr in nT)

        enddo ! i=1,nlonp1
!       write(6,"('transf: j=',i3,' sinalat=',e12.4)") j,sinalat
!       write(6,"('transf: j=',i3,' clm2(:,j)=',/,(6e12.4))") 
!    |    j,clm2(:,j)
!       write(6,"('transf: j=',i3,' sini(:,j)=',/,(6e12.4))") 
!    |    j,sini(:,j)
!       write(6,"('transf: j=',i3,' be3(:,j)=',/,(6e12.4))") 
!    |    j,be3(:,j)
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
        adotv(1,0,k,1) = (9.*adotv(1,1,k,1) - adotv(1,2,k,1))/8.
        adotv(1,nlatp1,k,1) =  (9.*adotv(1,nlat  ,k,1) -
     |                             adotv(1,nlat-1,k,1))/8.
        adotv(1,0,k,2) = (9.*adotv(1,1,k,2) - adotv(1,2,k,2))/8.
        adotv(1,nlatp1,k,2) = (9.*adotv(1,nlat  ,k,2) -
     |                            adotv(1,nlat-1,k,2))/8.
!
! Extend in longitude:
        do i = 2,nlon
          zz(i,0,k)             = zz(1,0,k)
          zz(i,nlatp1,k)        = zz(1,nlatp1,k)
          ssigma1(i,0,k)        = ssigma1(1,0,k)
          ssigma1(i,nlatp1,k)   = ssigma1(1,nlatp1,k)
          ssigma2(i,0,k)        = ssigma2(1,0,k)
          ssigma2(i,nlatp1,k)   = ssigma2(1,nlatp1,k)
          adotv(i,0,k,1)        = (9.*adotv(i,1,k,1) -
     |                                adotv(i,2,k,1))/8.
          adotv(i,nlatp1,k,1)   = (9.*adotv(i,nlat,k,1)-
     |                                adotv(i,nlat-1,k,1))/8.
          adotv(i,0,k,2)        = (9.*adotv(i,1,k,2) -
     |                                adotv(i,2,k,2))/8.
          adotv(i,nlatp1,k,2)   = (9.*adotv(i,nlat,k,2) -
     |                                adotv(i,nlat-1,k,2))/8.
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
        enddo ! k = -2,nlev
        zz(nlonp1,j,nlevp1)= zz(1,j,nlevp1)
        adota(nlonp1,j,1) = adota(1,j,1)
        adota(nlonp1,j,2) = adota(1,j,2)
        a1dta2(nlonp1,j)  = a1dta2(1,j)
        sini(nlonp1,j)    = sini(1,j)
        be3(nlonp1,j)     = be3(1,j)
      enddo ! j = 0,nlatp1
!
! Save to secondary histories:
!     real,dimension(nlonp1,0:nlatp1,-2:nlev) :: ssigma1,ssigma2
!     real,dimension(nlonp1,0:nlatp1,-2:nlev,2) :: adotv
!     real,dimension(nlonp4,nlevp1) :: sigma1_plt,...
!
!     sigma1_plt(:,:) = 0.
!     sigma2_plt(:,:) = 0.
!     zz_plt(:,:) = 0.
!     do j=1,nlat
!       do k=1,nlev
!         do i=1,nlon
!           sigma1_plt(i,k) = ssigma1(i,j,k)
!           sigma2_plt(i,k) = ssigma2(i,j,k)
!           zz_plt(i,k)     = zpoten(i,j,k)
!           adotv1_plt(i,k) = adotv(i,j,k,1)
!           adotv2_plt(i,k) = adotv(i,j,k,2)
!         enddo
!       enddo
!       call addfsech_ik('SIGMA1',' ',' ',sigma1_plt,
!    |    1,nlonp4,nlevp1,nlev,j)
!       call addfsech_ik('SIGMA2',' ',' ',sigma2_plt,
!    |    1,nlonp4,nlevp1,nlev,j)
!       call addfsech_ik('ZZ'    ,' ',' ',zz_plt,
!    |    1,nlonp4,nlevp1,nlev,j)
!       call addfsech_ik('ADOTV1',' ',' ',adotv1_plt,
!    |    1,nlonp4,nlevp1,nlev,j)
!       call addfsech_ik('ADOTV2',' ',' ',adotv2_plt,
!    |    1,nlonp4,nlevp1,nlev,j)
!     enddo ! j=1,nlat

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
      enddo ! k=-2,nlev
!
! Magnetic latitude loop:
      maglat_loop: do j=2,nmlat-1
        do k=-2,nlev
          call geo2mag(sigma_pedm(1,k),ssigma1(1,0,k),ig,jg,wt,nlonp1,
     |      nmlonp1,nmlon,nmlat,j)
          call geo2mag(sigma_hallm(1,k),ssigma2(1,0,k),ig,jg,wt,nlonp1,
     |      nmlonp1,nmlon,nmlat,j)
          do n=1,2
            call geo2mag(adotvm(1,k,n),adotv(1,0,k,n),ig,jg,wt,nlonp1,
     |        nmlonp1,nmlon,nmlat,j)
          enddo ! n=1,2
          sigma_pedm(nmlonp1,k) = sigma_pedm(1,k) ! periodic point
          sigma_hallm(nmlonp1,k) = sigma_hallm(1,k)
        enddo ! k=-2,nlev

!       call addfsech_ik('SIGMA1M'  ,' ',' ',sigma_pedm,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!       call addfsech_ik('SIGMA2M'  ,' ',' ',sigma_hallm,
!    |    1,nmlonp1,nmlev,nmlev-1,j)

        do n=1,2
          call geo2mag(adotam(1,n),adota(1,0,n),ig,jg,wt,nlonp1,
     |      nmlonp1,nmlon,nmlat,j)
        enddo ! n=1,2

        do k=-2,nlevp1
          call geo2mag(zpotenm(1,k),zz(1,0,k),ig,jg,wt,nlonp1,
     |      nmlonp1,nmlon,nmlat,j)
        enddo ! k=-2,nlevp1

        call geo2mag(a1a2m ,a1dta2,ig,jg,wt,nlonp1,nmlonp1,nmlon,nmlat,
     |    j)
        call geo2mag(siniam,sini  ,ig,jg,wt,nlonp1,nmlonp1,nmlon,nmlat,
     |    j)
        call geo2mag(bmodm ,be3   ,ig,jg,wt,nlonp1,nmlonp1,nmlon,nmlat,
     |    j)
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
          enddo ! i=1,nmlon
!
! Without calculation of k_(q,lam)
        else  ! icalkqlam /= 1
          do i=1,nmlon  
            do k=-2,nlev
              if (zpotenm(i,k) < z0) zpotenm(i,k) = z0
              zpotenm3d(i,j,k) = zpotenm(i,k)
            enddo ! k=-2,nlev
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
        enddo ! jj=1,2
      endif ! calkqlam==1
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
      j = nmlat/2+1
      do i = 1,nmlon
        zigm11(i,j)= .125*(zigm11(i,j-1)+ zigm11(i,j+1))
        zigm22(i,j)= .125*(zigm22(i,j-1)+ zigm22(i,j+1))
        zigmc(i,j) = .125*(zigmc(i,j-1) + zigmc(i,j+1))
        zigm2(i,j) = .06 *(zigm2(i,j-1) + zigm2(i,j+1))
        rim(i,j,1) = .06 *(rim(i,j-1,1) + rim(i,j+1,1))
        rim(i,j,2) = .06 *(rim(i,j-1,2) + rim(i,j+1,2))
      enddo ! i = 1,nmlon
!
! Include the boundary condition at the equator eq.(5.30) in
! Richmond (1995) Ionospheric Electrodynamics use. Mag. Apex Coord. 
!   J.Goemag.Geoelectr. 47,191-212
! Sig_phiphi/abs(sin Im) = 0.5*Sig_cowling/abs(sin Im)
!              = 0.5/abs(sin Im)*(Sig_phiphi - Sig_philam*sig_lamphi/Sig_lamlam)
!              = 0.5/abs(sin Im)*(Sig_phiphi + (Sig_h-sig_c)*(Sig_h+sig_c)/Sig_lamlam)
!  rim(1) / |sin I_m| = I_1 = R/2*(K_mphi - Sig_philam/Sig_lamlam*K_mlam)
!
      j = nmlath      ! nmlath = (nmlat+1)/2
      do i = 1,nmlon
        zigm11(i,j) = zigm11(i,j)+ (zigm2(i,j)-zigmc(i,j))*
     |                (zigm2(i,j)+zigmc(i,j))/zigm22(i,j)
        rim(i,j,1)  = rim(i,j,1) - (zigm2(i,j)-zigmc(i,j))/
     |                zigm22(i,j)*rim(i,j,2)
        zigm11(i,j) = 0.5*zigm11(i,j)
        rim(i,j,1)  = 0.5*rim(i,j,1)
      enddo

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
        zigm2 (nmlonp1,j) = zigm2 (1,j)
        zigm22(nmlonp1,j) = zigm22(1,j)
        rim(nmlonp1,j,:)  = rim(1,j,:)
      enddo ! j=1,nmlat
!
! the PDE is divided by 1/ DT0DTS
! Sigma_(phi phi) = zigm11/ rcos0s * dt0dts
! Sigma_(lam lam) = zigm22 * rcos0s / dt0dts
! Sigma_(phi lam) = +-(zigm2-zigmc)
! Sigma_(lam phi) = -+(zigm2+zigmc)
! K_(m phi)^D     =   rim(1) * dt0dts
! K_(m lam)^D     = +-rim(2) * rcos0s
!
! Polar values for Z:
      do k=-2,nlevp1
        zpotenm3d(1,1,k) = 
     |    (4.*sddot(nmlon,unitv,zpotenm3d(1,2,k))-
     |        sddot(nmlon,unitv,zpotenm3d(1,3,k)))/(3.*float(nmlon))
        zpotenm3d(1,nmlat,k) = 
     |    (4.*sddot(nmlon,unitv,zpotenm3d(1,nmlat-1,k))-
     |        sddot(nmlon,unitv,zpotenm3d(1,nmlat-2,k)))/
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
      zigm11(1,    1) = (4.*sddot(nmlon,unitv,zigm11(1,      2))-
     1  sddot(nmlon,unitv,zigm11(1,      3)))/(3.*float(nmlon))
      zigm11(1,nmlat) = (4.*sddot(nmlon,unitv,zigm11(1,nmlat-1))-
     1  sddot(nmlon,unitv,zigm11(1,nmlat-2)))/(3.*float(nmlon))
      zigmc(1,    1) = (4.*sddot(nmlon,unitv,zigmc(1,      2))-
     1  sddot(nmlon,unitv,zigmc(1,      3)))/(3.*float(nmlon))
      zigmc(1,nmlat) = (4.*sddot(nmlon,unitv,zigmc(1,nmlat-1))-
     1  sddot(nmlon,unitv,zigmc(1,nmlat-2)))/(3.*float(nmlon))
      zigm2(1,    1) = (4.*sddot(nmlon,unitv,zigm2(1,      2))-
     1  sddot(nmlon,unitv,zigm2(1,      3)))/(3.*float(nmlon))
      zigm2(1,nmlat) = (4.*sddot(nmlon,unitv,zigm2(1,nmlat-1))-
     1  sddot(nmlon,unitv,zigm2(1,nmlat-2)))/(3.*float(nmlon))
      zigm22(1,    1) = (4.*sddot(nmlon,unitv,zigm22(1,      2))-
     1  sddot(nmlon,unitv,zigm22(1,      3)))/(3.*float(nmlon))
      zigm22(1,nmlat) = (4.*sddot(nmlon,unitv,zigm22(1,nmlat-1))-
     1  sddot(nmlon,unitv,zigm22(1,nmlat-2)))/(3.*float(nmlon))
! 
! Extend over longitude                                        
      do i = 2,nmlon
        zigm11(i,    1)  = zigm11(1,    1)
        zigm11(i,nmlat)  = zigm11(1,nmlat)
        zigmc(i,     1)  = zigmc(1,     1)
        zigmc(i, nmlat)  = zigmc(1, nmlat)
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
        zigm2 (nmlonp1,j) = zigm2 (1,j)
        zigm22(nmlonp1,j) = zigm22(1,j)
        rim(nmlonp1,j,:)  = rim(1,j,:)
      enddo ! j=1,nmlat
!
! Save to secondary history:
!     real,dimension(nmlonp1,nmlat) ::
!    |  zigm11,  ! sigma11*cos(theta0)
!    |  zigmc,   ! sigmac
!    |  zigm2,   ! sigma2
!    |  zigm22   ! sigma22/cos(theta0)
!
!     do j=1,nmlat
!       do i=1,nmlonp1
!         zigm11_plt(i,:) = zigm11(i,j)
!         zigmc_plt (i,:) = zigmc (i,j)
!         zigm2_plt (i,:) = zigm2 (i,j)
!         zigm22_plt(i,:) = zigm22(i,j)
!         rim1_plt  (i,:) = rim   (i,j,1)
!         rim2_plt  (i,:) = rim   (i,j,2)
!       enddo ! i=1,nmlonp1
!       call addfsech_ik('ZIGM11'  ,' ',' ',zigm11_plt,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!       call addfsech_ik('ZIGMC'   ,' ',' ',zigmc_plt,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!       call addfsech_ik('ZIGM2'   ,' ',' ',zigm2_plt,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!       call addfsech_ik('ZIGM22'  ,' ',' ',zigm22_plt,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!       call addfsech_ik('RIM1'    ,' ',' ',rim1_plt,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!       call addfsech_ik('RIM2'    ,' ',' ',rim2_plt,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!     enddo ! i=1,nmlat
! 
! Folding south on to northern hemisphere and calculation of RHS has 
!   been moved to subroutine dynamo.
!
      end subroutine transf
!-----------------------------------------------------------------------
      subroutine fieldline_integrals(latm)
      use cons_module,only: h0,r0,ylatm
!
! Perform approximate field line quadratures at geomagnetic grid.
!
! Args:
      integer,intent(in) :: latm
!
! Local:
      integer :: i,k
      real :: z0
      real :: sinlm,clm2,absinim,ra,sqomrra,sqrra,afac,htfac
      real :: rtrmn,rora,del,omdel,sig1,sig2,ue1,ue2
      real,dimension(nmlonp1) :: sindm, cosdm, ram, aam, cosiam,
     |  csthdam, rtadram
      real,dimension(nmlonp1,-2:nlevp1) :: rrm, sinidm, cosidm,
     |  costhdm, rtramrm
!
! Quantities needed for field-line integration:
      real,dimension(nmlonp1,-2:nlev) :: htfunc, htfunc2
!
      z0 = h0
      sinlm   = sin(ylatm(latm))             ! sin(lam_m)
      clm2    = 1. - sinlm*sinlm             ! cos^2(lam_m)
      absinim = abs(sinlm)/sqrt(1.-.75*clm2) ! | sin I_m |
      ra      = r0/clm2                      ! (R_E + H_0)/cos^2(lam_m)
      sqomrra = sqrt(1.-r0/ra)               ! sqrt(1/ R_0/R_A) = sin(lam_m)
      sqrra   = sqrt(r0/ra)                  ! sqrt(R_0/R_A)
      afac    = 2.*sqrt(ra-r0)               ! 2*sqrt(R_A-R_0)
      htfac   = sqrt(ra-.75*r0)              ! sqrt(R_A -3/4*R_0)    
                                             ! ds = (A dh) / (2*sqrt(h_A-h))
!
! 2*sqrt( h_A - h_0 )/ |sin I_m | w.r to reference height A(h_R)
      do i=1,nmlon
        aam(i) = afac/abs(siniam(i))
      enddo ! i=1,nmlon
!
! Evaluate 3d fields:
      do k=-2,nlevp1
        do i=1,nmlon
!
! rr = r0+z-z0 radius of magnetic point
          rrm(i,k) = r0+zpotenm(i,k)-z0
!
! rtramr = ra-r if +ive, zero otherwise
          rtramrm(i,k) = amax1(0.,ra-rrm(i,k))
          rtramrm(i,k) = sqrt(rtramrm(i,k))
        enddo ! i=1,nmlon
      enddo ! k=-2,nlevp1
!
! Interpolation to the middle of horizontal slices:
! htfunc = factor by which to multiply AAM(I)*d(sqrt(ra-r)) = ds

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
        rim(i,latm,1)  = 0.
        rim(i,latm,2)  = 0.
      enddo ! i=1,nmlon
!
! Compute integrals:
      do k=-2,nlev
        do i=1,nmlon
!
! (R_E+h)/(R_E+h_A) < 1 -> h_A > h
          rora = amin1(1.,rrm(i,k)/ra) 
!
! (lam_m - lam) / lam_m = 
! sqrt(1-r_0/r_A)sqrt(r/r_A) - sqrt(r_0/r_A)sqrt(1-r/r_A)
          del = (sqomrra*sqrt(rora)-sqrra*sqrt(1.-rora))/
     |      abs(ylatm(latm))
          omdel = 1. - del
!
! Interpolate conductivities and winds in latitude along field line, assuming
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
! height varying factors: ds = aam*htfunc
!    d_1^2/D   = 1/htfunc * adotam(i,1)
!    d_2^2/D   = htfunc   * adotam(i,2)
!    d_1*d_2/D = 1        * a1a2m(i)
!
! int (sigma_p*d_1^2/D) ds : d_1^2/D s. below loop 11
          zigm11(i,latm) = zigm11(i,latm) + sig1*rtramrm(i,k)
!
! int (sigma_p*d_2^2/D) ds : d_2^2/D s. below loop 11
          zigm22(i,latm) = zigm22(i,latm) + sig1*rtramrm(i,k)*
     |      htfunc2(i,k)
!
! Following uses Richmond convention for sign of zigmc 
! (= -zigmc of Ridley)
!
! zigmc: int (sigma_p*d_1*d_2/D) ds
! zigm2: int (sigma_h) ds
! rim(1): int [sigma_p*d_1^2/D u_e2+(sigma_h-(sigma_p*d_1*d_2)/D) u_e1] ds
! rim(2): int [(sigma_h+sigma_p*d_1*d_2/D) u_e2-sigma_p*d_2^2/D u_e1 ] ds
!
          zigmc(i,latm) = zigmc(i,latm) + sig1*rtramrm(i,k)*htfunc(i,k)
          zigm2(i,latm) = zigm2(i,latm) + sig2*rtramrm(i,k)*htfunc(i,k)
          rim(i,latm,1) = rim(i,latm,1) + (sig1*adotam(i,1)*ue2 +(sig2 -
     |      sig1*a1a2m(i))*htfunc(i,k)*ue1)*rtramrm(i,k)
          rim(i,latm,2) = rim(i,latm,2)+(sig1*adotam(i,2)*htfunc2(i,k)*
     |      ue1-(sig2 + sig1*a1a2m(i))*htfunc(i,k)*ue2)*rtramrm(i,k)
        enddo ! i=1,nmlon
      enddo ! k=-2,nlev
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
      do i = 1,nmlon
        zigm11(i,latm) = 1.e-2*zigm11(i,latm)*aam(i)*adotam(i,1)
        zigm22(i,latm) = 1.e-2*zigm22(i,latm)*aam(i)*adotam(i,2)
        zigmc(i,latm)  = 1.e-2*zigmc(i,latm)*aam(i)*a1a2m(i)
        zigm2(i,latm)  = 1.e-2*zigm2(i,latm)*aam(i)
        rim(i,latm,1)  = 1.e-2*rim(i,latm,1)*aam(i)*bmodm(i)
        rim(i,latm,2)  = 1.e-2*rim(i,latm,2)*aam(i)*bmodm(i)
      enddo ! i = 1,nmlon
      end subroutine fieldline_integrals
!-----------------------------------------------------------------------
      subroutine rhspde
      use cons_module,only: pi_dyn,dlatm,dlonm,r0
!
! am_02/02  calculate the RHS of the PDE from rim(1) and rim(2)
! this part was in transf.F before- seems to be more logical here
!
! Local:
      integer :: i,j,jj
      real :: pi
      real,dimension(nmlat) :: tint1,tint2,tint3
      real,dimension(-1:nmlat) :: tint33
      real,parameter :: unitv(nmlon)=1.
!
! External:
      real,external :: sddot ! in util.F
!
      pi = pi_dyn
!
! Calculate right hand side of PDE from rim:
      do j=1,nmlat
        tint1(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo ! j=1,nmlat
!
! nmlat=97, nmlath=(nmlat+1)/2=49
      do j=2,nmlath-1   !  2,48                              ! do 45
        jj = j+nmlath-1 ! 49,96
        do i=1,nmlon
          tint33(i) = rim(i,jj,1)         ! tint33(1:nmlon)
        enddo ! i=1,nmlon
        do i=1,2
          tint33(i-2) = tint33(i-2+nmlon) ! -1:0 <= nmlon-1:nmlon
          tint33(i+nmlon) = tint33(i)     ! nmlon+1:nmlon+2 <= 1,2
        enddo ! i=1,2
!
! Perform differentiation of rim(1) w.r.t. lamda:
! ( d K_(m phi)^D / d phi )/ (cos (lam_0) / (d lam_0/ d lam_m) =
! ( d K_(m phi)^D / d phi )/ cos(lam_m) / DT0DTS / RCOS0S =
! ( d K_(m phi)^D(0) / d phi )/ cos(theta_0)
!
        do i=1,nmlon                                         ! do 40
          rhs(i,j) = 1./(dlonm*tint1(nmlath+j-1))*
     |      .5*(tint33(i+1)-tint33(i-1))  ! tint33 2:nmlon+1 and 0:nmlon-1
        enddo ! i=1,nmlon
      enddo ! j=2,nmlath-1
!
! Perform differentiation of rim(2) w.r.t. lam_0:
!  +/- (d [ K_(m lam)^D * cos(lam_m)]/ d lam_0 ) /cos ( lam_0) =
!  + (d [ K_(m lam)^D * cos(lam_m)]/ d lam_m ) /cos ( lam_m) / (RCOS0S*DT0DTS) =
!  +/- (d [ K_(m lam)^D(0) * cos(lam_0)]/ d lam_0 ) /cos ( lam_0) =
!
      do j=nmlath+1,nmlat-1 ! 50,96                          ! do 41
        jj = j-nmlath+1     !  2,48
        do i=1,nmlon
          rhs(i,jj) = rhs(i,jj)+1./(dlatm*tint1(j))*
     |      .5*(rim(i,j+1,2)*tint1(j+1)-rim(i,j-1,2)*tint1(j-1))
        enddo ! i=1,nmlon
      enddo ! j=nmlath+1,nmlat-1
!
! Now deal with nmlat:
      rhs(1,nmlath) = -2./float(nmlon)*
     |  sddot(nmlon,unitv,rim(1,nmlat-1,2))/tint1(nmlat-1)
!
! Equator:
! Include the boundary condition at the equator.
! rhs(equator)/R = d (K_mphi^DT(0) - sig_philam/sig_lamlam*K_mlam^DT(0)) / d phi_m
!                + d (cos lam_0 * K_mlam^DT(0))/ d lam_0
! from Cicely's notes:
! I_1 = 0.5*(K_(m phi)^DT(0) - Sig_(phi lam)/Sig_(lam lam)*K_(ml am)^DT(0))
! I_2 = K_(m lam)^DT(0)
! differentiate
! rhs = (I_1(i+1/2,j)-I_1(i-1/2,j))/dlonm +
!       (2*cos(lam_0)_(j+1/2)*I_2(i,j+1/2))/dlat_0
!
      i = 1
      rhs(i,1) = 0.5/dlonm*(rim(i+1,nmlath,1)-rim(nmlon,nmlath,1))
      rhs(i,1) = rhs(i,1)+1./dlatm*(tint1(nmlath)*rim(i,nmlath,2)+
     |                            tint1(nmlath+1)*rim(i,nmlath+1,2))
      do i = 2,nmlon-1
        rhs(i,1) = 0.5/dlonm*(rim(i+1,nmlath,1)-rim(i-1,nmlath,1))
        rhs(i,1) = rhs(i,1)+1./dlatm*(tint1(nmlath)*rim(i,nmlath,2)+
     |                              tint1(nmlath+1)*rim(i,nmlath+1,2))
      enddo ! i = 2,nmlon-1
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
! scale (multiply by earth radius in meter  = R0*1.E-2)
![( d K_(m phi)^D / d phi /(cos(theta_m)?) +
! (d [ K_(m lam)^D * cos(lam_m)]/ d lam_m ) /cos ( lam_m) ] * R / (RCOS0S*DT0DTS)
! ~ J_(Mr)*r^2*cos(theta_m)/cos(theta_0)/DT0DTS
! theta_m = theta_s
!
      do j=1,nmlath
        do i=1,nmlonp1
          rhs(i,j) = rhs(i,j)*r0*1.e-2
        enddo ! i=1,nmlonp1
      enddo ! j=1,nmlath
!
! Save rhs to secondary histories (redundant in vertical and
!  zero in north hem):
!     real,dimension(nmlonp1,nmlath)  :: rhs ! right-hand side 
!     real,dimension(nmlonp1,-2:nlev) :: rhs_plt ! diag
!
!     rhs_plt = 0.
!     do j=1,nmlat
!       if (j <= nmlath) then
!         do i=1,nmlonp1
!           rhs_plt(i,:) = rhs(i,j)
!         enddo
!       endif
!       call addfsech_ik('RHS',' ',' ',rhs_plt,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!     enddo ! j=1,nmlat
      end subroutine rhspde
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
      subroutine stencil(zigm,nlon0,nlat0,sym,cee,array,ncoef)
!
! Calculate contribution fo 3 by 3 stencil from coefficient zigm
! at each grid point and level.
!
! Args:
      integer,intent(in) :: 
     |  nlon0, ! longitude dimension of finest grid level
     |  nlat0, ! latitude dimension of finest grid level
     |  ncoef   ! integer identifier of coefficient
      real,intent(in) :: 
     |  zigm(nlon0,nlat0), ! coefficients (nlon0+1/2,(nlat0+1)/2) 
     |  sym                  !  1. if zigm symmetric w.r.t. equator 
                             ! -1. for antisymmetry
      real,intent(inout) ::
     |  array(nlon0+32,1)   ! work space (-15:nlon0+16,-15:nlat0+2)
      real,intent(inout) :: 
     |  cee(*)  ! output stencil array consisting of c0,c1,c2,c3,c4
!
! Local:
      integer :: nc,nlon,nlat,n
!
! Perform half-way interpolation and extend zigm in array:
!
      call htrpex(zigm,array,nlon0,nlat0,sym)
!
! Calculate contribution to stencil for each grid point and level:
      nc = 1
      nlon = nlon0
      nlat = nlat0
      do n=1,5 ! 5 levels of resolution
        call cnm(array,nlon0,nlat0,nlon,nlat,cee(nc),ncoef)
        nc = nc+9*nlon*nlat
        if (n==1) nc = nc+nlon*nlat
        nlon = (nlon+1)/2
        nlat = (nlat+1)/2
      enddo ! n=1,5
      end subroutine stencil
!-----------------------------------------------------------------------
      subroutine stencmd(zigm,cs,nlon0,nlat0,sym,cee,array,ncoef)
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
     |  array(nlon0+32,1)   ! work space (-15:nlon0+16,-15:nlat0+2)
      real,intent(inout) :: 
     |  cee(*)  ! output stencil array consisting of c0,c1,c2,c3,c4
!
! Local:
      integer :: nc,nlon,nlat,n
!
! Perform half-way interpolation and extend zigm in array:
!
      call htrpex(zigm,array,nlon0,nlat0,sym)
!
! Calculate contribution to stencil for each grid point and level:
!
      nc = 1
      nlon = nlon0
      nlat = nlat0
!
! Calculate modified and unmodified stencil on finest grid
!
      call cnmmod(array,nlon0,nlat0,nlon,nlat,cee(nc),ncoef)
!
! Stencils on other grid levels remain the same.
      nc = nc+10*nlon*nlat 
      nlon = (nlon+1)/2
      nlat = (nlat+1)/2
!
      do n=2,5
        call cnm(array,nlon0,nlat0,nlon,nlat,cee(nc),ncoef)
        nc = nc+9*nlon*nlat
        if (n==1) nc = nc+nlon*nlat
        nlon = (nlon+1)/2
        nlat = (nlat+1)/2
      enddo ! n=1,5
      end subroutine stencmd
!-----------------------------------------------------------------------
      subroutine htrpex(coeff,array,nmlon0,nmlat0,sym)
!
! Perform half-way interpolation on array coeff and extend over 16 grid
! points. Result returned in array.
!
! Args:
      integer,intent(in) :: nmlon0,nmlat0
      real,intent(in) :: coeff(nmlon0,nmlat0),sym
      real,intent(out) :: array(nmlon0+32,*)
!
! Local:
      integer :: i,j,jj
!
! Copy coeff into positions in array:
      do j=1,nmlat0
        jj = nmlat0-j+1
        do i=1,nmlon0
          array(i,j) = sym*coeff(i,jj)
        enddo ! i=1,nmlon0
      enddo ! j=1,nmlat0
!
! Extend over 32 grid spaces to allow for a total of 5 grid levels:
      do i=1,16
        do j=1,nmlat0
          array(1-i,j) = array(nmlon0-i,j) 
          array(nmlon0+i,j) = array(1+i,j)
        enddo ! j=1,nmlat0
      enddo ! i=1,16
      end subroutine htrpex
!-----------------------------------------------------------------------
      subroutine cnm(array,nlon0,nlat0,nlon,nlat,c,ncoef)
!
! Compute contribution to stencil from zigm(ncoef) on grid nlon by nlat,
! Finest grid is nlon0 by nlat0.
!
! Args:
      integer,intent(in) :: 
     |  nlon0,nlat0, ! finest grid dimensions
     |  nlon,nlat    ! output grid dimensions
!
! ncoef: integer id of coefficient:
! ncoef = 1 for zigm11
! ncoef = 2 for zigm12 (=zigmc+zigm2)
! ncoef = 3 for zigm21 (=zigmc-zigm2)
! ncoef = 4 for zigm22
!
      integer,intent(in) :: ncoef

      real,intent(in) :: 
     |  array(nlon0+32,*) ! zigm on grid nlon0 x nlat0
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
! Scan array nlon x nlat calculating and adding contributions to stencil
! from zigm(ncoef)
      i0 = 1-nint
      j0 = 1-nint
!
! zigm11:
! am 2001-6-27 include boundary condition at equator
      if (ncoef==1) then
        do j = 1,nlat-1
          do i = 1,nlon
            c(i,j,1) = c(i,j,1)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i+1)*nint,j0+j*nint))
            c(i,j,5) = c(i,j,5)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
            c(i,j,9) = c(i,j,9)-.5*(array(i0+(i+1)*nint,j0+j*nint)+
     |        2.*array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigm12 (=zigmc+zigm2)
      elseif (ncoef==2) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,2) = c(i,j,2)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i+1)*nint,j0+j*nint))
            c(i,j,4) = c(i,j,4)-.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
            c(i,j,6) = c(i,j,6)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
            c(i,j,8) = c(i,j,8)-.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i+1)*nint,j0+j*nint))
            wk(i,1) = .5*(array(i0+(i+1)*nint,j0+j*nint)-
     |        array(i0+(i-1)*nint,j0+j*nint))
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
            c(i,j,2) = c(i,j,2)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+i*nint,j0+(j+1)*nint))
            c(i,j,4) = c(i,j,4)-.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+i*nint,j0+(j+1)*nint))
            c(i,j,6) = c(i,j,6)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+i*nint,j0+(j-1)*nint))
            c(i,j,8) = c(i,j,8)-.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+i*nint,j0+(j-1)*nint))
            wk(i,1) = .5*(array(i0+i*nint,j0+(j+1)*nint)-
     |        array(i0+i*nint,j0+(j-1)*nint))
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
          c(i,j,2) = c(i,j,2)+.5*(array(i0+i*nint,j0+j*nint)+
     |      array(i0+i*nint,j0+(j+1)*nint))
          c(i,j,4) = c(i,j,4)-.5*(array(i0+i*nint,j0+j*nint)+
     |      array(i0+i*nint,j0+(j+1)*nint))
          wk(i,1) = .5*(array(i0+i*nint,j0+j*nint)+
     |      array(i0+i*nint,j0+(j+1)*nint))
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
            c(i,j,3) = c(i,j,3)+.5*(array(i0+i*nint,j0+j*nint)
     |        +array(i0+i*nint,j0+(j+1)*nint))
            c(i,j,7) = c(i,j,7)+.5*(array(i0+i*nint,j0+j*nint)
     |        +array(i0+i*nint,j0+(j-1)*nint))
            c(i,j,9) = c(i,j,9)-.5*(array(i0+i*nint,j0+(j-1)*nint)
     |        +2.*array(i0+i*nint,j0+j*nint)
     |        +array(i0+i*nint,j0+(j+1)*nint))
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
        j = 1
        do i=1,nlon
          c(i,j,3) = c(i,j,3)+.5*(array(i0+i*nint,j0+j*nint)
     |      +array(i0+i*nint,j0+(j+1)*nint))
          c(i,j,9) = c(i,j,9)-.5*(array(i0+i*nint,j0+j*nint)+
     |                            array(i0+i*nint,j0+(j+1)*nint))
        enddo ! i=1,nlon
      endif ! ncoef
      end subroutine cnm
!-----------------------------------------------------------------------
      subroutine cnmmod(array,nlon0,nlat0,nlon,nlat,c,ncoef)
!
! Compute contribution to stencil from zigm(ncoef) on grid nlon by nlat,
! Finest grid is nlon0 by nlat0.
!
! Args:
      integer,intent(in) :: 
     |  nlon0,nlat0, ! finest grid dimensions
     |  nlon,nlat    ! output grid dimensions
!
! ncoef: integer id of coefficient:
! ncoef = 1 for zigm11
! ncoef = 2 for zigm12 (=zigmc+zigm2)
! ncoef = 3 for zigm21 (=zigmc-zigm2)
! ncoef = 4 for zigm22
!
      integer,intent(in) :: ncoef

      real,intent(in) :: 
     |  array(nlon0+32,*) ! zigm on grid nlon0 x nlat0
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
! Scan array nlon x nlat calculating and adding contributions to stencil
! from zigm(ncoef)
      i0 = 1-nint
      j0 = 1-nint
!
! zigm11:
! am 2001-6-27 include boundary condition at equator
      if (ncoef==1) then
        do j = 1,nlat-1
          do i = 1,nlon
            c(i,j,1) = c(i,j,1)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i+1)*nint,j0+j*nint))
            c(i,j,5) = c(i,j,5)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
            c(i,j,9) = c(i,j,9)-.5*(array(i0+(i+1)*nint,j0+j*nint)+
     |        2.*array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
!
! Unmodified:
            cofum(i,j,1) = cofum(i,j,1)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i+1)*nint,j0+j*nint))
            cofum(i,j,5) = cofum(i,j,5)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
            cofum(i,j,9) = cofum(i,j,9)-
     |        .5*(array(i0+(i+1)*nint,j0+j*nint)+
     |         2.*array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigm12 (=zigmc+zigm2)
      elseif (ncoef==2) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,2) = c(i,j,2)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i+1)*nint,j0+j*nint))
            c(i,j,4) = c(i,j,4)-.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
            c(i,j,6) = c(i,j,6)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i-1)*nint,j0+j*nint))
            c(i,j,8) = c(i,j,8)-.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+(i+1)*nint,j0+j*nint))
            wk(i,1) = .5*(array(i0+(i+1)*nint,j0+j*nint)-
     |        array(i0+(i-1)*nint,j0+j*nint))
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
            c(i,j,2) = c(i,j,2)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+i*nint,j0+(j+1)*nint))
            c(i,j,4) = c(i,j,4)-.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+i*nint,j0+(j+1)*nint))
            c(i,j,6) = c(i,j,6)+.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+i*nint,j0+(j-1)*nint))
            c(i,j,8) = c(i,j,8)-.5*(array(i0+i*nint,j0+j*nint)+
     |        array(i0+i*nint,j0+(j-1)*nint))
            wk(i,1) = .5*(array(i0+i*nint,j0+(j+1)*nint)-
     |        array(i0+i*nint,j0+(j-1)*nint))
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
          c(i,j,2) = c(i,j,2)+.5*(array(i0+i*nint,j0+j*nint)+
     |      array(i0+i*nint,j0+(j+1)*nint))
          c(i,j,4) = c(i,j,4)-.5*(array(i0+i*nint,j0+j*nint)+
     |      array(i0+i*nint,j0+(j+1)*nint))
          wk(i,1) = .5*(array(i0+i*nint,j0+j*nint)+
     |      array(i0+i*nint,j0+(j+1)*nint))

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
            c(i,j,3) = c(i,j,3)+.5*(array(i0+i*nint,j0+j*nint)
     |        +array(i0+i*nint,j0+(j+1)*nint))
            c(i,j,7) = c(i,j,7)+.5*(array(i0+i*nint,j0+j*nint)
     |        +array(i0+i*nint,j0+(j-1)*nint))
            c(i,j,9) = c(i,j,9)-.5*(array(i0+i*nint,j0+(j-1)*nint)
     |        +2.*array(i0+i*nint,j0+j*nint)
     |        +array(i0+i*nint,j0+(j+1)*nint))
!
! Unmodified:
            cofum(i,j,3) = cofum(i,j,3)+.5*(array(i0+i*nint,j0+j*nint)
     |        +array(i0+i*nint,j0+(j+1)*nint))
            cofum(i,j,7) = cofum(i,j,7)+.5*(array(i0+i*nint,j0+j*nint)
     |        +array(i0+i*nint,j0+(j-1)*nint))
            cofum(i,j,9) = cofum(i,j,9)-.5*(array(i0+i*nint,j0+(j-1)*
     |        nint)+2.*array(i0+i*nint,j0+j*nint)+
     |        array(i0+i*nint,j0+(j+1)*nint))

          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
        j = 1
        do i=1,nlon
          c(i,j,3) = c(i,j,3)+.5*(array(i0+i*nint,j0+j*nint)
     |      +array(i0+i*nint,j0+(j+1)*nint))
          c(i,j,9) = c(i,j,9)-.5*(array(i0+i*nint,j0+j*nint)+
     |                            array(i0+i*nint,j0+(j+1)*nint))
          cofum(i,j,3) = cofum(i,j,3)+.5*(array(i0+i*nint,j0+j*nint)+
     |                                  array(i0+i*nint,j0+(j+1)*nint))
          cofum(i,j,9) = cofum(i,j,9)-.5*(array(i0+i*nint,j0+j*nint)+
     |                                  array(i0+i*nint,j0+(j+1)*nint))

        enddo ! i=1,nlon
      endif ! ncoef
      end subroutine cnmmod
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
      if (nint==1.and.igrid > 0) then
        do i = 1,nlon
          c(i,1,10) = c(i,1,10)/cs(1)
        enddo ! i = 1,nlon
      endif
      end subroutine divide
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
      subroutine ceee(cee,nx,ny,cf)
!
! Called from mudpack solvers to transfer coefficients.
!
      save
      integer,intent(in) :: nx,ny
      real,intent(in) :: cee(nx,ny,1)
      real,intent(out) :: cf(nx,ny,1)
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
     |  rcos0s,re,dt1dts
      use fields_module,only: phim3d ! phim3d(nmlonp1,nmlat,-2:nlevp1)
!
! Calculate 3-d dynamo potential phim3d in geomagnetic coordinates from
!   2-d potential phim output by the PDE solver mudpack.
! On input, phim(nmlonp1,nmlat) is defined from mud.
!
! Local:
      integer :: k,i,j,ip1f,ip2f,ip3f
      real :: cosltm,pi,phimn,phims,csth0,sym,ctm
      real,dimension(nmlonp1) :: thetam,pslot,qslot
      integer,dimension(nmlonp1) :: islot,jslot
      real,parameter :: unitv(nmlon)=1.
      real,dimension(nmlonp1,nmlat) :: ed1,ed2
      real,dimension(nmlonp1,nmlat,-2:nlevp1) :: ed13d,ed23d
!
! Externals:
      real,external :: sddot ! in util.F
!
      pi = pi_dyn
!
! Calculate ed1,ed2 components of electric field:
      do j=2,nmlat-1
        csth0 = cos(-pi/2.+(j-1)*dlatm)
        do i=2,nmlon
          ed1(i,j) = -(phim(i+1,j)-phim(i-1,j))/(2.*dlonm*csth0)*
     |      rcos0s(j)/(r0*1.e-2)
        enddo ! i=2,nmlon
        ed1(1,j) = -(phim(2,j)-phim(nmlon,j))/(2.*dlonm*csth0)*
     |      rcos0s(j)/(r0*1.e-2)
        ed1(nmlonp1,j) = ed1(1,j)
      enddo ! j=2,nmlat-1
!
! Southern hemisphere:
      do j=2,nmlath-1
        do i=1,nmlonp1
          ed2(i,j) = -(phim(i,j+1)-phim(i,j-1))/(2.*dlatm)*dt1dts(j)/
     |      (r0*1.e-2)
        enddo ! i=1,nmlonp1
      enddo ! j=2,nmlath-1
!
! Northern hemisphere:
      do j=nmlath+1,nmlat-1
        do i=1,nmlonp1
          ed2(i,j) = (phim(i,j+1)-phim(i,j-1))/(2.*dlatm)*dt1dts(j)/
     |      (r0*1.e-2)
        enddo ! i=1,nmlonp1
      enddo ! j=nmlath+1,nmlat-1
!
! Poles:
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
!
! Equator
       ed2(i,nmlath) = (4.*phim(i,nmlath+1)-phim(i,nmlath+2)
     |    -3.*phim(i,nmlath))/(2.*dlatm)/(R0*1.e-2)
       ed2(i,nmlath+1) = (4.*phim(i,nmlath+2)-phim(i,nmlath+3)
     |    -3.*phim(i,nmlath+1))/(2.*dlatm)/(R0*1.e-2)
       ed2(i,nmlath-1) = (4.*phim(i,nmlath-2)-phim(i,nmlath-3)
     |    -3.*phim(i,nmlath-1))/(2.*dlatm)/(R0*1.e-2)
      enddo ! i = 1,nmlonp1
!
! Comment from tgcm15 (threed.F):
! For each geomagnetic grid point over the 3d grid, determine the
!   location of the foot of the magnetic field passing through the
!   grid point. Assume that the field lines are dipole in shape. 
!   Recalling that the model assumes constant electric potential
!   along each field line, interpolation in the 2d array phim gives
!   the potential at the original grid point.
!
      do k=-2,nlevp1
        do j=2,nmlat-1
          cosltm = cos(ylatm(j))
          sym = 1.                   ! sym in new dynamo only
          if (j < nmlath) sym = -1.
          do i=1,nmlon
!
! Thetam(i) = thetas for each value of i:
! (note old dynamo incorrectly uses r00 instead of re here) 
!
            ctm = sqrt(r0/(re+zpotenm3d(i,j,k)))*cosltm
            thetam(i) = acos(ctm)
!
! Find position of thetam in table of theta0 vs thetas:
            pslot(i) = thetam(i)*180./pi+1.
            islot(i) = pslot(i)
            pslot(i) = pslot(i)-float(islot(i))
!
! Determine theta0 by interpolation in table:
! bf 9/7/02 bugfix: added parens at beginning of statement and 
!    before *sym
!
            thetam(i) = ((1.-pslot(i))*table(islot(i),2)+pslot(i)*
     |        table(islot(i)+1,2))*sym ! thetam negative for south hem
!
! Now identify geomagnetic grid element in which the point lies in phim:
            islot(i) = i
            pslot(i) = 0.
            qslot(i) = (thetam(i)+pi/2.)/dlatm+1.
            jslot(i) = qslot(i)
            qslot(i) = qslot(i)-float(jslot(i))
!
! Interpolate in phim using bilinear interpolation:
            phim3d(i,j,k) = (1.-qslot(i))*((1.-pslot(i))*
     |        phim(islot(i)  ,jslot(i))  +pslot(i)*
     |        phim(islot(i)+1,jslot(i))) +qslot(i)*((1.-pslot(i))*
     |        phim(islot(i)  ,jslot(i)+1)+pslot(i)*
     |        phim(islot(i)+1,jslot(i)+1))
!
            ed13d(i,j,k) = (1.-qslot(i))*((1.-pslot(i))*
     |        ed1(islot(i),jslot(i))+pslot(i)*
     |        ed1(islot(i)+1,jslot(i)))+qslot(i)*((1.-pslot(i))*
     |        ed1(islot(i),jslot(i)+1)+pslot(i)*
     |        ed1(islot(i)+1,jslot(i)+1))
            ed23d(i,j,k) = (1.-qslot(i))*((1.-pslot(i))*
     |        ed2(islot(i),jslot(i))+pslot(i)*
     |        ed2(islot(i)+1,jslot(i)))+qslot(i)*((1.-pslot(i))*
     |        ed2(islot(i),jslot(i)+1)+pslot(i)*
     |        ed2(islot(i)+1,jslot(i)+1))
          enddo ! i=1,nmlon
        enddo ! j=2,nmlat-1
      enddo ! k=-2,nlevp1
!
! Calculate value of phim3d at geomagnetic poles:
      phims = sddot(nmlon,unitv,phim(1,1))/float(nmlon)
      phimn = sddot(nmlon,unitv,phim(1,nmlat))/float(nmlon)
      do k = -2,nlevp1
        do i = 1,nmlon
          phim3d(i,1,k) = phims
          phim3d(i,nmlat,k) = phimn
!
          ed13d(i,1,k)     = ed1(i,1)
          ed13d(i,nmlat,k) = ed1(i,nmlat)
          ed23d(i,1,k)     = ed2(i,1)
          ed23d(i,nmlat,k) = ed2(i,nmlat)
        enddo ! i = 1,nmlon
      enddo ! k = -2,nlevp1
!
! Periodic point:
      do k=1,nlevp1
        do j=1,nmlat
          phim3d(nmlonp1,j,k) = phim3d(1,j,k)
          ed13d(nmlonp1,j,k)  = ed13d(1,j,k)
          ed23d(nmlonp1,j,k)  = ed23d(1,j,k)
        enddo ! j=1,nmlat
      enddo ! k=1,nlevp1
!
! Save 3d electric potential on magnetic grid to secondary history:
!     do j=1,nmlat
!       call addfsech_ik('PHIM3D_IK','PHIM3D_IK','[V]',phim3d(:,j,:),
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!
! Save electric field diagnostic:
!       call addfsech_ik('ED1M3D'  ,' ',' ',ed13d(:,j,:),
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!       call addfsech_ik('ED2M3D'  ,' ',' ',ed23d(:,j,:),
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!     enddo ! j=1,nmlat
!     call addfsech_ij('PHIM3D_IJ','PHIM3D_IJ','[V]',phim3d(:,:,1),
!    |   1,nmlonp1,1,nmlat)
      end subroutine threed
!-----------------------------------------------------------------------
      subroutine mag2geo(am,ag,im,jm,dim,djm,lg,lm,nlong,nlatg,nlonm,
     |  nlatm)
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
! This routine is similiar to sub geo2mag in dynamo.F.
!
! Args:
      integer,intent(in) :: nlonp1_geo,nlonp1_mag,nlon_mag,nlat_mag,lat
      integer,dimension(nlonp1_mag,nlat_mag),intent(in) :: long,latg
      real,intent(in) :: fgeo(nlonp1_geo,*),wght(4,nlonp1_mag,nlat_mag)
      real,intent(out) :: fmag(nlonp1_mag,*)
!     integer,intent(in) :: iprint
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
      subroutine mp_dynamo_gather
      use mpi_module,only: mytid,lon0,lon1,lat0,lat1,mxlat,mxlon,
     |  ntask,mpi_comm_world,mpi_real8,irstat,tasks,handle_mpi_err
!
! Local:
      integer,parameter :: nfld=7
      integer :: k,i,j,n,nn,ier,len,nlons,nlats,ifld,msgtag,ireqrecv,
     |  ireqsend,ilon0,ilon1,ilat0,ilat1
      integer :: idest = 0
      real :: fmin,fmax
      real,allocatable :: 
     |  rcvbuf(:,:,:,:,:), ! Receive buffer for root task (i,j,k,ifld,ntask)
     |  sndbuf(:,:,:,:)    ! Send buffer for slave tasks (i,j,k,ifld)
      character(len=8) :: fnames(nfld) =
     |  (/'SIGMAPED','SIGMAHAL','ZPOTEN  ','SCHEIGHT','UNVEL   ',
     |    'VNVEL   ','WNVEL   '/)
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
!-----------------------------------------------------------------------
      end module dynamo_module
