!
      module dynamo_old_module
!
! Module for electrodynamic dynamo (the "e" in tiegcm). 
!
! 4/02 btf: This is a rewrite of the dynamo in tgcm15, as originally 
! written by E. Cicely Ridley, probably for tgcm12, completed in
! about 1992. This does not include new dynamo code written by
! Astrid Maute and Art Richmond during 2001-2002.
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
     |  nlon,   ! number of geographic longitudes (e.g. at 5 deg, nlon==72)
     |  nlonp1, ! nlon+1
     |  nlonp4, ! nlon+4
     |  nlev,   ! number of pressure levels (e.g., 28)
     |  nlevp1, ! nlev+1
     |  nlat,   ! number of geographic latitudes (e.g. at 5 deg, nlat==36)
     |  nlatp1, ! nlat+1
     |  nlatp2, ! nlat+2
     |  nmlon,  ! number of geomagnetic grid longitudes
     |  nmlonp1,! nmlon+1
     |  nmlat,  ! number of geomagnetic grid latitudes
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
! Coefficients are stored in 1-d array cee0(ncee)
! cee transmits descretized dynamo PDE coefficients to the multi-grid 
!   mudpack solver. (cee was formerly in ceee.h)
! The common block /cee_com/ is retained from earlier versions because
!   of the equivalencing below of coefficient arrays c0, c1, etc.
!
      real :: cee0(ncee) 
      common/cee_com/ cee0
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
! nc(6) = ncee, the dimension of the entire cee0 array, containing
!   coefficients for all 5 levels of resolution.
! 
      integer :: nc_old(6)
!
      real :: 
     |  c0(nmlon0,nmlat0,10),
     |  c1(nmlon1,nmlat1,9),
     |  c2(nmlon2,nmlat2,9),
     |  c3(nmlon3,nmlat3,9),
     |  c4(nmlon4,nmlat4,9)
      equivalence 
     |  (cee0,c0), 
     |  (cee0(nc1),c1),
     |  (cee0(nc2),c2),
     |  (cee0(nc3),c3),
     |  (cee0(nc4),c4)
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
      real,dimension(nmlonp1,-2:nlev,2) :: adotvm,axvm
      real,dimension(nmlonp1,2)         :: bxam,adotam,azm 
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
! For dot products:
      real,parameter :: unitvm(nmlon)=1., unitv(nlon)=1.
!
! Electric potential from heelis or weimer:
      real,dimension(nmlonp1,nmlat0) ::
     |   pfrac, ! NH fraction of potential
     |   phihm  ! potential in magnetic
      contains
!-----------------------------------------------------------------------
      subroutine dynamo_old
!
! Driver for dynamo code. This is called from advance only if the
! user input parameter DYNAMO is > 0. As of 4/02, if this code is
! run without the dynamo, it will stop in settei with very high
! electron temperatures (> 5000 deg) after less than 6 hours model 
! time. This is executed by master task only. 
!
      use cons_module,only: pi_dyn,dlonm,dlatm,ylatm,rtd
      use magfield_module,only: im,jm,dim,djm
      use fields_module,only: dynpot,phim3d
      use input_module,only: potential_model
!
! Local:
      real :: array(-15:nmlon0+16,nmlat0),cs(nmlat0)
      integer :: n,ncc,nmaglon,nmaglat
      integer :: i,j,jj,jjj,j0,jntl,k
      real :: sym,wk(nmlon0,3)
      real,dimension(nmlonp1,nlevp1) :: phim_plt, phihm_plt
      real,dimension(nlonp1,0:nlatp1) :: phih   ! potential in geographic
!
!     write(6,"('Enter dynamo.')")
!
! Transform needed fields to geomagnetic coordinates, perform field-line
! integrations, and evaluate PDE coefficients and RHS.
!
      call transf
!
! Set index array nc and magnetic latitude cosine array:
      nc_old(1) = nc0
      nc_old(2) = nc1
      nc_old(3) = nc2
      nc_old(4) = nc3
      nc_old(5) = nc4
      nc_old(6) = ncee
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
! Set zigm11 to zero at megnetic poles to avoid floating exception 
! (values at poles are not used):
!
      do i = 1,nmlonp1
        zigm11(i,1) = 0.0
        zigm11(i,nmlat) = 0.0
      enddo
!
! Clear array for difference stencils at all levels:
      call clearcee(cee0,nmlon0,nmlat0)
!
! Calculate contribution to stencils from each PDE coefficient:
      array(1,1) = 1.
      sym = 1.
      call stencil(zigm11,nmlon0,nmlat0,sym,cee0,array(1,1),1)
!
      array(1,1) = 4.
      sym = 1.
      call stencil(zigm22,nmlon0,nmlat0,sym,cee0,array(1,1),4)
!
      array(1,1) = 2.
      sym = -1.
      call stencil(zigmc,nmlon0,nmlat0,sym,cee0,array(1,1),2)
!
      array(1,1) = 3.
      sym = -1.
      call stencil(zigm2,nmlon0,nmlat0,sym,cee0,array(1,1),3)
!
! Insert RHS in finest stencil (formerly sub rths):
      do j = 1,nmlat0
        jj = nmlath-nmlat0+j
        do i = 1,nmlon0
          c0(i,j,10) = rhs(i,jj)
        enddo ! i = 1,nmlon0
      enddo ! j = 1,nmlat0
!
! Set boundary condition at equator and pole:
      call edges(c0,nmlon0,nmlat0)
      call edges(c1,nmlon1,nmlat1)
      call edges(c2,nmlon2,nmlat2)
      call edges(c3,nmlon3,nmlat3)
      call edges(c4,nmlon4,nmlat4)
!
! Divide stencils by cos(theta):
      call divide(c0,nmlon0,nmlat0,nmlon0,nmlat0,cs)
      call divide(c1,nmlon1,nmlat1,nmlon0,nmlat0,cs)
      call divide(c2,nmlon2,nmlat2,nmlon0,nmlat0,cs)
      call divide(c3,nmlon3,nmlat3,nmlon0,nmlat0,cs)
      call divide(c4,nmlon4,nmlat4,nmlon0,nmlat0,cs)
!
! Set value of solution to 1. at pole:
      do i=1,nmlon0
        c0(i,nmlat0,10) = 1.
      enddo
!
! Save single-level heelis potential:
!     do j=1,nmlat
!       do i=1,nmlonp1
!         phihm_plt(i,:) = phihm(i,j)
!       enddo ! i=1,nmlonp1
!       call addfsech_ik('PHIHM'  ,' ',' ',phihm_plt,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!     enddo ! j=1,nmlat
!
! Modify stencils and RHS so that heelis potential is inserted at
! high latitude. The SH high lat potential will be added back later.
!  pfrac = fraction of dynamo in solution in the NH. = 1 low lat, = 0 hi lat
!    cons_module: crit(1)=15, crit(2)=30 deg colats, or hi-lat > 75 deg,
!      dynamo < 60 deg, and combination between 60-75 mag lat.
! The dynamo is symmetric about the magnetic equator, but the high latitude
!  is anti-symmetric in both hemispheres.  However, since Mudpack uses the
!  NH potential pattern, then the SH potential pattern must be added
!  back into the 2-D phim before the call threed, and before it is
!  transformed to geographic coordinates.
!
      ncc = 1
      nmaglon = nmlon0
      nmaglat = nmlat0
      do n=1,5
        call stenmod(nmaglon,nmaglat,cee0(ncc),phihm(1,nmlat0),pfrac)
        ncc = ncc+9*nmaglon*nmaglat
        if (n==1) ncc = ncc+nmaglon*nmaglat ! rhs is in 10th slot
        nmaglon = (nmaglon+1)/2
        nmaglat = (nmaglat+1)/2
      enddo ! n=1,5
!
! Call mudpack to solve PDE. Solution is returned in rim:
! (3rd arg to mud indicates to mud to use cee0 and nc_old)
      jntl = 0
      call mud(rim,jntl,-1)
!
! Copy output potential from rim to phim(nmlonp1,nmlat):
!  Correct the SH potential for the anti-symmetric imposed NH high lat poten
!  pfrac = fraction of dynamo in solution in the NH. = 1 low lat, = 0 hi lat
!    cons_module: crit(1)=15, crit(2)=30 deg colats, or hi-lat > 75 deg,
!      dynamo < 60 deg, and combination between 60-75 mag lat.
!  ie, need (1.-pfrac(i,nmlat0)) at phim(i,1), etc
! jn index for NH part of potential (nmlat down to ~nmlat0)
! jp index for NH pfrac (nmlat0 down to 1)
      do j=1,nmlat0
        jn = nmlat - j + 1
        jp = nmlat0 - j + 1
        do i=1,nmlonp1
          phim(i,j)=rim(i,j,1)+(1.-pfrac(i,jp))*(phihm(i,j)-phihm(i,jn))
        enddo ! i=1,nmlonp1
!       write(6,"('dynamo: j=',3i3,' rim phim s,n,frac=',/,(5e12.4))")
!    |    j,jn,jp,(rim(i,j,1),phim(i,j),phihm(i,j),phihm(i,jn),
!    |    pfrac(i,jp),i=1,nmlonp1,20)
      enddo ! j=1,nmlat0
      do j=nmlat0+1,nmlat
        jp = j - nmlat0
        do i=1,nmlonp1
          phim(i,j) = rim(i,j,1)
        enddo ! i=1,nmlonp1
!       write(6,"('dynamo: j=',2i3,' phim n,frac=',/,(6e12.4))")
!    |    j,jp,(phim(i,j),phihm(i,j),pfrac(i,jp),i=1,nmlonp1,20)
      enddo ! j=1,nmlat
!
! Save correct single-level potential to secondary histories (mag grid):
!     do j=1,nmlat
!       do i=1,nmlonp1
!         phim_plt(i,:) = phimsn(i,j)
!       enddo ! i=1,nmlonp1
!       call addfsech_ik('PHIM'  ,' ',' ',phim_plt,
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!     enddo ! j=1,nmlat
!
! Call threed to generate 3-d potential array in geomagnetic coordinates
!   from 2-d solver output phim, corrected for the SH potential.
!   phim3d(nmlonp1,nmlat,-2:nlevp1) is in fields.F.
!
      call threed
!
! Transform phim3d to geographic coorinates in dynpot (fields.F):
!   phim3d(nmlonp1,nmlat,-2:nlevp1)  ! 3d electric potential magnetic
!   dynpot(nlonp1,0:nlatp1,nlevp1),  ! 3d electric potential geographic
!
!     subroutine mag2geo(am,ag,im,jm,dim,djm,lg,lm,nlong,nlatg,nlonm,
!    |  nlatm)
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
!     do j=1,nlat
!       call addfsech_ik('DYNPOT',' ',' ',dynpot(:,j,:),
!    |    1,nlonp1,nlevp1,nlevp1-1,j)
!     enddo ! j=1,nlat
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
      end subroutine dynamo_old
!-----------------------------------------------------------------------
!     subroutine dynpot2mag
!
! 9/12/02: this did not work because of circular dependency when this
!          routine was use-associated into nchist, so will just do it
!          at end of dynamo above.
!
!     use fields_module,only: dynpot,phim3d
!     use magfield_module,only: im,jm,dim,djm
!
! Transform phim3d to geographic coorinates in dynpot (fields.F) prior
!   to writing a history. This is called from wrf4d (nchist.F).
!   phim3d(nmlonp1,nmlat,-2:nlevp1)  ! 3d electric potential magnetic
!   dynpot(nlonp1,0:nlatp1,nlevp1),  ! 3d electric potential geographic
!
!     subroutine mag2geo(am,ag,im,jm,dim,djm,lg,lm,nlong,nlatg,nlonm,
!    |  nlatm)
! Local:
!     integer :: k,j
!
!     do k=1,nlevp1
!       call mag2geo(phim3d(1,1,k),dynpot(1,0,k),im(1,0),jm(1,0),
!    |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
!     enddo ! k=1,nlevp1
!
! Periodic point:
!     do k=1,nlevp1
!       do j=0,nlatp1
!         dynpot(nlonp1,j,k) = dynpot(1,j,k)
!       enddo ! j=0,nlatp1
!     enddo ! k=1,nlevp1
!     end subroutine dynpot2mag
!-----------------------------------------------------------------------
      subroutine divide(c,nlon,nlat,nlon0,nlat0,cs)
!
! Divide stencil C by cos(theta(i,j))
!
! Args:
      integer,intent(in) :: nlon,nlat,nlon0,nlat0
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
      if (nint==1) then
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
      if (ncoef==1) then
        do j = 2,nlat-1
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
      subroutine prep_dynamo_old(tn,un,vn,w,z,barm,ped,hall,
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

      enddo ! lat=lat0,lat1
#ifdef MPI
!
! Gather dynamo input fields to the root task, defining module data 
! above at the global domain on the root task:
!
      call mp_dynamo_gather
#endif
      end subroutine prep_dynamo_old
!-----------------------------------------------------------------------
      subroutine geo2mag(fg,nglon,nglat, fm,nmlon,nmlat, ig,jg,wt,
     |  im0,im1,k0,k1,lat, iprint)
!
! Transform given field fg on the geographic grid to geomagnetic space.
! Field on the magnetic grid is returned in fm. This is called from a
!   magnetic lat loop in transf.
! Note this routine is different than the geo2mag in dynamo.F.
!
! Args:
      integer,intent(in) :: 
     |  nglon,nglat,     ! number of geographic lons,lats
     |  nmlon,nmlat,     ! number of magnetic lons,lats
     |  im0,im1,         ! beginning and ending longitude indices
     |  k0,k1,           ! beginning and ending pressure indices
     |  lat,             ! current magnetic latitude index
     |  ig(nmlon,nmlat), ! geog lon index at each geomag grid point
     |  jg(nmlon,nmlat)  ! geog lat index at each geomag grid point
      real,intent(in)  :: 
     |  wt(4,nmlon,nmlat),    ! weights at grid cell corners
     |  fg(nglon,nglat,k0:k1) ! 3d input field at geog grid
      real,intent(out) :: 
     |  fm(nmlon,k0:k1)       ! 2d output field at current maglat
      integer,optional :: iprint
!
! Local: 
      integer :: i,k
      real :: fgtmp(nglon,nglat,k0:k1) ! for debug
!
      do k=k0,k1
        do i=im0,im1
          fm(i,k) = 
     |      fg(ig(i,lat)  ,jg(i,lat)  ,k)*wt(1,i,lat)+
     |      fg(ig(i,lat)+1,jg(i,lat)  ,k)*wt(2,i,lat)+
     |      fg(ig(i,lat)+1,jg(i,lat)+1,k)*wt(3,i,lat)+
     |      fg(ig(i,lat)  ,jg(i,lat)+1,k)*wt(4,i,lat)
!
          if (present(iprint)) then
            if (iprint > 0) then
              write(6,"('geo2mag: i=',i3,' k=',i3,
     |          ' j=',i3,' ig=',i3,' jg=',i3,' fm=',e12.4)") i,k,lat,
     |          ig(i,lat),jg(i,lat),fm(i,k)
              write(6,"('fg(ig,jg)=',e12.4,' fg(ig+1,jg)=',e12.4,
     |          ' fg(ig+1,jg+1)=',e12.4,' fg(ig,jg+1)=',e12.4)")
     |          fg(ig(i,lat)  ,jg(i,lat)  ,k),
     |          fg(ig(i,lat)+1,jg(i,lat)  ,k),
     |          fg(ig(i,lat)+1,jg(i,lat)+1,k),
     |          fg(ig(i,lat)  ,jg(i,lat)+1,k)
            endif
          endif
        enddo ! i=im0,im1
      enddo ! k=k0,k1
      end subroutine geo2mag
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
      subroutine transf
!
! Define needed fields on the geographic grid, then transform to
! geomagnetic space (via sub geo2mag). 
!
      use cons_module,only: pi,dlatm,dlonm,r0
      use magfield_module,only: ! data read from mag data file (magfield.F)
     |  bmod,p,   ! (nlonp1,0:nlatp1)
     |  xb,yb,zb, ! (nlonp1,0:nlatp1)
     |  av,       ! (nlonp1,0:nlatp1,3,2)
     |  ig,       ! (nmlonp1,nmlat) geog lon at each geomag grid point
     |  jg,       ! (nmlonp1,nmlat) geog lat at each geomag grid point
     |  wt        ! (4,nmlonp1,nmlat) interpolation weights
!
! Local:
      integer :: k,kk,i,ii,j,jj,n,lat,iprint
      real :: fmin,fmax
      real,parameter ::
     |  rl1=5.e5,
     |  rl2=3.e5,
     |  z0=90.e5
!
! Fields to be transformed to geomagnetic space (formerly in transmag.h):
! (these will be input to the geographic to magnetic transformation). 
      real,dimension(nlonp1,0:nlatp1,-2:nlev) ::
     |  ssigma1,  ! pedersen conductivity
     |  ssigma2,  ! hall conductivity
     |  ww,       ! vertical motion (omega)
     |  vxb
      real,dimension(nlonp1,0:nlatp1,-2:nlevp1) :: zz
      real,dimension(nlonp1,0:nlatp1,-2:nlev,2) :: adotv,axv
      real,dimension(nlonp1,0:nlatp1,2)         :: bxa,adota
      real,dimension(nlonp1,0:nlatp1)           :: a1dta2,sini
!
! Local integrations:
      real,dimension(nmlat) :: tint1,tint2,tint3
      real,dimension(-1:nmlat) :: tint33
!
! Diagnostics for plotting:
      real,dimension(nmlonp1,nlevp1) ::
     |  zigm11_plt, zigmc_plt, zigm2_plt, zigm22_plt, rim1_plt, rim2_plt
!
! Externals:
      real,external :: sddot ! in util.F
!
! Pack inputs from 3->nlon+3 to 1->nlon+1, as in tgcm15:
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
! Check inputs from prep_dynamo:

!
! Polar points of magnetic quantities to be transformed to geomagnetic
! space. bmod, p, and av (excluding polar points) were read from magnetic 
! data file by sub nc_rdmag (magfield.F).
!
      bmod(1,0) = (9.*sddot(nlon,unitv,bmod(1,1))-
     |  sddot(nlon,unitv,bmod(1,2)))/(8.*float(nlon))

      bmod(1,nlatp1) = (9.*sddot(nlon,unitv,bmod(1,nlat))-
     |  sddot(nlon,unitv,bmod(1,nlat-1)))/(8.*float(nlon))

      p(1,0) = (9.*sddot(nlon,unitv,p(1,1))-
     |  sddot(nlon,unitv,p(1,2)))/(8.*float(nlon))

      p(1,nlatp1) = (9.*sddot(nlon,unitv,p(1,nlat))-
     |  sddot(nlon,unitv,p(1,nlat-1)))/(8.*float(nlon))

      av(1,0,3,1) = (9.*sddot(nlon,unitv,av(1,1,3,1))-
     |  sddot(nlon,unitv,av(1,2,3,1)))/(8.*float(nlon))

      av(1,nlatp1,3,1) = (9.*sddot(nlon,unitv,av(1,nlat,3,1))-
     |  sddot(nlon,unitv,av(1,nlat-1,3,1)))/(8.*float(nlon))

      av(1,0,3,2) = (9.*sddot(nlon,unitv,av(1,1,3,2))-
     |  sddot(nlon,unitv,av(1,2,3,2)))/(8.*float(nlon))

      av(1,nlatp1,3,2) = (9.*sddot(nlon,unitv,av(1,nlat,3,2))-
     |  sddot(nlon,unitv,av(1,nlat-1,3,2)))/(8.*float(nlon))

      do i=2,nlonp1
        bmod(i,0) = bmod(1,0)
        bmod(i,nlatp1) = bmod(1,nlatp1)
        p(i,0) = p(1,0)
        p(i,nlatp1) = p(1,nlatp1)
        av(i,0,3,1) = av(1,0,3,1)
        av(i,nlatp1,3,1) = av(1,nlatp1,3,1)
        av(i,0,3,2) = av(1,0,3,2)
        av(i,nlatp1,3,2) = av(1,nlatp1,3,2)
      enddo ! i=2,nlonp1
!
! Calculate quantities to be transformed to geomagnetic space:
      do k=1,nlev        
        do j=1,nlat
          do i=1,nlonp1
            ssigma1(i,j,k) = sigma_ped (i,j,k)
            ssigma2(i,j,k) = sigma_hall(i,j,k)
            zz(i,j,k) = zpoten(i,j,k)
            ww(i,j,k) = wnvel (i,j,k)
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
      do k=-2,0
        do j=1,nlat
          do i=1,nlonp1
            zz(i,j,k) = z0+float(k+2)*(zz(i,j,1)-z0)/3.
            ww(i,j,k) = ww(i,j,1)
          enddo ! i=1,nlonp1
        enddo ! j=1,nlat
      enddo ! k=-2,0
      do k=-2,0
        do j=1,nlat
          do i=1,nlonp1
            ssigma1(i,j,k) = ssigma1(i,j,1)*exp((zz(i,j,k)+zz(i,j,k+1)-
     |        zz(i,j,1)-zz(i,j,2))/(2.*rl1))
            ssigma2(i,j,k) = ssigma2(i,j,1)*exp((zz(i,j,k)+zz(i,j,k+1)-
     |        zz(i,j,1)-zz(i,j,2))/(2.*rl2))
          enddo ! i=1,nlonp1
        enddo ! j=1,nlat
      enddo ! k=-2,0
!
! ADOTV = AV(X)*U + AV(Y)*V
! V X A = U*AV(Y) - V*AV(X)
!
      do n=1,2
        do k=-2,nlev
          kk = k
          if (kk < 1) kk = 1
          do j=1,nlat
            do i=1,nlonp1
              adotv(i,j,k,n) = av(i,j,1,n)*unvel(i,j,kk)+
     |                         av(i,j,2,n)*vnvel(i,j,kk)
              axv(i,j,k,n) =   av(i,j,1,n)*vnvel(i,j,kk)-
     |                         av(i,j,2,n)*unvel(i,j,kk)
            enddo ! i=1,nlonp1
          enddo ! j=1,nlat
        enddo ! k=-2,nlev
      enddo ! n=1,2
!
! 2-d fields:
! TINT1 = SIN(D),   TINT2 = COS(D)
! VXB = U*COS(D)-V*SIN(D)
! BXA = B(X)*AV(Y) - B(Y)*AV(X)
! ADOTA = AV(X)**2 + AV(Y)**2
!
! xb,yb were read from magfield data file (magfield module).
      do j=1,nlat                                              ! do 31
        do i=1,nlonp1                                          ! do 32
          tint1(i) = sqrt(xb(i,j)**2+yb(i,j)**2) 
          tint2(i) = xb(i,j)/tint1(i)
          tint1(i) = yb(i,j)/tint1(i)
        enddo ! i=1,nlonp1
        do k=-2,nlev                                           ! do 34
          kk = k
          if (kk < 1) kk = 1 
          do i=1,nlonp1                                        ! do 34
            vxb(i,j,k) = unvel(i,j,kk)*tint2(i)-
     |                   vnvel(i,j,kk)*tint1(i)
          enddo ! i=1,nlonp1
        enddo ! k=-2,nlev
        do n=1,2                                               ! do 31
          do i=1,nlonp1
            bxa(i,j,n) = tint1(i)*av(i,j,2,n)-
     |                   tint2(i)*av(i,j,1,n) 
            adota(i,j,n) = av(i,j,1,n)**2+av(i,j,2,n)**2
          enddo ! i=1,nlonp1
        enddo ! n=1,2
      enddo ! j=1,nlat
!
! A1DTA2 = AV(X,1)*AV(X,2) + AV(Y,1)*AV(Y,2)
! SINI = ZB/BMOD
!
      do j=1,nlat                                              ! do 33
        do i=1,nlonp1
          a1dta2(i,j) = av(i,j,1,1)*av(i,j,1,2)+
     |                  av(i,j,2,1)*av(i,j,2,2)
          sini(i,j) = zb(i,j)/bmod(i,j)
        enddo ! i=1,nlonp1
      enddo ! j=1,nlat
!
! Calculate values at geographic poles:
! (in earlier versions, ssigma1 was over indexed through the common 
!  for all 11 fields)
!
! Dimensions from above, for reference:
!     real,dimension(nlonp1,0:nlatp1,-2:nlev) :: 
!    |  ssigma1, ssigma2, ww, vxb 
!     real,dimension(nlonp1,0:nlatp1,-2:nlevp1) :: zz
!     real,dimension(nlonp1,0:nlatp1,-2:nlev,2) :: adotv,axv
!
      do k=-2,nlev                                             ! do 5
        ssigma1(1,0,k) = 
     |    (9.*sddot(nlon,unitv,ssigma1(1,1,k))-
     |        sddot(nlon,unitv,ssigma1(1,2,k)))/(8.*float(nlon))
        ssigma1(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,ssigma1(1,nlat  ,k))-
     |        sddot(nlon,unitv,ssigma1(1,nlat-1,k)))/(8.*float(nlon))

!       write(6,"(/,'k=',i3,' nlon=',i3,' ssigma1(1,1,k)=',e12.4,
!    |    ' (1,2,k)=',e12.4)") k,nlon,ssigma1(1,1,k),ssigma1(1,2,k)
!       write(6,"('k=',i3,' nlat=',i3,' ssigma1(1,nlat,k)=',e12.4,
!    |    ' (1,nlat-1,k)=',e12.4)") k,nlat,
!    |    ssigma1(1,nlat,k),ssigma1(1,nlat-1,k)
!       write(6,"('ssigma1(1,0,k)=',e12.4,' ssigma1(1,nlatp1,k)=',
!    |    e12.4)") ssigma1(1,0,k),ssigma1(1,nlatp1,k)

        ssigma2(1,0,k) = 
     |    (9.*sddot(nlon,unitv,ssigma2(1,1,k))-
     |        sddot(nlon,unitv,ssigma2(1,2,k)))/(8.*float(nlon))
        ssigma2(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,ssigma2(1,nlat  ,k))-
     |        sddot(nlon,unitv,ssigma2(1,nlat-1,k)))/(8.*float(nlon))

        ww(1,0,k) = 
     |    (9.*sddot(nlon,unitv,ww(1,1,k))-
     |        sddot(nlon,unitv,ww(1,2,k)))/(8.*float(nlon))
        ww(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,ww(1,nlat  ,k))-
     |        sddot(nlon,unitv,ww(1,nlat-1,k)))/(8.*float(nlon))

        zz(1,0,k) = 
     |    (9.*sddot(nlon,unitv,zz(1,1,k))-
     |        sddot(nlon,unitv,zz(1,2,k)))/(8.*float(nlon))
        zz(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,zz(1,nlat  ,k))-
     |        sddot(nlon,unitv,zz(1,nlat-1,k)))/(8.*float(nlon))

        vxb(1,0,k) = 
     |    (9.*sddot(nlon,unitv,vxb(1,1,k))-
     |        sddot(nlon,unitv,vxb(1,2,k)))/(8.*float(nlon))
        vxb(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,vxb(1,nlat  ,k))-
     |        sddot(nlon,unitv,vxb(1,nlat-1,k)))/(8.*float(nlon))

        do n=1,2
          adotv(1,0,k,n) = 
     |      (9.*sddot(nlon,unitv,adotv(1,1,k,n))-
     |          sddot(nlon,unitv,adotv(1,2,k,n)))/(8.*float(nlon))
          adotv(1,nlatp1,k,n) = 
     |      (9.*sddot(nlon,unitv,adotv(1,nlat  ,k,n))-
     |          sddot(nlon,unitv,adotv(1,nlat-1,k,n)))/(8.*float(nlon))

          axv(1,0,k,n) = 
     |      (9.*sddot(nlon,unitv,axv(1,1,k,n))-
     |          sddot(nlon,unitv,axv(1,2,k,n)))/(8.*float(nlon))
          axv(1,nlatp1,k,n) = 
     |      (9.*sddot(nlon,unitv,axv(1,nlat  ,k,n))-
     |          sddot(nlon,unitv,axv(1,nlat-1,k,n)))/(8.*float(nlon))
        enddo ! n=1,2
!
        do i=2,nlon  
          ssigma1(i,0,k) = ssigma1(1,0,k)
          ssigma1(i,nlatp1,k) = ssigma1(1,nlatp1,k)

          ssigma2(i,0,k) = ssigma2(1,0,k)
          ssigma2(i,nlatp1,k) = ssigma2(1,nlatp1,k)

          ww(i,0,k) = ww(1,0,k)
          ww(i,nlatp1,k) = ww(1,nlatp1,k)

          zz(i,0,k) = zz(1,0,k)
          zz(i,nlatp1,k) = zz(1,nlatp1,k)

          vxb(i,0,k) = vxb(1,0,k)
          vxb(i,nlatp1,k) = vxb(1,nlatp1,k)

          adotv(i,0,k,:) = adotv(1,0,k,:)
          adotv(i,nlatp1,k,:) = adotv(1,nlatp1,k,:)

          axv(i,0,k,:) = axv(1,0,k,:)
          axv(i,nlatp1,k,:) = axv(1,nlatp1,k,:)

        enddo ! i=2,nlon  
      enddo ! k=-2,nlev
!
! sini(nlonp1,0:nlatp1), a1dta2(nlonp1,0:nlatp1) (no k)
      sini(1,0) = 
     |  (9.*sddot(nlon,unitv,sini(1,1))-
     |      sddot(nlon,unitv,sini(1,2)))/(8.*float(nlon))
      sini(1,nlatp1) = 
     |  (9.*sddot(nlon,unitv,sini(1,nlat  ))-
     |      sddot(nlon,unitv,sini(1,nlat-1)))/(8.*float(nlon))
      a1dta2(1,0) = 
     |  (9.*sddot(nlon,unitv,a1dta2(1,1))-
     |      sddot(nlon,unitv,a1dta2(1,2)))/(8.*float(nlon))
      a1dta2(1,nlatp1) = 
     |  (9.*sddot(nlon,unitv,a1dta2(1,nlat  ))-
     |      sddot(nlon,unitv,a1dta2(1,nlat-1)))/(8.*float(nlon))
      do i=2,nlon  
        sini(i,0) = sini(1,0)
        sini(i,nlatp1) = sini(1,nlatp1)
        a1dta2(i,0) = a1dta2(1,0)
        a1dta2(i,nlatp1) = a1dta2(1,nlatp1)
      enddo
!
! Dimensions from above, for reference:
!     real,dimension(nlonp1,0:nlatp1,2)         :: bxa,adota
!
      do n=1,2
        bxa(1,0,n) = 
     |    (9.*sddot(nlon,unitv,bxa(1,1,n))-
     |        sddot(nlon,unitv,bxa(1,2,n)))/(8.*float(nlon))
        bxa(1,nlatp1,n) = 
     |    (9.*sddot(nlon,unitv,bxa(1,nlat  ,n))-
     |        sddot(nlon,unitv,bxa(1,nlat-1,n)))/(8.*float(nlon))

        adota(1,0,n) = 
     |    (9.*sddot(nlon,unitv,adota(1,1,n))-
     |        sddot(nlon,unitv,adota(1,2,n)))/(8.*float(nlon))
        adota(1,nlatp1,n) = 
     |    (9.*sddot(nlon,unitv,adota(1,nlat  ,n))-
     |        sddot(nlon,unitv,adota(1,nlat-1,n)))/(8.*float(nlon))
!
        do i=2,nlon
          bxa(i,0,n) = bxa(1,0,n)
          bxa(i,nlatp1,n) = bxa(1,nlatp1,n)

          adota(i,0,n) = adota(1,0,n)
          adota(i,nlatp1,n) = adota(1,nlatp1,n)
        enddo ! i=2,nlon
      enddo ! n=1,2
!
! Periodic points:
!
      ssigma1(nlonp1,:,:)   = ssigma1(1,:,:)
      ssigma2(nlonp1,:,:)   = ssigma2(1,:,:)
      ww     (nlonp1,:,:)   = ww     (1,:,:)
      zz     (nlonp1,:,:)   = zz     (1,:,:)
      vxb    (nlonp1,:,:)   = vxb    (1,:,:)
      adotv  (nlonp1,:,:,:) = adotv  (1,:,:,:)
      axv    (nlonp1,:,:,:) = axv    (1,:,:,:)
      a1dta2 (nlonp1,:)     = a1dta2 (1,:)
      sini   (nlonp1,:)     = sini   (1,:)
      bxa    (nlonp1,:,:)   = bxa    (1,:,:)
      adota  (nlonp1,:,:)   = adota  (1,:,:)
!
! Transform these fields to geomagnetic coordinate system, one magnetic 
!   latitude at a time. Output fields are at current magnetic latitude.
!
      maglat_loop: do j=2,nmlat-1                               ! do 7
        if (j==nmlat/2+1) cycle maglat_loop
!
! ssigma1(nlonp1,0:nlatp1,-2:nlev), sigma_pedm(nmlonp1,-2:nlevp1)
        call geo2mag(ssigma1,nlonp1,nlatp2, sigma_pedm,
     |    nmlonp1,nmlat, ig,jg,wt, 1,nmlon,-2,nlev,j) 
!
! ssigma2(nlonp1,0:nlatp1,-2:nlev), sigma_hallm(nmlonp1,-2:nlev)
        call geo2mag(ssigma2,nlonp1,nlatp2, sigma_hallm,nmlonp1,nmlat,
     |    ig,jg,wt, 1,nmlon,-2,nlev,j)
!
! ww(nlonp1,0:nlatp1,-2:nlev), wnvelm(nmlonp1,-2:nlev)
        call geo2mag(ww,nlonp1,nlatp2, wnvelm,nmlonp1,nmlat,
     |    ig,jg,wt, 1,nmlon,-2,nlev,j)
!
! vxb(nlonp1,0:nlatp1,-2:nlev), vxbm(nmlonp1,-2:nlev)
        call geo2mag(vxb,nlonp1,nlatp2, vxbm,nmlonp1,nmlat,
     |    ig,jg,wt, 1,nmlon,-2,nlev,j)
!
        do n=1,2
!
! adotv(nlonp1,0:nlatp1,-2:nlev,2), adotvm(nmlonp1,-2:nlev,2)
          call geo2mag(adotv(1,0,-2,n),nlonp1,nlatp2, adotvm(1,-2,n),
     |      nmlonp1,nmlat, ig,jg,wt, 1,nmlon,-2,nlev,j)
!
! axv(nlonp1,0:nlatp1,-2:nlev,2), axvm(nmlonp1,-2:nlev,2)
          call geo2mag(axv(1,0,-2,n),nlonp1,nlatp2, axvm(1,-2,n),
     |      nmlonp1,nmlat, ig,jg,wt, 1,nmlon,-2,nlev,j)
!
! bxa(nlonp1,0:nlatp1,2), bxam(nmlonp1,2)
          call geo2mag(bxa(1,0,n),nlonp1,nlatp2, bxam(1,n), ! no k
     |      nmlonp1,nmlat, ig,jg,wt, 1,nmlon,1,1,j)
!
! adota(nlonp1,0:nlatp1,2), adotam(nmlonp1,2)
          call geo2mag(adota(1,0,n),nlonp1,nlatp2, adotam(1,n), ! no k
     |      nmlonp1,nmlat, ig,jg,wt, 1,nmlon,1,1,j)
        enddo ! n=1,2
!
! zz(nlonp1,0:nlatp1,-2:nlevp1), zpotenm(nmlonp1,-2:nlevp1)
        call geo2mag(zz,nlonp1,nlatp2, zpotenm,nmlonp1,nmlat,
     |    ig,jg,wt, 1,nmlon,-2,nlevp1,j) ! to nlevp1
!
! Save geopotential in 3d array (both are module data above):
        do k=-2,nlevp1
          do i=1,nmlon
            zpotenm3d(i,j,k) = zpotenm(i,k)
          enddo ! i=1,nmlon
        enddo ! k=-2,nlevp1
!
! a1dta2(nlonp1,0:nlatp1), a1a2m(nmlonp1) (no k)
        call geo2mag(a1dta2,nlonp1,nlatp2, a1a2m,nmlonp1,nmlat,
     |    ig,jg,wt, 1,nmlon,1,1,j) ! no k
!
! sini(nlonp1,0:nlatp1), siniam(nmlonp1) (no k)
        call geo2mag(sini,nlonp1,nlatp2, siniam,nmlonp1,nmlat,
     |    ig,jg,wt, 1,nmlon,1,1,j)
!
! bmod(nlonp1,0:nlatp1), bmodm(nmlonp1) (bmod is from magfield)
        call geo2mag(bmod,nlonp1,nlatp2, bmodm,nmlonp1,nmlat,
     |    ig,jg,wt, 1,nmlon,1,1,j) ! no k
!
! p(nlonp1,0:nlatp1), pm(nmlonp1) (p is from magfield)
        call geo2mag(p,nlonp1,nlatp2, pm,nmlonp1,nmlat,
     |    ig,jg,wt, 1,nmlon,1,1,j) ! no k

        do n=1,2
!
! av(nlonp1,0:nlatp1,3,2), azm(nmlonp1,2)
          call geo2mag(av(1,0,3,n),nlonp1,nlatp2, azm(1,n), ! no k
     |      nmlonp1,nmlat, ig,jg,wt, 1,nmlon,1,1,j)
        enddo ! n=1,2
!
! Save fields on magnetic grid to secondary history:
!
! Periodic point:
        do k=-2,nlev
          sigma_pedm (nmlonp1,k)   = sigma_pedm (1,k)
          sigma_hallm(nmlonp1,k)   = sigma_hallm(1,k)
          wnvelm     (nmlonp1,k)   = wnvelm     (1,k)
          vxbm       (nmlonp1,k)   = vxbm       (1,k)
          zpotenm    (nmlonp1,k)   = zpotenm    (1,k)
          adotvm     (nmlonp1,k,1) = adotvm     (1,k,1)
          adotvm     (nmlonp1,k,2) = adotvm     (1,k,2)
          axvm       (nmlonp1,k,1) = axvm       (1,k,1)
          axvm       (nmlonp1,k,2) = axvm       (1,k,2)
        enddo ! k=-2,nlev
!
! Magnetic field line integrals:
! Sub fieldline_integrals calculates PDE coefficients
!   zigm11, zigmc, zigm2, zigm22, rim, rhs
!
        call fieldline_integrals(j)
!
! End magnetic latitude loop:
      enddo maglat_loop ! j=2,nmlat-1
!
! Equatorial values of Z (3d array) (nmlath=(nmlat+1)/2):
      do k=-2,nlevp1
        do i=1,nmlon
          zpotenm3d(i,nmlath,k) = .5*(zpotenm3d(i,nmlath-1,k)+
     |                                zpotenm3d(i,nmlath+1,k))
        enddo ! i=1,nmlon
!
! Polar values for Z:
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
! Save 3d geopotential to secondary history:
!     do j=1,nmlat
!       call addfsech_ik('ZM3D'  ,' ',' ',zpotenm3d(:,j,:),
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!     enddo
!
! Fold southern hemisphere over on to northern:
! (mpi message passing will be necessary here when dynamo is multitasked) 
! If nmlat=97, then nmlath=49. This assigns 96,50 <= 96,50 + 2,48
!
      do j=2,nmlath-1                                          ! do 36
        do i=1,nmlon
          zigm11(i,nmlat+1-j) = zigm11(i,nmlat+1-j)+zigm11(i,j)
          zigmc (i,nmlat+1-j) = zigmc (i,nmlat+1-j)+zigmc (i,j)
          zigm2 (i,nmlat+1-j) = zigm2 (i,nmlat+1-j)+zigm2 (i,j)
          zigm22(i,nmlat+1-j) = zigm22(i,nmlat+1-j)+zigm22(i,j)
          rim(i,nmlat+1-j,1) = rim(i,nmlat+1-j,1)+rim(i,j,1)
          rim(i,nmlat+1-j,2) = rim(i,nmlat+1-j,2)+rim(i,j,2)
        enddo ! i=1,nmlon
      enddo ! j=2,nmlath
!
! Pole:
      zigm11(1,nmlat) = (4.*sddot(nmlon,unitvm,zigm11(1,nmlat-1))-
     |  sddot(nmlon,unitvm,zigm11(1,nmlat-2)))/(3.*float(nmlon))
      zigmc (1,nmlat) = (4.*sddot(nmlon,unitvm,zigmc (1,nmlat-1))-
     |  sddot(nmlon,unitvm,zigmc (1,nmlat-2)))/(3.*float(nmlon))
      zigm2 (1,nmlat) = (4.*sddot(nmlon,unitvm,zigm2 (1,nmlat-1))-
     |  sddot(nmlon,unitvm,zigm2 (1,nmlat-2)))/(3.*float(nmlon))
      zigm22(1,nmlat) = (4.*sddot(nmlon,unitvm,zigm22(1,nmlat-1))-
     |  sddot(nmlon,unitvm,zigm22(1,nmlat-2)))/(3.*float(nmlon))
!
      do i=1,nmlon                    !  i=1,80               ! do 47
        ii = 1+mod(i-1+nmlon/2,nmlon) ! ii=41:80,1:40
        rim(i,nmlat,1) = .5*(rim(i,nmlat-1,1)-rim(ii,nmlat-1,1))
        rim(i,nmlat,2) = .5*(rim(i,nmlat-1,2)-rim(ii,nmlat-1,2))
      enddo ! i=1,nmlon
!
! Extend over longitude:
      do i=2,nmlon
        zigm11(i,nmlat) = zigm11(1,nmlat)
        zigmc (i,nmlat) = zigmc (1,nmlat)
        zigm2 (i,nmlat) = zigm2 (1,nmlat)
        zigm22(i,nmlat) = zigm22(1,nmlat)
!       rim(i,nmlat,1) = rim(1,nmlat,1)
!       rim(i,nmlat,2) = rim(1,nmlat,2)
      enddo ! i=2,nmlon
!
! Equator:
      do i=1,nmlon
        zigm11(i,nmlath) = 0.
        zigmc (i,nmlath) = 0.
        zigm2 (i,nmlath) = 0.
        zigm22(i,nmlath) = 0.
        rim(i,nmlath,1) = 0.
        rim(i,nmlath,2) = 0.
      enddo ! i=1,nmlon
!
! Periodic point:
      do j=1,nmlat
        zigm11(nmlonp1,j) = zigm11(1,j)
        zigmc (nmlonp1,j) = zigmc (1,j)
        zigm2 (nmlonp1,j) = zigm2 (1,j)
        zigm22(nmlonp1,j) = zigm22(1,j)
        rim(nmlonp1,j,1) = rim(1,j,1)
        rim(nmlonp1,j,2) = rim(1,j,2)
      enddo ! j=1,nmlat
!
! Save to secondary histories:
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
!
        do i=1,nmlon                                         ! do 40
          rhs(i,j) = 1./(dlonm*tint1(nmlath+j-1))*
     |      .5*(tint33(i+1)-tint33(i-1))  ! tint33 2:nmlon+1 and 0:nmlon-1
        enddo ! i=1,nmlon
      enddo ! j=2,nmlath-1
!
! Perform differentiation of rim(2) w.r.t. theta0:
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
     |  sddot(nmlon,unitvm,rim(1,nmlat-1,2))/tint1(nmlat-1)
!
! Equator:
      do i=1,nmlon
        rhs(i,1) = 1./dlatm*.5*tint1(nmlath+1)*rim(i,nmlath+1,2)
      enddo ! i=1,nmlon
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
! Scale (multiply by r0*1.e-8 for results in volts):
      do j=1,nmlath
        do i=1,nmlonp1
          rhs(i,j) = rhs(i,j)*r0*1.e-8
        enddo ! i=1,nmlonp1
      enddo ! j=1,nmlath
!
! Save RHS to secondary histories (redundant in vertical and
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
      end subroutine transf
!-----------------------------------------------------------------------
      subroutine fieldline_integrals(latm)
      use cons_module,only: r0,ylatm,rcos0s,dt0dts
!
! Args:
      integer,intent(in) :: latm ! current magnetic latitude index
!
! Local:
      real,parameter :: z0=90.e5
      integer :: i,k
      real,dimension(nmlonp1) :: sindm, cosdm, ram, aam, cosiam,
     |  csthdam, rtadram
      real,dimension(nmlonp1,-2:nlevp1) :: rrm, sinidm, cosidm,
     |  costhdm, rtramrm
      real,dimension(nmlat) :: tint1,tint2,tint3
!
! Diagnostics for plotting:
      real,dimension(nmlonp1,nlevp1) ::
     |  zigm11_plt, zigmc_plt, zigm2_plt, zigm22_plt, rim1_plt, rim2_plt
!
! Evaluate 2d fields:
! (siniam is from geo2mag call above)
      do i=1,nmlon
        ram(i) = r0/cos(ylatm(latm))**2
        cosiam(i) = sqrt(1.-siniam(i)**2)
        aam(i) = 4.*ram(i)*(ram(i)-r0)/((r0*siniam(i)/cosiam(i))**2+
     |           4.*(ram(i)-r0)**2)
        csthdam(i) = sqrt((ram(i)-aam(i)*(ram(i)-r0))/ram(i))
        rtadram(i) = sqrt(aam(i)/ram(i))
      enddo ! i=1,nmlon

!     write(6,"('fieldline_integrals: latm=',i3)") latm
!     write(6,"('latm=',i3,' ram='    ,/,(6e12.4))") latm,ram
!     write(6,"('latm=',i3,' siniam=' ,/,(6e12.4))") latm,siniam
!     write(6,"('latm=',i3,' cosiam=' ,/,(6e12.4))") latm,cosiam
!     write(6,"('latm=',i3,' aam='    ,/,(6e12.4))") latm,aam
!     write(6,"('latm=',i3,' csthdam=',/,(6e12.4))") latm,csthdam
!     write(6,"('latm=',i3,' rtadram=',/,(6e12.4))") latm,rtadram

      do k=-2,nlevp1                                          ! do 8
        do i=1,nmlon
!
! RR = R0+Z-Z0
          rrm(i,k) = r0+zpotenm(i,k)-z0
!
! RTRAMR = RA-R IF +IVE, ZERO OTHERWISE
          rtramrm(i,k) = ram(i)-rrm(i,k)
          if (rtramrm(i,k) < 0.) rtramrm(i,k) = 0.
!
! COSTHD = SQRT((RA-AA*(RA-R))/RA)
          costhdm(i,k) = sqrt((ram(i)-aam(i)*rtramrm(i,k))/ram(i))

!         if (latm >= 90) 
!    |    write(6,"('integrals: i=',i3,' k=',i3,' j=',i3,
!    |      ' ram=',e12.4,' aam=',e12.4,' rtramrm=',e12.4,' costhdm=',
!    |      e12.4)") i,k,latm,ram(i),aam(i),rtramrm(i,k),costhdm(i,k)

!
! RTRAMR = SQRT(RA-R)
          rtramrm(i,k) = sqrt(rtramrm(i,k))
!
! SINID = TAN(ID)
          sinidm(i,k) = sign(2./(rrm(i,k)*rtadram(i))*costhdm(i,k)*
     |      rtramrm(i,k),siniam(i))
!
! COSID = COS(ID)
          cosidm(i,k) = 1./sqrt(1.+sinidm(i,k)**2)
!
! SINID = SIN(ID)
          sinidm(i,k) = sinidm(i,k)*cosidm(i,k)
        enddo ! i=1,nmlon
      enddo ! k=-2,nlevp1
!
! Interpolate to middle of horizontal slices:
      do k=-2,nlev                                            ! do 6
        do i=1,nmlon
          rrm(i,k)     = .5*(rrm(i,k)+rrm(i,k+1))
          sinidm(i,k)  = .5*(sinidm(i,k)+sinidm(i,k+1)) 
          cosidm(i,k)  = .5*(cosidm(i,k)+cosidm(i,k+1)) 
          costhdm(i,k) = .5*(costhdm(i,k)+costhdm(i,k+1)) 
          rtramrm(i,k) = rtramrm(i,k)-rtramrm(i,k+1)
        enddo ! i=1,nmlon
!
! Periodic point:
!       rrm    (nmlonp1,k) = rrm    (1,k)
!       sinidm (nmlonp1,k) = sinidm (1,k)
!       cosidm (nmlonp1,k) = cosidm (1,k)
!       costhdm(nmlonp1,k) = costhdm(1,k)
!       rtramrm(nmlonp1,k) = rtramrm(1,k)
      enddo ! k=-2,nlev

!     call addfsech_ik('RRM'    ,' ',' ',rrm,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!     call addfsech_ik('SINIDM' ,' ',' ',sinidm,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!     call addfsech_ik('COSIDM' ,' ',' ',cosidm,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!     call addfsech_ik('COSTHDM',' ',' ',costhdm,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!     call addfsech_ik('RTRAMRM',' ',' ',rtramrm,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!
! Init:
      do i=1,nmlon
        tint1(i) = 0.
        tint2(i) = 0.
      enddo
!
! QUADRATURES FOR A1**2/P
      do k=-2,nlev                                            ! do 10
        do i=1,nmlon
!
! TINT1 = QUAD(SIG1*R*D(SQRT(RA-R))/COS(THETAD)**2)
          tint1(i) = tint1(i)+sigma_pedm(i,k)*rrm(i,k)/(costhdm(i,k)**2)
     |      *rtramrm(i,k)
!
! TINT2 = QUAD(SIG1*R*COS(ID)**2*D(SQRT(RA-R))/COS(THETAD)**2)
          tint2(i) = tint2(i)+sigma_pedm(i,k)*rrm(i,k)*(cosidm(i,k)/
     |      costhdm(i,k))**2*rtramrm(i,k)
        enddo ! i=1,nmlon

!
! diffs in sigma_pedm
!       write(6,"('integrals: k=',i3,' j=',i3,' sigma_pedm=',/,
!    |    (6e12.5))") k,latm,sigma_pedm(:,k)
!       write(6,"('integrals: k=',i3,' j=',i3,' rrm=',/,
!    |    (6e12.5))") k,latm,rrm(:,k)
!       write(6,"('integrals: k=',i3,' j=',i3,' costhdm=',/,
!    |    (6e12.5))") k,latm,costhdm(:,k)
!       write(6,"('integrals: k=',i3,' j=',i3,' cosidm=',/,
!    |    (6e12.5))") k,latm,cosidm(:,k)
!       write(6,"('integrals: k=',i3,' j=',i3,' rtramrm=',/,
!    |    (6e12.5))") k,latm,rtramrm(:,k)

      enddo ! k=-2,nlev
!
! Complete calculation and store in zigm11:
      do i=1,nmlon
        zigm11(i,latm) = bmodm(i)*csthdam(i)/cosiam(i)*rtadram(i)*
     |    (adotam(i,1)*tint1(i)+(azm(i,1)/cosiam(i))**2*tint2(i))/pm(i)
      enddo ! i=1,nmlon
!
! Init:
      do i=1,nmlon
        tint1(i) = 0.
        tint2(i) = 0.
      enddo
!
! QUADRATURES FOR A2**2/P                                      ! do 13
      do k=-2,nlev
        do i=1,nmlon
!
! TINT1 = QUAD(SIG1*R*TAN(ID)**2*D(SQRT(RA-R)))
          tint1(i) = tint1(i)+sigma_pedm(i,k)*rrm(i,k)*(sinidm(i,k)/
     |      cosidm(i,k))**2*rtramrm(i,k)
!
! TINT2 = QUAD(SIG1*R*D(SQRT(RA-R)))
          tint2(i) = tint2(i)+sigma_pedm(i,k)*rrm(i,k)*rtramrm(i,k)
        enddo ! i=1,nmlon
      enddo ! k=-2,nlev
!
! Complete calculation and store result in zigm22:
      do i=1,nmlon
        zigm22(i,latm) = bmodm(i)*cosiam(i)/csthdam(i)*rtadram(i)*
     |    (adotam(i,2)/siniam(i)**2*tint1(i)+(azm(i,2)/cosiam(i))**2*
     |    tint2(i))/pm(i)
      enddo ! i=1,nmlon
!
! Init:
      do i=1,nmlon
        tint1(i) = 0.
        tint2(i) = 0.
      enddo
!
! QUADRATURES FOR A1.A2/P
      do k=-2,nlev                                              ! do 16
        do i=1,nmlon
!
! TINT1 = QUAD(SIG1*R*TAN(ID)*D(SQRT(RA-R))/COS(THETAD))
          tint1(i) = tint1(i)+sigma_pedm(i,k)*rrm(i,k)*sinidm(i,k)/
     |      (cosidm(i,k)*costhdm(i,k))*rtramrm(i,k)
!
! TINT2 = QUAD(SIG1*R*COS(ID)*D(SQRT(RA-R))/COS(THETAD))
          tint2(i) = tint2(i)+sigma_pedm(i,k)*rrm(i,k)*cosidm(i,k)/
     |      costhdm(i,k)*rtramrm(i,k)
        enddo ! i=1,nmlon
      enddo ! k=-2,nlev
!
! Complete calculation and store in zigmc:
      do i=1,nmlon
        zigmc(i,latm) = bmodm(i)*rtadram(i)*(a1a2m(i)/siniam(i)*
     |    tint1(i)+azm(i,1)*azm(i,2)/(cosiam(i)**2)*tint2(i))/pm(i)
      enddo ! i=1,nmlon
!
! Init:
      do i=1,nmlon
        tint1(i) = 0.
      enddo
!
! QUADRATURE FOR B.(A2 X A1) = 1./B
      do k=-2,nlev                                              ! do 19
        do i=1,nmlon
          tint1(i) = tint1(i)+sigma_hallm(i,k)*rrm(i,k)/(cosidm(i,k)*
     |      costhdm(i,k))*rtramrm(i,k)
        enddo ! i=1,nmlon
      enddo ! k=-2,nlev
!
! Complete calculation and store result in zigm2:
      do i=1,nmlon
        zigm2(i,latm) = rtadram(i)*tint1(i)
      enddo ! i=1,nmlon
!
! Right-hand side terms:
!
! Init:
      do i=1,nmlon
        tint1(i) = 0.
        tint2(i) = 0.
        tint3(i) = 0.
      enddo
!
! A1.(SIGMA1*(V X B)+SIGMA2*V)
! A2.(SIGMA1*(V X B)+SIGMA2*V)
      do k=-2,nlev                                              ! do 22
        do i=1,nmlon
!
!  TINT1 = INTEGRAL(A1(X)*Q(X))
          tint1(i) = tint1(i)+(adotvm(i,k,1)*sigma_hallm(i,k)+
     |      azm(i,1)*wnvelm(i,k)*sigma_hallm(i,k)*cosidm(i,k)/cosiam(i)-
     |      axvm(i,k,1)*sinidm(i,k)*sigma_pedm(i,k)+ 
     |      bxam(i,1)*wnvelm(i,k)*cosidm(i,k)*sigma_pedm(i,k)+
     |     vxbm(i,k)*azm(i,1)*sigma_pedm(i,k)*cosidm(i,k)**2/cosiam(i))*
     |      rrm(i,k)**2/costhdm(i,k)*rtramrm(i,k)
!
          tint2(i)=tint2(i)+(adotvm(i,k,2)*sigma_hallm(i,k)*sinidm(i,k)/
     |        (cosidm(i,k)*siniam(i))+
     |      azm(i,2)*wnvelm(i,k)*sigma_hallm(i,k)/cosiam(i)-
     |      axvm(i,k,2)*sinidm(i,k)**2*sigma_pedm(i,k)/
     |        (siniam(i)*cosidm(i,k))+
     |      bxam(i,2)*wnvelm(i,k)*sigma_pedm(i,k)*sinidm(i,k)/siniam(i)+
     |      vxbm(i,k)*azm(i,2)*sigma_pedm(i,k)*cosidm(i,k)/cosiam(i))* 
     |      rrm(i,k)**2*rtramrm(i,k)
        enddo ! i=1,nmlon
      enddo ! k=-2,nlev
!
! EVALUATE RHS TERMS RIM(1), RIM(2)
      do i=1,nmlon
        rim(i,latm,1) = bmodm(i)**2/r0*rtadram(i)/cosiam(i)*tint1(i)
        rim(i,latm,2) = bmodm(i)**2/r0*rtadram(i)/csthdam(i)*tint2(i)
      enddo ! i=1,nmlon
!
! Scale coefficients and RHS:
      do i=1,nmlon                                              ! do 24
        zigm11(i,latm) = zigm11(i,latm)*rcos0s(latm)/(dt0dts(latm)*
     |    abs(siniam(i)))
        zigmc(i,latm) = zigmc(i,latm)/siniam(i)
        zigm22(i,latm) = zigm22(i,latm)*dt0dts(latm)/(rcos0s(latm)*
     |    abs(siniam(i)))
        rim(i,latm,1) = rim(i,latm,1)*sign(1.,siniam(i))/dt0dts(latm)
        rim(i,latm,2) = rim(i,latm,2)/rcos0s(latm)
      enddo ! i=1,nmlon
!
! Save to secondary histories (will be redundant in pressure):
!     real,dimension(nmlonp1,nlevp1) ::
!    |  zigm11_plt, zigmc_plt, zigm2_plt, zigm22_plt, rim1_plt, rim2_plt
!
!     do i=1,nmlon
!       zigm11_plt(i,:) = zigm11(i,latm)
!       zigmc_plt (i,:) = zigmc (i,latm)
!       zigm2_plt (i,:) = zigm2 (i,latm)
!       zigm22_plt(i,:) = zigm22(i,latm)
!       rim1_plt  (i,:) = rim   (i,latm,1)
!       rim2_plt  (i,:) = rim   (i,latm,2)
!     enddo ! i=1,nmlon
! Periodic point:
!     zigm11_plt(nmlonp1,:) = zigm11_plt(1,:)
!     zigmc_plt (nmlonp1,:) = zigmc_plt (1,:)
!     zigm2_plt (nmlonp1,:) = zigm2_plt (1,:)
!     zigm22_plt(nmlonp1,:) = zigm22_plt(1,:)
!     rim1_plt  (nmlonp1,:) = rim1_plt  (1,:)
!     rim2_plt  (nmlonp1,:) = rim2_plt  (1,:)

!     call addfsech_ik('ZIGM11'  ,' ',' ',zigm11_plt,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!     call addfsech_ik('ZIGMC'   ,' ',' ',zigmc_plt,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!     call addfsech_ik('ZIGM2'   ,' ',' ',zigm2_plt,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!     call addfsech_ik('ZIGM22'  ,' ',' ',zigm22_plt,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!     call addfsech_ik('RIM1'    ,' ',' ',rim1_plt,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)
!     call addfsech_ik('RIM2'    ,' ',' ',rim2_plt,
!    |  1,nmlonp1,nmlev,nmlev-1,latm)

      end subroutine fieldline_integrals
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
      subroutine threed
      use cons_module,only: ylatm,r0,r00,pi_dyn,table,dlatm
      use fields_module,only: phim3d ! phim3d(nmlonp1,nmlat,-2:nlevp1)
!
! Calculate 3-d dynamo potential phim3d in geomagnetic coordinates from
!   2-d potential phim output by the PDE solver mudpack.
! On input, phim(nmlonp1,nmlat) is defined from mud.
!
! Local:
      integer :: k,i,j
      real :: cosltm,pi,phimn,phims
      real,dimension(nmlonp1) :: thetam,pslot,qslot
      integer,dimension(nmlonp1) :: islot,jslot
!
! Externals:
      real,external :: sddot ! in util.F
!
      pi = pi_dyn
      do k=-2,nlevp1
        do j=2,nmlat-1
          cosltm = cos(ylatm(j))
          do i=1,nmlon
!
! Thetam(i) = thetas for each value of i:
            thetam(i) = acos(sqrt(r0/(r00+zpotenm3d(i,j,k)))*cosltm)
!
! Find position of thetam in table of theta0 vs thetas:
            pslot(i) = thetam(i)*180./pi+1.
            islot(i) = pslot(i)
            pslot(i) = pslot(i)-float(islot(i))
!
! Determine theta0 by interpolation in table:
            thetam(i) = (1.-pslot(i))*table(islot(i),2)+pslot(i)*
     |        table(islot(i)+1,2)
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
        enddo ! i = 1,nmlon
      enddo ! k = -2,nlevp1
!
! Periodic point:
      do k=1,nlevp1
        do j=1,nmlat
          phim3d(nmlonp1,j,k) = phim3d(1,j,k)
        enddo ! j=1,nmlat
      enddo ! k=1,nlevp1
!
! Save 3d electric potential on magnetic grid to secondary history:
!     do j=1,nmlat
!       call addfsech_ik('PHIM3D'  ,' ',' ',phim3d(:,j,:),
!    |    1,nmlonp1,nmlev,nmlev-1,j)
!     enddo ! j=1,nmlat
      end subroutine threed
!-----------------------------------------------------------------------
      subroutine mp_dynamo_gather
!
! Input fields to the dynamo have been prepared on the geographic grid
! at task subdomains (see sub prep_dynamo). Now these must be gathered
! to the root task before execution of the dynamo itself. This is called
! from prep_dynamo (dynamo.F).
!
! On input, these fields are defined at subdomains. 
! On output, they are defined globally at the root task.
!
      use mpi_module,only: mytid,lon0,lon1,lat0,lat1,mxlat,mxlon,
     |  ntask,mpi_comm_world,mpi_real8,irstat,tasks,handle_mpi_err
!
! Local:
      integer,parameter :: nfld=7
      integer :: k,i,j,n,ier,len,nlons,nlats,ifld,msgtag,ireqrecv,
     |  ireqsend,ilon0,ilon1,ilat0,ilat1
      integer :: idest = 0
      real :: fmin,fmax
      real,allocatable :: rcvbuf(:,:,:,:),sndbuf(:,:,:,:) ! (k,i,j,ifld)
      character(len=8) :: fnames(nfld) =
     |  (/'SIGMAPED','SIGMAHAL','ZPOTEN  ','SCHEIGHT','UNVEL   ',
     |    'VNVEL   ','WNVEL   '/)
!
! Allocate send and receive buffers:
! (mxlon,mxlat are max number of lons,lats held by all tasks)
      allocate(rcvbuf(mxlon,mxlat,nlevp1,nfld),stat=ier)
      if (ier /= 0)
     |  write(6,"('>>> mp_dynamo_gather: error allocating rcvbuf:',
     |    ' mxlon=',i3,' mxlat=',i3,' nlevp1=',i3,' nfld=',i3,
     |    ' ier=',i4)") mxlon,mxlat,nlevp1,nfld,ier
      allocate(sndbuf(mxlon,mxlat,nlevp1,nfld),stat=ier)
      if (ier /= 0)
     |  write(6,"('>>> mp_dynamo_gather: error allocating sndbuf:',
     |    ' mxlon=',i3,' mxlat=',i3,' nlevp1=',i3,' nfld=',i3,
     |    ' ier=',i4)") mxlon,mxlat,nlevp1,nfld,ier
      len = mxlon*mxlat*nlevp1*nfld ! buffer length
!
! Root receives from all other tasks:
      if (mytid==0) then ! root task
!
! Receive subdomains from slave tasks:
        do n=1,ntask-1   ! receive from slaves only
          msgtag = n
          call mpi_irecv(rcvbuf,len,MPI_REAL8,n,msgtag,MPI_COMM_WORLD,
     |      ireqrecv,ier)
          if (ier /= 0)
     |      call handle_mpi_err(ier,'mp_dynamo_gather irecv from n')
!
! Wait for receives to complete:
          call mpi_wait(ireqrecv,irstat,ier)
          if (ier /= 0) then
            call handle_mpi_err(ier,'mp_dynamo_gather wait for n')
          else
!           write(6,"('mp_dynamo_gather: root recvd from task ',i3)") n
          endif
!
! Transfer from receive buffer to fields:
!     real,dimension(nlonp4,nlat,nlevp1) ::
!
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
     |          rcvbuf(1:nlons,1:nlats,:,ifld)
            case(2)
              sigma_hall(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld)
            case(3)
              zpoten(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld)
            case(4)
              scheight(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld)
            case(5)
              unvel(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld)
            case(6)
              vnvel(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld)
            case(7)
              wnvel(ilon0:ilon1,ilat0:ilat1,:) =
     |          rcvbuf(1:nlons,1:nlats,:,ifld)
            end select
!           call fminmax(rcvbuf(1:nlons,1:nlats,:,ifld),
!    |        nlons*nlats*nlevp1,fmin,fmax)
!           write(6,"('mp_dynamo_gather: rcvbuf for field ',
!    |        a,': min,max=',2e12.4)") fnames(ifld),fmin,fmax
          enddo ! ifld=1,nfld
        enddo ! n=1,ntask-1 (receive from slaves only)
!
! Non-root tasks send to master:
      else  ! mytid /= 0
!
! Load send buffer:
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
!         call fminmax(sndbuf(1:nlons,1:nlats,:,ifld),
!    |      nlons*nlats*nlevp1,fmin,fmax)
!         write(6,"('mp_dynamo_gather: sndbuf for field ',
!    |      a,': min,max=',2e12.4)") fnames(ifld),fmin,fmax
        enddo
!
! Send to root:
        msgtag = mytid
        call mpi_isend(sndbuf,len,MPI_REAL8,idest,msgtag,MPI_COMM_WORLD,
     |    ireqsend,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_dynamo_gather isend to root')
!
! Wait for send to complete:
        call mpi_wait(ireqsend,irstat,ier)
        if (ier /= 0) then
          call handle_mpi_err(ier,
     |    'mp_dynamo_gather wait for send to root')
        endif
      endif ! mytid==0 (root or slave)
!
! Print global min,max from root:
!     if (mytid==0) then
!       do ifld=1,nfld
!         select case (ifld)
!         case(1)
!           call fminmax(sigma_ped,nlonp4*nlat*nlevp1,fmin,fmax)
!         case(2)
!           call fminmax(sigma_hall,nlonp4*nlat*nlevp1,fmin,fmax)
!         case(3)
!           call fminmax(zpoten,nlonp4*nlat*nlevp1,fmin,fmax)
!         case(4)
!           call fminmax(scheight,nlonp4*nlat*nlevp1,fmin,fmax)
!         case(5)
!           call fminmax(unvel,nlonp4*nlat*nlevp1,fmin,fmax)
!         case(6)
!           call fminmax(vnvel,nlonp4*nlat*nlevp1,fmin,fmax)
!         case(7)
!           call fminmax(wnvel,nlonp4*nlat*nlevp1,fmin,fmax)
!         end select
!         write(6,"('mp_dynamo_gather: field ',a,' global min,max=',
!    |      2e12.4)") fnames(ifld),fmin,fmax
!       enddo ! ifld=1,nfld
!     endif
!
! Release local buffer space:
      deallocate(rcvbuf)
      deallocate(sndbuf)
      end subroutine mp_dynamo_gather
!-----------------------------------------------------------------------
      subroutine print_cee
      integer :: n,j
!
      do n=1,10
        do j=1,nmlat0
          write(6,"('print_cee: j=',i3,' nmlat0=',i3,' n=',i2,
     |      ' c0(:,j,n)=',/,(6e12.4))") j,nmlat0,n,c0(:,j,n)
        enddo ! j=1,nmlat0
      enddo ! n=1,10
      do n=1,9
        do j=1,nmlat1
          write(6,"('print_cee: j=',i3,' nmlat1=',i3,' n=',i2,
     |      ' c1(:,j,n)=',/,(6e12.4))") j,nmlat1,n,c1(:,j,n)
        enddo ! j=1,nmlat1
      enddo ! n=1,9
      do n=1,9
        do j=1,nmlat2
          write(6,"('print_cee: j=',i3,' nmlat2=',i3,' n=',i2,
     |      ' c2(:,j,n)=',/,(6e12.4))") j,nmlat2,n,c2(:,j,n)
        enddo ! j=1,nmlat2
      enddo ! n=1,9
      do n=1,9
        do j=1,nmlat3
          write(6,"('print_cee: j=',i3,' nmlat3=',i3,' n=',i2,
     |      ' c3(:,j,n)=',/,(6e12.4))") j,nmlat3,n,c3(:,j,n)
        enddo ! j=1,nmlat3
      enddo ! n=1,9
      do n=1,9
        do j=1,nmlat4
          write(6,"('print_cee: j=',i3,' nmlat4=',i3,' n=',i2,
     |      ' c4(:,j,n)=',/,(6e12.4))") j,nmlat4,n,c4(:,j,n)
        enddo ! j=1,nmlat4
      enddo ! n=1,9

      end subroutine print_cee
      end module dynamo_old_module
