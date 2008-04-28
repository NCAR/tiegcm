!
      module nodynamo_module
      use dynamo_module
      implicit none
!
! These routines are called if namelist read parameter DYNAMO is 0.
!
      contains
!-----------------------------------------------------------------------
      subroutine nodynamo
!
! 12/02:  This routine solves for the 3-D electric potential and fields
!       when the dynamo flag is off (idynamo=0).  Hence, the low-latitude
!       potential is zero, while the high-latitude potential is set
!       with potential_model='HEELIS','WEIMER', or 'NONE'.
!
      use input_module,only: potential_model
      use cons_module,only: dlatm,dlonm,pi_dyn,ylatm,rtd
      use magfield_module,only: im,jm,dim,djm
      use fields_module,only: dynpot,phim3d
!
! Local:
      integer :: i,j,k
      real,dimension(nlonp1,0:nlatp1) :: phih   ! potential in geographic
!
!  Set phim=phihm of the chosen high-latitude model
      do i=1,nmlonp1
        do j=1,nmlat
          phim(i,j) = phihm(i,j)
        enddo  ! j=1,nmlat
      enddo  ! i=1,nmlonp1
!
!  transform geopotential height to mag. grid (used in subroutine threed)
      call transf_nodynamo
!
! Call threed to generate 3-d potential array in geomagnetic coordinates
!   from 2-d solver output phim, corrected for the SH potential.
!   phim3d(nmlonp1,nmlat,-2:nlevp1) is in fields.F.
!
      call threed
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
! dynpot(nlonp1,0:nlatp1,nlevp1),  ! 3d electric potential geographic
! dynpot_diag(nlonp4,nlat,nlevp1)  ! for addfld
!
!     do j=1,nlat
!       dynpot_diag(:,:,:) = spval
!       dynpot_diag(3:nlonp2,j,:) = dynpot(1:nlon,j,:)               ! 3,74 <= 1,72
!       dynpot_diag(1:2,j,:)      = dynpot_diag(nlonp1:nlonp2,j,:)   ! 1,2 <= 73,74
!       dynpot_diag(nlonp4-1:nlonp4,j,:) = dynpot_diag(3:4,j,:)      ! 75,76 <= 3,4
!       call addfld('DYNPOT','elect. potential (geo)','V',
!    |    dynpot_diag(:,j,:),'lon',1,nlonp4,'ilev',1,nlevp1,j)
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
!
! Save 2d heelis potential on geographic grid:
!      if (potential_model == 'WEIMER01') then
!       call addfld('PHIHM2D','2D WEIMER01 ELECTRIC POTENTIAL',
!    |    'V', phih,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!      elseif (potential_model == 'WEIMER05') then
!       call addfld('PHIHM2D','2D WEIMER05 ELECTRIC POTENTIAL',
!    |    'V', phih,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!      elseif (potential_model == 'HEELIS') then
!       call addfld('PHIHM2D','2D HEELIS ELECTRIC POTENTIAL',
!    |    'V',phih,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!      else
!       call addfld('PHIHM2D','2D ZERO ELECTRIC POTENTIAL',
!    |    'V',phih,'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!      endif
      end subroutine nodynamo
!-----------------------------------------------------------------------
      subroutine transf_nodynamo
!     
! am_11/03  used only if dynamo = 0
! need zpoten3d in magnetic coordinates for subroutine threed
!     
      use cons_module,only: h0,r0
      use magfield_module,only: 
     |  zb,       ! (nlonp1,0:nlatp1)
     |  ig,       ! (nmlonp1,nmlat) geog lon at each geomag grid point
     |  jg,       ! (nmlonp1,nmlat) geog lat at each geomag grid point
     |  wt        ! (4,nmlonp1,nmlat) interpolation weights
!
! Local:
      integer :: i,ii,k,kk,j,jj,lat,n
      real :: z0
!
! Fields to be transformed to geomagnetic space (formerly in transmag.h):
! (these will be input to the geographic to magnetic transformation). 
      real,dimension(nlonp1,0:nlatp1,-2:nlevp1) :: zz
!
      real,dimension(nlonp4,nlevp1) :: zz_plt
!
! External:
      real,external :: sddot ! in util.F
!
! Set constants:
!   z0 is lowest level for start of field line integration set in h0 
!     (h0 in cons module)
!
!
! Pack inputs from 3->nlon+3 to 1->nlon+1, as in tgcm15 (end of sub lamdas):
      do lat=1,nlat
        do i=1,nlon+1
          ii = i+2
          zpoten(i,lat,:) = zpoten(ii,lat,:)
        enddo ! i=1,nlonp4-3
      enddo ! lat=1,nlat
!
      z0 = h0
!
! Calculate quantities to be transformed to geomagnetic space:
      do k=1,nlev        
        do j=1,nlat
          do i=1,nlonp1
            zz(i,j,k) = zpoten(i,j,k)
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
          enddo ! i=1,nlonp1
        enddo ! j=1,nlat
      enddo ! k=0,-2,-1
!
! Values at poles:
      do k=-2,nlev
        zz(1,0,k) = 
     |    (9.*sddot(nlon,unitv,zz(1,1,k))-
     |        sddot(nlon,unitv,zz(1,2,k)))/(8.*float(nlon))
        zz(1,nlatp1,k) = 
     |    (9.*sddot(nlon,unitv,zz(1,nlat  ,k))-
     |        sddot(nlon,unitv,zz(1,nlat-1,k)))/(8.*float(nlon))
!
! Extend in longitude:
        do i = 2,nlon
          zz(i,0,k)             = zz(1,0,k)
          zz(i,nlatp1,k)        = zz(1,nlatp1,k)
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
!
! Extend in longitude:
      do i = 2,nlon
        zz(i,0,nlevp1)     = zz(1,0,nlevp1)
        zz(i,nlatp1,nlevp1)= zz(1,nlatp1,nlevp1)
      enddo ! i = 2,nlon
!
! Periodic points:
      do j = 0,nlatp1
        do k = -2,nlev
           zz(nlonp1,j,k)      = zz(1,j,k)
        enddo ! k = -2,nlev
        zz(nlonp1,j,nlevp1)= zz(1,j,nlevp1)
      enddo ! j = 0,nlatp1
!
! Transform needed fields to geomagnetic coordinate system, one latitude
! at a time. 
!
! Magnetic latitude loop:
      maglat_loop: do j=2,nmlat-1
        do k=-2,nlevp1
          call geo2mag(zpotenm(1,k),zz(1,0,k),ig,jg,wt,nlonp1,
     |      nmlonp1,nmlon,nmlat,j)
        enddo ! k=-2,nlevp1
! Without calculation of k_(q,lam)
        do i=1,nmlon  
          do k=-2,nlev
            if (zpotenm(i,k) < z0) zpotenm(i,k) = z0
            zpotenm3d(i,j,k) = zpotenm(i,k)
          enddo ! k=-2,nlev
          if (zpotenm(i,nlevp1) < z0) zpotenm(i,nlevp1) = z0
          zpotenm3d(i,j,nlevp1) = zpotenm(i,nlevp1)
        enddo ! i=1,nmlon 
      enddo maglat_loop ! j=2,nmlat-1 Main magnetic latitude loop
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
! Save 3d potential on magnetic grid to secondary history:
! PLEASE DO NOT COMMENT THIS OUT -- ZMAG is a mandatory mag field for sech
! real,dimension(nmlonp1,nmlat,-2:nlevp1) :: zpotenm3d
!
      do j=1,nmlat
        call addfld('ZMAG','GEOPOTENTIAL HEIGHT (MAGNETIC)','CM',
     |    zpotenm3d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nimlevp1,j)
      enddo ! j=1,nmlat
!
      end subroutine transf_nodynamo
!-----------------------------------------------------------------------
      subroutine prep_nodynamo(z,lev0,lev1,lon0,lon1,lat0,lat1)
      use cons_module,only: gask,grav
!      
! am_11/03: if dynamo = 0 need zpoten3d for subroutine threed
! Prepare geographic-grid fields for input to the dynamo, and gather them
! to the root task. This is executed by all tasks, and is called from 
! advance before the dynamo itself (which is executed by master task only).
!
      integer :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Input fields at task subdomains. These are from current (not updated)
! fields (itp):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  z     ! geopotential height
!     
! Local:
      integer :: k,i,lat,nk,lonbeg,lonend
      real :: fmin,fmax
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
            zpoten    (i,lat,k) = z   (k,i,lat)
          enddo ! k=lev0,lev1
          zpoten(i,lat,lev1) = z(lev1,i,lat)
        enddo ! i=lon0,lon1
!
      enddo ! lat=lat0,lat1

#ifdef MPI
!
! Gather dynamo input fields to the root task, defining module data 
! above at the global domain on the root task:
!
      call mp_dynamo_gather
#endif
      end subroutine prep_nodynamo
!-----------------------------------------------------------------------
      end module nodynamo_module

