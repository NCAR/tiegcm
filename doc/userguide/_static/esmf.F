      module my_esmf
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! This module contains calls to ESMF (Earth System Modeling Framework)
! to perform mag2geo and geo2mag regridding in an MPI environment for 
! the parallel dynamo code (pdynamo.F).
!
      use esmf
      use params_module,only: gmlon,gmlat,glat,glon,nlev,nlevp1,nlon,
     |  nlat,nmlon,nmlat,mlev0,mlev1
      use mpi_module,only: tasks,lon0,lon1,lat0,lat1,ntask,
     |  ntaski,ntaskj,nmagtaski,nmagtaskj,mlon0,mlon1,mlat0,mlat1,mytid
      use getapex_module,only: gdlatdeg,gdlondeg ! (nmlonp1,nmlat)
      implicit none
!
! Geographic and magnetic grids (as source and destination grids):
      type(ESMF_Grid),save :: 
     |  geo_src_grid, mag_src_grid, ! source grids (will not have periodic pts)
     |  geo_des_grid, mag_des_grid  ! destination grids (will have periodic pts)
!
! 3d (i,j,k) ESMF Fields on geographic subdomains:
      type(ESMF_Field),save :: geo_ped, geo_hal, geo_zpot, 
     |  geo_adotv1, geo_adotv2, geo_je1pg, geo_je2pg
      type(ESMF_Field),save :: esmf_3dfields(8) ! used in pdynamo.F
!
! 3d (i,j,k) ESMF fields regridded to magnetic subdomains:
      type(ESMF_Field),save :: mag_ped, mag_hal, mag_zpot, 
     |  mag_adotv1, mag_adotv2, mag_je1pg, mag_je2pg
!
! 2d (i,j) ESMF fields on geographic subdomains:
      type(ESMF_Field),save :: geo_sini, geo_adota1, geo_adota2,
     |  geo_a1dta2, geo_be3
!
! 2d (i,j) ESMF fields on magnetic subdomains:
      type(ESMF_Field),save :: mag_sini, mag_adota1, mag_adota2,
     |  mag_a1dta2, mag_be3
#if defined(INTERCOMM) || defined(CISMAH)
      type(ESMF_Field),save :: mag_zigm1, mag_zigm2, mag_nsrhs
      type(ESMF_Field),save :: geo_zigm1, geo_zigm2, geo_nsrhs
#endif
!
! 3d electric potential for mag to geo regridding:
      type(ESMF_Field),save :: mag_phi3d,mag_ephi3d,mag_elam3d,mag_emz3d
      type(ESMF_Field),save :: geo_phi3d,geo_ephi3d,geo_elam3d,geo_emz3d

      type(ESMF_RouteHandle),save :: ! route handles for Regridding
     |  routehandle_geo2mag, routehandle_mag2geo,
     |  routehandle_geo2mag_2d, routehandle_mag2geo_2d
!
      contains
!-----------------------------------------------------------------------
      subroutine esmf_init
      use cons_module,only: rtd
!
! Called once per run at model initialization.
! Create ESMF grids for geographic and magnetic, and create ESMF fields
!   as necessary on both grids. Define the 2d coordinates for each grid,
!   and save an ESMF routehandles for geo2mag and mag2geo regridding.
!
      integer :: rc,istat
      integer(ESMF_KIND_I4),pointer :: factorIndexList(:,:) ! esmf (2,nw)
      real(ESMF_KIND_R8),pointer :: factorList(:)           ! esmf weights (nw)
!
! Global geo indices for geo2mag with all 4 corners. nlat+2 includes poles
! (Used for diagnostics only)
      integer(ESMF_KIND_I4) :: geoIndex_g2m(nlon,nlat+2) 
      integer(ESMF_KIND_I4),pointer :: geoIndexSeq_g2m(:)  ! reshaped geoIndex
      integer(ESMF_KIND_I4),pointer :: magIndex_g2m(:,:)   ! nmlatlons,nmaglats
      integer(ESMF_KIND_I4),pointer :: magIndexSeq_g2m(:)  ! reshaped magIndex
!
! Global mag indices for mag2geo for diagnostic/debug:
      integer(ESMF_KIND_I4) :: magIndex_m2g(nmlon,nmlat) 
      integer(ESMF_KIND_I4),pointer :: magIndexSeq_m2g(:)  ! reshaped magIndex
      integer(ESMF_KIND_I4),pointer :: geoIndex_m2g(:,:)   ! nlonp4,nlat+2
      integer(ESMF_KIND_I4),pointer :: geoIndexSeq_m2g(:)  ! reshaped geoIndex

      real(ESMF_KIND_R8),pointer :: fptr(:,:,:)  ! for info only
      integer :: lbnd_destgeo(3),ubnd_destgeo(3) ! 3d bounds of destination geo grid
      integer :: lbnd_destmag(3),ubnd_destmag(3) ! 3d bounds of destination mag grid
      integer :: lbnd_srcgeo(3),ubnd_srcgeo(3)   ! 3d bounds of source geo grid
      integer :: lbnd_srcmag(3),ubnd_srcmag(3)   ! 3d bounds of source mag grid
!
! Open esmf:
      call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI,rc=rc)
!
! Make magnetic and geographic grids for geo2mag regridding:
      call create_geo_grid(geo_src_grid,'src')  ! geo source grid
      call create_mag_grid(mag_des_grid,'des')  ! mag destination grid
!
! Make grids for mag2geo regridding:
      call create_mag_grid(mag_src_grid,'src')  ! mag source grid
      call create_geo_grid(geo_des_grid,'des')  ! geo destination grid
!
! Create empty fields on geographic grid that will be transformed to
!   the magnetic grid and passed as input to the dynamo. This does not
!   assign any values.
!
! 3d fields on source geo grid (these exclude periodic points):
      call esmf_create_geofield(geo_ped,geo_src_grid,'PED     ' ,
     |  mlev0,mlev1)
      call esmf_create_geofield(geo_hal ,geo_src_grid,'HAL     ',
     |  mlev0,mlev1)
      call esmf_create_geofield(geo_zpot,geo_src_grid,'ZPOT    ',
     |  mlev0,mlev1)
      call esmf_create_geofield(geo_adotv1,geo_src_grid,'ADOTV1  ',
     |  mlev0,mlev1)
      call esmf_create_geofield(geo_adotv2,geo_src_grid,'ADOTV2  ',
     |  mlev0,mlev1)
      call esmf_create_geofield(geo_je1pg,geo_src_grid,'JE1PG   ',
     |  mlev0,mlev1)
      call esmf_create_geofield(geo_je2pg,geo_src_grid,'JE2PG   ',
     |  mlev0,mlev1)
!
! 2d fields on source geo grid (these exclude periodic points:
      call esmf_create_geofield(geo_sini  ,geo_src_grid,'SINI    ',0,0)
      call esmf_create_geofield(geo_adota1,geo_src_grid,'ADOTA1  ',0,0)
      call esmf_create_geofield(geo_adota2,geo_src_grid,'ADOTA2  ',0,0)
      call esmf_create_geofield(geo_a1dta2,geo_src_grid,'A1DTA2  ',0,0)
      call esmf_create_geofield(geo_be3   ,geo_src_grid,'BE3     ',0,0)
!
! 3d fields on destination mag grid (will include periodic point):
      call esmf_create_magfield(mag_ped ,mag_des_grid,'PED     ',
     |  mlev0,mlev1)
      call esmf_create_magfield(mag_hal ,mag_des_grid,'HAL     ',
     |  mlev0,mlev1)
      call esmf_create_magfield(mag_zpot,mag_des_grid,'ZPOT    ',
     |  mlev0,mlev1)
      call esmf_create_magfield(mag_adotv1,mag_des_grid,'ADOTV1  ' ,
     |  mlev0,mlev1)
      call esmf_create_magfield(mag_adotv2,mag_des_grid,'ADOTV2  ' ,
     |  mlev0,mlev1)
      call esmf_create_magfield(mag_je1pg,mag_des_grid,'JE1PG   ' ,
     |  mlev0,mlev1)
      call esmf_create_magfield(mag_je2pg,mag_des_grid,'JE2PG   ' ,
     |  mlev0,mlev1)
!
! Get 3d bounds of destination mag field:
      call ESMF_FieldGet(mag_ped,localDe=0,farrayPtr=fptr,
     |  computationalLBound=lbnd_destmag,
     |  computationalUBound=ubnd_destmag,rc=rc)
!     write(6,"('esmf_init: lbnd_destmag=',3i4,' ubnd_destmag=',3i4)")
!    |  lbnd_destmag,ubnd_destmag
!
! 2d fields on destination mag grid (will include periodic point):
      call esmf_create_magfield(mag_sini  ,mag_des_grid,'SINI    ',0,0)
      call esmf_create_magfield(mag_adota1,mag_des_grid,'ADOTA1  ',0,0)
      call esmf_create_magfield(mag_adota2,mag_des_grid,'ADOTA2  ',0,0)
      call esmf_create_magfield(mag_a1dta2,mag_des_grid,'A1DTA2  ',0,0)
      call esmf_create_magfield(mag_be3   ,mag_des_grid,'BE3     ',0,0)
!
! 2d source mag grid for CISM
#if defined(INTERCOMM) || defined(CISMAH)
      call esmf_create_magfield(mag_zigm1 ,mag_src_grid,'ZIGM1 ',0,0)
      call esmf_create_magfield(mag_zigm2 ,mag_src_grid,'ZIGM2 ',0,0)
      call esmf_create_magfield(mag_nsrhs ,mag_src_grid,'NSRHS ',0,0)
#endif
!
! 3d fields on source mag grid for mag2geo:
      call esmf_create_magfield(mag_phi3d ,mag_src_grid,'PHIM3D  ',
     |  mlev0,mlev1)
      call esmf_create_magfield(mag_ephi3d,mag_src_grid,'EPHI3D  ',
     |  mlev0,mlev1)
      call esmf_create_magfield(mag_elam3d,mag_src_grid,'ELAM3D  ',
     |  mlev0,mlev1)
      call esmf_create_magfield(mag_emz3d ,mag_src_grid,'EMZ3D   ',
     |  mlev0,mlev1)
!
! 2d fields on destination geo grid for mag2geo:
#if defined(INTERCOMM) || defined(CISMAH)
      call esmf_create_geofield(geo_zigm1 ,geo_des_grid,'ZIGM1 ',0,0)
      call esmf_create_geofield(geo_zigm2 ,geo_des_grid,'ZIGM2 ',0,0)
      call esmf_create_geofield(geo_nsrhs ,geo_des_grid,'NSRHS ',0,0)
#endif
!
! 3d fields on destination geo grid for mag2geo:
      call esmf_create_geofield(geo_phi3d ,geo_des_grid,'PHIG3D  ',
     |  mlev0,mlev1)
      call esmf_create_geofield(geo_ephi3d,geo_des_grid,'EPHI3D  ',
     |  mlev0,mlev1)
      call esmf_create_geofield(geo_elam3d,geo_des_grid,'ELAM3D  ',
     |  mlev0,mlev1)
      call esmf_create_geofield(geo_emz3d ,geo_des_grid,'EMZ3D   ',
     |  mlev0,mlev1)
!
! Get 3d bounds of source mag field:
      call ESMF_FieldGet(mag_phi3d,localDe=0,farrayPtr=fptr,
     |  computationalLBound=lbnd_srcmag,
     |  computationalUBound=ubnd_srcmag,rc=rc)

!     write(6,"('esmf_init: lon bnd_srcmag =',2i4,' gmlon=',2f9.3)")
!    |  lbnd_srcmag(1),ubnd_srcmag(1),
!    |  gmlon(lbnd_srcmag(1)),gmlon(ubnd_srcmag(1))
!     write(6,"('esmf_init: lat bnd_srcmag =',2i4,' gmlat=',2f9.3)")
!    |  lbnd_srcmag(2),ubnd_srcmag(2),
!    |  gmlat(lbnd_srcmag(2)),gmlat(ubnd_srcmag(2))
!
! Get 3d bounds of destination geo field:
      call ESMF_FieldGet(geo_phi3d,localDe=0,farrayPtr=fptr,
     |  computationalLBound=lbnd_destgeo,
     |  computationalUBound=ubnd_destgeo,rc=rc)

!     write(6,"('esmf_init: lon bnd_destgeo=',2i4,' glon=',2f9.3)")
!    |  lbnd_destgeo(1),ubnd_destgeo(1),
!    |  glon(lbnd_destgeo(1)),glon(ubnd_destgeo(1))
!     write(6,"('esmf_init: lat bnd_destgeo=',2i4,' glat=',2f9.3)")
!    |  lbnd_destgeo(2),ubnd_destgeo(2),
!    |  glat(lbnd_destgeo(2)),glat(ubnd_destgeo(2))
!
! Save route handles for grid transformations in both directions
! geo2mag and mag2geo. FieldRegridStore needs to be called only 
! once for each transformation before the timestep loop (src and
! dest fields are still required, so just use ped here). Once inside 
! the timestep loop, the same routehandle can be used for all fields 
! that are regridded in the given direction.
!
! These calls will leave *.vtk info files in execdir:
!     call ESMF_GridWriteVTK(geo_src_grid, 
!    |  staggerloc=ESMF_STAGGERLOC_CENTER, filename="geoSrcGrid",rc=rc)
!     call ESMF_GridWriteVTK(mag_des_grid, 
!    |  staggerloc=ESMF_STAGGERLOC_CENTER, filename="magDesGrid",rc=rc)
!
! Save route handle and get esmf indices and weights for geo2mag:
! Use polemethod=ESMF_POLEMETHOD_ALLAVG to define geographic poles:
!
      call ESMF_FieldRegridStore(srcField=geo_ped,dstField=mag_ped, 
     |  regridMethod=ESMF_REGRIDMETHOD_BILINEAR,
     |  polemethod=ESMF_POLEMETHOD_ALLAVG,
     |  routeHandle=routehandle_geo2mag,factorIndexList=factorIndexList,
     |  factorList=factorList,rc=rc)

      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> esmf_init: error return from ',
     |    'ESMF_FieldRegridStore for 3d geo2mag: rc=',i4)") rc
        call shutdown('esmf_FieldRegridStore for 3d geo2mag')
      endif
!
! Store route handle for geo2mag 3d fields:
!
! Use esmf index list and weights for 3d fields:
      call ESMF_FieldSMMStore(geo_ped,mag_ped,routehandle_geo2mag,
     |  factorList,factorIndexList,rc=rc) 

      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> esmf_init: error return from ',
     |    'ESMF_FieldSMMStore for 3d geo2mag: rc=',i4)") rc
        call shutdown('ESMF_FieldSMMStore')
      endif
!
! Store route handle geo2mag 2d fields:
!
! Use esmf index list and weights for 2d fields:
      call ESMF_FieldSMMStore(geo_sini,mag_sini,routehandle_geo2mag_2d,
     |  factorList,factorIndexList,rc=rc) 
 
      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> esmf_init: error return from ',
     |    'ESMF_FieldSMMStore for 2d geo2mag: rc=',i4)") rc
        call shutdown('esmf_FieldSMMStore')
      endif
!
! Save route handle and get esmf indices and weights for mag2geo:
! (this overwrites factorIndexList and factorList from geo2mag call above)
!
! These calls will leave *.vtk info files in execdir:
!     call ESMF_GridWriteVTK(mag_src_grid, 
!    |  staggerloc=ESMF_STAGGERLOC_CENTER, filename="magSrcGrid",rc=rc)
!     call ESMF_GridWriteVTK(geo_des_grid, 
!    |  staggerloc=ESMF_STAGGERLOC_CENTER, filename="geoDesGrid",rc=rc)

      call ESMF_FieldRegridStore(srcField=mag_phi3d,dstField=geo_phi3d, 
     |  regridMethod=ESMF_REGRIDMETHOD_BILINEAR,
     |  polemethod=ESMF_POLEMETHOD_ALLAVG,
     |  routeHandle=routehandle_mag2geo,factorIndexList=factorIndexList,
     |  factorList=factorList,rc=rc)

      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> esmf_init: error return from ',
     |    'ESMF_FieldRegridStore for 3d mag2geo: rc=',i4)") rc
        call shutdown('esmf_FieldRegridStore')
      endif
!
! mag2geo 3d fields:
!
! Use esmf index list and weights for 3d fields:
      call ESMF_FieldSMMStore(mag_phi3d,geo_phi3d,routehandle_mag2geo,
     |  factorList,factorIndexList,rc=rc) 

      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> esmf_init: error return from ',
     |    'ESMF_FieldSMMStore for 3d mag2geo: rc=',i4)") rc
        call shutdown('ESMF_FieldSMMStore')
      endif
!
! mag2geo 2d fields:
!
! Use esmf index list and weights for 3d fields:
#if defined(INTERCOMM) || defined(CISMAH)
      call ESMF_FieldRegridStore(srcField=mag_zigm1,dstField=geo_zigm1, 
     |  regridMethod=ESMF_REGRIDMETHOD_BILINEAR,
     |  polemethod=ESMF_POLEMETHOD_ALLAVG,
     |  routeHandle=routehandle_mag2geo_2d,
     |  factorIndexList=factorIndexList,
     |  factorList=factorList,rc=rc)

      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> esmf_init: error return from ',
     |    'ESMF_FieldRegridStore for 2d mag2geo: rc=',i4)") rc
        call shutdown('esmf_FieldRegridStore')
      endif

      call ESMF_FieldSMMStore(mag_zigm1,geo_zigm1,
     |  routehandle_mag2geo_2d, factorList,factorIndexList,rc=rc) 

      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> esmf_init: error return from ',
     |    'ESMF_FieldSMMStore for 2d mag2geo: rc=',i4)") rc
        call shutdown('ESMF_FieldSMMStore')
      endif
#endif

      end subroutine esmf_init
!-----------------------------------------------------------------------
      subroutine create_mag_grid(grid_out,srcdes)
!
! Create ESMF geomagnetic grid, w/ lon,lat coordinates.
! This is called from esmf_init during model initialization.
!
! Args:
      type(ESMF_Grid),intent(out) :: grid_out
      character(len=*),intent(in) :: srcdes
!
! Local:
      integer :: i,j,k,n,rc,istat
      integer :: lbnd(2),ubnd(2)
      real(ESMF_KIND_R8),pointer :: coordX(:,:),coordY(:,:)
      integer,allocatable :: nmlons_task(:) ! dim nmagtaski
      integer,allocatable :: nmlats_task(:) ! dim nmagtaskj
      type(ESMF_VM) :: vm
!
! We are creating either a source grid or a destination grid:
      if (srcdes /= 'src' .and. srcdes /= 'des') then
        write(6,"('>>> create_mag_grid: srcdes = ''',a,
     |    ''' but must be either ''src'' or ''des''')") srcdes
        call shutdown('srcdes')
      endif
!
! nmlons_task(nmagtaski) = number of mag lons per task in lon dim
      if (.not.allocated(nmlons_task)) allocate(nmlons_task(nmagtaski),
     |  stat=istat)
      if (istat /= 0) call shutdown('allocate nmlons_task')
      do i=1,nmagtaski
        loop: do n=0,ntask-1
          if (tasks(n)%magtidi==i-1) then
            nmlons_task(i) = tasks(n)%nmaglons
            exit loop
          endif
        enddo loop
      enddo
!
! Exclude periodic points (1 point fewer for mpi tasks at east end)
! for source grids (this overwrites above for eastern-most tasks):
!
      if (srcdes == 'src') then
        do n=0,ntask-1
          if (tasks(n)%magtidi==nmagtaski-1) then  ! east edge of proc matrix
            nmlons_task(tasks(n)%magtidi+1) = tasks(n)%nmaglons-1
          endif
        enddo
      endif
!
! nmlats_task(nmagtaskj) = number of mag lats per task in lat dim
      if (.not.allocated(nmlats_task)) allocate(nmlats_task(nmagtaskj),
     |  stat=istat)
      if (istat /= 0) call shutdown('allocate nmlats_task')
      do j=1,nmagtaskj
        loop1: do n=0,ntask-1
          if (tasks(n)%magtidj==j-1) then
            nmlats_task(j) = tasks(n)%nmaglats
            exit loop1
          endif
        enddo loop1
      enddo
!
! Create curvilinear magnetic grid (both coords depend
! on both dimensions, i.e., lon(i,j),lat(i,j)):
!
      grid_out = ESMF_GridCreate1PeriDim(
     |  countsPerDEDim1=nmlons_task, coordDep1=(/1,2/),
     |  countsPerDEDim2=nmlats_task, coordDep2=(/1,2/),
     |  indexflag=ESMF_INDEX_GLOBAL,rc=rc)

      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> create_mag_grid: error return from ',
     |    'ESMF_GridCreateShapeTile: rc=',i4)") rc
        call shutdown('create_mag_grid')
      endif
!
! Allocate coordinates:
!
      call ESMF_GridAddCoord(grid_out,
     |  staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)

      if (rc /=ESMF_SUCCESS) then
        write(6,"('>>> create_mag_grid: error return from ',
     |    'ESMF_GridAddCoord: rc=',i4)") rc
        call shutdown('create_mag_grid')
      endif
!
! Get pointer and set mag grid longitude coordinates:
!
      call ESMF_GridGetCoord(grid_out, coordDim=1, localDE=0,
     |  computationalLBound=lbnd, computationalUBound=ubnd,
     |  farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> create_mag_grid: error return from ',
     |    'ESMF_GridGetCoord for longitude coords: rc=',i4)") rc
        call shutdown('create_mag_grid')
      endif

! real,dimension(nmlonp1,nmlat) :: gdlatdeg,gdlondeg

      do j=lbnd(2),ubnd(2)
        do i=lbnd(1),ubnd(1)
          coordX(i,j) = gdlondeg(i,j)
        enddo
      enddo
!
! Get pointer and set mag grid latitude coordinates:
!
      call ESMF_GridGetCoord(grid_out, coordDim=2, localDE=0,
     |  computationalLBound=lbnd, computationalUBound=ubnd,
     |  farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> create_mag_grid: error return from ',
     |    'ESMF_GridGetCoord for latitude coords: rc=',i4)") rc
        call shutdown('create_mag_grid')
      endif

      do j=lbnd(2),ubnd(2)
        do i=lbnd(1),ubnd(1)
          coordY(i,j) = gdlatdeg(i,j)
        enddo
      enddo

!     write(6,"('Created ESMF ',a,' mag grid: lbnd,ubnd_lon=',2i4,
!    |  ' mlon0,1=',2i4,' lbnd,ubnd_lat=',2i4,' mlat0,1=',2i4)")
!    |  srcdes,lbnd(1),ubnd(1),mlon0,mlon1,
!    |         lbnd(2),ubnd(2),mlat0,mlat1

!     write(6,"(a,' mag grid: my nmlons=',i4,' coordX=',(8f9.2))") 
!    |  srcdes,ubnd(1)-lbnd(1)+1,coordX
!     write(6,"(a,' mag grid: my nmlats=',i4,' coordY=',(8f9.2))") 
!    |  srcdes,ubnd(2)-lbnd(2)+1,coordY

      end subroutine create_mag_grid
!-----------------------------------------------------------------------
      subroutine create_geo_grid(grid_out,srcdes)
!
! Create ESMF geographic grid, w/ lon,lat coordinates (excluding
! periodic points). This is called once per run from tgcm.F. 
!
! Args:
      type(ESMF_Grid),intent(out) :: grid_out
      character(len=*),intent(in) :: srcdes
!
! Local:
      integer :: i,j,k,n,rc,istat,lbnd_lat,ubnd_lat,lbnd_lon,ubnd_lon,
     |  lbnd(1),ubnd(1)
      real(ESMF_KIND_R8),pointer :: coordX(:),coordY(:)
      integer,allocatable :: nlons_task(:) ! dim ntaski
      integer,allocatable :: nlats_task(:) ! dim ntaskj
!
! We are creating either a source grid or a destination grid:
      if (srcdes /= 'src' .and. srcdes /= 'des') then
        write(6,"('>>> create_mag_grid: srcdes = ''',a,
     |    ''' but must be either ''src'' or ''des''')") srcdes
        call shutdown('srcdes')
      endif
!
! nlons_task(ntaski) = number of geo lons per task in lon dim
      if (.not.allocated(nlons_task)) allocate(nlons_task(ntaski),
     |  stat=istat)
      if (istat /= 0) call shutdown('allocate nlons_task')
      do i=1,ntaski
        loop: do n=0,ntask-1
          if (tasks(n)%mytidi==i-1) then
            nlons_task(i) = tasks(n)%nlons
            exit loop
          endif
        enddo loop
      enddo
!
! Exclude periodic points (2 points fewer for procs at each end)
! for source grids only:
!
      if (srcdes == 'src') then
        if (ntask==1) then ! for single task, remove 4 periodic points
          nlons_task(1) = tasks(0)%nlons-4
        else
          do n=0,ntask-1
! east or west edge of task table:
            if (tasks(n)%mytidi==ntaski-1.or.tasks(n)%mytidi==0) then
              nlons_task(tasks(n)%mytidi+1) = tasks(n)%nlons-2
            endif
          enddo
        endif ! ntask==1
      endif
!
! nlats_task(ntaskj) = number of geo lats per task in lat dim
      if (.not.allocated(nlats_task)) allocate(nlats_task(ntaskj),
     |  stat=istat)
      if (istat /= 0) call shutdown('allocate nlats_task')
      do j=1,ntaskj
        loop1: do n=0,ntask-1
          if (tasks(n)%mytidj==j-1) then
            nlats_task(j) = tasks(n)%nlats
            exit loop1
          endif
        enddo loop1
      enddo
!
! Add poles to geographic grid (destination grid only):
      if (srcdes == 'des') then
        do n=0,ntask-1
          if (tasks(n)%mytidj==ntaskj-1.or.tasks(n)%mytidj==0) then
! north or south edge of task table: add 1 lat for pole
            nlats_task(tasks(n)%mytidj+1) = tasks(n)%nlats+1
          endif
        enddo
      endif
!
! Create 2d geographic grid:
! (minimum global lat index is 0, to include the poles)
      if (srcdes == 'src') then
        grid_out = ESMF_GridCreate1PeriDim(
     |    countsPerDEDim1=nlons_task, coordDep1=(/1/), 
     |    countsPerDEDim2=nlats_task, coordDep2=(/2/),
     |    indexflag=ESMF_INDEX_GLOBAL,
     |    minIndex=(/1,1/),rc=rc)
      else ! destination grid has poles
        grid_out = ESMF_GridCreate1PeriDim(
     |    countsPerDEDim1=nlons_task, coordDep1=(/1/), 
     |    countsPerDEDim2=nlats_task, coordDep2=(/2/),
     |    indexflag=ESMF_INDEX_GLOBAL,
     |    minIndex=(/1,0/),rc=rc)
      endif

      if (rc /=ESMF_SUCCESS) then
        write(6,"(/,'>>> create_geo_grid: error return from ',
     |    'ESMF_GridCreate1PeriDim: rc=',i4)") rc
        call shutdown('create_geo_grid')
      endif
!
! Allocate coordinates:
!
      call ESMF_GridAddCoord(grid_out,
     |  staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)

      if (rc /=ESMF_SUCCESS) then
        write(6,"(/,'>>> create_geo_grid: error return from ',
     |    'ESMF_GridAddCoord')")
        call shutdown('create_geo_grid')
      endif
!
! Get pointer and set geo grid longitude coordinates:
!
      call ESMF_GridGetCoord(grid_out, coordDim=1, localDE=0,
     |  computationalLBound=lbnd, computationalUBound=ubnd,
     |  farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

      if (rc /=ESMF_SUCCESS) then
        write(6,"(/,'>>> create_geo_grid: error return from ',
     |    'ESMF_GridGetCoord for longitude coords')")
        call shutdown('create_geo_grid')
      endif

      lbnd_lon = lbnd(1) ; ubnd_lon = ubnd(1)
      if (srcdes == 'src') then ! source grid excludes periodic points
        do i=lbnd_lon,ubnd_lon
          coordX(i) = glon(i)          ! 1 -> 72
        enddo
      else                             ! destination grid includes periodic points
        do i=lbnd_lon,ubnd_lon
          if (i <= 2) then             ! comments assume dlon = 5 deg
            coordX(i) = glon(nlon-2+i) ! i=1,2 <- 73,74 == glon(nlon-1,nlon)
          elseif (i >= nlon+3) then
            coordX(i) = glon(i-nlon-2) ! i=75,76 <- 3,4 == glon(1,2)
          else
            coordX(i) = glon(i-2)      ! i=3,74 == glon(1,nlon)
          endif
        enddo
      endif
!
! Get pointer and set geo grid latitude coordinates, including poles:
!
      call ESMF_GridGetCoord(grid_out, coordDim=2, localDE=0,
     |  computationalLBound=lbnd, computationalUBound=ubnd,
     |  farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

      if (rc /=ESMF_SUCCESS) then
        write(6,"(/,'>>> create_geo_grid: error return from ',
     |    'ESMF_GridGetCoord for latitude coords')")
        call shutdown('create_geo_grid')
      endif

      lbnd_lat = lbnd(1) ; ubnd_lat = ubnd(1)
      if (srcdes == 'src') then 
        do j=lbnd_lat,ubnd_lat
          coordY(j) = glat(j)
        enddo
      else
        do j=lbnd_lat,ubnd_lat
          if (j==0) then
            coordY(j) = -90.
          elseif (j==nlat+1) then
            coordY(j) = +90.
          else
            coordY(j) = glat(j)
          endif
        enddo
      endif

!     write(6,"('Created ESMF ',a,' geo grid: lbnd,ubnd_lon=',
!    |  2i4,'  lon0,1=',2i4,' lbnd,ubnd_lat=',2i4,'  lat0,1=',2i4)") 
!    |  srcdes,lbnd_lon,ubnd_lon,lon0,lon1,lbnd_lat,ubnd_lat,
!    |  lat0,lat1

!     write(6,"(a,' geo grid: my nlons=',i4,' coordX=',/,(8f9.2))") 
!    |  srcdes,ubnd_lon-lbnd_lon+1,coordX
!     write(6,"(a,' geo grid: my nlats=',i4,' coordY=',/,(8f9.2))") 
!    |  srcdes,ubnd_lat-lbnd_lat+1,coordY

      end subroutine create_geo_grid
!-----------------------------------------------------------------------
      subroutine esmf_create_geofield(field,grid,name,mlev0,mlev1)
!
! Create ESMF field (2d or 3d) on geo grid (will exclude periodic points)
! If mlev1 == 0, field is 2d (i,j), otherwise field is 3d,
!   and 3rd dimension is ungridded. Because 3rd dim must be the same
!   for geo and mag fields, they are both dimensioned mlev0:mlev1 
!
! Args:
      integer,intent(in) :: mlev0,mlev1 ! if mlev0,1==0, field is 2d (i,j)
      type(ESMF_Grid),intent(in) :: grid
      character(len=*),intent(in) :: name
      type(ESMF_Field),intent(out) :: field
!
! Local:
      integer :: rc
      type(ESMF_ArraySpec) :: arrayspec
!
! Create 3d field (i,j,k), with non-distributed vertical dimension:
      if (mlev1 > 0) then
        call ESMF_ArraySpecSet(arrayspec,3,ESMF_TYPEKIND_R8,rc=rc)
        if (rc /= ESMF_SUCCESS) call shutdown('ESMF_ArraySpecSet')

        field = ESMF_FieldCreate(grid, arrayspec, 
     |    ungriddedLBound=(/mlev0/), ungriddedUBound=(/mlev1/),
     |    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

        if (rc /= ESMF_SUCCESS) call shutdown('ESMF_FieldCreate 3d')
!       write(6,"('Created 3d (i,j,k) esmf geo field ',a)") 
!    |    name
!
! Create 2d field (i,j):
      else                ! create 2d field
        call ESMF_ArraySpecSet(arrayspec,2,ESMF_TYPEKIND_R8,rc=rc)
        if (rc /= ESMF_SUCCESS) call shutdown('ESMF_ArraySpecSet')
        field = ESMF_FieldCreate(grid, arrayspec, 
     |    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
        if (rc /= ESMF_SUCCESS) call shutdown('ESMF_FieldCreate 2d')
!       write(6,"('Created 2d (i,j) esmf geo field ',a)") name
      endif
      end subroutine esmf_create_geofield
!-----------------------------------------------------------------------
      subroutine esmf_create_magfield(field,grid,name,mlev0,mlev1)
!
! Create ESMF field (2d or 3d) on mag grid. This will include the
!   mag periodic point, which will be zero after regridding.
! If mlev1 == 0, field is 2d (i,j), otherwise field is 3d,
!   and 3rd dimension is ungridded
!
! Args:
      integer,intent(in) :: mlev0,mlev1
      type(ESMF_Grid),intent(in) :: grid
      character(len=*),intent(in) :: name
      type(ESMF_Field),intent(out) :: field
!
! Local:
      integer :: rc
      type(ESMF_ArraySpec) :: arrayspec
      type(ESMF_Array) :: array3d,array2d
      type(ESMF_DistGrid) :: distgrid
!
! Get necessary information from the mag grid:
      call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CENTER,
     |  distgrid=distgrid,rc=rc)
      if (rc /= ESMF_SUCCESS) call shutdown('ESMF_GridGet')
!
! Create 3d mag field (i,j,k), with non-distributed vertical dimension:
! (add periodic point in longitude with computationalEdgeUWidth)
!
      if (mlev1 > 0) then
        call ESMF_ArraySpecSet(arrayspec,3,ESMF_TYPEKIND_R8,rc=rc)
        if (rc /= ESMF_SUCCESS)call shutdown('ESMF_ArraySpecSet 3d mag')

        array3d = ESMF_ArrayCreate(arrayspec=arrayspec,
     |    distgrid=distgrid,computationalEdgeUWidth=(/1,0/),
     |    undistLBound=(/mlev0/),undistUBound=(/mlev1/),
     |    indexflag=ESMF_INDEX_GLOBAL,rc=rc)
        if (rc /= ESMF_SUCCESS) call shutdown('ESMF_ArrayCreate 3d mag')

        field = ESMF_FieldCreate(grid, array3d, 
     |    ungriddedLBound=(/mlev0/), ungriddedUBound=(/mlev1/),
     |    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
        if (rc /= ESMF_SUCCESS) call shutdown('ESMF_FieldCreate 3d mag')
!
! Create 2d mag field (i,j):
! (add periodic point in longitude with computationalEdgeUWidth)
!
      else                ! create 2d field
        call ESMF_ArraySpecSet(arrayspec,2,ESMF_TYPEKIND_R8,rc=rc)
        if (rc /= ESMF_SUCCESS)call shutdown('ESMF_ArraySpecSet 2d mag')

        array2d = ESMF_ArrayCreate(arrayspec=arrayspec,
     |    distgrid=distgrid,computationalEdgeUWidth=(/1,0/),
     |    indexflag=ESMF_INDEX_GLOBAL,rc=rc)
        if (rc /= ESMF_SUCCESS) call shutdown('ESMF_ArrayCreate 2d mag')

        field = ESMF_FieldCreate(grid, array2d, 
     |    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
        if (rc /= ESMF_SUCCESS) call shutdown('ESMF_FieldCreate 2d mag')

      endif
      end subroutine esmf_create_magfield
!-----------------------------------------------------------------------
      subroutine esmf_set3d_geo(fields,fnames,f,nf,ilev0,ilev1,
     |  ilon0,ilon1,ilat0,ilat1)
!
! Set values of a 3d ESMF field on geographic source grid, prior to
!   geographic to magnetic grid transformation. Periodic points are 
!   excluded. ESMF will calculate geographic poles via the route handle.
! Note dimension order changes from input (k,i,j) to output (i,j,k).
!
! Args:
      integer,intent(in) :: nf
      type(ESMF_Field) ,intent(in) :: fields(nf) ! esmf fields on geo grid
      character(len=*) ,intent(in) :: fnames(nf) ! field names
!
! f is input data on model subdomains (including periodic points)
! (note esmf source field excludes periodic points)
      integer,intent(in) :: ilev0,ilev1,ilon0,ilon1,ilat0,ilat1
      real,intent(in) :: f(ilev0:ilev1,ilon0:ilon1,ilat0:ilat1,nf)
!
! Local:
      integer,parameter :: mxf=8 ! for call by dynamo_inputs
      integer :: i,j,k,rc,n,istat
      integer,save :: ncalls=0
      integer,save :: lbnd(3),ubnd(3) ! lower,upper bounds of 3d field
!
! fptr is esmf pointer (i,j,k) to 3d field, set by this subroutine
      real(ESMF_KIND_R8),pointer :: fptr(:,:,:)
!
! ftmp is allocated only once (see below), so the save attribute here 
! is important.
      real,allocatable,save :: ftmp(:,:,:,:) ! esmf bounds, plus nf

      if (nf > mxf) then
        write(6,"('>>> esmf_set3d_geo: nf cannot be greater than ',
     |    'mxf: nf=',i4,' mxf=',i4)") nf,mxf
        call shutdown('esmf_set3d_geo nf')
      endif
!
! This routine is called every timestep from dynamo_inputs for 8 fields,
! and is called once per run for a single field from geo2mag_3d 
! (called by define_phim3d, which is called from rdsource.F). 
! Here, ftmp is allocated once per run w/ mxf=8. This is a local array, 
! so it must have the save attribute, see above.
!
      if (.not.allocated(ftmp)) then
!
! Get array bounds:
        call ESMF_FieldGet(fields(1),localDe=0,farrayPtr=fptr,
     |    computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
        if (rc /= ESMF_SUCCESS) then
          write(6,"('>>> esmf_set3d_geo: error from ESMF_FieldGet: rc=',
     |      i4)") rc
          call shutdown('esmf_set3d_geo: ESMF_FieldGet')
        endif
!
! Do the allocation:
        allocate(ftmp(lbnd(1):ubnd(1),lbnd(2):ubnd(2),lbnd(3):ubnd(3),
     |    mxf),stat=istat)
        if (istat /= 0) then
          write(6,"('>>> esmf_set3d_geo: error allocating ftmp')")
          call shutdown('alloc ftmp')
        endif
!       write(6,"('esmf_set3d_geo allocated ftmp: mxf=',i3)") mxf
      endif
!
! Fields loop:
      do n=1,nf
        ftmp(:,:,:,n) = 0.
!
! Set interior latitudes (ftmp(i,j,k,n) <- f(k,i,j,n))
! ftmp excludes periodic points.
!
        do j=lbnd(2),ubnd(2)     ! lat
          if (j /= 0 .and. j /= nlat+1) then ! interior latitudes (not poles)
            do i=lbnd(1),ubnd(1) ! lon
              do k=ilev0,ilev1   ! lbnd(3)==mlev0==zpbot_dyn, so cannot loop from lbnc(3) here 
                ftmp(i,j,k,n) = f(k,i+2,j,n)  ! data is at i+2 (3->74)
              enddo ! lev
            enddo   ! lon
          endif     ! poles or interior
        enddo       ! lat
      enddo         ! n=1,nf
!
! Get and set pointer to the field:
      do n=1,nf
        call ESMF_FieldGet(fields(n),localDe=0,farrayPtr=fptr,
     |    computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
        if (rc /= ESMF_SUCCESS) then
          write(6,"('>>> esmf_set3d_geo: error from ESMF_FieldGet: rc=',
     |      i4)") rc
          call shutdown('esmf_set3d_geo:ESMF_FieldGet')
        endif
        fptr(:,:,:) = ftmp(:,:,:,n)
      enddo ! n=1,nf
      deallocate(ftmp)
      end subroutine esmf_set3d_geo
!-----------------------------------------------------------------------
      subroutine esmf_set2d_geo(field,grid,fname,f,ilon0,ilon1,
     |  ilat0,ilat1)
!
! Set values of a 2d ESMF field on geographic source grid, prior to
!   geographic to magnetic grid transformation. (Essentially the same
!   as esmf_set3d_geo, except for 2d fields instead of 3d)
! Periodic points are excluded, geographic poles are at j==0 and nlat+1.
!
! Args:
      type(ESMF_Field) ,intent(in) :: field
      type(ESMF_Grid)  ,intent(in) :: grid
      character(len=*) ,intent(in) :: fname  ! field name
      integer      ,intent(in) :: ilon0,ilon1,ilat0,ilat1
      real         ,intent(in) :: f(ilon0:ilon1,ilat0:ilat1)
!
! Local:
      integer :: i,j,rc
      real(ESMF_KIND_R8),pointer :: fptr(:,:) ! i,j
      integer :: lbnd(2),ubnd(2)
!
! Get pointer to the field:
      call ESMF_FieldGet(field,localDe=0,farrayPtr=fptr,
     |  computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
      if (rc /= ESMF_SUCCESS) then
        write(6,"('>>> esmf_set2d_geo: error from ESMF_FieldGet: rc=',
     |    i4)") rc
        call shutdown('esmf_set2d_geo:ESMF_FieldGet')
      endif
!
      fptr(:,:) = 0. ! init
!
! Set interior latitudes (excluding poles):
      do j=lbnd(2),ubnd(2)
        if (j /= 0 .and. j /= nlat+1) then
          do i=lbnd(1),ubnd(1)
            fptr(i,j) = f(i+2,j) ! data is at f(i+2,j) (3->74)
          enddo
        endif       ! interior latitudes only
      enddo
!
!     write(6,"('esmf_set2d_geo field ',a,': lon bnds=',2i4,
!    |  ' lat bnds=',2i4,' 2d mnmx=',2e12.4)") 
!    |  fname,lbnd(1),ubnd(1),lbnd(2),ubnd(2),
!    |  minval(fptr(:,:)),maxval(fptr(:,:))

      end subroutine esmf_set2d_geo
!-----------------------------------------------------------------------
      subroutine esmf_set3d_mag(fields,fnames,f,nf,ilev0,ilev1,
     |  ilon0,ilon1,ilat0,ilat1)
!
! Set values of a 3d ESMF field on magnetic grid, prior to magnetic to 
! geographic grid transformation.
!
! Args:
      integer,intent(in) :: nf
      type(ESMF_Field) ,intent(in) :: fields(nf) ! esmf fields on mag grid
      character(len=*) ,intent(in) :: fnames(nf) ! field names
!
! f is input data on model subdomains:
! (note esmf source field excludes periodic points)
      integer,intent(in) :: ilev0,ilev1,ilon0,ilon1,ilat0,ilat1
      real,intent(in) :: f(ilon0:ilon1,ilat0:ilat1,ilev0:ilev1,nf)
!
! Local:
      integer :: i,j,k,rc,n,istat
      integer,save :: ncalls=0,lbnd(3),ubnd(3) ! lower,upper bounds of 3d field
!
! fptr is esmf pointer (i,j,k) to 3d field, set by this subroutine
      real(ESMF_KIND_R8),pointer :: fptr(:,:,:)
!
! Fields loop:
      do n=1,nf
        call ESMF_FieldGet(fields(n),localDe=0,farrayPtr=fptr,
     |    computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
        fptr(:,:,:) = 0.
!
! Set ESMF pointer:
!
        do j=lbnd(2),ubnd(2)      ! lat
          do i=lbnd(1),ubnd(1)    ! lon
            do k=lbnd(3),ubnd(3)  ! lev
! 0-15:forrtl: severe (408): fort: (3): Subscript #3 of the array F has value -2 which is less than the lower bound of 1
              fptr(i,j,k) = f(i,j,k,n)
            enddo ! mlev
          enddo   ! mlon
        enddo     ! mlat
      enddo       ! n=1,nf
      end subroutine esmf_set3d_mag
!-----------------------------------------------------------------------
#if defined(INTERCOMM) || defined(CISMAH)
      subroutine esmf_set2d_mag(fields,fnames,f,nf,
     |  ilon0,ilon1,ilat0,ilat1)
!
! Set values of a 2d ESMF field on magnetic grid, prior to magnetic to 
! geographic grid transformation.
!
! Args:
      integer,intent(in) :: nf
      type(ESMF_Field) ,intent(in) :: fields(nf) ! esmf fields on mag grid
      character(len=*) ,intent(in) :: fnames(nf) ! field names
!
! f is input data on model subdomains:
! (note esmf source field excludes periodic points)
      integer,intent(in) :: ilon0,ilon1,ilat0,ilat1
      real,intent(in) :: f(ilon0:ilon1,ilat0:ilat1,nf)
!
! Local:
      integer :: i,j,rc,n,istat
      integer,save :: ncalls=0,lbnd(2),ubnd(2) ! lower,upper bounds of 2d field
!
! fptr is esmf pointer (i,j) to 2d field, set by this subroutine
      real(ESMF_KIND_R8),pointer :: fptr(:,:)
!
! Fields loop:
      do n=1,nf
        call ESMF_FieldGet(fields(n),localDe=0,farrayPtr=fptr,
     |    computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
        fptr(:,:) = 0.
!
! Set ESMF pointer:
!
        do j=lbnd(2),ubnd(2)      ! lat
          do i=lbnd(1),ubnd(1)    ! lon
            fptr(i,j) = f(i,j,n)
          enddo   ! mlon
        enddo     ! mlat
      enddo       ! n=1,nf
      end subroutine esmf_set2d_mag
#endif
!-----------------------------------------------------------------------
      subroutine esmf_get_3dfield(field, fptr, name)
!
! Get pointer to 3d esmf field (i,j,k):
! In pdynamo.F: call esmf_get_3dfield(mag_ped   ,ped_mag,   "PED     ")
!
! Args:
      type(ESMF_field),intent(in) :: field
      real,pointer,dimension(:,:,:),intent(out) :: fptr
      character(len=*),intent(in) :: name
!
! Local:
      integer :: rc,k,lbnd(3),ubnd(3)
      character(len=80) :: errmsg

      call ESMF_FieldGet(field,localDe=0,farrayPtr=fptr,
     |  computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
      if (rc /= ESMF_SUCCESS) then
        write(errmsg,"('esmf_get_field 3d field ',a)") trim(name)
        call shutdown(errmsg)
      endif

!     if (name=='PHIG3D  ') then
!       write(6,"('esmf_get_3dfield: ',a,' lbnd=',3i4,' ubnd=',3i4,
!    |    ' min,max=',2e12.4)") name,lbnd,ubnd,minval(fptr),maxval(fptr)
!     endif

      end subroutine esmf_get_3dfield
!-----------------------------------------------------------------------
      subroutine esmf_get_2dfield(field, fptr, name)
!
! Get pointer to 2d esmf field (i,j):
!
! Args:
      type(ESMF_field),intent(in) :: field
      real,pointer,dimension(:,:),intent(out) :: fptr
      character(len=*),intent(in) :: name
!
! Local:
      integer :: rc
      character(len=80) :: errmsg

!     fptr = 0. ! not allowed

      call ESMF_FieldGet(field,localDe=0,farrayPtr=fptr,rc=rc)
      if (rc /= ESMF_SUCCESS) then
        write(errmsg,"('esmf_get_field 2d field ',a)") trim(name)
        call shutdown(errmsg)
      endif

!     write(6,"('esmf_get_2dfield: ',a,' min,max=',2e12.4)") 
!    |  name,minval(fptr),maxval(fptr)

!     write(6,"('esmf_get_2dfield: ',a8,' mlon0,1=',2i4,' nmlons=',i4,
!    |  ' size(fptr,1:2)=',2i4,' fptr min,max=',2e12.4)") name,
!    |  mlon0,mlon1,mlon1-mlon0+1,size(fptr,1),size(fptr,2),
!    |  minval(fptr),maxval(fptr)

      end subroutine esmf_get_2dfield
!-----------------------------------------------------------------------
      subroutine esmf_regrid(srcfield,dstfield,direction,ndim)
!
! Args:
      integer :: ndim
      type(ESMF_Field),intent(inout) :: srcfield,dstfield
      character(len=*),intent(in) :: direction
!
! Local:
      integer :: rc
      type(ESMF_RouteHandle) :: routehandle
!
! Direction is either geo2mag or mag2geo. 
! Use corresponding route handle (module data)
!
      select case(trim(direction))
        case ('geo2mag')
          routehandle = routehandle_geo2mag
          if (ndim==2) then
            routehandle = routehandle_geo2mag_2d
!           call ESMF_FieldRegrid(srcfield,dstfield,routehandle,rc=rc)
!
! Do sparse matrix multiply for 2d geo2mag. Note routehandle_geo2mag
! may be made from apex or esmf indices and weights (see call 
! ESMF_FieldSMMStore in sub esmf_init)
!
            call ESMF_FieldSMM(srcfield,dstfield,routehandle,rc=rc)

            if (rc /= ESMF_SUCCESS) then
              write(6,"(/,'>>> esmf_regrid: error return from ',
     |          'ESMF_FieldSMM for 2d ',a,': rc=',i4)") trim(direction),
     |          rc
              call shutdown('ESMF_FieldRegrid')
            endif
          else
!
! Do sparse matrix multiply for 3d geo2mag. Note routehandle_geo2mag
! may be made from apex or esmf indices and weights (see call 
! ESMF_FieldSMMStore in sub esmf_init)
!
            routehandle = routehandle_geo2mag
            call ESMF_FieldSMM(srcfield,dstfield,routehandle,rc=rc)
            if (rc /= ESMF_SUCCESS) then
              write(6,"(/,'>>> esmf_regrid: error return from ',
     |          'ESMF_FieldSMM for 3d ',a,': rc=',i4)") trim(direction),
     |          rc
              call shutdown('ESMF_FieldSMM')
            endif
          endif
!
! Do sparse matrix multiply for 3d mag2geo. Note routehandle_mag2geo
! may be made from apex or esmf indices and weights (see call 
! ESMF_FieldSMMStore in sub esmf_init)
!
        case ('mag2geo')
          routehandle = routehandle_mag2geo
          if (ndim==2) then
            routehandle = routehandle_mag2geo_2d
          endif
          call ESMF_FieldSMM(srcfield,dstfield,routehandle,rc=rc)
          if (rc /= ESMF_SUCCESS) then
            write(6,"(/,'>>> esmf_regrid: error return from ',
     |        'ESMF_FieldSMM for 3d ',a,': rc=',i4)") trim(direction),
     |        rc
            call shutdown('ESMF_FieldSMM')
          endif
        case default
          write(6,"('>>> esmf_regrid: bad direction=',a)")
     |      trim(direction)
          call shutdown('esmf_regrid')
      end select
      end subroutine esmf_regrid
!-----------------------------------------------------------------------
      end module my_esmf
