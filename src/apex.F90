module apex
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! April, 2013: B. Foster (NCAR/HAO)
!
! This is a refactored version of the legacy apex code, originally written
! by Art Richmond and Roy Barnes, and others in the 1995-2000 timeframe.
! This new version is written in free-format fortran90. Subroutines and
! module data may be use-associated from this module. 
!
! Original reference for the legacy code:
!          Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
!          Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995. 
!
! This code should produce near-identical results as the legacy code, altho 
! the refactored version does not provide all subroutines and options available 
! in the old code, notably the ability to write and read-back an external file.
! 
! A typical calling sequence for a code calling this module is as follows:
!
! subroutine ggrid (legacy SUBROUTINE GGRID): 
!   Make a global lat,lon,alt grid for use in later calls (optional)
!
! subroutine apex_mka (legacy SUBROUTINE APXMKA): 
!   Make magnetic arrays x,y,z,v for use in later routines
!   (geographic lat,lon grid and altitudes are input) 
!   (This must be called before apex_mall and apex_q2g)
! 
! subroutine apex_mall (legacy ENTRY APXMALL): 
!   Calculate modified Apex coordinates and other magnetic field parameters
!   (usually called from lat,lon,alt nested loop)
! 
! subroutine apex_q2g (legacy ENTRY APXQ2G): 
!   Convert from quasi-dipole to geodetic coordinates
!   (usually called from lat,lon,alt nested loop)
!
  implicit none
  real,parameter :: re = 6371.2, eps = 1.e-5

  real,allocatable,save :: &
    xarray(:,:,:), & ! cos(quasi-dipole latitude)*cos(apex longitude)
    yarray(:,:,:), & ! cos(quasi-dipole latitude)*sin(apex longitude)
    zarray(:,:,:), & ! sin(quasi-dipole latitude)
    varray(:,:,:)    ! (VMP/VP)*((RE+ALT)/RE)**2
!
! This grid (geolat,geolon,geoalt is equivalent to gdlat,gdlon,gdalt, 
! as passed to apex_mka.
!
  integer :: nglat,nglon,ngalt
  real,allocatable,save :: geolat(:), geolon(:), geoalt(:)    

  integer,parameter :: nmax=13
  integer,parameter :: ncoef = nmax*nmax + 2*nmax + 1 ! 196
  real,dimension(ncoef) :: &
    gb, & ! Coefficients for magnetic field calculation
    gv    ! Coefficients for magnetic potential calculation
!
  real :: &
    rtd,  & ! radians to degrees
    dtr,  & ! degrees to radians
    pola    ! pole angle (deg); when geographic lat is poleward of pola,
            ! x,y,z,v arrays are forced to be constant (pola=89.995) 

  real,parameter ::        & ! Formerly common /APXCON/
    req  = 6378.160,       & ! Equatorial earth radius
    precise = 7.6e-11,     & ! Precision factor
    glatlim = 89.9,        & ! Limit above which gradients are recalculated
    xmiss = -32767.
!
! colat,elon,vp,ctp,stp were in two commons in legacy code:
! /APXDIPL/ and /DIPOLE/. Need to check if these need to be separated.
!
  real ::  & ! Formerly /APXDIPL/ and /DIPOLE/
    colat, & ! Geocentric colatitude of geomagnetic dipole north pole (deg)
    elon,  & ! East longitude of geomagnetic dipole north pole (deg)
    vp,    & ! Magnitude, in T.m, of dipole component of magnetic
             ! potential at geomagnetic pole and geocentric radius re
    ctp,stp
!
  real ::  & ! Formerly /FLDCOMD/
    bx,    & ! X comp. of field vector at the current tracing point (Gauss)
    by,    & ! Y comp. of field vector at the current tracing point (Gauss)
    bz,    & ! Z comp. of field vector at the current tracing point (Gauss)
    bb       ! Magnitude of field vector at the current tracing point (Gauss)

  real ::      & ! Formerly /APXIN/
    yapx(3,3)    ! Matrix of cartesian coordinates (loaded columnwise) 
!
! /ITRA/ was only in subs linapx and itrace, so can probably be removed from module data
!
  integer ::   & ! Formerly /ITRA/ 
    nstp         ! Step count. Incremented in sub linapx.
  real    ::   & 
    y(3),      & ! Array containing current tracing point cartesian coordinates.
    yp(3),     & ! Array containing previous tracing point cartesian coordinates.
    sgn,       & ! Determines direction of trace. Set in subprogram linapx
    ds           ! Step size (Km) Computed in subprogram linapx.

  real ::         & ! limits beyond which east-west gradients are computed 
    glatmn,glatmx   ! differently to avoid potential underflow (apex_mka)

contains
!-----------------------------------------------------------------------
subroutine ggrid(nvert,glatmin,glatmax,glonmin,glonmax,altmin,altmax, &
                 gplat,gplon,gpalt,mxlat,mxlon,mxalt,nlat,nlon,nalt)
!
! Given desired range of geographic latitude, longitude and altitude, 
! choose an appropriate grid that can be used in subsequent calls to 
! subs apex_mka, apex_mall, apex_q2g.
!
! Input args:
  integer,intent(in) :: nvert,mxlat,mxlon,mxalt
  real,intent(in) :: glatmin,glatmax,glonmin,glonmax,altmin,altmax
!
! Output args:
  integer,intent(out) :: nlat,nlon,nalt
  real,intent(out) :: gplat(mxlat),gplon(mxlon),gpalt(mxalt)
!
! Local:
  real :: dlon,dlat,diht,dnv,glonmaxx,x
  integer :: nlatmin,nlatmax,nlonmin,nlonmax,naltmin,naltmax
  integer :: i,j,k,kk
!
! Check inputs:
  if (glatmin > glatmax) then
    write(6,"('>>> ggrid: glatmin=',f9.2,' must be <= glatmax=',f9.2)") &
      glatmin,glatmax
    stop 'ggrid'
  endif
  if (glonmin > glonmax) then
    write(6,"('>>> ggrid: glonmin=',f9.2,' must be <= glonmax=',f9.2)") &
      glonmin,glonmax
    stop 'ggrid'
  endif
  if (altmin > altmax) then
    write(6,"('>>> ggrid: altmin=',f9.2,' must be <= altmax=',f9.2)") &
      altmin,altmax
    stop 'ggrid'
  endif
!
! Init outputs:
  nlat = 0 ; nlon = 0 ; nalt = 0
  gplat = 0. ; gplon = 0. ; gpalt = 0.
!
  dnv = float(nvert)
  dlon = 360. / (5.*dnv)
  dlat = 180. / (3.*dnv)
  diht = 1.   / dnv

  nlatmin = max(int((glatmin+90.)/dlat),0)
  nlatmax = min(int((glatmax+90.)/dlat+1.),3*nvert)
  nlonmin = max(int((glonmin+180.)/dlon),0)
 
  glonmaxx = min(glonmax,glonmin+360.)
  nlonmax = min(int((glonmaxx+180.)/dlon+1.),10*nvert)
    
  x = re/(re+altmax)/diht-eps
  naltmin = max(x,1.)
  naltmin = min(naltmin,nvert-1)
  x = re/(re+altmin)/diht+eps
  i = x + 1.
  naltmax = min(i,nvert)

  nlat = nlatmax - nlatmin + 1
  nlon = nlonmax - nlonmin + 1
  nlon = min(nlon,5*nvert+1)
  nalt = naltmax - naltmin + 1

  do j=1,nlat
    gplat(j) = dlat*float(nlatmin+j-1) - 90.
  enddo
  do i=1,nlon
    gplon(i) = dlon*float(nlonmin+i-1) - 180.
  enddo
  do k=1,nalt
    kk = naltmax - k +1
    gpalt(k) = re*(float(nvert-kk) - eps) / (float(kk)+eps)
  enddo
  if (gplon(nlon-1) >= glonmax) nlon = nlon-1
  gpalt(1) = max(gpalt(1),0.)

! write(6,"('ggrid: nlat=',i4,' gplat=',/,(6f9.2))") nlat,gplat
! write(6,"('ggrid: nlon=',i4,' gplon=',/,(6f9.2))") nlon,gplon
! write(6,"('ggrid: nalt=',i4,' gpalt=',/,(6f9.2))") nalt,gpalt

end subroutine ggrid
!-----------------------------------------------------------------------
subroutine apex_mka(date,gplat,gplon,gpalt,nlat,nlon,nalt,ier)
!
! Given a 3d lat,lon,altitude grid, calculate x,y,z,v arrays in module
! data above. These arrays are used later for calculating quantities
! involving gradients of Apex coordinates, such as base vectors in the
! Modified-Apex and Quasi-Dipole systems.
!
! This defines module 3d data xarray,yarray,zarray,varray
!
! Input args:
  real,intent(in) :: date              ! year and fraction
  integer,intent(in) :: nlat,nlon,nalt ! dimensions of 3d grid
  real,intent(inout) :: gplat(nlat),gplon(nlon),gpalt(nalt)
!
! Output args:
  integer,intent(out) :: ier
!
! Local:
  integer :: i,j,k,kpol,istat
  real :: reqore,rqorm1,cp,ct,st,sp,stmcpm,stmspm,ctm
  real :: aht,alat,phia,bmag,xmag,ymag,zdown,vmp ! apex_sub output
  real :: vnor,rp,reqam1,a,slp,clp,phiar
!
! Some parts of the legacy apex code use constants to set dtr,rtd,
! other parts use rtd=45./atan(1.), dtr=1./rtd. Differences are
! on the order of 1.e-18 to 1.e-14. Here, the atan method is used. 
!
!  rtd  = 5.72957795130823E1
!  dtr  = 1.745329251994330E-2
!
   rtd  = 45./atan(1.)
   dtr  = 1./rtd
!
! pola:
!   Pole angle (deg); when the geographic latitude is poleward of POLA, 
!   X,Y,Z,V are forced to be constant for all longitudes at each altitude.  
!   This makes POLA = 89.995
!
  pola = 90.-sqrt(precise)*rtd    ! Pole angle (deg)

  ier = 0
!
! Allocate 3d x,y,z,v arrays:
! These are not deallocated by this module. They can be deallocated
! by the calling program following the last call to the apex subs.
!
  if (.not.allocated(xarray)) then
    allocate(xarray(nlat,nlon,nalt),stat=istat)
    if (istat /= 0) stop 'allocate xarray'
    xarray = 0.
  endif
  if (.not.allocated(yarray)) then
    allocate(yarray(nlat,nlon,nalt),stat=istat)
    if (istat /= 0) stop 'allocate yarray'
    yarray = 0.
  endif
  if (.not.allocated(zarray)) then
    allocate(zarray(nlat,nlon,nalt),stat=istat)
    if (istat /= 0) stop 'allocate zarray'
    zarray = 0.
  endif
  if (.not.allocated(varray)) then
    allocate(varray(nlat,nlon,nalt),stat=istat)
    if (istat /= 0) stop 'allocate varray'
    varray = 0.
  endif
!
! Set geographic grid in module data for later reference:
! (these also are not deallocated by this module)
!
  nglon=nlon ; nglat=nlat ; ngalt=nalt
  allocate(geolat(nglat),stat=istat)
  allocate(geolon(nglon),stat=istat)
  allocate(geoalt(ngalt),stat=istat)
  geolat(:) = gplat(:)
  geolon(:) = gplon(:)
  geoalt(:) = gpalt(:)
!
! Set coefficients gb,gv (module data) for requested year:
!
  call cofrm(date)

! write(6,"('apex_mka after cofrm: ncoef=',i4,' gb=',/,(6f12.3))") ncoef,gb
! write(6,"('apex_mka after cofrm: ncoef=',i4,' gv=',/,(6f12.3))") ncoef,gv

  call dypol(colat,elon,vp)

  ctp = cos(colat*dtr)
  stp = sin(colat*dtr)

  reqore = req/re
  rqorm1 = reqore-1.

  do j=1,nlat
    ct = sin(gplat(j)*dtr)
    st = cos(gplat(j)*dtr)
    kpol = 0
    if (abs(gplat(j)) > pola) kpol = 1
    do i=1,nlon
      if (kpol==1.and.i > 1) then
        xarray(j,i,:) = xarray(j,1,:)
        yarray(j,i,:) = yarray(j,1,:)
        zarray(j,i,:) = zarray(j,1,:)
        varray(j,i,:) = varray(j,1,:)
        cycle
      endif  
      cp = cos((gplon(i)-elon)*dtr)
      sp = sin((gplon(i)-elon)*dtr)
!
!  ctm   is pseudodipole component of z
! -ctm   is pseudodipole component of v
!  stmcpm is pseudodipole component of x
!  stmspm is pseudodipole component of y
!
      ctm = ctp*ct + stp*st*cp
      stmcpm = st*ctp*cp - ct*stp
      stmspm = st*sp
      do k=1,nalt
        call apex_sub(date,gplat(j),gplon(i),gpalt(k),&
          aht,alat,phia,bmag,xmag,ymag,zdown,vmp)

        vnor = vmp/vp
        rp = 1. + gpalt(k)/re
        varray(j,i,k) = vnor*rp*rp + ctm
        reqam1 = req*(aht-1.)
        slp = sqrt(max(reqam1-gpalt(k),0.)/(reqam1+re))
!
! Reverse sign of slp in southern magnetic hemisphere
!
        if (zdown.lt.0.) slp = -slp
        clp = sqrt (rp/(reqore*aht-rqorm1))
        phiar = phia*dtr
        xarray(j,i,k) = clp*cos (phiar) - stmcpm
        yarray(j,i,k) = clp*sin (phiar) - stmspm
        zarray(j,i,k) = slp - ctm
      enddo ! k=1,nalt
    enddo ! i=1,nlon
  enddo ! j=1,nlat
!
! Establish for this grid polar latitude limits beyond which east-west
! gradients are computed differently to avoid potential underflow
! (glatmx,glatmn are in module data, glatlim is parameter constant)
!
  glatmx = max( glatlim,gplat(nlat-2))
  glatmn = min(-glatlim,gplat(2))

end subroutine apex_mka
!-----------------------------------------------------------------------
subroutine apex_mall(glat,glon,alt,hr, b,bhat,bmag,si,alon,xlatm,vmp,w,&
  d,be3,sim,d1,d2,d3,e1,e2,e3,xlatqd,f,f1,f2,ier)
!
! Compute Modified Apex coordinates, quasi-dipole coordinates,
! base vectors and other parameters by interpolation from
! precalculated arrays. Subroutine apex_mka must be called
! before calling this subroutine.
!
! Args:
  real,intent(in)  :: & ! Input
    glat             ,& ! Geographic (geodetic) latitude (deg)
    glon             ,& ! Geographic (geodetic) longitude (deg)
    alt              ,& ! Altitude (km)
    hr                  ! Reference altitude (km)

  real,intent(out) :: & ! Output
    b(3)             ,& ! Magnetic field components (east, north, up), in nT    
    bhat(3)          ,& ! components (east, north, up) of unit vector along 
                        ! geomagnetic field direction
    bmag             ,& ! Magnitude of magnetic field (nT)
    si               ,& ! sin(i)
    alon             ,& ! Apex longitude = modified apex longitude = 
                        ! quasi-dipole longitude (deg)
    xlatm            ,& ! Modified Apex latitude (deg)
    vmp              ,& ! Magnetic potential (T.m)
    w                ,& ! W of Richmond reference above, in km**2 /nT (i.e., 10**15 m**2 /T)
    d                ,& ! D of Richmond reference above
    be3              ,& ! B_e3 of reference above (= Bmag/D), in nT
    sim              ,& ! sin(I_m) described in Richmond reference above
    xlatqd           ,& ! Quasi-dipole latitude (deg)
    f                ,& ! F described in ref above for quasi-dipole coordinates
    f1(2),f2(2)         ! Components (east, north) of base vectors
!
  real,dimension(3),intent(out) :: d1,d2,d3,e1,e2,e3 ! Components of base vectors
  integer,intent(out) :: ier ! error return
!
! Local:
  real :: glonloc,cth,sth,glatx,clm,r3_2
  real :: fx,fy,fz,fv
  real :: dfxdth,dfydth,dfzdth,dfvdth, &
          dfxdln,dfydln,dfzdln,dfvdln, &
          dfxdh ,dfydh ,dfzdh ,dfvdh
  real,dimension(3) :: gradx,grady,gradz,gradv, grclm,clmgrp,rgrlp
  real ::                        & ! dummies for polar calls to intrp
    fxdum,fydum,fzdum,fvdum,     &
    dmxdth,dmydth,dmzdth,dmvdth, &
    dmxdh,dmydh,dmzdh,dmvdh
!
! Init:
!
  ier = 0
  glonloc = glon

  call intrp (glat,glonloc,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
             fx,fy,fz,fv,                                               &
             dfxdth,dfydth,dfzdth,dfvdth,                               &
             dfxdln,dfydln,dfzdln,dfvdln,                               &
             dfxdh ,dfydh ,dfzdh ,dfvdh, ier)

  if (ier /= 0) then
    call setmiss(xmiss,xlatm,alon,vmp,b,bmag,be3,sim,si,f,d,w, &
      bhat,d1,d2,d3,e1,e2,e3,f1,f2)
    write(6,"('apex_mall called setmiss: glat,glon,alt=',3f12.3)") &
      glat,glon,alt
    return
  endif

  call adpl(glat,glonloc,cth,sth,fx,fy,fz,fv, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)

  call gradxyzv(alt,cth,sth, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)
!
! If the point is very close to either the North or South
! geographic pole, recompute the east-west gradients after
! stepping a small distance from the pole.
!
  if (glat > glatmx .or. glat < glatmn) then
    glatx = glatmx
    if (glat < 0.) glatx = glatmn

    call intrp (glatx,glonloc,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
               fxdum,fydum,fzdum,fvdum,                                    &
               dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,           &
               dfvdln,dmxdh,dmydh,dmzdh,dmvdh, ier)

    call adpl(glatx,glonloc,cth,sth,fxdum,fydum,fzdum,fvdum, &
      dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,dfvdln)

    call grapxyzv(alt,cth,sth,dfxdln,dfydln,dfzdln,dfvdln, &
      gradx,grady,gradz,gradv)
  endif

  call gradlpv(hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
    xlatm,alon,vmp,grclm,clmgrp,xlatqd,rgrlp,b,clm,r3_2)

  call basevec(hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
               bmag,sim,si,f,d,w,bhat,d1,d2,d3,e1,e2,e3,f1,f2)

  be3 = bmag/d
  ier = 0

end subroutine apex_mall
!-----------------------------------------------------------------------
subroutine apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ier)
!
! Convert from quasi-dipole to geodetic coordinates. This subroutine
! (input magnetic, output geodetic) is the functional inverse of 
! subroutine apex_mall (input geodetic, output magnetic). Sub apex_mka
! must be called before this routine.
!
! Args:
  real,intent(in) ::  & ! inputs
    qdlat,            & ! quasi-dipole latitude (deg)
    qdlon,            & ! quasi-dipole longitude (deg)
    alt                 ! altitude (km)

  real,intent(out) :: & ! outputs
    gdlat,            & ! geodetic latitude (deg)
    gdlon               ! geodetic longitude (deg)
  integer,intent(out) :: ier ! error return
!
! Local:
  real :: x0,y0,z0,xnorm,xdif,ydif,zdif,dist2,hgrd2e,hgrd2n,hgrd2,&
    angdist,distlon,glatx,cal,sal,coslm,slm,cad,sad,slp,clm2,slm2,&
    sad2,cal2,clp2,clp,dylon
  real :: ylat,ylon ! first guess output by gm2gc, input to intrp
  integer :: iter
  integer,parameter :: niter=20
  real ::                        & ! output of sub intrp
    fx,fy,fz,fv,                 & ! interpolated values of x,y,z,v
    dfxdth,dfydth,dfzdth,dfvdth, & ! derivatives of x,y,z,v wrt colatitude
    dfxdln,dfydln,dfzdln,dfvdln, & ! derivatives of x,y,z,v wrt longitude
    dfxdh ,dfydh ,dfzdh ,dfvdh     ! derivatives of x,y,z,v wrt altitude
  real ::                        & ! dummies for polar calls to intrp
    fxdum,fydum,fzdum,fvdum,     &
    dmxdth,dmydth,dmzdth,dmvdth, &
    dmxdh,dmydh,dmzdh,dmvdh
  real :: cth,sth  ! output of adpl
  character(len=5) :: edge

  ier = 0 ; gdlat = 0. ; gdlon = 0.
!
! Determine quasi-cartesian coordinates on a unit sphere of the
! desired magnetic lat,lon in quasi-dipole coordinates.
!
  x0 = cos (qdlat*dtr) * cos (qdlon*dtr)
  y0 = cos (qdlat*dtr) * sin (qdlon*dtr)
  z0 = sin (qdlat*dtr)
!
! Initial guess:  use centered dipole, convert to geocentric coords
!
  call gm2gc (qdlat,qdlon,ylat,ylon)
!
! Iterate until (angular distance)**2 (units: radians) is within
! precise of location (qdlat,qdlon) on a unit sphere. 
! (precise is a parameter in module data)
!
  do iter=1,niter
!
! geolat,lon,alt and nglat,lon,alt are in module data (set by apex_mka)
!
    call intrp (ylat,ylon,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
               fx,fy,fz,fv,                                            &
               dfxdth,dfydth,dfzdth,dfvdth,                            &
               dfxdln,dfydln,dfzdln,dfvdln,                            &
               dfxdh ,dfydh ,dfzdh ,dfvdh, ier)
    if (ier /= 0) then
      write(6,"('>>> apex_q2g error from intrp')")
      stop 'qpex_q2g intrp'
    endif
!
!  Add-back of pseudodipole component to x,y,z,v and their derivatives.
!
    call adpl(ylat,ylon,cth,sth,fx,fy,fz,fv, &
      dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)
    distlon = cos(ylat*dtr)

    if (ylat > glatmx .or. ylat < glatmn) then ! glatmx,glatmn are module data
      glatx = glatmx
      if (ylat.lt.0.) glatx = glatmn
      distlon = cos (glatx*dtr)
      call intrp (glatx,ylon,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
                 fxdum,fydum,fzdum,fvdum,                                 &
                 dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,        &
                 dfvdln,dmxdh,dmydh,dmzdh,dmvdh, ier)
      if (ier /= 0) then
        write(6,"('>>> apex_q2g error from polar intrp')")
        stop 'qpex_q2g intrp'
      endif

      call adpl(glatx,ylon,cth,sth,fxdum,fydum,fzdum,fvdum, &
        dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,dfvdln)
    endif
!
! At this point, FX,FY,FZ are approximate quasi-cartesian
! coordinates on a unit sphere for the quasi-dipole coordinates
! corresponding to the geodetic coordinates YLAT, YLON.
! Normalize the vector length of (FX,FY,FZ) to unity using XNORM
! so that the resultant vector can be directly compared with the
! target vector (X0,Y0,Z0).
!
    xnorm = sqrt(fx*fx + fy*fy + fz*fz)
    xdif = fx/xnorm - x0
    ydif = fy/xnorm - y0
    zdif = fz/xnorm - z0
!
! dist2 = square of distance between normalized (fx,fy,fz) and x0,y0,z0.
!
    dist2 = xdif*xdif + ydif*ydif + zdif*zdif

    if (dist2 <= precise) then
      ier = 0
      gdlat = ylat
      gdlon = ylon
      return
    endif
!
! hgrd2* = one-half of east or north gradient of dist2 on unit sphere.
!
    hgrd2e =  (xdif*dfxdln + ydif*dfydln + zdif*dfzdln)/distlon
    hgrd2n = -(xdif*dfxdth + ydif*dfydth + zdif*dfzdth)
    hgrd2  = sqrt(hgrd2e*hgrd2e + hgrd2n*hgrd2n)
!
! angdist = magnitude of angular distance to be moved for new guess
!           of ylat, ylon.
!
    angdist = dist2/hgrd2
!
! Following spherical trigonometry moves ylat,ylon to new location,
! in direction of grad(dist2), by amount angdist.
!
    cal = -hgrd2n/hgrd2
    sal = -hgrd2e/hgrd2 
    coslm = cos(ylat*dtr)
    slm = sin(ylat*dtr)
    cad = cos(angdist) 
    sad = sin(angdist)
    slp = slm*cad + coslm*sad*cal

    clm2 = coslm*coslm
    slm2 = slm*slm
    sad2 = sad*sad
    cal2 = cal*caL
    clp2 = clm2 + slm2*sad2 - 2.*slm*cad*coslm*sad*cal -clm2*sad2*cal2
    clp = sqrt (max(0.,clp2))
    ylat = atan2(slp,clp)*rtd
!
! Restrict latitude iterations to stay within the interpolation grid
! limits, but let intrp find any longitude exceedence.  This is only
! an issue when the interpolation grid does not cover the entire
! magnetic pole region.
!
    ylat = min(ylat,geolat(nglat))
    ylat = max(ylat,geolat(1))
    dylon = atan2 (sad*sal,cad*coslm-sad*slm*cal)*rtd
    ylon = ylon + dylon
    if (ylon > geolon(nglon)) ylon = ylon - 360.
    if (ylon < geolon(1))     ylon = ylon + 360.

  enddo ! iter=1,niter

  write(6,"('>>> apex_q2g: ',i3,' iterations only reduced the angular')") niter
  write(6,"('              difference to ',f10.5,' degrees, where test criterion')") &
    sqrt(dist2)*rtd
  write(6,"('              is ',f10.5,' degrees.')") sqrt(precise)*rtd
  edge = '     '
  if (ylat == geolat(nglat)) edge = 'north'
  if (ylat == geolat(1))     edge = 'south'
  if (edge /= '     ') then
    write(6,"('Coordinates are on the ',a,' edge of the interpolation grid ')") edge
    write(6,"('and latitude is constrained to stay within grid limits when iterating.')") 
  endif
  ier = 1

end subroutine apex_q2g
!-----------------------------------------------------------------------
subroutine gradxyzv(alt,cth,sth, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)
!
! Calculates east,north,up components of gradients of x,y,z,v in
! geodetic coordinates.  All gradients are in inverse km.  Assumes
! flatness of 1/298.25 and equatorial radius (REQ) of 6378.16 km.
! 940803 A. D. Richmond
!
! Args:
  real,intent(in) :: alt,cth,sth
  real,dimension(3),intent(out) :: gradx,grady,gradz,gradv
  real,intent(in) ::             &
    dfxdth,dfydth,dfzdth,dfvdth, &
    dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh,dfydh,dfzdh,dfvdh
!
! Local:
  real :: d,d2,rho,dddthod,drhodth,dzetdth,ddisdth

!
!          40680925. = req**2 (rounded off)
!          272340.   = req**2 * E2, where E2 = (2. - 1./298.25)/298.25
!                      is square of eccentricity of ellipsoid.
!
  d2 = 40680925.e0 - 272340.e0*cth*cth
  d = sqrt(d2)
  rho = sth*(alt + 40680925.e0/d)
  dddthod = 272340.e0*cth*sth/d2
  drhodth = alt*cth + (40680925.e0/d)*(cth-sth*dddthod)
  dzetdth =-alt*sth - (40408585.e0/d)*(sth+cth*dddthod)
  ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)

  gradx(1) = dfxdln/rho
  grady(1) = dfydln/rho
  gradz(1) = dfzdln/rho
  gradv(1) = dfvdln/rho

  gradx(2) = -dfxdth/ddisdth
  grady(2) = -dfydth/ddisdth
  gradz(2) = -dfzdth/ddisdth
  gradv(2) = -dfvdth/ddisdth

  gradx(3) = dfxdh
  grady(3) = dfydh
  gradz(3) = dfzdh
  gradv(3) = dfvdh

end subroutine gradxyzv
!-----------------------------------------------------------------------
subroutine grapxyzv(alt,cth,sth, &
  dfxdln,dfydln,dfzdln,dfvdln,gradx,grady,gradz,gradv)
!
! Calculates east component of gradient near pole.
!
! Args:
  real,intent(in) :: alt,cth,sth
  real,intent(in) :: dfxdln,dfydln,dfzdln,dfvdln
  real,dimension(3),intent(out) :: gradx,grady,gradz,gradv
!
! Local:
  real :: d,d2,rho,dddthod,drhodth,dzetdth,ddisdth
!
!          40680925. = req**2 (rounded off)
!          272340.   = req**2 * E2, where E2 = (2. - 1./298.25)/298.25
!                      is square of eccentricity of ellipsoid.
!
  d2 = 40680925.e0 - 272340.e0*cth*cth
  d = sqrt(d2)
  rho = sth*(alt + 40680925.e0/d)
  dddthod = 272340.e0*cth*sth/d2
  drhodth = alt*cth + (40680925.e0/d)*(cth-sth*dddthod)
  dzetdth =-alt*sth - (40408585.e0/d)*(sth+cth*dddthod)
  ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)

  gradx(1) = dfxdln/rho
  grady(1) = dfydln/rho
  gradz(1) = dfzdln/rho
  gradv(1) = dfvdln/rho

end subroutine grapxyzv
!-----------------------------------------------------------------------
subroutine gradlpv(hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
  xlatm,xlonm,vmp,grclm,clmgrp,qdlat,rgrlp,b,clm,r3_2)
!
! Uses gradients of x,y,z,v to compute geomagnetic field and
! gradients of apex latitude, longitude.
!
! Args:
  real,intent(in) :: & ! scalar inputs
    hr,              & ! reference altitude (km)
    alt,             & ! altitude (km)
    fx,fy,fz,fv        ! interpolated values of x,y,z,v, plus 
                       ! pseudodipole component
  real,dimension(3),intent(in) :: & ! 3-component inputs
    gradx,grady,gradz,gradv ! interpolated gradients of x,y,z,v,
                            ! including pseudodipole components (east,north,up)
!
! Local:
  integer :: i
  real :: rr,r,rn,sqrror,cpm,spm,bo,rn2,x2py2,xnorm,xlp,slp,clp,grclp

  real,intent(out) :: & ! scalar outputs
    xlatm,  & !  modified apex latitude (lambda_m), degrees
    xlonm,  & !  apex longitude (phi_a), degrees
    vmp,    & !  magnetic potential, in T.m.
    qdlat,  & !  quasi-dipole latitude, degrees
    clm,    & !  cos(lambda_m)
    r3_2      !  ((re + alt)/(re + hr))**(3/2)

  real,dimension(3),intent(out) :: & ! 3-component outputs 
    grclm,   & ! grad(cos(lambda_m)), in km-1
    clmgrp,  & ! cos(lambda_m)*grad(phi_a), in km-1
    rgrlp,   & ! (re + alt)*grad(lambda')
    b          ! magnetic field, in nT

  xlatm=0. ; xlonm=0. ; vmp=0 ; grclm=0. ; clmgrp=0. ; rgrlp = 0. ; b=0.
  clm=0. ; r3_2=0. ; qdlat=0.

  rr = re + hr
  r  = re + alt
  rn = r/re
  sqrror = sqrt(rr/r)
  r3_2 = 1./sqrror/sqrror/sqrror
  xlonm = atan2(fy,fx)
  cpm = cos(xlonm)
  spm = sin(xlonm)
  xlonm = rtd*xlonm ! output
  bo = vp*1.e6 ! vp is module data; 1.e6 converts T to nT and km-1 to m-1
  rn2 = rn*rn
  vmp = vp*fv/rn2   ! output
  b(1) = -bo*gradv(1)/rn2
  b(2) = -bo*gradv(2)/rn2
  b(3) = -bo*(gradv(3)-2.*fv/r)/rn2

  x2py2 = fx*fx + fy*fy
  xnorm = sqrt(x2py2 + fz*fz)
  xlp = atan2(fz,sqrt(x2py2))
  slp = sin(xlp)
  clp = cos(xlp)
  qdlat = xlp*rtd   ! output
  clm = sqrror*clp  ! output
  if (clm > 1.) then
    write(6,"('>>> gradlpv: hr=',f12.3,' alt=',f12.3)") hr,alt
    write(6,"('    Point lies below field line that peaks at reference height.')")
    stop 'gradlpv'
  endif
  xlatm = rtd*acos(clm)
!
!  If southern magnetic hemisphere, reverse sign of xlatm
!
  if (slp < 0.) xlatm = -xlatm
  do i=1,3  
    grclp = cpm*gradx(i) + spm*grady(i)
    rgrlp(i) = r*(clp*gradz(i) - slp*grclp)
    grclm(i) = sqrror*grclp
    clmgrp(i) = sqrror*(cpm*grady(i)-spm*gradx(i))
  enddo
  grclm(3) = grclm(3) - sqrror*clp/(2.*r)

end subroutine gradlpv
!-----------------------------------------------------------------------
subroutine basevec(hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
                   bmag,sim,si,f,d,w,bhat,d1,d2,d3,e1,e2,e3,f1,f2)
!
! Computes base vectors and other parameters for apex coordinates.
! Vector components:  east, north, up
!
! Args:
  real,intent(in) :: & ! scalar inputs
    hr,      & ! reference altitude
    xlatm,   & ! modified apex latitude (deg)
    clm,     & ! cos(lambda_m)
    r3_2       ! ((re + altitude)/(re + hr))**(3/2)

  real,dimension(3),intent(in) :: & ! 3-component inputs
    grclm,   & ! grad(cos(lambda_m)), in km-1
    clmgrp,  & ! cos(lambda_m)*grad(phi_a), in km-1
    rgrlp,   & ! (re + altitude)*grad(lambda')
    b          ! ((re + altitude)/(re + hr))**(3/2)

  real,intent(out) :: & ! scalar output
    bmag,    & ! magnitude of magnetic field, in nT
    sim,     & ! sin(I_m) of Richmond reference
    si,      & ! sin(I)
    f,       & ! F of Richmond reference
    d,       & ! D of Richmond reference
    w          ! W of Richmond reference

  real,dimension(3),intent(out) :: & ! 3-component outputs
    bhat,             & ! unit vector along geomagnetic field direction
    d1,d2,d3,e1,e2,e3   ! base vectors of Richmond reference
  real,dimension(2),intent(out) :: & ! 2-component outputs
    f1(2),f2(2)         ! base vectors of Richmond reference
!
! Local:
  integer :: i
  real :: rr,simoslm,d1db,d2db

  rr = re + hr
  simoslm = 2./sqrt(4. - 3.*clm*clm)
  sim = simoslm*sin(xlatm*dtr)
  bmag = sqrt(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))
  d1db = 0.
  d2db = 0.
  do i=1,3
    bhat(i) = b(i)/bmag
    d1(i) = rr*clmgrp(i)
    d1db = d1db + d1(i)*bhat(i)
    d2(i) = rr*simoslm*grclm(i)
    d2db = d2db + d2(i)*bhat(i)
  enddo
!
! Ensure that d1,d2 are exactly perpendicular to B:
!
  do i=1,3
    d1(i) = d1(i) - d1db*bhat(i)
    d2(i) = d2(i) - d2db*bhat(i)
  enddo  
  e3(1) = d1(2)*d2(3) - d1(3)*d2(2)
  e3(2) = d1(3)*d2(1) - d1(1)*d2(3)
  e3(3) = d1(1)*d2(2) - d1(2)*d2(1) 
  d = bhat(1)*e3(1) + bhat(2)*e3(2) + bhat(3)*e3(3)
  do i=1,3
    d3(i) = bhat(i)/d
    e3(i) = bhat(i)*d ! Ensure that e3 lies along bhat.
  enddo
  e1(1) = d2(2)*d3(3) - d2(3)*d3(2)
  e1(2) = d2(3)*d3(1) - d2(1)*d3(3)
  e1(3) = d2(1)*d3(2) - d2(2)*d3(1)
  e2(1) = d3(2)*d1(3) - d3(3)*d1(2)
  e2(2) = d3(3)*d1(1) - d3(1)*d1(3)
  e2(3) = d3(1)*d1(2) - d3(2)*d1(1)
  w = rr*rr*clm*abs(sim)/(bmag*d)
  si = -bhat(3)
  f1(1) =  rgrlp(2)
  f1(2) = -rgrlp(1)
  f2(1) = -d1(2)*r3_2
  f2(2) =  d1(1)*r3_2
  f = f1(1)*f2(2) - f1(2)*f2(1)

end subroutine basevec
!-----------------------------------------------------------------------
subroutine dypol(colat,elon,vp)
!
! Output args:
  real,intent(out) :: &
    colat, & ! Geocentric colatitude of geomagnetic dipole north pole (deg)
    elon,  & ! East longitude of geomagnetic dipole north pole (deg)
    vp       ! Magnitude, in T.m, of dipole component of magnetic
             ! potential at geomagnetic pole and geocentric radius re
!
! Local:
  real :: gpl,ctp,stp
!
! Compute geographic colatitude and longitude of the north pole of
! earth centered dipole  
!
  gpl = sqrt( gb(2  )**2+ gb(3  )**2+ gb(4  )**2)
  ctp = gb(2  )/gpl
  stp = sqrt(1. - ctp*ctp)

  colat = (acos(ctp))*rtd
  elon = atan2( gb(4  ), gb(3  ))*rtd
!           
! Compute magnitude of magnetic potential at pole, radius Re.
!      .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through f0 in COFRM).
!
  vp = .2*gpl*re

end subroutine dypol
!-----------------------------------------------------------------------
subroutine apex_sub(date,dlat,dlon,alt,aht,alat,alon,bmag,xmag,ymag,zmag,vmp)
!
! Args:
  real,intent(in) :: date
  real,intent(inout) :: dlat,dlon,alt
  real,intent(out) :: aht,alat,alon,bmag,xmag,ymag,zmag,vmp
!
! Local:
  real :: clatp,polon,vpol,x,y,z,xre,yre,zre
  integer :: iflag

  call cofrm(date)
  call dypol(clatp,polon,vpol)
!
! colat,ctp,stp,elon,vp are in module data.
!
  colat = clatp
  ctp   = cos(clatp*dtr)
  stp   = sqrt(1.-ctp*ctp)

  elon  = polon
  vp    = vpol

  vmp = 0.
!
! Last 7 args of linapx are output:
!
  call linapx(dlat,dlon,alt, aht,alat,alon,xmag,ymag,zmag,bmag)

  xmag = xmag*1.e5
  ymag = ymag*1.e5
  zmag = zmag*1.e5
  bmag = bmag*1.e5
  call gd2cart (dlat,dlon,alt,x,y,z)
  iflag = 3
  xre = x/re ; yre = y/re ; zre = z/re
  call feldg(iflag,xre,yre,zre,bx,by,bz,vmp)

end subroutine apex_sub
!-----------------------------------------------------------------------
subroutine linapx(gdlat,glon,alt,aht,alat,alon,xmag,ymag,zmag,fmag)
!
! Input Args:
!
  real,intent(inout) :: & ! These may be changed by convrt, depending on iflag
    gdlat,              & ! latitude of starting point (deg)
    glon,               & ! longitude of starting point (deg)
    alt                   ! height of starting point (km)
!
! Output Args:
!
  real,intent(out) :: &
    aht,              & ! (Apex height+req)/req, where req is equatorial earth radius
    alat,             & ! Apex latitude (deg)
    alon,             & ! Apex longitude (deg)
    xmag,             & ! North component of magnetic field at starting point
    ymag,             & ! East component of magnetic field at starting point
    zmag,             & ! Down component of magnetic field at starting point
    fmag                ! Magnetic field magnitude at starting point
!
! Local:
!
  real :: gclat,r,singml,cgml2,rho,xlat,xlon,ht
  real :: bnrth,beast,bdown,babs,y1,y2,y3
  integer :: iflag,iapx
  integer,parameter :: maxs = 200
!
! Determine step size as a function of geomagnetic dipole
! coordinates of the starting point
!
  iflag = 2 ! gclat,r are returned
  call convrt(iflag,gdlat,alt,gclat,r,'linapx')

  singml = ctp*sin(gclat*dtr) + stp*cos(gclat*dtr)*cos((glon-elon)*dtr)
  cgml2 = max(0.25,1.-singml*singml)
  ds = .06*r/cgml2 - 370. ! ds is in module data

  yapx = 0. ! init (module data)
!
! Convert from geodetic to earth centered cartesian coordinates:
!
  call gd2cart(gdlat,glon,alt,y(1),y(2),y(3))
  nstp = 0
!
! Get magnetic field components to determine the direction for tracing field line:
!
  iflag = 1 
  call feldg(iflag,gdlat,glon,alt,xmag,ymag,zmag,fmag)

  sgn = sign(1.,-zmag)
!
! Use cartesian coordinates to get magnetic field components
! (from which gradients steer the tracing)
!
100 continue
  iflag = 2 ! module data bx,by,bz,bb are returned
  y1 = y(1)/re ; y2 = y(2)/re ; y3 = y(3)/re
  call feldg(iflag,y1,y2,y3,bx,by,bz,bb)
  nstp = nstp + 1
!
! Quit if too many steps.
!
  if (nstp >= maxs) then
    rho = sqrt(y(1)*y(1) + y(2)*y(2))
    iflag = 3 ! xlat and ht are returned
    call convrt(iflag,xlat,ht,rho,y(3),'linapx')
    xlon = rtd*atan2(y(2),y(1))
    iflag = 1
    call feldg(iflag,xlat,xlon,ht,bnrth,beast,bdown,babs)
    call dipapx(xlat,xlon,ht,bnrth,beast,bdown,aht,alon)
    alat = -sgn*rtd*acos(sqrt(1./aht))
    return
  endif
!
! Find next point using adams algorithm after 7 points
!
  call itrace(iapx)
  if (iapx == 1) goto 100
!
! Maximum radius just passed.  Find apex coords
!
  call fndapx(alt,zmag,aht,alat,alon)

end subroutine linapx
!-----------------------------------------------------------------------
subroutine convrt(iflag,gdlat,alt,x1,x2,caller)
!
! Convert space point from geodetic to geocentric or vice versa.
!
! iflag = 1: Convert from geodetic to cylindrical
!   Input:  gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!   Output: x1    = Distance from Earth's rotation axis (km)
!           x2    = Distance above (north of) Earth's equatorial plane (km)
!
! iflag = 2: Convert from geodetic to geocentric spherical
!   Input:  gdlat = Geodetic latitude (deg) 
!           alt   = Altitude above reference ellipsoid (km)
!   Output: x1    = Geocentric latitude (deg)
!           x2    = Geocentric distance (km)
!
! iflag = 3: Convert from cylindrical to geodetic
!   Input:  x1    = Distance from Earth's rotation axis (km)
!           x2    = Distance from Earth's equatorial plane (km)
!   Output: gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!
! iflag = 4: Convert from geocentric spherical to geodetic
!   Input:  x1    = Geocentric latitude (deg)
!           x2    = Geocentric distance (km)
!   Output: gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!
! Args:
  integer,intent(in) :: iflag
  real,intent(inout) :: gdlat,alt
  real,intent(inout) :: x1,x2
  character(len=*),intent(in) :: caller
!
! Local:
  real :: sinlat,coslat,d,z,rho,rkm,scl,gclat,ri,a2,a4,a6,a8,&
    ccl,s2cl,c2cl,s4cl,c4cl,s8cl,s6cl,dltcl,sgl
  real,parameter ::                                         &
    fltnvrs = 298.25                                      , &
    e2=(2.-1./fltnvrs)/fltnvrs                            , &
    e4=e2*e2, e6=e4*e2, e8=e4*e4                          , &
    ome2req = (1.-e2)*req                                 , &
    A21 =     (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024. , &
    A22 =     (                        E6 +     E8)/  32. , &
    A23 = -3.*(                     4.*E6 +  3.*E8)/ 256. , &
    A41 =    -(           64.*E4 + 48.*E6 + 35.*E8)/1024. , &
    A42 =     (            4.*E4 +  2.*E6 +     E8)/  16. , &
    A43 =                                   15.*E8 / 256. , &
    A44 =                                      -E8 /  16. , &
    A61 =  3.*(                     4.*E6 +  5.*E8)/1024. , &
    A62 = -3.*(                        E6 +     E8)/  32. , &
    A63 = 35.*(                     4.*E6 +  3.*E8)/ 768. , &
    A81 =                                   -5.*E8 /2048. , &
    A82 =                                   64.*E8 /2048. , &
    A83 =                                 -252.*E8 /2048. , &
    A84 =                                  320.*E8 /2048.

  if (iflag < 3) then ! geodetic to geocentric
!
! Compute rho,z
    sinlat = sin(gdlat*dtr)
    coslat = sqrt(1.-sinlat*sinlat)
    d      = sqrt(1.-e2*sinlat*sinlat)
    z      = (alt+ome2req/d)*sinlat
    rho    = (alt+req/d)*coslat
    x1 = rho
    x2 = z
    if (iflag == 1) return
!
! Compute gclat,rkm
    rkm   = sqrt(z*z+rho*rho)
    gclat = rtd*atan2(z,rho)
    x1 = gclat
    x2 = rkm
    return    ! iflag == 2
  endif ! iflag < 3
!
! Geocentric to geodetic:
  if (iflag == 3) then
    rho = x1
    z   = x2
    rkm = sqrt(z*z+rho*rho)
    scl = z/rkm
    gclat = asin(scl)*rtd
  elseif (iflag == 4) then
    gclat = x1
    rkm = x2
    scl = sin(gclat*dtr)
  else
    return
  endif
!
! iflag == 3 or 4:
!
  ri = req/rkm
  a2 = ri*(a21+ri*(a22+ri* a23))
  a4 = ri*(a41+ri*(a42+ri*(a43+ri*a44)))
  a6 = ri*(a61+ri*(a62+ri* a63))
  a8 = ri*(a81+ri*(a82+ri*(a83+ri*a84)))
  ccl = sqrt(1.-scl*scl)
  s2cl = 2.*scl*ccL
  c2cl = 2.*ccl*ccl-1.
  s4cl = 2.*s2cl*c2cl
  c4cl = 2.*c2cl*c2cl-1.
  s8cl = 2.*s4cl*c4cl
  s6cl = s2cl*c4cl+c2cl*s4cl
  dltcl = s2cl*a2+s4cl*a4+s6cl*a6+s8cl*a8
  gdlat = dltcl*rtd+gclat
  sgl = sin(gdlat*dtr)
  alt = rkm*cos(dltcl)-req*sqrt(1.-e2*sgl*sgl)

end subroutine convrt
!-----------------------------------------------------------------------
subroutine gd2cart(gdlat,glon,alt,x,y,z)
!
! Arg:
  real,intent(inout) :: gdlat,alt,z
  real,intent(in) :: glon
  real,intent(out) :: x,y
!
! Local:
  real :: ang,rho
  integer :: iflag

  iflag = 1 ! Convert from geodetic to cylindrical (rho,z are output)
  call convrt(iflag,gdlat,alt,rho,z,'gd2cart')

  ang = glon*dtr
  x = rho*cos(ang)
  y = rho*sin(ang)

end subroutine gd2cart
!-----------------------------------------------------------------------
subroutine feldg(iflag,glat,glon,alt,bnrth,beast,bdown,babs)
!
! Compute the DGRF/IGRF field components at the point glat,glon,alt.
! cofrm must be called to establish coefficients for correct date
! prior to calling FELDG.
!
! iflag = 1:
!   Inputs:
!     glat = Latitude of point (deg)
!     glon = Longitude of point (deg)
!     alt  = Height of point (km)
!   Outputs:
!     bnrth = North component of field vector (Gauss)
!     beast = East component of field vector (Gauss)
!     bdown = Downward component of field vector (Gauss)
!     babs  = Magnitude of field vector (Gauss)  
!
! iflag = 2:
!   Inputs:
!     glat = x coordinate (in units of earth radii 6371.2 km)
!     glon = y coordinate (in units of earth radii 6371.2 km)
!     alt  = z coordinate (in units of earth radii 6371.2 km)
!   Outputs:
!     bnrth = x component of field vector (Gauss)
!     beast = y component of field vector (Gauss)
!     bdown = z component of field vector (Gauss)
!     babs  = Magnitude of field vector (Gauss)
!
! iflag = 3:
!   Inputs:
!     glat = x coordinate (in units of earth radii 6371.2 km)
!     glon = y coordinate (in units of earth radii 6371.2 km)
!     alt  = z coordinate (in units of earth radii 6371.2 km)
!   Outputs:
!     bnrth = Dummy variable
!     beast = Dummy variable
!     babs  = Legacy code had "Dummy variable" here, but its
!             set at the end if iflag==3.
!
! Args:
  integer,intent(in) :: iflag
  real,intent(in)    :: glon
  real,intent(inout) :: glat,alt
  real,intent(out)   :: bnrth,beast,bdown,babs
!
! Local:
  integer :: i,is,ihmax,last,imax,mk,k,ih,m,il,ihm,ilm
  real :: rlat,ct,st,rlon,cp,sp,xxx,yyy,zzz,rq,f,x,y,z
  real :: xi(3),h(ncoef),g(ncoef)
  real :: s,t,bxxx,byyy,bzzz,brho

  if (iflag == 1) then
    is   = 1
    rlat = glat*dtr
    ct   = sin(rlat)
    st   = cos(rlat)
    rlon = glon*dtr
    cp   = cos(rlon)
    sp   = sin(rlon)
    call gd2cart(glat,glon,alt,xxx,yyy,zzz)
    xxx = xxx/re
    yyy = yyy/re
    zzz = zzz/re
  else
    is  = 2
    xxx = glat
    yyy = glon
    zzz = alt
  endif
  rq    = 1./(xxx**2+yyy**2+zzz**2)
  xi(1) = xxx*rq
  xi(2) = yyy*rq
  xi(3) = zzz*rq
  ihmax = nmax*nmax+1
  last  = ihmax+nmax+nmax
  imax  = nmax+nmax-1
!
! Legacy code checks here to see if iflag or last call to cofrm have changed.
! For now, just do it anyway:
!
  if (iflag /= 3) then
    do i=1,last
      g(i) = gb(i) ! gb is module data from cofrm
    enddo
  else
    do i=1,last
      g(i) = gv(i) ! gv is module data from cofrm
    enddo
  endif

  do i=ihmax,last
    h(i) = g(i)
  enddo

  mk = 3
  if (imax == 1) mk = 1

  do k=1,mk,2
    i = imax
    ih = ihmax

100 continue
    il = ih-i
    f = 2./float(i-k+2)
    x = xi(1)*f
    y = xi(2)*f
    z = xi(3)*(f+f)

    i = i-2
    if (i < 1) then
      h(il) = g(il) + z*h(ih) + 2.*(x*h(ih+1)+y*h(ih+2))
    elseif (i == 1) then
      h(il+2) = g(il+2) + z*h(ih+2) + x*h(ih+4) - y*(h(ih+3)+h(ih))
      h(il+1) = g(il+1) + z*h(ih+1) + y*h(ih+4) + x*(h(ih+3)-h(ih))
      h(il)   = g(il)   + z*h(ih)   + 2.*(x*h(ih+1)+y*h(ih+2))
    else
      do m=3,i,2
        ihm = ih+m
        ilm = il+m
        h(ilm+1) = g(ilm+1)+ z*h(ihm+1) + x*(h(ihm+3)-h(ihm-1))- &
                   y*(h(ihm+2)+h(ihm-2))
        h(ilm)   = g(ilm)  + z*h(ihm)   + x*(h(ihm+2)-h(ihm-2))+ &
                   y*(h(ihm+3)+h(ihm-1))
      enddo
      h(il+2) = g(il+2) + z*h(ih+2) + x*h(ih+4) - y*(h(ih+3)+h(ih))
      h(il+1) = g(il+1) + z*h(ih+1) + y*h(ih+4) + x*(h(ih+3)-h(ih))
      h(il)   = g(il)   + z*h(ih)   + 2.*(x*h(ih+1)+y*h(ih+2))
    endif

    ih = il
    if (i >= k) goto 100
  enddo ! k=1,mk,2

  s = .5*h(1)+2.*(h(2)*xi(3)+h(3)*xi(1)+h(4)*xi(2))
  t = (rq+rq)*sqrt(rq)
  bxxx = t*(h(3)-s*xxx)
  byyy = t*(h(4)-s*yyy)
  bzzz = t*(h(2)-s*zzz)
  babs = sqrt(bxxx**2+byyy**2+bzzz**2)
  if (is .eq. 1) then            ! (convert back to geodetic)
    beast = byyy*cp-bxxx*sp
    brho  = byyy*sp+bxxx*cp
    bnrth = bzzz*st-brho*ct
    bdown = -bzzz*ct-brho*st
  elseif (is .eq. 2) then        ! (leave in earth centered cartesian)
    bnrth = bxxx
    beast = byyy
    bdown = bzzz
  endif
!
! Magnetic potential computation makes use of the fact that the
! calculation of V is identical to that for r*Br, if coefficients
! in the latter calculation have been divided by (n+1) (coefficients
! GV).  Factor .1 converts km to m and gauss to tesla.
!
  if (iflag == 3) babs = (bxxx*xxx + byyy*yyy + bzzz*zzz)*re*.1

end subroutine feldg
!-----------------------------------------------------------------------
subroutine dipapx(gdlat,gdlon,alt,bnorth,beast,bdown,a,alon)
!
! Compute a, alon from local magnetic field using dipole and spherical approx.
! Reference from legacy code: 940501 A. D. Richmond
!
! Input:
!   gdlat  = geodetic latitude, degrees
!   gdlon  = geodetic longitude, degrees
!   alt    = altitude, km
!   bnorth = geodetic northward magnetic field component (any units)
!   beast  = eastward magnetic field component
!   bdown  = geodetic downward magnetic field component
! Output:
!   a      = apex radius, 1 + h_A/R_eq
!   alon   = apex longitude, degrees
!     
! Algorithm: 
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let G be point at GDLAT,GDLON.
!   Let E be point on sphere below apex of dipolar field line passing through G.
!   Let TD be dipole colatitude of point G, found by applying dipole formula
!     for dip angle to actual dip angle.
!   Let B be Pi plus local declination angle.  B is in the direction
!     from G to E.
!   Let TG be colatitude of G.
!   Let ANG be longitude angle from GM to G.
!   Let TE be colatitude of E.
!   Let TP be colatitude of GM.
!   Let A be longitude angle from G to E.
!   Let APANG = A + ANG
!   Let PA be geomagnetic longitude, i.e., Pi minus angle measured
!     counterclockwise from arc GM-E to arc GM-GP.
!   Let TF be arc length between GM and E.
!   Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.
!
  real,intent(in)  :: gdlat,gdlon,alt,bnorth,beast,bdown
  real,intent(out) :: a,alon
!
! Local:
  real :: bhor,std,ctd,sb,cb,ctg,stg,ang,sang,cang,cte,ste,sa,ca, &
          cottd,capang,sapang,stfcpa,stfspa,ha,r

  bhor = sqrt(bnorth*bnorth + beast*beast)
  if (bhor == 0.) then
    alon = 0.
    a = 1.e34
    return
  endif

  cottd = bdown*.5/bhor
  std = 1./sqrt(1.+cottd*cottd)
  ctd = cottd*std
  sb = -beast/bhor
  cb = -bnorth/bhor
  ctg = sin(gdlat*dtr)
  stg = cos(gdlat*dtr)
  ang = (gdlon-elon)*dtr
  sang = sin(ang)
  cang = cos(ang)
  cte = ctg*std + stg*ctd*cb
  ste = sqrt(1. - cte*cte)
  sa = sb*ctd/ste
  ca = (std*stg - ctd*ctg*cb)/ste
  capang = ca*cang - sa*sang
  sapang = ca*sang + sa*cang
  stfcpa = ste*ctp*capang - cte*stp
  stfspa = sapang*ste
  alon = atan2(stfspa,stfcpa)*rtd
  r = alt + re
  ha = alt + r*cottd*cottd
  a = 1. + ha/req

end subroutine dipapx
!-----------------------------------------------------------------------
subroutine itrace(iapx)
  save
!
! Uses 4-point ADAMS formula after initialization.
! First 7 iterations advance point by 3 steps.
!
! y(3), yp(3), yapx(3,3), sgn and nstp are in module data
! yploc(3,4) is local
!
! Arg:
  integer,intent(out) :: iapx
!
! Local:
  integer :: i,j
  real :: yploc(3,4) ! local yp (i.e., not module data yp)
  real :: term,d2,d6,d12,d24,rc,rp

  iapx = 1
!
! Field line is defined by the following differential equations
! in cartesian coordinates.
! (yapx,yp,y are module data)
!
  yploc(1,4) = sgn*bx/bb 
  yploc(2,4) = sgn*by/bb 
  yploc(3,4) = sgn*bz/bb 

  if (nstp > 7) then
    do i=1,3
      yapx(i,1) = yapx(i,2)
      yapx(i,2) = y(i)
      yp(i) = y(i)
      term = 55.*yploc(i,4)-59.*yploc(i,3)+37.*yploc(i,2)-9.*yploc(i,1)
      y(i) = yp(i) + d24*term
      yapx(i,3) = y(i)
      do j=1,3
        yploc(i,j) = yploc(i,j+1)
      enddo
    enddo
    rc = rdus ( y(1),  y(2),  y(3))
    rp = rdus (yp(1), yp(2), yp(3))
    if (rc < rp) iapx=2
    return
  endif

  do i=1,3
    select case (nstp)
      case (1)
        d2  = ds/2.
        d6  = ds/6.
        d12 = ds/12.
        d24 = ds/24.
        yploc(i,1)= yploc(i,4)
        yp(i)     = y(i)
        yapx(i,1) = y(i)
        y(i) = yp(i) + ds*yploc(i,1)
      case (2)
        yploc(i,2) = yploc(i,4)
        y(i) = yp(i) + d2*(yploc(i,2)+yploc(i,1))
      case (3)
        y(i) = yp(i) + d6*(2.*yploc(i,4)+yploc(i,2)+3.*yploc(i,1))
      case (4)
        yploc(i,2) = yploc(i,4)
        yapx(i,2)  = y(i)
        yp(i)      = y(i)
        y(i)       = yp(i) + d2*(3.*yploc(i,2)-yploc(i,1))
      case (5)
        y(i) = yp(i) + d12*(5.*yploc(i,4)+8.*yploc(i,2)-yploc(i,1))
      case (6)
        yploc(i,3) = yploc(i,4)
        yp(i)      = y(i)
        yapx(i,3)  = y(i)
        y(i)       = yp(i) + d12*(23.*yploc(i,3)-16.*yploc(i,2)+5.*yploc(i,1))
      case (7)
        yapx(i,1) = yapx(i,2)
        yapx(i,2) = yapx(i,3)
        y(i) = yp(i) + d24*(9.*yploc(i,4)+19.*yploc(i,3)-5.*yploc(i,2)+yploc(i,1))
        yapx(i,3) = y(i)
      case default
        write(6,"('>>> itrace: unresolved case nstp=',i4)") nstp
        stop 'itrace'
    end select
  enddo
!
! Signal if apex passed:
!
  if (nstp == 6 .or. nstp == 7) then
    rc = rdus( yapx(1,3), yapx(2,3), yapx(3,3))
    rp = rdus( yapx(1,2), yapx(2,2), yapx(3,2))
    if (rc < rp) iapx=2
  endif

end subroutine itrace
!-----------------------------------------------------------------------
real function rdus(d,e,f)
  real,intent(in) :: d,e,f
  rdus = sqrt(d**2 + e**2 + f**2)
end function rdus
!-----------------------------------------------------------------------
subroutine fndapx(alt,zmag,a,alat,alon)
!
! Find apex coords once tracing has signalled that the apex has been passed.
!
! Args:
  real,intent(in) :: alt,zmag
  real,intent(out) :: a,alat,alon
!
! Local:
  integer :: i,iflag_convrt, iflag_feldg
  real :: z(3),ht(3),yloc(3),gdlt,gdln,x,ydum,f,rho,xinter,rasq,xlon,ang,&
    cang,sang,r,cte,ste,stfcpa,stfspa
!
! Get geodetic field components.
!
  iflag_feldg = 1
  iflag_convrt = 3
  do i=1,3
    rho = sqrt(yapx(1,i)**2+yapx(2,i)**2)
    call convrt(iflag_convrt,gdlt,ht(i),rho,yapx(3,i),'fndapx')
    gdln = rtd*atan2(yapx(2,i),yapx(1,i))
    call feldg(iflag_feldg,gdlt,gdln,ht(i),x,ydum,z(i),f)
  enddo 
!
! Find cartesian coordinates at dip equator by interpolation
!
  do i=1,3
    call fint(z(1),z(2),z(3),yapx(i,1),yapx(i,2),yapx(i,3),0.,yloc(i))
  enddo
!
! Find apex height by interpolation
!
  call fint(z(1),z(2),z(3),ht(1),ht(2),ht(3),0.,xinter)
!
! Ensure that XINTER is not less than original starting altitude:
  xinter = max(alt,xinter)
  a = (req+xinter)/req
!
! Find apex coordinates , giving alat sign of dip at starting point.  
! Alon is the value of the geomagnetic longitude at the apex.
!
  if (a < 1.) then
    write(6,"('>>> fndapx: a=',e12.4,' < 1.')") a
    stop 'fndapx'
  endif

  rasq = rtd*acos(sqrt(1./a))
  alat = sign(rasq,zmag)
!
! Algorithm for ALON:
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let XLON be longitude of apex.
!   Let TE be colatitude of apex.
!   Let ANG be longitude angle from GM to apex.
!   Let TP be colatitude of GM.
!   Let TF be arc length between GM and apex.
!   Let PA = ALON be geomagnetic longitude, i.e., Pi minus angle measured 
!     counterclockwise from arc GM-apex to arc GM-GP.
!   Then, using notation C=cos, S=sin, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively. 
!
  xlon = atan2(yloc(2),yloc(1))
  ang  = xlon-elon*dtr
  cang = cos(ang)
  sang = sin(ang)
  r    = sqrt(yloc(1)**2+yloc(2)**2+yloc(3)**2)
  cte  = yloc(3)/r
  ste  = sqrt(1.-cte*cte)
  stfcpa = ste*ctp*cang - cte*stp
  stfspa = sang*ste
  alon = atan2(stfspa,stfcpa)*rtd

end subroutine fndapx
!-----------------------------------------------------------------------
subroutine fint(a1,a2,a3,a4,a5,a6,a7,result)
!
! Second degree interpolation
!
! Args:
  real,intent(in) :: a1,a2,a3,a4,a5,a6,a7
  real,intent(out) :: result

  result = ((a2-a3)*(a7-a2)*(a7-a3)*a4-(a1-a3)*(a7-a1)*(a7-a3)*a5+ &
    (a1-a2)*(a7-a1)*(a7-a2)*a6)/((a1-a2)*(a1-a3)*(a2-a3))

end subroutine fint
!-----------------------------------------------------------------------
subroutine gm2gc(gmlat,gmlon,gclat,gclon)
!
! Args:
  real,intent(in)  :: gmlat,gmlon
  real,intent(out) :: gclat,gclon
!
! Local:
  real :: ylat,ylon,stm,ctm,ctc

  stm = cos(gmlat*dtr)
  ctm = sin(gmlat*dtr)
  ctc = ctp*ctm - stp*stm*cos(gmlon*dtr) ! ctp,stp are module data
  ctc = min(ctc,1.)
  ctc = max(ctc,-1.)
  gclat = asin(ctc)*rtd
  gclon = atan2(stp*stm*sin(gmlon*dtr),ctm-ctp*ctc)
!
! elon is in module data, and was set by dypol (called from apex_mka)
!
  gclon = gclon*rtd + elon 
  if (gclon < -180.) gclon = gclon + 360.

end subroutine gm2gc
!-----------------------------------------------------------------------
subroutine intrp(glat,glon,alt, gplat,gplon,gpalt, nlat,nlon,nalt, &
                 fx,fy,fz,fv,                                      &
                 dfxdth,dfydth,dfzdth,dfvdth,                      &
                 dfxdln,dfydln,dfzdln,dfvdln,                      &
                 dfxdh ,dfydh ,dfzdh ,dfvdh, ier)
!
! Args:
!
  real,intent(in)    :: glat,glon,alt
  integer,intent(in) :: nlat,nlon,nalt
  real,intent(in)    :: gplat(nlat),gplon(nlon),gpalt(nalt)
  real,intent(out)   ::          &
    fx,fy,fz,fv,                 &
    dfxdth,dfydth,dfzdth,dfvdth, &
    dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh ,dfydh ,dfzdh ,dfvdh
  integer,intent(out) :: ier
!
! Local:
!
  integer :: i,j,k,i0,j0,k0
  real :: glonloc,xi,dlon,yj,dlat,hti,diht,zk,fac,omfac
  real :: dfxdn,dfxde,dfxdd, &
          dfydn,dfyde,dfydd, &
          dfzdn,dfzde,dfzdd, &
          dfvdn,dfvde,dfvdd, &
          dmf,dmdfdn,dmdfde,dmdfdd

  ier = 0
  glonloc = glon
  if (glonloc < gplon(1))    glonloc = glonloc + 360.
  if (glonloc > gplon(nlon)) glonloc = glonloc - 360.
!
  i0 = 0
  do i=1,nlat-1
    if (glat >= gplat(i).and.glat <= gplat(i+1)) then
      i0 = i
      dlat = gplat(i+1)-gplat(i)
      xi = (glat - gplat(i)) / dlat
      exit 
    endif
  enddo
  if (i0==0) then
    write(6,"('>>> intrp: could not bracket glat=',f9.3,' in gplat=',/,(6f9.2))") &
      glat,gplat
    ier = 1
    return 
  endif

  j0 = 0
  do j=1,nlon-1
    if (glon >= gplon(j).and.glon <= gplon(j+1)) then
      j0 = j
      dlon = gplon(j+1)-gplon(j)
      yj = (glon - gplon(j)) / dlon
      exit 
    endif
  enddo
  if (j0==0) then
    write(6,"('>>> intrp: could not bracket glon=',f9.3,' in gplon=',/,(6f9.2))") &
      glon,gplon
    ier = 1
    return 
  endif

  k0 = 0
  do k=1,nalt-1
    if (alt >= gpalt(k).and.alt <= gpalt(k+1)) then
      k0 = k
      hti = re/(re+alt)
      diht = re/(re+gpalt(k+1)) - re/(re+gpalt(k))
      zk = (hti - re/(re+gpalt(k))) / diht
      exit 
    endif
  enddo
  if (k0==0) then
    write(6,"('>>> intrp: could not bracket alt=',f12.3,' in gpalt=',/,(6f12.2))") &
      alt,gpalt
    ier = 1
    return 
  endif

  call trilin(xarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fx,dfxdn,dfxde,dfxdd)
  dfxdth = -dfxdn*rtd/dlat
  dfxdln =  dfxde*rtd/dlon
  dfxdh  = -hti*hti*dfxdd/(re*diht)

  call trilin(yarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fy,dfydn,dfyde,dfydd)
  dfydth = -dfydn*rtd/dlat
  dfydln =  dfyde*rtd/dlon
  dfydh  = -hti*hti*dfydd/(re*diht)

  call trilin(zarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fz,dfzdn,dfzde,dfzdd)
  dfzdth = -dfzdn*rtd/dlat
  dfzdln =  dfzde*rtd/dlon
  dfzdh  = -hti*hti*dfzdd/(re*diht)

  call trilin(varray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fv,dfvdn,dfvde,dfvdd)
  dfvdth = -dfvdn*rtd/dlat
  dfvdln =  dfvde*rtd/dlon
  dfvdh  = -hti*hti*dfvdd/(re*diht)

  if (nlat < 3) return
!
! Improve calculation of longitudinal derivatives near poles
!
  if (glat < dlat-90.) then
    fac = .5*xi
    omfac = 1. - fac
    xi = xi - 1.
    i0 = i0 + 1
    call trilin (xarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfxdln = dfxdln*omfac + fac*dmdfde*rtd/dlon
    call trilin (yarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfydln = dfydln*omfac + fac*dmdfde*rtd/dlon
    call trilin (varray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfvdln = dfvdln*omfac + fac*dmdfde*rtd/dlon
  endif

  if (glat > 90.-dlat) then
    fac = .5*(1.-xi)
    omfac = 1. - fac
    xi = xi + 1.
    i0 = i0 - 1
    call trilin (xarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfxdln = dfxdln*omfac + fac*dmdfde*rtd/dlon
    call trilin (yarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfydln = dfydln*omfac + fac*dmdfde*rtd/dlon
    call trilin (varray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfvdln = dfvdln*omfac + fac*dmdfde*rtd/dlon
  endif

end subroutine intrp
!-----------------------------------------------------------------------
subroutine trilin(u,nlat,nlon,xi,yj,zk,fu,dfudx,dfudy,dfudz)
!
! Args:
  integer,intent(in) :: &
    nlat,               & ! first dimension of u from calling routine
    nlon                  ! second dimension of u from calling routine
  real,intent(in)    :: &
    u(1:2,1:2,1:2),     & ! u(1,1,1) is address of lower corner of interpolation box
    xi,  & ! fractional distance across box in x direction
    yj,  & ! fractional distance across box in y direction
    zk     ! fractional distance across box in z direction
  real,intent(out)   :: &
    fu,                 & ! interpolated value of u
    dfudx,              & ! interpolated derivative of u with respect to i (x direction)
    dfudy,              & ! interpolated derivative of u with respect to j (y direction)
    dfudz                 ! interpolated derivative of u with respect to k (z direction)
!
! Local:
  real :: omxi,omyj,omzk

! write(6,"('Enter trilin: xi,yj,zk=',3e12.4)") xi,yj,zk
! write(6,"('Enter trilin: u(1,1,1),u(1,2,1),u(1,1,2),u(1,2,2)=',4e12.4)") &
!                          u(1,1,1),u(1,2,1),u(1,1,2),u(1,2,2)
! write(6,"('Enter trilin: u(2,1,1),u(2,2,1),u(2,1,2),u(2,2,2)=',4e12.4)") &
!                          u(2,1,1),u(2,2,1),u(2,1,2),u(2,2,2)

  omxi = 1. - xi
  omyj = 1. - yj
  omzk = 1. - zk

  fu = u(1,1,1)*omxi*omyj*omzk &
     + u(2,1,1)*xi*omyj*omzk   &
     + u(1,2,1)*omxi*yj*omzk   &
     + u(1,1,2)*omxi*omyj*zk   &
     + u(2,2,1)*xi*yj*omzk     &
     + u(2,1,2)*xi*omyj*zk     &
     + u(1,2,2)*omxi*yj*zk     &
     + u(2,2,2)*xi*yj*zk

  dfudx = (u(2,1,1)-u(1,1,1))*omyj*omzk &
        + (u(2,2,1)-u(1,2,1))*yj*omzk   &
        + (u(2,1,2)-u(1,1,2))*omyj*zk   &
        + (u(2,2,2)-u(1,2,2))*yj*zk
  dfudy = (u(1,2,1)-u(1,1,1))*omxi*omzk &
        + (u(2,2,1)-u(2,1,1))*xi*omzk   &
        + (u(1,2,2)-u(1,1,2))*omxi*zk   &
        + (u(2,2,2)-u(2,1,2))*xi*zk
  dfudz = (u(1,1,2)-u(1,1,1))*omxi*omyj &
        + (u(2,1,2)-u(2,1,1))*xi*omyj   &
        + (u(1,2,2)-u(1,2,1))*omxi*yj   &
        + (u(2,2,2)-u(2,2,1))*xi*yj

end subroutine trilin
!-----------------------------------------------------------------------
subroutine adpl(glat,glon,cth,sth,fx,fy,fz,fv, &
                dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)
!
!  Add-back of pseudodipole component to x,y,z,v and their derivatives.
!
! Args:
  real,intent(in)      :: glat,glon
  real,intent(out)     :: cth,sth
  real,intent(inout)   ::        &
     fx,fy,fz,fv,                &
    dfxdth,dfydth,dfzdth,dfvdth, &
    dfxdln,dfydln,dfzdln,dfvdln
!
! Local:
  real :: cph,sph,ctm

  cph = cos((glon-elon)*dtr)
  sph = sin((glon-elon)*dtr)
  cth = sin(glat*dtr)
  sth = cos(glat*dtr)
  ctm = ctp*cth + stp*sth*cph
  fx = fx + sth*ctp*cph - cth*stp
  fy = fy + sth*sph
  fz = fz + ctm
  fv = fv - ctm

  dfxdth = dfxdth + ctp*cth*cph + stp*sth
  dfydth = dfydth + cth*sph 
  dfzdth = dfzdth - ctp*sth + stp*cth*cph
  dfvdth = dfvdth + ctp*sth - stp*cth*cph

  dfxdln = dfxdln - ctp*sth*sph
  dfydln = dfydln + sth*cph
  dfzdln = dfzdln - stp*sth*sph
  dfvdln = dfvdln + stp*sth*sph

end subroutine adpl
!-----------------------------------------------------------------------
subroutine setmiss(xmiss,xlatm,alon,vmp,b,bmag,be3,sim,si,f,d,w, &
  bhat,d1,d2,d3,e1,e2,e3,f1,f2)
!
! Args:
  real,intent(in)  :: xmiss
  real,intent(out) :: xlatm,alon,vmp,bmag,be3,sim,si,f,d,w
  real,dimension(3),intent(out) :: bhat,d1,d2,d3,e1,e2,e3,b
  real,dimension(2),intent(out) :: f1,f2

  xlatm = xmiss
  alon  = xmiss
  vmp   = xmiss
  bmag  = xmiss
  be3   = xmiss
  sim   = xmiss
  si    = xmiss
  f     = xmiss
  d     = xmiss
  w     = xmiss
  bhat  = xmiss
  d1    = xmiss
  d2    = xmiss
  d3    = xmiss
  e1    = xmiss
  e2    = xmiss
  e3    = xmiss
  b     = xmiss
  f1    = xmiss
  f2    = xmiss

end subroutine setmiss
!-----------------------------------------------------------------------
subroutine cofrm(date)
  implicit none
!
! Input arg:
  real,intent(in) :: date
!
! Local:
  integer :: m,n,i,l,ll,lm,nmx,nc,kmx,k,mm,nn
  real :: t,one,tc,r,f,f0
  integer,parameter :: n1=120, n2=195, isv=0
  integer,parameter :: &
    ncn1=19, & ! number of coefficients dimensioned n1
    ncn2=6     ! number of coefficients dimensioned n2: increase with each IGRF update
  integer,parameter :: ngh = n1*ncn1 + n2*ncn2 + 1 ! not sure why the extra +1
  real,save :: g1(n1,ncn1), g2(n2,ncn2), gh(ngh)
  real,parameter :: alt = 0.

  if (date < 1900. .or. date > 2025.) then
    write(6,"('>>> cofrm: date=',f8.2,' Date must be >= 1900 and <= 2025')") date
    stop 'cofrm'
  endif
  if (date > 2020.) then
    write(6,"('>>> WARNING cofrm:')")
    write(6,"(/,'   This version of IGRF is intended for use up to ')")
    write(6,"('     2020. Values for ',f9.3,' will be computed but')") date
    write(6,"('     may be of reduced accuracy.',/)")
  endif

  g1(:,1) = (/                                               & ! 1900
   -31543.,-2298., 5922., -677., 2905.,-1061.,  924., 1121., &
     1022.,-1469., -330., 1256.,    3.,  572.,  523.,  876., &
      628.,  195.,  660.,  -69., -361., -210.,  134.,  -75., &
     -184.,  328., -210.,  264.,   53.,    5.,  -33.,  -86., &
     -124.,  -16.,    3.,   63.,   61.,   -9.,  -11.,   83., &
     -217.,    2.,  -58.,  -35.,   59.,   36.,  -90.,  -69., &
       70.,  -55.,  -45.,    0.,  -13.,   34.,  -10.,  -41., &
       -1.,  -21.,   28.,   18.,  -12.,    6.,  -22.,   11., &
        8.,    8.,   -4.,  -14.,   -9.,    7.,    1.,  -13., &
        2.,    5.,   -9.,   16.,    5.,   -5.,    8.,  -18., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,   -1., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    2.,    4.,    2.,    0.,    0.,   -6.  /) 

  g1(:,2) = (/                                               & ! 1905
   -31464.,-2298., 5909., -728., 2928.,-1086., 1041., 1065., &
     1037.,-1494., -357., 1239.,   34.,  635.,  480.,  880., &
      643.,  203.,  653.,  -77., -380., -201.,  146.,  -65., &
     -192.,  328., -193.,  259.,   56.,   -1.,  -32.,  -93., &
     -125.,  -26.,   11.,   62.,   60.,   -7.,  -11.,   86., &
     -221.,    4.,  -57.,  -32.,   57.,   32.,  -92.,  -67., &
       70.,  -54.,  -46.,    0.,  -14.,   33.,  -11.,  -41., &
        0.,  -20.,   28.,   18.,  -12.,    6.,  -22.,   11., &
        8.,    8.,   -4.,  -15.,   -9.,    7.,    1.,  -13., &
        2.,    5.,   -8.,   16.,    5.,   -5.,    8.,  -18., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    2.,    4.,    2.,    0.,    0.,   -6.  /)

  g1(:,3) = (/                                               & ! 1910
   -31354.,-2297., 5898., -769., 2948.,-1128., 1176., 1000., &
     1058.,-1524., -389., 1223.,   62.,  705.,  425.,  884., &
      660.,  211.,  644.,  -90., -400., -189.,  160.,  -55., &
     -201.,  327., -172.,  253.,   57.,   -9.,  -33., -102., &
     -126.,  -38.,   21.,   62.,   58.,   -5.,  -11.,   89., &
     -224.,    5.,  -54.,  -29.,   54.,   28.,  -95.,  -65., &
       71.,  -54.,  -47.,    1.,  -14.,   32.,  -12.,  -40., &
        1.,  -19.,   28.,   18.,  -13.,    6.,  -22.,   11., &
        8.,    8.,   -4.,  -15.,   -9.,    6.,    1.,  -13., &
        2.,    5.,   -8.,   16.,    5.,   -5.,    8.,  -18., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    2.,    4.,    2.,    0.,    0.,   -6.  /)

  g1(:,4) = (/                                               & ! 1915
   -31212.,-2306., 5875., -802., 2956.,-1191., 1309.,  917., &
     1084.,-1559., -421., 1212.,   84.,  778.,  360.,  887., &
      678.,  218.,  631., -109., -416., -173.,  178.,  -51., &
     -211.,  327., -148.,  245.,   58.,  -16.,  -34., -111., &
     -126.,  -51.,   32.,   61.,   57.,   -2.,  -10.,   93., &
     -228.,    8.,  -51.,  -26.,   49.,   23.,  -98.,  -62., &
       72.,  -54.,  -48.,    2.,  -14.,   31.,  -12.,  -38., &
        2.,  -18.,   28.,   19.,  -15.,    6.,  -22.,   11., &
        8.,    8.,   -4.,  -15.,   -9.,    6.,    2.,  -13., &
        3.,    5.,   -8.,   16.,    6.,   -5.,    8.,  -18., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    1.,    4.,    2.,    0.,    0.,   -6.  /)

  g1(:,5) = (/                                               & ! 1920
   -31060.,-2317., 5845., -839., 2959.,-1259., 1407.,  823., &
     1111.,-1600., -445., 1205.,  103.,  839.,  293.,  889., &
      695.,  220.,  616., -134., -424., -153.,  199.,  -57., &
     -221.,  326., -122.,  236.,   58.,  -23.,  -38., -119., &
     -125.,  -62.,   43.,   61.,   55.,    0.,  -10.,   96., &
     -233.,   11.,  -46.,  -22.,   44.,   18., -101.,  -57., &
       73.,  -54.,  -49.,    2.,  -14.,   29.,  -13.,  -37., &
        4.,  -16.,   28.,   19.,  -16.,    6.,  -22.,   11., &
        7.,    8.,   -3.,  -15.,   -9.,    6.,    2.,  -14., &
        4.,    5.,   -7.,   17.,    6.,   -5.,    8.,  -19., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    9.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    1.,    4.,    3.,    0.,    0.,   -6.  /)

  g1(:,6) = (/                                               & ! 1925
   -30926.,-2318., 5817., -893., 2969.,-1334., 1471.,  728., &
     1140.,-1645., -462., 1202.,  119.,  881.,  229.,  891., &
      711.,  216.,  601., -163., -426., -130.,  217.,  -70., &
     -230.,  326.,  -96.,  226.,   58.,  -28.,  -44., -125., &
     -122.,  -69.,   51.,   61.,   54.,    3.,   -9.,   99., &
     -238.,   14.,  -40.,  -18.,   39.,   13., -103.,  -52., &
       73.,  -54.,  -50.,    3.,  -14.,   27.,  -14.,  -35., &
        5.,  -14.,   29.,   19.,  -17.,    6.,  -21.,   11., &
        7.,    8.,   -3.,  -15.,   -9.,    6.,    2.,  -14., &
        4.,    5.,   -7.,   17.,    7.,   -5.,    8.,  -19., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    9.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    1.,    4.,    3.,    0.,    0.,   -6.  /)

  g1(:,7) = (/                                               & ! 1930
   -30805.,-2316., 5808., -951., 2980.,-1424., 1517.,  644., &
     1172.,-1692., -480., 1205.,  133.,  907.,  166.,  896., &
      727.,  205.,  584., -195., -422., -109.,  234.,  -90., &
     -237.,  327.,  -72.,  218.,   60.,  -32.,  -53., -131., &
     -118.,  -74.,   58.,   60.,   53.,    4.,   -9.,  102., &
     -242.,   19.,  -32.,  -16.,   32.,    8., -104.,  -46., &
       74.,  -54.,  -51.,    4.,  -15.,   25.,  -14.,  -34., &
        6.,  -12.,   29.,   18.,  -18.,    6.,  -20.,   11., &
        7.,    8.,   -3.,  -15.,   -9.,    5.,    2.,  -14., &
        5.,    5.,   -6.,   18.,    8.,   -5.,    8.,  -19., &
        8.,   10.,  -20.,    1.,   14.,  -12.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    9.,    3.,   10.,    0., &
       -2.,   -2.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    1.,    4.,    3.,    0.,    0.,   -6.  /)

  g1(:,8) = (/                                               & ! 1935
   -30715.,-2306., 5812.,-1018., 2984.,-1520., 1550.,  586., &
     1206.,-1740., -494., 1215.,  146.,  918.,  101.,  903., &
      744.,  188.,  565., -226., -415.,  -90.,  249., -114., &
     -241.,  329.,  -51.,  211.,   64.,  -33.,  -64., -136., &
     -115.,  -76.,   64.,   59.,   53.,    4.,   -8.,  104., &
     -246.,   25.,  -25.,  -15.,   25.,    4., -106.,  -40., &
       74.,  -53.,  -52.,    4.,  -17.,   23.,  -14.,  -33., &
        7.,  -11.,   29.,   18.,  -19.,    6.,  -19.,   11., &
        7.,    8.,   -3.,  -15.,   -9.,    5.,    1.,  -15., &
        6.,    5.,   -6.,   18.,    8.,   -5.,    7.,  -19., &
        8.,   10.,  -20.,    1.,   15.,  -12.,    5.,   11., &
       -3.,    1.,   -3.,   -2.,    9.,    3.,   11.,    0., &
       -2.,   -2.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6.  /)

  g1(:,9) = (/                                               & ! 1940
   -30654.,-2292., 5821.,-1106., 2981.,-1614., 1566.,  528., &
     1240.,-1790., -499., 1232.,  163.,  916.,   43.,  914., &
      762.,  169.,  550., -252., -405.,  -72.,  265., -141., &
     -241.,  334.,  -33.,  208.,   71.,  -33.,  -75., -141., &
     -113.,  -76.,   69.,   57.,   54.,    4.,   -7.,  105., &
     -249.,   33.,  -18.,  -15.,   18.,    0., -107.,  -33., &
       74.,  -53.,  -52.,    4.,  -18.,   20.,  -14.,  -31., &
        7.,   -9.,   29.,   17.,  -20.,    5.,  -19.,   11., &
        7.,    8.,   -3.,  -14.,  -10.,    5.,    1.,  -15., &
        6.,    5.,   -5.,   19.,    9.,   -5.,    7.,  -19., &
        8.,   10.,  -21.,    1.,   15.,  -12.,    5.,   11., &
       -3.,    1.,   -3.,   -2.,    9.,    3.,   11.,    1., &
       -2.,   -2.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6.  /)

  g1(:,10) = (/                                              & ! 1945
   -30594.,-2285., 5810.,-1244., 2990.,-1702., 1578.,  477., &
     1282.,-1834., -499., 1255.,  186.,  913.,  -11.,  944., &
      776.,  144.,  544., -276., -421.,  -55.,  304., -178., &
     -253.,  346.,  -12.,  194.,   95.,  -20.,  -67., -142., &
     -119.,  -82.,   82.,   59.,   57.,    6.,    6.,  100., &
     -246.,   16.,  -25.,   -9.,   21.,  -16., -104.,  -39., &
       70.,  -40.,  -45.,    0.,  -18.,    0.,    2.,  -29., &
        6.,  -10.,   28.,   15.,  -17.,   29.,  -22.,   13., &
        7.,   12.,   -8.,  -21.,   -5.,  -12.,    9.,   -7., &
        7.,    2.,  -10.,   18.,    7.,    3.,    2.,  -11., &
        5.,  -21.,  -27.,    1.,   17.,  -11.,   29.,    3., &
       -9.,   16.,    4.,   -3.,    9.,   -4.,    6.,   -3., &
        1.,   -4.,    8.,   -3.,   11.,    5.,    1.,    1., &
        2.,  -20.,   -5.,   -1.,   -1.,   -6.,    8.,    6., &
       -1.,   -4.,   -3.,   -2.,    5.,    0.,   -2.,   -2.  /)

  g1(:,11) = (/                                              & ! 1950
   -30554.,-2250., 5815.,-1341., 2998.,-1810., 1576.,  381., &
     1297.,-1889., -476., 1274.,  206.,  896.,  -46.,  954., &
      792.,  136.,  528., -278., -408.,  -37.,  303., -210., &
     -240.,  349.,    3.,  211.,  103.,  -20.,  -87., -147., &
     -122.,  -76.,   80.,   54.,   57.,   -1.,    4.,   99., &
     -247.,   33.,  -16.,  -12.,   12.,  -12., -105.,  -30., &
       65.,  -55.,  -35.,    2.,  -17.,    1.,    0.,  -40., &
       10.,   -7.,   36.,    5.,  -18.,   19.,  -16.,   22., &
       15.,    5.,   -4.,  -22.,   -1.,    0.,   11.,  -21., &
       15.,   -8.,  -13.,   17.,    5.,   -4.,   -1.,  -17., &
        3.,   -7.,  -24.,   -1.,   19.,  -25.,   12.,   10., &
        2.,    5.,    2.,   -5.,    8.,   -2.,    8.,    3., &
      -11.,    8.,   -7.,   -8.,    4.,   13.,   -1.,   -2., &
       13.,  -10.,   -4.,    2.,    4.,   -3.,   12.,    6., &
        3.,   -3.,    2.,    6.,   10.,   11.,    3.,    8.  /)

  g1(:,12) = (/                                              & ! 1955
   -30500.,-2215., 5820.,-1440., 3003.,-1898., 1581.,  291., &
     1302.,-1944., -462., 1288.,  216.,  882.,  -83.,  958., &
      796.,  133.,  510., -274., -397.,  -23.,  290., -230., &
     -229.,  360.,   15.,  230.,  110.,  -23.,  -98., -152., &
     -121.,  -69.,   78.,   47.,   57.,   -9.,    3.,   96., &
     -247.,   48.,   -8.,  -16.,    7.,  -12., -107.,  -24., &
       65.,  -56.,  -50.,    2.,  -24.,   10.,   -4.,  -32., &
        8.,  -11.,   28.,    9.,  -20.,   18.,  -18.,   11., &
        9.,   10.,   -6.,  -15.,  -14.,    5.,    6.,  -23., &
       10.,    3.,   -7.,   23.,    6.,   -4.,    9.,  -13., &
        4.,    9.,  -11.,   -4.,   12.,   -5.,    7.,    2., &
        6.,    4.,   -2.,    1.,   10.,    2.,    7.,    2., &
       -6.,    5.,    5.,   -3.,   -5.,   -4.,   -1.,    0., &
        2.,   -8.,   -3.,   -2.,    7.,   -4.,    4.,    1., &
       -2.,   -3.,    6.,    7.,   -2.,   -1.,    0.,   -3.  /)

  g1(:,13) = (/                                              & ! 1960
   -30421.,-2169., 5791.,-1555., 3002.,-1967., 1590.,  206., &
     1302.,-1992., -414., 1289.,  224.,  878., -130.,  957., &
      800.,  135.,  504., -278., -394.,    3.,  269., -255., &
     -222.,  362.,   16.,  242.,  125.,  -26., -117., -156., &
     -114.,  -63.,   81.,   46.,   58.,  -10.,    1.,   99., &
     -237.,   60.,   -1.,  -20.,   -2.,  -11., -113.,  -17., &
       67.,  -56.,  -55.,    5.,  -28.,   15.,   -6.,  -32., &
        7.,   -7.,   23.,   17.,  -18.,    8.,  -17.,   15., &
        6.,   11.,   -4.,  -14.,  -11.,    7.,    2.,  -18., &
       10.,    4.,   -5.,   23.,   10.,    1.,    8.,  -20., &
        4.,    6.,  -18.,    0.,   12.,   -9.,    2.,    1., &
        0.,    4.,   -3.,   -1.,    9.,   -2.,    8.,    3., &
        0.,   -1.,    5.,    1.,   -3.,    4.,    4.,    1., &
        0.,    0.,   -1.,    2.,    4.,   -5.,    6.,    1., &
        1.,   -1.,   -1.,    6.,    2.,    0.,    0.,   -7.  /)

  g1(:,14) = (/                                              & ! 1965
   -30334.,-2119., 5776.,-1662., 2997.,-2016., 1594.,  114., &
     1297.,-2038., -404., 1292.,  240.,  856., -165.,  957., &
      804.,  148.,  479., -269., -390.,   13.,  252., -269., &
     -219.,  358.,   19.,  254.,  128.,  -31., -126., -157., &
      -97.,  -62.,   81.,   45.,   61.,  -11.,    8.,  100., &
     -228.,   68.,    4.,  -32.,    1.,   -8., -111.,   -7., &
       75.,  -57.,  -61.,    4.,  -27.,   13.,   -2.,  -26., &
        6.,   -6.,   26.,   13.,  -23.,    1.,  -12.,   13., &
        5.,    7.,   -4.,  -12.,  -14.,    9.,    0.,  -16., &
        8.,    4.,   -1.,   24.,   11.,   -3.,    4.,  -17., &
        8.,   10.,  -22.,    2.,   15.,  -13.,    7.,   10., &
       -4.,   -1.,   -5.,   -1.,   10.,    5.,   10.,    1., &
       -4.,   -2.,    1.,   -2.,   -3.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    4.,   -4.,    4.,    0., &
        0.,   -2.,    2.,    3.,    2.,    0.,    0.,   -6.  /)

  g1(:,15) = (/                                              & ! 1970
   -30220.,-2068., 5737.,-1781., 3000.,-2047., 1611.,   25., &
     1287.,-2091., -366., 1278.,  251.,  838., -196.,  952., &
      800.,  167.,  461., -266., -395.,   26.,  234., -279., &
     -216.,  359.,   26.,  262.,  139.,  -42., -139., -160., &
      -91.,  -56.,   83.,   43.,   64.,  -12.,   15.,  100., &
     -212.,   72.,    2.,  -37.,    3.,   -6., -112.,    1., &
       72.,  -57.,  -70.,    1.,  -27.,   14.,   -4.,  -22., &
        8.,   -2.,   23.,   13.,  -23.,   -2.,  -11.,   14., &
        6.,    7.,   -2.,  -15.,  -13.,    6.,   -3.,  -17., &
        5.,    6.,    0.,   21.,   11.,   -6.,    3.,  -16., &
        8.,   10.,  -21.,    2.,   16.,  -12.,    6.,   10., &
       -4.,   -1.,   -5.,    0.,   10.,    3.,   11.,    1., &
       -2.,   -1.,    1.,   -3.,   -3.,    1.,    2.,    1., &
       -5.,    3.,   -1.,    4.,    6.,   -4.,    4.,    0., &
        1.,   -1.,    0.,    3.,    3.,    1.,   -1.,   -4.  /)

  g1(:,16) = (/                                              & ! 1975
   -30100.,-2013., 5675.,-1902., 3010.,-2067., 1632.,  -68., &
     1276.,-2144., -333., 1260.,  262.,  830., -223.,  946., &
      791.,  191.,  438., -265., -405.,   39.,  216., -288., &
     -218.,  356.,   31.,  264.,  148.,  -59., -152., -159., &
      -83.,  -49.,   88.,   45.,   66.,  -13.,   28.,   99., &
     -198.,   75.,    1.,  -41.,    6.,   -4., -111.,   11., &
       71.,  -56.,  -77.,    1.,  -26.,   16.,   -5.,  -14., &
       10.,    0.,   22.,   12.,  -23.,   -5.,  -12.,   14., &
        6.,    6.,   -1.,  -16.,  -12.,    4.,   -8.,  -19., &
        4.,    6.,    0.,   18.,   10.,  -10.,    1.,  -17., &
        7.,   10.,  -21.,    2.,   16.,  -12.,    7.,   10., &
       -4.,   -1.,   -5.,   -1.,   10.,    4.,   11.,    1., &
       -3.,   -2.,    1.,   -3.,   -3.,    1.,    2.,    1., &
       -5.,    3.,   -2.,    4.,    5.,   -4.,    4.,   -1., &
        1.,   -1.,    0.,    3.,    3.,    1.,   -1.,   -5.  /)

  g1(:,17) = (/                                              & ! 1980
   -29992.,-1956., 5604.,-1997., 3027.,-2129., 1663., -200., &
     1281.,-2180., -336., 1251.,  271.,  833., -252.,  938., &
      782.,  212.,  398., -257., -419.,   53.,  199., -297., &
     -218.,  357.,   46.,  261.,  150.,  -74., -151., -162., &
      -78.,  -48.,   92.,   48.,   66.,  -15.,   42.,   93., &
     -192.,   71.,    4.,  -43.,   14.,   -2., -108.,   17., &
       72.,  -59.,  -82.,    2.,  -27.,   21.,   -5.,  -12., &
       16.,    1.,   18.,   11.,  -23.,   -2.,  -10.,   18., &
        6.,    7.,    0.,  -18.,  -11.,    4.,   -7.,  -22., &
        4.,    9.,    3.,   16.,    6.,  -13.,   -1.,  -15., &
        5.,   10.,  -21.,    1.,   16.,  -12.,    9.,    9., &
       -5.,   -3.,   -6.,   -1.,    9.,    7.,   10.,    2., &
       -6.,   -5.,    2.,   -4.,   -4.,    1.,    2.,    0., &
       -5.,    3.,   -2.,    6.,    5.,   -4.,    3.,    0., &
        1.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6.  /)

  g1(:,18) = (/                                              & ! 1985
   -29873.,-1905., 5500.,-2072., 3044.,-2197., 1687., -306., &
     1296.,-2208., -310., 1247.,  284.,  829., -297.,  936., &
      780.,  232.,  361., -249., -424.,   69.,  170., -297., &
     -214.,  355.,   47.,  253.,  150.,  -93., -154., -164., &
      -75.,  -46.,   95.,   53.,   65.,  -16.,   51.,   88., &
     -185.,   69.,    4.,  -48.,   16.,   -1., -102.,   21., &
       74.,  -62.,  -83.,    3.,  -27.,   24.,   -2.,   -6., &
       20.,    4.,   17.,   10.,  -23.,    0.,   -7.,   21., &
        6.,    8.,    0.,  -19.,  -11.,    5.,   -9.,  -23., &
        4.,   11.,    4.,   14.,    4.,  -15.,   -4.,  -11., &
        5.,   10.,  -21.,    1.,   15.,  -12.,    9.,    9., &
       -6.,   -3.,   -6.,   -1.,    9.,    7.,    9.,    1., &
       -7.,   -5.,    2.,   -4.,   -4.,    1.,    3.,    0., &
       -5.,    3.,   -2.,    6.,    5.,   -4.,    3.,    0., &
        1.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6.  /)

  g1(:,19) = (/                                              & ! 1990
   -29775.,-1848., 5406.,-2131., 3059.,-2279., 1686., -373., &
     1314.,-2239., -284., 1248.,  293.,  802., -352.,  939., &
      780.,  247.,  325., -240., -423.,   84.,  141., -299., &
     -214.,  353.,   46.,  245.,  154., -109., -153., -165., &
      -69.,  -36.,   97.,   61.,   65.,  -16.,   59.,   82., &
     -178.,   69.,    3.,  -52.,   18.,    1.,  -96.,   24., &
       77.,  -64.,  -80.,    2.,  -26.,   26.,    0.,   -1., &
       21.,    5.,   17.,    9.,  -23.,    0.,   -4.,   23., &
        5.,   10.,   -1.,  -19.,  -10.,    6.,  -12.,  -22., &
        3.,   12.,    4.,   12.,    2.,  -16.,   -6.,  -10., &
        4.,    9.,  -20.,    1.,   15.,  -12.,   11.,    9., &
       -7.,   -4.,   -7.,   -2.,    9.,    7.,    8.,    1., &
       -7.,   -6.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    3.,   -2.,    6.,    4.,   -4.,    3.,    0., &
        1.,   -2.,    3.,    3.,    3.,   -1.,    0.,   -6.  /)

  g2(:,1) = (/                                               & ! 1995
   -29692.,-1784., 5306.,-2200., 3070.,-2366., 1681., -413., &
     1335.,-2267., -262., 1249.,  302.,  759., -427.,  940., &
      780.,  262.,  290., -236., -418.,   97.,  122., -306., &
     -214.,  352.,   46.,  235.,  165., -118., -143., -166., &
      -55.,  -17.,  107.,   68.,   67.,  -17.,   68.,   72., &
     -170.,   67.,   -1.,  -58.,   19.,    1.,  -93.,   36., &
       77.,  -72.,  -69.,    1.,  -25.,   28.,    4.,    5., &
       24.,    4.,   17.,    8.,  -24.,   -2.,   -6.,   25., &
        6.,   11.,   -6.,  -21.,   -9.,    8.,  -14.,  -23., &
        9.,   15.,    6.,   11.,   -5.,  -16.,   -7.,   -4., &
        4.,    9.,  -20.,    3.,   15.,  -10.,   12.,    8., &
       -6.,   -8.,   -8.,   -1.,    8.,   10.,    5.,   -2., &
       -8.,   -8.,    3.,   -3.,   -6.,    1.,    2.,    0., &
       -4.,    4.,   -1.,    5.,    4.,   -5.,    2.,   -1., &
        2.,   -2.,    5.,    1.,    1.,   -2.,    0.,   -7., &
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,  &
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,  &
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,  &
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)

  g2(:,2) = (/                                               & ! 2000
   -29619.4,-1728.2, 5186.1,-2267.7, 3068.4,-2481.6, 1670.9, &
     -458.0, 1339.6,-2288.0, -227.6, 1252.1,  293.4,  714.5, &
     -491.1,  932.3,  786.8,  272.6,  250.0, -231.9, -403.0, &
      119.8,  111.3, -303.8, -218.8,  351.4,   43.8,  222.3, &
      171.9, -130.4, -133.1, -168.6,  -39.3,  -12.9,  106.3, &
       72.3,   68.2,  -17.4,   74.2,   63.7, -160.9,   65.1, &
       -5.9,  -61.2,   16.9,    0.7,  -90.4,   43.8,   79.0, &
      -74.0,  -64.6,    0.0,  -24.2,   33.3,    6.2,    9.1, &
       24.0,    6.9,   14.8,    7.3,  -25.4,   -1.2,   -5.8, &
       24.4,    6.6,   11.9,   -9.2,  -21.5,   -7.9,    8.5, &
      -16.6,  -21.5,    9.1,   15.5,    7.0,    8.9,   -7.9, &
      -14.9,   -7.0,   -2.1,    5.0,    9.4,  -19.7,    3.0, &
       13.4,   -8.4,   12.5,    6.3,   -6.2,   -8.9,   -8.4, &
       -1.5,    8.4,    9.3,    3.8,   -4.3,   -8.2,   -8.2, &
        4.8,   -2.6,   -6.0,    1.7,    1.7,    0.0,   -3.1, &
        4.0,   -0.5,    4.9,    3.7,   -5.9,    1.0,   -1.2, &
        2.0,   -2.9,    4.2,    0.2,    0.3,   -2.2,   -1.1, &
       -7.4,    2.7,   -1.7,    0.1,   -1.9,    1.3,    1.5, &
       -0.9,   -0.1,   -2.6,    0.1,    0.9,   -0.7,   -0.7, &
        0.7,   -2.8,    1.7,   -0.9,    0.1,   -1.2,    1.2, &
       -1.9,    4.0,   -0.9,   -2.2,   -0.3,   -0.4,    0.2, &
        0.3,    0.9,    2.5,   -0.2,   -2.6,    0.9,    0.7, &
       -0.5,    0.3,    0.3,    0.0,   -0.3,    0.0,   -0.4, &
        0.3,   -0.1,   -0.9,   -0.2,   -0.4,   -0.4,    0.8, &
       -0.2,   -0.9,   -0.9,    0.3,    0.2,    0.1,    1.8, &
       -0.4,   -0.4,    1.3,   -1.0,   -0.4,   -0.1,    0.7, &
        0.7,   -0.4,    0.3,    0.3,    0.6,   -0.1,    0.3, &
        0.4,   -0.2,    0.0,   -0.5,    0.1,   -0.9 /)

  g2(:,3) = (/                                               & ! 2005
   -29554.63,-1669.05, 5077.99,-2337.24, 3047.69,-2594.50,   &
    1657.76, -515.43, 1336.30,-2305.83, -198.86, 1246.39,    &
     269.72,  672.51, -524.72,  920.55,  797.96,  282.07,    &
     210.65, -225.23, -379.86,  145.15,  100.00, -305.36,    &
    -227.00,  354.41,   42.72,  208.95,  180.25, -136.54,    &
    -123.45, -168.05,  -19.57,  -13.55,  103.85,   73.60,    &
      69.56,  -20.33,   76.74,   54.75, -151.34,   63.63,    &
     -14.58,  -63.53,   14.58,    0.24,  -86.36,   50.94,    &
      79.88,  -74.46,  -61.14,   -1.65,  -22.57,   38.73,    &
       6.82,   12.30,   25.35,    9.37,   10.93,    5.42,    &
     -26.32,    1.94,   -4.64,   24.80,    7.62,   11.20,    &
     -11.73,  -20.88,   -6.88,    9.83,  -18.11,  -19.71,    &
      10.17,   16.22,    9.36,    7.61,  -11.25,  -12.76,    &
      -4.87,   -0.06,    5.58,    9.76,  -20.11,    3.58,    &
      12.69,   -6.94,   12.67,    5.01,   -6.72,  -10.76,    &
      -8.16,   -1.25,    8.10,    8.76,    2.92,   -6.66,    &
      -7.73,   -9.22,    6.01,   -2.17,   -6.12,    2.19,    &
       1.42,    0.10,   -2.35,    4.46,   -0.15,    4.76,    &
       3.06,   -6.58,    0.29,   -1.01,    2.06,   -3.47,    &
       3.77,   -0.86,   -0.21,   -2.31,   -2.09,   -7.93,    &
       2.95,   -1.60,    0.26,   -1.88,    1.44,    1.44,    &
      -0.77,   -0.31,   -2.27,    0.29,    0.90,   -0.79,    &
      -0.58,    0.53,   -2.69,    1.80,   -1.08,    0.16,    &
      -1.58,    0.96,   -1.90,    3.99,   -1.39,   -2.15,    &
      -0.29,   -0.55,    0.21,    0.23,    0.89,    2.38,    &
      -0.38,   -2.63,    0.96,    0.61,   -0.30,    0.40,    &
       0.46,    0.01,   -0.35,    0.02,   -0.36,    0.28,    &
       0.08,   -0.87,   -0.49,   -0.34,   -0.08,    0.88,    &
      -0.16,   -0.88,   -0.76,    0.30,    0.33,    0.28,    &
       1.72,   -0.43,   -0.54,    1.18,   -1.07,   -0.37,    &
      -0.04,    0.75,    0.63,   -0.26,    0.21,    0.35,    &
       0.53,   -0.05,    0.38,    0.41,   -0.22,   -0.10,    &
      -0.57,   -0.18,   -0.82 /)

  g2(:,4) = (/                                                & ! 2010
     -29496.57,-1586.42, 4944.26,-2396.06, 3026.34,-2708.54,  &	
       1668.17, -575.73, 1339.85,-2326.54, -160.40, 1232.10,  &  	
    	251.75,  633.73, -537.03,  912.66,  808.97,  286.48,  &  	
    	166.58, -211.03, -356.83,  164.46,   89.40, -309.72,  &  	
       -230.87,  357.29,   44.58,  200.26,  189.01, -141.05,  &  	
       -118.06, -163.17,   -0.01,   -8.03,  101.04,   72.78,  &  	
    	 68.69,  -20.90,   75.92,   44.18, -141.40,   61.54,  &  	
    	-22.83,  -66.26,   13.10,    3.02,  -78.09,   55.40,  &  	
    	 80.44,  -75.00,  -57.80,   -4.55,  -21.20,   45.24,  &  	
    	  6.54,   14.00,   24.96,   10.46,    7.03,    1.64,  &  	
    	-27.61,    4.92,   -3.28,   24.41,    8.21,   10.84,  &  	
    	-14.50,  -20.03,   -5.59,   11.83,  -19.34,  -17.41,  &  	
    	 11.61,   16.71,   10.85,    6.96,  -14.05,  -10.74,  &  	
    	 -3.54,    1.64,    5.50,    9.45,  -20.54,    3.45,  &  	
    	 11.51,   -5.27,   12.75,    3.13,   -7.14,  -12.38,  &  	
    	 -7.42,   -0.76,    7.97,    8.43,    2.14,   -8.42,  &  	
    	 -6.08,  -10.08,    7.01,   -1.94,   -6.24,    2.73,  &  	
    	  0.89,   -0.10,   -1.07,    4.71,   -0.16,    4.44,  &  	
    	  2.45,   -7.22,   -0.33,   -0.96,    2.13,   -3.95,  &  	
    	  3.09,   -1.99,   -1.03,   -1.97,   -2.80,   -8.31,  &  	
    	  3.05,   -1.48,    0.13,   -2.03,    1.67,    1.65,  &  	
    	 -0.66,   -0.51,   -1.76,    0.54,    0.85,   -0.79,  &  	
    	 -0.39,    0.37,   -2.51,    1.79,   -1.27,    0.12,  &  	
    	 -2.11,    0.75,   -1.94,    3.75,   -1.86,   -2.12,  &  	
    	 -0.21,   -0.87,    0.30,    0.27,    1.04,    2.13,  &  	
    	 -0.63,   -2.49,    0.95,    0.49,   -0.11,    0.59,  &  	
    	  0.52,    0.00,   -0.39,    0.13,   -0.37,    0.27,  &  	
    	  0.21,   -0.86,   -0.77,   -0.23,    0.04,    0.87,  &  	
    	 -0.09,   -0.89,   -0.87,    0.31,    0.30,    0.42,  &  	
    	  1.66,   -0.45,   -0.59,    1.08,   -1.14,   -0.31,  &  	
    	 -0.07,    0.78,    0.54,   -0.18,    0.10,    0.38,  &  	
    	  0.49,    0.02,    0.44,    0.42,   -0.25,   -0.26,  &  	
    	 -0.53,   -0.26,   -0.79/)			  	
  
   g2(:,5) = (/                                               & ! 2015
    -29442.0,-1501.0, 4797.1,-2445.1, 3012.9,-2845.6, 1676.7, &  
      -641.9, 1350.7,-2352.3, -115.3, 1225.6,  244.9,  582.0, &   
      -538.4,  907.6,  813.7,  283.3,  120.4, -188.7, -334.9, &   
       180.9,	70.4, -329.5, -232.6,  360.1,	47.3,  192.4, &   
       197.0, -140.9, -119.3, -157.5,	16.0,	 4.1,  100.2, &   
    	70.0,	67.7,  -20.8,	72.7,	33.2, -129.9,	58.9, &   
       -28.9,  -66.7,	13.2,	 7.3,  -70.9,	62.6,	81.6, &   
       -76.1,  -54.1,	-6.8,  -19.5,	51.8,	 5.7,	15.0, &   
    	24.4,	 9.4,	 3.4,	-2.8,  -27.4,	 6.8,	-2.2, &   
    	24.2,	 8.8,	10.1,  -16.9,  -18.3,	-3.2,	13.3, &   
       -20.6,  -14.6,	13.4,	16.2,	11.7,	 5.7,  -15.9, &   
    	-9.1,	-2.0,	 2.1,	 5.4,	 8.8,  -21.6,	 3.1, &   
    	10.8,	-3.3,	11.8,	 0.7,	-6.8,  -13.3,	-6.9, &   
    	-0.1,	 7.8,	 8.7,	 1.0,	-9.1,	-4.0,  -10.5, &   
    	 8.4,	-1.9,	-6.3,	 3.2,	 0.1,	-0.4,	 0.5, &   
    	 4.6,	-0.5,	 4.4,	 1.8,	-7.9,	-0.7,	-0.6, &   
    	 2.1,	-4.2,	 2.4,	-2.8,	-1.8,	-1.2,	-3.6, &   
    	-8.7,	 3.1,	-1.5,	-0.1,	-2.3,	 2.0,	 2.0, &   
    	-0.7,	-0.8,	-1.1,	 0.6,	 0.8,	-0.7,	-0.2, &   
    	 0.2,	-2.2,	 1.7,	-1.4,	-0.2,	-2.5,	 0.4, &   
    	-2.0,	 3.5,	-2.4,	-1.9,	-0.2,	-1.1,	 0.4, &   
    	 0.4,	 1.2,	 1.9,	-0.8,	-2.2,	 0.9,	 0.3, &   
    	 0.1,	 0.7,	 0.5,	-0.1,	-0.3,	 0.3,	-0.4, &   
    	 0.2,	 0.2,	-0.9,	-0.9,	-0.1,	 0.0,	 0.7, &   
    	 0.0,	-0.9,	-0.9,	 0.4,	 0.4,	 0.5,	 1.6, &   
    	-0.5,	-0.5,	 1.0,	-1.2,	-0.2,	-0.1,	 0.8, &   
    	 0.4,	-0.1,	-0.1,	 0.3,	 0.4,	 0.1,	 0.5, &   
    	 0.5,	-0.3,	-0.4,	-0.4,	-0.3,	-0.8/)	         
!
    g2(1:80,6) = (/                                           & ! 2017
       10.3,   18.1,  -26.6,   -8.7,   -3.3,  -27.4,    2.1,  &  
      -14.1,	3.4,   -5.5,	8.2,   -0.7,   -0.4,  -10.1,  &   
     	1.8,   -0.7,	0.2,   -1.3,   -9.1,	5.3,	4.1,  &   
     	2.9,   -4.3,   -5.2,   -0.2,	0.5,	0.6,   -1.3,  &   
     	1.7,   -0.1,   -1.2,	1.4,	3.4,	3.9,	0.0,  &   
       -0.3,   -0.1,	0.0,   -0.7,   -2.1,	2.1,   -0.7,  &   
       -1.2,	0.2,	0.3,	0.9,	1.6,	1.0,	0.3,  &   
       -0.2,	0.8,   -0.5,	0.4,	1.3,   -0.2,	0.1,  &   
       -0.3,   -0.6,   -0.6,   -0.8,	0.1,	0.2,   -0.2,  &   
     	0.2,	0.0,   -0.3,   -0.6,	0.3,	0.5,	0.1,  &   
       -0.2,	0.5,	0.4,   -0.2,	0.1,   -0.3,   -0.4,  &   
     	0.3,	0.3,	0.0/)
    g2(81:n2,6) = 0.0
!
! Set gh from g1,g2:
!
  do n=1,ncn1
    i = (n-1)*n1
    gh(i+1:i+n1) = g1(:,n)
!   write(6,"('cofrm: n=',i3,' i+1:i+n1=',i4,':',i4)") n,i+1,i+n1
  enddo
  do n=1,ncn2
    i = n1*ncn1 + (n-1)*n2
    gh(i+1:i+n2) = g2(:,n)
!   write(6,"('cofrm: n=',i3,' i+1:i+n2=',i4,':',i4)") n,i+1,i+n2
  enddo
  gh(ngh) = 0. ! not sure why gh is dimensioned with the extra element, so set it to 0. 
  
  if (date < 2015.) then
    t   = 0.2*(date - 1900.)
    ll  = t
    one = ll
    t   = t - one
    if (date < 1995.) then
      nmx   = 10
      nc    = nmx*(nmx+2)
      ll    = nc*ll
      kmx   = (nmx+1)*(nmx+2)/2
    else
      nmx   = 13
      nc    = nmx*(nmx+2)
      ll    = 0.2*(date - 1995.0)
      ll    = 120*19 + nc*ll
      kmx   = (nmx+1)*(nmx+2)/2
    endif
    tc    = 1.0 - t
    if (isv.eq.1) then
      tc = -0.2
      t = 0.2
    endif
  else ! date >= 2015
    t     = date - 2015.0
    tc    = 1.0
    if (isv.eq.1) then
      t = 1.0
      tc = 0.0
    end if
    ll    = 2865
    nmx   = 13
    nc    = nmx*(nmx+2)
    kmx   = (nmx+1)*(nmx+2)/2
  endif ! date < 2015
  r = alt
  l = 1
  m = 1
  n = 0
!
! Set outputs gb(ncoef) and gv(ncoef)
! These are module data above.
! 
  gb(1) = 0.
  gv(1) = 0.
  f0 = -1.e-5
  do k=2,kmx
    if (n < m) then
      m = 0
      n = n+1
    endif ! n < m
    lm = ll + l
    if (m == 0) f0 = f0 * float(n)/2.
    if (m == 0) f  = f0 / sqrt(2.0)
    nn = n+1
    mm = 1      
      
    if (m /= 0) then
      f = f / sqrt(float(n-m+1) / float(n+m) )
      gb(l+1)  = (tc*gh(lm) + t*gh(lm+nc))* f
    else   
      gb(l+1)  = (tc*gh(lm) + t*gh(lm+nc))* f0
    endif  
    gv(l+1) = gb(l+1)/float(nn)
    if (m /= 0) then
      gb(l+2)  = (tc*gh(lm+1) + t*gh(lm+nc+1))*f
      gv(l+2) = gb(l+2)/float(nn)
      l = l+2
    else
      l = l+1
    endif
    m = m+1
  enddo

! write(6,"('cofrm: ncoef=',i4,' gb=',/,(6f12.3))") ncoef,gb
! write(6,"('cofrm: ncoef=',i4,' gv=',/,(6f12.3))") ncoef,gv

end subroutine cofrm
!-----------------------------------------------------------------------
subroutine subsol(iyr,iday,ihr,imn,sec,sbsllat,sbsllon)
!
! Find subsolar geographic latitude and longitude given the
! date and time (Universal Time).
!     
! This is based on formulas in Astronomical Almanac for the
! year 1996, p.  C24. (U.S.  Government Printing Office,
! 1994).  According to the Almanac, results are good to at
! least 0.01 degree latitude and 0.025 degree longitude
! between years 1950 and 2050.  Accuracy for other years has
! not been tested although the algorithm has been designed to
! accept input dates from 1601 to 2100.  Every day is assumed
! to have exactly 86400 seconds; thus leap seconds that
! sometimes occur on June 30 and December 31 are ignored:
! their effect is below the accuracy threshold of the algorithm.
!     
! 961026 A. D. Richmond, NCAR
!
! Input Args:
  integer,intent(in) :: &
    iyr,   & ! Year (e.g., 1994). IYR must be in the range: 1601 to 2100.
    iday,  & ! Day number of year (e.g., IDAY = 32 for Feb 1)
    ihr,   & ! Hour of day    (e.g., 13 for 13:49)
    imn      ! Minute of hour (e.g., 49 for 13:49)
  real,intent(in) :: sec ! Second and fraction after the hour/minute.
!
! Output Args:
  real,intent(out) :: &
    sbsllat, & ! geographic latitude of subsolar point (degrees)
    sbsllon    ! geographic longitude of subsolar point (-180 to +180)
!
! Local:
  integer,parameter :: minyear=1601, maxyear = 2100
  real,parameter :: & ! Use local params for compatability w/ legacy code,
                      ! but probably would be ok to use module data dtr,rtd
    d2r=0.0174532925199432957692369076847, &
    r2d=57.2957795130823208767981548147
  real :: yr,l0,g0,ut,df,lf,gf,l,g,grad,n,epsilon,epsrad,alpha,delta,&
    etdeg,aptime,lambda,lamrad,sinlam
  integer :: nleap,ncent,nrot

  sbsllat=0. ; sbsllon=0.

  yr = iyr-2000
!
! nleap (final) = number of leap days from (2000 January 1) to (IYR January 1)
!                 (negative if iyr is before 1997)
  nleap = (iyr-1601)/4
  nleap = nleap - 99
  if (iyr <= 1900) then
    if (iyr < minyear) then
      write(6,*) 'subsolr invalid before ',minyear,': input year = ',iyr
      stop 'subsolr'
    endif
    ncent = (iyr-minyear)/100
    ncent = 3 - ncent
    nleap = nleap + ncent
  endif
  if (iyr > maxyear) then
    write(6,*) 'subsolr invalid after ',maxyear,':  input year = ',iyr
    stop 'subsolr'
  endif
!
! L0 = Mean longitude of Sun at 12 UT on January 1 of IYR:
!     L0 = 280.461 + .9856474*(365*(YR-NLEAP) + 366*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = 280.461 + .9856474*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = (280.461 - 360.) + (.9856474*365 - 360.)*(YR-4*NLEAP)
!          + (.9856474*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
!
      l0 = -79.549 + (-.238699*(yr-4*nleap) + 3.08514e-2*nleap)
!                 
! G0 = Mean anomaly at 12 UT on January 1 of IYR:
!     G0 = 357.528 + .9856003*(365*(YR-NLEAP) + 366*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = 357.528 + .9856003*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP) 
!          - (ARBITRARY INTEGER)*360.
!        = (357.528 - 360.) + (.9856003*365 - 360.)*(YR-4*NLEAP)
!          + (.9856003*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
!
      g0 = -2.472 + (-.2558905*(yr-4*nleap) - 3.79617e-2*nleap)
!     
! Universal time in seconds:
      ut = float(ihr*3600 + imn*60) + sec
!
! Days (including fraction) since 12 UT on January 1 of IYR:
      df = (ut/86400. - 1.5) + iday
!
! Addition to Mean longitude of Sun since January 1 of IYR: 
      lf = .9856474*df
! 
! Addition to Mean anomaly since January 1 of IYR:
      gf = .9856003*df
! 
! Mean longitude of Sun:
      l = l0 + lf
! 
! Mean anomaly:
      g = g0 + gf
      grad = g*d2r
! 
! Ecliptic longitude:
      lambda = l + 1.915*sin(grad) + .020*sin(2.*grad)
      lamrad = lambda*d2r
      sinlam = sin(lamrad)
! 
! Days (including fraction) since 12 UT on January 1 of 2000:
      n = df + 365.*yr + float(nleap)
! 
! Obliquity of ecliptic: 
      epsilon = 23.439 - 4.e-7*n
      epsrad = epsilon*d2r
! 
! Right ascension:
      alpha = atan2(cos(epsrad)*sinlam,cos(lamrad))*r2d
! 
! Declination:
      delta = asin(sin(epsrad)*sinlam)*r2d
! 
! Subsolar latitude (output argument):
      sbsllat = delta
! 
! Equation of time (degrees):
      etdeg = l - alpha
      nrot = nint(etdeg/360.)
      etdeg = etdeg - float(360*nrot)
! 
! Apparent time (degrees):
! Earth rotates one degree every 240 s.
      aptime = ut/240. + etdeg
!
! Subsolar longitude (output argument):
      sbsllon = 180. - aptime
      nrot = nint(sbsllon/360.)
      sbsllon = sbsllon - float(360*nrot)

end subroutine subsol
!-----------------------------------------------------------------------
subroutine solgmlon(xlat,xlon,colat,elon,mlon)
!
! Compute geomagnetic longitude of the point with geocentric spherical
!  latitude and longitude of XLAT and XLON, respectively.
! 940719 A. D. Richmond, NCAR
!
! Algorithm:
!   Use spherical coordinates. 
!   Let GP be geographic pole. 
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let XLON be longitude of point P.
!   Let TE be colatitude of point P.
!   Let ANG be longitude angle from GM to P.
!   Let TP be colatitude of GM.
!   Let TF be arc length between GM and P.
!   Let PA = MLON be geomagnetic longitude, i.e., Pi minus angle measured
!     counterclockwise from arc GM-P to arc GM-GP. 
!   Then, using notation C=cos, S=sin, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.
!
! Input Args:
  real,intent(in)  :: xlat,xlon,colat,elon
! 
! Output Arg: Geomagnetic dipole longitude of the point (deg, -180. to 180.)
  real,intent(out) :: mlon 
!
! Local:
  real,parameter ::           &
    rtod=5.72957795130823e1,  &
    dtor=1.745329251994330e-2
  real :: ctp,stp,ang,cang,sang,cte,ste,stfcpa,stfspa

  ctp = cos(colat*dtor)
  stp = sqrt(1. - ctp*ctp)
  ang = (xlon-elon)*dtor
  cang = cos(ang)
  sang = sin(ang)
  cte = sin(xlat*dtor)
  ste = sqrt(1.-cte*cte)
  stfcpa = ste*ctp*cang - cte*stp
  stfspa = sang*ste
  mlon = atan2(stfspa,stfcpa)*rtod

end subroutine solgmlon
!-----------------------------------------------------------------------
end module apex
