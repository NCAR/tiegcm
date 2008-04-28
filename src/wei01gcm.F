!      
      module weimer_module
!      
! module used to calculate the Weimer 2001 model potential in both
!  hemispheres from Bz, By, Sws, Swd, and AL.
!
      use params_module,only: nmlat,nmlon,nmlonp1
      use dynamo_module,only: nmlat0,pfrac,phihm
      implicit none
!
! Coefficients for Weimer 2001 from reading coef file:
      integer,parameter :: MJ=3,ML=4,MM=3,MN=2,MO=2
      REAL  CS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
      REAL BCS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:MN)
      REAL  SS( 0:1 , 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
      REAL BSS( 0:1 , 0:1 , 0:MO, 0:1 , 0:MN)
!
! Coefficients for Weimer 2001 after model set with input, and f(rmlt)
      REAL Coef(0:1,0:5,0:5),BoundFit(0:1,0:5)

! blat = BoundaryLat(rmlt):
      real :: blat

! Flag to find the derivatives of the Legendre polynomials for E-field calcs
!  Set to .FALSE. if dynamo=1 or 2, .TRUE. if dynamo=-1 (not implemented)
      logical :: derivative

! Pi constants in single and double precision
      real :: pi
      double precision :: dpi

!  Integers which are set to be the same as MJ,ML,MM,MN,MO
!  Were integer*4 on Sun, but keep just integer on IBM
      integer :: MaxJ,MaxL,MaxM,MaxN,MaxO
!
      contains 
!-----------------------------------------------------------------------
!     subroutine calccloc(phihm)
      subroutine calccloc
!
!  NCAR Nov 02:  Calculate the convection center location (dskofc,offc),
!       radius (theta0), dayside entry for cusp location (phid, poten=0),
!       and nightside exit (phin, poten=0), and the associated auroral
!       radius (arad).  (Leave the aurora center as is, dskofa=0, offa=1.)
!
! Input phihm: High lat model potential in magnetic coordinates (single level).
!
! (dimension 2 is for south, north hemispheres)
!  Calculate crad, offc, dskofc, and phid and phin if possible
!   Use Fig 8 of Heelis et al. [JGR, 85, 3315-3324, 1980]
!     This shows:  arad = 18.7 deg, crad = 16.7 deg (so arad = crad + 2 deg)
!           offa = offc = 3 deg (so offa = offc)
!           dskofc = 2 deg, dskofa = -0.5 deg  (so dskofa = dskofc - 2.5 deg)
!   Parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                             and phin (phin(MLT)=23.50 +/- 0.15By - 12)
!   (In aurora_cons, phid=0., phin=180.*rtd)
!     (For zero By, should be phid=21.39MLT*15*rtd, phin=11.5*15*rtd)
! These are the dimensions and descriptions (corrected phid,n) from aurora.F:
!    |  theta0(2), ! convection reversal boundary in radians
!    |  offa(2),   ! offset of oval towards 0 MLT relative to magnetic pole (rad)
!    |  dskofa(2), ! offset of oval in radians towards 18 MLT (f(By))
!    |  phid(2),   ! dayside convection entrance in MLT-12 converted to radians (f(By))
!    |  phin(2),   ! night convection entrance in MLT-12 converted to radians (f(By))
!    |  rrad(2),   ! radius of auroral circle in radians
!    |  offc(2),   ! offset of convection towards 0 MLT relative to mag pole (rad)
!    |  dskofc(2)  ! offset of convection in radians towards 18 MLT (f(By))
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc)
!
! Additional auroral parameters (see sub aurora_cons):
      use aurora_module,only: theta0,offa,dskofa,phid,phin,rrad,offc,
     |                        dskofc
      use magfield_module,only: sunlons 
      use input_module,only:  ! from user input
     |  byimf    ! By component of IMF (nT)      (e.g., 0.)
      use cons_module,only: rtd,
     |  ylonm,ylatm  ! magnetic grid lons, lats in radians
      implicit none
!
! Args:
!     integer,intent(in) :: nmlat0,nmlat,nmlonp1
!     real,dimension(nmlonp1,nmlat),intent(in) :: phihm ! potential in magnetic
!
! Local:
      integer :: i,i1,i2,ih,j,j1,j2,k
      real :: vnx(2,2),hem,mltdn,mltdx,mltnn,mltnx,mltd,mltn
      integer :: jnx(2,2),inx(2,2),jinx(nmlonp1,2)
      real :: vinx(nmlonp1,2),latinx(nmlonp1,2),mltinx(nmlonp1,2)
      integer :: iflgnn,iflgdx,inm3,inp3,ixm3,ixp3,i06,i18
      real :: offcn,offcx,offcdeg,dskof,ofdc,crad,arad,crad0,craduse
      real :: cosl06,sinl06,colat06,cosm06,cosl18,sinl18,colat18,cosm18
!
!  Look at both hemispheres
      do ih=1,2
!
	if (ih .eq. 1) then
	  j1 = 1
	  j2 = nmlat0
	  hem = -1.
	else
	  j1 = nmlat0 + 1
	  j2 = nmlat
	  hem = 1.
	endif
! Print out un-revised values:
!       write (6,"(1x,'Original convection/oval params (hem,By,off,dsk',
!    |    ',rad,phid,n=',10f9.4)") hem,byimf,offc(ih)*rtd,offa(ih)*rtd,
!    |    dskofc(ih)*rtd,dskofa(ih)*rtd,theta0(ih)*rtd,rrad(ih)*rtd,
!    |    phid(ih)*rtd/15.+12.,phin(ih)*rtd/15.+12.
!  Find min/max
	vnx(ih,1) = 0.
	vnx(ih,2) = 0.
	do j=j1,j2
	  do i=1,nmlonp1-1
	    if (phihm(i,j) .gt. vnx(ih,2)) then
	      vnx(ih,2) = phihm(i,j)
	      jnx(ih,2) = j
	      inx(ih,2) = i
	    endif
	    if (phihm(i,j) .lt. vnx(ih,1)) then
	      vnx(ih,1) = phihm(i,j)
	      jnx(ih,1) = j
	      inx(ih,1) = i
	    endif
	  enddo  !  i=1,nmlonp1-1
	enddo  !  j=j1,j2
!       write (6,"(1x,'ih min/max pot,lat,mlt=',i2,3e12.4,2x,3e12.4))")
!    |    ih,(vnx(ih,i),ylatm(jnx(ih,i)),(ylonm(inx(ih,i))-sunlons(1))
!    |    * rtd/15.+12.,i=1,2)

!  Set default values
      do k=1,2
	do i=1,nmlonp1
	  jinx(i,k) = -99
	  vinx(i,k) = -999.
	  mltinx(i,k) = -999.
	  latinx(i,k) = -999.
	enddo  ! i=1,nmlonp1
      enddo   !  k=1,2

!  Find min/max +/-8 hrs (nmlonp1/3) from peaks and +/-4 lats away
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc) in rad
! Min:
      k = 1
      mltdn = -999.
      mltnn = -999.
      iflgnn = 0
	j1 = jnx(ih,k) - 4
	if (j1 .lt. 1) j1 = 1
	j2 = jnx(ih,k) + 4
	if (j2 .gt. nmlat) j2 = nmlat
	i1 = inx(ih,k) - nmlonp1/3
	if (i1 .lt. 1) i1=1
	i2 = inx(ih,k) + nmlonp1/3
	if (i2 .gt. nmlonp1) i2=nmlonp1
! Look at mid-point part
!       write (6,"(1x,'k j1,2 i1,2=',5i3)") k,j1,j2,i1,i2
	do i=i1,i2
	vinx(i,k) = 0.
	mltinx(i,k) = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
	if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  MLT is 0.3 MLT apart (24/80=0.3)
	if (abs(mltinx(i,k)-18.) .lt. 0.15) i18=i
	  do j=j1,j2
	    if (phihm(i,j) .lt. vinx(i,k)) then
	      vinx(i,k) = phihm(i,j)
	      jinx(i,k) = j
	      latinx(i,k) = ylatm(j) * rtd
	    endif
	  enddo   !  j=j1,j2
	  if (vinx(i,k) .ge. 0.) then
!  Look at vinx=0 for low values of i (decreasing time - phid)
	    if (mltinx(i,k) .gt. 4.5 .and. mltinx(i,k) .lt. 16.5) then
		 mltdn = mltinx(i,k)
	    else
!  Look at vinx=0 for high values of i (increasing time - phin)
	      if (iflgnn .eq. 0) then
		mltnn = mltinx(i,k)
		iflgnn = 1
	      endif
	    endif
	  endif
	enddo   !  i=i1,i2
!         write (6,"(1x,'knx i j v mlt lat =',3i3,3e12.4)") (k,i,
!    |     jinx(i,k),vinx(i,k),mltinx(i,k),latinx(i,k),i=i1,i2)
!  Now look at i<1 for dusk side:
	if (inx(ih,k) - nmlonp1/3 .lt. 1) then
	  i1 = inx(ih,k) - nmlonp1/3 + nmlonp1 - 1
	  i2 = nmlonp1
	  do i=i1,i2
	  vinx(i,k) = 0.
	  mltinx(i,k) = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
	  if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	  if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  MLT is 0.3 MLT apart (24/80=0.3)
	  if (abs(mltinx(i,k)-18.) .lt. 0.15) i18=i
	    do j=j1,j2
	      if (phihm(i,j) .lt. vinx(i,k)) then
		vinx(i,k) = phihm(i,j)
		jinx(i,k) = j
		latinx(i,k) = ylatm(j) * rtd
	      endif
	    enddo   !  j=j1,j2
	    if (vinx(i,k) .ge. 0.) then
!  Look at vinx=0 for low values of i (decreasing time - phid)
	      if (mltinx(i,k) .gt. 4.5 .and. mltinx(i,k) .lt. 16.5)
     |             mltdn = mltinx(i,k)
	    endif
	  enddo   !  i=i1,i2
!         write (6,"(1x,'knx i j v mlt lat =',3i3,3e12.4)") (k,i,
!    |     jinx(i,k),vinx(i,k),mltinx(i,k),latinx(i,k),i=i1,i2)
	endif
!  Now look at i>nmlonp1 for dusk side:
	if (inx(ih,k) + nmlonp1/3 .gt. nmlonp1) then
	  i2 = inx(ih,k) + nmlonp1/3 - nmlonp1 + 1
	  i1 = 1
	  do i=i1,i2
	  vinx(i,k) = 0.
	  mltinx(i,k) = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
	  if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	  if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  MLT is 0.3 MLT apart (24/80=0.3)
	  if (abs(mltinx(i,k)-18.) .lt. 0.15) i18=i
	    do j=j1,j2
	      if (phihm(i,j) .lt. vinx(i,k)) then
		vinx(i,k) = phihm(i,j)
		jinx(i,k) = j
		latinx(i,k) = ylatm(j) * rtd
	      endif
	    enddo   !  j=j1,j2
	  if (vinx(i,k) .ge. 0.) then
!  Look at vinx=0 for high values of i (increasing time - phin)
	    if (mltinx(i,k) .le. 4.5 .or. mltinx(i,k) .ge. 16.5) then
	      if (iflgnn .eq. 0) then
		mltnn = mltinx(i,k)
		iflgnn = 1
	      endif
	    endif
	  endif
	  enddo   !  i=i1,i2
!         write (6,"(1x,'knx i j v mlt lat =',3i3,3e12.4)") (k,i,
!    |     jinx(i,k),vinx(i,k),mltinx(i,k),latinx(i,k),i=i1,i2)
	endif
! Max:
      k = 2
      mltnx = -999.
      mltdx = -999.
      iflgdx = 0
	j1 = jnx(ih,k) - 4
	if (j1 .lt. 1) j1 = 1
	j2 = jnx(ih,k) + 4
	if (j2 .gt. nmlat) j2 = nmlat
	i1 = inx(ih,k) - nmlonp1/3
	if (i1 .lt. 1) i1=1
	i2 = inx(ih,k) + nmlonp1/3
	if (i2 .gt. nmlonp1) i2=nmlonp1
! Look at mid-point part
!       write (6,"(1x,'k j1,2 i1,2=',5i3)") k,j1,j2,i1,i2
	do i=i1,i2
	vinx(i,k) = 0.
	mltinx(i,k) = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
	if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  MLT is 0.3 MLT apart (24/80=0.3)
	if (abs(mltinx(i,k)-6.) .lt. 0.15) i06=i
	  do j=j1,j2
	    if (phihm(i,j) .gt. vinx(i,k)) then
	      vinx(i,k) = phihm(i,j)
	      jinx(i,k) = j
	      latinx(i,k) = ylatm(j) * rtd
	    endif
	  enddo   !  j=j1,j2
	  if (vinx(i,k) .le. 0.) then
!  Look at vinx=0 for low values of i (decreasing time - phin)
	    if (mltinx(i,k) .le. 4.5 .or. mltinx(i,k) .ge. 16.5) then
		 mltnx = mltinx(i,k)
	    else
!  Look at vinx=0 for high values of i (increasing time - phid)
	      if (iflgdx .eq. 0) then
		mltdx = mltinx(i,k)
		iflgdx = 1
	      endif
	    endif
	  endif
	enddo   !  i=i1,i2
!         write (6,"(1x,'knx i j v mlt lat =',3i3,3e12.4)") (k,i,
!    |     jinx(i,k),vinx(i,k),mltinx(i,k),latinx(i,k),i=i1,i2)
!  Now look at i<1 for dawn side:
	if (inx(ih,k) - nmlonp1/3 .lt. 1) then
	  i1 = inx(ih,k) - nmlonp1/3 + nmlonp1 - 1
	  i2 = nmlonp1
	  do i=i1,i2
	  vinx(i,k) = 0.
	  mltinx(i,k) = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
	  if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	  if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  MLT is 0.3 MLT apart (24/80=0.3)
	  if (abs(mltinx(i,k)-6.) .lt. 0.15) i06=i
	    do j=j1,j2
	      if (phihm(i,j) .gt. vinx(i,k)) then
		vinx(i,k) = phihm(i,j)
		jinx(i,k) = j
		latinx(i,k) = ylatm(j) * rtd
	      endif
	    enddo   !  j=j1,j2
!  Look at vinx=0 for low values of i (decreasing time - phin)
	  if (vinx(i,k) .le. 0.) then
	    if (mltinx(i,k) .le. 4.5 .or. mltinx(i,k) .ge. 16.5)
     |           mltnx = mltinx(i,k)
	  endif
	  enddo   !  i=i1,i2
!         write (6,"(1x,'knx i j v mlt lat =',3i3,3e12.4)") (k,i,
!    |     jinx(i,k),vinx(i,k),mltinx(i,k),latinx(i,k),i=i1,i2)
	endif
!  Now look at i>nmlonp1 for dawn side:
	if (inx(ih,k) + nmlonp1/3 .gt. nmlonp1) then
	  i2 = inx(ih,k) + nmlonp1/3 - nmlonp1 + 1
	  i1 = 1
	  do i=i1,i2
	  vinx(i,k) = 0.
	  mltinx(i,k) = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
	  if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	  if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  MLT is 0.3 MLT apart (24/80=0.3)
	  if (abs(mltinx(i,k)-6.) .lt. 0.15) i06=i
	    do j=j1,j2
	      if (phihm(i,j) .gt. vinx(i,k)) then
		vinx(i,k) = phihm(i,j)
		jinx(i,k) = j
		latinx(i,k) = ylatm(j) * rtd
	      endif
	    enddo   !  j=j1,j2
	  if (vinx(i,k) .le. 0.) then
!  Look at vinx=0 for high values of i (increasing time - phid)
	    if (mltinx(i,k) .gt. 4.5 .and. mltinx(i,k) .lt. 16.5) then
	      if (iflgdx .eq. 0) then
		mltdx = mltinx(i,k)
		iflgdx = 1
	      endif
	    endif
	  endif
	  enddo   !  i=i1,i2
!         write (6,"(1x,'knx i j v mlt lat =',3i3,3e12.4)") (k,i,
!    |     jinx(i,k),vinx(i,k),mltinx(i,k),latinx(i,k),i=i1,i2)
	endif
!  Now look at vinx=0 to find phid,n
!  Have mltdx,n from 4.5 to 16.5 MLT (or -999.)
	if (mltdx .ge. 0. .or. mltdn .ge. 0.) then
	  if (mltdx*mltdn .ge. 0.) mltd = 0.5*(mltdx+mltdn)
	  if (mltdx .ge. 0. .and. mltdn .lt. 0.) mltd = mltdx
	  if (mltdn .ge. 0. .and. mltdx .lt. 0.) mltd = mltdn
	else
!  Use parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                             and phin (phin(MLT)=23.50 +/- 0.15By - 12)
	  mltd = 9.39 - hem*0.21*byimf
	endif
	phid(ih) = (mltd-12.) * 15. / rtd
	if (mltnx .ge. 0. .or. mltnn .ge. 0.) then
!  Make mltnx,n from 16.5 to 28.5 MLT (or -999.)
	  if (mltnx .ge. 0. .and. mltnx .lt. 12.) mltnx = mltnx + 24.
	  if (mltnn .ge. 0. .and. mltnn .lt. 12.) mltnn = mltnn + 24.
	  if (mltnx*mltnn .ge. 0.) mltn = 0.5*(mltnx+mltnn)
	  if (mltnx .ge. 0. .and. mltnn .lt. 0.) mltn = mltnx
	  if (mltnn .ge. 0. .and. mltnx .lt. 0.) mltn = mltnn
	else
!  Use parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                             and phin (phin(MLT)=23.50 +/- 0.15By - 12)
	  mltn = 23.50 - hem*0.15*byimf
	endif
	phin(ih) = (mltn-12.) * 15. / rtd
!       write (6,"(1x,'mltd,n,x  mltn,n,x phid,n =',8e12.4)")
!    |   mltd,mltdn,mltdx,mltn,mltnn,mltnx,phid(ih),phin(ih)

!  Estimate dskofc from lat of peak at 6 and 18 MLT (colat(18-6), lat(6-18))
	dskof = abs(latinx(i06,2)) - abs(latinx(i18,1))
	dskofc(ih) = dskof / rtd
!  oval offset is 2.5 deg towards dawn (more neg dskof)
	dskofa(ih) = (dskof-2.5) / rtd
!       write (6,"(1x,'18,06 mlt,lat dskof,c,a=',7e12.4)")
!    |    mltinx(i18,1),latinx(i18,1),mltinx(i06,2),latinx(i06,2),
!    |    dskof,dskofc(ih),dskofa(ih)

!  Estimate offc from lat of peak +/-3 hrs from each maximum
!   (In colat, is nighside-dayside, but in lat is dayside-nightside)
	inm3 = inx(ih,1) - (nmlonp1-1)/8
	inp3 = inx(ih,1) + (nmlonp1-1)/8
	ixm3 = inx(ih,2) - (nmlonp1-1)/8
	ixp3 = inx(ih,2) + (nmlonp1-1)/8
	if (inm3 .lt. 1) inm3 = inm3 + nmlonp1 - 1
	if (inp3 .gt. nmlonp1) inp3 = inp3 - nmlonp1 + 1
	if (ixm3 .lt. 1) ixm3 = ixm3 + nmlonp1 - 1
	if (ixp3 .gt. nmlonp1) ixp3 = ixp3 - nmlonp1 + 1
	offcn = abs(latinx(inm3,1)) - abs(latinx(inp3,1))
	offcx = abs(latinx(ixp3,2)) - abs(latinx(ixm3,2))
!  07/03:  Correction to make sure offset is at least 0.5 deg towards 0 MLT
!       offcdeg = 0.5*(offcn+offcx)
	offcdeg = max(0.5,0.5*(offcn+offcx))
	offc(ih) = offcdeg / rtd
	offa(ih) = offcdeg / rtd
!       write (6,"(1x,'min/max latd3,n3 offc =',8e12.4)")
!    |    latinx(inm3,1),latinx(inp3,1),offcn,latinx(ixp3,2),
!    |    latinx(ixm3,2),offcx,offcdeg,offc(ih)

!  Estimate theta0 from 6-18 MLT line first
	crad0 = 90. - 0.5*abs(latinx(i18,1)+latinx(i06,2))
!  Estimate theta0 from 6-18 MLT line in 'convection circle coordinates'
! Transform to convection circle coordinates:
	ofdc = sqrt(offcdeg**2+dskof**2)
	sinl18 = sin(abs(latinx(i18,1))/rtd)
	cosl18 = cos(abs(latinx(i18,1))/rtd)
	cosm18 = cos(mltinx(i18,1)*15./rtd+asin(dskof/ofdc))
	colat18 = cos(ofdc/rtd)*sinl18-sin(ofdc/rtd)*cosl18*cosm18
	colat18 = acos(colat18)*rtd
!       write (6,"(1x,'18 sinl,cosl,cosm,colat asin=',5e12.4)")
!    |   sinl18,cosl18,cosm18,colat18,asin(dskof/ofdc)
	sinl06 = sin(abs(latinx(i06,2))/rtd)
	cosl06 = cos(abs(latinx(i06,2))/rtd)
	cosm06 = cos(mltinx(i06,2)*15./rtd+asin(dskof/ofdc))
	colat06 = cos(ofdc/rtd)*sinl06-sin(ofdc/rtd)*cosl06*cosm06
	colat06 = acos(colat06)*rtd
!       write (6,"(1x,'06 sinl,cosl,cosm,colat asin=',5e12.4)")
!    |   sinl06,cosl06,cosm06,colat06,asin(dskof/ofdc)
	crad = 0.5*(colat18+colat06)
!  Make sure crad is largest of crad and crad0
	craduse = max(crad,crad0)
	if (craduse .gt. crad) write (6,"(1x,'Used crad0 from 6-18',
     |   'instead of calc crad =',2e12.4)") crad0,crad
	arad = craduse + 2.
	theta0(ih) = craduse / rtd
	rrad(ih) = arad / rtd
!       write (6,"(1x,'radius: 0,18,06,c,a deg,rad=',7e12.4)")crad0,
!    |    colat18,colat06,crad,arad,theta0(ih),rrad(ih)
!
! Print out revised values:
!       write (6,"('Revised convection/oval params (hem,By,off,dsk,',
!    |    'rad,phid,n=',i2,9e12.4)")ih,byimf,offc(ih)*rtd,offa(ih)*rtd,
!    |    dskofc(ih)*rtd,dskofa(ih)*rtd,theta0(ih)*rtd,rrad(ih)*rtd,
!    |    phid(ih)*rtd/15.+12.,phin(ih)*rtd/15.+12.
      enddo  ! ih=1,2
      return
      end subroutine calccloc
!-----------------------------------------------------------------------
!
!     subroutine weimer01(pfrac,phihm,nmlat0,nmlat,nmlon,nmlonp1)
      subroutine weimer01
      use heelis_module,only: colath
      use input_module,only: weimer_ncfile
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc)
      use magfield_module,only: sunlons 
      use cons_module,only: rtd,
     |  ylonm,ylatm  ! magnetic grid lons, lats in radians
      use input_module,only:  ! from user input
     |  byimf,   ! By component of IMF (nT)      (e.g., 0.)
     |  bzimf,   ! Bz component of IMF (nT)      (e.g., 0.)
     |  swvel,   ! solar wind velocity (km/s)    (default 400 km/s)
     |  swden,   ! solar wind density (#/cm3)    (default 4 #/cm3)
     |  AL,      ! AL lower magnetic auroral activity index in nT
		 !   if present, ALUSE=true;  if absent, AL=-20, ALUSE=false
     |  ALUSE    ! locigal to use AL in Weimer model or not (false if no AL)
      use init_module,only: iyear,iday,uthr
!
! Weimer driver, called from sub advance.F
! These routines return pfrac, and phihm to the dynamo, and re-calculates
!  the values of theta0,offc,dskofc,arad,phid,phin,offa,dskofa for dynamics
!   (see argument descriptions below). 
! This is not a module because of difficulties with circular dependencies
!   with the dynamo. Instead, fields are passed through argument lists.
!
! Args:
! pfrac:  Fractional presence of dynamo equation in NH given critical
!           convection colatitudes crit(2) (15,30).
! phihm:  High lat model potential in magnetic coordinates (single level).
!
!     integer,intent(in) :: nmlat0,nmlon,nmlat,nmlonp1
!     real,dimension(nmlonp1,nmlat0)  :: pfrac  ! fraction of potential
!     real,dimension(nmlonp1,nmlat)   :: phihm  ! potential in magnetic
      implicit none
      LOGICAL USEAL
      real :: r2d,by,bz,sws,swd,alind,bt,angl,tilt,ut,hem,htilt,hangl,
     |        rmlt
     |        ,gmlt,rmlat,blats,blatn,phikvs,phikvn,ets,etn,eps,epn
     |        ,et1s,ep1s,et1n,ep1n
      integer :: i,j,iyr,nda,imo,ida
!
! First time through read the coeff file - assuming NOT set indef!!!
      if (MaxJ == 0) then
!
!  Print Copyright for the 2001 Weimer Model
       write (6,"('Copyright 2001 Dan Weimer/MHC potential model')")
       derivative = .FALSE.
       call rd_weimer_coeffs(weimer_ncfile)
       write (6,"(4x,'Read Weimer coefs, MaxJ =',i2)") MaxJ
      endif
!
! Calculate the weimer potential phihm in geomagnetic coordinates in each
!  hemisphere separately
!
      R2D=57.2957795130823208767981548147
!     BY    = 5.
!     BZ = -6.
!     SWS = 400.
!     SWD   =  6.0
!     AL    =  -65.
!     USEAL = .TRUE.
      BY    = byimf
      BZ = bzimf
      SWS = swvel
      SWD   =  swden
      alind = al
      USEAL = ALUSE
      BT    = SQRT (BY*BY + BZ*BZ)
      ANGL  = ATAN2 (BY,BZ)*R2D
      IYR   = iyear
!  Convert from day number to month and day
      nda = iday

!     write(6,"('weimer01: iyr=',i5,' nda=',i4)") iyr,nda
!     write(6,"('weimer01: by=',e12.4,' bz=',e12.4,' sws=',e12.4,
!    |  ' swd=',e12.4,' alind=',e12.4,' useal=',l1)") 
!    |  by,bz,sws,swd,alind,useal

      call CVT2MD(6,IYR,NDA,IMO,IDA)
      ut = uthr
      TILT = GET_TILT (IYR,IMO,ida,ut)

! In epotval01, we find blat=BoundaryLat(rmlt), below which the potential
!  is constant.  EpotVal assumes positive maglat.

!  Do SH calculation
      HEM = -1.
      HTILT = HEM * TILT
      HANGL = HEM * ANGL
      CALL SETMODEL01 (HANGL,BT,HTILT,SWS,SWD,ALind,USEAL)
!     write(6,"('weimer01 after setmodel01: coef=',/,(6e12.4))") coef

      do j=1,nmlat0
        do i=1,nmlon
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc) in rad
          rmlt = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
!  Convert from kV to V
          phihm(i,j) = epotval01(abs(ylatm(j))*rtd,rmlt) * 1000.
       enddo ! i=1,nmlon

!      write(6,"('weimer01: j=',i3,' SH phihm(:,j)=',/,(6e12.4))") 
!    |   j,phihm(:,j)

      enddo ! j=1,nmlat0

!  Do NH calculation
      HEM = 1.
      HTILT = HEM * TILT
      HANGL = HEM * ANGL
      CALL SETMODEL01 (HANGL,BT,HTILT,SWS,SWD,ALind,USEAL)

      do j=nmlat0+1,nmlat
       do i=1,nmlon
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc) in rad
	 rmlt = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
!  Convert from kV to V
	 phihm(i,j) = epotval01(abs(ylatm(j))*rtd,rmlt) * 1000.
       enddo ! i=1,nmlon

!      write(6,"('weimer01: j=',i3,' NH phihm(:,j)=',/,(6e12.4))") 
!    |   j,phihm(:,j)

      enddo ! j=nmlat0+1,nmlat

!
! Periodic points:
      do j=1,nmlat
        phihm(nmlonp1,j) = phihm(1,j)
      enddo ! j=1,nmlat

!  Re-calculate the values of theta0,offc,dskofc,arad,phid,phin,offa,dskofa
!     call calccloc(phihm,nmlat,nmlonp1,nmlat0)
!     call calccloc(phihm)
      call calccloc
!
! Calculate pfrac fractional presence of dynamo equation using critical
!  convection colatitudes crit(2).  (crit is in cons module)
!
!     call colath(pfrac,nmlat0,nmlonp1)
      call colath
!
! From dynamo.F:
!     real,dimension(nmlonp1,nmlat0) ::
!    |   pfrac, ! NH fraction of potential
!    |   phihm  ! potential in magnetic
!
!     do j=1,nmlat0
!       write(6,"('weimer01: j=',i3,' pfrac=',/,(6e12.4))") j,pfrac(:,j)
!       write(6,"('weimer01: j=',i3,' phihm=',/,(6e12.4))") j,phihm(:,j)
!     enddo ! j=1,nmlat0
      end subroutine weimer01
!-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE CVT2MD(MSGUN,IYEAR,NDA,MON,DAY)
C     SUBROUTINE CVT2MD(MSGUN,IYEAR,NDA,MMDD)
C          This sub converts NDA, the day number of the year, IYEAR,
C          into the appropriate month and day of month (integers)
      implicit none
      integer :: msgun,iyear,nda,mon,miss,numd,i,mmdd
      INTEGER DAY
      INTEGER LMON(12)
      PARAMETER (MISS=-32767)
      SAVE LMON
      DATA LMON/31,28,31,30,31,30,31,31,30,31,30,31/
 
      LMON(2)=28
      IF(MOD(IYEAR,4) .EQ. 0)LMON(2)=29
 
      NUMD=0
      DO 100 I=1,12
      IF(NDA.GT.NUMD .AND. NDA.LE.NUMD+LMON(I))GO TO 200
      NUMD=NUMD+LMON(I)
  100 CONTINUE
      WRITE(MSGUN,'('' CVT2MD:  Unable to convert year & day of year'',
     +              I5,'','',I5,''to month & day of month'')')IYEAR,NDA
      MMDD = MISS
      MON = MISS
      DAY = MISS
      RETURN
  200 MON=I
      DAY=NDA-NUMD
      MMDD = MON*100 + DAY
 
      RETURN
!     END
	end subroutine cvt2md
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
*
* Subroutines to calculate the electric potentials from the "Weimer 2K" model of
* the polar cap ionospheric electric potentials described in the publication: 
* Weimer, D. R., An improved model of ionospheric electric potentials including
* substorm perturbations and application to the Geospace Environment Modeling
* November 24, 1996 event, Journal of Geophysical Research, Vol. 106, p. 407, 2001.
*
* To use, first call procedure SETMODEL01 with the specified input parameters:
*   angle: IMF Y-Z clock angle in degrees, 0=northward, 180=southward
*   Bt: Magnitude of IMF in Y-Z plane in nT
*   Tilt: dipole tilt angle in degrees.
*   SWVel: solar wind velocity in km/sec
*   SWDen: solar wind density in #/cc
*   ALindex: (optional) AL index in nT
*
* The function EPOTVAL01(gLAT,gMLT) can then be used repeatively to get the
* electric potential in kV at the desired location.
* Input coordinates assume use of 'altitude adjusted' corrected geomagnetic
* coordinates for R=1, also refered to as AACGM0.
*
* The function BOUNDARYLAT(gMLT) can be used to get the latitude of the boundary
*   where the potential goes to zero.  This boundary is a function of MLT, and
*   varies with the SETMODEL01 parameters.  The potential is zero everywhere below
*   this boundary.
*
* Two data files are provided:
*	'w2klittle.dat' for LITTLE_ENDIAN machines.
*	'w2kbig.dat'    for    BIG_ENDIAN machines.
* You must copy or rename the correct one to the file 'w2k.dat'
*
* This code is protected by copyright and is distributed
* for research or educational use only.
* Commerical use without written permission from Dan Weimer/MRC is prohibited.

CNCARGCM wei01gcm.f 11/02 B. Emery Add to tiegcm1 in dynamics.F
CNCARGCM  The electric fields are not needed so eliminate??
CNCAR      Revisions for use at NCAR:
C            (1) Change behavior at minimum magnetic latitude.  When approaching
C                the model equatorial edge (which varies with MLT) the electric
C                potential returned used to go to zero discontinuously; although
C                intended as a flag, it created artificial gradients in the
C                electric field calculation.  Now the potential returned is
C                constant (that of the minimum latitude) for any latitude at
C                or equatorward of the minimum.
C            (2) Accomodate running simultaneously 1996 and 2001 versions.  To
C                avoid name collisions this required: (i) revising names (e.g.,
C                adding '01') for differing subprograms, and (ii) relocating
C                common routines into another file (weicom.f).
C            (3) Pass the coefficients file name and unit number into READCOEF01
C                rather than using hard coded values.
C            (4) Add wrapper subroutines for non-ANSI trig functions which
C                input angles in degrees.
C eliminate? (5) Add electric field routine (GECMP01) to deterine the electric
C                potential gradient.
C eliminate? (6) Add wrapper routine (WEIEPOT01) for use with AMIE; this is a
C                substitute for calling SETMODEL01 and EPOTVAL01.
C            (7) Remove blanks in some statements to fit in 72 columns to
C                adhere to ANSI Fortran 77.
CNCAR      NCAR changes are delimited by "CNCAR"
*
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
!NCAR   Nov 02: Do not need to change call statement if MaxL,MaxM in module
!       SUBROUTINE DLEGENDRE(x,Plm,dPlm,derivative)
        SUBROUTINE DLEGENDRE(x,lmax,mmax,Plm,dPlm,derivative)
	implicit none
	integer :: lmax,mmax,l,l2,m,mlimit
* compute Double Precision Associate Legendre Function P_l^m(x)
* for all l up to lmax and all m up to mmax.
* Returns results in array Plm.
* If the LOGICAL flag "derivative" is set, then the first derivatives are also
* computed, and put into array dPlm.
* The recursion formulas keep a count of the exponents of the factor SQRT(1-x^2)
*  in both the numerator and denominator, which may cancel each other in the
*  final evaluation, particularly with the derivatives.  This prevents infinities
*  at x=-1 or +1. 
* If X is out of range ( abs(x)>1 ) then value is returns as if x=1.
	DOUBLE PRECISION x,xx,Plm(0:10,0:10),P(0:10,0:10,0:1),fact,sfact
	DOUBLE PRECISION dPlm(0:10,0:10),dP(0:10,0:10,0:2),anum,term
	LOGICAL derivative

	DO l=0,lmax
	    DO m=0,mmax
		  Plm(l,m)=0.D0
		  P(l,m,0)=0.D0
		  P(l,m,1)=0.D0
	    ENDDO
	ENDDO
	IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN
	  Print *,'Bad arguments to DLegendre'
	  RETURN
	ENDIF

* Copy x to xx, and make sure it is in range of -1. to +1.
	xx=MIN(x,1.D0)
	xx=MAX(xx,-1.D0)

	P(0,0,1)=1.D0
	IF(lmax.GT.0) P(1,0,1)=xx
	IF(lmax.GT.1)THEN
	   DO L=2,lmax
	    P(L,0,1)=( (2.D0*L-1)*xx*P(L-1,0,1) - (L-1)*P(L-2,0,1) ) / L
	   ENDDO
	ENDIF

	fact=1.D0-xx**2
	sfact=DSQRT(fact)

	IF(mmax .GT. 0)THEN
		DO M=1,mmax
		  DO L=M,lmax
			L2=MAX( L-2 ,  0 )
			P(L,M,1)= P(L2,M,1) -(2*L-1)*P(L-1,M-1,0)*fact
			P(L,M,0)= P(L2,M,0) -(2*L-1)*P(L-1,M-1,1)
		  ENDDO
	    ENDDO
	ENDIF

	IF(Derivative)Then
* First zero arrays
		DO l=0,lmax
			DO m=0,mmax
				dPlm(l,m)=0.D0
				dP(l,m,0)=0.D0
				dP(l,m,1)=0.D0
				dP(l,m,2)=0.D0
			ENDDO
		ENDDO

		IF(lmax .GT. 0) dP(1,0,1)=1.D0

		IF(lmax .GT. 1)THEN
			DO L=2,lmax  
				dP(L,0,1)=( (2*L-1)*P(L-1,0,1) + 
     $                  (2*L-1)*xx*dP(L-1,0,1) - 
     $                  (L-1)*dP(L-2,0,1) ) / L
			ENDDO
		ENDIF

		IF(mmax .GT. 0)THEN
		  DO M=1,mmax  
		    DO L=M,lmax
		      L2=MAX( L-2 ,  0 )
		      dP(L,M,1)= dP(L2,M,1) - (2*L-1)*fact*dP(L-1,M-1,0)
     $                 - (2*L-1)*dP(L-1,M-1,2) + (2*L-1)*xx*P(L-1,M-1,0)
		      dP(L,M,0)= dP(L2,M,0) - (2*L-1)*dP(L-1,M-1,1)
		      dP(L,M,2)=dP(L2,M,2) +(2*L-1)*xx*P(L-1,M-1,1)
		    ENDDO
		  ENDDO
		ENDIF

		DO L=0,lmax  
	      mlimit=MIN(mmax,L)
		  DO M=0,mlimit
* Prevent a divide by zero
		    anum=dP(L,M,2) !numerator
			IF(sfact.NE.0.)Then !denominator is OK
			  term=anum/sfact 
			ELSE !denominator is zero
			  IF(DABS(anum).LT.1.D-7)THEN
				term=0.D0 !return 0 in cases where numerator is near zero
			  ELSE !return nearly infinity with same sign as numerator
				term=DSIGN(1.D36,anum) 
			  ENDIF
			ENDIF
			dPlm(L,M)=dP(L,M,1) + dP(L,M,0)*sfact + term
		  ENDDO
		ENDDO

	ENDIF !End doing derivative

	DO L=0,lmax
	    mlimit=MIN(mmax,L)
	    DO M=0,mlimit
		  Plm(L,M)=P(L,M,1) + P(L,M,0)*sfact
	    ENDDO
	ENDDO	

	RETURN

!       END
	end subroutine dlegendre
!
!-----------------------------------------------------------------------
! NCAR  Nov 02:  Write 2001 coefficient file in netCDF format for generic
!                use.  Replace READCOEF01 with rd_weimer_coeffs
!-----------------------------------------------------------------------
!
      subroutine rd_weimer_coeffs(ncfile)
      use nchist_module,only: handle_ncerr
      implicit none
!
#ifdef SUN
#include <netcdf3.inc>
#else
#include <netcdf.inc>
#endif
!
! Read netcdf file containing Weimer 2001 coefficients:
!
! Args:
      character(len=*),intent(in) :: ncfile
!
! Local:
      real :: step,stpr,stpd,stp2,cstp,sstp,alamx,d2r,r2d
      integer :: istat,ncid,mjp1,mlp1,mmp1,mnp1,mop1
      integer :: id_mjp1,id_mlp1,id_mmp1,id_mnp1,id_mop1,id_two,
     |  idims5(5),idims6(6)
      integer :: start_5d(5),count_5d(5),start_6d(6),count_6d(6)
      integer :: idv_cs,idv_bcs,idv_ss,idv_bss
      character(len=120) :: msg,desc
      character(len=80) :: dskfile

CNCAR      Feb 01:  Initialize constants used in GECMP01
C          Sep 01:  Omit unneeded min lat variables because of switch to
C                   constant potential at and below min lat (hardcoded
C                   in EPOTVAL01).
      COMMON /CECMP/ ALAMX,STPD,STP2,CSTP,SSTP
C            ALAMX = Absolute max latitude (deg) for normal gradient calc.
C            STPD  = Angular dist (deg) of step @ 300km above earth (r=6371km)
C            STP2  = Denominator in gradient calc
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
!
! Acquire mss or disk file:
      dskfile = ' '
      call getfile(ncfile,dskfile)
!
! Open netcdf file:
      istat = nf_open(dskfile,NF_NOWRITE,ncid)
      if (istat /= NF_NOERR) then
        write(msg,"('Error opening netcdf file ',a)") trim(dskfile)
        call handle_ncerr(istat,msg)
        call shutdown('rd_weimer_coeffs')
      else
        write(6,"('Opened weimer coefficients file ',a)") trim(dskfile)
      endif
!
! Get dimensions and check against dimensions declared in module data:
!
! MJp1:
      istat = nf_inq_dimid(ncid,'MJp1',id_mjp1)
      istat = nf_inq_dimlen(ncid,id_mjp1,mjp1)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting MJp1 dimension')
      if (mjp1 /= mj+1) then
        write(6,"(/,'>>> rd_weimer_coeffs: mjp1 /= mj+1: mjp1=',i3,
     |    ' mj+1=',i3)") mjp1,mj+1 
        call shutdown('mjp1')
      endif
!
! MLp1:
      istat = nf_inq_dimid(ncid,'MLp1',id_mlp1)
      istat = nf_inq_dimlen(ncid,id_mlp1,mlp1)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting MLp1 dimension')
      if (mlp1 /= ml+1) then
        write(6,"(/,'>>> rd_weimer_coeffs: mlp1 /= ml+1: mlp1=',i3,
     |    ' ml+1=',i3)") mlp1,ml+1 
        call shutdown('mlp1')
      endif
!
! MMp1:
      istat = nf_inq_dimid(ncid,'MMp1',id_mmp1)
      istat = nf_inq_dimlen(ncid,id_mmp1,mmp1)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting MMp1 dimension')
      if (mmp1 /= mm+1) then
        write(6,"(/,'>>> rd_weimer_coeffs: mmp1 /= mm+1: mmp1=',i3,
     |    ' mm+1=',i3)") mmp1,mm+1 
        call shutdown('mmp1')
      endif
!
! MNp1:
      istat = nf_inq_dimid(ncid,'MNp1',id_mnp1)
      istat = nf_inq_dimlen(ncid,id_mnp1,mnp1)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting MNp1 dimension')
      if (mnp1 /= mn+1) then
        write(6,"(/,'>>> rd_weimer_coeffs: mnp1 /= mn+1: mnp1=',i3,
     |    ' mn+1=',i3)") mnp1,mn+1 
        call shutdown('mnp1')
      endif
!
! MOp1:
      istat = nf_inq_dimid(ncid,'MOp1',id_mop1)
      istat = nf_inq_dimlen(ncid,id_mop1,mop1)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting MOp1 dimension')
      if (mop1 /= mo+1) then
        write(6,"(/,'>>> rd_weimer_coeffs: mop1 /= mo+1: mop1=',i3,
     |    ' mo+1=',i3)") mop1,mo+1 
        call shutdown('mop1')
      endif
!
! Get coefficient fields:
!
! CS:
      istat = nf_inq_varid(ncid,'CS',idv_cs)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting CS var id')
      istat = nf_get_var_double(ncid,idv_cs,cs) ! use with 8-byte reals
!     istat = nf_get_var_real(ncid,idv_cs,cs)   ! use with 4-byte reals
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting variable CS')
!     write(6,"('rd_weimer_coeffs: cs=',/,(6e12.4))") cs
!
! BCS:
      istat = nf_inq_varid(ncid,'BCS',idv_bcs)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting BCS var id')
      istat = nf_get_var_double(ncid,idv_bcs,bcs) ! use with 8-byte reals
!     istat = nf_get_var_real(ncid,idv_bcs,bcs)   ! use with 4-byte reals
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting variable BCS')
!     write(6,"('rd_weimer_coeffs: bcs=',/,(6e12.4))") bcs
!
! SS:
      istat = nf_inq_varid(ncid,'SS',idv_ss)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting SS var id')
      istat = nf_get_var_double(ncid,idv_ss,ss) ! use with 8-byte reals
!     istat = nf_get_var_real(ncid,idv_ss,ss)   ! use with 4-byte reals
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting variable SS')
!     write(6,"('rd_weimer_coeffs: ss=',/,(6e12.4))") ss
!
! BSS:
      istat = nf_inq_varid(ncid,'BSS',idv_bss)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting BSS var id')
      istat = nf_get_var_double(ncid,idv_bss,bss) ! use with 8-byte reals
!     istat = nf_get_var_real(ncid,idv_bss,bss)   ! use with 4-byte reals
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_weimer_coeffs: Error getting variable BSS')
!     write(6,"('rd_weimer_coeffs: bss=',/,(6e12.4))") bss
!
      istat = nf_close(ncid)

!  Assign common integer values in INTVALS
      MaxJ = MJp1 - 1
      MaxL = MLp1 - 1
      MaxM = MMp1 - 1
      MaxN = MNp1 - 1
      MaxO = MOp1 - 1

!  Find pi, single and double precision
      pi=2.*ASIN(1.)
      dpi=2.D0*DASIN(1.D0)

!  Find constants for stepping in lat/lon to find E field from potential
!   in GECMP01 (Probably will not be used in GCM code.)
!   STEP = 10 km gives some large E fields  -- use 100 km instead
      STEP = 10.
      STEP = 100.
      STPR = STEP/6671.
      STPD = STPR*R2D
      STP2 = 2.*STEP
      CSTP = COS (STPR)
      SSTP = SQRT (1. - CSTP*CSTP)
      ALAMX = 90. - STPD

      end subroutine rd_weimer_coeffs
!-----------------------------------------------------------------------

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
CNCAR      Sep 01: Change name from SetModel to be distinct

!NCAR      Nov 02: Change COMMONS to correlate with Fortran 90 module
!                  structure with routine rd_weimer_coeffs to read the coefficients
!                  on a netCDF file.

	SUBROUTINE SetModel01 (angle,Bt,Tilt,SWVel,SWDen,ALindex,UseAL)
C       SUBROUTINE SetModel   (angle,Bt,Tilt,SWVel,SWDen,ALindex,UseAL)
	implicit none
	real :: SinTilt,omega,acoef
	integer :: klimit,mlimit
CNCAR
*
* Calculate the complete set of spherical harmonic coeficients,
* given an aribitrary IMF angle (degrees from northward toward +Y),
* magnitude Bt (nT), dipole tilt angle (degrees), 
* solar wind velocity (km/sec), SWDen (#/cc),
* ALindex (nT), and Logical flag to use optional AL index.
*
* Sets the value of Coef and Boundfit in the common block SetW2kCoef.
*
	REAL angle,Bt,Tilt,SWVel,SWDen,ALindex
C       LOGICAL First,UseAL
	LOGICAL UseAL

CNCAR      Feb 01: Move logic block to subroutine RdCoef01
C       DATA First/.TRUE./
C       SAVE First
C       INTEGER unit
C       CHARACTER*15 cfile
C       CHARACTER*30 Copyright
!       PARAMETER (MJ=3,ML=4,MM=3,MN=2,MO=2)
!       REAL  CS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
!       REAL BCS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:MN)
!       REAL  SS( 0:1 , 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
!       REAL BSS( 0:1 , 0:1 , 0:MO, 0:1 , 0:MN)
	REAL XA(0:MJ),XB(0:MJ),FSC(0:1,0:4),PSS(0:1)
!       REAL Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi
!       DOUBLE PRECISION dpi
!!      INTEGER*4 i,j,k,l,m,n,o
	INTEGER   i,j,k,l,m,n,o
!       INTEGER*4 maxj,MaxL,MaxM,MaxN,MaxO
!       COMMON /AllW2kCoefs/MaxJ,MaxO,CS,BCS,SS,BSS
!       COMMON /SetW2kCoef/MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi

CNCAR      Feb 01: Move logic block to subroutine RdCoef01
C  All First=TRUE in new subroutine ReadCoef
C       If(First)Then
C       cfile='w2k.dat' !make sure correct ENDIAN type file is used.
C       unit=99
C       OPEN(UNIT=unit,FILE=cfile,STATUS='OLD',form='UNFORMATTED')
C       READ(unit) Copyright
C       PRINT *,Copyright
C       READ(unit) Maxj,MaxL,MaxM,MaxN,MaxO
C       If(maxj.NE.MJ .OR. MaxL.NE.ML .OR. MaxM.NE.MM .OR.
C    $   MaxN.NE.MN .OR. MaxO.NE.MO)Then
C               PRINT *,'Data File Error'
C               STOP !Data file did not match sixe expected for arrays
C       Endif
C       READ(unit) CS
C       READ(unit) BCS
C       READ(unit) SS
C       READ(unit) BSS
C       CLOSE(unit)
C       pi=2.*ASIN(1.)
C       dpi=2.D0*DASIN(1.D0)
C       First=.FALSE.
C       Endif
CNCAR

	SinTilt=SIND(Tilt)
	omega=angle*pi/180.
	XA(0)=1.
	XA(1)=Bt**(2./3.) *SWvel
	XA(2)=SinTilt
	XA(3)=SWvel**2 *SWDen
	XB(0)=1.
	XB(1)=Bt
	XB(2)=SinTilt
	XB(3)=SWvel**2 *SWDen

	DO l=0,MaxL
	    mlimit=MIN(l,MaxM)
	    DO m=0,mlimit
		  klimit=MIN(m,1)
		  DO k=0,klimit
			acoef=0. !rezero for summation
			DO j=0,MaxJ
			    DO o=0,MaxO
				  DO i=0,1
					FSC(i,o)=CS(j,i,o,k,l,m)
				  ENDDO
			    ENDDO
     			    acoef=acoef+ XA(j)*FSVAL(omega,MaxO,fsc)
			ENDDO
			IF(UseAL)THEN
			    DO j=0,1
				  DO o=0,MaxO
					DO i=0,1
					    FSC(i,o)=SS(j,i,o,k,l,m)
					ENDDO
				  ENDDO
				  PSS(j)=FSVAL(omega,MaxO,fsc)
			    ENDDO
			    acoef=acoef + PSS(0) + PSS(1)*ALindex
			ENDIF
			Coef(k,l,m)=acoef
! 5/7/04 btf: Coef is NaNs and zeroes.
!                       write(6,"('setmodel01: k=',i3,' l=',i3,' m=',i3,
!    |                    ' coef=',e12.4)") k,l,m,coef(k,l,m)
		  ENDDO
	    ENDDO
	ENDDO

	DO n=0,MaxN
	    klimit=MIN(n,1)
	    DO k=0,klimit
		  acoef=0. !rezero for summation
		  DO j=0,MaxJ
			DO o=0,MaxO
			    DO i=0,1
				  FSC(i,o)=BCS(j,i,o,k,n)
			    ENDDO
			ENDDO
     			acoef=acoef+ XB(j)*FSVAL(omega,MaxO,fsc)
		  ENDDO
		  IF(UseAL)THEN
			DO j=0,1
			    DO o=0,MaxO
				  DO i=0,1
					FSC(i,o)=BSS(j,i,o,k,n)
				  ENDDO
			    ENDDO
			    PSS(j)=FSVAL(omega,MaxO,fsc)
			ENDDO
			acoef=acoef + PSS(0) + PSS(1)*ALindex
		  ENDIF
		  BoundFit(k,n)=acoef
	    ENDDO
	ENDDO
!       write(6,"('end setmodel01: Coef=',/,(6e12.4))") Coef
	end subroutine setmodel01
****************** Copyright 1996, 2001, Dan Weimer/MRC ***********************
	FUNCTION BoundaryLat(gmlt)
	implicit none
	real :: BoundaryLat
	REAL gmlt
!NCAR      Nov 02: Change COMMONS to correlate with Fortran 90 module
!                  structure with routine rd_weimer_coeffs to read the coefficients
!                  on a netCDF file.
!       REAL Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi
!       DOUBLE PRECISION dpi
!       INTEGER MaxL,MaxM,MaxN
!       COMMON /SetW2kCoef/MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi
	BoundaryLat=FSVal(gmlt*pi/12.,MaxN,BoundFit)
	RETURN
!       END
	end function boundarylat
****************** Copyright 1996, 2001, Dan Weimer/MRC ***********************
CNCAR      Sep 01: Change name from EpotVal to be distinct

	FUNCTION EPOTVAL01 (gLAT,gMLT)
C       FUNCTION EpotVal   (gLAT,gMLT)
	implicit none
	real :: epotval01,glatlim,colat,bcolat
	integer :: l,m,mlimit
CNCAR

* Return the value of the electric potential in kV at
* corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).
*
* Must first call ReadCoef and SetModel to set up the model coeficients for
* the desired values of Bt, IMF clock angle, Dipole tilt angle, and SW Vel.
*
      REAL gLAT,gMLT
      DOUBLE PRECISION Phi,Z,O,x,ct,Phim,DC
      DOUBLE PRECISION Plm(0:10,0:10), OPlm(0:10,0:10)

!NCAR      Nov 02: Change COMMONS to correlate with Fortran 90 module
!                  structure with routine rd_weimer_coeffs to read the coefficients
!                  on a netCDF file.  Add dPlm and dOPlm.
!     REAL Coef(0:1,0:5,0:5), BoundFit(0:1,0:5),pi
!     DOUBLE PRECISION dpi
!     COMMON /SetW2kCoef/ MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi
      DOUBLE PRECISION dPlm(0:10,0:10), dOPlm(0:10,0:10)
!     logical :: derivative

      blat = BoundaryLat (gMLT)
CNCAR      Feb 01:  For latitudes at and equatorward of the model minimum
C          (a function of MLT), limit the latitude used to never be less
C          than the model minimum, thus returning a constant potential for
C          points at and equatorward of BoundaryLat (gmlt).
C     IF (glat .GT. blat) THEN
      glatlim = max (glat,blat)
CNCAR

      Phi = DBLE (gMLT)*dpi/12.D0

CNCAR
      colat = 90.-glatlim
C     colat = 90.-glat
CNCAR

      bcolat  = 90.-blat
      x  = DBLE (colat)*dpi/DBLE(bcolat)
      DC = DBLE (Coef(0,0,0))
      Z  = DC
      O  = DC
      ct = DCOS(x)

!     write(6,"('epotval01: init z=',e12.4,' o=',e12.4)") z,o

!NCAR   Nov 02:  change 0.D0 to dPlm and dOPlm and .FALSE. to derivative
      CALL DLegendre (ct,   MaxL,MaxM,Plm ,dPlm,derivative)
!     CALL DLegendre (ct,Plm ,dPlm,.FALSE.)
!     CALL DLegendre (ct,Plm ,dPlm,derivative)
C          Also find value at outer boundary at angle Pi, or cos(pi)=-1.
!     CALL DLegendre (-1.D0,MaxL,MaxM,OPlm ,0.D0,.FALSE.)
      CALL DLegendre (-1.D0,MaxL,MaxM,OPlm ,dOPlm,derivative)
!     CALL DLegendre (-1.D0,OPlm ,dOPlm,derivative)
      DO l=1,MaxL

! btf 5/7/04: coef(0,l,0) is NaN
!       write(6,"('epotval01: l=',i4,' plm(l,0)=',e12.4,' coef(0,l,0)=',
!    |    e12.4)") l,plm(l,0),coef(0,l,0)

	Z = Z +  Plm(l,0)*DBLE(Coef(0,l,0))
	O = O + OPlm(l,0)*DBLE(Coef(0,l,0))
	mlimit = MIN(l,MaxM)
	DO m=1,mlimit
	  phim = phi*m
	  Z = Z + Plm(l,m)*(DBLE(Coef(0,l,m))*DCOS(Phim) +
     $                      DBLE(Coef(1,l,m))*DSIN(Phim) )
	  O = O +OPlm(l,m)*DBLE(Coef(0,l,m))

!         write(6,"('epotval01: l=',i4,' m=',i3,' z=',e12.4,' o=',
!    |      e12.4)") l,m,z,o

	ENDDO
      ENDDO
CNCAR
      EPOTVAL01 = SNGL(Z-O)
C     EpotVal   = SNGL(Z-O)
C     ELSE
C       EpotVal=0.
C     ENDIF
CNCAR
      RETURN
!     END
	end function epotval01

CNCAR
      SUBROUTINE GECMP01 (AMLA,RMLT,ET,EP)
!NCAR   Nov 02:  GECMP01 can give large E fields when STP2 is small, use 200km
      implicit none
      real :: amla,rmlt,et,ep,d2r,r2d,alamx,stpd,stp2,cstp,sstp,
     |        xmlt,xmlt1,amla1,p1,p2,dphi
      integer :: kpol
C          Get Electric field components for the 2001 Weimer electrostatic
C          potential model.  Before use, first load coefficients (CALL
C          READCOEF01) and initialize model conditions (CALL SETMODEL01).
C          The electric field, the gradient of the electic potential,
C          is determined using finite differences over a distance of
C          STEP km, defined in READCOEF01 and accessed here via common
C          block CECMP.
C
C          It has been useful to modify STEP.  One application discovered
C          the Weimer electrostatic potential can have minor unrealistic
C          perturbations over short distances.  To avoid these ripples
C          STEP was increased to 5 degrees arc (582 km at 300 km altitude,
C          R=6671 km).
C
C          INPUTS:
C            AMLA = Absolute value of magnetic latitude (deg)
C            RMLT = Magnetic local time (hours).
C          RETURNS:
C            ET = Etheta (magnetic equatorward*) E field component (V/m)
C            EP = Ephi   (magnetic eastward)     E field component (V/m)
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., when ET > 0, the direction is
C              southward when in northern magnetic hemisphere
C              northward when in southern magnetic hemisphere
C          Since the Weimer model uses only NH (positive) latitudes, the
C          sign of ET for SH should be changed outside of this routine.
C
C          HISTORY:
C          Jan 97:  Initial implementation at NCAR for the 1996 version.
C          Feb 01:  Error corrected in determining DPHI.  Old version
C          used wrong spherical right triangle formula (inadvertently
C          computing along a latitude line); more importantly, it
C          neglected to convert mlt from degrees to hour angle input
C          to EPOTVAL01.
C          Sep 01:  Revised equatorward boundary logic to the scheme Barbara
C          uses for AMIE:  As absolute magnetic latitude decreases to the
C          model minimum, the electric potential from EPOTVAL01 now goes to
C          a non-zero constant (rather than zero), thus obviating the need
C          for special handling here.  Therefore, special logic for determining
C          gradients at lower limit has been removed.

      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)

C          CECMP contains constants initialized in READCOEF01
      COMMON /CECMP/ ALAMX,STPD,STP2,CSTP,SSTP

      ET = -99999.
      EP = -99999.
      IF (AMLA .LT. 0.) GO TO 100

C          Calculate -(latitude gradient) by stepping along the magnetic
C          latitude line in each direction (flipping coordinates when
C          going over pole to keep lat <= 90).
      KPOL  = 0
      XMLT  = RMLT
   10 XMLT1 = XMLT
      AMLA1 = AMLA + STPD
      IF (AMLA1 .GT. 90.) THEN
	AMLA1 = 180. - AMLA1
	XMLT1 = XMLT1 + 12.
      ENDIF
      P1 = EPOTVAL01 (AMLA1    ,XMLT1)
      P2 = EPOTVAL01 (AMLA-STPD,XMLT )
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STP2

C          Calculate -(lon gradient).  Step along the magnetic meridion
C          in both directions to obtain the electric potential
      IF (AMLA .LT. ALAMX) THEN
	AMLA1 = ASIN(SIN(AMLA*D2R)*CSTP)
	DPHI  = ASIN (SSTP/COS(AMLA1))*R2D/15.      ! 15 converts from degrees to hours
	AMLA1 = AMLA1*R2D
	P1 = EPOTVAL01 (AMLA1,XMLT+DPHI)
	P2 = EPOTVAL01 (AMLA1,XMLT-DPHI)
      ELSE
	AMLA = 90.
	XMLT = XMLT + 6.
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STP2
      IF (KPOL .EQ. 1) EP = -EP

  100 RETURN
!     END
	end subroutine gecmp01

      SUBROUTINE WEIEPOT01 (IYR,IMO,IDA,IHR,IMN,SWS,SWD,BY,BZ,AL,USEAL,
     +                     IHEM,ISETM,RMLA,RMLT, ET,EP,EPOT)
      implicit none
      real :: sws,swd,by,bz,al,rmla,rmlt,et,ep,epot,r2d,angl,tilt,h,
     |        bt,hr,hangl,htilt,amla
      integer :: iyr,imo,ida,ihr,imn,ihem,isetm
C          Interface to Weimer-01 (a.k.a w2k) for AMIE's combined electrostatic
C          potential models.  This replaces two calls (SETMODEL01 and
C          EPOTVAL01) and automates sign changes for southern hemisphere.
C          INPUTS:
C            IYR   = UT year ((4-digit)
C            IMO   = month of IYR
C            IDA   = day of IMO
C            IHR   = hour of day
C            IMN   = min of hour
C            SWS   = Solar wind speed (km/s)
C            SWD   = Solar wind density (#/cm3)
C            BY    = IMF By component in GSM coordinates (nt)
C            BZ    = IMF Bz component in GSM coordinates (nt)
C            AL    = AL index (nT)
C            USEAL = logical (T or F) to use AL index or not
C            IHEM  = Hemisphere flag: (-1) southern, (1) northern
C            ISETM = Model conditions change flag: (0) no-change
C                                                  (1) time, IMF or SW changed
C            RMLA  = Magnetic latitude  (deg) of point to determine Potential
C                    (RMLA should be NH positive latitudes only, since SH
C                    values found by changing sign of By (or ANGL) and tilt.)
C            RMLT  = Magnetic longitude (hrs) of point to determine Potential
C          RETURNS:
C            ET   = Etheta (magnetic equatorward*) E field component (V/m)
C            EP   = Ephi   (magnetic eastward)     E field component (V/m)
C            EPOT = Electrostatic potential (kV)
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., when ET > 0, the direction is
C              southward when in northern magnetic hemisphere
C              northward when in southern magnetic hemisphere
C
C          Since the AMIE model assumes a NH solution although the
C          latitudes are negative for SH data, the sign of ET
C          for the SH should be changed outside of this routine.

C          HISTORY:
C          Jan 97:  Initial implementation. B. Emery.
C          Feb 01:  Remove bug (HR was not defined)
C          Sep 01:  Add common block s.t. computed values are available to
C          the calling routine without changing the argument list
      COMMON /WEIAT/ ANGL, TILT
      PARAMETER (R2D=57.2957795130823208767981548147)
      LOGICAL USEAL

      H = REAL (IHEM)

      IF (ISETM .EQ. 1) THEN
	ANGL = ATAN2 (BY,BZ)*R2D
	BT   = SQRT (BY*BY + BZ*BZ)
	HR = FLOAT(IHR) + FLOAT(IMN)/60.
	TILT = GET_TILT (IYR,IMO,IDA,HR)

	HANGL = H * ANGL
	HTILT = H * TILT
	CALL SETMODEL01 (HANGL,BT,HTILT,SWS,SWD,AL,USEAL)
C     WRITE (6,'('' WEIMER-2001: '',I4,5I3,8F7.2,L)') IYR,IMO,IDA,IHR,
C    +                    IMN,IHEM,SWS,SWD,AL,BY,BZ,BT,HANGL,HTILT,USEAL
      ENDIF

C          NH assumed latitudes only
      AMLA  = ABS (RMLA)
      EPOT = EPOTVAL01 (AMLA,RMLT)
      CALL GECMP01 (AMLA,RMLT,ET,EP)

      RETURN
!     END
	end subroutine weiepot01
CNCAR
C          Routines common to both 1996 and 2001 versions of Weimer's
C          Ionospheric Electrostatic Potential Model

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
	FUNCTION FSVal(omega,MaxN,FSC)
	implicit none
	real :: FSVal
* Fourier Series Value
* Return value of Sine/Cosine Fourier series for N terms up to MaxN
* at angle omega, given the F.S. coeficients in array FSC
	REAL omega,FSC(0:1,0:*)
	INTEGER MaxN,n
	REAL Y,theta
	Y=0.
	DO n=0,MaxN
	  theta=omega*n
	  Y=Y + FSC(0,n)*COS(theta) + FSC(1,n)*SIN(theta)
	ENDDO
	FSVal=Y
	RETURN
!       END
	end function fsval
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
* COORDINATE TRANSFORMATION UTILITIES

CNCAR      Feb 01:  Changed TRANS to GET_TILT s.t. the dipole tilt angle is
C          returned.

	FUNCTION GET_TILT (YEAR,MONTH,DAY,HOUR)
C       SUBROUTINE TRANS(YEAR,MONTH,DAY,HOUR,IDBUG)
	implicit none
	real :: GET_TILT,B3,B32,B33
	integer :: IYR,JD,MJD,I,J,K
CNCAR

	INTEGER YEAR,MONTH,DAY,IDBUG
	REAL HOUR
C         
C      THIS SUBROUTINE DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
C      TRANSFORMATIONS, IDENTIFIED BY K.
C          K=1 TRANSFORMS GSE to GEO
C          K=2     "      GEO to MAG
C          K=3     "      GSE to MAG
C          K=4     "      GSE to GSM
C          K=5     "      GEO to GSM
C          K=6     "      GSM to MAG
C          K=7     "      GSE to GEI
C          K=8     "      GEI to GEO
C          K=9     "      GSM to SM 
C	   K=10    "      GEO to SM 
C	   K=11    "      MAG to SM 
C
C      IF IDBUG IS NOT 0, THEN OUTPUTS DIAGNOSTIC INFORMATION TO
C      FILE UNIT=IDBUG
C	
	INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
	INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
	INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
	INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
	INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG

	PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
	PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
	PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
	PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
	PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
	PARAMETER (MAGSM =11,SMMAG =-11)
C
C      The formal names of the coordinate systems are:
C	GSE - Geocentric Solar Ecliptic
C	GEO - Geographic
C	MAG - Geomagnetic
C	GSM - Geocentric Solar Magnetospheric
C	SM  - Solar Magnetic
C	
C      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
C      ST(I) AND CT(I) ARE SINES & COSINES.       
C
C      Program author:  D. R. Weimer
C
C      Some of this code has been copied from subroutines which had been
C      obtained from D. Stern, NASA/GSFC.  Other formulas are from "Space 
C      Physics Coordinate Transformations: A User Guide" by M. Hapgood (1991).
C
C      The formulas for the calculation of Greenwich mean sidereal time (GMST)
C      and the sun's location are from "Almanac for Computers 1990",
C      U.S. Naval Observatory.
C
	REAL UT,T0,GMSTD,GMSTH,ECLIP,MA,LAMD,SUNLON

CNCAR      Feb 01:  Eliminate unused routines from translib.for: ROTATE,
C          ROTATEV, FROMCART, TOCART, MLT, MAGLONG, SUNLOC.  Remaining
C          are ADJUST and JULDAY
CNCAR      Nov 02: Commons MFIELD and TRANSDAT now only in TRANS (GET_TILT)
C                  So eliminate them as commons.  For Fortran 90, eliminate
C                  the DATA statement for assignments (not block_data)
C       COMMON/MFIELD/EPOCH,TH0,PH0,DIPOLE
C       COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
C         
C       DATA EPOCH,TH0,PH0,DIPOLE/1980.,11.19,-70.76,.30574/
	REAL EPOCH,TH0,PH0,DIPOLE
	REAL CX(9),ST(6),CT(6),AM(3,3,11)
C         
C  TH0 = geog co-lat of NH magnetic pole
C  PH0 = geog longitude of NH magnetic pole
C  DIPOLE = magnitude of the B field in gauss at the equator
	EPOCH = 1980.
	TH0 = 11.19
	PH0 = -70.76
	DIPOLE = 0.30574
CNCAR

CNCAR      Feb 01:  Prevent debug printing to a file
	IDBUG = 0
CNCAR

	IF(YEAR.LT.1900)THEN
	  IYR=1900+YEAR
	ELSE
	  IYR=YEAR
	ENDIF
	UT=HOUR
	JD=JULDAY(MONTH,DAY,IYR)
	MJD=JD-2400001
	T0=(FLOAT(MJD)-51544.5)/36525.0
	GMSTD=100.4606184 + 36000.770*T0 + 3.87933E-4*T0*T0 +
     $        15.0410686*UT
	CALL ADJUST(GMSTD)
	GMSTH=GMSTD*24./360.
	ECLIP=23.439 - 0.013*T0
	MA=357.528 + 35999.050*T0 + 0.041066678*UT
	CALL ADJUST(MA)
	LAMD=280.460 + 36000.772*T0 + 0.041068642*UT
	CALL ADJUST(LAMD)
	SUNLON=LAMD + (1.915-0.0048*T0)*SIND(MA) + 0.020*SIND(2.*MA)
	CALL ADJUST(SUNLON)
	IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) YEAR,MONTH,DAY,HOUR
	  WRITE(IDBUG,*) 'MJD=',MJD
	  WRITE(IDBUG,*) 'T0=',T0
	  WRITE(IDBUG,*) 'GMSTH=',GMSTH
	  WRITE(IDBUG,*) 'ECLIPTIC OBLIQUITY=',ECLIP
	  WRITE(IDBUG,*) 'MEAN ANOMALY=',MA
	  WRITE(IDBUG,*) 'MEAN LONGITUDE=',LAMD
	  WRITE(IDBUG,*) 'TRUE LONGITUDE=',SUNLON
	ENDIF

	CX(1)= GMSTD
	CX(2) = ECLIP
	CX(3) = SUNLON
	CX(4) = TH0
	CX(5) = PH0
c Derived later:
c       CX(6) = Dipole tilt angle  
c       CX(7) = Angle between sun and magnetic pole
c       CX(8) = Subsolar point latitude
c       CX(9) = Subsolar point longitude

	DO I=1,5
	  ST(I) = SIND(CX(I))
	  CT(I) = COSD(CX(I))
	ENDDO
C         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0.         
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)      
C         
      AM(1,1,GEIGEO) = CT(1)      
      AM(1,2,GEIGEO) = ST(1)      
      AM(1,3,GEIGEO) = 0.         
      AM(2,1,GEIGEO) = -ST(1)     
      AM(2,2,GEIGEO) = CT(1)      
      AM(2,3,GEIGEO) = 0.         
      AM(3,1,GEIGEO) = 0.         
      AM(3,2,GEIGEO) = 0.         
      AM(3,3,GEIGEO) = 1.         
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) +
     $    AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +  AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
C         
      AM(1,1,GEOMAG) = CT(4)*CT(5) 
      AM(1,2,GEOMAG) = CT(4)*ST(5) 
      AM(1,3,GEOMAG) =-ST(4)       
      AM(2,1,GEOMAG) =-ST(5)       
      AM(2,2,GEOMAG) = CT(5)       
      AM(2,3,GEOMAG) = 0.
      AM(3,1,GEOMAG) = ST(4)*CT(5) 
      AM(3,2,GEOMAG) = ST(4)*ST(5) 
      AM(3,3,GEOMAG) = CT(4)       
C         
      DO I=1,3   
      DO J=1,3   
       AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) +
     $   AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +  AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
C         
      B32 = AM(3,2,GSEMAG)         
      B33 = AM(3,3,GSEMAG)         
      B3  = SQRT(B32*B32+B33*B33)       
      IF (B33.LE.0.) B3 = -B3    
C         
      AM(2,2,GSEGSM) = B33/B3      
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)   
      AM(3,2,GSEGSM) = B32/B3      
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)   
      AM(1,1,GSEGSM) = 1.
      AM(1,2,GSEGSM) = 0.
      AM(1,3,GSEGSM) = 0.
      AM(2,1,GSEGSM) = 0.
      AM(3,1,GSEGSM) = 0.
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) +
     $    AM(I,2,GSEGSM)*AM(J,2,GSEGEO) + AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) +
     $   AM(I,2,GEOMAG)*AM(J,2,GEOGSM) + AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO
C
	ST(6) = AM(3,1,GSEMAG)       
	CT(6) = SQRT(1.-ST(6)*ST(6))      
	CX(6) = ASIND(ST(6))     

        AM(1,1,GSMSM) = CT(6)
        AM(1,2,GSMSM) = 0.
        AM(1,3,GSMSM) = -ST(6)
        AM(2,1,GSMSM) = 0.
        AM(2,2,GSMSM) = 1.
        AM(2,3,GSMSM) = 0.
        AM(3,1,GSMSM) = ST(6)
        AM(3,2,GSMSM) = 0.
        AM(3,3,GSMSM) = CT(6)
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOSM) = AM(I,1,GSMSM)*AM(1,J,GEOGSM) +
     $    AM(I,2,GSMSM)*AM(2,J,GEOGSM) +  AM(I,3,GSMSM)*AM(3,J,GEOGSM)
      ENDDO
      ENDDO
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,MAGSM) = AM(I,1,GSMSM)*AM(J,1,GSMMAG) +
     $   AM(I,2,GSMSM)*AM(J,2,GSMMAG) + AM(I,3,GSMSM)*AM(J,3,GSMMAG)
      ENDDO
      ENDDO
C
      CX(7)=ATAN2D( AM(2,1,11) , AM(1,1,11) )
      CX(8)=ASIND( AM(3,1,1) )
      CX(9)=ATAN2D( AM(2,1,1) , AM(1,1,1) )

      IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) 'Dipole tilt angle=',CX(6)
	  WRITE(IDBUG,*) 'Angle between sun and magnetic pole=',CX(7)
	  WRITE(IDBUG,*) 'Subsolar point latitude=',CX(8)
	  WRITE(IDBUG,*) 'Subsolar point longitude=',CX(9)

        DO K=1,11
         WRITE(IDBUG,1001) K
         DO I=1,3
           WRITE(IDBUG,1002) (AM(I,J,K),J=1,3)
         ENDDO
        ENDDO
 1001   FORMAT(' ROTATION MATRIX ',I2)
 1002   FORMAT(3F9.5)
      ENDIF

CNCAR      Mar 96: return the dipole tilt from this function call.
      GET_TILT = CX(6)
CNCAR

      RETURN
!     END
	end function get_tilt
******************************************************************************
CNCAR      Feb 01:  Eliminate unused routines from translib.for: ROTATE,
C          ROTATEV, FROMCART, TOCART, MLT, MAGLONG, SUNLOC.  Remaining
C          are ADJUST and JULDAY
CNCAR
	SUBROUTINE ADJUST(ANGLE)
	implicit none
	real :: angle
C	ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
 10	CONTINUE
	IF(ANGLE.LT.0.)THEN
	  ANGLE=ANGLE+360.
	  GOTO 10
	ENDIF
 20	CONTINUE
 	IF(ANGLE.GE.360.)THEN
	  ANGLE=ANGLE-360.
	  GOTO 20
	ENDIF
	RETURN
!       END
	end subroutine adjust
******************************************************************************
      FUNCTION JULDAY(MM,ID,IYYY)
      implicit none
      integer :: igreg, iyyy, mm, id, jy, jm, ja, julday
      PARAMETER (IGREG=15+31*(10+12*1582))
      IF (IYYY.EQ.0) PAUSE 'There is no Year Zero.'
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF
      RETURN
!     END
	end function julday

CNCAR      Routines added to work around non-ANSI trig functions which
C          input degrees instead of radians:  SIND, COSD, ASIND, ATAN2D

      FUNCTION SIND (DEG)
      implicit none
      real :: sind,d2r,r2d,deg
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      SIND = SIN (DEG * D2R)
      RETURN
!     END
	end function sind

      FUNCTION COSD (DEG)
      implicit none
      real :: cosd,d2r,r2d,deg
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)

      COSD = COS (DEG * D2R)
      RETURN
!     END
	end function cosd

      FUNCTION ASIND (RNUM)
      implicit none
      real :: asind,d2r,r2d,rnum
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      ASIND = R2D * ASIN (RNUM)
      RETURN
!     END
	end function asind

      FUNCTION ATAN2D (RNUM1,RNUM2)
      implicit none
      real :: atan2d,d2r,r2d,rnum1,rnum2
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      ATAN2D = R2D * ATAN2 (RNUM1,RNUM2)
      RETURN
!     END
	end function atan2d
CNCAR

!-----------------------------------------------------------------------
      end module weimer_module
