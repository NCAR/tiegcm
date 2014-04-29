       subroutine calccloc (istep,avmincrad,model_ctpoten)
! ... Change log .......................................................

! Nov 02:  Calculate the convection center location (dskofc,offc),
!       radius (theta0), dayside entry for cusp location (phid, poten=0),
!       and nightside exit (phin, poten=0), and the associated auroral
!       radius (arad).  (Leave the aurora center as is, dskofa=0, offa=1.)
! Jan 11:  bae:  calccloc has a problem with Bz>0, |Bz/By|>1 conditions where
!       multiple cells are possible, so for these conditions, set defaults of:
!       theta0 = 10 deg, offc = 4.2 deg, and dskofc = 0 deg (offc and dskofc from 2005 model)
! Dec 11: bae: Revise for CMIT tests using MIX potS(27,181) potN(27,181)
!       27 co-lats in radians or from 90 (pole), 89.9,.. to 51.662, 50.199 variable deg
!       181 theta from 0 to 2pi (2deg) where 0=noon MLT for NH and 0=midnight MLT for SH
! Jan 12: bae: Move calccloc from wei01gcm.F to util.F; divide dskofc by 2; 
!       check if 'nogood' and use defaults if: MLT(max)>12, MLT(min)<12, dsckofc<-10 or >+10,
!       offc<-5 or >+10;  set defaults from 2005 Weimer model, but defaults not based on IMF
! Aug 12: bae: Revise for ifbad (0 if good; 1, 2, or 3 for various failures) each hem
!       and put in opposite hem (flipped for By) if only 1 bad, or default if both bad.
!       NOTE:  These revisions for ifbad.ne.0 can be revised outside this routine as desired
! Nov 12 Dec 12 Jan 13: Major cleanup.  See SVN commit log for
!       details. Most of the above comments are now irellevant.
! May 13 bae: revisions for ifbad=0 or 1 outside calccloc and in revcloc.
!     Added more printout for istep and fixed problems:
!      (1) kmlt_min,max integer not real, 
!      (2) i1,i2 for knx1 bad from j_min not i_min,  
!      (3) inx used for inm,p ixm,p not i_min,max (deleted inx) 
!      (4) j1 and j2 not revised for NH (ih=2) so always using SH phihm to calc theta0,
!           and when j1,j2 correct for phihm, change mlatdeg(nmlat) from -90 to +90 mlat
! Aug 13 bae: found ifbad=2 for NH mostly since checked j2 for nmlat_subset instead of nmlat
!
! ... Description ......................................................
!       ih=1,2 for SH and NH
! Calculate crad, offc, dskofc if possible
! Use Fig 8 of Heelis et al. [JGR, 85, 3315-3324, 1980]
!     This shows:  arad = 18.7 deg, crad = 16.7 deg (so arad = crad + 2 deg)
!           offa = offc = 3 deg (so offa = offc)
!           dskofc = 2 deg, dskofa = -0.5 deg  (so dskofa = dskofc - 2.5 deg)
!
! ... Use Association ..................................................
      use aurora_module,only:  ! dimension (2) is for south, north hemispheres
     |     theta0, ! theta0(2), ! convection reversal boundary in radians
     |     offc,   ! offc(2),   ! offset of convection towards 0 MLT relative to mag pole (rad)
     |     dskofc  ! dskofc(2)  ! offset of convection in radians towards 18 MLT (f(By))
      use magfield_module,only: sunlons  ! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc)
      use cons_module,only: 
     |     ylonm,ylatm, ! magnetic grid longitudes (nmlonp1) and latitudes (nmlat) in radians
     |     rtd, dtr
      use params_module,only: 
     |     nmlat,nmlonp1 ! phihm dimensions
      use dynamo_module,only: phihm ! Input phihm(potS,N): High lat model potential in magnetic coordinates (single level).

      implicit none
! ... Arguments ........................................................  
      integer,intent(in) :: istep  
      real,intent(in) :: avmincrad    
      real,intent(out) :: model_ctpoten(2)
! ... Constants ........................................................      
      integer,parameter :: ifwr=0 ! ifwr=1,2 is a print flag for extra output
      integer,parameter :: nmlat_subset=nmlat/4-1
! ... Local variables ..................................................     
! 05/13 bae:  Use the TIEGCM distorted magnetic grid of ylatm,ylonm with phihm 
!             Do not use the real mlat/mlt grids of AMIE or CMIT etc
!      real,dimension(nmlonp1) ::
!     |     mltsh, ! CMIT/MIX: hours from noon clockwise PM to AM to noon, 
!                  ! AMIE: hours from 0 to 24 MLT
!     |     mltnh  ! CMIT/MIX: hours from midnight clockwise AM to PM to midnight 
!                  ! AMIE: hours from 0 to 24 MLT!
!      real,dimension(nmlat_subset) :: 
!     |     ylatm_deg_nh, ! Magnetic latitude in Northern Hemisphere (degrees) equ to +90
!     |     ylatm_deg_sh  ! Magnetic latitude in Southern Hemisphere (degrees) -90 to equ

      integer :: i,i1,i2,ih,j,j1,j2,k
      real :: phihm_min(2), phihm_max(2), kmlt_min(2), kmlt_max(2)
      integer j_min(2), i_min(2)
      integer j_max(2), i_max(2)  
      integer :: ifbad(2)    

      integer :: jinx(nmlonp1,2)
      real :: vinx(nmlonp1,2),latinx(nmlonp1,2),mltinx(nmlonp1,2)
      integer :: inm3,inp3,ixm3,ixp3,i06,i18
      real :: offcn,offcx,offcdeg,dskof,ofdc,crad,crad0,craduse
      real :: asind,dmlthalf,latnm3,latnp3,latxm3,latxp3
      real :: mlatdeg(nmlat),mltarr(nmlonp1)		
      real :: dskofcen(2),offcen(2),radcen(2)
      real :: rhon,rhox,thetandeg,thetaxdeg,xn,yn,xx,yx,xcen,ycen
      real :: x06,y06,rho06, cosl06,sinl06,colat06,cosm06, ! 06 => dawn
     |     x18,y18,rho18, cosl18,sinl18,colat18,cosm18, ! 18 => dusk
     |     cradcoord

! ... Begin ............................................................

! Calculate magnetic local time.
! Find MLT from sunlons (sunlons(1-nlat) are all the same value so use the first value)
      do i=1,nmlonp1
       mltarr(i) = (ylonm(i)-sunlons(1)) * rtd / 15. + 12.
	if (mltarr(i) .gt. 24.) mltarr(i) = mltarr(i) - 24.
	if (mltarr(i) .lt.  0.) mltarr(i) = mltarr(i) + 24.
      enddo

! !  Get points from poles to about +/-35 mlat or use nmlat_subset=nmlat/4-1
! nmlat_subset=23 from -90 about every 2 deg to -70, every 3 to -30:
!  23 magS = -0.9000E+02 -0.8812E+02 -0.8624E+02 -0.8433E+02 -0.8240E+02 
!            -0.8043E+02 -0.7841E+02 -0.7633E+02 -0.7419E+02 -0.7197E+02
!            -0.6967E+02 -0.6728E+02 -0.6480E+02 -0.6222E+02 -0.5955E+02 
!            -0.5678E+02 -0.5393E+02 -0.5100E+02 -0.4801E+02 -0.4498E+02
!            -0.4192E+02 -0.3885E+02 -0.3582E+02

      do j=1,nmlat
       mlatdeg(j) = ylatm(j)*rtd
      enddo
      if (istep<2) write (6,"(1x,'mlatdeg =',10f6.2)") mlatdeg
      if (istep<2) write (6,"(1x,'mltarr =',10f6.2)") mltarr
    
!
!  Look at both hemispheres (ih=1 SH, ih=2 NH)
      do ih=1,2
 
	if (ih .eq. 1) then
!  SH (ih=1) j1 and j2 are used often, so must redefine them for each ih
          j1 = 1
          j2 = nmlat_subset
	else
!  NH (ih=2)
          j1 = nmlat + 1 - nmlat_subset
          j2 = nmlat
        endif
!
! Print out un-revised values:
        if (ifwr .ge. 1 .or. istep<3) write (6,"(1x,
     | 'Original convection params istep,ih offc,dskofc,crad =',2i5,
     |    3f9.4)") istep,ih,offc(ih)*rtd,dskofc(ih)*rtd,
     |    theta0(ih)*rtd

!  Find min/max values & location of potential
! ih = hemisphere selection
	phihm_min(ih) =  99999999.0
	phihm_max(ih) = -99999999.0
	do j=j1,j2
	  do i=1,nmlonp1-1
	    if (phihm(i,j) .gt. phihm_max(ih)) then
	      phihm_max(ih) = phihm(i,j)
	      j_max(ih) = j
	      i_max(ih) = i
              kmlt_max(ih) = mltarr(i)
              if (abs(mlatdeg(j)) .gt. 89.99) kmlt_max(ih) = 6. !FIXME: Magic numbers
	    endif
	    if (phihm(i,j) .lt. phihm_min(ih)) then
	      phihm_min(ih) = phihm(i,j)
	      j_min(ih) = j
	      i_min(ih) = i
              kmlt_min(ih) = mltarr(i)
              if (abs(mlatdeg(j)) .gt. 89.99) kmlt_min(ih) = 18. !FIXME: Magic numbers
	    endif
	  enddo  !  i=1,nmlonp1-1
	enddo  !  j=j1,j2

! 02/10: Calculate the model ctpoten in kV from model (max - min) in V
	model_ctpoten(ih) = 0.001 * (phihm_max(ih) - phihm_min(ih))

       if (ifwr .ge. 1 .or. istep<3) 
     | write (6,"(1x,'istep ih i,j,min,mlt,mlat i,j,max,mlt,mlat CP =',
     |  i5,1x,3i4,3f8.2,1x,2i4,3f8.2,1x,f8.2)")
     |  istep,ih,i_min(ih),j_min(ih),phihm_min(ih)*0.001,kmlt_min(ih),
     |  mlatdeg(j_min(ih)),i_max(ih),j_max(ih),
     |  phihm_max(ih)*0.001,kmlt_max(ih),mlatdeg(j_max(ih)),
     |  model_ctpoten(ih)

!  Feb 2012 bae:  Find cartesian coordinates for min/max, and assume center of circle fit is the midpoint
!   theta = atan2(y,x), rho = sqrt(x*x+y*y);  x=rho*cos(theta), y=rho*sin(theta)
!   theta = (MLT-6)*360/24, rho=colat
      rhon = 90.-abs(mlatdeg(j_min(ih)))
      rhox = 90.-abs(mlatdeg(j_max(ih)))
      thetandeg = (mltarr(i_min(ih)) - 6.)*360./24.
      thetaxdeg = (mltarr(i_max(ih)) - 6.)*360./24.
      xn = rhon*cos(thetandeg*dtr)
      yn = rhon*sin(thetandeg*dtr)
      xx = rhox*cos(thetaxdeg*dtr)
      yx = rhox*sin(thetaxdeg*dtr)
      dskofcen(ih) = -(xx - 0.5*(xx-xn))
      offcen(ih) = -0.5*(yx+yn)
!  Center in x,y coordinates is y=-offcen(ih) and x=-dskofcen(ih)
      xcen = -dskofcen(ih)
      ycen = -offcen(ih)
      radcen(ih) = sqrt( (xx-xcen)*(xx-xcen) + (yx-ycen)*(yx-ycen) )
       if (ifwr .ge. 1 .or. istep<3) 
     | write (6,"(1x,'ifbad rhon,x thetan,x x,yn x,yx off,dsk,',
     |  'radc =',2i2,11f6.1)") ifbad,rhon,rhox,thetandeg,
     |  thetaxdeg,xn,yn,xx,yx,offcen(ih),dskofcen(ih),radcen(ih)
 
!  Feb 2012 bae:  Calculate ifbad for test of MLT min>12. and MLT max<12. (bad min<12., max>12.)
        ifbad(ih) = 0
        if (kmlt_min(ih) .lt. 12. .or. kmlt_max(ih) .gt. 12.) then
          ifbad(ih) = 1
!  Skip the rest of the calculation but put defaults in at the end for this or other ifbad cases
	  cycle
        endif

!  When ifcen=1, skip calculation of dskoffc and offc (set by revcloc to 0 and 1.1), and use radcen as crad
!  NOTE:  This is a BAD idea since radcen is usually worse than crad so ifcen is deleted
!!        if (ifcen==1) cycle

! Feb and Aug 2012:  Need to set ifbad(ih) > 1 for other problems with undef dskof or offcdeg w then bad crad when use the -999 instead of something good!

!  Set default values
      do k=1,2
	do i=1,nmlonp1
	  jinx(i,k) = -999.
	  vinx(i,k) = -999.
!	  mltinx(i,k) = -999.	! mltinx is always defined later on, so don't check on it alone
	  latinx(i,k) = -999.
	enddo  ! i=1,nmlonp1
      enddo   !  k=1,2
!  Find ~half-3/4 d(MLT) spacing for testing near edge (MLT spaced 24./(nmlonp1-1) hours apart)
      dmlthalf = 0.75 * 24./(nmlonp1-1)
      if (istep<3) write (6,"(1x,'ifbad dmlthalf =',2i2,f8.2)") 
     |  ifbad,dmlthalf

!  Find min/max +/-8 hrs (nmlonp1/3) from peaks and +/-4 lats away
!  Feb 2012:  If have small cell (less than 6 MLT wide), then latinx remains -999 
!    if go to the sign of the other cell.
!    Therefore, need to find 3 hours or less away from min/max to where latinx is not -999
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc) in rad - NOT USED!
! Min:
      k = 1
      i06 = -999
      i18 = -999
	j1 = j_min(ih) - 4
	if (j1 .lt. 1) j1 = 1		! do not go beyond S Pole at j=1
	j2 = j_min(ih) + 4
	if (j2 .gt. nmlat) j2 = nmlat	! do not go beyond N Pole at j=nmlat
	i1 = i_min(ih) - nmlonp1/3
	if (i1 .lt. 1) i1=1
	i2 = i_min(ih) + nmlonp1/3
	if (i2 .gt. nmlonp1) i2=nmlonp1
! Look at mid-point part
        if (ifwr .ge. 1 .or. istep<3) 
     |    write (6,"(1x,'ih k j1,2 i1,2=',2i2,2i3,2i4)") 
     |     ih,k,j1,j2,i1,i2
	do i=i1,i2
	vinx(i,k) = 0.
	mltinx(i,k) = mltarr(i)
	if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  dmlthalf*2=d(MLT) - find edge ~18MLT
	if (abs(mltinx(i,k)-18.) .lt. dmlthalf) i18=i
	  do j=j1,j2
	    if (phihm(i,j) .lt. vinx(i,k)) then
	      vinx(i,k) = phihm(i,j)
	      jinx(i,k) = j
              latinx(i,k) = mlatdeg(j)
	    endif
	  enddo   !  j=j1,j2
	enddo   !  i=i1,i2
      if (ifwr .ge. 2 .or. istep<2) write (6,"(1x,'knx1 i18 =',i4)")i18
       if (ifwr .ge. 2 .or. istep<2)
     |   write (6,"(1x,'knx1 i j v mlt dmlt18 lat =',3i4,4e12.4)") 
     |     (k,i,jinx(i,k),vinx(i,k),mltinx(i,k),abs(mltinx(i,k)-18.),
     |     latinx(i,k),i=i1,i2)

!  Now look at i<1 for dusk side:
	if (i_min(ih) - nmlonp1/3 .lt. 1) then
	  i1 = i_min(ih) - nmlonp1/3 + nmlonp1 - 1
	  i2 = nmlonp1
       if (ifwr .ge. 1 .or. istep<3) 
     |   write (6,"(1x,'i<1 dusk: i1,2=',2i4)") i1,i2
	  do i=i1,i2
	  vinx(i,k) = 0.
	mltinx(i,k) = mltarr(i)
	  if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	  if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  dmlthalf*2=d(MLT) - find edge ~18MLT
	if (abs(mltinx(i,k)-18.) .lt. dmlthalf) i18=i
	    do j=j1,j2
	      if (phihm(i,j) .lt. vinx(i,k)) then
		vinx(i,k) = phihm(i,j)
		jinx(i,k) = j
              latinx(i,k) = mlatdeg(j)
	      endif
	    enddo   !  j=j1,j2
	  enddo   !  i=i1,i2
      if (ifwr .ge. 2 .or. istep<2) write (6,"(1x,'knx2 i18 =',i4)")i18
        if (ifwr .ge. 2 .or. istep<2)
     |   write (6,"(1x,'knx2 i j v mlt dmlt18 lat =',3i4,4e12.4)") 
     |     (k,i,jinx(i,k),vinx(i,k),mltinx(i,k),abs(mltinx(i,k)-18.),
     |     latinx(i,k),i=i1,i2)
         endif	! if (i_min(ih) - nmlonp1/3 .lt. 1) then
!  Now look at i>nmlonp1 for dusk side:
      if (i_min(ih) + nmlonp1/3 .gt. nmlonp1) then
	  i2 = i_min(ih) + nmlonp1/3 - nmlonp1 + 1
	  i1 = 1
      if (ifwr .ge. 1 .or. istep<3) 
     |    write (6,"(1x,'i>nmlonp1 dusk: i1,2=',2i4)") i1,i2
	  do i=i1,i2
	  vinx(i,k) = 0.
	mltinx(i,k) = mltarr(i)
	  if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	  if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  dmlthalf*2=d(MLT) - find edge ~18MLT
	if (abs(mltinx(i,k)-18.) .lt. dmlthalf) i18=i
	    do j=j1,j2
	      if (phihm(i,j) .lt. vinx(i,k)) then
		vinx(i,k) = phihm(i,j)
		jinx(i,k) = j
                latinx(i,k) = mlatdeg(j)
	      endif
	    enddo   !  j=j1,j2
	  enddo   !  i=i1,i2
      if (ifwr .ge. 2 .or. istep<2) write (6,"(1x,'knx3 i18 =',i4)")i18
        if (ifwr .ge. 2 .or. istep<2)
     |   write (6,"(1x,'knx3 i j v mlt dmlt18 lat =',3i4,4e12.4)") 
     |     (k,i,jinx(i,k),vinx(i,k),mltinx(i,k),abs(mltinx(i,k)-18.),
     |     latinx(i,k),i=i1,i2)
      endif	! if (inx(ih,k) + nmlonp1/3 .gt. nmlonp1)
! Max:
      k = 2
	j1 = j_max(ih) - 4
	if (j1 .lt. 1) j1 = 1		! do not go beyond S Pole at j=1
	j2 = j_max(ih) + 4
	if (j2 .gt. nmlat) j2 = nmlat	! do not go beyond N Pole at j=nmlat
	i1 = i_max(ih) - nmlonp1/3
	if (i1 .lt. 1) i1=1
	i2 = i_max(ih) + nmlonp1/3
	if (i2 .gt. nmlonp1) i2=nmlonp1
! Look at mid-point part
       if (ifwr .ge. 1 .or. istep<3) 
     |   write(6,"(1x,'ih k j1,2 i1,2=',6i3)") ih,k,j1,j2,i1,i2
	do i=i1,i2
	vinx(i,k) = 0.
	mltinx(i,k) = mltarr(i)
	if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  dmlthalf*2=d(MLT) - find edge ~6MLT 
	if (abs(mltinx(i,k)-6.) .lt. dmlthalf) i06=i
	  do j=j1,j2
	    if (phihm(i,j) .gt. vinx(i,k)) then
	      vinx(i,k) = phihm(i,j)
	      jinx(i,k) = j
              latinx(i,k) = mlatdeg(j)
	    endif
	  enddo   !  j=j1,j2
	enddo   !  i=i1,i2
      if (ifwr .ge. 2 .or. istep<2) write (6,"(1x,'knx4 i06 =',i4)")i06
       if (ifwr .ge. 2 .or. istep<2)
     |   write (6,"(1x,'knx4 i j v mlt dmlt06 lat =',3i4,4e12.4)") 
     |     (k,i,jinx(i,k),vinx(i,k),mltinx(i,k),abs(mltinx(i,k)-6.),
     |     latinx(i,k),i=i1,i2)
!  Now look at i<1 for dawn side:
	if (i_max(ih) - nmlonp1/3 .lt. 1) then
	  i1 = i_max(ih) - nmlonp1/3 + nmlonp1 - 1
	  i2 = nmlonp1
       if (ifwr .ge. 1 .or. istep<3) 
     |   write (6,"(1x,'i<1 dawn: i1,2=',2i4)") i1,i2
	  do i=i1,i2
	  vinx(i,k) = 0.
	mltinx(i,k) = mltarr(i)
	  if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	  if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  dmlthalf*2=d(MLT) - find edge ~6MLT
	if (abs(mltinx(i,k)-6.) .lt. dmlthalf) i06=i
	    do j=j1,j2
	      if (phihm(i,j) .gt. vinx(i,k)) then
		vinx(i,k) = phihm(i,j)
		jinx(i,k) = j
              latinx(i,k) = mlatdeg(j)
	      endif
	    enddo   !  j=j1,j2
!  Look at vinx=0 for low values of i (decreasing time - phin)
	  enddo   !  i=i1,i2
      if (ifwr .ge. 2 .or. istep<2) write (6,"(1x,'knx5 i06 =',i4)")i06
        if (ifwr .ge. 2 .or. istep<2)
     |   write (6,"(1x,'knx5 i j v mlt dmlt06 lat =',3i4,4e12.4)") 
     |     (k,i,jinx(i,k),vinx(i,k),mltinx(i,k),abs(mltinx(i,k)-6.),
     |     latinx(i,k),i=i1,i2)
	endif
!  Now look at i>nmlonp1 for dawn side:
      if (i_max(ih) + nmlonp1/3 .gt. nmlonp1) then
	  i2 = i_max(ih) + nmlonp1/3 - nmlonp1 + 1
	  i1 = 1
      if (ifwr .ge. 1 .or. istep<3) 
     |    write (6,"(1x,'i>nmlonp1 dawn: i1,2=',2i4)") i1,i2
	  do i=i1,i2
	  vinx(i,k) = 0.
	mltinx(i,k) = mltarr(i)
	  if (mltinx(i,k) .gt. 24.) mltinx(i,k) = mltinx(i,k) - 24.
	  if (mltinx(i,k) .lt.  0.) mltinx(i,k) = mltinx(i,k) + 24.
!  dmlthalf*2=d(MLT) - find edge ~6MLT
	if (abs(mltinx(i,k)-6.) .lt. dmlthalf) i06=i
	    do j=j1,j2
	      if (phihm(i,j) .gt. vinx(i,k)) then
		vinx(i,k) = phihm(i,j)
		jinx(i,k) = j
              latinx(i,k) = mlatdeg(j)
	      endif
	    enddo   !  j=j1,j2
	  enddo   !  i=i1,i2
      if (ifwr .ge. 2 .or. istep<2) write (6,"(1x,'knx6 i06 =',i4)")i06
       if (ifwr .ge. 2 .or. istep<2)
     |   write (6,"(1x,'knx6 i j v mlt dmlt06 lat =',3i4,4e12.4)") 
     |     (k,i,jinx(i,k),vinx(i,k),mltinx(i,k),abs(mltinx(i,k)-6.),
     |     latinx(i,k),i=i1,i2)
      endif	! if (i_max(ih) + nmlonp1/3 .gt. nmlonp1) 
       if (i06 .lt. 1 .or. i18 .lt. 1) then
	  ifbad(ih) = 2
	  cycle
       endif
       if (ifwr .ge. 1 .or. istep<3) 
     |   write (6,"(1x,'ih lat_i18(1),_i06(2) =',i2,2f6.1)")
     |    ih,latinx(i18,1),latinx(i06,2)
        if (latinx(i06,2) .lt. -990. .or. latinx(i18,1) .lt. -990.) then
          ifbad(ih) = 2
          cycle
       endif
!  Estimate dskofc from lat of peak at 6 and 18 MLT (colat(18-6), lat(6-18))
!	dskof = abs(latinx(i06,2)) - abs(latinx(i18,1))
! Dec 2011 should be HALF this difference
	dskof = (abs(latinx(i06,2)) - abs(latinx(i18,1)))/2.
       if (ifwr .ge. 1 .or. istep<3) 
     |   write (6,"(1x,'ih dskof =',i2,f6.1)") ih,dskof

!  Estimate offc from lat of peak +/-3 hrs (nmlonp1-1)/8 from each maximum
!   (In colat, is nightside-dayside, but in lat is dayside-nightside)
	inm3 = i_min(ih) - (nmlonp1-1)/8
	inp3 = i_min(ih) + (nmlonp1-1)/8
	ixm3 = i_max(ih) - (nmlonp1-1)/8
	ixp3 = i_max(ih) + (nmlonp1-1)/8
	if (inm3 .lt. 1) inm3 = inm3 + nmlonp1 - 1
	if (inp3 .gt. nmlonp1) inp3 = inp3 - nmlonp1 + 1
	if (ixm3 .lt. 1) ixm3 = ixm3 + nmlonp1 - 1
	if (ixp3 .gt. nmlonp1) ixp3 = ixp3 - nmlonp1 + 1
	latnm3 = latinx(inm3,1)
! Feb 2012:  3 hours is too long a time if the cell is small so find the lesser of 
!   3 hr or when latinx is not -999 for cases when the cell is less than 6 hours wide
	if (latnm3 .lt. -990.) then
          do i=inm3+1,i_min(ih)
	    if (latinx(i,1) .gt. -990. .and. latinx(i-1,1) .lt. -990.) 
     |        latnm3 = latinx(i,1)
	  enddo
        endif
        latnp3 = latinx(inp3,1)
	if (latnp3 .lt. -990.) then
          do i=i_min(ih)+1,inp3
	    if (latinx(i-1,1) .gt. -990. .and. latinx(i,1) .lt. -990.) 
     |        latnp3 = latinx(i-1,1)
	  enddo
        endif
        latxm3 = latinx(ixm3,2)
	if (latxm3 .lt. -990.) then
          do i=ixm3+1,i_max(ih)
	    if (latinx(i,2) .gt. -990. .and. latinx(i-1,2) .lt. -990.) 
     |        latxm3 = latinx(i,2)
	  enddo
        endif
	latxp3 = latinx(ixp3,2)
	if (latxp3 .lt. -990.) then
          do i=i_max(ih)+1,ixp3
	    if (latinx(i-1,2) .gt. -990. .and. latinx(i,2) .lt. -990.) 
     |        latxp3 = latinx(i-1,2)
	  enddo
        endif
      if (ifwr .ge. 1 .or. istep<3) 
     |  write (6,"(1x,'ih imin,max inm,p3 ixm,p3 =',7i4)") 
     |    ih,i_min(ih),i_max(ih),inm3,inp3,ixm3,ixp3
       if (latnm3 .lt. -990. .or. latnp3 .lt. -990. .or. latxp3 .lt. 
     |   -990 .or. latxm3 .lt. -990.) then
!  Do not revise when ifbad(ih) already is 2 (should not see 1)
         ifbad(ih) = 3
	 cycle
       endif
! Feb 2012:  offc=0.5*(abs(lat_12MLT)-abs(lat_24MLT), but have to look away from noon-midnight
!   if min/max at 18/06MLT (is not), then +/-3h is 45 deg or 0.7071 in y 
!    so want 0.5*(offcn+offcx)*0.5/0.7071
!   Better if look at 17-19MLT and 7-5MLT y differences, and NOT lat diffs!
	offcn = abs(latnm3) - abs(latnp3)
	offcx = abs(latxp3) - abs(latxm3)
        offcdeg = 0.5*(offcn+offcx)*0.5/0.7071

       if (ifwr .ge. 1 .or. istep<3) 
     | write(6,"(1x,'ih lat,mlt_inm3,inp3,ixp3,ixm3 offcn,x,deg=',i2,
     | 11f6.1)") ih,latnm3,mltinx(inm3,1),latnp3,mltinx(inp3,1),latxp3,
     | mltinx(ixp3,2),latxm3,mltinx(ixm3,2),offcn,offcx,offcdeg

!  Estimate theta0 from 6-18 MLT line first
	crad0 = 90. - 0.5*abs(latinx(i18,1)+latinx(i06,2))
!  Estimate theta0 from 6-18 MLT line in 'convection circle coordinates'
! Transform to convection circle coordinates:
	ofdc = sqrt(offcdeg**2+dskof**2)
! Jan 2012: If ofdc=0 (center on magnetic pole), don't divide by 0 or get NaN!
        asind = 0.
        if (abs(ofdc) .gt. 1.e-5) asind = asin(dskof/ofdc) !FIXME: magic number. Define/use epsilon from cons module?
	sinl18 = sin(abs(latinx(i18,1))*dtr)
	cosl18 = cos(abs(latinx(i18,1))*dtr)
	cosm18 = cos(mltinx(i18,1)*15.*dtr+asind)
	colat18 = cos(ofdc*dtr)*sinl18-sin(ofdc*dtr)*cosl18*cosm18
	colat18 = acos(colat18)*rtd
       if (ifwr .ge. 1 .or. istep<3) 
     |   write (6,"(1x,'ih 18 sinl,cosl,cosm,colat asin=',i2,5e12.4)") 
     |     ih,sinl18,cosl18,cosm18,colat18,asind
!  Feb 2012 bae:  Find rho from ofdc midpoint where rho = sqrt(x*x+y*y)
        y18 = yn + offcdeg
        x18 = xn + dskof
        rho18 = sqrt(x18*x18+y18*y18)
       if (ifwr .ge. 1 .or. istep<3) 
     |  write (6,"(1x,'ih x18 y18 rho18 =',i2,3f6.1)") ih,x18,y18,rho18 
!   theta = atan2(y,x), rho = sqrt(x*x+y*y);  x=rho*cos(theta), y=rho*sin(theta)
!   theta = (MLT-6)*360/24, rho=colat
	sinl06 = sin(abs(latinx(i06,2))*dtr)
	cosl06 = cos(abs(latinx(i06,2))*dtr)
	cosm06 = cos(mltinx(i06,2)*15.*dtr+asind)
	colat06 = cos(ofdc*dtr)*sinl06-sin(ofdc*dtr)*cosl06*cosm06
	colat06 = acos(colat06)*rtd
        if (ifwr .ge. 1 .or. istep<3) 
     |   write (6,"(1x,'ih 06 sinl,cosl,cosm,colat asin=',i2,5e12.4)")
     |     ih,sinl06,cosl06,cosm06,colat06,asind
!  Feb 2012 bae:  Find rho from ofdc midpoint where rho = sqrt(x*x+y*y)
        y06 = yx + offcdeg
        x06 = xx + dskof
        rho06 = sqrt(x06*x06+y06*y06)
       if (ifwr .ge. 1 .or. istep<3) 
     |   write (6,"(1x,'ih x06 y06 rho06 =',i2,3f6.1)")ih,x06,y06,rho06 
	cradcoord = 0.5*(colat18+colat06)
        crad = 0.5*(rho06+rho18)

!  Make sure crad is largest of crad and crad0 (within 0.001 deg in Jan 2012 to avoid printout)
	craduse = max(crad,crad0-0.001)
	if (craduse .gt. crad+0.001) write (6,"(1x,'Used crad0 from 6-18',
     |   ' instead of calc crad =',2e12.4)") crad0,crad
	theta0(ih) = craduse / rtd
       if (ifwr .ge. 1 .or. istep<3) write (6,"(1x,'ih=',i2,' radius:',
     |   ' 0,18,06,c rho18,06,c,a deg=',8f6.1)") ih,crad0,colat18,
     |   colat06,cradcoord,rho18,rho06,crad
!
	offc(ih) = offcdeg / rtd
       if (ifwr .ge. 1 .or. istep<3)
     |  write(6,"(1x,'ih min/max latd3,n3 offc =',i2,8e12.4)") 
     |   ih,latinx(inm3,1),latinx(inp3,1),offcn,latinx(ixp3,2),
     |   latinx(ixm3,2),offcx,offcdeg,offc(ih)

!  oval offset is 2.5 deg towards dawn (more neg dskof)
        dskofc(ih) = dskof / rtd
       if (ifwr .ge. 1 .or. istep<3) 
     |  write (6,"(1x,'ih 18,06 mlt,lat dskof,c,a=',i2,7e12.4)")
     |    ih,mltinx(i18,1),latinx(i18,1),mltinx(i06,2),
     |    latinx(i06,2),dskof,dskofc(ih)

      if (ifwr .ge. 1 .or. istep<3) write (6,
     |  "('Revised convection/oval params hem,off,dsk,n=',
     |  i2,9e12.4)")ih,offc(ih)*rtd,dskofc(ih)*rtd,
     |  theta0(ih)*rtd

      enddo  ! ih=1,2 (and cycle point for 'cycle' of the do loop to continue ih=2 or end)

!  Print out revised convection parameters each timestep
      write (6,"(1x,'calccloc before revcloc (orig ifbad>0): istep ',
     | 'ifbad offcS,N dskofcS,N cradS,N radcenS,N =',i6,2i2,1x,4f7.2,
     | 1x,2f7.2,1x,2f7.2)")
     | istep,ifbad(1),ifbad(2),offc(1)*rtd,offc(2)*rtd,dskofc(1)*rtd,
     |  dskofc(2)*rtd,theta0(1)*rtd,theta0(2)*rtd,radcen(1),radcen(2)

! Revise center and theta0 depending of values of ifbad in each hemisphere
      call revcloc (ifbad,ifwr,istep,avmincrad,radcen,kmlt_min,
     |  kmlt_max,model_ctpoten)

      return
      end subroutine calccloc
!-----------------------------------------------------------------------
      subroutine revcloc (ifbad,ifwr,istep,avmincrad,radcen,kmlt_min,
     |  kmlt_max,model_ctpoten)

!  This subroutine sets defaults for the convection center and radius
!   if the calccloc routine has one or more bad estimates (ifbad not 0).
!  The defaults should be varied according to taste and the high-latitude
!   model used.

!  Aug 2012:  Check if bad convection center or radius in SH and NH.  
!   If both hems bad, can leave the defaults found in aurora_cons for the
!    center or the radius (crad in deg, theta0 in radians used in colath).  
!    crad could be radcen or a min crad of avmincrad from specific model
!  If 1 hemisphere only is bad (ihb) and one is good (ihg), could set
!    ihb to ihg with neg By, or set dskofc=0 and offc to some non-zero
!      number so the auroral coordinates are OK in aurora.F
!  The auroral center (offa,dskofa) and radius (arad) can be revised in aurora.F
!   from values of convection set here depending on flags like CMIT, AMIE, or irdpot
!   phin,phid can use the defaults set in aurora.F that depend on IMF By
!  
! ... Use Association ..................................................
      use aurora_module,only:  ! dimension (2) is for south, north hemispheres
     |     theta0, ! theta0(2), ! convection reversal boundary in radians
     |     offc,   ! offc(2),   ! offset of convection towards 0 MLT relative to mag pole (rad)
     |     dskofc  ! dskofc(2)  ! offset of convection in radians towards 18 MLT (f(By))
      use cons_module,only: rtd
      use input_module,only: ctpoten
      implicit none

! ... Arguments ........................................................ 
      integer,intent(in) :: ifbad(2),ifwr,istep 
      real,intent(in) :: avmincrad 
      real,intent(in) :: radcen(2)	! convection radius in degrees from min/max
      real,intent(in) :: kmlt_min(2),kmlt_max(2) ! min/max MLTs of neg/pos CP for print-out
      real,intent(in) :: model_ctpoten(2)	 ! CP in SH and NH 

! ... Local variables ..................................................  
      integer :: ih,ihb,ihg 
      real :: craduse,offcdeg,dskof

! ... Begin ............................................................

!  NOTE: revcloc assumes that the center is fixed at default_offcdeg for all cases since
!    the dawn-dusk dskofc is usually opposite signs in the hemispheres, and can be large
!   If this is not the case, then revcloc needs to be changed to use the calculated offc,dskofc

! Set ctpoten for output as average of CP in both hemispheres for ifrdpot=1
       ctpoten = 0.5*(model_ctpoten(1)+model_ctpoten(2))

! TEMP to fort.51 like advance.F 
!  NOTE:  With more than 1 processor, writing to fort.51 gives random results 
!  Thus for nproc=4, often model_ctpoten=0, but if this is written to unit 6, print-out is good
        if (istep==1) write (51,"('calccloc: ',
     |   'istep,ifbad,radcen,theta0,offc,dskofc,model_ctpoten'/
     |   '     revcloc: theta0,offc,dskofc'/
     |   ' advance: istep,modeltime(dayn,hr,min,sec),model_hpi,power,',
     |   'power_ave,model_ctpoten,ctpoten,model_qj,qj_ave')")
	write (51,"(1x,i10,2i2,10e12.4)")
     |  istep,ifbad,radcen,theta0*rtd,offc*rtd,dskofc*rtd,model_ctpoten
	write (6,"(1x,i10,2i2,10e12.4,' istep S,N: ifbad radcen theta0',
     |  ' offc dskofc CP')")
     |  istep,ifbad,radcen,theta0*rtd,offc*rtd,dskofc*rtd,model_ctpoten

! Set default center of convection if desired (note, both CANNOT be zero!)
! CANNOT have offc=dskofc=0 or divide by zero in potm in heelis.F
          offcdeg = 1.1		! Default center in deg towards midnight of 1.1 used for Heelis in aurora_cons
	  dskof = 0.

! Use defaults (radcen, offcdeg, dskof) if both ifbad are bad 
!  and check radcen(1,2) to see if they are <10, >30, or >4 deg apart
!   where radcen is the distance between min and max potentials (assuming mpole center)
      if (ifbad(1) .gt. 0 .and. ifbad(2) .gt. 0) then
       write(6,"(1x,'Use defaults, ifbad are both bad,radcen,kmlt,CP:',
     |   2i2,8f7.2)") ifbad,radcen,kmlt_min,kmlt_max,model_ctpoten
       do ih=1,2
          offc(ih) = offcdeg / rtd
          dskofc(ih) = dskof / rtd
       enddo	! ih=1,2
        if (radcen(1)<10. .or. radcen(2)<10. .or. radcen(1)>30. .or.
     |      radcen(2)>30. .or. abs(radcen(2)-radcen(1))>4.) then
!  6/13 bae found SWMF had SH~18deg and NH~28deg for ifbad=1,1 and died shortly thereafter in threed table
!    Set theta0=min(radcen(ih)) in this case 
          write (6,"(1x,'radcen problems for both bad hem at istep =',
     |     i10)") istep
	  craduse = max(avmincrad,min(radcen(1),radcen(2)))
          theta0(1) = craduse / rtd
          theta0(2) = craduse / rtd
        else
! 8/13 bae found CMIT 06348 radcen~10deg 0524, radcen~22deg 0530-0557 w CP~15-25kV
         if (radcen(1)>15. .and. radcen(2)>15. .and. model_ctpoten(1)
     |     <35. .and. model_ctpoten(2)<35.) then
          write (6,"(1x,'radcen big for low CP, set crad=avmincrad ',
     |     f5.2,' at istep=',i10)") avmincrad,istep
          craduse = avmincrad
          theta0(1) = craduse / rtd
          theta0(2) = craduse / rtd
         else
          do ih=1,2
	   craduse = max(avmincrad,radcen(ih))
	   theta0(ih) = craduse / rtd
          enddo	! ih=1,2  
         endif  
        endif
       go to 1000      
      endif

!  If the calculation of the convection location is good in both hemispheres,
!    can use the default center and keep the crad (as desired)
      if (ifbad(1)+ifbad(2) .eq. 0) then
       do ih=1,2
          offc(ih) = offcdeg / rtd
          dskofc(ih) = dskof / rtd
       enddo	! ih=1,2
       go to 1000    
      endif

! Now revise if only 1 hemisphere is bad
       if (ifbad(1)+ifbad(2)>0) then
       if (ifbad(1)>0) then
        ihb = 1
        ihg = 2 
       else
        ihb = 2
        ihg = 1
       endif
       theta0(ihb) = theta0(ihg)
! Set to opposite hem
!       offc(ihb) = offc(ihg)
!       dskofc(ihb) = -dskofc(ihg)
! Or set center to default
         do ih=1,2
	  offc(ih) = offcdeg / rtd
          dskofc(ih) = dskof / rtd
         enddo	! ih=1,2
       go to 1000    
       endif        ! for revisions for 1 bad hemispheres

!  Should not reach this point so stop
      write (6,"(1x,'revcloc: ifbad problem - stop',2i2)")
      stop

1000  continue

      if (ifwr .ge. 1 .or. istep<3) write (6,
     |  "('Revised convection/oval params SH,NH: off,dsk,rad=',
     |    6f8.2)") offc(1)*rtd,dskofc(1)*rtd,theta0(1)*rtd,
     |    offc(2)*rtd,dskofc(2)*rtd,theta0(2)*rtd

! TEMP to fort.51 like advance.F 
!  NOTE:  With more than 1 processor, writing to fort.51 gives random results 
!  Thus for nproc=4, often model_ctpoten=0, but if this is written to unit 6, print-out is good
	write (51,"(15x,24x,6e12.4)") theta0*rtd,offc*rtd,dskofc*rtd
	write (6,"(15x,24x,6e12.4,' revcloc S,N: theta0,offc,dskofc')")
     |    theta0*rtd,offc*rtd,dskofc*rtd

      return

      end subroutine revcloc
!-----------------------------------------------------------------------
