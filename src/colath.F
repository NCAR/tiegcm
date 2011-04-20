!-----------------------------------------------------------------------
      subroutine colath
!
! Calculate pfrac fractional presence of dynamo equation using critical
!  convection colatitudes crit(2).  (crit is in cons module)  All in radians.
!
! 01/11: Move colath from heelis.F to colath.F.  Calculate crit in colath from theta0 
!       found from ctpoten in aurora_cons, where offc=1.1 deg and dskofc = -0.08 deg 
!       (converted to radians).  In Weimer, theta0, offc, and dskofc are recalculated.
!       (For Weimer 2005, offc=4.2 deg, dskofc = 0 deg, and theta0 = crad = bndyfitr/2)
!       AMIE and CISM could: (1) use calccloc in wei01gcm.F to calculate theta0, offc, dskofc,
!        from their potential patterns (where if Bz>0 and Bz>|By|, set theta0=10), 
!        (2) let aurora_cons calculate theta0 from their ctpoten,
!	 (3) use the old default crit(1,2) = 15,30 deg (in radians) from cons.F by not
!        setting ctpoten, but using the NAMELIST CTPOTEN = 13. (or less) and BYIMF = 0 
!        (or do no set BYIMF because if missing, it is set to 0 in input.F) with
!        POTENTIAL_MODEL = 'HEELIS' or 'NONE' (The low CTPOTEN ensures theta0<10 deg, while
!	 byimf=0 ensures dskofa in both hemispheres is identical at -0.08 deg.),
!        (4) set theta0=10 default in aurora_cons as CISM does
! TEST: NONE, CTPOTEN=13., BYIMF=0 or no BYIMF - crit1,2 = 15,30, theta0=9.917, offc=1.1, dskofc=-0.08
! TEST: HEELIS, CTPOTEN=13., no BYIMF - crit1,2 = 15,30, theta0=9.917, offc=1.1, dskofc=-0.08
! TEST bz=+4.4, By=2.0, Vsw=484, 02080: W01 theta0=10 (if Bz>0,Bz>|By|), W05 theta0=13.36 deg
!
      use dynamo_module,only: nmlat0,pfrac
      use params_module,only: nmlonp1
      use aurora_module,only: theta0, dskofc, offc ! see aurora.F
      use magfield_module,only: sunlons 
      use cons_module,only: rtd,
     |  crit,        ! critical colatitudes crit(2)
     |  ylonm,ylatm  ! magnetic grid lons, lats
      implicit none
!
! Local:
      integer :: i,j
      real :: sinlat,coslat,aslonc,ofdc,cosofc,sinofc,crit1deg
      real,dimension(nmlonp1,nmlat0) :: colatc

! 01/11 bae:  Revise crit in rad so crit(1)=theta0(=crad in rad)+5deg, crit(2)=crit(1)+15deg
        crit1deg = max(15.,0.5*(theta0(1)+theta0(2))*rtd + 5.)
	crit1deg = min(30.,crit1deg) 
! To make the same as in cons.F, comment out next line
 	crit(1) = crit1deg/rtd
	crit(2) = crit(1) + 15./rtd
!
! offc(2), dskofc(2) are for northern hemisphere aurora (see aurora.F)
! 01/11 bae:  Revised so that use average of both hemispheres instead of NH only
!     ofdc = sqrt(offc(2)**2+dskofc(2)**2)
      ofdc = sqrt((0.5*(offc(1)+offc(2)))**2 +
     |  (0.5*(dskofc(1)+dskofc(2)))**2)
      cosofc = cos(ofdc)
      sinofc = sin(ofdc)
      aslonc = asin(0.5*(dskofc(1)+dskofc(2))/ofdc)
!     aslonc = asin(dskofc(2)/ofdc)
! TEMP
!      write (6,"(1x,'COLATH: crit1,2 theta0 offc dskofc =',8e12.4)")
!    |   crit(1)*rtd,crit(2)*rtd,theta0(1)*rtd,theta0(2)*rtd,
!    |   offc(1)*rtd,offc(2)*rtd,dskofc(1)*rtd,dskofc(2)*rtd
!
! Define colatc with northern convection circle coordinates
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc)
!
      do j=1,nmlat0
        sinlat = sin(abs(ylatm(j+nmlat0-1)))
        coslat = cos(    ylatm(j+nmlat0-1))
        do i=1,nmlonp1
          colatc(i,j) = cos(ylonm(i)-sunlons(1)+aslonc)
          colatc(i,j) = acos(cosofc*sinlat-sinofc*coslat*colatc(i,j))
        enddo ! i=1,nmlonp1

!       write(6,"('colath: j=',i3,' colatc(:,j)=',/,(6e12.4))")
!    |    j,colatc(:,j)*rtd

!
! Calculate fractional presence of dynamo equation at each northern
! hemisphere geomagnetic grid point. Output in pfrac(nmlonp1,nmlat0)
!
	do i=1,nmlonp1
          pfrac(i,j) = (colatc(i,j)-crit(1))/(crit(2)-crit(1))
          if (pfrac(i,j) < 0.) pfrac(i,j) = 0.
          if (pfrac(i,j) >= 1.) pfrac(i,j) = 1.
	enddo ! i=1,nmlonp1

!       write(6,"('colath: j=',i3,' pfrac(:,j)=',/,(6e12.4))")
!    |    j,pfrac(:,j)

      enddo ! j=1,nmlat0

      end subroutine colath
!-----------------------------------------------------------------------
