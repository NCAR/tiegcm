      subroutine colath
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate pfrac fractional presence of dynamo equation using critical
!  convection colatitudes crit(2).  (crit is in cons module)
!
      use aurora_module,only: theta0, dskofc, offc ! see aurora.F
      use magfield_module,only: sunlons 
      use cons_module,only: rtd,
     |  crit,        ! critical colatitudes crit(2)
     |  ylonm,ylatm  ! magnetic grid lons, lats
      use pdynamo_module,only: nmlat0,pfrac ! pfrac is output
      use params_module,only: nmlonp1,rp
      use mage_coupling_module,only: aurllbN,aurllbS ! colat in degrees
      
      implicit none
!
! Local:
      integer :: i,j
      real(rp) :: sinlat,coslat,aslonc,ofdc,cosofc,sinofc,crit1deg
      real(rp),dimension(nmlonp1,nmlat0) :: colatc

#ifdef GAMERA
! Lin, 09/2021: Use dynamic auroral boundary defined by the lowest lat where numflux<1e6
        !crit(1) = (max(aurllbN,aurllbS)+5._rp)/rtd
! Lin, 10/2021: using 5 deg wide merging zone tends to
! produce artificial eastward wind and drift bands at mid latitudes
        !crit(2) = crit(1)+10._rp/rtd
        crit1deg = max(aurllbN,aurllbS)
        ! at low colat, the merge zone is 10deg wide and is offset 5 deg from aurbc
        ! at high colat, the merge zone is 5 deg wide and starts at aurbc
        ! starts at 15 degree
        !crit(1) = crit1deg+max(7.5_rp+crit1deg*(-5._rp/30._rp),-5._rp)
        crit(1) = crit1deg + 5._rp
        !crit(1) = crit1deg
        crit(1) = max(crit(1),15._rp)/rtd
        ! starts at 25 degree
        !crit(2) = crit1deg+max(20._rp+crit1deg*(-10._rp/30._rp),5._rp)
        !crit(2) = max(crit(2),25._rp)/rtd
        crit(2) = crit(1) + 10._rp/rtd
!        crit1deg = max(aurllbN,aurllbS)
!        if (crit1deg .le. 45._rp) then
!            crit(1) = crit1deg+max(7.5_rp+
!     |                crit1deg*(-5._rp/30._rp),-5._rp)
!            crit(2) = crit(1) + 10._rp/rtd
!        else
!            crit(1) = 45._rp/rtd
!            crit(2) = crit(1) + (crit1deg*0.5_rp-12.5_rp)/rtd
!        endif
#else
! 01/11 bae:  Revise crit in rad so crit(1)=theta0(=crad in rad)+5deg, crit(2)=crit(1)+15deg
        crit1deg = max(15._rp,0.5_rp*(theta0(1)+theta0(2))*rtd + 5._rp)
        crit1deg = min(30._rp,crit1deg)
! To make the same as in cons.F, comment out next line
        crit(1) = crit1deg/rtd
        crit(2) = crit(1) + 15._rp/rtd
#endif
!
! offc(2), dskofc(2) are for northern hemisphere aurora (see aurora.F)
! 01/11 bae:  Revised so that use average of both hemispheres instead of NH only
!     ofdc = sqrt(offc(2)**2+dskofc(2)**2)
      ofdc = sqrt((0.5_rp*(offc(1)+offc(2)))**2 +
     |  (0.5_rp*(dskofc(1)+dskofc(2)))**2)
      cosofc = cos(ofdc)
      sinofc = sin(ofdc)
      aslonc = asin(0.5_rp*(dskofc(1)+dskofc(2))/ofdc)
! TEMP
!       write (6,"(1x,'COLATH: crit1,2 dskofc offc deg=',6e12.4)")
!    |   crit(1)*rtd,crit(2)*rtd,dskofc(1)*rtd,offc(1)*rtd,
!    |   dskofc(2)*rtd,offc(2)*rtd
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
          if (pfrac(i,j) < 0._rp) pfrac(i,j) = 0._rp
          if (pfrac(i,j) >= 1._rp) pfrac(i,j) = 1._rp
	enddo ! i=1,nmlonp1

!       write(6,"('colath: j=',i3,' colatc(:,j)=',/,(6e12.4))")
!    |    j,colatc(:,j)
!       write(6,"('colath: j=',i3,' pfrac(:,j)=',/,(6e12.4))")
!    |    j,pfrac(:,j)

      enddo ! j=1,nmlat0
      end subroutine colath
