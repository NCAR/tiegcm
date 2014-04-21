!      
      module weimer_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!      
! module used to calculate the Weimer 2001 model potential in both
!  hemispheres from Bz, By, Sws, Swd, and AL.
!
      use params_module,only: nmlat,nmlon,nmlonp1
      use dynamo_module,only: nmlat0,phihm
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
! 01/10 bae:  Have ctpoten from both hemispheres from Weimer 2001
      real :: weictpoten(2)

      contains 

!-----------------------------------------------------------------------
!
      subroutine weimer01

! 01/10:  Added ctpoten to store average weictpoten from calccloc
      use input_module,only: ctpoten
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
!           convection colatitudes crit(2) (15,30 in cons.F) re-calc in colath
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

!  Re-calculate SH values of offa, dskofa, arad, and Heelis phid and phin from
!   Weimer 2001 boundarylat values of theta0 and 2005 setboundary offc, dskofc
	call wei01loc (1)

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

!  Re-calculate NH values of offa, dskofa, arad, and Heelis phid and phin from
!   Weimer 2001 boundarylat values of theta0 and 2005 setboundary offc, dskofc
	call wei01loc (2)

!
! Periodic points:
      do j=1,nmlat
        phihm(nmlonp1,j) = phihm(1,j)
      enddo ! j=1,nmlat

!  Re-calculate the values of theta0,offc,dskofc,arad,phid,phin,offa,dskofa
! 01/11 bae:  Comment out call to calccloc (although it works) since wei01loc is used instead
!      call calccloc (weictpoten)
! 01/10 bae:  Return the average of the weictpoten from the SH and NH in ctpoten
      ctpoten = 0.5*(weictpoten(1)+weictpoten(2))
!      write (6,"(1x,'weimer01: CPS,N ctpoten=',3f8.1)") weictpoten,
!     |    ctpoten
!
! Calculate pfrac fractional presence of dynamo equation using critical
!  convection colatitudes crit(2).  (crit is in cons module, re-calc in colath)
!
      call colath
!
! From dynamo.F:
!     real,dimension(nmlonp1,nmlat0) ::
!    |   pfrac, ! NH fraction of potential
!    |   phihm  ! potential in magnetic
!
      end subroutine weimer01
!-----------------------------------------------------------------------
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
      use input_module,only: mxlen_filename
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
      character(len=mxlen_filename) :: dskfile

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
      IF (IYYY.EQ.0) call shutdown('There is no Year Zero.')
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
      subroutine wei01loc (ih)
! ih=1,2 for SH,NH called from weimer01

!  Jan 2012 bae: Made wei01loc from similar wei05loc (crad=bndyfitr/2.)
!    where 18 of 24 calls to BoundaryLat(MLT) for MLT=15.5-8.5 are averaged
!    into avblat to find crad=(90.-avblat)/2.
!   Added byloc to limit size of byimf (probably incorrect assymmetry in hems)
!  NCAR July 08 bae:  Write wei05loc to use values from Weimer 2005
!       setboundary for convection offc=4.2 deg, dskofc=0, theta0=bndyfitr
!       where deg should be in radians.  wei05loc replaces calccloc05.
!
! (dimension 2 is for south, north hemispheres)
!  Calculate offa, dskofa, rrad, phid, and phin from Weimer 2005 offc, dskofc, theta0
!   Use Fig 8 of Heelis et al. [JGR, 85, 3315-3324, 1980]
!     This shows:  arad = 18.7 deg, crad = 16.7 deg (so arad = crad + 2 deg)
!           offa = offc = 3 deg (so offa = offc)
!           dskofc = 2 deg, dskofa = -0.5 deg  (so dskofa = dskofc - 2.5 deg)
!   Parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                             and phin (phin(MLT)=23.50 +/- 0.15By - 12)
!   (In aurora_cons, phid=0., phin=180.*rtd)
!     (For zero By, should be phid=21.39MLT*15*rtd, phin=11.5*15*rtd)
!  05/08:  But formulae for ra-rc using IMF CP between 5-7 deg (not 2) so use
!          difference of ra(max IMF CP or HP)-rc(IMF CP) as in aurora.F
! These are the dimensions and descriptions (corrected phid,n) from aurora.F:
!    |  theta0(2), ! convection reversal boundary in radians
!    |  offa(2),   ! offset of oval towards 0 MLT relative to magnetic pole (rad)
!    |  dskofa(2), ! offset of oval in radians towards 18 MLT (f(By))
!    |  phid(2),   ! dayside convection entrance in MLT-12 converted to radians (f(By)),
!                      phid is the MLT-12 location of the cusp on the dayside
!    |  phin(2),   ! night convection entrance in MLT-12 converted to radians (f(By))
!    |  rrad(2),   ! radius of auroral circle in radians
!    |  offc(2),   ! offset of convection towards 0 MLT relative to mag pole (rad)
!    |  dskofc(2)  ! offset of convection in radians towards 18 MLT (f(By))
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc)
!
      use dynamo_module,only: nmlat0,phihm
! Additional auroral parameters (see sub aurora_cons):
      use aurora_module,only: theta0,offa,dskofa,phid,phin,rrad,offc,
     |                        dskofc
      use magfield_module,only: sunlons 
      use input_module,only:  ! from user input
     |  byimf    ! By component of IMF (nT)      (e.g., 0.)
! 05/08:  Added bzimf,ctpoten,power
     | ,bzimf    ! Bz component of IMF (nT)      (e.g., 0., not used except in print)
     | ,ctpoten  ! cross-cap potential (kV)      (e.g., 45.)
     | ,power    ! hemispheric power   (GW)      (e.g., 16.)
      use cons_module,only: rtd,
     |  ylonm,ylatm  ! magnetic grid lons, lats in radians
      implicit none
!
! Args:
      integer,intent(in) :: ih
!
! Local:
! 05/08:  Added rccp etc
      real :: rccp,racp,rahp,ramx,diffrac,plevel,tmltmin,tmltmax
      real :: diffrac2
      real :: offcdegp(2),dskofcp(2),radegp(2)
      integer :: i,i1,i2,j,j1,j2
      real :: vnx(2,2),hem,mltd,mltn
      integer :: jnx(2,2),inx(2,2)
      real :: offcdeg,dskof,arad,crad
      real :: cosl06,sinl06,colat06,cosm06,cosl18,sinl18,colat18,cosm18
! 01/11:  Added blatmlt(24), ml and avblat for wei01loc (different from wei05loc)
      integer :: ml
      real :: blatmlt(24),avblat
      real :: byloc

! Limit size of byimf in phin and phid calculations (as in aurora.F) 
!  NOTE:  This byloc is assymetric in hemisphere, which is probably not correct
      byloc = byimf
      if (byloc .gt. 7.) byloc = 7.
      if (byloc .lt. -11.) byloc = -11.
!
!  ih=1 is SH, ih=2 is NH
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
! 05/08: Calculate weictpoten in kV from Weimer model min/max in V
	weictpoten(ih) = 0.001 * (vnx(ih,2) - vnx(ih,1))
	tmltmin = (ylonm(inx(ih,1))-sunlons(1)) * rtd/15. + 12.
	if (tmltmin .gt. 24.) tmltmin = tmltmin - 24.
	tmltmax = (ylonm(inx(ih,2))-sunlons(1)) * rtd/15. + 12.
	if (tmltmax .gt. 24.) tmltmax = tmltmax - 24.
!       write (6,"('ih Bz By Hp ctpoten,wei min/max potV,lat,mlt=',i2,
!    |    5f8.2,2x,e12.4,2f8.2,2x,e12.4,2f8.2))") ih,bzimf,byimf,power,
!    |    ctpoten,weictpoten(ih),
!    |    vnx(ih,1),ylatm(jnx(ih,1))*rtd,tmltmin,
!    |    vnx(ih,2),ylatm(jnx(ih,2))*rtd,tmltmax
! 05/08: From aurora_cons, calculate convection and aurora radii using IMF convection
!   and power (plevel);  racp (DMSP/NOAA) - rccp (AMIE) = 5.32 (Bz>0) to 6.62 (Bz<0) deg
!  Heelis et al [1980, JGR, 85, pp 3315-3324] Fig 8: ra=rc+2deg, and is 2.5 deg to dusk
	rccp = -3.80+8.48*(weictpoten(ih)**0.1875)
	racp = -0.43+9.69*(weictpoten(ih)**0.1875)
	plevel = 0.
	if (power >=1.00) plevel = 2.09*alog(power)
	rahp = 14.20 + 0.96*plevel
	ramx = max(racp,rahp)
	diffrac = ramx - rccp

!  Set default values
!  Use parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                             and phin (phin(MLT)=23.50 +/- 0.15By - 12)
	mltd = 9.39 - hem*0.21*byloc
	mltn = 23.50 - hem*0.15*byloc
	phid(ih) = (mltd-12.) * 15. / rtd
	phin(ih) = (mltn-12.) * 15. / rtd
! 05/18/08:  Note that phid,phin are only for Heelis and are irrelevant for Weimer
!       write (6,"(1x,'mltd mltn phid,n =',4f8.2)")
!    |   mltd,mltn,phid(ih)*rtd/15.,phin(ih)*rtd/15.
!  Use default constant value of offcdegp from setboundary in Weimer 2005
	offcdeg = 4.2
	offcdegp(ih) = offcdeg
	offc(ih) = offcdegp(ih) / rtd
	offa(ih) = offcdegp(ih) / rtd
!       write (6,"(1x,'offcdeg,rad =',2e12.4)") offcdeg,offc(ih)
	dskof = 0.
	dskofcp(ih) = dskof
	dskofc(ih) = dskof / rtd
!  oval offset is 2.5 deg towards dawn (more neg dskof)
	dskofa(ih) = (dskof-2.5) / rtd
!       write (6,"(1x,'dskof,c,a=',3f8.2)")
!    |    dskof,dskofc(ih)*rtd,dskofa(ih)*rtd
! Set crad from bndyfitr/2 of setboundary of Weimer 2005
!	crad = bndyfitr/2.
        avblat = 0.
	do ml=1,24
         blatmlt(ml) = BoundaryLat(ml-0.5)
! average 9x from 0.5 to 8.5 MLT and 9x from 15.5 to 23.5 MLT (skip 6x 9.5 to 14.5 MLT)
	 if (ml.lt.10 .or. ml.gt.15) avblat = avblat + blatmlt(ml)
        enddo 
! Convert from average boundary latitude to crad       
	avblat = avblat/18.
        crad = (90. - avblat)/2.
!       write (6,"(1x,'wei01loc: ih,bz,y,avblat,crad =',i2,4f8.2)") 
!    |   ih,bzimf,byimf,avblat,crad
!       write (6,"(1x,'blatmlt avblat crad =',26f8.2)") blatmlt,avblat,crad
!  Fig 8 Heelis et al [1980]: ra=rc+2deg, and shifted 2.5 deg to dusk
	arad = crad + 2.
! 05/08:  Make ra=rc+diffrac(=ramx-rccp) - same difference as in aurora.F
! Choose to have arad=crad(Weimer) + diffrac(same diff as in aurora.F)
	arad = crad + diffrac
! 08/08: OR make ra=ramx=max(racp,rahp) so diffrac=arad-crad
	diffrac2 = ramx - crad
! Choose to have arad=ramx (same as in aurora.F as determined by P/CP)
!       arad = ramx
	theta0(ih) = crad / rtd
	radegp(ih) = arad
	rrad(ih) = arad / rtd
!       write (6,"(1x,'radius: crad,rccp,racp,rahp diffa-c',
!    |   '(aurF,ramx-Weic) ramx,Weic+d,arad deg=',9f8.2)") crad,rccp,
!    |   racp,rahp,diffrac,diffrac2,ramx,crad+diffrac,arad

! Print out revised values (revised 05/08):
!       write (6,"(1x,'Revised convection/oval params (off,dsk,',
!    |    'rad,phid,n=',8f9.4)")offc(ih)*rtd,offa(ih)*rtd,
!    |    dskofc(ih)*rtd,dskofa(ih)*rtd,theta0(ih)*rtd,rrad(ih)*rtd,
!    |    phid(ih)*rtd/15.+12.,phin(ih)*rtd/15.+12.

      return
      end subroutine wei01loc
!-----------------------------------------------------------------------
      end module weimer_module
