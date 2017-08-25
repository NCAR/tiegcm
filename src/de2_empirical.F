
      module de2_empirical_module
      
      use params_module,only: nmlon,nmlonp1,nmlat
      use heelis_module,only: phihm
      use pdynamo_module,only: Ed1SD_mag,Ed2SD_mag,pfrac
      use cons_module,only: crit, pi, dtr
      use input_module,only: ctpoten
      use magfield_module,only: sunlons 
      
      implicit none
      
      logical,save :: first_time
      logical,save :: first_time_SD
      
      contains
!--------------------------------------------------------------------------      
      subroutine de2_empirical_poten
         use cons_module, only: ylonm, ylatm,re,rtd
         use params_module, only: nmlat, nmlon, nmlonp1
	 use hist_module, only: modeltime
	 use input_module, only: byimf, bzimf,efieldvar
         use addfld_module

         implicit none
	 
         real :: poten_calc_model_value
	 real, dimension(nmlonp1) :: lons
	 real, dimension(nmlat)   :: lats
	 real :: mlat,lon_deg, MLT, Btrans, imf_angleNorth, 
     &           imf_angleSouth,sinTNorth, sinTSouth,
     &           value_poynt, value_poten,Ed1_value,Ed2_value,
     &           Ed1_glb,Ed2_glb
         real :: get_mlt,get_sinT,valx,valm,de2_ctpoten(2),
     |     value_Sbg
         real :: Stot(nmlonp1,nmlat),Sbg(nmlonp1,nmlat),
     |     ED1_mag(nmlonp1,nmlat),ED2_mag(nmlonp1,nmlat)
	 integer :: i, j, nmlat0
	 
	 
         sinTNorth = get_sinT(0)
         sinTSouth = get_sinT(183)   ! I.e.: pretend that we are 183 days later
	 Btrans = sqrt(byimf ** 2 + bzimf ** 2)
	 imf_angleNorth = atan2(byimf, bzimf)  ! / 3.1416 * 180.
         if (imf_angleNorth .lt. 0) imf_angleNorth = imf_angleNorth + 
     |         2*3.1416 ! + 360
	 imf_angleSouth = atan2(-byimf, bzimf)  ! / 3.1416 * 180.
         if (imf_angleSouth .lt. 0) imf_angleSouth = imf_angleSouth + 
     |         2*3.1416
	 nmlat0=(nmlat+1)/2
!	 write(6,*) 'in_de2_conditions',sinTNorth,sinTSouth,Btrans,
!     |     imf_angleNorth,imf_angleSouth
	 
	 if(first_time) first_time_SD = .true.
! loop over all points	 	 
	 do i=1, nmlonp1            
	    do j=1, nmlat
	       mlat    = ylatm(j)*rtd  ! from radian to degree
	       lon_deg = ylonm(i)*rtd  ! from radian to degree
               !mlt     = get_mlt(lon_deg)
               mlt = (ylonm(i)-sunlons(1))*rtd/15.+12.
               !if(i.eq.1.and.j.eq.1) then
	       !  call calc_pot_S(16.77,75.77,0.,imf_angleNorth,
     &         !     5.,first_time,value_poynt,value_poten) 
               !   write(6,*) 'de_specific',value_poten
	       !  stop
               ! endif
	       if (mlat.gt.45.) then
                 call calc_pot_S(mlt,mlat,sinTNorth,imf_angleNorth,
     &            Btrans,first_time,value_poynt,value_poten,value_Sbg,
     &            Ed1_glb,Ed2_glb) 
     
                  phihm(i, j) = value_poten !* 1e-3
                  Stot(i, j)  = value_poynt
                  Sbg(i, j)   = value_Sbg
                  ED1_mag(i, j) = Ed1_glb
                  ED2_mag(i, j) = Ed2_glb
                
		    if (efieldvar .ne. 0) then
		      call calc_SDEd(mlt,mlat,sinTNorth,
     &  		imf_angleNorth,Btrans,first_time_SD,Ed1_value,
     &  		Ed2_value)
! change the unit from mV/m to v/rad, since emx,emy are under v/rad ====
! dS_lon = R * cos(lat).dlon 
! dS_lat = R * dlat 	  
		      Ed1SD_mag(i,j) = Ed1_value*1e-3  !*1e-3*re*1e-2 ! *
!     &  	 cos(mlat *dtr)*1.
		      Ed2SD_mag(i,j) = Ed2_value*1e-3  !*1e-3*re *1e-2*1. 
		     first_time_SD = .false.
		  endif

               else if (mlat.lt.-50) then    ! Notice the -mlat
                 call  calc_pot_S(mlt,-mlat,sinTSouth,imf_angleSouth,
     &            Btrans,first_time,value_poynt,value_poten,value_Sbg,
     &            Ed1_glb,Ed2_glb) 
     
                  phihm(i, j) = value_poten !* 1e-3
                  Stot(i, j)  = value_poynt
                  Sbg(i, j)   = value_Sbg
                  ED1_mag(i, j) = Ed1_glb
                  ED2_mag(i, j) = Ed2_glb
		  
		  if (efieldvar .ne. 0) then
		     call calc_SDEd(mlt,-mlat,sinTSouth,
     &  		imf_angleSouth,Btrans,first_time_SD,Ed1_value,
     &  		Ed2_value)
! change the unit from mV/m to v/rad, since emx,emy are under v/rad ====
! dS_lon = R * cos(lat).dlon 
! dS_lat = R * dlat 
! convert from mV/m to V/m
		     Ed1SD_mag(i,j) = Ed1_value*1e-3 !*re*1e-2  ! *
!     &  	 cos(mlat *dtr)*1.
		     Ed2SD_mag(i,j) = Ed2_value*1e-3  !*re *1e-2*1.
		     first_time_SD = .false.
		  endif
!
               else
	          phihm(i, j) = 0
	       endif
	       !write(6,'(a10,3(x,f10.5))') 'de2_pts ',mlt,mlat,phihm(i, j)*1e-3
               
	    end do
	 end do
	 
      valx = maxval(phihm(:,1:nmlath-1))	 
      valm = minval(phihm(:,1:nmlath-1))	
      de2_ctpoten(1) = valx-valm
!      write(6,*) 'de_ctpoten_sh',valx,valm,de2_ctpoten(1)
      valx = maxval(phihm(:,nmlath:))	 
      valm = minval(phihm(:,nmlath:))	
      de2_ctpoten(2) = valx-valm
!      write(6,*) 'de_ctpoten_nh',valx,valm,de2_ctpoten(2)
      ctpoten = 0.5*(de2_ctpoten(2)+de2_ctpoten(1))
      
      call de2loc(1)  ! test to see if I get an aurora
      call de2loc(2)
         
      call addfld('potde2','potde2','V',phihm,'mlon',1,nmlonp1,
     |  'mlat',1,nmlat,0)
      call addfld('SDEd1','SDEd1','V/m',Ed1SD_mag,'mlon',1,nmlonp1,
     |  'mlat',1,nmlat,0)
      call addfld('SDEd2','SDEd2','V/m',Ed2SD_mag,'mlon',1,nmlonp1,
     |  'mlat',1,nmlat,0)
      call addfld('Ed1DE2','Ed1DE2','V/m',Ed1_mag,'mlon',1,nmlonp1,
     |  'mlat',1,nmlat,0)
      call addfld('Ed2DE2','Ed2DE2','V/m',Ed2_mag,'mlon',1,nmlonp1,
     |  'mlat',1,nmlat,0)
      call addfld('PoynTot','PoynTot','W/m2',Stot,'mlon',1,nmlonp1,
     |  'mlat',1,nmlat,0)
      call addfld('PoynBg','PoynBg','W/m2',Sbg,'mlon',1,nmlonp1,
     |  'mlat',1,nmlat,0)
      
      end subroutine de2_empirical_poten
!--------------------------------------------------------------------------       
      end module
!--------------------------------------------------------------------------  
!--------------------------------------------------------------------------------
! Inspired pretty directly from Christine's routine get_sinT in fit_data_btrans
! See the comment on dayShift in get_subsol
      real function get_sinT(dayShift)
         use apex, only: colat,elon, vp,dtr
	 
	 implicit none 
	 
	 real, parameter :: pi = 3.14159
	 real :: subsol_colat, subsol_lon
	 real :: theta_prad, phi_prad, subsol_colatrad, subsol_lonrad
         integer :: dayShift

         call get_subsol(subsol_colat, subsol_lon, dayShift)

         ! Change all angles from degrees to radians:

         theta_prad = dtr*colat ! COLAT
	 phi_prad   = dtr*elon ! ELON
	 subsol_colatrad = dtr*subsol_colat ! SBSLLAT
	 subsol_lonrad   = dtr*subsol_lon ! SBSLLON

         ! Calculate dipole tilt angle:

	 get_sinT = cos(theta_prad)*sin(subsol_colatrad) 
     &   	 + sin(theta_prad)*cos(subsol_colatrad)*
     &             cos(phi_prad-subsol_lonrad)
     
      end function get_sinT 
!--------------------------------------------------------------------------------     
! Both parameters are out, and in degrees
! get_subsol returns the position of the subsolar point, for the time defined in
! hist_module::modeltime
! dayShift was added because when it's summer in the northern hemisphere,
! it's winter in the other one, and the data is fitted for the 
! northern hemisphere. So, for the North, dayShift = 0. For the South, = 183
      subroutine get_subsol(colat_out, lon_out, dayShift)
         use hist_module, only: modeltime
         use apex,only: subsol
         use init_module,only: iyear

	 implicit none 
	 
         integer :: dayShift
	 
      	 integer :: day, hour, minute
         real :: colat_out, lon_out, second
	 
	 day = modulo(modeltime(1) + dayShift, 365)
         hour = modeltime(2)
	 minute = modeltime(3)
	 second = modeltime(4)
         call subsol(iyear,day,hour,minute,second,
     &             colat_out, lon_out)
	 
      end subroutine get_subsol
!---------------------------------------------------------------------------
      real function get_sinIm(mlat_in) 
      use cons_module, only: pi
      
      implicit none 
	 
      real :: mlat_in

      get_sinIm = sin(mlat_in*pi/180.)/(1-0.75*cos(mlat_in*pi/180.))**0.5

      end function get_sinIm
!---------------------------------------------------------------------------
      real function get_beta(bt_in) 
      implicit none 
      real :: bt_in
      real, parameter :: beta_inf = 16.

      get_beta = bt_in /sqrt(1.+(bt_in*bt_in)/(beta_inf*beta_inf))

      end function get_beta
!-------------------------------------------------------------------------------- 

!-----------------------------------------------------------------------
      subroutine de2loc (ih)
! ih=1,2 for SH,NH called from weimer05
! These are the dimensions and descriptions (corrected phid,n) from aurora.F:
!    |  theta0(2), ! convection reversal boundary in radians
!    |  offa(2),   ! offset of oval towards 0 MLT relative to magnetic pole (rad)
!    |  dskofa(2), ! offset of oval in radians towards 18 MLT (f(By))
!    |  phid(2),   ! dayside convection entrance in MLT-12 converted to radians (f(By))
!                      phid is the MLT-12 location of the cusp on the dayside
!    |  phin(2),   ! night convection entrance in MLT-12 converted to radians (f(By))
!    |  rrad(2),   ! radius of auroral circle in radians
!    |  offc(2),   ! offset of convection towards 0 MLT relative to mag pole (rad)
!    |  dskofc(2)  ! offset of convection in radians towards 18 MLT (f(By))
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc)
!
      use params_module,only: nmlat,nmlon,nmlonp1
      use pdynamo_module,only: nmlat0,phihm
! Additional auroral parameters (see sub aurora_cons):
      use aurora_module,only: theta0,offa,dskofa,phid,phin,rrad,offc,
     |                        dskofc
      use magfield_module,only: sunlons 
      use input_module,only:  ! from user input
     |  byimf,    ! By component of IMF (nT)      (e.g., 0.)
! 05/08:  Added bzimf,ctpoten,power
     | bzimf,    ! Bz component of IMF (nT)      (e.g., 0.)
     | swvel,    ! solar wind velocity (km/s)    (e.g., 400.)
     | swden,    ! solar wind density (#/cm3)    (e.g., 4.)
     | ctpoten,  ! cross-cap potential (kV)      (e.g., 45.)
     | power    ! hemispheric power   (GW)      (e.g., 16.)
      use cons_module,only: rtd,
     |  ylonm,ylatm,  ! magnetic grid lons, lats in radians
     |  crit
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
      real :: byloc,de2ctpoten(2)
      real :: angle,tilt,bt,bndyfitr,get_sinT
!

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
	
! from Weimer05
        tilt  = get_sinT(0)
	bt    = sqrt(byimf**2 + bzimf**2)
	angle = atan2(byimf, bzimf)
        if (angle.lt.0) angle=angle+ 2*3.1416 
        tilt  = hem * tilt
        angle = hem * angle
	write(6,*) 'de2loc_b',angle,bt,tilt,swvel,swden,hem
        call setboundary_de2(angle,bt,tilt,swvel,swden,bndyfitr)
	write(6,*) 'de2loc_a',bndyfitr
	
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
	de2ctpoten(ih) = 0.001 * (vnx(ih,2) - vnx(ih,1))
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
	rccp =  -2.43+9.69*(de2ctpoten(ih)**0.1875)    ! -3.80+8.48*(de2ctpoten(ih)**0.1875)
	racp = rccp ! -0.43+9.69*(de2ctpoten(ih)**0.1875)
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
	offcdeg = 0.01 ! 1.1 ! 1.1 ! -4.2 ! 4.2 ! 0.05  ! 4.2 am_change 8/2017
	offcdegp(ih) = offcdeg
	offc(ih) = offcdegp(ih) / rtd
	offa(ih) = offcdegp(ih) / rtd
!       write (6,"(1x,'offcdeg,rad =',2e12.4)") offcdeg,offc(ih)
	dskof = 0.      ! 4.2 am_change 8/2017
	dskofcp(ih) = dskof
	dskofc(ih) = dskof / rtd
!  oval offset is 2.5 deg towards dawn (more neg dskof)
	dskofa(ih) = (dskof-8.) / rtd
!	dskofa(ih) = (dskof-2.5) / rtd
!	dskofa(ih) =dskof-(6.5*abs(byloc)/7)  ! 4.2 am_change 8/2017
!	dskofa(ih) = (dskofa(ih)) / rtd  ! am_change  8/2017
!       write (6,"(1x,'dskof,c,a=',3f8.2)")
!    |    dskof,dskofc(ih)*rtd,dskofa(ih)*rtd
! Set crad from bndyfitr/2 of setboundary of Weimer 2005
        crad = -2.43+9.69*(de2ctpoten(ih)**0.1875)  ! from Heelis
!	crad = bndyfitr/2.
!	crad = crit(1)*rtd
!      write (6,"(1x,'wei05loc: ih,bz,y,crad =',i2,3f8.2)") 
!    |    ih,bzimf,byimf,crad
!  Fig 8 Heelis et al [1980]: ra=rc+2deg, and shifted 2.5 deg to dusk
!	arad = crad + 2.
! 05/08:  Make ra=rc+diffrac(=ramx-rccp) - same difference as in aurora.F
! Choose to have arad=crad(Weimer) + diffrac(same diff as in aurora.F)
	arad = crad + 2. ! + diffrac
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
       write (6,"(1x,'Revised convection/oval params (off,dsk,',
     |    'rad,phid,n=',8f9.4)")offc(ih)*rtd,offa(ih)*rtd,
     |    dskofc(ih)*rtd,dskofa(ih)*rtd,theta0(ih)*rtd,rrad(ih)*rtd,
     |    phid(ih)*rtd/15.+12.,phin(ih)*rtd/15.+12.

      return
      end subroutine de2loc
!-----------------------------------------------------------------------
      subroutine setboundary_de2(angle,bt,tilt,swvel,swden,bndyfitr)
!
! Sets the coefficients that define the low-latitude boundary model,
!   given the IMF and solar wind values.
      use cons_module, only: dtr
!
      implicit none
!
! Args:
      real,intent(in) :: angle,bt,tilt,swvel,swden
      real,intent(out):: bndyfitr
!
! Local:
      integer :: i
      integer, parameter :: na=6
      real :: swp,xc,theta,ct,st,tilt2,cosa,btx,x(na),c(na) ! na = 6
      real :: tmat(3,3),ttmat(3,3)    
      real :: ex_bndy(2),bndya(6)
!
! Calculate the transformation matrix to the coordinate system
! of the offset pole.
!
      xc = 4.2
      theta = xc*(dtr)
      ct = cos(theta)
      st = sin(theta)
      
      ! the following is from wei05sc.nc files
      ex_bndy =(/0.1595479995, 0.2477539927/)
      bndya   =(/16.65162249,  -0.8299262114, 
     |           1.754533395,  -2.740607467, 
     |           0.02176603485, 0.2679466709/)
!
      tmat(1,:) = (/ ct, 0., st/) 
      tmat(2,:) = (/ 0., 1., 0./) 
      tmat(3,:) = (/-st, 0., ct/)
!
      ttmat(1,:) = (/ct, 0.,-st/)
      ttmat(2,:) = (/ 0.,1., 0./)
      ttmat(3,:) = (/st, 0., ct/)
!
      swp = swden*swvel**2*1.6726e-6 ! pressure
      tilt2 = tilt**2
      cosa = cos(angle*dtr)
      btx = 1.-exp(-bt*ex_bndy(1))
      if (bt > 1.) then
        btx = btx*bt**ex_bndy(2)
      else
        cosa = 1.+bt*(cosa-1.) ! remove angle dependency for IMF under 1 nT
      endif
      x = (/1., cosa, btx, btx*cosa, swvel, swp/)
      c = bndya
      bndyfitr = 0.
      do i=1,na
        bndyfitr = bndyfitr+x(i)*c(i)

!       write(6,"('setboundry: i=',i3,' bndyfitr=',e12.4)") i,bndyfitr

      enddo
      end subroutine setboundary_de2
!-----------------------------------------------------------------------
