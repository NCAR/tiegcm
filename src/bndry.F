!
      module bndry_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Lower boundary conditions: diurnal, semi-diurnal, annual are
! lower boundary is only set in these subroutines if no input file
! from the GSWM is provided. Therefore a check is included for
! the tidal contribution (diurnal, semi-diurnal)
!
      use cons_module,only: t0,pi,atm_amu,gask,grav,freq_semidi,
     |  re,dlamda,tbound,tgrad,cs,cor,tn
      use params_module,only: nlat,nlonp4,nlevp1,dz,zibot
!
! GSWM boundary perturbation (see gswm.F):
      use gswm_module,only:  
     |  gswm_mi_di_z, gswm_mi_sdi_z, gswm_nm_di_z, gswm_nm_sdi_z,
     |  gswm_mi_di_t, gswm_mi_sdi_t, gswm_nm_di_t, gswm_nm_sdi_t,
     |  gswm_mi_di_u, gswm_mi_sdi_u, gswm_nm_di_u, gswm_nm_sdi_u,
     |  gswm_mi_di_v, gswm_mi_sdi_v, gswm_nm_di_v, gswm_nm_sdi_v
      implicit none
!
      complex,parameter :: ci=(0.,1.)
      complex ::
     |  zb(nlat),zb2(nlat),
     |  tb(nlat),tb2(nlat),
     |  ub(nlat),ub2(nlat),
     |  vb(nlat),vb2(nlat),
     |  bnd(nlonp4),bnd2(nlonp4)
!
      real :: b(nlonp4,2,2),fb(nlonp4,2)            ! for bndry_comp
!
! OX bottom boundary (set in sub comp_ox, and used in sub bndry_comp):
!     real :: xoxlb(nlonp4,nlat)
!
      contains
!-----------------------------------------------------------------------
      subroutine lowbound
!
      use init_module,only: igswm_mi_di,igswm_mi_sdi,	! GSWM input flags
     |   iday,iter		! for dday calculation
      use input_module,only: step
!     use chemrates_module,only: co2mix
!     use cco2gr_module,only: cco2gr,set_cco2_data
!
! Call the lower boundary routines. This is called once per model run
! from main tgcm.F.
!
! Local:
      integer :: k
      real :: dday
!
! tbound = value of T at lower boundary, and tgrad = gradient
!
      t0(1) = tbound
      t0(2) = tbound+dz*tgrad
      if(igswm_mi_sdi == 0) call bndry_semidiurnal  ! check for GSWM semidiurnal
      if(igswm_mi_di  == 0) call bndry_diurnal      ! check for GSWM diurnal
      dday = float(iday)+amod(float(iter)*float(step),86400.)/86400.
!
! t0==0:
      do k=1,nlevp1
        t0(K)=0.
      enddo
 !     
      end subroutine lowbound
!-----------------------------------------------------------------------
      subroutine bndry_semidiurnal
!
! Lower boundary conditions for semi-diurnal tide:
! 2/00: 1998 spherepack lib code (sphpac.f) replaces old lib 
!       alfpac.f for legendre polynomials and Hough functions.
! This routine calculates complex ZB,TB,UB,VB.
!
      use input_module,only: tide
!
! Local:
      integer,parameter :: nalf=19, malf=2
      real :: p(nlat,nalf,malf),hough(nlat,5,malf),cp(nalf/2+1)
      complex :: dzb(nlat)
      real :: B(5,19),RL(5),BHOUR(5),rlamda,xdot(19),ydot(19),
     |  ptscal,theta,ptjm(2*nlat+1)
      integer :: n,jm,l,lm1,m,mm1,j,ld,i,nm1
!
      COMPLEX ZEE(5),CL(5),EXPDLM
      DATA B/
     2 0.969152, 0.0     , 0.216046, 0.0     , 0.093838,
     3 0.0     , 0.909763, 0.0     , 0.342113, 0.0     ,
     4-0.245226, 0.0     , 0.798445, 0.0     , 0.421218,
     5 0.0     ,-0.408934, 0.0     , 0.645517, 0.0     ,
     6 0.024633, 0.0     ,-0.543993, 0.0     , 0.464159,
     7 0.0     , 0.071127, 0.0     ,-0.643189, 0.0     ,
     8-0.001292, 0.0     , 0.139613, 0.0     ,-0.699495,
     9 0.0     ,-0.006673, 0.0     , 0.225090, 0.0     ,
     1 0.000042, 0.0     ,-0.019654, 0.0     , 0.320141,
     1 0.0     , 0.000394, 0.0     ,-0.043345, 0.0     ,
     2-0.000001, 0.0     , 0.001772, 0.0     ,-0.079831,
     3 0.0     ,-0.000016, 0.0     , 0.005401, 0.0     ,
     4 0.0     , 0.0     ,-0.000112, 0.0     , 0.012932,
     5 0.0     , 0.0     , 0.0     ,-0.000476, 0.0     ,
     6 0.0     , 0.0     , 0.000005, 0.0     ,-0.001490,
     7 0.0     , 0.0     , 0.0     , 0.000031, 0.0     ,
     8 0.0     , 0.0     , 0.0     , 0.0     , 0.000129,
     9 0.0     , 0.0     , 0.0     ,-0.000002, 0.0     ,
     1 0.0     , 0.0     , 0.0     , 0.0     ,-0.000009/
      DATA RL/7.8519E5, 3.6665E5, 2.1098E5, 1.3671E5, 0.9565E5/
      real,external :: sddot ! in util.F
!
      bhour = tide(6:10)
      if (all(tide==0.)) goto 13
      DO N=1,5
        ZEE(N)=TIDE(N)*CEXP(CI*pi*BHOUR(N)/6.)
        CL(N)=CSQRT(CMPLX(gask/(atm_amu*grav*RL(N))*
     |    (T0(1)*2./7.+(T0(2)-T0(1))/dz)-.25))-.5*CI
      enddo
      JM=2*nlat+1
!
C     ****     SET UP HOUGH FUNCTIONS
! 
! Using new (1998) spherepack (sphpac.f):
      do n=2,nalf+1
        nm1 = n-1
        do m=2,malf+1
          mm1=m-1
          call alfk(n,m,cp)
          do j=1,jm
            theta = float(j-1)*pi/float(jm-1)
            call lfpt(n,m,theta,cp,ptscal)
            ptjm(j) = ptscal
          enddo
          do j=1,nlat
            p(j,nm1,mm1) = ptjm(2*(nlat+1-j))
          enddo
        enddo
        DO J=1,nlat
          P(J,nm1,2)=SQRT(FLOAT(n*(n+1)-6))*P(J,nm1,2)-2.*TN(J)*
     |               P(J,nm1,1)
        enddo
      enddo
!
! util.F: real function sddot(n,x,y)
      DO L=1,5
        DO LD=1,2
          DO J=1,nlat
            xdot(:) = p(j,:,ld)
            ydot(:) = b(l,:)
            HOUGH(J,L,LD)=sddot(19,xdot,ydot)
          enddo
        enddo
      enddo
C     ****      GENERATE ZB, UB, VB AND TB
      DO 5 J=1,nlat
        TB(J)=0.
        ZB(J)=0.
        DZB(J)=0.
    5 CONTINUE
      DO L=1,5
        DO J=1,nlat
          ZB(J)=ZB(J)+ZEE(L)*HOUGH(J,L,1)
          DZB(J)=DZB(J)+ZEE(L)*HOUGH(J,L,2)
          TB(J)=TB(J)+CI*atm_amu*grav/gask*ZEE(L)*CL(L)*HOUGH(J,L,1)
        enddo
      enddo
      DO 7 J=1,nlat
        UB(J)=freq_semidi*re*(1.-(COR(J)/freq_semidi)**2)
        VB(J)=CI*grav*(DZB(J)-2.*COR(J)/(freq_semidi*CS(J))*ZB(J))/UB(J)
        UB(J)=grav*(COR(J)/freq_semidi*DZB(J)-2./CS(J)*ZB(J))/UB(J)
    7 CONTINUE
      GO TO 11
   13 CONTINUE
C     ****      ZERO BOUNDARY CONDITION
      DO 12 J=1,nlat
        ZB(J)=0.
        TB(J)=0.
        UB(J)=0.
        VB(J)=0.
   12 CONTINUE
   11 CONTINUE
C     ****      CALCULATE LONGITUDINAL STRUCTURE
      RLAMDA = -2.*dlamda
      BND(1)=CEXP(CI*2.*RLAMDA)
      EXPDLM=CEXP(CI*2.*dlamda)
      DO 9 I=2,nlonp4
        BND(I)=BND(I-1)*EXPDLM
    9 CONTINUE
      end subroutine bndry_semidiurnal
!-----------------------------------------------------------------------
      subroutine bndry_diurnal
!
C     ****     TIDAL BOUNDARY CONDITION FOR DIURNAL GRAVITATIONAL MODE
C     ****     (1,1)
! 2/00: 1998 spherepack lib code (sphpac.f) replaces old lib 
!       alfpac.f for legendre polynomials and Hough functions.
! This routine calculates complex ZB2,TB2,UB2,VB2.
!
      use input_module,only: tide2
!
! Local:
      integer,parameter :: nalf=19, malf=2
      real :: p(nlat,nalf,malf),hough(nlat,5,malf),cp(nalf/2+1)
      complex :: dzb(nlat)
      real :: B(1,19),RL (1),BHOUR2(1),rlamda,xdot(19),ydot(19),
     |  ptscal,theta,ptjm(2*nlat+1),pik
      integer :: l,m,j,n,jm,ld,i
      COMPLEX ZEE(1),CL(1),EXPDLM
!
      DATA B/
     2 0.282710,
     3 0.0     ,
     4-0.638229,
     5 0.0     ,
     6 0.620521,
     7 0.0     ,
     8-0.336408,
     9 0.0     ,
     1 0.117021,
     1 0.0     ,
     2-0.028332,
     3 0.0     ,
     4 0.005042,
     5 0.0     ,
     6-0.000686,
     7 0.0     ,
     8 0.000074,
     9 0.0     ,
     1-0.000006/
      DATA RL/0.6909E5/
      real,external :: sddot ! in util.F
!
      bhour2(1) = tide2(2)
!
! pik is as per tgcm15 (4*atan(1)), without -lmass (pi with and without
! -lmass differs in 16-17 decimal places, apparently causing diffs in zee, 
! and subsequently zb2). This does not seem to be a problem with semidiurnal.
!
      pik = 3.14159265358979312 
      if (all(tide2==0.)) goto 13
      DO 1 N=1,1
        ZEE(N)=TIDE2(N)*CEXP(CI*pik*BHOUR2(N)/12.)
        CL(N)=CSQRT(CMPLX(gask/(atm_amu*grav*RL(N))*
     |  (T0(1)*2./7.+(T0(2)-T0(1))/dz)-.25))-.5*CI
    1 CONTINUE
      JM=2*nlat+1
!
C     ****     SET UP HOUGH FUNCTIONS
!
! Using new (1998) spherepack (sphpac.f):
      do n=1,19
        do m=1,2
          call alfk(n,m,cp)
          do j=1,jm
            theta = float(j-1)*pi/float(jm-1)
            call lfpt(n,m,theta,cp,ptscal)
            ptjm(j)=ptscal
          enddo
          do j=1,nlat
            p(j,n,m) = ptjm(2*(nlat+1-j))
          enddo
        enddo
        DO J=1,nlat
          P(J,n,2)=SQRT(FLOAT(n*(n+1)-2))*P(J,n,2)-TN(J)*P(J,n,1)
        enddo
      enddo
!
! util.F: real function sddot(n,x,y)
      DO L=1,1
        DO LD=1,2
          DO J=1,nlat
            xdot(:) = p(j,:,ld)
            ydot(:) = b(l,:)
            HOUGH(J,L,LD)=sddot(19,xdot,ydot)
          enddo
        enddo
      enddo
C     ****      GENERATE ZB2, UB2, VB2 AND TB2
      DO 5 J=1,nlat
      TB2(J)=0.
      ZB2(J)=0.
      DZB(J)=0.
    5 CONTINUE
      DO 6 L=1,1
      DO 6 J=1,nlat
      ZB2(J)=ZB2(J)+ZEE(L)*HOUGH(J,L,1)
      DZB(J)=DZB(J)+ZEE(L)*HOUGH(J,L,2)
      TB2(J)=TB2(J)+CI*atm_amu*grav/gask*ZEE(L)*CL(L)*HOUGH(J,L,1)
    6 CONTINUE
      DO 7 J=1,nlat
      UB2(J)=.5*freq_semidi*re*(1.-(COR(J)/(.5*freq_semidi))**2)
      VB2(J)=CI*grav*(DZB(J)-COR(J)/(.5*freq_semidi*CS(J))*ZB2(J))/
     |  UB2(J)
      UB2(J)=grav*(COR(J)/(.5*freq_semidi)*DZB(J)-1./CS(J)*ZB2(J))/
     |  UB2(J)
    7 CONTINUE
      GO TO 11
   13 CONTINUE
C     ****      ZERO BOUNDARY CONDITION
      DO 12 J=1,nlat
      ZB2(J)=0.
      TB2(J)=0.
      UB2(J)=0.
      VB2(J)=0.
   12 CONTINUE
   11 CONTINUE
C     ****      CALCULATE LONGITUDINAL STRUCTURE
      RLAMDA = -2.*dlamda
      BND2(1)=CEXP(CI*RLAMDA)
      EXPDLM=CEXP(CI*dlamda)
      DO 9 I=2,nlonp4
      BND2(I)=BND2(I-1)*EXPDLM
    9 CONTINUE
      end subroutine bndry_diurnal
!-----------------------------------------------------------------------
      subroutine bndcmp
C     ****
C     ****     CALCULATE MATRICES B(nlonp4,2,2) AND VECTORS FB(nlonp4,2)
C     ****       REPRESENTING THE LOWER BOUBNDARY CONDITION IN COMP,
C     ****       WHERE PSI1 AND PSI2 ARE CALCULATED:
C     ****
C     ****         PSI(K=-1/2) = B * PSI(K=1/2) + FB
C     ****
C     ****     BNDCMP CALLS THE SUBROUTINE BNDEF TO DEFINE THE 2 X 2
C     ****       MATRICES E, F AND THE 2 VECTOR G IN THE GENERAL
C     ****       LOWER BOUNDARY CONDITION:
C     ****
C     ****         E * D(PSI)/DS + F * PSI + G = 0.
C     ****
C     ****         WHERE:
C     ****           PSI = |PSI1| AND THE BOUNDARY CONDITION IS APPLIED
C     ****                 |    |
C     ****                 |PSI2|
C     ****
C     ****           AT LEVEL ZERO
C     ****
C     ****     THIS SUBROUTINE THEN EVALUATES B AND FB FROM:
C     ****
C     ****       B = (E/DS - F/2.)**(-1) * (E/DS + F/2.)
C     ****
C     ****       FB = (E/DS - F/2.)**(-1) * G
C     ****
!
! Local:
      real :: ee(nlonp4,2,2),ff(nlonp4,2,2),gg(nlonp4,2),
     |  wm1(nlonp4,2,2),wm2(nlonp4,2,2),wm3(nlonp4,2,2),
     |  ws1(nlonp4)
      integer :: l,m,i
C     ****
C     ****     CALL BNDEF TO DEFINE E, F AND G IN S1, S2 AND S3
C     ****
      call bndef(ee,ff,gg)
C     ****
C     ****     WM1 = (E/DS - F/2.)
C     ****
C     ****     WM1 = (E/DS + F/2.)
C     ****
      do l = 1,2
        do m = 1,2
          do i = 1,nlonp4
            wm1(i,l,m) = ee(i,l,m)/dz-ff(i,l,m)/2.
            wm2(i,l,m) = ee(i,l,m)/dz+ff(i,l,m)/2.
          enddo
        enddo
      enddo
C     ****
C     ****     WM3 = WM1**(-1)
C     ****
C     ****       WS1 = DET(WM1)
C     ****
      do i = 1,nlonp4
        ws1(i) = wm1(i,1,1)*wm1(i,2,2)-wm1(i,1,2)*wm1(i,2,1)
      enddo
C     ****
C     ****     NOW INVERSE OF WM1 IN WM3
C     ****
      do i = 1,nlonp4
        wm3(i,1,1) =  wm1(i,2,2)/ws1(i)
        wm3(i,1,2) = -wm1(i,1,2)/ws1(i)
        wm3(i,2,1) = -wm1(i,2,1)/ws1(i)
        wm3(i,2,2) =  wm1(i,1,1)/ws1(i)
      enddo
C     ****
C     ****     B = WM3 * WM2
C     ****
! b and fb are bndry_module module data.
      do l = 1,2
        do m = 1,2
          do i = 1,nlonp4
            b(i,l,m) = wm3(i,l,1)*wm2(i,1,m)+wm3(i,l,2)*wm2(i,2,m)
          enddo
        enddo
      enddo
C     ****
C     ****     FB = WM3 * G
C     ****
      do l = 1,2
        do i = 1,nlonp4
          fb(i,l) = wm3(i,l,1)*gg(i,1)+wm3(i,l,2)*gg(i,2)
        enddo
      enddo
      end subroutine bndcmp
!-----------------------------------------------------------------------
      SUBROUTINE BNDEF(ee,ff,gg)
C     ****
C     ****     BNDEF DEFINES THE LOWER BOUNDARY CONDITION FOR THIS
C     ****       VERSION OF THE MODEL
C     ****
C     ****     THE LOWER BOUNDARY CONDITION FOR COMP IS:
C     ****
C     ****       E * D(PSI)/DS + F * PSI +G = 0.
C     ****
C     ****         WHERE:
C     ****           PSI = VECTOR(PSI1,PSI2)
C     ****           E AND F ARE 2 X 2 MATRICES
C     ****           G = VECTOR(G1,G2)
C     ****
C     ****           E, F AND G MAY BE FUNCTIONS OF LATITUDE & LONGITUDE
C     ****
C     ****       THIS SUBROUTINE DEFINES E, F AND G for BNDCMP
C     ****
!
! Args:
      real,intent(out) :: EE(nlonp4,2,2),FF(nlonp4,2,2),GG(nlonp4,2)
!
! Local:
      real :: alfa
      integer :: i
C     ****
C     ****     IN TIGCM AND TIEGCM:
C     ****
C     ****       E = |0.  0.|
C     ****           |      |
C     ****           |0.  1.|
C     ****
C     ****       F = |1.  1.|
C     ****           |      |
C     ****           |0. -1.|
C     ****
C     ****       G = |-ALFA|
C     ****           |     |
C     ****           |   0.|
C     ****
C     ****       WHERE:
C     ****         ALFA = 0.22 + 0.014 = 0.234
C     ****
      DATA ALFA/0.234/
      DO 1 I = 1,nlonp4
        EE(I,1,1) = 0.
        EE(I,1,2) = 0.
        EE(I,2,1) = 0.
        EE(I,2,2) = 1.
        FF(I,1,1) = 1.
        FF(I,1,2) = 1.
        FF(I,2,1) = 0.
        FF(I,2,2) = -1.
        GG(I,1) = -ALFA
        GG(I,2) = 0.
    1 CONTINUE
      RETURN
      end subroutine bndef
!
!-----------------------------------------------------------------------
! Begin GSWM subroutines (formerly in bndry_gswm.F):
!
! Get gswm tn lbc perturbations for sub dt.
!-----------------------------------------------------------------------
!
      subroutine lbc_gswm_dt(tnlbc,lon0,lon1,lat0,lat1)
      
      use input_module,only: step
      use init_module,only:  iter,igswm_mi_di,igswm_mi_sdi,
     |   igswm_nm_di,igswm_nm_sdi
      use cons_module,only:  freq_semidi,tbound
      
      implicit none
      
      
 ! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1      ! dimension for tnlbc    
      real,intent(out) :: tnlbc(lon0:lon1,lat0:lat1) ! lower boundary condition

! Local:
      integer :: i,lat
      real :: rstep
      complex :: expt,expt2,expta
!     
! Calculate exponentials
!      
      rstep = float(step)
      expt  = cexp(ci*freq_semidi*rstep*iter)
      expt2 = cexp(ci*.5*freq_semidi*rstep*iter)
      expta = 1.
!
! GSWM migrating diurnal and semi-diurnal:
      if(igswm_mi_di == 1.and.igswm_mi_sdi == 1) then
	do lat=lat0,lat1
	  do i=lon0,lon1
 	    tnlbc(i,lat) = gswm_mi_sdi_t(i,lat)+tbound	     ! semidiurnal tide
	    tnlbc(i,lat) = tnlbc(i,lat)+ gswm_mi_di_t(i,lat) ! diurnal tide
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
!      
! GSWM migrating semi-diurnal:
      elseif(igswm_mi_di == 0.and.igswm_mi_sdi == 1) then
	do lat=lat0,lat1
	  do i=lon0,lon1
 	    tnlbc(i,lat) = gswm_mi_sdi_t(i,lat)+tbound    ! semidiurnal tide
	    tnlbc(i,lat) = tnlbc(i,lat)+
     |        real(tb2(lat)*bnd2(i)*expt2)	       ! diurnal tide
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
!      
! GSWM migrating diurnal:
      elseif(igswm_mi_di == 1.and.igswm_mi_sdi == 0) then
	do lat=lat0,lat1
	  do i=lon0,lon1
 	    tnlbc(i,lat) = real(tb(lat)*bnd(i)*expt)+tbound  ! semidiurnal tide
	    tnlbc(i,lat) = tnlbc(i,lat)+ gswm_mi_di_t(i,lat)    ! diurnal tide
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
      else
! 
! No gswm:
        do lat=lat0,lat1
          do i=lon0,lon1
             tnlbc(i,lat) = real(tb(lat)*bnd(i)*expt)+tbound ! semidiurnal tide
            tnlbc(i,lat) = tnlbc(i,lat)+
     |        real(tb2(lat)*bnd2(i)*expt2) 		    ! diurnal tide
          enddo   ! i=lon0,lon1
        enddo	  ! j=lat0,lat1
      endif
!
! GSWM non-migrating diurnal:
      if(igswm_nm_di == 1) then  ! nonmigrating diurnal tide
	do lat=lat0,lat1
	  do i=lon0,lon1
	    tnlbc(i,lat) = tnlbc(i,lat)+ gswm_nm_di_t(i,lat)
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
      endif
!
! GSWM non-migrating semi-diurnal:
      if(igswm_nm_sdi == 1) then  ! nonmigrating semidiurnal tide
	do lat=lat0,lat1
	  do i=lon0,lon1
	    tnlbc(i,lat) = tnlbc(i,lat)+ gswm_nm_sdi_t(i,lat)
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
      endif
!      
      end subroutine lbc_gswm_dt
!
!-----------------------------------------------------------------------
! calculate geopotential height
! z(1) = zb
!
      subroutine lbc_gswm_addiag(z,lon0,lon1,lat0,lat1)
!      
      use input_module,only: step
      use init_module,only:  iter,igswm_mi_di,igswm_mi_sdi,
     |   igswm_nm_di,igswm_nm_sdi
      use cons_module,only:  freq_semidi,dt
!      
      implicit none
!     
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1     ! dimension for tnlbc    
      real,intent(out) :: z(lon0:lon1,lat0:lat1) ! lower boundary condition

! Local:
      integer :: i,j
      real :: rstep
      complex :: expt,expt2,expta
!
! Calculate exponentials
!      
      expt  = cexp(ci*freq_semidi*dt*iter)
      expt2 = cexp(ci*.5*freq_semidi*dt*iter)
      expta = 1.
!
! GSWM migrating diurnal and semi-diurnal:
      if(igswm_mi_di == 1.and.igswm_mi_sdi == 1) then
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = gswm_mi_sdi_z(i,j)		  ! semidiurnal tide
	    z(i,j) = z(i,j)+ gswm_mi_di_z(i,j)	    ! diurnal tide
	  enddo
	enddo
!
! GSWM migrating semi-diurnal:
      elseif(igswm_mi_di == 0.and.igswm_mi_sdi == 1) then
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = gswm_mi_sdi_z(i,j)		  ! semidiurnal tide
	    z(i,j) = z(i,j)+real(zb2(j)*bnd2(i)*expt2) ! diurnal tide
	  enddo
	enddo
!
! GSWM migrating diurnal:
      elseif(igswm_mi_di == 1.and.igswm_mi_sdi == 0) then
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = real(zb(j)*bnd(i)*expt)	       ! semidiurnal tide
	    z(i,j) = z(i,j)+ gswm_mi_di_z(i,j)	       ! diurnal tide
	  enddo
	enddo
!
! No gswm:
      else
        do j = lat0,lat1
          do i=lon0,lon1
            z(i,j) = real(zb(j)*bnd(i)*expt)	       ! semidiurnal tide
            z(i,j) = z(i,j)+real(zb2(j)*bnd2(i)*expt2) ! diurnal tide
          enddo
        enddo
      endif
!     
! GSWM non-migrating diurnal:
      if(igswm_nm_di == 1) then	! nonmigrating diurnal tide
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = z(i,j)+ gswm_nm_di_z(i,j)
	  enddo
	enddo
      endif
!     
! GSWM non-migrating semi-diurnal:
      if(igswm_nm_sdi == 1) then	! nonmigrating semidiurnal tide
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = z(i,j)+ gswm_nm_sdi_z(i,j)
	  enddo
	enddo
      endif
      end subroutine lbc_gswm_addiag
!-----------------------------------------------------------------------
      subroutine lbc_gswm_duv(unlbc,vnlbc,unlbc_diag,vnlbc_diag,
     |    lon0,lon1,lat,lev0,lev1,expt,expt2,expta)
      use init_module,only:  igswm_mi_di,igswm_mi_sdi,
     |   igswm_nm_di,igswm_nm_sdi
!     
      implicit none
! 
! Args:
      integer,intent(in) :: lon0,lon1,lat,lev0,lev1  ! dimension for lbc    
      complex,intent(in) :: expt,expt2,expta 
      real,intent(out) :: unlbc(lon0:lon1),	! lower boundary condition
     |                    vnlbc(lon0:lon1), 
     |                    unlbc_diag(lev0:lev1,lon0:lon1), ! for output
     |                    vnlbc_diag(lev0:lev1,lon0:lon1)
! Local:
      integer :: i
!     
! GSWM migrating diurnal and semi-diurnal:
      if(igswm_mi_di == 1.and.igswm_mi_sdi == 1) then
        do i=lon0,lon1
          unlbc(i) = gswm_mi_sdi_u(i,lat)  	       ! semidiurnal tide
          vnlbc(i) = gswm_mi_sdi_v(i,lat)
          unlbc(i) = unlbc(i) + gswm_mi_di_u(i,lat)    ! diurnal tide
          vnlbc(i) = vnlbc(i) + gswm_mi_di_v(i,lat)
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
!
! GSWM migrating semi-diurnal:
      elseif(igswm_mi_di == 0.and.igswm_mi_sdi == 1) then
        do i=lon0,lon1
          unlbc(i) = gswm_mi_sdi_u(i,lat)  	         ! semidiurnal tide
          vnlbc(i) = gswm_mi_sdi_v(i,lat)
          unlbc(i) = unlbc(i) + real(ub2(lat)*bnd2(i)*expt2) ! diurnal tide
          vnlbc(i) = vnlbc(i) + real(vb2(lat)*bnd2(i)*expt2)
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
!
! GSWM migrating diurnal:
      elseif(igswm_mi_di == 1.and.igswm_mi_sdi == 0) then
        do i=lon0,lon1
          unlbc(i) =  real(ub(lat)*bnd(i)*expt)         ! semidiurnal tide
          vnlbc(i) =  real(vb(lat)*bnd(i)*expt)
          unlbc(i) = unlbc(i) + gswm_mi_di_u(i,lat)            ! diurnal tide
          vnlbc(i) = vnlbc(i) + gswm_mi_di_v(i,lat)
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
!
! No gswm:
      else
        do i=lon0,lon1
          unlbc(i) = real(ub(lat)*bnd(i)*expt)          ! Semidiurnal tide
          vnlbc(i) = real(vb(lat)*bnd(i)*expt)
          unlbc(i) = unlbc(i) + real(ub2(lat)*bnd2(i)*expt2)! diurnal tide
          vnlbc(i) = vnlbc(i) + real(vb2(lat)*bnd2(i)*expt2)
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
      endif
!    
! GSWM non-migrating diurnal:
      if(igswm_nm_di == 1) then      ! nonmigrating diurnal tide
        do i=lon0,lon1
          unlbc(i) = unlbc(i) + gswm_nm_di_u(i,lat)
          vnlbc(i) = vnlbc(i) + gswm_nm_di_v(i,lat)
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
      endif
!    
! GSWM non-migrating semi-diurnal:
      if(igswm_nm_sdi == 1) then      ! nonmigrating semidiurnal tide
        do i=lon0,lon1
          unlbc(i) = unlbc(i) + gswm_nm_sdi_u(i,lat)
          vnlbc(i) = vnlbc(i) + gswm_nm_sdi_v(i,lat)
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
      endif
      end subroutine lbc_gswm_duv
!-----------------------------------------------------------------------
      end module bndry_module
