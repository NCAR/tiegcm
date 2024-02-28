!
      module calculate_ecf_module
!
! Coded by Tong Dang, August, 2017
! This program is to calculate the rest euv flux during eclipse period
      use params_module,only: nlon,nlonp4,nlat
      use eclipse_module, only: eclip
!
      implicit none
!
      contains
!-----------------------------------------------------------------------
      subroutine calculate_ecf(eclipfactor,lon0,lon1,lev0,lev1,lat
     | ,newmag_factor)
!
! coded by Dang, 2017
      use params_module,only: glon,glat
      use hist_module, only: modeltime
      use chapman_module, only: chi
      use init_module, only: iyear, iday
! Args:
      integer,intent(in) :: lat,lon0,lon1,lev0,lev1
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: eclipfactor
     | ,newmag_factor
! Local:
      real :: longitude,latitude,height,ut
      real*8 :: newmag
      real :: t1_SLT,t4_SLT,tm_SLT,SLT,mag,EclipseFactor
      integer :: k,i

! Currently height is fixed at 300km for all altitudes
      height=300.
      ut=modeltime(2)+modeltime(3)/60.+modeltime(4)/3600.

!      if(modeltime(1) .gt. 233 .or. ut .gt. 21.2 .or. ut .lt. 15.6 .or.
!     . modeltime(1) .lt. 233) then
!         eclipfactor(lev0:lev1,lon0:lon1)=1.0
!         newmag_factor(lev0:lev1,lon0:lon1)=0.0
!      else

      do i=lon0,lon1

          if(i .ne. 1 .and. i .ne. 2 .and. i .ne. nlon+3 .and. 
     . i .ne. nlon+4) then
           longitude=glon(i-2)
           latitude=glat(lat)
           SLT=modeltime(2)+modeltime(3)/60.+modeltime(4)/3600.
     . + longitude/15.
           if(SLT .LT. 0) then
              SLT = SLT + 24.
           endif
           if(SLT .GT. 24) then
              SLT = SLT - 24.
           endif
         
          if(chi(i,lat)/3.14159*180. .le. (180.-asind(6371.0/(
     .6371.0+height)))) then 
           call eclip(ut,SLT,longitude,latitude,height,newmag, 
     .           t1_SLT,t4_SLT,tm_SLT,mag, iyear, iday)
           call Eclipse(newmag,EclipseFactor)
           eclipfactor(lev0:lev1,i)=EclipseFactor
           newmag_factor(lev0:lev1,i)=newmag
           else
             eclipfactor(lev0:lev1,i)=1.0
             newmag_factor(lev0:lev1,i)=0.0
           endif

          endif
       enddo
   
      if(lon0 .eq. 1) then 
        eclipfactor(lev0:lev1,1)=eclipfactor(lev0:lev1,3)
        eclipfactor(lev0:lev1,2)=eclipfactor(lev0:lev1,3)
      endif
      if(lon1 .eq. nlon+4)then
        eclipfactor(lev0:lev1,nlon+3)=eclipfactor(lev0:lev1,nlon+2)
        eclipfactor(lev0:lev1,nlon+4)=eclipfactor(lev0:lev1,nlon+2)
      endif

!      endif
      end subroutine calculate_ecf
!-----------------------------------------------------------------------
	Subroutine Eclipse(ecf,EclipseFactor)
CC    calculation of the ratio of φTi / φT0 = (12.3Spi+Sci)/(12.3sp0 + Sc0)
CC    calculation of the unmasked area for photosphere for P and C regions
CC    during the solar eclipse.
CC    P region is an inner central disk over the "photosphere";
CC    C region is an outer ring between the "coronal" and "phosphere".
CC    reference to J. J Curto et al.(2006), Modeling the geomagnetic effects 
CC    caused by the solar eclipse of 11 August 1999.

CC    initialize 
        implicit real*8 (A-Z)
        Rp=1.0; Rc=2.1; 
	  pi=3.1415926;
	  Sp0=pi*Rp*Rp;
	  Sc0=pi*(Rc*Rc-Rp*Rp);
CC      dl should be smaller than Rp+Rc; that is dl <= 3.1
C         v=3.1*2/Duration !! for duration=240 mins v= 3.1/120 mins =0.025833
C 	                   !! for duration=180 mins v= 3.1/90 mins =0.034444
C         dl=dabs(3.1-v*eclipseTime);  
        dl=2.0-ecf*2;

        if(dl<=Rc-Rp) then 
	      A1=0;
            E=2*Rp*Rp*acos(dl/2/Rp) - dl*sqrt(Rp*Rp - dl*dl*0.25);
	  end if
	  if(dl>Rc-Rp .AND. dl<2*Rp) then
	      E=2*Rp*Rp*acos(dl/2/Rp) - dl*sqrt(Rp*Rp - dl*dl*0.25);
	      Xd=(Rc*Rc - Rp*Rp + dl*dl)/2/dl;
            A2=Rc*Rc*acos(Xd/Rc) - Xd*sqrt(Rc*Rc -Xd*Xd);
	      A3=Rp*Rp*acos((dl-Xd)/Rp) - (dl-Xd)*
     & 	     sqrt(Rp*Rp-(dl-Xd)*(dl-Xd));
            A1=pi*Rp*Rp - A2 -A3;
	  endif
	  if(dl>=2*Rp .AND. dl<=3.1) then
	      E=0;
	      Xd=(Rc*Rc - Rp*Rp + dl*dl)/2/dl;
            A2=Rc*Rc*acos(Xd/Rc) - Xd*sqrt(Rc*Rc -Xd*Xd);
	      A3=Rp*Rp*acos((dl-Xd)/Rp) - (dl-Xd)*
     & 	     sqrt(Rp*Rp-(dl-Xd)*(dl-Xd));
            A1=pi*Rp*Rp - A2 -A3;
        endif

	  Sc=pi*(Rc*Rc - 2*Rp*Rp) + E + A1;
        Sp=pi*Rp*Rp - E;
	  EclipseFactor=(12.3*Sp+Sc)/(12.3*Sp0+Sc0);
	  if(ecf<0.0001) EclipseFactor=1.0
	End Subroutine Eclipse
!-----------------------------------------------------------------------
      end module calculate_ecf_module
