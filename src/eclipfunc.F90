!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        module eclipse_module
        implicit none
        real*8 :: obsvconst(0:5)
        real*8 :: elements(0:27)
        real*8 :: c1(0:39),c2(0:39),c3(0:39),c4(0:39),mid(0:39)
        real*8 :: duration, V,P, azi, alt, t,cover
        integer :: year,month,day,typepe
        logical :: doFill = .true.

        real*8,dimension(:,:),allocatable :: eclipse_list

  contains
  !!***********************************************************
  subroutine get_ymd(iyear, idoy, month, day)

!  USE init_module, only: iyear, iday
  integer, intent(in):: iyear, idoy
  integer, intent(out):: month, day
  integer, dimension(12):: month_leapyear=[31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  integer, dimension(12):: month_notleap=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  integer :: i, doy, year

  doy = idoy
  year = iyear
  if (MOD(year, 4)==0 .and. (MOD(year, 100) /= 0 .or. MOD(year, 400)==0)) then
    do i = 1, 12
       if (doy>month_leapyear(i)) then
          doy = doy-month_leapyear(i)
          cycle
       end if
       if (doy <= month_leapyear(i)) then
          month = i
          day = doy
          exit
       end if
    end do
  else
     do i = 1, 12
        if (doy > month_notleap(i)) then
           doy = doy-month_notleap(i)
           cycle
        end if

        if (doy <= month_notleap(i)) then
           month = i
           day = doy
           exit
        end if
     end do
  end if

  end subroutine get_ymd

  subroutine eclip(UT,SLT,longitude,latitude,height,newmag, &
                  & t1_SLT,t4_SLT,tm_SLT,mag, iyear, idoy)
   real, parameter  :: PI =3.14151926
   real*8 :: tmp,mag,ratio,vv(1:5),pp(1:5),azimuth(1:5),altitude(1:5),timee(1:5)
   integer :: i,j, eclitype, iyear, idoy  !! 0 - none, 1 - partial, 2 - annular, 3 - total
   real :: UT,sUT,longitude,latitude,height,t1,t4,tm
   real*8 :: newmag
   real :: t1_SLT,t4_SLT,tm_SLT,SLT,x;
! Dang, 2017
   real :: R,M,D,S1 
    x=longitude

!! Observer constants -
!! (0) North Latitude (degree)
!! (1) East Longitude (degree)
!! (2) Altitude (metres)
!! (3) West time zone (hours)
!! (4) rho dsin O'
!! (5) rho dcos O'
!! (6) index into the elements array for the eclipse in question
!!
    ! Dang, 2017

    call get_ymd(iyear, idoy, month, day)
    year = iyear
!    year = 2021
!    month= 6
!    day = 10
!    print *, 'label_cxt', year, month, day ! for debug
!   year=2021
!   month=6
!   day=10
!! Note that correcting for refraction will involve creating a "virtual" altitude
!! for each contact, and hence a different value of rho and O' for each contact!
   obsvconst(0)= latitude   !! degree
   obsvconst(1)= - longitude;   !! degree
   obsvconst(2)=height*1000.0; !! meter
   obsvconst(3)=0.0;  
!!! Eclipse Elements------------------------------
!
! First line -
!  (0) Julian date of maximum eclipse
!  (1) t0 - the TDT hour at which t=0
!  (2) tmin - the lowest allowed value of t
!  (3) tmax - the highest allowed value of t
!  (4) dUTC - the difference between the civilian "GMT" timescale and TDT
!  (5) dT - the difference between UT and TDT
! Second line -
!  (6) X0, X1, X2, X3 - X elements
! Third line -
! (10) Y0, Y1, Y2, Y3 - Y elements
! Fourth line -
! (14) D0, D1, D2 - D elements
! Fifth line -
! (17) M0, M1, M2 - mu elements
! Sixth line -
! (20) L10, L11, L12 - L1 elements
! Seventh line -
! (23) L20, L21, L22 - L2 elements
! Eigth line -
! (26) tan f1
! (27) tan f2------------------------------------------------------------
! Dang, 2017
      call findEclipse(year,month,day)
      tmp=atan(0.99664719*tan(obsvconst(0)*Pi/180.0));
      obsvconst(4)=0.99664719*dsin(tmp)+(obsvconst(2)/6378140.0)*dsin(obsvconst(0)*Pi/180.0);
      obsvconst(5)=dcos(tmp)+(obsvconst(2)/6378140.0*dcos(obsvconst(0)*Pi/180.0));
!! initial the values 
      timee=0; altitude=0; azimuth=0; pp=0; vv=0; eclitype=0; mag=0; ratio=0; duration=0; cover=0;
	  t1=0; t4=0; tm=0
      call getall()
!!==============================================================================================

!! Eclipse circumstances
!!  (0) Event type (C1=-2, C2=-1, Mid=0, C3=1, C4=2)
!!  (1) t
!! -- time-only dependent circumstances (and their per-hour derivatives) follow --
!!  (2) x  (3) y  (4)  d (5) dsin d  (6) dcos d  (7) mu  (8) l1  (9) l2  (10) dx
!! (11) dy (12) dd  (13) dmu  (14) dl1  (15) dl2
!! -- time and location dependent circumstances follow --
!! (16) h  (17) dsin h  (18) dcos h  (19) xi  (20) eta  (21) zeta  (22) dxi  (23) deta  (24) u
!! (25) v  (26) a  (27) b  (28) l1' (29) l2' (30) n^2
!! -- observational circumstances follow --
!! (31) p  (32) alt  (33) q  (34) v  (35) azi
!! (36) m (mid eclipse only) or limb correction applied (where available!)
!! (37) magnitude (mid eclipse only)
!! (38) moon/sun (mid eclipse only)
!! (39) calculated local event type for a transparent earth (0 = none, 1 = partial, 2 = annular, 3 = total)
	  if (mid(39)==0) then
	     newmag=0.0  !! there is no eclipse
		 goto 30
	  endif
!! displaymid--------------------------------------------------  
      call gettime(mid);  timee(5)=t; tm=t;
      call getalt(mid);   altitude(5)=alt;
      call getazi(mid);   azimuth(5)=azi;
      call getp(mid);     pp(5)=p;
      call getv(mid);     vv(5)=v;
!! Display the information about 1st and 4th contact--------------
      call gettime(c1); timee(1)=t; t1=t;
      call getalt(c1);  altitude(1)=alt;
      call getazi(c1);  azimuth(1)=azi;
      call getp(c1);    pp(1)=p;
      call getv(c1);    vv(1)=v;
      call gettime(c4); timee(4)=t; t4=t;
      call getalt(c4);  altitude(4)=alt;
      call getazi(c4);  azimuth(4)=azi;
      call getp(c4);    pp(4)=p;
      call getv(c4);    vv(4)=v;
      mag=floor(100000.0*mid(37)+0.5)/100000.0;
      ratio= floor(100000.0*mid(38)+0.5)/100000.0;
!	  t1_SLT=mod(t1+longitude/15.0+24.0,24.0);
!	  t4_SLT=mod(t4+longitude/15.0+24.0,24.0);
!	  tm_SLT=mod(tm+longitude/15.0+24.0,24.0);
          t1_SLT=mod(t1+x/15.0+24.0,24.0);
          t4_SLT=mod(t4+x/15.0+24.0,24.0);
          tm_SLT=mod(tm+x/15.0+24.0,24.0);

!	  if ((SLT<t1_SLT .OR. SLT>t4_SLT) .AND. (SLT>4 .AND. SLT<23) ) then
          if ((t1_SLT<t4_SLT .AND. (SLT<t1_SLT .OR. SLT>t4_SLT)) .OR. (t1_SLT>=t4_SLT .AND. (SLT<t1_SLT .AND. SLT>t4_SLT))) then
!              if(latitude .gt. 85 .and. SLT .gt. 22)
!              print *,'SLT:',SLT,t1_SLT,t4_SLT
!              endif
	      newmag=0;
		  goto 30
	  endif

!! linear interpolation     
!	  if (SLT<tm_SLT) then
!		  newmag=mag*(SLT-t1_SLT)/(tm_SLT-t1_SLT);
!	  else
!		  newmag=mag*(SLT-t4_SLT)/(tm_SLT-t4_SLT);
!	  endif

! Calculate the actual ratio of the obscured diameter of the solar photosphere   
! Dang, 2017

      R=0.5
      M=mag   
      D=sqrt(4*M*R-M*M)
      if (UT<=timee(5)) then 
        d=(timee(5)-UT)/(timee(5)-timee(1))*D
      else
        d=(timee(5)-UT)/(timee(5)-timee(4))*D
      endif
      S1=sqrt(d*d+(2*R-M)*(2*R-M))
      newmag=2*R-S1
30    return
  End subroutine eclip

!!-----------------------------Calculate a limb correction------------------------------------
   subroutine limbcorrectionc2(p, q)

      integer ::  n
	  real*8 :: i,C2limb2003May(0:317),p,q
	  real*8, parameter  :: PI =3.14151926

! C2 limb corrections for the 2003 May 31 annular eclipse in seconds--------------
!
! The first 3 elements of the array tell us that the remaining data starts at contact angle
! 197.37 degrees in 0.4 degree increments for 315 data points, which means that the last
! element is for angle 322.97 degrees
!
! These limb corrections were calculated by Fred Espenak, NASA/GSFC---------------
     C2limb2003May = [197.37, 0.4, 315.0, &
 &  9.85, 9.35, 8.85, 8.33, 7.80, 7.26, 6.71, 6.14, 5.56, 5.22, 5.06, &
 &  4.89, 4.70, 4.51, 4.31, 4.10, 3.88, 3.82, 3.91, 4.00, 4.07, 4.14, &
 &  4.20, 4.26, 4.31, 4.35, 4.38, 4.41, 4.43, 4.44, 4.44, 4.44, 4.43, &
 &  4.42, 4.39, 4.36, 4.33, 4.28, 4.23, 4.17, 4.11, 4.04, 3.96, 3.87, &
 &  3.84, 3.85, 3.86, 3.86, 3.86, 3.85, 3.83, 3.81, 3.78, 3.74, 3.70, &
 &  3.65, 3.59, 3.53, 3.47, 3.39, 3.31, 3.22, 3.13, 3.03, 2.93, 2.81, &
 &  2.70, 2.57, 2.44, 2.31, 2.16, 2.01, 1.86, 1.70, 1.59, 1.64, 1.69, &
 &  1.72, 1.76, 1.78, 1.81, 1.87, 1.92, 1.97, 2.01, 2.05, 2.08, 2.11, &
 &  2.13, 2.15, 2.16, 2.23, 2.37, 2.51, 2.64, 2.77, 2.89, 3.01, 3.12, &
 &  3.23, 3.33, 3.43, 3.52, 3.61, 3.69, 3.77, 3.84, 3.91, 3.97, 4.03, &
 &  4.08, 4.13, 4.17, 4.21, 4.24, 4.27, 4.29, 4.31, 4.32, 4.33, 4.33, &
 &  4.33, 4.32, 4.31, 4.29, 4.27, 4.24, 4.20, 4.17, 4.12, 4.08, 4.02, &
 &  3.96, 3.90, 3.83, 3.76, 3.68, 3.60, 3.51, 3.42, 3.32, 3.22, 3.11, &
 &  3.00, 2.88, 2.76, 2.63, 2.53, 2.50, 2.46, 2.41, 2.36, 2.31, 2.25, &
 &  2.19, 2.12, 2.04, 1.96, 1.88, 1.79, 1.70, 1.60, 1.49, 1.39, 1.27, &
 &  1.15, 1.03, 0.90, 0.77, 0.63, 0.49, 0.34, 0.19, 0.17, 0.21, 0.24, &
 &  0.26, 0.29, 0.30, 0.31, 0.32, 0.32, 0.31, 0.30, 0.29, 0.26, 0.24, &
 &  0.21, 0.17, 0.13, 0.09, 0.05, 0.18, 0.31, 0.43, 0.55, 0.66, 0.77, &
 &  0.87, 0.97, 1.05, 1.14, 1.22, 1.29, 1.36, 1.42, 1.48, 1.53, 1.57, &
 &  1.61, 1.65, 1.68, 1.70, 1.72, 1.73, 1.73, 1.73, 1.73, 1.72, 1.70, &
 &  1.69, 1.75, 1.81, 1.87, 1.92, 1.96, 2.00, 2.03, 2.06, 2.08, 2.10, &
 &  2.11, 2.12, 2.12, 2.12, 2.13, 2.16, 2.18, 2.19, 2.20, 2.21, 2.25, &
 &  2.35, 2.44, 2.52, 2.59, 2.66, 2.71, 2.77, 2.81, 2.85, 2.89, 2.91, &
 &  2.93, 2.94, 2.95, 2.94, 3.07, 3.22, 3.37, 3.51, 3.64, 3.77, 3.88, &
 &  3.99, 4.09, 4.18, 4.27, 4.34, 4.41, 4.47, 4.52, 4.57, 4.60, 4.63, &
 &  4.65, 4.66, 4.66, 4.66, 4.64, 4.62, 4.59, 4.55, 4.50, 4.45, 4.39, &
 &  4.31, 4.23, 4.15, 4.05, 3.94, 3.83, 3.71, 3.58, 3.44, 3.29, 3.14, &
 &  2.97, 2.80, 2.62, 2.88, 3.47, 4.05, 4.62, 5.16, 5.70, 6.21, 6.71, &
 &  7.19, 7.66, 8.10, 8.54, 8.95, 9.35, 9.73 ];

      i = (p * 180 / PI) - C2limb2003May(0)
      do while (i >= 360.0)
	    i=i- 360.0
	  end do
      do while (i < 0.0) 
	    i=i+360.0
	  end do
      i = i / C2limb2003May(1)
      if (i >= (C2limb2003May(2) - 1)) then 
         q = 999.0;
		 goto 10
      end if
      n = floor(i)
      q = (((C2limb2003May(n+4) - C2limb2003May(n+3)) * (i - n)) + C2limb2003May(n+3))  
10   end subroutine limbcorrectionc2
!!------------------------------------------------------------------
   subroutine limbcorrectionc3(p, q)

      integer ::  n
	  real*8 :: i,C3limb2003May(0:342),p,q
	  real*8, parameter  :: PI =3.14151926

! C3 limb corrections for the 2003 May 31 annular eclipse in seconds
!
! The first 3 elements of the array tell us that the remaining data starts at contact angle
! 353.37 degrees in 0.4 degree increments for 340 data points, which means that the last
! element is for angle 128.97 degrees
!
! These limb corrections were calculated by Fred Espenak, NASA/GSFC

 C3limb2003May = [353.37, 0.4, 340.0, &  
 &  -8.67, -6.49, -6.14, -6.36, -6.62, -6.87, -7.09, -7.29, -7.47, -7.64, -7.78, &  
 &  -7.90, -8.01, -8.09, -8.15, -8.20, -8.22, -8.23, -8.22, -8.18, -8.13, -8.06, &  
 &  -7.97, -7.86, -7.74, -7.59, -7.43, -7.25, -7.05, -6.83, -6.60, -6.34, -6.07, &  
 &  -5.78, -5.48, -5.15, -4.81, -4.79, -4.82, -4.84, -4.85, -4.84, -4.83, -4.80, &  
 &  -4.76, -4.70, -4.64, -4.56, -4.48, -4.38, -4.26, -4.14, -4.01, -3.86, -3.81, &  
 &  -3.79, -3.78, -3.92, -4.05, -4.17, -4.29, -4.40, -4.50, -4.59, -4.67, -4.74, &  
 &  -4.81, -4.87, -4.92, -4.96, -5.00, -5.02, -5.04, -5.05, -5.05, -5.04, -5.03, &  
 &  -5.00, -4.97, -4.93, -4.89, -4.83, -4.77, -4.70, -4.62, -4.53, -4.44, -4.33, &  
 &  -4.22, -4.20, -4.26, -4.32, -4.37, -4.42, -4.45, -4.49, -4.51, -4.53, -4.54, &  
 &  -4.55, -4.54, -4.54, -4.52, -4.50, -4.47, -4.44, -4.39, -4.35, -4.29, -4.23, &  
 &  -4.16, -4.09, -4.01, -3.92, -3.82, -3.72, -3.62, -3.50, -3.38, -3.37, -3.40, &  
 &  -3.42, -3.44, -3.46, -3.47, -3.47, -3.46, -3.45, -3.44, -3.42, -3.39, -3.36, &  
 &  -3.32, -3.27, -3.22, -3.18, -3.12, -3.06, -3.00, -2.93, -2.91, -2.90, -2.87, &  
 &  -2.85, -2.81, -2.77, -2.73, -2.68, -2.62, -2.57, -2.52, -2.63, -2.75, -2.87, &  
 &  -2.99, -3.10, -3.21, -3.31, -3.40, -3.49, -3.58, -3.66, -3.73, -3.81, -3.87, &  
 &  -3.93, -3.99, -4.04, -4.08, -4.12, -4.16, -4.19, -4.21, -4.23, -4.25, -4.26, &  
 &  -4.26, -4.26, -4.25, -4.24, -4.23, -4.23, -4.23, -4.22, -4.21, -4.20, -4.18, &  
 &  -4.15, -4.12, -4.09, -4.05, -4.04, -4.09, -4.14, -4.19, -4.23, -4.26, -4.31, &  
 &  -4.37, -4.44, -4.49, -4.54, -4.59, -4.63, -4.67, -4.76, -4.88, -4.99, -5.11, &  
 &  -5.21, -5.31, -5.41, -5.50, -5.58, -5.66, -5.74, -5.81, -5.88, -5.94, -5.99, &  
 &  -6.04, -6.09, -6.13, -6.16, -6.20, -6.22, -6.24, -6.26, -6.27, -6.27, -6.28, &  
 &  -6.28, -6.27, -6.26, -6.24, -6.22, -6.19, -6.16, -6.12, -6.08, -6.03, -5.98, &  
 &  -5.92, -5.86, -5.79, -5.72, -5.64, -5.56, -5.47, -5.47, -5.45, -5.44, -5.41, &  
 &  -5.39, -5.35, -5.31, -5.27, -5.22, -5.17, -5.11, -5.05, -4.98, -4.90, -4.82, &  
 &  -4.74, -4.67, -4.60, -4.53, -4.46, -4.43, -4.61, -4.78, -4.95, -5.11, -5.27, &  
 &  -5.42, -5.57, -5.70, -5.84, -5.96, -6.08, -6.20, -6.30, -6.40, -6.50, -6.59, &  
 &  -6.67, -6.75, -6.82, -6.88, -6.96, -7.13, -7.29, -7.44, -7.58, -7.72, -7.86, &  
 &  -7.98, -8.10, -8.22, -8.32, -8.42, -8.52, -8.60, -8.68, -8.75, -8.82, -8.88, &  
 &  -8.93, -8.98, -9.02, -9.06, -9.08, -9.11, -9.12, -9.13, -9.13, -9.12, -9.11, &  
 &  -9.09, -9.06, -9.03, -8.98, -8.94, -8.90, -8.87, -8.83, -8.79, -8.74, -8.68, &  
 &  -8.62, -8.54, -8.47, -8.38, -8.29, -8.18, -8.09, -8.59, -9.08, -9.56 ];
      i = (p * 180 / PI) - C3limb2003May(0)
      do while (i >= 360.0)
	    i=i- 360.0
	  end do
      do while (i < 0.0) 
	    i=i+360.0
	  end do
      i = i / C3limb2003May(1)
      if (i >= (C3limb2003May(2) - 1)) then 
         q = 999.0; 
		 goto 10
      end if
      n = floor(i)
      q = (((C3limb2003May(n+4) - C3limb2003May(n+3)) * (i - n)) + C3limb2003May(n+3))  
 10  end subroutine limbcorrectionc3
!!==============================================================================================

!! Eclipse circumstances
!!  (0) Event type (C1=-2, C2=-1, Mid=0, C3=1, C4=2)
!!  (1) t
!! -- time-only dependent circumstances (and their per-hour derivatives) follow --
!!  (2) x  (3) y  (4)  d (5) dsin d  (6) dcos d  (7) mu  (8) l1  (9) l2  (10) dx
!! (11) dy (12) dd  (13) dmu  (14) dl1  (15) dl2
!! -- time and location dependent circumstances follow --
!! (16) h  (17) dsin h  (18) dcos h  (19) xi  (20) eta  (21) zeta  (22) dxi  (23) deta  (24) u
!! (25) v  (26) a  (27) b  (28) l1' (29) l2' (30) n^2
!! -- observational circumstances follow --
!! (31) p  (32) alt  (33) q  (34) v  (35) azi
!! (36) m (mid eclipse only) or limb correction applied (where available!)
!! (37) magnitude (mid eclipse only)
!! (38) moon/sun (mid eclipse only)
!! (39) calculated local event type for a transparent earth (0 = none, 1 = partial, 2 = annular, 3 = total)

!!------------------------------------------------------------------------
!! Populate the circumstances array with the time-only dependent circumstances (x, y, d, m, ...)
   subroutine timedependent(circumstances) 
	 real*8 circumstances(0:39)
	 real*8, parameter  :: PI =3.14151926
     real*8 :: t, ans

  t = circumstances(1)
  !! Calculate x
  ans = elements(9) * t + elements(8)
  ans = ans * t + elements(7)
  ans = ans * t + elements(6)
  circumstances(2) = ans
  !! Calculate dx
  ans = 3.0 * elements(9) * t + 2.0 * elements(8)
  ans = ans * t + elements(7)
  circumstances(10) = ans
  !! Calculate y
  ans = elements(13) * t + elements(12)
  ans = ans * t + elements(11)
  ans = ans * t + elements(10)
  circumstances(3) = ans
  !! Calculate dy
  ans = 3.0 * elements(13) * t + 2.0 * elements(12)
  ans = ans * t + elements(11)
  circumstances(11) = ans
  !! Calculate d
  ans = elements(16) * t + elements(15)
  ans = ans * t + elements(14)
  ans = ans * PI / 180.0
  circumstances(4) = ans
  !! dsin d and dcos d
  circumstances(5) = dsin(ans)
  circumstances(6) = dcos(ans)
  !! Calculate dd
  ans = 2.0 * elements(16) * t + elements(15)
  ans = ans * PI / 180.0
  circumstances(12) = ans
  !! Calculate m
  ans = elements(19) * t + elements(18)
  ans = ans * t + elements(17)
  if (ans >= 360.0) then 
    ans = ans - 360.0
  end if
  ans = ans * PI / 180.0
  circumstances(7) = ans
  !! Calculate dm
  ans = 2.0 * elements(19) * t + elements(18)
  ans = ans * PI / 180.0
  circumstances(13) = ans
  !! Calculate l1 and dl1
  typepe = circumstances(0)
  if ((typepe == -2) .OR. (typepe == 0) .OR. (typepe == 2)) then
    ans = elements(22) * t + elements(21)
    ans = ans * t + elements(20)
    circumstances(8) = ans
    circumstances(14) = 2.0 * elements(22) * t + elements(21)
  end if
  !! Calculate l2 and dl2
  if ((typepe == -1) .OR. (typepe == 0) .OR. (typepe == 1)) then
    ans = elements(25) * t + elements(24)
    ans = ans * t + elements(23)
    circumstances(9) = ans
    circumstances(15) = 2.0 * elements(25) * t + elements(24)
  end if
   end subroutine timedependent
!!====================================================================
!!------------------------------------------------------------------------
!! Populate the circumstances array with the time and location dependent circumstances
  subroutine timelocdependent(circumstances)
	 real*8, parameter  :: PI =3.14151926
	 real*8 :: circumstances(0:39)

  call timedependent(circumstances)
  !! Calculate h, dsin h, dcos h
  circumstances(16) = circumstances(7) - obsvconst(1)*PI/180.0 - (elements(5) / 13713.44)
  circumstances(17) = dsin(circumstances(16))
  circumstances(18) = dcos(circumstances(16))
  !! Calculate xi
  circumstances(19) = obsvconst(5) * circumstances(17)
  !! Calculate eta
  circumstances(20) = obsvconst(4) * circumstances(6) - obsvconst(5) * circumstances(18) * circumstances(5)
  !! Calculate zeta
  circumstances(21) = obsvconst(4) * circumstances(5) + obsvconst(5) * circumstances(18) * circumstances(6)
  !! Calculate dxi
  circumstances(22) = circumstances(13) * obsvconst(5) * circumstances(18)
  !! Calculate deta
  circumstances(23) = circumstances(13) * circumstances(19) * circumstances(5) - circumstances(21) * circumstances(12)
  !! Calculate u
  circumstances(24) = circumstances(2) - circumstances(19)
  !! Calculate v
  circumstances(25) = circumstances(3) - circumstances(20)
  !! Calculate a
  circumstances(26) = circumstances(10) - circumstances(22)
  !! Calculate b
  circumstances(27) = circumstances(11) - circumstances(23)
  !! Calculate l1'
  typepe = circumstances(0)
  if ((typepe == -2) .OR. (typepe == 0) .OR. (typepe == 2)) then
    circumstances(28) = circumstances(8) - circumstances(21) * elements(26)
  endif
  !! Calculate l2'
  if ((typepe == -1) .OR. (typepe == 0) .OR. (typepe == 1)) then
    circumstances(29) = circumstances(9) - circumstances(21) * elements(27)
  endif
  !! Calculate n^2
  circumstances(30) = circumstances(26) * circumstances(26) + circumstances(27) * circumstances(27)

 end subroutine timelocdependent
 !!==========================================================================================

!! Iterate on C1 or C4-----------------------------------------------------------------------
  subroutine c1c4iterate(circumstances) 
  real*8 :: circumstances(0:39) 
  real*8 :: sign, tmp, n,tt
  integer :: iter

  call timelocdependent(circumstances)
  if (circumstances(0) < 0) then
    sign=-1.0
  else 
    sign=1.0
  end if
  tmp=1.0
  iter=0;
  do while (((tmp > 0.000001) .OR. (tmp < -0.000001)) .AND. (iter < 50)) 
    n = sqrt(circumstances(30))
    tmp = circumstances(26) * circumstances(25) - circumstances(24) * circumstances(27)
    tmp = tmp / n / circumstances(28)
	tt=1.0 - tmp * tmp;
	if(tt<0) tt=0;
    tmp = sign * sqrt(tt) * circumstances(28) / n
    tmp = (circumstances(24) * circumstances(26) + circumstances(25) * circumstances(27)) / circumstances(30) - tmp
    circumstances(1) = circumstances(1) - tmp
    call timelocdependent(circumstances)
    iter = iter +1
  end do
  end subroutine c1c4iterate
!!======================================================================================
!! Get C1 and C4 data--------------------------------------------------------------------
!!   Entry conditions -
!!   1. The mid array must be populated
!!   2. The magnitude at mid eclipse must be > 0.0
  subroutine getc1c4() 
  real*8 :: tmp, n,tt
  n = sqrt(mid(30))
  tmp = mid(26) * mid(25) - mid(24) * mid(27)
  tmp = tmp / n / mid(28)
  tt=1.0 - tmp * tmp;
  if(tt<0) tt=0;
  tmp = sqrt(tt) * mid(28) / n
  c1(0) = -2;
  c4(0) = 2;
  c1(1) = mid(1) - tmp;
  c4(1) = mid(1) + tmp;
  call c1c4iterate(c1)
  call c1c4iterate(c4)
  end subroutine getc1c4

!!===============================================================

!! Iterate on C2 or C3--------------------------------------------------------
  subroutine c2c3iterate(circumstances)
    real*8 :: sign,  tmp, n,circumstances(0:39),tt
	integer :: iter

  call timelocdependent(circumstances)
  if (circumstances(0) < 0) then
    sign=-1.0
   else 
    sign=1.0
  endif
  if (mid(29) < 0.0) then
    sign = -sign
  endif
  tmp=1.0
  iter=0
  do while (((tmp > 0.000001) .OR. (tmp < -0.000001)) .AND. (iter < 50)) 
    n = sqrt(circumstances(30))
    tmp = circumstances(26) * circumstances(25) - circumstances(24) * circumstances(27)
    tmp = tmp / n / circumstances(29)
    tt=1.0 - tmp * tmp;
	if(tt<0) tt=0;
    tmp = sign * sqrt(tt) * circumstances(29) / n
    tmp = (circumstances(24) * circumstances(26) + circumstances(25) * circumstances(27)) / circumstances(30) - tmp
    circumstances(1) = circumstances(1) - tmp
    call timelocdependent(circumstances)
    iter= iter+1
  end do
end subroutine c2c3iterate
!!===============================================================================
!!--------------------------------------------------------------------------------
!! Get C2 and C3 data
!!   Entry conditions -
!!   1. The mid array must be populated
!!   2. There mut be either a total or annular eclipse at the location!
  subroutine getc2c3() 
  real*8 :: tmp, n, tt

  n = sqrt(mid(30))
  tmp = mid(26) * mid(25) - mid(24) * mid(27)
  tmp = tmp / n / mid(29)
  tt=1.0 - tmp * tmp;
  if(tt<0) tt=0;
  tmp = sqrt(tt) * mid(29) / n
  c2(0) = -1
  c3(0) = 1
  if (mid(29) < 0.0) then
    c2(1) = mid(1) + tmp
    c3(1) = mid(1) - tmp
   else 
    c2(1) = mid(1) - tmp
    c3(1) = mid(1) + tmp
  endif
  call c2c3iterate(c2)
  call c2c3iterate(c3)
  end subroutine getc2c3
!==============================================================================================
!----------------------------------------------------------------------------------------------
!! Get the observational circumstances
  subroutine observational(circumstances) 
  real*8, parameter  :: PI =3.14151926
  real*8 :: contacttype, dcoslat, dsinlat
  real*8 :: circumstances(0:39)
  !! We are looking at an "external" contact UNLESS this is a total eclipse AND we are looking at
  !! c2 or c3, in which case it is an INTERNAL contact! Note that if we are looking at mid eclipse,
  !! then we may not have determined the type of eclipse (mid[39]) just yet!
  if (circumstances(0) == 0) then
    contacttype = 1.0
   else if ((mid(39) == 3) .AND. ((circumstances(0) == -1) .OR. (circumstances(0) == 1))) then
      contacttype = -1.0
       else 
      contacttype = 1.0
   endif
  !! Calculate p
  circumstances(31) = atan2(contacttype*circumstances(24), contacttype*circumstances(25))
  !! Calculate alt
  dsinlat = dsin(obsvconst(0)*Pi/180.0)
  dcoslat = dcos(obsvconst(0)*Pi/180.0)
  circumstances(32) = dasin(circumstances(5) * dsinlat + circumstances(6) * dcoslat * circumstances(18))
  !! Calculate q
  circumstances(33) = dasin(dcoslat * circumstances(17) / dcos(circumstances(32)))
  if (circumstances(20) < 0.0) then
    circumstances(33) = PI - circumstances(33)
  endif
  !! Calculate v
  circumstances(34) = circumstances(31) - circumstances(33)
  !! Calculate azi
  circumstances(35) = atan2(-1.0*circumstances(17)*circumstances(6), circumstances(5)*dcoslat &
    &- circumstances(18)*dsinlat*circumstances(6))
  end subroutine observational
!!=======================================================================================
!!----------------------------------------------------------------------------------------
!! Calculate mid eclipse
   subroutine getmid() 
    real*8 :: tmp
    integer :: iter
    mid(0) = 0
    mid(1) = 0.0
    iter = 0
    tmp = 1.0
    call timelocdependent(mid)
    do while (((tmp > 0.000001) .OR. (tmp < -0.000001)) .AND. (iter < 50)) 
       tmp = (mid(24) * mid(26) + mid(25) * mid(27)) / mid(30)
       mid(1) = mid(1) - tmp
       iter = iter + 1
       call timelocdependent(mid)
    end do
  end subroutine getmid
!!====================================================================================
!! Populate the c1, c2, mid, c3 and c4 arrays
subroutine getall() 
  real*8 :: C2limb2003May,C3limb2003May
  call getmid()
  call observational(mid)
  !! Calculate m, magnitude and moon/sun
  mid(36) = sqrt(mid(24)*mid(24) + mid(25)*mid(25))
  mid(37) = (mid(28) - mid(36)) / (mid(28) + mid(29))
  mid(38) = (mid(28) - mid(29)) / (mid(28) + mid(29))
  if (mid(37) > 0.0) then
    call getc1c4()
    if ((mid(36) < mid(29)) .OR. (mid(36) < -mid(29))) then
      call getc2c3()
      if (mid(29) < 0.0) then
        mid(39) = 3 !! Total eclipse
       else 
        mid(39) = 2 !! Annular eclipse
      end if
      call observational(c2)
      call observational(c3)
      !! 2003 May 31 eclipse limb corrections -
      if (year==2003 .and. month==5 .and. day==31) then
        call limbcorrectionc2(c2(31), c2(36))
        call limbcorrectionc3(c3(31), c3(36))
        if (c2(36) < 990.0) then
          c2(1) = c2(1) + c2(36) / 3600.0
        endif
        if (c3(36) < 990.0) then
          c3(1) = c3(1) + c3(36) / 3600.0
        endif
      else 
        c2(36) = 999.9
        c3(36) = 999.9
      endif
    else 
      mid(39) = 1 !! Partial eclipse
    end if
    call observational(c1)
    call observational(c4)
 else 
    mid(39) = 0 !! No eclipse
 end if
end subroutine getall
!!===========================================================================================
!!----------------------Get the local time of an event---------------------------------------
  subroutine gettime(circumstances)
     real*8 :: circumstances(0:39)
!! Calculate the local time. Add 0.05 seconds, as we will be rounding up to the nearest 0.1 sec
      t = circumstances(1) + elements(1) - obsvconst(3) - (elements(4) - 0.05) / 3600.0
      if (t < 0.0) then
       t = t + 24.0
      endif
      if (t >= 24.0) then
       t = t - 24.0
      endif
 end subroutine gettime
!!===========================================================================================
!!----------------------------- Get the altitude----------------------------------------------
subroutine getalt(circumstances) 
  real*8 :: circumstances(0:39)
  real*8, parameter  :: PI =3.14151926
  alt = circumstances(32) * 180.0 / PI
end subroutine getalt
!!==========================================================================================
!!------------------------------Get the azimuth--------------------------------------------
subroutine getazi(circumstances) 
  real*8 :: circumstances(0:39)
  real*8, parameter  :: PI =3.14151926

  azi = circumstances(35) * 180.0 / PI
  if (azi < 0.0) then
    azi = azi + 360.0
  endif
  if (azi >= 360.0) then
    azi = azi - 360.0
  endif
end subroutine getazi
!!============================================================================================
!!------------------------------Get the P-----------------------------------------------------
subroutine getp(circumstances) 
  real*8 :: circumstances(0:39)
  real*8, parameter  :: PI =3.14151926
  p = circumstances(31) * 180.0 / PI
  if (p < 0.0) then
    p = p + 360.0
  endif
  if (p >= 360.0) then
    p = p - 360.0
  endif
end subroutine getp
!!===========================================================================================
!!---------------------------------------Get V-----------------------------------------------
subroutine getv(circumstances) 
  real*8 :: circumstances(0:39)
  real*8, parameter  :: PI =3.14151926
  V = floor(120.5 - circumstances(34) * 60.0 / PI) / 10.0
  if (V > 13.0) then
    V = V - 12.0
  endif
  if (V > 13.0) then
    V = V - 12.0
  endif
  if (V < 1.0) then
    V = V + 12.0
  endif
end subroutine getV
!!============================================================================================
!!-----------------------Get duration------------------------------------------------------
subroutine getduration() 
  duration=c3(1)-c2(1);
  if (duration<0.0) then
    duration=duration+24.0
  else if (duration >= 24.0) then
    duration=duration-24.0
  endif
end subroutine getduration
!!========================================================================================

!!----------------------------- Get the coverage-------------------------------------------
subroutine getcoverage() 
  real*8, parameter  :: PI =3.14151926
  real*8 :: a, b, c

  if (mid(37) <= 0.0) then
    cover=0
  end if
  if (mid(37) >= 1.0) then
    cover=1.0
  end if
  if (mid(39) == 2) then
    c = mid(38) * mid(38)
   else 
    c = dacos((mid(28)*mid(28) + mid(29)*mid(29) - 2.0*mid(36)*mid(36)) / (mid(28)*mid(28) - mid(29)*mid(29)))
    b = dacos((mid(28)*mid(29) + mid(36)*mid(36))/mid(36)/(mid(28)+mid(29)))
    a = PI - b - c
    c = ((mid(38)*mid(38)*a + b) - mid(38)*dsin(c))/PI
  end if
  cover = c+0.5/1000.0
end subroutine getcoverage
!--------------------------------------------------------------------------
  subroutine findEclipse(syear,smonth,sday)
	 integer :: iyear,imonth,iday,syear,smonth,sday,i
    iyear=1900;

     if (doFill) call fillEclipse()

     i = 1
     do while(iyear<=syear .and. i <= 155)
        iyear = int(eclipse_list(i,1))
        imonth = int(eclipse_list(i,2))
        iday = int(eclipse_list(i,3))
        elements(0:27) = eclipse_list(i,4:31)

	    if (iyear==syear .and. imonth==smonth .and. iday==sday) then
	       goto 20
        end if

        i = i+1
      end do
10    close(9)
      write(*,*)'there is no elipse at this time'
      return
!! Get the observer's geocentric position
20    close(9)
  end subroutine findEclipse
!--------------------------------------------------------------------------
  subroutine fillEclipse()
    integer :: i

    if (.not. allocated(eclipse_list)) allocate(eclipse_list(155,31))

    i = 1
    eclipse_list(i,:) = (/ 1970.,3.,7.,2440653.250,18.0,-3.0,3.0,40.4,40.4,-0.03259033,0.50545257,-0.2268E-04,-0.8157E-05,0.49146292,0.27499732,-0.7894E-06,-0.4643E-05,-5.23251677,0.01554100,0.1004E-05,87.22360229,15.00374413,0.1334E-05,0.53947592,0.00006010,-0.1271E-04,-0.00664825,0.00005980,-0.1264E-04,0.00471047,0.00468701 /)
    i = i + 1
    eclipse_list(i,:) = [ 1970.,8.,31.,2440830.500,22.0,-4.0,3.0,40.8,40.8,-0.21531861,0.44769949,-0.1163E-04,-0.4932E-05,-0.49274859,-0.23453267,-0.2417E-04,0.2747E-05,8.54803467,-0.01448508,-0.1905E-05,149.93501282,15.00432968,0.1094E-05,0.56760037,-0.00001617,-0.9730E-05,0.02133636,-0.00001609,-0.9681E-05,0.00463379,0.00461071 ]
    i = i + 1
    eclipse_list(i,:) = [  1971.,2.,25.,2441008.000,10.0,-3.0,2.0,41.3,41.3,-0.32151362,0.51789355,-0.4075E-05,-0.8654E-05,1.09238625,0.26550746,-0.5096E-05,-0.4655E-05,-9.26003456,0.01480197,0.2003E-05,326.68786621,15.00288391,0.2096E-05,0.53718776,-0.00004126,-0.1303E-04,-0.00892498,-0.00004105,-0.1297E-04,0.00472292,0.00469940 ]
    i = i + 1
    eclipse_list(i,:) = [ 1971.,7.,22.,2441155.000,10.0,-2.0,1.0,41.8,41.8,0.68554115,0.50543481,-0.8301E-04,-0.6387E-05,1.37126946,-0.15726262,-0.1810E-03,0.2156E-05,20.36965942,-0.00779889,-0.4491E-05,328.41119385,15.00082302,0.1827E-05,0.55454224,0.00010073,-0.1054E-04,0.00834334,0.00010023,-0.1048E-04,0.00460218,0.00457926 ]
    i = i + 1
    eclipse_list(i,:) = [ 1971.,8.,20.,2441184.500,23.0,-3.0,2.0,41.8,41.8,-0.38978973,0.46392810,-0.2562E-04,-0.5405E-05,-1.21715128,-0.22285873,-0.1140E-04,0.2742E-05,12.42899895,-0.01321412,-0.2800E-05,164.15605164,15.00361347,0.1640E-05,0.56232190,0.00007557,-0.1005E-04,0.01608425,0.00007519,-0.1000E-04,0.00462199,0.00459897 ]
    i = i + 1
    eclipse_list(i,:) = [ 1972.,1.,16.,2441333.000,11.0,-3.0,3.0,42.2,42.3,0.23215397,0.52077079,-0.3423E-04,-0.7120E-05,-0.90778792,0.15136743,0.1880E-03,-0.2203E-05,-21.05788994,0.00739399,0.5356E-05,342.60723877,14.99772739,0.2309E-05,0.55650145,-0.00012142,-0.1138E-04,0.01029241,-0.00012082,-0.1132E-04,0.00475380,0.00473012 ]
    i = i + 1
    eclipse_list(i,:) = [ 1972.,7.,10.,2441509.250,20.0,-3.0,3.0,43.2,42.8,0.28180188,0.54694635,-0.6469E-04,-0.8288E-05,0.63913912,-0.13100556,-0.1790E-03,0.2150E-05,22.14701271,-0.00513984,-0.5069E-05,118.65975952,14.99997902,0.1392E-05,0.53948885,0.00010210,-0.1186E-04,-0.00663504,0.00010160,-0.1180E-04,0.00459882,0.00457592 ]
    i = i + 1
    eclipse_list(i,:) = [ 1973.,1.,4.,2441687.250,16.0,-4.0,3.0,44.2,43.4,0.16819540,0.50076556,-0.2606E-04,-0.5809E-05,-0.23490599,0.10499598,0.1607E-03,-0.1339E-05,-22.68856621,0.00422709,0.6098E-05,58.72224426,14.99662018,0.1436E-05,0.57128686,-0.00008145,-0.1019E-04,0.02500416,-0.00008104,-0.1014E-04,0.00475638,0.00473270 ]
    i = i + 1
    eclipse_list(i,:) = [ 1973.,6.,30.,2441864.000,12.0,-3.0,3.0,44.2,44.0,0.19178969,0.57546759,-0.4074E-04,-0.9741E-05,-0.11111309,-0.09466382,-0.1580E-03,0.1747E-05,23.16779137,-0.00232936,-0.5477E-05,359.10876465,14.99943066,0.7414E-06,0.53061301,0.00001900,-0.1277E-04,-0.01546669,0.00001890,-0.1270E-04,0.00459825,0.00457535 ]
    i = i + 1
    eclipse_list(i,:) = [ 1973.,12.,24.,2442041.250,15.0,-3.0,4.0,44.2,44.5,-0.07312115,0.49988997,-0.2011E-04,-0.5556E-05,0.41127804,0.06085132,0.1394E-03,-0.7894E-06,-23.41661263,0.00070772,0.6460E-05,45.08232880,14.99615288,0.1331E-06,0.57518399,0.00002672,-0.9898E-05,0.02888189,0.00002659,-0.9849E-05,0.00475528,0.00473160 ]
    i = i + 1
    eclipse_list(i,:) = [ 1974.,6.,20.,2442218.750,5.0,-3.0,3.0,45.2,45.0,0.04381487,0.57430005,-0.5640E-05,-0.9404E-05,-0.83066511,-0.04926458,-0.1205E-03,0.9398E-06,23.43549919,0.00056747,-0.5681E-05,254.65780640,14.99921703,-0.2617E-07,0.53341234,-0.00007637,-0.1248E-04,-0.01268134,-0.00007599,-0.1242E-04,0.00460061,0.00457770 ]
    i = i + 1
    eclipse_list(i,:) = [ 1974.,12.,13.,2442395.250,16.0,-2.0,3.0,45.2,45.4,-0.15130030,0.52346957,-0.1601E-04,-0.6510E-05,1.07530773,0.01739977,0.1145E-03,-0.3063E-06,-23.15630150,-0.00282684,0.6370E-05,61.45851517,14.99647713,-0.1185E-05,0.56478947,0.00011429,-0.1066E-04,0.01853920,0.00011372,-0.1060E-04,0.00475049,0.00472683 ]
    i = i + 1
    eclipse_list(i,:) = [ 1975.,5.,11.,2442543.750,7.0,-2.0,3.0,46.2,45.9,-0.33642676,0.52215034,0.4953E-04,-0.6700E-05,1.02199304,0.09154182,-0.1504E-03,-0.1122E-05,17.73483467,0.01059320,-0.3929E-05,285.90759277,15.00166416,-0.1912E-05,0.55563033,-0.00010870,-0.1063E-04,0.00942592,-0.00010816,-0.1058E-04,0.00462997,0.00460691 ]
    i = i + 1
    eclipse_list(i,:) = [ 1975.,11.,3.,2442720.000,13.0,-2.0,3.0,46.2,46.3,-0.37404397,0.55952728,0.3073E-04,-0.9079E-05,-0.96607953,-0.12643839,0.1574E-03,0.1949E-05,-14.96477413,-0.01276184,0.3133E-05,19.09655762,15.00145340,-0.2801E-05,0.54039729,0.00008898,-0.1262E-04,-0.00573143,0.00008854,-0.1256E-04,0.00471334,0.00468987 ]
    i = i + 1
    eclipse_list(i,:) = [ 1976.,4.,29.,2442898.000,10.0,-3.0,4.0,47.2,46.8,-0.27503377,0.49430403,0.3464E-04,-0.5586E-05,0.28391948,0.11243641,-0.9296E-04,-0.1207E-05,14.56363678,0.01261925,-0.3242E-05,330.66961670,15.00265312,-0.1773E-05,0.56741798,-0.00003436,-0.9799E-05,0.02115488,-0.00003419,-0.9750E-05,0.00464281,0.00461969 ]
    i = i + 1
    eclipse_list(i,:) = [ 1976.,10.,23.,2443074.750,5.0,-3.0,3.0,47.2,47.3,-0.21339460,0.56328297,0.3373E-04,-0.9577E-05,-0.28150886,-0.15026081,0.9024E-04,0.2470E-05,-11.44277287,-0.01423271,0.2278E-05,258.90756226,15.00283051,-0.2332E-05,0.53576642,-0.00001092,-0.1302E-04,-0.01033920,-0.00001087,-0.1295E-04,0.00469991,0.00467650 ]
    i = i + 1
    eclipse_list(i,:) = [ 1977.,4.,18.,2443252.000,11.0,-4.0,3.0,48.2,47.8,0.33826455,0.49262398,-0.4068E-05,-0.5675E-05,-0.32184276,0.13328393,-0.5049E-04,-0.1474E-05,10.86802959,0.01420036,-0.2466E-05,345.16491699,15.00354385,-0.1384E-05,0.56645668,0.00004991,-0.9956E-05,0.02019831,0.00004967,-0.9906E-05,0.00465662,0.00463343 ]
    i = i + 1
    eclipse_list(i,:) = [ 1977.,10.,12.,2443429.250,20.0,-3.0,4.0,48.2,48.3,-0.13762422,0.54020298,0.3324E-04,-0.8306E-05,0.44132689,-0.16101514,0.1801E-04,0.2389E-05,-7.58560658,-0.01522500,0.1461E-05,123.38939667,15.00387192,-0.1741E-05,0.54333234,-0.00009004,-0.1214E-04,-0.00281097,-0.00008959,-0.1207E-04,0.00468615,0.00466281 ]
    i = i + 1
    eclipse_list(i,:) = [ 1978.,4.,7.,2443606.250,15.0,-2.0,3.0,49.2,48.9,0.28864408,0.51485884,-0.1756E-04,-0.6887E-05,-1.07040870,0.15584747,0.3422E-05,-0.2035E-05,6.84699249,0.01528305,-0.1684E-05,44.45822525,15.00419140,-0.8395E-06,0.55398393,0.00011724,-0.1103E-04,0.00778766,0.00011666,-0.1098E-04,0.00467059,0.00464733 ]
    i = i + 1
    eclipse_list(i,:) = [ 1978.,10.,2.,2443783.750,6.0,-2.0,3.0,49.2,49.3,0.11055008,0.50530893,0.1458E-04,-0.6505E-05,1.18390203,-0.16086034,-0.3930E-04,0.1985E-05,-3.44974113,-0.01575957,0.6480E-06,272.61972046,15.00454426,-0.1041E-05,0.55768645,-0.00009952,-0.1074E-04,0.01147172,-0.00009902,-0.1069E-04,0.00467218,0.00464891 ]
    i = i + 1
    eclipse_list(i,:) = [ 1979.,2.,26.,2443931.250,17.0,-3.0,3.0,50.2,49.8,-0.20600221,0.55702728,-0.1022E-04,-0.9348E-05,0.87545109,0.16251618,0.3907E-05,-0.2634E-05,-8.75745583,0.01517580,0.1864E-05,71.75366211,15.00310802,0.1961E-05,0.53790218,0.00004136,-0.1295E-04,-0.00821413,0.00004115,-0.1289E-04,0.00472157,0.00469805 ]
    i = i + 1
    eclipse_list(i,:) = [ 1979.,8.,22.,2444108.250,17.0,-3.0,3.0,50.2,50.2,-0.42897618,0.48865560,-0.3122E-05,-0.5468E-05,-0.88319886,-0.12830278,-0.2930E-04,0.1372E-05,11.82268143,-0.01367850,-0.2628E-05,74.26847076,15.00381565,0.1589E-05,0.56708407,0.00002143,-0.9707E-05,0.02082269,0.00002132,-0.9658E-05,0.00462399,0.00460096 ]
    i = i + 1
    eclipse_list(i,:) = [ 1980.,2.,16.,2444285.750,9.0,-3.0,3.0,51.2,50.7,0.00103041,0.56041682,-0.1246E-04,-0.9307E-05,0.22978656,0.14254689,0.7345E-04,-0.2280E-05,-12.56627846,0.01400729,0.2736E-05,311.46118164,15.00194931,0.2509E-05,0.53934973,-0.00006218,-0.1289E-04,-0.00677381,-0.00006187,-0.1283E-04,0.00473271,0.00470914 ]
    i = i + 1
    eclipse_list(i,:) = [ 1980.,8.,10.,2444462.250,19.0,-3.0,4.0,51.2,51.0,-0.14559500,0.51036757,-0.2789E-04,-0.6268E-05,-0.16445617,-0.11048350,-0.7825E-04,0.1288E-05,15.35354233,-0.01204436,-0.3360E-05,103.69969177,15.00288391,0.1913E-05,0.55838150,0.00010527,-0.1031E-04,0.01216348,0.00010474,-0.1026E-04,0.00461406,0.00459108 ]
    i = i + 1
    eclipse_list(i,:) = [ 1981.,2.,4.,2444640.500,22.0,-3.0,3.0,51.2,51.5,0.01344806,0.53927886,-0.5695E-05,-0.7751E-05,-0.49118471,0.11136021,0.1351E-03,-0.1507E-05,-16.03494835,0.01230621,0.3656E-05,146.50033569,15.00048733,0.2838E-05,0.55204868,-0.00011762,-0.1171E-04,0.00586185,-0.00011704,-0.1165E-04,0.00474236,0.00471874 ]
    i = i + 1
    eclipse_list(i,:) = [ 1981.,7.,31.,2444816.750,4.0,-3.0,3.0,52.2,51.8,0.21525431,0.54793638,-0.4561E-04,-0.8034E-05,0.55176628,-0.08917374,-0.1411E-03,0.1241E-05,18.30328941,-0.01003576,-0.3987E-05,238.42449951,15.00187397,0.2028E-05,0.54301703,0.00011136,-0.1158E-04,-0.00312449,0.00011080,-0.1152E-04,0.00460618,0.00458324 ]
    i = i + 1
    eclipse_list(i,:) = [ 1982.,1.,25.,2444994.750,5.0,-3.0,2.0,52.2,52.2,0.32500023,0.51230949,-0.1698E-04,-0.6172E-05,-1.19656229,0.07533986,0.1772E-03,-0.8139E-06,-19.03972244,0.01001413,0.4602E-05,251.94006348,14.99889755,0.2778E-05,0.56794775,-0.00010186,-0.1040E-04,0.02168172,-0.00010136,-0.1035E-04,0.00475014,0.00472648 ]
    i = i + 1
    eclipse_list(i,:) = [ 1982.,6.,21.,2445142.000,12.0,-2.0,2.0,52.2,52.6,0.07337533,0.58126396,0.1214E-04,-0.9886E-05,-1.20875442,0.05667380,-0.9915E-04,-0.1094E-05,23.44344330,-0.00006997,-0.5732E-05,359.59475708,14.99923611,0.1250E-06,0.53016573,0.00000105,-0.1282E-04,-0.01591175,0.00000105,-0.1275E-04,0.00460010,0.00457719 ]
    i = i + 1
    eclipse_list(i,:) = [ 1982.,7.,20.,2445171.250,19.0,-2.0,2.0,53.2,52.6,0.27990681,0.57789016,-0.3696E-04,-0.9681E-05,1.26646936,-0.06000159,-0.2132E-03,0.9256E-06,20.61553574,-0.00772004,-0.4495E-05,103.43074036,15.00088024,0.1878E-05,0.53173220,0.00004363,-0.1266E-04,-0.01435306,0.00004341,-0.1259E-04,0.00460088,0.00457797 ]
    i = i + 1
    eclipse_list(i,:) = [ 1982.,12.,15.,2445319.000,10.0,-3.0,2.0,53.2,52.9,0.40209046,0.50094628,0.4707E-06,-0.5655E-05,1.08119714,-0.07597227,0.1116E-03,0.9831E-06,-23.26393318,-0.00201910,0.6436E-05,331.25390625,14.99632454,-0.8990E-06,0.57330143,0.00004606,-0.1001E-04,0.02700871,0.00004583,-0.9964E-05,0.00475184,0.00472817 ]
    i = i + 1
    eclipse_list(i,:) = [ 1983.,6.,11.,2445496.750,5.0,-3.0,3.0,53.2,53.4,0.24025983,0.56096524,0.2993E-04,-0.8859E-05,-0.45982304,0.09958333,-0.1340E-03,-0.1707E-05,23.04493523,0.00283416,-0.5538E-05,255.15377808,14.99935913,-0.6813E-06,0.53647351,-0.00009957,-0.1220E-04,-0.00963538,-0.00009907,-0.1214E-04,0.00460455,0.00458162 ]
    i = i + 1
    eclipse_list(i,:) = [ 1983.,12.,4.,2445673.000,13.0,-4.0,3.0,54.2,53.7,0.34178901,0.51896471,0.3853E-05,-0.6775E-05,0.33128408,-0.12362552,0.1372E-03,0.1766E-05,-22.20571709,-0.00537030,0.5880E-05,17.49466324,14.99718857,-0.2011E-05,0.56005812,0.00010999,-0.1102E-04,0.01383140,0.00010944,-0.1097E-04,0.00474486,0.00472123 ]
    i = i + 1
    eclipse_list(i,:) = [ 1984.,5.,30.,2445851.250,17.0,-4.0,3.0,54.2,54.0,0.05601003,0.52088910,0.5115E-04,-0.6906E-05,0.29867074,0.13301232,-0.1506E-03,-0.1918E-05,21.86890030,0.00569501,-0.5138E-05,75.62094116,14.99985313,-0.1315E-05,0.55110812,-0.00012498,-0.1088E-04,0.00492632,-0.00012436,-0.1083E-04,0.00461163,0.00458866 ]
    i = i + 1
    eclipse_list(i,:) = [ 1984.,11.,22.,2446027.500,23.0,-3.0,3.0,54.2,54.3,-0.04324950,0.54060394,0.3248E-04,-0.8425E-05,-0.31485939,-0.17163493,0.1592E-03,0.2844E-05,-20.32706833,-0.00828826,0.5025E-05,168.42231750,14.99853134,-0.2661E-05,0.54431510,0.00009636,-0.1240E-04,-0.00183322,0.00009588,-0.1234E-04,0.00473529,0.00471171 ]
    i = i + 1
    eclipse_list(i,:) = [ 1985.,5.,19.,2446205.500,21.0,-2.0,3.0,54.2,54.6,-0.57545817,0.48482391,0.6909E-04,-0.5532E-05,0.93886423,0.16021264,-0.1524E-03,-0.1978E-05,19.89410591,0.00840439,-0.4529E-05,135.89306641,15.00068855,-0.1761E-05,0.56449181,-0.00005251,-0.9881E-05,0.01824332,-0.00005225,-0.9832E-05,0.00462127,0.00459825 ]
    i = i + 1
    eclipse_list(i,:) = [ 1985.,11.,12.,2446382.000,14.0,-2.0,3.0,55.2,54.8,-0.45631406,0.54434079,0.7062E-04,-0.9204E-05,-0.87380719,-0.21001159,0.1700E-03,0.3742E-05,-17.77673531,-0.01071540,0.4029E-05,33.95781708,15.00004387,-0.2838E-05,0.53679305,0.00001187,-0.1309E-04,-0.00931778,0.00001181,-0.1302E-04,0.00472405,0.00470053 ]
    i = i + 1
    eclipse_list(i,:) = [ 1986.,4.,9.,2446529.750,6.0,-2.0,3.0,55.2,55.0,0.34347934,0.45629576,-0.6408E-05,-0.5336E-05,-1.04254282,0.24137661,-0.4147E-07,-0.3006E-05,7.48054028,0.01493897,-0.1812E-05,269.57363892,15.00397682,-0.9154E-06,0.56416452,0.00008974,-0.1020E-04,0.01791754,0.00008930,-0.1014E-04,0.00466905,0.00464580 ]
    i = i + 1
    eclipse_list(i,:) = [ 1986.,10.,3.,2446707.250,19.0,-3.0,3.0,55.2,55.2,0.42664278,0.48638743,0.5560E-05,-0.7032E-05,0.89859253,-0.26658231,-0.2943E-04,0.4025E-05,-4.06904554,-0.01545765,0.7723E-06,107.75020599,15.00441647,-0.1153E-05,0.54660767,-0.00011005,-0.1172E-04,0.00044811,-0.00010951,-0.1166E-04,0.00467360,0.00465033 ]
    i = i + 1
    eclipse_list(i,:) = [ 1987.,3.,29.,2446884.000,13.0,-3.0,3.0,55.2,55.5,0.22913754,0.47956255,-0.1677E-04,-0.6666E-05,-0.22226425,0.26463550,-0.1341E-04,-0.3878E-05,3.29555345,0.01560403,-0.8768E-06,13.77047348,15.00431061,-0.1970E-06,0.55027473,0.00011283,-0.1146E-04,0.00409688,0.00011227,-0.1140E-04,0.00468337,0.00466005 ]
    i = i + 1
    eclipse_list(i,:) = [ 1987.,9.,23.,2447061.750,3.0,-3.0,4.0,55.2,55.7,0.04201295,0.45559466,0.1015E-04,-0.5497E-05,0.29572821,-0.25436282,-0.1888E-04,0.3243E-05,0.17404810,-0.01560810,-0.1801E-06,226.83702087,15.00475883,-0.3808E-06,0.56109357,-0.00009196,-0.1039E-04,0.01486192,-0.00009151,-0.1034E-04,0.00465969,0.00463649 ]
    i = i + 1
    eclipse_list(i,:) = [ 1988.,3.,18.,2447238.500,2.0,-3.0,3.0,56.2,55.9,-0.19468223,0.50395757,-0.6478E-06,-0.8187E-05,0.37093231,0.28034675,-0.1931E-04,-0.4760E-05,-0.91812879,0.01579716,0.6044E-07,207.96788025,15.00425816,0.6059E-06,0.53825742,0.00006378,-0.1273E-04,-0.00786061,0.00006347,-0.1267E-04,0.00469703,0.00467364 ]
    i = i + 1
    eclipse_list(i,:) = [ 1988.,9.,11.,2447415.750,5.0,-4.0,3.0,56.2,56.1,-0.11086277,0.44259563,-0.5833E-05,-0.4860E-05,-0.47318095,-0.24311568,-0.4106E-05,0.2836E-05,4.48111725,-0.01527447,-0.1095E-05,255.84449768,15.00472450,0.4079E-06,0.56850761,-0.00001323,-0.9737E-05,0.02223908,-0.00001316,-0.9688E-05,0.00464598,0.00462284 ]
    i = i + 1
    eclipse_list(i,:) = [ 1989.,3.,7.,2447593.250,18.0,-2.0,3.0,56.2,56.4,-0.59675008,0.51151329,0.2298E-04,-0.8529E-05,0.92568552,0.27690777,-0.2172E-04,-0.4837E-05,-5.08024645,0.01555584,0.1036E-05,87.25234222,15.00377750,0.1383E-05,0.53674406,-0.00003365,-0.1298E-04,-0.00936647,-0.00003348,-0.1292E-04,0.00471024,0.00468678 ]
    i = i + 1
    eclipse_list(i,:) = [ 1989.,8.,31.,2447769.750,6.0,-3.0,2.0,56.2,56.6,-0.33567983,0.45761248,-0.1956E-04,-0.5356E-05,-1.17003071,-0.23831227,0.8664E-05,0.2941E-05,8.64246845,-0.01444713,-0.2002E-05,269.90679932,15.00434589,0.1092E-05,0.56235063,0.00007804,-0.1012E-04,0.01611279,0.00007765,-0.1007E-04,0.00463289,0.00460981 ]
    i = i + 1
    eclipse_list(i,:) = [ 1990.,1.,26.,2447918.250,20.0,-3.0,2.0,57.2,56.9,0.56778508,0.50742018,-0.5084E-04,-0.6860E-05,-0.79894918,0.18643744,0.1613E-03,-0.2664E-05,-18.62072372,0.01014584,0.4480E-05,116.84642029,14.99911308,0.2737E-05,0.55697238,-0.00013388,-0.1130E-04,0.01076099,-0.00013321,-0.1124E-04,0.00474887,0.00472522 ]
    i = i + 1
    eclipse_list(i,:) = [ 1990.,7.,22.,2448094.750,3.0,-3.0,3.0,57.2,57.2,0.19995058,0.53828216,-0.6685E-04,-0.8231E-05,0.73352063,-0.16929638,-0.1664E-03,0.2766E-05,20.34538078,-0.00779849,-0.4521E-05,223.39602661,15.00091457,0.1856E-05,0.53880906,0.00010617,-0.1195E-04,-0.00731149,0.00010565,-0.1189E-04,0.00460177,0.00457885 ]
    i = i + 1
    eclipse_list(i,:) = [ 1991.,1.,15.,2448272.500,0.0,-4.0,3.0,58.2,57.6,0.12669675,0.49010640,-0.3299E-04,-0.5625E-05,-0.24713154,0.14324532,0.1457E-03,-0.1781E-05,-21.07115364,0.00739146,0.5462E-05,177.62878418,14.99764156,0.2325E-05,0.57169992,-0.00007536,-0.1015E-04,0.02541518,-0.00007498,-0.1010E-04,0.00475431,0.00473064 ]
    i = i + 1
    eclipse_list(i,:) = [ 1991.,7.,11.,2448449.250,19.0,-3.0,3.0,58.2,58.0,-0.06725452,0.56713754,-0.3692E-04,-0.9612E-05,0.01212629,-0.13797723,-0.1529E-03,0.2496E-05,22.09743500,-0.00522018,-0.5113E-05,103.63541412,15.00005627,0.1494E-05,0.53043801,0.00002524,-0.1279E-04,-0.01564084,0.00002511,-0.1272E-04,0.00459868,0.00457578 ]
    i = i + 1
    eclipse_list(i,:) = [ 1992.,1.,4.,2448626.500,23.0,-3.0,4.0,58.2,58.3,-0.12921846,0.49348110,-0.3039E-04,-0.5479E-05,0.39099392,0.10225838,0.1321E-03,-0.1265E-05,-22.72167015,0.00413210,0.6190E-05,163.76997375,14.99658298,0.1347E-05,0.57505029,0.00003213,-0.9914E-05,0.02874886,0.00003197,-0.9865E-05,0.00475634,0.00473266 ]
    i = i + 1
    eclipse_list(i,:) = [ 1992.,6.,30.,2448804.000,12.0,-3.0,3.0,58.2,58.7,-0.23039281,0.56738216,-0.3849E-05,-0.9207E-05,-0.72325015,-0.09407896,-0.1228E-03,0.1672E-05,23.14090157,-0.00244503,-0.5551E-05,359.07891846,14.99942684,0.8676E-06,0.53387910,-0.00007127,-0.1242E-04,-0.01221687,-0.00007092,-0.1235E-04,0.00459843,0.00457553 ]
    i = i + 1
    eclipse_list(i,:) = [ 1992.,12.,24.,2448980.500,1.0,-3.0,2.0,59.2,59.1,0.12051979,0.52160698,-0.3970E-04,-0.6541E-05,1.09276879,0.06171783,0.1141E-03,-0.8767E-06,-23.41921997,0.00064265,0.6515E-05,195.11341858,14.99621201,0.9969E-07,0.56435484,0.00010116,-0.1073E-04,0.01810671,0.00010066,-0.1067E-04,0.00475461,0.00473093 ]
    i = i + 1
    eclipse_list(i,:) = [ 1993.,5.,21.,2449129.000,14.0,-2.0,3.0,59.2,59.5,-0.30413148,0.52511007,0.4287E-04,-0.6671E-05,1.11016607,0.05911719,-0.1739E-03,-0.6939E-06,20.25684738,0.00822132,-0.4589E-05,30.85806656,15.00062275,-0.1776E-05,0.55588204,-0.00010481,-0.1054E-04,0.00967644,-0.00010429,-0.1049E-04,0.00461953,0.00459652 ]
    i = i + 1
    eclipse_list(i,:) = [ 1993.,11.,13.,2449305.500,22.0,-3.0,2.0,60.2,59.9,-0.03945032,0.56676370,0.1336E-04,-0.9248E-05,-1.04917336,-0.09547210,0.1879E-03,0.1457E-05,-18.14668655,-0.01068719,0.4127E-05,153.89949036,14.99989796,-0.2895E-05,0.54058981,0.00007179,-0.1270E-04,-0.00553991,0.00007143,-0.1264E-04,0.00472565,0.00470211 ]
    i = i + 1
    eclipse_list(i,:) = [ 1994.,5.,10.,2449483.250,17.0,-3.0,4.0,60.2,60.3,-0.17341146,0.49906293,0.2956E-04,-0.5629E-05,0.38365141,0.08693927,-0.1183E-03,-0.9145E-06,17.68613052,0.01064192,-0.3994E-05,75.90597534,15.00162029,-0.1957E-05,0.56693691,-0.00003184,-0.9762E-05,0.02067623,-0.00003169,-0.9714E-05,0.00463085,0.00460779 ]
    i = i + 1
    eclipse_list(i,:) = [ 1994.,11.,3.,2449660.000,14.0,-3.0,3.0,61.2,60.6,0.11256027,0.56878257,0.2071E-04,-0.9657E-05,-0.38557220,-0.12578036,0.1233E-03,0.2049E-05,-15.10090446,-0.01268595,0.3262E-05,34.10152817,15.00142097,-0.2801E-05,0.53662097,-0.00003135,-0.1303E-04,-0.00948892,-0.00003119,-0.1296E-04,0.00471339,0.00468992 ]
    i = i + 1
    eclipse_list(i,:) = [ 1995.,4.,29.,2449837.250,18.0,-4.0,3.0,61.2,61.1,0.29701722,0.49877352,-0.7098E-06,-0.5786E-05,-0.27902138,0.11401000,-0.7489E-04,-0.1259E-05,14.48232937,0.01265913,-0.3271E-05,90.65514374,15.00268173,-0.1808E-05,0.56509030,0.00005580,-0.9978E-05,0.01883878,0.00005552,-0.9929E-05,0.00464336,0.00462024 ]
    i = i + 1
    eclipse_list(i,:) = [ 1995.,10.,24.,2450014.750,5.0,-4.0,3.0,61.2,61.5,0.33003256,0.54305160,0.1812E-04,-0.8268E-05,0.27634686,-0.14410999,0.5246E-04,0.2111E-05,-11.58056641,-0.01419584,0.2412E-05,258.92681885,15.00273991,-0.2360E-05,0.54477739,-0.00011489,-0.1209E-04,-0.00137312,-0.00011432,-0.1203E-04,0.00470036,0.00467695 ]
    i = i + 1
    eclipse_list(i,:) = [ 1996.,4.,17.,2450191.500,23.0,-3.0,2.0,62.2,61.8,0.46805340,0.52094173,-0.1941E-04,-0.7060E-05,-0.96886700,0.14228961,-0.2448E-04,-0.1877E-05,10.83629036,0.01419081,-0.2523E-05,165.15081787,15.00361347,-0.1409E-05,0.55229938,0.00010841,-0.1110E-04,0.00611154,0.00010787,-0.1104E-04,0.00465648,0.00463329 ]
    i = i + 1
    eclipse_list(i,:) = [ 1996.,10.,12.,2450369.000,14.0,-3.0,3.0,62.2,62.1,0.29616871,0.50603628,0.1452E-04,-0.6435E-05,1.08322012,-0.15152179,-0.1020E-04,0.1843E-05,-7.63950300,-0.01523369,0.1547E-05,33.40553284,15.00378227,-0.1770E-05,0.55937213,-0.00010672,-0.1069E-04,0.01314900,-0.00010619,-0.1064E-04,0.00468661,0.00466327 ]
    i = i + 1
    eclipse_list(i,:) = [ 1997.,3.,9.,2450516.500,1.0,-2.0,3.0,62.2,62.4,-0.50502318,0.55430734,0.1337E-04,-0.9343E-05,0.80381435,0.17435855,-0.2039E-04,-0.2846E-05,-4.55016375,0.01586202,0.8992E-06,192.33589172,15.00396252,0.1259E-05,0.53691423,0.00004861,-0.1296E-04,-0.00919715,0.00004836,-0.1289E-04,0.00470880,0.00468535 ]
    i = i + 1
    eclipse_list(i,:) = [ 1997.,9.,2.,2450693.500,0.0,-3.0,3.0,63.2,62.7,-0.33188948,0.48463467,-0.2565E-05,-0.5425E-05,-0.98140365,-0.14306034,-0.3121E-05,0.1538E-05,7.98437691,-0.01483813,-0.1845E-05,180.05313110,15.00450420,0.1001E-05,0.56766528,0.00002239,-0.9729E-05,0.02140094,0.00002228,-0.9681E-05,0.00463496,0.00461188 ]
    i = i + 1
    eclipse_list(i,:) = [ 1998.,2.,26.,2450871.250,17.0,-3.0,4.0,63.2,63.1,-0.33883452,0.55449575,0.1027E-04,-0.9172E-05,0.15050644,0.16113159,0.4941E-04,-0.2582E-05,-8.60943222,0.01521536,0.1724E-05,71.77923584,15.00313091,0.1953E-05,0.53918231,-0.00005214,-0.1282E-04,-0.00694038,-0.00005188,-0.1276E-04,0.00472161,0.00469809 ]
    i = i + 1
    eclipse_list(i,:) = [ 1998.,8.,22.,2451047.500,2.0,-3.0,4.0,63.2,63.3,-0.12741752,0.50677437,-0.2696E-04,-0.6278E-05,-0.24000636,-0.13235684,-0.5198E-04,0.1570E-05,11.89284897,-0.01363701,-0.2592E-05,209.24638367,15.00384808,0.1520E-05,0.55812109,0.00010754,-0.1038E-04,0.01190435,0.00010700,-0.1033E-04,0.00462302,0.00460000 ]
    i = i + 1
    eclipse_list(i,:) = [ 1999.,2.,16.,2451225.750,7.0,-4.0,3.0,64.2,63.5,0.34196952,0.53182095,-0.1825E-04,-0.7574E-05,-0.40032700,0.13612878,0.1049E-03,-0.1848E-05,-12.45960808,0.01406560,0.2627E-05,281.47988892,15.00190926,0.2509E-05,0.55227423,-0.00013244,-0.1161E-04,0.00608628,-0.00013178,-0.1156E-04,0.00473334,0.00470976 ]
    i = i + 1
    eclipse_list(i,:) = [ 1999.,8.,11.,2451402.000,11.0,-3.0,3.0,64.2,63.7,0.07005581,0.54430360,-0.4061E-04,-0.8069E-05,0.50283802,-0.11849286,-0.1158E-03,0.1690E-05,15.32734108,-0.01203498,-0.3257E-05,343.69030762,15.00298214,0.1905E-05,0.54249305,0.00011683,-0.1168E-04,-0.00364585,0.00011625,-0.1162E-04,0.00461301,0.00459004 ]
    i = i + 1
    eclipse_list(i,:) = [ 2000.,2.,5.,2451580.000,13.0,-3.0,2.0,64.2,63.9,0.32956558,0.50568569,-0.2072E-04,-0.6043E-05,-1.18088663,0.10520779,0.1530E-03,-0.1170E-05,-16.02865028,0.01232359,0.3612E-05,11.51423740,15.00039768,0.2821E-05,0.56820714,-0.00009628,-0.1034E-04,0.02193982,-0.00009580,-0.1029E-04,0.00474337,0.00471975 ]
    i = i + 1
    eclipse_list(i,:) = [ 2000.,7.,1.,2451727.250,20.0,-2.0,2.0,64.2,64.0,0.28074276,0.58385366,-0.1043E-04,-0.9957E-05,-1.27723312,0.01061091,-0.9229E-04,-0.2900E-06,23.04426193,-0.00308376,-0.5589E-05,119.01594543,14.99953175,0.9498E-06,0.53008103,-0.00001814,-0.1281E-04,-0.01599601,-0.00001805,-0.1275E-04,0.00459821,0.00457531 ]
    i = i + 1
    eclipse_list(i,:) = [ 2000.,7.,31.,2451756.500,2.0,-2.0,2.0,64.2,64.0,0.06577264,0.57378113,-0.3147E-04,-0.9690E-05,1.22254825,-0.09595570,-0.1900E-03,0.1545E-05,18.21852875,-0.01008576,-0.3844E-05,208.41378784,15.00197792,0.2064E-05,0.53161967,0.00005064,-0.1271E-04,-0.01446506,0.00005038,-0.1264E-04,0.00460562,0.00458268 ]
    i = i + 1
    eclipse_list(i,:) = [ 2000.,12.,25.,2451904.250,18.0,-3.0,2.0,64.2,64.1,0.27835506,0.50631285,-0.8421E-05,-0.5756E-05,1.12071812,-0.03355852,0.1120E-03,0.4911E-06,-23.37237930,0.00144918,0.6511E-05,89.91230774,14.99619579,0.3845E-06,0.57311630,0.00005135,-0.1005E-04,0.02682453,0.00005110,-0.9999E-05,0.00475540,0.00473171 ]
    i = i + 1
    eclipse_list(i,:) = [ 2001.,6.,21.,2452082.000,12.0,-3.0,3.0,64.2,64.2,0.01035507,0.56538600,0.2920E-04,-0.8862E-05,-0.57182974,0.05512570,-0.1339E-03,-0.9818E-06,23.43976593,-0.00018302,-0.5651E-05,359.56451416,14.99919319,0.2060E-06,0.53699583,-0.00009439,-0.1212E-04,-0.00911568,-0.00009392,-0.1206E-04,0.00460053,0.00457762 ]
    i = i + 1
    eclipse_list(i,:) = [ 2001.,12.,14.,2452258.250,21.0,-3.0,3.0,64.2,64.3,0.12366933,0.52894205,0.1518E-05,-0.6983E-05,0.39464471,-0.08130247,0.1455E-03,0.1207E-05,-23.25179863,-0.00205529,0.6342E-05,136.28274536,14.99641228,-0.9329E-06,0.55955315,0.00011852,-0.1110E-04,0.01332891,0.00011793,-0.1105E-04,0.00475089,0.00472723 ]
    i = i + 1
    eclipse_list(i,:) = [ 2002.,6.,10.,2452436.500,0.0,-4.0,3.0,64.2,64.4,0.09391917,0.52745426,0.4029E-04,-0.6937E-05,0.21901830,0.09319092,-0.1594E-03,-0.1364E-05,23.05546570,0.00279431,-0.5490E-05,180.13783264,14.99926090,-0.6142E-06,0.55168074,-0.00012314,-0.1079E-04,0.00549607,-0.00012253,-0.1074E-04,0.00460544,0.00458251 ]
    i = i + 1
    eclipse_list(i,:) = [ 2002.,12.,4.,2452612.750,8.0,-4.0,3.0,64.2,64.4,0.18614703,0.55325019,0.1548E-04,-0.8705E-05,-0.35446295,-0.13091308,0.1807E-03,0.2219E-05,-22.22649956,-0.00529656,0.5781E-05,302.48504639,14.99727154,-0.1955E-05,0.54418993,0.00008340,-0.1248E-04,-0.00195779,0.00008299,-0.1242E-04,0.00474379,0.00472017 ]
    i = i + 1
    eclipse_list(i,:) = [ 2003.,5.,31.,2452790.750,4.0,-3.0,3.0,64.2,64.7,-0.32385439,0.49390221,0.5450E-04,-0.5600E-05,0.94521183,0.12632880,-0.1691E-03,-0.1600E-05,21.84598541,0.00574630,-0.5100E-05,240.61499023,14.99978828,-0.1000E-05,0.56440508,-0.00005300,-0.9800E-05,0.01815700,-0.00005270,-0.9800E-05,0.00461260,0.00458960 ]
    i = i + 1
    eclipse_list(i,:) = [ 2003.,11.,23.,2452967.500,23.0,-3.0,2.0,64.2,64.8,-0.19795258,0.55689740,0.5705E-04,-0.9431E-05,-0.94790143,-0.17392465,0.1990E-03,0.3137E-05,-20.40538597,-0.00817705,0.4966E-05,168.39866638,14.99853039,-0.2544E-05,0.53735191,-0.00000323,-0.1312E-04,-0.00876173,-0.00000322,-0.1305E-04,0.00473484,0.00471126 ]
    i = i + 1
    eclipse_list(i,:) = [ 2004.,4.,19.,2453115.000,14.0,-3.0,2.0,65.0,65.0,0.69319844,0.46464908,-0.1011E-04,-0.5494E-05,-0.92219394,0.22851951,-0.2593E-04,-0.2882E-05,11.41963577,0.01377343,-0.2658E-05,30.25052071,15.00335503,-0.1480E-05,0.56265110,0.00007811,-0.1023E-04,0.01641171,0.00007773,-0.1018E-04,0.00465501,0.00463182 ]
    i = i + 1
    eclipse_list(i,:) = [ 2004.,10.,14.,2453292.750,3.0,-3.0,3.0,65.2,65.2,0.47698101,0.48947230,0.1588E-04,-0.7005E-05,0.91833633,-0.25634456,-0.6181E-05,0.3839E-05,-8.23845959,-0.01485915,0.1689E-05,228.50190735,15.00359821,-0.1873E-05,0.54827005,-0.00011405,-0.1167E-04,0.00210218,-0.00011349,-0.1161E-04,0.00468820,0.00466486 ]
    i = i + 1
    eclipse_list(i,:) = [ 2005.,4.,8.,2453469.250,21.0,-4.0,3.0,65.4,65.4,0.35023588,0.48577625,-0.1011E-04,-0.6846E-05,-0.20743079,0.25757581,-0.3666E-04,-0.3827E-05,7.48629189,0.01489852,-0.1748E-05,134.57572937,15.00404167,-0.8759E-06,0.54857194,0.00010693,-0.1152E-04,0.00240261,0.00010639,-0.1146E-04,0.00466870,0.00464545 ]
    i = i + 1
    eclipse_list(i,:) = [ 2005.,10.,3.,2453647.000,11.0,-4.0,3.0,65.6,65.6,0.36604807,0.45547906,0.9929E-05,-0.5440E-05,0.17581950,-0.25080654,0.7474E-05,0.3169E-05,-4.09223843,-0.01549049,0.6973E-06,347.75479126,15.00434017,-0.1129E-05,0.56268519,-0.00010220,-0.1035E-04,0.01644556,-0.00010169,-0.1029E-04,0.00467412,0.00465084 ]
    i = i + 1
    eclipse_list(i,:) = [ 2006.,3.,29.,2453824.000,10.0,-3.0,3.0,65.8,65.8,-0.28991744,0.50608802,0.1814E-04,-0.8282E-05,0.27903894,0.27899188,-0.3863E-04,-0.4769E-05,3.39884138,0.01555752,-0.8382E-06,328.79370117,15.00436878,-0.1321E-06,0.53702563,0.00006441,-0.1275E-04,-0.00908622,0.00006409,-0.1268E-04,0.00468267,0.00465935 ]
    i = i + 1
    eclipse_list(i,:) = [ 2006.,9.,22.,2454001.000,12.0,-4.0,3.0,66.0,66.0,-0.06064556,0.44051006,0.2814E-05,-0.4828E-05,-0.43147260,-0.24612810,0.1519E-04,0.2867E-05,0.26179335,-0.01562140,-0.2616E-06,1.81183231,15.00473213,-0.3423E-06,0.56944436,-0.00000793,-0.9749E-05,0.02317113,-0.00000789,-0.9700E-05,0.00465938,0.00463618 ]
    i = i + 1
    eclipse_list(i,:) = [ 2007.,3.,19.,2454178.500,3.0,-3.0,2.0,66.2,66.2,-0.28987253,0.50867426,0.2184E-04,-0.8458E-05,1.06556320,0.28123668,-0.5353E-04,-0.4898E-05,-0.74650520,0.01579541,0.8996E-07,223.00160217,15.00427628,0.6587E-06,0.53624809,-0.00005461,-0.1293E-04,-0.00985993,-0.00005434,-0.1286E-04,0.00469673,0.00467334 ]
    i = i + 1
    eclipse_list(i,:) = [ 2007.,9.,11.,2454355.000,13.0,-3.0,2.0,66.4,66.4,-0.33108848,0.45388386,-0.9881E-05,-0.5347E-05,-1.10171282,-0.24807695,0.2776E-04,0.3078E-05,4.58239222,-0.01524352,-0.1191E-05,15.81616974,15.00475311,0.4169E-06,0.56241173,0.00008297,-0.1019E-04,0.01617353,0.00008256,-0.1014E-04,0.00464494,0.00462181 ]
    i = i + 1
    eclipse_list(i,:) = [ 2008.,2.,7.,2454503.750,4.0,-3.0,3.0,66.6,66.6,0.41246781,0.49431282,-0.4200E-04,-0.6605E-05,-0.86424482,0.21449882,0.1395E-03,-0.3021E-05,-15.51263428,0.01237704,0.3479E-05,236.47483826,15.00059891,0.2725E-05,0.55744088,-0.00012466,-0.1122E-04,0.01122719,-0.00012404,-0.1116E-04,0.00474181,0.00471819 ]
    i = i + 1
    eclipse_list(i,:) = [ 2008.,8.,1.,2454680.000,10.0,-2.0,3.0,66.8,66.8,0.10179438,0.52857709,-0.6343E-04,-0.8129E-05,0.85061806,-0.20252265,-0.1512E-03,0.3309E-05,17.86754036,-0.01012085,-0.3848E-05,328.42254639,15.00201130,0.1962E-05,0.53825587,0.00011109,-0.1203E-04,-0.00786193,0.00011054,-0.1197E-04,0.00460662,0.00458368 ]
    i = i + 1
    eclipse_list(i,:) = [ 2009.,1.,26.,2454857.750,8.0,-4.0,4.0,67.0,67.0,0.09962367,0.47819689,-0.3525E-04,-0.5442E-05,-0.26380020,0.17630750,0.1255E-03,-0.2159E-05,-18.64789963,0.01014513,0.4588E-05,296.85946655,14.99901676,0.2751E-05,0.57193238,-0.00006966,-0.1011E-04,0.02564650,-0.00006932,-0.1005E-04,0.00474988,0.00472623 ]
    i = i + 1
    eclipse_list(i,:) = [ 2009.,7.,22.,2455034.500,3.0,-4.0,3.0,67.2,67.2,0.23998830,0.55639547,-0.5762E-04,-0.9430E-05,-0.00328382,-0.17745794,-0.1344E-03,0.3176E-05,20.26424217,-0.00787358,-0.4559E-05,223.38414001,15.00100231,0.1928E-05,0.53045028,0.00000631,-0.1281E-04,-0.01562861,0.00000628,-0.1274E-04,0.00460142,0.00457850 ]
    i = i + 1
    eclipse_list(i,:) = [ 2010.,1.,15.,2455211.750,7.0,-3.0,4.0,67.4,67.4,-0.17324369,0.48452041,-0.3710E-04,-0.5374E-05,0.36640391,0.14049202,0.1170E-03,-0.1704E-05,-21.12929916,0.00730746,0.5566E-05,282.66671753,14.99759007,0.2250E-05,0.57469970,0.00003723,-0.9929E-05,0.02840006,0.00003704,-0.9880E-05,0.00475458,0.00473090 ]
    i = i + 1
    eclipse_list(i,:) = [ 2010.,7.,11.,2455389.250,20.0,-3.0,2.0,67.6,67.6,0.07409976,0.55725145,-0.2761E-04,-0.8967E-05,-0.71703005,-0.13665789,-0.1121E-03,0.2357E-05,22.03570175,-0.00534120,-0.5182E-05,118.60980225,15.00006771,0.1570E-05,0.53444636,-0.00009080,-0.1236E-04,-0.01165245,-0.00009035,-0.1230E-04,0.00459887,0.00457597 ]
    i = i + 1
    eclipse_list(i,:) = [ 2011.,1.,4.,2455565.750,9.0,-3.0,3.0,67.8,67.8,-0.14059769,0.51627612,-0.4174E-04,-0.6514E-05,1.05582511,0.10513873,0.1064E-03,-0.1449E-05,-22.74122429,0.00406206,0.6256E-05,313.80661011,14.99663258,0.1312E-05,0.56361586,0.00010976,-0.1079E-04,0.01737137,0.00010922,-0.1074E-04,0.00475578,0.00473209 ]
    i = i + 1
    eclipse_list(i,:) = [ 2011.,6.,1.,2455714.500,21.0,-2.0,3.0,68.0,68.0,-0.20262891,0.52627444,0.2975E-04,-0.6643E-05,1.20555782,0.02219590,-0.1931E-03,-0.2178E-06,22.09197044,0.00551956,-0.5104E-05,135.53350830,14.99974346,-0.1319E-05,0.55623817,-0.00010256,-0.1046E-04,0.01003081,-0.00010205,-0.1041E-04,0.00461123,0.00458827 ]
    i = i + 1
    eclipse_list(i,:) = [ 2011.,7.,1.,2455743.750,9.0,-2.0,1.0,68.0,68.0,-0.05926928,0.53551126,-0.3978E-05,-0.7370E-05,-1.50203753,-0.08818307,-0.8331E-04,0.1354E-05,23.12097740,-0.00254033,-0.5636E-05,314.04870605,14.99935722,0.9142E-06,0.54779315,-0.00012719,-0.1108E-04,0.00162788,-0.00012656,-0.1102E-04,0.00459904,0.00457613 ]
    i = i + 1
    eclipse_list(i,:) = [ 2011.,11.,25.,2455890.750,6.0,-2.0,3.0,68.2,68.2,-0.31001288,0.57290632,0.2052E-04,-0.9391E-05,-1.02766418,-0.05788106,0.2105E-03,0.8317E-06,-20.67681122,-0.00810005,0.5032E-05,273.29544067,14.99839687,-0.2637E-05,0.54062742,0.00008218,-0.1277E-04,-0.00550252,0.00008177,-0.1271E-04,0.00473611,0.00471252 ]
    i = i + 1
    eclipse_list(i,:) = [ 2012.,5.,20.,2456068.500,0.0,-4.0,3.0,68.4,68.4,-0.00223737,0.50318277,0.1833E-04,-0.5658E-05,0.48552877,0.05605376,-0.1411E-03,-0.5633E-06,20.22055626,0.00827127,-0.4662E-05,180.85285950,15.00057697,-0.1794E-05,0.56651109,-0.00003124,-0.9729E-05,0.02025256,-0.00003108,-0.9681E-05,0.00462052,0.00459751 ]
    i = i + 1
    eclipse_list(i,:) = [ 2012.,11.,13.,2456245.500,22.0,-3.0,3.0,68.6,68.6,-0.18370187,0.57434052,0.3505E-04,-0.9731E-05,-0.34676909,-0.09407004,0.1484E-03,0.1507E-05,-18.24755478,-0.01060200,0.4254E-05,153.89453125,14.99986744,-0.2871E-05,0.53746939,-0.00002289,-0.1303E-04,-0.00864481,-0.00002277,-0.1297E-04,0.00472560,0.00470207 ]
    i = i + 1
    eclipse_list(i,:) = [ 2013.,5.,10.,2456422.500,0.0,-3.0,4.0,68.8,68.8,-0.17514551,0.50528872,0.1445E-04,-0.5908E-05,-0.30429307,0.08888988,-0.9584E-04,-0.9696E-06,17.60548210,0.01070052,-0.4027E-05,180.90058899,15.00166130,-0.1998E-05,0.56369555,0.00007880,-0.1001E-04,0.01745101,0.00007840,-0.9955E-05,0.00463136,0.00460829 ]
    i = i + 1
    eclipse_list(i,:) = [ 2013.,11.,3.,2456600.000,13.0,-3.0,3.0,69.0,69.0,0.18323015,0.54694951,0.2823E-04,-0.8258E-05,0.29471111,-0.12007563,0.7904E-04,0.1729E-05,-15.20966434,-0.01263576,0.3402E-05,19.10831261,15.00132561,-0.2794E-05,0.54630607,-0.00011210,-0.1204E-04,0.00014792,-0.00011154,-0.1198E-04,0.00471376,0.00469029 ]
    i = i + 1
    eclipse_list(i,:) = [ 2014.,4.,29.,2456776.750,6.0,-3.0,3.0,69.2,69.2,0.18519846,0.52826846,-0.4991E-05,-0.7244E-05,-0.98351729,0.12211305,-0.4732E-04,-0.1619E-05,14.44979954,0.01265839,-0.3327E-05,270.65426636,15.00275612,-0.1836E-05,0.55056149,0.00011865,-0.1115E-04,0.00438229,0.00011806,-0.1110E-04,0.00464335,0.00462023 ]
    i = i + 1
    eclipse_list(i,:) = [ 2014.,10.,23.,2456954.500,22.0,-3.0,2.0,69.4,69.4,0.40257463,0.50848579,0.1603E-04,-0.6393E-05,1.02156377,-0.13553271,0.1798E-04,0.1623E-05,-11.61852455,-0.01419850,0.2503E-05,153.92974854,15.00265121,-0.2382E-05,0.56102097,-0.00011083,-0.1064E-04,0.01478960,-0.00011028,-0.1059E-04,0.00470062,0.00467721 ]
    i = i + 1
    eclipse_list(i,:) = [ 2015.,3.,20.,2457102.000,10.0,-3.0,2.0,69.6,69.6,-0.16825180,0.55374151,0.5801E-05,-0.9359E-05,0.93907130,0.17864808,-0.5413E-04,-0.2927E-05,-0.21265797,0.01603605,-0.2087E-07,328.10589600,15.00441456,0.5250E-06,0.53595513,0.00002650,-0.1295E-04,-0.01015145,0.00002636,-0.1289E-04,0.00469512,0.00467173 ]
    i = i + 1
    eclipse_list(i,:) = [ 2015.,9.,13.,2457278.750,7.0,-3.0,3.0,69.8,69.8,-0.29285854,0.48214203,0.1666E-05,-0.5405E-05,-1.06144142,-0.15181828,0.2264E-04,0.1638E-05,3.89024282,-0.01556351,-0.1047E-05,285.97631836,15.00485039,0.3000E-06,0.56828785,0.00002582,-0.9758E-05,0.02202037,0.00002569,-0.9709E-05,0.00464739,0.00462424 ]
    i = i + 1
    eclipse_list(i,:) = [ 2016.,3.,9.,2457456.500,2.0,-3.0,3.0,70.0,70.0,-0.06248009,0.55027688,0.4649E-05,-0.9054E-05,0.25385040,0.17212328,0.1713E-04,-0.2749E-05,-4.37971497,0.01588556,0.7605E-06,207.37260437,15.00397110,0.1244E-05,0.53889132,-0.00007038,-0.1275E-04,-0.00722991,-0.00007003,-0.1269E-04,0.00470875,0.00468530 ]
    i = i + 1
    eclipse_list(i,:) = [ 2016.,9.,1.,2457633.000,9.0,-3.0,4.0,70.2,70.2,-0.16135593,0.50406349,-0.2136E-04,-0.6303E-05,-0.29966924,-0.14815211,-0.2585E-04,0.1781E-05,8.06329823,-0.01480240,-0.1811E-05,315.03265381,15.00454521,0.9563E-06,0.55795211,0.00011155,-0.1047E-04,0.01173616,0.00011099,-0.1042E-04,0.00463402,0.00461095 ]
    i = i + 1
    eclipse_list(i,:) = [ 2017.,2.,26.,2457811.000,15.0,-3.0,3.0,70.4,70.4,0.17598145,0.52535623,-0.6196E-05,-0.7411E-05,-0.42554048,0.15325409,0.7922E-04,-0.2075E-05,-8.49164200,0.01526081,0.1616E-05,41.80067444,15.00308609,0.1926E-05,0.55249369,-0.00012569,-0.1152E-04,0.00630473,-0.00012506,-0.1146E-04,0.00472198,0.00469847 ]
    i = i + 1
    eclipse_list(i,:) = [ 2017.,8.,21.,2457987.250,18.0,-3.0,4.0,70.6,70.6,-0.12952891,0.54064256,-0.2931E-04,-0.8087E-05,0.48540425,-0.14163995,-0.9049E-04,0.2051E-05,11.86695576,-0.01362157,-0.2494E-05,89.24772644,15.00393772,0.1507E-05,0.54211676,0.00012407,-0.1178E-04,-0.00402028,0.00012345,-0.1172E-04,0.00462230,0.00459928 ]
    i = i + 1
    eclipse_list(i,:) = [ 2018.,2.,15.,2458165.250,21.0,-3.0,2.0,70.8,70.8,0.36366239,0.49905220,-0.2121E-04,-0.5920E-05,-1.15753388,0.12833364,0.1268E-03,-0.1440E-05,-12.46403790,0.01408048,0.2579E-05,131.48355103,15.00182343,0.2494E-05,0.56828094,-0.00009227,-0.1028E-04,0.02201329,-0.00009181,-0.1022E-04,0.00473410,0.00471052 ]
    i = i + 1
    eclipse_list(i,:) = [ 2018.,7.,13.,2458312.750,3.0,-2.0,2.0,71.0,71.0,-0.09923430,0.58281463,-0.1318E-05,-0.9927E-05,-1.35077310,-0.03329333,-0.7703E-04,0.4696E-06,21.84531212,-0.00593743,-0.5210E-05,223.57409668,15.00024128,0.1668E-05,0.53019166,-0.00001183,-0.1280E-04,-0.01588594,-0.00001177,-0.1274E-04,0.00459887,0.00457596 ]
    i = i + 1
    eclipse_list(i,:) = [ 2018.,8.,11.,2458342.000,10.0,-2.0,2.0,71.0,71.0,0.36755145,0.56849587,-0.4773E-04,-0.9621E-05,1.09391034,-0.12629345,-0.1598E-03,0.2063E-05,15.21672916,-0.01207630,-0.3122E-05,328.69937134,15.00307846,0.1878E-05,0.53172153,0.00003379,-0.1275E-04,-0.01436375,0.00003362,-0.1269E-04,0.00461273,0.00458976 ]
    i = i + 1
    eclipse_list(i,:) = [ 2019.,1.,6.,2458489.500,2.0,-3.0,2.0,71.2,71.2,0.12841544,0.50823843,-0.1620E-04,-0.5820E-05,1.14402676,0.00842363,0.1036E-03,-0.1089E-08,-22.54492188,0.00484782,0.6184E-05,208.61903381,14.99673653,0.1554E-05,0.57272607,0.00005746,-0.1008E-04,0.02643623,0.00005717,-0.1003E-04,0.00475629,0.00473260 ]
    i = i + 1
    eclipse_list(i,:) = [ 2019.,7.,2.,2458667.250,19.0,-3.0,3.0,71.4,71.4,-0.21557993,0.56620866,0.2740E-04,-0.8788E-05,-0.65071219,0.01063996,-0.1272E-03,-0.2684E-06,23.01294899,-0.00318702,-0.5496E-05,103.98381042,14.99950504,0.1079E-05,0.53765464,-0.00008982,-0.1203E-04,-0.00846013,-0.00008937,-0.1198E-04,0.00459846,0.00457556 ]
    i = i + 1
    eclipse_list(i,:) = [ 2019.,12.,26.,2458843.750,5.0,-3.0,4.0,71.6,71.6,-0.14036173,0.53561032,-0.1459E-05,-0.7150E-05,0.42407426,-0.03665511,0.1458E-03,0.6043E-06,-23.37347412,0.00140653,0.6424E-05,254.94107056,14.99627018,0.3278E-06,0.55891114,0.00012839,-0.1118E-04,0.01269010,0.00012775,-0.1112E-04,0.00475485,0.00473117 ]
    i = i + 1
    eclipse_list(i,:) = [ 2020.,6.,21.,2459021.750,7.0,-4.0,3.0,71.8,71.8,0.15431480,0.53115451,0.2586E-04,-0.6921E-05,0.13640723,0.05138710,-0.1610E-03,-0.7911E-06,23.43567085,-0.00023279,-0.5590E-05,284.53991699,14.99910927,0.2419E-06,0.55234158,-0.00012230,-0.1071E-04,0.00615367,-0.00012169,-0.1066E-04,0.00460097,0.00457806 ]
    i = i + 1
    eclipse_list(i,:) = [ 2020.,12.,14.,2459198.250,16.0,-3.0,3.0,72.0,72.0,-0.18176793,0.56335676,0.2165E-04,-0.8956E-05,-0.26964971,-0.08581217,0.1884E-03,0.1504E-05,-23.25776482,-0.00198600,0.6238E-05,61.27030945,14.99649620,-0.8593E-06,0.54388565,0.00009702,-0.1256E-04,-0.00226052,0.00009654,-0.1249E-04,0.00475029,0.00472663 ]
    i = i + 1
    eclipse_list(i,:) = [ 2021.,6.,10.,2459376.000,11.0,-3.0,3.0,72.2,72.2,-0.01865130,0.50122893,0.3420E-04,-0.5706E-05,0.92611074,0.08877646,-0.1797E-03,-0.1132E-05,23.04228592,0.00284114,-0.5435E-05,345.13113403,14.99919891,-0.6637E-06,0.56440377,-0.00005512,-0.9795E-05,0.01815576,-0.00005484,-0.9747E-05,0.00460603,0.00458309 ]
    i = i + 1
    eclipse_list(i,:) = [ 2021.,12.,4.,2459552.750,8.0,-3.0,2.0,72.4,72.4,0.02526544,0.56830269,0.3910E-04,-0.9655E-05,-0.98367131,-0.13151418,0.2213E-03,0.2409E-05,-22.27472115,-0.00517795,0.5701E-05,302.45623779,14.99727917,-0.1860E-05,0.53782868,-0.00001604,-0.1313E-04,-0.00828731,-0.00001596,-0.1307E-04,0.00474351,0.00471989 ]
    i = i + 1
    eclipse_list(i,:) = [ 2022.,4.,30.,2459700.250,21.0,-3.0,2.0,72.6,72.6,0.61813062,0.47531474,-0.1515E-05,-0.5684E-05,-1.02807009,0.20964050,-0.4317E-04,-0.2682E-05,14.97104359,0.01216706,-0.3476E-05,135.70941162,15.00246811,-0.1861E-05,0.56109691,0.00008474,-0.1027E-04,0.01486530,0.00008432,-0.1022E-04,0.00464210,0.00461898 ]
    i = i + 1
    eclipse_list(i,:) = [ 2022.,10.,25.,2459878.000,11.0,-3.0,3.0,72.8,72.8,0.45485225,0.49554944,0.2763E-04,-0.7023E-05,0.96874976,-0.23958749,0.1676E-04,0.3562E-05,-12.17347527,-0.01374551,0.2660E-05,348.98565674,15.00242996,-0.2461E-05,0.54990345,-0.00011520,-0.1162E-04,0.00372742,-0.00011463,-0.1156E-04,0.00470198,0.00467856 ]
    i = i + 1
    eclipse_list(i,:) = [ 2023.,4.,20.,2460054.750,4.0,-3.0,4.0,73.0,73.0,0.02690837,0.49501815,0.1358E-04,-0.7053E-05,-0.42734402,0.24419928,-0.4941E-04,-0.3671E-05,11.41178703,0.01374089,-0.2582E-05,240.24565125,15.00341797,-0.1467E-05,0.54682785,0.00012161,-0.1157E-04,0.00066722,0.00012101,-0.1151E-04,0.00465503,0.00463184 ]
    i = i + 1
    eclipse_list(i,:) = [ 2023.,10.,14.,2460232.250,18.0,-3.0,3.0,73.3,73.3,0.16971114,0.45855319,0.2780E-04,-0.5428E-05,0.33483642,-0.24136706,0.2405E-04,0.3030E-05,-8.24418926,-0.01488811,0.1615E-05,93.50385284,15.00352955,-0.1849E-05,0.56433600,-0.00008907,-0.1031E-04,0.01808814,-0.00008863,-0.1026E-04,0.00468827,0.00466492 ]
    i = i + 1
    eclipse_list(i,:) = [ 2024.,4.,8.,2460409.250,18.0,-3.0,3.0,73.5,73.5,-0.31818941,0.51171160,0.3261E-04,-0.8416E-05,0.21979441,0.27095896,-0.5944E-04,-0.4659E-05,7.58619833,0.01484435,-0.1703E-05,89.59259033,15.00408363,-0.8164E-06,0.53583723,0.00006179,-0.1276E-04,-0.01026868,0.00006148,-0.1270E-04,0.00466833,0.00464509 ]
    i = i + 1
    eclipse_list(i,:) = [ 2024.,10.,2.,2460586.250,19.0,-4.0,3.0,73.7,73.7,-0.06798750,0.44161698,0.1361E-04,-0.4830E-05,-0.36319631,-0.24356307,0.3393E-04,0.2835E-05,-3.98724580,-0.01551079,0.6083E-06,107.73171997,15.00433159,-0.1104E-05,0.57037348,-0.00000017,-0.9765E-05,0.02409558,-0.00000017,-0.9717E-05,0.00467347,0.00465020 ]
    i = i + 1
    eclipse_list(i,:) = [ 2025.,3.,29.,2460764.000,11.0,-3.0,2.0,73.9,73.9,-0.40281191,0.50941223,0.4146E-04,-0.8439E-05,0.96572608,0.27883485,-0.7227E-04,-0.4841E-05,3.56601572,0.01553849,-0.8087E-06,343.83151245,15.00436687,-0.7590E-07,0.53578967,-0.00005328,-0.1286E-04,-0.01031602,-0.00005302,-0.1280E-04,0.00468239,0.00465907 ]
    i = i + 1
    eclipse_list(i,:) = [ 2025.,9.,21.,2460940.250,20.0,-3.0,2.0,74.1,74.1,-0.39001107,0.45315927,0.3144E-05,-0.5387E-05,-1.00185990,-0.25216335,0.4563E-04,0.3155E-05,0.36471686,-0.01559954,-0.3545E-06,121.78099823,15.00477219,-0.3196E-06,0.56251699,0.00009090,-0.1026E-04,0.01627823,0.00009045,-0.1021E-04,0.00465839,0.00463519 ]
    i = i + 1
    eclipse_list(i,:) = [ 2026.,2.,17.,2461089.000,12.0,-3.0,3.0,74.3,74.3,0.32201946,0.48272234,-0.3147E-04,-0.6367E-05,-0.92694736,0.23553929,0.1168E-03,-0.3267E-05,-11.87930012,0.01404900,0.2440E-05,356.51254272,15.00198364,0.2377E-05,0.55774337,-0.00011808,-0.1112E-04,0.01152824,-0.00011749,-0.1107E-04,0.00473212,0.00470856 ]
    i = i + 1
    eclipse_list(i,:) = [ 2026.,8.,12.,2461265.250,18.0,-3.0,3.0,74.5,74.5,0.47558114,0.51892501,-0.7730E-04,-0.8041E-05,0.77115941,-0.23016810,-0.1246E-03,0.3768E-05,14.79667187,-0.01206486,-0.3098E-05,88.74527740,15.00308990,0.1778E-05,0.53797895,0.00009394,-0.1212E-04,-0.00813748,0.00009347,-0.1206E-04,0.00461413,0.00459115 ]
    i = i + 1
    eclipse_list(i,:) = [ 2027.,2.,6.,2461443.250,16.0,-4.0,4.0,74.7,74.7,0.11174255,0.46649510,-0.3371E-04,-0.5268E-05,-0.27327067,0.20318554,0.1025E-03,-0.2456E-05,-15.54794216,0.01238295,0.3584E-05,56.48994446,15.00050926,0.2747E-05,0.57195204,-0.00006532,-0.1006E-04,0.02566613,-0.00006500,-0.1001E-04,0.00474266,0.00471904 ]
    i = i + 1
    eclipse_list(i,:) = [ 2027.,8.,2.,2461620.000,10.0,-3.0,3.0,74.9,74.9,-0.01969843,0.54471231,-0.4463E-04,-0.9221E-05,0.16003945,-0.21115834,-0.1217E-03,0.3757E-05,17.76247406,-0.01018110,-0.3875E-05,328.41891479,15.00209618,0.2041E-05,0.53062016,0.00001381,-0.1283E-04,-0.01545961,0.00001374,-0.1277E-04,0.00460646,0.00458352 ]
    i = i + 1
    eclipse_list(i,:) = [ 2028.,1.,26.,2461797.250,15.0,-3.0,4.0,75.1,75.1,-0.20521516,0.47425699,-0.3901E-04,-0.5264E-05,0.34029821,0.17385875,0.9684E-04,-0.2087E-05,-18.72824669,0.01007411,0.4695E-05,41.88721848,14.99896240,0.2710E-05,0.57414103,0.00004203,-0.9944E-05,0.02784414,0.00004182,-0.9894E-05,0.00475017,0.00472651 ]
    i = i + 1
    eclipse_list(i,:) = [ 2028.,7.,22.,2461974.500,3.0,-3.0,3.0,75.3,75.3,-0.15433459,0.54498911,-0.2147E-04,-0.8679E-05,-0.58644259,-0.17460850,-0.1021E-03,0.2954E-05,20.18231201,-0.00797363,-0.4621E-05,223.37434387,15.00101662,0.1991E-05,0.53526050,-0.00008587,-0.1230E-04,-0.01084236,-0.00008544,-0.1224E-04,0.00460163,0.00457871 ]
    i = i + 1
    eclipse_list(i,:) = [ 2029.,1.,14.,2462151.250,17.0,-2.0,3.0,75.5,75.5,-0.40737015,0.50815260,-0.3929E-04,-0.6453E-05,0.98107040,0.14552829,0.9214E-04,-0.1988E-05,-21.16300774,0.00724073,0.5641E-05,72.68831635,14.99763107,0.2225E-05,0.56268984,0.00011889,-0.1085E-04,0.01645001,0.00011830,-0.1080E-04,0.00475413,0.00473045 ]
    i = i + 1
    eclipse_list(i,:) = [ 2029.,6.,12.,2462299.750,4.0,-2.0,2.0,75.7,75.7,-0.01072016,0.52476060,0.1041E-04,-0.6536E-05,1.29541421,-0.01763651,-0.2057E-03,0.2900E-06,23.15931892,0.00259083,-0.5408E-05,240.03134155,14.99919891,-0.5841E-06,0.55668622,-0.00010270,-0.1038E-04,0.01047663,-0.00010219,-0.1033E-04,0.00460491,0.00458197 ]
    i = i + 1
    eclipse_list(i,:) = [ 2029.,7.,11.,2462329.250,16.0,-2.0,1.0,75.8,75.8,-0.13726678,0.52526319,-0.9467E-05,-0.7030E-05,-1.42716122,-0.12804167,-0.7690E-04,0.1876E-05,22.00244713,-0.00542267,-0.5259E-05,58.59796906,15.00000668,0.1633E-05,0.54877949,-0.00012690,-0.1100E-04,0.00260931,-0.00012627,-0.1094E-04,0.00459944,0.00457653 ]
    i = i + 1
    eclipse_list(i,:) = [ 2029.,12.,5.,2462476.250,15.0,-2.0,3.0,76.0,76.0,-0.06375255,0.57663524,-0.2660E-05,-0.9501E-05,-1.05967152,-0.01401654,0.2295E-03,0.1004E-06,-22.44545174,-0.00505369,0.5749E-05,47.30568695,14.99717236,-0.1912E-05,0.54066575,0.00006986,-0.1283E-04,-0.00546439,0.00006952,-0.1277E-04,0.00474464,0.00472101 ]
    i = i + 1
    eclipse_list(i,:) = [ 2030.,6.,1.,2462653.750,6.0,-3.0,4.0,76.2,76.2,-0.26931521,0.50563711,0.1815E-04,-0.5681E-05,0.55198264,0.02101492,-0.1585E-03,-0.1570E-06,22.06129646,0.00558135,-0.5180E-05,270.53576660,14.99970150,-0.1362E-05,0.56617349,-0.00001297,-0.9705E-05,0.01991664,-0.00001290,-0.9657E-05,0.00461206,0.00458909 ]
    i = i + 1
    eclipse_list(i,:) = [ 2030.,11.,25.,2462830.750,7.0,-3.0,3.0,76.4,76.4,0.04423760,0.57877976,0.1766E-04,-0.9774E-05,-0.39266980,-0.05518899,0.1743E-03,0.8359E-06,-20.76099968,-0.00798900,0.5154E-05,288.27084351,14.99836063,-0.2586E-05,0.53823727,-0.00003787,-0.1303E-04,-0.00788074,-0.00003768,-0.1297E-04,0.00473614,0.00471255 ]
    i = i + 1
    eclipse_list(i,:) = [ 2031.,5.,21.,2463007.750,7.0,-3.0,4.0,76.6,76.6,-0.11470663,0.51123929,0.7228E-05,-0.6025E-05,-0.21123368,0.05793294,-0.1182E-03,-0.6058E-06,20.15914917,0.00833890,-0.4697E-05,285.84768677,15.00061989,-0.1871E-05,0.56242877,0.00008064,-0.1004E-04,0.01619054,0.00008024,-0.9987E-05,0.00462090,0.00459789 ]
    i = i + 1
    eclipse_list(i,:) = [ 2031.,11.,14.,2463185.500,21.0,-3.0,3.0,76.8,76.8,-0.01978426,0.55094385,0.3659E-04,-0.8228E-05,0.31495243,-0.08906509,0.1046E-03,0.1244E-05,-18.33680916,-0.01053444,0.4392E-05,138.89105225,14.99976349,-0.2869E-05,0.54779804,-0.00010683,-0.1199E-04,0.00163239,-0.00010630,-0.1193E-04,0.00472611,0.00470257 ]
    i = i + 1
    eclipse_list(i,:) = [ 2032.,5.,9.,2463362.000,13.0,-2.0,3.0,77.0,77.0,-0.07427618,0.53595459,0.5256E-05,-0.7426E-05,-0.96542984,0.09540578,-0.7024E-04,-0.1259E-05,17.59290504,0.01069434,-0.4083E-05,15.88661385,15.00173855,-0.2017E-05,0.54887635,0.00012718,-0.1121E-04,0.00270559,0.00012654,-0.1116E-04,0.00463105,0.00460799 ]
    i = i + 1
    eclipse_list(i,:) = [ 2032.,11.,3.,2463539.750,6.0,-3.0,2.0,77.2,77.2,0.44932297,0.51201904,0.1702E-04,-0.6382E-05,0.99080956,-0.11286831,0.4517E-04,0.1327E-05,-15.23991680,-0.01263290,0.3496E-05,274.11730957,15.00123024,-0.2807E-05,0.56262904,-0.00011268,-0.1060E-04,0.01638958,-0.00011212,-0.1054E-04,0.00471418,0.00469071 ]
    i = i + 1
    eclipse_list(i,:) = [ 2033.,3.,30.,2463687.250,18.0,-3.0,3.0,77.4,77.4,-0.31876293,0.55542439,0.2270E-04,-0.9424E-05,0.92470294,0.17566095,-0.8011E-04,-0.2889E-05,4.09368324,0.01571873,-0.8907E-06,88.92649841,15.00445461,-0.2408E-06,0.53496599,0.00002761,-0.1294E-04,-0.01113563,0.00002747,-0.1288E-04,0.00468079,0.00465748 ]
    i = i + 1
    eclipse_list(i,:) = [ 2033.,9.,23.,2463864.000,14.0,-3.0,3.0,77.6,77.6,-0.30991709,0.48154482,0.8664E-05,-0.5408E-05,-1.11703670,-0.15454406,0.4778E-04,0.1670E-05,-0.33981988,-0.01584517,-0.2248E-06,31.94151878,15.00480556,-0.4559E-06,0.56892198,0.00003177,-0.9790E-05,0.02265132,0.00003162,-0.9742E-05,0.00466078,0.00463757 ]
    i = i + 1
    eclipse_list(i,:) = [ 2034.,3.,20.,2464042.000,10.0,-3.0,3.0,77.8,77.8,-0.25951615,0.54816294,0.2337E-04,-0.8964E-05,0.22078310,0.17557891,-0.7998E-05,-0.2790E-05,-0.05512699,0.01604236,-0.1583E-06,328.13885498,15.00440121,0.5047E-06,0.53865451,-0.00006649,-0.1267E-04,-0.00746552,-0.00006616,-0.1261E-04,0.00469521,0.00467183 ]
    i = i + 1
    eclipse_list(i,:) = [ 2034.,9.,12.,2464218.250,16.0,-3.0,4.0,78.0,78.0,-0.28082183,0.50283420,-0.1063E-04,-0.6348E-05,-0.32437146,-0.15778460,-0.8466E-06,0.1917E-05,3.97190952,-0.01553381,-0.1011E-05,60.95006561,15.00490379,0.2652E-06,0.55782509,0.00011884,-0.1055E-04,0.01160974,0.00011825,-0.1050E-04,0.00464622,0.00462308 ]
    i = i + 1
    eclipse_list(i,:) = [ 2035.,3.,9.,2464396.500,23.0,-3.0,3.0,78.2,78.2,0.07955908,0.52057379,0.4932E-05,-0.7276E-05,-0.43279579,0.16309454,0.5322E-04,-0.2195E-05,-4.27334070,0.01591991,0.6503E-06,162.39717102,15.00390434,0.1224E-05,0.55264676,-0.00012195,-0.1142E-04,0.00645701,-0.00012134,-0.1137E-04,0.00470957,0.00468612 ]
    i = i + 1
    eclipse_list(i,:) = [ 2035.,9.,2.,2464572.500,2.0,-3.0,3.0,78.4,78.4,0.13437596,0.53777355,-0.3596E-04,-0.8123E-05,0.34897482,-0.15846506,-0.5955E-04,0.2323E-05,8.01771450,-0.01478295,-0.1713E-05,210.03163147,15.00464153,0.9310E-06,0.54194397,0.00011035,-0.1188E-04,-0.00419227,0.00010980,-0.1182E-04,0.00463289,0.00460982 ]
    i = i + 1
    eclipse_list(i,:) = [ 2036.,2.,27.,2464750.750,5.0,-3.0,2.0,78.6,78.6,0.44412693,0.49340099,-0.2008E-04,-0.5811E-05,-1.11428463,0.14454030,0.9968E-04,-0.1622E-05,-8.49969006,0.01528089,0.1568E-05,251.81062317,15.00299644,0.1913E-05,0.56821597,-0.00009059,-0.1021E-04,0.02194867,-0.00009013,-0.1016E-04,0.00472313,0.00469961 ]
    i = i + 1
    eclipse_list(i,:) = [ 2036.,7.,23.,2464898.000,11.0,-2.0,1.0,78.8,78.8,0.09012426,0.57882214,-0.1798E-04,-0.9837E-05,-1.44783211,-0.07336815,-0.5466E-04,0.1159E-05,19.89421082,-0.00853668,-0.4645E-05,343.36456299,15.00124264,0.2058E-05,0.53045923,-0.00003060,-0.1280E-04,-0.01561971,-0.00003044,-0.1273E-04,0.00460193,0.00457901 ]
    i = i + 1
    eclipse_list(i,:) = [ 2036.,8.,21.,2464927.250,17.0,-2.0,3.0,78.8,78.8,0.03654715,0.56328863,-0.2802E-04,-0.9556E-05,1.11035466,-0.14969717,-0.1354E-03,0.2466E-05,11.74119091,-0.01364626,-0.2357E-05,74.26192474,15.00402737,0.1488E-05,0.53193229,0.00004448,-0.1280E-04,-0.01415403,0.00004426,-0.1274E-04,0.00462183,0.00459881 ]
    i = i + 1
    eclipse_list(i,:) = [ 2037.,1.,16.,2465075.000,10.0,-3.0,2.0,79.0,79.0,-0.01334805,0.50710255,-0.2147E-04,-0.5836E-05,1.15150464,0.04756258,0.8750E-04,-0.4619E-06,-20.83011246,0.00796920,0.5515E-05,327.55361938,14.99783230,0.2418E-05,0.57210594,0.00006330,-0.1012E-04,0.02581923,0.00006299,-0.1007E-04,0.00475407,0.00473039 ]
    i = i + 1
    eclipse_list(i,:) = [ 2037.,7.,13.,2465252.500,3.0,-3.0,3.0,79.2,79.2,0.14160752,0.56359959,0.4717E-07,-0.8689E-05,-0.73372155,-0.03182168,-0.1131E-03,0.4051E-06,21.78242683,-0.00604630,-0.5105E-05,223.55377197,15.00022411,0.1687E-05,0.53840691,-0.00011015,-0.1197E-04,-0.00771159,-0.00010960,-0.1191E-04,0.00459938,0.00457648 ]
    i = i + 1
    eclipse_list(i,:) = [ 2038.,1.,5.,2465429.000,14.0,-4.0,3.0,79.4,79.4,0.10901411,0.53854102,-0.2483E-04,-0.7270E-05,0.41855288,0.00799331,0.1379E-03,-0.1441E-08,-22.55480576,0.00481185,0.6102E-05,28.64703369,14.99681282,0.1541E-05,0.55818731,0.00011648,-0.1126E-04,0.01196985,0.00011590,-0.1120E-04,0.00475543,0.00473174 ]
    i = i + 1
    eclipse_list(i,:) = [ 2038.,7.,2.,2465607.000,14.0,-4.0,3.0,79.6,79.6,0.23929454,0.53154081,0.9912E-05,-0.6866E-05,0.04415986,0.00971130,-0.1544E-03,-0.2305E-06,22.99406433,-0.00324049,-0.5430E-05,28.96660233,14.99942017,0.1076E-05,0.55315948,-0.00012261,-0.1063E-04,0.00696745,-0.00012200,-0.1058E-04,0.00459941,0.00457651 ]
    i = i + 1
    eclipse_list(i,:) = [ 2038.,12.,26.,2465783.500,1.0,-3.0,3.0,79.8,79.8,-0.02056916,0.56985623,-0.1691E-06,-0.9127E-05,-0.28739566,-0.03797106,0.1915E-03,0.7289E-06,-23.36258125,0.00148085,0.6306E-05,194.92645264,14.99637413,0.3991E-06,0.54353130,0.00008670,-0.1262E-04,-0.00261315,0.00008627,-0.1256E-04,0.00475380,0.00473012 ]
    i = i + 1
    eclipse_list(i,:) = [ 2039.,6.,21.,2465961.250,17.0,-3.0,3.0,80.0,80.0,-0.18950747,0.50586927,0.2760E-04,-0.5741E-05,0.81659442,0.04954139,-0.1800E-03,-0.6708E-06,23.43237877,-0.00017776,-0.5544E-05,74.54014587,14.99903870,0.2092E-06,0.56454962,-0.00003876,-0.9760E-05,0.01830092,-0.00003856,-0.9711E-05,0.00460181,0.00457889 ]
    i = i + 1
    eclipse_list(i,:) = [ 2039.,12.,15.,2466138.250,16.0,-2.0,3.0,80.2,80.2,-0.36587262,0.57692868,0.4728E-04,-0.9809E-05,-0.90212750,-0.08494870,0.2295E-03,0.1594E-05,-23.27405167,-0.00186174,0.6139E-05,61.23015976,14.99652195,-0.7869E-06,0.53823078,-0.00000043,-0.1314E-04,-0.00788728,-0.00000043,-0.1307E-04,0.00474994,0.00472628 ]

    doFill = .false.

    end subroutine fillEclipse

end module eclipse_module
